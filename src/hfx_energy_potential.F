!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate HFX energy and potential 
!> \par History
!>      11.2006 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE hfx_energy_potential 
  
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE bibliography,                    ONLY: cite_reference,&
                                             guidon2008
  USE cell_types,                      ONLY: cell_type,&
                                             pbc,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE cp_files,                        ONLY: get_unit_number
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE erf_fn,                          ONLY: erfc
  USE gamma,                           ONLY: init_md_ftable
  USE hfx_communication,               ONLY: distribute_ks_matrix,&
                                             get_full_density
  USE hfx_compression_methods,         ONLY: hfx_add_mult_cache_elements,&
                                             hfx_add_single_cache_element,&
                                             hfx_decompress_first_cache,&
                                             hfx_flush_last_cache,&
                                             hfx_get_mult_cache_elements,&
                                             hfx_get_single_cache_element,&
                                             hfx_reset_cache_and_container
  USE hfx_helpers,                     ONLY: FNV64,&
                                             get_occupation
  USE hfx_libint_interface,            ONLY: evaluate_eri,&
                                             set_eps_cutoff
  USE hfx_libint_wrapper,              ONLY: initialize_libint,&
                                             lib_int
  USE hfx_load_balance_methods,        ONLY: load_balance_energy,&
                                             update_load_balance_energy
  USE hfx_screening_methods,           ONLY: comp_screen_mat,&
                                             comp_screen_mat_per,&
                                             screen_farfield,&
                                             screen_farfield_per,&
                                             update_pmax_mat
  USE hfx_types,                       ONLY: &
       hfx_basis_info_type, hfx_basis_type, hfx_box_scheme, hfx_cache_type, &
       hfx_cell_type, hfx_container_type, hfx_general_type, &
       hfx_init_container, hfx_load_balance_type, hfx_memory_type, &
       hfx_potential_type, hfx_screening_type, hfx_type
  USE input_constants,                 ONLY: do_hfx_potential_coulomb,&
                                             do_hfx_potential_gaussian,&
                                             do_hfx_potential_id,&
                                             do_hfx_potential_mix_cl,&
                                             do_hfx_potential_mix_lg,&
                                             do_hfx_potential_short,&
                                             do_hfx_potential_truncated
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: m_walltime
  USE mathconstants
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum,&
                                             mp_sum_i8,&
                                             mp_sync
  USE orbital_pointers
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE sparse_matrix_types,             ONLY: first_block_node,&
                                             get_block_node,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type
  USE t_c_g0,                          ONLY: init
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC  integrate_four_center
  TYPE(lib_int),   ALLOCATABLE, DIMENSION(:), SAVE      :: lib
  INTEGER, SAVE                                         :: init_t_c_g0_lmax = -1

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_energy_potential'

!***
  
  CONTAINS

! *****************************************************************************
!> \brief computes four center integrals for a full basis set and updates the
!>      Kohn-Sham-Matrix and energy. Uses all 8 eri symmetries
!> \param qs_env  
!> \param ks_matri x
!> \param energ y
!> \param rho density matrix
!> \param para_en v
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2007 created [Manuel Guidon]
!>      08.2007 optimized load balance [Manuel Guidon] 
!>      09.2007 new parallelization [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE integrate_four_center(qs_env,ks_matrix,energy,rho,hfx_section,para_env,& 
                                   geometry_did_change,irep,distribute_fock_matrix,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: hfx_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL                                  :: geometry_did_change
    INTEGER                                  :: irep
    LOGICAL, INTENT(IN)                      :: distribute_fock_matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_four_center', &
      routineP = moduleN//':'//routineN

    INTEGER :: bin, buffer_left, buffer_size, buffer_start, cache_size, &
      current_counter, handle, handle_dist_ks, handle_getP, handle_main, i, &
      i_thread, iatom, ikind, iset, iw, jatom, jkind, jset, katom, kkind, &
      kset, latom, lkind, lset, ma, max_am, max_set, maxsgf, mb, mc, md, &
      my_bin_size, my_istart, my_jstart, my_kstart, my_lstart, n_threads, &
      natom, nbits, nco_max, ncoa, ncob, ncoc, ncod, nints, nseta, nsetb, &
      nsetc, nsetd, nsgfa_total, nsgfb_total, nsgfc_total, nsgfd_total, &
      nspins, p_index, pa, pa1, pa2, pb, pb1, pb2, pc, pc1, pc2, pd, pd1, &
      pd2, pgf_radius_box, radius1, radius2, sgfa, sgfb, sgfc, sgfd, stat
    INTEGER(int_8) :: box_idx1, box_idx2, box_idx3, box_idx4, &
      estimate_to_store_int, mem_compression_counter, &
      mem_compression_counter_disk, mem_eris, mem_eris_disk, mem_maxvals, &
      my_current_counter, my_original_cpu_id, n_processes, ncpu, neris_disk, &
      neris_incore, neris_onthefly, neris_total, stor_count_int_disk, &
      storage_counter_integrals
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of, last_sgf_global
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, npgfa, npgfb, npgfc, npgfd, nsgfa, nsgfb, &
      nsgfc, nsgfd
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfc, first_sgfd
    LOGICAL :: buffer_overflow, do_correct_p_screening, do_p_screening, &
      do_periodic, failure, first_j_loop, first_k_loop, first_l_loop
    REAL(dp) :: actual_pmax_entry, bintime_start, bintime_stop, &
      cartesian_estimate, center1(3), center2(3), center2_folded(3), &
      compression_factor, compression_factor_disk, current_max_error, &
      ene_x_aa, ene_x_bb, eps_box, eps_schwarz, eps_storage, fac, &
      hf_fraction, Kbox, ks_bc, ks_bd, ks_beta_bc, ks_beta_bd, &
      max_contraction_val, max_val1, max_val2, omega2Inv, omega_box, p_bc, &
      p_bd, p_beta_bc, p_beta_bd, p_screening_correction_factor, &
      pbc_center(3), pmax_entry, R1, R2, ra(3), rb(3), Rbox, rc(3), rd(3), &
      screen_max, spherical_estimate, symm_fac, temp(3), tmp
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: primitive_integrals, T1
    REAL(dp), DIMENSION(:, :), POINTER :: full_density, full_density_beta, &
      full_ks, full_ks_beta, max_contraction, screen_atom, sparse_block, &
      sphi_a, sphi_b, sphi_c, sphi_d, zeta, zetb, zetc, zetd
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: actual_pmax_set, screen_set
    REAL(dp), DIMENSION(:, :, :, :, :), &
      POINTER                                :: screen_periodic
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(hfx_basis_info_type), POINTER       :: basis_info
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_cache_type), DIMENSION(:), &
      POINTER                                :: integral_caches, &
                                                integral_caches_disk
    TYPE(hfx_cache_type), POINTER            :: maxval_cache, &
                                                maxval_cache_disk
    TYPE(hfx_container_type), DIMENSION(:), &
      POINTER                                :: integral_containers, &
                                                integral_containers_disk
    TYPE(hfx_container_type), POINTER        :: maxval_container, &
                                                maxval_container_disk
    TYPE(hfx_general_type)                   :: general_parameter
    TYPE(hfx_load_balance_type)              :: load_balance_parameter
    TYPE(hfx_memory_type), POINTER           :: memory_parameter
    TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(hfx_screening_type)                 :: screening_parameter
    TYPE(hfx_type), POINTER                  :: actual_x_data
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_block_node_type), POINTER      :: block_node

!$  INTEGER :: omp_get_max_threads,omp_get_thread_num
    REAL(dp), DIMENSION(:, :, :, :), &
              POINTER, SAVE                  :: shm_initial_p
    REAL(dp), DIMENSION(:, :, :, :), &
              POINTER, SAVE                  :: shm_actual_pmax_set
    TYPE(hfx_box_scheme), DIMENSION(:,:), &
                           POINTER, SAVE     :: shm_far_field_boxes
    INTEGER, SAVE                            :: shm_number_of_p_entries
    INTEGER, DIMENSION(:,:), POINTER, SAVE   :: shm_is_assoc_atomic_block
    INTEGER(int_8)                           :: shm_neris_total, shm_neris_onthefly, &
                                                shm_storage_counter_integrals, shm_neris_incore, &
                                                shm_neris_disk, shm_stor_count_int_disk
    LOGICAL                                  :: do_disk_storage, use_disk_storage
    INTEGER(int_8)                           :: counter, neris_tmp, nprim_ints, shm_nprim_ints,&
                                                mb_size_p, mb_size_f, mb_size_s, mb_size_buffers
    INTEGER                                  :: nkind, l_max, istat
    REAL(dp)                                 :: ks_occ
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL                                  :: exist
    CHARACTER(LEN=512)                       :: error_msg
    INTEGER                                  :: unit_id

    CALL mp_sync(para_env%group)
    CALL timeset(routineN,handle)

    CALL cite_reference(Guidon2008)

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, cell=cell, error=error)

    !! Get occupation of KS-matrix
    IF( geometry_did_change) THEN
      ks_occ = get_occupation(ks_matrix(1)%matrix,para_env,error)
      IF( ks_occ < 100.0_dp ) THEN
        WRITE(error_msg,'(A,A,A)') "The Kohn Sham matrix is not 100% occupied. This "//&
              "may result in uncorrect Hartree-Fock results. Try to decrease EPS_PGF_ORB "//&
              "in the QS section. "
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                         error_msg,&
                         error,&
                         only_ionode=.TRUE.)
      END IF
    END IF

    !! Calculate l_max used in fgamma , because init_md_ftable is definitely not thread safe
    nkind = SIZE(atomic_kind_set,1)
    l_max = 0
    DO ikind=1,nkind
      l_max = MAX(l_max,MAXVAL(qs_env%x_data(1,1)%basis_parameter(ikind)%lmax))
    ENDDO
    l_max = 4*l_max 
    CALL init_md_ftable(l_max)

    IF(qs_env%x_data(1,1)%potential_parameter%potential_type == do_hfx_potential_truncated) THEN
      IF(l_max>init_t_c_g0_lmax) THEN
        unit_id = get_unit_number()
        OPEN(unit_id,FILE=qs_env%x_data(1,1)%potential_parameter%filename)
        CALL init(l_max,unit_id)
        CLOSE(unit_id)
        init_t_c_g0_lmax = l_max
      END IF
    CALL set_eps_cutoff(qs_env%x_data(1,1)%screening_parameter%eps_schwarz)
    END IF


    n_threads = 1
!$  n_threads = omp_get_max_threads()
    ALLOCATE(lib(n_threads), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    shm_neris_total = 0 
    shm_nprim_ints = 0
    shm_neris_onthefly = 0
    shm_storage_counter_integrals = 0
    shm_stor_count_int_disk = 0
    shm_neris_incore = 0
    shm_neris_disk = 0

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(qs_env,&
!$OMP                                  ks_matrix,&
!$OMP                                  energy,&
!$OMP                                  rho,&
!$OMP                                  hfx_section,&
!$OMP                                  para_env,&
!$OMP                                  geometry_did_change,&
!$OMP                                  irep,&
!$OMP                                  distribute_fock_matrix,&
!$OMP                                  error,&
!$OMP                                  ncoset,&
!$OMP                                  lib,&
!$OMP                                  full_ks,&
!$OMP                                  full_ks_beta,&
!$OMP                                  n_threads,&
!$OMP                                  full_density,&
!$OMP                                  full_density_beta,&
!$OMP                                  shm_initial_p,&
!$OMP                                  shm_actual_pmax_set,&
!$OMP                                  screen_max,&
!$OMP                                  screen_atom,&
!$OMP                                  screen_set,&
!$OMP                                  screen_periodic,&
!$OMP                                  shm_far_field_boxes,&
!$OMP                                  shm_is_assoc_atomic_block,&
!$OMP                                  shm_number_of_p_entries,&
!$OMP                                  shm_neris_total,&
!$OMP                                  shm_neris_onthefly,&
!$OMP                                  shm_storage_counter_integrals,&
!$OMP                                  shm_stor_count_int_disk,&
!$OMP                                  shm_neris_incore,&
!$OMP                                  shm_neris_disk,&
!$OMP                                  shm_nprim_ints,&
!$OMP                                  cell)


    i_thread = 0
!$  i_thread = omp_get_thread_num()
    
    current_max_error = 0.0_dp
    actual_x_data => qs_env%x_data(irep, i_thread + 1)
    do_periodic = actual_x_data%periodic_parameter%do_periodic

    screening_parameter    = actual_x_data%screening_parameter
    potential_parameter    = actual_x_data%potential_parameter
    IF(potential_parameter%omega /= 0.0_dp) THEN
      omega2Inv = 1.0_dp/(potential_parameter%omega**2)
    ELSE
      omega2Inv = 0.0_dp
    END IF

    general_parameter      = actual_x_data%general_parameter
    load_balance_parameter = actual_x_data%load_balance_parameter
    memory_parameter       => actual_x_data%memory_parameter

    maxval_container       => actual_x_data%maxval_container
    maxval_cache           => actual_x_data%maxval_cache

    integral_containers    => actual_x_data%integral_containers
    integral_caches        => actual_x_data%integral_caches

    maxval_container_disk       => actual_x_data%maxval_container_disk
    maxval_cache_disk           => actual_x_data%maxval_cache_disk

    integral_containers_disk    => actual_x_data%integral_containers_disk
    integral_caches_disk        => actual_x_data%integral_caches_disk


    basis_parameter   => actual_x_data%basis_parameter
    basis_info        => actual_x_data%basis_info

    !! initalize some counters
    neris_total = 0_int_8
    neris_incore = 0_int_8
    neris_disk = 0_int_8
    neris_onthefly = 0_int_8
    mem_eris = 0_int_8
    mem_eris_disk = 0_int_8
    mem_maxvals = 0_int_8
    compression_factor = 0.0_dp
    compression_factor_disk = 0.0_dp
    nprim_ints = 0_int_8


    use_disk_storage = .FALSE.
    counter = 0_int_8
    do_disk_storage = memory_parameter%do_disk_storage
    IF(geometry_did_change) THEN
     memory_parameter%ram_counter = HUGE(memory_parameter%ram_counter) 
    END IF
   
    !! We do not need the containers if MAX_MEM = 0
    IF(.NOT. memory_parameter%do_all_on_the_fly) THEN 
      !! IF new md step -> reinitialize containers
      IF(geometry_did_change) THEN
        CALL hfx_init_container(6, maxval_container, memory_parameter%actual_memory_usage, .FALSE., error)
        DO i=1,64
          CALL hfx_init_container(6, integral_containers(i),memory_parameter%actual_memory_usage, .FALSE., error)
        END DO
      END IF

      !! Decompress the first cache for maxvals and integrals
      IF( .NOT. geometry_did_change) THEN
        CALL hfx_decompress_first_cache(6, maxval_cache, maxval_container, memory_parameter%actual_memory_usage, .FALSE.)
        DO i=1,64
          CALL hfx_decompress_first_cache(i, integral_caches(i), integral_containers(i), &
                                           memory_parameter%actual_memory_usage,.FALSE.)
        END DO
      END IF
    END IF
    
    !! Since the I/O routines are no thread-safe, i.e. the procedure to get the unit number, put a lock here
!$OMP CRITICAL
    !! If we do disk storage, we have to initialize the containers/caches
    IF(do_disk_storage) THEN
      !! IF new md step -> reinitialize containers
      IF(geometry_did_change) THEN
        CALL hfx_init_container(6, maxval_container_disk, memory_parameter%actual_memory_usage_disk, do_disk_storage, error)
        DO i=1,64
          CALL hfx_init_container(i, integral_containers_disk(i),memory_parameter%actual_memory_usage_disk, do_disk_storage, error)
        END DO
      END IF
      !! Decompress the first cache for maxvals and integrals
      IF( .NOT. geometry_did_change) THEN
        CALL hfx_decompress_first_cache(6, maxval_cache_disk, maxval_container_disk, &
                                         memory_parameter%actual_memory_usage_disk, .TRUE.)
        DO i=1,64
          CALL hfx_decompress_first_cache(i, integral_caches_disk(i), integral_containers_disk(i), &
                                           memory_parameter%actual_memory_usage_disk, .TRUE.)
        END DO
      END IF
    END IF
!$OMP END CRITICAL

    !! Get screening parameter
    eps_schwarz = screening_parameter%eps_schwarz
    eps_box     = screening_parameter%eps_box
    !! get storage epsilon
    eps_storage = eps_schwarz*memory_parameter%eps_storage_scaling
    cache_size = memory_parameter%cache_size

    do_p_screening = screening_parameter%do_initial_p_screening
    do_correct_p_screening = screening_parameter%do_init_p_screen_corr
    p_screening_correction_factor = screening_parameter%initial_p_screening_correction

    !! If we have a hybrid functional, we may need only a fraction of exact exchange
    hf_fraction = general_parameter%fraction 
   
    !! The number of integrals that fit into the given MAX_MEMORY

    !! Parameters related to the potential 1/r, erf(wr)/r, erfc(wr/r)
    potential_parameter = actual_x_data%potential_parameter
    
    !! Variable to check if we calculate the integrals in-core or on the fly
    !! If TRUE -> on the fly
    IF(.NOT. memory_parameter%do_all_on_the_fly) THEN
      buffer_overflow = .FALSE. 
    ELSE
      buffer_overflow = .TRUE.
    END IF
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
   
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    error=error)
 
    natom = SIZE(particle_set,1)
    ALLOCATE(kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             kind_of=kind_of)

    maxsgf = basis_info%max_sgf
    max_am = basis_info%max_am
    max_set =basis_info%max_set
    !! Initialize memory used by libint
    CALL initialize_libint(lib(i_thread+1), max_am, error)

    nspins = qs_env%dft_control%nspins
    !! Helper array to map local basis function indeces to global ones 
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO iatom=1,natom
      ikind = kind_of(iatom)
      last_sgf_global(iatom) = last_sgf_global(iatom-1)+ basis_parameter(ikind)%nsgf_total
    END DO
!$OMP BARRIER
!$OMP MASTER
    !! Calculate helper array that stores if a certain atomic pair is associated in the KS matrix
    IF( geometry_did_change ) THEN
      actual_x_data%is_assoc_atomic_block = 0
      actual_x_data%number_of_p_entries = 0
      DO iatom=1,ks_matrix(1)%matrix%nblock_row
        block_node => first_block_node(ks_matrix(1)%matrix,iatom)
        ikind = kind_of(iatom)
        DO WHILE (ASSOCIATED(block_node))
          CALL get_block_node(block_node=block_node,&
                              block_col=jatom,&
                              BLOCK=sparse_block)
          jkind = kind_of(jatom)
          actual_x_data%is_assoc_atomic_block(iatom,jatom) = para_env%mepos + 1
          actual_x_data%is_assoc_atomic_block(jatom,iatom) = para_env%mepos + 1
          actual_x_data%number_of_p_entries =  actual_x_data%number_of_p_entries + &
                       basis_parameter(ikind)%nsgf_total*basis_parameter(jkind)%nsgf_total
          block_node => next_block_node(block_node)
        END DO
      END DO
      CALL mp_sum(actual_x_data%is_assoc_atomic_block,para_env%group)
      CALL mp_max(actual_x_data%number_of_p_entries, para_env%group)
      shm_number_of_p_entries = actual_x_data%number_of_p_entries
      shm_is_assoc_atomic_block => actual_x_data%is_assoc_atomic_block
    END IF

    !! Let master thread get the density (avoid problems with MPI)
    !! Get the full density from all the processors
    NULLIFY(full_density)
    NULLIFY(full_density_beta)
    CALL mp_sync(para_env%group) ! to clarify timings
    CALL timeset(routineN//"_getP",handle_getP)
    CALL get_full_density(para_env, full_density, rho%rho_ao(1)%matrix, shm_number_of_p_entries,&
                          shm_is_assoc_atomic_block, natom, last_sgf_global, &
                          kind_of, basis_parameter, error)   
    IF(nspins == 2) THEN
      CALL get_full_density(para_env, full_density_beta, rho%rho_ao(2)%matrix, shm_number_of_p_entries,&
                            shm_is_assoc_atomic_block, natom, last_sgf_global, &
                            kind_of, basis_parameter, error)
    END IF

    CALL symmetrize_matrix(full_density,"upper_to_lower")
    IF(nspins == 2) CALL symmetrize_matrix(full_density_beta,"upper_to_lower")


    CALL mp_sync(para_env%group) ! to clarify timings
    CALL timestop(handle_getP)

    !! Calculate the max values of the density matrix actual_pmax stores the data from the actual density matrix
    !! and x_data%initial_p stores the same for the initial guess. The initial guess is updated only in the case of
    !! a changed geometry 
    NULLIFY(shm_initial_p) 
    shm_initial_p => actual_x_data%initial_p
    IF( do_p_screening) THEN
      IF( do_correct_p_screening) THEN
        ALLOCATE(actual_pmax_set(max_set,max_set,natom,natom), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        NULLIFY(actual_pmax_set)
        shm_actual_pmax_set => actual_pmax_set
        CALL update_pmax_mat(actual_pmax_set, full_density, full_density_beta, &
                             particle_set,max_set, kind_of, basis_parameter, error)
      END IF
      IF(geometry_did_change ) THEN
        CALL update_pmax_mat(actual_x_data%initial_p, &
                             full_density,full_density_beta,&
                             particle_set,max_set, kind_of, basis_parameter, error)
      END IF
    ELSE
      IF(geometry_did_change) THEN
        actual_x_data%initial_p = 1.0_dp
      END IF  
    END IF
!$OMP END MASTER
!$OMP BARRIER

    ALLOCATE(max_contraction(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    max_contraction=0.0_dp
    DO jatom=1,natom
      jkind = kind_of(jatom)
      lb_max => basis_parameter(jkind)%lmax
      nsetb = basis_parameter(jkind)%nset
      npgfb => basis_parameter(jkind)%npgf
      first_sgfb => basis_parameter(jkind)%first_sgf
      sphi_b => basis_parameter(jkind)%sphi
      nsgfb => basis_parameter(jkind)%nsgf
      DO jset = 1,nsetb
        ! takes the primitive to contracted transformation into account
        ncob = npgfb(jset)*ncoset(lb_max(jset))
        sgfb = first_sgfb(1,jset)
        ! if the primitives are assumed to be all of max_val2, max_val2*p2s_b becomes
        ! the maximum value after multiplication with sphi_b
        max_contraction(jset,jatom) = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,i))),i=sgfb,sgfb+nsgfb(jset)-1)/))
      ENDDO
    ENDDO

!$OMP MASTER    
    NULLIFY(full_ks, full_ks_beta)
    ALLOCATE(actual_x_data%full_ks_alpha(SIZE(full_density,1),SIZE(full_density,2)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    full_ks => actual_x_data%full_ks_alpha
    full_ks = 0.0_dp
    IF(nspins==2) THEN
      ALLOCATE(actual_x_data%full_ks_beta(SIZE(full_density,1),SIZE(full_density,2)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      full_ks_beta => actual_x_data%full_ks_beta
      full_ks_beta = 0.0_dp
    END IF

    !! Initialize schwarz screening matrices for near field estimates and boxing screening matrices
    !! for far field estimates. The update is only performed if the geomtry of the system changed.
    !! If the system is periodic, then the corresponding routines are called and some variables 
    !! are initialized
    IF( do_periodic ) THEN
      ALLOCATE(screen_periodic(SIZE(actual_x_data%neighbor_cells),max_set,max_set,natom,natom),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL comp_screen_mat_per(lib(i_thread+1), screen_set,screen_atom,screen_max,screen_periodic,particle_set,max_set,&
                               potential_parameter,.NOT.geometry_did_change,&
                               kind_of, basis_parameter,&
                               actual_x_data%x_screening_buffer_periodic,&
                               actual_x_data%neighbor_cells, i_thread, cell, error)
      IF(geometry_did_change) THEN
        NULLIFY(shm_far_field_boxes)
        CALL screen_farfield_per(actual_x_data%far_field_boxes, &
                                 particle_set, max_set, &
                                 actual_x_data%neighbor_cells, kind_of, basis_parameter, &
                                 cell, error)
        shm_far_field_boxes => actual_x_data%far_field_boxes
       END IF

    ELSE
      CALL comp_screen_mat(lib(i_thread+1),screen_set,screen_atom,screen_max,particle_set,max_set,&
                           potential_parameter,.NOT.geometry_did_change,kind_of, basis_parameter, i_thread,&
                           actual_x_data%x_screening_buffer, error)
      IF(geometry_did_change) THEN
        NULLIFY(shm_far_field_boxes)
        CALL screen_farfield(actual_x_data%far_field_boxes,&
                             particle_set, max_set,&
                             kind_of, basis_parameter, error)
        shm_far_field_boxes => actual_x_data%far_field_boxes
      END IF
    END IF
    CALL timeset(routineN//"_main",handle_main)
!$OMP END MASTER
!$OMP BARRIER


    !! Initialize a prefactor depending on the fraction and the number of spins
    IF(nspins == 1) THEN
      fac = 0.5_dp * hf_fraction
    ELSE
      fac = 1.0_dp * hf_fraction
    END IF
    !! Call routines that distribute the load on all processes. If we want to screen on a initial density matrix, there is
    !! an optional parameter. Of course, this is only done if the geometry did change
    IF( geometry_did_change ) THEN
      IF( actual_x_data%b_first_load_balance_energy ) THEN
        CALL load_balance_energy(actual_x_data,eps_schwarz,eps_box,particle_set,max_set,para_env,&
                                 potential_parameter,screen_set,screen_atom,screen_max,&
                                 shm_is_assoc_atomic_block,do_periodic,load_balance_parameter, &
                                 kind_of, basis_parameter, shm_initial_p, i_thread, n_threads, &
                                 shm_far_field_boxes, cell, error)
        actual_x_data%b_first_load_balance_energy = .FALSE.
      ELSE
        
        CALL update_load_balance_energy(actual_x_data,para_env, &
                                        load_balance_parameter, &
                                        i_thread, n_threads, error)
      END IF
    END IF
    !! precompute maximum nco and allocate scratch
    nco_max=0

    DO iatom=1,natom
      ikind = kind_of(iatom)
      nseta = basis_parameter(ikind)%nset
      npgfa => basis_parameter(ikind)%npgf
      la_max => basis_parameter(ikind)%lmax
      DO iset = 1, nseta
         nco_max = MAX(nco_max,npgfa(iset)*ncoset(la_max(iset)))
      ENDDO
    ENDDO

    !! Allocate the arrays for the integrals. T1 is only used as a work array in the contraction using the DGEMM calls
    ALLOCATE(primitive_integrals(nco_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    primitive_integrals=0.0_dp
    ALLOCATE(T1(nco_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    T1=0.0_dp

    !! Start caluclating integrals of the form (ab|cd) or (ij|kl)
    !! In order to do so, there is a main four-loop structre that takes into account the two symmetries
    !!
    !!   (ab|cd) = (ba|cd) = (ab|dc) = (ba|dc)
    !!
    !! by iterating in the following way
    !!
    !! DO iatom=1,natom               and       DO katom=1,natom
    !!   DO jatom=iatom,natom                     DO latom=katom,natom
    !!
    !! The third symmetry
    !!  
    !!  (ab|cd) = (cd|ab)
    !!
    !! is taken into account by the following criterion:
    !!
    !! IF(katom+latom<=iatom+jatom)  THEN
    !!   IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
    !!
    !! Depending on the degeneracy of an integral the exchange contribution is multiplied by a corresponding
    !! factor ( symm_fac ).
    !!
    !! If one quartet does not pass the screening we CYCLE on the outer most possible loop. Thats why we use
    !! different hierarchies of short range screening matrices.
    !!
    !! If we do a parallel run, each process owns a unique array of workloads. Here, a workload is
    !! defined as :
    !!
    !! istart, jstart, kstart, lstart, number_of_integrals, initial_cpu_id
    !!
    !! This tells the process where to start the main loops and how many bunches of integrals it has to
    !! calculate. The original parallelization is a simple modulo distribution that is binned and 
    !! optimized in the load_balance routines. Since the Monte Carlo routines can swap processors, 
    !! we need to know which was the inital cpu_id.
    !! Furthermore, the indices iatom, jatom, katom, latom have to be set to istart, jstart, kstart and
    !! lstart only the first time the loop is executed. All subsequent loops have to start with one or
    !! iatom and katom respectively. Therefore, we use flags like first_j_loop etc.
    
    ncpu = para_env%num_pe 
    n_processes = ncpu * n_threads

    my_bin_size = SIZE(actual_x_data%distribution_energy)
    ! reset timings for the new SCF round
    IF (geometry_did_change) THEN
       DO bin=1,my_bin_size
          actual_x_data%distribution_energy(bin)%time_first_scf=0.0_dp
          actual_x_data%distribution_energy(bin)%time_other_scf=0.0_dp
          actual_x_data%distribution_energy(bin)%time_forces=0.0_dp
       ENDDO
    ENDIF
    DO bin=1,my_bin_size
      bintime_start=m_walltime()
      !! Get information on workload
      first_j_loop = .TRUE.
      first_l_loop = .TRUE.
      first_k_loop = .TRUE.
      my_istart = actual_x_data%distribution_energy(bin)%istart
      my_jstart = actual_x_data%distribution_energy(bin)%jstart
      my_kstart = actual_x_data%distribution_energy(bin)%kstart
      my_lstart = actual_x_data%distribution_energy(bin)%lstart
      my_original_cpu_id = actual_x_data%distribution_energy(bin)%original_cpu_id 
      my_current_counter = 0
      !! It could happen, that one bin is empty. In that case just skip the loop (but still time it)
      IF(my_istart == 0 .OR. my_jstart == 0 .OR. my_kstart == 0 .OR. my_lstart == 0 .OR. &
         actual_x_data%distribution_energy(bin)%number_of_integrals == 0) THEN
         my_istart = natom + 1 ; my_jstart = natom + 1 ; my_kstart = natom + 1 ; my_lstart = natom + 1
      ENDIF
      DO iatom=my_istart,natom
        ikind = kind_of(iatom)
        la_max => basis_parameter(ikind)%lmax
        la_min => basis_parameter(ikind)%lmin
        npgfa => basis_parameter(ikind)%npgf
        nseta = basis_parameter(ikind)%nset
        zeta => basis_parameter(ikind)%zet
        nsgfa => basis_parameter(ikind)%nsgf
        first_sgfa => basis_parameter(ikind)%first_sgf
        sphi_a => basis_parameter(ikind)%sphi
        nsgfa_total = basis_parameter(ikind)%nsgf_total
        ra=particle_set(iatom)%r(:)
        pa1 = last_sgf_global(iatom) - nsgfa_total
 
        IF(first_j_loop) jatom=my_jstart-1
        IF(.NOT. first_j_loop) jatom = iatom-1
        DO WHILE(jatom+1<=natom)
          jatom = jatom + 1
          !! Be consistent with the modulo distribution of the load balance
          box_idx1 = (iatom-1)*natom+jatom
          IF(MODULO(FNV64(box_idx1),n_processes)/=my_original_cpu_id-1) CYCLE
          jkind = kind_of(jatom)
          lb_max => basis_parameter(jkind)%lmax
          lb_min => basis_parameter(jkind)%lmin
          npgfb => basis_parameter(jkind)%npgf
          nsetb = basis_parameter(jkind)%nset
          zetb => basis_parameter(jkind)%zet
          nsgfb => basis_parameter(jkind)%nsgf
          first_sgfb => basis_parameter(jkind)%first_sgf
          sphi_b => basis_parameter(jkind)%sphi
          nsgfb_total = basis_parameter(jkind)%nsgf_total

          rb=particle_set(jatom)%r(:)
          pb1 = last_sgf_global(jatom) - nsgfb_total
          IF(first_k_loop) katom=my_kstart-1
          IF(.NOT. first_k_loop) katom = 0
          DO WHILE(katom+1<=natom)
            katom = katom + 1
            kkind = kind_of(katom)
            lc_max => basis_parameter(kkind)%lmax
            lc_min => basis_parameter(kkind)%lmin
            npgfc => basis_parameter(kkind)%npgf
            nsetc = basis_parameter(kkind)%nset
            zetc => basis_parameter(kkind)%zet
            nsgfc => basis_parameter(kkind)%nsgf
            first_sgfc => basis_parameter(kkind)%first_sgf
            sphi_c => basis_parameter(kkind)%sphi
            nsgfc_total = basis_parameter(kkind)%nsgf_total

            rc=particle_set(katom)%r(:)
            pc1 = last_sgf_global(katom) - nsgfc_total

            IF(first_l_loop) THEN
              latom=my_lstart-1
            ELSE
              latom = katom-1
            END IF
            DO WHILE(latom+1<=natom)
              latom = latom + 1
              !! third permutational symmetry
              IF(katom+latom<=iatom+jatom)  THEN
                IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                !! if we pass through here, we have to update the workload counter
                my_current_counter = my_current_counter + 1
                !! if we reach the counter in the distribution we are finished with the curren bunch
                !! of integrals. We also have to reset the flags for the main loops
                IF(my_current_counter <= actual_x_data%distribution_energy(bin)%number_of_integrals) THEN
                  
                  IF (screen_atom(jatom,iatom) * screen_max < eps_schwarz) CYCLE
  
                  !! screen near field atomic part
                  IF (screen_atom(jatom,iatom) * screen_atom(latom,katom) < eps_schwarz) CYCLE
  
                  !! we want to be consistent with the KS matrix. If none of the atomic indices
                  !! is associated cycle
                  IF(.NOT. (shm_is_assoc_atomic_block(latom,iatom)>=1 .AND. &
                            shm_is_assoc_atomic_block(katom,iatom)>=1 .AND. &
                            shm_is_assoc_atomic_block(katom,jatom)>=1 .AND. &
                            shm_is_assoc_atomic_block(latom,jatom)>=1 ) ) CYCLE
  
                  !! calculate symmetry_factor accordin to degeneracy of atomic quartet
                  symm_fac = 0.5_dp
                  IF(iatom==jatom) symm_fac = symm_fac*2.0_dp
                  IF(katom==latom) symm_fac = symm_fac*2.0_dp
                  IF(iatom==katom .AND. jatom==latom .AND. iatom/=jatom .AND. katom/=latom) symm_fac = symm_fac*2.0_dp
                  IF(iatom==katom .AND. iatom==jatom .AND. katom==latom) symm_fac = symm_fac*2.0_dp
                  symm_fac = 1.0_dp / symm_fac
  
                  lkind = kind_of(latom)
                  ld_max => basis_parameter(lkind)%lmax
                  ld_min => basis_parameter(lkind)%lmin
                  npgfd => basis_parameter(lkind)%npgf
                  nsetd = basis_parameter(lkind)%nset
                  zetd => basis_parameter(lkind)%zet
                  nsgfd => basis_parameter(lkind)%nsgf
                  first_sgfd => basis_parameter(lkind)%first_sgf
                  sphi_d => basis_parameter(lkind)%sphi
                  nsgfd_total = basis_parameter(lkind)%nsgf_total
  
                  rd=particle_set(latom)%r(:)
                  pd1 = last_sgf_global(latom) - nsgfd_total
                  box_idx3 = (katom-1)*natom+latom
                  
                  !! Loop over all sets of the current atomic quartet
    
                  !! At this stage, check for memory used in compression

                  IF(geometry_did_change ) THEN
                    IF(.NOT.memory_parameter%do_all_on_the_fly) THEN
                      mem_compression_counter = memory_parameter%actual_memory_usage *&
                                                memory_parameter%cache_size
                      IF(mem_compression_counter > memory_parameter%max_compression_counter) THEN
                        buffer_overflow = .TRUE.
                        memory_parameter%ram_counter = counter
                      ELSE
                        counter = counter + 1
                        buffer_overflow = .FALSE.
                      END IF 
                    END IF
                  ELSE
                    IF(.NOT.memory_parameter%do_all_on_the_fly) THEN
                      IF(memory_parameter%ram_counter == counter) THEN
                        buffer_overflow = .TRUE.
                      ELSE
                        counter = counter + 1
                       buffer_overflow = .FALSE.
                      END IF
                    END IF
                  END IF

                  IF( buffer_overflow .AND. do_disk_storage ) THEN
                    use_disk_storage = .TRUE.
                    buffer_overflow = .FALSE.
                  END IF
                  
                  IF( use_disk_storage ) THEN
                    mem_compression_counter_disk = memory_parameter%actual_memory_usage_disk *&
                                                   memory_parameter%cache_size
                    IF(mem_compression_counter_disk > memory_parameter%max_compression_counter_disk) THEN
                      buffer_overflow = .TRUE.
                      use_disk_storage = .FALSE.
                    END IF
                  END IF
                 
                  DO iset = 1, nseta
                    ncoa = npgfa(iset)*ncoset(la_max(iset))
                    sgfa = first_sgfa(1,iset)
                    DO jset = 1,nsetb
                      ncob = npgfb(jset)*ncoset(lb_max(jset))
                      sgfb = first_sgfb(1,jset)
                      max_val1 = screen_set(jset,iset,jatom,iatom)
                      !! Near field screening
                      IF (max_val1 * screen_atom(latom,katom) < eps_schwarz) CYCLE
                      box_idx2 = (iset-1)*max_set+jset
                      center1 = shm_far_field_boxes(box_idx2, box_idx1)%center
                      R1 = shm_far_field_boxes(box_idx2, box_idx1)%R
                      radius1 = shm_far_field_boxes(box_idx2, box_idx1)%radius
                      DO kset = 1,nsetc
                        ncoc = npgfc(kset)*ncoset(lc_max(kset))
                        sgfc = first_sgfc(1,kset)
                        DO lset = 1,nsetd
                          ncod = npgfd(lset)*ncoset(ld_max(lset)) 
                          sgfd = first_sgfd(1,lset)
                          max_val2 = max_val1 * screen_set(lset,kset,latom,katom)
                          !! Near field screening
                          IF(max_val2<eps_schwarz) CYCLE
                          box_idx4 = (kset-1)*max_set+lset
                          IF(.NOT. do_periodic) THEN
                            center2_folded = shm_far_field_boxes(box_idx4, box_idx3)%center
                          ELSE
                            center2 = shm_far_field_boxes(box_idx4, box_idx3)%center
                            temp = center2-center1
                            pbc_center = pbc(temp,cell)
                            center2_folded = center1 + pbc_center
                          END IF
  
                          R2 = shm_far_field_boxes(box_idx4, box_idx3)%R
                          radius2 = shm_far_field_boxes(box_idx4, box_idx3)%radius
                          !! get max_vals if we screen on initial density
                          pmax_entry = MAX(shm_initial_p(kset,iset,katom,iatom), &
                                           shm_initial_p(lset,jset,latom,jatom), &
                                           shm_initial_p(lset,iset,latom,iatom), &
                                           shm_initial_p(kset,jset,katom,jatom))
                          IF( do_correct_p_screening ) THEN
                            actual_pmax_entry = MAX(shm_actual_pmax_set(kset,iset,katom,iatom), &
                                                    shm_actual_pmax_set(lset,jset,latom,jatom), &
                                                    shm_actual_pmax_set(lset,iset,latom,iatom), &
                                                    shm_actual_pmax_set(kset,jset,katom,jatom))
                          END IF
                          max_val2 = max_val2 * pmax_entry
                          !! If wanted, correct the previous screening by recalculating the omitted
                          !! integrals on the fly.
                          IF(max_val2<eps_schwarz) THEN
                            IF( do_correct_p_screening ) THEN
                              IF(pmax_entry*p_screening_correction_factor < actual_pmax_entry) THEN
                                max_contraction_val =  max_contraction(iset,iatom) * &
                                                       max_contraction(jset,jatom) * &
                                                       max_contraction(kset,katom) * &
                                                       max_contraction(lset,latom) * pmax_entry
                                CALL on_the_fly(do_periodic, ra, rb, rc, rd, zeta, zetb, zetc, zetd, &
                                                la_min, la_max, lb_min, lb_max, lc_min, lc_max, ld_min, ld_max, &
                                                npgfa, npgfb, npgfc, npgfd, ncoa, ncob, ncoc, ncod, &
                                                nsgfa, nsgfb, nsgfc, nsgfd, sgfa, sgfb, sgfc, sgfd, &
                                                primitive_integrals, potential_parameter, &
                                                screen_periodic, eps_schwarz, sphi_a, sphi_b, sphi_c, sphi_d, T1,&
                                                full_ks, full_ks_beta, full_density, full_density_beta, fac, symm_fac, &
                                                iset, jset, kset, lset, iatom, jatom, katom, latom, pa1, pb1, pc1, pd1,&
                                                nspins, max_contraction_val, actual_x_data, i_thread, cell, neris_tmp,&
                                                error)
                                nprim_ints = nprim_ints + neris_tmp
                              END IF
                            END IF
                            CYCLE
                          END IF
  
                          !! Far field box screening also taking into account screening on initial density matrix
                          Rbox = (center1(1)-center2_folded(1))**2 + &
                                 (center1(2)-center2_folded(2))**2 + &
                                 (center1(3)-center2_folded(3))**2
                          Rbox = SQRT(Rbox) - R1 - R2 
                          pgf_radius_box = radius1 + radius2 + 2
                          IF(Rbox>pgf_radius_box) THEN
                            max_val2 = 0.0_dp
                            SELECT CASE(potential_parameter%potential_type)
                              CASE(do_hfx_potential_truncated)
                              Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                     shm_far_field_boxes(box_idx4, box_idx3)%K
                              max_val2 = 1.0_dp/Rbox*Kbox
                              IF( max_val2 < eps_box ) CYCLE
                              max_val2 = max_val2 * pmax_entry
                              CASE(do_hfx_potential_coulomb)
                                Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                       shm_far_field_boxes(box_idx4, box_idx3)%K
                                max_val2 = 1.0_dp/Rbox*Kbox
                                IF(max_val2 < eps_box ) CYCLE
                                max_val2 = max_val2 * pmax_entry
                              CASE(do_hfx_potential_short)
                                omega_box = shm_far_field_boxes(box_idx2, box_idx1)%alphaInv + &
                                            shm_far_field_boxes(box_idx4, box_idx3)%alphaInv + &
                                            omega2Inv 
                                omega_box = 1.0_dp/omega_box
                                Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                       shm_far_field_boxes(box_idx4, box_idx3)%K
                                max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox
                                IF(max_val2 < eps_box ) CYCLE
                                max_val2 = max_val2 * pmax_entry
                              CASE(do_hfx_potential_mix_cl)
                                Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                       shm_far_field_boxes(box_idx4, box_idx3)%K
                                max_val2 = 1.0_dp/Rbox*Kbox
                                IF(max_val2 < eps_box ) CYCLE
                                max_val2 = max_val2 * pmax_entry
                              CASE (do_hfx_potential_gaussian, do_hfx_potential_mix_lg, do_hfx_potential_id)
                                max_val2 = HUGE(max_val2)
                            END SELECT
                            IF(max_val2 < eps_box) THEN
                              IF( do_correct_p_screening ) THEN
                                IF(pmax_entry*p_screening_correction_factor < actual_pmax_entry) THEN
                                  max_contraction_val =  max_contraction(iset,iatom) * &
                                                         max_contraction(jset,jatom) * &
                                                         max_contraction(kset,katom) * &
                                                         max_contraction(lset,latom) * pmax_entry
                                  CALL on_the_fly(do_periodic, ra, rb, rc, rd, zeta, zetb, zetc, zetd, &
                                                  la_min, la_max, lb_min, lb_max, lc_min, lc_max, ld_min, ld_max, &
                                                  npgfa, npgfb, npgfc, npgfd, ncoa, ncob, ncoc, ncod, &
                                                  nsgfa, nsgfb, nsgfc, nsgfd, sgfa, sgfb, sgfc, sgfd, &
                                                  primitive_integrals, potential_parameter, &
                                                  screen_periodic, eps_schwarz, sphi_a, sphi_b, sphi_c, sphi_d, T1, &
                                                  full_ks, full_ks_beta, full_density, full_density_beta, fac, symm_fac, &
                                                  iset, jset, kset, lset, iatom, jatom, katom, latom, pa1, pb1, pc1, pd1,&
                                                  nspins, max_contraction_val, actual_x_data, i_thread, cell, neris_tmp,&
                                                  error)
                                   nprim_ints = nprim_ints + neris_tmp
                                END IF 
                              END IF
                              CYCLE
                            END IF
                          END IF
  
                          !! store current number of integrals, update total number and number of integrals in buffer
                          current_counter = nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                          IF(buffer_overflow) THEN
                            neris_onthefly = neris_onthefly + current_counter 
                          END IF
  
                          !! Get integrals from buffer and update Kohn-Sham matrix
                          IF(.NOT.buffer_overflow .AND. .NOT.geometry_did_change) THEN
                            nints = current_counter
                            IF(.NOT. use_disk_storage) THEN
                              CALL hfx_get_single_cache_element(estimate_to_store_int, 6,&
                                                       maxval_cache, maxval_container, memory_parameter%actual_memory_usage, &
                                                       use_disk_storage)
                            ELSE
                               CALL hfx_get_single_cache_element(estimate_to_store_int, 6,&
                                                       maxval_cache_disk, maxval_container_disk, &
                                                       memory_parameter%actual_memory_usage_disk, &
                                                       use_disk_storage)
                            END IF
                            spherical_estimate = SET_EXPONENT(1.0_dp,estimate_to_store_int+1)
                            IF(spherical_estimate*pmax_entry<eps_schwarz) CYCLE
                            nbits = EXPONENT(spherical_estimate*pmax_entry/eps_storage) + 1                         
                            buffer_left = nints
                            buffer_start = 1
                            DO WHILE (buffer_left > 0 )
                              buffer_size = MIN(buffer_left, cache_size)
                              IF(.NOT. use_disk_storage) THEN
                                CALL hfx_get_mult_cache_elements(T1(buffer_start), &
                                                                  buffer_size, nbits, &
                                                                  integral_caches(nbits), &
                                                                  integral_containers(nbits), &
                                                                  eps_storage, pmax_entry, &
                                                                  memory_parameter%actual_memory_usage, &
                                                                  use_disk_storage, error)
                              ELSE
                                CALL hfx_get_mult_cache_elements(T1(buffer_start), &
                                                                  buffer_size, nbits, &
                                                                  integral_caches_disk(nbits), &
                                                                  integral_containers_disk(nbits), &
                                                                  eps_storage, pmax_entry, &
                                                                  memory_parameter%actual_memory_usage_disk, &
                                                                  use_disk_storage, error)
                              END IF
                              buffer_left = buffer_left - buffer_size
                              buffer_start = buffer_start + buffer_size
                            END DO 
                            i=1
                            IF(nspins == 1) THEN
                              DO md = 1,nsgfd(lset)
                                pd2 = md+sgfd-1
                                pd = pd1 + pd2
                                DO mc = 1,nsgfc(kset)
                                  pc2 = mc+sgfc-1
                                  pc = pc1 + pc2
                                  DO mb = 1,nsgfb(jset)
                                    pb2 = mb+sgfb-1
                                    pb = pb1 + pb2
                                    ks_bd = 0.0_dp
                                    ks_bc = 0.0_dp
                                    p_bd = full_density(pb,pd)
                                    p_bc = full_density(pb,pc)
                                    DO ma = 1,nsgfa(iset)
                                      pa2 = ma+sgfa-1
                                      pa = pa1 + pa2
                                      tmp = fac*symm_fac*T1(i)
                                      i=i+1
                                      ks_bc = ks_bc + tmp*full_density(pa,pd)
                                      ks_bd = ks_bd + tmp*full_density(pa,pc)
                                      !$OMP ATOMIC
                                      full_ks(pa,pd) = full_ks(pa,pd)-tmp*p_bc
                                      !$OMP ATOMIC
                                      full_ks(pa,pc) = full_ks(pa,pc)-tmp*p_bd
                                    END DO
                                    !$OMP ATOMIC 
                                    full_ks(pb,pd) = full_ks(pb,pd) - ks_bd
                                    !$OMP ATOMIC
                                    full_ks(pb,pc) = full_ks(pb,pc) - ks_bc
                                  END DO
                                END DO
                              END DO
                            ELSE
                              DO md = 1,nsgfd(lset)
                                pd2 = md+sgfd-1
                                pd = pd1 + pd2
                                DO mc = 1,nsgfc(kset)
                                  pc2 = mc+sgfc-1
                                  pc = pc1 + pc2
                                  DO mb = 1,nsgfb(jset)
                                    pb2 = mb+sgfb-1
                                    pb = pb1 + pb2
                                    ks_bd = 0.0_dp
                                    ks_bc = 0.0_dp
                                    p_bd = full_density(pb,pd)
                                    p_bc = full_density(pb,pc)
                                    ks_beta_bd = 0.0_dp
                                    ks_beta_bc = 0.0_dp
                                    p_beta_bd = full_density_beta(pb,pd)
                                    p_beta_bc = full_density_beta(pb,pc)
                                    DO ma = 1,nsgfa(iset)
                                      pa2 = ma+sgfa-1
                                      pa = pa1 + pa2
                                      tmp = fac*symm_fac*T1(i)
                                      i=i+1
                                      ks_bc = ks_bc + tmp*full_density(pa,pd)
                                      ks_bd = ks_bd + tmp*full_density(pa,pc)
                                      !$OMP ATOMIC
                                      full_ks(pa,pd) = full_ks(pa,pd)-tmp*p_bc
                                      !$OMP ATOMIC
                                      full_ks(pa,pc) = full_ks(pa,pc)-tmp*p_bd
                                      ks_beta_bc = ks_beta_bc + tmp*full_density_beta(pa,pd)
                                      ks_beta_bd = ks_beta_bd + tmp*full_density_beta(pa,pc)
                                      !$OMP ATOMIC
                                      full_ks_beta(pa,pd) = full_ks_beta(pa,pd)-tmp*p_beta_bc
                                      !$OMP ATOMIC
                                      full_ks_beta(pa,pc) = full_ks_beta(pa,pc)-tmp*p_beta_bd
                                    END DO
                                    !$OMP ATOMIC
                                    full_ks(pb,pd) = full_ks(pb,pd) - ks_bd
                                    !$OMP ATOMIC
                                    full_ks(pb,pc) = full_ks(pb,pc) - ks_bc
                                    !$OMP ATOMIC
                                    full_ks_beta(pb,pd) = full_ks_beta(pb,pd) - ks_beta_bd
                                    !$OMP ATOMIC
                                    full_ks_beta(pb,pc) = full_ks_beta(pb,pc) - ks_beta_bc
                                  END DO
                                END DO
                              END DO
                            END IF
                          END IF
    
                          !! Calculate integrals if we run out of buffer or the geometry did change
                          IF(geometry_did_change .OR. buffer_overflow) THEN
                            max_contraction_val =  max_contraction(iset,iatom) * &
                                                   max_contraction(jset,jatom) * &
                                                   max_contraction(kset,katom) * &
                                                   max_contraction(lset,latom) * pmax_entry
  
                            IF( .NOT.do_periodic ) THEN 
                              CALL coulomb4(ra, rb, rc, rd,&
                                            zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset), &
                                            zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                                            la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                            lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                            npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                            ncoa,ncob,ncoc,ncod,&
                                            primitive_integrals(1),&
                                            potential_parameter, eps_schwarz, max_contraction_val, &
                                            cartesian_estimate, i_thread, neris_tmp)
                            ELSE
                              CALL coulomb4_per(ra, rb, rc, rd,&
                                                zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset), &
                                                zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                                                la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                                lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                                npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                                ncoa,ncob,ncoc,ncod,&
                                                primitive_integrals(1),&
                                                potential_parameter, &
                                                actual_x_data%neighbor_cells, screen_periodic(:,iset,jset,iatom,jatom),&
                                                screen_periodic(:,kset,lset,katom,latom), eps_schwarz, &
                                                max_contraction_val, cartesian_estimate, i_thread, cell, neris_tmp,&
                                                error)
                           
                            END IF
                            
                            nints = nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                            neris_total = neris_total + nints
                            nprim_ints = nprim_ints + neris_tmp
                            IF(cartesian_estimate == 0.0_dp) cartesian_estimate = TINY(cartesian_estimate)
                            estimate_to_store_int = EXPONENT(cartesian_estimate)
                            estimate_to_store_int = MAX(estimate_to_store_int,-15_int_8)
                            cartesian_estimate = SET_EXPONENT(1.0_dp,estimate_to_store_int+1)
                            IF(.NOT.buffer_overflow .AND. geometry_did_change) THEN
                              IF(cartesian_estimate < eps_schwarz ) THEN
                                IF(.NOT. use_disk_storage) THEN
                                  CALL hfx_add_single_cache_element(estimate_to_store_int, 6,&
                                                  maxval_cache, maxval_container, memory_parameter%actual_memory_usage, &
                                                  use_disk_storage, error)
                                ELSE
                                  CALL hfx_add_single_cache_element(estimate_to_store_int, 6,&
                                                  maxval_cache_disk, maxval_container_disk, &
                                                  memory_parameter%actual_memory_usage_disk, &
                                                  use_disk_storage, error)
                                END IF
                              END IF
                            END IF

                            IF(cartesian_estimate<eps_schwarz) CYCLE 
                            CALL dgemm("T","N",ncob*ncoc*ncod,nsgfa(iset),ncoa,&
                                       1.0_dp, primitive_integrals(1),ncoa,&
                                       sphi_a(1,sgfa), SIZE(sphi_a,1),&
                                       0.0_dp, T1(1),ncob*ncoc*ncod)
    
                            CALL dgemm("T","N",nsgfa(iset)*ncoc*ncod,nsgfb(jset),ncob,&
                                       1.0_dp, T1(1),ncob,&
                                       sphi_b(1,sgfb), SIZE(sphi_b,1),&
                                       0.0_dp, primitive_integrals(1),nsgfa(iset)*ncoc*ncod)
    
                            CALL dgemm("T","N",nsgfa(iset)*nsgfb(jset)*ncod,nsgfc(kset),ncoc,&
                                       1.0_dp, primitive_integrals(1),ncoc,&
                                       sphi_c(1,sgfc), SIZE(sphi_c,1),&
                                       0.0_dp, T1(1),nsgfa(iset)*nsgfb(jset)*ncod)
    
                            CALL dgemm("T","N",nsgfa(iset)*nsgfb(jset)*nsgfc(kset),nsgfd(lset),ncod,&
                                       1.0_dp, T1(1),ncod,&
                                       sphi_d(1,sgfd), SIZE(sphi_d,1),&
                                       0.0_dp, primitive_integrals(1),nsgfa(iset)*nsgfb(jset)*nsgfc(kset))
                            !! Compress the array for storage
                            spherical_estimate = 0.0_dp
                            DO i=1,nints
                              spherical_estimate = MAX(spherical_estimate,ABS(primitive_integrals(i)))
                            END DO
                            IF(spherical_estimate == 0.0_dp) spherical_estimate = TINY(spherical_estimate)
                            estimate_to_store_int = EXPONENT(spherical_estimate)
                            estimate_to_store_int = MAX(estimate_to_store_int, -15_int_8)
                            
                            IF(.NOT. buffer_overflow .AND. geometry_did_change) THEN
                              IF(.NOT.use_disk_storage) THEN
                                CALL hfx_add_single_cache_element(estimate_to_store_int, 6, &
                                                                   maxval_cache, maxval_container, &
                                                                   memory_parameter%actual_memory_usage, &
                                                                   use_disk_storage, error)
                              ELSE
                                CALL hfx_add_single_cache_element(estimate_to_store_int, 6, &
                                                                   maxval_cache_disk, maxval_container_disk, &
                                                                   memory_parameter%actual_memory_usage_disk, &
                                                                   use_disk_storage, error)
                              END IF
                            END IF
                            spherical_estimate = SET_EXPONENT(1.0_dp,estimate_to_store_int+1)
                            IF(spherical_estimate*pmax_entry<eps_schwarz) CYCLE 
                            IF(.NOT.buffer_overflow) THEN
                              nbits = EXPONENT(spherical_estimate*pmax_entry/eps_storage) + 1
                              buffer_left = nints
                              buffer_start = 1
                              DO WHILE (buffer_left > 0 )
                                buffer_size = MIN(buffer_left, CACHE_SIZE)
                                IF(.NOT. use_disk_storage) THEN
                                  neris_incore = neris_incore + nints
                                  CALL hfx_add_mult_cache_elements(primitive_integrals(buffer_start),&
                                                                    buffer_size, nbits, &
                                                                    integral_caches(nbits), &
                                                                    integral_containers(nbits), &
                                                                    eps_storage, pmax_entry, &
                                                                    memory_parameter%actual_memory_usage, &
                                                                    use_disk_storage, error)
                                ELSE
                                  neris_disk = neris_disk + nints
                                  CALL hfx_add_mult_cache_elements(primitive_integrals(buffer_start),&
                                                                    buffer_size, nbits, &
                                                                    integral_caches_disk(nbits), &
                                                                    integral_containers_disk(nbits), &
                                                                    eps_storage, pmax_entry, &
                                                                    memory_parameter%actual_memory_usage_disk, &
                                                                    use_disk_storage, error)
                                END IF
                                buffer_left = buffer_left - buffer_size
                                buffer_start = buffer_start + buffer_size
                              END DO
                            ELSE
                             !! In order to be consistent with in-core part, round all the eris wrt. eps_schwarz
                             DO i=1,nints
                               primitive_integrals(i) = primitive_integrals(i)*pmax_entry
                               IF(ABS(primitive_integrals(i)) > eps_storage) THEN
                                 primitive_integrals(i) = ANINT(primitive_integrals(i)/eps_storage, dp)*eps_storage/pmax_entry
                               ELSE
                                 primitive_integrals(i) = 0.0_dp
                               END IF
                             END DO
                            END IF
                            !! Update Kohn-Sham matrix
                            p_index=0 
                            IF(nspins == 1) THEN
                              DO md = 1,nsgfd(lset)
                                pd2 = md+sgfd-1
                                pd = pd1 + pd2
                                DO mc = 1,nsgfc(kset)
                                  pc2 = mc+sgfc-1
                                  pc = pc1 + pc2
                                  DO mb = 1,nsgfb(jset)
                                    pb2 = mb+sgfb-1
                                    pb = pb1 + pb2
                                    ks_bd = 0.0_dp
                                    ks_bc = 0.0_dp
                                    p_bd = full_density(pb,pd)
                                    p_bc = full_density(pb,pc)
                                    DO ma = 1,nsgfa(iset)
                                      pa2 = ma+sgfa-1
                                      pa = pa1 + pa2
                                      p_index=p_index+1
                                      tmp = fac*symm_fac*primitive_integrals(p_index)
                                      ks_bc = ks_bc + tmp*full_density(pa,pd)
                                      ks_bd = ks_bd + tmp*full_density(pa,pc)
                                      !$OMP ATOMIC 
                                      full_ks(pa,pd) = full_ks(pa,pd)-tmp*p_bc
                                      !$OMP ATOMIC
                                      full_ks(pa,pc) = full_ks(pa,pc)-tmp*p_bd
                                    END DO
                                    !$OMP ATOMIC 
                                    full_ks(pb,pd) = full_ks(pb,pd) - ks_bd
                                    !$OMP ATOMIC
                                    full_ks(pb,pc) = full_ks(pb,pc) - ks_bc
                                  END DO
                                END DO
                              END DO
                            ELSE
                              DO md = 1,nsgfd(lset)
                                pd2 = md+sgfd-1
                                pd = pd1 + pd2
                                DO mc = 1,nsgfc(kset)
                                  pc2 = mc+sgfc-1
                                  pc = pc1 + pc2
                                  DO mb = 1,nsgfb(jset)
                                    pb2 = mb+sgfb-1
                                    pb = pb1 + pb2
                                    ks_bd = 0.0_dp
                                    ks_bc = 0.0_dp
                                    p_bd = full_density(pb,pd)
                                    p_bc = full_density(pb,pc)
                                    ks_beta_bd = 0.0_dp
                                    ks_beta_bc = 0.0_dp
                                    p_beta_bd = full_density_beta(pb,pd)
                                    p_beta_bc = full_density_beta(pb,pc)
                                    DO ma = 1,nsgfa(iset)
                                      pa2 = ma+sgfa-1
                                      pa = pa1 + pa2
                                      p_index=p_index+1
                                      tmp = fac*symm_fac*primitive_integrals(p_index)
                                      ks_bc = ks_bc + tmp*full_density(pa,pd)
                                      ks_bd = ks_bd + tmp*full_density(pa,pc)
                                      !$OMP ATOMIC
                                      full_ks(pa,pd) = full_ks(pa,pd)-tmp*p_bc
                                      !$OMP ATOMIC
                                      full_ks(pa,pc) = full_ks(pa,pc)-tmp*p_bd
                                      ks_beta_bc = ks_beta_bc + tmp*full_density_beta(pa,pd)
                                      ks_beta_bd = ks_beta_bd + tmp*full_density_beta(pa,pc)
                                      !$OMP ATOMIC
                                      full_ks_beta(pa,pd) = full_ks_beta(pa,pd)-tmp*p_beta_bc
                                      !$OMP ATOMIC
                                      full_ks_beta(pa,pc) = full_ks_beta(pa,pc)-tmp*p_beta_bd
                                    END DO
                                    !$OMP ATOMIC
                                    full_ks(pb,pd) = full_ks(pb,pd) - ks_bd
                                    !$OMP ATOMIC
                                    full_ks(pb,pc) = full_ks(pb,pc) - ks_bc
                                    !$OMP ATOMIC
                                    full_ks_beta(pb,pd) = full_ks_beta(pb,pd) - ks_beta_bd
                                    !$OMP ATOMIC
                                    full_ks_beta(pb,pc) = full_ks_beta(pb,pc) - ks_beta_bc
                                  END DO
                                END DO
                              END DO
                            END IF
                          END IF
                        END DO !lset
                      END DO !kset
                    END DO !jset
                  END DO !iset
                  IF( do_disk_storage ) THEN
                    buffer_overflow = .TRUE.
                  END IF
                END IF
              END IF !symmetry cycle 
            END DO !latom
            first_l_loop = .FALSE.
          END DO !katom
          first_k_loop = .FALSE.
        END DO !jatom
        first_j_loop = .FALSE.
      END DO !iatom

      bintime_stop = m_walltime()
      IF (geometry_did_change) THEN
         actual_x_data%distribution_energy(bin)%time_first_scf=bintime_stop-bintime_start 
      ELSE
         actual_x_data%distribution_energy(bin)%time_other_scf= &
                actual_x_data%distribution_energy(bin)%time_other_scf+bintime_stop-bintime_start 
      ENDIF
    END DO !bin
!$OMP BARRIER

    DEALLOCATE(primitive_integrals,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(T1,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP BARRIER
    !! Get some number about ERIS
!$OMP ATOMIC
    shm_neris_total = shm_neris_total + neris_total
!$OMP ATOMIC
    shm_neris_onthefly = shm_neris_onthefly + neris_onthefly
!$OMP ATOMIC
    shm_nprim_ints = shm_nprim_ints + nprim_ints
    storage_counter_integrals = memory_parameter%actual_memory_usage * &
                                memory_parameter%cache_size
    stor_count_int_disk = memory_parameter%actual_memory_usage_disk * &
                          memory_parameter%cache_size
!$OMP ATOMIC
    shm_storage_counter_integrals = shm_storage_counter_integrals + storage_counter_integrals
!$OMP ATOMIC
    shm_stor_count_int_disk =  shm_stor_count_int_disk + stor_count_int_disk
!$OMP ATOMIC
    shm_neris_incore = shm_neris_incore + neris_incore
!$OMP ATOMIC
    shm_neris_disk = shm_neris_disk + neris_disk
!$OMP BARRIER

!$OMP MASTER   
    !! Calculate the exchange energies from the Kohn-Sham matrix. Before we can go on, we have to symmetrize.
    ene_x_aa = 0.0_dp
    ene_x_bb = 0.0_dp
   
    mb_size_p = SIZE(full_density,1)**2/1024/128
    mb_size_f = SIZE(full_ks,1)**2/1024/128
    IF( nspins == 2) THEN
      mb_size_f = mb_size_f * 2
      mb_size_p = mb_size_p * 2
    END IF
    !! size of primitive_integrals and T1 buffer
    mb_size_buffers = 2*INT(nco_max,int_8)**4
    !! screening buffers
    IF( do_periodic ) THEN
      mb_size_buffers = mb_size_buffers + SIZE(actual_x_data%neighbor_cells)*INT(max_set,int_8)**2*INT(natom,int_8)**2
    ENDIF
    !! non periodic screening buffers
    mb_size_buffers = mb_size_buffers + INT(max_set,int_8)**2*INT(natom,int_8)**2 + INT(natom,int_8)**2
    !! mb
    mb_size_buffers = mb_size_buffers /1024/128

    DO pa=1,SIZE(full_ks,1)
      DO pd=pa,SIZE(full_ks,2)
        IF(pa/=pd) THEN
          full_ks(pa,pd) = (full_ks(pa,pd) + full_ks(pd,pa))*0.5_dp
          full_ks(pd,pa) = full_ks(pa,pd)
          ene_x_aa = ene_x_aa + 2.0_dp * full_density(pa,pd)*full_ks(pa,pd)
          IF(nspins==2) THEN
            full_ks_beta(pa,pd) = (full_ks_beta(pa,pd) + &
                                   full_ks_beta(pd,pa)) * 0.5_dp
            full_ks_beta(pd,pa) = full_ks_beta(pa,pd) 
            ene_x_bb = ene_x_bb + 2.0_dp * full_density_beta(pa,pd)*full_ks_beta(pa,pd)
          END IF
        ELSE
          ene_x_aa = ene_x_aa + full_density(pa,pd)*full_ks(pa,pd)
          IF(nspins==2) ene_x_bb = ene_x_bb + full_density_beta(pa,pd)*full_ks_beta(pa,pd)
        END IF
      END DO
    END DO
    
    IF( distribute_fock_matrix) THEN
      !! Distribute the current KS-matrix to all the processes
      CALL mp_sync(para_env%group) ! to clarify timings 
      CALL timeset(routineN//"_dist_KS",handle_dist_ks)
      CALL distribute_ks_matrix(para_env, full_ks, ks_matrix(1)%matrix, shm_number_of_p_entries, &
                                shm_is_assoc_atomic_block,&
                                natom, last_sgf_global, kind_of, basis_parameter, error)

      NULLIFY(full_ks)
      DEALLOCATE(actual_x_data%full_ks_alpha, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      IF(nspins == 2) THEN
        CALL distribute_ks_matrix(para_env, full_ks_beta, ks_matrix(2)%matrix, shm_number_of_p_entries, &
                                  shm_is_assoc_atomic_block,&
                                  natom, last_sgf_global, kind_of, basis_parameter, error)
        NULLIFY(full_ks_beta)
        DEALLOCATE(actual_x_data%full_ks_beta, STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END IF
      CALL mp_sync(para_env%group) ! to clarify timings 
      CALL timestop(handle_dist_ks)
    END IF

    CALL mp_sum(ene_x_aa, para_env%group)
    IF(nspins==2) CALL mp_sum(ene_x_bb,para_env%group)

    !! Update energy type
    IF( distribute_fock_matrix) THEN
      energy%ex = 0.5_dp*(ene_x_aa+ene_x_bb) + energy%ex
    ELSE
      energy%ex = 0.5_dp*(ene_x_aa+ene_x_bb)
    END IF

    CALL timestop(handle_main)

    !! Print some memeory information if this is the first step
    IF( geometry_did_change ) THEN
      CALL mp_sum_i8(shm_storage_counter_integrals, para_env%group)
      CALL mp_sum_i8(shm_neris_onthefly, para_env%group)
      CALL mp_sum_i8(shm_neris_incore, para_env%group)
      CALL mp_sum_i8(shm_neris_disk, para_env%group)
      CALL mp_sum_i8(shm_neris_total, para_env%group)
      CALL mp_sum_i8(shm_stor_count_int_disk, para_env%group)
      CALL mp_sum_i8(shm_nprim_ints, para_env%group)
      mem_eris = (shm_storage_counter_integrals+128*1024-1)/1024/128
      compression_factor = REAL(shm_neris_incore,dp)/REAL(shm_storage_counter_integrals,dp)
      mem_eris_disk = (shm_stor_count_int_disk+128*1024-1)/1024/128
      compression_factor_disk = REAL(shm_neris_disk,dp)/REAL(shm_stor_count_int_disk,dp)
     
      IF( shm_neris_incore == 0 ) THEN 
        mem_eris = 0
        compression_factor = 0.0_dp
      END IF
      IF( shm_neris_disk == 0 ) THEN
        mem_eris_disk = 0
        compression_factor_disk = 0.0_dp
      END IF
        
      iw = cp_print_key_unit_nr(logger,hfx_section,"HF_INFO",&
           extension=".scfLog",error=error)
      IF (iw>0) THEN
        WRITE (UNIT=iw,FMT="(/,(T3,A,T60,I20))")&
                 "HFX_MEM_INFO| Number of cart. primitive ERI's calculated:      ",shm_nprim_ints

        WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
                 "HFX_MEM_INFO| Number of sph. ERI's calculated:           ",shm_neris_total+shm_neris_onthefly

        WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
                 "HFX_MEM_INFO| Number of sph. ERI's stored in-core:        ",shm_neris_incore

        WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
                 "HFX_MEM_INFO| Number of sph. ERI's stored on disk:        ",shm_neris_disk
  
        WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
                 "HFX_MEM_INFO| Number of sph. ERI's calculated on the fly: ",shm_neris_onthefly
  
        WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
                 "HFX_MEM_INFO| Total memory consumption ERI's RAM [MB's]:            ",mem_eris
      
        WRITE (UNIT=iw,FMT="((T3,A,T60,F20.2))")&
                 "HFX_MEM_INFO| Total compression factor ERI's RAM:                  ",compression_factor
             
        WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
                 "HFX_MEM_INFO| Total memory consumption ERI's disk [MB's]:       ",mem_eris_disk

        WRITE (UNIT=iw,FMT="((T3,A,T60,F20.2))")&
                 "HFX_MEM_INFO| Total compression factor ERI's disk:             ",compression_factor_disk

        WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
                 "HFX_MEM_INFO| Size of density matrix [MB's]:             ", mb_size_p

        WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
                 "HFX_MEM_INFO| Size of fock matrix [MB's]:             ", mb_size_f
       
        IF( do_periodic) THEN
          WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
                 "HFX_MEM_INFO| Size of buffers [MB's]:             ", mb_size_buffers
          WRITE (UNIT=iw,FMT="((T3,A,T60,I20),/)")&
                 "HFX_MEM_INFO| Number of periodic image cells considered: ", SIZE(actual_x_data%neighbor_cells)
        ELSE
           WRITE (UNIT=iw,FMT="((T3,A,T60,I20),/)")&
                 "HFX_MEM_INFO| Size of buffers [MB's]:             ", mb_size_buffers
        END IF

      END IF
  
      CALL cp_print_key_finished_output(iw,logger,hfx_section,&
                                      "HF_INFO", error=error)
    END IF
!$OMP END MASTER

    !! flush caches if the geometry changed 
    IF(geometry_did_change ) THEN
      CALL hfx_flush_last_cache(6, maxval_cache, maxval_container, memory_parameter%actual_memory_usage, .FALSE., error)
      DO i=1,64
        CALL hfx_flush_last_cache(i, integral_caches(i), integral_containers(i), &
                                   memory_parameter%actual_memory_usage, .FALSE., error)
      END DO
    END IF
    !! reset all caches except we calculate all on the fly
    IF(.NOT.memory_parameter%do_all_on_the_fly) THEN
      CALL hfx_reset_cache_and_container(maxval_cache, maxval_container, memory_parameter%actual_memory_usage, .FALSE.)
      DO i=1,64
        CALL hfx_reset_cache_and_container(integral_caches(i), integral_containers(i), memory_parameter%actual_memory_usage, &
                                            .FALSE.)
      END DO
    END IF

    !! Since the I/O routines are no thread-safe, i.e. the procedure to get the unit number, put a lock here
!$OMP CRITICAL
    IF(do_disk_storage ) THEN
      !! flush caches if the geometry changed 
      IF(geometry_did_change ) THEN
        CALL hfx_flush_last_cache(6, maxval_cache_disk, maxval_container_disk, &
                                   memory_parameter%actual_memory_usage_disk, .TRUE., error)
        DO i=1,64
          CALL hfx_flush_last_cache(i, integral_caches_disk(i), integral_containers_disk(i), &
                                     memory_parameter%actual_memory_usage_disk, .TRUE., error)
        END DO
      END IF
      !! reset all caches except we calculate all on the fly
      CALL hfx_reset_cache_and_container(maxval_cache_disk, maxval_container_disk, memory_parameter%actual_memory_usage_disk, &
                                          do_disk_storage)
      DO i=1,64
        CALL hfx_reset_cache_and_container(integral_caches_disk(i), integral_containers_disk(i),&
                                            memory_parameter%actual_memory_usage_disk, do_disk_storage)
      END DO
    END IF   
!$OMP END CRITICAL
!$OMP BARRIER
    !! Clean up
    DEALLOCATE(last_sgf_global,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP MASTER
    IF( do_periodic ) THEN
      DEALLOCATE(screen_periodic, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(screen_set, screen_atom, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(full_density,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(do_correct_p_screening) THEN
      DEALLOCATE(actual_pmax_set,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(nspins==2) THEN
      DEALLOCATE(full_density_beta,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
!$OMP END MASTER
    DEALLOCATE(max_contraction, kind_of, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !
    ! this wraps to free_libint, but currently causes a segfault
    ! as a result, we don't call it, but some memory remains allocated
    !  
    ! CALL terminate_libint(lib(i_thread))
    
! many timings somewhat more clear
!$OMP BARRIER
!$OMP MASTER
    CALL mp_sync(para_env%group)
!$OMP END MASTER
!$OMP END PARALLEL
    DEALLOCATE(lib, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL mp_sync(para_env%group)
    CALL timestop(handle)
  END SUBROUTINE integrate_four_center

! *****************************************************************************
!> \brief calculates two-electron integrals of a quartet/shell using the library 
!>      lib_int
!> \param r position
!> \param z zeta
!> \param l angular momentum
!> \param npgf number of primitive cartesian gaussian in actual shell
!> \param primitive_integrals array of primitive_integrals
!> \param potential_parameter contains info for libint 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2006 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE coulomb4(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                      la_min, la_max, lb_min, lb_max,&
                      lc_min, lc_max, ld_min, ld_max,&
                      npgfa, npgfb, npgfc, npgfd,&
                      ncoa,ncob,ncoc,ncod,&
                      primitive_integrals,&
                      potential_parameter, eps_schwarz, max_contraction_val, &
                      cartesian_estimate, i_thread, neris_tmp)
  
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    INTEGER, INTENT(IN)                      :: npgfa, ld_max, ld_min, &
                                                lc_max, lc_min, lb_max, &
                                                lb_min, la_max, la_min, &
                                                npgfd, npgfc, npgfb
    REAL(dp), DIMENSION(1:npgfa), INTENT(IN) :: zeta
    REAL(dp), DIMENSION(1:npgfb), INTENT(IN) :: zetb
    REAL(dp), DIMENSION(1:npgfc), INTENT(IN) :: zetc
    REAL(dp), DIMENSION(1:npgfd), INTENT(IN) :: zetd
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncoc, ncod
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod)      :: primitive_integrals
    TYPE(hfx_potential_type)                 :: potential_parameter
    REAL(dp), INTENT(IN)                     :: eps_schwarz, &
                                                max_contraction_val
    REAL(dp), INTENT(INOUT)                  :: cartesian_estimate
    INTEGER, INTENT(IN)                      :: i_thread
    INTEGER(int_8)                           :: neris_tmp

    INTEGER :: ipgf, jpgf, kpgf, la, lb, lc, ld, lpgf, offset_a, offset_a1, &
      offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1
    REAL(dp)                                 :: tmp_max

    neris_tmp = 0
    primitive_integrals = 0.0_dp
    cartesian_estimate = 0.0_dp
    DO lpgf = 1,npgfd
      offset_d1 = (lpgf-1)*ncoset(ld_max)
      DO ld = ld_min, ld_max
        offset_d = offset_d1 + ncoset(ld-1)
        DO kpgf = 1,npgfc
          offset_c1 = (kpgf-1)*ncoset(lc_max)
          DO lc = lc_min, lc_max
            offset_c = offset_c1 + ncoset(lc-1)
            DO jpgf = 1,npgfb
              offset_b1 = (jpgf-1)*ncoset(lb_max)
              DO lb = lb_min, lb_max
                offset_b = offset_b1 + ncoset(lb-1)
                DO ipgf = 1,npgfa
                  offset_a1 = (ipgf-1)*ncoset(la_max)
                  DO la = la_min,la_max
                    offset_a = offset_a1 + ncoset(la-1)
                    !Build primitives
                    tmp_max = 0.0_dp
                    CALL evaluate_eri(lib(i_thread+1), ra, rb, rc, rd,&
                                          zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                          la, lb, lc ,ld,&
                                          ncoa,ncob,ncoc,ncod,&
                                          offset_a,offset_b,offset_c,offset_d, &
                                          primitive_integrals, &
                                          potential_parameter, max_contraction_val, tmp_max, eps_schwarz,&
                                          neris_tmp)
                    cartesian_estimate = MAX(cartesian_estimate, tmp_max)                                     
                  END DO !la
                END DO !ipgf
              END DO !lb
            END DO !jpgf
          END DO !lc
        END DO !kpgf
      END DO !ld
    END DO !lpgf
  END SUBROUTINE coulomb4

! *****************************************************************************
!> \brief calculates two-electron integrals of a quartet/shell using the library 
!>      lib_int in the case of periodic boundary conditions
!> \param r position
!> \param z zeta
!> \param l angular momentum
!> \param npgf number of primitive cartesian gaussian in actual shell
!> \param primitive_integrals array of primitive_integrals
!> \param potential_parameter contains info for libint 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      The minimum image convention (mic) is applied in the following way:
!>        1. Fix A and shift B according to mic
!>        2. Fix C and shift D according to mic
!>        3. Calculate interaction centers P and Q
!>        4. Fix P (i.e. A and B) and shift Q (i.e. C and D) according to mic
!> \par History
!>      06.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE coulomb4_per(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                          la_min, la_max, lb_min, lb_max,&
                          lc_min, lc_max, ld_min, ld_max,&
                          npgfa, npgfb, npgfc, npgfd,&
                          ncoa,ncob,ncoc,ncod,&
                          primitive_integrals,&
                          potential_parameter, neighbor_cells, &
                          screen1, screen2, eps_schwarz, max_contraction_val, &
                          cart_estimate, i_thread, cell, neris_tmp, error)

    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    INTEGER, INTENT(IN)                      :: npgfa, ld_max, ld_min, &
                                                lc_max, lc_min, lb_max, &
                                                lb_min, la_max, la_min, &
                                                npgfd, npgfc, npgfb
    REAL(dp), DIMENSION(1:npgfa), INTENT(IN) :: zeta
    REAL(dp), DIMENSION(1:npgfb), INTENT(IN) :: zetb
    REAL(dp), DIMENSION(1:npgfc), INTENT(IN) :: zetc
    REAL(dp), DIMENSION(1:npgfd), INTENT(IN) :: zetd
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncoc, ncod
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod)      :: primitive_integrals
    TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(hfx_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(dp), DIMENSION(:)                   :: screen1, screen2
    REAL(dp), INTENT(IN)                     :: eps_schwarz, &
                                                max_contraction_val
    REAL(dp), INTENT(INOUT)                  :: cart_estimate
    INTEGER, INTENT(IN)                      :: i_thread
    TYPE(cell_type), POINTER                 :: cell
    INTEGER(int_8)                           :: neris_tmp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'coulomb4_per', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, ipgf, j, jpgf, k, kpgf, la, lb, lc, ld, lpgf, offset_a, &
      offset_a1, offset_b, offset_b1, offset_c, offset_c1, offset_d, &
      offset_d1, stat
    LOGICAL                                  :: failure
    REAL(dp) :: B11(3), C11(3), D11(3), Eta, im_B(3), im_D(3), P(3), &
      pbc_B(3), pbc_D(3), PQ(3), Q(3), s(3), s1(3), s2(3), s3(3), s4(3), &
      shift(3), temp(3), tmp_D(3), tmp_max, Zeta1, Zeta_A, Zeta_B, Zeta_C, &
      Zeta_D
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: tmp_cart_estimate

    failure = .FALSE.
    cart_estimate = 0.0_dp
    neris_tmp = 0
    ALLOCATE(tmp_cart_estimate(SIZE(neighbor_cells), SIZE(neighbor_cells)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    tmp_cart_estimate = 0.0_dp
    primitive_integrals = 0.0_dp
    temp = rb -ra
    pbc_B = pbc(temp,cell)
    B11 = ra + pbc_B
    CALL real_to_scaled(s1,B11,cell)
    temp = rd-rc
    pbc_D = pbc(temp,cell)
    D11 = rc + pbc_D
    CALL real_to_scaled(s2,D11,cell)
    DO i = 1,SIZE(neighbor_cells)
      s = s1 + neighbor_cells(i)%cell
      CALL scaled_to_real(im_B,s,cell)
      DO j = 1,SIZE(neighbor_cells)
        IF(screen1(i)*screen2(j)<eps_schwarz) CYCLE
        s = s2 + neighbor_cells(j)%cell
        CALL scaled_to_real(im_D,s,cell)
        DO lpgf = 1,npgfd
          offset_d1 = (lpgf-1)*ncoset(ld_max)
          Zeta_D = zetd(lpgf)
          DO ld = ld_min, ld_max
            offset_d = offset_d1 + ncoset(ld-1)
            DO kpgf = 1,npgfc
              offset_c1 = (kpgf-1)*ncoset(lc_max)
              Zeta_C = zetc(kpgf)
              Eta  =  Zeta_C + Zeta_D
              Q = (Zeta_C*rc + Zeta_D*im_D)/Eta
              DO lc = lc_min, lc_max
                offset_c = offset_c1 + ncoset(lc-1)
                DO jpgf = 1,npgfb
                  offset_b1 = (jpgf-1)*ncoset(lb_max)
                  Zeta_B = zetb(jpgf)
                  DO lb = lb_min, lb_max
                    offset_b = offset_b1 + ncoset(lb-1)
                    DO ipgf = 1,npgfa
                      offset_a1 = (ipgf-1)*ncoset(la_max)
                      Zeta_A = zeta(ipgf)
                      Zeta1 = Zeta_A + Zeta_B
                      P = (Zeta_A*ra + Zeta_B*im_B)/Zeta1
                      temp = P-Q
                      PQ = pbc(temp,cell)
                      shift = -PQ + temp
                      temp = rc+shift
                      CALL real_to_scaled(s3,temp,cell)
                      temp = im_D+shift
                      CALL real_to_scaled(s4,temp,cell)
                      DO k = 1,SIZE(neighbor_cells)
                        s = s3 + neighbor_cells(k)%cell
                        CALL scaled_to_real(C11,s,cell)
                        s = s4 + neighbor_cells(k)%cell
                        CALL scaled_to_real(tmp_D,s,cell)
                        DO la = la_min,la_max
                          offset_a = offset_a1 + ncoset(la-1)
                          !Build primitives
                          tmp_max = 0.0_dp
                          CALL evaluate_eri(lib(i_thread+1), ra, im_B, C11, tmp_D,&
                                            zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                            la, lb, lc ,ld,&
                                            ncoa,ncob,ncoc,ncod,&
                                            offset_a,offset_b,offset_c,offset_d, &
                                            primitive_integrals, &
                                            potential_parameter, max_contraction_val, tmp_max, eps_schwarz,&
                                            neris_tmp)
                          tmp_cart_estimate(i,j) = MAX(tmp_cart_estimate(i,j),tmp_max)
                        END DO !la
                      END DO !k-shell
                    END DO !ipgf
                  END DO !lb
                END DO !jpgf
              END DO !lc
            END DO !kpgf
          END DO !ld
        END DO !lpgf
      END DO !j-shell
    END DO !i-shell
    DO i=1,SIZE(neighbor_cells)
      DO j=1,SIZE(neighbor_cells)
        cart_estimate = MAX(cart_estimate, tmp_cart_estimate(i,j))
      END DO
    END DO

    DEALLOCATE(tmp_cart_estimate, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE coulomb4_per

! *****************************************************************************
!> \brief This routine does the correction to the screening on an in inital
!>      density matrix. It is basically a copy from the out-core part in
!>      integrate_four_center.
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE on_the_fly(do_periodic, ra, rb, rc, rd, zeta, zetb, zetc, zetd, &
                        la_min, la_max, lb_min, lb_max, lc_min, lc_max, ld_min, ld_max, &
                        npgfa, npgfb, npgfc, npgfd, ncoa, ncob, ncoc, ncod, &
                        nsgfa, nsgfb, nsgfc, nsgfd, sgfa, sgfb, sgfc, sgfd, &
                        primitive_integrals, potential_parameter, screen_periodic, &
                        eps_schwarz, sphi_a, sphi_b, sphi_c, sphi_d, T1,& 
                        full_ks, full_ks_beta, full_density, full_density_beta, fac, symm_fac, &
                        iset, jset, kset, lset, iatom, jatom, katom, latom, pa1, pb1, pc1, pd1, &
                        nspins, max_contraction_val, actual_x_data, i_thread, cell, neris_tmp, error)

    LOGICAL, INTENT(IN)                      :: do_periodic
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: zeta, zetb, zetc, zetd
    INTEGER, DIMENSION(:), INTENT(IN)        :: la_min, la_max, lb_min, &
                                                lb_max, lc_min, lc_max, &
                                                ld_min, ld_max, npgfa, npgfb, &
                                                npgfc, npgfd
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncoc, ncod
    INTEGER, DIMENSION(:), INTENT(IN)        :: nsgfa, nsgfb, nsgfc, nsgfd
    INTEGER, INTENT(IN)                      :: sgfa, sgfb, sgfc, sgfd
    REAL(dp), DIMENSION(:)                   :: primitive_integrals
    TYPE(hfx_potential_type)                 :: potential_parameter
    REAL(dp), DIMENSION(:, :, :, :, :), &
      POINTER                                :: screen_periodic
    REAL(dp), INTENT(IN)                     :: eps_schwarz
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: sphi_a, sphi_b, sphi_c, sphi_d
    REAL(dp), DIMENSION(:)                   :: T1
    REAL(dp), DIMENSION(:, :)                :: full_ks, full_ks_beta, &
                                                full_density, &
                                                full_density_beta
    REAL(dp)                                 :: fac, symm_fac
    INTEGER, INTENT(IN)                      :: iset, jset, kset, lset, &
                                                iatom, jatom, katom, latom, &
                                                pa1, pb1, pc1, pd1, nspins
    REAL(dp), INTENT(IN)                     :: max_contraction_val
    TYPE(hfx_type), POINTER                  :: actual_x_data
    INTEGER                                  :: i_thread
    TYPE(cell_type), POINTER                 :: cell
    INTEGER(int_8)                           :: neris_tmp
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ma, mb, mc, md, p_index, pa, &
                                                pa2, pb, pb2, pc, pc2, pd, pd2
    REAL(dp)                                 :: cartesian_estimate, tmp

    IF( .NOT.do_periodic ) THEN 
      CALL coulomb4(ra, rb, rc, rd,&
                    zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset), &
                    zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                    la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                    lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                    npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                    ncoa,ncob,ncoc,ncod,&
                    primitive_integrals,&
                    potential_parameter, eps_schwarz, max_contraction_val, &
                    cartesian_estimate, i_thread, neris_tmp)
    ELSE
      CALL coulomb4_per(ra, rb, rc, rd,&
                        zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset), &
                        zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                        la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                        lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                        npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                        ncoa,ncob,ncoc,ncod,&
                        primitive_integrals,&
                        potential_parameter, &
                        actual_x_data%neighbor_cells, screen_periodic(:,iset,jset,iatom,jatom),&
                        screen_periodic(:,kset,lset,katom,latom), eps_schwarz, &
                        max_contraction_val, cartesian_estimate, i_thread, cell, neris_tmp,&
                        error)
    END IF

    IF(cartesian_estimate < eps_schwarz) RETURN

    CALL dgemm("T","N",ncob*ncoc*ncod,nsgfa(iset),ncoa,&
               1.0_dp, primitive_integrals(1),ncoa,&
               sphi_a(1,sgfa), SIZE(sphi_a,1),&
               0.0_dp, T1(1),ncob*ncoc*ncod)
  
    CALL dgemm("T","N",nsgfa(iset)*ncoc*ncod,nsgfb(jset),ncob,&
               1.0_dp, T1(1),ncob,&
               sphi_b(1,sgfb), SIZE(sphi_b,1),&
               0.0_dp, primitive_integrals(1),nsgfa(iset)*ncoc*ncod)
  
    CALL dgemm("T","N",nsgfa(iset)*nsgfb(jset)*ncod,nsgfc(kset),ncoc,&
               1.0_dp, primitive_integrals(1),ncoc,&
               sphi_c(1,sgfc), SIZE(sphi_c,1),&
               0.0_dp, T1(1),nsgfa(iset)*nsgfb(jset)*ncod)
  
    CALL dgemm("T","N",nsgfa(iset)*nsgfb(jset)*nsgfc(kset),nsgfd(lset),ncod,&
               1.0_dp, T1(1),ncod,&
               sphi_d(1,sgfd), SIZE(sphi_d,1),&
               0.0_dp, primitive_integrals(1),nsgfa(iset)*nsgfb(jset)*nsgfc(kset))

    p_index=0 
    DO md = 1,nsgfd(lset)
      pd2 = md+sgfd-1
      pd = pd1 + pd2
      DO mc = 1,nsgfc(kset)
        pc2 = mc+sgfc-1
        pc = pc1 + pc2
        DO mb = 1,nsgfb(jset)
          pb2 = mb+sgfb-1
          pb = pb1 + pb2
          DO ma = 1,nsgfa(iset)
            pa2 = ma+sgfa-1
            pa = pa1 + pa2
            p_index=p_index+1
            tmp = fac*symm_fac*primitive_integrals(p_index)
            full_ks(pa,pd) = full_ks(pa,pd)-tmp*full_density(pb,pc)
            full_ks(pa,pc) = full_ks(pa,pc)-tmp*full_density(pb,pd)
            full_ks(pb,pc) = full_ks(pb,pc)-tmp*full_density(pa,pd)
            full_ks(pb,pd) = full_ks(pb,pd)-tmp*full_density(pa,pc)
            IF( nspins == 2 ) THEN
              full_ks_beta(pa,pd) = full_ks_beta(pa,pd)-tmp*full_density_beta(pb,pc)
              full_ks_beta(pa,pc) = full_ks_beta(pa,pc)-tmp*full_density_beta(pb,pd)
              full_ks_beta(pb,pc) = full_ks_beta(pb,pc)-tmp*full_density_beta(pa,pd)
              full_ks_beta(pb,pd) = full_ks_beta(pb,pd)-tmp*full_density_beta(pa,pc)
            END IF
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE on_the_fly

END MODULE hfx_energy_potential
