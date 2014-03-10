!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines to calculate derivatives with respect to basis function origin
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE hfx_derivatives

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_matrix_type
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_files,                        ONLY: get_unit_number
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_types,                     ONLY: dbcsr_type_antisymmetric
  USE gamma,                           ONLY: init_md_ftable
  USE hfx_communication,               ONLY: get_full_density
  USE hfx_compression_methods,         ONLY: hfx_add_mult_cache_elements,&
                                             hfx_add_single_cache_element,&
                                             hfx_decompress_first_cache,&
                                             hfx_flush_last_cache,&
                                             hfx_get_mult_cache_elements,&
                                             hfx_get_single_cache_element,&
                                             hfx_reset_cache_and_container
  USE hfx_libint_interface,            ONLY: evaluate_deriv_eri
  USE hfx_libint_wrapper_types,        ONLY: lib_deriv
  USE hfx_load_balance_methods,        ONLY: collect_load_balance_info,&
                                             hfx_load_balance,&
                                             hfx_update_load_balance
  USE hfx_pair_list_methods,           ONLY: build_atomic_pair_list,&
                                             build_pair_list,&
                                             build_pair_list_pgf,&
                                             build_pgf_product_list
  USE hfx_screening_methods,           ONLY: update_pmax_mat
  USE hfx_types,                       ONLY: &
       alloc_containers, dealloc_containers, hfx_basis_info_type, &
       hfx_basis_type, hfx_cache_type, hfx_cell_type, hfx_container_type, &
       hfx_distribution, hfx_general_type, hfx_init_container, &
       hfx_load_balance_type, hfx_memory_type, hfx_p_kind, hfx_pgf_list, &
       hfx_pgf_product_list, hfx_potential_type, hfx_screen_coeff_type, &
       hfx_screening_type, hfx_task_list_type, hfx_type, init_t_c_g0_lmax, &
       log_zero, pair_list_type, pair_set_list_type
  USE input_constants,                 ONLY: do_hfx_potential_mix_cl_trunc,&
                                             do_hfx_potential_truncated,&
                                             hfx_do_eval_forces
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: m_walltime
  USE mathconstants
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE t_c_g0,                          ONLY: init
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  PUBLIC  derivatives_four_center, forces4

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_derivatives'



!***

  CONTAINS

! *****************************************************************************
!> \brief computes four center derivatives for a full basis set and updates the
!>      forces%fock_4c arrays. Uses all 8 eri symmetries
!> \param qs_env
!> \param ks_matrix
!> \param energy
!> \param rho density matrix
!> \param hfx_section HFX input section
!> \param para_env para_env
!> \param irep ID of HFX replica
!> \param adiabatic_rescale_factor parameter used for MCY3 hybrid
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2007 created [Manuel Guidon]
!>      08.2007 optimized load balance [Manuel Guidon]
!>      02.2009 completely rewritten screening part [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE derivatives_four_center(qs_env,rho_ao, hfx_section, para_env, &
                                     irep, use_virial, adiabatic_rescale_factor, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao
    TYPE(section_vals_type), POINTER         :: hfx_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: irep
    LOGICAL                                  :: use_virial
    REAL(dp), INTENT(IN), OPTIONAL           :: adiabatic_rescale_factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derivatives_four_center', &
      routineP = moduleN//':'//routineN

    INTEGER :: bin, coord, forces_map(4,2), handle, handle_bin, handle_getP, &
      handle_load, handle_main, i, i_atom, i_thread, iatom, ikind, iset, &
      j_atom, jatom, jkind, jset, k_atom, katom, kkind, kset, l_atom, latom, &
      lkind, lset, max_am, max_set, my_bin_size, n_threads, natom, ncob, &
      ncos_max, nseta, nsetb, nspins, sgfb, stat
    INTEGER(int_8)                           :: my_current_counter, &
                                                my_istart, n_processes, ncpu
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                last_sgf_global
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, npgfa, npgfb, npgfc, npgfd, nsgfa, nsgfb, &
      nsgfc, nsgfd
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: do_periodic, failure, &
                                                screen_pmat_forces
    REAL(dp) :: bintime_start, bintime_stop, cartesian_estimate, &
      cartesian_estimate_virial, eps_schwarz, fac, hf_fraction, &
      max_contraction_val, max_val1, max_val2, my_adiabatic_rescale_factor, &
      pmax_entry, ra(3), rb(3), rc(3), rd(3), spherical_estimate_virial, &
      symm_fac
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: work_forces
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      TARGET                                 :: primitive_forces
    REAL(dp), DIMENSION(:), POINTER          :: full_density, &
                                                full_density_beta, &
                                                full_density_mp2, T2
    REAL(dp), DIMENSION(:, :), POINTER       :: max_contraction, sphi_b, &
                                                zeta, zetb, zetc, zetd
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(hfx_basis_info_type), POINTER       :: basis_info
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_general_type)                   :: general_parameter
    TYPE(hfx_load_balance_type), POINTER     :: load_balance_parameter
    TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(hfx_screening_type)                 :: screening_parameter
    TYPE(hfx_type), POINTER                  :: actual_x_data
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force

!$  INTEGER :: omp_get_max_threads,omp_get_thread_num
    INTEGER, SAVE                            :: shm_number_of_p_entries
    INTEGER, DIMENSION(:,:), POINTER, SAVE   :: shm_is_assoc_atomic_block
    INTEGER                                  :: nkind, l_max
    TYPE(cell_type), POINTER                 :: cell
    INTEGER                                  :: unit_id
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:,:,:,:,:,:), POINTER        :: screen_coeffs_pgf, radii_pgf
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:,:,:,:), POINTER            :: screen_coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:,:), POINTER                :: screen_coeffs_kind, tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2
    REAL(dp)                                 :: ln_10, rab2, rcd2, log10_eps_schwarz,&
                                                log10_pmax, max_val2_set, log_2
    INTEGER                                  :: nsgf_max
    REAL(dp), DIMENSION(:), ALLOCATABLE      :: pbd_buf, pbc_buf, pad_buf, pac_buf, &
                                                pbd_buf_mp2, pbc_buf_mp2, pad_buf_mp2, pac_buf_mp2, &
                                                pbd_buf_beta, pbc_buf_beta, pad_buf_beta, pac_buf_beta
    INTEGER(int_8)                           :: atom_block, tmp_block, nblocks
    INTEGER                                  :: iatom_block, jatom_block, katom_block, latom_block
    INTEGER                                  :: i_list_ij,i_list_kl, &
                                                i_set_list_ij, i_set_list_ij_start, i_set_list_ij_stop,&
                                                i_set_list_kl, i_set_list_kl_start, i_set_list_kl_stop
    TYPE(pair_list_type)                     :: list_ij, list_kl
    REAL(KIND=dp)                            :: coeffs_kind_max0
    TYPE(pair_set_list_type), DIMENSION(:), ALLOCATABLE :: set_list_ij, set_list_kl
    TYPE(hfx_p_kind), DIMENSION(:), POINTER  :: shm_initial_p
    REAL(dp)                                 :: pmax_tmp
    INTEGER                                  :: kind_kind_idx, swap_id
    REAL(dp), DIMENSION(:,:), POINTER        :: ptr_p_1, ptr_p_2, ptr_p_3, ptr_p_4, shm_pmax_block
    INTEGER, DIMENSION(:), POINTER           :: shm_block_offset
    INTEGER, DIMENSION(:,:), POINTER         :: shm_atomic_block_offset
    INTEGER, DIMENSION(:,:,:,:), POINTER     :: shm_set_offset
    INTEGER                                  :: atomic_offset_bd, atomic_offset_bc, atomic_offset_ad, &
                                                atomic_offset_ac, nneighbors, max_pgf, &
                                                atomic_offset_ab, atomic_offset_cd
    INTEGER, DIMENSION(:,:), POINTER         :: offset_bd_set, offset_bc_set, offset_ad_set, offset_ac_set
    INTEGER, DIMENSION(:,:), POINTER         :: offset_ab_set, offset_cd_set
    TYPE(hfx_pgf_list), DIMENSION(:), &
                            ALLOCATABLE      :: pgf_list_ij, pgf_list_kl
    TYPE(hfx_pgf_product_list), DIMENSION(:), &
                        ALLOCATABLE          :: pgf_product_list
    REAL(dp), DIMENSION(:,:,:,:), POINTER &
                                             :: sphi_a_ext, sphi_b_ext, sphi_c_ext, sphi_d_ext
    REAL(dp), DIMENSION(:,:,:), POINTER &
                                             :: sphi_a_ext_set, sphi_b_ext_set, sphi_c_ext_set, sphi_d_ext_set
    INTEGER, DIMENSION(:, :), POINTER        :: nsgfl_a, nsgfl_b, &
                                                nsgfl_c, nsgfl_d
    REAL(dp), DIMENSION(:), ALLOCATABLE      :: ede_work, ede_work2, ede_work_forces,  &
                                                ede_buffer1, ede_buffer2, &
                                                ede_primitives_tmp
    TYPE(lib_deriv)                          :: private_deriv
    LOGICAL                                  :: do_dynamic_load_balancing, bins_left, do_it
    TYPE(hfx_distribution), POINTER          :: distribution_forces
    INTEGER                                  :: shm_total_bins, shm_task_counter, my_thread_id, my_bin_id ,&
                                                storage_id
    TYPE(hfx_task_list_type), DIMENSION(:), &
                              POINTER        :: shm_task_list, tmp_task_list
    INTEGER(int_8), DIMENSION(:), &
                    ALLOCATABLE              :: tmp_task_list_cost
    INTEGER, DIMENSION(:), ALLOCATABLE       :: tmp_index
    INTEGER(int_8)                           :: current_ram_counter
    INTEGER, DIMENSION(:), ALLOCATABLE       :: nimages
    INTEGER                                  :: sphi_d_u1,sphi_d_u2,sphi_d_u3
    INTEGER                                  :: sphi_c_u1,sphi_c_u2,sphi_c_u3
    INTEGER                                  :: sphi_b_u1,sphi_b_u2,sphi_b_u3
    INTEGER                                  :: sphi_a_u1,sphi_a_u2,sphi_a_u3
    INTEGER                                  :: iatom_start, iatom_end, jatom_start, jatom_end, katom_start,&
                                                katom_end, latom_start, latom_end, my_process_id
    TYPE(hfx_cache_type), DIMENSION(:), &
      POINTER                                :: integral_caches
    TYPE(hfx_cache_type), POINTER            :: maxval_cache
    TYPE(hfx_container_type), DIMENSION(:), &
      POINTER                                :: integral_containers
    TYPE(hfx_container_type), POINTER        :: maxval_container
    INTEGER                                  :: cache_size, bits_max_val, current_counter, nints, &
                                                nbits, buffer_left, buffer_start, buffer_size, iw
    TYPE(hfx_memory_type), POINTER           :: memory_parameter
    INTEGER(int_8)                           :: mem_compression_counter, counter, neris_onthefly, &
                                                estimate_to_store_int, neris_total, nprim_ints, neris_tmp,&
                                                max_val_memory, neris_incore
    LOGICAL                                  :: buffer_overflow, use_disk_storage
    REAL(dp)                                 :: spherical_estimate, eps_storage
    INTEGER(int_8)                           :: shm_neris_total, shm_neris_onthefly, &
                                                shm_storage_counter_integrals, shm_neris_incore, &
                                                shm_neris_disk, shm_stor_count_int_disk, shm_stor_count_max_val,&
                                                tmp_i8(6), mem_eris, mem_max_val, shm_nprim_ints, stor_count_max_val,&
                                                storage_counter_integrals
    REAL(dp)                                 :: compression_factor, cartesian_estimate_orig
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL                                  :: treat_forces_in_core
    REAL(dp), DIMENSION(:,:), POINTER        :: shm_pmax_atom

    TYPE(cp_dbcsr_p_type), POINTER, DIMENSION(:) :: rho_ao_mp2
    INTEGER::ii
    REAL(dp), DIMENSION(:), ALLOCATABLE      :: ede_work_virial, ede_work2_virial,&
                                                ede_primitives_tmp_virial
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      TARGET                                 :: primitive_forces_virial
    TYPE(virial_type), POINTER               :: virial
    REAL(dp)                                 :: tmp_virial(3,3), pmax_blocks, pmax_atom
    INTEGER                                  :: k,j
    LOGICAL, DIMENSION(:,:), POINTER         :: shm_atomic_pair_list
    LOGICAL                                  :: do_print_load_balance_info
    LOGICAL                                  :: with_mp2_density, is_anti_symmetric

    with_mp2_density=.FALSE.
    IF(ASSOCIATED(qs_env%mp2_env)) with_mp2_density=.TRUE.

    IF(with_mp2_density) rho_ao_mp2=>qs_env%matrix_p_mp2
    is_anti_symmetric=cp_dbcsr_get_matrix_type(rho_ao(1)%matrix).EQ.dbcsr_type_antisymmetric

    CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, atomic_kind_set=atomic_kind_set, &
                    cell=cell, error=error)
    !! One atom systems have no contribution to forces
    IF( SIZE(particle_set,1) == 1) THEN
      IF( .NOT. use_virial ) RETURN
    END IF

    CALL timeset(routineN,handle)
    !! Calculate l_max used in fgamma , because init_md_ftable is definitely not thread safe
    nkind = SIZE(atomic_kind_set,1)
    l_max = 0
    DO ikind=1,nkind
      l_max = MAX(l_max,MAXVAL(qs_env%x_data(1,1)%basis_parameter(ikind)%lmax))
    ENDDO
    l_max = 4*l_max + 1
    CALL init_md_ftable(l_max)

    IF(qs_env%x_data(1,1)%potential_parameter%potential_type == do_hfx_potential_truncated .OR. &
       qs_env%x_data(1,1)%potential_parameter%potential_type == do_hfx_potential_mix_cl_trunc ) THEN
      IF(l_max>init_t_c_g0_lmax) THEN
        unit_id = get_unit_number()
        OPEN(unit_id,FILE=qs_env%x_data(1,1)%potential_parameter%filename)
        CALL init(l_max, unit_id, para_env%mepos, para_env%group)
        CLOSE(unit_id)
        init_t_c_g0_lmax = l_max
      END IF
    END IF


    n_threads = 1
!$  n_threads = omp_get_max_threads()

    shm_neris_total = 0
    shm_nprim_ints = 0
    shm_neris_onthefly = 0
    shm_storage_counter_integrals = 0
    shm_neris_incore = 0
    shm_stor_count_max_val = 0


    !! get force array
    CALL get_qs_env(qs_env=qs_env, force=force, virial=virial, error=error)

    my_adiabatic_rescale_factor = 1.0_dp
    IF(PRESENT(adiabatic_rescale_factor)) THEN
      my_adiabatic_rescale_factor = adiabatic_rescale_factor
    END IF


!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(qs_env,&
!$OMP                                  rho_ao,&
!$OMP                                  rho_ao_mp2,&
!$OMP                                  with_mp2_density,&
!$OMP                                  hfx_section,&
!$OMP                                  para_env,&
!$OMP                                  irep,&
!$OMP                                  error,&
!$OMP                                  ncoset,&
!$OMP                                  nco,&
!$OMP                                  nso,&
!$OMP                                  n_threads,&
!$OMP                                  full_density,&
!$OMP                                  full_density_mp2,&
!$OMP                                  full_density_beta,&
!$OMP                                  shm_initial_p,&
!$OMP                                  shm_is_assoc_atomic_block,&
!$OMP                                  shm_number_of_p_entries,&
!$OMP                                  force,&
!$OMP                                  virial,&
!$OMP                                  my_adiabatic_rescale_factor,&
!$OMP                                  shm_neris_total,&
!$OMP                                  shm_neris_onthefly,&
!$OMP                                  shm_storage_counter_integrals,&
!$OMP                                  shm_neris_incore,&
!$OMP                                  shm_nprim_ints,&
!$OMP                                  shm_stor_count_max_val,&
!$OMP                                  cell,&
!$OMP                                  screen_coeffs_set,&
!$OMP                                  screen_coeffs_kind,&
!$OMP                                  screen_coeffs_pgf,&
!$OMP                                  nkind,&
!$OMP                                  radii_pgf,&
!$OMP                                  shm_atomic_block_offset,&
!$OMP                                  shm_set_offset,&
!$OMP                                  shm_block_offset,&
!$OMP                                  shm_task_counter,&
!$OMP                                  shm_task_list,&
!$OMP                                  shm_total_bins,&
!$OMP                                  shm_pmax_block,&
!$OMP                                  shm_atomic_pair_list,&
!$OMP                                  shm_pmax_atom,&
!$OMP                                  use_virial,&
!$OMP                                  do_print_load_balance_info)

    i_thread = 0
!$  i_thread = omp_get_thread_num()

    ln_10 = LOG(10.0_dp)
    log_2 = LOG10(2.0_dp)
    actual_x_data => qs_env%x_data(irep, i_thread + 1)
    do_periodic = actual_x_data%periodic_parameter%do_periodic

    screening_parameter    = actual_x_data%screening_parameter
    general_parameter      = actual_x_data%general_parameter
    potential_parameter    = actual_x_data%potential_parameter
    basis_info             => actual_x_data%basis_info

    load_balance_parameter => actual_x_data%load_balance_parameter
    basis_parameter        => actual_x_data%basis_parameter

    ! IF(with_mp2_density) THEN
    !   ! here we also do a copy of the original load_balance_parameter
    !   ! since if MP2 forces are required after the calculation of the HFX
    !   ! derivatives we need to calculate again the SCF energy.
    !   ALLOCATE(load_balance_parameter_energy)
    !   load_balance_parameter_energy = actual_x_data%load_balance_parameter
    ! END IF

    memory_parameter       => actual_x_data%memory_parameter

    cache_size = memory_parameter%cache_size
    bits_max_val = memory_parameter%bits_max_val

    use_disk_storage = .FALSE.

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    error=error)

    max_set = basis_info%max_set
    max_am = basis_info%max_am
    natom = SIZE(particle_set,1)

    treat_forces_in_core = memory_parameter%treat_forces_in_core

    hf_fraction = general_parameter%fraction
    hf_fraction = hf_fraction * my_adiabatic_rescale_factor
    eps_schwarz = screening_parameter%eps_schwarz_forces
    IF( eps_schwarz < 0.0_dp) THEN
      log10_eps_schwarz = log_zero
    ELSE
      !! ** make eps_schwarz tighter in case of a stress calculation
      IF( use_virial ) eps_schwarz = eps_schwarz / actual_x_data%periodic_parameter%R_max_stress
      log10_eps_schwarz = LOG10(eps_schwarz)
    END IF
    screen_pmat_forces =  screening_parameter%do_p_screening_forces
    eps_storage = eps_schwarz*memory_parameter%eps_storage_scaling

    counter = 0_int_8

    !! Copy the libint_guy
    private_deriv = actual_x_data%lib_deriv

    !! Get screening parameter

    !! Allocate work-space for forces
    ALLOCATE(work_forces(((max_am+1)*(max_am+2)/2)**4*12),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind,&
                             kind_of=kind_of)

    nspins = qs_env%dft_control%nspins

    !! Create helper arrray for mapping local basis functions to global ones
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO iatom=1,natom
      ikind = kind_of(iatom)
      last_sgf_global(iatom) = last_sgf_global(iatom-1) + basis_parameter(ikind)%nsgf_total
    END DO

    ALLOCATE(max_contraction(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    max_contraction=0.0_dp
    max_pgf = 0
    DO jatom=1,natom
      jkind = kind_of(jatom)
      lb_max => basis_parameter(jkind)%lmax
      npgfb => basis_parameter(jkind)%npgf
      nsetb = basis_parameter(jkind)%nset
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
        max_pgf = MAX(max_pgf,npgfb(jset))
      ENDDO
    ENDDO

    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads

    !! initalize some counters
    neris_total = 0_int_8
    neris_incore = 0_int_8
    neris_onthefly = 0_int_8
    mem_eris = 0_int_8
    mem_max_val = 0_int_8
    compression_factor = 0.0_dp
    nprim_ints = 0_int_8
    neris_tmp = 0_int_8
    max_val_memory = 1_int_8



!$OMP MASTER
    !! Set pointer for is_assoc helper
    shm_is_assoc_atomic_block => actual_x_data%is_assoc_atomic_block
    shm_number_of_p_entries = actual_x_data%number_of_p_entries
    shm_atomic_block_offset => actual_x_data%atomic_block_offset
    shm_set_offset => actual_x_data%set_offset
    shm_block_offset => actual_x_data%block_offset

    NULLIFY(full_density)
    NULLIFY(full_density_beta)
    ALLOCATE(full_density(shm_block_offset(ncpu+1)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    !! Get the full density from all the processors
    CALL timeset(routineN//"_getP",handle_getP)
    CALL get_full_density(para_env, full_density, rho_ao(1)%matrix, shm_number_of_p_entries,&
                          shm_block_offset, natom, &
                          kind_of, basis_parameter, get_max_vals_spin=.FALSE.,&
                          antisymmetric=is_anti_symmetric,error=error)
    ! for now only closed shell case
    IF(with_mp2_density) THEN
      NULLIFY(full_density_mp2)
      ALLOCATE(full_density_mp2(shm_block_offset(ncpu+1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL get_full_density(para_env, full_density_mp2, rho_ao_mp2(1)%matrix, shm_number_of_p_entries,&
                            shm_block_offset, natom, &
                            kind_of, basis_parameter, get_max_vals_spin=.FALSE.,&
                            antisymmetric=is_anti_symmetric,error=error)
    END IF

    IF(nspins == 2) THEN
      ALLOCATE(full_density_beta(shm_block_offset(ncpu+1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL get_full_density(para_env, full_density_beta, rho_ao(2)%matrix, shm_number_of_p_entries,&
                            shm_block_offset, natom, &
                            kind_of, basis_parameter, get_max_vals_spin=.FALSE.,&
                            antisymmetric=is_anti_symmetric,error=error)
    END IF
    CALL timestop(handle_getP)

    !! Calculate max entries for screening on actual density. If screen_p_mat_forces = FALSE, the
    !! matrix is initialized to 1.0
    IF ( screen_pmat_forces ) THEN
      NULLIFY(shm_initial_p)
      shm_initial_p => actual_x_data%initial_p
      shm_pmax_atom => actual_x_data%pmax_atom

      IF( memory_parameter%treat_forces_in_core) THEN
        shm_initial_p => actual_x_data%initial_p_forces
        shm_pmax_atom => actual_x_data%pmax_atom_forces
      END IF
      IF( memory_parameter%recalc_forces ) THEN
        CALL update_pmax_mat(shm_initial_p, &
                             actual_x_data%map_atom_to_kind_atom, &
                             actual_x_data%set_offset, &
                             actual_x_data%atomic_block_offset, &
                             shm_pmax_atom,&
                             full_density,full_density_beta,&
                             natom, kind_of, basis_parameter, &
                             nkind, actual_x_data%is_assoc_atomic_block)
      END IF
    END IF

    ! restore as full density the HF density
    ! maybe in the future
    IF(with_mp2_density) THEN
      full_density=full_density-full_density_mp2
      ! full_density_mp2=full_density+full_density_mp2
    END IF

    screen_coeffs_set   => actual_x_data%screen_funct_coeffs_set
    screen_coeffs_kind  => actual_x_data%screen_funct_coeffs_kind
    screen_coeffs_pgf   => actual_x_data%screen_funct_coeffs_pgf
    radii_pgf           => actual_x_data%pair_dist_radii_pgf
    shm_pmax_block      => actual_x_data%pmax_block

!$OMP END MASTER
!$OMP BARRIER

!$OMP MASTER
    CALL timeset(routineN//"_load",handle_load)
!$OMP END MASTER
    !! Load balance the work
    IF( actual_x_data%memory_parameter%recalc_forces ) THEN
      IF( actual_x_data%b_first_load_balance_forces ) THEN
        CALL hfx_load_balance(actual_x_data,eps_schwarz,particle_set,max_set,para_env,&
                              screen_coeffs_set,screen_coeffs_kind,&
                              shm_is_assoc_atomic_block, do_periodic,&
                              load_balance_parameter, kind_of, basis_parameter,shm_initial_p, &
                              shm_pmax_atom, i_thread, n_threads,&
                              cell, screen_pmat_forces, actual_x_data%map_atom_to_kind_atom, &
                              nkind, hfx_do_eval_forces, shm_pmax_block, use_virial, error)
        actual_x_data%b_first_load_balance_forces = .FALSE.
      ELSE
        CALL hfx_update_load_balance(actual_x_data, para_env, &
                                     load_balance_parameter, &
                                     i_thread, n_threads, hfx_do_eval_forces, &
                                     error)
      END IF
    END IF
!$OMP MASTER
    CALL timestop(handle_load)
!$OMP END MASTER

    !! precompute maximum nco and allocate scratch
    ncos_max=0
    nsgf_max=0
    DO iatom=1,natom
      ikind = kind_of(iatom)
      nseta = basis_parameter(ikind)%nset
      npgfa => basis_parameter(ikind)%npgf
      la_max => basis_parameter(ikind)%lmax
      nsgfa => basis_parameter(ikind)%nsgf
      DO iset = 1, nseta
         ncos_max = MAX(ncos_max,ncoset(la_max(iset)))
         nsgf_max = MAX(nsgf_max, nsgfa(iset))
      ENDDO
    ENDDO

    !! Allocate work arrays
    ALLOCATE(primitive_forces(12*nsgf_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    primitive_forces=0.0_dp

    ! ** Allocate buffers for pgf_lists
    nneighbors = SIZE(actual_x_data%neighbor_cells)
    ALLOCATE(pgf_list_ij(max_pgf**2), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(pgf_list_kl(max_pgf**2), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(pgf_product_list(nneighbors**3), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(nimages(max_pgf**2), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    DO i=1,max_pgf**2
      ALLOCATE(pgf_list_ij(i)%image_list(nneighbors), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(pgf_list_kl(i)%image_list(nneighbors), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO


    ALLOCATE(pbd_buf(nsgf_max**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(pbc_buf(nsgf_max**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(pad_buf(nsgf_max**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(pac_buf(nsgf_max**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(with_mp2_density) THEN
      ALLOCATE(pbd_buf_mp2(nsgf_max**2),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(pbc_buf_mp2(nsgf_max**2),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(pad_buf_mp2(nsgf_max**2),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(pac_buf_mp2(nsgf_max**2),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF( nspins == 2 ) THEN
      ALLOCATE(pbd_buf_beta(nsgf_max**2),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(pbc_buf_beta(nsgf_max**2),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(pad_buf_beta(nsgf_max**2),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(pac_buf_beta(nsgf_max**2),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    ALLOCATE(ede_work(ncos_max**4*12),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ede_work2(ncos_max**4*12),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ede_work_forces(ncos_max**4*12),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ede_buffer1(ncos_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ede_buffer2(ncos_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ede_primitives_tmp(nsgf_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF( use_virial ) THEN
      ALLOCATE(primitive_forces_virial(12*nsgf_max**4*3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      primitive_forces_virial=0.0_dp
      ALLOCATE(ede_work_virial(ncos_max**4*12*3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(ede_work2_virial(ncos_max**4*12*3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(ede_primitives_tmp_virial(nsgf_max**4),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_virial = 0.0_dp
    ELSE
      ! dummy allocation
      ALLOCATE(primitive_forces_virial(1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      primitive_forces_virial=0.0_dp
      ALLOCATE(ede_work_virial(1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(ede_work2_virial(1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(ede_primitives_tmp_virial(1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF


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
    !! Furthermore, if iatom==jatom==katom==latom we cycle, because the derivatives are zero anyway.
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
    !! istart, jstart, kstart, lstart, number_of_atom_quartets, initial_cpu_id
    !!
    !! This tells the process where to start the main loops and how many bunches of integrals it has to
    !! calculate. The original parallelization is a simple modulo distribution that is binned and
    !! optimized in the load_balance routines. Since the Monte Carlo routines can swap processors,
    !! we need to know which was the inital cpu_id.
    !! Furthermore, the indices iatom, jatom, katom, latom have to be set to istart, jstart, kstart and
    !! lstart only the first time the loop is executed. All subsequent loops have to start with one or
    !! iatom and katom respectively. Therefore, we use flags like first_j_loop etc.

!$OMP BARRIER
!$OMP MASTER
    CALL timeset(routineN//"_main",handle_main)
!$OMP END MASTER
!$OMP BARRIER


    coeffs_kind_max0=MAXVAL(screen_coeffs_kind(:,:)%x(2))
    ALLOCATE(set_list_ij((max_set*load_balance_parameter%block_size)**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(set_list_kl((max_set*load_balance_parameter%block_size)**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

!$OMP BARRIER

    IF(actual_x_data%memory_parameter%recalc_forces) THEN
     actual_x_data%distribution_forces%ram_counter = HUGE(distribution_forces%ram_counter)
     memory_parameter%ram_counter_forces = HUGE(memory_parameter%ram_counter_forces)
    END IF

    IF(treat_forces_in_core) THEN
      buffer_overflow = .FALSE.
    ELSE
      buffer_overflow = .TRUE.
    END IF

    do_dynamic_load_balancing = .TRUE.
    IF( n_threads == 1 ) do_dynamic_load_balancing = .FALSE.

    IF( do_dynamic_load_balancing ) THEN
      my_bin_size = SIZE(actual_x_data%distribution_forces)
    ELSE
      my_bin_size = 1
    END IF
    !! We do not need the containers if MAX_MEM = 0
    IF( treat_forces_in_core ) THEN
      !! IF new md step -> reinitialize containers
      IF(actual_x_data%memory_parameter%recalc_forces) THEN
        CALL dealloc_containers(actual_x_data, hfx_do_eval_forces, error)
        CALL alloc_containers(actual_x_data, my_bin_size, hfx_do_eval_forces, error)

        DO bin = 1,my_bin_size
          maxval_container => actual_x_data%maxval_container_forces(bin)
          integral_containers => actual_x_data%integral_containers_forces(:,bin)
          CALL hfx_init_container(maxval_container, memory_parameter%actual_memory_usage, .FALSE., error)
          DO i=1,64
            CALL hfx_init_container(integral_containers(i),memory_parameter%actual_memory_usage, .FALSE., error)
          END DO
        END DO
      END IF

      !! Decompress the first cache for maxvals and integrals
      IF( .NOT. actual_x_data%memory_parameter%recalc_forces) THEN
        DO bin = 1,my_bin_size
          maxval_cache => actual_x_data%maxval_cache_forces(bin)
          maxval_container => actual_x_data%maxval_container_forces(bin)
          integral_caches => actual_x_data%integral_caches_forces(:,bin)
          integral_containers => actual_x_data%integral_containers_forces(:,bin)
          CALL hfx_decompress_first_cache(bits_max_val, maxval_cache, maxval_container, &
                                          memory_parameter%actual_memory_usage, .FALSE.)
          DO i=1,64
            CALL hfx_decompress_first_cache(i, integral_caches(i), integral_containers(i), &
                                            memory_parameter%actual_memory_usage,.FALSE.)
          END DO
        END DO
      END IF
    END IF

!$OMP BARRIER
!$OMP MASTER
    IF( do_dynamic_load_balancing ) THEN
      ! ** Lets contstruct the task list
      shm_total_bins = 0
      DO i = 1,n_threads
        shm_total_bins = shm_total_bins + SIZE(qs_env%x_data(irep, i)%distribution_forces)
      END DO
      ALLOCATE(tmp_task_list(shm_total_bins), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      shm_task_counter = 0
      DO i = 1,n_threads
        DO bin = 1,SIZE(qs_env%x_data(irep,i)%distribution_forces)
          shm_task_counter = shm_task_counter + 1
          tmp_task_list(shm_task_counter)%thread_id = i
          tmp_task_list(shm_task_counter)%bin_id = bin
          tmp_task_list(shm_task_counter)%cost = qs_env%x_data(irep,i)%distribution_forces(bin)%cost
        END DO
      END DO

      ! ** Now sort the task list
      ALLOCATE(tmp_task_list_cost(shm_total_bins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(tmp_index(shm_total_bins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO i = 1,shm_total_bins
        tmp_task_list_cost(i) = tmp_task_list(i)%cost
      END DO

      CALL sort(tmp_task_list_cost,shm_total_bins,tmp_index)

      ALLOCATE(actual_x_data%task_list(shm_total_bins), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO i = 1,shm_total_bins
        actual_x_data%task_list(i) = tmp_task_list(tmp_index(shm_total_bins - i +1))
      END DO

      shm_task_list => actual_x_data%task_list
      shm_task_counter = 0

      DEALLOCATE(tmp_task_list_cost, tmp_index, tmp_task_list, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    !! precalculate maximum density matrix elements in blocks
    shm_pmax_block => actual_x_data%pmax_block
    actual_x_data%pmax_block = 0.0_dp
    IF ( screen_pmat_forces ) THEN
      DO iatom_block = 1,SIZE(actual_x_data%blocks)
        iatom_start = actual_x_data%blocks(iatom_block)%istart
        iatom_end = actual_x_data%blocks(iatom_block)%iend
        DO jatom_block = 1,SIZE(actual_x_data%blocks)
          jatom_start = actual_x_data%blocks(jatom_block)%istart
          jatom_end = actual_x_data%blocks(jatom_block)%iend
          shm_pmax_block(iatom_block,jatom_block) = MAXVAL(shm_pmax_atom(iatom_start:iatom_end,jatom_start:jatom_end))
        END DO
      END DO
    END IF
    shm_atomic_pair_list => actual_x_data%atomic_pair_list_forces
    IF( actual_x_data%memory_parameter%recalc_forces ) THEN
      CALL build_atomic_pair_list(natom, shm_atomic_pair_list, kind_of,basis_parameter,particle_set, &
                                  do_periodic,screen_coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell,&
                                  actual_x_data%blocks)
    END IF

    my_bin_size = SIZE(actual_x_data%distribution_forces)
    DO bin=1,my_bin_size
       actual_x_data%distribution_forces(bin)%time_forces=0.0_dp
    ENDDO
!$OMP END MASTER
!$OMP BARRIER

    IF(treat_forces_in_core) THEN
      IF(my_bin_size > 0 ) THEN
        maxval_container       => actual_x_data%maxval_container_forces(1)
        maxval_cache           => actual_x_data%maxval_cache_forces(1)

        integral_containers    => actual_x_data%integral_containers_forces(:,1)
        integral_caches        => actual_x_data%integral_caches_forces(:,1)
      END IF
    END IF

    nblocks = load_balance_parameter%nblocks

    bins_left = .TRUE.
    do_it = .TRUE.
    bin = 0
    DO WHILE (bins_left)
      IF( .NOT. do_dynamic_load_balancing) THEN
        bin = bin + 1
        IF( bin > my_bin_size ) THEN
          do_it = .FALSE.
          bins_left = .FALSE.
        ELSE
          do_it = .TRUE.
          bins_left = .TRUE.
          distribution_forces => actual_x_data%distribution_forces(bin)
         END IF
      ELSE
!$OMP CRITICAL
        shm_task_counter = shm_task_counter + 1
        IF( shm_task_counter <= shm_total_bins) THEN
          my_thread_id = shm_task_list(shm_task_counter)%thread_id
          my_bin_id = shm_task_list(shm_task_counter)%bin_id
          IF( treat_forces_in_core) THEN
            maxval_cache => qs_env%x_data(irep,my_thread_id)%maxval_cache_forces(my_bin_id)
            maxval_container => qs_env%x_data(irep,my_thread_id)%maxval_container_forces(my_bin_id)
            integral_caches => qs_env%x_data(irep,my_thread_id)%integral_caches_forces(:,my_bin_id)
            integral_containers => qs_env%x_data(irep,my_thread_id)%integral_containers_forces(:,my_bin_id)
          END IF
          distribution_forces => qs_env%x_data(irep,my_thread_id)%distribution_forces(my_bin_id)
          do_it = .TRUE.
          bins_left = .TRUE.
          IF( actual_x_data%memory_parameter%recalc_forces) THEN
            distribution_forces%ram_counter = HUGE(distribution_forces%ram_counter)
          END IF
          counter = 0_Int_8
        ELSE
          do_it = .FALSE.
          bins_left = .FALSE.
        END IF
!$OMP END CRITICAL
      END IF

      IF( .NOT. do_it ) CYCLE
!$OMP MASTER
      CALL timeset(routineN//"_bin",handle_bin)
!$OMP END MASTER
      bintime_start=m_walltime()
      my_istart = distribution_forces%istart
      my_current_counter = 0
      IF(distribution_forces%number_of_atom_quartets == 0 .OR. &
         my_istart == -1_int_8) my_istart = nblocks **4
      atomic_blocks: DO atom_block = my_istart,nblocks**4-1,n_processes
        latom_block = MODULO(atom_block,nblocks)+1
        tmp_block = atom_block/nblocks
        katom_block = MODULO(tmp_block,nblocks)+1
        IF(latom_block<katom_block) CYCLE
        tmp_block = tmp_block/nblocks
        jatom_block = MODULO(tmp_block,nblocks)+1
        tmp_block = tmp_block/nblocks
        iatom_block = MODULO(tmp_block,nblocks)+1
        IF(jatom_block<iatom_block) CYCLE
        my_current_counter = my_current_counter + 1
        IF(my_current_counter > distribution_forces%number_of_atom_quartets) EXIT atomic_blocks

        iatom_start = actual_x_data%blocks(iatom_block)%istart
        iatom_end =   actual_x_data%blocks(iatom_block)%iend
        jatom_start = actual_x_data%blocks(jatom_block)%istart
        jatom_end =   actual_x_data%blocks(jatom_block)%iend
        katom_start = actual_x_data%blocks(katom_block)%istart
        katom_end =   actual_x_data%blocks(katom_block)%iend
        latom_start = actual_x_data%blocks(latom_block)%istart
        latom_end =   actual_x_data%blocks(latom_block)%iend

        pmax_blocks = MAX( shm_pmax_block(katom_block, iatom_block) + &
                           shm_pmax_block(latom_block, jatom_block), &
                           shm_pmax_block(latom_block, iatom_block) + &
                           shm_pmax_block(katom_block, jatom_block) )

        IF( 2.0_dp*coeffs_kind_max0 + pmax_blocks  < log10_eps_schwarz ) CYCLE

        CALL build_pair_list(natom,list_ij,set_list_ij,iatom_start, iatom_end, &
                             jatom_start, jatom_end, &
                             kind_of,basis_parameter,particle_set, &
                             do_periodic,screen_coeffs_set,screen_coeffs_kind,coeffs_kind_max0, &
                             log10_eps_schwarz,cell,pmax_blocks, shm_atomic_pair_list)

        CALL build_pair_list(natom,list_kl,set_list_kl, katom_start, katom_end, &
                             latom_start, latom_end, &
                             kind_of,basis_parameter,particle_set, &
                             do_periodic,screen_coeffs_set,screen_coeffs_kind,coeffs_kind_max0, &
                             log10_eps_schwarz,cell,pmax_blocks,shm_atomic_pair_list)


        DO i_list_ij=1,list_ij%n_element
          iatom=list_ij%elements(i_list_ij)%pair(1)
          jatom=list_ij%elements(i_list_ij)%pair(2)
          i_set_list_ij_start=list_ij%elements(i_list_ij)%set_bounds(1)
          i_set_list_ij_stop=list_ij%elements(i_list_ij)%set_bounds(2)
          ikind=list_ij%elements(i_list_ij)%kind_pair(1)
          jkind=list_ij%elements(i_list_ij)%kind_pair(2)
          ra=list_ij%elements(i_list_ij)%r1
          rb=list_ij%elements(i_list_ij)%r2
          rab2=list_ij%elements(i_list_ij)%dist2

          la_max => basis_parameter(ikind)%lmax
          la_min => basis_parameter(ikind)%lmin
          npgfa => basis_parameter(ikind)%npgf
          nseta = basis_parameter(ikind)%nset
          zeta => basis_parameter(ikind)%zet
          nsgfa => basis_parameter(ikind)%nsgf
          sphi_a_ext => basis_parameter(ikind)%sphi_ext(:,:,:,:)
          nsgfl_a => basis_parameter(ikind)%nsgfl
          sphi_a_u1=UBOUND(sphi_a_ext,1)
          sphi_a_u2=UBOUND(sphi_a_ext,2)
          sphi_a_u3=UBOUND(sphi_a_ext,3)


          lb_max => basis_parameter(jkind)%lmax
          lb_min => basis_parameter(jkind)%lmin
          npgfb => basis_parameter(jkind)%npgf
          nsetb = basis_parameter(jkind)%nset
          zetb => basis_parameter(jkind)%zet
          nsgfb => basis_parameter(jkind)%nsgf
          sphi_b_ext => basis_parameter(jkind)%sphi_ext(:,:,:,:)
          nsgfl_b => basis_parameter(jkind)%nsgfl
          sphi_b_u1=UBOUND(sphi_b_ext,1)
          sphi_b_u2=UBOUND(sphi_b_ext,2)
          sphi_b_u3=UBOUND(sphi_b_ext,3)

          i_atom = atom_of_kind(iatom)
          forces_map(1,1) = ikind
          forces_map(1,2) = i_atom
          j_atom = atom_of_kind(jatom)
          forces_map(2,1) = jkind
          forces_map(2,2) = j_atom

          DO i_list_kl=1,list_kl%n_element
            katom=list_kl%elements(i_list_kl)%pair(1)
            latom=list_kl%elements(i_list_kl)%pair(2)

            !All four centers equivalent => zero-contribution
!VIIRAL            IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE
            IF( .NOT. (katom+latom<=iatom+jatom) )  CYCLE
            IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE

            i_set_list_kl_start=list_kl%elements(i_list_kl)%set_bounds(1)
            i_set_list_kl_stop=list_kl%elements(i_list_kl)%set_bounds(2)
            kkind=list_kl%elements(i_list_kl)%kind_pair(1)
            lkind=list_kl%elements(i_list_kl)%kind_pair(2)
            rc=list_kl%elements(i_list_kl)%r1
            rd=list_kl%elements(i_list_kl)%r2
            rcd2=list_kl%elements(i_list_kl)%dist2

            IF( screen_pmat_forces) THEN
              pmax_atom = MAX(shm_pmax_atom(katom,iatom)+ &
                              shm_pmax_atom(latom,jatom), &
                              shm_pmax_atom(latom,iatom)+ &
                              shm_pmax_atom(katom,jatom) )
            ELSE
              pmax_atom = 0.0_dp
            END IF


            IF( (screen_coeffs_kind(jkind,ikind)%x(1)*rab2+&
                 screen_coeffs_kind(jkind,ikind)%x(2) ) + &
                 (screen_coeffs_kind(lkind,kkind)%x(1)*rcd2+&
                 screen_coeffs_kind(lkind,kkind)%x(2)) + pmax_atom < log10_eps_schwarz ) CYCLE

            IF(.NOT. (shm_is_assoc_atomic_block(latom,iatom)>=1 .AND. &
                      shm_is_assoc_atomic_block(katom,iatom)>=1 .AND. &
                      shm_is_assoc_atomic_block(katom,jatom)>=1 .AND. &
                      shm_is_assoc_atomic_block(latom,jatom)>=1 ) ) CYCLE
            k_atom = atom_of_kind(katom)
            forces_map(3,1) = kkind
            forces_map(3,2) = k_atom

            l_atom = atom_of_kind(latom)
            forces_map(4,1) = lkind
            forces_map(4,2) = l_atom

            IF(nspins == 1) THEN
              fac = 0.25_dp * hf_fraction
            ELSE
              fac = 0.5_dp * hf_fraction
            END IF
            !calculate symmetry_factor
            symm_fac = 0.25_dp
            IF(iatom==jatom) symm_fac = symm_fac*2.0_dp
            IF(katom==latom) symm_fac = symm_fac*2.0_dp
            IF(iatom==katom .AND. jatom==latom .AND. iatom/=jatom .AND. katom/=latom) symm_fac = symm_fac*2.0_dp
            IF(iatom==katom .AND. iatom==jatom .AND. katom==latom) symm_fac = symm_fac*2.0_dp

            symm_fac = 1.0_dp/symm_fac
            fac = fac * symm_fac

            lc_max => basis_parameter(kkind)%lmax
            lc_min => basis_parameter(kkind)%lmin
            npgfc => basis_parameter(kkind)%npgf
            zetc => basis_parameter(kkind)%zet
            nsgfc => basis_parameter(kkind)%nsgf
            sphi_c_ext => basis_parameter(kkind)%sphi_ext(:,:,:,:)
            nsgfl_c => basis_parameter(kkind)%nsgfl
            sphi_c_u1=UBOUND(sphi_c_ext,1)
            sphi_c_u2=UBOUND(sphi_c_ext,2)
            sphi_c_u3=UBOUND(sphi_c_ext,3)


            ld_max => basis_parameter(lkind)%lmax
            ld_min => basis_parameter(lkind)%lmin
            npgfd => basis_parameter(lkind)%npgf
            zetd => basis_parameter(lkind)%zet
            nsgfd => basis_parameter(lkind)%nsgf
            sphi_d_ext => basis_parameter(lkind)%sphi_ext(:,:,:,:)
            nsgfl_d => basis_parameter(lkind)%nsgfl
            sphi_d_u1=UBOUND(sphi_d_ext,1)
            sphi_d_u2=UBOUND(sphi_d_ext,2)
            sphi_d_u3=UBOUND(sphi_d_ext,3)

            atomic_offset_bd = shm_atomic_block_offset(jatom,latom)
            atomic_offset_bc = shm_atomic_block_offset(jatom,katom)
            atomic_offset_ad = shm_atomic_block_offset(iatom,latom)
            atomic_offset_ac = shm_atomic_block_offset(iatom,katom)

            IF( jatom < latom) THEN
              offset_bd_set => shm_set_offset(:,:,lkind,jkind)
            ELSE
              offset_bd_set => shm_set_offset(:,:,jkind,lkind)
            END IF
            IF( jatom < katom) THEN
              offset_bc_set => shm_set_offset(:,:,kkind,jkind)
            ELSE
              offset_bc_set => shm_set_offset(:,:,jkind,kkind)
            END IF
            IF( iatom < latom) THEN
              offset_ad_set => shm_set_offset(:,:,lkind,ikind)
            ELSE
              offset_ad_set => shm_set_offset(:,:,ikind,lkind)
            END IF
            IF( iatom < katom) THEN
              offset_ac_set => shm_set_offset(:,:,kkind,ikind)
            ELSE
              offset_ac_set => shm_set_offset(:,:,ikind,kkind)
            END IF

            IF( screen_pmat_forces ) THEN
              swap_id = 16
              kind_kind_idx = INT(get_1D_idx(kkind,ikind,INT(nkind,int_8)))
              IF( ikind >= kkind) THEN
                ptr_p_1 => shm_initial_p(kind_kind_idx)%p_kind(:,:, &
                                         actual_x_data%map_atom_to_kind_atom(katom), &
                                         actual_x_data%map_atom_to_kind_atom(iatom))
              ELSE
                ptr_p_1 => shm_initial_p(kind_kind_idx)%p_kind(:,:, &
                                         actual_x_data%map_atom_to_kind_atom(iatom), &
                                         actual_x_data%map_atom_to_kind_atom(katom))
                swap_id = swap_id + 1
              END IF
              kind_kind_idx = INT(get_1D_idx(lkind,jkind,INT(nkind,int_8)))
              IF( jkind>=lkind ) THEN
                ptr_p_2 => shm_initial_p(kind_kind_idx)%p_kind(:,:, &
                                         actual_x_data%map_atom_to_kind_atom(latom), &
                                         actual_x_data%map_atom_to_kind_atom(jatom))
              ELSE
                ptr_p_2 => shm_initial_p(kind_kind_idx)%p_kind(:,:, &
                                         actual_x_data%map_atom_to_kind_atom(jatom), &
                                         actual_x_data%map_atom_to_kind_atom(latom))
                swap_id = swap_id + 2
              END IF
              kind_kind_idx = INT(get_1D_idx(lkind,ikind,INT(nkind,int_8)))
              IF( ikind>=lkind ) THEN
                ptr_p_3 => shm_initial_p(kind_kind_idx)%p_kind(:,:, &
                                         actual_x_data%map_atom_to_kind_atom(latom), &
                                         actual_x_data%map_atom_to_kind_atom(iatom))
              ELSE
                ptr_p_3 => shm_initial_p(kind_kind_idx)%p_kind(:,:, &
                                         actual_x_data%map_atom_to_kind_atom(iatom), &
                                         actual_x_data%map_atom_to_kind_atom(latom))
                swap_id = swap_id + 4
              END IF
              kind_kind_idx = INT(get_1D_idx(kkind,jkind,INT(nkind,int_8)))
              IF( jkind>=kkind) THEN
                ptr_p_4 => shm_initial_p(kind_kind_idx)%p_kind(:,:, &
                                         actual_x_data%map_atom_to_kind_atom(katom), &
                                         actual_x_data%map_atom_to_kind_atom(jatom))
              ELSE
                ptr_p_4 => shm_initial_p(kind_kind_idx)%p_kind(:,:, &
                                         actual_x_data%map_atom_to_kind_atom(jatom), &
                                         actual_x_data%map_atom_to_kind_atom(katom))
                swap_id = swap_id + 8
              END IF
            END IF

            !! At this stage, check for memory used in compression

            IF(actual_x_data%memory_parameter%recalc_forces ) THEN
              IF(treat_forces_in_core) THEN
                ! ** We know the maximum amount of integrals that we can store per MPI-process
                ! ** Now we need to sum the current memory usage among all openMP threads
                ! ** We can just read what is currently stored on the corresponding x_data type
                ! ** If this is thread i and it tries to read the data from thread j, that is
                ! ** currently writing that data, we just dont care, because the possible error
                ! ** is of the order of CACHE_SIZE
                mem_compression_counter = 0
                DO i=1,n_threads
                  mem_compression_counter = mem_compression_counter + &
                                            qs_env%x_data(irep,i)%memory_parameter%actual_memory_usage *&
                                            memory_parameter%cache_size
                END DO
                IF(mem_compression_counter + memory_parameter%final_comp_counter_energy &
                                      > memory_parameter%max_compression_counter) THEN
                  buffer_overflow = .TRUE.
                  IF( do_dynamic_load_balancing ) THEN
                    distribution_forces%ram_counter = counter
                  ELSE
                    memory_parameter%ram_counter_forces = counter
                  END IF
                ELSE
                  counter = counter + 1
                  buffer_overflow = .FALSE.
                END IF
              END IF
            ELSE
              IF(treat_forces_in_core ) THEN
                IF( do_dynamic_load_balancing) THEN
                  IF(distribution_forces%ram_counter == counter) THEN
                    buffer_overflow = .TRUE.
                  ELSE
                    counter = counter + 1
                    buffer_overflow = .FALSE.
                  END IF
                ELSE
                  IF(memory_parameter%ram_counter_forces == counter) THEN
                    buffer_overflow = .TRUE.
                  ELSE
                    counter = counter + 1
                    buffer_overflow = .FALSE.
                  END IF
                END IF
              END IF
            END IF


            DO i_set_list_ij=i_set_list_ij_start, i_set_list_ij_stop
              iset=set_list_ij(i_set_list_ij)%pair(1)
              jset=set_list_ij(i_set_list_ij)%pair(2)

              ncob = npgfb(jset)*ncoset(lb_max(jset))
              max_val1 = screen_coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + &
                         screen_coeffs_set(jset,iset,jkind,ikind)%x(2)
              !! Near field screening
              IF (max_val1 + (screen_coeffs_kind(lkind,kkind)%x(1)*rcd2 + &
                              screen_coeffs_kind(lkind,kkind)%x(2) ) + pmax_atom < log10_eps_schwarz) CYCLE
              sphi_a_ext_set => sphi_a_ext(:,:,:,iset)
              sphi_b_ext_set => sphi_b_ext(:,:,:,jset)


              DO i_set_list_kl=i_set_list_kl_start, i_set_list_kl_stop
                kset=set_list_kl(i_set_list_kl)%pair(1)
                lset=set_list_kl(i_set_list_kl)%pair(2)

                max_val2_set = (screen_coeffs_set(lset,kset,lkind,kkind)%x(1)*rcd2 + &
                                       screen_coeffs_set(lset,kset,lkind,kkind)%x(2) )
                max_val2 = max_val1 + max_val2_set

                !! Near field screening
                IF(max_val2+ pmax_atom <log10_eps_schwarz) CYCLE
                sphi_c_ext_set => sphi_c_ext(:,:,:,kset)
                sphi_d_ext_set => sphi_d_ext(:,:,:,lset)

                IF( screen_pmat_forces ) THEN
                  CALL get_pmax_val(ptr_p_1, ptr_p_2, ptr_p_3, ptr_p_4, &
                                    iset, jset, kset, lset, &
                                    pmax_tmp, swap_id)

                  log10_pmax = log_2 + pmax_tmp
                ELSE
                  log10_pmax = 0.0_dp
                END IF

                max_val2 = max_val2 + log10_pmax
                IF(max_val2<log10_eps_schwarz)  CYCLE

                pmax_entry = EXP(log10_pmax*ln_10)
                !! store current number of integrals, update total number and number of integrals in buffer
                current_counter = nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)*12
                IF(buffer_overflow) THEN
                  neris_onthefly = neris_onthefly + current_counter
                END IF

                !! Get integrals from buffer and update Kohn-Sham matrix
                IF(.NOT.buffer_overflow .AND. .NOT.actual_x_data%memory_parameter%recalc_forces) THEN
                  nints = current_counter
                  CALL hfx_get_single_cache_element(estimate_to_store_int, 6,&
                                             maxval_cache, maxval_container, memory_parameter%actual_memory_usage, &
                                             use_disk_storage)
                  spherical_estimate = SET_EXPONENT(1.0_dp,estimate_to_store_int+1)
                  IF(spherical_estimate*pmax_entry<eps_schwarz) CYCLE
                  nbits = EXPONENT(ANINT(spherical_estimate*pmax_entry/eps_storage)) + 1
                  buffer_left = nints
                  buffer_start = 1
                  DO WHILE (buffer_left > 0 )
                    buffer_size = MIN(buffer_left, cache_size)
                    CALL hfx_get_mult_cache_elements(primitive_forces(buffer_start), &
                                                     buffer_size, nbits, &
                                                     integral_caches(nbits), &
                                                     integral_containers(nbits), &
                                                     eps_storage, pmax_entry, &
                                                     memory_parameter%actual_memory_usage, &
                                                     use_disk_storage)
                    buffer_left = buffer_left - buffer_size
                    buffer_start = buffer_start + buffer_size
                  END DO
                END IF
                !! Calculate integrals if we run out of buffer or the geometry did change
                IF(actual_x_data%memory_parameter%recalc_forces .OR. buffer_overflow) THEN
                  max_contraction_val =  max_contraction(iset,iatom) * &
                                         max_contraction(jset,jatom) * &
                                         max_contraction(kset,katom) * &
                                         max_contraction(lset,latom) * &
                                         pmax_entry
                  tmp_R_1 => radii_pgf(:,:,jset,iset,jkind,ikind)
                  tmp_R_2 => radii_pgf(:,:,lset,kset,lkind,kkind)
                  tmp_screen_pgf1 => screen_coeffs_pgf(:,:,jset,iset,jkind,ikind)
                  tmp_screen_pgf2 => screen_coeffs_pgf(:,:,lset,kset,lkind,kkind)

                  CALL forces4(private_deriv, ra, rb, rc, rd, npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset), &
                               la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                               lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                               nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset),&
                               sphi_a_u1,sphi_a_u2,sphi_a_u3,&
                               sphi_b_u1,sphi_b_u2,sphi_b_u3,&
                               sphi_c_u1,sphi_c_u2,sphi_c_u3,&
                               sphi_d_u1,sphi_d_u2,sphi_d_u3,&
                               zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset),&
                               zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                               primitive_forces,&
                               potential_parameter, &
                               actual_x_data%neighbor_cells, screen_coeffs_set(jset,iset,jkind,ikind)%x,&
                               screen_coeffs_set(lset,kset,lkind,kkind)%x, eps_schwarz, &
                               max_contraction_val, cartesian_estimate,  cell, neris_tmp, &
                               log10_pmax, log10_eps_schwarz, &
                               tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2,&
                               pgf_list_ij, pgf_list_kl, pgf_product_list, &
                               nsgfl_a(:,iset), nsgfl_b(:,jset),&
                               nsgfl_c(:,kset), nsgfl_d(:,lset),&
                               sphi_a_ext_set, sphi_b_ext_set, sphi_c_ext_set, sphi_d_ext_set,&
                               ede_work,ede_work2,ede_work_forces, &
                               ede_buffer1,ede_buffer2,ede_primitives_tmp, &
                               nimages,do_periodic,use_virial, ede_work_virial, ede_work2_virial, &
                               ede_primitives_tmp_virial, primitive_forces_virial, cartesian_estimate_virial)

                  nints = nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)*12
                  neris_total = neris_total + nints
                  nprim_ints = nprim_ints + neris_tmp
                  IF(cartesian_estimate == 0.0_dp) cartesian_estimate = TINY(cartesian_estimate)
                  estimate_to_store_int = EXPONENT(cartesian_estimate)
                  estimate_to_store_int = MAX(estimate_to_store_int,-15_int_8)
                  cartesian_estimate = SET_EXPONENT(1.0_dp,estimate_to_store_int+1)
                  IF(.NOT.buffer_overflow .AND. actual_x_data%memory_parameter%recalc_forces) THEN
                    IF(cartesian_estimate < eps_schwarz ) THEN
                      CALL hfx_add_single_cache_element(estimate_to_store_int, 6,&
                                        maxval_cache, maxval_container, memory_parameter%actual_memory_usage, &
                                        use_disk_storage, max_val_memory, error)
                    END IF
                  END IF

                  IF( .NOT. use_virial ) THEN
                    IF(cartesian_estimate<eps_schwarz) CYCLE
                  END IF

                  !! Compress the array for storage
                  spherical_estimate = 0.0_dp
                  DO i=1,nints
                    spherical_estimate = MAX(spherical_estimate,ABS(primitive_forces(i)))
                  END DO

                  IF( use_virial ) THEN
                    spherical_estimate_virial = 0.0_dp
                    DO i=1,3*nints
                      spherical_estimate_virial = MAX(spherical_estimate_virial, ABS(primitive_forces_virial(i)))
                    END DO
                  END IF

                  IF(spherical_estimate == 0.0_dp) spherical_estimate = TINY(spherical_estimate)
                  estimate_to_store_int = EXPONENT(spherical_estimate)
                  estimate_to_store_int = MAX(estimate_to_store_int, -15_int_8)
                  IF(.NOT. buffer_overflow .AND. actual_x_data%memory_parameter%recalc_forces) THEN
                    CALL hfx_add_single_cache_element(estimate_to_store_int, 6, &
                                                      maxval_cache, maxval_container, &
                                                      memory_parameter%actual_memory_usage, &
                                                      use_disk_storage, max_val_memory, error)
                  END IF
                  spherical_estimate = SET_EXPONENT(1.0_dp,estimate_to_store_int+1)
                  IF( .NOT. use_virial ) THEN
                    IF(spherical_estimate*pmax_entry<eps_schwarz) CYCLE
                  END IF
                  IF(.NOT.buffer_overflow) THEN
                    nbits = EXPONENT(ANINT(spherical_estimate*pmax_entry/eps_storage)) + 1
                    buffer_left = nints
                    buffer_start = 1
                    neris_incore = neris_incore + nints

                    DO WHILE (buffer_left > 0 )
                      buffer_size = MIN(buffer_left, CACHE_SIZE)
                      CALL hfx_add_mult_cache_elements(primitive_forces(buffer_start),&
                                                       buffer_size, nbits, &
                                                       integral_caches(nbits), &
                                                       integral_containers(nbits), &
                                                       eps_storage, pmax_entry, &
                                                       memory_parameter%actual_memory_usage, &
                                                       use_disk_storage, error)
                      buffer_left = buffer_left - buffer_size
                      buffer_start = buffer_start + buffer_size
                    END DO
                  ELSE
                    !! In order to be consistent with in-core part, round all the eris wrt. eps_schwarz
                    !! but only if we treat forces in-core
                    IF( treat_forces_in_core ) THEN
                      DO i=1,nints
                        primitive_forces(i) = primitive_forces(i)*pmax_entry
                        IF(ABS(primitive_forces(i)) > eps_storage) THEN
                          primitive_forces(i) = ANINT(primitive_forces(i)/eps_storage, dp)*eps_storage/pmax_entry
                        ELSE
                          primitive_forces(i) = 0.0_dp
                        END IF
                      END DO
                    END IF
                  END IF
                END IF
                IF(with_mp2_density) THEN
                  CALL prefetch_density_matrix(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                               full_density, pbd_buf, pbc_buf, pad_buf, pac_buf, &
                                               iatom, jatom, katom, latom, &
                                               iset, jset, kset, lset, offset_bd_set, offset_bc_set,&
                                               offset_ad_set, offset_ac_set, atomic_offset_bd, atomic_offset_bc,&
                                               atomic_offset_ad, atomic_offset_ac,is_anti_symmetric)
                  CALL prefetch_density_matrix(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                               full_density_mp2, pbd_buf_mp2, pbc_buf_mp2, pad_buf_mp2, pac_buf_mp2, &
                                               iatom, jatom, katom, latom, &
                                               iset, jset, kset, lset, offset_bd_set, offset_bc_set,&
                                               offset_ad_set, offset_ac_set, atomic_offset_bd, atomic_offset_bc,&
                                               atomic_offset_ad, atomic_offset_ac,is_anti_symmetric)
                ELSE
                  CALL prefetch_density_matrix(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                               full_density, pbd_buf, pbc_buf, pad_buf, pac_buf, &
                                               iatom, jatom, katom, latom, &
                                               iset, jset, kset, lset, offset_bd_set, offset_bc_set,&
                                               offset_ad_set, offset_ac_set, atomic_offset_bd, atomic_offset_bc,&
                                               atomic_offset_ad, atomic_offset_ac,is_anti_symmetric)
                END IF
                IF( nspins == 2 ) THEN
                  CALL prefetch_density_matrix(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                               full_density_beta, pbd_buf_beta, pbc_buf_beta, pad_buf_beta, pac_buf_beta, &
                                               iatom, jatom, katom, latom, &
                                               iset, jset, kset, lset, offset_bd_set, offset_bc_set,&
                                               offset_ad_set, offset_ac_set, atomic_offset_bd, atomic_offset_bc,&
                                               atomic_offset_ad, atomic_offset_ac,is_anti_symmetric)
                END IF
                IF( spherical_estimate*pmax_entry >= eps_schwarz ) THEN
                  DO coord = 1,12
                    T2 => primitive_forces((coord-1)*nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)+1:&
                                            coord*nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset))

                    IF(with_mp2_density) THEN
                      CALL update_forces(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                         pbd_buf, pbc_buf, pad_buf, pac_buf, fac, &
                                         T2, force, forces_map, coord,&
                                         pbd_buf_mp2, pbc_buf_mp2, pad_buf_mp2, pac_buf_mp2)
                    ELSE
                      CALL update_forces(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                         pbd_buf, pbc_buf, pad_buf, pac_buf, fac, &
                                         T2, force, forces_map, coord)
                    END IF
                    IF( nspins == 2 ) THEN
                      CALL update_forces(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                         pbd_buf_beta, pbc_buf_beta, pad_buf_beta, pac_buf_beta, fac, &
                                         T2, force, forces_map, coord)
                    END IF
                  END DO !coord
                END IF
                IF (use_virial .AND. spherical_estimate_virial*pmax_entry >=eps_schwarz ) THEN
                  DO coord = 1,12
                    DO i=1,3
                      T2 => primitive_forces_virial(((i-1)*12+coord-1)*nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)+1:&
                                              ((i-1)*12+coord)*nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset))
                      IF(with_mp2_density) THEN
                        CALL update_virial(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                           pbd_buf, pbc_buf, pad_buf, pac_buf, fac, &
                                           T2, tmp_virial, coord, i, &
                                           pbd_buf_mp2, pbc_buf_mp2, pad_buf_mp2, pac_buf_mp2)
                      ELSE
                        CALL update_virial(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                           pbd_buf, pbc_buf, pad_buf, pac_buf, fac, &
                                           T2, tmp_virial, coord, i)
                      END IF
                      IF( nspins == 2 ) THEN
                        CALL update_virial(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                           pbd_buf_beta, pbc_buf_beta, pad_buf_beta, pac_buf_beta, fac, &
                                           T2, tmp_virial, coord, i)
                      END IF
                    END DO
                  END DO !coord
                END IF

              END DO !i_set_list_kl
            END DO !i_set_list_ij
          END DO !i_list_kl
        END DO !i_list_ij
      END DO atomic_blocks
      bintime_stop=m_walltime()
      distribution_forces%time_forces=bintime_stop-bintime_start
!$OMP MASTER
      CALL timestop(handle_bin)
!$OMP END MASTER
    END DO  !bin

!$OMP MASTER
    logger => cp_error_get_logger(error)
    do_print_load_balance_info = .FALSE.
    do_print_load_balance_info = BTEST(cp_print_key_should_output(logger%iter_info,hfx_section,&
               "LOAD_BALANCE%PRINT/LOAD_BALANCE_INFO",error=error),cp_p_file)
!$OMP END MASTER
!$OMP BARRIER
    IF(do_print_load_balance_info) THEN
      iw = -1
!$OMP MASTER
      iw = cp_print_key_unit_nr(logger,hfx_section,"LOAD_BALANCE%PRINT/LOAD_BALANCE_INFO",&
                              extension=".scfLog",error=error)
!$OMP END MASTER

      CALL collect_load_balance_info(para_env, actual_x_data, iw, n_threads, i_thread, &
                                     hfx_do_eval_forces, error)

!$OMP MASTER
      CALL cp_print_key_finished_output(iw,logger,hfx_section,&
                                       "LOAD_BALANCE%PRINT/LOAD_BALANCE_INFO", error=error)
!$OMP END MASTER
    END IF


    IF ( use_virial ) THEN
       DO i = 1, 3
          DO j = 1, 3
             DO k = 1, 3
                !$OMP ATOMIC
                virial%pv_fock_4c(i,j) = virial%pv_fock_4c(i,j) + tmp_virial(i,k)*cell%hmat(j,k)
             END DO
          END DO
       END DO
    END IF

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
    stor_count_max_val = max_val_memory * memory_parameter%cache_size
!$OMP ATOMIC
    shm_storage_counter_integrals = shm_storage_counter_integrals + storage_counter_integrals
!$OMP ATOMIC
    shm_neris_incore = shm_neris_incore + neris_incore
!$OMP ATOMIC
    shm_stor_count_max_val = shm_stor_count_max_val + stor_count_max_val
!$OMP BARRIER

    ! IF(with_mp2_density) THEN
    !   ! restore original load_balance_parameter
    !   actual_x_data%load_balance_parameter = load_balance_parameter_energy
    !   DEALLOCATE(load_balance_parameter_energy)
    ! END IF

!$OMP MASTER
    !! Print some memeory information if this is the first step
    IF( actual_x_data%memory_parameter%recalc_forces ) THEN
      tmp_i8(1:6)=(/shm_storage_counter_integrals,shm_neris_onthefly,shm_neris_incore,&
                    shm_neris_total,shm_nprim_ints,shm_stor_count_max_val/)
      CALL mp_sum(tmp_i8,para_env%group)
      shm_storage_counter_integrals=tmp_i8(1)
      shm_neris_onthefly=tmp_i8(2)
      shm_neris_incore=tmp_i8(3)
      shm_neris_total=tmp_i8(4)
      shm_nprim_ints=tmp_i8(5)
      shm_stor_count_max_val=tmp_i8(6)
      mem_eris = (shm_storage_counter_integrals+128*1024-1)/1024/128
      compression_factor = REAL(shm_neris_incore,dp)/REAL(shm_storage_counter_integrals,dp)
      mem_max_val = (shm_stor_count_max_val+128*1024-1)/1024/128

      IF( shm_neris_incore == 0 ) THEN
        mem_eris = 0
        compression_factor = 0.0_dp
      END IF

      logger => cp_error_get_logger(error)

      iw = cp_print_key_unit_nr(logger,hfx_section,"HF_INFO",&
           extension=".scfLog",error=error)
      IF (iw>0) THEN
        WRITE (UNIT=iw,FMT="(/,(T3,A,T65,I16))")&
                 "HFX_MEM_INFO| Number of cart. primitive DERIV's calculated:      ",shm_nprim_ints

        WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Number of sph. DERIV's calculated:           ",shm_neris_total

        WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Number of sph. DERIV's stored in-core:        ",shm_neris_incore

        WRITE (UNIT=iw,FMT="((T3,A,T62,I19))")&
                 "HFX_MEM_INFO| Number of sph. DERIV's calculated on the fly: ",shm_neris_onthefly

        WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Total memory consumption DERIV's RAM [MB's]:            ",mem_eris

        WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Whereof max-vals [MB's]:            ",mem_max_val

        WRITE (UNIT=iw,FMT="((T3,A,T60,F21.2),/)")&
                 "HFX_MEM_INFO| Total compression factor DERIV's RAM:                  ",compression_factor
      END IF

      CALL cp_print_key_finished_output(iw,logger,hfx_section,&
                                      "HF_INFO", error=error)
    END IF
!$OMP END MASTER

    !! flush caches if the geometry changed
    IF( do_dynamic_load_balancing ) THEN
      my_bin_size = SIZE(actual_x_data%distribution_forces)
    ELSE
      my_bin_size = 1
    END IF

    IF(actual_x_data%memory_parameter%recalc_forces ) THEN
      IF(treat_forces_in_core) THEN
        DO bin = 1,my_bin_size
          maxval_cache => actual_x_data%maxval_cache_forces(bin)
          maxval_container => actual_x_data%maxval_container_forces(bin)
          integral_caches => actual_x_data%integral_caches_forces(:,bin)
          integral_containers => actual_x_data%integral_containers_forces(:,bin)
          CALL hfx_flush_last_cache(bits_max_val, maxval_cache, maxval_container, memory_parameter%actual_memory_usage, &
                                    .FALSE., error)
          DO i=1,64
            CALL hfx_flush_last_cache(i, integral_caches(i), integral_containers(i), &
                                      memory_parameter%actual_memory_usage, .FALSE., error)
          END DO
        END DO
      END IF
    END IF
    !! reset all caches except we calculate all on the fly
    IF(treat_forces_in_core) THEN
      DO bin = 1,my_bin_size
        maxval_cache => actual_x_data%maxval_cache_forces(bin)
        maxval_container => actual_x_data%maxval_container_forces(bin)
        integral_caches => actual_x_data%integral_caches_forces(:,bin)
        integral_containers => actual_x_data%integral_containers_forces(:,bin)

        CALL hfx_reset_cache_and_container(maxval_cache, maxval_container, memory_parameter%actual_memory_usage, .FALSE.)
        DO i=1,64
          CALL hfx_reset_cache_and_container(integral_caches(i), integral_containers(i), &
                                             memory_parameter%actual_memory_usage, &
                                             .FALSE.)
        END DO
      END DO
    END IF

    IF( actual_x_data%memory_parameter%recalc_forces ) THEN
      actual_x_data%memory_parameter%recalc_forces = .FALSE.
    END IF

!$OMP BARRIER
!$OMP MASTER
    CALL timestop(handle_main)
!$OMP END MASTER
!$OMP BARRIER
    DEALLOCATE(last_sgf_global, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP MASTER
    DEALLOCATE(full_density,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(with_mp2_density) THEN
      DEALLOCATE(full_density_mp2,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(nspins==2) THEN
      DEALLOCATE(full_density_beta,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

!$OMP END MASTER
    DEALLOCATE(primitive_forces, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(work_forces, atom_of_kind,kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(max_contraction,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(pbd_buf, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(pbc_buf, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(pad_buf, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(pac_buf, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF(with_mp2_density) THEN
      DEALLOCATE(pbd_buf_mp2, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(pbc_buf_mp2, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(pad_buf_mp2, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(pac_buf_mp2, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    DO i=1,max_pgf**2
      DEALLOCATE(pgf_list_ij(i)%image_list, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(pgf_list_kl(i)%image_list, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO

    DEALLOCATE(pgf_list_ij, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(pgf_list_kl, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(pgf_product_list, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(set_list_ij, set_list_kl, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(primitive_forces_virial, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(ede_work, ede_work2, ede_work_forces, ede_buffer1, ede_buffer2, ede_primitives_tmp, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(ede_work_virial, ede_work2_virial, ede_primitives_tmp_virial, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(nimages, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF( nspins == 2 ) THEN
      DEALLOCATE(pbd_buf_beta,pbc_buf_beta,pad_buf_beta,pac_buf_beta, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    !
    ! this wraps to free_libint, but currently causes a segfault
    ! as a result, we don't call it, but some memory remains allocated
    !
    ! CALL terminate_libderiv(private_deriv)

!$OMP END PARALLEL

    CALL timestop(handle)

  END SUBROUTINE derivatives_four_center

  SUBROUTINE forces4(deriv, ra, rb, rc, rd, npgfa, npgfb, npgfc, npgfd, &
                     la_min, la_max, lb_min, lb_max,&
                     lc_min, lc_max, ld_min, ld_max, nsgfa, nsgfb, nsgfc, nsgfd,&
                     sphi_a_u1,sphi_a_u2,sphi_a_u3,&
                     sphi_b_u1,sphi_b_u2,sphi_b_u3,&
                     sphi_c_u1,sphi_c_u2,sphi_c_u3,&
                     sphi_d_u1,sphi_d_u2,sphi_d_u3,&
                     zeta, zetb, zetc, zetd,&
                     primitive_integrals,&
                     potential_parameter, neighbor_cells,&
                     screen1, screen2, eps_schwarz, max_contraction_val,&
                     cart_estimate, cell, neris_tmp, log10_pmax,&
                     log10_eps_schwarz,R1_pgf, R2_pgf, pgf1, pgf2,&
                     pgf_list_ij, pgf_list_kl,&
                     pgf_product_list,&
                     nsgfl_a, nsgfl_b, nsgfl_c,&
                     nsgfl_d,&
                     sphi_a_ext, sphi_b_ext, sphi_c_ext, sphi_d_ext,&
                     ede_work,ede_work2,ede_work_forces, &
                     ede_buffer1,ede_buffer2,ede_primitives_tmp, &
                     nimages, do_periodic, use_virial, ede_work_virial, &
                     ede_work2_virial, ede_primitives_tmp_virial, primitive_integrals_virial,&
                     cart_estimate_virial)

    TYPE(lib_deriv)                          :: deriv
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    INTEGER, INTENT(IN) :: npgfa, npgfb, npgfc, npgfd, la_min, la_max, &
      lb_min, lb_max, lc_min, lc_max, ld_min, ld_max, nsgfa, nsgfb, nsgfc, &
      nsgfd, sphi_a_u1, sphi_a_u2, sphi_a_u3, sphi_b_u1, sphi_b_u2, &
      sphi_b_u3, sphi_c_u1, sphi_c_u2, sphi_c_u3, sphi_d_u1, sphi_d_u2, &
      sphi_d_u3
    REAL(dp), DIMENSION(1:npgfa), INTENT(IN) :: zeta
    REAL(dp), DIMENSION(1:npgfb), INTENT(IN) :: zetb
    REAL(dp), DIMENSION(1:npgfc), INTENT(IN) :: zetc
    REAL(dp), DIMENSION(1:npgfd), INTENT(IN) :: zetd
    REAL(dp), DIMENSION(nsgfa, nsgfb, nsgfc,&
      nsgfd, 12)                             :: primitive_integrals
    TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(hfx_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(dp), INTENT(IN)                     :: screen1(2), screen2(2), &
                                                eps_schwarz, &
                                                max_contraction_val
    REAL(dp)                                 :: cart_estimate
    TYPE(cell_type), POINTER                 :: cell
    INTEGER(int_8)                           :: neris_tmp
    REAL(dp), INTENT(IN)                     :: log10_pmax, log10_eps_schwarz
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: R1_pgf, R2_pgf, pgf1, pgf2
    TYPE(hfx_pgf_list), DIMENSION(*)         :: pgf_list_ij, pgf_list_kl
    TYPE(hfx_pgf_product_list), DIMENSION(*) :: pgf_product_list
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nsgfl_a, nsgfl_b, nsgfl_c, &
                                                nsgfl_d
    REAL(dp), INTENT(IN) :: sphi_a_ext(sphi_a_u1,sphi_a_u2,sphi_a_u3), &
      sphi_b_ext(sphi_b_u1,sphi_b_u2,sphi_b_u3), &
      sphi_c_ext(sphi_c_u1,sphi_c_u2,sphi_c_u3), &
      sphi_d_ext(sphi_d_u1,sphi_d_u2,sphi_d_u3)
    REAL(dp), DIMENSION(*)                   :: ede_work, ede_work2, &
                                                ede_work_forces, ede_buffer1, &
                                                ede_buffer2, &
                                                ede_primitives_tmp
    INTEGER, DIMENSION(*)                    :: nimages
    LOGICAL, INTENT(IN)                      :: do_periodic, use_virial
    REAL(dp), DIMENSION(*)                   :: ede_work_virial, &
                                                ede_work2_virial, &
                                                ede_primitives_tmp_virial
    REAL(dp), DIMENSION(nsgfa, nsgfb, nsgfc,&
      nsgfd, 12, 3)                          :: primitive_integrals_virial
    REAL(dp)                                 :: cart_estimate_virial

    INTEGER :: ipgf, jpgf, kpgf, la, lb, lc, ld, list_ij, list_kl, lpgf, &
      max_l, ncoa, ncob, ncoc, ncod, nelements_ij, nelements_kl, nneighbors3, &
      nproducts, nsgfla, nsgflb, nsgflc, nsgfld, nsoa, nsob, nsoc, nsod, &
      s_offset_a, s_offset_b, s_offset_c, s_offset_d
    REAL(dp)                                 :: EtaInv, tmp_max, &
                                                tmp_max_virial, Zeta_A, &
                                                Zeta_B, Zeta_C, Zeta_D, &
                                                ZetaInv

    CALL build_pair_list_pgf(npgfa, npgfb, pgf_list_ij, zeta, zetb, screen1, screen2,&
                             pgf1, R1_pgf, log10_pmax, log10_eps_schwarz, ra, rb,&
                             nelements_ij, &
                             neighbor_cells, nimages, do_periodic)
    CALL build_pair_list_pgf(npgfc, npgfd, pgf_list_kl, zetc, zetd, screen2, screen1, &
                             pgf2, R2_pgf, log10_pmax, log10_eps_schwarz, rc, rd,&
                             nelements_kl, &
                             neighbor_cells, nimages, do_periodic)

    cart_estimate = 0.0_dp
    primitive_integrals = 0.0_dp
    IF( use_virial ) THEN
      primitive_integrals_virial = 0.0_dp
      cart_estimate_virial = 0.0_dp
    END IF
    neris_tmp = 0_int_8
    max_l = la_max + lb_max + lc_max + ld_max + 1
    DO list_ij = 1,nelements_ij
      Zeta_A = pgf_list_ij(list_ij)%zeta
      Zeta_B = pgf_list_ij(list_ij)%zetb
      ZetaInv = pgf_list_ij(list_ij)%ZetaInv
      ipgf = pgf_list_ij(list_ij)%ipgf
      jpgf = pgf_list_ij(list_ij)%jpgf

      DO list_kl = 1,nelements_kl
        Zeta_C = pgf_list_kl(list_kl)%zeta
        Zeta_D = pgf_list_kl(list_kl)%zetb
        EtaInv = pgf_list_kl(list_kl)%ZetaInv
        kpgf = pgf_list_kl(list_kl)%ipgf
        lpgf = pgf_list_kl(list_kl)%jpgf

        nneighbors3 = SIZE(neighbor_cells)*pgf_list_ij(list_ij)%nimages*pgf_list_kl(list_kl)%nimages
        CALL build_pgf_product_list(nneighbors3, pgf_list_ij(list_ij), pgf_list_kl(list_kl), pgf_product_list, &
                                    nproducts, log10_pmax, log10_eps_schwarz, neighbor_cells, cell, &
                                    potential_parameter, max_l, do_periodic)
        s_offset_a = 0
        DO la = la_min,la_max
          s_offset_b = 0
          ncoa = nco(la)
          nsgfla = nsgfl_a(la)
          nsoa = nso(la)

          DO lb = lb_min, lb_max
            s_offset_c = 0
            ncob = nco(lb)
            nsgflb = nsgfl_b(lb)
            nsob = nso(lb)

            DO lc = lc_min, lc_max
              s_offset_d = 0
              ncoc = nco(lc)
              nsgflc = nsgfl_c(lc)
              nsoc = nso(lc)

              DO ld = ld_min, ld_max
                ncod = nco(ld)
                nsgfld = nsgfl_d(ld)
                nsod = nso(ld)
                tmp_max = 0.0_dp
                tmp_max_virial = 0.0_dp
                CALL evaluate_deriv_eri(deriv, nproducts, pgf_product_list,  &
                                        la, lb, lc ,ld,&
                                        ncoa,ncob,ncoc,ncod,&
                                        nsgfa, nsgfb, nsgfc, nsgfd,&
                                        primitive_integrals, &
                                        max_contraction_val, tmp_max, eps_schwarz, neris_tmp,&
                                        Zeta_A, Zeta_B, Zeta_C, Zeta_D, ZetaInv,EtaInv,&
                                        s_offset_a, s_offset_b, s_offset_c, s_offset_d,&
                                        nsgfla, nsgflb, nsgflc, nsgfld, nsoa, nsob, nsoc, nsod, &
                                        sphi_a_ext(1,la+1,ipgf), &
                                        sphi_b_ext(1,lb+1,jpgf), &
                                        sphi_c_ext(1,lc+1,kpgf), &
                                        sphi_d_ext(1,ld+1,lpgf), &
                                        ede_work,ede_work2,ede_work_forces, &
                                        ede_buffer1,ede_buffer2,ede_primitives_tmp, use_virial,&
                                        ede_work_virial, ede_work2_virial, ede_primitives_tmp_virial, &
                                        primitive_integrals_virial, cell, tmp_max_virial)
                cart_estimate = MAX(tmp_max,cart_estimate)
                IF( use_virial ) cart_estimate_virial = MAX(tmp_max_virial, cart_estimate_virial)
                s_offset_d = s_offset_d + nsod*nsgfld
              END DO !ld
              s_offset_c = s_offset_c + nsoc*nsgflc
            END DO !lc
            s_offset_b = s_offset_b + nsob*nsgflb
          END DO !lb
          s_offset_a = s_offset_a + nsoa*nsgfla
        END DO !la
      END DO
    END DO

  END SUBROUTINE forces4

! *****************************************************************************
!> \brief Given a 2d index pair, this function returns a 1d index pair for
!>        a symmetric upper triangle NxN matrix
!>        The compiler should inline this function, therefore it appears in
!>        several modules
!> \param i,j 2d index
!> \param N matrix size
!> \par History
!>      03.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  PURE FUNCTION get_1D_idx(i,j,N)
    INTEGER, INTENT(IN)                      :: i, j
    INTEGER(int_8), INTENT(IN)               :: N
    INTEGER(int_8)                           :: get_1D_idx

    INTEGER(int_8)                           :: min_ij

    min_ij = MIN(i,j)
    get_1D_idx = min_ij*N + MAX(i,j) - (min_ij-1)*min_ij/2 - N

  END FUNCTION get_1D_idx

! *****************************************************************************
!> \brief This routine prefetches density matrix elements, i.e. reshuffles the
!>        data in a way that they can be accessed later on in a cache friendly
!>        way
!> \param pa1,pb1,pc1,pd1 offsets
!> \param ma_max, mb_max, mc_max, md_max: Size of matrix blocks
!> \param sgfa, sgfb, sgfc, sgfd offsets
!> \param pbd, pbc, pad, pac buffers that will contain P(b,d), P(b,c) ...
!> \par History
!>      03.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE prefetch_density_matrix(ma_max, mb_max, mc_max, md_max, &
                                     density, pbd, pbc, pad, pac, &
                                     iatom, jatom, katom, latom, &
                                     iset, jset, kset, lset, offset_bd_set, offset_bc_set,&
                                     offset_ad_set, offset_ac_set, atomic_offset_bd, atomic_offset_bc,&
                                     atomic_offset_ad, atomic_offset_ac,anti_symmetric)

    INTEGER, INTENT(IN)                      :: ma_max, mb_max, mc_max, md_max
    REAL(dp), DIMENSION(:), INTENT(IN)       :: density
    REAL(dp), DIMENSION(*), INTENT(INOUT)    :: pbd, pbc, pad, pac
    INTEGER, INTENT(IN)                      :: iatom, jatom, katom, latom, &
                                                iset, jset, kset, lset
    INTEGER, DIMENSION(:, :), POINTER        :: offset_bd_set, offset_bc_set, &
                                                offset_ad_set, offset_ac_set
    LOGICAL                                  :: anti_symmetric
    INTEGER, INTENT(IN)                      :: atomic_offset_bd, &
                                                atomic_offset_bc, &
                                                atomic_offset_ad, &
                                                atomic_offset_ac

    INTEGER                                  :: i, j, ma, mb, mc, md, &
                                                offset_ac, offset_ad, &
                                                offset_bc, offset_bd
    REAL(dp)                                 :: fac

    fac=1.0_dp
    IF(anti_symmetric)fac=-1.0_dp
    IF( jatom >= latom ) THEN
      i = 1
      offset_bd = offset_bd_set(jset,lset) + atomic_offset_bd - 1
      j = offset_bd
      DO md=1,md_max
        DO mb=1,mb_max
          pbd(i) = density(j)*fac
          i = i + 1
          j = j + 1
        END DO
      END DO
    ELSE
      i = 1
      offset_bd = offset_bd_set(lset,jset) + atomic_offset_bd - 1
      DO md=1,md_max
        j = offset_bd + md - 1
        DO mb=1,mb_max
          pbd(i) = density(j)
          i = i + 1
          j = j + md_max
        END DO
      END DO
    END IF
    IF ( jatom >= katom ) THEN
      i = 1
      offset_bc = offset_bc_set(jset,kset) + atomic_offset_bc - 1
      j = offset_bc
      DO mc=1,mc_max
        DO mb=1,mb_max
          pbc(i) = density(j)*fac
          i = i + 1
          j = j + 1
        END DO
      END DO
    ELSE
      i = 1
      offset_bc = offset_bc_set(kset,jset) + atomic_offset_bc - 1
      DO mc=1,mc_max
        j = offset_bc + mc -1
        DO mb=1,mb_max
          pbc(i) = density(j)
          i = i + 1
          j = j + mc_max
        END DO
      END DO
    END IF
    IF( iatom >= latom) THEN
      i = 1
      offset_ad = offset_ad_set(iset,lset) + atomic_offset_ad - 1
      j = offset_ad
      DO md=1,md_max
        DO ma=1,ma_max
          pad(i) = density(j)*fac
          i = i + 1
          j = j + 1
        END DO
      END DO
    ELSE
      i = 1
      offset_ad = offset_ad_set(lset,iset) + atomic_offset_ad - 1
      DO md=1,md_max
        j = offset_ad + md -1
        DO ma=1,ma_max
          pad(i) = density(j)
          i = i + 1
          j = j + md_max
        END DO
      END DO
    END IF
    IF( iatom >= katom ) THEN
      i = 1
      offset_ac = offset_ac_set(iset,kset) + atomic_offset_ac - 1
      j = offset_ac
      DO mc=1,mc_max
        DO ma=1,ma_max
          pac(i) = density(j)*fac
          i = i + 1
          j = j + 1
        END DO
      END DO
    ELSE
      i = 1
      offset_ac = offset_ac_set(kset,iset) + atomic_offset_ac - 1
      DO mc=1,mc_max
        j = offset_ac + mc - 1
        DO ma=1,ma_max
          pac(i) = density(j)
          i = i + 1
          j = j + mc_max
        END DO
      END DO
    END IF
  END SUBROUTINE prefetch_density_matrix

! *****************************************************************************
!> \brief This routine updates the forces using bufferd density matrices
!> \param ma_max, mb_max, mc_max, md_max: Size of matrix blocks
!> \param pbd, pbc, pad, pac buffers that will contain P(b,d), P(b,c) ...
!> \param fac mulitplication factor (spin, symmetry)
!> \param prim primitve forces
!> \param force storage loacation for forces
!> \param forces_map index table
!> \param coord which of the 12 coords to be updated
!> \par History
!>      03.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE update_forces(ma_max, mb_max, mc_max, md_max, &
                           pbd, pbc, pad, pac, fac, &
                           prim, force, forces_map, coord, &
                           pbd_mp2, pbc_mp2, pad_mp2, pac_mp2)

    INTEGER, INTENT(IN)                      :: ma_max, mb_max, mc_max, md_max
    REAL(dp), DIMENSION(*), INTENT(IN)       :: pbd, pbc, pad, pac
    REAL(dp), INTENT(IN)                     :: fac
    REAL(dp), DIMENSION(ma_max*mb_max*mc_max&
      *md_max), INTENT(IN)                   :: prim
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    INTEGER, INTENT(IN)                      :: forces_map(4,2), coord
    REAL(dp), DIMENSION(*), INTENT(IN), &
      OPTIONAL                               :: pbd_mp2, pbc_mp2, pad_mp2, &
                                                pac_mp2

    INTEGER                                  :: ma, mb, mc, md, p_index
    LOGICAL                                  :: with_mp2_density
    REAL(dp)                                 :: temp1, temp1_mp2, temp2, &
                                                temp3, temp3_mp2, temp4

    with_mp2_density=.FALSE.
    IF(PRESENT(pbd_mp2).AND.&
       PRESENT(pbc_mp2).AND.&
       PRESENT(pad_mp2).AND.&
       PRESENT(pac_mp2)) with_mp2_density=.TRUE.

    IF(with_mp2_density) THEN
      p_index = 0
      temp4 = 0.0_dp
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,mb_max
            temp1 = pbc((mc-1)*mb_max+mb) * fac
            temp3 = pbd((md-1)*mb_max+mb) * fac
            temp1_mp2 = pbc_mp2((mc-1)*mb_max+mb) * fac
            temp3_mp2 = pbd_mp2((md-1)*mb_max+mb) * fac
            DO ma = 1,ma_max
              p_index = p_index + 1
              ! HF-SCF
              temp2 = temp1 * pad((md-1)*ma_max+ma) + &
                      temp3 * pac((mc-1)*ma_max+ma)
              ! MP2+HF
              temp2 = temp2 + &
                      pac((mc-1)*ma_max+ma) * temp3_mp2 + &
                      pac_mp2((mc-1)*ma_max+ma) * temp3 + &
                      pad((md-1)*ma_max+ma) * temp1_mp2 + &
                      pad_mp2((md-1)*ma_max+ma) * temp1
              temp4 = temp4 + temp2 * prim(p_index) 
            END DO !ma
          END DO !mb
        END DO !mc
      END DO !md
    ELSE
      p_index = 0
      temp4 = 0.0_dp
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,mb_max
            temp1 = pbc((mc-1)*mb_max+mb) * fac
            temp3 = pbd((md-1)*mb_max+mb) * fac
            DO ma = 1,ma_max
              p_index = p_index + 1
              temp2 = temp1 * pad((md-1)*ma_max+ma) + &
                      temp3 * pac((mc-1)*ma_max+ma)
              temp4 = temp4 + temp2 * prim(p_index)
            END DO !ma
          END DO !mb
        END DO !mc
      END DO !md
    END IF

!$OMP ATOMIC
    force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1, &
          forces_map((coord-1)/3 + 1,2)) = &
    force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1, &
          forces_map((coord-1)/3 + 1,2)) - &
          temp4

  END SUBROUTINE update_forces

! *****************************************************************************
!> \brief This routine updates the virial using bufferd density matrices
!> \param ma_max, mb_max, mc_max, md_max: Size of matrix blocks
!> \param pbd, pbc, pad, pac buffers that will contain P(b,d), P(b,c) ...
!> \param fac mulitplication factor (spin, symmetry)
!> \param prim primitve forces
!> \param force storage loacation for forces
!> \param forces_map index table
!> \param coord which of the 12 coords to be updated
!> \par History
!>      03.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE update_virial(ma_max, mb_max, mc_max, md_max, &
                           pbd, pbc, pad, pac, fac, &
                           prim, tmp_virial, coord, l, &
                           pbd_mp2, pbc_mp2, pad_mp2, pac_mp2)

    INTEGER, INTENT(IN)                      :: ma_max, mb_max, mc_max, md_max
    REAL(dp), DIMENSION(*), INTENT(IN)       :: pbd, pbc, pad, pac
    REAL(dp), INTENT(IN)                     :: fac
    REAL(dp), DIMENSION(ma_max*mb_max*mc_max&
      *md_max), INTENT(IN)                   :: prim
    REAL(dp)                                 :: tmp_virial(3,3)
    INTEGER, INTENT(IN)                      :: coord, l
    REAL(dp), DIMENSION(*), INTENT(IN), &
      OPTIONAL                               :: pbd_mp2, pbc_mp2, pad_mp2, &
                                                pac_mp2

    INTEGER                                  :: i, j, ma, mb, mc, md, p_index
    LOGICAL                                  :: with_mp2_density
    REAL(dp)                                 :: temp1, temp1_mp2, temp2, &
                                                temp3, temp3_mp2, temp4

    with_mp2_density=.FALSE.
    IF(PRESENT(pbd_mp2).AND.&
       PRESENT(pbc_mp2).AND.&
       PRESENT(pad_mp2).AND.&
       PRESENT(pac_mp2)) with_mp2_density=.TRUE.

    IF(with_mp2_density) THEN
      p_index = 0
      temp4 = 0.0_dp
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,mb_max
            temp1 = pbc((mc-1)*mb_max+mb) * fac
            temp3 = pbd((md-1)*mb_max+mb) * fac
            temp1_mp2 = pbc_mp2((mc-1)*mb_max+mb) * fac
            temp3_mp2 = pbd_mp2((md-1)*mb_max+mb) * fac
            DO ma = 1,ma_max
              p_index = p_index + 1
              ! HF-SCF
              temp2 = temp1 * pad((md-1)*ma_max+ma) + &
                      temp3 * pac((mc-1)*ma_max+ma)
              ! MP2+HF
              temp2 = temp2 + &
                      pac((mc-1)*ma_max+ma) * temp3_mp2 + &
                      pac_mp2((mc-1)*ma_max+ma) * temp3 + &
                      pad((md-1)*ma_max+ma) * temp1_mp2 + &
                      pad_mp2((md-1)*ma_max+ma) * temp1
              temp4 = temp4 + temp2 * prim(p_index)
            END DO !ma
          END DO !mb
        END DO !mc
      END DO !md
    ELSE
      p_index = 0
      temp4 = 0.0_dp
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,mb_max
            temp1 = pbc((mc-1)*mb_max+mb) * fac
            temp3 = pbd((md-1)*mb_max+mb) * fac
            DO ma = 1,ma_max
              p_index = p_index + 1
              temp2 = temp1 * pad((md-1)*ma_max+ma) + &
                      temp3 * pac((mc-1)*ma_max+ma)
              temp4 = temp4 + temp2 * prim(p_index)
            END DO !ma
          END DO !mb
        END DO !mc
      END DO !md
    END IF

    j = l
    i = MOD(coord-1,3)+1
    tmp_virial(i,j) = tmp_virial(i,j) - temp4
  END SUBROUTINE update_virial

#include "hfx_get_pmax_val.f90"

END MODULE hfx_derivatives
