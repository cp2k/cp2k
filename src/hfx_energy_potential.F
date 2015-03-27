!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
                                             guidon2008,&
                                             guidon2009
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_get_matrix_type,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_trace,&
                                             dbcsr_type_antisymmetric
  USE cp_files,                        ONLY: get_unit_number
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE gamma,                           ONLY: init_md_ftable
  USE hfx_communication,               ONLY: distribute_ks_matrix,&
                                             get_atomic_block_maps,&
                                             get_full_density
  USE hfx_compression_methods,         ONLY: hfx_add_mult_cache_elements,&
                                             hfx_add_single_cache_element,&
                                             hfx_decompress_first_cache,&
                                             hfx_flush_last_cache,&
                                             hfx_get_mult_cache_elements,&
                                             hfx_get_single_cache_element,&
                                             hfx_reset_cache_and_container
  USE hfx_contract_block,              ONLY: contract_block
  USE hfx_libint_interface,            ONLY: evaluate_eri
  USE hfx_libint_wrapper_types,        ONLY: lib_int
  USE hfx_load_balance_methods,        ONLY: collect_load_balance_info,&
                                             hfx_load_balance,&
                                             hfx_update_load_balance
  USE hfx_pair_list_methods,           ONLY: build_atomic_pair_list,&
                                             build_pair_list,&
                                             build_pair_list_pgf,&
                                             build_pgf_product_list
  USE hfx_screening_methods,           ONLY: calc_pair_dist_radii,&
                                             calc_screening_functions,&
                                             update_pmax_mat
  USE hfx_types,                       ONLY: &
       alloc_containers, dealloc_containers, hfx_basis_info_type, &
       hfx_basis_type, hfx_cache_type, hfx_cell_type, hfx_container_type, &
       hfx_create_neighbor_cells, hfx_distribution, hfx_general_type, &
       hfx_init_container, hfx_load_balance_type, hfx_memory_type, &
       hfx_p_kind, hfx_pgf_list, hfx_pgf_product_list, hfx_potential_type, &
       hfx_reset_memory_usage_counter, hfx_screen_coeff_type, &
       hfx_screening_type, hfx_task_list_type, hfx_type, init_t_c_g0_lmax, &
       log_zero, pair_list_type, pair_set_list_type
  USE input_constants,                 ONLY: do_hfx_potential_mix_cl_trunc,&
                                             do_hfx_potential_truncated,&
                                             hfx_do_eval_energy
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: m_flush,&
                                             m_memory,&
                                             m_walltime
  USE mathconstants,                   ONLY: fac
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum,&
                                             mp_sync
  USE orbital_pointers,                ONLY: nco,&
                                             ncoset,&
                                             nso
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE t_c_g0,                          ONLY: init
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC ::  integrate_four_center, coulomb4

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_energy_potential'

!***

  CONTAINS

! *****************************************************************************
!> \brief computes four center integrals for a full basis set and updates the
!>      Kohn-Sham-Matrix and energy. Uses all 8 eri symmetries
!> \param qs_env ...
!> \param ks_matrix ...
!> \param energy ...
!> \param rho_ao density matrix in ao basis
!> \param hfx_section input_section HFX
!> \param para_env ...
!> \param geometry_did_change flag that indicates we have to recalc integrals
!> \param irep Index for the HFX replica
!> \param distribute_fock_matrix Flag that indicates whether to communicate the
!>        new fock matrix or not
!> \param ispin ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param do_exx ...
!> \par History
!>      06.2007 created [Manuel Guidon]
!>      08.2007 optimized load balance [Manuel Guidon]
!>      09.2007 new parallelization [Manuel Guidon]
!>      02.2009 completely rewritten screening part [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE integrate_four_center(qs_env,ks_matrix,energy,rho_ao,hfx_section,para_env,&
                                   geometry_did_change,irep,distribute_fock_matrix,&
                                   ispin, error, do_exx)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao
    TYPE(section_vals_type), POINTER         :: hfx_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL                                  :: geometry_did_change
    INTEGER                                  :: irep
    LOGICAL, INTENT(IN)                      :: distribute_fock_matrix
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL, OPTIONAL                        :: do_exx

    CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_four_center', &
      routineP = moduleN//':'//routineN

    INTEGER :: act_atomic_block_offset, act_set_offset, atomic_offset_ac, &
      atomic_offset_ad, atomic_offset_bc, atomic_offset_bd, bin, &
      bits_max_val, buffer_left, buffer_size, buffer_start, cache_size, &
      current_counter, handle, handle_bin, handle_dist_ks, handle_getP, &
      handle_load, handle_main, i, i_list_ij, i_list_kl, i_set_list_ij, &
      i_set_list_ij_start, i_set_list_ij_stop, i_set_list_kl, &
      i_set_list_kl_start, i_set_list_kl_stop, i_thread, iatom, iatom_block, &
      ikind, iset, iw, j, jatom, jatom_block, jkind, jset, katom, &
      katom_block, kind_kind_idx, kkind, kset, l_max, latom, latom_block, &
      lkind, lset, ma
    INTEGER :: max_am, max_pgf, max_set, mb, my_bin_size, n_threads, natom, &
      nbits, ncob, ncos_max, nints, nkind, nneighbors, nseta, nsetb, &
      nsgf_max, nspins, pa, sgfb, stat, swap_id, unit_id
    INTEGER(int_8) :: estimate_to_store_int, mem_compression_counter, &
      mem_compression_counter_disk, mem_eris, mem_eris_disk, mem_max_val, &
      my_current_counter, my_istart, n_processes, ncpu, neris_disk, &
      neris_incore, neris_onthefly, neris_total, shm_mem_compression_counter, &
      stor_count_int_disk, stor_count_max_val, storage_counter_integrals
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of, last_sgf_global
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, npgfa, npgfb, npgfc, npgfd, nsgfa, nsgfb, &
      nsgfc, nsgfd
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb, nsgfl_a, nsgfl_b, &
                                                nsgfl_c, nsgfl_d
    LOGICAL                                  :: buffer_overflow, &
                                                do_p_screening, do_periodic, &
                                                failure, my_do_exx
    REAL(dp) :: bintime_start, bintime_stop, cartesian_estimate, &
      compression_factor, compression_factor_disk, ene_x_aa, ene_x_bb, &
      eps_schwarz, eps_storage, fac, hf_fraction, max_contraction_val, &
      max_val1, max_val2, pmax_entry, ra(3), rb(3), rc(3), rd(3), &
      spherical_estimate, symm_fac
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: primitive_integrals
    REAL(dp), DIMENSION(:, :), POINTER       :: max_contraction, sphi_b, &
                                                zeta, zetb, zetc, zetd
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_fit_hfx
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
    TYPE(hfx_load_balance_type), POINTER     :: load_balance_parameter
    TYPE(hfx_memory_type), POINTER           :: memory_parameter
    TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(hfx_screening_type)                 :: screening_parameter
    TYPE(hfx_type), DIMENSION(:, :), POINTER :: x_data
    TYPE(hfx_type), POINTER                  :: actual_x_data, &
                                                shm_master_x_data
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!$  INTEGER :: omp_get_max_threads,omp_get_thread_num
    INTEGER, SAVE                            :: shm_number_of_p_entries
    INTEGER, DIMENSION(:,:), POINTER, SAVE   :: shm_is_assoc_atomic_block
    INTEGER(int_8)                           :: shm_neris_total, shm_neris_onthefly, &
                                                shm_storage_counter_integrals, shm_neris_incore, &
                                                shm_neris_disk, shm_stor_count_int_disk, shm_stor_count_max_val
    LOGICAL                                  :: do_disk_storage, use_disk_storage
    INTEGER(int_8)                           :: counter, neris_tmp, nprim_ints, shm_nprim_ints,&
                                                mb_size_p, mb_size_f, mb_size_buffers
    TYPE(cell_type), POINTER                 :: cell
    CHARACTER(LEN=512)                       :: error_msg
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:,:,:,:,:,:), POINTER        :: screen_coeffs_pgf, radii_pgf
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:,:,:,:), POINTER            :: screen_coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:,:), POINTER                :: screen_coeffs_kind
    REAL(dp)                                 :: rab2, rcd2, &
                                                log10_eps_schwarz, &
                                                log10_pmax, &
                                                ln_10
    INTEGER(int_8)                           :: subtr_size_mb

    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:,:), POINTER                :: tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2

    INTEGER(int_8)                           :: memsize_before,memsize_after,&
                                                max_val_memory, tmp_i8(8)
    REAL(dp)                                 :: max_val2_set
    TYPE(lib_int)                            :: private_lib
    REAL(dp), DIMENSION(:), POINTER          :: full_density, full_density_beta, full_ks,&
                                                full_ks_beta
    REAL(dp), DIMENSION(:), ALLOCATABLE      :: pbd_buf, pbc_buf, pad_buf, pac_buf,&
                                                kbd_buf, kbc_buf, kad_buf, kac_buf
    LOGICAL                                  :: treat_lsd_in_core, ks_fully_occ
    INTEGER(int_8)                           :: atom_block, tmp_block, nblocks
    TYPE(pair_list_type)                     :: list_ij, list_kl
    REAL(KIND=dp)                            :: coeffs_kind_max0
    TYPE(pair_set_list_type), DIMENSION(:), ALLOCATABLE :: set_list_ij, set_list_kl
    TYPE(hfx_p_kind), DIMENSION(:), POINTER  :: shm_initial_p
    REAL(dp), DIMENSION(:,:), POINTER        :: ptr_p_1, ptr_p_2, ptr_p_3, ptr_p_4
    INTEGER, DIMENSION(:), POINTER           :: shm_block_offset
    INTEGER, DIMENSION(:,:), POINTER         :: shm_atomic_block_offset
    INTEGER, DIMENSION(:,:,:,:), POINTER     :: shm_set_offset
    INTEGER, DIMENSION(:,:), POINTER         :: offset_bd_set, offset_bc_set, offset_ad_set, offset_ac_set
    REAL(dp)                                 :: ene_x_aa_diag, ene_x_bb_diag
    TYPE(hfx_pgf_list), DIMENSION(:), &
                            ALLOCATABLE      :: pgf_list_ij, pgf_list_kl
    TYPE(hfx_pgf_product_list), DIMENSION(:), &
                        ALLOCATABLE          :: pgf_product_list
    REAL(dp), DIMENSION(:,:,:,:), POINTER &
                                             :: sphi_a_ext, sphi_b_ext, sphi_c_ext, sphi_d_ext
    REAL(dp), DIMENSION(:,:,:), POINTER &
                                             :: sphi_a_ext_set, sphi_b_ext_set, sphi_c_ext_set, sphi_d_ext_set
    LOGICAL                                  :: do_dynamic_load_balancing, bins_left, do_it
    TYPE(hfx_distribution), POINTER          :: distribution_energy
    INTEGER                                  :: shm_total_bins, shm_task_counter, my_thread_id, my_bin_id
    TYPE(hfx_task_list_type), DIMENSION(:), &
                              POINTER        :: shm_task_list, tmp_task_list
    INTEGER(int_8), DIMENSION(:), &
                    ALLOCATABLE              :: tmp_task_list_cost
    INTEGER, DIMENSION(:), ALLOCATABLE       :: tmp_index
    REAL(dp), DIMENSION(:), ALLOCATABLE      :: ee_work, ee_work2, &
                                                ee_buffer1, ee_buffer2, &
                                                ee_primitives_tmp
    INTEGER, DIMENSION(:), ALLOCATABLE       :: nimages
    INTEGER                                  :: sphi_d_u1, sphi_d_u2, sphi_d_u3
    INTEGER                                  :: sphi_c_u1, sphi_c_u2, sphi_c_u3
    INTEGER                                  :: sphi_b_u1, sphi_b_u2, sphi_b_u3
    INTEGER                                  :: sphi_a_u1, sphi_a_u2, sphi_a_u3
    INTEGER                                  :: iatom_start, iatom_end, jatom_start, jatom_end, katom_start,&
                                                katom_end, latom_start, latom_end
    REAL(dp), DIMENSION(:,:), POINTER        :: shm_pmax_atom, shm_pmax_block
    REAL(dp)                                 :: pmax_atom, screen_kind_ij, screen_kind_kl, pmax_blocks,afac
    REAL(dp), DIMENSION(:), POINTER          :: p_work
    LOGICAL, DIMENSION(:,:), POINTER         :: shm_atomic_pair_list
    LOGICAL                                  :: do_print_load_balance_info,is_anti_symmetric
    TYPE(dft_control_type),  &
      POINTER                                :: dft_control

    NULLIFY(dft_control)

#if defined(__USE_PAT)
    INTEGER :: istat
#include "pat_apif.h"
#endif

    CALL timeset(routineN,handle)

#if defined(__USE_PAT)
  CALL PAT_sampling_state(PAT_STATE_ON,istat)
#endif

    CALL cite_reference(Guidon2008)
    CALL cite_reference(Guidon2009)

    ! ** This is not very clean, but effective. ispin can only be 2 if we do the beta spin part in core
    IF( ispin == 2 ) geometry_did_change = .FALSE.

    logger => cp_error_get_logger(error)

    is_anti_symmetric=cp_dbcsr_get_matrix_type(rho_ao(1)%matrix).EQ.dbcsr_type_antisymmetric

    IF( geometry_did_change ) THEN
      memsize_before=m_memory()
      CALL mp_max(memsize_before, para_env%group)
      iw = cp_print_key_unit_nr(logger,hfx_section,"HF_INFO",&
                                extension=".scfLog",error=error)
      IF (iw>0) THEN
        WRITE (UNIT=iw,FMT="(/,(T3,A,T60,I21))")&
              "HFX_MEM_INFO| Est. max. program size before HFX [MB's]:", memsize_before/(1024*1024)
        CALL m_flush(iw)
      END IF
      CALL cp_print_key_finished_output(iw,logger,hfx_section,&
                                       "HF_INFO", error=error)
    END IF

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                    matrix_ks_aux_fit_hfx=matrix_ks_aux_fit_hfx, error=error)

    my_do_exx=.FALSE.
    IF(PRESENT(do_exx)) my_do_exx=do_exx
    NULLIFY(x_data)
    IF(my_do_exx) THEN
      x_data=>qs_env%mp2_env%ri_rpa%x_data
    ELSE
      x_data=>qs_env%x_data
    END IF

    !! Calculate l_max used in fgamma , because init_md_ftable is definitely not thread safe
    nkind = SIZE(atomic_kind_set,1)
    l_max = 0
    DO ikind=1,nkind
      l_max = MAX(l_max,MAXVAL(x_data(1,1)%basis_parameter(ikind)%lmax))
    ENDDO
    l_max = 4*l_max
    CALL init_md_ftable(l_max)

    IF(x_data(1,1)%potential_parameter%potential_type == do_hfx_potential_truncated .OR. &
       x_data(1,1)%potential_parameter%potential_type == do_hfx_potential_mix_cl_trunc) THEN
      IF(l_max>init_t_c_g0_lmax) THEN
        IF( para_env%mepos == 0 ) THEN
          unit_id = get_unit_number()
          OPEN(unit_id,FILE=x_data(1,1)%potential_parameter%filename)
        END IF
        CALL init(l_max, unit_id, para_env%mepos, para_env%group)
        IF(para_env%mepos == 0 ) THEN
          CLOSE(unit_id)
        END IF
        init_t_c_g0_lmax = l_max
      END IF
    END IF


    n_threads = 1
!$  n_threads = omp_get_max_threads()

    shm_neris_total = 0
    shm_nprim_ints = 0
    shm_neris_onthefly = 0
    shm_storage_counter_integrals = 0
    shm_stor_count_int_disk = 0
    shm_neris_incore = 0
    shm_neris_disk = 0
    shm_stor_count_max_val = 0

!$OMP PARALLEL DEFAULT(NONE) SHARED(qs_env,&
!$OMP                                  x_data,&
!$OMP                                  ks_matrix,&
!$OMP                                  energy,&
!$OMP                                  rho_ao,&
!$OMP                                  matrix_ks_aux_fit_hfx,&
!$OMP                                  hfx_section,&
!$OMP                                  para_env,&
!$OMP                                  geometry_did_change,&
!$OMP                                  irep,&
!$OMP                                  distribute_fock_matrix,&
!$OMP                                  error,&
!$OMP                                  ncoset,&
!$OMP                                  nso,&
!$OMP                                  nco,&
!$OMP                                  full_ks,&
!$OMP                                  full_ks_beta,&
!$OMP                                  n_threads,&
!$OMP                                  full_density,&
!$OMP                                  full_density_beta,&
!$OMP                                  shm_initial_p,&
!$OMP                                  shm_is_assoc_atomic_block,&
!$OMP                                  shm_number_of_p_entries,&
!$OMP                                  shm_neris_total,&
!$OMP                                  shm_neris_onthefly,&
!$OMP                                  shm_storage_counter_integrals,&
!$OMP                                  shm_stor_count_int_disk,&
!$OMP                                  shm_neris_incore,&
!$OMP                                  shm_neris_disk,&
!$OMP                                  shm_nprim_ints,&
!$OMP                                  shm_stor_count_max_val,&
!$OMP                                  cell,&
!$OMP                                  screen_coeffs_set,&
!$OMP                                  screen_coeffs_kind,&
!$OMP                                  screen_coeffs_pgf,&
!$OMP                                  radii_pgf,&
!$OMP                                  nkind,&
!$OMP                                  ispin,&
!$OMP                                  is_anti_symmetric,&
!$OMP                                  shm_atomic_block_offset,&
!$OMP                                  shm_set_offset,&
!$OMP                                  shm_block_offset,&
!$OMP                                  shm_task_counter,&
!$OMP                                  shm_task_list,&
!$OMP                                  shm_total_bins,&
!$OMP                                  shm_master_x_data,&
!$OMP                                  shm_pmax_atom,&
!$OMP                                  shm_pmax_block,&
!$OMP                                  shm_atomic_pair_list,&
!$OMP                                  shm_mem_compression_counter,&
!$OMP                                  do_print_load_balance_info) &
!$OMP PRIVATE(ln_10,i_thread,actual_x_data,do_periodic,screening_parameter,potential_parameter,&
!$OMP         general_parameter,load_balance_parameter,memory_parameter,cache_size,bits_max_val,&
!$OMP         basis_parameter,basis_info,treat_lsd_in_core,ncpu,n_processes,neris_total,neris_incore,&
!$OMP         neris_disk,neris_onthefly,mem_eris,mem_eris_disk,mem_max_val,compression_factor,&
!$OMP         compression_factor_disk,nprim_ints,neris_tmp,max_val_memory,max_am,do_p_screening,&
!$OMP         max_set,particle_set,atomic_kind_set,natom,kind_of,stat,failure,ncos_max,nsgf_max,ikind,&
!$OMP         nseta,npgfa,la_max,nsgfa,primitive_integrals,pbd_buf,pbc_buf,pad_buf,pac_buf,kbd_buf,kbc_buf,&
!$OMP         kad_buf,kac_buf,ee_work,ee_work2,ee_buffer1,ee_buffer2,ee_primitives_tmp,nspins,max_contraction,&
!$OMP         max_pgf,jkind,lb_max,nsetb,npgfb,first_sgfb,sphi_b,nsgfb,ncob,sgfb,nneighbors,pgf_list_ij,pgf_list_kl,&
!$OMP         pgf_product_list,nimages,ks_fully_occ,error_msg,subtr_size_mb,use_disk_storage,counter,do_disk_storage,&
!$OMP         maxval_container_disk,maxval_cache_disk,integral_containers_disk,integral_caches_disk,eps_schwarz,&
!$OMP         log10_eps_schwarz,eps_storage,hf_fraction,buffer_overflow,logger,private_lib,last_sgf_global,handle_getp,&
!$OMP         p_work,fac,handle_load,do_dynamic_load_balancing,my_bin_size,maxval_container,integral_containers,maxval_cache,&
!$OMP         integral_caches,tmp_task_list,tmp_task_list_cost,tmp_index,handle_main,coeffs_kind_max0,set_list_ij,&
!$OMP         set_list_kl,iatom_start,iatom_end,jatom_start,jatom_end,nblocks,bins_left,do_it,distribution_energy,&
!$OMP         my_thread_id,my_bin_id,handle_bin,bintime_start,my_istart,my_current_counter,latom_block,tmp_block,&
!$OMP         katom_block,katom_start,katom_end,latom_start,latom_end,pmax_blocks,list_ij,list_kl,i_set_list_ij_start,&
!$OMP         i_set_list_ij_stop,ra,rb,rab2,la_min,zeta,sphi_a_ext,nsgfl_a,sphi_a_u1,sphi_a_u2,sphi_a_u3,&
!$OMP         lb_min,zetb,sphi_b_ext,nsgfl_b,sphi_b_u1,sphi_b_u2,sphi_b_u3,katom,latom,i_set_list_kl_start,i_set_list_kl_stop,&
!$OMP         kkind,lkind,rc,rd,rcd2,pmax_atom,screen_kind_ij,screen_kind_kl,symm_fac,lc_max,lc_min,npgfc,zetc,nsgfc,sphi_c_ext,&
!$OMP         nsgfl_c,sphi_c_u1,sphi_c_u2,sphi_c_u3,ld_max,ld_min,npgfd,zetd,nsgfd,sphi_d_ext,nsgfl_d,sphi_d_u1,sphi_d_u2,&
!$OMP         sphi_d_u3,atomic_offset_bd,atomic_offset_bc,atomic_offset_ad,atomic_offset_ac,offset_bd_set,offset_bc_set,&
!$OMP         offset_ad_set,offset_ac_set,swap_id,kind_kind_idx,ptr_p_1,ptr_p_2,ptr_p_3,ptr_p_4,mem_compression_counter,&
!$OMP         mem_compression_counter_disk,max_val1,sphi_a_ext_set,sphi_b_ext_set,kset,lset,max_val2_set,max_val2,&
!$OMP         sphi_c_ext_set,sphi_d_ext_set,pmax_entry,log10_pmax,current_counter,nints,estimate_to_store_int,&
!$OMP         spherical_estimate,nbits,buffer_left,buffer_start,buffer_size,max_contraction_val,tmp_r_1,tmp_r_2,&
!$OMP         tmp_screen_pgf1,tmp_screen_pgf2,cartesian_estimate,bintime_stop,iw,memsize_after,storage_counter_integrals,&
!$OMP         stor_count_int_disk,stor_count_max_val,ene_x_aa,ene_x_bb,mb_size_p,mb_size_f,mb_size_buffers,afac,ene_x_aa_diag,&
!$OMP         ene_x_bb_diag,act_atomic_block_offset,act_set_offset,j,handle_dist_ks,tmp_i8,dft_control,bin)


    ln_10 = LOG(10.0_dp)
    i_thread = 0
!$  i_thread = omp_get_thread_num()

    actual_x_data => x_data(irep, i_thread + 1)
!$OMP MASTER
    shm_master_x_data => x_data(irep,1)
!$OMP END MASTER
!$OMP BARRIER



    do_periodic = actual_x_data%periodic_parameter%do_periodic

    IF( do_periodic ) THEN
      ! ** Rebuild neighbor lists in case the cell has changed (i.e. NPT MD)
      actual_x_data%periodic_parameter%number_of_shells = actual_x_data%periodic_parameter%mode
      CALL hfx_create_neighbor_cells(actual_x_data, actual_x_data%periodic_parameter%number_of_shells_from_input,&
                                     cell, i_thread, &
                                     error)
    END IF


    screening_parameter    = actual_x_data%screening_parameter
    potential_parameter    = actual_x_data%potential_parameter

    general_parameter      = actual_x_data%general_parameter
    load_balance_parameter => actual_x_data%load_balance_parameter
    memory_parameter       => actual_x_data%memory_parameter



    cache_size = memory_parameter%cache_size
    bits_max_val = memory_parameter%bits_max_val

    basis_parameter   => actual_x_data%basis_parameter
    basis_info        => actual_x_data%basis_info

    treat_lsd_in_core = general_parameter%treat_lsd_in_core

    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads

    !! initalize some counters
    neris_total = 0_int_8
    neris_incore = 0_int_8
    neris_disk = 0_int_8
    neris_onthefly = 0_int_8
    mem_eris = 0_int_8
    mem_eris_disk = 0_int_8
    mem_max_val = 0_int_8
    compression_factor = 0.0_dp
    compression_factor_disk = 0.0_dp
    nprim_ints = 0_int_8
    neris_tmp = 0_int_8
    max_val_memory = 1_int_8

    max_am = basis_info%max_am

    do_p_screening = screening_parameter%do_initial_p_screening
    max_set = basis_info%max_set
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    dft_control=dft_control,&
                    error=error)

    natom = SIZE(particle_set,1)

    ALLOCATE(kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             kind_of=kind_of)

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
         nsgf_max = MAX(nsgf_max,nsgfa(iset))
      ENDDO
    ENDDO
    !! Allocate the arrays for the integrals.
    ALLOCATE(primitive_integrals(nsgf_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    primitive_integrals = 0.0_dp

    ALLOCATE(pbd_buf(nsgf_max**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(pbc_buf(nsgf_max**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(pad_buf(nsgf_max**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(pac_buf(nsgf_max**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(kbd_buf(nsgf_max**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(kbc_buf(nsgf_max**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(kad_buf(nsgf_max**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(kac_buf(nsgf_max**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ee_work(ncos_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ee_work2(ncos_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ee_buffer1(ncos_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ee_buffer2(ncos_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ee_primitives_tmp(nsgf_max**4),STAT=stat) ! XXXXX could be wrong
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    nspins = dft_control%nspins

    ALLOCATE(max_contraction(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    max_contraction=0.0_dp
    max_pgf = 0
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
        max_pgf = MAX(max_pgf,npgfb(jset))
      ENDDO
    ENDDO

    ! ** Allocate buffers for pgf_lists
    nneighbors = SIZE(actual_x_data%neighbor_cells)
    ALLOCATE(pgf_list_ij(max_pgf**2), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(pgf_list_kl(max_pgf**2), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(pgf_product_list(nneighbors**3), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(nimages(max_pgf**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO i=1,max_pgf**2
      ALLOCATE(pgf_list_ij(i)%image_list(nneighbors), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(pgf_list_kl(i)%image_list(nneighbors), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
!$OMP BARRIER
!$OMP MASTER
    !! Calculate helper array that stores if a certain atomic pair is associated in the KS matrix
    IF( geometry_did_change ) THEN
      CALL get_atomic_block_maps(ks_matrix(1)%matrix,basis_parameter,kind_of,&
                                 shm_master_x_data%is_assoc_atomic_block,&
                                 shm_master_x_data%number_of_p_entries ,&
                                 para_env, &
                                 shm_master_x_data%atomic_block_offset,&
                                 shm_master_x_data%set_offset,&
                                 shm_master_x_data%block_offset,&
                                 shm_master_x_data%map_atoms_to_cpus,&
                                 nkind,&
                                 error)

      shm_is_assoc_atomic_block => shm_master_x_data%is_assoc_atomic_block

      !! Get occupation of KS-matrix
      ks_fully_occ = .TRUE.
      outer: DO iatom=1,natom
        DO jatom=iatom,natom
          IF(shm_is_assoc_atomic_block(jatom,iatom) == 0 ) THEN
            ks_fully_occ = .FALSE.
            EXIT outer
          END IF
        END DO
      END DO outer

      IF ( .NOT. ks_fully_occ) THEN
        WRITE(error_msg,'(A,A,A)') "The Kohn Sham matrix is not 100% occupied. This "//&
             "may result in uncorrect Hartree-Fock results. Try to decrease EPS_PGF_ORB "//&
             "and EPS_FILTER_MATRIX in the QS section. "
        CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                       error_msg,&
                       error,&
                       only_ionode=.TRUE.)
      END IF
    END IF

    ! ** Set pointers
    shm_number_of_p_entries = shm_master_x_data%number_of_p_entries
    shm_is_assoc_atomic_block => shm_master_x_data%is_assoc_atomic_block
    shm_atomic_block_offset => shm_master_x_data%atomic_block_offset
    shm_set_offset => shm_master_x_data%set_offset
    shm_block_offset => shm_master_x_data%block_offset
!$OMP END MASTER
!$OMP BARRIER

    ! ** Reset storage counter given by MAX_MEMORY by subtracting all buffers
    ! ** Fock and density Matrices (shared)
    subtr_size_mb = 2_int_8*shm_block_offset(ncpu+1)
    ! ** if non restricted ==> alpha/beta spin
    IF (.NOT. treat_lsd_in_core ) THEN
      IF( nspins == 2) subtr_size_mb = subtr_size_mb * 2_int_8
    END IF
    ! ** Initial P only MAX(alpha,beta) (shared)
    IF(do_p_screening .OR. screening_parameter%do_p_screening_forces) THEN
      subtr_size_mb = subtr_size_mb + memory_parameter%size_p_screen
    END IF
    ! ** In core forces require their own initial P
    IF(screening_parameter%do_p_screening_forces) THEN
      IF( memory_parameter%treat_forces_in_core) THEN
        subtr_size_mb = subtr_size_mb + memory_parameter%size_p_screen
      END IF
    END IF
    ! ** primitive buffer (not shared by the threads)
    subtr_size_mb = subtr_size_mb + nsgf_max**4*n_threads
    ! ** density + fock buffers
    subtr_size_mb = subtr_size_mb + 8_int_8*nsgf_max**2*n_threads
    ! ** screening functions (shared)
    ! ** coeffs_pgf
    subtr_size_mb = subtr_size_mb + max_pgf**2*max_set**2*nkind**2
    ! ** coeffs_set
    subtr_size_mb = subtr_size_mb + max_set**2*nkind**2
    ! ** coeffs_kind
    subtr_size_mb = subtr_size_mb + nkind**2
    ! ** radii_pgf
    subtr_size_mb = subtr_size_mb + max_pgf**2*max_set**2*nkind**2

    ! ** is_assoc (shared)
    subtr_size_mb = subtr_size_mb + natom**2

    ! ** pmax_atom (shared)
    IF( do_p_screening) THEN
      subtr_size_mb = subtr_size_mb + natom**2
    END IF
    IF(screening_parameter%do_p_screening_forces) THEN
      IF( memory_parameter%treat_forces_in_core) THEN
        subtr_size_mb = subtr_size_mb + natom**2
      END IF
    END IF

    ! ** Convert into MB's
    subtr_size_mb = subtr_size_mb * 8_int_8/1024_int_8/1024_int_8

    ! ** Subtracting all these buffers from from MAX_MEMORY yields the amount
    ! ** of RAM that is left for the compressed integrals. When using threads
    ! ** all the available memory is shared among all n_threads. i.e. the faster
    ! ** ones can steal from the slower ones

    CALL hfx_reset_memory_usage_counter(memory_parameter, subtr_size_mb)

    use_disk_storage = .FALSE.
    counter = 0_int_8
    do_disk_storage = memory_parameter%do_disk_storage
    IF (do_disk_storage ) THEN
      maxval_container_disk       => actual_x_data%maxval_container_disk
      maxval_cache_disk           => actual_x_data%maxval_cache_disk

      integral_containers_disk    => actual_x_data%integral_containers_disk
      integral_caches_disk        => actual_x_data%integral_caches_disk
    END IF


    IF(geometry_did_change) THEN
     memory_parameter%ram_counter = HUGE(memory_parameter%ram_counter)
    END IF

    IF( geometry_did_change) THEN
      memory_parameter%recalc_forces = .TRUE.
    ELSE
      IF( .NOT. memory_parameter%treat_forces_in_core ) memory_parameter%recalc_forces = .TRUE.
    END IF

    !! Get screening parameter
    eps_schwarz = screening_parameter%eps_schwarz
    IF( eps_schwarz <= 0.0_dp) THEN
      log10_eps_schwarz = log_zero
    ELSE
      log10_eps_schwarz = LOG10(eps_schwarz)
    END IF
    !! get storage epsilon
    eps_storage = eps_schwarz*memory_parameter%eps_storage_scaling


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
    logger => cp_error_get_logger(error)


    private_lib = actual_x_data%lib

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
    !! Let master thread get the density (avoid problems with MPI)
    !! Get the full density from all the processors
    NULLIFY(full_density)
    NULLIFY(full_density_beta)
    ALLOCATE(full_density(shm_block_offset(ncpu+1)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF( .NOT. treat_lsd_in_core .OR. nspins == 1 ) THEN
      CALL timeset(routineN//"_getP",handle_getP)
      CALL get_full_density(para_env, full_density, rho_ao(ispin)%matrix, shm_number_of_p_entries,&
                            shm_master_x_data%block_offset, natom, &
                            kind_of, basis_parameter, get_max_vals_spin=.FALSE., antisymmetric=is_anti_symmetric,& 
                            error=error)

      IF(nspins == 2) THEN
        ALLOCATE(full_density_beta(shm_block_offset(ncpu+1)),STAT=stat)
        CALL get_full_density(para_env, full_density_beta, rho_ao(2)%matrix, shm_number_of_p_entries,&
                              shm_master_x_data%block_offset, natom,  &
                              kind_of, basis_parameter, get_max_vals_spin=.FALSE.,antisymmetric=is_anti_symmetric,&
                              error=error)
      END IF
      CALL timestop(handle_getP)

      !! Calculate the max values of the density matrix actual_pmax stores the data from the actual density matrix
      !! and x_data%initial_p stores the same for the initial guess. The initial guess is updated only in the case of
      !! a changed geometry
      NULLIFY(shm_initial_p)
      IF( do_p_screening) THEN
        shm_initial_p => shm_master_x_data%initial_p
        shm_pmax_atom => shm_master_x_data%pmax_atom
        IF(geometry_did_change ) THEN
          CALL update_pmax_mat(shm_master_x_data%initial_p, &
                               shm_master_x_data%map_atom_to_kind_atom, &
                               shm_master_x_data%set_offset,&
                               shm_master_x_data%atomic_block_offset,&
                               shm_pmax_atom,&
                               full_density,full_density_beta,&
                               natom, kind_of, basis_parameter, &
                               nkind, shm_master_x_data%is_assoc_atomic_block)
        END IF
      END IF
    ELSE
      IF( do_p_screening ) THEN
        CALL timeset(routineN//"_getP",handle_getP)
        CALL get_full_density(para_env, full_density, rho_ao(1)%matrix, shm_number_of_p_entries,&
                              shm_master_x_data%block_offset, natom,  &
                              kind_of, basis_parameter,  get_max_vals_spin=.TRUE., &
                              rho_beta=rho_ao(2)%matrix, antisymmetric=is_anti_symmetric, error=error)
        CALL timestop(handle_getP)

        !! Calculate the max values of the density matrix actual_pmax stores the data from the actual density matrix
        !! and x_data%initial_p stores the same for the initial guess. The initial guess is updated only in the case of
        !! a changed geometry
        NULLIFY(shm_initial_p)
        shm_initial_p => actual_x_data%initial_p
        shm_pmax_atom => shm_master_x_data%pmax_atom
        IF(geometry_did_change ) THEN
          CALL update_pmax_mat(shm_master_x_data%initial_p, &
                               shm_master_x_data%map_atom_to_kind_atom, &
                               shm_master_x_data%set_offset,&
                               shm_master_x_data%atomic_block_offset,&
                               shm_pmax_atom,&
                               full_density,full_density_beta,&
                               natom, kind_of, basis_parameter, &
                               nkind, shm_master_x_data%is_assoc_atomic_block)
        END IF
      END IF
      ! ** Now get the density(ispin)
      CALL get_full_density(para_env, full_density, rho_ao(ispin)%matrix, shm_number_of_p_entries,&
                            shm_master_x_data%block_offset, natom, &
                            kind_of, basis_parameter, get_max_vals_spin=.FALSE., antisymmetric=is_anti_symmetric, error=error)
    END IF

    NULLIFY(full_ks, full_ks_beta)
    ALLOCATE(shm_master_x_data%full_ks_alpha(shm_block_offset(ncpu+1)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    full_ks => shm_master_x_data%full_ks_alpha
    full_ks = 0.0_dp


    IF (.NOT. treat_lsd_in_core ) THEN
      IF(nspins==2) THEN
        ALLOCATE(shm_master_x_data%full_ks_beta(shm_block_offset(ncpu+1)),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        full_ks_beta => shm_master_x_data%full_ks_beta
        full_ks_beta = 0.0_dp
      END IF
    END IF

    !! Initialize schwarz screening matrices for near field estimates and boxing screening matrices
    !! for far field estimates. The update is only performed if the geomtry of the system changed.
    !! If the system is periodic, then the corresponding routines are called and some variables
    !! are initialized

!$OMP END MASTER
!$OMP BARRIER


    IF( .NOT. shm_master_x_data%screen_funct_is_initialized ) THEN
      CALL calc_pair_dist_radii(qs_env, basis_parameter,&
                                shm_master_x_data%pair_dist_radii_pgf, max_set, max_pgf, eps_schwarz,&
                                n_threads, i_thread, error)
!$OMP BARRIER
      CALL calc_screening_functions(qs_env, basis_parameter, private_lib, shm_master_x_data%potential_parameter,&
                                    shm_master_x_data%screen_funct_coeffs_set,&
                                    shm_master_x_data%screen_funct_coeffs_kind, &
                                    shm_master_x_data%screen_funct_coeffs_pgf, &
                                    shm_master_x_data%pair_dist_radii_pgf,&
                                    max_set, max_pgf, n_threads, i_thread, p_work, error)

!$OMP MASTER
      shm_master_x_data%screen_funct_is_initialized = .TRUE.
!$OMP END MASTER
   END IF
!$OMP BARRIER

!$OMP MASTER
    screen_coeffs_set   => shm_master_x_data%screen_funct_coeffs_set
    screen_coeffs_kind  => shm_master_x_data%screen_funct_coeffs_kind
    screen_coeffs_pgf   => shm_master_x_data%screen_funct_coeffs_pgf
    radii_pgf           => shm_master_x_data%pair_dist_radii_pgf
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
!$OMP BARRIER
!$OMP MASTER
    CALL timeset(routineN//"_load",handle_load)
!$OMP END MASTER
!$OMP BARRIER
    IF( geometry_did_change ) THEN
      IF( actual_x_data%b_first_load_balance_energy ) THEN
        CALL hfx_load_balance(actual_x_data,eps_schwarz,particle_set,max_set,para_env,&
                              screen_coeffs_set,screen_coeffs_kind,&
                              shm_is_assoc_atomic_block,do_periodic,load_balance_parameter, &
                              kind_of, basis_parameter, shm_initial_p, shm_pmax_atom, i_thread, n_threads, &
                              cell, do_p_screening, actual_x_data%map_atom_to_kind_atom, &
                              nkind, hfx_do_eval_energy, shm_pmax_block, use_virial=.FALSE., error=error)
        actual_x_data%b_first_load_balance_energy = .FALSE.
      ELSE
        CALL hfx_update_load_balance(actual_x_data,para_env, &
                                     load_balance_parameter, &
                                     i_thread, n_threads, hfx_do_eval_energy, &
                                     error)
      END IF
    END IF
!$OMP BARRIER
!$OMP MASTER
    CALL timestop(handle_load)
!$OMP END MASTER
!$OMP BARRIER

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
    !! istart, jstart, kstart, lstart, number_of_atom_quartets, initial_cpu_id
    !!
    !! This tells the process where to start the main loops and how many bunches of integrals it has to
    !! calculate. The original parallelization is a simple modulo distribution that is binned and
    !! optimized in the load_balance routines. Since the Monte Carlo routines can swap processors,
    !! we need to know which was the inital cpu_id.
    !! Furthermore, the indices iatom, jatom, katom, latom have to be set to istart, jstart, kstart and
    !! lstart only the first time the loop is executed. All subsequent loops have to start with one or
    !! iatom and katom respectively. Therefore, we use flags like first_j_loop etc.


    do_dynamic_load_balancing = .TRUE.

    IF( n_threads == 1 .OR. do_disk_storage ) do_dynamic_load_balancing = .FALSE.

    IF( do_dynamic_load_balancing ) THEN
      my_bin_size = SIZE(actual_x_data%distribution_energy)
    ELSE
      my_bin_size = 1
    END IF
    !! We do not need the containers if MAX_MEM = 0
    IF(.NOT. memory_parameter%do_all_on_the_fly) THEN
      !! IF new md step -> reinitialize containers
      IF(geometry_did_change) THEN
        CALL dealloc_containers(actual_x_data, hfx_do_eval_energy, error)
        CALL alloc_containers(actual_x_data, my_bin_size, hfx_do_eval_energy, error)

        DO bin = 1,my_bin_size
          maxval_container => actual_x_data%maxval_container(bin)
          integral_containers => actual_x_data%integral_containers(:,bin)
          CALL hfx_init_container(maxval_container, memory_parameter%actual_memory_usage, .FALSE., error)
          DO i=1,64
            CALL hfx_init_container(integral_containers(i),memory_parameter%actual_memory_usage, .FALSE., error)
          END DO
        END DO
      END IF

      !! Decompress the first cache for maxvals and integrals
      IF( .NOT. geometry_did_change) THEN
        DO bin = 1,my_bin_size
          maxval_cache => actual_x_data%maxval_cache(bin)
          maxval_container => actual_x_data%maxval_container(bin)
          integral_caches => actual_x_data%integral_caches(:,bin)
          integral_containers => actual_x_data%integral_containers(:,bin)
          CALL hfx_decompress_first_cache(bits_max_val, maxval_cache, maxval_container, &
                                          memory_parameter%actual_memory_usage, .FALSE.)
          DO i=1,64
            CALL hfx_decompress_first_cache(i, integral_caches(i), integral_containers(i), &
                                            memory_parameter%actual_memory_usage,.FALSE.)
          END DO
        END DO
      END IF
    END IF

    !! Since the I/O routines are no thread-safe, i.e. the procedure to get the unit number, put a lock here
!$OMP CRITICAL
    !! If we do disk storage, we have to initialize the containers/caches
    IF(do_disk_storage) THEN
      !! IF new md step -> reinitialize containers
      IF(geometry_did_change) THEN
        CALL hfx_init_container(maxval_container_disk, memory_parameter%actual_memory_usage_disk, do_disk_storage, error)
        DO i=1,64
          CALL hfx_init_container(integral_containers_disk(i),memory_parameter%actual_memory_usage_disk, do_disk_storage, error)
        END DO
      END IF
      !! Decompress the first cache for maxvals and integrals
      IF( .NOT. geometry_did_change) THEN
        CALL hfx_decompress_first_cache(bits_max_val, maxval_cache_disk, maxval_container_disk, &
                                         memory_parameter%actual_memory_usage_disk, .TRUE.)
        DO i=1,64
          CALL hfx_decompress_first_cache(i, integral_caches_disk(i), integral_containers_disk(i), &
                                           memory_parameter%actual_memory_usage_disk, .TRUE.)
        END DO
      END IF
    END IF
!$OMP END CRITICAL

!$OMP BARRIER
!$OMP MASTER

    IF( do_dynamic_load_balancing ) THEN
      ! ** Lets contstruct the task list
      shm_total_bins = 0
      DO i = 1,n_threads
        shm_total_bins = shm_total_bins + SIZE(x_data(irep, i)%distribution_energy)
      END DO
      ALLOCATE(tmp_task_list(shm_total_bins), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      shm_task_counter = 0
      DO i = 1,n_threads
        DO bin = 1,SIZE(x_data(irep,i)%distribution_energy)
          shm_task_counter = shm_task_counter + 1
          tmp_task_list(shm_task_counter)%thread_id = i
          tmp_task_list(shm_task_counter)%bin_id = bin
          tmp_task_list(shm_task_counter)%cost = x_data(irep,i)%distribution_energy(bin)%cost
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

      ALLOCATE(shm_master_x_data%task_list(shm_total_bins), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO i = 1,shm_total_bins
        shm_master_x_data%task_list(i) = tmp_task_list(tmp_index(shm_total_bins - i +1))
      END DO

      shm_task_list => shm_master_x_data%task_list
      shm_task_counter = 0

      DEALLOCATE(tmp_task_list_cost, tmp_index, tmp_task_list, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
!$OMP END MASTER
!$OMP BARRIER

    IF( my_bin_size > 0 ) THEN
      maxval_container       => actual_x_data%maxval_container(1)
      maxval_cache           => actual_x_data%maxval_cache(1)

      integral_containers    => actual_x_data%integral_containers(:,1)
      integral_caches        => actual_x_data%integral_caches(:,1)
    END IF

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
!$OMP MASTER

    !! precalculate maximum density matrix elements in blocks
    actual_x_data%pmax_block = 0.0_dp
    shm_pmax_block => actual_x_data%pmax_block
    IF( do_p_screening ) THEN
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
    shm_atomic_pair_list => actual_x_data%atomic_pair_list
    IF( geometry_did_change) THEN
      CALL build_atomic_pair_list(natom, shm_atomic_pair_list, kind_of,basis_parameter,particle_set, &
                                  do_periodic,screen_coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell,&
                                  actual_x_data%blocks)
    END IF

    my_bin_size = SIZE(actual_x_data%distribution_energy)
    ! reset timings for the new SCF round
    IF (geometry_did_change) THEN
      DO bin=1,my_bin_size
        actual_x_data%distribution_energy(bin)%time_first_scf=0.0_dp
        actual_x_data%distribution_energy(bin)%time_other_scf=0.0_dp
        actual_x_data%distribution_energy(bin)%time_forces=0.0_dp
       END DO
    ENDIF
!$OMP END MASTER
!$OMP BARRIER

    my_bin_size = SIZE(actual_x_data%distribution_energy)
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
          distribution_energy => actual_x_data%distribution_energy(bin)
        END IF
      ELSE
!$OMP CRITICAL
        shm_task_counter = shm_task_counter + 1
        IF( shm_task_counter <= shm_total_bins) THEN
          my_thread_id = shm_task_list(shm_task_counter)%thread_id
          my_bin_id = shm_task_list(shm_task_counter)%bin_id
          maxval_cache => x_data(irep,my_thread_id)%maxval_cache(my_bin_id)
          maxval_container => x_data(irep,my_thread_id)%maxval_container(my_bin_id)
          integral_caches => x_data(irep,my_thread_id)%integral_caches(:,my_bin_id)
          integral_containers => x_data(irep,my_thread_id)%integral_containers(:,my_bin_id)
          distribution_energy => x_data(irep,my_thread_id)%distribution_energy(my_bin_id)
          do_it = .TRUE.
          bins_left = .TRUE.
          IF( geometry_did_change) THEN
            distribution_energy%ram_counter = HUGE(distribution_energy%ram_counter)
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
      my_istart = distribution_energy%istart
      my_current_counter = 0
      IF(distribution_energy%number_of_atom_quartets == 0 .OR. &
         my_istart == -1_int_8) my_istart = nblocks **4
      atomic_blocks:  DO atom_block = my_istart,nblocks**4-1,n_processes
        latom_block = INT(MODULO(atom_block,nblocks))+1
        tmp_block = atom_block/nblocks
        katom_block = INT(MODULO(tmp_block,nblocks))+1
        IF(latom_block<katom_block) CYCLE
        tmp_block = tmp_block/nblocks
        jatom_block = INT(MODULO(tmp_block,nblocks))+1
        tmp_block = tmp_block/nblocks
        iatom_block = INT(MODULO(tmp_block,nblocks))+1
        IF(jatom_block<iatom_block) CYCLE
        my_current_counter = my_current_counter + 1
        IF(my_current_counter > distribution_energy%number_of_atom_quartets) EXIT atomic_blocks

        iatom_start = actual_x_data%blocks(iatom_block)%istart
        iatom_end =   actual_x_data%blocks(iatom_block)%iend
        jatom_start = actual_x_data%blocks(jatom_block)%istart
        jatom_end =   actual_x_data%blocks(jatom_block)%iend
        katom_start = actual_x_data%blocks(katom_block)%istart
        katom_end =   actual_x_data%blocks(katom_block)%iend
        latom_start = actual_x_data%blocks(latom_block)%istart
        latom_end =   actual_x_data%blocks(latom_block)%iend

        pmax_blocks = MAX( shm_pmax_block(katom_block, iatom_block), &
                           shm_pmax_block(latom_block, jatom_block), &
                           shm_pmax_block(latom_block, iatom_block), &
                           shm_pmax_block(katom_block, jatom_block) )

        IF( 2.0_dp*coeffs_kind_max0 + pmax_blocks  < log10_eps_schwarz ) CYCLE


        CALL build_pair_list(natom, list_ij,set_list_ij,iatom_start, iatom_end, &
                             jatom_start, jatom_end, &
                             kind_of,basis_parameter,particle_set, &
                             do_periodic,screen_coeffs_set,screen_coeffs_kind, &
                             coeffs_kind_max0,log10_eps_schwarz,cell, pmax_blocks, &
                             shm_atomic_pair_list)


        CALL build_pair_list(natom, list_kl,set_list_kl, katom_start, katom_end, &
                             latom_start, latom_end, &
                             kind_of,basis_parameter,particle_set, &
                             do_periodic,screen_coeffs_set,screen_coeffs_kind, &
                             coeffs_kind_max0,log10_eps_schwarz, cell, pmax_blocks, &
                             shm_atomic_pair_list)

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

          DO i_list_kl=1,list_kl%n_element
            katom=list_kl%elements(i_list_kl)%pair(1)
            latom=list_kl%elements(i_list_kl)%pair(2)

            IF (.NOT. (katom+latom<=iatom+jatom)) CYCLE
            IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
            i_set_list_kl_start=list_kl%elements(i_list_kl)%set_bounds(1)
            i_set_list_kl_stop=list_kl%elements(i_list_kl)%set_bounds(2)
            kkind=list_kl%elements(i_list_kl)%kind_pair(1)
            lkind=list_kl%elements(i_list_kl)%kind_pair(2)
            rc=list_kl%elements(i_list_kl)%r1
            rd=list_kl%elements(i_list_kl)%r2
            rcd2=list_kl%elements(i_list_kl)%dist2

            IF( do_p_screening) THEN
              pmax_atom = MAX(shm_pmax_atom(katom,iatom), &
                              shm_pmax_atom(latom,jatom), &
                              shm_pmax_atom(latom,iatom), &
                              shm_pmax_atom(katom,jatom) )
            ELSE
              pmax_atom = 0.0_dp
            END IF

            screen_kind_ij = screen_coeffs_kind(jkind,ikind)%x(1)*rab2+&
                             screen_coeffs_kind(jkind,ikind)%x(2)
            screen_kind_kl = screen_coeffs_kind(lkind,kkind)%x(1)*rcd2+&
                             screen_coeffs_kind(lkind,kkind)%x(2)

            IF( screen_kind_ij + screen_kind_kl + pmax_atom < log10_eps_schwarz ) CYCLE

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

            IF( do_p_screening) THEN
              swap_id = 0
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

            IF(geometry_did_change ) THEN
              IF(.NOT.memory_parameter%do_all_on_the_fly) THEN
                ! ** We know the maximum amount of integrals that we can store per MPI-process
                ! ** Now we need to sum the current memory usage among all openMP threads
                ! ** We can just read what is currently stored on the corresponding x_data type
                ! ** If this is thread i and it tries to read the data from thread j, that is
                ! ** currently writing that data, we just dont care, because the possible error
                ! ** is of the order of CACHE_SIZE
                mem_compression_counter = 0
                DO i=1,n_threads
                  mem_compression_counter = mem_compression_counter + &
                                            x_data(irep,i)%memory_parameter%actual_memory_usage *&
                                            memory_parameter%cache_size
                END DO
                IF(mem_compression_counter > memory_parameter%max_compression_counter) THEN
                  buffer_overflow = .TRUE.
                  IF( do_dynamic_load_balancing ) THEN
                    distribution_energy%ram_counter = counter
                  ELSE
                    memory_parameter%ram_counter = counter
                  END IF
                ELSE
                  counter = counter + 1
                  buffer_overflow = .FALSE.
                END IF
              END IF
            ELSE
              IF(.NOT.memory_parameter%do_all_on_the_fly) THEN
                IF( do_dynamic_load_balancing) THEN
                  IF(distribution_energy%ram_counter == counter) THEN
                    buffer_overflow = .TRUE.
                  ELSE
                    counter = counter + 1
                    buffer_overflow = .FALSE.
                  END IF

                ELSE
                  IF(memory_parameter%ram_counter == counter) THEN
                    buffer_overflow = .TRUE.
                  ELSE
                    counter = counter + 1
                    buffer_overflow = .FALSE.
                  END IF
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

            DO i_set_list_ij=i_set_list_ij_start, i_set_list_ij_stop
              iset=set_list_ij(i_set_list_ij)%pair(1)
              jset=set_list_ij(i_set_list_ij)%pair(2)

              ncob = npgfb(jset)*ncoset(lb_max(jset))
              max_val1 = screen_coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + &
                         screen_coeffs_set(jset,iset,jkind,ikind)%x(2)

              IF( max_val1 + screen_kind_kl + pmax_atom < log10_eps_schwarz) CYCLE

              sphi_a_ext_set => sphi_a_ext(:,:,:,iset)
              sphi_b_ext_set => sphi_b_ext(:,:,:,jset)
              DO i_set_list_kl=i_set_list_kl_start, i_set_list_kl_stop
                kset=set_list_kl(i_set_list_kl)%pair(1)
                lset=set_list_kl(i_set_list_kl)%pair(2)

                max_val2_set = (screen_coeffs_set(lset,kset,lkind,kkind)%x(1)*rcd2 + &
                                screen_coeffs_set(lset,kset,lkind,kkind)%x(2) )
                max_val2 = max_val1 + max_val2_set

                !! Near field screening
                IF(max_val2 + pmax_atom <log10_eps_schwarz) CYCLE
                sphi_c_ext_set => sphi_c_ext(:,:,:,kset)
                sphi_d_ext_set => sphi_d_ext(:,:,:,lset)
                !! get max_vals if we screen on initial density
                IF( do_p_screening) THEN
                  CALL get_pmax_val(ptr_p_1, ptr_p_2, ptr_p_3, ptr_p_4, &
                                    iset, jset, kset, lset, &
                                    pmax_entry, swap_id)
                ELSE
                  pmax_entry = 0.0_dp
                END IF
                log10_pmax = pmax_entry
                max_val2 = max_val2 + log10_pmax
                IF(max_val2<log10_eps_schwarz) CYCLE
                pmax_entry = EXP(log10_pmax*ln_10)

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
                  nbits = EXPONENT(ANINT(spherical_estimate*pmax_entry/eps_storage)) + 1
                  buffer_left = nints
                  buffer_start = 1
                  DO WHILE (buffer_left > 0 )
                    buffer_size = MIN(buffer_left, cache_size)
                    IF(.NOT. use_disk_storage) THEN
                      CALL hfx_get_mult_cache_elements(primitive_integrals(buffer_start), &
                                                       buffer_size, nbits, &
                                                       integral_caches(nbits), &
                                                       integral_containers(nbits), &
                                                       eps_storage, pmax_entry, &
                                                       memory_parameter%actual_memory_usage, &
                                                       use_disk_storage)
                    ELSE
                      CALL hfx_get_mult_cache_elements(primitive_integrals(buffer_start), &
                                                       buffer_size, nbits, &
                                                       integral_caches_disk(nbits), &
                                                       integral_containers_disk(nbits), &
                                                       eps_storage, pmax_entry, &
                                                       memory_parameter%actual_memory_usage_disk, &
                                                       use_disk_storage)
                    END IF
                    buffer_left = buffer_left - buffer_size
                    buffer_start = buffer_start + buffer_size
                  END DO
                END IF
                !! Calculate integrals if we run out of buffer or the geometry did change
                IF(geometry_did_change .OR. buffer_overflow) THEN
                  max_contraction_val =  max_contraction(iset,iatom) * &
                                         max_contraction(jset,jatom) * &
                                         max_contraction(kset,katom) * &
                                         max_contraction(lset,latom) * pmax_entry
                  tmp_R_1 => radii_pgf(:,:,jset,iset,jkind,ikind)
                  tmp_R_2 => radii_pgf(:,:,lset,kset,lkind,kkind)
                  tmp_screen_pgf1 => screen_coeffs_pgf(:,:,jset,iset,jkind,ikind)
                  tmp_screen_pgf2 => screen_coeffs_pgf(:,:,lset,kset,lkind,kkind)

                  CALL coulomb4(private_lib, ra, rb, rc, rd, npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset), &
                                la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset),&
                                sphi_a_u1,sphi_a_u2,sphi_a_u3,&
                                sphi_b_u1,sphi_b_u2,sphi_b_u3,&
                                sphi_c_u1,sphi_c_u2,sphi_c_u3,&
                                sphi_d_u1,sphi_d_u2,sphi_d_u3,&
                                zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset),&
                                zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                                primitive_integrals,&
                                potential_parameter, &
                                actual_x_data%neighbor_cells, screen_coeffs_set(jset,iset,jkind,ikind)%x,&
                                screen_coeffs_set(lset,kset,lkind,kkind)%x, eps_schwarz, &
                                max_contraction_val, cartesian_estimate,  cell, neris_tmp,&
                                log10_pmax, log10_eps_schwarz, &
                                tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2,&
                                pgf_list_ij, pgf_list_kl, pgf_product_list, &
                                nsgfl_a(:,iset), nsgfl_b(:,jset),&
                                nsgfl_c(:,kset), nsgfl_d(:,lset),&
                                sphi_a_ext_set,&
                                sphi_b_ext_set,&
                                sphi_c_ext_set,&
                                sphi_d_ext_set,&
                                ee_work,ee_work2,ee_buffer1,ee_buffer2, ee_primitives_tmp,&
                                nimages,do_periodic, p_work)

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
                                        use_disk_storage, max_val_memory, error)
                      ELSE
                        CALL hfx_add_single_cache_element(estimate_to_store_int, 6,&
                                        maxval_cache_disk, maxval_container_disk, &
                                        memory_parameter%actual_memory_usage_disk, &
                                        use_disk_storage, error=error)
                      END IF
                    END IF
                  END IF

                  IF(cartesian_estimate<eps_schwarz) CYCLE

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
                                                         use_disk_storage, max_val_memory, error)
                    ELSE
                      CALL hfx_add_single_cache_element(estimate_to_store_int, 6, &
                                                         maxval_cache_disk, maxval_container_disk, &
                                                         memory_parameter%actual_memory_usage_disk, &
                                                         use_disk_storage, error=error)
                    END IF
                  END IF
                  spherical_estimate = SET_EXPONENT(1.0_dp,estimate_to_store_int+1)
                  IF(spherical_estimate*pmax_entry<eps_schwarz) CYCLE
                  IF(.NOT.buffer_overflow) THEN
                    nbits = EXPONENT(ANINT(spherical_estimate*pmax_entry/eps_storage)) + 1
                    buffer_left = nints
                    buffer_start = 1
                    IF( .NOT. use_disk_storage) THEN
                      neris_incore = neris_incore + nints
                    ELSE
                      neris_disk = neris_disk + nints
                    END IF
                    DO WHILE (buffer_left > 0 )
                      buffer_size = MIN(buffer_left, CACHE_SIZE)
                      IF(.NOT. use_disk_storage) THEN
                        CALL hfx_add_mult_cache_elements(primitive_integrals(buffer_start),&
                                                          buffer_size, nbits, &
                                                          integral_caches(nbits), &
                                                          integral_containers(nbits), &
                                                          eps_storage, pmax_entry, &
                                                          memory_parameter%actual_memory_usage, &
                                                          use_disk_storage, error)
                      ELSE
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
                END IF
!!!  DEBUG, print out primitive integrals and indices. Only works serial no OMP  !!!
!               CALL print_integrals(iatom,jatom,katom,latom,shm_set_offset,shm_atomic_block_offset,&
!                                    iset, jset, kset, lset,nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset),primitive_integrals)
!!!  DEBUG END
                IF(.NOT.is_anti_symmetric)THEN
                  !! Update Kohn-Sham matrix
                  CALL update_fock_matrix(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                          fac, symm_fac, full_density, full_ks, &
                                          primitive_integrals, pbd_buf, pbc_buf, pad_buf, pac_buf, kbd_buf,&
                                          kbc_buf, kad_buf, kac_buf, iatom, jatom, katom, latom,&
                                          iset, jset, kset, lset, offset_bd_set, offset_bc_set, offset_ad_set, offset_ac_set,&
                                          atomic_offset_bd, atomic_offset_bc, atomic_offset_ad, atomic_offset_ac, is_anti_symmetric)
                  IF( .NOT. treat_lsd_in_core ) THEN
                    IF( nspins == 2 ) THEN
                      CALL update_fock_matrix(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                              fac, symm_fac, full_density_beta, full_ks_beta, &
                                              primitive_integrals, pbd_buf, pbc_buf, pad_buf, pac_buf, kbd_buf,&
                                              kbc_buf, kad_buf, kac_buf, iatom, jatom, katom, latom,&
                                              iset, jset, kset, lset, offset_bd_set, offset_bc_set, offset_ad_set,offset_ac_set,&
                                              atomic_offset_bd, atomic_offset_bc, atomic_offset_ad, atomic_offset_ac,&
                                              is_anti_symmetric)
                    END IF
                  END IF
                ELSE
                  !! Update Kohn-Sham matrix
                  CALL update_fock_matrix_as(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                          fac, symm_fac, full_density, full_ks, &
                                          primitive_integrals, pbd_buf, pbc_buf, pad_buf, pac_buf, kbd_buf,&
                                          kbc_buf, kad_buf, kac_buf, iatom, jatom, katom, latom,&
                                          iset, jset, kset, lset, offset_bd_set, offset_bc_set, offset_ad_set, offset_ac_set,&
                                          atomic_offset_bd, atomic_offset_bc, atomic_offset_ad, atomic_offset_ac, is_anti_symmetric)
                  IF( .NOT. treat_lsd_in_core ) THEN
                    IF( nspins == 2 ) THEN
                      CALL update_fock_matrix_as(nsgfa(iset), nsgfb(jset), nsgfc(kset), nsgfd(lset), &
                                              fac, symm_fac, full_density_beta, full_ks_beta, &
                                              primitive_integrals, pbd_buf, pbc_buf, pad_buf, pac_buf, kbd_buf,&
                                              kbc_buf, kad_buf, kac_buf, iatom, jatom, katom, latom,&
                                              iset, jset, kset, lset, offset_bd_set, offset_bc_set, offset_ad_set,offset_ac_set,&
                                              atomic_offset_bd, atomic_offset_bc, atomic_offset_ad, atomic_offset_ac,&
                                              is_anti_symmetric)
                    END IF
                  END IF           
                END IF
              END DO ! i_set_list_kl
            END DO ! i_set_list_ij
            IF( do_disk_storage ) THEN
              buffer_overflow = .TRUE.
            END IF
          END DO !i_list_ij
        END DO !i_list_kl
      END DO atomic_blocks
      bintime_stop = m_walltime()
      IF (geometry_did_change) THEN
        distribution_energy%time_first_scf=bintime_stop-bintime_start
      ELSE
        distribution_energy%time_other_scf= &
            distribution_energy%time_other_scf+bintime_stop-bintime_start
      ENDIF
!$OMP MASTER
      CALL timestop(handle_bin)
!$OMP END MASTER
    END DO !bin


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
                                     hfx_do_eval_energy, error)

!$OMP MASTER
      CALL cp_print_key_finished_output(iw,logger,hfx_section,&
                                       "LOAD_BALANCE%PRINT/LOAD_BALANCE_INFO", error=error)
!$OMP END MASTER
    END IF

!$OMP BARRIER
!$OMP MASTER
    memsize_after=m_memory()
!$OMP END MASTER
!$OMP BARRIER

    DEALLOCATE(primitive_integrals,STAT=stat)
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
    stor_count_max_val = max_val_memory * memory_parameter%cache_size
!$OMP ATOMIC
    shm_storage_counter_integrals = shm_storage_counter_integrals + storage_counter_integrals
!$OMP ATOMIC
    shm_stor_count_int_disk =  shm_stor_count_int_disk + stor_count_int_disk
!$OMP ATOMIC
    shm_neris_incore = shm_neris_incore + neris_incore
!$OMP ATOMIC
    shm_neris_disk = shm_neris_disk + neris_disk
!$OMP ATOMIC
    shm_stor_count_max_val = shm_stor_count_max_val + stor_count_max_val
!$OMP BARRIER


    ! ** Calculate how much memory has already been used (might be needed for in-core forces
!$OMP MASTER
    shm_mem_compression_counter = 0
    DO i=1,n_threads
      shm_mem_compression_counter = shm_mem_compression_counter + &
                                    x_data(irep,i)%memory_parameter%actual_memory_usage *&
                                    memory_parameter%cache_size
    END DO
!$OMP END MASTER
!$OMP BARRIER
    actual_x_data%memory_parameter%final_comp_counter_energy = shm_mem_compression_counter


!$OMP MASTER
    !! Calculate the exchange energies from the Kohn-Sham matrix. Before we can go on, we have to symmetrize.
    ene_x_aa = 0.0_dp
    ene_x_bb = 0.0_dp

    mb_size_p = shm_block_offset(ncpu+1) /1024/128
    mb_size_f = shm_block_offset(ncpu+1) /1024/128
    IF (.NOT. treat_lsd_in_core ) THEN
      IF( nspins == 2) THEN
        mb_size_f = mb_size_f * 2
        mb_size_p = mb_size_p * 2
      END IF
    END IF
    !! size of primitive_integrals(not shared)
    mb_size_buffers = INT(nsgf_max,int_8)**4 * n_threads
    !! fock density buffers (not shared)
    mb_size_buffers = mb_size_buffers + INT(nsgf_max,int_8)**2 * n_threads
    subtr_size_mb = subtr_size_mb + 8_int_8*nsgf_max**2*n_threads
    !! size of screening functions (shared)
    mb_size_buffers = mb_size_buffers + max_pgf**2*max_set**2*nkind**2 &
                                      + max_set**2*nkind**2 &
                                      + nkind**2 &
                                      + max_pgf**2*max_set**2*nkind**2
    !! is_assoc (shared)
    mb_size_buffers = mb_size_buffers+natom**2
    ! ** pmax_atom (shared)
    IF( do_p_screening) THEN
      mb_size_buffers = mb_size_buffers + natom**2
    END IF
    IF(screening_parameter%do_p_screening_forces) THEN
      IF( memory_parameter%treat_forces_in_core) THEN
        mb_size_buffers = mb_size_buffers + natom**2
      END IF
    END IF
    ! ** Initial P only MAX(alpha,beta) (shared)
    IF(do_p_screening .OR. screening_parameter%do_p_screening_forces) THEN
      mb_size_buffers = mb_size_buffers + memory_parameter%size_p_screen
    END IF
    ! ** In core forces require their own initial P
    IF(screening_parameter%do_p_screening_forces) THEN
      IF( memory_parameter%treat_forces_in_core) THEN
        mb_size_buffers = mb_size_buffers + memory_parameter%size_p_screen
      END IF
    END IF


    !! mb
    mb_size_buffers = mb_size_buffers /1024/128

    afac=1.0_dp
    IF(is_anti_symmetric) afac=-1.0_dp
    CALL timestop(handle_main)
    ene_x_aa_diag = 0.0_dp
    ene_x_bb_diag = 0.0_dp
    DO iatom = 1,natom
      ikind = kind_of(iatom)
      nseta = basis_parameter(ikind)%nset
      nsgfa => basis_parameter(ikind)%nsgf
      jatom = iatom
      jkind = kind_of(jatom)
      nsetb = basis_parameter(jkind)%nset
      nsgfb => basis_parameter(jkind)%nsgf
      act_atomic_block_offset = shm_atomic_block_offset(jatom,iatom)
      DO iset=1,nseta
        DO jset=1,nsetb
          act_set_offset = shm_set_offset(jset,iset,jkind,ikind)
          i = act_set_offset + act_atomic_block_offset - 1
          DO ma=1,nsgfa(iset)
            j = shm_set_offset(iset,jset,jkind,ikind) + act_atomic_block_offset - 1 + ma -1
            DO mb=1,nsgfb(jset)
              IF(i>j) THEN
                full_ks(i) = (full_ks(i) + full_ks(j)*afac)
                full_ks(j) = full_ks(i)*afac 
                IF( .NOT. treat_lsd_in_core .AND. nspins == 2 ) THEN
                  full_ks_beta(i) = (full_ks_beta(i) + full_ks_beta(j)*afac)
                  full_ks_beta(j) = full_ks_beta(i)*afac
                END IF
              END IF
              ! ** For adiabatically rescaled functionals we need the energy coming from the diagonal elements
              IF( i==j ) THEN
                ene_x_aa_diag = ene_x_aa_diag + full_ks(i)*full_density(i)
                IF( .NOT. treat_lsd_in_core .AND. nspins == 2 ) THEN
                  ene_x_bb_diag = ene_x_bb_diag + full_ks_beta(i)*full_density_beta(i)
                END IF
              END IF
              i = i+ 1
              j = j+nsgfa(iset)
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL mp_sync(para_env%group)
    afac=1.0_dp
    IF( is_anti_symmetric) afac=0._dp
    IF( distribute_fock_matrix) THEN
      !! Distribute the current KS-matrix to all the processes
      CALL timeset(routineN//"_dist_KS",handle_dist_ks)
      CALL distribute_ks_matrix(para_env, full_ks, ks_matrix(ispin)%matrix, shm_number_of_p_entries, &
                                shm_block_offset, natom, kind_of, basis_parameter, &
                                off_diag_fac=0.5_dp,diag_fac=afac,error=error)

      NULLIFY(full_ks)
      DEALLOCATE(shm_master_x_data%full_ks_alpha, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      IF( .NOT. treat_lsd_in_core ) THEN
        IF(nspins == 2) THEN
          CALL distribute_ks_matrix(para_env, full_ks_beta, ks_matrix(2)%matrix, shm_number_of_p_entries, &
                                    shm_block_offset, natom, kind_of, basis_parameter, &
                                    off_diag_fac=0.5_dp, diag_fac=afac, error=error)
          NULLIFY(full_ks_beta)
          DEALLOCATE(shm_master_x_data%full_ks_beta, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
      END IF
      CALL timestop(handle_dist_ks)
    END IF


    IF( distribute_fock_matrix) THEN
      !! ** Calculate the exchange energy
      CALL cp_dbcsr_trace(ks_matrix(ispin)%matrix, rho_ao(ispin)%matrix, &
                          ene_x_aa, error=error)
      !for ADMMS, we need the exchange matrix k(d) for both spins 
      IF( dft_control%do_admm ) THEN
        CALL cp_dbcsr_copy(matrix_ks_aux_fit_hfx(ispin)%matrix, ks_matrix(ispin)%matrix, &
                           name="HF exch. part of matrix_ks_aux_fit for ADMMS", error=error) 
      END IF

      IF( nspins == 2 .AND. .NOT. treat_lsd_in_core ) THEN
        CALL cp_dbcsr_trace(ks_matrix(2)%matrix, rho_ao(2)%matrix, &
                            ene_x_bb, error=error)
        !for ADMMS, we need the exchange matrix k(d) for both spins
        IF( dft_control%do_admm ) THEN
          CALL cp_dbcsr_copy(matrix_ks_aux_fit_hfx(2)%matrix, ks_matrix(2)%matrix, &
                             name="HF exch. part of matrix_ks_aux_fit for ADMMS", error=error)          
        END IF
      END IF

      !! Update energy type
      energy%ex = 0.5_dp*(ene_x_aa+ene_x_bb) + energy%ex
    ELSE
      ! ** It is easier to correct the following expression by the diagonal energy contribution,
      ! ** than explicitly going throuhg the diagonal elements
      DO pa=1,SIZE(full_ks)
        ene_x_aa = ene_x_aa + full_ks(pa)*full_density(pa)
      END DO
      ! ** Now correct
      ene_x_aa = (ene_x_aa + ene_x_aa_diag)*0.5_dp
      IF( nspins == 2 ) THEN
        DO pa=1,SIZE(full_ks)
          ene_x_bb = ene_x_bb + full_ks_beta(pa)*full_density_beta(pa)
        END DO
        ! ** Now correct
        ene_x_bb = (ene_x_bb + ene_x_bb_diag)*0.5_dp
      END IF
      CALL mp_sum(ene_x_aa, para_env%group)
      IF(nspins==2) CALL mp_sum(ene_x_bb,para_env%group)
      energy%ex = 0.5_dp*(ene_x_aa+ene_x_bb)
    END IF

    !! Print some memeory information if this is the first step
    IF( geometry_did_change ) THEN
      tmp_i8(1:8)=(/shm_storage_counter_integrals,shm_neris_onthefly,shm_neris_incore,shm_neris_disk, &
                    shm_neris_total,shm_stor_count_int_disk,shm_nprim_ints,shm_stor_count_max_val/)
      CALL mp_sum(tmp_i8,para_env%group)
      shm_storage_counter_integrals=tmp_i8(1)
      shm_neris_onthefly=tmp_i8(2)
      shm_neris_incore=tmp_i8(3)
      shm_neris_disk=tmp_i8(4)
      shm_neris_total=tmp_i8(5)
      shm_stor_count_int_disk=tmp_i8(6)
      shm_nprim_ints=tmp_i8(7)
      shm_stor_count_max_val=tmp_i8(8)
      CALL mp_max(memsize_after, para_env%group)
      mem_eris = (shm_storage_counter_integrals+128*1024-1)/1024/128
      compression_factor = REAL(shm_neris_incore,dp)/REAL(shm_storage_counter_integrals,dp)
      mem_eris_disk = (shm_stor_count_int_disk+128*1024-1)/1024/128
      compression_factor_disk = REAL(shm_neris_disk,dp)/REAL(shm_stor_count_int_disk,dp)
      mem_max_val = (shm_stor_count_max_val+128*1024-1)/1024/128

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
        WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Number of cart. primitive ERI's calculated:      ",shm_nprim_ints

        WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Number of sph. ERI's calculated:           ",shm_neris_total

        WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Number of sph. ERI's stored in-core:        ",shm_neris_incore

        WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Number of sph. ERI's stored on disk:        ",shm_neris_disk

        WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Number of sph. ERI's calculated on the fly: ",shm_neris_onthefly

        WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Total memory consumption ERI's RAM [MB's]:            ",mem_eris

        WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Whereof max-vals [MB's]:            ",mem_max_val

        WRITE (UNIT=iw,FMT="((T3,A,T60,F21.2))")&
                 "HFX_MEM_INFO| Total compression factor ERI's RAM:                  ",compression_factor

        WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Total memory consumption ERI's disk [MB's]:       ",mem_eris_disk

        WRITE (UNIT=iw,FMT="((T3,A,T60,F21.2))")&
                 "HFX_MEM_INFO| Total compression factor ERI's disk:             ",compression_factor_disk

        WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Size of density/Fock matrix [MB's]:             ", 2_int_8*mb_size_p

        IF( do_periodic) THEN
          WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Size of buffers [MB's]:             ", mb_size_buffers
          WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Number of periodic image cells considered: ", SIZE(shm_master_x_data%neighbor_cells)
        ELSE
           WRITE (UNIT=iw,FMT="((T3,A,T60,I21))")&
                 "HFX_MEM_INFO| Size of buffers [MB's]:             ", mb_size_buffers
        END IF
        WRITE (UNIT=iw,FMT="((T3,A,T60,I21),/)")&
                 "HFX_MEM_INFO| Est. max. program size after HFX  [MB's]:", memsize_after/(1024*1024)
        CALL m_flush(iw)
      END IF

      CALL cp_print_key_finished_output(iw,logger,hfx_section,&
                                      "HF_INFO", error=error)
    END IF
!$OMP END MASTER

    !! flush caches if the geometry changed
    IF( do_dynamic_load_balancing ) THEN
      my_bin_size = SIZE(actual_x_data%distribution_energy)
    ELSE
      my_bin_size = 1
    END IF

    IF(geometry_did_change ) THEN
      IF(.NOT.memory_parameter%do_all_on_the_fly) THEN
        DO bin = 1,my_bin_size
          maxval_cache => actual_x_data%maxval_cache(bin)
          maxval_container => actual_x_data%maxval_container(bin)
          integral_caches => actual_x_data%integral_caches(:,bin)
          integral_containers => actual_x_data%integral_containers(:,bin)
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
    IF(.NOT.memory_parameter%do_all_on_the_fly) THEN
      DO bin = 1,my_bin_size
        maxval_cache => actual_x_data%maxval_cache(bin)
        maxval_container => actual_x_data%maxval_container(bin)
        integral_caches => actual_x_data%integral_caches(:,bin)
        integral_containers => actual_x_data%integral_containers(:,bin)

        CALL hfx_reset_cache_and_container(maxval_cache, maxval_container, memory_parameter%actual_memory_usage, .FALSE.)
        DO i=1,64
          CALL hfx_reset_cache_and_container(integral_caches(i), integral_containers(i), &
                                             memory_parameter%actual_memory_usage, &
                                             .FALSE.)
        END DO
      END DO
    END IF

    !! Since the I/O routines are no thread-safe, i.e. the procedure to get the unit number, put a lock here
!$OMP CRITICAL
    IF(do_disk_storage ) THEN
      !! flush caches if the geometry changed
      IF(geometry_did_change ) THEN
        CALL hfx_flush_last_cache(bits_max_val, maxval_cache_disk, maxval_container_disk, &
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
    DEALLOCATE(full_density,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF( .NOT. treat_lsd_in_core ) THEN
      IF(nspins==2) THEN
        DEALLOCATE(full_density_beta,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    IF (do_dynamic_load_balancing) THEN
      DEALLOCATE(shm_master_x_data%task_list, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
!$OMP END MASTER
    DEALLOCATE(pbd_buf, pbc_buf, pad_buf, pac_buf, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(kbd_buf,kbc_buf,kad_buf,kac_buf, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(set_list_ij, set_list_kl, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

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

    DEALLOCATE(max_contraction, kind_of, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(ee_work, ee_work2, ee_buffer1, ee_buffer2, ee_primitives_tmp, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(nimages, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

!$OMP BARRIER
!$OMP END PARALLEL

#if defined(__USE_PAT)
  CALL PAT_sampling_state(PAT_STATE_OFF,istat)
#endif

    CALL timestop(handle)
  END SUBROUTINE integrate_four_center

! *****************************************************************************
!> \brief calculates two-electron integrals of a quartet/shell using the library
!>      lib_int in the periodic case
!> \param lib ...
!> \param ra ...
!> \param rb ...
!> \param rc ...
!> \param rd ...
!> \param npgfa ...
!> \param npgfb ...
!> \param npgfc ...
!> \param npgfd ...
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param lc_min ...
!> \param lc_max ...
!> \param ld_min ...
!> \param ld_max ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param nsgfd ...
!> \param sphi_a_u1 ...
!> \param sphi_a_u2 ...
!> \param sphi_a_u3 ...
!> \param sphi_b_u1 ...
!> \param sphi_b_u2 ...
!> \param sphi_b_u3 ...
!> \param sphi_c_u1 ...
!> \param sphi_c_u2 ...
!> \param sphi_c_u3 ...
!> \param sphi_d_u1 ...
!> \param sphi_d_u2 ...
!> \param sphi_d_u3 ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param zetd ...
!> \param primitive_integrals array of primitive_integrals
!> \param potential_parameter contains info for libint
!> \param neighbor_cells Periodic images
!> \param screen1 set based coefficients for near field screening
!> \param screen2 set based coefficients for near field screening
!> \param eps_schwarz threshold
!> \param max_contraction_val maximum multiplication factor for cart -> sph
!> \param cart_estimate maximum calculated integral value
!> \param cell cell
!> \param neris_tmp counter for calculated cart integrals
!> \param log10_pmax logarithm of initial p matrix max element
!> \param log10_eps_schwarz log of threshold
!> \param R1_pgf coefficients for radii of product distribution function
!> \param R2_pgf coefficients for radii of product distribution function
!> \param pgf1 schwarz coefficients pgf basid
!> \param pgf2 schwarz coefficients pgf basid
!> \param pgf_list_ij ...
!> \param pgf_list_kl ...
!> \param pgf_product_list ...
!> \param nsgfl_a ...
!> \param nsgfl_b ...
!> \param nsgfl_c ...
!> \param nsgfl_d ...
!> \param sphi_a_ext ...
!> \param sphi_b_ext ...
!> \param sphi_c_ext ...
!> \param sphi_d_ext ...
!> \param ee_work ...
!> \param ee_work2 ...
!> \param ee_buffer1 ...
!> \param ee_buffer2 ...
!> \param ee_primitives_tmp ...
!> \param nimages ...
!> \param do_periodic ...
!> \param p_work ...
!> \par History
!>      11.2006 created [Manuel Guidon]
!>      02.2009 completely rewritten screening part [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE coulomb4(lib, ra, rb, rc, rd, npgfa, npgfb, npgfc, npgfd, &
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
                      ee_work,ee_work2,ee_buffer1,ee_buffer2, ee_primitives_tmp,&
                      nimages, do_periodic, p_work)

    TYPE(lib_int)                            :: lib
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
    REAL(dp), &
      DIMENSION(nsgfa, nsgfb, nsgfc, nsgfd)  :: primitive_integrals
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
    REAL(dp), DIMENSION(*)                   :: ee_work, ee_work2, &
                                                ee_buffer1, ee_buffer2, &
                                                ee_primitives_tmp
    INTEGER, DIMENSION(*)                    :: nimages
    LOGICAL, INTENT(IN)                      :: do_periodic
    REAL(dp), DIMENSION(:), POINTER          :: p_work

    INTEGER :: ipgf, jpgf, kpgf, la, lb, lc, ld, list_ij, list_kl, lpgf, &
      max_l, ncoa, ncob, ncoc, ncod, nelements_ij, nelements_kl, nneighbors3, &
      nproducts, nsgfla, nsgflb, nsgflc, nsgfld, nsoa, nsob, nsoc, nsod, &
      s_offset_a, s_offset_b, s_offset_c, s_offset_d
    REAL(dp)                                 :: EtaInv, tmp_max, ZetaInv

    CALL build_pair_list_pgf(npgfa, npgfb, pgf_list_ij, zeta, zetb, screen1, screen2,&
                             pgf1, R1_pgf, log10_pmax, log10_eps_schwarz, ra, rb,&
                             nelements_ij, &
                             neighbor_cells, nimages, do_periodic)
    CALL build_pair_list_pgf(npgfc, npgfd, pgf_list_kl, zetc, zetd, screen2, screen1, &
                             pgf2, R2_pgf, log10_pmax, log10_eps_schwarz, rc, rd,&
                             nelements_kl, &
                             neighbor_cells, nimages, do_periodic)

    cart_estimate = 0.0_dp
    neris_tmp = 0
    primitive_integrals = 0.0_dp
    max_l = la_max + lb_max + lc_max + ld_max
    DO list_ij = 1,nelements_ij
      ZetaInv = pgf_list_ij(list_ij)%ZetaInv
      ipgf = pgf_list_ij(list_ij)%ipgf
      jpgf = pgf_list_ij(list_ij)%jpgf

      DO list_kl = 1,nelements_kl
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
                CALL evaluate_eri(lib, nproducts, pgf_product_list, &
                                  la, lb, lc ,ld,&
                                  ncoa,ncob,ncoc,ncod,&
                                  nsgfa, nsgfb, nsgfc, nsgfd,&
                                  primitive_integrals, &
                                  max_contraction_val, tmp_max, eps_schwarz,&
                                  neris_tmp,ZetaInv,EtaInv,&
                                  s_offset_a, s_offset_b, s_offset_c, s_offset_d,&
                                  nsgfla, nsgflb, nsgflc, nsgfld, nsoa, nsob, nsoc, nsod, &
                                  sphi_a_ext(1,la+1,ipgf), &
                                  sphi_b_ext(1,lb+1,jpgf), &
                                  sphi_c_ext(1,lc+1,kpgf), &
                                  sphi_d_ext(1,ld+1,lpgf), &
                                  ee_work,ee_work2,ee_buffer1,ee_buffer2, ee_primitives_tmp,&
                                  p_work)
                cart_estimate = MAX(tmp_max,cart_estimate)
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

  END SUBROUTINE coulomb4

! *****************************************************************************
!> \brief Given a 2d index pair, this function returns a 1d index pair for
!>        a symmetric upper triangle NxN matrix
!>        The compiler should inline this function, therefore it appears in
!>        several modules
!> \param i 2d index
!> \param j 2d index
!> \param N matrix size
!> \retval get_1D_idx ...
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
!> \brief This routine prefetches density/fock matrix elements and stores them
!>        in cache friendly arrays. These buffers are then used to update the
!>        fock matrix
!> \param ma_max Size of matrix blocks
!> \param mb_max Size of matrix blocks
!> \param mc_max Size of matrix blocks
!> \param md_max Size of matrix blocks
!> \param fac multiplication factor (spin)
!> \param symm_fac multiplication factor (symmetry)
!> \param density upper triangular density matrix
!> \param ks upper triangular fock matrix
!> \param prim primitive integrals
!> \param pbd buffer that will contain P(b,d)
!> \param pbc buffer that will contain P(b,c)
!> \param pad buffer that will contain P(a,d)
!> \param pac buffer that will contain P(a,c)
!> \param kbd buffer for KS(b,d)
!> \param kbc buffer for KS(b,c)
!> \param kad buffer for KS(a,d)
!> \param kac buffer for KS(a,c)
!> \param iatom ...
!> \param jatom ...
!> \param katom ...
!> \param latom ...
!> \param iset ...
!> \param jset ...
!> \param kset ...
!> \param lset ...
!> \param offset_bd_set ...
!> \param offset_bc_set ...
!> \param offset_ad_set ...
!> \param offset_ac_set ...
!> \param atomic_offset_bd ...
!> \param atomic_offset_bc ...
!> \param atomic_offset_ad ...
!> \param atomic_offset_ac ...
!> \param anti_symmetric ...
!> \par History
!>      03.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************

  SUBROUTINE update_fock_matrix(ma_max, mb_max, mc_max, md_max, &
                                fac, symm_fac, density, ks, prim, &
                                pbd, pbc, pad, pac, kbd, kbc, kad, kac, &
                                iatom, jatom, katom, latom, &
                                iset, jset, kset, lset,  offset_bd_set, offset_bc_set, offset_ad_set,&
                                offset_ac_set, atomic_offset_bd, atomic_offset_bc, atomic_offset_ad,&
                                atomic_offset_ac,anti_symmetric)

    INTEGER, INTENT(IN)                      :: ma_max, mb_max, mc_max, md_max
    REAL(dp), INTENT(IN)                     :: fac, symm_fac
    REAL(dp), DIMENSION(:), INTENT(IN)       :: density
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: ks
    REAL(dp), DIMENSION(ma_max*mb_max*mc_max&
      *md_max), INTENT(IN)                   :: prim
    REAL(dp), DIMENSION(*), INTENT(INOUT)    :: pbd, pbc, pad, pac, kbd, kbc, &
                                                kad, kac
    INTEGER, INTENT(IN)                      :: iatom, jatom, katom, latom, &
                                                iset, jset, kset, lset
    INTEGER, DIMENSION(:, :), POINTER        :: offset_bd_set, offset_bc_set, &
                                                offset_ad_set, offset_ac_set
    INTEGER, INTENT(IN)                      :: atomic_offset_bd, &
                                                atomic_offset_bc, &
                                                atomic_offset_ad, &
                                                atomic_offset_ac
    LOGICAL                                  :: anti_symmetric

    INTEGER                                  :: i, j, ma, mb, mc, md, &
                                                offset_ac, offset_ad, &
                                                offset_bc, offset_bd

    IF( jatom >= latom ) THEN
      i = 1
      offset_bd = offset_bd_set(jset,lset) + atomic_offset_bd - 1
      j = offset_bd
      DO md=1,md_max
        DO mb=1,mb_max
          pbd(i) = density(j)
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
          pbc(i) = density(j)
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
          pad(i) = density(j)
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
          pac(i) = density(j)
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

    CALL contract_block(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,fac*symm_fac)
    IF(jatom>=latom) THEN
      i = 1
      j = offset_bd
      DO md=1,md_max
        DO mb=1,mb_max
          !$OMP ATOMIC
          ks(j) = ks(j) + kbd(i)
          i = i + 1
          j = j + 1
        END DO
      END DO
    ELSE
      i = 1
      DO md=1,md_max
        j = offset_bd+md-1
        DO mb=1,mb_max
          !$OMP ATOMIC
          ks(j) = ks(j) + kbd(i)
          i = i + 1
          j = j + md_max
        END DO
      END DO
    END IF
    IF(jatom>=katom) THEN
      i = 1
      j = offset_bc
      DO mc=1,mc_max
        DO mb=1,mb_max
          !$OMP ATOMIC
          ks(j) = ks(j) + kbc(i)
          i = i + 1
          j = j + 1
        END DO
      END DO
    ELSE
      i = 1
      DO mc=1,mc_max
        j = offset_bc + mc -1
        DO mb=1,mb_max
          !$OMP ATOMIC
          ks(j) = ks(j) + kbc(i)
          i = i + 1
          j = j + mc_max
        END DO
      END DO
    END IF
    IF(iatom>=latom) THEN
      i = 1
      j = offset_ad
      DO md=1,md_max
        DO ma=1,ma_max
          !$OMP ATOMIC
          ks(j) = ks(j) + kad(i)
          i = i + 1
          j = j + 1
        END DO
      END DO
    ELSE
      i = 1
      DO md=1,md_max
        j = offset_ad + md - 1
        DO ma=1,ma_max
          !$OMP ATOMIC
          ks(j) = ks(j) + kad(i)
          i = i + 1
          j = j + md_max
        END DO
      END DO
    END IF
    IF(iatom>=katom) THEN
      i = 1
      j = offset_ac
      DO mc=1,mc_max
        DO ma=1,ma_max
          !$OMP ATOMIC
          ks(j) = ks(j) + kac(i)
          i = i + 1
          j = j + 1
        END DO
      END DO
    ELSE
      i = 1
      DO mc=1,mc_max
        j = offset_ac + mc -1
        DO ma=1,ma_max
          !$OMP ATOMIC
          ks(j) = ks(j) + kac(i)
          i = i + 1
          j = j + mc_max
        END DO
      END DO
    END IF
  END SUBROUTINE update_fock_matrix

! *****************************************************************************
!> \brief ...
!> \param ma_max ...
!> \param mb_max ...
!> \param mc_max ...
!> \param md_max ...
!> \param fac ...
!> \param symm_fac ...
!> \param density ...
!> \param ks ...
!> \param prim ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param iatom ...
!> \param jatom ...
!> \param katom ...
!> \param latom ...
!> \param iset ...
!> \param jset ...
!> \param kset ...
!> \param lset ...
!> \param offset_bd_set ...
!> \param offset_bc_set ...
!> \param offset_ad_set ...
!> \param offset_ac_set ...
!> \param atomic_offset_bd ...
!> \param atomic_offset_bc ...
!> \param atomic_offset_ad ...
!> \param atomic_offset_ac ...
!> \param anti_symmetric ...
! *****************************************************************************
  SUBROUTINE update_fock_matrix_as(ma_max, mb_max, mc_max, md_max, &
                                fac, symm_fac, density, ks, prim, &
                                pbd, pbc, pad, pac, kbd, kbc, kad, kac, &
                                iatom, jatom, katom, latom, &
                                iset, jset, kset, lset,  offset_bd_set, offset_bc_set, offset_ad_set,&
                                offset_ac_set, atomic_offset_bd, atomic_offset_bc, atomic_offset_ad,&
                                atomic_offset_ac,anti_symmetric)

    INTEGER, INTENT(IN)                      :: ma_max, mb_max, mc_max, md_max
    REAL(dp), INTENT(IN)                     :: fac, symm_fac
    REAL(dp), DIMENSION(:), INTENT(IN)       :: density
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: ks
    REAL(dp), DIMENSION(ma_max*mb_max*mc_max&
      *md_max), INTENT(IN)                   :: prim
    REAL(dp), DIMENSION(*), INTENT(INOUT)    :: pbd, pbc, pad, pac, kbd, kbc, &
                                                kad, kac
    INTEGER, INTENT(IN)                      :: iatom, jatom, katom, latom, &
                                                iset, jset, kset, lset
    INTEGER, DIMENSION(:, :), POINTER        :: offset_bd_set, offset_bc_set, &
                                                offset_ad_set, offset_ac_set
    INTEGER, INTENT(IN)                      :: atomic_offset_bd, &
                                                atomic_offset_bc, &
                                                atomic_offset_ad, &
                                                atomic_offset_ac
    LOGICAL                                  :: anti_symmetric

    INTEGER                                  :: i, j, ma, mb, mc, md, &
                                                offset_ac, offset_ad, &
                                                offset_bc, offset_bd

    IF( jatom >=  latom ) THEN
      i = 1
      offset_bd = offset_bd_set(jset,lset) + atomic_offset_bd - 1
      j = offset_bd
      DO md=1,md_max
        DO mb=1,mb_max
          pbd(i) = +density(j)
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
          pbd(i) = -density(j)
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
          pbc(i) = -density(j)
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
          pad(i) = -density(j)
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
          pac(i) = +density(j)
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
          pac(i) = -density(j)
          i = i + 1
          j = j + mc_max
        END DO
      END DO
    END IF

    CALL contract_block(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,fac*symm_fac)

    IF(jatom>=latom) THEN
      i = 1
      j = offset_bd
      DO md=1,md_max
        DO mb=1,mb_max
          !$OMP ATOMIC
          ks(j) = ks(j) + kbd(i)
          i = i + 1
          j = j + 1
        END DO
      END DO
    ELSE
      i = 1
      DO md=1,md_max
        j = offset_bd+md-1
        DO mb=1,mb_max
          !$OMP ATOMIC
          ks(j) = ks(j) - kbd(i)
          i = i + 1
          j = j + md_max
        END DO
      END DO
    END IF
    IF(jatom>=katom) THEN
      i = 1
      j = offset_bc
      DO mc=1,mc_max
        DO mb=1,mb_max
          !$OMP ATOMIC
          ks(j) = ks(j) - kbc(i)
          i = i + 1
          j = j + 1
        END DO
      END DO
    ELSE
      i = 1
      DO mc=1,mc_max
        j = offset_bc + mc -1
        DO mb=1,mb_max
          !$OMP ATOMIC
          ks(j) = ks(j) + kbc(i)
          i = i + 1
          j = j + mc_max
        END DO
      END DO
    END IF
    IF(iatom>=latom) THEN
      i = 1
      j = offset_ad
      DO md=1,md_max
        DO ma=1,ma_max
          !$OMP ATOMIC
          ks(j) = ks(j) - kad(i)
          i = i + 1
          j = j + 1
        END DO
      END DO
    ELSE
      i = 1
      DO md=1,md_max
        j = offset_ad + md - 1
        DO ma=1,ma_max
          !$OMP ATOMIC
          ks(j) = ks(j) + kad(i)
          i = i + 1
          j = j + md_max
        END DO
      END DO
    END IF
! XXXXXXXXXXXXXXXXXXXXXXXX
    IF(iatom>=katom) THEN
      i = 1
      j = offset_ac
      DO mc=1,mc_max
        DO ma=1,ma_max
          !$OMP ATOMIC
          ks(j) = ks(j) + kac(i)
          i = i + 1
          j = j + 1
        END DO
      END DO
    ELSE
      i = 1
      DO mc=1,mc_max
        j = offset_ac + mc -1
        DO ma=1,ma_max
          !$OMP ATOMIC
          ks(j) = ks(j) - kac(i)
          i = i + 1
          j = j + mc_max
        END DO
      END DO
    END IF
 
  END SUBROUTINE update_fock_matrix_as

! *****************************************************************************
!> \brief ...
!> \param i ...
!> \param j ...
!> \param k ...
!> \param l ...
!> \param set_offsets ...
!> \param atom_offsets ...
!> \param iset ...
!> \param jset ...
!> \param kset ...
!> \param lset ...
!> \param ma_max ...
!> \param mb_max ...
!> \param mc_max ...
!> \param md_max ...
!> \param prim ...
! *****************************************************************************
  SUBROUTINE print_integrals(i,j,k,l,set_offsets,atom_offsets,iset,jset,kset,lset,ma_max,mb_max,mc_max,md_max,prim)
    INTEGER                                  :: i, j, k, l
    INTEGER, DIMENSION(:, :, :, :), POINTER  :: set_offsets
    INTEGER, DIMENSION(:, :), POINTER        :: atom_offsets
    INTEGER                                  :: iset, jset, kset, lset, &
                                                ma_max, mb_max, mc_max, md_max
    REAL(dp), DIMENSION(ma_max*mb_max*mc_max&
      *md_max), INTENT(IN)                   :: prim

    INTEGER                                  :: iint, ma, mb, mc, md

    iint=0
    DO md = 1,md_max
      DO mc = 1,mc_max
        DO mb = 1,mb_max
          DO ma = 1,ma_max
             iint=iint+1
             IF(ABS(prim(iint)).gt.0.0000000000001)&
             WRITE(99,*)atom_offsets(i,1)+ma+set_offsets(iset,1,i,1)-1,&
                        atom_offsets(j,1)+ma+set_offsets(jset,1,j,1)-1,&
                        atom_offsets(k,1)+ma+set_offsets(kset,1,k,1)-1,&
                        atom_offsets(l,1)+ma+set_offsets(lset,1,l,1)-1,&
                        prim(iint)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE print_integrals

#include "hfx_get_pmax_val.f90"
END MODULE hfx_energy_potential
