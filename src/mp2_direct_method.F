!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate MP2 energy
!> \par History
!>      06.2011 created [Mauro Del Ben]
!> \author Mauro Del Ben
! *****************************************************************************
MODULE mp2_direct_method
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE gamma,                           ONLY: init_md_ftable
  USE hfx_energy_potential,            ONLY: coulomb4
  USE hfx_libint_wrapper_types,        ONLY: lib_int
  USE hfx_load_balance_methods,        ONLY: cost_model,&
                                             p1_energy,&
                                             p2_energy,&
                                             p3_energy
  USE hfx_pair_list_methods,           ONLY: build_pair_list_mp2
  USE hfx_screening_methods,           ONLY: calc_pair_dist_radii,&
                                             calc_screening_functions
  USE hfx_types,                       ONLY: &
       hfx_basis_info_type, hfx_basis_type, hfx_cell_loop_info_type, &
       hfx_create_neighbor_cells, hfx_general_type, hfx_load_balance_type, &
       hfx_memory_type, hfx_pgf_list, hfx_pgf_product_list, &
       hfx_potential_type, hfx_screen_coeff_type, hfx_screening_type, &
       hfx_type, log_zero, pair_set_list_type
  USE input_constants,                 ONLY: do_mp2_potential_TShPSC
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE mathconstants,                   ONLY: zero
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sendrecv,&
                                             mp_sum
  USE mp2_types,                       ONLY: init_TShPSC_lmax,&
                                             mp2_biel_type,&
                                             mp2_type,&
                                             pair_list_type_mp2
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE t_sh_p_s_c,                      ONLY: free_C0,&
                                             init
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC ::  mp2_canonical_direct_single_batch

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_direct_method'

!***

  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param Emp2 ...
!> \param Emp2_Cou ...
!> \param Emp2_ex ...
!> \param mp2_env ...
!> \param qs_env ...
!> \param rho density matrix
!> \param hfx_section input_section HFX
!> \param para_env ...
!> \param mp2_biel ...
!> \param dimen ...
!> \param C ...
!> \param Auto ...
!> \param i_batch_start ...
!> \param Ni_occupied ...
!> \param occupied ...
!> \param elements_ij_proc ...
!> \param ij_list_proc ...
!> \param Nj_occupied ...
!> \param j_batch_start ...
!> \param occupied_beta ...
!> \param C_beta ...
!> \param Auto_beta ...
!> \param Integ_MP2 ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2011 created [Mauro Del Ben]
!> \author Mauro Del Ben
! *****************************************************************************
  SUBROUTINE mp2_canonical_direct_single_batch(Emp2,Emp2_Cou,Emp2_ex,mp2_env,qs_env,rho,hfx_section,para_env,&
                                               mp2_biel,dimen,C,Auto,i_batch_start,Ni_occupied,&
                                               occupied,elements_ij_proc, ij_list_proc,Nj_occupied,j_batch_start,&
                                               occupied_beta,C_beta,Auto_beta,Integ_MP2,error)

    REAL(KIND=dp)                            :: Emp2, Emp2_Cou, Emp2_ex
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: hfx_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mp2_biel_type)                      :: mp2_biel
    INTEGER                                  :: dimen
    REAL(KIND=dp), DIMENSION(dimen, dimen)   :: C
    REAL(KIND=dp), DIMENSION(dimen)          :: Auto
    INTEGER                                  :: i_batch_start, Ni_occupied, &
                                                occupied, elements_ij_proc
    INTEGER, DIMENSION(elements_ij_proc, 2)  :: ij_list_proc
    INTEGER                                  :: Nj_occupied, j_batch_start
    INTEGER, OPTIONAL                        :: occupied_beta
    REAL(KIND=dp), DIMENSION(dimen, dimen), &
      OPTIONAL                               :: C_beta
    REAL(KIND=dp), DIMENSION(dimen), &
      OPTIONAL                               :: Auto_beta
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :), OPTIONAL        :: Integ_MP2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'mp2_canonical_direct_single_batch', &
      routineP = moduleN//':'//routineN

    INTEGER :: bits_max_val, cache_size, case_index, case_send_receive, &
      counter_proc, elements_ij_proc_rec, elements_kl_proc, global_counter, &
      handle, i, i_list_ij, i_list_kl, i_set_list_ij, i_set_list_ij_start, &
      i_set_list_ij_stop, i_set_list_kl, i_set_list_kl_start, &
      i_set_list_kl_stop, i_start, i_thread, iatom, iatom_end, iatom_start, &
      iiB, ij_elem_max, ikind, index_ij_rec, index_ij_send, index_kl, &
      index_proc_ij, index_proc_shift, irep, iset, jatom, jatom_end, &
      jatom_start, jjB, jkind, jset, katom, katom_end, katom_start, kkB, &
      kkind, kset, l_max, latom, latom_end, latom_start, lkind, llB, lset
    INTEGER :: max_am, max_num_call_sec_transf, max_pgf, max_set, multiple, &
      my_num_call_sec_transf, n_threads, natom, ncob, ncos_max, nints, nkind, &
      nneighbors, nseta, nsetb, nsgf_max, nspins, primitive_counter, &
      proc_num, proc_receive, proc_send, R_offset_rec, Rsize_rec, &
      S_offset_rec, same_size_kl_index, sgfb, sphi_a_u1, sphi_a_u2, &
      sphi_a_u3, sphi_b_u1, sphi_b_u2, sphi_b_u3, sphi_c_u1, sphi_c_u2, &
      sphi_c_u3, sphi_d_u1, sphi_d_u2, sphi_d_u3, Ssize_rec, stat, step_size, &
      total_num_RS_task, unit_id
    INTEGER(int_8) :: estimate_to_store_int, max_val_memory, mem_eris, &
      mem_eris_disk, mem_max_val, n_processes, ncpu, neris_disk, &
      neris_incore, neris_onthefly, neris_tmp, neris_total, nprim_ints, &
      shm_neris_disk, shm_neris_incore, shm_neris_onthefly, shm_neris_total, &
      shm_nprim_ints, shm_stor_count_int_disk, shm_stor_count_max_val, &
      shm_storage_counter_integrals
    INTEGER, ALLOCATABLE, DIMENSION(:) :: kind_of, last_sgf_global, nimages, &
      proc_map, proc_num_task, same_size_kl_elements_counter
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: kl_list_proc, &
                                                task_counter_RS, &
                                                task_counter_RS_temp
    INTEGER, DIMENSION(4)                    :: RS_counter_temp
    INTEGER, DIMENSION(5)                    :: size_parameter_rec, &
                                                size_parameter_send
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, npgfa, npgfb, npgfc, npgfd, nsgfa, nsgfb, &
      nsgfc, nsgfd, shm_block_offset
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb, nsgfl_a, nsgfl_b, &
                                                nsgfl_c, nsgfl_d, &
                                                shm_atomic_block_offset
    INTEGER, DIMENSION(:, :), POINTER, SAVE  :: shm_is_assoc_atomic_block
    INTEGER, DIMENSION(:, :, :, :), POINTER  :: shm_set_offset
    INTEGER, SAVE                            :: shm_number_of_p_entries
    LOGICAL                                  :: alpha_beta_case, &
                                                copy_integrals, do_periodic, &
                                                failure, treat_lsd_in_core
    LOGICAL, DIMENSION(:, :), POINTER        :: shm_atomic_pair_list
    REAL(KIND=dp) :: cartesian_estimate, coeffs_kind_max0, &
      compression_factor, compression_factor_disk, cost_tmp, eps_schwarz, &
      eps_storage, hf_fraction, ln_10, log10_eps_schwarz, log10_pmax, &
      max_contraction_val, max_val1, max_val2, max_val2_set, pmax_atom, &
      pmax_entry, ra(3), rab2, rb(3), rc(3), rcd2, rd(3), screen_kind_ij, &
      screen_kind_kl
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cost_RS, cost_RS_temp, &
      ee_buffer1, ee_buffer2, ee_primitives_tmp, ee_work, ee_work2, &
      primitive_integrals
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: BIb_RS_mat_rec, C_beta_T
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb, BIb_RS_mat_rec_big, &
                                                zero_mat_big
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: BI1, MNRS
    REAL(KIND=dp), DIMENSION(:), POINTER     :: p_work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: max_contraction, &
                                                shm_pmax_block, sphi_b, zeta, &
                                                zetb, zetc, zetd
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: sphi_a_ext_set, &
                                                sphi_b_ext_set, &
                                                sphi_c_ext_set, sphi_d_ext_set
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: sphi_a_ext, sphi_b_ext, &
                                                sphi_c_ext, sphi_d_ext
    REAL(KIND=dp), DIMENSION(dimen, 2)       :: zero_mat
    REAL(KIND=dp), DIMENSION(dimen, dimen)   :: C_T
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(hfx_basis_info_type), POINTER       :: basis_info
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_cell_loop_info_type)            :: cell_loop_info
    TYPE(hfx_general_type)                   :: general_parameter
    TYPE(hfx_load_balance_type), POINTER     :: load_balance_parameter
    TYPE(hfx_memory_type), POINTER           :: memory_parameter
    TYPE(hfx_pgf_list), ALLOCATABLE, &
      DIMENSION(:)                           :: pgf_list_ij, pgf_list_kl
    TYPE(hfx_pgf_product_list), &
      ALLOCATABLE, DIMENSION(:)              :: pgf_product_list
    TYPE(hfx_potential_type)                 :: mp2_potential_parameter, &
                                                potential_parameter
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: screen_coeffs_kind, tmp_R_1, &
                                                tmp_R_2, tmp_screen_pgf1, &
                                                tmp_screen_pgf2
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: screen_coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :, :, :), POINTER   :: radii_pgf, screen_coeffs_pgf
    TYPE(hfx_screening_type)                 :: screening_parameter
    TYPE(hfx_type), POINTER                  :: actual_x_data, &
                                                shm_master_x_data
    TYPE(lib_int)                            :: private_lib
    TYPE(pair_list_type_mp2)                 :: list_ij, list_kl
    TYPE(pair_set_list_type), ALLOCATABLE, &
      DIMENSION(:)                           :: set_list_ij, set_list_kl
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CALL timeset(routineN,handle)
    NULLIFY(dft_control)

    ! The Integ_MP2 will contain the (ia|jb) integrals, necessary for example
    ! for the RI-MP2 basis optimization. In this case the number of ij batches 
    ! has to be equal to 1 (all integrals over molecular orbitals are compituted
    ! in a sigle step).
    copy_integrals=.FALSE.
    IF(PRESENT(Integ_MP2)) copy_integrals=.TRUE.

    alpha_beta_case=.FALSE.

    irep=1

    CALL qs_rho_get(rho, rho_ao=rho_ao, error=error)

    logger => cp_error_get_logger(error)

    CALL get_qs_env(qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    error=error)

    ! kpoint info for coulomb4
    cell_loop_info%dokp = .FALSE.

    !! Calculate l_max used in fgamma , because init_md_ftable is definitely not thread safe
    nkind = SIZE(atomic_kind_set,1)
    l_max = 0
    DO ikind=1,nkind
      l_max = MAX(l_max,MAXVAL(qs_env%x_data(1,1)%basis_parameter(ikind)%lmax))
    ENDDO
    l_max = 4*l_max
    CALL init_md_ftable(l_max)

    IF(mp2_env%potential_parameter%potential_type == do_mp2_potential_TShPSC) THEN
      IF(l_max>init_TShPSC_lmax) THEN
        IF( para_env%mepos == 0 ) THEN
          CALL open_file(unit_number=unit_id,file_name=mp2_env%potential_parameter%filename)
        END IF
        CALL init(l_max, unit_id, para_env%mepos, para_env%group)
        IF(para_env%mepos == 0 ) THEN
          CALL close_file(unit_id)
        END IF
        init_TShPSC_lmax = l_max
      END IF
      mp2_potential_parameter%cutoff_radius=mp2_env%potential_parameter%truncation_radius/2.0_dp
    END IF

    mp2_potential_parameter%potential_type=mp2_env%potential_parameter%potential_type

    n_threads = 1

    shm_neris_total = 0
    shm_nprim_ints = 0
    shm_neris_onthefly = 0
    shm_storage_counter_integrals = 0
    shm_stor_count_int_disk = 0
    shm_neris_incore = 0
    shm_neris_disk = 0
    shm_stor_count_max_val = 0

    ln_10 = LOG(10.0_dp)
    i_thread = 0

    actual_x_data => qs_env%x_data(irep, i_thread + 1)

    shm_master_x_data => qs_env%x_data(irep,1)



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

    max_set = basis_info%max_set
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
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

    ! ** Set pointers
    shm_number_of_p_entries = shm_master_x_data%number_of_p_entries
    shm_is_assoc_atomic_block => shm_master_x_data%is_assoc_atomic_block
    shm_atomic_block_offset => shm_master_x_data%atomic_block_offset
    shm_set_offset => shm_master_x_data%set_offset
    shm_block_offset => shm_master_x_data%block_offset

    !!! Skipped part

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


    !!!!!!!! Missing part on the density matrix

    !! Initialize schwarz screening matrices for near field estimates and boxing screening matrices
    !! for far field estimates. The update is only performed if the geomtry of the system changed.
    !! If the system is periodic, then the corresponding routines are called and some variables
    !! are initialized


    IF( .NOT. shm_master_x_data%screen_funct_is_initialized ) THEN
      CALL calc_pair_dist_radii(qs_env, basis_parameter,&
                                shm_master_x_data%pair_dist_radii_pgf, max_set, max_pgf, eps_schwarz,&
                                n_threads, i_thread, error)
      CALL calc_screening_functions(qs_env, basis_parameter, private_lib, shm_master_x_data%potential_parameter,&
                                    shm_master_x_data%screen_funct_coeffs_set,&
                                    shm_master_x_data%screen_funct_coeffs_kind, &
                                    shm_master_x_data%screen_funct_coeffs_pgf, &
                                    shm_master_x_data%pair_dist_radii_pgf,&
                                    max_set, max_pgf, n_threads, i_thread, p_work, error)

      shm_master_x_data%screen_funct_is_initialized = .TRUE.
   END IF

    screen_coeffs_set   => shm_master_x_data%screen_funct_coeffs_set
    screen_coeffs_kind  => shm_master_x_data%screen_funct_coeffs_kind
    screen_coeffs_pgf   => shm_master_x_data%screen_funct_coeffs_pgf
    radii_pgf           => shm_master_x_data%pair_dist_radii_pgf

!!!!!!!!!
    ALLOCATE(list_ij%elements(natom**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(list_kl%elements(natom**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!!!!!!!!!

    coeffs_kind_max0=MAXVAL(screen_coeffs_kind(:,:)%x(2))
    ALLOCATE(set_list_ij((max_set*natom)**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(set_list_kl((max_set*natom)**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !! precalculate maximum density matrix elements in blocks
    actual_x_data%pmax_block = 0.0_dp
    shm_pmax_block => actual_x_data%pmax_block

    shm_atomic_pair_list => actual_x_data%atomic_pair_list

    iatom_start=1
    iatom_end=natom
    jatom_start=1
    jatom_end=natom
    katom_start=1
    katom_end=natom
    latom_start=1
    latom_end=natom

    CALL build_pair_list_mp2(natom, list_ij,set_list_ij,iatom_start, iatom_end, &
                     jatom_start, jatom_end, &
                     kind_of,basis_parameter,particle_set, &
                     do_periodic,screen_coeffs_set,screen_coeffs_kind, &
                     coeffs_kind_max0,log10_eps_schwarz,cell, 0.D+00, &
                     shm_atomic_pair_list)



    CALL build_pair_list_mp2(natom, list_kl,set_list_kl, katom_start, katom_end, &
                     latom_start, latom_end, &
                     kind_of,basis_parameter,particle_set, &
                     do_periodic,screen_coeffs_set,screen_coeffs_kind, &
                     coeffs_kind_max0,log10_eps_schwarz, cell, 0.D+00, &
                     shm_atomic_pair_list)

     ALLOCATE(BIb(dimen,dimen,elements_ij_proc))
     BIb=0.0D+00
     C_T=TRANSPOSE(C)

     IF(PRESENT(occupied_beta).AND.PRESENT(C_beta).AND.PRESENT(Auto_beta)) THEN
       ALLOCATE(C_beta_T(dimen,dimen))
       C_beta_T(:,:)=TRANSPOSE(C_beta)
       alpha_beta_case=.TRUE.
     END IF

     ij_elem_max=elements_ij_proc
     CALL mp_max(ij_elem_max,para_env%group)

     ! proc_map, vector that replicate the processor numbers also
     ! for negative and positive number > num_pe
     ! needed to know which is the processor, to respect to another one,
     ! for a given shift
     ALLOCATE(proc_map(-para_env%num_pe:2*para_env%num_pe-1))
     DO iiB=0,para_env%num_pe-1
       proc_map(iiB)=iiB
       proc_map(-iiB-1)=para_env%num_pe-iiB-1
       proc_map(para_env%num_pe+iiB)=iiB
     END DO

     ! calculate the minimum multiple of num_pe >= to Ni_occupied*occupied, in such a way
     ! that the i, j loop is performed exactly the same number of time for each procewssor
     multiple=0
     DO
       multiple=multiple+para_env%num_pe
       IF(multiple>=Ni_occupied*Nj_occupied) EXIT
     END DO


     ! proc_num_task contains the numer of time second occupied
     ! orbital transformation is called for each processor, needs for balancing
     ! the point to point send
     ALLOCATE(proc_num_task(0:para_env%num_pe-1),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     proc_num_task=0

     ! step_size=nsgf_max

     counter_proc=0
     !!!! DO iiB=nsgf_max**2+1, 0, -step_size
       DO i_list_ij=1,list_ij%n_element
         iatom=list_ij%elements(i_list_ij)%pair(1)
         jatom=list_ij%elements(i_list_ij)%pair(2)
         i_set_list_ij_start=list_ij%elements(i_list_ij)%set_bounds(1)
         i_set_list_ij_stop=list_ij%elements(i_list_ij)%set_bounds(2)
         ikind=list_ij%elements(i_list_ij)%kind_pair(1)
         jkind=list_ij%elements(i_list_ij)%kind_pair(2)

         nsgfb => basis_parameter(jkind)%nsgf
         nsgfa => basis_parameter(ikind)%nsgf

         ! proc_num=MOD(i_list_ij,para_env%num_pe)

         DO i_set_list_ij=i_set_list_ij_start, i_set_list_ij_stop
           iset=set_list_ij(i_set_list_ij)%pair(1)
           jset=set_list_ij(i_set_list_ij)%pair(2)
           IF(iatom==jatom .AND. jset<iset) CYCLE

           !!!! IF(nsgfa(iset)*nsgfb(jset)>iiB-step_size .AND. nsgfa(iset)*nsgfb(jset)<=iiB) THEN
             counter_proc=counter_proc+1
             proc_num=MOD(counter_proc,para_env%num_pe)

             proc_num_task(proc_num)=proc_num_task(proc_num)+1
           !!!! END IF

          END DO
       END DO
     !!!! END DO
     ! calculate the exact maximum  number of call to the second occupied
     ! orbital transformation
     ! max_num_call_sec_transf=MAXVAL(proc_num_task)

     ! distribute the RS pair over all processor
     ALLOCATE(kl_list_proc(proc_num_task(para_env%mepos),3),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     kl_list_proc=0

     !!!! ALLOCATE(same_size_kl_elements_counter((nsgf_max**2+1)/step_size+1))
     !!!! same_size_kl_elements_counter=0

     !!!! same_size_kl_index=0
     counter_proc=0
     elements_kl_proc=0
     !!!! DO iiB=nsgf_max**2+1, 0, -step_size
       DO i_list_ij=1,list_ij%n_element
         iatom=list_ij%elements(i_list_ij)%pair(1)
         jatom=list_ij%elements(i_list_ij)%pair(2)
         i_set_list_ij_start=list_ij%elements(i_list_ij)%set_bounds(1)
         i_set_list_ij_stop=list_ij%elements(i_list_ij)%set_bounds(2)
         ikind=list_ij%elements(i_list_ij)%kind_pair(1)
         jkind=list_ij%elements(i_list_ij)%kind_pair(2)

         nsgfb => basis_parameter(jkind)%nsgf
         nsgfa => basis_parameter(ikind)%nsgf

         ! proc_num=MOD(i_list_ij,para_env%num_pe)

         DO i_set_list_ij=i_set_list_ij_start, i_set_list_ij_stop
           iset=set_list_ij(i_set_list_ij)%pair(1)
           jset=set_list_ij(i_set_list_ij)%pair(2)
           IF(iatom==jatom .AND. jset<iset) CYCLE

           !!!! IF(nsgfa(iset)*nsgfb(jset)>iiB-step_size .AND. nsgfa(iset)*nsgfb(jset)<=iiB) THEN
             counter_proc=counter_proc+1
             proc_num=MOD(counter_proc,para_env%num_pe)

             IF(proc_num==para_env%mepos) THEN
               elements_kl_proc=elements_kl_proc+1
               kl_list_proc(elements_kl_proc,1)=i_list_ij
               kl_list_proc(elements_kl_proc,2)=i_set_list_ij
               kl_list_proc(elements_kl_proc,3)=counter_proc
             END IF
           !!!! END IF

          END DO
       END DO
       !!!! same_size_kl_index=same_size_kl_index+1
       !!!! same_size_kl_elements_counter(same_size_kl_index)=elements_kl_proc
     !!!! END DO

     total_num_RS_task=SUM(proc_num_task)
     ALLOCATE(task_counter_RS(total_num_RS_task,4),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     ALLOCATE(cost_RS(total_num_RS_task),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     task_counter_RS=0
     cost_RS=0.0_dp

     DO case_index=1, 2

       my_num_call_sec_transf=0
       DO index_kl=1, elements_kl_proc

              i_list_ij=kl_list_proc(index_kl,1)
              i_set_list_ij=kl_list_proc(index_kl,2)

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

       !!!!!     DO i_set_list_ij=i_set_list_ij_start, i_set_list_ij_stop
                iset=set_list_ij(i_set_list_ij)%pair(1)
                jset=set_list_ij(i_set_list_ij)%pair(2)

       !!!!!         IF(iatom==jatom .AND. jset<iset) CYCLE

                ncob = npgfb(jset)*ncoset(lb_max(jset))
                max_val1 = screen_coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + &
                           screen_coeffs_set(jset,iset,jkind,ikind)%x(2)


                sphi_a_ext_set => sphi_a_ext(:,:,:,iset)
                sphi_b_ext_set => sphi_b_ext(:,:,:,jset)

                IF(case_index==1) THEN
                  global_counter=kl_list_proc(index_kl,3)
                  task_counter_RS(global_counter,1)=i_list_ij
                  task_counter_RS(global_counter,2)=i_set_list_ij
                  task_counter_RS(global_counter,3)=nsgfb(jset)*nsgfa(iset)
                END IF

                IF(ALLOCATED(BI1)) DEALLOCATE(BI1)
                ALLOCATE(BI1(dimen,Ni_occupied,nsgfb(jset),nsgfa(iset)),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                BI1=0.D+00

            !!!t1=m_walltime()

            DO i_list_kl=1,list_kl%n_element
            !!!!!DO i_list_kl=1+para_env%mepos,list_kl%n_element,para_env%num_pe

              katom=list_kl%elements(i_list_kl)%pair(1)
              latom=list_kl%elements(i_list_kl)%pair(2)

              i_set_list_kl_start=list_kl%elements(i_list_kl)%set_bounds(1)
              i_set_list_kl_stop=list_kl%elements(i_list_kl)%set_bounds(2)
              kkind=list_kl%elements(i_list_kl)%kind_pair(1)
              lkind=list_kl%elements(i_list_kl)%kind_pair(2)
              rc=list_kl%elements(i_list_kl)%r1
              rd=list_kl%elements(i_list_kl)%r2
              rcd2=list_kl%elements(i_list_kl)%dist2

              pmax_atom = 0.0_dp

              screen_kind_ij = screen_coeffs_kind(jkind,ikind)%x(1)*rab2+&
                               screen_coeffs_kind(jkind,ikind)%x(2)
              screen_kind_kl = screen_coeffs_kind(lkind,kkind)%x(1)*rcd2+&
                               screen_coeffs_kind(lkind,kkind)%x(2)

              !!!!! Change the loop order
              IF( max_val1 + screen_kind_kl + pmax_atom < log10_eps_schwarz) CYCLE
              !!!!!
              IF( screen_kind_ij + screen_kind_kl + pmax_atom < log10_eps_schwarz ) CYCLE

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

                DO i_set_list_kl=i_set_list_kl_start, i_set_list_kl_stop
                  kset=set_list_kl(i_set_list_kl)%pair(1)
                  lset=set_list_kl(i_set_list_kl)%pair(2)

                  IF(katom==latom .AND. lset<kset) CYCLE

                  max_val2_set = (screen_coeffs_set(lset,kset,lkind,kkind)%x(1)*rcd2 + &
                                  screen_coeffs_set(lset,kset,lkind,kkind)%x(2) )
                  max_val2 = max_val1 + max_val2_set

                  !! Near field screening
                  IF(max_val2 + pmax_atom <log10_eps_schwarz) CYCLE
                  sphi_c_ext_set => sphi_c_ext(:,:,:,kset)
                  sphi_d_ext_set => sphi_d_ext(:,:,:,lset)
                  !! get max_vals if we screen on initial density
                  pmax_entry = 0.0_dp

                  log10_pmax = pmax_entry
                  max_val2 = max_val2 + log10_pmax
                  IF(max_val2<log10_eps_schwarz) CYCLE
                  pmax_entry = EXP(log10_pmax*ln_10)

                  IF(case_index==2) THEN
                    IF(ALLOCATED(MNRS)) DEALLOCATE(MNRS)
                    ALLOCATE(MNRS(nsgfd(lset),nsgfc(kset),nsgfb(jset),nsgfa(iset)),STAT=stat)
                    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                    MNRS=0.D+00

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
                                    mp2_potential_parameter, &
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
                                    nimages,do_periodic, p_work, cell_loop_info)
                  !END IF

                    nints = nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                    neris_total = neris_total + nints
                    nprim_ints = nprim_ints + neris_tmp
                    IF(cartesian_estimate == 0.0_dp) cartesian_estimate = TINY(cartesian_estimate)
                    estimate_to_store_int = EXPONENT(cartesian_estimate)
                    estimate_to_store_int = MAX(estimate_to_store_int,-15_int_8)
                    cartesian_estimate = SET_EXPONENT(1.0_dp,estimate_to_store_int+1)


                   IF(cartesian_estimate<eps_schwarz) CYCLE

                  !IF(case_index==2) THEN
                      primitive_counter=0
                      DO llB=1, nsgfd(lset)
                        DO kkB=1, nsgfc(kset)
                          DO jjB=1, nsgfb(jset)
                            DO iiB=1, nsgfa(iset)
                              primitive_counter=primitive_counter+1
                              MNRS(llB,kkB,jjB,iiB)=primitive_integrals(primitive_counter)
                            END DO
                          END DO
                        END DO
                      END DO

                      CALL transform_occupied_orbitals_first(dimen,iatom,jatom,katom,latom,&
                                                 iset,jset,kset,lset,&
                                                 nsgfa(iset),nsgfb(jset),nsgfc(kset),nsgfd(lset),&
                                                 i_batch_start,Ni_occupied,&
                                                 MNRS,C_T,mp2_biel,BI1,error)
                   ELSE
                      task_counter_RS(global_counter,4)=task_counter_RS(global_counter,4)+1

                      cost_tmp = 0.0_dp
                      cost_tmp = cost_model(nsgfd(lset),nsgfc(kset),nsgfb(jset),nsgfa(iset), &
                                                          npgfd(lset),npgfc(kset),npgfb(jset),npgfa(iset),&
                                                          max_val2/log10_eps_schwarz,&
                                                          p1_energy,p2_energy,p3_energy)
                      cost_RS(global_counter)=cost_RS(global_counter)+cost_tmp
                   END IF


                 END DO ! i_set_list_kl
              END DO ! i_list_kl

              IF(case_index==2) THEN
                my_num_call_sec_transf=my_num_call_sec_transf+1
                IF(.NOT.alpha_beta_case) THEN
                  IF(.NOT. mp2_env%direct_canonical%big_send) THEN
                    CALL transform_occupied_orbitals_second(dimen,iatom,jatom,iset,jset,&
                                                     nsgfa(iset),nsgfb(jset),Ni_occupied,Nj_occupied,j_batch_start,&
                                                     BI1,C_T,mp2_biel,para_env,elements_ij_proc,&
                                                     multiple,proc_map,BIb,error)
                  ELSE
                    CALL transform_occupied_orbitals_second_big(dimen,iatom,jatom,iset,jset,&
                                                     nsgfa(iset),nsgfb(jset),Ni_occupied,Nj_occupied,j_batch_start,&
                                                     ij_elem_max,BI1,C_T,mp2_biel,para_env,elements_ij_proc,&
                                                     proc_map,BIb,error)
                  END IF
                ELSE
                  IF(.NOT. mp2_env%direct_canonical%big_send) THEN
                    CALL transform_occupied_orbitals_second(dimen,iatom,jatom,iset,jset,&
                                                     nsgfa(iset),nsgfb(jset),Ni_occupied,Nj_occupied,j_batch_start,&
                                                     BI1,C_beta_T,mp2_biel,para_env,elements_ij_proc,&
                                                     multiple,proc_map,BIb,error)
                  ELSE
                    CALL transform_occupied_orbitals_second_big(dimen,iatom,jatom,iset,jset,&
                                                     nsgfa(iset),nsgfb(jset),Ni_occupied,Nj_occupied,j_batch_start,&
                                                     ij_elem_max,BI1,C_beta_T,mp2_biel,para_env,elements_ij_proc,&
                                                     proc_map,BIb,error)
                  END IF
                END IF
              END IF

            !!!!!!! END DO !i_set_list_ij
          END DO !i_list_ij

          IF(case_index==1) THEN
            CALL mp_sum(task_counter_RS,para_env%group)
            CALL mp_sum(cost_RS,para_env%group)
            ALLOCATE(task_counter_RS_temp(total_num_RS_task,4),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

            ALLOCATE(cost_RS_temp(total_num_RS_task),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

            step_size=1
            ALLOCATE(same_size_kl_elements_counter((nsgf_max**2+1)/step_size+1),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

            same_size_kl_elements_counter=0

            same_size_kl_index=0
            global_counter=0
            DO iiB=nsgf_max**2+1, 0, -step_size
              DO jjB=1,total_num_RS_task
                IF(task_counter_RS(jjB,3)>iiB-step_size .AND. task_counter_RS(jjB,3)<=iiB) THEN
                  global_counter=global_counter+1
                  task_counter_RS_temp(global_counter,1:4)=task_counter_RS(jjB,1:4)
                  cost_RS_temp(global_counter)=cost_RS(jjB)
                END IF
              END DO
              same_size_kl_index=same_size_kl_index+1
              same_size_kl_elements_counter(same_size_kl_index)=global_counter
            END DO

            DEALLOCATE(task_counter_RS)
            DEALLOCATE(cost_RS)

            i_start=1
            DO same_size_kl_index=1, SIZE(same_size_kl_elements_counter)
              DO iiB=i_start, same_size_kl_elements_counter(same_size_kl_index)
                DO jjB=iiB+1, same_size_kl_elements_counter(same_size_kl_index)

                  IF(cost_RS_temp(jjB)>=cost_RS_temp(iiB)) THEN
                    RS_counter_temp=task_counter_RS_temp(iiB,1:4)
                    task_counter_RS_temp(iiB,1:4)=task_counter_RS_temp(jjB,1:4)
                    task_counter_RS_temp(jjB,1:4)=RS_counter_temp

                    cost_tmp=cost_RS_temp(iiB)
                    cost_RS_temp(iiB)=cost_RS_temp(jjB)
                    cost_RS_temp(jjB)=cost_tmp
                  END IF
                END DO
              END DO
              i_start=same_size_kl_elements_counter(same_size_kl_index)+1
            END DO

            proc_num_task=0
            DO counter_proc=1,total_num_RS_task
              proc_num=MOD(counter_proc,para_env%num_pe)
              proc_num_task(proc_num)=proc_num_task(proc_num)+1
            END DO

            max_num_call_sec_transf=MAXVAL(proc_num_task)

            DEALLOCATE(kl_list_proc)
            ALLOCATE(kl_list_proc(proc_num_task(para_env%mepos),2),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

            kl_list_proc=0

            elements_kl_proc=0
            DO counter_proc=1,total_num_RS_task
              proc_num=MOD(counter_proc,para_env%num_pe)
              IF(proc_num==para_env%mepos) THEN
                elements_kl_proc=elements_kl_proc+1
                kl_list_proc(elements_kl_proc,1)=task_counter_RS_temp(counter_proc,1)
                kl_list_proc(elements_kl_proc,2)=task_counter_RS_temp(counter_proc,2)
              END IF
            END DO

            DEALLOCATE(task_counter_RS_temp)
            DEALLOCATE(cost_RS_temp)
          END IF
        END DO ! case_index

        size_parameter_send(1)=1
        size_parameter_send(2)=1
        size_parameter_send(3)=0
        size_parameter_send(4)=0
        size_parameter_send(5)=elements_ij_proc

        IF(mp2_env%direct_canonical%big_send) THEN
          ALLOCATE(zero_mat_big(dimen,2,ij_elem_max),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

        END IF

        DO iiB=my_num_call_sec_transf+1, max_num_call_sec_transf
          DO index_proc_shift=0, para_env%num_pe-1

            proc_send=proc_map(para_env%mepos+index_proc_shift)
            proc_receive=proc_map(para_env%mepos-index_proc_shift)

          !DO index_proc_shift=0, para_env%num_pe-1
          !  proc_send=proc_map(0-para_env%mepos+index_proc_shift)
          !  proc_receive=proc_send

            IF(proc_send/=para_env%mepos) THEN
              ! the processor starts to send (and receive?)

              CALL mp_sendrecv(size_parameter_send,proc_send,size_parameter_rec,proc_receive,para_env%group)

              Rsize_rec=size_parameter_rec(1)
              Ssize_rec=size_parameter_rec(2)
              R_offset_rec=size_parameter_rec(3)
              S_offset_rec=size_parameter_rec(4)
              elements_ij_proc_rec=size_parameter_rec(5)
              IF(.NOT. mp2_env%direct_canonical%big_send) THEN
                ALLOCATE(BIb_RS_mat_rec(dimen,Rsize_rec+Ssize_rec),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              ELSE
                ALLOCATE(BIb_RS_mat_rec_big(dimen,Rsize_rec+Ssize_rec,ij_elem_max),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              END IF
            ELSE
              elements_ij_proc_rec=elements_ij_proc
            END IF

            IF(.NOT. mp2_env%direct_canonical%big_send) THEN
              index_ij_send=0
              index_ij_rec=0
              DO index_proc_ij=proc_send+1,multiple,para_env%num_pe

                zero_mat=0.D+00
                case_send_receive=-1
                IF(proc_send==para_env%mepos) THEN
                  case_send_receive=0
                ELSE
                  case_send_receive=1
                END IF

                SELECT CASE(case_send_receive)
                  CASE(0)
                    ! do nothing
                  CASE(1)

                    CALL mp_sendrecv(zero_mat,proc_send,BIb_RS_mat_rec,proc_receive,para_env%group)

                    index_ij_rec=index_ij_rec+1
                    IF(index_ij_rec<=elements_ij_proc .AND. elements_ij_proc>0) THEN

                      BIb(1:dimen,R_offset_rec+1:R_offset_rec+Rsize_rec,index_ij_rec)= &
                      BIb(1:dimen,R_offset_rec+1:R_offset_rec+Rsize_rec,index_ij_rec)+ &
                      BIb_RS_mat_rec(1:dimen,1:Rsize_rec)

                      BIb(1:dimen,S_offset_rec+1:S_offset_rec+Ssize_rec,index_ij_rec)= &
                      BIb(1:dimen,S_offset_rec+1:S_offset_rec+Ssize_rec,index_ij_rec)+ &
                      BIb_RS_mat_rec(1:dimen,Rsize_rec+1:Rsize_rec+Ssize_rec)

                    END IF
                  CASE DEFAULT
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                END SELECT

              END DO
            ELSE
              zero_mat_big=0.D+00
              case_send_receive=-1
              IF(proc_send==para_env%mepos) THEN
                case_send_receive=0
              ELSE
                case_send_receive=1
              END IF

              SELECT CASE(case_send_receive)
                CASE(0)
                  ! do nothing
                CASE(1)

                  CALL mp_sendrecv(zero_mat_big,proc_send,BIb_RS_mat_rec_big,proc_receive,para_env%group)

                    BIb(1:dimen,R_offset_rec+1:R_offset_rec+Rsize_rec,1:elements_ij_proc)= &
                    BIb(1:dimen,R_offset_rec+1:R_offset_rec+Rsize_rec,1:elements_ij_proc)+ &
                    BIb_RS_mat_rec_big(1:dimen,1:Rsize_rec,1:elements_ij_proc)

                    BIb(1:dimen,S_offset_rec+1:S_offset_rec+Ssize_rec,1:elements_ij_proc)= &
                    BIb(1:dimen,S_offset_rec+1:S_offset_rec+Ssize_rec,1:elements_ij_proc)+ &
                    BIb_RS_mat_rec_big(1:dimen,Rsize_rec+1:Rsize_rec+Ssize_rec,1:elements_ij_proc)

                CASE DEFAULT
                  CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
              END SELECT
            END IF

            IF(proc_send/=para_env%mepos) THEN
              IF(.NOT. mp2_env%direct_canonical%big_send) THEN
                DEALLOCATE(BIb_RS_mat_rec)
              ELSE
                DEALLOCATE(BIb_RS_mat_rec_big)
              END IF
            END IF

          END DO
        END DO

        IF(mp2_env%direct_canonical%big_send) THEN
          DEALLOCATE(zero_mat_big)
        END IF

        logger => cp_error_get_logger(error)

    DEALLOCATE(primitive_integrals,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF(.NOT.alpha_beta_case) THEN
      CALL transform_virtual_orbitals_and_accumulate(dimen,occupied,dimen-occupied,i_batch_start,&
                                                     j_batch_start,BIb,C,Auto,elements_ij_proc,ij_list_proc,&
                                                     nspins,Emp2,Emp2_Cou,Emp2_ex,error=error)
    ELSE
      CALL transform_virtual_orbitals_and_accumulate_ABcase(dimen,occupied,occupied_beta,dimen-occupied,dimen-occupied_beta,&
                                                     i_batch_start,j_batch_start,&
                                                     BIb,C,C_beta,Auto,Auto_beta,&
                                                     elements_ij_proc,ij_list_proc,Emp2,Emp2_Cou,error=error)
      DEALLOCATE(C_beta_T)
    END IF

    IF(copy_integrals) THEN
      IF(.NOT.alpha_beta_case) THEN
        ALLOCATE(Integ_MP2(dimen-occupied,dimen-occupied,occupied,occupied),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        Integ_MP2=0.0_dp
        DO i=1, elements_ij_proc
          iiB=ij_list_proc(i,1)
          jjB=ij_list_proc(i,2)
          Integ_MP2(:,:,iiB+i_batch_start,jjB+j_batch_start)=BIb(1:dimen-occupied,1:dimen-occupied,i)
        END DO
      ELSE
        ALLOCATE(Integ_MP2(dimen-occupied,dimen-occupied_beta,occupied,occupied_beta),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        Integ_MP2=0.0_dp
        DO i=1, elements_ij_proc
          iiB=ij_list_proc(i,1)
          jjB=ij_list_proc(i,2)
          Integ_MP2(:,:,iiB+i_batch_start,jjB+j_batch_start)=BIb(1:dimen-occupied,1:dimen-occupied_beta,i)
        END DO
      END IF
    END IF
    DEALLOCATE(BIb)

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

    IF(mp2_env%potential_parameter%potential_type == do_mp2_potential_TShPSC) THEN
      init_TShPSC_lmax = -1
      CALL free_C0()
    END IF

    CALL timestop(handle)

  END SUBROUTINE mp2_canonical_direct_single_batch

! *****************************************************************************
!> \brief ...
!> \param dimen ...
!> \param latom ...
!> \param katom ...
!> \param jatom ...
!> \param iatom ...
!> \param lset ...
!> \param kset ...
!> \param jset ...
!> \param iset ...
!> \param Ssize ...
!> \param Rsize ...
!> \param Nsize ...
!> \param Msize ...
!> \param i_batch_start ...
!> \param Ni_occupied ...
!> \param MNRS ...
!> \param C_T ...
!> \param mp2_biel ...
!> \param BI1 ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE transform_occupied_orbitals_first(dimen,latom,katom,jatom,iatom,&
                                             lset,kset,jset,iset,&
                                             Ssize,Rsize,Nsize,Msize,&
                                             i_batch_start,Ni_occupied,&
                                             MNRS,C_T,mp2_biel,BI1,error)

    INTEGER :: dimen, latom, katom, jatom, iatom, lset, kset, jset, iset, &
      Ssize, Rsize, Nsize, Msize, i_batch_start, Ni_occupied
    REAL(KIND=dp), &
      DIMENSION(Msize, Nsize, Rsize, Ssize)  :: MNRS
    REAL(KIND=dp), DIMENSION(dimen, dimen)   :: C_T
    TYPE(mp2_biel_type)                      :: mp2_biel
    REAL(KIND=dp), DIMENSION(dimen, &
      Ni_occupied, Rsize, Ssize)             :: BI1
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'transform_occupied_orbitals_first', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i_global, m, M_global, &
                                                M_offset, M_start, n, &
                                                N_global, N_offset, r, &
                                                R_offset, R_start, s, S_offset
    REAL(KIND=dp)                            :: MNRS_element

    N_offset=mp2_biel%index_table(jatom,jset)-1
    M_offset=mp2_biel%index_table(iatom,iset)-1
    S_offset=mp2_biel%index_table(latom,lset)-1
    R_offset=mp2_biel%index_table(katom,kset)-1

    DO S=1, Ssize
      R_start=1
      IF(katom==latom .AND. kset==lset) R_start=S
      DO R=R_start, Rsize

        ! fast i don't know why
        DO N=1, Nsize
          N_global=N+N_offset
          M_start=1
          IF(iatom==jatom .AND. iset==jset) THEN
            M=N
            M_global=M+M_offset
            MNRS_element=MNRS(M,N,R,S)
            DO i=1, Ni_occupied
              i_global=i+i_batch_start
              BI1(N_global,i,R,S)=BI1(N_global,i,R,S)+C_T(i_global,M_global)*MNRS_element
            END DO
            M_start=N+1
          END IF

          DO M=M_start, Msize
            M_global=M+M_offset
            MNRS_element=MNRS(M,N,R,S)
            DO i=1, Ni_occupied
               i_global=i+i_batch_start
               BI1(N_global,i,R,S)=BI1(N_global,i,R,S)+C_T(i_global,M_global)*MNRS_element
               BI1(M_global,i,R,S)=BI1(M_global,i,R,S)+C_T(i_global,N_global)*MNRS_element
            END DO
          END DO
        END DO

       END DO
     END DO

  END SUBROUTINE transform_occupied_orbitals_first

! *****************************************************************************
!> \brief ...
!> \param dimen ...
!> \param latom ...
!> \param katom ...
!> \param lset ...
!> \param kset ...
!> \param Ssize ...
!> \param Rsize ...
!> \param Ni_occupied ...
!> \param Nj_occupied ...
!> \param j_batch_start ...
!> \param BI1 ...
!> \param C_T ...
!> \param mp2_biel ...
!> \param para_env ...
!> \param elements_ij_proc ...
!> \param multiple ...
!> \param proc_map ...
!> \param BIb ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE transform_occupied_orbitals_second(dimen,latom,katom,lset,kset,&
                                             Ssize,Rsize,Ni_occupied,Nj_occupied,j_batch_start,&
                                             BI1,C_T,mp2_biel,para_env,&
                                             elements_ij_proc,&
                                             multiple,proc_map,BIb,error)

    INTEGER                                  :: dimen, latom, katom, lset, &
                                                kset, Ssize, Rsize, &
                                                Ni_occupied, Nj_occupied, &
                                                j_batch_start
    REAL(KIND=dp), DIMENSION(dimen, &
      Ni_occupied, Rsize, Ssize)             :: BI1
    REAL(KIND=dp), DIMENSION(dimen, dimen)   :: C_T
    TYPE(mp2_biel_type)                      :: mp2_biel
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: elements_ij_proc, multiple
    INTEGER, DIMENSION((-para_env%num_pe):(2&
      *para_env%num_pe-1))                   :: proc_map
    REAL(KIND=dp), DIMENSION(dimen, dimen, &
      elements_ij_proc)                      :: BIb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'transform_occupied_orbitals_second', &
      routineP = moduleN//':'//routineN

    INTEGER :: case_send_receive, elements_ij_proc_rec, handle, i, &
      index_ij_rec, index_ij_send, index_proc_ij, index_proc_shift, j, n, &
      proc_receive, proc_send, r, R_global, R_offset, R_offset_rec, R_start, &
      Rsize_rec, s, S_global, S_offset, S_offset_rec, Ssize_rec, stat
    REAL(KIND=dp), &
      DIMENSION(dimen, Rsize+Ssize)          :: BIb_RS_mat
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: BIb_RS_mat_rec
    REAL(KIND=dp)                            :: C_T_R, C_T_S
    LOGICAL                                  :: failure
    INTEGER, DIMENSION(5)                    :: size_parameter_rec, &
                                                size_parameter_send

    CALL timeset(routineN,handle)
    failure=.FALSE.

    S_offset=mp2_biel%index_table(latom,lset)-1
    R_offset=mp2_biel%index_table(katom,kset)-1

    size_parameter_send(1)=Rsize
    size_parameter_send(2)=Ssize
    size_parameter_send(3)=R_offset
    size_parameter_send(4)=S_offset
    size_parameter_send(5)=elements_ij_proc

    DO index_proc_shift=0, para_env%num_pe-1

      proc_send=proc_map(para_env%mepos+index_proc_shift)
      proc_receive=proc_map(para_env%mepos-index_proc_shift)

      IF(proc_send/=para_env%mepos) THEN
        ! the processor starts to send (and receive?)

        CALL mp_sendrecv(size_parameter_send,proc_send,size_parameter_rec,proc_receive,para_env%group)

        Rsize_rec=size_parameter_rec(1)
        Ssize_rec=size_parameter_rec(2)
        R_offset_rec=size_parameter_rec(3)
        S_offset_rec=size_parameter_rec(4)
        elements_ij_proc_rec=size_parameter_rec(5)
        ALLOCATE(BIb_RS_mat_rec(dimen,Rsize_rec+Ssize_rec),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      ELSE
        elements_ij_proc_rec=elements_ij_proc
      END IF

      case_send_receive=-1
      IF(proc_send==para_env%mepos) THEN
        case_send_receive=0
      ELSE
        case_send_receive=1
      END IF

      index_ij_send=0
      index_ij_rec=0
      DO index_proc_ij=proc_send+1,multiple,para_env%num_pe

        BIb_RS_mat=zero
        IF(index_proc_ij<=Ni_occupied*Nj_occupied) THEN

          index_ij_send=index_ij_send+1

          i=(index_proc_ij-1)/Nj_occupied+1
          j=index_proc_ij-(i-1)*Nj_occupied+j_batch_start

          DO S=1, Ssize
            S_global=S+S_offset
            R_start=1
            IF(katom==latom .AND. kset==lset) R_start=S
            DO R=R_start, Rsize
               R_global=R+R_offset

                 IF(R_global/=S_global) THEN
                   C_T_R=C_T(j,R_global)
                   C_T_S=C_T(j,S_global)
                   DO N=1, dimen
                     BIb_RS_mat(N,R)=BIb_RS_mat(N,R)+C_T_S*BI1(N,i,R,S)
                   END DO
                   DO N=1, dimen
                     BIb_RS_mat(N,Rsize+S)=BIb_RS_mat(N,Rsize+S)+C_T_R*BI1(N,i,R,S)
                   END DO
                 ELSE
                   C_T_S=C_T(j,S_global)
                   DO N=1, dimen
                     BIb_RS_mat(N,R)=BIb_RS_mat(N,R)+C_T_S*BI1(N,i,R,S)
                   END DO
                 END IF

            END DO
          END DO

         END IF

         SELECT CASE(case_send_receive)
           CASE(0)
             ! the processor is the sender and receiver itself
             IF(index_ij_send<=elements_ij_proc .AND. elements_ij_proc>0) THEN

               BIb(1:dimen,R_offset+1:R_offset+Rsize,index_ij_send)= &
                           BIb(1:dimen,R_offset+1:R_offset+Rsize,index_ij_send)+BIb_RS_mat(1:dimen,1:Rsize)

               BIb(1:dimen,S_offset+1:S_offset+Ssize,index_ij_send)= &
                           BIb(1:dimen,S_offset+1:S_offset+Ssize,index_ij_send)+BIb_RS_mat(1:dimen,Rsize+1:Rsize+Ssize)

             END IF
           CASE(1)

             CALL mp_sendrecv(BIb_RS_mat,proc_send,BIb_RS_mat_rec,proc_receive,para_env%group)

             index_ij_rec=index_ij_rec+1
             IF(index_ij_rec<=elements_ij_proc .AND. elements_ij_proc>0) THEN

               BIb(1:dimen,R_offset_rec+1:R_offset_rec+Rsize_rec,index_ij_rec)= &
               BIb(1:dimen,R_offset_rec+1:R_offset_rec+Rsize_rec,index_ij_rec)+ &
               BIb_RS_mat_rec(1:dimen,1:Rsize_rec)

               BIb(1:dimen,S_offset_rec+1:S_offset_rec+Ssize_rec,index_ij_rec)= &
               BIb(1:dimen,S_offset_rec+1:S_offset_rec+Ssize_rec,index_ij_rec)+ &
               BIb_RS_mat_rec(1:dimen,Rsize_rec+1:Rsize_rec+Ssize_rec)

             END IF
           CASE DEFAULT
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
         END SELECT

      END DO ! loop over the ij of the processor

      IF(proc_send/=para_env%mepos) THEN
        DEALLOCATE(BIb_RS_mat_rec)
      END IF

    END DO ! loop over the processor starting from itself

    CALL timestop(handle)

  END SUBROUTINE transform_occupied_orbitals_second

! *****************************************************************************
!> \brief ...
!> \param dimen ...
!> \param latom ...
!> \param katom ...
!> \param lset ...
!> \param kset ...
!> \param Ssize ...
!> \param Rsize ...
!> \param Ni_occupied ...
!> \param Nj_occupied ...
!> \param j_batch_start ...
!> \param ij_elem_max ...
!> \param BI1 ...
!> \param C_T ...
!> \param mp2_biel ...
!> \param para_env ...
!> \param elements_ij_proc ...
!> \param proc_map ...
!> \param BIb ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE transform_occupied_orbitals_second_big(dimen,latom,katom,lset,kset,&
                                             Ssize,Rsize,Ni_occupied,Nj_occupied,j_batch_start,&
                                             ij_elem_max,BI1,C_T,mp2_biel,para_env,&
                                             elements_ij_proc,&
                                             proc_map,BIb,error)

    INTEGER :: dimen, latom, katom, lset, kset, Ssize, Rsize, Ni_occupied, &
      Nj_occupied, j_batch_start, ij_elem_max
    REAL(KIND=dp), DIMENSION(dimen, &
      Ni_occupied, Rsize, Ssize)             :: BI1
    REAL(KIND=dp), DIMENSION(dimen, dimen)   :: C_T
    TYPE(mp2_biel_type)                      :: mp2_biel
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: elements_ij_proc
    INTEGER, DIMENSION((-para_env%num_pe):(2&
      *para_env%num_pe-1))                   :: proc_map
    REAL(KIND=dp), DIMENSION(dimen, dimen, &
      elements_ij_proc)                      :: BIb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'transform_occupied_orbitals_second_big', &
      routineP = moduleN//':'//routineN

    INTEGER :: case_send_receive, elements_ij_proc_rec, handle, i, &
      index_ij_rec, index_ij_send, index_proc_ij, index_proc_shift, j, n, &
      proc_receive, proc_send, r, R_global, R_offset, R_offset_rec, R_start, &
      Rsize_rec, s, S_global, S_offset, S_offset_rec, Ssize_rec, stat
    REAL(KIND=dp), DIMENSION(dimen, &
      Rsize+Ssize, ij_elem_max)              :: BIb_RS_mat
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_RS_mat_rec
    REAL(KIND=dp)                            :: C_T_R, C_T_S
    LOGICAL                                  :: failure
    INTEGER, DIMENSION(5)                    :: size_parameter_rec, &
                                                size_parameter_send

    CALL timeset(routineN,handle)
    failure=.FALSE.

    S_offset=mp2_biel%index_table(latom,lset)-1
    R_offset=mp2_biel%index_table(katom,kset)-1

    size_parameter_send(1)=Rsize
    size_parameter_send(2)=Ssize
    size_parameter_send(3)=R_offset
    size_parameter_send(4)=S_offset
    size_parameter_send(5)=elements_ij_proc

    DO index_proc_shift=0, para_env%num_pe-1

      proc_send=proc_map(para_env%mepos+index_proc_shift)
      proc_receive=proc_map(para_env%mepos-index_proc_shift)

      IF(proc_send/=para_env%mepos) THEN
        ! the processor starts to send (and receive?)

        CALL mp_sendrecv(size_parameter_send,proc_send,size_parameter_rec,proc_receive,para_env%group)

        Rsize_rec=size_parameter_rec(1)
        Ssize_rec=size_parameter_rec(2)
        R_offset_rec=size_parameter_rec(3)
        S_offset_rec=size_parameter_rec(4)
        elements_ij_proc_rec=size_parameter_rec(5)
        ALLOCATE(BIb_RS_mat_rec(dimen,Rsize_rec+Ssize_rec,ij_elem_max),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ELSE
        elements_ij_proc_rec=elements_ij_proc
      END IF


      index_ij_send=0
      index_ij_rec=0
      BIb_RS_mat=zero

      case_send_receive=-1
      IF(proc_send==para_env%mepos) THEN
        case_send_receive=0
      ELSE
        case_send_receive=1
      END IF

      DO index_proc_ij=proc_send+1,Ni_occupied*Nj_occupied,para_env%num_pe

         index_ij_send=index_ij_send+1

         i=(index_proc_ij-1)/Nj_occupied+1
         j=index_proc_ij-(i-1)*Nj_occupied+j_batch_start

         DO S=1, Ssize
           S_global=S+S_offset
           R_start=1
           IF(katom==latom .AND. kset==lset) R_start=S
           DO R=R_start, Rsize
              R_global=R+R_offset

                IF(R_global/=S_global) THEN
                  C_T_R=C_T(j,R_global)
                  C_T_S=C_T(j,S_global)
                  DO N=1, dimen
                    BIb_RS_mat(N,R,index_ij_send)=BIb_RS_mat(N,R,index_ij_send)+C_T_S*BI1(N,i,R,S)
                  END DO
                  DO N=1, dimen
                    BIb_RS_mat(N,Rsize+S,index_ij_send)=BIb_RS_mat(N,Rsize+S,index_ij_send)+C_T_R*BI1(N,i,R,S)
                  END DO
                ELSE
                  C_T_S=C_T(j,S_global)
                  DO N=1, dimen
                    BIb_RS_mat(N,R,index_ij_send)=BIb_RS_mat(N,R,index_ij_send)+C_T_S*BI1(N,i,R,S)
                  END DO
                END IF

           END DO
         END DO

      END DO

      SELECT CASE(case_send_receive)
      CASE(0)
        ! the processor is the sender and receiver itself
          BIb(1:dimen,R_offset+1:R_offset+Rsize,1:elements_ij_proc)= &
          BIb(1:dimen,R_offset+1:R_offset+Rsize,1:elements_ij_proc)+ &
          BIb_RS_mat(1:dimen,1:Rsize,1:elements_ij_proc)

          BIb(1:dimen,S_offset+1:S_offset+Ssize,1:elements_ij_proc)= &
          BIb(1:dimen,S_offset+1:S_offset+Ssize,1:elements_ij_proc)+ &
          BIb_RS_mat(1:dimen,Rsize+1:Rsize+Ssize,1:elements_ij_proc)

      CASE(1)

        CALL mp_sendrecv(BIb_RS_mat,proc_send,BIb_RS_mat_rec,proc_receive,para_env%group)

          BIb(1:dimen,R_offset_rec+1:R_offset_rec+Rsize_rec,1:elements_ij_proc)= &
          BIb(1:dimen,R_offset_rec+1:R_offset_rec+Rsize_rec,1:elements_ij_proc)+ &
          BIb_RS_mat_rec(1:dimen,1:Rsize_rec,1:elements_ij_proc)

          BIb(1:dimen,S_offset_rec+1:S_offset_rec+Ssize_rec,1:elements_ij_proc)= &
          BIb(1:dimen,S_offset_rec+1:S_offset_rec+Ssize_rec,1:elements_ij_proc)+ &
          BIb_RS_mat_rec(1:dimen,Rsize_rec+1:Rsize_rec+Ssize_rec,1:elements_ij_proc)

      CASE DEFAULT
        CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT

      IF(proc_send/=para_env%mepos) THEN
        DEALLOCATE(BIb_RS_mat_rec)
      END IF

    END DO ! loop over the processor starting from itself

    CALL timestop(handle)

  END SUBROUTINE transform_occupied_orbitals_second_big

! *****************************************************************************
!> \brief ...
!> \param dimen ...
!> \param occupied ...
!> \param virtual ...
!> \param i_batch_start ...
!> \param j_batch_start ...
!> \param BIb ...
!> \param C ...
!> \param Auto ...
!> \param elements_ij_proc ...
!> \param ij_list_proc ...
!> \param nspins ...
!> \param Emp2 ...
!> \param Emp2_Cou ...
!> \param Emp2_ex ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE transform_virtual_orbitals_and_accumulate(dimen,occupied,virtual,i_batch_start,&
                                                   j_batch_start,BIb,C,Auto,elements_ij_proc,&
                                                   ij_list_proc,nspins,Emp2,Emp2_Cou,Emp2_ex,error)

    INTEGER                                  :: dimen, occupied, virtual, &
                                                i_batch_start, j_batch_start
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb
    REAL(KIND=dp), DIMENSION(dimen, dimen)   :: C
    REAL(KIND=dp), DIMENSION(dimen)          :: Auto
    INTEGER                                  :: elements_ij_proc
    INTEGER, DIMENSION(elements_ij_proc, 2)  :: ij_list_proc
    INTEGER                                  :: nspins
    REAL(KIND=dp)                            :: Emp2, Emp2_Cou, Emp2_ex
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'transform_virtual_orbitals_and_accumulate', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: zero = 0.0_dp

    INTEGER                                  :: a, a_global, b, b_global, &
                                                handle, i, i_global, &
                                                index_ij, j, j_global, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: iajb, ibja, parz, two
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: BIa

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ALLOCATE(BIa(dimen,virtual),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    BIa=zero
    DO index_ij=1, elements_ij_proc

      CALL DGEMM('T','N',dimen,virtual,dimen,1.0_dp,Bib(1,1,index_ij),dimen,C(1,occupied+1),dimen,0.0_dp,Bia(1,1),dimen)
      Bib(1:dimen,1:virtual,index_ij)=Bia(1:dimen,1:virtual)

    END DO

    DEALLOCATE(BIa)
    ALLOCATE(BIa(virtual,virtual),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    BIa=zero
    DO index_ij=1, elements_ij_proc

      CALL DGEMM('T','N',virtual,virtual,dimen,1.0_dp,Bib(1,1,index_ij),dimen,C(1,occupied+1),dimen,0.0_dp,&
                 BIa(1,1),virtual)
      BIb(1:virtual,1:virtual,index_ij)=BIa(1:virtual,1:virtual)

    END DO

    two=2.0_dp/nspins
    DO index_ij=1, elements_ij_proc
       i=ij_list_proc(index_ij,1)
       j=ij_list_proc(index_ij,2)
       i_global=i+i_batch_start
       j_global=j+j_batch_start
       DO a=1, virtual
         a_global=a+occupied
         DO b=1, virtual
           b_global=b+occupied
           iajb=BIb(a,b,index_ij)
           ibja=BIb(b,a,index_ij)
           parz=iajb/(Auto(i_global)+Auto(j_global)-Auto(a_global)-Auto(b_global))
           ! parz=parz*(two*iajb-ibja)   !Full
           ! parz=parz*(iajb)            !Coulomb
           ! parz=parz*(ibja)            !Coulomb
           ! Emp2=Emp2+parz/nspins
           Emp2_Cou=Emp2_Cou+parz*two*(iajb)/nspins
           Emp2_ex=Emp2_ex-parz*(ibja)/nspins
           Emp2=Emp2+parz*(two*iajb-ibja)/nspins
         END DO
       END DO
    END DO

    DEALLOCATE(BIa)

    CALL timestop(handle)

  END SUBROUTINE transform_virtual_orbitals_and_accumulate

! *****************************************************************************
!> \brief ...
!> \param dimen ...
!> \param occ_i ...
!> \param occ_j ...
!> \param virt_i ...
!> \param virt_j ...
!> \param i_batch_start ...
!> \param j_batch_start ...
!> \param BIb ...
!> \param C_i ...
!> \param C_j ...
!> \param Auto_i ...
!> \param Auto_j ...
!> \param elements_ij_proc ...
!> \param ij_list_proc ...
!> \param Emp2 ...
!> \param Emp2_Cou ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE transform_virtual_orbitals_and_accumulate_ABcase(dimen,occ_i,occ_j,virt_i,virt_j,i_batch_start,&
                                                   j_batch_start,BIb,C_i,C_j,Auto_i,Auto_j,elements_ij_proc,&
                                                   ij_list_proc,Emp2,Emp2_Cou,error)

    INTEGER                                  :: dimen, occ_i, occ_j, virt_i, &
                                                virt_j, i_batch_start, &
                                                j_batch_start
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb
    REAL(KIND=dp), DIMENSION(dimen, dimen)   :: C_i, C_j
    REAL(KIND=dp), DIMENSION(dimen)          :: Auto_i, Auto_j
    INTEGER                                  :: elements_ij_proc
    INTEGER, DIMENSION(elements_ij_proc, 2)  :: ij_list_proc
    REAL(KIND=dp)                            :: Emp2, Emp2_Cou
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'transform_virtual_orbitals_and_accumulate_ABcase', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: two = 2.D+00, zero = 0.D+00

    INTEGER                                  :: a, a_global, b, b_global, &
                                                handle, i, i_global, &
                                                index_ij, j, j_global, n, s, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: iajb, parz
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: BIa

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ALLOCATE(BIa(dimen,virt_i),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO index_ij=1, elements_ij_proc

      DO a=1, virt_i
        a_global=a+occ_i
        DO S=1, dimen
          parz=zero
          DO N=1, dimen
            parz=parz+C_i(N,a_global)*BIb(N,S,index_ij)
          END DO
          BIa(S,a)=parz
        END DO
      END DO

      BIb(1:dimen,1:virt_i,index_ij)=BIa

    END DO

    DEALLOCATE(BIa)
    ALLOCATE(BIa(virt_i,virt_j),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO index_ij=1, elements_ij_proc

      DO a=1, virt_i
        DO b=1, virt_j
          b_global=b+occ_j
          parz=zero
          DO S=1, dimen
            parz=parz+C_j(S,b_global)*BIb(S,a,index_ij)
          END DO
          BIa(a,b)=parz
        END DO
      END DO

      BIb(1:virt_i,1:virt_j,index_ij)=BIa

    END DO

    DO index_ij=1, elements_ij_proc
       i=ij_list_proc(index_ij,1)
       j=ij_list_proc(index_ij,2)
       i_global=i+i_batch_start
       j_global=j+j_batch_start
       DO a=1, virt_i
         a_global=a+occ_i
         DO b=1, virt_j
           b_global=b+occ_j
           iajb=BIb(a,b,index_ij)
           parz=iajb*iajb/(Auto_i(i_global)+Auto_j(j_global)-Auto_i(a_global)-Auto_j(b_global))
           ! Emp2=Emp2+parz/two
           Emp2_Cou=Emp2_Cou+parz/two
           Emp2=Emp2+parz/two
         END DO
       END DO
    END DO

    DEALLOCATE(BIa)

    CALL timestop(handle)

  END SUBROUTINE transform_virtual_orbitals_and_accumulate_ABcase

END MODULE mp2_direct_method
