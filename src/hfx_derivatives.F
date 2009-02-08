!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines to calculate dervivatives with respect to basis function origin 
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE hfx_derivatives 
  
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_files,                        ONLY: get_unit_number
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE gamma,                           ONLY: init_md_ftable
  USE hfx_communication,               ONLY: get_full_density
  USE hfx_helpers,                     ONLY: FNV64
  USE hfx_libint_interface,            ONLY: evaluate_deriv_eri
  USE hfx_libint_wrapper,              ONLY: initialize_libderiv,&
                                             initialize_libint,&
                                             lib_deriv,&
                                             lib_int
  USE hfx_load_balance_methods,        ONLY: load_balance_forces,&
                                             update_load_balance_forces
  USE hfx_screening_methods,           ONLY: update_pmax_mat
  USE hfx_types,                       ONLY: &
       hfx_basis_info_type, hfx_basis_type, hfx_cell_type, hfx_general_type, &
       hfx_load_balance_type, hfx_potential_type, hfx_screen_coeff_type, &
       hfx_screening_type, hfx_type, log_zero
  USE input_constants,                 ONLY: do_hfx_potential_truncated
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: m_walltime
  USE mathconstants
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_sync
  USE orbital_pointers
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE t_c_g0,                          ONLY: init
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC  derivatives_four_center
  TYPE(lib_int),   ALLOCATABLE, DIMENSION(:), SAVE      :: lib
  TYPE(lib_deriv), ALLOCATABLE, DIMENSION(:), SAVE      :: deriv

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_derivatives'


  INTEGER, SAVE                                         :: init_t_c_g0_lmax = -1

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
  SUBROUTINE derivatives_four_center(qs_env,ks_matrix,energy,rho,hfx_section,para_env,& 
                                     irep, adiabatic_rescale_factor, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: hfx_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: irep
    REAL(dp), INTENT(IN), OPTIONAL           :: adiabatic_rescale_factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derivatives_four_center', &
      routineP = moduleN//':'//routineN

    INTEGER :: bin, coord, forces_map(4,2), handle, i, i_atom, i_thread, &
      iatom, ikind, iset, j_atom, jatom, jkind, jset, k_atom, katom, kkind, &
      kset, l_atom, latom, lkind, lset, ma, max_am, max_set, maxsgf, mb, mc, &
      md, my_bin_size, my_istart, my_jstart, my_kstart, my_lstart, n_threads, &
      natom, nco_max, ncoa, ncob, ncoc, ncod, nseta, nsetb, nsetc, nsetd, &
      nsgfa_total, nsgfb_total, nsgfc_total, nsgfd_total, nspins, p_index, &
      pa, pa1, pa2, pb, pb1, pb2, pc, pc1, pc2, pd, pd1, pd2, sgfa, sgfb, &
      sgfc, sgfd, stat
    INTEGER(int_8)                           :: box_idx1, my_current_counter, &
                                                my_original_cpu_id, &
                                                n_processes, ncpu
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                last_sgf_global
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, npgfa, npgfb, npgfc, npgfd, nsgfa, nsgfb, &
      nsgfc, nsgfd
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfc, first_sgfd
    LOGICAL                                  :: do_periodic, failure, &
                                                first_j_loop, first_k_loop, &
                                                first_l_loop, &
                                                screen_pmat_forces
    REAL(dp) :: bintime_start, bintime_stop, cartesian_estimate, eps_schwarz, &
      fac, hf_fraction, max_contraction_val, max_val1, max_val2, &
      my_adiabatic_rescale_factor, omega2Inv, pmax_entry, ra(3), rb(3), &
      rc(3), rd(3), symm_fac, temp(3), temp1, temp1b, temp2, temp3, temp3b
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: work_forces
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      TARGET                                 :: primitive_forces
    REAL(dp), DIMENSION(:), POINTER          :: T1, T2
    REAL(dp), DIMENSION(:, :), POINTER :: full_density, full_density_beta, &
      max_contraction, sphi_a, sphi_b, sphi_c, sphi_d, zeta, zetb, zetc, zetd
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: pmax_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(hfx_basis_info_type), POINTER       :: basis_info
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_general_type)                   :: general_parameter
    TYPE(hfx_load_balance_type)              :: load_balance_parameter
    TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(hfx_screening_type)                 :: screening_parameter
    TYPE(hfx_type), POINTER                  :: actual_x_data
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force

!$  INTEGER :: omp_get_max_threads,omp_get_thread_num
    REAL(dp), DIMENSION(:, :, :, :), &
              POINTER, SAVE                  :: shm_initial_p
    INTEGER, SAVE                            :: shm_number_of_p_entries
    INTEGER, DIMENSION(:,:), POINTER, SAVE   :: shm_is_assoc_atomic_block
    INTEGER(int_8)                           :: shm_neris_total, shm_neris_onthefly, &
                                                shm_storage_counter_integrals, shm_neris_incore
    INTEGER                                  :: nkind, l_max
    TYPE(cell_type), POINTER                 :: cell
    INTEGER                                  :: unit_id
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:,:,:,:,:,:), POINTER        :: screen_coeffs_pgf, radii_pgf
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:,:,:,:), POINTER            :: screen_coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:,:), POINTER                :: screen_coeffs_kind, tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2
    REAL(dp)                                 :: ln_10, B11(3), rab2, D11(3), rcd2, log10_eps_schwarz,&
                                                log10_pmax, max_val2_set, pbc_B(3), pbc_D(3), log_2


    CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, atomic_kind_set=atomic_kind_set, &
                    cell=cell, error=error)
    !! One atom systems have no contribution to forces
    IF( SIZE(particle_set,1) == 1) THEN
      RETURN
    END IF

    CALL mp_sync(para_env%group)
    CALL timeset(routineN,handle)

    !! Calculate l_max used in fgamma , because init_md_ftable is definitely not thread safe
    nkind = SIZE(atomic_kind_set,1)
    l_max = 0
    DO ikind=1,nkind
      l_max = MAX(l_max,MAXVAL(qs_env%x_data(1,1)%basis_parameter(ikind)%lmax))
    ENDDO
    l_max = 4*l_max + 1
    CALL init_md_ftable(l_max)

    IF(qs_env%x_data(1,1)%potential_parameter%potential_type == do_hfx_potential_truncated) THEN
      IF(l_max>init_t_c_g0_lmax) THEN
        unit_id = get_unit_number()
        OPEN(unit_id,FILE=qs_env%x_data(1,1)%potential_parameter%filename)
        CALL init(l_max,unit_id)
        CLOSE(unit_id)
        init_t_c_g0_lmax = l_max
      END IF
      ! ** We use the same cutoff as for the energy
! **    CALL set_eps_cutoff(qs_env%x_data(1,1)%screening_parameter%eps_schwarz_forces)
    END IF


    n_threads = 1
!$  n_threads = omp_get_max_threads()
    ALLOCATE(deriv(n_threads), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(lib(n_threads), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    !! get force array
    CALL get_qs_env(qs_env=qs_env, force=force, error=error)

    my_adiabatic_rescale_factor = 1.0_dp
    IF(PRESENT(adiabatic_rescale_factor)) THEN
      my_adiabatic_rescale_factor = adiabatic_rescale_factor
    END IF

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(qs_env,&
!$OMP                                  ks_matrix,&
!$OMP                                  energy,&
!$OMP                                  rho,&
!$OMP                                  hfx_section,&
!$OMP                                  para_env,&
!$OMP                                  irep,&
!$OMP                                  error,&
!$OMP                                  ncoset,&
!$OMP                                  deriv,&
!$OMP                                  lib,&
!$OMP                                  n_threads,&
!$OMP                                  full_density,&
!$OMP                                  full_density_beta,&
!$OMP                                  shm_initial_p,&
!$OMP                                  pmax_set,&
!$OMP                                  shm_is_assoc_atomic_block,&
!$OMP                                  shm_number_of_p_entries,&
!$OMP                                  shm_neris_total,&
!$OMP                                  shm_neris_onthefly,&
!$OMP                                  shm_storage_counter_integrals,&
!$OMP                                  shm_neris_incore,&
!$OMP                                  force,&
!$OMP                                  my_adiabatic_rescale_factor,&
!$OMP                                  cell,&
!$OMP                                  screen_coeffs_set,&
!$OMP                                  screen_coeffs_kind,&
!$OMP                                  screen_coeffs_pgf,&
!$OMP                                  radii_pgf)

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

    IF(potential_parameter%omega /= 0.0_dp) THEN
      omega2Inv = 1.0_dp/(potential_parameter%omega**2)
    ELSE
      omega2Inv = 0.0_dp
    END IF

    load_balance_parameter = actual_x_data%load_balance_parameter
    basis_parameter        => actual_x_data%basis_parameter
   
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    error=error)
 
    maxsgf = basis_info%max_sgf
    max_set = basis_info%max_set
    max_am = basis_info%max_am 
    natom = SIZE(particle_set,1)


    hf_fraction = general_parameter%fraction
    hf_fraction = hf_fraction * my_adiabatic_rescale_factor
    eps_schwarz = screening_parameter%eps_schwarz_forces
    IF( eps_schwarz < 0.0_dp) THEN
      log10_eps_schwarz = log_zero
    ELSE
      log10_eps_schwarz = LOG10(eps_schwarz) 
    END IF
    screen_pmat_forces =  screening_parameter%do_p_screening_forces

    !! Initialize memory for libint
    CALL initialize_libderiv(deriv(i_thread+1), max_am, error)
    CALL initialize_libint(lib(i_thread+1), max_am, error)

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
      ENDDO
    ENDDO

!$OMP MASTER
    NULLIFY(full_density)
    NULLIFY(full_density_beta)
    !! Set pointer for is_assoc helper
    shm_is_assoc_atomic_block => actual_x_data%is_assoc_atomic_block 
    shm_number_of_p_entries = actual_x_data%number_of_p_entries
    !! Get the full density from all the processors
    CALL get_full_density(para_env, full_density, rho%rho_ao(1)%matrix, actual_x_data%number_of_p_entries,&
                          actual_x_data%is_assoc_atomic_block, natom, last_sgf_global, &
                          kind_of, basis_parameter, error)   
    IF(nspins == 2) THEN
      CALL get_full_density(para_env, full_density_beta, rho%rho_ao(2)%matrix, actual_x_data%number_of_p_entries,&
                           actual_x_data%is_assoc_atomic_block, natom, last_sgf_global, &
                           kind_of, basis_parameter, error)
    END IF
    CALL symmetrize_matrix(full_density,"upper_to_lower")
    IF(nspins == 2) CALL symmetrize_matrix(full_density_beta,"upper_to_lower")

    !! Calculate max entries for screening on actual density. If screen_p_mat_forces = FALSE, the
    !! matrix is initialized to 1.0
    IF ( screen_pmat_forces ) THEN
      NULLIFY(shm_initial_p)
      shm_initial_p => actual_x_data%initial_p
      CALL update_pmax_mat(shm_initial_p, full_density,full_density_beta,natom, &
                           max_set,  kind_of, basis_parameter, last_sgf_global, error)
    END IF 

    !! Calculatematrices for screening purpose taking into account periodicity. The farfield part is already
    !! in memory

    screen_coeffs_set   => actual_x_data%screen_funct_coeffs_set
    screen_coeffs_kind  => actual_x_data%screen_funct_coeffs_kind
    screen_coeffs_pgf   => actual_x_data%screen_funct_coeffs_pgf
    radii_pgf           => actual_x_data%pair_dist_radii_pgf

!$OMP END MASTER
!$OMP BARRIER

    !! Load balance the work
    IF( actual_x_data%b_first_load_balance_forces ) THEN
      CALL load_balance_forces(actual_x_data,eps_schwarz,particle_set,max_set,para_env,&
                               potential_parameter,screen_coeffs_set,screen_coeffs_kind,&
                               shm_initial_p,shm_is_assoc_atomic_block, do_periodic,&
                               load_balance_parameter, kind_of, basis_parameter,i_thread, n_threads,&
                               cell, screen_pmat_forces,error)
      actual_x_data%b_first_load_balance_forces = .FALSE.
    ELSE
      CALL update_load_balance_forces(actual_x_data, para_env, &
                                      load_balance_parameter, &
                                      i_thread, n_threads, error)
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

    !! Allocate work arrays
    ALLOCATE(primitive_forces(12*nco_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    primitive_forces=0.0_dp
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
    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads
!$OMP BARRIER
    my_bin_size = SIZE(actual_x_data%distribution_forces)
    DO bin=1,my_bin_size
       actual_x_data%distribution_forces(bin)%time_forces=0.0_dp
    ENDDO
    !! outer loop will be cycled if actual work_load done
    DO bin=1,my_bin_size
      bintime_start=m_walltime()
      !! Get information on workload
      first_j_loop = .TRUE.
      first_l_loop = .TRUE.
      first_k_loop = .TRUE.
      my_istart = actual_x_data%distribution_forces(bin)%istart
      my_jstart = actual_x_data%distribution_forces(bin)%jstart
      my_kstart = actual_x_data%distribution_forces(bin)%kstart
      my_lstart = actual_x_data%distribution_forces(bin)%lstart
      my_original_cpu_id = actual_x_data%distribution_forces(bin)%original_cpu_id 
      my_current_counter = 0
      !! It could happen, that one bin is empty. In that case just skip
      IF(my_istart == 0 .OR. my_jstart == 0 .OR. my_kstart == 0 .OR. my_lstart == 0 .OR. &
         actual_x_data%distribution_forces(bin)%number_of_atom_quartets == 0) THEN
         my_istart = natom + 1 ; my_jstart = natom + 1 ;my_kstart = natom + 1 ;my_lstart = natom + 1 
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
        i_atom = atom_of_kind(iatom)
        pa1 = last_sgf_global(iatom) - nsgfa_total
        forces_map(1,1) = ikind
        forces_map(1,2) = i_atom
        IF(first_j_loop) jatom=my_jstart-1
        IF(.NOT. first_j_loop) jatom = iatom-1
        DO WHILE(jatom+1<=natom)
          jatom = jatom + 1
          !! Be consistent to the modulo distribution of the load balance
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
          ! ** pbc for rb
          IF( do_periodic ) THEN
            temp = rb -ra
            pbc_B = pbc(temp,cell)
            B11 = ra + pbc_B
            rab2 = (ra(1)-B11(1))**2+(ra(2)-B11(2))**2+(ra(3)-B11(3))**2
          ELSE
            rab2 = (ra(1)-rb(1))**2+(ra(2)-rb(2))**2+(ra(3)-rb(3))**2
            B11 = ra - rb
          END IF
          
          j_atom = atom_of_kind(jatom)
          pb1 = last_sgf_global(jatom) - nsgfb_total
          forces_map(2,1) = jkind
          forces_map(2,2) = j_atom
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
            k_atom = atom_of_kind(katom)
            pc1 = last_sgf_global(katom) - nsgfc_total
            forces_map(3,1) = kkind
            forces_map(3,2) = k_atom
            IF(first_l_loop) THEN
              latom=my_lstart-1
            ELSE
              latom = katom-1
            END IF
            DO WHILE(latom+1<=natom)
              latom = latom + 1
              lkind = kind_of(latom)
              rd=particle_set(latom)%r(:)
             
              ! ** pbc for rd
              IF( do_periodic ) THEN
                temp = rd -rc
                pbc_D = pbc(temp,cell)
                D11 = rc + pbc_D
                rcd2 = (rc(1)-D11(1))**2+(rc(2)-D11(2))**2+(rc(3)-D11(3))**2
              ELSE
                rcd2 = (rc(1)-rd(1))**2+(rc(2)-rd(2))**2+(rc(3)-rd(3))**2
                D11 = rc - rd
              END IF

              !All four centers equivalent => zero-contribution
              IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE 
              IF(katom+latom<=iatom+jatom)  THEN
                IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE

                my_current_counter = my_current_counter + 1
                IF(my_current_counter <= actual_x_data%distribution_forces(bin)%number_of_atom_quartets) THEN
                  IF( (screen_coeffs_kind(jkind,ikind)%x(1)*rab2+&
                       screen_coeffs_kind(jkind,ikind)%x(2) ) + &
                      (screen_coeffs_kind(lkind,kkind)%x(1)*rcd2+&
                       screen_coeffs_kind(lkind,kkind)%x(2)) < log10_eps_schwarz ) CYCLE
 
                  IF(.NOT. (shm_is_assoc_atomic_block(latom,iatom)>=1 .AND. &
                            shm_is_assoc_atomic_block(katom,iatom)>=1 .AND. &
                            shm_is_assoc_atomic_block(katom,jatom)>=1 .AND. &
                            shm_is_assoc_atomic_block(latom,jatom)>=1 ) ) CYCLE
  
                  ld_max => basis_parameter(lkind)%lmax
                  ld_min => basis_parameter(lkind)%lmin
                  npgfd => basis_parameter(lkind)%npgf
                  nsetd = basis_parameter(lkind)%nset
                  zetd => basis_parameter(lkind)%zet
                  nsgfd => basis_parameter(lkind)%nsgf
                  first_sgfd => basis_parameter(lkind)%first_sgf
                  sphi_d => basis_parameter(lkind)%sphi
                  nsgfd_total = basis_parameter(lkind)%nsgf_total
  
 
                  l_atom = atom_of_kind(latom)
                  lkind = kind_of(latom)
                  pd1 = last_sgf_global(latom) - nsgfd_total
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
                  DO iset = 1, nseta
                    ncoa = npgfa(iset)*ncoset(la_max(iset))
                    sgfa = first_sgfa(1,iset)
                    DO jset = 1,nsetb
                      ncob = npgfb(jset)*ncoset(lb_max(jset))
                      sgfb = first_sgfb(1,jset)
                      max_val1 = screen_coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + &
                                 screen_coeffs_set(jset,iset,jkind,ikind)%x(2)
                      IF (max_val1 + (screen_coeffs_kind(lkind,kkind)%x(1)*rcd2 + &
                                      screen_coeffs_kind(lkind,kkind)%x(2) ) < log10_eps_schwarz) CYCLE

                      DO kset = 1,nsetc
                        ncoc = npgfc(kset)*ncoset(lc_max(kset))
                        sgfc = first_sgfc(1,kset)
                        DO lset = 1,nsetd

                          ncod = npgfd(lset)*ncoset(ld_max(lset)) 
                          sgfd = first_sgfd(1,lset)
                          max_val2_set = (screen_coeffs_set(lset,kset,lkind,kkind)%x(1)*rcd2 + &
                                          screen_coeffs_set(lset,kset,lkind,kkind)%x(2) )
                          max_val2 = max_val1 + max_val2_set
                          IF(max_val2<log10_eps_schwarz) CYCLE
                          IF( screen_pmat_forces ) THEN 
                            log10_pmax = log_2 + MAX(shm_initial_p(kset,iset,katom,iatom) + shm_initial_p(lset,jset,latom,jatom),&
                                                     shm_initial_p(lset,iset,latom,iatom) + shm_initial_p(kset,jset,katom,jatom))
                          ELSE
                            log10_pmax = 0.0_dp
                          END IF

                          max_val2 = max_val2 + log10_pmax
                          IF(max_val2<log10_eps_schwarz)  CYCLE

                          pmax_entry = EXP(log10_pmax*ln_10) 
                          max_contraction_val =  max_contraction(iset,iatom) * &
                                                   max_contraction(jset,jatom) * &
                                                   max_contraction(kset,katom) * &
                                                   max_contraction(lset,latom) * &
                                                   pmax_entry
                          tmp_R_1 => radii_pgf(:,:,jset,iset,jkind,ikind)
                          tmp_R_2 => radii_pgf(:,:,lset,kset,lkind,kkind)
                          tmp_screen_pgf1 => screen_coeffs_pgf(:,:,jset,iset,jkind,ikind)
                          tmp_screen_pgf2 => screen_coeffs_pgf(:,:,lset,kset,lkind,kkind)

                          IF( .NOT. do_periodic ) THEN
                            CALL forces4(ra, rb, rc, rd,&
                                         zeta(:,iset), zetb(:,jset), zetc(:,kset), zetd(:,lset),&
                                         la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                         lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                         npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                         ncoa, ncob, ncoc, ncod,&
                                         primitive_forces,&
                                         work_forces,&
                                         potential_parameter, eps_schwarz, max_contraction_val, cartesian_estimate,&
                                         i_thread, rab2, rcd2, max_val2_set, log10_eps_schwarz, log10_pmax,&
                                         tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2, B11, D11)
                          ELSE
                            CALL forces4_periodic(ra, rb, rc, rd,&
                                                  zeta(:,iset), zetb(:,jset), zetc(:,kset), zetd(:,lset),&
                                                  la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                                  lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                                  npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                                  ncoa, ncob, ncoc, ncod,&
                                                  primitive_forces,&
                                                  work_forces,&
                                                  potential_parameter, &
                                                  actual_x_data%neighbor_cells, screen_coeffs_set(jset,iset,jkind,ikind)%x,& 
                                                  screen_coeffs_set(lset,kset,lkind,kkind)%x, eps_schwarz,&
                                                  max_contraction_val, cartesian_estimate, i_thread, cell, &
                                                  log10_pmax, B11, D11, log10_eps_schwarz, &
                                                  tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2,error)
                          END IF 
  
                          IF(cartesian_estimate<eps_schwarz) CYCLE
                          DO coord = 1,12
                            T2 => primitive_forces((coord-1)*ncoa*ncob*ncoc*ncod+1:coord*ncoa*ncob*ncoc*ncod)
                            CALL dgemm("T","N",ncob*ncoc*ncod,nsgfa(iset),ncoa,&
                                       1.0_dp, T2(1), ncoa,&
                                       sphi_a(1,sgfa), SIZE(sphi_a,1),&
                                       0.0_dp, T1(1),ncob*ncoc*ncod)
  
                            CALL dgemm("T","N",nsgfa(iset)*ncoc*ncod,nsgfb(jset),ncob,&
                                       1.0_dp, T1(1),ncob,&
                                       sphi_b(1,sgfb), SIZE(sphi_b,1),&
                                       0.0_dp, T2(1),&
                                       nsgfa(iset)*ncoc*ncod)
    
                            CALL dgemm("T","N",nsgfa(iset)*nsgfb(jset)*ncod,nsgfc(kset),ncoc,&
                                       1.0_dp, T2(1),&
                                       ncoc,&
                                       sphi_c(1,sgfc), SIZE(sphi_c,1),&
                                       0.0_dp, T1(1),nsgfa(iset)*nsgfb(jset)*ncod)
    
                            CALL dgemm("T","N",nsgfa(iset)*nsgfb(jset)*nsgfc(kset),nsgfd(lset),ncod,&
                                       1.0_dp, T1(1),ncod,&
                                       sphi_d(1,sgfd), SIZE(sphi_d,1),&
                                       0.0_dp, T2(1),&
                                       nsgfa(iset)*nsgfb(jset)*nsgfc(kset))
                            p_index = 0
                            DO md = 1,nsgfd(lset)
                              pd2 = md+sgfd-1
                              pd = pd1 + pd2
                              DO mc = 1, nsgfc(kset)
                                pc2 = mc+sgfc-1
                                pc = pc1 + pc2
                                DO mb = 1, nsgfb(jset)
                                  pb2 = mb+sgfb-1
                                  pb = pb1 + pb2
                                  temp1 = full_density(pc,pb) * fac
                                  IF( nspins == 2 ) temp1b = full_density_beta(pc,pb) * fac
                                  DO ma = 1,nsgfa(iset) 
                                    pa2 = ma+sgfa-1
                                    pa = pa1 + pa2
                                    temp3 = full_density(pb,pd) * fac
                                    IF( nspins == 2 ) temp3b = full_density_beta(pb,pd) * fac
                                    p_index = p_index + 1 
                                    temp2 = temp1 * full_density(pa,pd) + temp3 * full_density(pa,pc)

                                    !$OMP ATOMIC
                                    force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1, &
                                          forces_map((coord-1)/3 + 1,2)) = &
                                    force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1, &
                                          forces_map((coord-1)/3 + 1,2)) - &
                                    temp2 * T2(p_index)
                                    IF(nspins == 2) THEN
                                      temp2 = temp1b*full_density_beta(pa,pd) + temp3b*full_density_beta(pa,pc)
                                      !$OMP ATOMIC
                                      force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1,&
                                   forces_map((coord-1)/3 + 1,2)) = &
                                      force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1,&
                                            forces_map((coord-1)/3 + 1,2)) - &
                                      temp2 * T2(p_index)
                                    END IF 
                                  END DO !ma
                                END DO !mb
                              END DO !mc
                            END DO !md
                          END DO !coord
                        END DO !lset
                      END DO !ksetm
                    END DO !jset
                  END DO !iset
                END IF
              END IF !cycle
            END DO !latom
            first_l_loop = .FALSE.
          END DO !katom
          first_k_loop = .FALSE.
        END DO !jatom
        first_j_loop = .FALSE.
      END DO !iatom
      bintime_stop=m_walltime()
      actual_x_data%distribution_forces(bin)%time_forces=bintime_stop-bintime_start
    END DO  !bin
 
!$OMP BARRIER
    DEALLOCATE(last_sgf_global, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP MASTER
    DEALLOCATE(full_density,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(nspins==2) THEN
      DEALLOCATE(full_density_beta,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

!$OMP END MASTER
    DEALLOCATE(T1,primitive_forces, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(work_forces, atom_of_kind,kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(max_contraction,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
 
    !
    ! this wraps to free_libint, but currently causes a segfault
    ! as a result, we don't call it, but some memory remains allocated
    !  
    ! CALL terminate_libderiv(deriv)

!$OMP END PARALLEL
    DEALLOCATE(lib, deriv, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL mp_sync(para_env%group)
    CALL timestop(handle)
  END SUBROUTINE derivatives_four_center


! *****************************************************************************
!> \brief calculates two-electron derivatived of a quartet/shell using the library 
!>      lib_deriv
!> \param r position
!> \param z zeta
!> \param l angular momentum
!> \param npgf number of primitive cartesian gaussian in actual shell
!> \param nco number of cartesian gaussians
!> \param primitive_forces array of primitive_forces
!> \param work_forces work_storage
!> \param potential_parameter Infos about potential
!> \param eps_schwarz threshold
!> \param max_contraction_val max multiplication factor cart -> sph
!> \param cartesian estimate max calculated integral
!> \param i_thread current thread ID
!> \param rab2, rcd2 distance (ab| and (cd|
!> \param max_val2_set schwarz estimate for product
!> \param log10_eps_schwarz logarithm of threshold
!> \param log10_pmax logarithm of max p value
!> \param R1_pgf, R2_pgf Radii functions of pair distributions
!> \param pgf1, pgf2 near field screening functions
!> \param B11, D11 ra-rb and rc-rd
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2007 created [Manuel Guidon]
!>      02.2009 completely rewritte screening part [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE forces4(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                    la_min, la_max, lb_min, lb_max,&
                    lc_min, lc_max, ld_min, ld_max,&
                    npgfa, npgfb, npgfc, npgfd,&
                    ncoa, ncob, ncoc, ncod,&
                    primitive_forces, work_forces,&
                    potential_parameter, eps_schwarz, max_contraction_val, &
                    cartesian_estimate, i_thread,rab2, rcd2, max_val2_set,&
                    log10_eps_schwarz, log10_pmax,R1_pgf, R2_pgf, pgf1, pgf2, B11, D11)

  
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb, zetc, zetd
    INTEGER, INTENT(IN) :: la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
      ld_min, ld_max, npgfa, npgfb, npgfc, npgfd, ncoa, ncob, ncoc, ncod
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod, 12)  :: primitive_forces
    REAL(dp), DIMENSION(*)                   :: work_forces
    TYPE(hfx_potential_type)                 :: potential_parameter
    REAL(dp), INTENT(IN)                     :: eps_schwarz, &
                                                max_contraction_val
    REAL(dp), INTENT(INOUT)                  :: cartesian_estimate
    INTEGER, INTENT(IN)                      :: i_thread
    REAL(dp)                                 :: rab2, rcd2, max_val2_set, &
                                                log10_eps_schwarz, log10_pmax
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: R1_pgf, R2_pgf, pgf1, pgf2
    REAL(dp), INTENT(IN)                     :: B11(3), D11(3)

    CHARACTER(len=*), PARAMETER :: routineN = 'forces4', &
      routineP = moduleN//':'//routineN

    INTEGER :: ipgf, jpgf, kpgf, la, lb, lc, ld, lpgf, offset_a, offset_a1, &
      offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1
    LOGICAL                                  :: failure
    REAL(dp) :: Eta, EtaInv, P(3), pgf_max_1, pgf_max_2, Q(3), R1, R2, Rho, &
      RhoInv, rpq2, S1234, S1234a, tmp_max, W(3), Zeta1, Zeta_A, Zeta_B, &
      Zeta_C, Zeta_D, ZetaInv, ZetapEtaInv

    failure = .FALSE.
    cartesian_estimate = 0.0_dp

    primitive_forces = 0.0_dp
    DO ipgf = 1,npgfa
      offset_a1 = (ipgf-1)*ncoset(la_max)
      Zeta_A = zeta(ipgf)
      DO jpgf = 1,npgfb
        pgf_max_1 = pgf1(jpgf,ipgf)%x(1)*rab2+pgf1(jpgf,ipgf)%x(2)
        IF( pgf_max_1 + max_val2_set + log10_pmax < log10_eps_schwarz) CYCLE
        R1 = R1_pgf(jpgf,ipgf)%x(1)*rab2 + &
             R1_pgf(jpgf,ipgf)%x(2)
        offset_b1 = (jpgf-1)*ncoset(lb_max)
        Zeta_B = zetb(jpgf)
        Zeta1 = Zeta_A + Zeta_B
        ZetaInv = 1.0_dp/Zeta1
        S1234a = (-Zeta_A*Zeta_B*ZetaInv*rab2)
        P = (Zeta_A*ra + Zeta_B*rb)*ZetaInv
        DO kpgf = 1,npgfc
          offset_c1 = (kpgf-1)*ncoset(lc_max)
          Zeta_C = zetc(kpgf)
          DO lpgf = 1,npgfd
            pgf_max_2 = pgf2(lpgf,kpgf)%x(1)*rcd2+pgf2(lpgf,kpgf)%x(2)
            IF (pgf_max_1 + pgf_max_2 + log10_pmax < log10_eps_schwarz ) CYCLE
            R2 = R2_pgf(lpgf,kpgf)%x(1)*rcd2 + &
                 R2_pgf(lpgf,kpgf)%x(2)
            offset_d1 = (lpgf-1)*ncoset(ld_max)
            Zeta_D = zetd(lpgf)
            Eta  =  Zeta_C + Zeta_D
            EtaInv = 1.0_dp/Eta
            ZetapEtaInv = Zeta1+Eta
            ZetapEtaInv = 1.0_dp/ZetapEtaInv
            Rho = Zeta1*Eta*ZetapEtaInv
            RhoInv = 1.0_dp/Rho
            S1234 = EXP(S1234a-Zeta_C*Zeta_D*EtaInv*rcd2)
            Q = (Zeta_C*rc + Zeta_D*rd)*EtaInv
            rpq2 = (P(1)-Q(1))**2+(P(2)-Q(2))**2+(P(3)-Q(3))**2
            W = (Zeta1*P+Eta*Q)*ZetapEtaInv
            DO la = la_min,la_max
              offset_a = offset_a1 + ncoset(la-1)
              DO lb = lb_min, lb_max
                offset_b = offset_b1 + ncoset(lb-1)
                DO lc = lc_min, lc_max
                  offset_c = offset_c1 + ncoset(lc-1)
                  DO ld = ld_min, ld_max
                    offset_d = offset_d1 + ncoset(ld-1) 
                    !Build primitives
                    tmp_max = 0.0_dp
                    CALL evaluate_deriv_eri(deriv(i_thread+1), ra, rb, rc, rd,&
                                            zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                            la, lb, lc ,ld,&
                                            work_forces,&
                                            ncoa, ncob, ncoc, ncod, &
                                            primitive_forces, &
                                            offset_a, offset_b, offset_c, offset_d,&
                                            potential_parameter, max_contraction_val, tmp_max, eps_schwarz,&
                                            Zeta1,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                                            P,Q,W,rpq2,B11,D11,&
                                            R1,R2)
                    cartesian_estimate = MAX(cartesian_estimate, tmp_max)
                  END DO !ld
                END DO !lc
              END DO !lb
            END DO !la
          END DO !lpgf
        END DO !kpgf
      END DO !jpgf
    END DO !ipgf
  
  END SUBROUTINE forces4
! *****************************************************************************
!> \brief calculates two-electron derivatived of a quartet/shell using the library 
!>      lib_deriv in the periodic case
!> \param r position
!> \param z zeta
!> \param l angular momentum
!> \param npgf number of primitive cartesian gaussian in actual shell
!> \param nco number of cartesian gaussians
!> \param primitive_forces array of primitive_forces
!> \param work_forces work_storage
!> \param potential_parameter Infos about potential
!> \param neighbor_cells periodic image cells
!> \param screen1, screen2 set based near field screening functions
!> \param eps_schwarz threshold
!> \param max_contraction_val max multiplication factor cart -> sph
!> \param cartesian estimate max calculated integral
!> \param i_thread current thread ID
!> \param cell cell
!> \param log10_pmax logarithm of max p value
!> \param B11, D11 ra-rb and rc-rd after first MIC
!> \param log10_eps_schwarz logarithm of threshold
!> \param R1_pgf, R2_pgf Radii functions of pair distributions
!> \param pgf1, pgf2 near field screening functions
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2007 created [Manuel Guidon]
!>      02.2009 completely rewritte screening part [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE forces4_periodic(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                              la_min, la_max, lb_min, lb_max,&
                              lc_min, lc_max, ld_min, ld_max,&
                              npgfa, npgfb, npgfc, npgfd,&
                              ncoa, ncob, ncoc, ncod,&
                              primitive_forces, work_forces,&
                              potential_parameter, neighbor_cells,&
                              screen1, screen2, eps_schwarz, max_contraction_val,&
                              cart_estimate, i_thread, cell, log10_pmax,&
                              B11,D11,log10_eps_schwarz,R1_pgf, R2_pgf, pgf1, pgf2,error)
  
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb, zetc, zetd
    INTEGER, INTENT(IN) :: la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
      ld_min, ld_max, npgfa, npgfb, npgfc, npgfd, ncoa, ncob, ncoc, ncod
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod, 12)  :: primitive_forces
    REAL(dp), DIMENSION(*)                   :: work_forces
    TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(hfx_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(dp), INTENT(IN)                     :: screen1(2), screen2(2), &
                                                eps_schwarz, &
                                                max_contraction_val
    REAL(dp), INTENT(INOUT)                  :: cart_estimate
    INTEGER, INTENT(IN)                      :: i_thread
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), INTENT(IN)                     :: log10_pmax, B11(3), D11(3), &
                                                log10_eps_schwarz
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: R1_pgf, R2_pgf, pgf1, pgf2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'forces4_periodic', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, ipgf, j, jpgf, k, kpgf, la, lb, lc, ld, lpgf, offset_a, &
      offset_a1, offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1
    LOGICAL                                  :: failure
    REAL(dp) :: AB(3), C11(3), CD(3), Eta, EtaInv, im_B(3), im_D(3), &
      max_val1, max_val2, P(3), pgf_max_1, pgf_max_2, PQ(3), Q(3), R1, R2, &
      rab2, rcd2, Rho, RhoInv, rpq2, S1234, S1234a, shift(3), temp(3), &
      temp_CC(3), temp_DD(3), tmp_D(3), tmp_max, W(3), Zeta1, Zeta_A, Zeta_B, &
      Zeta_C, Zeta_D, ZetaInv, ZetapEtaInv

    failure = .FALSE.
    cart_estimate = 0.0_dp

    primitive_forces = 0.0_dp

    DO i = 1,SIZE(neighbor_cells)
      im_B = B11 + neighbor_cells(i)%cell_r(:)
      AB = ra-im_B
      rab2 = AB(1)**2+AB(2)**2+AB(3)**2
      max_val1 = screen1(1)*rab2+screen1(2)
      IF( max_val1 + screen2(2) + log10_pmax < log10_eps_schwarz ) CYCLE
      DO j = 1,SIZE(neighbor_cells)
        im_D = D11 + neighbor_cells(j)%cell_r(:)
        rcd2 = (im_D(1)-rc(1))**2+(im_D(2)-rc(2))**2+(im_D(3)-rc(3))**2
        max_val2 = screen2(1)*rcd2+screen2(2)
        IF( max_val1 + max_val2 + log10_pmax < log10_eps_schwarz ) CYCLE
        DO ipgf = 1,npgfa
          Zeta_A = zeta(ipgf)
          offset_a1 = (ipgf-1)*ncoset(la_max)
          DO jpgf = 1,npgfb
            pgf_max_1 = pgf1(jpgf,ipgf)%x(1)*rab2+pgf1(jpgf,ipgf)%x(2)
            IF( pgf_max_1 + max_val2 + log10_pmax < log10_eps_schwarz) CYCLE
            R1 = R1_pgf(jpgf,ipgf)%x(1)*rab2 + &
                 R1_pgf(jpgf,ipgf)%x(2)
            Zeta_B = zetb(jpgf)
            offset_b1 = (jpgf-1)*ncoset(lb_max)
            Zeta1 = Zeta_A + Zeta_B
            ZetaInv = 1.0_dp/Zeta1
            S1234a = (-Zeta_A*Zeta_B*ZetaInv*rab2)
            P = (Zeta_A*ra + Zeta_B*im_B)*ZetaInv
            DO kpgf = 1,npgfc
              Zeta_C = zetc(kpgf)
              offset_c1 = (kpgf-1)*ncoset(lc_max)
              DO lpgf = 1,npgfd
                pgf_max_2 = pgf2(lpgf,kpgf)%x(1)*rcd2+pgf2(lpgf,kpgf)%x(2)
                IF (pgf_max_1 + pgf_max_2 + log10_pmax < log10_eps_schwarz ) CYCLE
                R2 = R2_pgf(lpgf,kpgf)%x(1)*rcd2 +&
                     R2_pgf(lpgf,kpgf)%x(2)
                offset_d1 = (lpgf-1)*ncoset(ld_max)
                Zeta_D = zetd(lpgf)
                Eta  =  Zeta_C + Zeta_D
                EtaInv = 1.0_dp/Eta
                ZetapEtaInv = Zeta1+Eta
                ZetapEtaInv = 1.0_dp/ZetapEtaInv
                Rho = Zeta1*Eta*ZetapEtaInv
                RhoInv = 1.0_dp/Rho
                S1234 = EXP(S1234a-Zeta_C*Zeta_D*EtaInv*rcd2)
                Q = (Zeta_C*rc + Zeta_D*im_D)*EtaInv
                temp = P-Q
                PQ = pbc(temp,cell)
                shift = -PQ + temp
                temp_CC = rc+shift
                temp_DD = im_D+shift
                DO k = 1,SIZE(neighbor_cells)
                  C11 = temp_CC + neighbor_cells(k)%cell_r(:)
                  tmp_D = temp_DD + neighbor_cells(k)%cell_r(:)
                  Q = (Zeta_C*C11 + Zeta_D*tmp_D)*EtaInv
                  W = (Zeta1*P+Eta*Q)*ZetapEtaInv
                  CD = C11-tmp_D
                  rpq2 = (P(1)-Q(1))**2+(P(2)-Q(2))**2+(P(3)-Q(3))**2
                  DO la = la_min,la_max
                    offset_a = offset_a1 + ncoset(la-1)
                    DO lb = lb_min, lb_max
                      offset_b = offset_b1 + ncoset(lb-1)
                      DO lc = lc_min, lc_max
                        offset_c = offset_c1 + ncoset(lc-1)
                        DO ld = ld_min, ld_max
                          offset_d = offset_d1 + ncoset(ld-1)
                          !Build primitives
                          tmp_max = 0.0_dp
                          CALL evaluate_deriv_eri(deriv(i_thread+1), ra, im_B, C11, tmp_D,&
                                                  zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                                  la, lb, lc ,ld,&
                                                  work_forces,&
                                                  ncoa, ncob, ncoc, ncod, &
                                                  primitive_forces, &
                                                  offset_a, offset_b, offset_c, offset_d,&
                                                  potential_parameter, max_contraction_val, tmp_max, eps_schwarz,&
                                                  Zeta1,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                                                  P,Q,W,rpq2,AB,CD,&
                                                  R1,R2)
                          cart_estimate = MAX(tmp_max,cart_estimate)
                        END DO !ld
                      END DO !lc
                    END DO !lb
                  END DO !la
                END DO !k
              END DO !lpgf
            END DO !kpgf
          END DO !jpgf
        END DO !ipgf
      END DO !j
    END DO !i
  END SUBROUTINE forces4_periodic

END MODULE hfx_derivatives
