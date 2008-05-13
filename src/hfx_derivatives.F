!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
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
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE erf_fn,                          ONLY: erfc
  USE gamma,                           ONLY: init_md_ftable
  USE hfx_communication,               ONLY: get_full_density
  USE hfx_helpers,                     ONLY: FNV64
  USE hfx_libint_interface,            ONLY: evaluate_deriv_eri,&
                                             pbc
  USE hfx_libint_wrapper,              ONLY: initialize_libderiv,&
                                             initialize_libint,&
                                             lib_deriv,&
                                             lib_int
  USE hfx_load_balance_methods,        ONLY: load_balance_forces,&
                                             update_load_balance_forces
  USE hfx_screening_methods,           ONLY: comp_pmax_mat,&
                                             comp_screen_mat,&
                                             comp_screen_mat_per
  USE hfx_types,                       ONLY: &
       hfx_basis_info_type, hfx_basis_type, hfx_box_scheme, hfx_cell_type, &
       hfx_general_type, hfx_load_balance_type, hfx_potential_type, &
       hfx_screening_type, hfx_type
  USE input_constants,                 ONLY: do_hfx_potential_coulomb,&
                                             do_hfx_potential_gaussian,&
                                             do_hfx_potential_id,&
                                             do_hfx_potential_mix_cl,&
                                             do_hfx_potential_mix_lg,&
                                             do_hfx_potential_short
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
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC  derivatives_four_center
  TYPE(lib_int),   ALLOCATABLE, DIMENSION(:), SAVE      :: lib
  TYPE(lib_deriv), ALLOCATABLE, DIMENSION(:), SAVE      :: deriv

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_derivatives'

!***
  
  CONTAINS


! *****************************************************************************
!> \brief calculates two-electron derivatived of a quartet/shell using the library 
!>      lib_deriv
!> \param r position
!> \param z zeta
!> \param l angular momentum
!> \param npgf number of primitive cartesian gaussian in actual shell
!> \param primitive_forces array of primitive_forces
!> \param work_forces work_storage
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE forces4(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                    la_min, la_max, lb_min, lb_max,&
                    lc_min, lc_max, ld_min, ld_max,&
                    npgfa, npgfb, npgfc, npgfd,&
                    ncoa, ncob, ncoc, ncod,&
                    primitive_forces, work_forces,&
                    potential_parameter, eps_schwarz, max_contraction_val, &
                    cartesian_estimate, i_thread)
  
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

    CHARACTER(len=*), PARAMETER :: routineN = 'forces4', &
      routineP = moduleN//':'//routineN

    INTEGER :: ipgf, jpgf, kpgf, la, lb, lc, ld, lpgf, offset_a, offset_a1, &
      offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1
    LOGICAL                                  :: failure
    REAL(dp)                                 :: tmp_max_all

    failure = .FALSE.
    cartesian_estimate = 0.0_dp

    primitive_forces = 0.0_dp
    DO ipgf = 1,npgfa
      offset_a1 = (ipgf-1)*ncoset(la_max)
      DO jpgf = 1,npgfb
        offset_b1 = (jpgf-1)*ncoset(lb_max)
        DO kpgf = 1,npgfc
          offset_c1 = (kpgf-1)*ncoset(lc_max)
          DO lpgf = 1,npgfd
            offset_d1 = (lpgf-1)*ncoset(ld_max)
            DO la = la_min,la_max
              offset_a = offset_a1 + ncoset(la-1)
              DO lb = lb_min, lb_max
                offset_b = offset_b1 + ncoset(lb-1)
                DO lc = lc_min, lc_max
                  offset_c = offset_c1 + ncoset(lc-1)
                  DO ld = ld_min, ld_max
                    offset_d = offset_d1 + ncoset(ld-1) 
                    !Build primitives
                    tmp_max_all = 0.0_dp
                    CALL evaluate_deriv_eri(deriv(i_thread+1), ra, rb, rc, rd,&
                                           zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                           la, lb, lc ,ld,&
                                           work_forces,&
                                           ncoa, ncob, ncoc, ncod, &
                                           primitive_forces, &
                                           offset_a, offset_b, offset_c, offset_d,&
                                           potential_parameter, max_contraction_val, tmp_max_all, eps_schwarz)
                    cartesian_estimate = MAX(cartesian_estimate, tmp_max_all)
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
!> \brief computes four center derivatives for a full basis set and updates the
!>      forces%fock_4c arrays. Uses all 8 eri symmetries
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
      pa, pa1, pa2, pb, pb1, pb2, pc, pc1, pc2, pd, pd1, pd2, pgf_radius_box, &
      radius1, radius2, sgfa, sgfb, sgfc, sgfd, stat
    INTEGER(int_8) :: box_idx1, box_idx2, box_idx3, box_idx4, &
      my_current_counter, my_original_cpu_id, n_processes, ncpu
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
    REAL(dp) :: bintime_start, bintime_stop, box_size(3), box_size_inv(3), &
      cartesian_estimate, center1(3), center2(3), center2_folded(3), eps_box, &
      eps_schwarz, fac, hf_fraction, Kbox, max_contraction_val, max_val1, &
      max_val2, my_adiabatic_rescale_factor, omega2Inv, omega_box, &
      pbc_center(3), pmax_entry, pmax_max, R1, R2, ra(3), rb(3), Rbox, rc(3), &
      rd(3), screen_max, symm_fac, temp(3), temp1, temp1b, temp2, temp3, &
      temp3b
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: work_forces
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      TARGET                                 :: primitive_forces
    REAL(dp), DIMENSION(:), POINTER          :: T1, T2
    REAL(dp), DIMENSION(:, :), POINTER :: full_density, full_density_beta, &
      max_contraction, pmax_atom, screen_atom, sphi_a, sphi_b, sphi_c, &
      sphi_d, zeta, zetb, zetc, zetd
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: pmax_set, screen_set
    REAL(dp), DIMENSION(:, :, :, :, :), &
      POINTER                                :: screen_periodic
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
    TYPE(hfx_box_scheme), DIMENSION(:,:), &
                           POINTER, SAVE     :: shm_far_field_boxes
    INTEGER, SAVE                            :: shm_number_of_p_entries
    INTEGER, DIMENSION(:,:), POINTER, SAVE   :: shm_is_assoc_atomic_block
    INTEGER(int_8)                           :: shm_neris_total, shm_neris_onthefly, &
                                                shm_storage_counter_integrals, shm_neris_incore
    INTEGER                                  :: nkind, l_max


    CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, atomic_kind_set=atomic_kind_set, error=error)
    !! One atom systems have no contribution to forces
    IF( SIZE(particle_set,1) == 1) THEN
      RETURN
    END IF

    CALL mp_sync(para_env%group)
    CALL timeset(routineN,"I",'',handle)

    !! Calculate l_max used in fgamma , because init_md_ftable is definitely not thread safe
    nkind = SIZE(atomic_kind_set,1)
    l_max = 0
    DO ikind=1,nkind
      l_max = MAX(l_max,MAXVAL(qs_env%x_data(1,1)%basis_parameter(ikind)%lmax))
    ENDDO
    l_max = 4*l_max + 1
    CALL init_md_ftable(l_max)


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
!$OMP                                  shm_neris_incore,&
!$OMP                                  box_size,&
!$OMP                                  box_size_inv,&
!$OMP                                  force,&
!$OMP                                  my_adiabatic_rescale_factor)

    i_thread = 0
!$  i_thread = omp_get_thread_num()


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
    eps_box = screening_parameter%eps_box_forces
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
    CALL comp_pmax_mat(pmax_set,pmax_atom,pmax_max, &
                       full_density,full_density_beta,screen_pmat_forces,particle_set,max_set,&
                       kind_of, basis_parameter, error)
    !! Calculatematrices for screening purpose taking into account periodicity. The farfield part is already
    !! in memory
    IF( do_periodic ) THEN
      box_size(1) = qs_env%cell%hmat(1,1)
      box_size(2) = qs_env%cell%hmat(2,2)
      box_size(3) = qs_env%cell%hmat(3,3)
      box_size_inv(1) = 1.0_dp/box_size(1)
      box_size_inv(2) = 1.0_dp/box_size(2)
      box_size_inv(3) = 1.0_dp/box_size(3)

      ALLOCATE(screen_periodic(SIZE(actual_x_data%neighbor_cells),max_set,max_set,natom,natom),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL comp_screen_mat_per(lib(i_thread+1), screen_set,screen_atom,screen_max,screen_periodic,particle_set,max_set,&
                               potential_parameter,.FALSE.,&
                               kind_of, basis_parameter,&
                               neighbor_cells=actual_x_data%neighbor_cells,&
                               box_size=box_size, box_size_inv=box_size_inv, &
                               i_thread=i_thread, error=error)

      !! Farfield screening matrix already there
      shm_far_field_boxes => actual_x_data%far_field_boxes
    ELSE

      CALL comp_screen_mat(lib(i_thread+1),screen_set,screen_atom,screen_max,particle_set,max_set,&
                           potential_parameter,.FALSE., kind_of, basis_parameter, i_thread, &
                           error=error)
      !! Farfield screening matrix already there
    
      shm_far_field_boxes => actual_x_data%far_field_boxes
    END IF

!$OMP END MASTER
!$OMP BARRIER

    !! Load balance the work
    IF( actual_x_data%b_first_load_balance_forces ) THEN
      CALL load_balance_forces(actual_x_data,eps_schwarz,eps_box,particle_set,max_set,para_env,&
                               potential_parameter,screen_set,screen_atom,screen_max,&
                               pmax_set,shm_is_assoc_atomic_block, do_periodic,&
                               load_balance_parameter, box_size, box_size_inv, kind_of, basis_parameter,i_thread, n_threads,&
                               shm_far_field_boxes, error)
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
         actual_x_data%distribution_forces(bin)%number_of_integrals == 0) THEN
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
 
              !All four centers equivalent => zero-contribution
              IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE 
              IF(katom+latom<=iatom+jatom)  THEN
                IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE

                my_current_counter = my_current_counter + 1
                IF(my_current_counter <= actual_x_data%distribution_forces(bin)%number_of_integrals) THEN
               
                  IF (screen_atom(jatom,iatom) * screen_max < eps_schwarz) CYCLE
                  IF (screen_atom(jatom,iatom) * screen_atom(latom,katom) < eps_schwarz) CYCLE
                  IF(.NOT. (shm_is_assoc_atomic_block(latom,iatom)>=1 .AND. &
                            shm_is_assoc_atomic_block(katom,iatom)>=1 .AND. &
                            shm_is_assoc_atomic_block(katom,jatom)>=1 .AND. &
                            shm_is_assoc_atomic_block(latom,jatom)>=1 ) ) CYCLE
  
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
              
                  l_atom = atom_of_kind(latom)
                  lkind = kind_of(latom)
                  pd1 = last_sgf_global(latom) - nsgfd_total
                  forces_map(4,1) = lkind
                  forces_map(4,2) = l_atom
                  box_idx3 = (katom-1)*natom+latom
  
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
                      max_val1 = screen_set(jset,iset,jatom,iatom)
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
                          pmax_entry = 2.0_dp *MAX(pmax_set(kset,iset,katom,iatom) * pmax_set(lset,jset,latom,jatom),&
                                                   pmax_set(lset,iset,latom,iatom) * pmax_set(kset,jset,katom,jatom))
                          max_val2 = max_val1 * pmax_entry
                          IF(max_val2<eps_schwarz)  CYCLE
                          box_idx4 = (kset-1)*max_set+lset
                          IF(.NOT. do_periodic) THEN
                            center2_folded = shm_far_field_boxes(box_idx4, box_idx3)%center
                          ELSE
                            center2 = shm_far_field_boxes(box_idx4, box_idx3)%center
                            temp = center2-center1
                            CALL pbc(temp, box_size, box_size_inv, pbc_center)
                            center2_folded = center1 + pbc_center
                          END IF
                          R2 = shm_far_field_boxes(box_idx4, box_idx3)%R
                          radius2 = shm_far_field_boxes(box_idx4, box_idx3)%radius
                          Rbox = (center1(1)-center2_folded(1))**2 + &
                                 (center1(2)-center2_folded(2))**2 + &
                                 (center1(3)-center2_folded(3))**2
                          Rbox = SQRT(Rbox) - R1 - R2
                          pgf_radius_box = radius1 + radius2 + 2
  
                          IF(Rbox>pgf_radius_box) THEN
                            max_val2 = 0.0_dp
                            SELECT CASE(potential_parameter%potential_type)
                              CASE(do_hfx_potential_coulomb)
                                Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                       shm_far_field_boxes(box_idx4, box_idx3)%K
                                max_val2 = 1.0_dp/Rbox*Kbox* pmax_entry
                              CASE(do_hfx_potential_short)
                                omega_box = shm_far_field_boxes(box_idx2, box_idx1)%alphaInv + &
                                            shm_far_field_boxes(box_idx4, box_idx3)%alphaInv + &
                                            omega2Inv 
                                omega_box = 1.0_dp/omega_box
                                Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                       shm_far_field_boxes(box_idx4, box_idx3)%K
                                max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox *  pmax_entry
                              CASE(do_hfx_potential_mix_cl)
                                Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                       shm_far_field_boxes(box_idx4, box_idx3)%K
                                max_val2 = 1.0_dp/Rbox*Kbox* pmax_entry
                              CASE (do_hfx_potential_gaussian,do_hfx_potential_mix_lg, do_hfx_potential_id) 
                                max_val2 = HUGE(max_val2)
                            END SELECT
                            IF(max_val2<eps_box) CYCLE
                          END IF
                          IF( .NOT. do_periodic ) THEN
                            max_contraction_val =  max_contraction(iset,iatom) * &
                                                   max_contraction(jset,jatom) * &
                                                   max_contraction(kset,katom) * &
                                                   max_contraction(lset,latom) * &
                                                   pmax_entry
                            CALL forces4(ra, rb, rc, rd,&
                                         zeta(:,iset), zetb(:,jset), zetc(:,kset), zetd(:,lset),&
                                         la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                         lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                         npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                         ncoa, ncob, ncoc, ncod,&
                                         primitive_forces,&
                                         work_forces,&
                                         potential_parameter, eps_schwarz, max_contraction_val, cartesian_estimate,&
                                         i_thread)
                          ELSE
                            max_contraction_val =  max_contraction(iset,iatom) * &
                                                   max_contraction(jset,jatom) * &
                                                   max_contraction(kset,katom) * &
                                                   max_contraction(lset,latom) * &
                                                   pmax_entry
                            CALL forces4_periodic(ra, rb, rc, rd,&
                                                  zeta(:,iset), zetb(:,jset), zetc(:,kset), zetd(:,lset),&
                                                  la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                                  lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                                  npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                                  ncoa, ncob, ncoc, ncod,&
                                                  primitive_forces,&
                                                  work_forces,&
                                                  potential_parameter, box_size, box_size_inv, &
                                                  actual_x_data%neighbor_cells, screen_periodic(:,iset,jset,iatom,jatom),&
                                                  screen_periodic(:,kset,lset,katom,latom), eps_schwarz, &
                                                  max_contraction_val, cartesian_estimate, i_thread, error)
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
    IF( do_periodic ) THEN
      DEALLOCATE(screen_periodic, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(screen_set, screen_atom, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(full_density,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(nspins==2) THEN
      DEALLOCATE(full_density_beta,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(pmax_set,pmax_atom,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

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
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE derivatives_four_center

! *****************************************************************************
!> \brief calculates two-electron derivatived of a quartet/shell using the library 
!>      lib_deriv wrt pbc
!> \param r position
!> \param z zeta
!> \param l angular momentum
!> \param npgf number of primitive cartesian gaussian in actual shell
!> \param primitive_forces array of primitive_forces
!> \param work_forces work_storage
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE forces4_periodic(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                              la_min, la_max, lb_min, lb_max,&
                              lc_min, lc_max, ld_min, ld_max,&
                              npgfa, npgfb, npgfc, npgfd,&
                              ncoa, ncob, ncoc, ncod,&
                              primitive_forces, work_forces,&
                              potential_parameter, box_size, box_size_inv, neighbor_cells, screen1, screen2,&
                              eps_schwarz, max_contraction_val, cartesian_estimate, i_thread, error)
  
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb, zetc, zetd
    INTEGER, INTENT(IN) :: la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
      ld_min, ld_max, npgfa, npgfb, npgfc, npgfd, ncoa, ncob, ncoc, ncod
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod, 12)  :: primitive_forces
    REAL(dp), DIMENSION(*)                   :: work_forces
    TYPE(hfx_potential_type)                 :: potential_parameter
    REAL(dp), INTENT(IN)                     :: box_size(3), box_size_inv(3)
    TYPE(hfx_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(dp), DIMENSION(:)                   :: screen1, screen2
    REAL(dp), INTENT(IN)                     :: eps_schwarz, &
                                                max_contraction_val
    REAL(dp), INTENT(INOUT)                  :: cartesian_estimate
    INTEGER, INTENT(IN)                      :: i_thread
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'forces4_periodic', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, ipgf, j, jpgf, kpgf, la, lb, lc, ld, lpgf, offset_a, &
      offset_a1, offset_b, offset_b1, offset_c, offset_c1, offset_d, &
      offset_d1, stat
    LOGICAL                                  :: failure
    REAL(dp) :: B11(3), C11(3), D11(3), Eta, im_B(3), im_D(3), P(3), &
      pbc_B(3), pbc_D(3), PQ(3), Q(3), shift(3), temp(3), tmp_D(3), &
      tmp_max_all, Zeta1, Zeta_A, Zeta_B, Zeta_C, Zeta_D
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: tmp_cart_estimate

    failure = .FALSE.
    cartesian_estimate = 0.0_dp
    ALLOCATE(tmp_cart_estimate(SIZE(neighbor_cells), SIZE(neighbor_cells)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    tmp_cart_estimate = 0.0_dp 

    primitive_forces = 0.0_dp
    temp = rb -ra
    CALL pbc(temp, box_size, box_size_inv, pbc_B)
    B11 = ra + pbc_B
    temp = rd-rc
    CALL pbc(temp, box_size, box_size_inv, pbc_D)
    D11 = rc + pbc_D
    DO i = 1,SIZE(neighbor_cells)
      im_B(1) = B11(1) + neighbor_cells(i)%cell(1) * box_size(1)
      im_B(2) = B11(2) + neighbor_cells(i)%cell(2) * box_size(2)
      im_B(3) = B11(3) + neighbor_cells(i)%cell(3) * box_size(3)
      DO j = 1,SIZE(neighbor_cells)
        IF(screen1(i)*screen2(j)<eps_schwarz) CYCLE
        im_D(1) = D11(1) + neighbor_cells(j)%cell(1) * box_size(1)
        im_D(2) = D11(2) + neighbor_cells(j)%cell(2) * box_size(2)
        im_D(3) = D11(3) + neighbor_cells(j)%cell(3) * box_size(3)
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
                      CALL pbc(temp,box_size, box_size_inv, PQ)
                      shift = -PQ + temp
                      C11 = rc + shift
                      tmp_D = im_D + shift
                      DO la = la_min,la_max
                        offset_a = offset_a1 + ncoset(la-1)
                        !Build primitives
                        tmp_max_all = 0.0_dp
                        CALL evaluate_deriv_eri(deriv(i_thread+1), ra, im_B, C11, tmp_D,&
                                                zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                                la, lb, lc ,ld,&
                                                work_forces,&
                                                ncoa, ncob, ncoc, ncod, &
                                                primitive_forces, &
                                                offset_a, offset_b, offset_c, offset_d,&
                                                potential_parameter, max_contraction_val, tmp_max_all, eps_schwarz)
                        tmp_cart_estimate(i,j) = MAX(tmp_cart_estimate(i,j),tmp_max_all)
                      END DO !la
                    END DO !ipgf
                  END DO !lb
                END DO !jpgf
              END DO !lc
            END DO !kpgf
          END DO !ld
        END DO !lpgf
      END DO
    END DO
    DO i=1,SIZE(neighbor_cells)
      DO j=1,SIZE(neighbor_cells)
        cartesian_estimate = MAX(cartesian_estimate, tmp_cart_estimate(i,j))
      END DO
    END DO
    DEALLOCATE(tmp_cart_estimate, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE forces4_periodic

END MODULE hfx_derivatives 
