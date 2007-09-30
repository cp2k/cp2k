!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!***** cp2k/four_center_integrals [1.0] *
!!
!!   NAME
!!     four_center_integrals
!!
!!   FUNCTION
!!     routines to calculate Hartree-Fock Exchange and the corresponding energy
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon 
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!****************************************************************************

MODULE four_center_integrals
  
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE erf_fn,                          ONLY: erfc
  USE input_constants,                 ONLY: do_hf_potential_coulomb,&
                                             do_hf_potential_short
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_8
  USE libint_interface,                ONLY: evaluate_deriv_eri,&
                                             evaluate_eri,&
                                             evaluate_eri_screen,&
                                             initialize_libderiv,&
                                             initialize_libint,&
                                             lib_deriv,&
                                             lib_int,&
                                             pbc
  USE mathconstants
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sum_i8,&
                                             mp_sync
  USE orbital_pointers
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE sparse_matrix_types,             ONLY: copy_local_sm_to_replicated_fm,&
                                             first_block_node,&
                                             get_block_node,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE x_hf,                            ONLY: x_hf_box_scheme,&
                                             x_hf_cell_type,&
                                             x_hf_distribution,&
                                             x_hf_potential_type,&
                                             x_hf_set_distr_energy,&
                                             x_hf_set_distr_forces,&
                                             x_hf_type, &
                                             x_hf_screening_type, &
                                             x_hf_potential_type, &
                                             x_hf_general_type
  USE x_hf_distribution_mc,            ONLY: x_hf_optimize_distribution
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC  integrate_four_center, derivatives_four_center

  TYPE(lib_int), SAVE      :: lib
  TYPE(lib_deriv), SAVE    :: deriv
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'four_center_integrals'
  REAL(KIND=dp), PARAMETER :: epsilon = 0.00001_dp


!***
!****************************************************************************
  
  CONTAINS

!!***** four_center_integrals/integrate_four_center [1.0] *
!!
!!   NAME
!!     integrate_four_center
!!
!!   FUNCTION
!!     computes four center integrals for a full basis set and updates the
!!     Kohn-Sham-Matrix and energy. Uses all 8 eri symmetries
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env 
!!     - ks_matrix
!!     - energy
!!     - rho: density matrix
!!     - para_env
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     06.2007 created [Manuel Guidon]
!!     08.2007 optimized load balance [Manuel Guidon] 
!!     09.2007 new parallelization [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE integrate_four_center(qs_env,ks_matrix,energy,rho,hf_section,para_env,& 
                                   geometry_did_change,do_periodic,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: hf_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL                                  :: geometry_did_change, &
                                                do_periodic
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_four_center', &
      routineP = moduleN//':'//routineN

    INTEGER :: bin, box_idx1, box_idx2, box_idx3, box_idx4, current_counter, &
      handle, handle_dist_ks, handle_getP, handle_main, iatom, &
      integral_counter, iset, iw, jatom, jset, katom, kset, latom, lset, ma, &
      max_am, max_k, max_set, maxsgf, mb, mc, md, my_bin_size, my_istart, &
      my_jstart, my_kstart, my_lstart, my_original_cpu_id, natom, nco_max, &
      ncoa, ncob, ncoc, ncod, nseta, nsetb, nsetc, nsetd, nsgfa_total, &
      nsgfb_total, nsgfc_total, nsgfd_total, nspins, &
      number_of_precalc_integrals, p_index, pa, pa1, pa2, pb, pb1, pb2, pc, &
      pc1, pc2, pd, pd1, pd2, pgf_radius_box, sgfa, sgfb, sgfc, sgfd, stat, &
      total_integrals_to_store
    INTEGER :: total_integrals_to_store_aver, total_integrals_to_store_max
    INTEGER(KIND=int_8)                      :: my_current_counter
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: last_sgf_global
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, npgfa, npgfb, npgfc, npgfd, nsgfa, nsgfb, &
      nsgfc, nsgfd
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfc, first_sgfd, &
                                                is_assoc_atomic_block
    LOGICAL :: buffer_overflow, do_correct_p_screening, do_p_screening, &
      failure, first_j_loop, first_k_loop, first_l_loop
    REAL(dp) :: actual_pmax_entry, box_size(3), ene_x_aa, ene_x_bb, &
      eps_schwarz, fac, hf_fraction, K, Kbox, max_val1, max_val2, &
      memory_usage, omega_box, p_screening_correction_factor, pmax_entry, &
      ra(3), rb(3), Rbox, rc(3), rd(3), screen_max, symm_fac, tmp
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: primitive_integrals, T1
    REAL(dp), DIMENSION(:, :), POINTER       :: full_ks, full_ks_beta, &
                                                screen_atom, sphi_a, sphi_b, &
                                                sphi_c, sphi_d, zeta, zetb, &
                                                zetc, zetd
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: actual_pmax_set, screen_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: full_density, &
                                                full_density_beta, &
                                                sparse_block, &
                                                sparse_block_beta
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: screen_periodic
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a, orb_basis_b, &
                                                orb_basis_c, orb_basis_d
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_block_node_type), POINTER      :: block_node, block_node_beta
    TYPE(x_hf_potential_type)                :: potential_parameter
    TYPE(x_hf_screening_type)                :: screening_parameter
    TYPE(x_hf_general_type)                  :: general_parameter
    REAL(dp)                                 :: eps_box


    CALL mp_sync(para_env%group)
    CALL timeset(routineN,"I",'',handle)

    screening_parameter = qs_env%x_data%screening_parameter
    potential_parameter = qs_env%x_data%potential_parameter
    general_parameter   = qs_env%x_data%general_parameter

    !! Get screening paramter
    eps_schwarz = screening_parameter%eps_schwarz
    eps_box     = screening_parameter%eps_box
    do_p_screening = screening_parameter%do_initial_p_screening
    do_correct_p_screening = screening_parameter%do_init_p_screen_corr
    p_screening_correction_factor = screening_parameter%initial_p_screening_correction

    !! If we have a hybrid functional, we may need only a fraction of exact exchange
    hf_fraction = general_parameter%fraction 
   
    !! The number of integrals that fit into the given MAX_MEMORY
    number_of_precalc_integrals = qs_env%x_data%number_of_elements

    !! Parameters related to the potential 1/r, erf(wr)/r, erfc(wr/r)
    potential_parameter = qs_env%x_data%potential_parameter
    
    !! Initialize counter for the integrals stored in memory
    integral_counter = 1
    !! Initialize memory_usage and total number of integrals to be stored on this process
    memory_usage = 0.0_dp
    total_integrals_to_store=0

    !! Variable to check if we calculate the integrals in-core or on the fly
    !! If TRUE -> on the fly
    buffer_overflow = .FALSE.

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
   
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    error=error)
 
    natom = SIZE(particle_set,1)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxsgf=maxsgf,&
                             maxnset=max_set,&
                             maxlgto=max_am)

    !! Initialize memory used by libint
    CALL initialize_libint(lib, max_am, error)
    
    nspins = qs_env%dft_control%nspins

    !! Helper array to map local basis function indeces to global ones 
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO iatom=1,natom
      atom_kind => particle_set(iatom)%atomic_kind
 
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
                             nsgf=nsgfa_total)
      last_sgf_global(iatom) = last_sgf_global(iatom-1)+nsgfa_total
    END DO
 
    !! Get the full density from all the processors
    CALL timeset(routineN//"_getP",handle_getP)
    ALLOCATE(full_density(1,1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL copy_local_sm_to_replicated_fm(rho%rho_ao(1)%matrix,full_density)
    CALL symmetrize_matrix(full_density,"upper_to_lower")
    CALL mp_sum(full_density,para_env%group)
    IF(nspins == 2) THEN
      ALLOCATE(full_density_beta(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL copy_local_sm_to_replicated_fm(rho%rho_ao(2)%matrix,full_density_beta)
      CALL symmetrize_matrix(full_density_beta,"upper_to_lower")
      CALL mp_sum(full_density_beta,para_env%group)
    END IF
    CALL timestop(handle_getP)

    
    !! Calculate the max values of the density matrix actual_pmax stores the data from the actual density matrix
    !! and x_data%initial_p stores the same for the initial guess. The initial guess is updated only in the case of
    !! a changed geometry 
    IF( do_p_screening) THEN
      pmax_entry = 0.0_dp
      actual_pmax_entry = 0.0_dp
      ALLOCATE(actual_pmax_set(max_set,max_set,natom,natom), STAT=stat)
      actual_pmax_set = 0.0_dp
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL update_pmax_matrix(actual_pmax_set, full_density, full_density_beta, &
                              particle_set,atomic_kind_set,error)
      IF(geometry_did_change ) THEN
        qs_env%x_data%initial_p = 0.0_dp
        CALL update_pmax_matrix(qs_env%x_data%initial_p, &
                    full_density,full_density_beta,&
                    particle_set,atomic_kind_set,error)
      END IF
    END IF
    
    !! Get the full ks_matrix from all processors
    ALLOCATE(full_ks(SIZE(full_density,1),SIZE(full_density,2)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    full_ks = 0.0_dp
    IF(nspins==2) THEN
      ALLOCATE(full_ks_beta(SIZE(full_density,1),SIZE(full_density,2)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      full_ks_beta = 0.0_dp
    END IF

    !! Initialize schwarz screening matrices for near field estimates and boxing screening matrices
    !! for far field estimates. The update is only performed if the geomtry of the system changed.
    !! If the system is periodic, then the corresponding routines are called and some variables such as
    !! box_size are initialized

    max_k = qs_env%x_data%number_of_k_loops
    IF( do_periodic ) THEN
      box_size(1) = qs_env%cell%hmat(1,1)
      box_size(2) = qs_env%cell%hmat(2,2)
      box_size(3) = qs_env%cell%hmat(3,3)
      max_k = MIN(max_k, SIZE(qs_env%x_data%neighbor_cells))
      ALLOCATE(screen_periodic(natom,natom,max_set,max_set,SIZE(qs_env%x_data%neighbor_cells)))
      CALL comp_screen_mat_per(screen_set,screen_atom,screen_max,screen_periodic,particle_set,atomic_kind_set,&
                               potential_parameter,.NOT.geometry_did_change,&
                               qs_env%x_data%x_screening_buffer_periodic,&
                               qs_env%x_data%neighbor_cells,box_size,max_k,error)
      IF(geometry_did_change) THEN
        CALL screen_farfield_per(qs_env%x_data%far_field_boxes, &
                                 particle_set, atomic_kind_set, potential_parameter, &
                                 qs_env%x_data%neighbor_cells, box_size, error)
       END IF

    ELSE
      CALL compute_screening_matrices(screen_set,screen_atom,screen_max,particle_set,atomic_kind_set,&
                                    potential_parameter,.NOT.geometry_did_change,qs_env%x_data%x_screening_buffer, error)
      IF(geometry_did_change) THEN
        CALL screen_farfield(qs_env%x_data%far_field_boxes,&
                            particle_set, atomic_kind_set, potential_parameter,&
                            error)
      END IF
    END IF
    CALL timeset(routineN//"_main",handle_main)

    !! Initialize a prefactor depending on the fraction and the number of spins
    IF(nspins == 1) THEN
      fac = 0.5_dp * hf_fraction
    ELSE
      fac = 1.0_dp * hf_fraction
    END IF
    
    !! Calculate helper array that stores if a certain atomic pair is associated in the KS matrix
    ALLOCATE(is_assoc_atomic_block(natom,natom))
    is_assoc_atomic_block = 0
    DO iatom=1,ks_matrix(1)%matrix%nblock_row
      block_node => first_block_node(ks_matrix(1)%matrix,iatom)
      DO WHILE (ASSOCIATED(block_node))
        CALL get_block_node(block_node=block_node,&
                            block_col=jatom,&
                            BLOCK=sparse_block)

        is_assoc_atomic_block(iatom,jatom) = 1
        is_assoc_atomic_block(jatom,iatom) = 1
        block_node => next_block_node(block_node)
      END DO
    END DO
    CALL mp_sum(is_assoc_atomic_block,para_env%group)

    !! Call routines that distribute the load on all processes. If we want to screen on a initial density matrix, there is
    !! an optional parameter. Of course, this is only done if the geometry did change
    IF( geometry_did_change ) THEN
      IF( do_p_screening) THEN
        CALL load_balance_energy(qs_env%x_data,eps_schwarz,eps_box,particle_set,atomic_kind_set,para_env,&
                                 potential_parameter,screen_set,screen_atom,screen_max,&
                                 is_assoc_atomic_block,ks_matrix,do_periodic,pmax_set=qs_env%x_data%initial_p,error=error)
      ELSE
        CALL load_balance_energy(qs_env%x_data,eps_schwarz,eps_box,particle_set,atomic_kind_set,para_env,&
                                 potential_parameter,screen_set,screen_atom,screen_max,&
                                 is_assoc_atomic_block,ks_matrix,do_periodic,error=error)
      END IF
    END IF

    !! precompute maximum nco and allocate scratch
    nco_max=0
    DO iatom=1,natom
      atom_kind => particle_set(iatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a, lmax=la_max, lmin=la_min,&
                             npgf=npgfa, nset=nseta, zet=zeta, nsgf_set=nsgfa,&
                             first_sgf=first_sgfa, sphi=sphi_a, nsgf=nsgfa_total)
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
    
    
    my_bin_size = SIZE(qs_env%x_data%distribution_energy)
    first_j_loop = .TRUE.
    first_l_loop = .TRUE.
    first_k_loop = .TRUE.
    !! outer loop will be cycled if actual work_load done
    outer: DO bin=1,my_bin_size
      !! Get information on workload
      my_istart = qs_env%x_data%distribution_energy(bin)%istart
      my_jstart = qs_env%x_data%distribution_energy(bin)%jstart
      my_kstart = qs_env%x_data%distribution_energy(bin)%kstart
      my_lstart = qs_env%x_data%distribution_energy(bin)%lstart
      my_original_cpu_id = qs_env%x_data%distribution_energy(bin)%original_cpu_id 
      my_current_counter = 0
      !! It could happen, that one bin is empty. In that case just cycle
      IF(my_istart == 0 .OR. my_jstart == 0 .OR. my_kstart == 0 .OR. my_lstart == 0 .OR. &
         qs_env%x_data%distribution_energy(bin)%number_of_integrals == 0) CYCLE
      DO iatom=my_istart,natom
        !! Be consistent to the modulo distribution of the load balance
        IF(MOD(iatom,para_env%num_pe) == my_original_cpu_id-1) THEN
          atom_kind => particle_set(iatom)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atom_kind,&
                               orb_basis_set=orb_basis_a)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_a, lmax=la_max, lmin=la_min,&
                                 npgf=npgfa, nset=nseta, zet=zeta, nsgf_set=nsgfa,&
                                 first_sgf=first_sgfa, sphi=sphi_a, nsgf=nsgfa_total)
          ra=particle_set(iatom)%r(:)
          pa1 = last_sgf_global(iatom) - nsgfa_total
   
          IF(first_j_loop) jatom=my_jstart-1
          IF(.NOT. first_j_loop) jatom = iatom-1
          DO WHILE(jatom+1<=natom)
            jatom = jatom + 1
            IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
            atom_kind => particle_set(jatom)%atomic_kind
            CALL get_atomic_kind(atomic_kind=atom_kind,&
                                 orb_basis_set=orb_basis_b)
            CALL get_gto_basis_set(gto_basis_set=orb_basis_b, lmax=lb_max, lmin=lb_min,&
                                   npgf=npgfb, nset=nsetb, zet=zetb, nsgf_set=nsgfb,&
                                   first_sgf=first_sgfb, sphi=sphi_b, nsgf=nsgfb_total)
            rb=particle_set(jatom)%r(:)
            pb1 = last_sgf_global(jatom) - nsgfb_total

        
            IF(first_k_loop) katom=my_kstart-1
            IF(.NOT. first_k_loop) katom = 0
            DO WHILE(katom+1<=natom)
              katom = katom + 1
              atom_kind => particle_set(katom)%atomic_kind
              CALL get_atomic_kind(atomic_kind=atom_kind,&
                                   orb_basis_set=orb_basis_c)
              CALL get_gto_basis_set(gto_basis_set=orb_basis_c, lmax=lc_max, lmin=lc_min,&
                                     npgf=npgfc, nset=nsetc, zet=zetc, nsgf_set=nsgfc,&
                                   first_sgf=first_sgfc, sphi=sphi_c, nsgf=nsgfc_total)
              rc=particle_set(katom)%r(:)
              pc1 = last_sgf_global(katom) - nsgfc_total

              IF(first_l_loop) THEN
                latom=my_lstart-1
              ELSE
                latom = katom-1
              END IF
              DO WHILE(latom+1<=natom)
                latom = latom + 1
                !! we want to be consisten with the KS matrix. If none of the atomic indices
                !! is associated cycle
                IF(.NOT. (is_assoc_atomic_block(iatom,latom)>=1 .AND. &
                          is_assoc_atomic_block(iatom,katom)>=1 .AND. &
                          is_assoc_atomic_block(jatom,katom)>=1 .AND. &
                          is_assoc_atomic_block(jatom,latom)>=1 ) ) CYCLE

                !! third permutational symmetry
                IF(katom+latom<=iatom+jatom)  THEN
                  IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                  
                  !! screen near field atomic part
                  IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
                  
                  !! if we pass through here, we have to update the workload counter
                  my_current_counter = my_current_counter + 1
                  !! if we reach the counter in the distribution we are finished with the curren bunch
                  !! of integrals. We also have to reset the flags for the main loops
                  IF(my_current_counter > qs_env%x_data%distribution_energy(bin)%number_of_integrals) THEN
                    first_j_loop = .TRUE.
                    first_l_loop = .TRUE.
                    first_k_loop = .TRUE.
                    CYCLE outer
                  END IF

                  !! calculate symmetry_factor accordin to degeneracy of atomic quartet
                  symm_fac = 0.5_dp
                  IF(iatom==jatom) symm_fac = symm_fac*2.0_dp
                  IF(katom==latom) symm_fac = symm_fac*2.0_dp
                  IF(iatom==katom .AND. jatom==latom .AND. iatom/=jatom .AND. katom/=latom) symm_fac = symm_fac*2.0_dp
                  IF(iatom==katom .AND. iatom==jatom .AND. katom==latom) symm_fac = symm_fac*2.0_dp
                  symm_fac = 1.0_dp / symm_fac

                  atom_kind => particle_set(latom)%atomic_kind
                  CALL get_atomic_kind(atomic_kind=atom_kind,&
                                       orb_basis_set=orb_basis_d)
                  CALL get_gto_basis_set(gto_basis_set=orb_basis_d, lmax=ld_max, lmin=ld_min,&
                                         npgf=npgfd, nset=nsetd, zet=zetd, nsgf_set=nsgfd,&
                                         first_sgf=first_sgfd, sphi=sphi_d, nsgf=nsgfd_total)
                  rd=particle_set(latom)%r(:)
                  pd1 = last_sgf_global(latom) - nsgfd_total
             
                  !! Loop over all sets of the current atomic quartet
                  DO iset = 1, nseta
                    ncoa = npgfa(iset)*ncoset(la_max(iset))
                    sgfa = first_sgfa(1,iset)
                    DO jset = 1,nsetb
                      ncob = npgfb(jset)*ncoset(lb_max(jset))
                      sgfb = first_sgfb(1,jset)
                      max_val1 = screen_set(iatom,jatom,iset,jset)
                      !! Near field screening
                      IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                      DO kset = 1,nsetc
                        ncoc = npgfc(kset)*ncoset(lc_max(kset))
                        sgfc = first_sgfc(1,kset)
                        DO lset = 1,nsetd
                          ncod = npgfd(lset)*ncoset(ld_max(lset)) 
                          sgfd = first_sgfd(1,lset)
                          max_val2 = max_val1 * screen_set(katom,latom,kset,lset)
                          !! Near field screening
                          IF(max_val2<eps_schwarz) CYCLE

                          !! get max_vals if we screen on initial density
                          IF( do_p_screening) THEN
                            pmax_entry = MAX(qs_env%x_data%initial_p(iset,kset,iatom,katom), &
                                             qs_env%x_data%initial_p(jset,lset,jatom,latom), &
                                             qs_env%x_data%initial_p(iset,lset,iatom,latom), &
                                             qs_env%x_data%initial_p(kset,jset,katom,jatom))
                            actual_pmax_entry = MAX(actual_pmax_set(iset,kset,iatom,katom), &
                                                    actual_pmax_set(jset,lset,jatom,latom), &
                                                    actual_pmax_set(iset,lset,iatom,latom), &
                                                    actual_pmax_set(kset,jset,katom,jatom))
                            max_val2 = max_val2 * pmax_entry
                          END IF
                          !! If wanted, correct the previous screening by recalculating the omitted
                          !! integrals on the fly.
                          IF(max_val2<eps_schwarz) THEN
                            IF(do_p_screening) THEN
                              IF( do_correct_p_screening ) THEN
                                IF(pmax_entry*p_screening_correction_factor < actual_pmax_entry) THEN
                                  CALL on_the_fly(do_periodic, ra, rb, rc, rd, zeta, zetb, zetc, zetd, &
                                                  la_min, la_max, lb_min, lb_max, lc_min, lc_max, ld_min, ld_max, &
                                                  npgfa, npgfb, npgfc, npgfd, ncoa, ncob, ncoc, ncod, &
                                                  nsgfa, nsgfb, nsgfc, nsgfd, sgfa, sgfb, sgfc, sgfd, &
                                                  primitive_integrals, potential_parameter, box_size, qs_env, screen_periodic, &
                                                  eps_schwarz, max_k, sphi_a, sphi_b, sphi_c, sphi_d, T1, &
                                                  full_ks, full_ks_beta, full_density, full_density_beta, fac, symm_fac, &
                                                  iset, jset, kset, lset, iatom, jatom, katom, latom, pa1, pb1, pc1, pd1, nspins)
                                END IF
                              ENDIF
                            END IF
                            CYCLE
                          END IF

                          !! Far field box screening also taking into account screening on initial density matrix
                          box_idx1 = (iatom-1)*natom+jatom
                          box_idx2 = (iset-1)*max_set+jset
                          box_idx3 = (katom-1)*natom+latom
                          box_idx4 = (kset-1)*max_set+lset
                          Rbox = (qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%center(1) - &
                                  qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%center(1))**2 + &
                                 (qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%center(2) - &
                                  qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%center(2))**2 + &
                                 (qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%center(3) - &
                                  qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%center(3))**2
                          Rbox = SQRT(Rbox)
                          Rbox = Rbox - qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                        qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%R
                          pgf_radius_box = INT(1.0_dp/ &
                                (SQRT(2.0_dp*qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                          pgf_radius_box = pgf_radius_box + INT(1.0_dp/ & 
                                (SQRT(2.0_dp*qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2
                          IF(Rbox>pgf_radius_box) THEN
                            max_val2 = 0.0_dp
                            SELECT CASE(potential_parameter%potential_type)
                              CASE(do_hf_potential_coulomb)
                                Kbox = qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = 1.0_dp/Rbox*Kbox
                                IF(max_val2 < eps_box ) CYCLE
                                IF( do_p_screening) THEN
                                  max_val2 = max_val2 * pmax_entry
                                END IF
                              CASE(do_hf_potential_short)
                                omega_box = 1.0_dp/qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                            1.0_dp/qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                            1.0_dp/(potential_parameter%omega**2)
                                omega_box = 1.0_dp/omega_box
                                Kbox = qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox
                                IF(max_val2 < eps_box ) CYCLE
                                IF( do_p_screening) THEN
                                  max_val2 = max_val2 * pmax_entry
                                END IF
                            END SELECT
                            IF(max_val2 < eps_box) THEN
                              IF(do_p_screening) THEN
                                IF( do_correct_p_screening ) THEN
                                  IF(pmax_entry*p_screening_correction_factor < actual_pmax_entry) THEN
                                    CALL on_the_fly(do_periodic, ra, rb, rc, rd, zeta, zetb, zetc, zetd, &
                                                    la_min, la_max, lb_min, lb_max, lc_min, lc_max, ld_min, ld_max, &
                                                    npgfa, npgfb, npgfc, npgfd, ncoa, ncob, ncoc, ncod, &
                                                    nsgfa, nsgfb, nsgfc, nsgfd, sgfa, sgfb, sgfc, sgfd, &
                                                    primitive_integrals, potential_parameter, box_size, qs_env, screen_periodic, &
                                                    eps_schwarz, max_k, sphi_a, sphi_b, sphi_c, sphi_d, T1, &
                                                    full_ks, full_ks_beta, full_density, full_density_beta, fac, symm_fac, &
                                                    iset, jset, kset, lset, iatom, jatom, katom, latom, pa1, pb1, pc1, pd1, nspins)
                                  END IF 
                                END IF
                              END IF
                              CYCLE
                            END IF
                          END IF

                          !! store current number of integrals, update total number and number of integrals in buffer
                          current_counter = nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                        
                          total_integrals_to_store=total_integrals_to_store+current_counter
                          buffer_overflow = integral_counter + current_counter - 1 > number_of_precalc_integrals
 
                          !! Get integrals from buffer and update Kohn-Sham matrix
                          IF(.NOT.buffer_overflow .AND. .NOT.geometry_did_change) THEN 
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
                                    tmp = fac*symm_fac*qs_env%x_data%x_buffer(integral_counter)
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
                                    integral_counter = integral_counter + 1
                                  END DO
                                END DO
                              END DO
                            END DO
                          END IF
    
                          !! Calculate integrals if we run out of buffer or the geometry did change
                          IF(geometry_did_change .OR. buffer_overflow) THEN
                            IF( .NOT.do_periodic ) THEN 
                              CALL coulomb4(ra, rb, rc, rd,&
                                            zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset), &
                                            zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                                            la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                            lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                            npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                            ncoa,ncob,ncoc,ncod,&
                                            primitive_integrals(1),&
                                            potential_parameter)
                            ELSE
                              CALL coulomb4_periodic(ra, rb, rc, rd,&
                                                     zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset), &
                                                     zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                                                     la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                                     lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                                     npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                                     ncoa,ncob,ncoc,ncod,&
                                                     primitive_integrals(1),&
                                                     potential_parameter, box_size, &
                                                     qs_env%x_data%neighbor_cells, screen_periodic(iatom,jatom,iset,jset,:),&
                                                     screen_periodic(katom,latom,kset,lset,:), eps_schwarz, max_k)
                            END IF
  
                            !! contract primitive cartesian integrals to contracted spherical integrals
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


                            !! Update Kohn-Sham matrix
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
                                    IF(geometry_did_change .AND. (.NOT.buffer_overflow)) THEN
                                      qs_env%x_data%x_buffer(integral_counter) = &
                                      primitive_integrals(p_index)
                                      integral_counter = integral_counter + 1
                                    END IF
                                  END DO
                                END DO
                              END DO
                            END DO
                          END IF
                        END DO !lset
                      END DO !kset
                    END DO !jset
                  END DO !iset

                END IF !symmetry cycle 
              END DO !latom
              first_l_loop = .FALSE.
            END DO !katom
            first_k_loop = .FALSE.
          END DO !jatom
          first_j_loop = .FALSE.
        ENDIF
      END DO !iatom
      !! Because of the Monte Carlo optimization it is not guaranteed that there are no more bunches to calculate.
      !! Therefore we have to reset the flags
      first_j_loop = .TRUE.
      first_k_loop = .TRUE.
      first_l_loop = .TRUE.
    END DO outer !bin

    DEALLOCATE(primitive_integrals,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(T1,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    


    !! Calculate the exchange energies from the Kohn-Sham matrix. Before we can go on, we have to symmetrize.
    ene_x_aa = 0.0_dp
    ene_x_bb = 0.0_dp
   
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
    
    CALL mp_sum(full_ks,para_env%group)
    IF(nspins==2) CALL mp_sum(full_ks_beta,para_env%group)
    CALL mp_sum(ene_x_aa, para_env%group)
    IF(nspins==2) CALL mp_sum(ene_x_bb,para_env%group)


    !! Distribute the current KS-matrix to all the processes
    CALL timeset(routineN//"_dist_KS",handle_dist_ks)
    DO iatom=1,ks_matrix(1)%matrix%nblock_row
      block_node => first_block_node(ks_matrix(1)%matrix,iatom)
      IF(nspins == 2) block_node_beta => first_block_node(ks_matrix(2)%matrix,iatom)
      ! parallel distribution might mean that there won't be an associated latom
      IF (.NOT. ASSOCIATED(block_node)) CYCLE
      atom_kind => particle_set(iatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
                             nsgf=nsgfa_total,nset=nseta,&
                             first_sgf=first_sgfa,&
                             nsgf_set=nsgfa)
      pa1 = last_sgf_global(iatom) - nsgfa_total
      block_node => first_block_node(ks_matrix(1)%matrix,iatom)
      IF(nspins == 2) block_node_beta => first_block_node(ks_matrix(2)%matrix,iatom)
      DO WHILE (ASSOCIATED(block_node))
        CALL get_block_node(block_node=block_node,&
                            block_col=latom,&
                            BLOCK=sparse_block)
        block_node => next_block_node(block_node)
        IF(nspins == 2) THEN
          CALL get_block_node(block_node=block_node_beta,&
                              block_col=latom,&
                              BLOCK=sparse_block_beta)
          block_node_beta => next_block_node(block_node_beta)
        END IF
        atom_kind => particle_set(latom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_d)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_d,&
                               nsgf=nsgfd_total,nset=nsetd,&
                               first_sgf=first_sgfd,&
                               nsgf_set=nsgfd)
        pd1 = last_sgf_global(latom) - nsgfd_total
        DO lset=1,nsetd
          sgfd = first_sgfd(1,lset)
          DO iset=1,nseta
            sgfa = first_sgfa(1,iset)
            DO md = 1,nsgfd(lset)
              pd2 = md+sgfd-1
              pd = pd1 + pd2
              DO ma = 1,nsgfa(iset)
                pa2 = ma+sgfa-1
                pa = pa1 + pa2
                sparse_block(pa2,pd2) = full_ks(pa,pd)
                IF(nspins==2) sparse_block_beta(pa2,pd2) = full_ks_beta(pa,pd)
              END DO 
            END DO
          END DO    
        END DO
      END DO
    END DO
    
    CALL timestop(handle_dist_ks)
    !! Update energy type
    energy%ex = 0.5_dp*(ene_x_aa+ene_x_bb)

    CALL timestop(handle_main)
   

    !! Print some information 
    iw = cp_print_key_unit_nr(logger,hf_section,"HF_INFO",&
             extension=".scfLog",error=error)
    total_integrals_to_store_max=total_integrals_to_store
    total_integrals_to_store_aver=(total_integrals_to_store+para_env%num_pe-1)/(para_env%num_pe)
    CALL mp_max(total_integrals_to_store_max,para_env%group)
    CALL mp_sum(total_integrals_to_store_aver,para_env%group)

    IF (iw>0) THEN
      IF(memory_usage < integral_counter * dp_size /1024_dp/1024_dp) THEN
        memory_usage = integral_counter * dp_size /1024_dp/1024_dp
      END IF
      WRITE (UNIT=iw,FMT="(/,(T3,A,T60,F20.10))")&
               "Memory used for x_buffer [MB]:                 ",memory_usage
      WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
               "Number of integrals stored in buffer:          ",integral_counter-1

      WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
               "Maximum number of elements to store:          ",total_integrals_to_store_max

      WRITE (UNIT=iw,FMT="((T3,A,T60,I20),/)")&
               "Average number of elements to store:          ",total_integrals_to_store_aver
    END IF

    CALL cp_print_key_finished_output(iw,logger,hf_section,&
                                      "HF_INFO", error=error)


    !! Clean up
    IF( do_periodic ) THEN
      DEALLOCATE(screen_periodic, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    DEALLOCATE(screen_set, screen_atom, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(full_ks, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(last_sgf_global,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(full_density,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(is_assoc_atomic_block,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(nspins==2) THEN
      DEALLOCATE(full_density_beta,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(full_ks_beta, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(do_p_screening) THEN
      DEALLOCATE(actual_pmax_set,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    !
    ! this wraps to free_libint, but currently causes a segfault
    ! as a result, we don't call it, but some memory remains allocated
    !  
    ! CALL terminate_libint(lib)

    ! many timings somewhat more clear
    CALL mp_sync(para_env%group)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE integrate_four_center
!***************************************************************************


!!***** four_center_integrals/coulomb4 [1.0] *
!!
!!   NAME
!!     coulomb4
!!
!!   FUNCTION
!!     calculates two-electron integrals of a quartet/shell using the library 
!!     lib_int
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - r: position
!!     - z: zeta
!!     - l: angular momentum
!!     - npgf: number of primitive cartesian gaussian in actual shell
!!     - primitive_integrals: array of primitive_integrals
!!     - potential_parameter: contains info for libint 
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE coulomb4(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                      la_min, la_max, lb_min, lb_max,&
                      lc_min, lc_max, ld_min, ld_max,&
                      npgfa, npgfb, npgfc, npgfd,&
                      ncoa,ncob,ncoc,ncod,&
                      primitive_integrals,&
                      potential_parameter)
  
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
    TYPE(x_hf_potential_type)                :: potential_parameter

    INTEGER :: ipgf, jpgf, kpgf, la, lb, lc, ld, lpgf, offset_a, offset_a1, &
      offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1

    primitive_integrals = 0.0_dp

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
                    CALL evaluate_eri(lib, ra, rb, rc, rd,&
                                      zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                      la, lb, lc ,ld,&
                                      ncoa,ncob,ncoc,ncod,&
                                      offset_a,offset_b,offset_c,offset_d, &
                                      primitive_integrals, &
                                      potential_parameter)
                  END DO !la
                END DO !ipgf
              END DO !lb
            END DO !jpgf
          END DO !lc
        END DO !kpgf
      END DO !ld
    END DO !lpgf
  END SUBROUTINE coulomb4
!***************************************************************************

!!***** four_center_integrals/screen4 [1.0] *
!!
!!   NAME
!!     screen4
!!
!!   FUNCTION
!!     calculates max values of two-electron integrals in a quartet/shell
!!     w.r.t. different zetas using the library lib_int
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - r: position
!!     - z: zeta
!!     - l: angular momentum
!!     - npgf: number of primitive cartesian gaussian in actual shell
!!     - max_val: schwarz screening value
!!     - potential_parameter: contains info for libint
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!    03.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE screen4(ra, rb, zeta, zetb, &
                      la_min, la_max, lb_min, lb_max,&
                      npgfa, npgfb, &
                      max_val, potential_parameter)

    REAL(dp), INTENT(IN)                     :: ra(3), rb(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb
    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, &
                                                lb_max, npgfa, npgfb
    REAL(dp), INTENT(INOUT)                  :: max_val
    TYPE(x_hf_potential_type)                :: potential_parameter

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'screen4'

    INTEGER                                  :: ipgf, jpgf, la, lb
    REAL(dp)                                 :: max_val_temp

    max_val_temp = max_val
    DO ipgf = 1,npgfa
      DO jpgf = 1,npgfb
        DO la = la_min,la_max
          DO lb = lb_min, lb_max
            !Build primitives
            CALL evaluate_eri_screen(lib, ra, rb, ra, rb,&
                                     zeta(ipgf), zetb(jpgf), zeta(ipgf), zetb(jpgf),&
                                     la, lb, la ,lb,&
                                     max_val_temp, potential_parameter)
            max_val = MAX(max_val, max_val_temp)
          END DO !lb
        END DO !la
      END DO !jpgf
    END DO !ipgf

  END SUBROUTINE screen4


!!***** four_center_integrals/compute_screening_matrices [1.0] *
!!
!!   NAME
!!     compute_screening_matrices
!!
!!   FUNCTION
!!     calculates the screening matrices 
!!
!!   NOTES
!!     - already computes the square root
!!     - also takes the primitive to sgf transformation into account
!!     - in this way all stored contracted integrals are smaller or equal to their estimate
!!   INPUTS
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!    05.2007 refactored [Joost VandeVondele]
!!
!!****************************************************************************

  SUBROUTINE compute_screening_matrices(screen_set,screen_atom,screen_max,particle_set,&
                                        atomic_kind_set,potential_parameter,read_from_memory,buffer,error)

    REAL(kind=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: screen_set
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: screen_atom
    REAL(KIND=dp)                            :: screen_max
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(x_hf_potential_type)                :: potential_parameter
    LOGICAL                                  :: read_from_memory
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      OPTIONAL, POINTER                      :: buffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_screening_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, jatom, jset, k, &
                                                katom, kset, max_set, natom, &
                                                ncob, nsetb, nsetc, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: lb_max, lb_min, lc_max, &
                                                lc_min, npgfb, npgfc, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: contraction_max
    REAL(dp), DIMENSION(:, :), POINTER       :: sphi_b, zetb, zetc
    REAL(KIND=dp)                            :: max_val, p2s_b, p2s_c, rb(3), &
                                                rc(3)
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_b, orb_basis_c

    CALL timeset(routineN,handle)

    failure = .FALSE.

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set)
    ALLOCATE(screen_set(natom,natom,max_set,max_set),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(screen_atom(natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(contraction_max(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    screen_set=0.0_dp
    screen_atom=0.0_dp
    contraction_max=0.0_dp

    DO jatom=1,natom

      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             first_sgf=first_sgfb,&
                             sphi=sphi_b,&
                             nsgf_set=nsgfb)
      DO jset = 1,nsetb
        ! takes the primitive to contracted transformation into account
        ncob = npgfb(jset)*ncoset(lb_max(jset))
        sgfb = first_sgfb(1,jset)
        ! if the primitives are assumed to be all of max_val2, max_val2*p2s_b becomes
        ! the maximum value after multiplication with sphi_b
        contraction_max(jset,jatom) = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,k))),k=sgfb,sgfb+nsgfb(jset)-1)/))
      ENDDO
    ENDDO

    DO jatom=1,natom

      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             lmin=lb_min,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             zet=zetb)
      rb=particle_set(jatom)%r(:)
      DO katom = 1,natom
        atom_kind => particle_set(katom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_c)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                               lmax=lc_max,&
                               lmin=lc_min,&
                               npgf=npgfc,&
                               nset=nsetc,&
                               zet=zetc)
        rc=particle_set(katom)%r(:)
        DO jset = 1,nsetb
          p2s_b = contraction_max(jset,jatom)
          DO kset = 1,nsetc
            p2s_c = contraction_max(kset,katom)
            max_val = 0.0_dp
            IF( read_from_memory ) THEN
              screen_set(jatom,katom,jset,kset) = buffer(jatom,katom,jset,kset)
            ELSE
              CALL screen4(rb, rc, &
                            zetb(:,jset), zetc(:,kset), &
                            lb_min(jset), lb_max(jset), lc_min(kset), lc_max(kset),&
                            npgfb(jset), npgfc(kset), &
                            max_val, potential_parameter)
              screen_set(jatom,katom,jset,kset) = SQRT(max_val)*p2s_b*p2s_c
              IF( PRESENT(buffer) ) buffer(jatom,katom,jset,kset) = screen_set(jatom,katom,jset,kset)
            END IF
          END DO
        END DO
        screen_atom(jatom,katom)=MAXVAL(screen_set(jatom,katom,1:nsetb,1:nsetc))
      END DO
    END DO
    screen_max=MAXVAL(screen_atom)
    DEALLOCATE(contraction_max, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
  END SUBROUTINE compute_screening_matrices

!!***** four_center_integrals/compute_pmax_matrices [1.0] *
!!
!!   NAME
!!     compute_pmax_matrices
!!
!!   FUNCTION
!!     calculates the maximum of the density matrix in compressed form for screening purposes 
!!
!!   NOTES
!!     - computes for each pair of shells the maximum absolute value of p
!!
!!   INPUTS
!!     - screen_pmat_forces: if false just returns a matrix of 1.0_dp
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!    05.2007 created [Joost VandeVondele]
!!
!****************************************************************************
  
  SUBROUTINE compute_pmax_matrices(pmax_set,pmax_atom,pmax_max, &
              full_density,full_density_beta,screen_pmat_forces,particle_set,atomic_kind_set,error)

    REAL(kind=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: pmax_set
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: pmax_atom
    REAL(KIND=dp)                            :: pmax_max
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: full_density, &
                                                full_density_beta
    LOGICAL, INTENT(IN)                      :: screen_pmat_forces
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_pmax_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, jatom, jset, katom, kset, max_set, mb, mc, natom, &
      nsetb, nsetc, nsgfb_total, nsgfc_total, pb, pb1, pb2, pc, pc1, pc2, &
      sgfb, sgfc, stat
    INTEGER, DIMENSION(:), POINTER           :: last_sgf_global, nsgfb, nsgfc
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb, first_sgfc
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: pmax_tmp
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_b, orb_basis_c

    CALL timeset(routineN,handle)
    failure = .FALSE.

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set)

    ! create a mapping array of easy access to the density matrix elements
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO jatom=1,natom
      atom_kind => particle_set(jatom)%atomic_kind
 
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             nsgf=nsgfb_total)
      last_sgf_global(jatom) = last_sgf_global(jatom-1)+nsgfb_total
    END DO

    ! allocate and compute pmax
    ALLOCATE(pmax_set(max_set,max_set,natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(pmax_atom(natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (screen_pmat_forces) THEN

       pmax_set=0.0_dp
       pmax_atom=0.0_dp

       DO jatom=1,natom

         atom_kind => particle_set(jatom)%atomic_kind
         CALL get_atomic_kind(atomic_kind=atom_kind,&
                              orb_basis_set=orb_basis_b)
         CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                                nsgf=nsgfb_total,&
                                nset=nsetb,& 
                                first_sgf=first_sgfb,&
                                nsgf_set=nsgfb)
         pb1 = last_sgf_global(jatom) - nsgfb_total

         DO katom = 1,natom
           atom_kind => particle_set(katom)%atomic_kind
           CALL get_atomic_kind(atomic_kind=atom_kind,&
                                orb_basis_set=orb_basis_c)
           CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                                  nsgf=nsgfc_total,&
                                  nset=nsetc,&
                                  first_sgf=first_sgfc,&
                                  nsgf_set=nsgfc)
           pc1 = last_sgf_global(katom) - nsgfc_total

           DO jset = 1,nsetb
             sgfb = first_sgfb(1,jset)
             DO kset = 1,nsetc
               sgfc = first_sgfc(1,kset)
               pmax_tmp = 0.0_dp
               DO mc = 1,nsgfc(kset)
                  pc2 = mc+sgfc-1
                  pc = pc1 + pc2
                  DO mb = 1,nsgfb(jset)
                     pb2 = mb+sgfb-1
                     pb = pb1 + pb2
                     pmax_tmp=MAX(pmax_tmp,ABS(full_density(pb,pc)))
                     IF( ASSOCIATED(full_density_beta) ) pmax_tmp=MAX(pmax_tmp,ABS(full_density_beta(pb,pc)))
                  ENDDO
               ENDDO
               pmax_set(jset,kset,jatom,katom) = pmax_tmp
             END DO
           END DO
           pmax_atom(jatom,katom)=MAXVAL(pmax_set(1:nsetb,1:nsetc,jatom,katom))
         END DO
       END DO

    ELSE

       pmax_set=1.0_dp
       pmax_atom=1.0_dp

    ENDIF

    pmax_max=MAXVAL(pmax_atom)

    DEALLOCATE(last_sgf_global,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE compute_pmax_matrices
!***************************************************************************

!!***** four_center_integrals/forces4 [1.0] *
!!
!!   NAME
!!     forces4
!!
!!   FUNCTION
!!     calculates two-electron derivatived of a quartet/shell using the library 
!!     lib_deriv
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - r: position
!!     - z: zeta
!!     - l: angular momentum
!!     - npgf: number of primitive cartesian gaussian in actual shell
!!     - primitive_forces: array of primitive_forces
!!     - work_forces: work_storage
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     03.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE forces4(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                    la_min, la_max, lb_min, lb_max,&
                    lc_min, lc_max, ld_min, ld_max,&
                    npgfa, npgfb, npgfc, npgfd,&
                    ncoa, ncob, ncoc, ncod,&
                    primitive_forces, work_forces, eps_schwarz,&
                    potential_parameter )
  
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb, zetc, zetd
    INTEGER, INTENT(IN) :: la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
      ld_min, ld_max, npgfa, npgfb, npgfc, npgfd, ncoa, ncob, ncoc, ncod
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod, 12)  :: primitive_forces
    REAL(dp), DIMENSION(:, :), POINTER       :: work_forces
    REAL(dp), INTENT(IN)                     :: eps_schwarz
    TYPE(x_hf_potential_type)                :: potential_parameter

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'forces4'

    INTEGER :: ipgf, jpgf, kpgf, la, lb, lc, ld, lpgf, offset_a, offset_a1, &
      offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1

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
                    CALL evaluate_deriv_eri(deriv, ra, rb, rc, rd,&
                                           zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                           la, lb, lc ,ld,&
                                           work_forces,&
                                           ncoa, ncob, ncoc, ncod, &
                                           primitive_forces, &
                                           offset_a, offset_b, offset_c, offset_d,&
                                           potential_parameter)

                  END DO !ld
                END DO !lc
              END DO !lb
            END DO !la

          END DO !lpgf
        END DO !kpgf
      END DO !jpgf
    END DO !ipgf
  
  END SUBROUTINE forces4
!***************************************************************************

!!***** four_center_integrals/derivatives_four_center [1.0] *
!!
!!   NAME
!!     derivatives_four_center
!!
!!   FUNCTION
!!     computes four center derivatives for a full basis set and updates the
!!     forces%fock_4c arrays. Uses all 8 eri symmetries
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env 
!!     - ks_matrix
!!     - energy
!!     - rho: density matrix
!!     - para_env
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     06.2007 created [Manuel Guidon]
!!     08.2007 optimized load balance [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE derivatives_four_center(qs_env,ks_matrix,energy,rho,hf_section,para_env,& 
                                    do_periodic, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: hf_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: do_periodic
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derivatives_four_center', &
      routineP = moduleN//':'//routineN

    INTEGER :: bin, box_idx1, box_idx2, box_idx3, box_idx4, coord, &
      forces_map(4,2), handle, i_atom, iatom, ikind, iset, j_atom, jatom, &
      jkind, jset, k_atom, katom, kkind, kset, l_atom, latom, lkind, lset, &
      ma, max_am, max_k, max_set, maxsgf, mb, mc, md, my_bin_size, my_istart, &
      my_jstart, my_kstart, my_lstart, natom, nco_max, ncoa, ncob, ncoc, &
      ncod, nseta, nsetb, nsetc, nsetd, nsgfa_total, nsgfb_total, &
      nsgfc_total, nsgfd_total, nspins, p_index, pa, pa1, pa2, pb, pb1, pb2, &
      pc, pc1, pc2, pd, pd1, pd2, pgf_radius_box, sgfa, sgfb, sgfc, sgfd, stat, &
      my_original_cpu_id
    INTEGER(KIND=int_8)                      :: my_current_counter
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                last_sgf_global
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, ncgfa, ncgfb, ncgfc, ncgfd, npgfa, npgfb, &
      npgfc, npgfd, nsgfa, nsgfb, nsgfc, nsgfd
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfc, first_sgfd, &
                                                is_assoc_atomic_block
    LOGICAL                                  :: failure, first_j_loop, &
                                                first_k_loop, first_l_loop, &
                                                screen_pmat_forces
    REAL(dp)                                 :: Kbox, omega_box, Rbox
    REAL(KIND=dp) :: box_size(3), eps_schwarz, fac, hf_fraction, max_val1, &
      max_val2, pmax_max, ra(3), rb(3), rc(3), rd(3), screen_max, symm_fac, &
      temp1, temp1b, temp2, temp3, temp3b
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: primitive_forces
    REAL(KIND=dp), DIMENSION(:), POINTER     :: T1, T2
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: full_density, &
      full_density_beta, pmax_atom, screen_atom, sparse_block, sphi_a, &
      sphi_b, sphi_c, sphi_d, work_forces, zeta, zetb, zetc, zetd
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: pmax_set, screen_set
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: screen_periodic
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a, orb_basis_b, &
                                                orb_basis_c, orb_basis_d
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_block_node_type), POINTER      :: block_node
    TYPE(x_hf_potential_type)                :: potential_parameter
    TYPE(x_hf_screening_type)                :: screening_parameter
    TYPE(x_hf_general_type)                  :: general_parameter
    REAL(dp)                                 :: eps_box


    CALL mp_sync(para_env%group)
    CALL timeset(routineN,"I",'',handle)


    screening_parameter = qs_env%x_data%screening_parameter
    general_parameter = qs_env%x_data%general_parameter
    potential_parameter = qs_env%x_data%potential_parameter

    
    !! initalize force array to zero 
    CALL get_qs_env(qs_env=qs_env, force=force, error=error)
    DO ikind = 1,SIZE(force)
      force(ikind)%fock_4c(:,:) = 0.0_dp
    END DO

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    error=error)
 
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxsgf=maxsgf,&
                             maxnset=max_set,&
                             maxlgto=max_am)
   
    natom = SIZE(particle_set,1)

    !! One atom systems have no contribution to forces
    IF( natom == 1) THEN
      CALL timestop(0.0_dp,handle)
      RETURN
    END IF


    hf_fraction = general_parameter%fraction
    eps_schwarz = screening_parameter%eps_schwarz_forces
    eps_box = screening_parameter%eps_box_forces
    screen_pmat_forces =  screening_parameter%do_p_screening_forces

    !! Initialize memory for libint
    CALL initialize_libderiv(deriv, max_am, error)

    !! Get screening parameter

    !! Allocate work-space for forces
    ALLOCATE(work_forces(((max_am+1)*(max_am+2)/2)**4,12),STAT=stat)
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
      atom_kind => particle_set(iatom)%atomic_kind
 
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
                             nsgf=nsgfa_total)
      last_sgf_global(iatom) = last_sgf_global(iatom-1)+nsgfa_total
    END DO

    !! Get full density from all processes  
    ALLOCATE(full_density(1,1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL copy_local_sm_to_replicated_fm(rho%rho_ao(1)%matrix,full_density)
    CALL symmetrize_matrix(full_density,"upper_to_lower")
    CALL mp_sum(full_density,para_env%group)
    IF(nspins == 2) THEN
      ALLOCATE(full_density_beta(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL copy_local_sm_to_replicated_fm(rho%rho_ao(2)%matrix,full_density_beta)
      CALL symmetrize_matrix(full_density_beta,"upper_to_lower")
      CALL mp_sum(full_density_beta,para_env%group)
    ELSE
      NULLIFY(full_density_beta)
    END IF

    !! Calculate max entries for screening on actual density. If screen_p_mat_forces = FALSE, the
    !! matrix is initialized to 1.0
    CALL compute_pmax_matrices(pmax_set,pmax_atom,pmax_max, &
                  full_density,full_density_beta,screen_pmat_forces,particle_set,atomic_kind_set,error)


    !! Calculatematrices for screening purpose taking into account periodicity. The farfield part is already
    !! in memory
    max_k = qs_env%x_data%number_of_k_loops
    IF( do_periodic ) THEN
      box_size(1) = qs_env%cell%hmat(1,1)
      box_size(2) = qs_env%cell%hmat(2,2)
      box_size(3) = qs_env%cell%hmat(3,3)
      max_k = MIN(max_k, SIZE(qs_env%x_data%neighbor_cells))
      ALLOCATE(screen_periodic(natom,natom,max_set,max_set,SIZE(qs_env%x_data%neighbor_cells)))
      CALL comp_screen_mat_per(screen_set,screen_atom,screen_max,screen_periodic,particle_set,atomic_kind_set,&
                                               potential_parameter,.FALSE.,&
                                               neighbor_cells=qs_env%x_data%neighbor_cells,&
                                               box_size=box_size,&
                                               max_k=max_k,error=error)
      !! Farfield screening matrix already there
    ELSE
      CALL compute_screening_matrices(screen_set,screen_atom,screen_max,particle_set,atomic_kind_set,&
                                      potential_parameter,.FALSE.,error=error)
      !! Farfield screening matrix already there
    END IF

    !! create helper array that stores the sparsity of the KS-matrix
    ALLOCATE(is_assoc_atomic_block(natom,natom))
    is_assoc_atomic_block = 0
    DO iatom=1,ks_matrix(1)%matrix%nblock_row
      block_node => first_block_node(ks_matrix(1)%matrix,iatom)
      DO WHILE (ASSOCIATED(block_node))
        CALL get_block_node(block_node=block_node,&
                            block_col=jatom,&
                            BLOCK=sparse_block)

        is_assoc_atomic_block(iatom,jatom) = 1
        is_assoc_atomic_block(jatom,iatom) = 1
        block_node => next_block_node(block_node)
      END DO
    END DO

    CALL mp_sum(is_assoc_atomic_block,para_env%group)


    !! Load balance the work
    CALL load_balance_forces(qs_env%x_data,eps_schwarz,eps_box,particle_set,atomic_kind_set,para_env,&
                             potential_parameter,screen_set,screen_atom,screen_max,&
                             pmax_set,is_assoc_atomic_block,ks_matrix,do_periodic,error)

    !! precompute maximum nco and allocate scratch
    nco_max=0
    DO iatom=1,natom
      atom_kind => particle_set(iatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a, lmax=la_max, lmin=la_min,&
                             npgf=npgfa, nset=nseta, zet=zeta, nsgf_set=nsgfa,&
                             first_sgf=first_sgfa, sphi=sphi_a, nsgf=nsgfa_total)
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

    my_bin_size = SIZE(qs_env%x_data%distribution_forces)
    first_j_loop = .TRUE.
    first_l_loop = .TRUE.
    first_k_loop = .TRUE.
    !! outer loop will be cycled if actual work_load done
    outer: DO bin=1,my_bin_size
      !! Get information on workload
      my_istart = qs_env%x_data%distribution_forces(bin)%istart
      my_jstart = qs_env%x_data%distribution_forces(bin)%jstart
      my_kstart = qs_env%x_data%distribution_forces(bin)%kstart
      my_lstart = qs_env%x_data%distribution_forces(bin)%lstart
      my_original_cpu_id = qs_env%x_data%distribution_forces(bin)%original_cpu_id 
      my_current_counter = 0
      !! It could happen, that one bin is empty. In that case just cycle
      IF(my_istart == 0 .OR. my_jstart == 0 .OR. my_kstart == 0 .OR. my_lstart == 0 .OR. &
         qs_env%x_data%distribution_forces(bin)%number_of_integrals == 0) CYCLE
      DO iatom=my_istart,natom
        !! Be consistent to the modulo distribution of the load balance
        IF(MOD(iatom,para_env%num_pe) == my_original_cpu_id-1) THEN

          atom_kind => particle_set(iatom)%atomic_kind

          CALL get_atomic_kind(atomic_kind=atom_kind,&
                               orb_basis_set=orb_basis_a)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
                                 lmax=la_max,&
                                 lmin=la_min,&
                                 npgf=npgfa,&
                                 nset=nseta,&
                                 zet=zeta,&
                                 nsgf_set=nsgfa,&
                                 ncgf_set=ncgfa,&
                                 first_sgf=first_sgfa,&
                                 sphi=sphi_a,&
                                 nsgf=nsgfa_total)
          ra=particle_set(iatom)%r(:)
          i_atom = atom_of_kind(iatom)
          ikind = kind_of(iatom)
          pa1 = last_sgf_global(iatom) - nsgfa_total
          forces_map(1,1) = ikind
          forces_map(1,2) = i_atom
          IF(first_j_loop) jatom=my_jstart-1
          IF(.NOT. first_j_loop) jatom = iatom-1
          DO WHILE(jatom+1<=natom)
            jatom = jatom + 1
  
            IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
            atom_kind => particle_set(jatom)%atomic_kind
            CALL get_atomic_kind(atomic_kind=atom_kind,&
                                 orb_basis_set=orb_basis_b)
            CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                                   lmax=lb_max,&
                                   lmin=lb_min,&
                                   npgf=npgfb,&
                                   nset=nsetb,&
                                   zet=zetb,&
                                   nsgf_set=nsgfb,&
                                   ncgf_set=ncgfb,&
                                   first_sgf=first_sgfb,&
                                   sphi=sphi_b,&
                                   nsgf=nsgfb_total)
            rb=particle_set(jatom)%r(:)
            j_atom = atom_of_kind(jatom)
            jkind = kind_of(jatom)
            pb1 = last_sgf_global(jatom) - nsgfb_total
            forces_map(2,1) = jkind
            forces_map(2,2) = j_atom
            IF(first_k_loop) katom=my_kstart-1
              IF(.NOT. first_k_loop) katom = 0
              DO WHILE(katom+1<=natom)
                katom = katom + 1
              atom_kind => particle_set(katom)%atomic_kind
              CALL get_atomic_kind(atomic_kind=atom_kind,&
                                   orb_basis_set=orb_basis_c)
              CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                                     lmax=lc_max,&
                                     lmin=lc_min,&
                                     npgf=npgfc,&
                                     nset=nsetc,&
                                     zet=zetc,&
                                     nsgf_set=nsgfc,&
                                     ncgf_set=ncgfc,&
                                     first_sgf=first_sgfc,&
                                     sphi=sphi_c,&
                                     nsgf=nsgfc_total)
              rc=particle_set(katom)%r(:)
              k_atom = atom_of_kind(katom)
              kkind = kind_of(katom)
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
                IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
                IF(.NOT. (is_assoc_atomic_block(iatom,latom)>=1 .AND. &
                        is_assoc_atomic_block(iatom,katom)>=1 .AND. &
                        is_assoc_atomic_block(jatom,katom)>=1 .AND. &
                        is_assoc_atomic_block(jatom,latom)>=1 ) ) CYCLE
    
                IF(katom+latom<=iatom+jatom)  THEN
                  IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
    
                  my_current_counter = my_current_counter + 1
                     IF(my_current_counter > qs_env%x_data%distribution_forces(bin)%number_of_integrals) THEN
                       first_j_loop = .TRUE.
                       first_l_loop = .TRUE.
                       first_k_loop = .TRUE.
                       CYCLE outer
                     END IF
    
                  atom_kind => particle_set(latom)%atomic_kind
                  CALL get_atomic_kind(atomic_kind=atom_kind,&
                                       orb_basis_set=orb_basis_d)
                  CALL get_gto_basis_set(gto_basis_set=orb_basis_d,&
                                         lmax=ld_max,&
                                         lmin=ld_min,&
                                         npgf=npgfd,&
                                         nset=nsetd,&
                                         zet=zetd,&
                                         nsgf_set=nsgfd,&
                                         ncgf_set=ncgfd,&
                                         first_sgf=first_sgfd,&
                                         sphi=sphi_d,&
                                         nsgf=nsgfd_total)
 
                  rd=particle_set(latom)%r(:)
                
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
                      max_val1 = screen_set(iatom,jatom,iset,jset)
                      IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                      DO kset = 1,nsetc
                        ncoc = npgfc(kset)*ncoset(lc_max(kset))
                        sgfc = first_sgfc(1,kset)
                        DO lset = 1,nsetd
                          ncod = npgfd(lset)*ncoset(ld_max(lset)) 
                          sgfd = first_sgfd(1,lset)
     
                          max_val2 = max_val1 * screen_set(katom,latom,kset,lset) * &
                                     2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                                pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                          IF(max_val2<eps_schwarz)  CYCLE
                          box_idx1 = (iatom-1)*natom+jatom
                          box_idx2 = (iset-1)*max_set+jset
                          box_idx3 = (katom-1)*natom+latom
                          box_idx4 = (kset-1)*max_set+lset
                          Rbox = (qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%center(1) - &
                                  qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%center(1))**2 + &
                                 (qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%center(2) - &
                                  qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%center(2))**2 + &
                                 (qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%center(3) - &
                                  qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%center(3))**2
                          Rbox = SQRT(Rbox)
                          Rbox = Rbox - qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                        qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%R
                          pgf_radius_box = INT(1.0_dp/ &
                                (SQRT(2.0_dp*qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                          pgf_radius_box = pgf_radius_box + INT(1.0_dp/ &
                                (SQRT(2.0_dp*qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2
                          IF(Rbox>pgf_radius_box) THEN
                            max_val2 = 0.0_dp
                            SELECT CASE(potential_parameter%potential_type)
                              CASE(do_hf_potential_coulomb)
                                Kbox = qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = 1.0_dp/Rbox*Kbox* &
                                           2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                           pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                             CASE(do_hf_potential_short)
                                omega_box = 1.0_dp/qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                            1.0_dp/qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                            1.0_dp/(potential_parameter%omega**2)
                                omega_box = 1.0_dp/omega_box
                                Kbox = qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox * &
                                           2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                           pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                            END SELECT
                            IF(max_val2<eps_box) CYCLE
                          END IF
                          IF( .NOT. do_periodic ) THEN
                            CALL forces4(ra, rb, rc, rd,&
                                         zeta(:,iset), zetb(:,jset), zetc(:,kset), zetd(:,lset),&
                                         la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                         lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                         npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                         ncoa, ncob, ncoc, ncod,&
                                         primitive_forces,&
                                         work_forces,&
                                         eps_schwarz,&
                                         potential_parameter)
                          ELSE
                            CALL forces4_periodic(ra, rb, rc, rd,&
                                                  zeta(:,iset), zetb(:,jset), zetc(:,kset), zetd(:,lset),&
                                                  la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                                  lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                                  npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                                  ncoa, ncob, ncoc, ncod,&
                                                  primitive_forces,&
                                                  work_forces,&
                                                  potential_parameter, box_size, &
                                                  qs_env%x_data%neighbor_cells, screen_periodic(iatom,jatom,iset,jset,:),&
                                                  screen_periodic(katom,latom,kset,lset,:), eps_schwarz, max_k)
                          END IF
   
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
                                    force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1, &
                                          forces_map((coord-1)/3 + 1,2)) = &
                                    force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1, &
                                          forces_map((coord-1)/3 + 1,2)) - &
                                    temp2 * T2(p_index) !primitive_forces((coord-1)*ncoa*ncob*ncoc*ncod+p_index)
                                    IF(nspins == 2) THEN
                                      temp2 = temp1b*full_density_beta(pa,pd) + temp3b*full_density_beta(pa,pc)
                                      force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1,&
                                            forces_map((coord-1)/3 + 1,2)) = &
                                      force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1,&
                                            forces_map((coord-1)/3 + 1,2)) - &
                                      temp2 * T2(p_index)!primitive_forces((coord-1)*ncoa*ncob*ncoc*ncod+p_index)
                                    END IF 
                                  END DO
                                END DO !ma
                              END DO !mb
                            END DO !mc
                          END DO !md
      
                        END DO !lset
                      END DO !ksetm
                    END DO !jset
                  END DO !iset
                END IF !cycle
              END DO !latom
              first_l_loop = .FALSE.
            END DO !katom
            first_k_loop = .FALSE.
          END DO !jatom
          first_j_loop = .FALSE.
        END IF
      END DO !iatom
      !! Because of the Monte Carlo optimization it is not guaranteed that there are no more bunches to calculate.
      !! Therefore we have to reset the flags
      first_j_loop = .TRUE.
      first_k_loop = .TRUE.
      first_l_loop = .TRUE.
    END DO outer !bin
   
    IF( do_periodic ) THEN
      DEALLOCATE(screen_periodic, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(T1,primitive_forces, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(last_sgf_global,screen_set, screen_atom, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(full_density,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(is_assoc_atomic_block,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(nspins==2) THEN
      DEALLOCATE(full_density_beta,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(work_forces, atom_of_kind,kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(pmax_set,pmax_atom,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    !
    ! this wraps to free_libint, but currently causes a segfault
    ! as a result, we don't call it, but some memory remains allocated
    !  
    ! CALL terminate_libderiv(deriv)
    
    CALL mp_sync(para_env%group)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE derivatives_four_center

!!***** four_center_integrals/coulomb4_periodic [1.0] *
!!
!!   NAME
!!     coulomb4_periodic
!!
!!   FUNCTION
!!     calculates two-electron integrals of a quartet/shell using the library 
!!     lib_int in the case of periodic boundary conditions
!!
!!   NOTES
!!     The minimum image convention (mic) is applied in the following way:
!!       1. Fix A and shift B according to mic
!!       2. Fix C and shift D according to mic
!!       3. Calculate interaction centers P and Q
!!       4. Fix P (i.e. A and B) and shift Q (i.e. C and D) according to mic
!!
!!   INPUTS
!!     - r: position
!!     - z: zeta
!!     - l: angular momentum
!!     - npgf: number of primitive cartesian gaussian in actual shell
!!     - primitive_integrals: array of primitive_integrals
!!     - box_size, neighbor_cells: contains info for pbc
!!     - potential_parameter: contains info for libint 
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     06.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE coulomb4_periodic(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                               la_min, la_max, lb_min, lb_max,&
                               lc_min, lc_max, ld_min, ld_max,&
                               npgfa, npgfb, npgfc, npgfd,&
                               ncoa,ncob,ncoc,ncod,&
                               primitive_integrals,&
                               potential_parameter, box_size, neighbor_cells, &
                               screen1, screen2, eps_schwarz, max_k)
  
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
    TYPE(x_hf_potential_type)                :: potential_parameter
    REAL(dp), INTENT(IN)                     :: box_size(3)
    TYPE(x_hf_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(dp), DIMENSION(:)                   :: screen1, screen2
    REAL(dp), INTENT(IN)                     :: eps_schwarz
    INTEGER, INTENT(IN)                      :: max_k

    INTEGER :: i, ipgf, j, jpgf, k, kpgf, la, lb, lc, ld, lpgf, offset_a, &
      offset_a1, offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1
    REAL(dp) :: B11(3), C11(3), D11(3), Eta, im_B(3), im_D(3), P(3), &
      pbc_B(3), pbc_D(3), PQ(3), Q(3), shift(3), shift1(3), temp(3), &
      tmp_D(3), Zeta1, Zeta_A, Zeta_B, Zeta_C, Zeta_D

    primitive_integrals = 0.0_dp
    temp = rb -ra
    CALL pbc(temp, box_size, pbc_B)
    B11 = ra + pbc_B
    temp = rd-rc
    CALL pbc(temp, box_size, pbc_D)
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
                      CALL pbc(temp,box_size,PQ)
                      shift = -PQ + temp
                      DO k=1,max_k!size(neighbor_cells)
                        shift1(1) =  neighbor_cells(k)%cell(1) * box_size(1)
                        shift1(2) =  neighbor_cells(k)%cell(2) * box_size(2)
                        shift1(3) =  neighbor_cells(k)%cell(3) * box_size(3)
                        C11 = rc + shift + shift1
                        tmp_D = im_D + shift + shift1
                        DO la = la_min,la_max
                          offset_a = offset_a1 + ncoset(la-1)
                          !Build primitives
                          CALL evaluate_eri(lib, ra, im_B, C11, tmp_D,&
                                            zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                            la, lb, lc ,ld,&
                                            ncoa,ncob,ncoc,ncod,&
                                            offset_a,offset_b,offset_c,offset_d, &
                                            primitive_integrals, &
                                            potential_parameter)
                        END DO !la
                      END DO
                    END DO !ipgf
                  END DO !lb
                END DO !jpgf
              END DO !lc
            END DO !kpgf
          END DO !ld
        END DO !lpgf
      END DO
    END DO
  END SUBROUTINE coulomb4_periodic
!***************************************************************************

!!***** four_center_integrals/comp_screen_mat_per [1.0] *
!!
!!   NAME
!!     comp_screen_mat_per
!!
!!   FUNCTION
!!     calculates the screening matrices in the case of periodic boundary
!!     conditions
!!
!!   NOTES
!!     - already computes the square root
!!     - also takes the primitive to sgf transformation into account
!!     - in this way all stored contracted integrals are smaller or equal to their estimate
!!
!!   INPUTS
!!     - screen_xxx: screening matrices
!!     - read_from_memory: Flag to control reading (when true, data come from buffer)
!!     - buffer: stores the matrices
!!     - neighbor_cells, box_size: contain information fort pbc
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!    06.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE comp_screen_mat_per(screen_set, screen_atom, screen_max,screen_periodic,particle_set,&
                                 atomic_kind_set,potential_parameter, read_from_memory, buffer,&
                                 neighbor_cells, box_size, max_k, error)
    REAL(kind=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: screen_set
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: screen_atom
    REAL(KIND=dp)                            :: screen_max
    REAL(kind=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: screen_periodic
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(x_hf_potential_type)                :: potential_parameter
    LOGICAL                                  :: read_from_memory
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), OPTIONAL, &
      POINTER                                :: buffer
    TYPE(x_hf_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(KIND=dp)                            :: box_size(3)
    INTEGER, INTENT(IN)                      :: max_k
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'comp_screen_mat_per', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, jatom, jset, k, &
                                                katom, kset, max_set, natom, &
                                                ncob, nsetb, nsetc, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: lb_max, lb_min, lc_max, &
                                                lc_min, npgfb, npgfc, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: max_val
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: contraction_max
    REAL(dp), DIMENSION(:, :), POINTER       :: sphi_b, zetb, zetc
    REAL(KIND=dp)                            :: p2s_b, p2s_c, rb(3), rc(3)
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_b, orb_basis_c

    CALL timeset(routineN,handle)

    failure = .FALSE.

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set)
    ALLOCATE(screen_set(natom,natom,max_set,max_set),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(screen_atom(natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(contraction_max(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    screen_set=0.0_dp
    screen_atom=0.0_dp
    contraction_max=0.0_dp
    screen_periodic = 0.0_dp
    ALLOCATE(max_val(SIZE(neighbor_cells)))
    max_val = 0.0_dp
    
    DO jatom=1,natom

      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             first_sgf=first_sgfb,&
                             sphi=sphi_b,&
                             nsgf_set=nsgfb)
      DO jset = 1,nsetb
        ! takes the primitive to contracted transformation into account
        ncob = npgfb(jset)*ncoset(lb_max(jset))
        sgfb = first_sgfb(1,jset)
        ! if the primitives are assumed to be all of max_val2, max_val2*p2s_b becomes
        ! the maximum value after multiplication with sphi_b
        contraction_max(jset,jatom) = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,k))),k=sgfb,sgfb+nsgfb(jset)-1)/))
      ENDDO
    ENDDO

    DO jatom=1,natom

      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             lmin=lb_min,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             zet=zetb)
      rb=particle_set(jatom)%r(:)
      DO katom = 1,natom
        atom_kind => particle_set(katom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_c)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                               lmax=lc_max,&
                               lmin=lc_min,&
                               npgf=npgfc,&
                               nset=nsetc,&
                               zet=zetc)
        rc=particle_set(katom)%r(:)
        DO jset = 1,nsetb
          p2s_b = contraction_max(jset,jatom)
          DO kset = 1,nsetc
            p2s_c = contraction_max(kset,katom)
            IF( read_from_memory ) THEN
              screen_periodic(jatom,katom,jset,kset,:) = buffer(jatom,katom,jset,kset,:)
            ELSE
              max_val = 0.0_dp
              CALL screen4_periodic(rb, rc, &
                                    zetb(:,jset), zetc(:,kset), &
                                    lb_min(jset), lb_max(jset), lc_min(kset), lc_max(kset),&
                                    npgfb(jset), npgfc(kset), &
                                    max_val, potential_parameter, neighbor_cells, box_size, max_k)
                screen_periodic(jatom,katom,jset,kset,:) = SQRT(max_val)*p2s_b*p2s_c
                IF( PRESENT(buffer) ) buffer(jatom,katom,jset,kset,:) = screen_periodic(jatom,katom,jset,kset,:)
              END IF
              screen_set(jatom,katom,jset,kset) = MAXVAL(screen_periodic(jatom,katom,jset,kset,:))
          END DO
        END DO
        screen_atom(jatom,katom)=MAXVAL(screen_set(jatom,katom,1:nsetb,1:nsetc))
      END DO
    END DO

    screen_max=MAXVAL(screen_atom)
    DEALLOCATE(contraction_max,max_val, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)

  END SUBROUTINE comp_screen_mat_per
!***************************************************************************

!!***** four_center_integrals/screen4_periodic [1.0] *
!!
!!   NAME
!!     screen4_periodic
!!
!!   FUNCTION
!!     calculates max values of two-electron integrals in a quartet/shell
!!     w.r.t. different zetas using the library lib_int in the case of 
!!     periodic boundary conditions
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - r: position
!!     - z: zeta
!!     - l: angular momentum
!!     - npgf: number of primitive cartesian gaussian in actual shell
!!     - max_val: schwarz screening value
!!     - potential_parameter: contains info for libint
!!     - neighbor_cells: Array containing the considered cells
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!    06.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE screen4_periodic(ra, rb, zeta, zetb, &
                             la_min, la_max, lb_min, lb_max,&
                             npgfa, npgfb, &
                             max_val, potential_parameter, neighbor_cells, &
                             box_size, max_k)

    REAL(dp), INTENT(IN)                     :: ra(3), rb(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb
    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, &
                                                lb_max, npgfa, npgfb
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: max_val
    TYPE(x_hf_potential_type)                :: potential_parameter
    TYPE(x_hf_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(KIND=dp), INTENT(IN)                :: box_size(3)
    INTEGER, INTENT(IN)                      :: max_k

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'screen4_periodic'

    INTEGER                                  :: i, ipgf, jpgf, k, la, lb
    REAL(dp)                                 :: B11(3), im_B(3), &
                                                max_val_temp, pbc_B(3), &
                                                shift1(3), temp(3)

    temp = rb -ra
    CALL pbc(temp, box_size, pbc_B)
    B11 = ra + pbc_B
    DO i = 1,SIZE(neighbor_cells)
      im_B(1) = B11(1) + neighbor_cells(i)%cell(1) * box_size(1)
      im_B(2) = B11(2) + neighbor_cells(i)%cell(2) * box_size(2)
      im_B(3) = B11(3) + neighbor_cells(i)%cell(3) * box_size(3)
      DO k=1,max_k!size(neighbor_cells)
        shift1(1) =  neighbor_cells(k)%cell(1) * box_size(1)
        shift1(2) =  neighbor_cells(k)%cell(2) * box_size(2)
        shift1(3) =  neighbor_cells(k)%cell(3) * box_size(3)
        DO ipgf = 1,npgfa
          DO jpgf = 1,npgfb
            DO la = la_min,la_max
              DO lb = lb_min, lb_max
                !Build primitives
                 max_val_temp = 0.0_dp
                 CALL evaluate_eri_screen(lib, ra, im_B, ra, im_B,&
                                          zeta(ipgf), zetb(jpgf), zeta(ipgf), zetb(jpgf),&
                                          la, lb, la ,lb,&
                                           max_val_temp, potential_parameter)
                 max_val(i) = MAX(max_val(i), max_val_temp)
              END DO !lb
            END DO !la
          END DO !jpgf
        END DO !ipgf
      END DO
    END DO
  END SUBROUTINE screen4_periodic
!***************************************************************************

!!***** four_center_integrals/forces4_periodic [1.0] *
!!
!!   NAME
!!     forces4_periodic
!!
!!   FUNCTION
!!     calculates two-electron derivatived of a quartet/shell using the library 
!!     lib_deriv wrt pbc
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - r: position
!!     - z: zeta
!!     - l: angular momentum
!!     - npgf: number of primitive cartesian gaussian in actual shell
!!     - primitive_forces: array of primitive_forces
!!     - work_forces: work_storage
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     07.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE forces4_periodic(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                              la_min, la_max, lb_min, lb_max,&
                              lc_min, lc_max, ld_min, ld_max,&
                              npgfa, npgfb, npgfc, npgfd,&
                              ncoa, ncob, ncoc, ncod,&
                              primitive_forces, work_forces,&
                              potential_parameter, box_size, neighbor_cells, screen1, screen2, eps_schwarz, max_k)
  
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb, zetc, zetd
    INTEGER, INTENT(IN) :: la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
      ld_min, ld_max, npgfa, npgfb, npgfc, npgfd, ncoa, ncob, ncoc, ncod
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod, 12)  :: primitive_forces
    REAL(dp), DIMENSION(:, :), POINTER       :: work_forces
    TYPE(x_hf_potential_type)                :: potential_parameter
    REAL(KIND=dp), INTENT(IN)                :: box_size(3)
    TYPE(x_hf_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(dp), DIMENSION(:)                   :: screen1, screen2
    REAL(dp), INTENT(IN)                     :: eps_schwarz
    INTEGER, INTENT(IN)                      :: max_k

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'forces4_periodic'

    INTEGER :: i, ipgf, j, jpgf, k, kpgf, la, lb, lc, ld, lpgf, offset_a, &
      offset_a1, offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1
    REAL(dp) :: B11(3), C11(3), D11(3), Eta, im_B(3), im_D(3), P(3), &
      pbc_B(3), pbc_D(3), PQ(3), Q(3), shift(3), shift1(3), temp(3), &
      tmp_D(3), Zeta1, Zeta_A, Zeta_B, Zeta_C, Zeta_D

    primitive_forces = 0.0_dp
    temp = rb -ra
    CALL pbc(temp, box_size, pbc_B)
    B11 = ra + pbc_B
    temp = rd-rc
    CALL pbc(temp, box_size, pbc_D)
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
                      CALL pbc(temp,box_size,PQ)
                      shift = -PQ + temp
                      DO k=1,max_k!size(neighbor_cells)
                        shift1(1) =  neighbor_cells(k)%cell(1) * box_size(1)
                        shift1(2) =  neighbor_cells(k)%cell(2) * box_size(2)
                        shift1(3) =  neighbor_cells(k)%cell(3) * box_size(3)
                     
                        C11 = rc + shift + shift1
                        tmp_D = im_D + shift + shift1
                        DO la = la_min,la_max
                          offset_a = offset_a1 + ncoset(la-1)
                          !Build primitives
                          CALL evaluate_deriv_eri(deriv, ra, im_B, C11, tmp_D,&
                                                  zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                                  la, lb, lc ,ld,&
                                                  work_forces,&
                                                  ncoa, ncob, ncoc, ncod, &
                                                  primitive_forces, &
                                                  offset_a, offset_b, offset_c, offset_d,&
                                                  potential_parameter)
                        END DO !la
                      END DO
                    END DO !ipgf
                  END DO !lb
                END DO !jpgf
              END DO !lc
            END DO !kpgf
          END DO !ld
        END DO !lpgf
      END DO
    END DO
  END SUBROUTINE forces4_periodic
!***************************************************************************

!!***** four_center_integrals/screen_farfield [1.0] *
!!
!!   NAME
!!     screen_farfield
!!
!!   FUNCTION
!!     Estimates the contribution of an ERI taking 1/R dependence into account.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - boxes: Buffer to store screening matrices
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     08.2007 created [Manuel Guidon]
!!
!!****************************************************************************
 
  SUBROUTINE screen_farfield(boxes, particle_set, atomic_kind_set, potential_parameter, error)
    TYPE(x_hf_box_scheme), DIMENSION(:, :), &
      POINTER                                :: boxes
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(x_hf_potential_type)                :: potential_parameter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'screen_farfield', &
      routineP = moduleN//':'//routineN

    INTEGER :: box_idx1, box_idx2, handle, jatom, jpgf, jset, k, katom, kpgf, &
      kset, max_set, natom, ncob, nsetb, nsetc, pgf_radius, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: lb_max, lb_min, lc_max, &
                                                lc_min, npgfb, npgfc, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: failure
    REAL(dp)                                 :: px_max, px_min, py_max, &
                                                py_min, pz_max, pz_min, rbc2
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: contraction_max
    REAL(dp), DIMENSION(:, :), POINTER       :: rpgfb, rpgfc, sphi_b, zetb, &
                                                zetc
    REAL(KIND=dp) :: alpha, K_max, Kbc, p2s_b, p2s_c, R_max, rb(3), Rbc, &
      rc(3), rp(3), system_size(3), x_max, x_min, y_max, y_min, z_max, z_min, &
      zet1, zet2
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_b, orb_basis_c

    CALL timeset(routineN,handle)

    failure = .FALSE.

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set)
   
    ALLOCATE(contraction_max(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    contraction_max=0.0_dp

    x_min = HUGE(x_min)
    x_max =-HUGE(x_max)
    y_min = HUGE(y_min)
    y_max = -HUGE(y_max)
    z_min = HUGE(z_min)
    z_max = -HUGE(z_max)
    DO jatom=1,natom
      rb=particle_set(jatom)%r(:)
      IF(rb(1)>x_max) x_max=rb(1)
      IF(rb(1)<x_min) x_min=rb(1)
      IF(rb(2)>y_max) y_max=rb(2)
      IF(rb(2)<y_min) y_min=rb(2)
      IF(rb(3)>z_max) z_max=rb(3)
      IF(rb(3)<z_min) z_min=rb(3)
    END DO

    system_size(1) = x_max - x_min
    system_size(2) = y_max - y_min
    system_size(3) = z_max - z_min
    DEALLOCATE(boxes,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(boxes(natom**2,max_set**2), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    boxes(:,:)%center(1) = 0.0_dp
    boxes(:,:)%center(2) = 0.0_dp
    boxes(:,:)%center(3) = 0.0_dp
    boxes(:,:)%K = -HUGE(Kbc)
    boxes(:,:)%R = -HUGE(Rbc)
    boxes(:,:)%alpha = HUGE(alpha)
    DO jatom=1,natom
      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             first_sgf=first_sgfb,&
                             sphi=sphi_b,&
                             nsgf_set=nsgfb)
      DO jset = 1,nsetb
        ! takes the primitive to contracted transformation into account
        ncob = npgfb(jset)*ncoset(lb_max(jset))
        sgfb = first_sgfb(1,jset)
        ! if the primitives are assumed to be all of max_val2, max_val2*p2s_b becomes
        ! the maximum value after multiplication with sphi_b
        contraction_max(jset,jatom) = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,k))),k=sgfb,sgfb+nsgfb(jset)-1)/))
      ENDDO
    ENDDO

    DO jatom=1,natom
      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             lmin=lb_min,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             pgf_radius=rpgfb,&
                             zet=zetb)
      rb=particle_set(jatom)%r(:)
      DO katom = 1,natom
        atom_kind => particle_set(katom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_c)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                               lmax=lc_max,&
                               lmin=lc_min,&
                               npgf=npgfc,&
                               nset=nsetc,&
                               pgf_radius=rpgfc,&
                               zet=zetc)
        rc=particle_set(katom)%r(:)
        rbc2 = (rb(1)-rc(1))**2+(rb(2)-rc(2))**2+(rb(3)-rc(3))**2
        box_idx1 = (jatom-1)*natom + katom 
        DO jset = 1,nsetb
          p2s_b = contraction_max(jset,jatom)
          DO kset = 1,nsetc
            p2s_c = contraction_max(kset,katom)
            K_max = TINY(K_max)
            R_max = TINY(R_max)
            px_min = HUGE(px_min)
            px_max = -HUGE(px_max)
            py_min = HUGE(py_min)
            py_max = -HUGE(py_max)
            pz_min = HUGE(pz_min)
            pz_max = -HUGE(pz_max)
            box_idx2 = (jset-1)* max_set + kset
            DO jpgf = 1,npgfb(jset)
              zet1 = zetb(jpgf,jset)
              DO kpgf = 1,npgfc(kset)
                zet2 = zetc(kpgf,kset)
                alpha = zet1+zet2
                rp = (zet1*rb+zet2*rc)/alpha
                IF(rp(1)<px_min) px_min = rp(1)
                IF(rp(2)<py_min) py_min = rp(2)
                IF(rp(3)<pz_min) pz_min = rp(3)
                IF(rp(1)>px_max) px_max = rp(1)
                IF(rp(2)>py_max) py_max = rp(2)
                IF(rp(3)>pz_max) pz_max = rp(3)
                Kbc = 1.0_dp/alpha*EXP(-zet1*zet2/alpha*(rbc2))
                IF( alpha < boxes(box_idx1, box_idx2)%alpha ) boxes(box_idx1, box_idx2)%alpha = alpha
                IF( Kbc > boxes(box_idx1, box_idx2)%K ) boxes(box_idx1, box_idx2)%K = Kbc*p2s_b*p2s_c*&
                    pi**1.25_dp*SQRT(2.0_dp)
              END DO
            END DO
            boxes(box_idx1, box_idx2)%center(1) = (px_max + px_min) / 2.0_dp
            boxes(box_idx1, box_idx2)%center(2) = (py_max + py_min) / 2.0_dp
            boxes(box_idx1, box_idx2)%center(3) = (pz_max + pz_min) / 2.0_dp
            R_max = TINY(R_max)
            DO jpgf = 1,npgfb(jset)
              zet1 = zetb(jpgf,jset)
              DO kpgf = 1,npgfc(kset)
                zet2 = zetc(kpgf,kset)
                alpha = zet1+zet2
                rp = (zet1*rb+zet2*rc)/alpha
                R_max = SQRT((rp(1)-boxes(box_idx1, box_idx2)%center(1))**2 + &
                             (rp(2)-boxes(box_idx1, box_idx2)%center(2))**2 + &
                             (rp(3)-boxes(box_idx1, box_idx2)%center(3))**2)
                IF( R_max > boxes(box_idx1, box_idx2)%R ) boxes(box_idx1, box_idx2)%R = R_max
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(contraction_max, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
  END SUBROUTINE screen_farfield
!***************************************************************************

!!***** four_center_integrals/screen_farfield_per [1.0] *
!!
!!   NAME
!!     screen_farfield_per
!!
!!   FUNCTION
!!     Estimates the contribution of an ERI taking 1/R dependence into account.
!!     Periodic case
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - boxes: Buffer to store screening matrices
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************
 
  SUBROUTINE screen_farfield_per(boxes, particle_set, atomic_kind_set, potential_parameter, &
                                 neighbor_cells, box_size, error)
    TYPE(x_hf_box_scheme), DIMENSION(:, :), &
      POINTER                                :: boxes
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(x_hf_potential_type)                :: potential_parameter
    TYPE(x_hf_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(KIND=dp), INTENT(IN)                :: box_size(3)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'screen_farfield_per', &
      routineP = moduleN//':'//routineN

    INTEGER :: box_idx1, box_idx2, handle, i, jatom, jpgf, jset, k, katom, &
      kpgf, kset, max_set, natom, ncob, nsetb, nsetc, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: lb_max, lb_min, lc_max, &
                                                lc_min, npgfb, npgfc, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: failure
    REAL(dp)                                 :: C11(3), im_C(3), pbc_C(3), &
                                                px_max, px_min, py_max, &
                                                py_min, pz_max, pz_min, rbc2, &
                                                temp(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: contraction_max
    REAL(dp), DIMENSION(:, :), POINTER       :: rpgfb, rpgfc, sphi_b, zetb, &
                                                zetc
    REAL(KIND=dp)                            :: alpha, K_max, Kbc, p2s_b, &
                                                p2s_c, R_max, rb(3), Rbc, &
                                                rc(3), rp(3), zet1, zet2
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_b, orb_basis_c

    CALL timeset(routineN,handle)

    failure = .FALSE.

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set)
   
    ALLOCATE(contraction_max(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    contraction_max=0.0_dp

    DEALLOCATE(boxes,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(boxes(natom**2,max_set**2), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    boxes(:,:)%center(1) = 0.0_dp
    boxes(:,:)%center(2) = 0.0_dp
    boxes(:,:)%center(3) = 0.0_dp
    boxes(:,:)%K = -HUGE(Kbc)
    boxes(:,:)%R = -HUGE(Rbc)
    boxes(:,:)%alpha = HUGE(alpha)
    DO jatom=1,natom
      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             first_sgf=first_sgfb,&
                             sphi=sphi_b,&
                             nsgf_set=nsgfb)
      DO jset = 1,nsetb
        ! takes the primitive to contracted transformation into account
        ncob = npgfb(jset)*ncoset(lb_max(jset))
        sgfb = first_sgfb(1,jset)
        ! if the primitives are assumed to be all of max_val2, max_val2*p2s_b becomes
        ! the maximum value after multiplication with sphi_b
        contraction_max(jset,jatom) = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,k))),k=sgfb,sgfb+nsgfb(jset)-1)/))
      ENDDO
    ENDDO

    DO jatom=1,natom
      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             lmin=lb_min,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             pgf_radius=rpgfb,&
                             zet=zetb)
      rb=particle_set(jatom)%r(:)
      DO katom = 1,natom
        atom_kind => particle_set(katom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_c)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                               lmax=lc_max,&
                               lmin=lc_min,&
                               npgf=npgfc,&
                               nset=nsetc,&
                               pgf_radius=rpgfc,&
                               zet=zetc)
        rc=particle_set(katom)%r(:)
        box_idx1 = (jatom-1)*natom + katom 
        temp = rc-rb
        CALL pbc(temp, box_size, pbc_C)
        C11 = rb + pbc_C
        DO jset = 1,nsetb
          p2s_b = contraction_max(jset,jatom)
          DO kset = 1,nsetc
            p2s_c = contraction_max(kset,katom)
            K_max = TINY(K_max)
            R_max = TINY(R_max)
            px_min = HUGE(px_min)
            px_max = -HUGE(px_max)
            py_min = HUGE(py_min)
            py_max = -HUGE(py_max)
            pz_min = HUGE(pz_min)
            pz_max = -HUGE(pz_max)
            box_idx2 = (jset-1)* max_set + kset
            DO i = 1,SIZE(neighbor_cells)
              im_C(1) = C11(1) + neighbor_cells(i)%cell(1) * box_size(1)
              im_C(2) = C11(2) + neighbor_cells(i)%cell(2) * box_size(2)
              im_C(3) = C11(3) + neighbor_cells(i)%cell(3) * box_size(3)
              rbc2 = (rb(1)-im_C(1))**2+(rb(2)-im_C(2))**2+(rb(3)-im_C(3))**2
              DO jpgf = 1,npgfb(jset)
                zet1 = zetb(jpgf,jset)
                DO kpgf = 1,npgfc(kset)
                  zet2 = zetc(kpgf,kset)
                  alpha = zet1+zet2
                  rp = (zet1*rb+zet2*im_C)/alpha
                  IF(rp(1)<px_min) px_min = rp(1)
                  IF(rp(2)<py_min) py_min = rp(2)
                  IF(rp(3)<pz_min) pz_min = rp(3)
                  IF(rp(1)>px_max) px_max = rp(1)
                  IF(rp(2)>py_max) py_max = rp(2)
                  IF(rp(3)>pz_max) pz_max = rp(3)
                  Kbc = 1.0_dp/alpha*EXP(-zet1*zet2/alpha*(rbc2))
                  IF( alpha < boxes(box_idx1, box_idx2)%alpha ) boxes(box_idx1, box_idx2)%alpha = alpha
                  IF( Kbc > boxes(box_idx1, box_idx2)%K ) boxes(box_idx1, box_idx2)%K = Kbc*p2s_b*p2s_c*&
                    pi**1.25_dp*SQRT(2.0_dp)
                END DO
              END DO
            END DO
            boxes(box_idx1, box_idx2)%center(1) = (px_max + px_min) / 2.0_dp
            boxes(box_idx1, box_idx2)%center(2) = (py_max + py_min) / 2.0_dp
            boxes(box_idx1, box_idx2)%center(3) = (pz_max + pz_min) / 2.0_dp
            R_max = TINY(R_max)
            DO i = 1,SIZE(neighbor_cells)
              im_C(1) = C11(1) + neighbor_cells(i)%cell(1) * box_size(1)
              im_C(2) = C11(2) + neighbor_cells(i)%cell(2) * box_size(2)
              im_C(3) = C11(3) + neighbor_cells(i)%cell(3) * box_size(3)
              rbc2 = (rb(1)-im_C(1))**2+(rb(2)-im_C(2))**2+(rb(3)-im_C(3))**2
              DO jpgf = 1,npgfb(jset)
                zet1 = zetb(jpgf,jset)
                DO kpgf = 1,npgfc(kset)
                  zet2 = zetc(kpgf,kset)
                  alpha = zet1+zet2
                  rp = (zet1*rb+zet2*im_C)/alpha
                  R_max = SQRT((rp(1)-boxes(box_idx1, box_idx2)%center(1))**2 + &
                               (rp(2)-boxes(box_idx1, box_idx2)%center(2))**2 + &
                               (rp(3)-boxes(box_idx1, box_idx2)%center(3))**2)
                  IF( R_max > boxes(box_idx1, box_idx2)%R ) boxes(box_idx1, box_idx2)%R = R_max
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(contraction_max, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
  END SUBROUTINE screen_farfield_per
!***************************************************************************

!!***** four_center_integrals/update_pmax_matrix [1.0] *
!!
!!   NAME
!!     update_pmax_matrix
!!
!!   FUNCTION
!!     updates the maximum of the density matrix in compressed form for screening purposes 
!!
!!   NOTES
!!     - updates for each pair of shells the maximum absolute value of p
!!
!!   INPUTS
!!    -  pmax_set: buffer to store matrix
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!    09.2007 created [Manuel Guidon]
!!
!****************************************************************************
  
  SUBROUTINE update_pmax_matrix(pmax_set, &
              full_density,full_density_beta,particle_set,atomic_kind_set,error)

    REAL(kind=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: pmax_set
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: full_density, &
                                                full_density_beta
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_pmax_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, jatom, jset, katom, kset, max_set, mb, mc, natom, &
      nsetb, nsetc, nsgfb_total, nsgfc_total, pb, pb1, pb2, pc, pc1, pc2, &
      sgfb, sgfc, stat
    INTEGER, DIMENSION(:), POINTER           :: last_sgf_global, nsgfb, nsgfc
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb, first_sgfc
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: pmax_tmp
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_b, orb_basis_c

    CALL timeset(routineN,handle)
    failure = .FALSE.

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set)

    ! create a mapping array of easy access to the density matrix elements
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO jatom=1,natom
      atom_kind => particle_set(jatom)%atomic_kind
 
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             nsgf=nsgfb_total)
      last_sgf_global(jatom) = last_sgf_global(jatom-1)+nsgfb_total
    END DO

    pmax_set=0.0_dp

    DO jatom=1,natom

      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             nsgf=nsgfb_total,&
                             nset=nsetb,& 
                             first_sgf=first_sgfb,&
                             nsgf_set=nsgfb)
      pb1 = last_sgf_global(jatom) - nsgfb_total

      DO katom = 1,natom
        atom_kind => particle_set(katom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_c)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                               nsgf=nsgfc_total,&
                               nset=nsetc,&
                               first_sgf=first_sgfc,&
                               nsgf_set=nsgfc)
        pc1 = last_sgf_global(katom) - nsgfc_total

        DO jset = 1,nsetb
          sgfb = first_sgfb(1,jset)
          DO kset = 1,nsetc
            sgfc = first_sgfc(1,kset)
            pmax_tmp = 0.0_dp
            DO mc = 1,nsgfc(kset)
               pc2 = mc+sgfc-1
               pc = pc1 + pc2
               DO mb = 1,nsgfb(jset)
                  pb2 = mb+sgfb-1
                  pb = pb1 + pb2
                  pmax_tmp=MAX(pmax_tmp,ABS(full_density(pb,pc)))
                  IF( ASSOCIATED(full_density_beta) ) pmax_tmp=MAX(pmax_tmp,ABS(full_density_beta(pb,pc)))
               ENDDO
            ENDDO
            pmax_set(jset,kset,jatom,katom) = pmax_tmp
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE(last_sgf_global,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE update_pmax_matrix
!***************************************************************************

!!***** four_center_integrals/on_the_fly [1.0] *
!!
!!   NAME
!!     on_the_fly
!!
!!   FUNCTION
!!     This routine does the correction to the screening on an in inital
!!     density matrix. It is basically a copy from the out-core part in
!!     integrate_four_center.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     -

!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE on_the_fly(do_periodic, ra, rb, rc, rd, zeta, zetb, zetc, zetd, &
                        la_min, la_max, lb_min, lb_max, lc_min, lc_max, ld_min, ld_max, &
                        npgfa, npgfb, npgfc, npgfd, ncoa, ncob, ncoc, ncod, &
                        nsgfa, nsgfb, nsgfc, nsgfd, sgfa, sgfb, sgfc, sgfd, &
                        primitive_integrals, potential_parameter, box_size, qs_env, screen_periodic, &
                        eps_schwarz, max_k, sphi_a, sphi_b, sphi_c, sphi_d, T1, & 
                        full_ks, full_ks_beta, full_density, full_density_beta, fac, symm_fac, &
                        iset, jset, kset, lset, iatom, jatom, katom, latom, pa1, pb1, pc1, pd1, nspins)

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
    TYPE(x_hf_potential_type)                :: potential_parameter
    REAL(dp)                                 :: box_size(3)
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: screen_periodic
    REAL(dp), INTENT(IN)                     :: eps_schwarz
    INTEGER, INTENT(IN)                      :: max_k
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: sphi_a, sphi_b, sphi_c, sphi_d
    REAL(dp), DIMENSION(:)                   :: T1
    REAL(dp), DIMENSION(:, :)                :: full_ks, full_ks_beta, &
                                                full_density, &
                                                full_density_beta
    REAL(dp)                                 :: fac, symm_fac
    INTEGER, INTENT(IN)                      :: iset, jset, kset, lset, &
                                                iatom, jatom, katom, latom, &
                                                pa1, pb1, pc1, pd1, nspins

    INTEGER                                  :: ma, mb, mc, md, p_index, pa, &
                                                pa2, pb, pb2, pc, pc2, pd, pd2
    REAL(dp)                                 :: tmp

    IF( .NOT.do_periodic ) THEN 
      CALL coulomb4(ra, rb, rc, rd,&
                    zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset), &
                    zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                    la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                    lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                    npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                    ncoa,ncob,ncoc,ncod,&
                    primitive_integrals,&
                    potential_parameter)
    ELSE
      CALL coulomb4_periodic(ra, rb, rc, rd,&
                             zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset), &
                             zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                             la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                             lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                             npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                             ncoa,ncob,ncoc,ncod,&
                             primitive_integrals,&
                             potential_parameter, box_size, &
                             qs_env%x_data%neighbor_cells, screen_periodic(iatom,jatom,iset,jset,:),&
                             screen_periodic(katom,latom,kset,lset,:), eps_schwarz, max_k)
    END IF

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
!***************************************************************************

!!***** four_center_integrals/load_balance_energy [1.0] *
!!
!!   NAME
!!     load_balance_energy 
!!
!!   FUNCTION
!!     Distributes the computation of eri's to all available processes.
!!
!!   NOTES
!!     The optimization is done via a binning procedure followed by simple 
!!     Monte Carlo procedure:
!!     In a first step the total amount of integrals in the system is calculated,
!!     taking into account the sparsity of the KS-matrix , the screening based
!!     on near/farfield approximations and if desired the screening on an initial
!!     density matrix. In a second step, these integrals are distributed to bins
!!     of a given size. Finally, a Monte Carlo procedure optimizes the assignment
!!     of the different loads to each process
!!     At the end each process owns an unique array of inidices-ranges 
!!     that are used to decide whether a process has to calculate a certain
!!     bunch of integrals or not
!!     
!!
!!   INPUTS
!!     - x_data: Object that stores the indices array
!!     - eps_schwarz: screening parameter
!!     - particle_set, atomic_kind_set, para_env
!!     - potential_parameter: short/longrange screening information 
!!     - screen_set, screen_atom, screen_max: screening matrices
!!     - is_assoc_atomic_kind: KS-matrix sparsity
!!     - ks_matrix: used for parallel computation of total number of integrals
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     06.2007 created [Manuel Guidon]
!!     08.2007 new parallel scheme [Manuel Guidon]
!!     09.2007 new 'modulo' parellel scheme and Monte Carlo step [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE load_balance_energy(x_data, eps_schwarz, eps_box, particle_set,atomic_kind_set,para_env,&
                                 potential_parameter, screen_set, screen_atom, screen_max,&
                                 is_assoc_atomic_block_global,ks_matrix,do_periodic,error, pmax_set)
    TYPE(x_hf_type), POINTER                 :: x_data
    REAL(KIND=dp), INTENT(IN)                :: eps_schwarz, eps_box
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(x_hf_potential_type), INTENT(IN)    :: potential_parameter
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: screen_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: screen_atom
    REAL(KIND=dp), INTENT(IN)                :: screen_max
    INTEGER, DIMENSION(:, :), POINTER        :: is_assoc_atomic_block_global
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    LOGICAL                                  :: do_periodic
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(dp), DIMENSION(:, :, :, :), &
      OPTIONAL, POINTER                      :: pmax_set

    CHARACTER(LEN=*), PARAMETER :: routineN = 'load_balance_energy', &
      routineP = moduleN//':'//routineN

    INTEGER :: box_idx1, box_idx2, box_idx3, box_idx4, handle, i, iatom, &
      ibin, icpu, ipos, iset, jatom, jset, katom, kset, latom, lset, max_set, &
      mepos, natom, nbins, ncpu, nseta, nsetb, nsetc, nsetd, pgf_radius_box, &
      stat, steps
    INTEGER(int_8)                           :: optimal_cost, tmp_counter1
    INTEGER(int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: optimal_costs
    INTEGER(int_8), DIMENSION(:), POINTER    :: cost_matrix
    INTEGER(KIND=int_8) :: current_counter, distribution_counter_end, &
      distribution_counter_start, int_per_cpu, restore_counter, &
      total_integrals
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_pos
    INTEGER, DIMENSION(:), POINTER           :: distribution_vector, npgfa, &
                                                npgfb, npgfc, npgfd, nsgfa, &
                                                nsgfb, nsgfc, nsgfd
    LOGICAL                                  :: do_p_screening, failure
    REAL(dp)                                 :: initial_temperature, K, Kbox, &
                                                omega_box, ra(3), rb(3), &
                                                Rbox, rc(3), rd(3)
    REAL(KIND=dp)                            :: max_val1, max_val2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zeta, zetb, zetc, zetd
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(x_hf_distribution), DIMENSION(:), &
      POINTER                                :: ptr_to_tmp_dist
    TYPE(x_hf_distribution), &
      DIMENSION(:, :), POINTER               :: binned_dist, full_dist, &
                                                tmp_dist

    do_p_screening = .FALSE.
    IF(PRESENT(pmax_set)) do_p_screening = .TRUE.
    
    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set) 
    CALL timeset(routineN,handle)
    failure=.FALSE.

    ncpu = para_env%num_pe
    !! If there is only 1 cpu skip the binning
    IF(ncpu == 1) THEN
      ALLOCATE(tmp_dist(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1,1)%number_of_integrals = HUGE(tmp_dist(1,1)%number_of_integrals)
      tmp_dist(1,1)%istart = 1
      tmp_dist(1,1)%jstart = 1
      tmp_dist(1,1)%kstart = 1
      tmp_dist(1,1)%lstart = 1
      tmp_dist(1,1)%original_cpu_id = 1
      ptr_to_tmp_dist => tmp_dist(1,:)
      CALL x_hf_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)
      DEALLOCATE(tmp_dist,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
      !! Calculate total numbers of integrals that have to be calculated (wrt screening and symmetry)
      total_integrals = 0
      DO iatom=1,natom
        IF(MOD(iatom,ncpu) == para_env%mepos) THEN
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, orb_basis_set=orb_basis)
          ra=particle_set(iatom)%r(:)
          CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nseta, nsgf_set=nsgfa, zet=zeta, npgf=npgfa)
          DO jatom=iatom,natom
            IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
            CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, orb_basis_set=orb_basis)
            rb=particle_set(jatom)%r(:)
            CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetb, nsgf_set=nsgfb, zet=zetb, npgf=npgfb)
            DO katom=1,natom
              CALL get_atomic_kind(atomic_kind=particle_set(katom)%atomic_kind, orb_basis_set=orb_basis)
              rc=particle_set(katom)%r(:)
              CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetc, nsgf_set=nsgfc, zet=zetc, npgf=npgfc)
              DO latom=katom,natom
                IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
                IF(.NOT. (is_assoc_atomic_block_global(iatom,latom)  >=1 .AND. &
                          is_assoc_atomic_block_global(iatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,latom)  >=1 ) ) CYCLE 
 
                IF((katom+latom<=iatom+jatom))  THEN
                  IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                 
                  CALL get_atomic_kind(atomic_kind=particle_set(latom)%atomic_kind, orb_basis_set=orb_basis)
                  rd=particle_set(latom)%r(:)
                  CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetd, nsgf_set=nsgfd, zet=zetd, npgf=npgfd)

                  DO iset = 1, nseta
                    DO jset = 1,nsetb
                      max_val1 = screen_set(iatom,jatom,iset,jset)
                      IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                      DO kset = 1,nsetc
                        DO lset = 1,nsetd
                          max_val2 = max_val1 * screen_set(katom,latom,kset,lset)
                          IF( do_p_screening) THEN
                            max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                                  pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                          END IF
                          IF(max_val2<eps_schwarz) CYCLE
                          box_idx1 = (iatom-1)*natom+jatom
                          box_idx2 = (iset-1)*max_set+jset
                          box_idx3 = (katom-1)*natom+latom
                          box_idx4 = (kset-1)*max_set+lset
                          Rbox = (x_data%far_field_boxes(box_idx1,box_idx2)%center(1) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(1))**2 + &
                                 (x_data%far_field_boxes(box_idx1,box_idx2)%center(2) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(2))**2 + &
                                 (x_data%far_field_boxes(box_idx1,box_idx2)%center(3) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(3))**2
                          Rbox = SQRT(Rbox)
                          Rbox = Rbox - x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                        x_data%far_field_boxes(box_idx3, box_idx4)%R
                          pgf_radius_box = INT(1.0_dp/ &
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                          pgf_radius_box = pgf_radius_box + INT(1.0_dp/ &
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2
                          IF(Rbox>pgf_radius_box) THEN
                            max_val2 = 0.0_dp
                            SELECT CASE(potential_parameter%potential_type)
                              CASE(do_hf_potential_coulomb)
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = 1.0_dp/Rbox*Kbox
                                IF( do_p_screening) THEN
                                  max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                             pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                                END IF
                              CASE(do_hf_potential_short)
                                omega_box = 1.0_dp/x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                            1.0_dp/x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                            1.0_dp/(potential_parameter%omega**2)
                                omega_box = 1.0_dp/omega_box
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox
                                IF( do_p_screening) THEN
                                  max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                                 pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                                END IF
                            END SELECT
                            IF(max_val2<eps_box) CYCLE
                          END IF
                          total_integrals = total_integrals + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                END IF
              ENDDO
            ENDDO
          ENDDO
        END IF
      ENDDO

      !! Bin the total amount of integrals on each cpu
      nbins = 16
      !floor value
      int_per_cpu = (total_integrals + nbins -1 ) / (nbins) 

      ALLOCATE(binned_dist(1,nbins), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      binned_dist(:,:)%istart = 0
      binned_dist(:,:)%jstart = 0
      binned_dist(:,:)%kstart = 0
      binned_dist(:,:)%lstart = 0
      binned_dist(:,:)%number_of_integrals = 0
      binned_dist(:,:)%cost = 0
      binned_dist(:,:)%original_cpu_id = 0
 
      tmp_counter1=0
      current_counter = 0
      restore_counter = 0
      ipos = 1
      mepos = 1
      distribution_counter_start = 1
      distribution_counter_end = 0
      ibin = 1
      DO iatom=1,natom
        IF(MOD(iatom,ncpu) == para_env%mepos) THEN
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, orb_basis_set=orb_basis)
          ra=particle_set(iatom)%r(:)
          CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nseta, nsgf_set=nsgfa, zet=zeta, npgf=npgfa)
          DO jatom=iatom,natom
            IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
            CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, orb_basis_set=orb_basis)
            rb=particle_set(jatom)%r(:)
            CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetb, nsgf_set=nsgfb, zet=zetb, npgf=npgfb)
            DO katom=1,natom
              CALL get_atomic_kind(atomic_kind=particle_set(katom)%atomic_kind, orb_basis_set=orb_basis)
              rc=particle_set(katom)%r(:)
              CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetc, nsgf_set=nsgfc, zet=zetc, npgf=npgfc)
              DO latom=katom,natom
                IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
                IF(.NOT. (is_assoc_atomic_block_global(iatom,latom)  >=1 .AND. &
                          is_assoc_atomic_block_global(iatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,latom)  >=1 ) ) CYCLE 
 
                IF((katom+latom<=iatom+jatom))  THEN
                  IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                  distribution_counter_end = distribution_counter_end + 1

                  IF(binned_dist(1,ibin)%istart == 0 ) binned_dist(1,ibin)%istart = iatom
                  IF(binned_dist(1,ibin)%jstart == 0 ) binned_dist(1,ibin)%jstart = jatom
                  IF(binned_dist(1,ibin)%kstart == 0 ) binned_dist(1,ibin)%kstart = katom
                  IF(binned_dist(1,ibin)%lstart == 0 ) binned_dist(1,ibin)%lstart = latom

                  CALL get_atomic_kind(atomic_kind=particle_set(latom)%atomic_kind, orb_basis_set=orb_basis)
                  rd=particle_set(latom)%r(:)
                  CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetd, nsgf_set=nsgfd, zet=zetd, npgf=npgfd)

                  
                  DO iset = 1, nseta
                    DO jset = 1,nsetb
                      max_val1 = screen_set(iatom,jatom,iset,jset)
                      IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                      DO kset = 1,nsetc
                        DO lset = 1,nsetd
                          max_val2 = max_val1 * screen_set(katom,latom,kset,lset)
                          IF( do_p_screening) THEN
                            max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                                  pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                          END IF
                          IF(max_val2<eps_schwarz) CYCLE
                          box_idx1 = (iatom-1)*natom+jatom
                          box_idx2 = (iset-1)*max_set+jset
                          box_idx3 = (katom-1)*natom+latom
                          box_idx4 = (kset-1)*max_set+lset
                          Rbox = (x_data%far_field_boxes(box_idx1,box_idx2)%center(1) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(1))**2 + &
                                 (x_data%far_field_boxes(box_idx1,box_idx2)%center(2) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(2))**2 + &
                                 (x_data%far_field_boxes(box_idx1,box_idx2)%center(3) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(3))**2
                          Rbox = SQRT(Rbox)
                          Rbox = Rbox - x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                        x_data%far_field_boxes(box_idx3, box_idx4)%R
                          pgf_radius_box = INT(1.0_dp/ &
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                          pgf_radius_box = pgf_radius_box + INT(1.0_dp/ &
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2
                          IF(Rbox>pgf_radius_box) THEN
                            max_val2 = 0.0_dp
                            SELECT CASE(potential_parameter%potential_type)
                              CASE(do_hf_potential_coulomb)
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = 1.0_dp/Rbox*Kbox
                                IF( do_p_screening) THEN
                                  max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                             pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                                END IF
                              CASE(do_hf_potential_short)
                                omega_box = 1.0_dp/x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                            1.0_dp/x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                            1.0_dp/(potential_parameter%omega**2)
                                omega_box = 1.0_dp/omega_box
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox
                                IF( do_p_screening) THEN
                                  max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                                 pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                                END IF
                            END SELECT
                            IF(max_val2<eps_box) CYCLE
                          END IF
                          current_counter = current_counter + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO

                  IF( current_counter >= int_per_cpu) THEN
                    binned_dist(1,ibin)%number_of_integrals = distribution_counter_end - distribution_counter_start +1
                    binned_dist(1,ibin)%cost = binned_dist(1,ibin)%cost + current_counter
                    ibin = ibin + 1
                    distribution_counter_start = distribution_counter_end + 1
                    restore_counter = current_counter
                    current_counter = 0
                  END IF
                END IF
              ENDDO
            ENDDO
          ENDDO
        END IF
      ENDDO

      !! Fill the last bin if necessary
      IF( current_counter /= 0 ) THEN  
        binned_dist(1,ibin)%cost = binned_dist(1,ibin)%cost + current_counter
        binned_dist(1,ibin)%number_of_integrals = distribution_counter_end-distribution_counter_start+1
      END IF

      !! store the original cpu_di
      binned_dist(1,:)%original_cpu_id = para_env%mepos + 1
     
      !! store all local results in a big cost matrix
      ALLOCATE(cost_matrix(ncpu*nbins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      cost_matrix = 0
      DO icpu=1,ncpu
        IF(para_env%mepos+1 == icpu ) THEN
          DO i=1,nbins
            cost_matrix((icpu-1)*nbins + i) = binned_dist(1,i)%cost
          END DO
        END IF
      END DO
      !! Communicate the results
      CALL mp_sum_i8(cost_matrix(:), para_env%group) 
      
      !! Find an optimal distribution i.e. assign each element of the cost matrix to a certain process
      ALLOCATE(distribution_vector(ncpu*nbins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      steps = 10000
      initial_temperature = 1.0_dp
      ALLOCATE(optimal_costs(ncpu), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      optimal_costs = 0
      optimal_cost = 0
      !! Optimize load on each node independently
      CALL x_hf_optimize_distribution(para_env, ncpu*nbins, ncpu, cost_matrix, distribution_vector, &
                                      steps, initial_temperature, optimal_cost, error)
      
      !! Find optimal load balance wrt the cost
      optimal_costs(para_env%mepos+1) = optimal_cost
      CALL mp_sum_i8(optimal_costs(1:ncpu), para_env%group)
      icpu = 0
      optimal_cost = HUGE(optimal_cost)
      !! Get process id of the best distribution
      DO i = 1,ncpu
        IF(optimal_costs(i) < optimal_cost ) THEN
          optimal_cost = optimal_costs(i)
          icpu = i
        END IF
      END DO

      !! Broadcast the best distribution_vector
      CALL mp_bcast(distribution_vector,icpu-1,para_env%group) 


      !! Collect local data to global array
      ALLOCATE(full_dist(ncpu,nbins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      full_dist(:,:)%istart = 0
      full_dist(:,:)%jstart = 0
      full_dist(:,:)%kstart = 0
      full_dist(:,:)%lstart = 0
      full_dist(:,:)%number_of_integrals = 0
      full_dist(:,:)%cost = 0
      full_dist(:,:)%original_cpu_id = 0
      
      DO mepos = 1,ncpu
        IF(para_env%mepos == mepos-1) THEN
          full_dist(mepos,:) = binned_dist(1,:)
        END IF
      END DO

      DO mepos=1,ncpu
        DO i=1,nbins
          CALL mp_sum(full_dist(mepos,i)%istart, para_env%group)
          CALL mp_sum(full_dist(mepos,i)%jstart, para_env%group)
          CALL mp_sum(full_dist(mepos,i)%kstart, para_env%group)
          CALL mp_sum(full_dist(mepos,i)%lstart, para_env%group)
          CALL mp_sum_i8(full_dist(mepos,i)%number_of_integrals, para_env%group)
          CALL mp_sum_i8(full_dist(mepos,i)%cost, para_env%group)
          CALL mp_sum(full_dist(mepos,i)%original_cpu_id, para_env%group)
        END DO
      END DO

      !! reorder the distribution according to the distribution vector
      ALLOCATE(tmp_pos(ncpu),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_pos = 1
      ALLOCATE(tmp_dist(ncpu,nbins*ncpu),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      
      DO mepos = 1,ncpu
        DO icpu=1,ncpu
          DO i=1,nbins
            IF(distribution_vector((icpu-1)*nbins + i) == mepos) THEN
              tmp_dist(mepos,tmp_pos(mepos)) = full_dist(icpu,i)
              tmp_pos(mepos) = tmp_pos(mepos) + 1
            END IF
          END DO
        END DO
      END DO
      
      !! Assign the load to each process
      NULLIFY(ptr_to_tmp_dist) 
      DO icpu=1,ncpu
        IF(para_env%mepos == icpu -1) THEN
          ptr_to_tmp_dist => tmp_dist(icpu,1:tmp_pos(icpu)-1)
        END IF
      END DO
     
      CALL x_hf_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)
      DEALLOCATE(tmp_dist, full_dist, tmp_pos, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(optimal_costs, distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(cost_matrix,binned_dist, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)
  END SUBROUTINE load_balance_energy
!***************************************************************************

!!***** four_center_integrals/load_balance_forces [1.0] *
!!
!!   NAME
!!     load_balance_forces
!!
!!   FUNCTION
!!     Distributes the computation of eri-derivatives to all available processes.
!!
!!   NOTES
!!     The optimization is done via a simple binning procedure:
!!     In a first step the total amount of integrals in the system is calculated,
!!     taking into account the sparsity of the KS-matrix, the screening based
!!     on the Schwarz-Inequality and the density-matrix. Then, the procedure uses 
!!     different bin-sizes in order to optimize the load balance.
!!     At the end each process owns an unique array of inidices-ranges 
!!     that are used to decide whether a process has to calculate a certain
!!     bunch of integrals.
!!     
!!
!!   INPUTS
!!     - x_data: Object that stores the indices array
!!     - eps_schwarz: screening parameter
!!     - particle_set, atomic_kind_set, para_env
!!     - potential_parameter: short/longrange screening information 
!!     - screen_set, screen_atom, screen_max, pmax_set: screening matrices
!!     - is_assoc_atomic_kind: KS-matrix sparsity
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     06.2007 created [Manuel Guidon]
!!     08.2007 new parallel scheme [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE load_balance_forces(x_data, eps_schwarz, eps_box, particle_set,atomic_kind_set,para_env,&
                                 potential_parameter, screen_set, screen_atom, screen_max,&
                                 pmax_set, is_assoc_atomic_block_global, ks_matrix, do_periodic, &
                                 error)
    TYPE(x_hf_type), POINTER                 :: x_data
    REAL(KIND=dp), INTENT(IN)                :: eps_schwarz, eps_box
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(x_hf_potential_type), INTENT(IN)    :: potential_parameter
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: screen_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: screen_atom
    REAL(KIND=dp), INTENT(IN)                :: screen_max
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: pmax_set
    INTEGER, DIMENSION(:, :), POINTER        :: is_assoc_atomic_block_global
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    LOGICAL                                  :: do_periodic
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'load_balance_forces', &
      routineP = moduleN//':'//routineN

    INTEGER :: box_idx1, box_idx2, box_idx3, box_idx4, first_from, first_to, &
      handle, i, iatom, icpu, ipos, iset, jatom, jset, katom, kset, latom, &
      lset, max_set, mepos, natom, nseta, nsetb, nsetc, nsetd, &
      pgf_radius_box, restore_iatom, restore_jatom, restore_katom, &
      restore_latom, stat, ncpu, ibin, steps, nbins
    INTEGER(KIND=int_8) :: current_counter, distribution_counter_end, &
      distribution_counter_start, int_per_cpu, max_dev, min_dev, &
      opt_max, opt_nbin, restore_counter, total_integrals
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: cost_per_cpu
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: cost_array
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_pos
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: first_dist
    INTEGER, DIMENSION(:), POINTER           :: nsgfa, nsgfb, nsgfc, nsgfd, distribution_vector
    LOGICAL                                  :: failure
    REAL(dp)                                 :: Kbox, omega_box, Rbox
    REAL(KIND=dp)                            :: max_val1, max_val2
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(x_hf_distribution), DIMENSION(:), &
      POINTER                                :: ptr_to_tmp_dist
    TYPE(x_hf_distribution), &
      DIMENSION(:, :), POINTER               :: tmp_dist, binned_dist, full_dist
    INTEGER(int_8)                           :: optimal_cost, tmp_counter1
    INTEGER(int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: optimal_costs
    INTEGER(int_8), DIMENSION(:), POINTER    :: cost_matrix
    REAL(dp)                                 :: initial_temperature

    natom=SIZE(particle_set)

    CALL timeset(routineN,handle)
    failure=.FALSE.
    ncpu = para_env%num_pe
    IF(ncpu == 1) THEN
      ALLOCATE(tmp_dist(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1,1)%number_of_integrals = HUGE(tmp_dist(1,1)%number_of_integrals)
      tmp_dist(1,1)%istart = 1
      tmp_dist(1,1)%jstart = 1
      tmp_dist(1,1)%kstart = 1
      tmp_dist(1,1)%lstart = 1
      tmp_dist(1,1)%original_cpu_id = 1
      ptr_to_tmp_dist => tmp_dist(1,:)
      CALL x_hf_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)
      DEALLOCATE(tmp_dist,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxnset=max_set)
      total_integrals = 0
      DO iatom=1,natom
        IF(MOD(iatom,ncpu) == para_env%mepos) THEN
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, orb_basis_set=orb_basis)
          CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nseta, nsgf_set=nsgfa)
          DO jatom=iatom,natom 
            IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
            CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, orb_basis_set=orb_basis)
            CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetb, nsgf_set=nsgfb)
            DO katom=1,natom
              CALL get_atomic_kind(atomic_kind=particle_set(katom)%atomic_kind, orb_basis_set=orb_basis)
              CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetc, nsgf_set=nsgfc)
              DO latom=katom,natom
                !All four centers equivalent => zero-contribution
                IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE
  
                IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
                IF(.NOT. (is_assoc_atomic_block_global(iatom,latom)  >=1 .AND. &
                          is_assoc_atomic_block_global(iatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,latom)  >=1 ) ) CYCLE
  
                IF((katom+latom<=iatom+jatom))  THEN
                  IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                
                  CALL get_atomic_kind(atomic_kind=particle_set(latom)%atomic_kind, orb_basis_set=orb_basis)
                  CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetd, nsgf_set=nsgfd)
  
                  DO iset = 1, nseta
                    DO jset = 1,nsetb
                      max_val1 = screen_set(iatom,jatom,iset,jset)
                      IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                      DO kset = 1,nsetc
                        DO lset = 1,nsetd
                          max_val2 = max_val1 * screen_set(katom,latom,kset,lset) * &
                                     2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                                pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                          IF(max_val2<eps_schwarz) CYCLE

                          box_idx1 = (iatom-1)*natom+jatom
                          box_idx2 = (iset-1)*max_set+jset
                          box_idx3 = (katom-1)*natom+latom
                          box_idx4 = (kset-1)*max_set+lset
                          Rbox = (x_data%far_field_boxes(box_idx1,box_idx2)%center(1) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(1))**2 + &
                                 (x_data%far_field_boxes(box_idx1,box_idx2)%center(2) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(2))**2 + &
                                 (x_data%far_field_boxes(box_idx1,box_idx2)%center(3) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(3))**2
                          Rbox = SQRT(Rbox)
                          Rbox = Rbox - x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                        x_data%far_field_boxes(box_idx3, box_idx4)%R
                          pgf_radius_box = INT(1.0_dp/ &
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                          pgf_radius_box = pgf_radius_box + INT(1.0_dp/ &
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2

                          IF(Rbox>pgf_radius_box) THEN
                            max_val2 = 0.0_dp
                            SELECT CASE(potential_parameter%potential_type)
                              CASE(do_hf_potential_coulomb)
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = 1.0_dp/Rbox*Kbox*&
                                           2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                           pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                              CASE(do_hf_potential_short)
                                omega_box = 1.0_dp/x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                            1.0_dp/x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                            1.0_dp/(potential_parameter%omega**2)
                                omega_box = 1.0_dp/omega_box
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox* &
                                           2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                           pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                            END SELECT
                            IF(max_val2<eps_box) CYCLE
                          END IF
                          total_integrals = total_integrals + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                END IF
              ENDDO
            ENDDO
          ENDDO
       END IF
      ENDDO
      
      nbins = 16
      int_per_cpu = (total_integrals + nbins -1 ) / nbins
      ALLOCATE(binned_dist(1,nbins), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      binned_dist(:,:)%istart = 0
      binned_dist(:,:)%jstart = 0
      binned_dist(:,:)%kstart = 0
      binned_dist(:,:)%lstart = 0
      binned_dist(:,:)%number_of_integrals = 0
      binned_dist(:,:)%cost = 0
      binned_dist(:,:)%original_cpu_id = 0
 
      tmp_counter1=0
      current_counter = 0
      restore_counter = 0
      ipos = 1
      mepos = 1
      distribution_counter_start = 1
      distribution_counter_end = 0
      ibin = 1
      DO iatom=1,natom
        IF(MOD(iatom,ncpu) == para_env%mepos) THEN
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, orb_basis_set=orb_basis)
          CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nseta, nsgf_set=nsgfa)
          DO jatom=iatom,natom
            IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
            CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, orb_basis_set=orb_basis)
            CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetb, nsgf_set=nsgfb)
            DO katom=1,natom
              CALL get_atomic_kind(atomic_kind=particle_set(katom)%atomic_kind, orb_basis_set=orb_basis)
              CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetc, nsgf_set=nsgfc)
              DO latom=katom,natom
                !All four centers equivalent => zero-contribution
                IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE
   
                IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
                IF(.NOT. (is_assoc_atomic_block_global(iatom,latom)  >=1 .AND. &
                          is_assoc_atomic_block_global(iatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,latom)  >=1 ) ) CYCLE
  
                IF((katom+latom<=iatom+jatom))  THEN
                  IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                  distribution_counter_end = distribution_counter_end + 1
              
                  IF(binned_dist(1,ibin)%istart == 0 ) binned_dist(1,ibin)%istart = iatom
                  IF(binned_dist(1,ibin)%jstart == 0 ) binned_dist(1,ibin)%jstart = jatom
                  IF(binned_dist(1,ibin)%kstart == 0 ) binned_dist(1,ibin)%kstart = katom
                  IF(binned_dist(1,ibin)%lstart == 0 ) binned_dist(1,ibin)%lstart = latom

                  CALL get_atomic_kind(atomic_kind=particle_set(latom)%atomic_kind, orb_basis_set=orb_basis)
                  CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetd, nsgf_set=nsgfd)
  
                  DO iset = 1, nseta
                    DO jset = 1,nsetb
                     max_val1 = screen_set(iatom,jatom,iset,jset)
                      IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                      DO kset = 1,nsetc
                       DO lset = 1,nsetd
                          max_val2 = max_val1 * screen_set(katom,latom,kset,lset) * &
                                     2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                                pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                          IF(max_val2<eps_schwarz) CYCLE
                          box_idx1 = (iatom-1)*natom+jatom
                          box_idx2 = (iset-1)*max_set+jset
                          box_idx3 = (katom-1)*natom+latom
                          box_idx4 = (kset-1)*max_set+lset
                          Rbox = (x_data%far_field_boxes(box_idx1,box_idx2)%center(1) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(1))**2 + &
                                 (x_data%far_field_boxes(box_idx1,box_idx2)%center(2) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(2))**2 + &
                                 (x_data%far_field_boxes(box_idx1,box_idx2)%center(3) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(3))**2
                          Rbox = SQRT(Rbox)
                          Rbox = Rbox - x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                        x_data%far_field_boxes(box_idx3, box_idx4)%R
                          pgf_radius_box = INT(1.0_dp/ &
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                          pgf_radius_box = pgf_radius_box + INT(1.0_dp/ &
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2
                          IF(Rbox>pgf_radius_box) THEN
                            max_val2 = 0.0_dp
                            SELECT CASE(potential_parameter%potential_type)
                              CASE(do_hf_potential_coulomb)
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = 1.0_dp/Rbox*Kbox * &
                                           2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                           pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                              CASE(do_hf_potential_short)
                                omega_box = 1.0_dp/x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                            1.0_dp/x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                            1.0_dp/(potential_parameter%omega**2)
                                omega_box = 1.0_dp/omega_box
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox * &
                                           2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                           pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                            END SELECT
                            IF(max_val2<eps_box) CYCLE
                          END IF
                          current_counter = current_counter + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                  IF( current_counter >= int_per_cpu) THEN
                    binned_dist(1,ibin)%number_of_integrals = distribution_counter_end - distribution_counter_start +1
                    binned_dist(1,ibin)%cost = binned_dist(1,ibin)%cost + current_counter
                    ibin = ibin + 1
                    distribution_counter_start = distribution_counter_end + 1
                    restore_counter = current_counter
                    current_counter = 0
                  END IF
                END IF
              ENDDO
            ENDDO
          ENDDO
        END IF
      ENDDO
      
      !! Fill the last bin if necessary
      IF( current_counter /= 0 ) THEN  
        binned_dist(1,ibin)%cost = binned_dist(1,ibin)%cost + current_counter
        binned_dist(1,ibin)%number_of_integrals = distribution_counter_end-distribution_counter_start+1
      END IF

      !! store the original cpu_di
      binned_dist(1,:)%original_cpu_id = para_env%mepos + 1
     
      !! store all local results in a big cost matrix
      ALLOCATE(cost_matrix(ncpu*nbins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      cost_matrix = 0
      DO icpu=1,ncpu
        IF(para_env%mepos+1 == icpu ) THEN
          DO i=1,nbins
            cost_matrix((icpu-1)*nbins + i) = binned_dist(1,i)%cost
          END DO
        END IF
      END DO


      !! Communicate the results
      CALL mp_sum_i8(cost_matrix(:), para_env%group) 

      !! Find an optimal distribution i.e. assign each element of the cost matrix to a certain process
      ALLOCATE(distribution_vector(ncpu*nbins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      steps = 10000
      initial_temperature = 1.0_dp
      ALLOCATE(optimal_costs(ncpu), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      optimal_costs = 0
      optimal_cost = 0
      !! Optimize load on each node independently
      CALL x_hf_optimize_distribution(para_env, ncpu*nbins, ncpu, cost_matrix, distribution_vector, &
                                      steps, initial_temperature, optimal_cost, error)
      
      !! Find optimal load balance wrt the cost
      optimal_costs(para_env%mepos+1) = optimal_cost
      CALL mp_sum_i8(optimal_costs(1:ncpu), para_env%group)
      icpu = 0
      optimal_cost = HUGE(optimal_cost)
      !! Get process id of the best distribution
      DO i = 1,ncpu
        IF(optimal_costs(i) < optimal_cost ) THEN
          optimal_cost = optimal_costs(i)
          icpu = i
        END IF
      END DO

      !! Broadcast the best distribution_vector
      CALL mp_bcast(distribution_vector,icpu-1,para_env%group) 


      !! Collect local data to global array
      ALLOCATE(full_dist(ncpu,nbins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      full_dist(:,:)%istart = 0
      full_dist(:,:)%jstart = 0
      full_dist(:,:)%kstart = 0
      full_dist(:,:)%lstart = 0
      full_dist(:,:)%number_of_integrals = 0
      full_dist(:,:)%cost = 0
      full_dist(:,:)%original_cpu_id = 0
      
      DO mepos = 1,ncpu
        IF(para_env%mepos == mepos-1) THEN
          full_dist(mepos,:) = binned_dist(1,:)
        END IF
      END DO

      DO mepos=1,ncpu
        DO i=1,nbins
          CALL mp_sum(full_dist(mepos,i)%istart, para_env%group)
          CALL mp_sum(full_dist(mepos,i)%jstart, para_env%group)
          CALL mp_sum(full_dist(mepos,i)%kstart, para_env%group)
          CALL mp_sum(full_dist(mepos,i)%lstart, para_env%group)
          CALL mp_sum_i8(full_dist(mepos,i)%number_of_integrals, para_env%group)
          CALL mp_sum_i8(full_dist(mepos,i)%cost, para_env%group)
          CALL mp_sum(full_dist(mepos,i)%original_cpu_id, para_env%group)
        END DO
      END DO

      !! reorder the distribution according to the distribution vector
      ALLOCATE(tmp_pos(ncpu),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_pos = 1
      ALLOCATE(tmp_dist(ncpu,nbins*ncpu),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      
      DO mepos = 1,ncpu
        DO icpu=1,ncpu
          DO i=1,nbins
            IF(distribution_vector((icpu-1)*nbins + i) == mepos) THEN
              tmp_dist(mepos,tmp_pos(mepos)) = full_dist(icpu,i)
              tmp_pos(mepos) = tmp_pos(mepos) + 1
            END IF
          END DO
        END DO
      END DO

      
      !! Assign the load to each process
      NULLIFY(ptr_to_tmp_dist) 
      DO icpu=1,ncpu
        IF(para_env%mepos == icpu -1) THEN
          ptr_to_tmp_dist => tmp_dist(icpu,1:tmp_pos(icpu)-1)
        END IF
      END DO
      
      CALL x_hf_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)
      DEALLOCATE(tmp_dist, full_dist, tmp_pos, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(optimal_costs, distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(cost_matrix,binned_dist, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE load_balance_forces

END MODULE four_center_integrals
