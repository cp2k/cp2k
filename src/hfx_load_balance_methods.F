!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines for optimizing load balance between processes in HFX calculations
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE hfx_load_balance_methods
  USE cell_types,                      ONLY: cell_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE hfx_pair_list_methods,           ONLY: build_atomic_pair_list,&
                                             build_pair_list
  USE hfx_types,                       ONLY: &
       hfx_basis_type, hfx_block_range_type, hfx_distribution, &
       hfx_load_balance_type, hfx_p_kind, hfx_screen_coeff_type, &
       hfx_set_distr_energy, hfx_set_distr_forces, hfx_type, pair_list_type, &
       pair_set_list_type
  USE input_constants,                 ONLY: hfx_do_eval_energy,&
                                             hfx_do_eval_forces
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_gatherv,&
                                             mp_isendrecv,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync,&
                                             mp_waitall
  USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: hfx_load_balance, &
         hfx_update_load_balance,&
         collect_load_balance_info, cost_model, p1_energy,p2_energy,p3_energy

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_load_balance_methods'

  REAL(KIND=dp), PARAMETER :: p1_energy(12)=(/2.9461408209700424_dp,       1.0624718662999657_dp, &
                                             -1.91570128356921242E-002_dp,-1.6668495454436603_dp, &
                                              1.7512639006523709_dp,      -9.76074323945336081E-002_dp, &
                                              2.6230786127311889_dp,      -0.31870737623014189_dp, &
                                              7.9588203912690973_dp,       1.8331423413134813_dp, &
                                             -0.15427618665346299_dp,      0.19749436090711650_dp/)
  REAL(KIND=dp), PARAMETER :: p2_energy(12)=(/ 2.3104682960662593_dp,      1.8744052737304417_dp, &
                                              -9.36564055598656797E-002_dp,0.64284973765086939_dp, &
                                               1.0137565430060556_dp,     -6.80088178288954567E-003_dp, &
                                               1.1692629207374552_dp,     -2.6314710080507573_dp, &
                                              19.237814781880786_dp,       1.0505934173661349_dp, &
                                               0.80382371955699250_dp,     0.49903401991818103_dp/)
  REAL(KIND=dp), PARAMETER :: p3_energy(2) =(/7.82336287670072350E-002_dp, 0.38073304105744837_dp/)
  REAL(KIND=dp), PARAMETER :: p1_forces(12)=(/2.5746279948798874_dp,       1.3420575378609276_dp, &
                                             -9.41673106447732111E-002_dp, 0.94568006899317825_dp, &
                                             -1.4511897117448544_dp,       0.59178934677316952_dp, &
                                              2.7291149361757236_dp,      -0.50555512044800210_dp, &
                                              8.3508180969609871_dp,       1.6829982496141809_dp, &
                                             -0.74895370472152600_dp,      0.43801726744197500_dp/)
  REAL(KIND=dp), PARAMETER :: p2_forces(12)=(/2.6398568961569020_dp,       2.3024918834564101_dp, &
                                              5.33216585432061581E-003_dp, 0.45572145697283628_dp, &
                                              1.8119743851500618_dp,      -0.12533918548421166_dp, &
                                             -1.4040312084552751_dp,      -4.5331650463917859_dp, &
                                             12.593431549069477_dp,        1.1311978374487595_dp, &
                                              1.4245996087624646_dp,       1.1425350529853495_dp/)
  REAL(KIND=dp), PARAMETER :: p3_forces(2) =(/0.12051930516830946_dp,      1.3828051586144336_dp/)

!***

  CONTAINS

! *****************************************************************************
!> \brief Distributes the computation of eri's to all available processes.
!> \param x_data Object that stores the indices array
!> \param eps_schwarz screening parameter
!> \param particle_set , atomic_kind_set, para_env
!> \param max_set Maximum number of set to be considered
!> \param para_env para_env
!> \param coeffs_set screening functions
!> \param coeffs_kind screening functions
!> \param is_assoc_atomic_block_global KS-matrix sparsity
!> \param do_periodic flag for periodicity
!> \param load_balance_parameter Paramters for Monte-Carlo routines
!> \param kind_of helper array for mapping
!> \param basis_parameter Basis set parameters
!> \param pmax_set Initial screening matrix
!> \param pmax_atom ...
!> \param i_thread Process ID of current Thread
!> \param n_threads Total Number of Threads
!> \param cell cell
!> \param do_p_screening Flag for initial p screening
!> \param map_atom_to_kind_atom ...
!> \param nkind ...
!> \param eval_type ...
!> \param pmax_block ...
!> \param use_virial ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2007 created [Manuel Guidon]
!>      08.2007 new parallel scheme [Manuel Guidon]
!>      09.2007 new 'modulo' parellel scheme and Monte Carlo step [Manuel Guidon]
!>      11.2007 parallelize load balance on box_idx1 [Manuel Guidon]
!>      02.2009 completely refactored [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      The optimization is done via a binning procedure followed by simple
!>      Monte Carlo procedure:
!>      In a first step the total amount of integrals in the system is calculated,
!>      taking into account the sparsity of the KS-matrix , the screening based
!>      on near/farfield approximations and if desired the screening on an initial
!>      density matrix.
!>      In a second step, bins are generate that contain approximately the same number
!>      of integrals, and a cost for these bins is estimated (currently the number of integrals)
!>      In a third step, a Monte Carlo procedure optimizes the assignment
!>      of the different loads to each process
!>      At the end each process owns an unique array of *atomic* indices-ranges
!>      that are used to decide whether a process has to calculate a certain
!>      bunch of integrals or not
! *****************************************************************************
  SUBROUTINE hfx_load_balance(x_data, eps_schwarz, particle_set,max_set,para_env,&
                          coeffs_set, coeffs_kind, &
                          is_assoc_atomic_block_global, do_periodic,&
                          load_balance_parameter, kind_of, basis_parameter, pmax_set, &
                          pmax_atom, i_thread, n_threads, cell, &
                          do_p_screening, map_atom_to_kind_atom, nkind, eval_type, &
                          pmax_block, use_virial, error)
    TYPE(hfx_type), POINTER                  :: x_data
    REAL(dp), INTENT(IN)                     :: eps_schwarz
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: max_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: coeffs_kind
    INTEGER, DIMENSION(:, :)                 :: is_assoc_atomic_block_global
    LOGICAL                                  :: do_periodic
    TYPE(hfx_load_balance_type), POINTER     :: load_balance_parameter
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_p_kind), DIMENSION(:), POINTER  :: pmax_set
    REAL(dp), DIMENSION(:, :), POINTER       :: pmax_atom
    INTEGER, INTENT(IN)                      :: i_thread, n_threads
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(IN)                      :: do_p_screening
    INTEGER, DIMENSION(:), POINTER           :: map_atom_to_kind_atom
    INTEGER, INTENT(IN)                      :: nkind, eval_type
    REAL(dp), DIMENSION(:, :), POINTER       :: pmax_block
    LOGICAL, INTENT(IN)                      :: use_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_load_balance', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=512)                       :: error_msg
    INTEGER :: block_size, current_block_id, data_from, dest, handle, &
      handle_inner, handle_range, i, iatom_block, iatom_end, iatom_start, &
      ibin, icpu, j, jatom_block, jatom_end, jatom_start, katom_block, &
      katom_end, katom_start, latom_block, latom_end, latom_start, mepos, &
      my_process_id, n_processes, natom, nbins, nblocks, ncpu, new_iatom_end, &
      new_iatom_start, new_jatom_end, new_jatom_start, non_empty_blocks, &
      objective_block_size, objective_nblocks, req(2), source, stat, &
      total_blocks
    INTEGER(int_8) :: atom_block, cost_per_bin, cost_per_core, current_cost, &
      distribution_counter_end, distribution_counter_start, &
      global_quartet_counter, local_quartet_counter, self_cost_per_block, &
      tmp_block, total_block_self_cost
    INTEGER(int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: buffer_in, buffer_out
    INTEGER(int_8), DIMENSION(:), POINTER    :: local_cost_matrix, recbuffer, &
                                                sendbuffer, swapbuffer
    INTEGER(int_8), DIMENSION(:), POINTER, &
      SAVE                                   :: cost_matrix
    INTEGER(int_8), SAVE                     :: shm_global_quartet_counter, &
                                                shm_local_quartet_counter
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount, rdispl, tmp_index, &
                                                tmp_pos, to_be_sorted
    INTEGER, DIMENSION(:), POINTER, SAVE     :: shm_distribution_vector
    INTEGER, SAVE                            :: shm_nblocks
    LOGICAL                                  :: changed, failure, &
                                                last_bin_needs_to_be_filled, &
                                                optimized
    LOGICAL, DIMENSION(:, :), POINTER, SAVE  :: atomic_pair_list
    REAL(dp)                                 :: coeffs_kind_max0, &
                                                log10_eps_schwarz, log_2, &
                                                pmax_blocks
    TYPE(hfx_block_range_type), &
      DIMENSION(:), POINTER                  :: blocks_guess, tmp_blocks, &
                                                tmp_blocks2
    TYPE(hfx_block_range_type), &
      DIMENSION(:), POINTER, SAVE            :: shm_blocks
    TYPE(hfx_distribution), DIMENSION(:), &
      POINTER                                :: binned_dist, ptr_to_tmp_dist, &
                                                tmp_dist
    TYPE(hfx_distribution), &
      DIMENSION(:, :), POINTER, SAVE         :: full_dist
    TYPE(pair_list_type)                     :: list_ij, list_kl
    TYPE(pair_set_list_type), ALLOCATABLE, &
      DIMENSION(:)                           :: set_list_ij, set_list_kl

!$OMP BARRIER
!$OMP MASTER
    CALL timeset(routineN,handle)
!$OMP END MASTER
!$OMP BARRIER

    failure=.FALSE.

    log10_eps_schwarz = LOG10(eps_schwarz)
    log_2 = LOG10(2.0_dp)
    coeffs_kind_max0=MAXVAL(coeffs_kind(:,:)%x(2))
    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads
    natom=SIZE(particle_set)

    block_size = load_balance_parameter%block_size
    ALLOCATE(set_list_ij((max_set*block_size)**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(set_list_kl((max_set*block_size)**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    IF( .NOT. load_balance_parameter%blocks_initialized ) THEN
!$OMP BARRIER
!$OMP MASTER
      CALL timeset(routineN//"_range",handle_range)

      nblocks = MAX((natom+block_size-1)/block_size,1)
      ALLOCATE(blocks_guess(nblocks),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(tmp_blocks(natom),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(tmp_blocks2(natom),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      pmax_blocks = 0.0_dp
      SELECT CASE(eval_type)
      CASE (hfx_do_eval_energy)
        atomic_pair_list => x_data%atomic_pair_list
      CASE (hfx_do_eval_forces)
        atomic_pair_list => x_data%atomic_pair_list_forces
      END SELECT
      atomic_pair_list = .TRUE.
      CALL init_blocks(nkind, para_env, natom,block_size,nblocks,blocks_guess,&
                       list_ij, list_kl, set_list_ij, set_list_kl, &
                       particle_set,&
                       coeffs_set, coeffs_kind, &
                       is_assoc_atomic_block_global, do_periodic,&
                       kind_of, basis_parameter, pmax_set, pmax_atom, &
                       pmax_blocks, cell, &
                       do_p_screening, map_atom_to_kind_atom, eval_type,&
                       log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, atomic_pair_list)

      total_block_self_cost = 0

      DO i=1,nblocks
         total_block_self_cost = total_block_self_cost + blocks_guess(i)%cost
      END DO

      CALL mp_sum(total_block_self_cost, para_env%group)

      objective_block_size = load_balance_parameter%block_size
      objective_nblocks = MAX((natom+objective_block_size-1)/objective_block_size,1)

      self_cost_per_block =  (total_block_self_cost + objective_nblocks -1 ) / (objective_nblocks)

      DO i = 1,nblocks
        tmp_blocks2(i) = blocks_guess(i)
      END DO

      optimized = .FALSE.
      i = 0
      DO WHILE( .NOT. optimized)
        i = i + 1
        current_block_id = 0
        changed = .FALSE.
        DO atom_block = 1,nblocks
          current_block_id = current_block_id + 1
          iatom_start = tmp_blocks2(atom_block)%istart
          iatom_end = tmp_blocks2(atom_block)%iend
          IF(tmp_blocks2(atom_block)%cost > 1.5_dp * self_cost_per_block .AND. iatom_end-iatom_start > 0) THEN
            changed = .TRUE.
            new_iatom_start = iatom_start
            new_iatom_end  = (iatom_end - iatom_start+1) /2 + iatom_start - 1
            new_jatom_start = new_iatom_end + 1
            new_jatom_end = iatom_end
            tmp_blocks(current_block_id)%istart = new_iatom_start
            tmp_blocks(current_block_id)%iend = new_iatom_end
            tmp_blocks(current_block_id)%cost = estimate_block_cost(natom, nkind, list_ij, list_kl, set_list_ij, set_list_kl, &
                                                new_iatom_start, new_iatom_end, new_iatom_start, new_iatom_end, &
                                                new_iatom_start, new_iatom_end, new_iatom_start, new_iatom_end, &
                                                particle_set,&
                                                coeffs_set, coeffs_kind, &
                                                is_assoc_atomic_block_global, do_periodic,&
                                                kind_of, basis_parameter, pmax_set, pmax_atom, pmax_blocks, &
                                                cell, &
                                                do_p_screening, map_atom_to_kind_atom, eval_type,&
                                                log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, atomic_pair_list)
            current_block_id = current_block_id + 1
            tmp_blocks(current_block_id)%istart = new_jatom_start
            tmp_blocks(current_block_id)%iend = new_jatom_end
            tmp_blocks(current_block_id)%cost = estimate_block_cost(natom, nkind, list_ij, list_kl, set_list_ij, set_list_kl, &
                                                new_jatom_start, new_jatom_end, new_jatom_start, new_jatom_end, &
                                                new_jatom_start, new_jatom_end, new_jatom_start, new_jatom_end, &
                                                particle_set,&
                                                coeffs_set, coeffs_kind, &
                                                is_assoc_atomic_block_global, do_periodic,&
                                                kind_of, basis_parameter, pmax_set, pmax_atom, pmax_blocks, &
                                                cell, &
                                                do_p_screening, map_atom_to_kind_atom, eval_type,&
                                                log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, atomic_pair_list)
          ELSE
            tmp_blocks(current_block_id)%istart = iatom_start
            tmp_blocks(current_block_id)%iend = iatom_end
            tmp_blocks(current_block_id)%cost = tmp_blocks2(atom_block)%cost
          END IF
        END DO
        IF( .NOT. changed) optimized = .TRUE.
        IF( i> 20) optimized = .TRUE.
        nblocks = current_block_id
        DO atom_block = 1,nblocks
          tmp_blocks2(atom_block) = tmp_blocks(atom_block)
        END DO
      END DO

      DEALLOCATE(tmp_blocks2, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      ! ** count number of non empty blocks on each node
      non_empty_blocks = 0
      DO atom_block = 1,nblocks
        IF(tmp_blocks(atom_block)%istart == 0) CYCLE
        non_empty_blocks = non_empty_blocks + 1
      END DO

      ALLOCATE(rcount(ncpu), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      rcount = 0
      rcount(para_env%mepos+1) = non_empty_blocks
      CALL mp_sum(rcount, para_env%group)

      ! ** sum all non_empty_blocks
      total_blocks = 0
      DO i=1,ncpu
        total_blocks = total_blocks + rcount(i)
      END DO

      ! ** calculate offsets
      ALLOCATE(rdispl(ncpu), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      rcount = rcount * 3
      rdispl(1) = 0
      DO i = 2,ncpu
        rdispl(i) = rdispl(i-1) + rcount(i-1)
      END DO

      ALLOCATE(buffer_in(3*non_empty_blocks), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


      non_empty_blocks = 0
      DO atom_block = 1,nblocks
        IF(tmp_blocks(atom_block)%istart == 0) CYCLE
        buffer_in(non_empty_blocks*3+1) = tmp_blocks(atom_block)%istart
        buffer_in(non_empty_blocks*3+2) = tmp_blocks(atom_block)%iend
        buffer_in(non_empty_blocks*3+3) = tmp_blocks(atom_block)%cost
        non_empty_blocks = non_empty_blocks + 1
      END DO

      nblocks = total_blocks

      ALLOCATE(tmp_blocks2(nblocks), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


      ALLOCATE(buffer_out(3*nblocks), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      ! ** Gather all three arrays
      CALL mp_allgather(buffer_in,buffer_out,rcount,rdispl,para_env%group)

      DO i=1,nblocks
        tmp_blocks2(i)%istart = INT(buffer_out((i-1)*3+1))
        tmp_blocks2(i)%iend =   INT(buffer_out((i-1)*3+2))
        tmp_blocks2(i)%cost = buffer_out((i-1)*3+3)
      END DO

      ! ** Now we sort the blocks
      ALLOCATE(to_be_sorted(nblocks),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(tmp_index(nblocks),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO atom_block = 1,nblocks
        to_be_sorted(atom_block) = tmp_blocks2(atom_block)%istart
      END DO

      CALL sort(to_be_sorted,nblocks,tmp_index)

      ALLOCATE(x_data%blocks(nblocks), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO atom_block = 1,nblocks
        x_data%blocks(atom_block) = tmp_blocks2(tmp_index(atom_block))
      END DO

      shm_blocks => x_data%blocks
      shm_nblocks = nblocks

      ! ** Set nblocks in structure
      load_balance_parameter%nblocks = nblocks

      DEALLOCATE(blocks_guess, tmp_blocks, tmp_blocks2, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DEALLOCATE(rcount, rdispl, buffer_in, buffer_out, to_be_sorted, tmp_index, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      load_balance_parameter%blocks_initialized = .TRUE.

      x_data%blocks = shm_blocks
      load_balance_parameter%nblocks = shm_nblocks
      load_balance_parameter%blocks_initialized = .TRUE.

      ALLOCATE(x_data%pmax_block(shm_nblocks,shm_nblocks),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      x_data%pmax_block = 0.0_dp
      pmax_block => x_data%pmax_block
      CALL timestop(handle_range)
!$OMP END MASTER
!$OMP BARRIER

      IF( .NOT. load_balance_parameter%blocks_initialized ) THEN
        ALLOCATE(x_data%blocks(shm_nblocks), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        x_data%blocks = shm_blocks
        load_balance_parameter%nblocks = shm_nblocks
        load_balance_parameter%blocks_initialized = .TRUE.
      END IF
      !! ** precalculate maximum density matrix elements in blocks
!$OMP BARRIER
    END IF

!$OMP BARRIER
!$OMP MASTER
    pmax_block => x_data%pmax_block
    pmax_block = 0.0_dp
    IF( do_p_screening ) THEN
      DO iatom_block = 1,shm_nblocks
        iatom_start = x_data%blocks(iatom_block)%istart
        iatom_end = x_data%blocks(iatom_block)%iend
        DO jatom_block = 1,shm_nblocks
          jatom_start = x_data%blocks(jatom_block)%istart
          jatom_end = x_data%blocks(jatom_block)%iend
          pmax_block(iatom_block,jatom_block) = MAXVAL(pmax_atom(iatom_start:iatom_end,jatom_start:jatom_end))
        END DO
      END DO
    END IF

    SELECT CASE(eval_type)
      CASE (hfx_do_eval_energy)
        atomic_pair_list => x_data%atomic_pair_list
      CASE (hfx_do_eval_forces)
        atomic_pair_list => x_data%atomic_pair_list_forces
    END SELECT
    CALL build_atomic_pair_list(natom, atomic_pair_list, kind_of,basis_parameter,particle_set, &
                                do_periodic,coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell,&
                                x_data%blocks)

!$OMP END MASTER
!$OMP BARRIER

    !! If there is only 1 cpu skip the binning
    IF(n_processes == 1) THEN
      ALLOCATE(tmp_dist(1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1)%number_of_atom_quartets = HUGE(tmp_dist(1)%number_of_atom_quartets)
      tmp_dist(1)%istart = 0_int_8
      ptr_to_tmp_dist => tmp_dist(:)
      SELECT CASE( eval_type)
      CASE (hfx_do_eval_energy)
        CALL hfx_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)
      CASE (hfx_do_eval_forces)
        CALL hfx_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)
      END SELECT
      DEALLOCATE(tmp_dist,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      !! Calculate total numbers of integrals that have to be calculated (wrt screening and symmetry)
!$OMP BARRIER
!$OMP MASTER
    CALL timeset(routineN//"_count",handle_inner)
!$OMP END MASTER
!$OMP BARRIER

      cost_per_core = 0_int_8
      my_process_id = para_env%mepos*n_threads + i_thread
      nblocks = load_balance_parameter%nblocks

      DO atom_block = my_process_id,INT(nblocks,KIND=int_8)**4-1,n_processes

        latom_block = INT(MODULO(atom_block,INT(nblocks,KIND=int_8)))+1
        tmp_block = atom_block/nblocks
        katom_block = INT(MODULO(tmp_block,INT(nblocks,KIND=int_8)))+1
        IF(latom_block<katom_block) CYCLE
        tmp_block = tmp_block/nblocks
        jatom_block = INT(MODULO(tmp_block,INT(nblocks,KIND=int_8)))+1
        tmp_block = tmp_block/nblocks
        iatom_block = INT(MODULO(tmp_block,INT(nblocks,KIND=int_8)))+1
        IF(jatom_block<iatom_block) CYCLE

        iatom_start = x_data%blocks(iatom_block)%istart
        iatom_end =   x_data%blocks(iatom_block)%iend
        jatom_start = x_data%blocks(jatom_block)%istart
        jatom_end =   x_data%blocks(jatom_block)%iend
        katom_start = x_data%blocks(katom_block)%istart
        katom_end =   x_data%blocks(katom_block)%iend
        latom_start = x_data%blocks(latom_block)%istart
        latom_end =   x_data%blocks(latom_block)%iend

        SELECT CASE(eval_type)
        CASE (hfx_do_eval_energy)
          pmax_blocks = MAX( pmax_block(katom_block, iatom_block), &
                             pmax_block(latom_block, jatom_block), &
                             pmax_block(latom_block, iatom_block), &
                             pmax_block(katom_block, jatom_block) )
        CASE (hfx_do_eval_forces)
          pmax_blocks = MAX( pmax_block(katom_block, iatom_block) + &
                             pmax_block(latom_block, jatom_block), &
                             pmax_block(latom_block, iatom_block) + &
                             pmax_block(katom_block, jatom_block) )
        END SELECT

        IF( 2.0_dp*coeffs_kind_max0 + pmax_blocks  < log10_eps_schwarz ) CYCLE

        cost_per_core = cost_per_core &
                        + estimate_block_cost(natom, nkind, list_ij, list_kl, set_list_ij, set_list_kl, &
                                              iatom_start, iatom_end, jatom_start, jatom_end, &
                                              katom_start, katom_end, latom_start, latom_end, &
                                              particle_set,&
                                              coeffs_set, coeffs_kind, &
                                              is_assoc_atomic_block_global, do_periodic,&
                                              kind_of, basis_parameter, pmax_set, pmax_atom, pmax_blocks, &
                                              cell, &
                                              do_p_screening, map_atom_to_kind_atom, eval_type,&
                                              log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, atomic_pair_list)

      END DO ! atom_block

      nbins = load_balance_parameter%nbins
      cost_per_bin = (cost_per_core + nbins -1 ) / (nbins)

!$OMP BARRIER
!$OMP MASTER
    CALL timestop(handle_inner)
!$OMP END MASTER
!$OMP BARRIER

! new load balancing test
IF ( .FALSE. ) THEN
    CALL hfx_recursive_load_balance( n_processes, my_process_id, nblocks, &
                                 natom, nkind, list_ij, list_kl, set_list_ij, set_list_kl, &
                                 particle_set,&
                                 coeffs_set, coeffs_kind, &
                                 is_assoc_atomic_block_global, do_periodic,&
                                 kind_of, basis_parameter, pmax_set, pmax_atom, pmax_blocks, &
                                 cell, x_data, para_env, pmax_block,&
                                 do_p_screening, map_atom_to_kind_atom, eval_type,&
                                 log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, atomic_pair_list)
END IF

!$OMP BARRIER
!$OMP MASTER
    CALL timeset(routineN//"_bin",handle_inner)
!$OMP END MASTER
!$OMP BARRIER

      ALLOCATE(binned_dist(nbins), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      binned_dist(:)%istart = -1_int_8
      binned_dist(:)%number_of_atom_quartets = 0_int_8
      binned_dist(:)%cost = 0_int_8
      binned_dist(:)%time_first_scf = 0.0_dp
      binned_dist(:)%time_other_scf =0.0_dp
      binned_dist(:)%time_forces = 0.0_dp

      current_cost = 0
      mepos = 1
      distribution_counter_start = 1
      distribution_counter_end = 0
      ibin = 1

      global_quartet_counter = 0
      local_quartet_counter = 0
      last_bin_needs_to_be_filled = .FALSE.
      DO atom_block = my_process_id,INT(nblocks,KIND=int_8)**4-1,n_processes
        latom_block = INT(MODULO(atom_block,INT(nblocks,KIND=int_8)))+1
        tmp_block = atom_block/nblocks
        katom_block = INT(MODULO(tmp_block,INT(nblocks,KIND=int_8)))+1
        IF(latom_block<katom_block) CYCLE
        tmp_block = tmp_block/nblocks
        jatom_block = INT(MODULO(tmp_block,INT(nblocks,KIND=int_8)))+1
        tmp_block = tmp_block/nblocks
        iatom_block = INT(MODULO(tmp_block,INT(nblocks,KIND=int_8)))+1
        IF(jatom_block<iatom_block) CYCLE

        distribution_counter_end = distribution_counter_end + 1
        global_quartet_counter = global_quartet_counter + 1
        last_bin_needs_to_be_filled = .TRUE.

        IF(binned_dist(ibin)%istart == -1_int_8 ) binned_dist(ibin)%istart = atom_block

        iatom_start = x_data%blocks(iatom_block)%istart
        iatom_end =   x_data%blocks(iatom_block)%iend
        jatom_start = x_data%blocks(jatom_block)%istart
        jatom_end =   x_data%blocks(jatom_block)%iend
        katom_start = x_data%blocks(katom_block)%istart
        katom_end =   x_data%blocks(katom_block)%iend
        latom_start = x_data%blocks(latom_block)%istart
        latom_end =   x_data%blocks(latom_block)%iend

        SELECT CASE(eval_type)
        CASE (hfx_do_eval_energy)
          pmax_blocks = MAX( pmax_block(katom_block, iatom_block), &
                             pmax_block(latom_block, jatom_block), &
                             pmax_block(latom_block, iatom_block), &
                             pmax_block(katom_block, jatom_block) )
        CASE (hfx_do_eval_forces)
          pmax_blocks = MAX( pmax_block(katom_block, iatom_block) + &
                             pmax_block(latom_block, jatom_block), &
                             pmax_block(latom_block, iatom_block) + &
                             pmax_block(katom_block, jatom_block) )
        END SELECT


        IF( 2.0_dp*coeffs_kind_max0 + pmax_blocks  < log10_eps_schwarz ) CYCLE

        current_cost = current_cost &
                        + estimate_block_cost(natom, nkind, list_ij, list_kl, set_list_ij, set_list_kl, &
                                              iatom_start, iatom_end, jatom_start, jatom_end, &
                                              katom_start, katom_end, latom_start, latom_end, &
                                              particle_set,&
                                              coeffs_set, coeffs_kind, &
                                              is_assoc_atomic_block_global, do_periodic,&
                                              kind_of, basis_parameter, pmax_set, pmax_atom, pmax_blocks, &
                                              cell, &
                                              do_p_screening, map_atom_to_kind_atom, eval_type,&
                                              log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, atomic_pair_list)


        IF( current_cost >= cost_per_bin) THEN
          IF( ibin==nbins) THEN
            binned_dist(ibin)%number_of_atom_quartets = binned_dist(ibin)%number_of_atom_quartets + &
                               distribution_counter_end - distribution_counter_start + 1
          ELSE
            binned_dist(ibin)%number_of_atom_quartets = distribution_counter_end - distribution_counter_start + 1
          END IF
          binned_dist(ibin)%cost = binned_dist(ibin)%cost + current_cost
          ibin = MIN(ibin + 1,nbins)
          distribution_counter_start = distribution_counter_end + 1
          current_cost = 0
          last_bin_needs_to_be_filled = .FALSE.
        END IF
      ENDDO

!$OMP BARRIER
!$OMP MASTER
    CALL timestop(handle_inner)
    CALL timeset(routineN//"_dist",handle_inner)
!$OMP END MASTER
!$OMP BARRIER
      !! Fill the last bin if necessary
      IF( last_bin_needs_to_be_filled ) THEN
        binned_dist(ibin)%cost = binned_dist(ibin)%cost + current_cost
        IF( ibin == nbins ) THEN
          binned_dist(ibin)%number_of_atom_quartets = binned_dist(ibin)%number_of_atom_quartets + &
                                                   distribution_counter_end-distribution_counter_start+1
        ELSE
          binned_dist(ibin)%number_of_atom_quartets = distribution_counter_end-distribution_counter_start+1
        END IF
      END IF

      !! Sanity-Check
      DO ibin=1,nbins
        local_quartet_counter = local_quartet_counter + binned_dist(ibin)%number_of_atom_quartets
      END DO
!$OMP BARRIER
!$OMP MASTER
      shm_local_quartet_counter = 0
      shm_global_quartet_counter = 0
!$OMP END MASTER
!$OMP BARRIER
!$OMP ATOMIC
      shm_local_quartet_counter = shm_local_quartet_counter + local_quartet_counter
!$OMP ATOMIC
      shm_global_quartet_counter = shm_global_quartet_counter + global_quartet_counter

!$OMP BARRIER
!$OMP MASTER
      CALL mp_sum(shm_local_quartet_counter, para_env%group)
      CALL mp_sum(shm_global_quartet_counter, para_env%group)
      IF(para_env%mepos == 0 ) THEN
        IF( shm_local_quartet_counter /= shm_global_quartet_counter ) THEN
          WRITE(error_msg,'(A,I0,A,I0,A)') "HFX Sanity check for parallel distribution failed. "//&
                "Number of local quartets (", shm_local_quartet_counter, &
                ") and number of global quartets (", shm_global_quartet_counter,&
                ") are different. Please send in a bug report."
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                         error_msg,&
                         error,&
                         only_ionode=.TRUE.)
        END IF
      END IF
!$OMP END MASTER

!$OMP BARRIER
!$OMP MASTER
      ALLOCATE(cost_matrix(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      cost_matrix = 0
!$OMP END MASTER
!$OMP BARRIER
      icpu = para_env%mepos+1
      DO i=1,nbins
        cost_matrix((icpu-1)*nbins*n_threads +i_thread*nbins + i) = binned_dist(i)%cost
      END DO
      mepos = para_env%mepos
!$OMP BARRIER

!$OMP MASTER
      ! sync before/after ring of isendrecv
      CALL mp_sync(para_env%group)

      ALLOCATE(sendbuffer(nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(recbuffer(nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      sendbuffer = cost_matrix(mepos*nbins*n_threads+1:mepos*nbins*n_threads+nbins*n_threads)

      dest  =MODULO(mepos+1,ncpu)
      source=MODULO(mepos-1,ncpu)
      DO icpu = 0, ncpu-1
        IF (icpu .NE. ncpu-1) THEN
          CALL mp_isendrecv(sendbuffer,dest,recbuffer,source, &
                            para_env%group,req(1),req(2),13)
        ENDIF
        data_from = MODULO(mepos - icpu,ncpu)
        cost_matrix(data_from*nbins*n_threads+1:data_from*nbins*n_threads+nbins*n_threads) = sendbuffer
        IF (icpu .NE. ncpu -1) THEN
          CALL mp_waitall(req)
        ENDIF
        swapbuffer=>sendbuffer
        sendbuffer=>recbuffer
        recbuffer=>swapbuffer
      ENDDO
      DEALLOCATE(recbuffer, sendbuffer, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER

!$OMP BARRIER
!$OMP MASTER
    CALL timestop(handle_inner)
    CALL timeset(routineN//"_opt",handle_inner)
!$OMP END MASTER
!$OMP BARRIER

      !! Find an optimal distribution i.e. assign each element of the cost matrix to a certain process
!$OMP BARRIER
      ALLOCATE(local_cost_matrix(SIZE(cost_matrix,1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      local_cost_matrix = cost_matrix
!$OMP MASTER
      ALLOCATE(shm_distribution_vector(ncpu*nbins*n_threads), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      CALL optimize_distribution(ncpu*nbins*n_threads, ncpu*n_threads, local_cost_matrix, &
                                 shm_distribution_vector, x_data%load_balance_parameter%do_randomize, &
                                 error)

      CALL timestop(handle_inner)
      CALL timeset(routineN//"_redist",handle_inner)
      !! Collect local data to global array
      ALLOCATE(full_dist(ncpu*n_threads,nbins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      full_dist(:,:)%istart = 0_int_8
      full_dist(:,:)%number_of_atom_quartets = 0_int_8
      full_dist(:,:)%cost = 0_int_8
      full_dist(:,:)%time_first_scf = 0.0_dp
      full_dist(:,:)%time_other_scf = 0.0_dp
      full_dist(:,:)%time_forces = 0.0_dp
!$OMP END MASTER
!$OMP BARRIER
      mepos = para_env%mepos + 1
      full_dist((mepos-1)*n_threads+i_thread+1,:) = binned_dist(:)

!$OMP BARRIER
!$OMP MASTER
      ALLOCATE(sendbuffer(3*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(recbuffer(3*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      mepos = para_env%mepos
      DO j=1,n_threads
        DO i=1,nbins
          sendbuffer((j-1)*3*nbins+(i-1)*3+1) = full_dist(mepos*n_threads+j,i)%istart
          sendbuffer((j-1)*3*nbins+(i-1)*3+2) = full_dist(mepos*n_threads+j,i)%number_of_atom_quartets
          sendbuffer((j-1)*3*nbins+(i-1)*3+3) = full_dist(mepos*n_threads+j,i)%cost
        END DO
      END DO

      ! sync before/after ring of isendrecv
      CALL mp_sync(para_env%group)
      dest   = MODULO(mepos+1,ncpu)
      source = MODULO(mepos-1,ncpu)
      DO icpu = 0, ncpu-1
        IF (icpu .NE. ncpu-1) THEN
          CALL mp_isendrecv(sendbuffer,dest,recbuffer,source, &
                            para_env%group,req(1),req(2),13)
        ENDIF
        data_from = MODULO(mepos - icpu,ncpu)
        DO j=1,n_threads
          DO i=1,nbins
            full_dist(data_from*n_threads+j,i)%istart = sendbuffer((j-1)*3*nbins+(i-1)*3+1)
            full_dist(data_from*n_threads+j,i)%number_of_atom_quartets = sendbuffer((j-1)*3*nbins+(i-1)*3+2)
            full_dist(data_from*n_threads+j,i)%cost = sendbuffer((j-1)*3*nbins+(i-1)*3+3)
          END DO
        END DO

        IF (icpu .NE. ncpu -1) THEN
          CALL mp_waitall(req)
        ENDIF
        swapbuffer=>sendbuffer
        sendbuffer=>recbuffer
        recbuffer=>swapbuffer
      ENDDO
      DEALLOCATE(recbuffer, sendbuffer, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      ! sync before/after ring of isendrecv
      CALL mp_sync(para_env%group)
!$OMP END MASTER
!$OMP BARRIER
      !! reorder the distribution according to the distribution vector
      ALLOCATE(tmp_pos(ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_pos = 1
      ALLOCATE(tmp_dist(nbins*ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      tmp_dist(:)%istart = 0_int_8
      tmp_dist(:)%number_of_atom_quartets = 0_int_8
      tmp_dist(:)%cost = 0_int_8
      tmp_dist(:)%time_first_scf = 0.0_dp
      tmp_dist(:)%time_other_scf = 0.0_dp
      tmp_dist(:)%time_forces = 0.0_dp

      DO icpu= 1,n_processes
        DO i=1,nbins
            mepos = my_process_id + 1
            IF(shm_distribution_vector((icpu-1)*nbins + i) == mepos) THEN
              tmp_dist(tmp_pos(mepos)) = full_dist(icpu,i)
              tmp_pos(mepos) = tmp_pos(mepos) + 1
            END IF
        END DO
      END DO

      !! Assign the load to each process
      NULLIFY(ptr_to_tmp_dist)
      mepos = my_process_id+1
      ptr_to_tmp_dist => tmp_dist(1:tmp_pos(mepos)-1)
      SELECT CASE (eval_type)
      CASE (hfx_do_eval_energy)
        CALL hfx_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)
      CASE (hfx_do_eval_forces)
        CALL hfx_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)
      END SELECT

!$OMP BARRIER
!$OMP MASTER
      DEALLOCATE(full_dist,cost_matrix,shm_distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER
      DEALLOCATE(tmp_dist, tmp_pos, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(binned_dist, local_cost_matrix, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(set_list_ij, set_list_kl, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

!$OMP BARRIER
!$OMP MASTER
    CALL timestop(handle_inner)
!$OMP END MASTER
!$OMP BARRIER
    END IF
!$OMP BARRIER
!$OMP MASTER
    CALL timestop(handle)
!$OMP END MASTER
!$OMP BARRIER
  END SUBROUTINE hfx_load_balance

! *****************************************************************************
!> \brief Reference implementation of new recursive load balancing routine
!>        Computes a local list of atom_blocks (p_atom_blocks,q_atom_blocks) for
!>        each process in a P-Q grid such that every process has more or less the
!>        same amount of work. Has no output at the moment (not used) but writes
!>        its computed load balance values into a file. Possible output is ready
!>        to use in the two arrays p_atom_blocks & q_atom_blocks
!> \param n_processes ...
!> \param my_process_id ...
!> \param nblocks ...
!> \param natom ...
!> \param nkind ...
!> \param list_ij ...
!> \param list_kl ...
!> \param set_list_ij ...
!> \param set_list_kl ...
!> \param particle_set ...
!> \param coeffs_set ...
!> \param coeffs_kind ...
!> \param is_assoc_atomic_block_global ...
!> \param do_periodic ...
!> \param kind_of ...
!> \param basis_parameter ...
!> \param pmax_set ...
!> \param pmax_atom ...
!> \param pmax_blocks ...
!> \param cell ...
!> \param x_data ...
!> \param para_env ...
!> \param pmax_block ...
!> \param do_p_screening ...
!> \param map_atom_to_kind_atom ...
!> \param eval_type ...
!> \param log10_eps_schwarz ...
!> \param log_2 ...
!> \param coeffs_kind_max0 ...
!> \param use_virial ...
!> \param atomic_pair_list ...
!> \par History
!>      03.2011 created [Michael Steinlechner]
!> \author Michael Steinlechner
! *****************************************************************************

SUBROUTINE hfx_recursive_load_balance( n_processes, my_process_id, nblocks, &
                             natom, nkind, list_ij, list_kl, set_list_ij, set_list_kl, &
                             particle_set,&
                             coeffs_set, coeffs_kind, &
                             is_assoc_atomic_block_global, do_periodic,&
                             kind_of, basis_parameter, pmax_set, pmax_atom, pmax_blocks, &
                             cell, x_data, para_env, pmax_block,&
                             do_p_screening, map_atom_to_kind_atom, eval_type,&
                             log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, atomic_pair_list)

! input variables:
    INTEGER, INTENT(IN)                      :: n_processes, my_process_id, &
                                                nblocks, natom, nkind
    TYPE(pair_list_type), INTENT(IN)         :: list_ij, list_kl
    TYPE(pair_set_list_type), ALLOCATABLE, &
      DIMENSION(:), INTENT(IN)               :: set_list_ij, set_list_kl
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN), POINTER                    :: particle_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), INTENT(IN), &
      POINTER                                :: coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), INTENT(IN), POINTER   :: coeffs_kind
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: is_assoc_atomic_block_global
    LOGICAL, INTENT(IN)                      :: do_periodic
    INTEGER, INTENT(IN)                      :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      INTENT(IN), POINTER                    :: basis_parameter
    TYPE(hfx_p_kind), DIMENSION(:), &
      INTENT(IN), POINTER                    :: pmax_set
    REAL(dp), DIMENSION(:, :), INTENT(IN), &
      POINTER                                :: pmax_atom
    REAL(dp)                                 :: pmax_blocks
    TYPE(cell_type), INTENT(IN), POINTER     :: cell
    TYPE(hfx_type), INTENT(IN), POINTER      :: x_data
    TYPE(cp_para_env_type), INTENT(IN), &
      POINTER                                :: para_env
    REAL(dp), DIMENSION(:, :), INTENT(IN), &
      POINTER                                :: pmax_block
    LOGICAL, INTENT(IN)                      :: do_p_screening
    INTEGER, DIMENSION(:), INTENT(IN), &
      POINTER                                :: map_atom_to_kind_atom
    INTEGER, INTENT(IN)                      :: eval_type
    REAL(dp), INTENT(IN)                     :: log10_eps_schwarz, log_2, &
                                                coeffs_kind_max0
    LOGICAL, INTENT(IN)                      :: use_virial
    LOGICAL, DIMENSION(:, :), INTENT(IN), &
      POINTER                                :: atomic_pair_list

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_recursive_load_balance', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom_block, iatom_end, iatom_start, j, &
      jatom_block, jatom_end, jatom_start, katom_block, katom_end, &
      katom_start, latom_block, latom_end, latom_start, nP, nQ, numBins, p, &
      q, sizeP, sizeQ, unit_nr
    INTEGER(int_8)                           :: local_cost, pidx, qidx, sumP, &
                                                sumQ
    INTEGER(int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: local_cost_vector
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: blocksize, p_atom_blocks, &
                                                permute, q_atom_blocks
    REAL(dp)                                 :: maximum, mean

! internal variables:

!$OMP BARRIER
!$OMP MASTER
    CALL timeset(routineN,handle)
!$OMP END MASTER
!$OMP BARRIER


    ! calculate best p/q distribution grid for the n_processes
    CALL hfx_calculate_PQ(p, q, numBins, n_processes)

    ALLOCATE(blocksize(numBins))
    ALLOCATE(permute(nblocks**2))
    DO i=1,nblocks**2
        permute(i) = i
    END DO

    ! call the main recursive permutation routine.
    ! Output:
    !   blocksize :: vector (size numBins) with the sizes for each column/row block
    !   permute   :: permutation vector
    CALL hfx_recursive_permute( blocksize, 1, nblocks**2, numBins,&
                                permute, 1,&
                                my_process_id, n_processes, nblocks, &
                                natom, nkind, list_ij, list_kl, set_list_ij, set_list_kl, &
                                particle_set,&
                                coeffs_set, coeffs_kind, &
                                is_assoc_atomic_block_global, do_periodic,&
                                kind_of, basis_parameter, pmax_set, pmax_atom, pmax_blocks, &
                                cell, x_data, para_env, pmax_block,&
                                do_p_screening, map_atom_to_kind_atom, eval_type,&
                                log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, atomic_pair_list)

    ! number of blocks per processor in p-direction (vertical)
    nP = numBins/p
    ! number of blocks per processor in q-direction (horizontal)
    nQ = numBins/q

    ! calc own position in P-Q-processor grid (PQ-grid is column-major)
    pidx = MODULO(INT(my_process_id),INT(p)) + 1
    qidx = my_process_id / p + 1

    sizeP = SUM(blocksize( (nP*(pidx-1)+1):(nP*pidx) ))
    sizeQ = SUM(blocksize( (nQ*(qidx-1)+1):(nQ*qidx) ))

    sumP = SUM(blocksize( 1:(nP*(pidx-1)) ))
    sumQ = SUM(blocksize( 1:(nQ*(qidx-1)) ))

    ALLOCATE(p_atom_blocks(sizeP))
    ALLOCATE(q_atom_blocks(sizeQ))

    p_atom_blocks(:) = permute( (sumP+1):(sumP+sizeP) )
    q_atom_blocks(:) = permute( (sumQ+1):(sumQ+sizeQ) )

    ! from here on, we are actually finished, each process has been
    ! assigned a (p_atom_blocks,q_atom_blocks) pair list.
    ! what follows is just a small routine to calculate the local cost
    ! for each processor which is then written to a file.

    ! calculate local cost for each processor!
    ! ****************************************
    local_cost = 0
    DO i = 1,sizeP
        DO j= 1,sizeQ

    !       get corresponding 4D block indices out of our own P-Q-block
            latom_block = MODULO(q_atom_blocks(j),nblocks)
            iatom_block = q_atom_blocks(j)/nblocks + 1
            jatom_block = MODULO(p_atom_blocks(i),nblocks)
            katom_block = p_atom_blocks(i)/nblocks + 1

    !       symmetry checks.
            IF(latom_block<katom_block) CYCLE
            IF(jatom_block<iatom_block) CYCLE

            iatom_start = x_data%blocks(iatom_block)%istart
            iatom_end =   x_data%blocks(iatom_block)%iend
            jatom_start = x_data%blocks(jatom_block)%istart
            jatom_end =   x_data%blocks(jatom_block)%iend
            katom_start = x_data%blocks(katom_block)%istart
            katom_end =   x_data%blocks(katom_block)%iend
            latom_start = x_data%blocks(latom_block)%istart
            latom_end =   x_data%blocks(latom_block)%iend

    !       whatever.
            SELECT CASE(eval_type)
            CASE (hfx_do_eval_energy)
              pmax_blocks = MAX( pmax_block(katom_block, iatom_block), &
                                 pmax_block(latom_block, jatom_block), &
                                 pmax_block(latom_block, iatom_block), &
                                 pmax_block(katom_block, jatom_block) )
            CASE (hfx_do_eval_forces)
              pmax_blocks = MAX( pmax_block(katom_block, iatom_block) + &
                                 pmax_block(latom_block, jatom_block), &
                                 pmax_block(latom_block, iatom_block) + &
                                 pmax_block(katom_block, jatom_block) )
            END SELECT

    !       screening.
            IF( 2.0_dp*coeffs_kind_max0 + pmax_blocks  < log10_eps_schwarz ) CYCLE

    !       estimate the cost of this atom_block.
            local_cost = local_cost + estimate_block_cost( natom, nkind, list_ij, list_kl, set_list_ij, &
                                                  set_list_kl, &
                                                  iatom_start, iatom_end, jatom_start, jatom_end, &
                                                  katom_start, katom_end, latom_start, latom_end, &
                                                  particle_set,&
                                                  coeffs_set, coeffs_kind, &
                                                  is_assoc_atomic_block_global, do_periodic,&
                                                  kind_of, basis_parameter, pmax_set, pmax_atom, pmax_blocks, &
                                                  cell, &
                                                  do_p_screening, map_atom_to_kind_atom, eval_type,&
                                                  log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, atomic_pair_list)
          END DO
    END DO

    ALLOCATE(local_cost_vector(n_processes))
    local_cost_vector = 0
    local_cost_vector(my_process_id+1) = local_cost
    CALL mp_sum(local_cost_vector, para_env%group)

    mean = SUM(local_cost_vector) / n_processes
    maximum = MAXVAL(local_cost_vector)

    !$OMP BARRIER
    !$OMP MASTER
        ! only output once
        IF( my_process_id == 0) THEN
            CALL open_file(unit_number=unit_nr,file_name="loads.dat")
            WRITE(unit_nr,*) 'maximum cost:', maximum
            WRITE(unit_nr,*) 'mean cost:', mean
            WRITE(unit_nr,*) 'load balance ratio max/mean: ', maximum/mean
            WRITE(unit_nr,*) '-------- detailed per-process costs ---------'
            DO i=1,n_processes
                WRITE(unit_nr,*) local_cost_vector(i)
            END DO
            CALL close_file(unit_nr)
        END IF
    !$OMP END MASTER
    !$OMP BARRIER

    DEALLOCATE(local_cost_vector)
    DEALLOCATE(p_atom_blocks,q_atom_blocks)
    DEALLOCATE(blocksize,permute)

!$OMP BARRIER
!$OMP MASTER
    CALL timestop(handle)
!$OMP END MASTER
!$OMP BARRIER

END SUBROUTINE hfx_recursive_load_balance


! *****************************************************************************
!> \brief Small routine to calculate the optimal P-Q-processor grid distribution
!>        for a given number of processors N
!>        and the corresponding number of Bins for the load balancing routine
!> \param p     number of rows on P-Q process grid (output)
!> \param q     number of columns on P-Q process grid (output)
!> \param nBins number of Bins (output)
!> \param N     number of processes (input)
!> \par History
!>      03.2011 created [Michael Steinlechner]
!> \author Michael Steinlechner
! *****************************************************************************
  SUBROUTINE hfx_calculate_PQ(p, q, nBins, N)

    INTEGER, INTENT(OUT)                     :: p, q, nBins
    INTEGER, INTENT(IN)                      :: N

    INTEGER                                  :: a, b, k
    REAL(dp)                                 :: sqN

    k=2
    sqN = SQRT(REAL(N,KIND=dp))
    p = 1

    DO WHILE (REAL(k,KIND=dp) <= sqN)
        IF (MODULO(N,k) == 0) THEN
            p = k
        END IF
        k = k+1
    END DO
    q=N/p

    ! now compute the least common multiple of p & q to get the number of necessary bins
    ! compute using the relation LCM(p,q) = abs(p*q) / GCD(p,q)
    ! and use euclid's algorithm for GCD computation.
    a = p
    b = q

    DO WHILE( b .NE. 0)
        IF (a > b) THEN
            a = a-b
        ELSE
            b = b-a
        END IF
    END DO
    ! gcd(p,q) is now saved in a

    nBins = p*q / a

  END SUBROUTINE


! *****************************************************************************
!> \brief Recursive permutation routine for the load balancing of the integral
!>       computation
!> \param blocksize     vector of blocksizes, size(nProc), which contains for
!>                      each process the local blocksize (OUTPUT)
!> \param blockstart    starting row/column idx of the block which is to be examined
!>                      at this point (INPUT)
!> \param blockend      ending row/column idx of the block which is to be examined
!>                      (INPUT)
!> \param nProc_in      number of bins into which the current block has to be divided
!>                      (INPUT)
!> \param permute       permutation vector which balances column/row cost
!>                      size(nblocks^2). (OUTPUT)
!> \param step ...
!> \param my_process_id ...
!> \param n_processes ...
!> \param nblocks ...
!> \param natom ...
!> \param nkind ...
!> \param list_ij ...
!> \param list_kl ...
!> \param set_list_ij ...
!> \param set_list_kl ...
!> \param particle_set ...
!> \param coeffs_set ...
!> \param coeffs_kind ...
!> \param is_assoc_atomic_block_global ...
!> \param do_periodic ...
!> \param kind_of ...
!> \param basis_parameter ...
!> \param pmax_set ...
!> \param pmax_atom ...
!> \param pmax_blocks ...
!> \param cell ...
!> \param x_data ...
!> \param para_env ...
!> \param pmax_block ...
!> \param do_p_screening ...
!> \param map_atom_to_kind_atom ...
!> \param eval_type ...
!> \param log10_eps_schwarz ...
!> \param log_2 ...
!> \param coeffs_kind_max0 ...
!> \param use_virial ...
!> \param atomic_pair_list ...
!> \par History
!>      03.2011 created [Michael Steinlechner]
!> \author Michael Steinlechner
! *****************************************************************************
  RECURSIVE SUBROUTINE hfx_recursive_permute( blocksize, blockstart, blockend, nProc_in,&
                                      permute, step,&
                                      my_process_id, n_processes, nblocks, &
                                      natom, nkind, list_ij, list_kl, set_list_ij, set_list_kl, &
                                      particle_set,&
                                      coeffs_set, coeffs_kind, &
                                      is_assoc_atomic_block_global, do_periodic,&
                                      kind_of, basis_parameter, pmax_set, pmax_atom, pmax_blocks, &
                                      cell, x_data, para_env, pmax_block,&
                                      do_p_screening, map_atom_to_kind_atom, eval_type,&
                                      log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, atomic_pair_list)


    INTEGER                                  :: nProc_in, blockend, blockstart
    INTEGER, DIMENSION(nProc_in)             :: blocksize
    INTEGER                                  :: nblocks, n_processes, &
                                                my_process_id
    INTEGER, INTENT(IN)                      :: step
    INTEGER, DIMENSION(nblocks*nblocks)      :: permute
    INTEGER                                  :: natom
    INTEGER, INTENT(IN)                      :: nkind
    TYPE(pair_list_type)                     :: list_ij, list_kl
    TYPE(pair_set_list_type), ALLOCATABLE, &
      DIMENSION(:)                           :: set_list_ij, set_list_kl
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: coeffs_kind
    INTEGER, DIMENSION(:, :)                 :: is_assoc_atomic_block_global
    LOGICAL                                  :: do_periodic
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_p_kind), DIMENSION(:), POINTER  :: pmax_set
    REAL(dp), DIMENSION(:, :), POINTER       :: pmax_atom
    REAL(dp)                                 :: pmax_blocks
    TYPE(cell_type), POINTER                 :: cell
    TYPE(hfx_type), POINTER                  :: x_data
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(dp), DIMENSION(:, :), POINTER       :: pmax_block
    LOGICAL, INTENT(IN)                      :: do_p_screening
    INTEGER, DIMENSION(:), POINTER           :: map_atom_to_kind_atom
    INTEGER, INTENT(IN)                      :: eval_type
    REAL(dp)                                 :: log10_eps_schwarz, log_2, &
                                                coeffs_kind_max0
    LOGICAL, INTENT(IN)                      :: use_virial
    LOGICAL, DIMENSION(:, :), POINTER        :: atomic_pair_list

    INTEGER :: col, endoffset, i, iatom_block, iatom_end, iatom_start, idx, &
      inv_perm, jatom_block, jatom_end, jatom_start, katom_block, katom_end, &
      katom_start, latom_block, latom_end, latom_start, nbins, nProc, row, &
      startoffset
    INTEGER(int_8)                           :: atom_block, tmp_block
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ithblocksize, localblocksize
    INTEGER, &
      DIMENSION(blockend-blockstart+1)       :: bin_perm, tmp_perm
    REAL(dp)                                 :: partialcost
    REAL(dp), DIMENSION(nblocks*nblocks)     :: cost_vector

    nProc = nProc_in
    cost_vector = 0.0_dp

!   loop over local atom_blocks.
    DO atom_block = my_process_id,INT(nblocks,KIND=int_8)**4-1,n_processes

!       get corresponding 4D block indices
        latom_block = INT(MODULO(atom_block,INT(nblocks,KIND=int_8)))+1
        tmp_block = atom_block/nblocks
        katom_block = INT(MODULO(tmp_block,INT(nblocks,KIND=int_8)))+1
        IF(latom_block<katom_block) CYCLE
        tmp_block = tmp_block/nblocks
        jatom_block = INT(MODULO(tmp_block,INT(nblocks,KIND=int_8)))+1
        tmp_block = tmp_block/nblocks
        iatom_block = INT(MODULO(tmp_block,INT(nblocks,KIND=int_8)))+1
        IF(jatom_block<iatom_block) CYCLE

!       get 2D indices of this atom_block (with permutation applied)
!       for this, we need to invert the permutation, this means
!       find position in permutation vector where value==idx

        row = (katom_block-1)*nblocks + jatom_block
        inv_perm = 1
        DO WHILE (permute(inv_perm) .NE. row)
            inv_perm = inv_perm + 1
        END DO
        row = inv_perm

        col = (iatom_block-1)*nblocks + latom_block
        inv_perm = 1
        DO WHILE (permute(inv_perm) .NE. col)
            inv_perm = inv_perm + 1
        END DO
        col = inv_perm

!       if row/col outside our current diagonal block, skip calculation.
        IF (col < blockstart .OR. col > blockend) CYCLE
        IF (row < blockstart .OR. row > blockend) CYCLE

        iatom_start = x_data%blocks(iatom_block)%istart
        iatom_end =   x_data%blocks(iatom_block)%iend
        jatom_start = x_data%blocks(jatom_block)%istart
        jatom_end =   x_data%blocks(jatom_block)%iend
        katom_start = x_data%blocks(katom_block)%istart
        katom_end =   x_data%blocks(katom_block)%iend
        latom_start = x_data%blocks(latom_block)%istart
        latom_end =   x_data%blocks(latom_block)%iend

!       whatever.
        SELECT CASE(eval_type)
        CASE (hfx_do_eval_energy)
          pmax_blocks = MAX( pmax_block(katom_block, iatom_block), &
                             pmax_block(latom_block, jatom_block), &
                             pmax_block(latom_block, iatom_block), &
                             pmax_block(katom_block, jatom_block) )
        CASE (hfx_do_eval_forces)
          pmax_blocks = MAX( pmax_block(katom_block, iatom_block) + &
                             pmax_block(latom_block, jatom_block), &
                             pmax_block(latom_block, iatom_block) + &
                             pmax_block(katom_block, jatom_block) )
        END SELECT

!       screening.
        IF( 2.0_dp*coeffs_kind_max0 + pmax_blocks  < log10_eps_schwarz ) CYCLE

!       every second recursion step, compute row sum instead of column sum

        IF( MODULO(step,2) .EQ. 0) THEN
            idx = row
        ELSE
            idx = col
        END IF

!       estimate the cost of this atom_block.
        partialcost = estimate_block_cost( natom, nkind, list_ij, list_kl, set_list_ij, &
                                              set_list_kl, &
                                              iatom_start, iatom_end, jatom_start, jatom_end, &
                                              katom_start, katom_end, latom_start, latom_end, &
                                              particle_set,&
                                              coeffs_set, coeffs_kind, &
                                              is_assoc_atomic_block_global, do_periodic,&
                                              kind_of, basis_parameter, pmax_set, pmax_atom, pmax_blocks, &
                                              cell, &
                                              do_p_screening, map_atom_to_kind_atom, eval_type,&
                                              log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, atomic_pair_list)

        cost_vector(idx) = cost_vector(idx) + partialcost
      END DO ! atom_block

!   sum costvector over all processes
    CALL mp_sum(cost_vector, para_env%group)

!   calculate next prime factor of nProc
    nBins = 2
    DO WHILE (MODULO(INT(nProc),INT(nBins)) .NE. 0)
        nBins = nBins + 1
    END DO

    nProc = nProc/nBins

! ... do the binning...

    ALLOCATE(localblocksize(nBins))
    CALL hfx_permute_binning( nBins, cost_vector(blockstart:blockend), blockend-blockstart+1, bin_perm, localblocksize )

!... and update the permutation vector

    tmp_perm = permute(blockstart:blockend)
    permute(blockstart:blockend) = tmp_perm(bin_perm)

!   split recursion into the nBins Bins
    IF (nProc > 1) THEN
        ALLOCATE(ithblocksize(nProc))
        DO i=1,nBins
            startoffset = SUM(localblocksize(1:(i-1)))
            endoffset = SUM(localblocksize(1:i))-1

            CALL hfx_recursive_permute( ithblocksize, blockstart+startoffset, blockstart+endoffset, nProc,&
                                      permute, step+1, &
                                      my_process_id, n_processes, nblocks, &
                                      natom, nkind, list_ij, list_kl, set_list_ij, set_list_kl, &
                                      particle_set, &
                                      coeffs_set, coeffs_kind, &
                                      is_assoc_atomic_block_global, do_periodic,&
                                      kind_of, basis_parameter, pmax_set, pmax_atom, pmax_blocks, &
                                      cell, x_data,para_env, pmax_block,&
                                      do_p_screening, map_atom_to_kind_atom, eval_type,&
                                      log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, atomic_pair_list)
            blocksize(((i-1)*nProc+1):(i*nProc)) = ithblocksize
        END DO
        DEALLOCATE(ithblocksize)
    ELSE
        DO i=1,nBins
            blocksize(i) = localblocksize(i)
        END DO
    END IF

    DEALLOCATE(localblocksize)

END SUBROUTINE hfx_recursive_permute



! *****************************************************************************
!> \brief small binning routine for the recursive load balancing
!>
!> \param nBins         number of Bins (INPUT)
!> \param costvector    vector of current row/column costs which have to be binned (INPUT)
!> \param maxbinsize    upper bound for bin size (INPUT)
!> \param perm          resulting permutation due to be binning routine (OUTPUT)
!> \param block_count   vector of size(nbins) which contains the size of each bin (OUTPUT)
!> \par History
!>      03.2011 created [Michael Steinlechner]
!> \author Michael Steinlechner
! *****************************************************************************
SUBROUTINE hfx_permute_binning ( nBins, costvector, maxbinsize, perm, block_count)

    INTEGER, INTENT(IN)                      :: nBins, maxbinsize
    REAL(dp), DIMENSION(maxbinsize), &
      INTENT(IN)                             :: costvector
    INTEGER, DIMENSION(maxbinsize), &
      INTENT(OUT)                            :: perm
    INTEGER, DIMENSION(nBins), INTENT(OUT)   :: block_count

    INTEGER                                  :: i, j, mod_idx, offset
    INTEGER, DIMENSION(nBins, maxbinsize)    :: bin
    INTEGER, DIMENSION(nBins)                :: bin_idx
    INTEGER, DIMENSION(maxbinsize)           :: idx
    REAL(dp), DIMENSION(maxbinsize)          :: vec
    REAL(dp), DIMENSION(nBins)               :: bincosts

! be careful not to change costvector (copy it!)

    vec = costvector
    block_count = 0
    bincosts = 0

    !sort the array (ascending)
    CALL sort(vec, maxbinsize, idx)

    ! count the loop down to distribute the largest cols/rows first
    DO i = maxbinsize, 1, -1
        IF (vec(i) == 0) THEN
            ! spread zero-cost col/rows evenly among procs
            mod_idx = MODULO(i,nBins) + 1 !(note the fortran offset by one!)
            block_count(mod_idx) = block_count(mod_idx) + 1
            bin(mod_idx, block_count(mod_idx) ) = idx(i)
        ELSE
            ! sort the bins so that the one with the lowest cost is at the
            ! first place, where we then assign the current col/row
            CALL sort(bincosts, nBins, bin_idx)
            block_count = block_count(bin_idx)
            bin = bin(bin_idx,:)

            bincosts(1) = bincosts(1) + vec(i)
            block_count(1) = block_count(1) + 1
            bin(1, block_count(1) ) = idx(i)
        END IF
    END DO

    ! construct permutation vector from the binning
    offset = 0
    DO i = 1, nBins
        DO j = 1, block_count(i)
            perm(offset+j) = bin(i,j)
        END DO
        offset = offset + block_count(i)
    END DO


END SUBROUTINE hfx_permute_binning



! *****************************************************************************
!> \brief Cheap way of redistributing the eri's
!> \param x_data Object that stores the indices array
!> \param para_env para_env
!> \param load_balance_parameter contains parmameter for Monte-Carlo routines
!> \param i_thread current thread ID
!> \param n_threads Total Number of threads
!> \param eval_type ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2007 created [Manuel Guidon]
!>      02.2009 optimize Memory Usage [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      The cost matrix is given by the walltime for each bin that is measured
!>      during the calculation
! *****************************************************************************
  SUBROUTINE hfx_update_load_balance(x_data, para_env, &
                                load_balance_parameter, &
                                i_thread, n_threads, eval_type, error)

    TYPE(hfx_type), POINTER                  :: x_data
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(hfx_load_balance_type)              :: load_balance_parameter
    INTEGER, INTENT(IN)                      :: i_thread, n_threads, eval_type
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_update_load_balance', &
      routineP = moduleN//':'//routineN

    INTEGER :: data_from, dest, end_idx, handle, i, ibin, icpu, iprocess, j, &
      mepos, my_bin_size, my_global_start_idx, my_process_id, n_processes, &
      nbins, ncpu, req(2), source, start_idx, stat
    INTEGER(int_8), DIMENSION(:), POINTER    :: local_cost_matrix, recbuffer, &
                                                sendbuffer, swapbuffer
    INTEGER(int_8), DIMENSION(:), POINTER, &
      SAVE                                   :: cost_matrix
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_pos
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: bins_per_rank
    INTEGER, ALLOCATABLE, DIMENSION(:, :), &
      SAVE                                   :: bin_histogram
    INTEGER, DIMENSION(:), POINTER, SAVE     :: shm_distribution_vector
    INTEGER, SAVE                            :: max_bin_size
    LOGICAL                                  :: failure
    TYPE(hfx_distribution), DIMENSION(:), &
      POINTER                                :: binned_dist, ptr_to_tmp_dist, &
                                                tmp_dist
    TYPE(hfx_distribution), &
      DIMENSION(:, :), POINTER, SAVE         :: full_dist

    failure=.FALSE.

!$OMP BARRIER
!$OMP MASTER
    CALL timeset(routineN,handle)
!$OMP END MASTER
!$OMP BARRIER

    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads
    !! If there is only 1 cpu skip the binning
    IF(n_processes == 1) THEN
      ALLOCATE(tmp_dist(1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1)%number_of_atom_quartets = HUGE(tmp_dist(1)%number_of_atom_quartets)
      tmp_dist(1)%istart = 0_int_8
      ptr_to_tmp_dist => tmp_dist(:)
      SELECT CASE (eval_type)
      CASE (hfx_do_eval_energy)
        CALL hfx_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)
      CASE (hfx_do_eval_forces)
        CALL hfx_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)
      END SELECT
      DEALLOCATE(tmp_dist,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      mepos = para_env%mepos
      my_process_id = para_env%mepos*n_threads + i_thread
      nbins = load_balance_parameter%nbins
!$OMP MASTER
      ALLOCATE(bin_histogram(n_processes,2),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      bin_histogram = 0
!$OMP END MASTER
!$OMP BARRIER
      SELECT CASE (eval_type)
      CASE (hfx_do_eval_energy)
        my_bin_size = SIZE(x_data%distribution_energy)
      CASE (hfx_do_eval_forces)
        my_bin_size = SIZE(x_data%distribution_forces)
      END SELECT
      bin_histogram(my_process_id+1,1) = my_bin_size
!$OMP BARRIER
!$OMP MASTER
      CALL mp_sum(bin_histogram(:,1), para_env%group)
      bin_histogram(1,2) = bin_histogram(1,1)
      DO iprocess=2,n_processes
        bin_histogram(iprocess,2) = bin_histogram(iprocess-1,2) + bin_histogram(iprocess,1)
      END DO

      max_bin_size = MAXVAL(bin_histogram(para_env%mepos*n_threads+1:para_env%mepos*n_threads+n_threads,1))
      CALL mp_max(max_bin_size, para_env%group)
!$OMP END MASTER
!$OMP BARRIER
      ALLOCATE(binned_dist(my_bin_size), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      !! Use old binned_dist, but with timings cost
      SELECT CASE (eval_type)
      CASE (hfx_do_eval_energy)
        binned_dist = x_data%distribution_energy
      CASE (hfx_do_eval_forces)
        binned_dist = x_data%distribution_forces
      END SELECT

      DO ibin=1,my_bin_size
        IF(binned_dist(ibin)%number_of_atom_quartets == 0 ) THEN
          binned_dist(ibin)%cost = 0
        ELSE
          SELECT CASE (eval_type)
          CASE (hfx_do_eval_energy)
            IF(.NOT. load_balance_parameter%rtp_redistribute)THEN
               binned_dist(ibin)%cost = INT((binned_dist(ibin)%time_first_scf + &
                                        binned_dist(ibin)%time_other_scf)*10000.0_dp, int_8)
            ELSE
               binned_dist(ibin)%cost = INT((binned_dist(ibin)%time_other_scf)*10000.0_dp, int_8)
            END IF
          CASE (hfx_do_eval_forces)
            binned_dist(ibin)%cost = INT((binned_dist(ibin)%time_forces)*10000.0_dp,int_8)
          END SELECT
        END IF
      END DO
!$OMP BARRIER
!$OMP MASTER
      !! store all local results in a big cost matrix
      ALLOCATE(cost_matrix(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      cost_matrix = 0
      ALLOCATE(sendbuffer(max_bin_size*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(recbuffer(max_bin_size*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER
      my_global_start_idx = bin_histogram(my_process_id+1,2)- my_bin_size
      icpu = para_env%mepos+1
      DO i=1,my_bin_size
        cost_matrix(my_global_start_idx + i) = binned_dist(i)%cost
      END DO

      mepos = para_env%mepos
!$OMP BARRIER
!$OMP MASTER
      ALLOCATE(bins_per_rank(ncpu), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      bins_per_rank = 0
      DO icpu=1,ncpu
        bins_per_rank(icpu) = SUM(bin_histogram((icpu-1)*n_threads+1:(icpu-1)*n_threads+n_threads,1))
      END DO
      sendbuffer(1:bins_per_rank(para_env%mepos+1)) = &
        cost_matrix(my_global_start_idx+1:my_global_start_idx+bins_per_rank(para_env%mepos+1))

      dest  =MODULO(mepos+1,ncpu)
      source=MODULO(mepos-1,ncpu)
      ! sync before/after ring of isendrecv
      CALL mp_sync(para_env%group)
      DO icpu = 0, ncpu-1
        IF (icpu .NE. ncpu-1) THEN
          CALL mp_isendrecv(sendbuffer,dest,recbuffer,source, &
                            para_env%group,req(1),req(2),13)
        ENDIF
        data_from = MODULO(mepos - icpu,ncpu)
        start_idx = SUM(bins_per_rank(1:data_from+1)) - bins_per_rank(data_from+1) + 1
        end_idx = start_idx + bins_per_rank(data_from+1) - 1
        cost_matrix(start_idx:end_idx) = sendbuffer(1:end_idx-start_idx+1)

        IF (icpu .NE. ncpu -1) THEN
          CALL mp_waitall(req)
        ENDIF
        swapbuffer=>sendbuffer
        sendbuffer=>recbuffer
        recbuffer=>swapbuffer
      ENDDO
      DEALLOCATE(recbuffer, sendbuffer, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ! sync before/after ring of isendrecv
      CALL mp_sync(para_env%group)
!$OMP END MASTER
!$OMP BARRIER
      ALLOCATE(local_cost_matrix(SIZE(cost_matrix,1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      local_cost_matrix = cost_matrix
!$OMP MASTER
      ALLOCATE(shm_distribution_vector(ncpu*nbins*n_threads), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL optimize_distribution(ncpu*nbins*n_threads, ncpu*n_threads, local_cost_matrix, &
                                 shm_distribution_vector, x_data%load_balance_parameter%do_randomize, &
                                 error)

      ALLOCATE(full_dist(ncpu*n_threads,max_bin_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      full_dist(:,:)%istart = 0_int_8
      full_dist(:,:)%number_of_atom_quartets = 0_int_8
      full_dist(:,:)%cost = 0_int_8
      full_dist(:,:)%time_first_scf = 0.0_dp
      full_dist(:,:)%time_other_scf = 0.0_dp
      full_dist(:,:)%time_forces = 0.0_dp
!$OMP END MASTER

!$OMP BARRIER
      mepos = para_env%mepos + 1
      full_dist((mepos-1)*n_threads+i_thread+1,1:my_bin_size) = binned_dist(1:my_bin_size)
!$OMP BARRIER
!$OMP MASTER
      ALLOCATE(sendbuffer(3*max_bin_size*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(recbuffer(3*max_bin_size*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      mepos = para_env%mepos
      DO j=1,n_threads
        DO i=1,max_bin_size
          sendbuffer((j-1)*3*max_bin_size+(i-1)*3+1) = full_dist(mepos*n_threads+j,i)%istart
          sendbuffer((j-1)*3*max_bin_size+(i-1)*3+2) = full_dist(mepos*n_threads+j,i)%number_of_atom_quartets
          sendbuffer((j-1)*3*max_bin_size+(i-1)*3+3) = full_dist(mepos*n_threads+j,i)%cost
        END DO
      END DO
      dest   = MODULO(mepos+1,ncpu)
      source = MODULO(mepos-1,ncpu)
      ! sync before/after ring of isendrecv
      CALL mp_sync(para_env%group)
      DO icpu = 0, ncpu-1
        IF (icpu .NE. ncpu-1) THEN
          CALL mp_isendrecv(sendbuffer,dest,recbuffer,source, &
                            para_env%group,req(1),req(2),13)
        ENDIF
        data_from = MODULO(mepos - icpu,ncpu)
        DO j=1,n_threads
          DO i=1,max_bin_size
            full_dist(data_from*n_threads+j,i)%istart = sendbuffer((j-1)*3*max_bin_size+(i-1)*3+1)
            full_dist(data_from*n_threads+j,i)%number_of_atom_quartets = sendbuffer((j-1)*3*max_bin_size+(i-1)*3+2)
            full_dist(data_from*n_threads+j,i)%cost = sendbuffer((j-1)*3*max_bin_size+(i-1)*3+3)
          END DO
        END DO

        IF (icpu .NE. ncpu -1) THEN
          CALL mp_waitall(req)
        ENDIF
        swapbuffer=>sendbuffer
        sendbuffer=>recbuffer
        recbuffer=>swapbuffer
      ENDDO
      ! sync before/after ring of isendrecv
      DEALLOCATE(recbuffer, sendbuffer, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL mp_sync(para_env%group)
!$OMP END MASTER
!$OMP BARRIER
      !! reorder the distribution according to the distribution vector
      ALLOCATE(tmp_pos(ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_pos = 1
      ALLOCATE(tmp_dist(nbins*ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      tmp_dist(:)%istart = 0_int_8
      tmp_dist(:)%number_of_atom_quartets = 0_int_8
      tmp_dist(:)%cost = 0_int_8
      tmp_dist(:)%time_first_scf = 0.0_dp
      tmp_dist(:)%time_other_scf = 0.0_dp
      tmp_dist(:)%time_forces = 0.0_dp

      mepos = my_process_id + 1
      DO icpu= 1,n_processes
        DO i = 1,bin_histogram(icpu,1)
          IF(shm_distribution_vector(bin_histogram(icpu,2)- bin_histogram(icpu,1) + i) == mepos) THEN
               tmp_dist(tmp_pos(mepos)) = full_dist(icpu,i)
               tmp_pos(mepos) = tmp_pos(mepos) + 1
          END IF
        END DO
      END DO

      !! Assign the load to each process
      NULLIFY(ptr_to_tmp_dist)
      mepos = my_process_id + 1
      ptr_to_tmp_dist => tmp_dist(1:tmp_pos(mepos)-1)
      SELECT CASE (eval_type)
      CASE (hfx_do_eval_energy)
        CALL hfx_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)
      CASE (hfx_do_eval_forces)
        CALL hfx_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)
      END SELECT

!$OMP BARRIER
!$OMP MASTER
      DEALLOCATE(full_dist, cost_matrix, shm_distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(bins_per_rank,bin_histogram, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER
      DEALLOCATE(tmp_dist,  tmp_pos, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(binned_dist, local_cost_matrix, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
!$OMP BARRIER
!$OMP MASTER
    CALL timestop(handle)
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE hfx_update_load_balance

! *****************************************************************************
!> \brief estimates the cost of a set quartet with info available at load balance time
!>        i.e. without much info on the primitives primitives
!> \param nsa ...
!> \param nsb ...
!> \param nsc ...
!> \param nsd ...
!> \param npgfa ...
!> \param npgfb ...
!> \param npgfc ...
!> \param npgfd ...
!> \param ratio ...
!> \param p1 ...
!> \param p2 ...
!> \param p3 ...
!> \retval res ...
!> \par History
!>      08.2009 created Joost VandeVondele
!> \author Joost VandeVondele
! *****************************************************************************
  FUNCTION cost_model(nsa,nsb,nsc,nsd,npgfa,npgfb,npgfc,npgfd,ratio,p1,p2,p3) RESULT(res)
    IMPLICIT NONE
    REAL(KIND=dp) :: estimate1,estimate2,estimate,ratio,switch,mu,sigma
    INTEGER(KIND=int_8) :: res
    REAL(KIND=dp), INTENT(IN) :: p1(12),p2(12),p3(2)

    INTEGER   :: nsa,nsb,nsc,nsd,npgfa,npgfb,npgfc,npgfd

    estimate1=estimate_basic(p1)
    estimate2=estimate_basic(p2)
    mu=LOG(ABS(1.0E6_dp*p3(1))+1)
    sigma=p3(2)*0.1_dp*mu
    switch=1.0_dp/(1.0_dp+EXP((LOG(estimate1)-mu)/sigma))
    estimate=estimate1*(1.0_dp-switch)+estimate2*switch
    res = INT(estimate*0.001_dp,KIND=int_8)+1

  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param p ...
!> \retval res ...
! *****************************************************************************
    REAL(KIND=dp) FUNCTION estimate_basic(p) RESULT(res)
      REAL(KIND=dp) :: p(12)
      REAL(KIND=dp) :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12
       p1=p(1) ; p2=p(2) ; p3=p(3) ; p4=p(4)
       p5=p(5) ; p6=p(6) ; p7=p(7) ; p8=p(8)
       p9=p(9) ; p10=p(10) ; p11=p(11) ; p12=p(12)
       res=poly2(nsa,p1,p2,p3)*poly2(nsb,p1,p2,p3)*poly2(nsc,p1,p2,p3)*poly2(nsd,p1,p2,p3) * &
           poly2(npgfa,p4,p5,p6)*poly2(npgfb,p4,p5,p6)*poly2(npgfc,p4,p5,p6)*poly2(npgfd,p4,p5,p6)*EXP(-p7*ratio+p8*ratio**2)+&
           1000.0_dp*p9+poly2(nsa,p10,p11,p12)*poly2(nsb,p10,p11,p12)*poly2(nsc,p10,p11,p12)*poly2(nsd,p10,p11,p12)
       res=1+ABS(res)
    END FUNCTION estimate_basic

! *****************************************************************************
!> \brief ...
!> \param x ...
!> \param a0 ...
!> \param a1 ...
!> \param a2 ...
!> \retval poly2 ...
! *****************************************************************************
    REAL(KIND=dp) FUNCTION poly2(x,a0,a1,a2)
       INTEGER :: x
       REAL(KIND=dp) :: a0,a1,a2
       poly2=a0+a1*x+a2*x*x
    END FUNCTION poly2

  END FUNCTION cost_model
! *****************************************************************************
!> \brief Minimizes the maximum cost per cpu by shuffling around all bins
!> \param total_number_of_bins ...
!> \param number_of_processes ...
!> \param bin_costs costs per bin
!> \param distribution_vector will contain the final distribution
!> \param do_randomize ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2009 created from a hack by Joost [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE optimize_distribution(total_number_of_bins, number_of_processes, bin_costs, &
                                   distribution_vector, do_randomize, error)
    INTEGER                                  :: total_number_of_bins, &
                                                number_of_processes
    INTEGER(int_8), DIMENSION(:), POINTER    :: bin_costs
    INTEGER, DIMENSION(:), POINTER           :: distribution_vector
    LOGICAL, INTENT(IN)                      :: do_randomize
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'optimize_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, itmp, j, nstep, stat
    INTEGER(int_8), DIMENSION(:), POINTER    :: my_cost_cpu, tmp_cost, &
                                                tmp_cpu_cost
    INTEGER, DIMENSION(:), POINTER           :: tmp_cpu_index, tmp_index
    LOGICAL                                  :: failure
    TYPE(rng_stream_type), POINTER           :: rng_stream

    failure = .FALSE.

    nstep=MAX(1,INT(number_of_processes)/2)

    ALLOCATE(tmp_cost(total_number_of_bins),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tmp_index(total_number_of_bins),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tmp_cpu_cost(number_of_processes),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tmp_cpu_index(number_of_processes),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(my_cost_cpu(number_of_processes),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    tmp_cost=bin_costs

    CALL sort(tmp_cost,total_number_of_bins,tmp_index)
    my_cost_cpu=0
    !
    ! assign the largest remaining bin to the CPU with the smallest load
    ! gives near perfect distributions for a sufficient number of bins ...
    ! doing this in chunks of nstep (where nstep ~ number_of_processes) makes this n log n and gives
    ! each cpu a similar number of tasks.
    ! it also avoids degenerate cases where thousands of zero sized tasks
    ! are assigned to the same (least loaded) cpu
    !
    IF ( do_randomize ) THEN
      NULLIFY(rng_stream)
      CALL create_rng_stream(rng_stream=rng_stream,&
                             name="uniform_rng",&
                             distribution_type=UNIFORM,error=error)
    END IF

    DO i=total_number_of_bins,1,-nstep
       tmp_cpu_cost=my_cost_cpu
       CALL sort(tmp_cpu_cost,INT(number_of_processes),tmp_cpu_index)
       IF( do_randomize ) THEN
         CALL reshuffle(MIN(i,nstep), tmp_cpu_index(1:MIN(i,nstep)),rng_stream, error)
       END IF
       DO j=1,MIN(i,nstep)
          itmp=tmp_cpu_index(j)
          distribution_vector(tmp_index(i-j+1))=itmp
          my_cost_cpu(itmp)=my_cost_cpu(itmp)+bin_costs(tmp_index(i-j+1))
       ENDDO
    ENDDO

    IF( do_randomize ) THEN
      CALL delete_rng_stream(rng_stream,error=error)
    END IF

    DEALLOCATE(tmp_cost,tmp_index,tmp_cpu_cost,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tmp_cpu_index,my_cost_cpu,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE optimize_distribution

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
!> \brief ...
!> \param natom ...
!> \param nkind ...
!> \param list_ij ...
!> \param list_kl ...
!> \param set_list_ij ...
!> \param set_list_kl ...
!> \param iatom_start ...
!> \param iatom_end ...
!> \param jatom_start ...
!> \param jatom_end ...
!> \param katom_start ...
!> \param katom_end ...
!> \param latom_start ...
!> \param latom_end ...
!> \param particle_set ...
!> \param coeffs_set ...
!> \param coeffs_kind ...
!> \param is_assoc_atomic_block_global ...
!> \param do_periodic ...
!> \param kind_of ...
!> \param basis_parameter ...
!> \param pmax_set ...
!> \param pmax_atom ...
!> \param pmax_blocks ...
!> \param cell ...
!> \param do_p_screening ...
!> \param map_atom_to_kind_atom ...
!> \param eval_type ...
!> \param log10_eps_schwarz ...
!> \param log_2 ...
!> \param coeffs_kind_max0 ...
!> \param use_virial ...
!> \param atomic_pair_list ...
!> \retval estimate_block_cost ...
! *****************************************************************************
  FUNCTION estimate_block_cost(natom, nkind, list_ij, list_kl, set_list_ij, set_list_kl, &
                               iatom_start, iatom_end, jatom_start, jatom_end, &
                               katom_start, katom_end, latom_start, latom_end, &
                               particle_set,&
                               coeffs_set, coeffs_kind, &
                               is_assoc_atomic_block_global, do_periodic,&
                               kind_of, basis_parameter, pmax_set, pmax_atom, pmax_blocks, &
                               cell, &
                               do_p_screening, map_atom_to_kind_atom, eval_type,&
                               log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, &
                               atomic_pair_list)

    INTEGER, INTENT(IN)                      :: natom, nkind
    TYPE(pair_list_type)                     :: list_ij, list_kl
    TYPE(pair_set_list_type), DIMENSION(:)   :: set_list_ij, set_list_kl
    INTEGER, INTENT(IN) :: iatom_start, iatom_end, jatom_start, jatom_end, &
      katom_start, katom_end, latom_start, latom_end
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(nkind, nkind)                :: coeffs_kind
    INTEGER, DIMENSION(:, :)                 :: is_assoc_atomic_block_global
    LOGICAL                                  :: do_periodic
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_p_kind), DIMENSION(:), POINTER  :: pmax_set
    REAL(dp), DIMENSION(:, :), POINTER       :: pmax_atom
    REAL(dp)                                 :: pmax_blocks
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(IN)                      :: do_p_screening
    INTEGER, DIMENSION(:), POINTER           :: map_atom_to_kind_atom
    INTEGER, INTENT(IN)                      :: eval_type
    REAL(dp)                                 :: log10_eps_schwarz, log_2, &
                                                coeffs_kind_max0
    LOGICAL, INTENT(IN)                      :: use_virial
    LOGICAL, DIMENSION(natom, natom)         :: atomic_pair_list
    INTEGER(int_8)                           :: estimate_block_cost

    INTEGER :: i_list_ij, i_list_kl, i_set_list_ij, i_set_list_ij_start, &
      i_set_list_ij_stop, i_set_list_kl, i_set_list_kl_start, &
      i_set_list_kl_stop, iatom, ikind, iset, jatom, jkind, jset, katom, &
      kind_kind_idx, kkind, kset, latom, lkind, lset, swap_id
    INTEGER, DIMENSION(:), POINTER           :: npgfa, npgfb, npgfc, npgfd, &
                                                nsgfa, nsgfb, nsgfc, nsgfd
    REAL(dp) :: actual_pmax_atom, cost_tmp, max_val1, max_val2, pmax_entry, &
      rab2, rcd2, screen_kind_ij, screen_kind_kl
    REAL(dp), DIMENSION(:, :), POINTER       :: ptr_p_1, ptr_p_2, ptr_p_3, &
                                                ptr_p_4

    estimate_block_cost = 0_int_8

    CALL build_pair_list(natom,list_ij, set_list_ij, iatom_start, iatom_end, jatom_start, jatom_end,&
                         kind_of, basis_parameter, particle_set, &
                         do_periodic, coeffs_set, coeffs_kind, coeffs_kind_max0, &
                         log10_eps_schwarz, cell, pmax_blocks, atomic_pair_list)


    CALL build_pair_list(natom,list_kl, set_list_kl, katom_start, katom_end, latom_start, latom_end,&
                         kind_of, basis_parameter, particle_set, &
                         do_periodic, coeffs_set, coeffs_kind, coeffs_kind_max0, &
                         log10_eps_schwarz,cell, pmax_blocks, atomic_pair_list)


    DO i_list_ij=1,list_ij%n_element
      iatom=list_ij%elements(i_list_ij)%pair(1)
      jatom=list_ij%elements(i_list_ij)%pair(2)
      i_set_list_ij_start=list_ij%elements(i_list_ij)%set_bounds(1)
      i_set_list_ij_stop=list_ij%elements(i_list_ij)%set_bounds(2)
      ikind=list_ij%elements(i_list_ij)%kind_pair(1)
      jkind=list_ij%elements(i_list_ij)%kind_pair(2)
      rab2=list_ij%elements(i_list_ij)%dist2

      nsgfa => basis_parameter(ikind)%nsgf
      nsgfb => basis_parameter(jkind)%nsgf
      npgfa => basis_parameter(ikind)%npgf
      npgfb => basis_parameter(jkind)%npgf

      DO i_list_kl=1,list_kl%n_element

        katom=list_kl%elements(i_list_kl)%pair(1)
        latom=list_kl%elements(i_list_kl)%pair(2)

        IF( .NOT. (katom+latom<=iatom+jatom))  CYCLE
        IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE

        IF( eval_type == hfx_do_eval_forces) THEN
          IF( .NOT. use_virial ) THEN
            IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE
          END IF
        END IF

        i_set_list_kl_start=list_kl%elements(i_list_kl)%set_bounds(1)
        i_set_list_kl_stop=list_kl%elements(i_list_kl)%set_bounds(2)
        kkind=list_kl%elements(i_list_kl)%kind_pair(1)
        lkind=list_kl%elements(i_list_kl)%kind_pair(2)
        rcd2=list_kl%elements(i_list_kl)%dist2

        nsgfc => basis_parameter(kkind)%nsgf
        nsgfd => basis_parameter(lkind)%nsgf
        npgfc => basis_parameter(kkind)%npgf
        npgfd => basis_parameter(lkind)%npgf

        IF( do_p_screening) THEN
          actual_pmax_atom = MAX(pmax_atom(katom,iatom), &
                                 pmax_atom(latom,jatom), &
                                 pmax_atom(latom,iatom), &
                                 pmax_atom(katom,jatom) )
        ELSE
          actual_pmax_atom = 0.0_dp
        END IF

        screen_kind_ij = coeffs_kind(jkind,ikind)%x(1)*rab2+&
                         coeffs_kind(jkind,ikind)%x(2)
        screen_kind_kl = coeffs_kind(lkind,kkind)%x(1)*rcd2+&
                         coeffs_kind(lkind,kkind)%x(2)
        IF( screen_kind_ij + screen_kind_kl + actual_pmax_atom < log10_eps_schwarz) CYCLE

        IF(.NOT. (is_assoc_atomic_block_global(latom,iatom)  >=1 .AND. &
                  is_assoc_atomic_block_global(katom,iatom)  >=1 .AND. &
                  is_assoc_atomic_block_global(katom,jatom)  >=1 .AND. &
                  is_assoc_atomic_block_global(latom,jatom)  >=1 ) ) CYCLE

        IF( do_p_screening ) THEN
          SELECT CASE( eval_type)
          CASE (hfx_do_eval_energy)
            swap_id = 0
            kind_kind_idx = INT(get_1D_idx(kkind,ikind,INT(nkind,int_8)))
            IF( ikind >= kkind) THEN
              ptr_p_1 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(katom), &
                                  map_atom_to_kind_atom(iatom))
            ELSE
              ptr_p_1 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(iatom), &
                                  map_atom_to_kind_atom(katom))
              swap_id = swap_id + 1
            END IF
            kind_kind_idx = INT(get_1D_idx(lkind,jkind,INT(nkind,int_8)))
            IF( jkind>=lkind ) THEN
              ptr_p_2 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(latom), &
                                  map_atom_to_kind_atom(jatom))
            ELSE
              ptr_p_2 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(jatom), &
                                  map_atom_to_kind_atom(latom))
              swap_id = swap_id + 2
            END IF
            kind_kind_idx = INT(get_1D_idx(lkind,ikind,INT(nkind,int_8)))
            IF( ikind>=lkind ) THEN
              ptr_p_3 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(latom), &
                                  map_atom_to_kind_atom(iatom))
            ELSE
              ptr_p_3 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(iatom), &
                                  map_atom_to_kind_atom(latom))
              swap_id = swap_id + 4
            END IF
            kind_kind_idx = INT(get_1D_idx(kkind,jkind,INT(nkind,int_8)))
            IF( jkind>=kkind) THEN
              ptr_p_4 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(katom), &
                                  map_atom_to_kind_atom(jatom))
            ELSE
              ptr_p_4 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(jatom), &
                                  map_atom_to_kind_atom(katom))
              swap_id = swap_id + 8
            END IF
          CASE (hfx_do_eval_forces)
            swap_id = 16
            kind_kind_idx = INT(get_1D_idx(kkind,ikind,INT(nkind,int_8)))
            IF( ikind >= kkind) THEN
              ptr_p_1 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(katom), &
                                  map_atom_to_kind_atom(iatom))
            ELSE
              ptr_p_1 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(iatom), &
                                  map_atom_to_kind_atom(katom))
              swap_id = swap_id + 1
            END IF
            kind_kind_idx = INT(get_1D_idx(lkind,jkind,INT(nkind,int_8)))
            IF( jkind>=lkind ) THEN
              ptr_p_2 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(latom), &
                                  map_atom_to_kind_atom(jatom))
            ELSE
              ptr_p_2 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(jatom), &
                                  map_atom_to_kind_atom(latom))
              swap_id = swap_id + 2
            END IF
            kind_kind_idx = INT(get_1D_idx(lkind,ikind,INT(nkind,int_8)))
            IF( ikind>=lkind ) THEN
              ptr_p_3 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(latom), &
                                  map_atom_to_kind_atom(iatom))
            ELSE
              ptr_p_3 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(iatom), &
                                  map_atom_to_kind_atom(latom))
              swap_id = swap_id + 4
            END IF
            kind_kind_idx = INT(get_1D_idx(kkind,jkind,INT(nkind,int_8)))
            IF( jkind>=kkind) THEN
              ptr_p_4 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(katom), &
                                  map_atom_to_kind_atom(jatom))
            ELSE
              ptr_p_4 => pmax_set(kind_kind_idx)%p_kind(:,:, &
                                  map_atom_to_kind_atom(jatom), &
                                  map_atom_to_kind_atom(katom))
              swap_id = swap_id + 8
            END IF
          END SELECT
        END IF

        DO i_set_list_ij=i_set_list_ij_start, i_set_list_ij_stop
          iset=set_list_ij(i_set_list_ij)%pair(1)
          jset=set_list_ij(i_set_list_ij)%pair(2)

          max_val1 = coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + &
                     coeffs_set(jset,iset,jkind,ikind)%x(2)

          IF( max_val1 + screen_kind_kl + actual_pmax_atom < log10_eps_schwarz) CYCLE
          DO i_set_list_kl=i_set_list_kl_start, i_set_list_kl_stop
            kset=set_list_kl(i_set_list_kl)%pair(1)
            lset=set_list_kl(i_set_list_kl)%pair(2)

            max_val2 = max_val1 + (coeffs_set(lset,kset,lkind,kkind)%x(1)*rcd2 + &
                                   coeffs_set(lset,kset,lkind,kkind)%x(2) )

            IF(max_val2 + actual_pmax_atom <log10_eps_schwarz ) CYCLE
            IF( do_p_screening ) THEN
              CALL get_pmax_val(ptr_p_1, ptr_p_2, ptr_p_3, ptr_p_4, &
                                iset, jset, kset, lset, &
                                pmax_entry, swap_id)
              IF (eval_type == hfx_do_eval_forces) THEN
                pmax_entry = log_2 + pmax_entry
              END IF
            ELSE
              pmax_entry = 0.0_dp
            END IF
            max_val2= max_val2 + pmax_entry
            IF( max_val2< log10_eps_schwarz) CYCLE
            SELECT CASE (eval_type)
            CASE (hfx_do_eval_energy)
              cost_tmp = cost_model(nsgfa(iset),nsgfb(jset),nsgfc(kset),nsgfd(lset), &
                                                          npgfa(iset),npgfb(jset),npgfc(kset),npgfd(lset), &
                                                          max_val2/log10_eps_schwarz,&
                                                          p1_energy,p2_energy,p3_energy)
              estimate_block_cost = estimate_block_cost + INT(cost_tmp,KIND=int_8)
            CASE (hfx_do_eval_forces)
              cost_tmp = cost_model(nsgfa(iset),nsgfb(jset),nsgfc(kset),nsgfd(lset), &
                                                          npgfa(iset),npgfb(jset),npgfc(kset),npgfd(lset), &
                                                          max_val2/log10_eps_schwarz,&
                                                          p1_forces,p2_forces,p3_forces)
              estimate_block_cost = estimate_block_cost + INT(cost_tmp,KIND=int_8)
            END SELECT
          ENDDO ! i_set_list_kl
        ENDDO ! i_set_list_ij
      ENDDO ! i_list_kl
    ENDDO ! i_list_ij

  END FUNCTION estimate_block_cost

! *****************************************************************************
!> \brief ...
!> \param nkind ...
!> \param para_env ...
!> \param natom ...
!> \param block_size ...
!> \param nblock ...
!> \param blocks ...
!> \param list_ij ...
!> \param list_kl ...
!> \param set_list_ij ...
!> \param set_list_kl ...
!> \param particle_set ...
!> \param coeffs_set ...
!> \param coeffs_kind ...
!> \param is_assoc_atomic_block_global ...
!> \param do_periodic ...
!> \param kind_of ...
!> \param basis_parameter ...
!> \param pmax_set ...
!> \param pmax_atom ...
!> \param pmax_blocks ...
!> \param cell ...
!> \param do_p_screening ...
!> \param map_atom_to_kind_atom ...
!> \param eval_type ...
!> \param log10_eps_schwarz ...
!> \param log_2 ...
!> \param coeffs_kind_max0 ...
!> \param use_virial ...
!> \param atomic_pair_list ...
! *****************************************************************************
  SUBROUTINE init_blocks(nkind, para_env, natom, block_size, nblock, blocks, &
                         list_ij, list_kl, set_list_ij, set_list_kl, &
                         particle_set,&
                         coeffs_set, coeffs_kind, &
                         is_assoc_atomic_block_global, do_periodic,&
                         kind_of, basis_parameter, pmax_set, pmax_atom, &
                         pmax_blocks, cell, &
                         do_p_screening, map_atom_to_kind_atom, eval_type,&
                         log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial,&
                         atomic_pair_list)

    INTEGER, INTENT(IN)                      :: nkind
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: natom, block_size, nblock
    TYPE(hfx_block_range_type), &
      DIMENSION(1:nblock)                    :: blocks
    TYPE(pair_list_type)                     :: list_ij, list_kl
    TYPE(pair_set_list_type), DIMENSION(:)   :: set_list_ij, set_list_kl
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: coeffs_kind
    INTEGER, DIMENSION(:, :)                 :: is_assoc_atomic_block_global
    LOGICAL                                  :: do_periodic
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_p_kind), DIMENSION(:), POINTER  :: pmax_set
    REAL(dp), DIMENSION(:, :), POINTER       :: pmax_atom
    REAL(dp)                                 :: pmax_blocks
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(IN)                      :: do_p_screening
    INTEGER, DIMENSION(:), POINTER           :: map_atom_to_kind_atom
    INTEGER, INTENT(IN)                      :: eval_type
    REAL(dp)                                 :: log10_eps_schwarz, log_2, &
                                                coeffs_kind_max0
    LOGICAL, INTENT(IN)                      :: use_virial
    LOGICAL, DIMENSION(natom, natom)         :: atomic_pair_list

    INTEGER                                  :: atom_block, i, iatom_block, &
                                                iatom_end, iatom_start, &
                                                my_cpu_rank, ncpus

    DO atom_block = 0,nblock-1
      iatom_block = MODULO(atom_block,nblock)+1
      iatom_start = (iatom_block-1)*block_size + 1
      iatom_end = MIN(iatom_block*block_size,natom)
      blocks(atom_block+1)%istart = iatom_start
      blocks(atom_block+1)%iend   = iatom_end
      blocks(atom_block+1)%cost   = 0_int_8
    END DO

    ncpus = para_env%num_pe
    my_cpu_rank = para_env%mepos
    DO i=1,nblock
      IF( MODULO(i,ncpus) /= my_cpu_rank)  THEN
        blocks(i)%istart = 0
        blocks(i)%iend = 0
        CYCLE
      END IF
      iatom_start = blocks(i)%istart
      iatom_end = blocks(i)%iend
      blocks(i)%cost = estimate_block_cost(natom, nkind, list_ij, list_kl, set_list_ij, set_list_kl, &
                                           iatom_start, iatom_end, iatom_start, iatom_end, &
                                           iatom_start, iatom_end, iatom_start, iatom_end, &
                                           particle_set,&
                                           coeffs_set, coeffs_kind, &
                                           is_assoc_atomic_block_global, do_periodic,&
                                           kind_of, basis_parameter, pmax_set, pmax_atom, pmax_blocks, &
                                           cell, &
                                           do_p_screening, map_atom_to_kind_atom, eval_type,&
                                           log10_eps_schwarz, log_2, coeffs_kind_max0, use_virial, atomic_pair_list)

    END DO
  END SUBROUTINE init_blocks


! *****************************************************************************
!> \brief ...
!> \param para_env ...
!> \param x_data ...
!> \param iw ...
!> \param n_threads ...
!> \param i_thread ...
!> \param eval_type ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE  collect_load_balance_info(para_env, x_data, iw, n_threads, i_thread, &
                                        eval_type, error)

    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(hfx_type), POINTER                  :: x_data
    INTEGER, INTENT(IN)                      :: iw, n_threads, i_thread, &
                                                eval_type
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'collect_load_balance_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, my_rank, nbins, &
                                                nranks, stat, total_bins
    INTEGER(int_8)                           :: avg_bin, avg_rank, max_bin, &
                                                max_rank, min_bin, min_rank, &
                                                sum_bin, sum_rank
    INTEGER(int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: buffer, buffer_in, &
                                                buffer_out, summary
    INTEGER(int_8), ALLOCATABLE, &
      DIMENSION(:), SAVE                     :: shm_cost_vector
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bins_per_rank, rdispl, &
                                                sort_idx
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: shm_bins_per_rank, shm_displ
    LOGICAL                                  :: failure

    failure = .FALSE.


    SELECT CASE (eval_type)
    CASE (hfx_do_eval_energy)
      nbins = SIZE(x_data%distribution_energy)
    CASE (hfx_do_eval_forces)
      nbins = SIZE(x_data%distribution_forces)
    END SELECT

!$OMP MASTER
    ALLOCATE(shm_bins_per_rank(n_threads),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(shm_displ(n_threads+1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER

    shm_bins_per_rank(i_thread+1) = nbins
!$OMP BARRIER
    nbins = 0
    DO i = 1,n_threads
      nbins = nbins + shm_bins_per_rank(i)
    END DO
    my_rank = para_env%mepos
    nranks = para_env%num_pe

!$OMP BARRIER
!$OMP MASTER
    ALLOCATE(bins_per_rank(nranks),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    bins_per_rank = 0

    bins_per_rank(my_rank+1) = nbins

    CALL mp_sum(bins_per_rank, para_env%group)

    total_bins = 0
    DO i = 1,nranks
      total_bins = total_bins + bins_per_rank(i)
    END DO

    ALLOCATE(shm_cost_vector(2*total_bins),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    shm_cost_vector = -1_int_8
    shm_displ(1) = 1
    DO i = 2,n_threads
      shm_displ(i) = shm_displ(i-1) + shm_bins_per_rank(i-1)
    END DO
    shm_displ(n_threads+1) = nbins+1
!$OMP END MASTER
!$OMP BARRIER
    j = 0
    SELECT CASE (eval_type)
    CASE (hfx_do_eval_energy)
      DO i = shm_displ(i_thread+1), shm_displ(i_thread+2)-1
        j = j + 1
        shm_cost_vector(2*(i-1)+1) = x_data%distribution_energy(j)%cost
        shm_cost_vector(2*i) = INT(x_data%distribution_energy(j)%time_first_scf*10000.0_dp,KIND=int_8)
      END DO
    CASE (hfx_do_eval_forces)
       DO i = shm_displ(i_thread+1), shm_displ(i_thread+2)-1
        j = j + 1
        shm_cost_vector(2*(i-1)+1) = x_data%distribution_forces(j)%cost
        shm_cost_vector(2*i) = INT(x_data%distribution_forces(j)%time_forces*10000.0_dp,KIND=int_8)
      END DO
    END SELECT
!$OMP BARRIER
!$OMP MASTER
    ! ** calculate offsets
    ALLOCATE(rdispl(nranks), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    bins_per_rank = bins_per_rank * 2
    rdispl(1) = 0
    DO i = 2,nranks
      rdispl(i) = rdispl(i-1) + bins_per_rank(i-1)
    END DO

    ALLOCATE(buffer_in(2*nbins), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(buffer_out(2*total_bins), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO i=1,nbins
      buffer_in(2*(i-1)+1) = shm_cost_vector(2*(i-1)+1)
      buffer_in(2*i) = shm_cost_vector(2*i)
    END DO

    CALL mp_gatherv(buffer_in,buffer_out,bins_per_rank,rdispl,para_env%source,para_env%group)

    IF( iw>0 ) THEN

      ALLOCATE(summary(2*nranks),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      summary = 0_int_8

      WRITE(iw,'( /, 1X, 79("-") )' )
      WRITE(iw,'( " -", 77X, "-" )' )
      SELECT CASE (eval_type)
      CASE (hfx_do_eval_energy)
        WRITE(iw,'( " -", 20X, A, 19X, "-" )' ) ' HFX LOAD BALANCE INFORMATION - ENERGY '
      CASE (hfx_do_eval_forces)
        WRITE(iw,'( " -", 20X, A, 19X, "-" )' ) ' HFX LOAD BALANCE INFORMATION - FORCES '
      END SELECT
      WRITE(iw,'( " -", 77X, "-" )' )
      WRITE(iw,'( 1X, 79("-") )' )

      WRITE(iw,FMT="(T3,A,T15,A,T35,A,T55,A)") "MPI RANK", "BIN #", "EST cost", "Processing time [s]"
      WRITE(iw, '( 1X, 79("-"), / )' )
      k = 0
      DO i=1,nranks
        DO j = 1, bins_per_rank(i)/2
          k = k + 1
          WRITE(iw,FMT="(T6,I5,T15,I5,T27,I16,T55,F19.8)") i-1, j, buffer_out(2*(k-1)+1), REAL(buffer_out(2*k),dp)/10000.0_dp
          summary(2*(i-1)+1) = summary(2*(i-1)+1) + buffer_out(2*(k-1)+1)
          summary(2*i) = summary(2*i) + buffer_out(2*k)
        END DO
      END DO

      !** Summary
      max_bin = 0_int_8
      min_bin = HUGE(min_bin)
      sum_bin = 0_int_8
      DO i=1,total_bins
        sum_bin = sum_bin + buffer_out(2*i)
        max_bin = MAX(max_bin,buffer_out(2*i))
        min_bin = MIN(min_bin,buffer_out(2*i))
      END DO
      avg_bin = sum_bin/total_bins

      max_rank = 0_int_8
      min_rank = HUGE(min_rank)
      sum_rank = 0_int_8
      DO i=1,nranks
        sum_rank = sum_rank + summary(2*i)
        max_rank = MAX(max_rank,summary(2*i))
        min_rank = MIN(min_rank,summary(2*i))
      END DO
      avg_rank = sum_rank/nranks

      WRITE(iw,FMT='(/,T3,A,/)') "SUMMARY:"
      WRITE(iw,FMT="(T3,A,T35,F19.8)") "Max bin", REAL(max_bin,dp)/10000.0_dp
      WRITE(iw,FMT="(T3,A,T35,F19.8)") "Min bin", REAL(min_bin,dp)/10000.0_dp
      WRITE(iw,FMT="(T3,A,T35,F19.8)") "Sum bin", REAL(sum_bin,dp)/10000.0_dp
      WRITE(iw,FMT="(T3,A,T35,F19.8,/)") "Avg bin", REAL(avg_bin,dp)/10000.0_dp
      WRITE(iw,FMT="(T3,A,T35,F19.8)") "Max rank", REAL(max_rank,dp)/10000.0_dp
      WRITE(iw,FMT="(T3,A,T35,F19.8)") "Min rank", REAL(min_rank,dp)/10000.0_dp
      WRITE(iw,FMT="(T3,A,T35,F19.8)") "Sum rank", REAL(sum_rank,dp)/10000.0_dp
      WRITE(iw,FMT="(T3,A,T35,F19.8,/)") "Avg rank", REAL(avg_rank,dp)/10000.0_dp


      ALLOCATE(buffer(nranks),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(sort_idx(nranks),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO i=1,nranks
        buffer(i) = summary(2*i)
      END DO

      CALL sort(buffer,nranks,sort_idx)

      WRITE(iw,FMT="(T3,A,T35,A,T55,A,/)") "MPI RANK", "EST cost", "Processing time [s]"
      DO i=nranks,1,-1
        WRITE(iw,FMT="(T6,I5,T27,I16,T55,F19.8)") sort_idx(i)-1, summary(2*(sort_idx(i)-1)+1),REAL(buffer(i),dp)/10000.0_dp
      END DO

      DEALLOCATE(summary, buffer, sort_idx, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    END IF

    DEALLOCATE(buffer_in, buffer_out, rdispl, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL mp_sync(para_env%group)

    DEALLOCATE(shm_bins_per_rank, shm_displ, shm_cost_vector, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE collect_load_balance_info

! *****************************************************************************
!> \brief ...
!> \param size ...
!> \param array ...
!> \param rng_stream ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE reshuffle(size, array, rng_stream, error)
    INTEGER                                  :: size
    INTEGER, DIMENSION(size)                 :: array
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, idx1, idx2, tmp
    REAL(dp)                                 :: x

    DO i = 1,size*10
      x = next_random_number(rng_stream,error=error)
      idx1 = INT(x*(size+1-1))+1
      x = next_random_number(rng_stream,error=error)
      idx2 = INT(x*(size+1-1))+1


      tmp = array(idx1)
      array(idx1) = array(idx2)
      array(idx2) = tmp
    END DO

  END SUBROUTINE

#include "hfx_get_pmax_val.f90"

END MODULE hfx_load_balance_methods
