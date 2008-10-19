!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines for optimizing load balance between processes in HFX calculations 
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE hfx_load_balance_methods 
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE erf_fn,                          ONLY: erfc
  USE f77_blas
  USE hfx_distribution_mc,             ONLY: optimize_distribution
  USE hfx_helpers,                     ONLY: FNV64
  USE hfx_libint_interface,            ONLY: pbc
  USE hfx_types,                       ONLY: hfx_basis_type,&
                                             hfx_box_scheme,&
                                             hfx_distribution,&
                                             hfx_load_balance_type,&
                                             hfx_potential_type,&
                                             hfx_set_distr_energy,&
                                             hfx_set_distr_forces,&
                                             hfx_type
  USE input_constants,                 ONLY: do_hfx_potential_coulomb,&
                                             do_hfx_potential_gaussian,&
                                             do_hfx_potential_id,&
                                             do_hfx_potential_mix_cl,&
                                             do_hfx_potential_mix_lg,&
                                             do_hfx_potential_short
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE mathconstants
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_isendrecv,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sum_i8,&
                                             mp_waitall
  USE orbital_pointers
  USE particle_types,                  ONLY: particle_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC load_balance_energy, load_balance_forces, update_load_balance_energy, update_load_balance_forces 

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_load_balance_methods'

!***
  
  CONTAINS

! *****************************************************************************
!> \brief Distributes the computation of eri's to all available processes.
!> \param x_data Object that stores the indices array
!> \param eps_schwarz screening parameter
!> \param particle_set , atomic_kind_set, para_env
!> \param potential_parameter short/longrange screening information 
!> \param screen_set , screen_atom, screen_max: screening matrices
!> \param is_assoc_atomic_kind KS-matrix sparsity
!> \param ks_matrix used for parallel computation of total number of integrals
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      The optimization is done via a binning procedure followed by simple 
!>      Monte Carlo procedure:
!>      In a first step the total amount of integrals in the system is calculated,
!>      taking into account the sparsity of the KS-matrix , the screening based
!>      on near/farfield approximations and if desired the screening on an initial
!>      density matrix. In a second step, these integrals are distributed to bins
!>      of a given size. Finally, a Monte Carlo procedure optimizes the assignment
!>      of the different loads to each process
!>      At the end each process owns an unique array of inidices-ranges 
!>      that are used to decide whether a process has to calculate a certain
!>      bunch of integrals or not
!> \par History
!>      06.2007 created [Manuel Guidon]
!>      08.2007 new parallel scheme [Manuel Guidon]
!>      09.2007 new 'modulo' parellel scheme and Monte Carlo step [Manuel Guidon]
!>      11.2007 parallelize load balance on box_idx1 [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE load_balance_energy(x_data, eps_schwarz, eps_box, particle_set,max_set,para_env,&
                                 potential_parameter, screen_set, screen_atom, screen_max,&
                                 is_assoc_atomic_block_global, do_periodic,&
                                 load_balance_parameter, box_size, box_size_inv, kind_of, basis_parameter, pmax_set, &
                                 i_thread, n_threads, shm_far_field_boxes, error)
    TYPE(hfx_type), POINTER                  :: x_data
    REAL(dp), INTENT(IN)                     :: eps_schwarz, eps_box
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: max_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(hfx_potential_type), INTENT(IN)     :: potential_parameter
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: screen_set
    REAL(dp), DIMENSION(:, :), POINTER       :: screen_atom
    REAL(dp), INTENT(IN)                     :: screen_max
    INTEGER, DIMENSION(:, :)                 :: is_assoc_atomic_block_global
    LOGICAL                                  :: do_periodic
    TYPE(hfx_load_balance_type)              :: load_balance_parameter
    REAL(dp), INTENT(IN)                     :: box_size(3), box_size_inv(3)
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: pmax_set
    INTEGER, INTENT(IN)                      :: i_thread, n_threads
    TYPE(hfx_box_scheme), DIMENSION(:, :), &
      POINTER                                :: shm_far_field_boxes
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'load_balance_energy', &
      routineP = moduleN//':'//routineN

    INTEGER :: data_from, dest, i, iatom, ibin, icpu, ikind, ipos, iset, j, &
      jatom, jkind, jset, katom, kkind, kset, latom, lkind, lset, &
      my_process_id, natom, nbins, nseta, nsetb, nsetc, nsetd, &
      pgf_radius_box, radius1, radius2, req(2), source, stat
    INTEGER(int_8) :: box_idx1, box_idx2, box_idx3, box_idx4, &
      current_counter, distribution_counter_end, distribution_counter_start, &
      int_per_cpu, mepos, n_processes, ncpu, optimal_cost, restore_counter, &
      tmp_counter1, total_integrals
    INTEGER(int_8), ALLOCATABLE, &
      DIMENSION(:), SAVE                     :: optimal_costs
    INTEGER(int_8), DIMENSION(:), POINTER    :: local_cost_matrix, recbuffer, &
                                                sendbuffer, swapbuffer
    INTEGER(int_8), DIMENSION(:), POINTER, &
      SAVE                                   :: cost_matrix
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_pos
    INTEGER, DIMENSION(:), POINTER           :: distribution_vector, nsgfa, &
                                                nsgfb, nsgfc, nsgfd
    INTEGER, DIMENSION(:), POINTER, SAVE     :: shm_distribution_vector
    INTEGER, SAVE                            :: rank_id_of_best_dist, &
                                                thread_id_of_best_dist
    LOGICAL                                  :: failure
    REAL(dp) :: center1(3), center2(3), center2_folded(3), Kbox, max_val1, &
      max_val2, omega2Inv, omega_box, pbc_center(3), pmax_entry, R1, R2, &
      Rbox, temp(3)
    TYPE(hfx_distribution), DIMENSION(:), &
      POINTER                                :: ptr_to_tmp_dist
    TYPE(hfx_distribution), &
      DIMENSION(:, :), POINTER               :: binned_dist, tmp_dist
    TYPE(hfx_distribution), &
      DIMENSION(:, :), POINTER, SAVE         :: full_dist

    IF(potential_parameter%omega /= 0.0_dp) THEN
      omega2Inv = 1.0_dp/(potential_parameter%omega**2)
    ELSE
      omega2Inv = 0.0_dp
    END IF
    natom=SIZE(particle_set)
    failure=.FALSE.

    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads
    !! If there is only 1 cpu skip the binning
    IF(n_processes == 1) THEN
      ALLOCATE(tmp_dist(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1,1)%number_of_integrals = HUGE(tmp_dist(1,1)%number_of_integrals)
      tmp_dist(1,1)%istart = 1
      tmp_dist(1,1)%jstart = 1
      tmp_dist(1,1)%kstart = 1
      tmp_dist(1,1)%lstart = 1
      tmp_dist(1,1)%original_cpu_id = 1
      ptr_to_tmp_dist => tmp_dist(1,:)
      CALL hfx_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)
      DEALLOCATE(tmp_dist,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      !! Calculate total numbers of integrals that have to be calculated (wrt screening and symmetry)
      total_integrals = 0
      my_process_id = para_env%mepos*n_threads + i_thread
     
      DO iatom=1,natom
        ikind = kind_of(iatom)
        nseta =  basis_parameter(ikind)%nset
        nsgfa => basis_parameter(ikind)%nsgf
        DO jatom=iatom,natom
          box_idx1 = (iatom-1)*natom+jatom
          IF(MODULO(FNV64(box_idx1),n_processes) /= my_process_id) CYCLE 
          IF (screen_atom(jatom,iatom) * screen_max < eps_schwarz) CYCLE
          jkind = kind_of(jatom)
          nsetb =  basis_parameter(jkind)%nset
          nsgfb => basis_parameter(jkind)%nsgf
          DO katom=1,natom
            kkind = kind_of(katom)
            nsetc =  basis_parameter(kkind)%nset
            nsgfc => basis_parameter(kkind)%nsgf
            DO latom=katom,natom
              IF (screen_atom(jatom,iatom) * screen_atom(latom,katom) < eps_schwarz) CYCLE
              IF(.NOT. (is_assoc_atomic_block_global(latom,iatom)  >=1 .AND. &
                        is_assoc_atomic_block_global(katom,iatom)  >=1 .AND. &
                        is_assoc_atomic_block_global(katom,jatom)  >=1 .AND. &
                        is_assoc_atomic_block_global(latom,jatom)  >=1 ) ) CYCLE 

              IF((katom+latom<=iatom+jatom))  THEN
                IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
               
                box_idx3 = (katom-1)*natom+latom
                lkind = kind_of(latom)
                nsetd =  basis_parameter(lkind)%nset
                nsgfd => basis_parameter(lkind)%nsgf
                DO iset = 1, nseta
                  DO jset = 1,nsetb
                    max_val1 = screen_set(iset,jset,iatom,jatom)
                    IF (max_val1 * screen_atom(latom,katom) < eps_schwarz) CYCLE
                    box_idx2 = (iset-1)*max_set+jset
                    center1 = shm_far_field_boxes(box_idx2, box_idx1)%center
                    R1 = shm_far_field_boxes(box_idx2, box_idx1)%R
                    radius1 = shm_far_field_boxes(box_idx2, box_idx1)%radius
                    DO kset = 1,nsetc
                      DO lset = 1,nsetd
                        max_val2 = max_val1 * screen_set(lset,kset,latom,katom)
                        pmax_entry = MAX(pmax_set(kset,iset,katom,iatom), &
                                         pmax_set(lset,jset,latom,jatom), &
                                         pmax_set(lset,iset,latom,iatom), &
                                         pmax_set(kset,jset,katom,jatom))

                        max_val2 = max_val2 * pmax_entry
                        IF(max_val2<eps_schwarz) CYCLE
                        box_idx4 = (kset-1)*max_set+lset
                        IF(.NOT. do_periodic) THEN
                          center2_folded = shm_far_field_boxes(box_idx4, box_idx3)%center
                        ELSE
                          center2 = shm_far_field_boxes(box_idx4, box_idx3)%center
                          temp = center2-center1
                          CALL pbc(temp, box_size, box_size_inv ,pbc_center)
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
                              max_val2 = 1.0_dp/Rbox*Kbox
                              max_val2 = max_val2 * pmax_entry
                            CASE(do_hfx_potential_short)
                              omega_box = shm_far_field_boxes(box_idx2, box_idx1)%alphaInv + &
                                          shm_far_field_boxes(box_idx4, box_idx3)%alphaInv + &
                                          omega2Inv 
                              omega_box = 1.0_dp/omega_box
                              Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                     shm_far_field_boxes(box_idx4, box_idx3)%K
                              max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox
                              max_val2 = max_val2 * pmax_entry
                            CASE(do_hfx_potential_mix_cl)
                              Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                     shm_far_field_boxes(box_idx4, box_idx3)%K
                              max_val2 = 1.0_dp/Rbox*Kbox
                              max_val2 = max_val2 * pmax_entry
                            CASE(do_hfx_potential_gaussian,do_hfx_potential_mix_lg, do_hfx_potential_id) 
                              max_val2 = HUGE(max_val2)
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
      ENDDO
      !! Bin the total amount of integrals on each cpu
      nbins = load_balance_parameter%nbins
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
      binned_dist(:,:)%time_first_scf = 0.0_dp
      binned_dist(:,:)%time_other_scf =0.0_dp
      binned_dist(:,:)%time_forces = 0.0_dp
 
      tmp_counter1=0
      current_counter = 0
      restore_counter = 0
      ipos = 1
      mepos = 1
      distribution_counter_start = 1
      distribution_counter_end = 0
      ibin = 1
      DO iatom=1,natom
        ikind = kind_of(iatom)
        nseta =  basis_parameter(ikind)%nset
        nsgfa => basis_parameter(ikind)%nsgf

        DO jatom=iatom,natom
          box_idx1 = (iatom-1)*natom+jatom
          IF(MODULO(FNV64(box_idx1),n_processes) /= my_process_id) CYCLE
          jkind = kind_of(jatom)
          nsetb =  basis_parameter(jkind)%nset
          nsgfb => basis_parameter(jkind)%nsgf
          DO katom=1,natom
            kkind = kind_of(katom)
            nsetc =  basis_parameter(kkind)%nset
            nsgfc => basis_parameter(kkind)%nsgf
            DO latom=katom,natom
              IF((katom+latom<=iatom+jatom))  THEN
                IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                distribution_counter_end = distribution_counter_end + 1

                IF(binned_dist(1,ibin)%istart == 0 ) binned_dist(1,ibin)%istart = iatom
                IF(binned_dist(1,ibin)%jstart == 0 ) binned_dist(1,ibin)%jstart = jatom
                IF(binned_dist(1,ibin)%kstart == 0 ) binned_dist(1,ibin)%kstart = katom
                IF(binned_dist(1,ibin)%lstart == 0 ) binned_dist(1,ibin)%lstart = latom

                IF (screen_atom(jatom,iatom) * screen_max < eps_schwarz) CYCLE
                IF (screen_atom(jatom,iatom) * screen_atom(latom,katom) < eps_schwarz) CYCLE
                IF(.NOT. (is_assoc_atomic_block_global(latom,iatom)  >=1 .AND. &
                          is_assoc_atomic_block_global(katom,iatom)  >=1 .AND. &
                          is_assoc_atomic_block_global(katom,jatom)  >=1 .AND. &
                          is_assoc_atomic_block_global(latom,jatom)  >=1 ) ) CYCLE

                box_idx3 = (katom-1)*natom+latom
                lkind = kind_of(latom)
                nsetd =  basis_parameter(lkind)%nset
                nsgfd => basis_parameter(lkind)%nsgf
                DO iset = 1, nseta
                  DO jset = 1,nsetb
                    max_val1 = screen_set(jset,iset,jatom,iatom)
                    IF (max_val1 * screen_atom(latom,katom) < eps_schwarz) CYCLE
                    box_idx2 = (iset-1)*max_set+jset
                    center1 = shm_far_field_boxes(box_idx2, box_idx1)%center
                    R1 = shm_far_field_boxes(box_idx2, box_idx1)%R
                    radius1 = shm_far_field_boxes(box_idx2, box_idx1)%radius
                    DO kset = 1,nsetc
                      DO lset = 1,nsetd
                        max_val2 = max_val1 * screen_set(lset,kset,latom,katom)
                        pmax_entry = MAX(pmax_set(kset,iset,katom,iatom), &
                                         pmax_set(lset,jset,latom,jatom), &
                                         pmax_set(lset,iset,latom,iatom), &
                                         pmax_set(kset,jset,katom,jatom))
                        max_val2 = max_val2 * pmax_entry
                        IF(max_val2<eps_schwarz) CYCLE
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
                              max_val2 = 1.0_dp/Rbox*Kbox
                              max_val2 = max_val2 * pmax_entry
                            CASE(do_hfx_potential_short)
                              omega_box = shm_far_field_boxes(box_idx2, box_idx1)%alphaInv + &
                                          shm_far_field_boxes(box_idx4, box_idx3)%alphaInv + &
                                          omega2Inv 
                              omega_box = 1.0_dp/omega_box
                              Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                     shm_far_field_boxes(box_idx4, box_idx3)%K
                              max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox
                              max_val2 = max_val2 * pmax_entry
                            CASE(do_hfx_potential_mix_cl)
                              Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                     shm_far_field_boxes(box_idx4, box_idx3)%K
                              max_val2 = 1.0_dp/Rbox*Kbox
                              max_val2 = max_val2 * pmax_entry
                            CASE(do_hfx_potential_gaussian,do_hfx_potential_mix_lg, do_hfx_potential_id) 
                              max_val2 = HUGE(max_val2)
                          END SELECT
                          IF(max_val2<eps_box) CYCLE
                        END IF
                        current_counter = current_counter + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO

                IF( current_counter >= int_per_cpu) THEN
                  binned_dist(1,ibin)%number_of_integrals = distribution_counter_end - distribution_counter_start + 1
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
      ENDDO

      !! Fill the last bin if necessary
      IF( current_counter /= 0 ) THEN  
        binned_dist(1,ibin)%cost = binned_dist(1,ibin)%cost + current_counter
        binned_dist(1,ibin)%number_of_integrals = distribution_counter_end-distribution_counter_start+1
      END IF

      !! store the original cpu_id
      binned_dist(1,:)%original_cpu_id = my_process_id + 1 !para_env%mepos + 1
!$OMP BARRIER     
      !! store all local results in a big cost matrix
!$OMP MASTER
      ALLOCATE(cost_matrix(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      cost_matrix = 0
!$OMP END MASTER
!$OMP BARRIER
      DO icpu=1,ncpu
        IF(para_env%mepos+1 == icpu ) THEN
          DO i=1,nbins
            cost_matrix((icpu-1)*nbins*n_threads +i_thread*nbins + i) = binned_dist(1,i)%cost 
          END DO
        END IF
      END DO

      mepos = para_env%mepos
!$OMP BARRIER
!$OMP MASTER
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
      !! Find an optimal distribution i.e. assign each element of the cost matrix to a certain process
      ALLOCATE(distribution_vector(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP MASTER
      ALLOCATE(optimal_costs(ncpu*n_threads), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      optimal_costs = 0
!$OMP END MASTER
!$OMP BARRIER
      optimal_cost = 0
      !! Optimize load on each node independently
      ALLOCATE(local_cost_matrix(SIZE(cost_matrix,1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      local_cost_matrix = cost_matrix
      CALL optimize_distribution(para_env, ncpu*nbins*n_threads, ncpu*n_threads, local_cost_matrix, distribution_vector, &
                                      load_balance_parameter, optimal_cost, my_process_id, i_thread, error)
      
      !! Find optimal load balance wrt the cost
      optimal_costs(my_process_id+1) = optimal_cost
!$OMP BARRIER
!$OMP MASTER
      CALL mp_sum_i8(optimal_costs(1:ncpu*n_threads), para_env%group)
      icpu = 0
      optimal_cost = HUGE(optimal_cost)
      !! Get process id of the best distribution
      DO i = 1,n_processes
        IF(optimal_costs(i) < optimal_cost ) THEN
          optimal_cost = optimal_costs(i)
          thread_id_of_best_dist = i
        END IF
      END DO
      rank_id_of_best_dist = (thread_id_of_best_dist-1)/n_threads
      ALLOCATE(shm_distribution_vector(ncpu*nbins*n_threads), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER
      IF(i_thread == MOD(thread_id_of_best_dist-1,n_threads)) THEN
        shm_distribution_vector = distribution_vector
        !! Broadcast the best distribution_vector
        CALL mp_bcast(shm_distribution_vector,rank_id_of_best_dist,para_env%group) 
      END IF      
!$OMP BARRIER
!$OMP MASTER
      !! Collect local data to global array
      ALLOCATE(full_dist(ncpu*n_threads,nbins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      full_dist(:,:)%istart = 0
      full_dist(:,:)%jstart = 0
      full_dist(:,:)%kstart = 0
      full_dist(:,:)%lstart = 0
      full_dist(:,:)%number_of_integrals = 0
      full_dist(:,:)%cost = 0
      full_dist(:,:)%original_cpu_id = 0
      full_dist(:,:)%time_first_scf = 0.0_dp
      full_dist(:,:)%time_other_scf = 0.0_dp
      full_dist(:,:)%time_forces = 0.0_dp
!$OMP END MASTER
!$OMP BARRIER
      DO mepos = 1,ncpu
        IF(para_env%mepos == mepos-1) THEN
          full_dist((mepos-1)*n_threads+i_thread+1,:) = binned_dist(1,:)
        END IF
      END DO
     
!$OMP BARRIER
!$OMP MASTER
      ALLOCATE(sendbuffer(7*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(recbuffer(7*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      mepos = para_env%mepos
      DO j=1,n_threads
        DO i=1,nbins
          sendbuffer((j-1)*7*nbins+(i-1)*7+1) = full_dist(mepos*n_threads+j,i)%istart
          sendbuffer((j-1)*7*nbins+(i-1)*7+2) = full_dist(mepos*n_threads+j,i)%jstart
          sendbuffer((j-1)*7*nbins+(i-1)*7+3) = full_dist(mepos*n_threads+j,i)%kstart
          sendbuffer((j-1)*7*nbins+(i-1)*7+4) = full_dist(mepos*n_threads+j,i)%lstart
          sendbuffer((j-1)*7*nbins+(i-1)*7+5) = full_dist(mepos*n_threads+j,i)%number_of_integrals
          sendbuffer((j-1)*7*nbins+(i-1)*7+6) = full_dist(mepos*n_threads+j,i)%cost
          sendbuffer((j-1)*7*nbins+(i-1)*7+7) = full_dist(mepos*n_threads+j,i)%original_cpu_id
        END DO
      END DO

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
            full_dist(data_from*n_threads+j,i)%istart = sendbuffer((j-1)*7*nbins+(i-1)*7+1)
            full_dist(data_from*n_threads+j,i)%jstart = sendbuffer((j-1)*7*nbins+(i-1)*7+2)
            full_dist(data_from*n_threads+j,i)%kstart = sendbuffer((j-1)*7*nbins+(i-1)*7+3)
            full_dist(data_from*n_threads+j,i)%lstart = sendbuffer((j-1)*7*nbins+(i-1)*7+4)
            full_dist(data_from*n_threads+j,i)%number_of_integrals = sendbuffer((j-1)*7*nbins+(i-1)*7+5)
            full_dist(data_from*n_threads+j,i)%cost = sendbuffer((j-1)*7*nbins+(i-1)*7+6)
            full_dist(data_from*n_threads+j,i)%original_cpu_id = sendbuffer((j-1)*7*nbins+(i-1)*7+7)
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
!$OMP END MASTER
!$OMP BARRIER
      !! reorder the distribution according to the distribution vector
      ALLOCATE(tmp_pos(ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_pos = 1
      ALLOCATE(tmp_dist(ncpu*n_threads,nbins*ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      tmp_dist(:,:)%istart = 0
      tmp_dist(:,:)%jstart = 0
      tmp_dist(:,:)%kstart = 0
      tmp_dist(:,:)%lstart = 0
      tmp_dist(:,:)%number_of_integrals = 0
      tmp_dist(:,:)%cost = 0
      tmp_dist(:,:)%original_cpu_id = 0
      tmp_dist(:,:)%time_first_scf = 0.0_dp
      tmp_dist(:,:)%time_other_scf = 0.0_dp
      tmp_dist(:,:)%time_forces = 0.0_dp

      
      DO mepos = 1,n_processes
        DO icpu= 1,n_processes
          DO i=1,nbins
            IF(shm_distribution_vector((icpu-1)*nbins + i) == mepos) THEN
              tmp_dist(mepos,tmp_pos(mepos)) = full_dist(icpu,i)
              tmp_pos(mepos) = tmp_pos(mepos) + 1
            END IF
          END DO
        END DO
      END DO

      !! Assign the load to each process
      NULLIFY(ptr_to_tmp_dist) 
      DO mepos=1,n_processes
        IF(my_process_id == mepos -1) THEN
          ptr_to_tmp_dist => tmp_dist(mepos,1:tmp_pos(mepos)-1)
        END IF
      END DO
      CALL hfx_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)
!$OMP BARRIER
!$OMP MASTER
      DEALLOCATE(full_dist,optimal_costs,cost_matrix,shm_distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER
      DEALLOCATE(tmp_dist, tmp_pos, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(binned_dist, local_cost_matrix, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE load_balance_energy

! *****************************************************************************
!> \brief Distributes the computation of eri-derivatives to all available processes.
!> \param x_data Object that stores the indices array
!> \param eps_schwarz screening parameter
!> \param particle_set , atomic_kind_set, para_env
!> \param potential_parameter short/longrange screening information 
!> \param screen_set , screen_atom, screen_max, pmax_set: screening matrices
!> \param is_assoc_atomic_kind KS-matrix sparsity
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      The optimization is done via a simple binning procedure:
!>      In a first step the total amount of integrals in the system is calculated,
!>      taking into account the sparsity of the KS-matrix, the screening based
!>      on the Schwarz-Inequality and the density-matrix. Then, the procedure uses 
!>      different bin-sizes in order to optimize the load balance.
!>      At the end each process owns an unique array of inidices-ranges 
!>      that are used to decide whether a process has to calculate a certain
!>      bunch of integrals.
!> \par History
!>      06.2007 created [Manuel Guidon]
!>      08.2007 new parallel scheme [Manuel Guidon]
!>      11.2007 parallelize load balance on box_idx1 [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE load_balance_forces(x_data, eps_schwarz, eps_box, particle_set,max_set,para_env,&
                                 potential_parameter, screen_set, screen_atom, screen_max,&
                                 pmax_set, is_assoc_atomic_block_global, do_periodic, &
                                 load_balance_parameter, box_size, box_size_inv, kind_of, basis_parameter,&
                                 i_thread, n_threads, shm_far_field_boxes, error)
    TYPE(hfx_type), POINTER                  :: x_data
    REAL(dp), INTENT(IN)                     :: eps_schwarz, eps_box
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: max_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(hfx_potential_type), INTENT(IN)     :: potential_parameter
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: screen_set
    REAL(dp), DIMENSION(:, :), POINTER       :: screen_atom
    REAL(dp), INTENT(IN)                     :: screen_max
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: pmax_set
    INTEGER, DIMENSION(:, :)                 :: is_assoc_atomic_block_global
    LOGICAL                                  :: do_periodic
    TYPE(hfx_load_balance_type)              :: load_balance_parameter
    REAL(dp), INTENT(IN)                     :: box_size(3), box_size_inv(3)
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    INTEGER, INTENT(IN)                      :: i_thread, n_threads
    TYPE(hfx_box_scheme), DIMENSION(:, :), &
      POINTER                                :: shm_far_field_boxes
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'load_balance_forces', &
      routineP = moduleN//':'//routineN

    INTEGER :: data_from, dest, i, iatom, ibin, icpu, ikind, ipos, iset, j, &
      jatom, jkind, jset, katom, kkind, kset, latom, lkind, lset, &
      my_process_id, natom, nbins, nseta, nsetb, nsetc, nsetd, &
      pgf_radius_box, radius1, radius2, req(2), source, stat
    INTEGER(int_8) :: box_idx1, box_idx2, box_idx3, box_idx4, &
      current_counter, distribution_counter_end, distribution_counter_start, &
      int_per_cpu, mepos, n_processes, ncpu, optimal_cost, restore_counter, &
      tmp_counter1, total_integrals
    INTEGER(int_8), ALLOCATABLE, &
      DIMENSION(:), SAVE                     :: optimal_costs
    INTEGER(int_8), DIMENSION(:), POINTER    :: local_cost_matrix, recbuffer, &
                                                sendbuffer, swapbuffer
    INTEGER(int_8), DIMENSION(:), POINTER, &
      SAVE                                   :: cost_matrix
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_pos
    INTEGER, DIMENSION(:), POINTER           :: distribution_vector, nsgfa, &
                                                nsgfb, nsgfc, nsgfd
    INTEGER, DIMENSION(:), POINTER, SAVE     :: shm_distribution_vector
    INTEGER, SAVE                            :: rank_id_of_best_dist, &
                                                thread_id_of_best_dist
    LOGICAL                                  :: failure
    REAL(dp) :: center1(3), center2(3), center2_folded(3), Kbox, max_val1, &
      max_val2, omega2Inv, omega_box, pbc_center(3), pmax_entry, R1, R2, &
      Rbox, temp(3)
    TYPE(hfx_distribution), DIMENSION(:), &
      POINTER                                :: ptr_to_tmp_dist
    TYPE(hfx_distribution), &
      DIMENSION(:, :), POINTER               :: binned_dist, tmp_dist
    TYPE(hfx_distribution), &
      DIMENSION(:, :), POINTER, SAVE         :: full_dist

    IF(potential_parameter%omega /= 0.0_dp) THEN
      omega2Inv = 1.0_dp/(potential_parameter%omega**2)
    ELSE
      omega2Inv = 0.0_dp
    END IF

    natom=SIZE(particle_set)

    failure=.FALSE.
    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads
    IF(n_processes == 1) THEN
      ALLOCATE(tmp_dist(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1,1)%number_of_integrals = HUGE(tmp_dist(1,1)%number_of_integrals)
      tmp_dist(1,1)%istart = 1
      tmp_dist(1,1)%jstart = 1
      tmp_dist(1,1)%kstart = 1
      tmp_dist(1,1)%lstart = 1
      tmp_dist(1,1)%original_cpu_id = 1
      ptr_to_tmp_dist => tmp_dist(1,:)
      CALL hfx_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)
      DEALLOCATE(tmp_dist,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      total_integrals = 0
      my_process_id = para_env%mepos*n_threads + i_thread
      DO iatom=1,natom
        ikind = kind_of(iatom)
        nseta =  basis_parameter(ikind)%nset
        nsgfa => basis_parameter(ikind)%nsgf
        DO jatom=iatom,natom 
          box_idx1 = (iatom-1)*natom+jatom
          IF(MODULO(FNV64(box_idx1),n_processes) /= my_process_id) CYCLE
          IF (screen_atom(jatom,iatom) * screen_max < eps_schwarz) CYCLE
          jkind = kind_of(jatom)
          nsetb =  basis_parameter(jkind)%nset
          nsgfb => basis_parameter(jkind)%nsgf
          box_idx1 = (iatom-1)*natom+jatom
          DO katom=1,natom
            kkind = kind_of(katom)
            nsetc =  basis_parameter(kkind)%nset
            nsgfc => basis_parameter(kkind)%nsgf
            DO latom=katom,natom
              !All four centers equivalent => zero-contribution
              IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE

              IF (screen_atom(jatom,iatom) * screen_atom(latom,katom) < eps_schwarz) CYCLE
              IF(.NOT. (is_assoc_atomic_block_global(latom,iatom)  >=1 .AND. &
                        is_assoc_atomic_block_global(katom,iatom)  >=1 .AND. &
                        is_assoc_atomic_block_global(katom,jatom)  >=1 .AND. &
                        is_assoc_atomic_block_global(latom,jatom)  >=1 ) ) CYCLE
              IF((katom+latom<=iatom+jatom))  THEN
                IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                lkind = kind_of(latom)
                nsetd =  basis_parameter(lkind)%nset
                nsgfd => basis_parameter(lkind)%nsgf
                box_idx3 = (katom-1)*natom+latom

                DO iset = 1, nseta
                  DO jset = 1,nsetb
                    max_val1 = screen_set(jset,iset,jatom,iatom)
                    IF (max_val1 * screen_atom(latom,katom) < eps_schwarz) CYCLE
                    box_idx2 = (iset-1)*max_set+jset
                    center1 = shm_far_field_boxes(box_idx2, box_idx1)%center
                    R1 = shm_far_field_boxes(box_idx2, box_idx1)%R
                    radius1 = shm_far_field_boxes(box_idx2, box_idx1)%radius
                    DO kset = 1,nsetc
                      DO lset = 1,nsetd
                        pmax_entry = 2.0_dp*MAX(pmax_set(kset,iset,katom,iatom) * pmax_set(lset,jset,latom,jatom),&
                                                pmax_set(lset,iset,latom,iatom) * pmax_set(kset,jset,katom,jatom))
                        max_val2 = max_val1 * pmax_entry
                        IF(max_val2<eps_schwarz) CYCLE
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
                              max_val2 = 1.0_dp/Rbox*Kbox * pmax_entry 
                            CASE(do_hfx_potential_short)
                              omega_box = shm_far_field_boxes(box_idx2, box_idx1)%alphaInv + &
                                          shm_far_field_boxes(box_idx4, box_idx3)%alphaInv + &
                                          omega2Inv 
                              omega_box = 1.0_dp/omega_box
                              Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                     shm_far_field_boxes(box_idx4, box_idx3)%K
                              max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox * pmax_entry 
                            CASE(do_hfx_potential_mix_cl)
                              Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                     shm_far_field_boxes(box_idx4, box_idx3)%K
                              max_val2 = 1.0_dp/Rbox*Kbox * pmax_entry
                            CASE(do_hfx_potential_gaussian,do_hfx_potential_mix_lg, do_hfx_potential_id) 
                              max_val2 = HUGE(max_val2)
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
      ENDDO

      nbins = load_balance_parameter%nbins 
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
        ikind = kind_of(iatom)
        nseta =  basis_parameter(ikind)%nset
        nsgfa => basis_parameter(ikind)%nsgf
        DO jatom=iatom,natom
          box_idx1 = (iatom-1)*natom+jatom
          IF(MODULO(FNV64(box_idx1),n_processes) /= my_process_id) CYCLE
          jkind = kind_of(jatom)
          nsetb =  basis_parameter(jkind)%nset
          nsgfb => basis_parameter(jkind)%nsgf
          box_idx1 = (iatom-1)*natom+jatom
          DO katom=1,natom
            kkind = kind_of(katom)
            nsetc =  basis_parameter(kkind)%nset
            nsgfc => basis_parameter(kkind)%nsgf
            DO latom=katom,natom
              !All four centers equivalent => zero-contribution
              IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE
              IF((katom+latom<=iatom+jatom))  THEN
                IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                distribution_counter_end = distribution_counter_end + 1
            
                IF(binned_dist(1,ibin)%istart == 0 ) binned_dist(1,ibin)%istart = iatom
                IF(binned_dist(1,ibin)%jstart == 0 ) binned_dist(1,ibin)%jstart = jatom
                IF(binned_dist(1,ibin)%kstart == 0 ) binned_dist(1,ibin)%kstart = katom
                IF(binned_dist(1,ibin)%lstart == 0 ) binned_dist(1,ibin)%lstart = latom

                IF (screen_atom(jatom,iatom) * screen_max < eps_schwarz) CYCLE
                IF (screen_atom(jatom,iatom) * screen_atom(latom,katom) < eps_schwarz) CYCLE
                IF(.NOT. (is_assoc_atomic_block_global(latom,iatom)  >=1 .AND. &
                          is_assoc_atomic_block_global(katom,iatom)  >=1 .AND. &
                          is_assoc_atomic_block_global(katom,jatom)  >=1 .AND. &
                          is_assoc_atomic_block_global(latom,jatom)  >=1 ) ) CYCLE
 
                lkind = kind_of(latom)
                nsetd =  basis_parameter(lkind)%nset
                nsgfd => basis_parameter(lkind)%nsgf
                box_idx3 = (katom-1)*natom+latom
                DO iset = 1, nseta
                  DO jset = 1,nsetb
                   max_val1 = screen_set(jset,iset,jatom,iatom)
                    IF (max_val1 * screen_atom(latom,katom) < eps_schwarz) CYCLE
                    box_idx2 = (iset-1)*max_set+jset
                    center1 = shm_far_field_boxes(box_idx2, box_idx1)%center
                    R1 = shm_far_field_boxes(box_idx2, box_idx1)%R
                    radius1 = shm_far_field_boxes(box_idx2, box_idx1)%radius
                    DO kset = 1,nsetc
                     DO lset = 1,nsetd
                        pmax_entry = 2.0_dp*MAX(pmax_set(kset,iset,katom,iatom) * pmax_set(lset,jset,latom,jatom),&
                                                pmax_set(lset,iset,latom,iatom) * pmax_set(kset,jset,katom,jatom))

                        max_val2 = max_val1 * pmax_entry
                        IF(max_val2<eps_schwarz) CYCLE
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
                              max_val2 = 1.0_dp/Rbox*Kbox * pmax_entry
                            CASE(do_hfx_potential_short)
                              omega_box = shm_far_field_boxes(box_idx2, box_idx1)%alphaInv + &
                                          shm_far_field_boxes(box_idx4, box_idx3)%alphaInv + &
                                          omega2Inv
                              omega_box = 1.0_dp/omega_box
                              Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                     shm_far_field_boxes(box_idx4, box_idx3)%K
                              max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox * pmax_entry
                            CASE(do_hfx_potential_mix_cl)
                              Kbox = shm_far_field_boxes(box_idx2, box_idx1)%K * &
                                     shm_far_field_boxes(box_idx4, box_idx3)%K
                              max_val2 = 1.0_dp/Rbox*Kbox * pmax_entry
                            CASE(do_hfx_potential_gaussian,do_hfx_potential_mix_lg, do_hfx_potential_id) 
                              max_val2 = HUGE(max_val2)
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
      ENDDO
      !! Fill the last bin if necessary
      IF( current_counter /= 0 ) THEN  
        binned_dist(1,ibin)%cost = binned_dist(1,ibin)%cost + current_counter
        binned_dist(1,ibin)%number_of_integrals = distribution_counter_end-distribution_counter_start+1
      END IF

      !! store the original cpu_di
      binned_dist(1,:)%original_cpu_id = my_process_id + 1 !para_env%mepos + 1
    
      !! store all local results in a big cost matrix
!$OMP BARRIER
!$OMP MASTER
      ALLOCATE(cost_matrix(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      cost_matrix = 0
!$OMP END MASTER
!$OMP BARRIER
      DO icpu=1,ncpu
        IF(para_env%mepos+1 == icpu ) THEN
          DO i=1,nbins
            cost_matrix((icpu-1)*nbins*n_threads +i_thread*nbins + i) = binned_dist(1,i)%cost
          END DO
        END IF
      END DO

      mepos = para_env%mepos
!$OMP BARRIER
!$OMP MASTER
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
      !! Find an optimal distribution i.e. assign each element of the cost matrix to a certain process
      ALLOCATE(distribution_vector(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP MASTER
      ALLOCATE(optimal_costs(ncpu*n_threads), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      optimal_costs = 0
!$OMP END MASTER

!$OMP BARRIER
      optimal_cost = 0
      !! Optimize load on each node independently
      ALLOCATE(local_cost_matrix(SIZE(cost_matrix,1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      local_cost_matrix = cost_matrix

      CALL optimize_distribution(para_env, ncpu*nbins*n_threads, ncpu*n_threads, local_cost_matrix, distribution_vector, &
                                 load_balance_parameter, optimal_cost, my_process_id, i_thread, error)
      !! Find optimal load balance wrt the cost
      optimal_costs(my_process_id+1) = optimal_cost
!$OMP BARRIER
!$OMP MASTER
      CALL mp_sum_i8(optimal_costs(1:ncpu*n_threads), para_env%group)
      icpu = 0
      optimal_cost = HUGE(optimal_cost)
      !! Get process id of the best distribution
      DO i = 1,n_processes
        IF(optimal_costs(i) < optimal_cost ) THEN
          optimal_cost = optimal_costs(i)
          thread_id_of_best_dist = i
        END IF
      END DO
      rank_id_of_best_dist = (thread_id_of_best_dist-1)/n_threads
      ALLOCATE(shm_distribution_vector(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER
      IF(i_thread == MOD(thread_id_of_best_dist-1,n_threads)) THEN
        shm_distribution_vector = distribution_vector
        !! Broadcast the best distribution_vector
        CALL mp_bcast(shm_distribution_vector,rank_id_of_best_dist,para_env%group)
      END IF
!$OMP BARRIER
!$OMP MASTER
      !! Collect local data to global array
      ALLOCATE(full_dist(ncpu*n_threads,nbins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      full_dist(:,:)%istart = 0
      full_dist(:,:)%jstart = 0
      full_dist(:,:)%kstart = 0
      full_dist(:,:)%lstart = 0
      full_dist(:,:)%number_of_integrals = 0
      full_dist(:,:)%cost = 0
      full_dist(:,:)%original_cpu_id = 0
      full_dist(:,:)%time_first_scf = 0.0_dp
      full_dist(:,:)%time_other_scf = 0.0_dp
      full_dist(:,:)%time_forces = 0.0_dp
!$OMP END MASTER
!$OMP BARRIER
      DO mepos = 1,ncpu
        IF(para_env%mepos == mepos-1) THEN
          full_dist((mepos-1)*n_threads+i_thread+1,:) = binned_dist(1,:)
        END IF
      END DO

!$OMP BARRIER
!$OMP MASTER
      ALLOCATE(sendbuffer(7*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(recbuffer(7*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      mepos = para_env%mepos
      DO j=1,n_threads
        DO i=1,nbins
          sendbuffer((j-1)*7*nbins+(i-1)*7+1) = full_dist(mepos*n_threads+j,i)%istart
          sendbuffer((j-1)*7*nbins+(i-1)*7+2) = full_dist(mepos*n_threads+j,i)%jstart
          sendbuffer((j-1)*7*nbins+(i-1)*7+3) = full_dist(mepos*n_threads+j,i)%kstart
          sendbuffer((j-1)*7*nbins+(i-1)*7+4) = full_dist(mepos*n_threads+j,i)%lstart
          sendbuffer((j-1)*7*nbins+(i-1)*7+5) = full_dist(mepos*n_threads+j,i)%number_of_integrals
          sendbuffer((j-1)*7*nbins+(i-1)*7+6) = full_dist(mepos*n_threads+j,i)%cost
          sendbuffer((j-1)*7*nbins+(i-1)*7+7) = full_dist(mepos*n_threads+j,i)%original_cpu_id
        END DO
      END DO

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
            full_dist(data_from*n_threads+j,i)%istart = sendbuffer((j-1)*7*nbins+(i-1)*7+1)
            full_dist(data_from*n_threads+j,i)%jstart = sendbuffer((j-1)*7*nbins+(i-1)*7+2)
            full_dist(data_from*n_threads+j,i)%kstart = sendbuffer((j-1)*7*nbins+(i-1)*7+3)
            full_dist(data_from*n_threads+j,i)%lstart = sendbuffer((j-1)*7*nbins+(i-1)*7+4)
            full_dist(data_from*n_threads+j,i)%number_of_integrals = sendbuffer((j-1)*7*nbins+(i-1)*7+5)
            full_dist(data_from*n_threads+j,i)%cost = sendbuffer((j-1)*7*nbins+(i-1)*7+6)
            full_dist(data_from*n_threads+j,i)%original_cpu_id = sendbuffer((j-1)*7*nbins+(i-1)*7+7)
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
!$OMP END MASTER
!$OMP BARRIER
      !! reorder the distribution according to the distribution vector
      ALLOCATE(tmp_pos(ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_pos = 1
      ALLOCATE(tmp_dist(ncpu*n_threads,nbins*ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      tmp_dist(:,:)%istart = 0
      tmp_dist(:,:)%jstart = 0
      tmp_dist(:,:)%kstart = 0
      tmp_dist(:,:)%lstart = 0
      tmp_dist(:,:)%number_of_integrals = 0
      tmp_dist(:,:)%cost = 0
      tmp_dist(:,:)%original_cpu_id = 0
      tmp_dist(:,:)%time_first_scf = 0.0_dp
      tmp_dist(:,:)%time_other_scf = 0.0_dp
      tmp_dist(:,:)%time_forces = 0.0_dp
      DO mepos = 1,n_processes
        DO icpu= 1,n_processes
          DO i=1,nbins
            IF(shm_distribution_vector((icpu-1)*nbins + i) == mepos) THEN
              tmp_dist(mepos,tmp_pos(mepos)) = full_dist(icpu,i)
              tmp_pos(mepos) = tmp_pos(mepos) + 1
            END IF
          END DO
        END DO
      END DO

      !! Assign the load to each process
      NULLIFY(ptr_to_tmp_dist)
      DO mepos=1,n_processes
        IF(my_process_id == mepos -1) THEN
          ptr_to_tmp_dist => tmp_dist(mepos,1:tmp_pos(mepos)-1)
        END IF
      END DO
      CALL hfx_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)
!$OMP BARRIER
!$OMP MASTER
      DEALLOCATE(full_dist,optimal_costs,cost_matrix,shm_distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER
      DEALLOCATE(tmp_dist, tmp_pos, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(binned_dist, local_cost_matrix, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE load_balance_forces

! *****************************************************************************
!> \brief Distributes the computation of eri's to all available processes.
!> \param x_data Object that stores the indices array
!> \param potential_parameter short/longrange screening information 
!> \param screen_set , screen_atom, screen_max: screening matrices
!> \param is_assoc_atomic_kind KS-matrix sparsity
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      The cost matrix is given by the walltime for each bin that is measured during the calculation
!> \par History
!>      12.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE update_load_balance_energy(x_data, para_env, &
                                        load_balance_parameter, &
                                        i_thread, n_threads, error)

    TYPE(hfx_type), POINTER                  :: x_data
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(hfx_load_balance_type)              :: load_balance_parameter
    INTEGER, INTENT(IN)                      :: i_thread, n_threads
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_load_balance_energy', &
      routineP = moduleN//':'//routineN

    INTEGER :: data_from, dest, end_idx, i, ibin, icpu, iprocess, j, &
      my_bin_size, my_global_start_idx, my_process_id, nbins, req(2), source, &
      start_idx, stat
    INTEGER(int_8)                           :: mepos, n_processes, ncpu, &
                                                optimal_cost
    INTEGER(int_8), ALLOCATABLE, &
      DIMENSION(:), SAVE                     :: optimal_costs
    INTEGER(int_8), DIMENSION(:), POINTER    :: local_cost_matrix, recbuffer, &
                                                sendbuffer, swapbuffer
    INTEGER(int_8), DIMENSION(:), POINTER, &
      SAVE                                   :: cost_matrix
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_pos
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: bins_per_rank
    INTEGER, ALLOCATABLE, DIMENSION(:, :), &
      SAVE                                   :: bin_histogram
    INTEGER, DIMENSION(:), POINTER           :: distribution_vector
    INTEGER, DIMENSION(:), POINTER, SAVE     :: shm_distribution_vector
    INTEGER, SAVE                            :: max_bin_size, &
                                                rank_id_of_best_dist, &
                                                thread_id_of_best_dist
    LOGICAL                                  :: failure
    TYPE(hfx_distribution), DIMENSION(:), &
      POINTER                                :: binned_dist, ptr_to_tmp_dist
    TYPE(hfx_distribution), &
      DIMENSION(:, :), POINTER               :: tmp_dist
    TYPE(hfx_distribution), &
      DIMENSION(:, :), POINTER, SAVE         :: full_dist

    failure=.FALSE.

    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads
    !! If there is only 1 cpu skip the binning
    IF(n_processes == 1) THEN
      ALLOCATE(tmp_dist(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1,1)%number_of_integrals = HUGE(tmp_dist(1,1)%number_of_integrals)
      tmp_dist(1,1)%istart = 1
      tmp_dist(1,1)%jstart = 1
      tmp_dist(1,1)%kstart = 1
      tmp_dist(1,1)%lstart = 1
      tmp_dist(1,1)%original_cpu_id = 1
      ptr_to_tmp_dist => tmp_dist(1,:)
      CALL hfx_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)
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
      my_bin_size = SIZE(x_data%distribution_energy)
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
      binned_dist = x_data%distribution_energy
      DO ibin=1,my_bin_size
        IF(binned_dist(ibin)%number_of_integrals == 0 ) THEN
          binned_dist(ibin)%cost = 0
        ELSE
          binned_dist(ibin)%cost = INT((binned_dist(ibin)%time_first_scf + &
                                   binned_dist(ibin)%time_other_scf)*10000.0_dp, int_8)
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
      DO icpu=1,ncpu
        IF(para_env%mepos+1 == icpu ) THEN
          DO i=1,my_bin_size
            cost_matrix(my_global_start_idx + i) = binned_dist(i)%cost
          END DO
        END IF
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
!$OMP END MASTER
!$OMP BARRIER
      !! Find an optimal distribution i.e. assign each element of the cost matrix to a certain process
      ALLOCATE(distribution_vector(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP MASTER
      ALLOCATE(optimal_costs(ncpu*n_threads), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      optimal_costs = 0
!$OMP END MASTER
!$OMP BARRIER
      optimal_cost = 0
      !! Optimize load on each node independently
      ALLOCATE(local_cost_matrix(SIZE(cost_matrix,1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      local_cost_matrix = cost_matrix
      CALL optimize_distribution(para_env, ncpu*nbins*n_threads, ncpu*n_threads, local_cost_matrix, distribution_vector, &
                                 load_balance_parameter, optimal_cost, my_process_id, i_thread, error)
      
      !! Find optimal load balance wrt the cost
      optimal_costs(my_process_id+1) = optimal_cost
!$OMP BARRIER
!$OMP MASTER
      CALL mp_sum_i8(optimal_costs(1:ncpu*n_threads), para_env%group)
      icpu = 0
      optimal_cost = HUGE(optimal_cost)
      !! Get process id of the best distribution
      DO i = 1,n_processes
        IF(optimal_costs(i) < optimal_cost ) THEN
          optimal_cost = optimal_costs(i)
          thread_id_of_best_dist = i
        END IF
      END DO
      rank_id_of_best_dist = (thread_id_of_best_dist-1)/n_threads
      ALLOCATE(shm_distribution_vector(ncpu*nbins*n_threads), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER

!$OMP BARRIER
      IF(i_thread == MOD(thread_id_of_best_dist-1,n_threads)) THEN
        shm_distribution_vector = distribution_vector
        !! Broadcast the best distribution_vector
        CALL mp_bcast(shm_distribution_vector,rank_id_of_best_dist,para_env%group)
      END IF
!$OMP BARRIER

!$OMP MASTER
      !! Collect local data to global array
      ALLOCATE(full_dist(ncpu*n_threads,max_bin_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      full_dist(:,:)%istart = 0
      full_dist(:,:)%jstart = 0
      full_dist(:,:)%kstart = 0
      full_dist(:,:)%lstart = 0
      full_dist(:,:)%number_of_integrals = 0
      full_dist(:,:)%cost = 0
      full_dist(:,:)%original_cpu_id = 0
      full_dist(:,:)%time_first_scf = 0.0_dp
      full_dist(:,:)%time_other_scf = 0.0_dp
      full_dist(:,:)%time_forces = 0.0_dp
!$OMP END MASTER

!$OMP BARRIER
      DO mepos = 1,ncpu
        IF(para_env%mepos == mepos-1) THEN
           full_dist((mepos-1)*n_threads+i_thread+1,1:my_bin_size) = binned_dist(1:my_bin_size)
        END IF
      END DO
!$OMP BARRIER
!$OMP MASTER
      ALLOCATE(sendbuffer(7*max_bin_size*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(recbuffer(7*max_bin_size*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      mepos = para_env%mepos
      DO j=1,n_threads
        DO i=1,max_bin_size
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+1) = full_dist(mepos*n_threads+j,i)%istart
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+2) = full_dist(mepos*n_threads+j,i)%jstart
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+3) = full_dist(mepos*n_threads+j,i)%kstart
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+4) = full_dist(mepos*n_threads+j,i)%lstart
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+5) = full_dist(mepos*n_threads+j,i)%number_of_integrals
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+6) = full_dist(mepos*n_threads+j,i)%cost
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+7) = full_dist(mepos*n_threads+j,i)%original_cpu_id
        END DO
      END DO
      dest   = MODULO(mepos+1,ncpu)
      source = MODULO(mepos-1,ncpu)
      DO icpu = 0, ncpu-1
        IF (icpu .NE. ncpu-1) THEN
          CALL mp_isendrecv(sendbuffer,dest,recbuffer,source, &
                            para_env%group,req(1),req(2),13)
        ENDIF
        data_from = MODULO(mepos - icpu,ncpu)
        DO j=1,n_threads
          DO i=1,max_bin_size
            full_dist(data_from*n_threads+j,i)%istart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+1)
            full_dist(data_from*n_threads+j,i)%jstart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+2)
            full_dist(data_from*n_threads+j,i)%kstart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+3)
            full_dist(data_from*n_threads+j,i)%lstart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+4)
            full_dist(data_from*n_threads+j,i)%number_of_integrals = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+5)
            full_dist(data_from*n_threads+j,i)%cost = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+6)
            full_dist(data_from*n_threads+j,i)%original_cpu_id = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+7)
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
!$OMP END MASTER
!$OMP BARRIER
      !! reorder the distribution according to the distribution vector
      ALLOCATE(tmp_pos(ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_pos = 1
      ALLOCATE(tmp_dist(ncpu*n_threads,nbins*ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(:,:)%istart = 0
      tmp_dist(:,:)%jstart = 0
      tmp_dist(:,:)%kstart = 0
      tmp_dist(:,:)%lstart = 0
      tmp_dist(:,:)%number_of_integrals = 0
      tmp_dist(:,:)%cost = 0
      tmp_dist(:,:)%original_cpu_id = 0
      tmp_dist(:,:)%time_first_scf = 0.0_dp
      tmp_dist(:,:)%time_other_scf = 0.0_dp
      tmp_dist(:,:)%time_forces = 0.0_dp

      DO mepos = 1,n_processes
        DO icpu= 1,n_processes
          DO i = 1,bin_histogram(icpu,1)
            IF(shm_distribution_vector(bin_histogram(icpu,2)- bin_histogram(icpu,1) + i) == mepos) THEN
              tmp_dist(mepos,tmp_pos(mepos)) = full_dist(icpu,i)
              tmp_pos(mepos) = tmp_pos(mepos) + 1
            END IF
          END DO
        END DO
      END DO
      
      !! Assign the load to each process
      NULLIFY(ptr_to_tmp_dist) 
      DO mepos=1,n_processes
        IF(my_process_id == mepos -1) THEN
          ptr_to_tmp_dist => tmp_dist(mepos,1:tmp_pos(mepos)-1)
        END IF
      END DO

      CALL hfx_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)

!$OMP BARRIER
!$OMP MASTER
      DEALLOCATE(full_dist, optimal_costs, cost_matrix, shm_distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(bins_per_rank,bin_histogram, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER
      DEALLOCATE(tmp_dist,  tmp_pos, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(binned_dist, local_cost_matrix, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE update_load_balance_energy

! *****************************************************************************
!> \brief Distributes the computation of eri's to all available processes.
!> \param x_data Object that stores the indices array
!> \param potential_parameter short/longrange screening information 
!> \param screen_set , screen_atom, screen_max: screening matrices
!> \param is_assoc_atomic_kind KS-matrix sparsity
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      The cost matrix is given by the walltime for each bin that is measured during the calculation
!> \par History
!>      12.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE update_load_balance_forces(x_data, para_env,&
                                        load_balance_parameter,i_thread, n_threads, error)

    TYPE(hfx_type), POINTER                  :: x_data
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(hfx_load_balance_type)              :: load_balance_parameter
    INTEGER, INTENT(IN)                      :: i_thread, n_threads
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_load_balance_forces', &
      routineP = moduleN//':'//routineN

    INTEGER :: data_from, dest, end_idx, i, ibin, icpu, iprocess, j, &
      my_bin_size, my_global_start_idx, my_process_id, nbins, req(2), source, &
      start_idx, stat
    INTEGER(int_8)                           :: mepos, n_processes, ncpu, &
                                                optimal_cost
    INTEGER(int_8), ALLOCATABLE, &
      DIMENSION(:), SAVE                     :: optimal_costs
    INTEGER(int_8), DIMENSION(:), POINTER    :: local_cost_matrix, recbuffer, &
                                                sendbuffer, swapbuffer
    INTEGER(int_8), DIMENSION(:), POINTER, &
      SAVE                                   :: cost_matrix
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_pos
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: bins_per_rank
    INTEGER, ALLOCATABLE, DIMENSION(:, :), &
      SAVE                                   :: bin_histogram
    INTEGER, DIMENSION(:), POINTER           :: distribution_vector
    INTEGER, DIMENSION(:), POINTER, SAVE     :: shm_distribution_vector
    INTEGER, SAVE                            :: max_bin_size, &
                                                rank_id_of_best_dist, &
                                                thread_id_of_best_dist
    LOGICAL                                  :: failure
    TYPE(hfx_distribution), DIMENSION(:), &
      POINTER                                :: binned_dist, ptr_to_tmp_dist
    TYPE(hfx_distribution), &
      DIMENSION(:, :), POINTER               :: tmp_dist
    TYPE(hfx_distribution), &
      DIMENSION(:, :), POINTER, SAVE         :: full_dist

    failure=.FALSE.

    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads
    !! If there is only 1 cpu skip the binning
    IF(n_processes == 1) THEN
      ALLOCATE(tmp_dist(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1,1)%number_of_integrals = HUGE(tmp_dist(1,1)%number_of_integrals)
      tmp_dist(1,1)%istart = 1
      tmp_dist(1,1)%jstart = 1
      tmp_dist(1,1)%kstart = 1
      tmp_dist(1,1)%lstart = 1
      tmp_dist(1,1)%original_cpu_id = 1
      ptr_to_tmp_dist => tmp_dist(1,:)
      CALL hfx_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)
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
      my_bin_size = SIZE(x_data%distribution_forces)
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
      binned_dist = x_data%distribution_forces
      DO ibin=1,my_bin_size
        IF(binned_dist(ibin)%number_of_integrals == 0 ) THEN
          binned_dist(ibin)%cost = 0
        ELSE
          binned_dist(ibin)%cost = INT((binned_dist(ibin)%time_forces)*10000.0_dp, int_8)
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
      DO icpu=1,ncpu
        IF(para_env%mepos+1 == icpu ) THEN
          DO i=1,my_bin_size
            cost_matrix(my_global_start_idx + i) = binned_dist(i)%cost
          END DO
        END IF
      END DO

      mepos = para_env%mepos
!$OMP BARRIER
!$OMP MASTER
      ALLOCATE(bins_per_rank(ncpu),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      bins_per_rank = 0
      DO icpu=1,ncpu
        bins_per_rank(icpu) = SUM(bin_histogram((icpu-1)*n_threads+1:(icpu-1)*n_threads+n_threads,1))
      END DO
      sendbuffer(1:bins_per_rank(para_env%mepos+1)) = &
        cost_matrix(my_global_start_idx+1:my_global_start_idx+bins_per_rank(para_env%mepos+1))

      dest  =MODULO(mepos+1,ncpu)
      source=MODULO(mepos-1,ncpu)
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
!$OMP END MASTER
!$OMP BARRIER
      !! Find an optimal distribution i.e. assign each element of the cost matrix to a certain process
      ALLOCATE(distribution_vector(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP MASTER
      ALLOCATE(optimal_costs(ncpu*n_threads), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      optimal_costs = 0
!$OMP END MASTER
!$OMP BARRIER
      optimal_cost = 0
      !! Optimize load on each node independently
      ALLOCATE(local_cost_matrix(SIZE(cost_matrix,1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      local_cost_matrix = cost_matrix
      CALL optimize_distribution(para_env, ncpu*nbins*n_threads, ncpu*n_threads, local_cost_matrix, distribution_vector, &
                                 load_balance_parameter, optimal_cost, my_process_id, i_thread, error)
      
      !! Find optimal load balance wrt the cost
      optimal_costs(my_process_id+1) = optimal_cost
!$OMP BARRIER
!$OMP MASTER
      CALL mp_sum_i8(optimal_costs(1:ncpu*n_threads), para_env%group)
      icpu = 0
      optimal_cost = HUGE(optimal_cost)
      !! Get process id of the best distribution
      DO i = 1,n_processes
        IF(optimal_costs(i) < optimal_cost ) THEN
          optimal_cost = optimal_costs(i)
          thread_id_of_best_dist = i
        END IF
      END DO
      rank_id_of_best_dist = (thread_id_of_best_dist-1)/n_threads
      ALLOCATE(shm_distribution_vector(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER

!$OMP BARRIER
      IF(i_thread == MOD(thread_id_of_best_dist-1,n_threads)) THEN
        shm_distribution_vector = distribution_vector
        !! Broadcast the best distribution_vector
        CALL mp_bcast(shm_distribution_vector,rank_id_of_best_dist,para_env%group)
      END IF
!$OMP BARRIER

!$OMP MASTER
      !! Collect local data to global array
      ALLOCATE(full_dist(ncpu*n_threads,max_bin_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      full_dist(:,:)%istart = 0
      full_dist(:,:)%jstart = 0
      full_dist(:,:)%kstart = 0
      full_dist(:,:)%lstart = 0
      full_dist(:,:)%number_of_integrals = 0
      full_dist(:,:)%cost = 0
      full_dist(:,:)%original_cpu_id = 0
      full_dist(:,:)%time_first_scf = 0.0_dp
      full_dist(:,:)%time_other_scf = 0.0_dp
      full_dist(:,:)%time_forces = 0.0_dp
!$OMP END MASTER

!$OMP BARRIER
      DO mepos = 1,ncpu
        IF(para_env%mepos == mepos-1) THEN
           full_dist((mepos-1)*n_threads+i_thread+1,1:my_bin_size) = binned_dist(1:my_bin_size)
        END IF
      END DO
!$OMP BARRIER
!$OMP MASTER
      ALLOCATE(sendbuffer(7*max_bin_size*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(recbuffer(7*max_bin_size*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      mepos = para_env%mepos
      DO j=1,n_threads
        DO i=1,max_bin_size
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+1) = full_dist(mepos*n_threads+j,i)%istart
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+2) = full_dist(mepos*n_threads+j,i)%jstart
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+3) = full_dist(mepos*n_threads+j,i)%kstart
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+4) = full_dist(mepos*n_threads+j,i)%lstart
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+5) = full_dist(mepos*n_threads+j,i)%number_of_integrals
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+6) = full_dist(mepos*n_threads+j,i)%cost
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+7) = full_dist(mepos*n_threads+j,i)%original_cpu_id
        END DO
      END DO
      dest   = MODULO(mepos+1,ncpu)
      source = MODULO(mepos-1,ncpu)
      DO icpu = 0, ncpu-1
        IF (icpu .NE. ncpu-1) THEN
          CALL mp_isendrecv(sendbuffer,dest,recbuffer,source, &
                            para_env%group,req(1),req(2),13)
        ENDIF
        data_from = MODULO(mepos - icpu,ncpu)
        DO j=1,n_threads
          DO i=1,max_bin_size
            full_dist(data_from*n_threads+j,i)%istart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+1)
            full_dist(data_from*n_threads+j,i)%jstart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+2)
            full_dist(data_from*n_threads+j,i)%kstart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+3)
            full_dist(data_from*n_threads+j,i)%lstart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+4)
            full_dist(data_from*n_threads+j,i)%number_of_integrals = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+5)
            full_dist(data_from*n_threads+j,i)%cost = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+6)
            full_dist(data_from*n_threads+j,i)%original_cpu_id = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+7)
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
!$OMP END MASTER
!$OMP BARRIER
      !! reorder the distribution according to the distribution vector
      ALLOCATE(tmp_pos(ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_pos = 1
      ALLOCATE(tmp_dist(ncpu*n_threads,nbins*ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(:,:)%istart = 0
      tmp_dist(:,:)%jstart = 0
      tmp_dist(:,:)%kstart = 0
      tmp_dist(:,:)%lstart = 0
      tmp_dist(:,:)%number_of_integrals = 0
      tmp_dist(:,:)%cost = 0
      tmp_dist(:,:)%original_cpu_id = 0
      tmp_dist(:,:)%time_first_scf = 0.0_dp
      tmp_dist(:,:)%time_other_scf = 0.0_dp
      tmp_dist(:,:)%time_forces = 0.0_dp

      DO mepos = 1,n_processes
        DO icpu= 1,n_processes
          DO i = 1,bin_histogram(icpu,1)
            IF(shm_distribution_vector(bin_histogram(icpu,2)- bin_histogram(icpu,1) + i) == mepos) THEN
              tmp_dist(mepos,tmp_pos(mepos)) = full_dist(icpu,i)
              tmp_pos(mepos) = tmp_pos(mepos) + 1
            END IF
          END DO
        END DO
      END DO
      
      !! Assign the load to each process
      NULLIFY(ptr_to_tmp_dist) 
      DO mepos=1,n_processes
        IF(my_process_id == mepos -1) THEN
          ptr_to_tmp_dist => tmp_dist(mepos,1:tmp_pos(mepos)-1)
        END IF
      END DO

      CALL hfx_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)

!$OMP BARRIER
!$OMP MASTER
      DEALLOCATE(full_dist, optimal_costs, cost_matrix, shm_distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(bins_per_rank,bin_histogram, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER
      DEALLOCATE(tmp_dist,  tmp_pos, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(binned_dist, local_cost_matrix, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE update_load_balance_forces

END MODULE hfx_load_balance_methods 
