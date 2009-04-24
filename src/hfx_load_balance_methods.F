!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines for optimizing load balance between processes in HFX calculations 
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE hfx_load_balance_methods 
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE hfx_types,                       ONLY: &
       hfx_basis_type, hfx_distribution, hfx_load_balance_type, &
       hfx_potential_type, hfx_screen_coeff_type, hfx_set_distr_energy, &
       hfx_set_distr_forces, hfx_type, pair_list_type, pair_set_list_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE mathconstants
  USE message_passing,                 ONLY: mp_isendrecv,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync,&
                                             mp_waitall
  USE orbital_pointers
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC load_balance_energy, load_balance_forces, update_load_balance_energy, update_load_balance_forces, build_pair_list

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_load_balance_methods'

!***
  
  CONTAINS

  SUBROUTINE build_pair_list(list,set_list,i_start,i_end,j_start,j_end,kind_of,basis_parameter,particle_set, &
                             do_periodic,coeffs_set,coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell)
  
    TYPE(pair_list_type), INTENT(OUT)        :: list
    TYPE(pair_set_list_type), DIMENSION(*), &
      INTENT(OUT)                            :: set_list
    INTEGER, INTENT(IN)                      :: i_start, i_end, j_start, j_end
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    LOGICAL, INTENT(IN)                      :: do_periodic
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: coeffs_kind
    REAL(KIND=dp), INTENT(IN)                :: coeffs_kind_max0, &
                                                log10_eps_schwarz
    TYPE(cell_type), POINTER                 :: cell

    INTEGER                                  :: iatom, ikind, iset, jatom, &
                                                jkind, jset, n_element, &
                                                nset_ij, nseta, nsetb
    REAL(KIND=dp)                            :: rab2
    REAL(KIND=dp), DIMENSION(3)              :: B11, pbc_B, ra, rb, temp

     n_element=0
     nset_ij=0
  
     DO iatom=i_start,i_end
       ikind = kind_of(iatom)
       nseta =  basis_parameter(ikind)%nset
       ra=particle_set(iatom)%r(:)
       DO jatom=j_start,j_end
         IF(jatom<iatom) CYCLE
         jkind = kind_of(jatom)
         nsetb =  basis_parameter(jkind)%nset
         rb=particle_set(jatom)%r(:)

         IF( do_periodic ) THEN
           temp = rb -ra
           pbc_B = pbc(temp,cell)
           B11 = ra + pbc_B
           rab2 = (ra(1)-B11(1))**2+(ra(2)-B11(2))**2+(ra(3)-B11(3))**2
         ELSE
           rab2 = (ra(1)-rb(1))**2+(ra(2)-rb(2))**2+(ra(3)-rb(3))**2
           B11 = ra - rb
         END IF

         IF( (coeffs_kind(jkind,ikind)%x(1)*rab2+&
              coeffs_kind(jkind,ikind)%x(2) ) + coeffs_kind_max0 < log10_eps_schwarz ) CYCLE

         n_element=n_element+1 
         list%elements(n_element)%pair(1)=iatom
         list%elements(n_element)%pair(2)=jatom
         list%elements(n_element)%kind_pair(1)=ikind
         list%elements(n_element)%kind_pair(2)=jkind
         list%elements(n_element)%r1   =ra
         list%elements(n_element)%r2   =rb
         list%elements(n_element)%vec  =B11
         list%elements(n_element)%dist2=rab2

         ! build a list of guaranteed overlapping sets
         list%elements(n_element)%set_bounds(1)=nset_ij+1
         DO iset = 1, nseta
            DO jset = 1,nsetb
               IF (coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + coeffs_set(jset,iset,jkind,ikind)%x(2) + &
                   coeffs_kind_max0 < log10_eps_schwarz ) CYCLE
               nset_ij=nset_ij+1
               set_list(nset_ij)%pair(1)=iset
               set_list(nset_ij)%pair(2)=jset
            ENDDO
         ENDDO
         list%elements(n_element)%set_bounds(2)=nset_ij

       ENDDO
     ENDDO
  
     list%n_element=n_element
  
  END SUBROUTINE build_pair_list

! *****************************************************************************
!> \brief Distributes the computation of eri's to all available processes.
!> \param x_data Object that stores the indices array
!> \param eps_schwarz screening parameter
!> \param particle_set , atomic_kind_set, para_env
!> \param max_set Maximum number of set to be considered
!> \param para_env para_env
!> \param potential_parameter short/longrange screening information 
!> \param coeffs_set, coeffs_kind: screening functions
!> \param is_assoc_atomic_kind KS-matrix sparsity
!> \param do_periodic flag for periodicity
!> \param load_balance_parameter Paramters for Monte-Carlo routines
!> \param kind_of helper array for mapping
!> \param basis_parameter Basis set parameters
!> \param pmax_set Initial screening matrix
!> \param i_thread Process ID of current Thread
!> \param n_thread Total Number of Threads
!> \param cell cell
!> \param do_p_screening Flag for initial p screening
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
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
!> \par History
!>      06.2007 created [Manuel Guidon]
!>      08.2007 new parallel scheme [Manuel Guidon]
!>      09.2007 new 'modulo' parellel scheme and Monte Carlo step [Manuel Guidon]
!>      11.2007 parallelize load balance on box_idx1 [Manuel Guidon]
!>      02.2009 completely refactored [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE load_balance_energy(x_data, eps_schwarz, particle_set,max_set,para_env,&
                                 potential_parameter, coeffs_set, coeffs_kind, &
                                 is_assoc_atomic_block_global, do_periodic,&
                                 load_balance_parameter, kind_of, basis_parameter, pmax_set, &
                                 i_thread, n_threads, cell, &
                                 do_p_screening, error)
    TYPE(hfx_type), POINTER                  :: x_data
    REAL(dp), INTENT(IN)                     :: eps_schwarz
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: max_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(hfx_potential_type), INTENT(IN)     :: potential_parameter
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: coeffs_kind
    INTEGER, DIMENSION(:, :)                 :: is_assoc_atomic_block_global
    LOGICAL                                  :: do_periodic
    TYPE(hfx_load_balance_type)              :: load_balance_parameter
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: pmax_set
    INTEGER, INTENT(IN)                      :: i_thread, n_threads
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(IN)                      :: do_p_screening
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'load_balance_energy', &
      routineP = moduleN//':'//routineN

    INTEGER :: block_size, data_from, dest, handle, handle_inner, i, &
      i_list_ij, i_list_kl, i_set_list_ij, i_set_list_ij_start, &
      i_set_list_ij_stop, i_set_list_kl, i_set_list_kl_start, &
      i_set_list_kl_stop, iatom, iatom_block, ibin, icpu, ikind, ipos, iset, &
      j, jatom, jatom_block, jkind, jset, katom, katom_block, kkind, kset, &
      latom, latom_block, lkind, lset, my_process_id, natom, nbins, req(2), &
      source, stat
    INTEGER(int_8) :: atom_block, current_counter, distribution_counter_end, &
      distribution_counter_start, int_per_cpu, mepos, n_processes, nblocks, &
      ncpu, restore_counter, tmp_block, total_integrals
    INTEGER(int_8), DIMENSION(:), POINTER    :: local_cost_matrix, recbuffer, &
                                                sendbuffer, swapbuffer
    INTEGER(int_8), DIMENSION(:), POINTER, &
      SAVE                                   :: cost_matrix
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_pos
    INTEGER, DIMENSION(:), POINTER           :: distribution_vector, nsgfa, &
                                                nsgfb, nsgfc, nsgfd
    INTEGER, DIMENSION(:), POINTER, SAVE     :: shm_distribution_vector
    LOGICAL                                  :: failure
    REAL(dp) :: coeffs_kind_max0, log10_eps_schwarz, max_val1, max_val2, &
      omega2Inv, pmax_entry, ra(3), rab2, rb(3), rc(3), rcd2, rd(3)
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

    IF(potential_parameter%omega /= 0.0_dp) THEN
      omega2Inv = 1.0_dp/(potential_parameter%omega**2)
    ELSE
      omega2Inv = 0.0_dp
    END IF
    natom=SIZE(particle_set)
    failure=.FALSE.
    log10_eps_schwarz = LOG10(eps_schwarz)
    coeffs_kind_max0=MAXVAL(coeffs_kind(:,:)%x(2))
    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads

    ALLOCATE(set_list_ij((max_set*load_balance_parameter%block_size)**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(set_list_kl((max_set*load_balance_parameter%block_size)**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !! If there is only 1 cpu skip the binning
    IF(n_processes == 1) THEN
      ALLOCATE(tmp_dist(1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1)%number_of_atom_quartets = HUGE(tmp_dist(1)%number_of_atom_quartets)
      tmp_dist(1)%istart = 0_int_8
      tmp_dist(1)%jstart = 1_int_8
      tmp_dist(1)%kstart = 1_int_8
      tmp_dist(1)%lstart = 1_int_8
      tmp_dist(1)%original_cpu_id = 1_int_8
      ptr_to_tmp_dist => tmp_dist(:)
      CALL hfx_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)
      DEALLOCATE(tmp_dist,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      !! Calculate total numbers of integrals that have to be calculated (wrt screening and symmetry)
!$OMP BARRIER
!$OMP MASTER
    CALL timeset(routineN//"_count",handle_inner)
!$OMP END MASTER
!$OMP BARRIER

      total_integrals = 0
      my_process_id = para_env%mepos*n_threads + i_thread
      block_size = load_balance_parameter%block_size
      nblocks = MAX((natom+block_size-1)/block_size,1)
      DO atom_block = my_process_id,nblocks**4-1,n_processes
        latom_block = MODULO(atom_block,nblocks)+1
        tmp_block = atom_block/nblocks
        katom_block = MODULO(tmp_block,nblocks)+1
        IF(latom_block<katom_block) CYCLE
        tmp_block = tmp_block/nblocks
        jatom_block = MODULO(tmp_block,nblocks)+1
        tmp_block = tmp_block/nblocks
        iatom_block = MODULO(tmp_block,nblocks)+1
        IF(jatom_block<iatom_block) CYCLE

        CALL build_pair_list(list_ij,set_list_ij,(iatom_block-1)*block_size +1,MIN(iatom_block*block_size,natom), &
                             (jatom_block-1)*block_size +1, MIN(jatom_block*block_size,natom),&
                             kind_of,basis_parameter,particle_set, &
                             do_periodic,coeffs_set,coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell)     

        CALL build_pair_list(list_kl,set_list_kl,(katom_block-1)*block_size +1,MIN(katom_block*block_size,natom), &
                             (latom_block-1)*block_size +1, MIN(latom_block*block_size,natom),&
                             kind_of,basis_parameter,particle_set, &
                             do_periodic,coeffs_set,coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell)     

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

          nsgfa => basis_parameter(ikind)%nsgf
          nsgfb => basis_parameter(jkind)%nsgf

          DO i_list_kl=1,list_kl%n_element

            katom=list_kl%elements(i_list_kl)%pair(1)
            latom=list_kl%elements(i_list_kl)%pair(2)

            IF( .NOT. (katom+latom<=iatom+jatom))  CYCLE
            IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE

            i_set_list_kl_start=list_kl%elements(i_list_kl)%set_bounds(1)
            i_set_list_kl_stop=list_kl%elements(i_list_kl)%set_bounds(2)
            kkind=list_kl%elements(i_list_kl)%kind_pair(1)
            lkind=list_kl%elements(i_list_kl)%kind_pair(2)
            rc=list_kl%elements(i_list_kl)%r1
            rd=list_kl%elements(i_list_kl)%r2
            rcd2=list_kl%elements(i_list_kl)%dist2

            nsgfc => basis_parameter(kkind)%nsgf
            nsgfd => basis_parameter(lkind)%nsgf

            IF( (coeffs_kind(jkind,ikind)%x(1)*rab2+&
                 coeffs_kind(jkind,ikind)%x(2) ) + &
                (coeffs_kind(lkind,kkind)%x(1)*rcd2+&
                 coeffs_kind(lkind,kkind)%x(2)) < log10_eps_schwarz ) CYCLE

            IF(.NOT. (is_assoc_atomic_block_global(latom,iatom)  >=1 .AND. &
                      is_assoc_atomic_block_global(katom,iatom)  >=1 .AND. &
                      is_assoc_atomic_block_global(katom,jatom)  >=1 .AND. &
                      is_assoc_atomic_block_global(latom,jatom)  >=1 ) ) CYCLE 

            DO i_set_list_ij=i_set_list_ij_start, i_set_list_ij_stop
              iset=set_list_ij(i_set_list_ij)%pair(1)
              jset=set_list_ij(i_set_list_ij)%pair(2)

              max_val1 = coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + &
                         coeffs_set(jset,iset,jkind,ikind)%x(2)
              IF( max_val1 + coeffs_kind(lkind,kkind)%x(1)*rcd2 + &
                             coeffs_kind(lkind,kkind)%x(2)  <log10_eps_schwarz) CYCLE

              DO i_set_list_kl=i_set_list_kl_start, i_set_list_kl_stop
                kset=set_list_kl(i_set_list_kl)%pair(1)
                lset=set_list_kl(i_set_list_kl)%pair(2)

                max_val2 = max_val1 + (coeffs_set(lset,kset,lkind,kkind)%x(1)*rcd2 + &
                                       coeffs_set(lset,kset,lkind,kkind)%x(2) )

                IF(max_val2 <log10_eps_schwarz ) CYCLE
                IF( do_p_screening ) THEN
                  pmax_entry = MAX(pmax_set(kset,iset,katom,iatom), &
                                   pmax_set(lset,jset,latom,jatom), &
                                   pmax_set(lset,iset,latom,iatom), &
                                   pmax_set(kset,jset,katom,jatom))
                ELSE
                  pmax_entry = 0.0_dp
                END IF
                max_val2= max_val2 + pmax_entry
                IF( max_val2< log10_eps_schwarz) CYCLE
                total_integrals = total_integrals + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
              ENDDO ! i_set_list_kl
            ENDDO ! i_set_list_ij
          ENDDO ! i_list_kl
        ENDDO ! i_list_ij
      ENDDO ! atom_block

      !! Bin the total amount of integrals on each cpu
      nbins = load_balance_parameter%nbins
      !floor value
      int_per_cpu = (total_integrals + nbins -1 ) / (nbins) 

!$OMP BARRIER
!$OMP MASTER
    CALL timestop(handle_inner)
    CALL timeset(routineN//"_bin",handle_inner)
!$OMP END MASTER
!$OMP BARRIER

      ALLOCATE(binned_dist(nbins), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      binned_dist(:)%istart = -1_int_8
      binned_dist(:)%jstart = 0_int_8
      binned_dist(:)%kstart = 0_int_8
      binned_dist(:)%lstart = 0_int_8
      binned_dist(:)%number_of_atom_quartets = 0_int_8
      binned_dist(:)%cost = 0_int_8

      binned_dist(:)%original_cpu_id = 0_int_8
      binned_dist(:)%time_first_scf = 0.0_dp
      binned_dist(:)%time_other_scf =0.0_dp
      binned_dist(:)%time_forces = 0.0_dp
 
      current_counter = 0
      restore_counter = 0
      ipos = 1
      mepos = 1
      distribution_counter_start = 1
      distribution_counter_end = 0
      ibin = 1

      DO atom_block = my_process_id,nblocks**4-1,n_processes
        latom_block = MODULO(atom_block,nblocks)+1
        tmp_block = atom_block/nblocks
        katom_block = MODULO(tmp_block,nblocks)+1
        IF(latom_block<katom_block) CYCLE
        tmp_block = tmp_block/nblocks
        jatom_block = MODULO(tmp_block,nblocks)+1
        tmp_block = tmp_block/nblocks
        iatom_block = MODULO(tmp_block,nblocks)+1
        IF(jatom_block<iatom_block) CYCLE

        distribution_counter_end = distribution_counter_end + 1

        IF(binned_dist(ibin)%istart == -1_int_8 ) binned_dist(ibin)%istart = atom_block
        IF(binned_dist(ibin)%jstart == 0_int_8 ) binned_dist(ibin)%jstart = jatom_block
        IF(binned_dist(ibin)%kstart == 0_int_8 ) binned_dist(ibin)%kstart = katom_block
        IF(binned_dist(ibin)%lstart == 0_int_8 ) binned_dist(ibin)%lstart = latom_block

        CALL build_pair_list(list_ij,set_list_ij,(iatom_block-1)*block_size +1,MIN(iatom_block*block_size,natom), &
                             (jatom_block-1)*block_size +1, MIN(jatom_block*block_size,natom),&
                             kind_of,basis_parameter,particle_set, &
                             do_periodic,coeffs_set,coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell)

        CALL build_pair_list(list_kl,set_list_kl,(katom_block-1)*block_size +1,MIN(katom_block*block_size,natom), &
                             (latom_block-1)*block_size +1, MIN(latom_block*block_size,natom),&
                             kind_of,basis_parameter,particle_set, &
                             do_periodic,coeffs_set,coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell)

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

          nsgfa => basis_parameter(ikind)%nsgf
          nsgfb => basis_parameter(jkind)%nsgf

          DO i_list_kl=1,list_kl%n_element

            katom=list_kl%elements(i_list_kl)%pair(1)
            latom=list_kl%elements(i_list_kl)%pair(2)

            IF( .NOT. (katom+latom<=iatom+jatom))  CYCLE
            IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE

            i_set_list_kl_start=list_kl%elements(i_list_kl)%set_bounds(1)
            i_set_list_kl_stop=list_kl%elements(i_list_kl)%set_bounds(2)
            kkind=list_kl%elements(i_list_kl)%kind_pair(1)
            lkind=list_kl%elements(i_list_kl)%kind_pair(2)
            rc=list_kl%elements(i_list_kl)%r1
            rd=list_kl%elements(i_list_kl)%r2
            rcd2=list_kl%elements(i_list_kl)%dist2

            nsgfc => basis_parameter(kkind)%nsgf
            nsgfd => basis_parameter(lkind)%nsgf

            IF( (coeffs_kind(jkind,ikind)%x(1)*rab2+&
                 coeffs_kind(jkind,ikind)%x(2) ) + &
                (coeffs_kind(lkind,kkind)%x(1)*rcd2+&
                 coeffs_kind(lkind,kkind)%x(2)) < log10_eps_schwarz ) CYCLE

            IF(.NOT. (is_assoc_atomic_block_global(latom,iatom)  >=1 .AND. &
                      is_assoc_atomic_block_global(katom,iatom)  >=1 .AND. &
                      is_assoc_atomic_block_global(katom,jatom)  >=1 .AND. &
                      is_assoc_atomic_block_global(latom,jatom)  >=1 ) ) CYCLE

            DO i_set_list_ij=i_set_list_ij_start, i_set_list_ij_stop
              iset=set_list_ij(i_set_list_ij)%pair(1)
              jset=set_list_ij(i_set_list_ij)%pair(2)

              max_val1 = coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + &
                         coeffs_set(jset,iset,jkind,ikind)%x(2)
              IF( max_val1 +( coeffs_kind(lkind,kkind)%x(1)*rcd2 + &
                              coeffs_kind(lkind,kkind)%x(2) ) < log10_eps_schwarz) CYCLE

              DO i_set_list_kl=i_set_list_kl_start, i_set_list_kl_stop
                kset=set_list_kl(i_set_list_kl)%pair(1)
                lset=set_list_kl(i_set_list_kl)%pair(2)

                max_val2 = max_val1 + (coeffs_set(lset,kset,lkind,kkind)%x(1)*rcd2 + &
                                       coeffs_set(lset,kset,lkind,kkind)%x(2) )
                IF(max_val2<log10_eps_schwarz) CYCLE
                IF( do_p_screening) THEN
                  pmax_entry = MAX(pmax_set(kset,iset,katom,iatom), &
                                   pmax_set(lset,jset,latom,jatom), &
                                   pmax_set(lset,iset,latom,iatom), &
                                   pmax_set(kset,jset,katom,jatom))
                ELSE
                  pmax_entry = 0.0_dp
                END IF
                max_val2 = max_val2 + pmax_entry
                IF(max_val2<log10_eps_schwarz) CYCLE
                current_counter = current_counter + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
              ENDDO ! i_set_list_kl
            ENDDO ! i_set_list_ij
          ENDDO ! i_list_kl
        ENDDO ! i_list_ij

        IF( current_counter >= int_per_cpu) THEN
          binned_dist(ibin)%number_of_atom_quartets = distribution_counter_end - distribution_counter_start + 1
          binned_dist(ibin)%cost = binned_dist(ibin)%cost + current_counter
          ! protect against int_per_cpu==0 and possibly roundoff
          ibin = MIN(ibin + 1,nbins)
          distribution_counter_start = distribution_counter_end + 1
          restore_counter = current_counter
          current_counter = 0
        END IF
      ENDDO

!$OMP BARRIER
!$OMP MASTER
    CALL timestop(handle_inner)
    CALL timeset(routineN//"_dist",handle_inner)
!$OMP END MASTER
!$OMP BARRIER
      !! Fill the last bin if necessary
      IF( current_counter /= 0 ) THEN  
        binned_dist(ibin)%cost = binned_dist(ibin)%cost + current_counter
        binned_dist(ibin)%number_of_atom_quartets = distribution_counter_end-distribution_counter_start+1
      END IF

      !! store the original cpu_id
      binned_dist(:)%original_cpu_id = my_process_id + 1 !para_env%mepos + 1
!$OMP BARRIER     
      !! store all local results in a big cost matrix
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
      ALLOCATE(distribution_vector(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP BARRIER
      ALLOCATE(local_cost_matrix(SIZE(cost_matrix,1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      local_cost_matrix = cost_matrix
      CALL optimize_distribution(para_env, ncpu*nbins*n_threads, ncpu*n_threads, local_cost_matrix, &
                                 distribution_vector, error)

!$OMP BARRIER
!$OMP MASTER
      CALL timestop(handle_inner)
      CALL timeset(routineN//"_redist",handle_inner)
      ALLOCATE(shm_distribution_vector(ncpu*nbins*n_threads), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      shm_distribution_vector = distribution_vector
      !! Collect local data to global array
      ALLOCATE(full_dist(ncpu*n_threads,nbins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      full_dist(:,:)%istart = 0_int_8
      full_dist(:,:)%jstart = 0_int_8
      full_dist(:,:)%kstart = 0_int_8
      full_dist(:,:)%lstart = 0_int_8
      full_dist(:,:)%number_of_atom_quartets = 0_int_8
      full_dist(:,:)%cost = 0_int_8
      full_dist(:,:)%original_cpu_id = 0_int_8
      full_dist(:,:)%time_first_scf = 0.0_dp
      full_dist(:,:)%time_other_scf = 0.0_dp
      full_dist(:,:)%time_forces = 0.0_dp
!$OMP END MASTER
!$OMP BARRIER
      mepos = para_env%mepos + 1
      full_dist((mepos-1)*n_threads+i_thread+1,:) = binned_dist(:)
     
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
          sendbuffer((j-1)*7*nbins+(i-1)*7+5) = full_dist(mepos*n_threads+j,i)%number_of_atom_quartets
          sendbuffer((j-1)*7*nbins+(i-1)*7+6) = full_dist(mepos*n_threads+j,i)%cost
          sendbuffer((j-1)*7*nbins+(i-1)*7+7) = full_dist(mepos*n_threads+j,i)%original_cpu_id
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
            full_dist(data_from*n_threads+j,i)%istart = sendbuffer((j-1)*7*nbins+(i-1)*7+1)
            full_dist(data_from*n_threads+j,i)%jstart = sendbuffer((j-1)*7*nbins+(i-1)*7+2)
            full_dist(data_from*n_threads+j,i)%kstart = sendbuffer((j-1)*7*nbins+(i-1)*7+3)
            full_dist(data_from*n_threads+j,i)%lstart = sendbuffer((j-1)*7*nbins+(i-1)*7+4)
            full_dist(data_from*n_threads+j,i)%number_of_atom_quartets = sendbuffer((j-1)*7*nbins+(i-1)*7+5)
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
      tmp_dist(:)%jstart = 0_int_8
      tmp_dist(:)%kstart = 0_int_8
      tmp_dist(:)%lstart = 0_int_8
      tmp_dist(:)%number_of_atom_quartets = 0_int_8
      tmp_dist(:)%cost = 0_int_8
      tmp_dist(:)%original_cpu_id = 0_int_8
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
      CALL hfx_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)

!$OMP BARRIER
!$OMP MASTER
      DEALLOCATE(full_dist,cost_matrix,shm_distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER
      DEALLOCATE(tmp_dist, tmp_pos, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(distribution_vector, STAT=stat)
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
  END SUBROUTINE load_balance_energy


! *****************************************************************************
!> \brief Distributes the computation of derivative-eri's to all available processes.
!> \param x_data Object that stores the indices array
!> \param eps_schwarz screening parameter
!> \param particle_set , atomic_kind_set, para_env
!> \param max_set Maximum number of set to be considered
!> \param para_env para_env
!> \param potential_parameter short/longrange screening information 
!> \param coeffs_set, coeffs_kind: screening functions
!> \param pmax_set Initial screening matrix
!> \param is_assoc_atomic_kind KS-matrix sparsity
!> \param do_periodic flag for periodicity
!> \param load_balance_parameter Paramters for Monte-Carlo routines
!> \param kind_of helper array for mapping
!> \param basis_parameter Basis set parameters
!> \param i_thread Process ID of current Thread
!> \param n_thread Total Number of Threads
!> \param cell cell
!> \param screen_pmat_forces Flag for p screening
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
!>      02.2009 completely refactored [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE load_balance_forces(x_data, eps_schwarz, particle_set,max_set,para_env,&
                                 potential_parameter, coeffs_set, coeffs_kind,&
                                 pmax_set, is_assoc_atomic_block_global, do_periodic, &
                                 load_balance_parameter, kind_of, basis_parameter,&
                                 i_thread, n_threads, cell, screen_pmat_forces, error)
    TYPE(hfx_type), POINTER                  :: x_data
    REAL(dp), INTENT(IN)                     :: eps_schwarz
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: max_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(hfx_potential_type), INTENT(IN)     :: potential_parameter
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: coeffs_kind
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: pmax_set
    INTEGER, DIMENSION(:, :)                 :: is_assoc_atomic_block_global
    LOGICAL                                  :: do_periodic
    TYPE(hfx_load_balance_type)              :: load_balance_parameter
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    INTEGER, INTENT(IN)                      :: i_thread, n_threads
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(IN)                      :: screen_pmat_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'load_balance_forces', &
      routineP = moduleN//':'//routineN

    INTEGER :: block_size, data_from, dest, i, i_list_ij, i_list_kl, &
      i_set_list_ij, i_set_list_ij_start, i_set_list_ij_stop, i_set_list_kl, &
      i_set_list_kl_start, i_set_list_kl_stop, iatom, iatom_block, ibin, &
      icpu, ikind, ipos, iset, j, jatom, jatom_block, jkind, jset, katom, &
      katom_block, kkind, kset, latom, latom_block, lkind, lset, &
      my_process_id, natom, nbins, req(2), source, stat
    INTEGER(int_8) :: atom_block, current_counter, distribution_counter_end, &
      distribution_counter_start, int_per_cpu, mepos, n_processes, nblocks, &
      ncpu, restore_counter, tmp_block, total_integrals
    INTEGER(int_8), DIMENSION(:), POINTER    :: local_cost_matrix, recbuffer, &
                                                sendbuffer, swapbuffer
    INTEGER(int_8), DIMENSION(:), POINTER, &
      SAVE                                   :: cost_matrix
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_pos
    INTEGER, DIMENSION(:), POINTER           :: distribution_vector, nsgfa, &
                                                nsgfb, nsgfc, nsgfd
    INTEGER, DIMENSION(:), POINTER, SAVE     :: shm_distribution_vector
    LOGICAL                                  :: failure
    REAL(dp) :: coeffs_kind_max0, log10_eps_schwarz, log10_pmax, log_2, &
      max_val1, max_val2, omega2Inv, ra(3), rab2, rb(3), rc(3), rcd2, rd(3)
    TYPE(hfx_distribution), DIMENSION(:), &
      POINTER                                :: binned_dist, ptr_to_tmp_dist, &
                                                tmp_dist
    TYPE(hfx_distribution), &
      DIMENSION(:, :), POINTER, SAVE         :: full_dist
    TYPE(pair_list_type)                     :: list_ij, list_kl
    TYPE(pair_set_list_type), ALLOCATABLE, &
      DIMENSION(:)                           :: set_list_ij, set_list_kl

    log_2 = LOG10(2.0_dp)
     
    IF(potential_parameter%omega /= 0.0_dp) THEN
      omega2Inv = 1.0_dp/(potential_parameter%omega**2)
    ELSE
      omega2Inv = 0.0_dp
    END IF

    natom=SIZE(particle_set)
    log10_eps_schwarz = LOG10(eps_schwarz)
    failure=.FALSE.
    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads
    coeffs_kind_max0=MAXVAL(coeffs_kind(:,:)%x(2))
    ALLOCATE(set_list_ij((max_set*load_balance_parameter%block_size)**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(set_list_kl((max_set*load_balance_parameter%block_size)**2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF(n_processes == 1) THEN
      ALLOCATE(tmp_dist(1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1)%number_of_atom_quartets = HUGE(tmp_dist(1)%number_of_atom_quartets)
      tmp_dist(1)%istart = 0_int_8
      tmp_dist(1)%jstart = 1_int_8
      tmp_dist(1)%kstart = 1_int_8
      tmp_dist(1)%lstart = 1_int_8
      tmp_dist(1)%original_cpu_id = 1_int_8
      ptr_to_tmp_dist => tmp_dist(:)
      CALL hfx_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)
      DEALLOCATE(tmp_dist,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      total_integrals = 0
      my_process_id = para_env%mepos*n_threads + i_thread
      block_size = load_balance_parameter%block_size
      nblocks = MAX((natom+block_size-1)/block_size,1)
      DO atom_block = my_process_id,nblocks**4-1,n_processes
        latom_block = MODULO(atom_block,nblocks)+1
        tmp_block = atom_block/nblocks
        katom_block = MODULO(tmp_block,nblocks)+1
        IF(latom_block<katom_block) CYCLE
        tmp_block = tmp_block/nblocks
        jatom_block = MODULO(tmp_block,nblocks)+1
        tmp_block = tmp_block/nblocks
        iatom_block = MODULO(tmp_block,nblocks)+1
        IF(jatom_block<iatom_block) CYCLE

        CALL build_pair_list(list_ij,set_list_ij,(iatom_block-1)*block_size +1,MIN(iatom_block*block_size,natom), &
                             (jatom_block-1)*block_size +1, MIN(jatom_block*block_size,natom),&
                             kind_of,basis_parameter,particle_set, &
                             do_periodic,coeffs_set,coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell)

        CALL build_pair_list(list_kl,set_list_kl,(katom_block-1)*block_size +1,MIN(katom_block*block_size,natom), &
                             (latom_block-1)*block_size +1, MIN(latom_block*block_size,natom),&
                             kind_of,basis_parameter,particle_set, &
                             do_periodic,coeffs_set,coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell)
        
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

          nsgfa => basis_parameter(ikind)%nsgf
          nsgfb => basis_parameter(jkind)%nsgf

          DO i_list_kl=1,list_kl%n_element

            katom=list_kl%elements(i_list_kl)%pair(1)
            latom=list_kl%elements(i_list_kl)%pair(2)

            IF( .NOT. (katom+latom<=iatom+jatom))  CYCLE
            IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
            IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE

            i_set_list_kl_start=list_kl%elements(i_list_kl)%set_bounds(1)
            i_set_list_kl_stop=list_kl%elements(i_list_kl)%set_bounds(2)
            kkind=list_kl%elements(i_list_kl)%kind_pair(1)
            lkind=list_kl%elements(i_list_kl)%kind_pair(2)
            rc=list_kl%elements(i_list_kl)%r1
            rd=list_kl%elements(i_list_kl)%r2
            rcd2=list_kl%elements(i_list_kl)%dist2

            nsgfc => basis_parameter(kkind)%nsgf
            nsgfd => basis_parameter(lkind)%nsgf

            IF( (coeffs_kind(jkind,ikind)%x(1)*rab2+&
                 coeffs_kind(jkind,ikind)%x(2) ) + &
                 (coeffs_kind(lkind,kkind)%x(1)*rcd2+&
                 coeffs_kind(lkind,kkind)%x(2)) < log10_eps_schwarz ) CYCLE

            IF(.NOT. (is_assoc_atomic_block_global(latom,iatom)  >=1 .AND. &
                      is_assoc_atomic_block_global(katom,iatom)  >=1 .AND. &
                      is_assoc_atomic_block_global(katom,jatom)  >=1 .AND. &
                      is_assoc_atomic_block_global(latom,jatom)  >=1 ) ) CYCLE

            DO i_set_list_ij=i_set_list_ij_start, i_set_list_ij_stop
              iset=set_list_ij(i_set_list_ij)%pair(1)
              jset=set_list_ij(i_set_list_ij)%pair(2)

              max_val1 = coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + &
                         coeffs_set(jset,iset,jkind,ikind)%x(2)
              IF( max_val1 + coeffs_kind(lkind,kkind)%x(1)*rcd2 + &
                             coeffs_kind(lkind,kkind)%x(2)  <log10_eps_schwarz) CYCLE

              DO i_set_list_kl=i_set_list_kl_start, i_set_list_kl_stop
                kset=set_list_kl(i_set_list_kl)%pair(1)
                lset=set_list_kl(i_set_list_kl)%pair(2)

                max_val2 = max_val1 + (coeffs_set(lset,kset,lkind,kkind)%x(1)*rcd2 + &
                                       coeffs_set(lset,kset,lkind,kkind)%x(2) )

                IF(max_val2 <log10_eps_schwarz ) CYCLE
                IF( screen_pmat_forces) THEN
                  log10_pmax = log_2 + MAX(pmax_set(kset,iset,katom,iatom) + pmax_set(lset,jset,latom,jatom),&
                                           pmax_set(lset,iset,latom,iatom) + pmax_set(kset,jset,katom,jatom))
                ELSE
                  log10_pmax = 0.0_dp
                END IF
                max_val2= max_val2 + log10_pmax
                IF( max_val2< log10_eps_schwarz) CYCLE
                total_integrals = total_integrals + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
              ENDDO ! i_set_list_kl
            ENDDO ! i_set_list_ij
          ENDDO ! i_list_kl
        ENDDO ! i_list_ij
      ENDDO ! atom_block

      nbins = load_balance_parameter%nbins 
      int_per_cpu = (total_integrals + nbins -1 ) / nbins
      ALLOCATE(binned_dist(nbins), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      binned_dist(:)%istart = -1_int_8
      binned_dist(:)%jstart = 0_int_8
      binned_dist(:)%kstart = 0_int_8
      binned_dist(:)%lstart = 0_int_8
      binned_dist(:)%number_of_atom_quartets = 0_int_8
      binned_dist(:)%cost = 0_int_8
      binned_dist(:)%original_cpu_id = 0_int_8
 
      current_counter = 0
      restore_counter = 0
      ipos = 1
      mepos = 1
      distribution_counter_start = 1
      distribution_counter_end = 0
      ibin = 1
      DO atom_block = my_process_id,nblocks**4-1,n_processes
        latom_block = MODULO(atom_block,nblocks)+1
        tmp_block = atom_block/nblocks
        katom_block = MODULO(tmp_block,nblocks)+1
        IF(latom_block<katom_block) CYCLE
        tmp_block = tmp_block/nblocks
        jatom_block = MODULO(tmp_block,nblocks)+1
        tmp_block = tmp_block/nblocks
        iatom_block = MODULO(tmp_block,nblocks)+1
        IF(jatom_block<iatom_block) CYCLE

        distribution_counter_end = distribution_counter_end + 1

        IF(binned_dist(ibin)%istart == -1_int_8 ) binned_dist(ibin)%istart = atom_block
        IF(binned_dist(ibin)%jstart == 0_int_8 ) binned_dist(ibin)%jstart = jatom_block
        IF(binned_dist(ibin)%kstart == 0_int_8 ) binned_dist(ibin)%kstart = katom_block
        IF(binned_dist(ibin)%lstart == 0_int_8 ) binned_dist(ibin)%lstart = latom_block

        CALL build_pair_list(list_ij,set_list_ij,(iatom_block-1)*block_size +1,MIN(iatom_block*block_size,natom), &
                             (jatom_block-1)*block_size +1, MIN(jatom_block*block_size,natom),&
                             kind_of,basis_parameter,particle_set, &
                             do_periodic,coeffs_set,coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell)

        CALL build_pair_list(list_kl,set_list_kl,(katom_block-1)*block_size +1,MIN(katom_block*block_size,natom), &
                             (latom_block-1)*block_size +1, MIN(latom_block*block_size,natom),&
                             kind_of,basis_parameter,particle_set, &
                             do_periodic,coeffs_set,coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell)

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

          nsgfa => basis_parameter(ikind)%nsgf
          nsgfb => basis_parameter(jkind)%nsgf

          DO i_list_kl=1,list_kl%n_element
            katom=list_kl%elements(i_list_kl)%pair(1)
            latom=list_kl%elements(i_list_kl)%pair(2)

            IF( .NOT. (katom+latom<=iatom+jatom))  CYCLE
            IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
            IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE

            i_set_list_kl_start=list_kl%elements(i_list_kl)%set_bounds(1)
            i_set_list_kl_stop=list_kl%elements(i_list_kl)%set_bounds(2)
            kkind=list_kl%elements(i_list_kl)%kind_pair(1)
            lkind=list_kl%elements(i_list_kl)%kind_pair(2)
            rc=list_kl%elements(i_list_kl)%r1
            rd=list_kl%elements(i_list_kl)%r2
            rcd2=list_kl%elements(i_list_kl)%dist2

            nsgfc => basis_parameter(kkind)%nsgf
            nsgfd => basis_parameter(lkind)%nsgf

            IF( (coeffs_kind(jkind,ikind)%x(1)*rab2+&
                 coeffs_kind(jkind,ikind)%x(2) ) + &
                (coeffs_kind(lkind,kkind)%x(1)*rcd2+&
                 coeffs_kind(lkind,kkind)%x(2)) < log10_eps_schwarz ) CYCLE
           
            IF(.NOT. (is_assoc_atomic_block_global(latom,iatom)  >=1 .AND. &
                      is_assoc_atomic_block_global(katom,iatom)  >=1 .AND. &
                      is_assoc_atomic_block_global(katom,jatom)  >=1 .AND. &
                      is_assoc_atomic_block_global(latom,jatom)  >=1 ) ) CYCLE

            DO i_set_list_ij=i_set_list_ij_start, i_set_list_ij_stop
              iset=set_list_ij(i_set_list_ij)%pair(1)
              jset=set_list_ij(i_set_list_ij)%pair(2)

              max_val1 = coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + &
                         coeffs_set(jset,iset,jkind,ikind)%x(2)
              IF( max_val1 +( coeffs_kind(lkind,kkind)%x(1)*rcd2 + &
                              coeffs_kind(lkind,kkind)%x(2) ) < log10_eps_schwarz) CYCLE

              DO i_set_list_kl=i_set_list_kl_start, i_set_list_kl_stop
                kset=set_list_kl(i_set_list_kl)%pair(1)
                lset=set_list_kl(i_set_list_kl)%pair(2)

                max_val2 = max_val1 + (coeffs_set(lset,kset,lkind,kkind)%x(1)*rcd2 + &
                                       coeffs_set(lset,kset,lkind,kkind)%x(2) )
                IF(max_val2<log10_eps_schwarz) CYCLE
                IF( screen_pmat_forces) THEN
                  log10_pmax = log_2 + MAX(pmax_set(kset,iset,katom,iatom) + pmax_set(lset,jset,latom,jatom),&
                                          pmax_set(lset,iset,latom,iatom) + pmax_set(kset,jset,katom,jatom))
                ELSE
                  log10_pmax = 0.0_dp
                END IF
                max_val2 = max_val2 + log10_pmax
                IF(max_val2<log10_eps_schwarz) CYCLE
                current_counter = current_counter + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IF( current_counter >= int_per_cpu) THEN
          binned_dist(ibin)%number_of_atom_quartets = distribution_counter_end - distribution_counter_start + 1
          binned_dist(ibin)%cost = binned_dist(ibin)%cost + current_counter
          ! protect against int_per_cpu==0 and possibly roundoff
          ibin = MIN(ibin + 1,nbins)
          distribution_counter_start = distribution_counter_end + 1
          restore_counter = current_counter
          current_counter = 0
        END IF
      ENDDO

      !! Fill the last bin if necessary
      IF( current_counter /= 0 ) THEN  
        binned_dist(ibin)%cost = binned_dist(ibin)%cost + current_counter
        binned_dist(ibin)%number_of_atom_quartets = distribution_counter_end-distribution_counter_start+1
      END IF

      !! store the original cpu_di
      binned_dist(:)%original_cpu_id = my_process_id + 1 !para_env%mepos + 1
    
      !! store all local results in a big cost matrix
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
      ALLOCATE(sendbuffer(nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(recbuffer(nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      sendbuffer = cost_matrix(mepos*nbins*n_threads+1:mepos*nbins*n_threads+nbins*n_threads)

      ! sync before/after ring of isendrecv
      CALL mp_sync(para_env%group)

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
      ! sync before/after ring of isendrecv
      CALL mp_sync(para_env%group)
!$OMP END MASTER

!$OMP BARRIER
      !! Find an optimal distribution i.e. assign each element of the cost matrix to a certain process
      ALLOCATE(distribution_vector(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

!$OMP BARRIER
      ALLOCATE(local_cost_matrix(SIZE(cost_matrix,1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      local_cost_matrix = cost_matrix
      CALL optimize_distribution(para_env, ncpu*nbins*n_threads, ncpu*n_threads, local_cost_matrix, &
                                 distribution_vector, error)
!$OMP BARRIER
!$OMP MASTER
      ALLOCATE(shm_distribution_vector(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      shm_distribution_vector = distribution_vector
      !! Collect local data to global array
      ALLOCATE(full_dist(ncpu*n_threads,nbins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      full_dist(:,:)%istart = 0_int_8
      full_dist(:,:)%jstart = 0_int_8
      full_dist(:,:)%kstart = 0_int_8
      full_dist(:,:)%lstart = 0_int_8
      full_dist(:,:)%number_of_atom_quartets = 0_int_8
      full_dist(:,:)%cost = 0_int_8
      full_dist(:,:)%original_cpu_id = 0_int_8
      full_dist(:,:)%time_first_scf = 0.0_dp
      full_dist(:,:)%time_other_scf = 0.0_dp
      full_dist(:,:)%time_forces = 0.0_dp
!$OMP END MASTER
!$OMP BARRIER
      mepos=para_env%mepos+1
      full_dist((mepos-1)*n_threads+i_thread+1,:) = binned_dist(:)

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
          sendbuffer((j-1)*7*nbins+(i-1)*7+5) = full_dist(mepos*n_threads+j,i)%number_of_atom_quartets
          sendbuffer((j-1)*7*nbins+(i-1)*7+6) = full_dist(mepos*n_threads+j,i)%cost
          sendbuffer((j-1)*7*nbins+(i-1)*7+7) = full_dist(mepos*n_threads+j,i)%original_cpu_id
        END DO
      END DO

      ! sync before/after a ring of isendrecv
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
            full_dist(data_from*n_threads+j,i)%istart = sendbuffer((j-1)*7*nbins+(i-1)*7+1)
            full_dist(data_from*n_threads+j,i)%jstart = sendbuffer((j-1)*7*nbins+(i-1)*7+2)
            full_dist(data_from*n_threads+j,i)%kstart = sendbuffer((j-1)*7*nbins+(i-1)*7+3)
            full_dist(data_from*n_threads+j,i)%lstart = sendbuffer((j-1)*7*nbins+(i-1)*7+4)
            full_dist(data_from*n_threads+j,i)%number_of_atom_quartets = sendbuffer((j-1)*7*nbins+(i-1)*7+5)
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
      ! sync before/after ring of isendrecv
      ! this one is a magic sync that can speedup the above message passing by a factor of 100 on an XT5
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
      tmp_dist(:)%jstart = 0_int_8
      tmp_dist(:)%kstart = 0_int_8
      tmp_dist(:)%lstart = 0_int_8
      tmp_dist(:)%number_of_atom_quartets = 0_int_8
      tmp_dist(:)%cost = 0_int_8
      tmp_dist(:)%original_cpu_id = 0_int_8
      tmp_dist(:)%time_first_scf = 0.0_dp
      tmp_dist(:)%time_other_scf = 0.0_dp
      tmp_dist(:)%time_forces = 0.0_dp
      mepos = my_process_id + 1
      DO icpu= 1,n_processes
        DO i=1,nbins
          IF(shm_distribution_vector((icpu-1)*nbins + i) == mepos) THEN
              tmp_dist(tmp_pos(mepos)) = full_dist(icpu,i)
              tmp_pos(mepos) = tmp_pos(mepos) + 1
          END IF
        END DO
      END DO

      !! Assign the load to each process
      NULLIFY(ptr_to_tmp_dist)
      mepos = my_process_id + 1
      ptr_to_tmp_dist => tmp_dist(1:tmp_pos(mepos)-1)

      CALL hfx_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)
!$OMP BARRIER
!$OMP MASTER
      DEALLOCATE(full_dist,cost_matrix,shm_distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER
      DEALLOCATE(tmp_dist, tmp_pos, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(binned_dist, local_cost_matrix, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(set_list_ij, set_list_kl, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE load_balance_forces

! *****************************************************************************
!> \brief Cheap way of redistributing the eri's
!> \param x_data Object that stores the indices array
!> \param para_env para_env
!> \param load_balance_parameter contains parmameter for Monte-Carlo routines
!> \param i_thread current thread ID
!> \param n_threads Total Number of threads
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      The cost matrix is given by the walltime for each bin that is measured 
!>      during the calculation
!> \par History
!>      12.2007 created [Manuel Guidon]
!>      02.2009 optimize Memory Usage [Manuel Guidon]
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
    INTEGER(int_8)                           :: mepos, n_processes, ncpu
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
    INTEGER, SAVE                            :: max_bin_size
    LOGICAL                                  :: failure
    TYPE(hfx_distribution), DIMENSION(:), &
      POINTER                                :: binned_dist, ptr_to_tmp_dist, &
                                                tmp_dist
    TYPE(hfx_distribution), &
      DIMENSION(:, :), POINTER, SAVE         :: full_dist

    failure=.FALSE.

    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads
    !! If there is only 1 cpu skip the binning
    IF(n_processes == 1) THEN
      ALLOCATE(tmp_dist(1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1)%number_of_atom_quartets = HUGE(tmp_dist(1)%number_of_atom_quartets)
      tmp_dist(1)%istart = 0_int_8
      tmp_dist(1)%jstart = 1_int_8
      tmp_dist(1)%kstart = 1_int_8
      tmp_dist(1)%lstart = 1_int_8
      tmp_dist(1)%original_cpu_id = 1_int_8
      ptr_to_tmp_dist => tmp_dist(:)
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
        IF(binned_dist(ibin)%number_of_atom_quartets == 0 ) THEN
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
      !! Find an optimal distribution i.e. assign each element of the cost matrix to a certain process
      ALLOCATE(distribution_vector(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP BARRIER
      ALLOCATE(local_cost_matrix(SIZE(cost_matrix,1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      local_cost_matrix = cost_matrix
      CALL optimize_distribution(para_env, ncpu*nbins*n_threads, ncpu*n_threads, local_cost_matrix, &
                                 distribution_vector, error)

!$OMP BARRIER
!$OMP MASTER
      ALLOCATE(shm_distribution_vector(ncpu*nbins*n_threads), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      shm_distribution_vector = distribution_vector
      ALLOCATE(full_dist(ncpu*n_threads,max_bin_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      full_dist(:,:)%istart = 0_int_8
      full_dist(:,:)%jstart = 0_int_8
      full_dist(:,:)%kstart = 0_int_8
      full_dist(:,:)%lstart = 0_int_8
      full_dist(:,:)%number_of_atom_quartets = 0_int_8
      full_dist(:,:)%cost = 0_int_8
      full_dist(:,:)%original_cpu_id = 0_int_8
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
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+5) = full_dist(mepos*n_threads+j,i)%number_of_atom_quartets
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+6) = full_dist(mepos*n_threads+j,i)%cost
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+7) = full_dist(mepos*n_threads+j,i)%original_cpu_id
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
            full_dist(data_from*n_threads+j,i)%istart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+1)
            full_dist(data_from*n_threads+j,i)%jstart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+2)
            full_dist(data_from*n_threads+j,i)%kstart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+3)
            full_dist(data_from*n_threads+j,i)%lstart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+4)
            full_dist(data_from*n_threads+j,i)%number_of_atom_quartets = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+5)
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
      ! sync before/after ring of isendrecv
      CALL mp_sync(para_env%group)
      DEALLOCATE(recbuffer, sendbuffer, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER
      !! reorder the distribution according to the distribution vector
      ALLOCATE(tmp_pos(ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_pos = 1
      ALLOCATE(tmp_dist(nbins*ncpu*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      tmp_dist(:)%istart = 0_int_8
      tmp_dist(:)%jstart = 0_int_8
      tmp_dist(:)%kstart = 0_int_8
      tmp_dist(:)%lstart = 0_int_8
      tmp_dist(:)%number_of_atom_quartets = 0_int_8
      tmp_dist(:)%cost = 0_int_8
      tmp_dist(:)%original_cpu_id = 0_int_8
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

      CALL hfx_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)

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
      DEALLOCATE(distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(binned_dist, local_cost_matrix, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE update_load_balance_energy

! *****************************************************************************
!> \brief Cheap way of redistributing the derivative eri's
!> \param x_data Object that stores the indices array
!> \param para_env para_env
!> \param load_balance_parameter contains parmameter for Monte-Carlo routines
!> \param i_thread current thread ID
!> \param n_threads Total Number of threads
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      The cost matrix is given by the walltime for each bin that is measured 
!>      during the calculation
!> \par History
!>      12.2007 created [Manuel Guidon]
!>      02.2009 optimize Memory Usage [Manuel Guidon]
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
    INTEGER(int_8)                           :: mepos, n_processes, ncpu
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
    INTEGER, SAVE                            :: max_bin_size
    LOGICAL                                  :: failure
    TYPE(hfx_distribution), DIMENSION(:), &
      POINTER                                :: binned_dist, ptr_to_tmp_dist, &
                                                tmp_dist
    TYPE(hfx_distribution), &
      DIMENSION(:, :), POINTER, SAVE         :: full_dist

    failure=.FALSE.

    ncpu = para_env%num_pe
    n_processes = ncpu * n_threads
    !! If there is only 1 cpu skip the binning
    IF(n_processes == 1) THEN
      ALLOCATE(tmp_dist(1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1)%number_of_atom_quartets = HUGE(tmp_dist(1)%number_of_atom_quartets)
      tmp_dist(1)%istart = 0_int_8
      tmp_dist(1)%jstart = 1_int_8
      tmp_dist(1)%kstart = 1_int_8
      tmp_dist(1)%lstart = 1_int_8
      tmp_dist(1)%original_cpu_id = 1_int_8
      ptr_to_tmp_dist => tmp_dist(:)
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
        IF(binned_dist(ibin)%number_of_atom_quartets == 0 ) THEN
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
      icpu = para_env%mepos+1
      DO i=1,my_bin_size
        cost_matrix(my_global_start_idx + i) = binned_dist(i)%cost
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
      ! sync before/after ring of isendrecv
      CALL mp_sync(para_env%group)
      DEALLOCATE(recbuffer, sendbuffer, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END MASTER
!$OMP BARRIER
      !! Find an optimal distribution i.e. assign each element of the cost matrix to a certain process
      ALLOCATE(distribution_vector(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP BARRIER
      ALLOCATE(local_cost_matrix(SIZE(cost_matrix,1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      local_cost_matrix = cost_matrix
      CALL optimize_distribution(para_env, ncpu*nbins*n_threads, ncpu*n_threads, local_cost_matrix, &
                                 distribution_vector, error)
      
!$OMP BARRIER
!$OMP MASTER
      ALLOCATE(shm_distribution_vector(ncpu*nbins*n_threads),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      shm_distribution_vector = distribution_vector
      !! Collect local data to global array
      ALLOCATE(full_dist(ncpu*n_threads,max_bin_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      full_dist(:,:)%istart = 0_int_8
      full_dist(:,:)%jstart = 0_int_8
      full_dist(:,:)%kstart = 0_int_8
      full_dist(:,:)%lstart = 0_int_8
      full_dist(:,:)%number_of_atom_quartets = 0_int_8
      full_dist(:,:)%cost = 0_int_8
      full_dist(:,:)%original_cpu_id = 0_int_8
      full_dist(:,:)%time_first_scf = 0.0_dp
      full_dist(:,:)%time_other_scf = 0.0_dp
      full_dist(:,:)%time_forces = 0.0_dp
!$OMP END MASTER

!$OMP BARRIER
      mepos = para_env%mepos + 1
      full_dist((mepos-1)*n_threads+i_thread+1,1:my_bin_size) = binned_dist(1:my_bin_size)
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
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+5) = full_dist(mepos*n_threads+j,i)%number_of_atom_quartets
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+6) = full_dist(mepos*n_threads+j,i)%cost
          sendbuffer((j-1)*7*max_bin_size+(i-1)*7+7) = full_dist(mepos*n_threads+j,i)%original_cpu_id
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
            full_dist(data_from*n_threads+j,i)%istart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+1)
            full_dist(data_from*n_threads+j,i)%jstart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+2)
            full_dist(data_from*n_threads+j,i)%kstart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+3)
            full_dist(data_from*n_threads+j,i)%lstart = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+4)
            full_dist(data_from*n_threads+j,i)%number_of_atom_quartets = sendbuffer((j-1)*7*max_bin_size+(i-1)*7+5)
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
      tmp_dist(:)%jstart = 0_int_8
      tmp_dist(:)%kstart = 0_int_8
      tmp_dist(:)%lstart = 0_int_8
      tmp_dist(:)%number_of_atom_quartets = 0_int_8
      tmp_dist(:)%cost = 0_int_8
      tmp_dist(:)%original_cpu_id = 0_int_8
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

      CALL hfx_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)

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
      DEALLOCATE(distribution_vector, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(binned_dist, local_cost_matrix, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE update_load_balance_forces

! *****************************************************************************
!> \brief Minimizes the maximum cost per cpu by shuffling around all bins
!> \param para_env para_env
!> \param total_number_of_bins 
!> \param number_of_processes
!> \param bin_costs costs per bin
!> \param distribution_vector will contain the final distribution
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2009 created from a hack by Joost [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE optimize_distribution(para_env, total_number_of_bins, number_of_processes, bin_costs, &
                                   distribution_vector, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER(int_8)                           :: total_number_of_bins, &
                                                number_of_processes
    INTEGER(int_8), DIMENSION(:), POINTER    :: bin_costs
    INTEGER, DIMENSION(:), POINTER           :: distribution_vector
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'optimize_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, itmp, j, nstep, stat
    INTEGER(int_8), DIMENSION(:), POINTER    :: my_cost_cpu, tmp_cost, &
                                                tmp_cpu_cost
    INTEGER, DIMENSION(:), POINTER           :: tmp_cpu_index, tmp_index
    LOGICAL                                  :: failure

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

    CALL sort(tmp_cost,INT(total_number_of_bins),tmp_index)
    my_cost_cpu=0
!
! assign the largest remaining bin to the CPU with the smallest load
! gives near perfect distributions for a sufficient number of bins ...
! doing this in chunks of nstep (where nstep ~ number_of_processes) makes this n log n and gives
! each cpu a similar number of tasks.
! it also avoids degenerate cases where thousands of zero sized tasks
! are assigned to the same (least loaded) cpu
!
    DO i=total_number_of_bins,1,-nstep
       tmp_cpu_cost=my_cost_cpu
       CALL sort(tmp_cpu_cost,INT(number_of_processes),tmp_cpu_index)
       DO j=1,MIN(i,nstep)
          itmp=tmp_cpu_index(j)
          distribution_vector(tmp_index(i-j+1))=itmp
          my_cost_cpu(itmp)=my_cost_cpu(itmp)+bin_costs(tmp_index(i-j+1))
       ENDDO
    ENDDO

    DEALLOCATE(tmp_cost,tmp_index,tmp_cpu_cost,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tmp_cpu_index,my_cost_cpu,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE optimize_distribution

END MODULE hfx_load_balance_methods
