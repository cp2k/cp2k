!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Third layer of the dbcsr matrix-matrix multiplication.
!> 	   It collects the full matrix blocks, which need to be multiplied,
!>         and stores their parameters in various stacks.
!>         After a certain amount of parameters is collected it dispatches 
!>         the filled stacks to either the CPU or the GPU.
!>
!> \author  Urban Borstnik
!>
!> <b>Modification history:</b>
!>  - 2010-02-23 Moved from dbcsr_operations
!>  - 2011-11    Moved parameter-stack processing routines to
!>               dbcsr_mm_methods.
!>  - 2013-01    reorganized code (Ole Schuett)
! *****************************************************************************

MODULE dbcsr_mm_csr

  USE array_types,                     ONLY: array_data
  USE dbcsr_config,                    ONLY: dbcsr_get_conf_nstacks,&
                                             mm_async,&
                                             mm_driver,&
                                             mm_driver_cuda,&
                                             mm_host_driver,&
                                             mm_stack_size,&
                                             mm_thread_workshare,&
                                             use_comm_thread
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: dp,&
                                             int_1,&
                                             int_4,&
                                             int_8,&
                                             real_8,&
                                             sp
  USE dbcsr_machine,                   ONLY: m_walltime
  USE dbcsr_message_passing,           ONLY: mp_testany
  USE dbcsr_mm_driver,                 ONLY: dbcsr_mm_driver_inner_init,&
                                             dbcsr_mm_driver_outer_finalize,&
                                             dbcsr_mm_driver_outer_init
  USE dbcsr_mm_stack,                  ONLY: &
       enqueue_ps_group, get_stack_or_process_queue, process_queue_mine, &
       process_queue_others, process_queue_preempt, t_calc_step, t_dev_idle, &
       t_dev_sync, t_process_stack
  USE dbcsr_mm_types,                  ONLY: carrier_type,&
                                             ele_type,&
                                             hash_table_type
  USE dbcsr_pq_methods,                ONLY: &
       dbcsr_pq_create, dbcsr_pq_destroy, dbcsr_pq_flush_level_chg, &
       dbcsr_pq_handoff_level_chg, dbcsr_pq_mark_done, dbcsr_pq_refresh, &
       dbcsr_pq_share, dbcsr_pq_unshare, dbcsr_ps_set_advance, &
       dbcsr_ps_set_create, dbcsr_ps_set_destroy, dbcsr_ps_set_get_group_p, &
       dbcsr_ps_set_get_n_working, dbcsr_ps_target_add_data, &
       dbcsr_ps_target_new, dbcsr_ps_target_release, dbcsr_psg_add_data_ab, &
       dbcsr_psg_get_state, dbcsr_psg_rm_data_ab, dbcsr_psg_view_close, &
       dbcsr_psg_view_open
  USE dbcsr_pq_types,                  ONLY: &
       dbcsr_pq_type, dbcsr_ps_group_type, dbcsr_ps_set_type, &
       dbcsr_ps_state_queued, dbcsr_ps_target_obj, dbcsr_ps_width, &
       mult_timers, p_a_first, p_b_first, p_c_blk, p_c_first, p_k, p_m, p_n, &
       tmr_mult_finish, tmr_mult_finish1, tmr_mult_finish2, tmr_mult_finish3, &
       tmr_sync
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
  USE dbcsr_toollib,                   ONLY: sort
  USE dbcsr_types,                     ONLY: dbcsr_memory_CUDA_host_pinned,&
                                             dbcsr_memory_default,&
                                             dbcsr_type
  USE dbcsr_util,                      ONLY: map_most_common,&
                                             xtime_set,&
                                             xtime_start,&
                                             xtime_stop

  !$ USE OMP_LIB
        
  IMPLICIT NONE
 
  PRIVATE
 
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_csr'
  
  LOGICAL, PARAMETER :: debug_mod  = .FALSE.
  LOGICAL, PARAMETER :: careful_mod = .FALSE.

  
  !> \var max_stack_block_size  The maximal block size to be specially
  !>                            treated.
  INTEGER, PARAMETER :: max_stack_block_size = HUGE (INT (0))
  
  REAL, PARAMETER                      :: default_resize_factor = 1.618034
  
  REAL(kind=dp)  :: index_time
  LOGICAL, PRIVATE :: do_index_time = .FALSE.

  
  
  PUBLIC :: csr_multiply_outer_init, csr_multiply_outer_finalize
  PUBLIC :: csr_multiply_inner_init, csr_multiply_inner_finalize
  PUBLIC :: csr_multiply_unwrap
  PUBLIC :: index_time
  
  
  CONTAINS

! *****************************************************************************
!> \brief Performs multiplication of smaller submatrices.
!>
!> This routine is used for setting up and calling csr_multiply with
!> as few argument indirections as possible.
! *****************************************************************************
  SUBROUTINE csr_multiply_unwrap(mi, mf, ni, nf, ki, kf,&
       ai, af, a_index, bi, bf, b_index, &
       carrier, error)
    INTEGER, INTENT(IN)                      :: mi, mf, ni, nf, ki, kf, ai, af
    INTEGER, DIMENSION(1:3, 1:af), &
      INTENT(IN)                             :: a_index
    INTEGER, INTENT(IN)                      :: bi, bf
    INTEGER, DIMENSION(1:3, 1:bf), &
      INTENT(IN)                             :: b_index
    TYPE(carrier_type), INTENT(INOUT)        :: carrier
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'csr_multiply_unwrap', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: ithread, max_new_nblks, &
                                                n_a_norms, n_b_norms, &
                                                nblks_new, nstacks
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: a_row_p, b_row_p
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: a_blk_info, b_blk_info
    LOGICAL                                  :: advance_memreg
    REAL(KIND=sp), DIMENSION(1:af-ai+1)      :: csr_a_norms
    REAL(KIND=sp), DIMENSION(1:bf-bi+1)      :: csr_b_norms
    TYPE(dbcsr_ps_group_type), POINTER       :: ps_group

!   ---------------------------------------------------------------------------

    IF (dbg) THEN
       WRITE(*,'(I7,1X,5(A,2(1X,I7)))')0,"uwr", mi, mf,",",ni,nf,",",ki,kf,"/",ai,af,",",bi,bf
    ENDIF

    IF (af-ai+1 .GT. 0 .AND. bf-bi+1 .GT. 0) THEN

       ! the maximum number of blocks can be safely estimated by considering both the rowxcol,
       ! but also the blocks the latter can never be larger than norec**2, which is a 'small' constant
       max_new_nblks = MIN(INT(mf-mi+1,int_8) * INT(nf-ni+1,int_8), &
                           INT(af-ai+1,int_8) * INT(bf-bi+1,int_8))
       nblks_new = carrier%lastblk + max_new_nblks

       CALL ensure_array_size(carrier%my_wm%row_i, ub=nblks_new,&
            factor=default_resize_factor, error=error)
       CALL ensure_array_size(carrier%my_wm%col_i, ub=nblks_new,&
            factor=default_resize_factor, error=error)
       CALL ensure_array_size(carrier%my_wm%blk_p, ub=nblks_new,&
            factor=default_resize_factor, error=error)

       ALLOCATE (a_row_p(mi:mf+1))
       ALLOCATE (b_row_p(ki:kf+1))
       ALLOCATE (a_blk_info(2,af-ai+1))
       ALLOCATE (b_blk_info(2,bf-bi+1))
       !
       IF (carrier%use_eps) THEN
          n_a_norms = af-ai+1
          n_b_norms = bf-bi+1
       ELSE
          n_a_norms = 0
          n_b_norms = 0
       ENDIF
       !
       ! Build the indices
       CALL build_csr_index (mi,mf,ai,af,a_row_p, a_blk_info, a_index,&
            n_a_norms, csr_a_norms, carrier%a_norms)
       CALL build_csr_index (ki,kf,bi,bf,b_row_p, b_blk_info, b_index,&
            n_b_norms, csr_b_norms, carrier%b_norms)
       !
       ! This is called here to get the ps_group, a parameter stack.
       ! The memory region (of the parameter set) it belongs to must
       ! remain unchanged throughout the csr_multiply call while the
       ! buffer is allowed to change.  The parameter stack groups of
       ! all the buffers of a memory region share one "super parameter
       ! stack" (\see dbcsr_ps_set_create).  This super parameter
       ! stack is passed directly to the csr_multiply routine
       ! (ps_group%master%s%parameters) for efficiency reasons.
       CALL get_stack_or_process_queue (carrier%param_sets, carrier%queue,&
            ps_group, driver=mm_host_driver, error=error)
       !
       nstacks = SIZE(ps_group%stacks)
       CALL csr_multiply (mi, mf, ni, nf, ki, kf,&
            ai, af, a_row_p, a_blk_info, bi, bf, b_row_p, b_blk_info,&
            carrier%c_hashes,&
            carrier%my_wm%row_i, carrier%my_wm%col_i, carrier%my_wm%blk_p,&
            carrier%lastblk, carrier%datasize,&
            carrier%m_sizes, carrier%n_sizes, carrier%k_sizes,&
            carrier%c_local_rows, carrier%c_local_cols,&
            carrier%c_has_symmetry, carrier%keep_sparsity, carrier%use_eps,&
            carrier%param_sets,&
            carrier%queue,&
            ps_group%master%s%parameters, advance_memreg,&
            carrier%m_size_maps, carrier%n_size_maps, carrier%k_size_maps,&
            carrier%m_size_maps_size, carrier%n_size_maps_size, carrier%k_size_maps_size,&
            carrier%nm_stacks, carrier%nn_stacks, carrier%nk_stacks, &
            carrier%stack_map,&
            nstacks, carrier%default_stack, &
            carrier%row_max_epss, csr_a_norms, csr_b_norms,&
            carrier%flop, carrier%t_gemm, error)
       !
       ! If a stack group was filled up during the csr_multiply call then
       ! we roll around to a new memory region.
       IF (advance_memreg) THEN
            CALL dbcsr_ps_set_advance (carrier%param_sets,&
                 advance_memreg = .TRUE.,&
                 error=error)
            ! And make some progress in MPI
            ithread = 0
!$          ithread = omp_get_thread_num()
            IF (use_comm_thread .AND. (ithread .EQ. 0)) &
              CALL progress_comms(carrier)
       END IF
       !
       DEALLOCATE (a_row_p)
       DEALLOCATE (b_row_p)
       DEALLOCATE (a_blk_info)
       DEALLOCATE (b_blk_info)
    ENDIF
  END SUBROUTINE csr_multiply_unwrap

! *****************************************************************************
!> \brief Call in MPI to progrss any outstanding communications
!> TODO: copied from dbcsr_internal_operations. MPI does not belong on this level.
! *****************************************************************************
  SUBROUTINE progress_comms(carrier)
    TYPE(carrier_type), INTENT(inout)        :: carrier

    CALL mp_testany(carrier%right_data_sr)
    CALL mp_testany(carrier%right_data_rr)
    CALL mp_testany(carrier%left_data_sr)
    CALL mp_testany(carrier%left_data_rr)
    CALL mp_testany(carrier%right_index_sr)
    CALL mp_testany(carrier%right_index_rr)
    CALL mp_testany(carrier%left_index_sr)
    CALL mp_testany(carrier%left_index_rr)
  END SUBROUTINE
  

! *****************************************************************************
!> \brief Moved out of dbcsr_multrec_init
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE csr_multiply_outer_init(left, right, product, carrier, fill_guess, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(carrier_type), INTENT(INOUT)        :: carrier
    REAL(KIND=real_8)                        :: fill_guess
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'csr_multiply_outer_init', &
      routineP = moduleN//':'//routineN

    INTEGER :: block_estimate, default_stack, istack, ithread, k_map, k_size, &
      m_map, m_size, mem_type, n_map, n_size, n_stack_buffers, &
      n_stack_mem_regions, nstacks, nthreads, ps_buffer, ps_g, ps_memreg
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: flop_index, flop_list, &
                                                most_common_k, most_common_m, &
                                                most_common_n
    INTEGER, DIMENSION(3)                    :: nxstacks
    LOGICAL                                  :: local_cols, local_indexing, &
                                                local_rows
    TYPE(dbcsr_ps_group_type), POINTER       :: ps_group
    TYPE(dbcsr_ps_target_obj)                :: product_target

    ithread = 0 ; nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM () ; nthreads = OMP_GET_NUM_THREADS ()     
    
    ! find out if we have local_indexing
    local_cols = right%local_indexing
    local_rows = left%local_indexing
    CALL dbcsr_assert (local_cols, "EQV", local_rows,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Local index useage must be consistent.", __LINE__, error=error)
    local_indexing = local_rows
    
    
    ! Setup the hash tables if needed
    IF (local_indexing) THEN
       ALLOCATE (carrier%c_hashes (product%nblkrows_local))
    ELSE
       ALLOCATE (carrier%c_hashes (product%nblkrows_total))
    ENDIF
    block_estimate=MAX(product%nblks,left%nblks,right%nblks)/nthreads
    IF (local_indexing) THEN
       CALL fill_hash_tables (carrier%c_hashes, product,block_estimate,&
            row_map=array_data(product%global_rows),&
            col_map=array_data(product%global_cols),&
            error=error)
    ELSE
       CALL fill_hash_tables (carrier%c_hashes, product,block_estimate,&
            error=error)
    ENDIF
    !
    ! Setup the MM stack
    CALL dbcsr_get_conf_nstacks (nxstacks, n_stack_buffers, n_stack_mem_regions,&
         error)
    carrier%nm_stacks = nxstacks(1)
    carrier%nn_stacks = nxstacks(2)
    carrier%nk_stacks = nxstacks(3)
    nstacks = nxstacks(1) * nxstacks(2) * nxstacks(3) + 1
    CALL dbcsr_assert (nstacks, "LE", INT (HUGE (carrier%stack_map)),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Too many stacks requested (global/dbcsr/n_size_*_stacks in input)",&
         __LINE__, error=error)
    !
    CALL dbcsr_ps_target_new (product_target, error)
    CALL dbcsr_ps_target_add_data (product_target%t,&
         product%wms(ithread+1)%data_area, carrier%datasize, error=error)
    
    IF (mm_driver .EQ. mm_driver_cuda) THEN
       mem_type = dbcsr_memory_CUDA_host_pinned
    ELSE
       mem_type = dbcsr_memory_default
    ENDIF
    
   ! Here the first driver init part used to be
    
    CALL dbcsr_ps_set_create (carrier%param_sets, mm_stack_size,&
         product_target, nstacks, n_stack_buffers, n_stack_mem_regions,&
         mem_type,error=error)
    
    !moved this to the end as product_target is need for driver_init    
    !CALL dbcsr_ps_target_release (product_target, error)
    
    ! Setup the block sizes mappings.  The carrier%*_sizes arrays are
    ! mappings from block sizes to an index from the number of most
    ! common block sizes.  I.e., if nxstacks(1) = 2 and the most
    ! common blocks sizes are 13 and 5, then the values in
    ! carrier%m_sizes will be 3 for all elements except 13 and 5;
    ! these two elemnts will be either 1 or 2.
    !
    ALLOCATE (most_common_m(nxstacks(1)))
    ALLOCATE (most_common_n(nxstacks(2)))
    ALLOCATE (most_common_k(nxstacks(3)))
    CALL map_most_common (carrier%m_sizes, carrier%m_size_maps, nxstacks(1),&
         most_common_m,&
         max_stack_block_size, carrier%max_m)
    carrier%m_size_maps_size = SIZE (carrier%m_size_maps)
    CALL map_most_common (carrier%n_sizes, carrier%n_size_maps, nxstacks(2),&
         most_common_n,&
         max_stack_block_size, carrier%max_n)
    carrier%n_size_maps_size = SIZE (carrier%n_size_maps)
    CALL map_most_common (array_data(right%row_blk_size),&
         carrier%k_size_maps, nxstacks(3), &
         most_common_k,&
         max_stack_block_size, carrier%max_k)
    carrier%k_size_maps_size = SIZE (carrier%k_size_maps)
    !
    ! Creates the stack map--a mapping from (mapped) stack block sizes
    ! (carrier%*_sizes) to a stack number.  Triples with even one
    ! uncommon size will be mapped to a general, non-size-specific
    ! stack.
    ALLOCATE (carrier%stack_map(nxstacks(2)+1, nxstacks(3)+1, nxstacks(1)+1))
    default_stack = nstacks
    carrier%default_stack = default_stack
    DO m_map = 1, nxstacks(1)+1
       IF (m_map .LE. nxstacks(1)) THEN
          m_size = most_common_m(m_map)
       ELSE
          m_size = 777
       ENDIF
       DO k_map = 1, nxstacks(3)+1
          IF (k_map .LE. nxstacks(3)) THEN
             k_size = most_common_k(k_map)
          ELSE
             k_size = 888
          ENDIF
          DO n_map = 1, nxstacks(2)+1
             IF (n_map .LE. nxstacks(2)) THEN
                n_size = most_common_n(n_map)
             ELSE
                n_size = 999
             ENDIF
             IF (       m_map .LE. nxstacks(1)&
                  .AND. k_map .LE. nxstacks(3)&
                  .AND. n_map .LE. nxstacks(2)) THEN
                ! This is the case when m, n, and k are all defined.
                ps_g = (m_map-1)*nxstacks(2)*nxstacks(3) +&
                       (k_map-1)*nxstacks(2) +&
                       n_map
                ps_g = nstacks-ps_g
                carrier%stack_map(n_map, k_map, m_map) = ps_g
                ! Also take care of the stack m, n, k descriptors
                DO ps_memreg = 1, n_stack_mem_regions
                   DO ps_buffer = 1, n_stack_buffers
                      ps_group => dbcsr_ps_set_get_group_p (carrier%param_sets,&
                           ps_buffer, ps_memreg, error=error)
                      ps_group%stacks(ps_g)%s%m     = m_size
                      ps_group%stacks(ps_g)%s%n     = n_size
                      ps_group%stacks(ps_g)%s%k     = k_size
                      ps_group%stacks(ps_g)%s%max_m = m_size
                      ps_group%stacks(ps_g)%s%max_n = n_size
                      ps_group%stacks(ps_g)%s%max_k = k_size
                      ps_group%stacks(ps_g)%s%defined_mnk = .TRUE.
                      ! The bins are quasi-sorted according from those
                      ! with largest sizes (best for accelerators) to
                      ! those with smallest (worst for accelerators).
                      ps_group%stacks(ps_g)%s%size_bin = ps_g
                   ENDDO
                ENDDO
             ELSE
                ! This is the case when at least one of m, n, or k is
                ! undefined.
                ps_g = default_stack
                carrier%stack_map(n_map, k_map, m_map) = default_stack
                ! Also take care of the stack m, n, k descriptors
                DO ps_memreg = 1, n_stack_mem_regions
                   DO ps_buffer = 1, n_stack_buffers
                      ps_group => dbcsr_ps_set_get_group_p (carrier%param_sets,&
                           ps_buffer, ps_memreg, error=error)
                      ps_group%stacks(ps_g)%s%m     = 0
                      ps_group%stacks(ps_g)%s%n     = 0
                      ps_group%stacks(ps_g)%s%k     = 0
                      ps_group%stacks(ps_g)%s%max_m = carrier%max_m
                      ps_group%stacks(ps_g)%s%max_n = carrier%max_n
                      ps_group%stacks(ps_g)%s%max_k = carrier%max_k
                      ps_group%stacks(ps_g)%s%defined_mnk = .FALSE.
                      ps_group%stacks(ps_g)%s%size_bin = default_stack
                   ENDDO
                ENDDO
             END IF
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE (most_common_m)
    DEALLOCATE (most_common_n)
    DEALLOCATE (most_common_k)

    ! sort to make the order fixed... all defined stacks first, default stack
    ! last. Next, sort according to flops, first stack lots of flops, last
    ! stack, few flops
    IF (nstacks>1) THEN
       ps_group => dbcsr_ps_set_get_group_p (carrier%param_sets, 1, 1, error=error)
       ALLOCATE(flop_list(nstacks-1),flop_index(nstacks-1))
       DO istack=1,nstacks-1
          flop_list(istack)=-2*ps_group%stacks(istack)%s%m*ps_group%stacks(istack)%s%n*ps_group%stacks(istack)%s%k
       ENDDO
       CALL sort(flop_list,nstacks-1,flop_index)
       DO ps_memreg = 1, n_stack_mem_regions
          DO ps_buffer = 1, n_stack_buffers
             ps_group => dbcsr_ps_set_get_group_p (carrier%param_sets,&
                           ps_buffer, ps_memreg, error=error)
             DO istack=1,nstacks-1
                ps_group%stacks(istack)%s%size_bin = flop_index(istack)
             ENDDO
          ENDDO
       ENDDO
       DEALLOCATE(flop_list,flop_index)
    ENDIF

    CALL dbcsr_mm_driver_outer_init(left, right, product, product_target, carrier, fill_guess, error=error)
    
    CALL dbcsr_ps_target_release (product_target, error)
    
    ALLOCATE (carrier%queue)
    CALL dbcsr_pq_create (carrier%queue, nstacks, error)
    CALL dbcsr_pq_share (carrier%queue, error)
    !
  END SUBROUTINE csr_multiply_outer_init
  

! *****************************************************************************
!> \brief Fills row hashtable from an existing matrix.
!> \param[in] block_estimate guess for the number of blocks in the product matrix, can be zero
! *****************************************************************************
  SUBROUTINE fill_hash_tables(hashes, matrix, block_estimate, row_map, col_map, error)
    TYPE(hash_table_type), DIMENSION(:), &
      INTENT(inout)                          :: hashes
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER                                  :: block_estimate
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: row_map, col_map
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fill_hash_tables', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, error_handler, i, imat, &
                                                n_rows, row

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    imat = 1
    !$ imat = OMP_GET_THREAD_NUM() + 1
    IF (PRESENT (row_map)) THEN
       n_rows = matrix%nblkrows_local
       CALL dbcsr_assert (SIZE(hashes), "EQ", n_rows,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Local row count mismatch", __LINE__, error=error)
    ELSE
       n_rows = matrix%nblkrows_total
       CALL dbcsr_assert (SIZE(hashes), "EQ", n_rows,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Global row count mismatch", __LINE__, error=error)
    ENDIF
    DO row = 1, n_rows
       ! create the hash table row with a reasonable initial size
       CALL hash_table_create (hashes(row), &
            MAX(8,(3*block_estimate)/MAX(1,n_rows)))
    ENDDO
    ! We avoid using the iterator because we will use the existing
    ! work matrix instead of the BCSR index.
    DO i = 1, matrix%wms(imat)%lastblk
       row = matrix%wms(imat)%row_i(i)
       col = matrix%wms(imat)%col_i(i)
       IF (PRESENT (row_map)) row = row_map(row)
       IF (PRESENT (col_map)) col = col_map(col)
       CALL hash_table_add(hashes(row), col, i, error=error)
    ENDDO
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE fill_hash_tables

 
! *****************************************************************************
!> \brief Moved out of dbcsr_multrec_finalize
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE csr_multiply_outer_finalize(carrier, product, error)
    TYPE(carrier_type), INTENT(INOUT)        :: carrier
    TYPE(dbcsr_type), INTENT(inout)          :: product
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'csr_multiply_outer_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    CALL dbcsr_mm_driver_outer_finalize(carrier, product, error)
    
! Clear hash tables

    DO i = 1, SIZE(carrier%c_hashes)
       CALL hash_table_release (carrier%c_hashes (i))
    ENDDO
    DEALLOCATE (carrier%c_hashes)
                

    
    ! Clear MM stack
    !CALL dbcsr_psg_rm_data_c (&
    !     dbcsr_ps_set_get_group(carrier%param_sets, 1, 1, error), error=error)

    ! maybe a call to dbcsr_psg_rm_data_ab in a critical section would be enough
    !$OMP CRITICAL (crit_data)
    CALL dbcsr_ps_set_destroy (carrier%param_sets, error=error)
    !$OMP END CRITICAL (crit_data)

    
     
    CALL dbcsr_pq_unshare (carrier%queue, error)
    CALL dbcsr_pq_destroy (carrier%queue, error)
    DEALLOCATE (carrier%queue)            
                
  END SUBROUTINE csr_multiply_outer_finalize
! *****************************************************************************
!> \brief Moved out of dbcsr_multrec_sparse
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE csr_multiply_inner_init(left, right, carrier, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(carrier_type), INTENT(INOUT)        :: carrier
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: ps_buffer, ps_memreg
    TYPE(dbcsr_ps_group_type), POINTER       :: stack_group

    !$OMP MASTER
    t_calc_step = 0.0_dp
    t_process_stack = 0.0_dp
    t_dev_sync = 0.0_dp
    t_dev_idle = 0.0_dp
    !$OMP END MASTER
    !$OMP BARRIER
    
! Setup the carrier data that changes in each multiplication step.

    DO ps_memreg = 1, carrier%param_sets%nmemregs
       DO ps_buffer = 1, carrier%param_sets%nbuffers
          stack_group => dbcsr_ps_set_get_group_p (carrier%param_sets,&
               ps_buffer, ps_memreg, error=error)
          !$OMP CRITICAL (crit_data)
          CALL dbcsr_psg_add_data_ab (stack_group,&
               left%data_area, right%data_area, error=error)
          !$OMP END CRITICAL (crit_data)
       ENDDO
    ENDDO
    !
    ! Quick reset of the queue to avoid heavyweight de/allocations.
    CALL dbcsr_pq_refresh (carrier%queue,&
         flush_level=(carrier%param_sets%nbuffers*carrier%param_sets%group_size),&
         error=error)
    
    CALL dbcsr_mm_driver_inner_init(left, right, carrier, error)
    !$OMP BARRIER
         
  END SUBROUTINE csr_multiply_inner_init
  
! *****************************************************************************
!> \brief Moved out of dbcsr_multrec_sparse
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE csr_multiply_inner_finalize(carrier, error)
    TYPE(carrier_type), INTENT(INOUT)        :: carrier
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: ps_buffer, ps_memreg
    LOGICAL                                  :: i_interface
    TYPE(dbcsr_ps_group_type), POINTER       :: stack_group

    CALL xtime_set (mult_timers(tmr_mult_finish))
    ! Mark myself as done with indexing.
    CALL dbcsr_pq_mark_done (carrier%queue, error)
    i_interface = .FALSE.
    IF (mm_thread_workshare) THEN
       !$OMP MASTER
       !
       ! I am the thread that interfaces to the accelerator.
       i_interface = mm_async
       !$OMP END MASTER
    ELSE
       i_interface = mm_async
    ENDIF
    !
    ! Queue all outstanding stacks.
    DO ps_memreg = 1, carrier%param_sets%nmemregs
       DO ps_buffer = 1, carrier%param_sets%nbuffers
          stack_group => dbcsr_ps_set_get_group_p (carrier%param_sets,&
               wait=.FALSE., error=error)
          ! Only enqueue the stack group if it is not already queued
          ! or being processed.
          IF (dbcsr_psg_get_state (stack_group, error) &
               .LT. dbcsr_ps_state_queued) THEN
             CALL enqueue_ps_group (carrier%queue,&
                  stack_group, error)
          ENDIF
       ENDDO
    ENDDO
    !
    ! Process my queues now.
    IF (.NOT. i_interface) THEN
       ! Non-master threads just do CPU calculations.
       CALL xtime_set (mult_timers(tmr_mult_finish1))                           
       CALL process_queue_mine(carrier%queue, driver=mm_host_driver,&
            error=error)
       CALL xtime_stop (mult_timers(tmr_mult_finish1))
    ELSE
       ! The master thread tries to fill up the accelerator.
       CALL xtime_set (mult_timers(tmr_mult_finish1))
       CALL process_queue_mine (carrier%queue, driver=mm_driver,&
            error=error)
       CALL xtime_stop (mult_timers(tmr_mult_finish1))

       IF (mm_thread_workshare) THEN
          CALL xtime_set (mult_timers(tmr_mult_finish2))
          CALL process_queue_others (carrier%queue, forever=.FALSE.,&
               driver=mm_driver, error=error)
          CALL xtime_stop (mult_timers(tmr_mult_finish2))
       ENDIF

       CALL xtime_set (mult_timers(tmr_mult_finish1))
       CALL process_queue_mine (carrier%queue,&
            driver=mm_host_driver, error=error)
       CALL xtime_stop (mult_timers(tmr_mult_finish1))

       IF (mm_thread_workshare) THEN
          CALL xtime_set (mult_timers(tmr_mult_finish3))
          CALL process_queue_others (carrier%queue, forever=.TRUE.,&
               driver=mm_driver, error=error)
          CALL xtime_stop (mult_timers(tmr_mult_finish3))
       ENDIF
    ENDIF
    CALL xtime_stop (mult_timers(tmr_mult_finish))
    
    IF (mm_thread_workshare) THEN
       CALL xtime_start (mult_timers(tmr_sync))
       !$OMP BARRIER
       CALL xtime_stop (mult_timers(tmr_sync))
    ENDIF
    
    !TODO: why is this only called on one stack_group see also csr_multiply_init
    !
    !$OMP CRITICAL (crit_data)
    CALL dbcsr_psg_rm_data_ab (stack_group, error=error)
    !$OMP END CRITICAL (crit_data)
  
  END SUBROUTINE csr_multiply_inner_finalize
  

! *****************************************************************************
!> \brief Performs multiplication of smaller submatrices.
!> \note The queue parameter could be INTENT(INOUT) when OpenMP is not used.
! *****************************************************************************
  SUBROUTINE csr_multiply(mi, mf, ni, nf, ki, kf,&
       ai, af, a_row_p, a_blk_info, bi, bf, b_row_p, b_blk_info,&
       c_hashes, c_row_i, c_col_i, c_blk_p, lastblk, datasize,&
       m_sizes, n_sizes, k_sizes,&
       c_local_rows, c_local_cols,&
       c_has_symmetry, keep_sparsity, use_eps,&
       param_sets, queue,&
       params_array, advance_memreg, &
       row_size_maps, col_size_maps, k_size_maps,&
       row_size_maps_size, col_size_maps_size, k_size_maps_size,&
       nm_stacks, nn_stacks, nk_stacks, stack_map,&
       nstacks, default_stack,&
       row_max_epss, left_norms, right_norms,&
       flop, t_gemm, error)
    INTEGER, INTENT(IN)                      :: mi, mf, ni, nf, ki, kf, ai, af
    INTEGER, DIMENSION(mi:mf+1), INTENT(IN)  :: a_row_p
    INTEGER, DIMENSION(2, 1:af-ai+1), &
      INTENT(IN)                             :: a_blk_info
    INTEGER, INTENT(IN)                      :: bi, bf
    INTEGER, DIMENSION(ki:kf+1), INTENT(IN)  :: b_row_p
    INTEGER, DIMENSION(2, 1:bf-bi+1), &
      INTENT(IN)                             :: b_blk_info
    TYPE(hash_table_type), DIMENSION(:), &
      INTENT(INOUT)                          :: c_hashes
    INTEGER, DIMENSION(*), INTENT(INOUT)     :: c_row_i, c_col_i, c_blk_p
    INTEGER, INTENT(INOUT)                   :: lastblk, datasize
    INTEGER, DIMENSION(*), INTENT(IN)        :: m_sizes, n_sizes, k_sizes, &
                                                c_local_rows, c_local_cols
    LOGICAL, INTENT(IN)                      :: c_has_symmetry, &
                                                keep_sparsity, use_eps
    TYPE(dbcsr_ps_set_type), INTENT(INOUT)   :: param_sets
    TYPE(dbcsr_pq_type), POINTER             :: queue
    INTEGER, DIMENSION(*), INTENT(INOUT)     :: params_array
    LOGICAL, INTENT(OUT)                     :: advance_memreg
    INTEGER, INTENT(IN)                      :: row_size_maps_size, &
                                                k_size_maps_size, &
                                                col_size_maps_size
    INTEGER(KIND=int_4), &
      DIMENSION(0:row_size_maps_size-1), &
      INTENT(IN)                             :: row_size_maps
    INTEGER(KIND=int_4), &
      DIMENSION(0:col_size_maps_size-1), &
      INTENT(IN)                             :: col_size_maps
    INTEGER(KIND=int_4), &
      DIMENSION(0:k_size_maps_size-1), &
      INTENT(IN)                             :: k_size_maps
    INTEGER, INTENT(IN)                      :: nm_stacks, nn_stacks, &
                                                nk_stacks
    INTEGER(KIND=int_1), DIMENSION(&
      nn_stacks+1, nk_stacks+1, nm_stacks+1)&
      , INTENT(IN)                           :: stack_map
    INTEGER, INTENT(IN)                      :: nstacks, default_stack
    REAL(kind=sp), DIMENSION(*)              :: row_max_epss, left_norms, &
                                                right_norms
    INTEGER(KIND=int_8), INTENT(INOUT)       :: flop
    REAL(KIND=dp), INTENT(INOUT)             :: t_gemm
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'csr_multiply', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE., &
                                                local_timing = .FALSE.

    INTEGER :: a_blk, a_col_l, a_row_l, b_blk, b_col_l, c, c_blk_id, &
      c_blk_pt, c_col_logical, c_nze, c_row_logical, ithread, k_size, m_size, &
      mapped_col_size, mapped_k_size, mapped_row_size, n_size, new_blk, &
      nworking, s_dp, ws, zero_first, zero_last
    INTEGER, DIMENSION(nstacks)              :: param_starts, stack_p
    INTEGER, DIMENSION(:), POINTER           :: params_array_direct
    INTEGER(KIND=int_4)                      :: offset
    LOGICAL                                  :: block_exists, do_preempt, &
                                                flush_stack, i_interface
    REAL(kind=dp)                            :: t_gemm_me
    REAL(kind=sp)                            :: a_norm, a_row_eps, b_norm
    TYPE(dbcsr_ps_group_type), POINTER       :: ps_group

!   ---------------------------------------------------------------------------

    IF (do_index_time) &
         index_time = -m_walltime()
    ithread = 0
    !$ ithread = omp_get_thread_num()
    t_gemm_me = REAL(0, KIND(t_gemm_me))
    !
    advance_memreg = .FALSE.
    i_interface = .FALSE.
    IF (mm_thread_workshare) THEN
       !$OMP MASTER
       ! I am the thread that interfaces to the accelerator.
       i_interface = mm_async
       !$OMP END MASTER
    ELSE
       i_interface = mm_async
    ENDIF
    !
    ! Each thread always shares the same zero_first, zero_last, and
    ! new_blk.  All of the parameter stacks should have the same
    ! value.
    ps_group => dbcsr_ps_set_get_group_p (param_sets, error=error)
    CALL dbcsr_psg_view_open (ps_group, params_array_direct, stack_p,&
         zero_first, zero_last, new_blk, param_starts, error=error)
    !
    ! New data blocks will be put into the data area starting at
    ! the c_blk_pt position.
    new_blk = lastblk
    c_blk_pt = datasize + 1
    ws = 1
    !
    a_row_cycle: DO a_row_l = mi, mf
       m_size = m_sizes(a_row_l)

       a_row_eps = row_max_epss (a_row_l)
       mapped_row_size = row_size_maps(m_size)

       a_blk_cycle: DO a_blk = a_row_p(a_row_l)+1, a_row_p(a_row_l+1)
          a_col_l = a_blk_info(1, a_blk)
          IF (debug_mod) WRITE(*,*)ithread,routineN//" A col", a_col_l,";",a_row_l
          k_size = k_sizes (a_col_l)
          mapped_k_size = k_size_maps(k_size)

          a_norm = left_norms(a_blk)
          b_blk_cycle: DO b_blk = b_row_p(a_col_l)+1, b_row_p(a_col_l+1)
             IF (dbg) THEN
                WRITE(*,'(1X,A,3(1X,I7),1X,A,1X,I16)')routineN//" trying B",&
                     a_row_l, b_blk_info(1,b_blk), a_col_l, "at", b_blk_info(2,b_blk)
             ENDIF
             b_norm = right_norms(b_blk)
             IF (a_norm * b_norm .LT. a_row_eps) THEN
                CYCLE
             ENDIF
             b_col_l = b_blk_info(1,b_blk)
             ! Don't calculate symmetric blocks.
             symmetric_product: IF (c_has_symmetry) THEN
                c_row_logical = c_local_rows (a_row_l)
                c_col_logical = c_local_cols (b_col_l)
                IF (c_row_logical .NE. c_col_logical&
                     .AND. my_checker_tr (c_row_logical, c_col_logical)) THEN
                   IF (dbg) THEN
                      WRITE(*,*)"Skipping symmetric block!", c_row_logical,&
                           c_col_logical
                   ENDIF
                   CYCLE
                ENDIF
             ENDIF symmetric_product

             c_blk_id = hash_table_get (c_hashes(a_row_l), b_col_l)
             IF (.FALSE.) THEN
                WRITE(*,'(1X,A,3(1X,I7),1X,A,1X,I16)')routineN//" coor",&
                     a_row_l, a_col_l, b_col_l,"c blk", c_blk_id
             ENDIF
             block_exists = c_blk_id .GT. 0

             n_size = n_sizes(b_col_l)
             c_nze = m_size * n_size
             !
             new_block_case: IF (.NOT. block_exists) THEN
                sparsity_enforcement: IF (keep_sparsity) THEN
                   CYCLE
                ENDIF sparsity_enforcement
                offset = c_blk_pt
                new_blk = new_blk+1
                c_blk_id = new_blk
                IF (dbg) WRITE(*,*)routineN//" new block offset, nze", offset, c_nze
                CALL hash_table_add(c_hashes(a_row_l),&
                     b_col_l, c_blk_id, error=error)
                !
                ! We still keep the linear index because it's
                ! easier than getting the values out of the
                ! hashtable in the end.
                c_row_i(new_blk) = a_row_l
                c_col_i(new_blk) = b_col_l
                c_blk_p(new_blk) = offset
                !
                c_blk_pt = c_blk_pt + c_nze
                datasize = datasize + c_nze
                zero_last = c_blk_pt - 1
             ELSE
                offset = c_blk_p(c_blk_id)
             ENDIF new_block_case
             !
             ! We should not call certain MM routines (netlib BLAS)
             ! with zero LDs; however, we still need to get to here
             ! to get new blocks.
             IF (careful_mod) THEN
                IF (c_nze .EQ. 0 .OR. k_size .EQ. 0) THEN
                   CALL dbcsr_assert (.FALSE.,&
                        dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                        "Can not call MM with LDx=0.", __LINE__, error=error)
                   CYCLE
                ENDIF
             ENDIF
             !
             mapped_col_size = col_size_maps (n_size)
             ws = stack_map (mapped_col_size, mapped_k_size, mapped_row_size)
             s_dp = (param_starts(ws)-1 + stack_p(ws)) * dbcsr_ps_width
             stack_p(ws) = stack_p(ws) + 1
             !
             params_array(s_dp+p_m) = m_size
             params_array(s_dp+p_n) = n_size
             params_array(s_dp+p_k) = k_size
             !
             params_array(s_dp+p_a_first) = a_blk_info(2, a_blk)
             !
             params_array(s_dp+p_b_first) = b_blk_info(2, b_blk)
             !
             params_array(s_dp+p_c_first) = offset
             params_array(s_dp+p_c_blk) = c_blk_id
             !
             flop = flop + INT(2*c_nze, int_8) * INT(k_size, int_8)
             !
             flush_stack = (stack_p(ws) .GE. mm_stack_size-1)
             IF (flush_stack) THEN
                advance_memreg = .TRUE.
                IF (do_index_time) index_time = index_time + m_walltime()
                CALL dbcsr_psg_view_close (ps_group, params_array_direct,&
                     stack_p, zero_first, zero_last, new_blk, error=error)
                do_preempt = mm_async
                IF (mm_async) THEN
                   nworking = dbcsr_ps_set_get_n_working (param_sets, error)
                   queue%nworking = nworking
                   IF (nworking .EQ. 0) THEN
                      !> If no stack in currently being processed then it
                      !> is assumed the accelerator can handle the stacks
                      !> on its own.  The queue's flush_level is thus
                      !> increased, making it less favorable for the host
                      !> CPU to process stacks from the queue.
                      !> \see process_queue_preempt
                      !> \see get_stack_or_process_queue
                      CALL dbcsr_pq_flush_level_chg (queue,&
                           MAX(1,param_sets%nbuffers/1), param_sets)
                      CALL dbcsr_pq_handoff_level_chg (queue,&
                           MAX(1,param_sets%group_size/2), param_sets)
                      do_preempt = .FALSE.
                   ENDIF
                ENDIF
                CALL enqueue_ps_group (queue, ps_group, error=error)
                !
                IF (i_interface) THEN
                   CALL process_queue_mine (queue, &
                        driver=mm_driver, error=error)
                   IF (mm_thread_workshare) THEN
                      CALL process_queue_others (queue, forever=.FALSE.,&
                           driver=mm_driver, error=error)
                   ENDIF
                ENDIF
                IF (do_preempt) THEN
                   ! Do some CPU stack calculations when needed.
                   CALL process_queue_preempt (queue, param_sets,&
                        driver=mm_host_driver, error=error)
                ENDIF
                !
                CALL dbcsr_ps_set_advance (param_sets,&
                     advance_memreg = .FALSE.,&
                     error=error)
                !
                CALL get_stack_or_process_queue (param_sets,&
                     queue, ps_group, driver=mm_host_driver, error=error)
                !
                CALL dbcsr_psg_view_open (ps_group, params_array_direct,&
                     stack_p, zero_first, zero_last, new_blk,&
                     param_starts, error=error)
                IF (do_index_time) index_time = index_time - m_walltime()
             ENDIF
          ENDDO b_blk_cycle ! b
       ENDDO a_blk_cycle ! a_col
    ENDDO a_row_cycle ! a_row
    !
    lastblk = new_blk
    datasize = c_blk_pt -1
    !
    CALL dbcsr_psg_view_close (ps_group, params_array_direct,&
         stack_p, zero_first, zero_last, new_blk, error=error)

    !
    t_gemm = t_gemm + t_gemm_me
    !
    IF (do_index_time) index_time = index_time + m_walltime()
  END SUBROUTINE csr_multiply


! *****************************************************************************
!> \brief  Builds and sorts a CSR index from a list index.
!> \author JV
!> <b>Modification history:</b>
!> - 2011-02-15 [UB] Adapted to use DBCSR-type CSR indexing
! *****************************************************************************
!  PURE SUBROUTINE build_csr_index(mi,mf,ai,af, row_p, blk_info, list_index)
  SUBROUTINE build_csr_index(mi,mf,ai,af, row_p, blk_info, list_index,&
       nnorms, csr_norms, list_norms)
    INTEGER, INTENT(IN)                      :: mi, mf, ai, af
    INTEGER, DIMENSION(mi:mf+1), INTENT(OUT) :: row_p
    INTEGER, DIMENSION(2, 1:af-ai+1), &
      INTENT(OUT)                            :: blk_info
    INTEGER, DIMENSION(3, 1:af), INTENT(IN)  :: list_index
    INTEGER, INTENT(IN)                      :: nnorms
    REAL(KIND=sp), DIMENSION(1:af-ai+1), &
      INTENT(OUT)                            :: csr_norms
    REAL(KIND=sp), DIMENSION(:), INTENT(IN)  :: list_norms

    CHARACTER(len=*), PARAMETER :: routineN = 'build_csr_index', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = .FALSE., &
                                                dbg = .FALSE.

    INTEGER                                  :: i, row
    INTEGER, DIMENSION(mi:mf)                :: counts
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------
! Counts blocks per row and calculates the offsets.

    IF (dbg) THEN
       WRITE(*,'(I7,1X,5(A,2(1X,I7)))')0,"bci", mi, mf,";",ai,af
       !write(*,'(3(I7))')list_index(:,ai:af)
    ENDIF

    counts(:) = 0
    DO i = ai, af
       IF (careful) THEN
          CALL dbcsr_assert (list_index(1,i), "GE", mi,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Out of range", __LINE__, error=error)
          CALL dbcsr_assert (list_index(1,i), "LE", mf,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Out of range", __LINE__, error=error)
       ENDIF
       counts(list_index(1,i)) = counts(list_index(1,i))+1
    ENDDO
    row_p(mi) = 0
    DO i = mi+1, mf+1
       row_p(i) = row_p(i-1) + counts(i-1)
    ENDDO
    ! Adds every block to its corresponding row.
    counts(:) = 0
    DO i = ai, af
       row = list_index(1,i)
       counts(row) = counts(row)+1
       IF (careful) THEN
          CALL dbcsr_assert (row_p(row) + counts(row), "LE", af-ai+1,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Out of range", __LINE__, error=error)
          CALL dbcsr_assert (row_p(row) + counts(row), "GE", 1,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Out of range", __LINE__, error=error)
       ENDIF
       blk_info(1, row_p(row) + counts(row)) = list_index(2,i)
       blk_info(2, row_p(row) + counts(row)) = list_index(3,i)
       IF (nnorms .GT. 0) THEN
          csr_norms(row_p(row) + counts(row)) = list_norms(i)
       ENDIF
    ENDDO
    IF (nnorms .EQ. 0) THEN
       csr_norms(:) = 0.0_sp
    ENDIF
  END SUBROUTINE build_csr_index
 
! *****************************************************************************
!> \brief Determines whether a transpose must be applied
!> \par Source
!> This function is copied from dbcsr_dist_operations for speed reasons.
!> \param[in] row   The absolute matrix row.
!> \param[in] column          The absolute matrix column.
! *****************************************************************************
  ELEMENTAL FUNCTION my_checker_tr(row, column) RESULT(transpose)
    INTEGER, INTENT(IN)                      :: row, column
    LOGICAL                                  :: transpose

    transpose = BTEST(column+row, 0) .EQV. column.GE.row

  END FUNCTION my_checker_tr
  

! -----------------------------------------------------------------------------
! Beginning of hashtable
  ! finds a prime equal or larger than i
  FUNCTION matching_prime(i) RESULT(res)
    INTEGER, INTENT(IN)                      :: i
    INTEGER                                  :: res

    INTEGER                                  :: j

    res=i
    j=0
    DO WHILE (j<res)
      DO j=2,res-1
         IF (MOD(res,j)==0) THEN
            res=res+1
            EXIT
         ENDIF
      ENDDO
    ENDDO
  END FUNCTION

  SUBROUTINE hash_table_create(hash_table,table_size)
    TYPE(hash_table_type)                    :: hash_table
    INTEGER, INTENT(IN)                      :: table_size

    INTEGER                                  :: j

! guarantee a minimal hash table size (8), so that expansion works

   j=3
   DO WHILE(2**j-1<table_size)
      j=j+1
   ENDDO
   hash_table%nmax=2**j-1
   hash_table%prime=matching_prime(hash_table%nmax)
   hash_table%nele=0
   ALLOCATE(hash_table%table(0:hash_table%nmax))
  END SUBROUTINE hash_table_create

  SUBROUTINE hash_table_release(hash_table)
    TYPE(hash_table_type)                    :: hash_table

   hash_table%nmax=0
   hash_table%nele=0
   DEALLOCATE(hash_table%table)

  END SUBROUTINE hash_table_release

  RECURSIVE SUBROUTINE hash_table_add(hash_table,c,p, error)
    TYPE(hash_table_type), INTENT(INOUT)     :: hash_table
    INTEGER, INTENT(IN)                      :: c, p
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    REAL(KIND=real_8), PARAMETER :: hash_table_expand = 1.5_real_8, &
      inv_hash_table_fill = 2.5_real_8

    INTEGER                                  :: i, j
    TYPE(ele_type), ALLOCATABLE, &
      DIMENSION(:)                           :: tmp_hash

! if too small, make a copy and rehash in a larger table

    IF (hash_table%nele*inv_hash_table_fill>hash_table%nmax) THEN
       ALLOCATE(tmp_hash(LBOUND(hash_table%table,1):UBOUND(hash_table%table,1)))
       tmp_hash=hash_table%table
       CALL hash_table_release(hash_table)
       CALL hash_table_create(hash_table,INT((UBOUND(tmp_hash,1)+8)*hash_table_expand))
       DO i=LBOUND(tmp_hash,1),UBOUND(tmp_hash,1)
          IF (tmp_hash(i)%c.NE.0) THEN
             CALL hash_table_add(hash_table,tmp_hash(i)%c,tmp_hash(i)%p,error)
          ENDIF
       ENDDO
       DEALLOCATE(tmp_hash)
    ENDIF

   hash_table%nele=hash_table%nele+1
   i=IAND(c*hash_table%prime,hash_table%nmax)

   DO j=i,hash_table%nmax
      IF (hash_table%table(j)%c==0 .OR. hash_table%table(j)%c==c) THEN
         hash_table%table(j)%c=c
         hash_table%table(j)%p=p
         RETURN
      ENDIF
   ENDDO
   DO j=0,i-1
      IF (hash_table%table(j)%c==0 .OR. hash_table%table(j)%c==c) THEN
         hash_table%table(j)%c=c
         hash_table%table(j)%p=p
         RETURN
      ENDIF
   ENDDO
  END SUBROUTINE hash_table_add

  PURE FUNCTION hash_table_get(hash_table,c) RESULT(p)
    TYPE(hash_table_type), INTENT(IN)        :: hash_table
    INTEGER, INTENT(IN)                      :: c
    INTEGER                                  :: p

    INTEGER                                  :: i, j

   i=IAND(c*hash_table%prime,hash_table%nmax)

   ! catch the likely case first
   IF (hash_table%table(i)%c==c) THEN
      p=hash_table%table(i)%p
      RETURN
   ENDIF

   DO j=i,hash_table%nmax
      IF (hash_table%table(j)%c==0 .OR. hash_table%table(j)%c==c) THEN
         p=hash_table%table(j)%p
         RETURN
      ENDIF
   ENDDO
   DO j=0,i-1
      IF (hash_table%table(j)%c==0 .OR. hash_table%table(j)%c==c) THEN
         p=hash_table%table(j)%p
         RETURN
      ENDIF
   ENDDO
  END FUNCTION hash_table_get

! End of hashtable
! -----------------------------------------------------------------------------

END MODULE dbcsr_mm_csr
