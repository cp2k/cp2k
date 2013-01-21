!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Fourth layer of the dbcsr matrix-matrix multiplication.
!> 	   It hides the differences between performing calculations on the 
!>         GPU or the CPU.	
!>
!> \author  Urban Borstnik
!>	
!> <b>Modification history:</b>
!>  - 2010-02-23 Moved from dbcsr_operations
!>  - 2011-11    Moved parameter-stack processing routines to
!>               dbcsr_mm_methods.
!>  - 2013-01    reorganized code (Ole Schuett)
! *****************************************************************************


MODULE dbcsr_mm_driver

  USE dbcsr_block_operations,          ONLY: block_add
  USE dbcsr_config,                    ONLY: mm_driver,&
                                             mm_driver_cuda,&
                                             mm_stack_size,&
                                             mm_thread_workshare
  USE dbcsr_cuda_device,               ONLY: dbcsr_cuda_create_streams,&
                                             dbcsr_cuda_destroy_streams
  USE dbcsr_cuda_memory,               ONLY: dbcsr_cuda_dev_mem_alloc,&
                                             dbcsr_cuda_dev_mem_hold,&
                                             dbcsr_cuda_dev_mem_new,&
                                             dbcsr_cuda_dev_mem_release,&
                                             dbcsr_cuda_dev_mem_zero
  USE dbcsr_cuda_methods,              ONLY: dbcsr_cuda_dev_mem_get_alloc
  USE dbcsr_cuda_operations,           ONLY: dbcsr_cuda_cp_dev_to_host,&
                                             dbcsr_cuda_cp_host_to_dev
  USE dbcsr_cuda_types,                ONLY: dbcsr_cuda_mem_type
  USE dbcsr_data_methods,              ONLY: dbcsr_data_ensure_size,&
                                             dbcsr_data_get_size,&
                                             dbcsr_data_get_type,&
                                             dbcsr_data_init,&
                                             dbcsr_data_new,&
                                             dbcsr_data_release
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: int_4,&
                                             int_8,&
                                             real_8
  USE dbcsr_mm_types,                  ONLY: carrier_type
  USE dbcsr_pq_methods,                ONLY: dbcsr_ps_set_get_group_p,&
                                             dbcsr_ps_target_add_data_cuda,&
                                             dbcsr_psg_add_data_cuda_ab,&
                                             dbcsr_psg_rm_data_cuda_ab
  USE dbcsr_pq_types,                  ONLY: dbcsr_ps_group_type,&
                                             dbcsr_ps_state_empty,&
                                             dbcsr_ps_target_obj,&
                                             dbcsr_ps_target_type,&
                                             dbcsr_ps_width
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_type,&
                                             dbcsr_type_int_4

 !$ USE OMP_LIB
        
  IMPLICIT NONE
 
  PRIVATE
 
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_driver'
            
  PUBLIC :: dbcsr_mm_driver_outer_init, dbcsr_mm_driver_outer_finalize 
  PUBLIC :: dbcsr_mm_driver_inner_init
  
  CONTAINS
  
  
   
! *****************************************************************************
!> \brief Moved out of internal_operations
!>        This loads the data for matrizes A and B onto the card.        
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_driver_inner_init(left, right, carrier, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(carrier_type), INTENT(INOUT)        :: carrier
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_driver_inner_init', &
      routineP = moduleN//':'//routineN

    IF (mm_driver .EQ. mm_driver_cuda) THEN
       !$OMP MASTER
       CALL dbcsr_cuda_cp_host_to_dev(&            
            left%data_area,&
            carrier%param_sets%groups(1,1)%master%s%left_data_cuda,&
            async = .TRUE.,&
            error=error)
       CALL dbcsr_cuda_cp_host_to_dev(&
            right%data_area,&
            carrier%param_sets%groups(1,1)%master%s%right_data_cuda,&
            async = .TRUE.,&
            error=error)
       !$OMP END MASTER
    ENDIF
    
  END SUBROUTINE dbcsr_mm_driver_inner_init
  
! *****************************************************************************
!> \brief Moved out of dbcsr_multrec_init (actually csr_multiply_outer_init)
!>        This allocates the memory for the matrizes A and B on the card and
!>        and leaves pointers to them in every stack.         
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_driver_outer_init(left, right, product, product_target, carrier, fill_guess, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(dbcsr_ps_target_obj), INTENT(INOUT) :: product_target
    TYPE(carrier_type), INTENT(INOUT)        :: carrier
    REAL(KIND=real_8)                        :: fill_guess
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_driver_outer_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_type, nthreads, &
                                                ps_buffer, ps_memreg
    TYPE(dbcsr_cuda_mem_type), POINTER, SAVE :: a_dev, b_dev
    TYPE(dbcsr_ps_group_type), POINTER       :: ps_group

    nthreads = 1
    !$ nthreads = OMP_GET_NUM_THREADS ()
                
               
    ! first driver init part --------------------------------------------------- 

    IF (mm_driver .EQ. mm_driver_cuda) THEN
       !$OMP MASTER
       !
       ! Setup threads<->streams
       CALL dbcsr_cuda_create_streams (nthreads, error)
       !
       data_type = dbcsr_data_get_type (left%data_area)
       CALL dbcsr_cuda_dev_mem_new (a_dev, data_type, error=error)
       CALL dbcsr_cuda_dev_mem_new (b_dev, data_type, error=error)
       CALL dbcsr_cuda_dev_mem_alloc (a_dev,&
            dbcsr_data_get_size (left%data_area), error=error)
       CALL dbcsr_cuda_dev_mem_alloc (b_dev,&
            dbcsr_data_get_size (right%data_area), error=error)
       !$OMP END MASTER
       CALL init_card_c (product, product_target%t,fill_guess, error)
       !$OMP BARRIER
       IF (.FALSE.) THEN
          ! Left in for mental accounting
          !$OMP CRITICAL
          CALL dbcsr_cuda_dev_mem_hold (a_dev, error=error)
          CALL dbcsr_cuda_dev_mem_hold (b_dev, error=error)
          !$OMP END CRITICAL
       ENDIF
    ELSE
       !$OMP MASTER
       NULLIFY (a_dev, b_dev)
       !$OMP END MASTER
       !$OMP BARRIER
    ENDIF
    ! end of first part of driver init ----------------------------------------- 
    
    ! second part of driver init -----------------------------------------------
     IF (mm_driver .EQ. mm_driver_cuda) THEN
       DO ps_memreg = 1, carrier%param_sets%nmemregs
          DO ps_buffer = 1, carrier%param_sets%nbuffers
             ps_group => dbcsr_ps_set_get_group_p (carrier%param_sets,&
                  ps_buffer, ps_memreg, error=error)
             !$OMP CRITICAL (crit_data_card)
             CALL dbcsr_psg_add_data_cuda_ab (ps_group,&
                  a_dev, b_dev, error=error)
             !$OMP END CRITICAL (crit_data_card)
          ENDDO
       ENDDO
       IF (.FALSE.) THEN
          ! Left in for mental accounting
          !$OMP CRITICAL
          CALL dbcsr_cuda_dev_mem_release (a_dev, error=error)
          CALL dbcsr_cuda_dev_mem_release (b_dev, error=error)
          !$OMP END CRITICAL
       ENDIF
       !$OMP BARRIER
       !$OMP MASTER
       CALL dbcsr_cuda_dev_mem_release (a_dev, error=error)
       CALL dbcsr_cuda_dev_mem_release (b_dev, error=error)
       !$OMP END MASTER
    ENDIF
  
  END SUBROUTINE dbcsr_mm_driver_outer_init

   
! *****************************************************************************
!> \brief Moved out of dbcsr_multrec_finalize (actually csr_multiply_outer_finalize)
!>        Copies the result matrix C back from the card,
!>        deallocates the memory on the card, and destroys the streams.         
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_driver_outer_finalize(carrier, product, error)
    TYPE(carrier_type), INTENT(INOUT)        :: carrier
    TYPE(dbcsr_type), INTENT(inout)          :: product
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: data_size, data_type, &
                                                ithread, ps_buffer, ps_memreg
    TYPE(dbcsr_data_obj)                     :: tmp_data
    TYPE(dbcsr_ps_group_type), POINTER       :: ps_group

   IF (mm_driver .EQ. mm_driver_cuda) THEN
       ! this part is from dbcsr_internal_operations ---------------------------
                      
       ithread = 0
       !$ ithread = omp_get_thread_num()
       
       ! Transfer data from GPU to host.
       CALL dbcsr_data_init (tmp_data)
       CALL dbcsr_data_new (tmp_data,&
            data_type = dbcsr_data_get_type(product%wms(ithread+1)%data_area),&
            data_size = carrier%datasize)
       CALL dbcsr_cuda_cp_dev_to_host(&
            carrier%param_sets%groups(1,1)%master%s%t%t%product_data_cuda,&
            !product%wms(1)%data_area,&
            tmp_data,&
            error=error)
       CALL dbcsr_data_ensure_size(product%wms(ithread+1)%data_area,&
            carrier%datasize, error=error)
       CALL block_add (product%wms(ithread+1)%data_area, tmp_data,&
            len=carrier%datasize, error=error)
       CALL dbcsr_data_release (tmp_data)
       ! end of part from dbcsr_internal_operations ----------------------------
    
       ! this part is from dbcsr_multrec_finalize ------------------------------
       DO ps_memreg = 1, carrier%param_sets%nmemregs
          DO ps_buffer = 1, carrier%param_sets%nbuffers
             ps_group => dbcsr_ps_set_get_group_p (carrier%param_sets,&
                  ps_buffer, ps_memreg, error=error)
             !$OMP CRITICAL (crit_data_card)
             CALL dbcsr_psg_rm_data_cuda_ab (ps_group, error=error)
             !$OMP END CRITICAL (crit_data_card)
          ENDDO
       ENDDO
    ENDIF

    !$OMP BARRIER
    
    !$OMP MASTER
    IF (mm_driver .EQ. mm_driver_cuda) THEN
       CALL dbcsr_cuda_destroy_streams (error=error)
    ENDIF
    !$OMP END MASTER
    
  END SUBROUTINE dbcsr_mm_driver_outer_finalize
 
 
 
! *****************************************************************************
!> \brief Perform allocations and setup needed for multiplication on an
!>        accelerator card.
! *****************************************************************************
  SUBROUTINE init_card_c (host_matrix, product_target, fill_guess, error)
    TYPE(dbcsr_type), INTENT(IN)             :: host_matrix
    TYPE(dbcsr_ps_target_type), &
      INTENT(INOUT)                          :: product_target
    REAL(KIND=real_8)                        :: fill_guess
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_card_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: clocks_size, data_size, &
                                                data_type, error_handle, &
                                                ithread
    INTEGER, POINTER                         :: state_tmp
    INTEGER, TARGET                          :: state_tmp_tgt
    LOGICAL                                  :: do_master_work
    TYPE(dbcsr_cuda_mem_type), POINTER       :: c_dev, c_locks_dev, &
                                                params_dev_l
    TYPE(dbcsr_cuda_mem_type), POINTER, SAVE :: params_dev, stack_state_dev

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handle, error)

    ! Only the master thread does common work in when threads share work.
    do_master_work = .NOT. mm_thread_workshare
    do_master_work = .TRUE.
    !$OMP MASTER
    do_master_work = .TRUE.
    !$OMP END MASTER

    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM ()
    data_type = dbcsr_data_get_type (host_matrix%wms(ithread+1)%data_area)
    data_size = dbcsr_data_get_size (host_matrix%wms(ithread+1)%data_area)
    !
    ! Allocate space for product data on the card.
    CALL dbcsr_cuda_dev_mem_new (c_dev, data_type, error)
    CALL dbcsr_cuda_dev_mem_alloc (c_dev, data_size, error=error)
    !
    ! To finish allocations before the zeroing is launched.
    !
    !$OMP BARRIER
    !
    ! Zero C data on card.  This assumes that data will be
    ! summed with the host data after the multiplication is done.
    CALL dbcsr_cuda_dev_mem_zero(c_dev, first=1,&
         last=data_size, error=error)
    !
    ! Attach C data on card to the product target descriptor
    CALL dbcsr_ps_target_add_data_cuda (product_target, c_dev, error)
    CALL dbcsr_cuda_dev_mem_release (c_dev, error=error)
    !
    ! Allocate C locks on GPU.  Locks should be cleared as soon as
    ! allocated.
    CALL dbcsr_cuda_dev_mem_new (c_locks_dev,&
         dbcsr_type_int_4, error)
    clocks_size = INT(host_matrix%nblkcols_local,KIND=int_8)* &
                  INT(host_matrix%nblkrows_local,KIND=int_8)* &
                  fill_guess
    clocks_size = MAX(1,INT(clocks_size, KIND=int_4))
    CALL dbcsr_cuda_dev_mem_alloc (c_locks_dev,INT(clocks_size,KIND=int_4), error=error)
    CALL dbcsr_cuda_dev_mem_zero (c_locks_dev,&
         first=1,last=dbcsr_cuda_dev_mem_get_alloc(c_locks_dev),&
         error=error)
    product_target%c_locks_dev =>  c_locks_dev
    !
    ! Have a source from which to update states and setup the on-card
    ! parameter stacks.
    !
    !$OMP MASTER
    CALL dbcsr_cuda_dev_mem_new (stack_state_dev, dbcsr_type_int_4,&
         error=error)
    CALL dbcsr_cuda_dev_mem_alloc (stack_state_dev, 1, error=error)
    state_tmp_tgt = dbcsr_ps_state_empty
    state_tmp => state_tmp_tgt
    CALL dbcsr_cuda_cp_host_to_dev (state_tmp, stack_state_dev%d_i, error=error)
    !$OMP END MASTER
    IF (do_master_work) THEN
       ! Allocate the on-card parameter stack.
       CALL dbcsr_cuda_dev_mem_new (params_dev_l, dbcsr_type_int_4, error)
       CALL dbcsr_cuda_dev_mem_alloc (params_dev_l, mm_stack_size*dbcsr_ps_width,&
            error=error)
       !$OMP MASTER
       params_dev => params_dev_l
       !$OMP END MASTER
    ENDIF
    !$OMP BARRIER
    !$OMP CRITICAL (crit_data_card)
    !
    IF (.FALSE.) THEN
       ! Left in for mental accounting
       CALL dbcsr_cuda_dev_mem_hold (stack_state_dev, error=error)
    ENDIF
    ! all threads share the master's stack_state_dev
    product_target%stack_state_dev => stack_state_dev
    IF (do_master_work) THEN
       product_target%params_dev => params_dev_l
    ELSE
       product_target%params_dev => params_dev
    ENDIF
    CALL dbcsr_cuda_dev_mem_hold (product_target%params_dev, error=error)
    IF (.FALSE.) THEN
       ! Left in for mental accounting
       CALL dbcsr_cuda_dev_mem_release (stack_state_dev, error=error)
    ENDIF
    !
    !$OMP END CRITICAL (crit_data_card)
    !$OMP BARRIER
    IF (do_master_work) THEN
       !$OMP CRITICAL (crit_data_card)
       CALL dbcsr_cuda_dev_mem_release (params_dev_l, error=error)
       !$OMP END CRITICAL (crit_data_card)
    ENDIF
    !
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE init_card_c
  
  
END MODULE dbcsr_mm_driver
