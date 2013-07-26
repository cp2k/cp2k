!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Third layer of the dbcsr matrix-matrix multiplication.
!>         It collects the full matrix blocks, which need to be multiplied,
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
!>  - 2013-01    extensive refactoring (Ole Schuett)
! *****************************************************************************

MODULE dbcsr_mm_csr

  USE array_types,                     ONLY: array_data
  USE dbcsr_config,                    ONLY: dbcsr_get_conf_nstacks,&
                                             mm_driver,&
                                             mm_host_driver
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: int_1,&
                                             int_4,&
                                             int_8,&
                                             real_8,&
                                             sp
  USE dbcsr_mm_driver,                 ONLY: dbcsr_mm_driver_inner_finalize,&
                                             dbcsr_mm_driver_inner_init,&
                                             dbcsr_mm_driver_lib_finalize,&
                                             dbcsr_mm_driver_lib_init,&
                                             dbcsr_mm_driver_outer_finalize,&
                                             dbcsr_mm_driver_outer_init,&
                                             dbcsr_mm_driver_process_stack
  USE dbcsr_mm_types,                  ONLY: &
       csr_carrier_type, dbcsr_ps_state_empty, dbcsr_ps_state_working, &
       ele_type, hash_table_type, p_a_first, p_b_first, p_c_blk, p_c_first, &
       p_k, p_m, p_n, stack_descriptor_type
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
  USE dbcsr_toollib,                   ONLY: sort
  USE dbcsr_types,                     ONLY: dbcsr_type
  USE dbcsr_util,                      ONLY: map_most_common

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_csr'
  LOGICAL, PARAMETER :: debug_mod  = .FALSE.
  LOGICAL, PARAMETER :: careful_mod = .FALSE.
  REAL, PARAMETER                      :: default_resize_factor = 1.618034

  !> \var max_stack_block_size  The maximal block size to be specially
  !>                            treated.
  INTEGER, PARAMETER :: max_stack_block_size = HUGE (INT (0))


  PUBLIC :: dbcsr_mm_csr_lib_init,   dbcsr_mm_csr_lib_finalize
  PUBLIC :: dbcsr_mm_csr_outer_init, dbcsr_mm_csr_outer_finalize
  PUBLIC :: dbcsr_mm_csr_inner_init, dbcsr_mm_csr_inner_finalize
  PUBLIC :: dbcsr_mm_csr_multiply


  CONTAINS

! *****************************************************************************
! \brief Initialize the library
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_csr_lib_init(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

     CALL dbcsr_mm_driver_lib_init(error)
  END SUBROUTINE


! *****************************************************************************
! \brief Finalize the library
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_csr_lib_finalize(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

     CALL dbcsr_mm_driver_lib_finalize(error)
  END SUBROUTINE


! *****************************************************************************
! \brief A wrapper around dbcsr_mm_csr_multiply_low to avoid expensive dereferencings.
! \author Ole Schuett
! *****************************************************************************
 SUBROUTINE dbcsr_mm_csr_multiply(mi, mf, ni, nf, ki, kf,&
       ai, af,&
       bi, bf,&
       m_sizes, n_sizes, k_sizes,&
       c_local_rows, c_local_cols,&
       c_has_symmetry, keep_sparsity, use_eps,&
       row_max_epss,&
       flop,&
       a_index, b_index, a_norms, b_norms,&
       product, csr_carrier, error)
    INTEGER, INTENT(IN)                      :: mi, mf, ni, nf, ki, kf, ai, &
                                                af, bi, bf
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: m_sizes, n_sizes, k_sizes, &
                                                c_local_rows, c_local_cols
    LOGICAL, INTENT(INOUT)                   :: c_has_symmetry, &
                                                keep_sparsity, use_eps
    REAL(kind=sp), DIMENSION(:)              :: row_max_epss
    INTEGER(KIND=int_8), INTENT(INOUT)       :: flop
    INTEGER, DIMENSION(1:3, 1:af), &
      INTENT(IN)                             :: a_index
    INTEGER, DIMENSION(1:3, 1:bf), &
      INTENT(IN)                             :: b_index
    REAL(KIND=sp), DIMENSION(:), POINTER     :: a_norms, b_norms
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(csr_carrier_type), INTENT(INOUT)    :: csr_carrier
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: ithread, max_new_nblks, &
                                                nblks_new

    ithread = 0
    !$ ithread = omp_get_thread_num()

    ! This has to be done here because ensure_array_size() expects a pointer.
    ! the maximum number of blocks can be safely estimated by considering both the rowxcol,
    ! but also the blocks the latter can never be larger than norec**2, which is a 'small' constant
    max_new_nblks = MIN(INT(mf-mi+1,int_8) * INT(nf-ni+1,int_8), &
                        INT(af-ai+1,int_8) * INT(bf-bi+1,int_8))

    nblks_new = product%wms(ithread+1)%lastblk + max_new_nblks

    CALL ensure_array_size(product%wms(ithread+1)%row_i, ub=nblks_new,&
         factor=default_resize_factor, error=error)
    CALL ensure_array_size(product%wms(ithread+1)%col_i, ub=nblks_new,&
         factor=default_resize_factor, error=error)
    CALL ensure_array_size(product%wms(ithread+1)%blk_p, ub=nblks_new,&
         factor=default_resize_factor, error=error)

    CALL dbcsr_mm_csr_multiply_low(mi=mi, mf=mf,ni=ni, nf=nf, ki=ki, kf=kf,&
            ai=ai, af=af,&
            bi=bi, bf=bf,&
            c_row_i=product%wms(ithread+1)%row_i,&
            c_col_i=product%wms(ithread+1)%col_i,&
            c_blk_p=product%wms(ithread+1)%blk_p,&
            lastblk =product%wms(ithread+1)%lastblk, &
            datasize=product%wms(ithread+1)%datasize,&
            m_sizes=m_sizes, n_sizes=n_sizes, k_sizes=k_sizes,&
            c_local_rows=c_local_rows, c_local_cols=c_local_cols,&
            c_has_symmetry=c_has_symmetry, keep_sparsity=keep_sparsity,&
            use_eps=use_eps,&
            row_max_epss=row_max_epss,&
            flop=flop,&
            row_size_maps=csr_carrier%m_size_maps,& 
            col_size_maps=csr_carrier%n_size_maps,& 
            k_size_maps=csr_carrier%k_size_maps,&
            row_size_maps_size=csr_carrier%m_size_maps_size,& 
            col_size_maps_size=csr_carrier%n_size_maps_size,& 
            k_size_maps_size=csr_carrier%k_size_maps_size,&
            nm_stacks=csr_carrier%nm_stacks, nn_stacks=csr_carrier%nn_stacks,&
            nk_stacks=csr_carrier%nk_stacks, &
            stack_map=csr_carrier%stack_map,&
            c_hashes=csr_carrier%c_hashes,&
            a_index=a_index, b_index=b_index,&
            a_norms=a_norms, b_norms=b_norms,&
            csr_carrier=csr_carrier,&
            error=error)

   END SUBROUTINE dbcsr_mm_csr_multiply


! *****************************************************************************
!> \brief Performs multiplication of smaller submatrices.
! *****************************************************************************
  SUBROUTINE dbcsr_mm_csr_multiply_low(mi, mf, ni, nf, ki, kf,&
       ai, af,&
       bi, bf,&
       c_hashes, c_row_i, c_col_i, c_blk_p, lastblk, datasize,&
       m_sizes, n_sizes, k_sizes,&
       c_local_rows, c_local_cols,&
       c_has_symmetry, keep_sparsity, use_eps,&
       row_size_maps, col_size_maps, k_size_maps,&
       row_size_maps_size, col_size_maps_size, k_size_maps_size,&
       nm_stacks, nn_stacks, nk_stacks, stack_map,&
       row_max_epss, &
       flop, a_index, b_index,a_norms, b_norms, csr_carrier, error)
    INTEGER, INTENT(IN)                      :: mi, mf, ni, nf, ki, kf, ai, &
                                                af, bi, bf
    TYPE(hash_table_type), DIMENSION(:), &
      INTENT(INOUT)                          :: c_hashes
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: c_row_i, c_col_i, c_blk_p
    INTEGER, INTENT(INOUT)                   :: lastblk, datasize
    INTEGER, DIMENSION(:), INTENT(IN)        :: m_sizes, n_sizes, k_sizes, &
                                                c_local_rows, c_local_cols
    LOGICAL, INTENT(IN)                      :: c_has_symmetry, &
                                                keep_sparsity, use_eps
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
    REAL(kind=sp), DIMENSION(:)              :: row_max_epss
    INTEGER(KIND=int_8), INTENT(INOUT)       :: flop
    INTEGER, DIMENSION(1:3, 1:af), &
      INTENT(IN)                             :: a_index
    INTEGER, DIMENSION(1:3, 1:bf), &
      INTENT(IN)                             :: b_index
    REAL(KIND=sp), DIMENSION(:), POINTER     :: a_norms, b_norms
    TYPE(csr_carrier_type), INTENT(INOUT)    :: csr_carrier
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_csr_multiply_low', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE., &
                                                local_timing = .FALSE.

    INTEGER :: a_blk, a_col_l, a_row_l, b_blk, b_col_l, c, c_blk_id, &
      c_col_logical, c_nze, c_row_logical, ithread, k_size, m_size, &
      mapped_col_size, mapped_k_size, mapped_row_size, n_a_norms, n_b_norms, &
      n_size, nstacks, s_dp, ws
    INTEGER, DIMENSION(mi:mf+1)              :: a_row_p
    INTEGER, DIMENSION(ki:kf+1)              :: b_row_p
    INTEGER, DIMENSION(:, :, :), POINTER     :: stackcol_data
    INTEGER, DIMENSION(:), POINTER           :: stackcol_fillcount
    INTEGER, DIMENSION(2, bf-bi+1)           :: b_blk_info
    INTEGER, DIMENSION(2, af-ai+1)           :: a_blk_info
    INTEGER(KIND=int_4)                      :: offset
    LOGICAL                                  :: block_exists, flush_stack
    REAL(kind=sp)                            :: a_norm, a_row_eps, b_norm
    REAL(KIND=sp), DIMENSION(1:af-ai+1)      :: left_norms
    REAL(KIND=sp), DIMENSION(1:bf-bi+1)      :: right_norms

!   ---------------------------------------------------------------------------

    ithread = 0
    !$ ithread = omp_get_thread_num()

    nstacks = SIZE(csr_carrier%stackgrid_data, 3)

    IF (use_eps) THEN
       n_a_norms = af-ai+1
       n_b_norms = bf-bi+1
    ELSE
       n_a_norms = 0
       n_b_norms = 0
    ENDIF


    !
    ! Build the indices
    CALL build_csr_index (mi,mf,ai,af,a_row_p, a_blk_info, a_index,&
         n_a_norms, left_norms, a_norms)
    CALL build_csr_index (ki,kf,bi,bf,b_row_p, b_blk_info, b_index,&
         n_b_norms, right_norms, b_norms)


    CALL get_free_stack_column(csr_carrier, stackcol_data, stackcol_fillcount, error)


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
             IF (block_exists) THEN
                offset = c_blk_p(c_blk_id)    
             ELSE
                IF (keep_sparsity) CYCLE

                offset = datasize + 1
                lastblk = lastblk+1
                datasize = datasize + c_nze
                c_blk_id = lastblk ! assign a new c-block-id

                IF (dbg) WRITE(*,*)routineN//" new block offset, nze", offset, c_nze
                CALL hash_table_add(c_hashes(a_row_l),&
                     b_col_l, c_blk_id, error=error)

                ! We still keep the linear index because it's
                ! easier than getting the values out of the
                ! hashtable in the end.
                c_row_i(lastblk) = a_row_l
                c_col_i(lastblk) = b_col_l
                c_blk_p(lastblk) = offset
             ENDIF

             ! TODO: this is only called withcareful_mod 
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

             mapped_col_size = col_size_maps (n_size)
             ws = stack_map (mapped_col_size, mapped_k_size, mapped_row_size)
             stackcol_fillcount(ws) = stackcol_fillcount(ws) + 1
             s_dp = stackcol_fillcount(ws)

             stackcol_data(p_m, s_dp, ws) = m_size
             stackcol_data(p_n, s_dp, ws) = n_size
             stackcol_data(p_k, s_dp, ws) = k_size
             stackcol_data(p_a_first, s_dp, ws) = a_blk_info(2, a_blk)
             stackcol_data(p_b_first, s_dp, ws) = b_blk_info(2, b_blk)
             stackcol_data(p_c_first, s_dp, ws) = offset
             stackcol_data(p_c_blk,   s_dp, ws) = c_blk_id

             flop = flop + INT(2*c_nze, int_8) * INT(k_size, int_8)

             flush_stack = stackcol_fillcount(ws) >= SIZE(stackcol_data, 2)
             IF (flush_stack) &
                      CALL get_free_stack_column(csr_carrier, stackcol_data, stackcol_fillcount, error)

          ENDDO b_blk_cycle ! b
       ENDDO a_blk_cycle ! a_col
    ENDDO a_row_cycle ! a_row


  END SUBROUTINE dbcsr_mm_csr_multiply_low

! *****************************************************************************
!> \brief Initializes a multiplication cycle for new set of C-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_csr_outer_init(csr_carrier, left, right, product, fill_guess,&
    m_sizes, n_sizes, k_sizes, error)
    TYPE(csr_carrier_type), INTENT(INOUT)    :: csr_carrier
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    REAL(KIND=real_8)                        :: fill_guess
    INTEGER, DIMENSION(:), POINTER           :: m_sizes, n_sizes, k_sizes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_csr_outer_init', &
      routineP = moduleN//':'//routineN

    INTEGER :: block_estimate, default_stack, istack, ithread, k_map, k_size, &
      m_map, m_size, n_map, n_size, n_stack_buffers, nstacks, nthreads, ps_g
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: flop_index, flop_list, &
                                                most_common_k, most_common_m, &
                                                most_common_n
    INTEGER, DIMENSION(3)                    :: nxstacks
    LOGICAL                                  :: local_indexing
    TYPE(stack_descriptor_type), &
      ALLOCATABLE, DIMENSION(:)              :: tmp_descr

    ithread = 0 ; nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM () ; nthreads = OMP_GET_NUM_THREADS ()

    ! find out if we have local_indexing
    CALL dbcsr_assert (right%local_indexing, "EQV", left%local_indexing,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Local index useage must be consistent.", __LINE__, error=error)
    local_indexing = left%local_indexing

    ! Setup the hash tables if needed
    block_estimate=MAX(product%nblks,left%nblks,right%nblks)/nthreads
    IF (local_indexing) THEN
       ALLOCATE (csr_carrier%c_hashes (product%nblkrows_local))
       CALL fill_hash_tables (csr_carrier%c_hashes, product,block_estimate,&
            row_map=array_data(product%global_rows),&
            col_map=array_data(product%global_cols),&
            error=error)
    ELSE
       ALLOCATE (csr_carrier%c_hashes (product%nblkrows_total))
       CALL fill_hash_tables (csr_carrier%c_hashes, product,block_estimate,&
            error=error)
    ENDIF

    ! Setup the MM stack
    CALL dbcsr_get_conf_nstacks (nxstacks, n_stack_buffers, error)
    csr_carrier%nm_stacks = nxstacks(1)
    csr_carrier%nn_stacks = nxstacks(2)
    csr_carrier%nk_stacks = nxstacks(3)
    nstacks = nxstacks(1) * nxstacks(2) * nxstacks(3) + 1
    CALL dbcsr_assert (nstacks, "LE", INT (HUGE (csr_carrier%stack_map)),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Too many stacks requested (global/dbcsr/n_size_*_stacks in input)",&
         __LINE__, error=error)


    ALLOCATE(csr_carrier%stackrow_descr(nstacks))

    DO istack=1, nstacks
       csr_carrier%stackrow_descr(istack)%total_entries = 0
    END DO

    ALLOCATE (most_common_m(nxstacks(1)))
    ALLOCATE (most_common_n(nxstacks(2)))
    ALLOCATE (most_common_k(nxstacks(3)))
    CALL map_most_common (m_sizes, csr_carrier%m_size_maps, nxstacks(1),&
         most_common_m,&
         max_stack_block_size, csr_carrier%max_m)
    csr_carrier%m_size_maps_size = SIZE (csr_carrier%m_size_maps)
    CALL map_most_common (n_sizes, csr_carrier%n_size_maps, nxstacks(2),&
         most_common_n,&
         max_stack_block_size, csr_carrier%max_n)
    csr_carrier%n_size_maps_size = SIZE (csr_carrier%n_size_maps)
    CALL map_most_common (array_data(right%row_blk_size),&
         csr_carrier%k_size_maps, nxstacks(3), &
         most_common_k,&
         max_stack_block_size, csr_carrier%max_k)
    csr_carrier%k_size_maps_size = SIZE (csr_carrier%k_size_maps)

    ! Creates the stack map--a mapping from (mapped) stack block sizes
    ! (carrier%*_sizes) to a stack number.  Triples with even one
    ! uncommon size will be mapped to a general, non-size-specific
    ! stack.
    ALLOCATE (csr_carrier%stack_map(nxstacks(2)+1, nxstacks(3)+1, nxstacks(1)+1))
    default_stack = nstacks

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
                       (k_map-1)*nxstacks(2) + n_map
                ps_g = nstacks-ps_g
                csr_carrier%stack_map(n_map, k_map, m_map) = ps_g
                ! Also take care of the stack m, n, k descriptors
                csr_carrier%stackrow_descr(ps_g)%m     = m_size
                csr_carrier%stackrow_descr(ps_g)%n     = n_size
                csr_carrier%stackrow_descr(ps_g)%k     = k_size
                csr_carrier%stackrow_descr(ps_g)%max_m = m_size
                csr_carrier%stackrow_descr(ps_g)%max_n = n_size
                csr_carrier%stackrow_descr(ps_g)%max_k = k_size
                csr_carrier%stackrow_descr(ps_g)%defined_mnk = .TRUE.
             ELSE
                ! This is the case when at least one of m, n, or k is
                ! undefined.
                ps_g = default_stack
                csr_carrier%stack_map(n_map, k_map, m_map) = default_stack
                ! Also take care of the stack m, n, k descriptors
                csr_carrier%stackrow_descr(ps_g)%m     = 0
                csr_carrier%stackrow_descr(ps_g)%n     = 0
                csr_carrier%stackrow_descr(ps_g)%k     = 0
                csr_carrier%stackrow_descr(ps_g)%max_m = csr_carrier%max_m
                csr_carrier%stackrow_descr(ps_g)%max_n = csr_carrier%max_n
                csr_carrier%stackrow_descr(ps_g)%max_k = csr_carrier%max_k
                csr_carrier%stackrow_descr(ps_g)%defined_mnk = .FALSE.
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
    ! The default stack shall remain at the end of the gridcolumn
    ALLOCATE(flop_list(nstacks-1),flop_index(nstacks-1), tmp_descr(nstacks))
    DO istack=1,nstacks-1
       flop_list(istack) = -2 * csr_carrier%stackrow_descr(istack)%m& 
                              * csr_carrier%stackrow_descr(istack)%n&
                              * csr_carrier%stackrow_descr(istack)%k
    ENDDO

    CALL sort(flop_list, nstacks-1, flop_index)
    tmp_descr = csr_carrier%stackrow_descr
    DO istack=1,nstacks-1
       csr_carrier%stackrow_descr(istack) = tmp_descr(flop_index(istack)) 
    ENDDO

    DO m_map = 1, SIZE(csr_carrier%stack_map, 1)
      DO k_map = 1, SIZE(csr_carrier%stack_map, 2)
        map_loop: DO n_map = 1, SIZE(csr_carrier%stack_map, 1)
          DO istack=1,nstacks-1
            IF(csr_carrier%stack_map(m_map, k_map, n_map) == flop_index(istack)) THEN
               csr_carrier%stack_map(m_map, k_map, n_map) = istack
               CYCLE map_loop
            END IF
          ENDDO
        ENDDO map_loop
      ENDDO
    ENDDO
    DEALLOCATE(flop_list,flop_index,tmp_descr)




    ! Driver allocates stackgrid_data and stackgrid_state, if needed as host-pinned.
    NULLIFY(csr_carrier%stackgrid_data, csr_carrier%stackgrid_state)
    CALL dbcsr_mm_driver_outer_init(left=left, right=right, product=product,&
           drv_carrier=csr_carrier%drv_carrier, fill_guess=fill_guess,&
           stackgrid_data=csr_carrier%stackgrid_data, &
           stackgrid_state=csr_carrier%stackgrid_state, &
           nstacks=nstacks, n_stack_buffers=n_stack_buffers, error=error)

     IF(.NOT. ASSOCIATED(csr_carrier%stackgrid_data)) &
        STOP "dbcsr_mm_csr_outer_init: stackgrid_data not associated"

     IF(.NOT. ASSOCIATED(csr_carrier%stackgrid_state)) &
        STOP "dbcsr_mm_csr_outer_init: stackgrid_state not associated"

    ALLOCATE(csr_carrier%stackgrid_fillcount(nstacks, n_stack_buffers))
    csr_carrier%stackgrid_fillcount = 0
    csr_carrier%stackgrid_state = dbcsr_ps_state_empty


  END SUBROUTINE dbcsr_mm_csr_outer_init


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
!> \brief Finalizes a multiplication cycle for a set of C-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_csr_outer_finalize(csr_carrier, product, error)
    TYPE(csr_carrier_type), INTENT(INOUT)    :: csr_carrier
    TYPE(dbcsr_type), INTENT(inout)          :: product
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_csr_outer_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

! The driver deallocates stackgrid_data and stackgrid_state as it might be host-pinned.

    CALL dbcsr_mm_driver_outer_finalize(drv_carrier=csr_carrier%drv_carrier,&
       product=product, stackgrid_data=csr_carrier%stackgrid_data,&
       stackgrid_state=csr_carrier%stackgrid_state,&
       error=error)

    IF(ASSOCIATED(csr_carrier%stackgrid_data)) & 
        STOP "dbcsr_mm_csr_outer_finalize: stackgrid_data still associated"
    IF(ASSOCIATED(csr_carrier%stackgrid_state)) & 
        STOP "dbcsr_mm_csr_outer_finalize: stackgrid_state still associated"

    ! Clear hash tables
    DO i = 1, SIZE(csr_carrier%c_hashes)
       CALL hash_table_release (csr_carrier%c_hashes(i))
    ENDDO
    DEALLOCATE (csr_carrier%c_hashes)

    DEALLOCATE(csr_carrier%stackrow_descr)
    DEALLOCATE(csr_carrier%stackgrid_fillcount)
    DEALLOCATE(csr_carrier%stack_map)
    DEALLOCATE(csr_carrier%m_size_maps)
    DEALLOCATE(csr_carrier%n_size_maps)
    DEALLOCATE(csr_carrier%k_size_maps)

  END SUBROUTINE dbcsr_mm_csr_outer_finalize


! *****************************************************************************
!> \brief Initializes a multiplication cycle for a new set of A- and B-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_csr_inner_init(left, right, csr_carrier, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(csr_carrier_type), INTENT(INOUT)    :: csr_carrier
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CALL dbcsr_mm_driver_inner_init(left, right, csr_carrier%drv_carrier, error)
  END SUBROUTINE dbcsr_mm_csr_inner_init

! *****************************************************************************
!> \brief Finalizes a multiplication cycle for a set of A- and B-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_csr_inner_finalize(csr_carrier, error)
    TYPE(csr_carrier_type), INTENT(INOUT)    :: csr_carrier
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

! submit all left over stacks

    CALL process_stacks(csr_carrier, driver=mm_driver, min_fill=0, error=error)
    CALL dbcsr_mm_driver_inner_finalize(csr_carrier%drv_carrier, error=error)
  END SUBROUTINE dbcsr_mm_csr_inner_finalize


! *****************************************************************************
!> \brief Processes stacks in the stackgrid which have min_fill or more entries.
!> \param driver      Which driver (smm, blas, cuda,...) should be used
!> \param min_fill    Only process stacks that contain at least min_fill items.
!> \param max_stacks  Process only given number of stacks. Not present = no limit.
!> \param start_below If set, the stackgrid will be search bottom-up.
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE process_stacks(csr_carrier, driver, min_fill, max_stacks, start_below, error)
    TYPE(csr_carrier_type), INTENT(INOUT)    :: csr_carrier
    INTEGER, INTENT(IN)                      :: driver, min_fill
    INTEGER, INTENT(IN), OPTIONAL            :: max_stacks
    LOGICAL, INTENT(IN), OPTIONAL            :: start_below
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i, j, j_end, j_incr, j_start, &
                                                n_stack_buffers, nstacks, &
                                                num_processed
    INTEGER, DIMENSION(:, :), POINTER        :: stack_data
    INTEGER, POINTER                         :: stack_fillcount, stack_state
    TYPE(stack_descriptor_type), POINTER     :: stack_descr

    num_processed = 0
    nstacks = SIZE(csr_carrier%stackgrid_data, 3)
    n_stack_buffers = SIZE(csr_carrier%stackgrid_data, 4)

    j_start=1; j_end=nstacks; j_incr=1
    IF(PRESENT(start_below)) THEN
       IF(start_below) THEN
          j_start=nstacks; j_end=1; j_incr=-1
       ENDIF
    ENDIF

    DO i=1, n_stack_buffers
      DO j=j_start, j_end, j_incr
         IF(csr_carrier%stackgrid_fillcount(j,i) > min_fill) THEN
            stack_data=>csr_carrier%stackgrid_data(:,:,j,i)
            stack_fillcount=>csr_carrier%stackgrid_fillcount(j,i)
            stack_state=>csr_carrier%stackgrid_state(j,i)
            stack_descr=>csr_carrier%stackrow_descr(j)

            CALL dbcsr_mm_driver_process_stack(&
                       drv_carrier=csr_carrier%drv_carrier,&
                       driver=driver,&
                       stack_data=stack_data,&
                       stack_fillcount=stack_fillcount,&
                       stack_state=stack_state,&
                       stack_descr=stack_descr,&
                       stack_id=j + (i-1)*nstacks,&
                       error=error)

            num_processed = num_processed + 1
            IF(PRESENT(max_stacks)) THEN
               IF(num_processed >= max_stacks) RETURN
            END IF
         END IF
      END DO
    END DO

END SUBROUTINE process_stacks


! *****************************************************************************
!> \brief Returns a free column of the stackgrid.
!>        If also processes some of the stacks in the stackgrid.
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE get_free_stack_column(csr_carrier, stackcol_data, stackcol_fillcount, error)
    TYPE(csr_carrier_type), INTENT(INOUT)    :: csr_carrier
    INTEGER, DIMENSION(:, :, :), POINTER     :: stackcol_data
    INTEGER, DIMENSION(:), POINTER           :: stackcol_fillcount
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i, j, mm_stack_size, &
                                                n_stack_buffers, nstacks

    n_stack_buffers = SIZE(csr_carrier%stackgrid_data, 4)
    nstacks = SIZE(csr_carrier%stackgrid_data, 3)
    mm_stack_size = SIZE(csr_carrier%stackgrid_data, 2)

    CALL process_stacks(csr_carrier, driver=mm_driver, min_fill=mm_stack_size/2, error=error)

    ! checking for unkown stack states
    DO i=1, n_stack_buffers
      DO j=1, nstacks
        IF((csr_carrier%stackgrid_state(j,i) /= dbcsr_ps_state_empty) &
            .AND. (csr_carrier%stackgrid_state(j,i) /= dbcsr_ps_state_working)) THEN
            WRITE (*,*) "odd stack state: ", csr_carrier%stackgrid_state(j,i)
            STOP
        END IF
      END DO
    END DO

    ! looking for a free stackgrid column
    outer_loop: DO WHILE(.TRUE.)
       DO i=1, n_stack_buffers
          IF(ALL(csr_carrier%stackgrid_state(:,i) == dbcsr_ps_state_empty)) THEN
             stackcol_data => csr_carrier%stackgrid_data(:,:,:,i)
             stackcol_fillcount => csr_carrier%stackgrid_fillcount(:,i)
             EXIT outer_loop
          END IF
       END DO

       ! Waiting for free stackgrid column, doing some cpu work in the meanwhile.
       CALL process_stacks(csr_carrier, driver=mm_host_driver, &
              min_fill=0, max_stacks=1, start_below=.TRUE., error=error)
    END DO outer_loop


  END SUBROUTINE get_free_stack_column


! *****************************************************************************
!> \brief Prints a nice tabular summary of the stackgrid's current state.
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE print_stackgrid(csr_carrier)
    TYPE(csr_carrier_type), INTENT(INOUT)    :: csr_carrier

    INTEGER                                  :: i, j, n_stack_buffers, nstacks

    n_stack_buffers = SIZE(csr_carrier%stackgrid_data, 4)
    nstacks = SIZE(csr_carrier%stackgrid_data, 3)

    ! print table header
    WRITE(*,"(14X,A)",advance='no') "|     total "
    DO i=1, n_stack_buffers
       WRITE(*,"(A,I3,A)",advance='no') "|", i, ". grp "
    END DO
    WRITE(*,"(A)",advance='yes') "|"

    ! print separation line below table header
    WRITE(*,"(A)",advance='no') "--------------+-----------+"
    DO j=1, n_stack_buffers
       WRITE(*,"(A)",advance='no') "---------+"
    END DO
    WRITE(*,"(A)",advance='yes') ""

    DO i=1, nstacks
      IF(csr_carrier%stackrow_descr(i)%defined_mnk) THEN
                  WRITE(*,"(A,I2,A,I2,A,I2,A)",advance='no') "|", & 
                      csr_carrier%stackrow_descr(i)%m, " x ", &
                      csr_carrier%stackrow_descr(i)%n, " x ", &
                      csr_carrier%stackrow_descr(i)%k, " "
       ELSE
                  WRITE(*,"(A)",advance='no') "|       mixed "
       END IF

       WRITE(*,"(A,I9,X)",advance='no') "| ", csr_carrier%stackrow_descr(i)%total_entries
       DO j=1, n_stack_buffers
          IF(csr_carrier%stackgrid_state(i,j)==dbcsr_ps_state_empty) THEN 
             WRITE(*,"(A,I7,A)",advance='no') "| ",csr_carrier%stackgrid_fillcount(i,j)," "
          ELSE
             WRITE(*,"(A,I1,A)",advance='no') "| busy(",csr_carrier%stackgrid_state(i,j),") "
          END IF
       END DO
       WRITE(*,"(A)",advance='yes') "|"
    END DO



END SUBROUTINE print_stackgrid



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
