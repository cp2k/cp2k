!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
                                             default_resize_factor,&
                                             mm_stack_size
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: int_1,&
                                             int_4,&
                                             int_8,&
                                             real_8,&
                                             sp
  USE dbcsr_mm_sched,                  ONLY: &
       dbcsr_mm_sched_barrier, dbcsr_mm_sched_begin_burst, &
       dbcsr_mm_sched_end_burst, dbcsr_mm_sched_finalize, &
       dbcsr_mm_sched_init, dbcsr_mm_sched_lib_finalize, &
       dbcsr_mm_sched_lib_init, dbcsr_mm_sched_phaseout, &
       dbcsr_mm_sched_process, dbcsr_mm_sched_type
  USE dbcsr_mm_types,                  ONLY: &
       dbcsr_ps_width, p_a_first, p_b_first, p_c_blk, p_c_first, p_k, p_m, &
       p_n, stack_descriptor_type
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
  USE dbcsr_toollib,                   ONLY: sort
  USE dbcsr_types,                     ONLY: dbcsr_type,&
                                             dbcsr_work_type
  USE dbcsr_util,                      ONLY: map_most_common

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_csr'
  LOGICAL, PARAMETER :: debug_mod  = .FALSE.
  LOGICAL, PARAMETER :: careful_mod = .FALSE.

  !> \var max_stack_block_size  The maximal block size to be specially
  !>                            treated.
  INTEGER, PARAMETER :: max_stack_block_size = HUGE (INT (0))


 ! *****************************************************************************
  TYPE dbcsr_mm_csr_type
     PRIVATE
     TYPE(hash_table_type), DIMENSION(:), POINTER  :: c_hashes => Null()
     INTEGER                        :: nm_stacks, nn_stacks, nk_stacks
     INTEGER(KIND=int_4), DIMENSION(:), POINTER :: m_size_maps => Null()
     INTEGER(KIND=int_4), DIMENSION(:), POINTER :: n_size_maps => Null()
     INTEGER(KIND=int_4), DIMENSION(:), POINTER :: k_size_maps => Null()
     INTEGER                        :: max_m, max_n, max_k
     INTEGER                        :: m_size_maps_size,&
                                       n_size_maps_size,&
                                       k_size_maps_size
     INTEGER(KIND=int_1), DIMENSION(:,:,:), POINTER :: stack_map => Null()
     TYPE(stack_descriptor_type), DIMENSION(:), POINTER  :: stacks_descr => Null()
     TYPE(dbcsr_work_type), POINTER           :: product_wm => Null()
     INTEGER, DIMENSION(:,:,:), POINTER       :: stacks_data  => Null()
     INTEGER, DIMENSION(:), POINTER           :: stacks_fillcount => Null()
     TYPE(dbcsr_mm_sched_type)                      :: sched
  END TYPE dbcsr_mm_csr_type


! *****************************************************************************
!> \brief Types needed for the hashtable.
! *****************************************************************************
  TYPE ele_type
     INTEGER :: c=0
     INTEGER :: p=0
  END TYPE ele_type

  TYPE hash_table_type
     TYPE(ele_type), DIMENSION(:), POINTER :: table
     INTEGER :: nele=0
     INTEGER :: nmax=0
     INTEGER :: prime=0
  END TYPE hash_table_type

 ! *****************************************************************************
  PUBLIC :: dbcsr_mm_csr_type
  PUBLIC :: dbcsr_mm_csr_lib_init,   dbcsr_mm_csr_lib_finalize
  PUBLIC :: dbcsr_mm_csr_phaseout
  PUBLIC :: dbcsr_mm_csr_init, dbcsr_mm_csr_finalize
  PUBLIC :: dbcsr_mm_csr_multiply, dbcsr_mm_csr_purge_stacks

  CONTAINS

! *****************************************************************************
! \brief Initialize the library
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_csr_lib_init(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

     CALL dbcsr_mm_sched_lib_init(error)
  END SUBROUTINE


! *****************************************************************************
! \brief Finalize the library
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_csr_lib_finalize(group, output_unit, error)
    INTEGER, INTENT(IN)                      :: group, output_unit
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

     CALL dbcsr_mm_sched_lib_finalize(group, output_unit, error)
  END SUBROUTINE


! *****************************************************************************
! \brief A wrapper around dbcsr_mm_csr_multiply_low to avoid expensive dereferencings.
! \author Ole Schuett
! *****************************************************************************
 SUBROUTINE dbcsr_mm_csr_multiply(this, left, right, mi, mf, ni, nf, ki, kf,&
       ai, af,&
       bi, bf,&
       m_sizes, n_sizes, k_sizes,&
       c_local_rows, c_local_cols,&
       c_has_symmetry, keep_sparsity, use_eps,&
       row_max_epss,&
       flop,&
       a_index, b_index, a_norms, b_norms,&
       error)
    TYPE(dbcsr_mm_csr_type), INTENT(INOUT)   :: this
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
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

    nblks_new = this%product_wm%lastblk + max_new_nblks

    CALL ensure_array_size(this%product_wm%row_i, ub=nblks_new,&
         factor=default_resize_factor, error=error)
    CALL ensure_array_size(this%product_wm%col_i, ub=nblks_new,&
         factor=default_resize_factor, error=error)
    CALL ensure_array_size(this%product_wm%blk_p, ub=nblks_new,&
         factor=default_resize_factor, error=error)

    CALL dbcsr_mm_csr_multiply_low(this, left=left, right=right,&
            mi=mi, mf=mf,ni=ni, nf=nf, ki=ki, kf=kf,&
            ai=ai, af=af,&
            bi=bi, bf=bf,&
            c_row_i=this%product_wm%row_i,&
            c_col_i=this%product_wm%col_i,&
            c_blk_p=this%product_wm%blk_p,&
            lastblk =this%product_wm%lastblk, &
            datasize=this%product_wm%datasize,&
            m_sizes=m_sizes, n_sizes=n_sizes, k_sizes=k_sizes,&
            c_local_rows=c_local_rows, c_local_cols=c_local_cols,&
            c_has_symmetry=c_has_symmetry, keep_sparsity=keep_sparsity,&
            use_eps=use_eps,&
            row_max_epss=row_max_epss,&
            flop=flop,&
            row_size_maps=this%m_size_maps,& 
            col_size_maps=this%n_size_maps,& 
            k_size_maps=this%k_size_maps,&
            row_size_maps_size=this%m_size_maps_size,& 
            col_size_maps_size=this%n_size_maps_size,& 
            k_size_maps_size=this%k_size_maps_size,&
            nm_stacks=this%nm_stacks, nn_stacks=this%nn_stacks,&
            nk_stacks=this%nk_stacks, &
            stack_map=this%stack_map,&
            stacks_data=this%stacks_data,&
            stacks_fillcount=this%stacks_fillcount,&
            c_hashes=this%c_hashes,&
            a_index=a_index, b_index=b_index,&
            a_norms=a_norms, b_norms=b_norms,&
            error=error)

   END SUBROUTINE dbcsr_mm_csr_multiply


! *****************************************************************************
!> \brief Performs multiplication of smaller submatrices.
! *****************************************************************************
  SUBROUTINE dbcsr_mm_csr_multiply_low(this, left, right, mi, mf, ni, nf, ki, kf,&
       ai, af, bi, bf,&
       c_row_i, c_col_i, c_blk_p, lastblk, datasize,&
       m_sizes, n_sizes, k_sizes,&
       c_local_rows, c_local_cols,&
       c_has_symmetry, keep_sparsity, use_eps,&
       row_max_epss, flop,&
       row_size_maps, col_size_maps, k_size_maps,&
       row_size_maps_size, col_size_maps_size, k_size_maps_size,&
       nm_stacks, nn_stacks, nk_stacks, stack_map,&
       stacks_data, stacks_fillcount, c_hashes,&
       a_index, b_index,a_norms, b_norms,&
       error)
    TYPE(dbcsr_mm_csr_type), INTENT(INOUT)   :: this
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    INTEGER, INTENT(IN)                      :: mi, mf, ni, nf, ki, kf, ai, &
                                                af, bi, bf
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: c_row_i, c_col_i, c_blk_p
    INTEGER, INTENT(INOUT)                   :: lastblk, datasize
    INTEGER, DIMENSION(:), INTENT(IN)        :: m_sizes, n_sizes, k_sizes, &
                                                c_local_rows, c_local_cols
    LOGICAL, INTENT(IN)                      :: c_has_symmetry, &
                                                keep_sparsity, use_eps
    REAL(kind=sp), DIMENSION(:)              :: row_max_epss
    INTEGER(KIND=int_8), INTENT(INOUT)       :: flop
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
    INTEGER, DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: stacks_data
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: stacks_fillcount
    TYPE(hash_table_type), DIMENSION(:), &
      INTENT(INOUT)                          :: c_hashes
    INTEGER, DIMENSION(1:3, 1:af), &
      INTENT(IN)                             :: a_index
    INTEGER, DIMENSION(1:3, 1:bf), &
      INTENT(IN)                             :: b_index
    REAL(KIND=sp), DIMENSION(:), POINTER     :: a_norms, b_norms
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
    INTEGER, DIMENSION(2, bf-bi+1)           :: b_blk_info
    INTEGER, DIMENSION(2, af-ai+1)           :: a_blk_info
    INTEGER(KIND=int_4)                      :: offset
    LOGICAL                                  :: block_exists
    REAL(kind=sp)                            :: a_norm, a_row_eps, b_norm
    REAL(KIND=sp), DIMENSION(1:af-ai+1)      :: left_norms
    REAL(KIND=sp), DIMENSION(1:bf-bi+1)      :: right_norms

!   ---------------------------------------------------------------------------

    ithread = 0
    !$ ithread = omp_get_thread_num()

    nstacks = SIZE(this%stacks_data, 3)

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

             ! TODO: this is only called with careful_mod
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
             stacks_fillcount(ws) = stacks_fillcount(ws) + 1
             s_dp = stacks_fillcount(ws)

             stacks_data(p_m, s_dp, ws) = m_size
             stacks_data(p_n, s_dp, ws) = n_size
             stacks_data(p_k, s_dp, ws) = k_size
             stacks_data(p_a_first, s_dp, ws) = a_blk_info(2, a_blk)
             stacks_data(p_b_first, s_dp, ws) = b_blk_info(2, b_blk)
             stacks_data(p_c_first, s_dp, ws) = offset
             stacks_data(p_c_blk,   s_dp, ws) = c_blk_id

             flop = flop + INT(2*c_nze, int_8) * INT(k_size, int_8)

             IF(stacks_fillcount(ws) >= SIZE(stacks_data, 2))&
                CALL flush_stacks(this, left=left, right=right, error=error)

          ENDDO b_blk_cycle ! b
       ENDDO a_blk_cycle ! a_col
    ENDDO a_row_cycle ! a_row

  END SUBROUTINE dbcsr_mm_csr_multiply_low

! *****************************************************************************
!> \brief Initializes a multiplication cycle for new set of C-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_csr_init(this, left, right, product, &
    m_sizes, n_sizes, k_sizes, error)
    TYPE(dbcsr_mm_csr_type), INTENT(INOUT)   :: this
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    INTEGER, DIMENSION(:), POINTER           :: m_sizes, n_sizes, k_sizes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_csr_init', &
      routineP = moduleN//':'//routineN

    INTEGER :: block_estimate, default_stack, error_handler, istack, ithread, &
      k_map, k_size, m_map, m_size, n_map, n_size, nstacks, nthreads, ps_g
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: flop_index, flop_list, &
                                                most_common_k, most_common_m, &
                                                most_common_n
    INTEGER, DIMENSION(3)                    :: nxstacks
    LOGICAL                                  :: local_indexing
    TYPE(stack_descriptor_type), &
      ALLOCATABLE, DIMENSION(:)              :: tmp_descr

    CALL dbcsr_error_set(routineN, error_handler, error)

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
       ALLOCATE (this%c_hashes (product%nblkrows_local))
       CALL fill_hash_tables (this%c_hashes, product,block_estimate,&
            row_map=array_data(product%global_rows),&
            col_map=array_data(product%global_cols),&
            error=error)
    ELSE
       ALLOCATE (this%c_hashes (product%nblkrows_total))
       CALL fill_hash_tables (this%c_hashes, product,block_estimate,&
            error=error)
    ENDIF

    ! Setup the MM stack
    CALL dbcsr_get_conf_nstacks (nxstacks, error)
    this%nm_stacks = nxstacks(1)
    this%nn_stacks = nxstacks(2)
    this%nk_stacks = nxstacks(3)
    nstacks = nxstacks(1) * nxstacks(2) * nxstacks(3) + 1
    CALL dbcsr_assert (nstacks, "LE", INT (HUGE (this%stack_map)),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Too many stacks requested (global/dbcsr/n_size_*_stacks in input)",&
         __LINE__, error=error)


    ALLOCATE(this%stacks_descr(nstacks))
    ALLOCATE(this%stacks_data(dbcsr_ps_width, mm_stack_size, nstacks))
    ALLOCATE(this%stacks_fillcount(nstacks))
    this%stacks_fillcount(:) = 0

    ALLOCATE (most_common_m(nxstacks(1)))
    ALLOCATE (most_common_n(nxstacks(2)))
    ALLOCATE (most_common_k(nxstacks(3)))
    CALL map_most_common (m_sizes, this%m_size_maps, nxstacks(1),&
         most_common_m,&
         max_stack_block_size, this%max_m)
    this%m_size_maps_size = SIZE (this%m_size_maps)
    CALL map_most_common (n_sizes, this%n_size_maps, nxstacks(2),&
         most_common_n,&
         max_stack_block_size, this%max_n)
    this%n_size_maps_size = SIZE (this%n_size_maps)
    CALL map_most_common (array_data(right%row_blk_size),&
         this%k_size_maps, nxstacks(3), &
         most_common_k,&
         max_stack_block_size, this%max_k)
    this%k_size_maps_size = SIZE (this%k_size_maps)

    ! Creates the stack map--a mapping from (mapped) stack block sizes
    ! (carrier%*_sizes) to a stack number.  Triples with even one
    ! uncommon size will be mapped to a general, non-size-specific
    ! stack.
    ALLOCATE (this%stack_map(nxstacks(2)+1, nxstacks(3)+1, nxstacks(1)+1))
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
                this%stack_map(n_map, k_map, m_map) = ps_g
                ! Also take care of the stack m, n, k descriptors
                this%stacks_descr(ps_g)%m     = m_size
                this%stacks_descr(ps_g)%n     = n_size
                this%stacks_descr(ps_g)%k     = k_size
                this%stacks_descr(ps_g)%max_m = m_size
                this%stacks_descr(ps_g)%max_n = n_size
                this%stacks_descr(ps_g)%max_k = k_size
                this%stacks_descr(ps_g)%defined_mnk = .TRUE.
             ELSE
                ! This is the case when at least one of m, n, or k is
                ! undefined.
                ps_g = default_stack
                this%stack_map(n_map, k_map, m_map) = default_stack
                ! Also take care of the stack m, n, k descriptors
                this%stacks_descr(ps_g)%m     = 0
                this%stacks_descr(ps_g)%n     = 0
                this%stacks_descr(ps_g)%k     = 0
                this%stacks_descr(ps_g)%max_m = this%max_m
                this%stacks_descr(ps_g)%max_n = this%max_n
                this%stacks_descr(ps_g)%max_k = this%max_k
                this%stacks_descr(ps_g)%defined_mnk = .FALSE.
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
       flop_list(istack) = -2 * this%stacks_descr(istack)%m& 
                              * this%stacks_descr(istack)%n&
                              * this%stacks_descr(istack)%k
    ENDDO

    CALL sort(flop_list, nstacks-1, flop_index)
    tmp_descr = this%stacks_descr
    DO istack=1,nstacks-1
       this%stacks_descr(istack) = tmp_descr(flop_index(istack)) 
    ENDDO

    DO m_map = 1, SIZE(this%stack_map, 1)
      DO k_map = 1, SIZE(this%stack_map, 2)
        map_loop: DO n_map = 1, SIZE(this%stack_map, 1)
          DO istack=1,nstacks-1
            IF(this%stack_map(m_map, k_map, n_map) == flop_index(istack)) THEN
               this%stack_map(m_map, k_map, n_map) = istack
               CYCLE map_loop
            END IF
          ENDDO
        ENDDO map_loop
      ENDDO
    ENDDO
    DEALLOCATE(flop_list,flop_index,tmp_descr)

    this%product_wm => product%wms(ithread+1)
    CALL dbcsr_mm_sched_init(this%sched, left=left, right=right,&
                product_wm=this%product_wm,error=error)

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_mm_csr_init


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
! \brief Signal approaching end of multiplication
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_csr_phaseout(this, error)
    TYPE(dbcsr_mm_csr_type), INTENT(INOUT)   :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CALL dbcsr_mm_sched_phaseout(this%sched, error)
  END SUBROUTINE dbcsr_mm_csr_phaseout


! *****************************************************************************
!> \brief Finalizes a multiplication cycle for a set of C-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_csr_finalize(this, error)
    TYPE(dbcsr_mm_csr_type), INTENT(INOUT)   :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_csr_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    CALL dbcsr_mm_sched_finalize(this%sched, error=error)

    ! Clear hash tables
    DO i = 1, SIZE(this%c_hashes)
       CALL hash_table_release (this%c_hashes(i))
    ENDDO
    DEALLOCATE(this%c_hashes)
    DEALLOCATE(this%stacks_descr)
    DEALLOCATE(this%stack_map)
    DEALLOCATE(this%m_size_maps)
    DEALLOCATE(this%n_size_maps)
    DEALLOCATE(this%k_size_maps)
    DEALLOCATE(this%stacks_fillcount)
    DEALLOCATE(this%stacks_data)
  END SUBROUTINE dbcsr_mm_csr_finalize



! *****************************************************************************
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE dbcsr_mm_csr_purge_stacks(this, left, right, error)
    TYPE(dbcsr_mm_csr_type), INTENT(INOUT)   :: this
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CALL flush_stacks(this, left, right, purge=.TRUE., error=error)
    CALL dbcsr_mm_sched_barrier(this%sched, error)
END SUBROUTINE dbcsr_mm_csr_purge_stacks

! *****************************************************************************
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE flush_stacks(this, left, right, purge, error)
    TYPE(dbcsr_mm_csr_type), INTENT(INOUT)   :: this
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    LOGICAL, INTENT(IN), OPTIONAL            :: purge
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i, min_fill, n_stacks
    INTEGER, DIMENSION(:, :), POINTER        :: stack_data
    INTEGER, POINTER                         :: stack_fillcount
    TYPE(stack_descriptor_type)              :: stack_descr

    n_stacks = SIZE(this%stacks_data, 3)
    min_fill = SIZE(this%stacks_data, 2) * 3 / 4 !TODO: play with this
    IF(PRESENT(purge)) THEN
      IF(purge) min_fill = 0
    ENDIF

    CALL dbcsr_mm_sched_begin_burst(this%sched, error)

    DO i=1, n_stacks
       IF(this%stacks_fillcount(i) > min_fill) THEN
          stack_data=>this%stacks_data(:,:,i)
          stack_fillcount=>this%stacks_fillcount(i)
          stack_descr = this%stacks_descr(i)

          CALL dbcsr_mm_sched_process(this%sched,&
                        left,right,&
                        stack_data=stack_data,&
                        stack_fillcount=stack_fillcount,&
                        stack_descr=stack_descr,&
                        error=error)

          stack_fillcount = 0
       ENDIF
    END DO

    CALL dbcsr_mm_sched_end_burst(this%sched, error)
END SUBROUTINE flush_stacks


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
