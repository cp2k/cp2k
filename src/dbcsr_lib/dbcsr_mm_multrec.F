!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Second layer of the dbcsr matrix-matrix multiplication.
!> 	   It divides the multiplication in a cache-oblivious manner.
!>
!> \author  Urban Borstnik
!>
!> <b>Modification history:</b>
!>  - 2010-02-23 Moved from dbcsr_operations
!>  - 2011-11    Moved parameter-stack processing routines to
!>               dbcsr_mm_methods.
!>  - 2013-01    reorganized code (Ole Schuett)
! *****************************************************************************
	
MODULE dbcsr_mm_multrec

  USE array_types,                     ONLY: array_data,&
                                             array_equality
  USE dbcsr_config,                    ONLY: mm_thread_workshare
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: int_8,&
                                             real_8,&
                                             sp
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_col_dist,&
                                             dbcsr_distribution_has_threads,&
                                             dbcsr_distribution_local_cols,&
                                             dbcsr_distribution_local_rows,&
                                             dbcsr_distribution_row_dist,&
                                             dbcsr_distribution_thread_dist
  USE dbcsr_mm_csr,                    ONLY: csr_multiply_inner_finalize,&
                                             csr_multiply_inner_init,&
                                             csr_multiply_outer_finalize,&
                                             csr_multiply_outer_init,&
                                             csr_multiply_unwrap
  USE dbcsr_mm_types,                  ONLY: carrier_type
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
  USE dbcsr_types,                     ONLY: dbcsr_type

  !$ USE OMP_LIB
        
  IMPLICIT NONE
 
  PRIVATE
 
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_multrec'
  LOGICAL, PARAMETER :: careful_mod = .FALSE.

  PUBLIC :: dbcsr_multrec_init
  PUBLIC :: dbcsr_multrec_sparse                                                                   
  PUBLIC :: dbcsr_multrec_finalize
  
CONTAINS

! *****************************************************************************
!> \brief Sets up recursive multiplication
!>
!>
!> \param[in] left, right     left and right DBCSR matrices
!> \param[in,out] product     resulting DBCSR product matrix
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix, default is no
!> \param[in] eps             (optional) on-the-fly filtering epsilon
! *****************************************************************************
  SUBROUTINE dbcsr_multrec_init(left, right, product, carrier,&
       right_data_sr, right_data_rr, left_data_sr, left_data_rr, &
       right_index_sr, right_index_rr, left_index_sr, left_index_rr, &
       keep_sparsity, eps, row_max_epss, fill_guess, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(carrier_type), INTENT(out)          :: carrier
    INTEGER, DIMENSION(:, :), POINTER :: right_data_sr, right_data_rr, &
      left_data_sr, left_data_rr, right_index_sr, right_index_rr, &
      left_index_sr, left_index_rr
    LOGICAL, INTENT(IN)                      :: keep_sparsity
    REAL(kind=real_8), INTENT(in), OPTIONAL  :: eps
    REAL(kind=sp), DIMENSION(:), INTENT(IN)  :: row_max_epss
    REAL(KIND=real_8)                        :: fill_guess
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multrec_init', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: c_nblkcols_local, &
                                                c_nblkrows_local, &
                                                error_handler, i, ithread
    INTEGER, DIMENSION(:), POINTER           :: c_local_cols, c_local_rows, &
                                                product_thread_dist
    LOGICAL                                  :: local_cols, local_indexing, &
                                                local_rows

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    
    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM ()
    carrier%id = ithread
    !
    ! Ensures that the index is correctly defined.
    CALL dbcsr_assert (left%list_indexing,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Must use list indexing for this routine.", __LINE__, error=error)
    !
    CALL dbcsr_assert ("NOT", left%bcsc,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Wrong routine for BCSC matrices.", __LINE__, error=error)
    CALL dbcsr_assert ("NOT", right%bcsc,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Wrong routine for BCSC matrices.", __LINE__, error=error)
    local_cols = right%local_indexing
    local_rows = left%local_indexing
    CALL dbcsr_assert (local_cols, "EQV", local_rows,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Local index useage must be consistent.", __LINE__, error=error)
    local_indexing = local_rows
    IF (local_cols) THEN
       CALL dbcsr_assert (left%local_indexing,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Wrong left format for local_cols.", __LINE__, error=error)
       CALL dbcsr_assert (right%local_indexing,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Wrong right format for local_cols.", __LINE__, error=error)
    ELSE
       CALL dbcsr_assert ("NOT",left%local_indexing,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Wrong left format for not local_cols.", __LINE__, error=error)
       CALL dbcsr_assert ("NOT",right%local_indexing,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Wrong right format for not local_cols.", __LINE__, error=error)
    ENDIF
    !
    ! Fill carrier data structure.
    carrier%local_indexing = local_indexing
    carrier%keep_sparsity = keep_sparsity
    carrier%c_has_symmetry = product%symmetry
    carrier%use_eps = PRESENT (eps)
    carrier%my_wm = product%wms(ithread+1)
    carrier%lastblk = product%wms(ithread+1)%lastblk
    carrier%original_lastblk = carrier%lastblk
    carrier%datasize = product%wms(ithread+1)%datasize
    carrier%flop = INT(0, int_8)
    carrier%right_data_sr => right_data_sr
    carrier%right_data_rr => right_data_rr
    carrier%right_index_sr => right_index_sr
    carrier%right_index_rr => right_index_rr
    carrier%left_data_sr => left_data_sr
    carrier%left_data_rr => left_data_rr
    carrier%left_index_sr => left_index_sr
    carrier%left_index_rr => left_index_rr
    IF (PRESENT (eps)) THEN
       carrier%eps = eps
    ELSE
       carrier%eps = 0.0_real_8
    ENDIF
    !
    !
    !$ NULLIFY (product_thread_dist)
    !$ CALL dbcsr_assert (dbcsr_distribution_has_threads (product%dist),&
    !$      dbcsr_fatal_level, dbcsr_internal_error, routineN,&
    !$      "Missing thread distribution.", __LINE__, error=error)
    !$ product_thread_dist => array_data (&
    !$      dbcsr_distribution_thread_dist (product%dist))
    !
    ! Find out the C/A rows and C/B columns and sizes.
    c_nblkrows_local = product%nblkrows_local
    c_local_rows => array_data (product%local_rows)
    c_nblkcols_local = product%nblkcols_local
    c_local_cols => array_data (product%local_cols)
    IF (local_indexing) THEN
       carrier%c_local_rows => c_local_rows
       carrier%c_local_cols => c_local_cols
    ELSE
       ALLOCATE (carrier%c_local_rows (product%nblkrows_total))
       ALLOCATE (carrier%c_local_cols (product%nblkcols_total))
       FORALL (i = 1 : product%nblkrows_total)
          carrier%c_local_rows(i) = i
       END FORALL
       FORALL (i = 1 : product%nblkcols_total)
          carrier%c_local_cols(i) = i
       END FORALL
    ENDIF
    IF (dbg) WRITE(*,*)"setting up for product", product%name
    IF (careful_mod) THEN
       IF (.NOT. array_equality (dbcsr_distribution_local_rows (product%dist),&
                                 product%local_rows)) THEN
          WRITE(*,*)"row dist", array_data(dbcsr_distribution_row_dist(product%dist))
          WRITE(*,*)"dist local rows", array_data(dbcsr_distribution_local_rows (product%dist))
          WRITE(*,*)" mat local rows", array_data(product%local_rows)
          CALL dbcsr_assert (.FALSE., &
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Array mismatch.", __LINE__, error=error)
       ENDIF
       IF (.NOT. array_equality (dbcsr_distribution_local_cols (product%dist),&
                                 product%local_cols)) THEN
          WRITE(*,*)"col dist", array_data(dbcsr_distribution_col_dist(product%dist))
          WRITE(*,*)"dist local cols", array_data(dbcsr_distribution_local_cols (product%dist))
          WRITE(*,*)" mat local cols", array_data(product%local_cols)
          CALL dbcsr_assert (.FALSE., &
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Array mismatch.", __LINE__, error=error)
       ENDIF
       CALL dbcsr_assert (SIZE(c_local_rows), "EQ", c_nblkrows_local,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Row count mismatch.", __LINE__, error=error)
       CALL dbcsr_assert (SIZE(c_local_cols), "EQ", c_nblkcols_local,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Column count mismatch.", __LINE__, error=error)
    ENDIF
    !
    ! And the k epsilons
    IF (local_indexing) THEN
       ALLOCATE (carrier%row_max_epss(c_nblkrows_local))
    ELSE
       ALLOCATE (carrier%row_max_epss(product%nblkrows_total))
    ENDIF
    IF (carrier%use_eps) THEN
       IF (local_indexing) THEN
          CALL local_filter_sp(row_max_epss, c_nblkrows_local, c_local_rows,&
               carrier%row_max_epss)
       ELSE
          IF (careful_mod) &
               CALL dbcsr_assert (SIZE(row_max_epss) .EQ. SIZE(carrier%row_max_epss),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "max epss local/global mismatch.", __LINE__, error=error)
          carrier%row_max_epss(:) = row_max_epss(:)
       ENDIF
    ELSE
       carrier%row_max_epss(:) = -HUGE(0.0_sp)
    ENDIF
    !
    IF (local_indexing) THEN
       ALLOCATE(carrier%m_sizes(c_nblkrows_local))
       CALL local_filter(array_data (product%row_blk_size), SIZE(c_local_rows),&
            c_local_rows, carrier%m_sizes)
       ALLOCATE(carrier%n_sizes(c_nblkcols_local))
       CALL local_filter(array_data (product%col_blk_size), SIZE(c_local_cols),&
            c_local_cols, carrier%n_sizes)
    ELSE
       ALLOCATE(carrier%m_sizes(product%nblkrows_total))
       carrier%m_sizes(:) = array_data (product%row_blk_size)
       ALLOCATE(carrier%n_sizes(product%nblkcols_total))
       carrier%n_sizes(:) = array_data (product%col_blk_size)
    ENDIF
    !
    NULLIFY (carrier%k_locals)
    NULLIFY (carrier%k_sizes)
    IF (.NOT. local_indexing) THEN
       ALLOCATE (carrier%k_locals(right%nblkrows_total))
       FORALL (i = 1:right%nblkrows_total)
          carrier%k_locals(i) = i
       END FORALL
       carrier%k_sizes => array_data (right%row_blk_size)
    ENDIF
    
    CALL csr_multiply_outer_init(left, right, product, carrier, fill_guess, error=error)
    
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_multrec_init


! *****************************************************************************
!> \brief Multiplies two DBCSR matrices using recursive algorithm
!>
!> This routine sets up the multiplication.  Specifically, it <ul>
!> <li> verifies input sanity
!> <li> converts everything into "local indexing"
!> </ul>
!>
!> \param[in] left, right     left and right DBCSR matrices
!> \param[in,out] product     resulting DBCSR product matrix
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix, default is no
!> \param[in] a_norms         (optional) norms of left-matrix blocks
!> \param[in] b_norms         (optional) norms of right-matrix blocks
!> \param[in] eps             (optional) on-the-fly filtering epsilon
!> \param[in] row_max_epss  (optional) on-the-fly filtering epsilon per block
!> \param[out] flop           (optional) number of effective double-precision
!>                            floating point operations performed
! *****************************************************************************
  SUBROUTINE dbcsr_multrec_sparse(left, right, product, carrier, flop,&
       keep_sparsity, a_norms, b_norms, eps, row_max_epss, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(carrier_type), INTENT(inout)        :: carrier
    INTEGER(KIND=int_8), INTENT(OUT)         :: flop
    LOGICAL, INTENT(IN)                      :: keep_sparsity
    REAL(kind=sp), DIMENSION(:), &
      INTENT(in), TARGET                     :: a_norms, b_norms
    REAL(kind=real_8), INTENT(in), OPTIONAL  :: eps
    REAL(kind=sp), DIMENSION(:), &
      INTENT(in), TARGET                     :: row_max_epss
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multrec_sparse', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: ithread, nthreads, t_a_f, &
                                                t_a_l, t_b_f, t_b_l
    INTEGER, DIMENSION(:), POINTER           :: k_locals
    LOGICAL                                  :: local_cols, local_indexing, &
                                                local_rows

!   ---------------------------------------------------------------------------
!

    ithread = 0 ; nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM () ; nthreads = OMP_GET_NUM_THREADS ()
    carrier%flop = 0
   
    CALL csr_multiply_inner_init(left, right, carrier, error=error)
    
    local_cols = right%local_indexing
    local_rows = left%local_indexing
    local_indexing = local_rows
    
    ! Find out the local A columns / B rows and sizes
    ! The right%local_rows is setup by the communication engine.
    IF (local_indexing) THEN
       k_locals => array_data (right%local_rows)
       carrier%k_locals => k_locals
       CALL ensure_array_size (carrier%k_sizes, ub=SIZE(k_locals), error=error)
       CALL local_filter(array_data(right%row_blk_size), SIZE(k_locals),&
            k_locals, carrier%k_sizes)
    ELSE
       k_locals => carrier%k_locals
    ENDIF
    ! Setup the block norms
    carrier%a_norms => a_norms
    carrier%b_norms => b_norms
    !
    IF (mm_thread_workshare) THEN
       !$OMP BARRIER
    ENDIF
    !
    ! Start local multiplication
    IF (.TRUE.) THEN
       t_a_f = 1
       t_a_l = left%nblks
       t_b_f = 1
       t_b_l = right%nblks
       !$ ithread = OMP_GET_THREAD_NUM()
       !$ t_a_f = left%thr_c(ithread+1)+1
       !$ t_a_l = left%thr_c(ithread+2)
       IF (left%local_indexing) THEN
          CALL sparse_multrec(&
               1, left%nblkrows_local,&
               1, right%nblkcols_local,&
               1, SIZE(k_locals),&
               t_a_f, t_a_l, left%coo_l,&
               t_b_f, t_b_l, right%coo_l,&
               carrier, error, 0)
       ELSE
          CALL sparse_multrec(&
               1, left%nblkrows_total,&
               1, right%nblkcols_total,&
               1, SIZE(k_locals),&
               t_a_f, t_a_l, left%coo_l,&
               t_b_f, t_b_l, right%coo_l,&
               carrier, error, 0)
       ENDIF
    ELSE
       ! For debugging.  This avoids the recursion but it might not
       ! work anymore.
       IF (left%local_indexing) THEN
          CALL csr_multiply_unwrap(&
               1, left%nblkrows_local,&
               1, right%nblkrows_local,&
               1, SIZE(k_locals),&
               1, left%nblks, left%row_p,&
               1, right%nblks, right%row_p,&
               carrier, error)
       ELSE
          CALL csr_multiply_unwrap(&
               1, left%nblkrows_total,&
               1, right%nblkrows_total,&
               1, SIZE(k_locals),&
               1, left%nblks, left%row_p,&
               1, right%nblks, right%row_p,&
               carrier, error)
       ENDIF
    ENDIF
    !
    CALL csr_multiply_inner_finalize(carrier, error=error)
    flop = carrier%flop
    !
    IF (ASSOCIATED (carrier%k_sizes) .AND. carrier%local_indexing) &
         DEALLOCATE (carrier%k_sizes)
    !
    product%wms(ithread+1)%lastblk = carrier%lastblk
    product%wms(ithread+1)%datasize = carrier%datasize
    !
  END SUBROUTINE dbcsr_multrec_sparse            
  

! *****************************************************************************
!> \brief Sets up recursive multiplication
!>
!>
!> \param[in] left, right     left and right DBCSR matrices
!> \param[in,out] product     resulting DBCSR product matrix
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix, default is no
!> \param[in] eps             (optional) on-the-fly filtering epsilon
! *****************************************************************************
  SUBROUTINE dbcsr_multrec_finalize(carrier, product, error)
    TYPE(carrier_type), INTENT(inout)        :: carrier
    TYPE(dbcsr_type), INTENT(inout)          :: product
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multrec_finalize', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler, ithread

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL csr_multiply_outer_finalize(carrier, product, error=error)
    
    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM()
    
    ! Release the carrier
    DEALLOCATE (carrier%m_sizes)
    DEALLOCATE (carrier%n_sizes)
    DEALLOCATE (carrier%row_max_epss)
    
    IF (carrier%local_indexing) THEN
       CALL remap_local2global(carrier%my_wm%row_i, carrier%my_wm%col_i,&
            carrier%c_local_rows, carrier%c_local_cols,&
            carrier%original_lastblk+1, carrier%lastblk)
    ENDIF
    !
    DEALLOCATE (carrier%m_size_maps)
    DEALLOCATE (carrier%n_size_maps)
    DEALLOCATE (carrier%k_size_maps)
    DEALLOCATE (carrier%stack_map)
    IF (.NOT. carrier%local_indexing) THEN
       DEALLOCATE (carrier%c_local_rows)
       DEALLOCATE (carrier%c_local_cols)
       DEALLOCATE (carrier%k_locals)
       NULLIFY (carrier%k_sizes)
    ENDIF
    !
    ! Reinstate WM
    carrier%my_wm%lastblk = carrier%lastblk
    carrier%my_wm%datasize = carrier%datasize
    product%wms(ithread+1) = carrier%my_wm

    NULLIFY(carrier%right_data_sr)
    NULLIFY(carrier%right_data_rr)
    NULLIFY(carrier%right_index_sr)
    NULLIFY(carrier%right_index_rr)
    NULLIFY(carrier%left_data_sr)
    NULLIFY(carrier%left_data_rr)
    NULLIFY(carrier%left_index_sr)
    NULLIFY(carrier%left_index_rr)
    !
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_multrec_finalize
  
  
  
! *****************************************************************************
!> \brief Performs recursive multiplication
!> \author Joost VandeVondele
! *****************************************************************************
  RECURSIVE SUBROUTINE sparse_multrec(mi, mf, ni, nf, ki, kf,&
       ai, af, a_index, bi, bf, b_index, &
       carrier, error, d)
    INTEGER, INTENT(IN)                      :: mi, mf, ni, nf, ki, kf, ai, af
    INTEGER, DIMENSION(3, 1:af), INTENT(IN)  :: a_index
    INTEGER, INTENT(IN)                      :: bi, bf
    INTEGER, DIMENSION(3, 1:bf), INTENT(IN)  :: b_index
    TYPE(carrier_type), INTENT(INOUT)        :: carrier
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER, INTENT(IN)                      :: d

    CHARACTER(len=*), PARAMETER :: routineN = 'sparse_multrec', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: norec = 512
    LOGICAL, PARAMETER                       :: careful = careful_mod, &
                                                dbg = .FALSE.

    INTEGER                                  :: acut, bcut, cut, K, M, N, s1

!   ---------------------------------------------------------------------------

    IF (dbg) THEN
       WRITE(*,'(I7,1X,5(A,2(1X,I7)))')d," rm", mi, mf,",",ni,nf,",",ki,kf,"/",ai,af,",",bi,bf
    ENDIF
    IF (.TRUE.) THEN
       IF (af .LT. ai .OR. bf .LT. bi .OR. mf .LT. mi .OR. nf .LT. ni .OR. kf .LT. ki) THEN
          IF (dbg) WRITE(*,*)"Empty"
          RETURN
       ENDIF
    ENDIF
    IF (af-ai+1 <= norec .AND. bf-bi+1 <= norec) THEN
       CALL csr_multiply_unwrap(&
            mi, mf, ni, nf, ki, kf,&
            ai, af, a_index,&
            bi, bf, b_index,&
            carrier, error)
       RETURN
    ENDIF

    M = mf-mi + 1
    N = nf-ni + 1
    K = kf-ki + 1
    IF (dbg) THEN
       WRITE(*,*)'m,k,n',M,K,N
    ENDIF
    IF (M >= MAX(N, K)) cut = 1
    IF (K >= MAX(N, M)) cut = 2
    IF (N >= MAX(M, K)) cut = 3
    SELECT CASE(cut)
    CASE(1)
       s1=M/2
       acut = find_cut_row(ai,af,a_index,mi+s1-1)
       CALL sparse_multrec(mi,mi+s1-1, ni,nf, ki,kf,&
            ai,acut-1,a_index, bi,bf,b_index, carrier, error,d+1)
       CALL sparse_multrec(mi+s1,mf, ni,nf, ki,kf,&
            acut,af,a_index, bi,bf,b_index, carrier, error,d+1)
    CASE(2)
       s1=K/2
       acut = find_cut_col(ai,af,a_index,ki+s1-1)
       IF (dbg) THEN
          WRITE(*,*)N,s1,ni+s1-1,"/",ai,af,acut
          WRITE(*,'(3(I7))')a_index
       ENDIF
       bcut = find_cut_row(bi,bf,b_index,ki+s1-1)
       IF (dbg) THEN
          WRITE(*,*)N,s1,ni+s1-1,"/",bi,bf,bcut
          WRITE(*,'(3(I7))')b_index
       ENDIF
       CALL sparse_multrec(mi,mf, ni,nf, ki,ki+s1-1,&
            ai,acut-1,a_index, bi,bcut-1,b_index, carrier, error,d+1)
       CALL sparse_multrec(mi,mf, ni,nf, ki+s1,kf,&
            acut,af,a_index, bcut,bf,b_index, carrier, error,d+1)
    CASE(3)
       s1=N/2
       bcut = find_cut_col(bi,bf,b_index,ni+s1-1)
       IF (dbg) THEN
          WRITE(*,*)N,s1,ni+s1-1,"/",bi,bf,bcut
          WRITE(*,'(3(I7))')b_index
       ENDIF
       CALL sparse_multrec(mi,mf, ni,ni+s1-1, ki,kf,&
            ai,af,a_index, bi,bcut-1,b_index, carrier, error,d+1)
       CALL sparse_multrec(mi,mf, ni+s1,nf, ki,kf,&
            ai,af,a_index, bcut,bf,b_index, carrier, error,d+1)
    END SELECT
  END SUBROUTINE sparse_multrec

              
! ******************************************************************************
!> \brief 
!> \author JV
! ******************************************************************************
  PURE FUNCTION find_cut_row(ai,af,a,val) RESULT(res)
    INTEGER, INTENT(IN)                      :: ai, af
    INTEGER, DIMENSION(3, 1:af), INTENT(IN)  :: a
    INTEGER, INTENT(IN)                      :: val
    INTEGER                                  :: res

    INTEGER                                  :: i, ihigh, ilow

! do a log(N) search along the ordered index

    ilow = ai
    IF (a(1,ilow) > val)  THEN
       res = ilow
       RETURN
    ENDIF

    ihigh = af
    IF (a(1,ihigh) <= val)  THEN
       res = ihigh+1
       RETURN
    ENDIF

    DO
       IF (ihigh-ilow == 1) EXIT
       i = (ilow + ihigh)/2
       IF (a(1,i)>val) THEN
          ihigh=i
       ELSE
          ilow=i
       ENDIF
    ENDDO
    res=ihigh

    ! the linear search version
    ! DO i=ai,af
    !    IF (a(i)%r>val) EXIT
    !ENDDO
    !res=i
  END FUNCTION find_cut_row                   

  
! ******************************************************************************
!> \brief 
!> \author JV
! ****************************************************************************** 
  PURE FUNCTION find_cut_col(ai,af,a,val) RESULT(res)
    INTEGER, INTENT(IN)                      :: ai, af
    INTEGER, DIMENSION(3, 1:af), INTENT(IN)  :: a
    INTEGER, INTENT(IN)                      :: val
    INTEGER                                  :: res

    INTEGER                                  :: i, ihigh, ilow

! do a log(N) search along the ordered index

    ilow = ai
    IF (a(2,ilow) > val)  THEN
       res = ilow
       RETURN
    ENDIF

    ihigh = af
    IF (a(2,ihigh) <= val)  THEN
       res = ihigh+1
       RETURN
    ENDIF

    DO
       IF (ihigh-ilow == 1) EXIT
       i = (ilow + ihigh)/2
       IF (a(2,i) > val) THEN
          ihigh = i
       ELSE
          ilow = i
       ENDIF
    ENDDO
    res = ihigh

    ! the linear search version
    ! DO i=ai,af
    !    IF (a(i)%c>val) EXIT
    !ENDDO
    !res=i
  END FUNCTION find_cut_col

  
  
! ******************************************************************************
!> \brief Packs a globally-indexed array into a locally-indexed array.
! ****************************************************************************** 
  PURE SUBROUTINE remap_local2global(row_i, col_i, local_rows, local_cols,&
       first, last)
    INTEGER, INTENT(in)                      :: last, first
    INTEGER, DIMENSION(:), INTENT(in)        :: local_cols, local_rows
    INTEGER, DIMENSION(1:last), &
      INTENT(inout)                          :: col_i, row_i

    INTEGER                                  :: i

    FORALL (i = first : last)
       row_i(i) = local_rows(row_i(i))
       col_i(i) = local_cols(col_i(i))
    END FORALL
  END SUBROUTINE remap_local2global

  !> \brief Maps between locally-indexed arrays with different local indexing.
  PURE SUBROUTINE remap_local2local(n_local_src, local_remap, local_src, &
       n_globals, global_dst)
    INTEGER, INTENT(IN)                      :: n_local_src
    INTEGER, DIMENSION(1:n_local_src), &
      INTENT(OUT)                            :: local_remap
    INTEGER, DIMENSION(1:n_local_src), &
      INTENT(IN)                             :: local_src
    INTEGER, INTENT(IN)                      :: n_globals
    INTEGER, DIMENSION(1:n_globals), &
      INTENT(IN)                             :: global_dst

    INTEGER                                  :: i

    FORALL (i = 1 : n_local_src)
       local_remap(i) = global_dst (local_src (i))
    END FORALL
  END SUBROUTINE remap_local2local


! *****************************************************************************
!> \brief Gathers the local elements from all data (full_data)
!>
!> \param[in] full_data       All elements
!> \param[in] nle             Number of local elements
!> \param[in] local_elements  List of local elements
!> \param[out] local_data     Local elements obtained from all elements
! *****************************************************************************
  PURE SUBROUTINE local_filter (full_data, nle, local_elements, local_data)
    INTEGER, DIMENSION(:), INTENT(IN)        :: full_data
    INTEGER, INTENT(IN)                      :: nle
    INTEGER, DIMENSION(1:nle), INTENT(IN)    :: local_elements
    INTEGER, DIMENSION(1:nle), INTENT(OUT)   :: local_data

    INTEGER                                  :: l

    FORALL (l = 1 : nle)
       local_data(l) = full_data(local_elements(l))
    END FORALL
  END SUBROUTINE local_filter

! *****************************************************************************
!> \brief Gathers the local elements from all data (full_data) for
!>        single precision elements.
!> \see local_filter
! *****************************************************************************
  PURE SUBROUTINE local_filter_sp (full_data, nle, local_elements, local_data)
    REAL(KIND=sp), DIMENSION(:), INTENT(IN)  :: full_data
    INTEGER, INTENT(IN)                      :: nle
    INTEGER, DIMENSION(1:nle), INTENT(IN)    :: local_elements
    REAL(KIND=sp), DIMENSION(1:nle), &
      INTENT(OUT)                            :: local_data

    INTEGER                                  :: l

    FORALL (l = 1 : SIZE(local_data))
       local_data(l) = full_data(local_elements(l))
    END FORALL
  END SUBROUTINE local_filter_sp
  


END MODULE dbcsr_mm_multrec
