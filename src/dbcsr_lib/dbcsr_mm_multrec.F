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
!>  - 2013-01    extensive refactoring (Ole Schuett)
! *****************************************************************************
	
MODULE dbcsr_mm_multrec

  USE array_types,                     ONLY: array_data,&
                                             array_equality
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
  USE dbcsr_mm_csr,                    ONLY: dbcsr_mm_csr_inner_finalize,&
                                             dbcsr_mm_csr_inner_init,&
                                             dbcsr_mm_csr_lib_finalize,&
                                             dbcsr_mm_csr_lib_init,&
                                             dbcsr_mm_csr_multiply,&
                                             dbcsr_mm_csr_outer_finalize,&
                                             dbcsr_mm_csr_outer_init
  USE dbcsr_mm_types,                  ONLY: multrec_carrier_type
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
  USE dbcsr_types,                     ONLY: dbcsr_type

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_multrec'
  LOGICAL, PARAMETER :: careful_mod = .FALSE.

  PUBLIC :: dbcsr_mm_multrec_lib_init, dbcsr_mm_multrec_lib_finalize
  PUBLIC :: dbcsr_mm_multrec_init, dbcsr_mm_multrec_finalize
  PUBLIC :: dbcsr_mm_multrec_multiply

  CONTAINS

! *****************************************************************************
! \brief Initialize the library
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_multrec_lib_init(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

     CALL dbcsr_mm_csr_lib_init(error)
  END SUBROUTINE


! *****************************************************************************
! \brief Finalize the library
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_multrec_lib_finalize(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

     CALL dbcsr_mm_csr_lib_finalize(error)
  END SUBROUTINE


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
  SUBROUTINE dbcsr_mm_multrec_init(left, right, product, multrec_carrier,&
       keep_sparsity, eps, row_max_epss, fill_guess, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(multrec_carrier_type), INTENT(out)  :: multrec_carrier
    LOGICAL, INTENT(IN)                      :: keep_sparsity
    REAL(kind=real_8), INTENT(in), OPTIONAL  :: eps
    REAL(kind=sp), DIMENSION(:), INTENT(IN)  :: row_max_epss
    REAL(KIND=real_8)                        :: fill_guess
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_multrec_init', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: c_nblkcols_local, &
                                                c_nblkrows_local, &
                                                error_handler, i, ithread
    INTEGER, DIMENSION(:), POINTER           :: c_local_cols, c_local_rows, &
                                                product_thread_dist

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM ()
    !
    ! Ensures that the index is correctly defined.
    CALL dbcsr_assert (left%list_indexing,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Must use list indexing for this routine.", __LINE__, error=error)

    CALL dbcsr_assert ("NOT", left%bcsc,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Wrong routine for BCSC matrices.", __LINE__, error=error)
    CALL dbcsr_assert ("NOT", right%bcsc,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Wrong routine for BCSC matrices.", __LINE__, error=error)
    CALL dbcsr_assert (right%local_indexing, "EQV", left%local_indexing,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Local index useage must be consistent.", __LINE__, error=error)

    IF (right%local_indexing) THEN
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
    ! Fill multrec_carrier data structure.
    multrec_carrier%local_indexing = left%local_indexing 
    multrec_carrier%keep_sparsity = keep_sparsity
    multrec_carrier%c_has_symmetry = product%symmetry
    multrec_carrier%use_eps = PRESENT (eps)
    multrec_carrier%original_lastblk = product%wms(ithread+1)%lastblk
    multrec_carrier%flop = INT(0, int_8)

    IF (PRESENT (eps)) THEN
       multrec_carrier%eps = eps
    ELSE
       multrec_carrier%eps = 0.0_real_8
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
    IF (multrec_carrier%local_indexing) THEN
       multrec_carrier%c_local_rows => c_local_rows
       multrec_carrier%c_local_cols => c_local_cols
    ELSE
       ALLOCATE (multrec_carrier%c_local_rows (product%nblkrows_total))
       ALLOCATE (multrec_carrier%c_local_cols (product%nblkcols_total))
       FORALL (i = 1 : product%nblkrows_total)
          multrec_carrier%c_local_rows(i) = i
       END FORALL
       FORALL (i = 1 : product%nblkcols_total)
          multrec_carrier%c_local_cols(i) = i
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
    IF (multrec_carrier%local_indexing) THEN
       ALLOCATE (multrec_carrier%row_max_epss(c_nblkrows_local))
    ELSE
       ALLOCATE (multrec_carrier%row_max_epss(product%nblkrows_total))
    ENDIF
    IF (multrec_carrier%use_eps) THEN
       IF (multrec_carrier%local_indexing) THEN
          CALL local_filter_sp(row_max_epss, c_nblkrows_local, c_local_rows,&
               multrec_carrier%row_max_epss)
       ELSE
          IF (careful_mod) &
               CALL dbcsr_assert (SIZE(row_max_epss) .EQ. SIZE(multrec_carrier%row_max_epss),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "max epss local/global mismatch.", __LINE__, error=error)
          multrec_carrier%row_max_epss(:) = row_max_epss(:)
       ENDIF
    ELSE
       multrec_carrier%row_max_epss(:) = -HUGE(0.0_sp)
    ENDIF
    !
    IF (multrec_carrier%local_indexing) THEN
       ALLOCATE(multrec_carrier%m_sizes(c_nblkrows_local))
       CALL local_filter(array_data (product%row_blk_size), SIZE(c_local_rows),&
            c_local_rows, multrec_carrier%m_sizes)
       ALLOCATE(multrec_carrier%n_sizes(c_nblkcols_local))
       CALL local_filter(array_data (product%col_blk_size), SIZE(c_local_cols),&
            c_local_cols, multrec_carrier%n_sizes)
    ELSE
       ALLOCATE(multrec_carrier%m_sizes(product%nblkrows_total))
       multrec_carrier%m_sizes(:) = array_data (product%row_blk_size)
       ALLOCATE(multrec_carrier%n_sizes(product%nblkcols_total))
       multrec_carrier%n_sizes(:) = array_data (product%col_blk_size)
    ENDIF
    !
    NULLIFY (multrec_carrier%k_locals)
    NULLIFY (multrec_carrier%k_sizes)
    IF (.NOT. multrec_carrier%local_indexing) THEN
       ALLOCATE (multrec_carrier%k_locals(right%nblkrows_total))
       FORALL (i = 1:right%nblkrows_total)
          multrec_carrier%k_locals(i) = i
       END FORALL
       multrec_carrier%k_sizes => array_data (right%row_blk_size)
    ENDIF

    CALL dbcsr_mm_csr_outer_init(csr_carrier=multrec_carrier%csr_carrier, &
       left=left, right=right, product=product,fill_guess=fill_guess,&
       m_sizes=multrec_carrier%m_sizes, n_sizes=multrec_carrier%n_sizes,& 
       k_sizes=multrec_carrier%k_sizes, error=error)

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_mm_multrec_init


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
  SUBROUTINE dbcsr_mm_multrec_multiply(left, right, product, multrec_carrier, flop,&
       keep_sparsity, a_norms, b_norms, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(multrec_carrier_type), &
      INTENT(inout)                          :: multrec_carrier
    INTEGER(KIND=int_8), INTENT(OUT)         :: flop
    LOGICAL, INTENT(IN)                      :: keep_sparsity
    REAL(kind=sp), DIMENSION(:), &
      INTENT(in), TARGET                     :: a_norms, b_norms
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_multrec_multiply', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: ithread, nthreads, t_a_f, &
                                                t_a_l, t_b_f, t_b_l
    INTEGER, DIMENSION(:), POINTER           :: k_locals

!   ---------------------------------------------------------------------------

    ithread = 0 ; nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM () ; nthreads = OMP_GET_NUM_THREADS ()
    multrec_carrier%flop = 0

    CALL dbcsr_mm_csr_inner_init(left=left, right=right, &
              csr_carrier=multrec_carrier%csr_carrier, error=error)

    ! Find out the local A columns / B rows and sizes
    ! The right%local_rows is setup by the communication engine.
    IF (multrec_carrier%local_indexing) THEN
       k_locals => array_data (right%local_rows)
       multrec_carrier%k_locals => k_locals
       CALL ensure_array_size (multrec_carrier%k_sizes, ub=SIZE(k_locals), error=error)
       CALL local_filter(array_data(right%row_blk_size), SIZE(k_locals),&
            k_locals, multrec_carrier%k_sizes)
    ELSE
       k_locals => multrec_carrier%k_locals
    ENDIF
    ! Setup the block norms
    multrec_carrier%a_norms => a_norms
    multrec_carrier%b_norms => b_norms


    ! Start local multiplication
    t_a_f = 1
    t_a_l = left%nblks
    t_b_f = 1
    t_b_l = right%nblks
    !$ ithread = OMP_GET_THREAD_NUM()
    !$ t_a_f = left%thr_c(ithread+1)+1
    !$ t_a_l = left%thr_c(ithread+2)
    IF (multrec_carrier%local_indexing) THEN
       CALL sparse_multrec(&
            1, left%nblkrows_local,&
            1, right%nblkcols_local,&
            1, SIZE(k_locals),&
            t_a_f, t_a_l, left%coo_l,&
            t_b_f, t_b_l, right%coo_l,&
            product, multrec_carrier, error, 0)
    ELSE
       CALL sparse_multrec(&
            1, left%nblkrows_total,&
            1, right%nblkcols_total,&
            1, SIZE(k_locals),&
            t_a_f, t_a_l, left%coo_l,&
            t_b_f, t_b_l, right%coo_l,&
            product, multrec_carrier, error, 0)
    ENDIF

    CALL dbcsr_mm_csr_inner_finalize(csr_carrier=multrec_carrier%csr_carrier, error=error)
    flop = multrec_carrier%flop
    !
    IF (ASSOCIATED (multrec_carrier%k_sizes) .AND. multrec_carrier%local_indexing) &
         DEALLOCATE (multrec_carrier%k_sizes)

  END SUBROUTINE dbcsr_mm_multrec_multiply


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
  SUBROUTINE dbcsr_mm_multrec_finalize(multrec_carrier, product, error)
    TYPE(multrec_carrier_type), &
      INTENT(inout)                          :: multrec_carrier
    TYPE(dbcsr_type), INTENT(inout)          :: product
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_multrec_finalize', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler, ithread

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_mm_csr_outer_finalize(csr_carrier=multrec_carrier%csr_carrier,&
          product=product, error=error)

    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM()

    ! Release the carrier
    DEALLOCATE (multrec_carrier%m_sizes)
    DEALLOCATE (multrec_carrier%n_sizes)
    DEALLOCATE (multrec_carrier%row_max_epss)

    IF (multrec_carrier%local_indexing) THEN
       CALL remap_local2global(product%wms(ithread+1)%row_i, &
                               product%wms(ithread+1)%col_i, &
                  multrec_carrier%c_local_rows, multrec_carrier%c_local_cols,&
                  multrec_carrier%original_lastblk+1, product%wms(ithread+1)%lastblk)
    ENDIF
    !

    IF (.NOT. multrec_carrier%local_indexing) THEN
       DEALLOCATE (multrec_carrier%c_local_rows)
       DEALLOCATE (multrec_carrier%c_local_cols)
       DEALLOCATE (multrec_carrier%k_locals)
       NULLIFY (multrec_carrier%k_sizes)
    ENDIF

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_mm_multrec_finalize




! *****************************************************************************
!> \brief Performs recursive multiplication
!> \author Joost VandeVondele
! *****************************************************************************
  RECURSIVE SUBROUTINE sparse_multrec(mi, mf, ni, nf, ki, kf,&
       ai, af, a_index, bi, bf, b_index, &
       product, carrier, error, d)
    INTEGER, INTENT(IN)                      :: mi, mf, ni, nf, ki, kf, ai, af
    INTEGER, DIMENSION(3, 1:af), INTENT(IN)  :: a_index
    INTEGER, INTENT(IN)                      :: bi, bf
    INTEGER, DIMENSION(3, 1:bf), INTENT(IN)  :: b_index
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(multrec_carrier_type), &
      INTENT(INOUT)                          :: carrier
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
        IF (af-ai+1 .GT. 0 .AND. bf-bi+1 .GT. 0) &
        CALL dbcsr_mm_csr_multiply(mi=mi, mf=mf,ni=ni, nf=nf, ki=ki, kf=kf,&
           ai=ai, af=af, &
           bi=bi, bf=bf,&
           m_sizes=carrier%m_sizes, n_sizes=carrier%n_sizes, k_sizes=carrier%k_sizes,&
           c_local_rows=carrier%c_local_rows, c_local_cols=carrier%c_local_cols,&
           c_has_symmetry=carrier%c_has_symmetry, keep_sparsity=carrier%keep_sparsity,&
           use_eps=carrier%use_eps, row_max_epss=carrier%row_max_epss,&
           flop=carrier%flop,&
           a_index=a_index, b_index=b_index,&
           a_norms=carrier%a_norms, b_norms=carrier%b_norms,&  
           csr_carrier=carrier%csr_carrier, &
           product=product,&
           error=error)
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
            ai,acut-1,a_index, bi,bf,b_index, product, carrier, error,d+1)
       CALL sparse_multrec(mi+s1,mf, ni,nf, ki,kf,&
            acut,af,a_index, bi,bf,b_index,  product, carrier, error,d+1)
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
            ai,acut-1,a_index, bi,bcut-1,b_index, product, carrier, error,d+1)
       CALL sparse_multrec(mi,mf, ni,nf, ki+s1,kf,&
            acut,af,a_index, bcut,bf,b_index, product, carrier, error,d+1)
    CASE(3)
       s1=N/2
       bcut = find_cut_col(bi,bf,b_index,ni+s1-1)
       IF (dbg) THEN
          WRITE(*,*)N,s1,ni+s1-1,"/",bi,bf,bcut
          WRITE(*,'(3(I7))')b_index
       ENDIF
       CALL sparse_multrec(mi,mf, ni,ni+s1-1, ki,kf,&
            ai,af,a_index, bi,bcut-1,b_index, product, carrier, error,d+1)
       CALL sparse_multrec(mi,mf, ni+s1,nf, ki,kf,&
            ai,af,a_index, bcut,bf,b_index,  product, carrier, error,d+1)
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
