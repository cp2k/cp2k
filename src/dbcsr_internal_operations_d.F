
! *****************************************************************************
!> \brief Issues actual DGEMM calls.
!>
!> \param[in] params           Stack of DGEMM parameters
!> \param[in] n                Number of parameters
!> \param[in] left_data_a      Left-matrix data
!> \param[in] right_data_a     Right-matrix data
!> \param[in,out] product_data_area  Data for results
!> \param[out] lflop           (optional) Number of FLOPs used by DGEMM
! *****************************************************************************
  SUBROUTINE process_dgemm_stack_d(params,&
       params_a, params_b, params_c,&
       stack_size, stack_size_a, stack_size_b, stack_size_c,&
       left_data_a, right_data_a, product_data_area, lflop, ltime, error)
    INTEGER, INTENT(IN)                       :: stack_size, stack_size_a,&
                                                 stack_size_b, stack_size_c
    TYPE(dgemm_join), DIMENSION(1:stack_size), &
      INTENT(IN)                              :: params
    TYPE(block_parameters), DIMENSION(1:stack_size_a), INTENT(IN) :: params_a
    TYPE(block_parameters), DIMENSION(1:stack_size_b), INTENT(IN) :: params_b
    TYPE(block_parameters), DIMENSION(1:stack_size_c), INTENT(IN) :: params_c
    REAL(kind=real_8), DIMENSION(:), INTENT(IN)         :: left_data_a, &
                                                 right_data_a
    TYPE(dbcsr_data_obj), INTENT(INOUT)       :: product_data_area
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                                :: lflop
    REAL(KIND=real_8), INTENT(OUT), OPTIONAL  :: ltime
    TYPE(dbcsr_error_type), INTENT(inout)     :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'process_dgemm_stack_d', &
      routineP = moduleN//':'//routineN
    REAL, PARAMETER                           :: resize_factor = 1.618034

    INTEGER                                   :: c, lda, ldb, ldc, maxs, r, sp,&
                                                 poff
    INTEGER                                   :: sp_c, sp_a, sp_b
    LOGICAL                                   :: do_resize

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY

    IF (PRESENT (lflop)) lflop = INT(0, int_8)
    IF (PRESENT (ltime)) ltime = 0.0_real_8
    ! Increase product data area size if necessary.
    do_resize = .FALSE.
    maxs = dbcsr_get_data_size(product_data_area)
    DO sp = 1, stack_size
       sp_c = params(sp)%p_c
       IF (params_c(sp_c)%offset+params_c(sp_c)%nze-1 .GT. maxs) THEN
          maxs = params_c(sp_c)%offset+params_c(sp_c)%nze-1
          do_resize = .TRUE.
       ENDIF
    ENDDO
    IF (do_resize) THEN
       CALL dbcsr_data_ensure_size (product_data_area,&
            maxs, factor=resize_factor, error=error)
    ENDIF
    ! Setup encapsulated data area
    DO sp = 1, stack_size
       sp_c = params(sp)%p_c
       sp_a = params(sp)%p_a
       sp_b = params(sp)%p_b
       IF (params(sp)%last_k .EQ. 0) CYCLE
       IF (params_c(sp_c)%nze .EQ. 0&
            .or.params_a(sp_a)%nze .EQ. 0&
            .or.params_b(sp_b)%nze .EQ. 0) THEN
          CYCLE
       ENDIF
       !
       poff = params_c(sp_c)%offset
       !
       ! Avoid nagging compilers
       IF (params(sp)%beta%r_dp .EQ. 0.0_real_8) THEN
          IF (.NOT. params_c(sp_c)%tr &
               .OR. params(sp)%last_n .EQ. params_c(sp_c)%logical_cols) THEN
             CALL dbcsr_data_clear (product_data_area,&
                  lb=poff,&
                  ub=poff-1+params_c(sp_c)%logical_rows*params(sp)%last_n)
          ELSE
             ! It must be transposed.
             !### optimize out the inner loop (turn into a range)
             DO r = 1, params_c(sp_c)%logical_rows
                DO c = 1, params(sp)%last_n
                   CALL dbcsr_data_clear(product_data_area,&
                        lb=poff-1+(r-1)*params_c(sp_c)%logical_cols+c,&
                        ub=poff-1+(r-1)*params_c(sp_c)%logical_cols+c)
                ENDDO
             ENDDO
             !FORALL (r = 1:params(sp)%lrows_c, c = 1:params(sp)%last_n)
             !   product_data_d((r-1)*params(sp)%lcols_c+c) = 0.0_dp
             !END FORALL
          ENDIF
       ENDIF
       ! Now prepare and issue the multiplies.
       lda = select_n_or_t (params_a(sp_a)%tr,&
            params_a(sp_a)%logical_rows, params_a(sp_a)%logical_cols)
       ldb = select_n_or_t (params_b(sp_b)%tr,&
            params_b(sp_b)%logical_rows, params_b(sp_b)%logical_cols)
       ldc = select_n_or_t (params_c(sp_c)%tr,&
            params_c(sp_c)%logical_rows, params_c(sp_c)%logical_cols)
       !write(*,*)"LDx", lda, ldb, ldc
       IF (.NOT. params_c(sp_c)%tr) THEN
          IF (PRESENT (ltime)) ltime = ltime + m_walltime()
          CALL DGEMM(&
               blas_mat_type(params_a(sp_a)%tr),&
               blas_mat_type(params_b(sp_b)%tr),&
               params_c(sp_c)%logical_rows, params(sp)%last_n,& !m, n
               params(sp)%last_k,& ! k
               params(sp)%alpha%r_dp,&
               left_data_a(params_a(sp_a)%offset), lda,&
               right_data_a(params_b(sp_b)%offset), ldb,&
               params(sp)%beta%r_dp,&
               product_data_area%d%r_dp(poff), ldc)
          IF (PRESENT (ltime)) ltime = m_walltime() - ltime
          IF (PRESENT (lflop)) &
               lflop = lflop + INT(2, int_8) * &
               INT(params_c(sp_c)%logical_rows &
               * params(sp)%last_n * params(sp)%last_k, int_8)
       ELSE
          IF (PRESENT (ltime)) ltime = ltime + m_walltime()
          CALL DGEMM(&
               blas_mat_type(.NOT.params_b(sp_b)%tr),&
               blas_mat_type(.NOT.params_a(sp_a)%tr),&
               params_c(sp_c)%logical_cols, params_c(sp_c)%logical_rows,& !m, n (switched)
               params(sp)%last_k,& ! k
               params(sp)%alpha%r_dp,&
               right_data_a(params_b(sp_b)%offset), ldb,&
               left_data_a(params_a(sp_a)%offset), lda,&
               params(sp)%beta%r_dp,&
               product_data_area%d%r_dp(poff), ldc)
          IF (PRESENT (ltime)) ltime = m_walltime() - ltime
          IF (PRESENT (lflop)) &
               lflop = lflop + INT(2, int_8) * &
               INT(params(sp)%last_n &
               * params_c(sp_c)%logical_cols * params(sp)%last_k, int_8)
       ENDIF
    ENDDO
  END SUBROUTINE process_dgemm_stack_d


