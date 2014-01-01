!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Processes MM stack and issues BLAS xGEMM calls
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[in,out] error        error
! *****************************************************************************
  SUBROUTINE blas_process_mm_stack_z(params,&
       stack_size,&
       a_data, b_data, c_data, error)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data
    TYPE(dbcsr_error_type), INTENT(inout)     :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'blas_process_mm_stack_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL ZGEMM('N',&
            'N',&
            params(p_m,sp), params(p_n,sp),& !m, n
            params(p_k,sp),& ! k
            CMPLX(1.0, 0.0, real_8),& ! alpha
            a_data(params(p_a_first,sp)),& ! A
            params(p_m,sp),& !lda
            b_data(params(p_b_first,sp)),& ! B
            params(p_k,sp),& !ldb
            CMPLX(1.0, 0.0, real_8),& ! beta
            c_data(params(p_c_first,sp)), params(p_m,sp))
    ENDDO
  END SUBROUTINE blas_process_mm_stack_z

! *****************************************************************************
!> \brief Processes MM stack and issues internal MM calls.
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[in,out] error        error
! *****************************************************************************
  SUBROUTINE internal_process_mm_stack_z(params, stack_size,&
       a_data, b_data, c_data, error)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data
    TYPE(dbcsr_error_type), INTENT(inout)     :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'internal_process_mm_stack_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL internal_mm_z_nn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO
  END SUBROUTINE internal_process_mm_stack_z


! *****************************************************************************
!> \brief Processes MM stack and issues SMM library calls
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[in,out] error        error
! *****************************************************************************
  SUBROUTINE smm_process_mm_stack_z(stack_descr, params,&
       stack_size,&
       a_data, b_data, c_data, error)
    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data
    TYPE(dbcsr_error_type), INTENT(inout)     :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'smm_process_mm_stack_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------
#if defined(__HAS_smm_vec) && defined(__HAS_smm_znn)

    IF(stack_descr%defined_mnk) THEN
       CALL smm_vec_znn (stack_descr%m, stack_descr%n, stack_descr%k, &
         a_data, b_data, c_data, stack_size, &
         dbcsr_ps_width, params, p_a_first, p_b_first, p_c_first)
       RETURN
    ENDIF

#endif

    DO sp = 1, stack_size
       CALL smm_znn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO

  END SUBROUTINE smm_process_mm_stack_z


! *****************************************************************************
!> \brief Processes MM stack and issues Plasma xGEMM calls.
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[in,out] error        error
! *****************************************************************************
  SUBROUTINE plasma_process_mm_stack_z(params, stack_size,&
       a_data, b_data, c_data, error)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data
    TYPE(dbcsr_error_type), INTENT(inout)     :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'plasma_process_mm_stack_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------
#ifdef __PLASMA
    INCLUDE 'plasmaf.h'
#else
    CALL dbcsr_assert(.FALSE.,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "PLASMA support not compiled.", __LINE__, error=error)
#endif
    !
    DO sp = 1, stack_size
#ifdef __PLASMA
       CALL plasma_ZGEMM(&
            'N',&
            'N',&
            params(p_m,sp), params(p_n,sp),& !m, n
            params(p_k,sp),& ! k
            CMPLX(1.0, 0.0, real_8),& ! alpha
            a_data(params(p_a_first,sp)),& ! A
            params(p_m,sp),& !lda
            b_data(params(p_b_first,sp)),& ! B
            params(p_k,sp),& !ldb
            CMPLX(1.0, 0.0, real_8),& ! beta
            c_data(params(p_c_first,sp)), params(p_m,sp),& !c, ldc
            plasma_info)
       CALL dbcsr_assert( plasma_info, "EQ", 0, dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "plasma_gemm failed", __LINE__, error=error)
#else
       CALL dbcsr_assert( .FALSE., dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "plasma badly set", __LINE__, error=error)
#endif
    ENDDO
  END SUBROUTINE plasma_process_mm_stack_z



  PURE SUBROUTINE internal_mm_z_nn(&
       M,N,K,A,B,C)
    INTEGER, INTENT(IN)                      :: M, N, K
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: C(M,N)
    COMPLEX(kind=real_8), INTENT(IN)                      :: B(K,N)
    COMPLEX(kind=real_8), INTENT(IN)                      :: A(M,K)
    C(:,:) = C(:,:) + MATMUL (A, B)
  END SUBROUTINE internal_mm_z_nn
