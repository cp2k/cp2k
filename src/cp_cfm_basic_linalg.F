!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_cfm_basic_linalg [1.0] *
!!
!!   NAME
!!     cp_cfm_basic_linalg
!!
!!   FUNCTION
!!     basic linear algebra operations for complex full matrixes
!!
!!   NOTES
!!     - not all functionality implemented
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!     Nearly literal copy of Fawzi's routines
!!
!!   SOURCE
!****************************************************************************
MODULE cp_cfm_basic_linalg
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_equivalent,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_retain,&
                                             cp_fm_struct_type
  USE cp_cfm_types,                     ONLY: cp_cfm_create,&
                                             cp_cfm_get_info,&
                                             cp_cfm_init_random,&
                                             cp_cfm_release,&
                                             cp_cfm_retain,&
                                             cp_cfm_set_all,&
                                             cp_cfm_set_element,&
                                             cp_cfm_get_element,&
                                             cp_cfm_p_type,&
                                             cp_cfm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_from_globenv,&
                                             cp_para_env_release,&
                                             cp_para_env_retain
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: wp=>dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync
  USE parallel,                        ONLY: cp2k_is_parallel
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_cfm_basic_linalg'

  PUBLIC :: cp_cfm_add, cp_cfm_gemm, cp_cfm_lu_decompose
!!***
!****************************************************************************
CONTAINS

  SUBROUTINE cp_cfm_add(alpha,matrix_a,beta,matrix_b,error)

    !   Purpose: Scale and add two BLACS matrices (a <- alpha*a + beta*b).

    !   History: - Creation (11.06.2001, Matthias Krack)

    !   ***************************************************************************

    TYPE(cp_cfm_type), POINTER           :: matrix_a
    COMPLEX(wp), INTENT(IN)                         :: alpha
    COMPLEX(wp), INTENT(in), OPTIONAL               :: beta
    TYPE(cp_cfm_type), POINTER, OPTIONAL :: matrix_b
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    !   *** Local variables ***

    INTEGER :: handle,mypcol,myprow,npcol,nprow,source

    COMPLEX(wp), DIMENSION(:,:), POINTER :: a,b
    CHARACTER(len=*), PARAMETER :: routineN='cp_cfm_add',&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    COMPLEX(wp) :: my_beta

    !   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    my_beta=CMPLX(0.0_wp,0.0_wp)
    IF(PRESENT(beta)) my_beta=beta
    NULLIFY(a,b)

    CPPrecondition(ASSOCIATED(matrix_a),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix_a%ref_count>0,cp_failure_level,routineP,error,failure)
    ! to do: use dscal,dcopy,daxp
    myprow=matrix_a%matrix_struct%context%mepos(1)
    mypcol=matrix_a%matrix_struct%context%mepos(2)

    a => matrix_a%local_data

    IF (my_beta==0.0_wp) THEN

       IF (alpha==CMPLX(0.0_wp,0.0_wp)) THEN
          a(:,:)=CMPLX(0.0_wp,0.0_wp)
       ELSE IF (alpha == CMPLX(1.0_wp,0.0_wp)) THEN
          RETURN
       ELSE
          a(:,:) = alpha*a(:,:)
       END IF

    ELSE
       CPPrecondition(PRESENT(matrix_b),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(matrix_b),cp_failure_level,routineP,error,failure)
       CPPrecondition(matrix_b%ref_count>0,cp_failure_level,routineP,error,failure)
       CALL cp_assert(matrix_a%matrix_struct%context%group==&
            matrix_b%matrix_struct%context%group,cp_failure_level,&
            cp_assertion_failed,fromWhere=routineP,&
            message="matrixes must be in the same blacs context"//&
            CPSourceFileRef,&
            error=error,failure=failure)

       IF (cp_fm_struct_equivalent(matrix_a%matrix_struct,&
            matrix_b%matrix_struct,error=error)) THEN

          b => matrix_b%local_data

          IF (alpha == CMPLX(0.0_wp,0.0_wp)) THEN
             IF (my_beta == CMPLX(1.0_wp,0.0_wp)) THEN
                a(:,:) = b(:,:)
             ELSE
                a(:,:) = my_beta*b(:,:)
             END IF
          ELSE IF (alpha == CMPLX(1.0_wp,0.0_wp)) THEN
             IF (my_beta == CMPLX(1.0_wp,0.0_wp)) THEN
                a(:,:) = a(:,:) + b(:,:)
             ELSE
                a(:,:) = a(:,:) + my_beta*b(:,:)
             END IF
          ELSE
             a(:,:) = alpha*a(:,:) + my_beta*b(:,:)
          END IF
       ELSE
#ifdef __SCALAPACK
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="to do (pdscal,pdcopy,pdaxpy)", error=error)
#else
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
#endif
       END IF
    END IF
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE cp_cfm_add

! *****************************************************************************

! *****************************************************************************
! computes the LU decomposition of a given matrix
! the actual purpose right now is to compute the determinant of a given matrix
! which is most efficiently done this way, but, indeed, destroys the matrix
! SERIOUS WARNING (KNOWN BUG) : the sign of the determinant depends on ipivot
! one should be able to find out if ipivot is an even or an odd permutation...
! *****************************************************************************
  SUBROUTINE cp_cfm_lu_decompose(matrix_a,almost_determinant)
    TYPE(cp_cfm_type), POINTER :: matrix_a
    COMPLEX(wp), INTENT(OUT)              :: almost_determinant

! *** locals ***
    INTEGER :: handle,lda,i,n,istat,info,npe
    INTEGER, DIMENSION(9) :: desca
    COMPLEX(wp), DIMENSION(:,:), POINTER :: a
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ipivot
    CHARACTER(LEN=*), PARAMETER :: routine =&
         "SUBROUTINE cp_cfm_lu_decompose"
    COMPLEX(wp), DIMENSION(:), POINTER :: diag
    COMPLEX(wp) :: determinant

    CALL timeset("cp_cfm_lu_decompose","I","",handle)

    a => matrix_a%local_data
    n = matrix_a%matrix_struct%nrow_global
    ALLOCATE(ipivot(n))

#if defined(__SCALAPACK)
    desca(:) = matrix_a%matrix_struct%descriptor(:)
    CALL pzgetrf(n,n,a(1,1),1,1,desca,ipivot,info)
    ALLOCATE(diag(n))
    diag(:)=CMPLX(0.0_wp,0.0_wp)
    DO i=1,n
       CALL cp_cfm_get_element(matrix_a,i,i,diag(i)) !  not completely optimal in speed i would say
    ENDDO
    determinant=CMPLX(1.0_wp,0.0_wp)
    do i=1,n
       determinant=determinant*diag(i)
    enddo
    DEALLOCATE(diag)
#else
    lda=size(a,1)
    CALL zgetrf(n,n,a(1,1),lda,ipivot,info)
    determinant=CMPLX(1.0_wp,0.0_wp)
    do i=1,n
       determinant=determinant*a(i,i)
    enddo
#endif
    ! info is allowed to be zero
    ! this does just signal a zero diagonal element
    DEALLOCATE(ipivot)
    almost_determinant=determinant ! notice that the sign is random
    CALL timestop(0.0_wp,handle)
  END SUBROUTINE


! *****************************************************************************

  SUBROUTINE cp_cfm_gemm(transa,transb,m,n,k,alpha,matrix_a,matrix_b,beta,&
       matrix_c,error,b_first_col)

!   Purpose: BLACS interface to the BLAS routine dgemm.

!   History: - Creation (07.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_cfm_type), POINTER:: matrix_a,matrix_b,matrix_c
    CHARACTER(LEN=1), INTENT(IN)              :: transa,transb
    COMPLEX(wp), INTENT(IN)                      :: alpha,beta
    INTEGER, INTENT(IN)                       :: k,m,n
    INTEGER, OPTIONAL, INTENT(IN)             :: b_first_col
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: handle,j,lda,ldb,ldc

    INTEGER, DIMENSION(9) :: desca,descb,descc

    COMPLEX(wp), DIMENSION(:,:), POINTER :: a,b,c
    REAL(wp) :: flops

!   ---------------------------------------------------------------------------

    CALL timeset("cp_cfm_gemm","I","Gflops",handle)

    a => matrix_a%local_data
    b => matrix_b%local_data
    c => matrix_c%local_data

    IF (PRESENT(b_first_col)) THEN
      j = b_first_col
    ELSE
      j = 1
    END IF

#if defined(__SCALAPACK)

    desca(:) = matrix_a%matrix_struct%descriptor(:)
    descb(:) = matrix_b%matrix_struct%descriptor(:)
    descc(:) = matrix_c%matrix_struct%descriptor(:)

    flops=4*2.0_wp*(1E-9_wp*m*n*k)/(matrix_a%matrix_struct%context%num_pe(1)*&
         matrix_a%matrix_struct%context%num_pe(2))
    CALL pzgemm(transa,transb,m,n,k,alpha,a(1,1),1,1,desca,b(1,1),1,j,&
                descb,beta,c(1,1),1,j,descc)
#else

    lda = SIZE(a,1)
    ldb = SIZE(b,1)
    ldc = SIZE(c,1)

    flops=4*2.0_wp*(1E-9_wp*m*n*k)
    CALL zgemm(transa,transb,m,n,k,alpha,a(1,1),lda,b(1,j),ldb,beta,c(1,j),ldc)

#endif
    CALL timestop(flops,handle)

  END SUBROUTINE cp_cfm_gemm

END MODULE cp_cfm_basic_linalg
