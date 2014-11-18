!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Collection of simple mathematical functions and subroutines
!> \par History
!>      FUNCTION angle updated and FUNCTION dihedral angle added; cleaned
!>      (13.03.2004,MK)
!> \author MK (15.11.1998)
! *****************************************************************************
MODULE mathlib

  
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: euler,&
                                             fac
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mathlib'
  REAL(KIND=dp), PARAMETER             :: eps_geo = 1.0E-6_dp

  ! Public subroutines

  PUBLIC :: build_rotmat,&
            jacobi,&
            diamat_all,&
            invmat,&
            invert_matrix,&
            symmetrize_matrix,&
            unit_matrix, diag

  ! Public functions

  PUBLIC :: angle,&
            binomial,&
            det_3x3,&
            dihedral_angle,&
            gcd,&
            inv_3x3,&
            lcm,&
            vector_product,&
            matmul_3x3,&
            matvec_3x3,&
            rotate_vector,&
            reflect_vector,&
            dotprod_3d,&
            transpose_3d,&
            expint

  INTERFACE det_3x3
    MODULE PROCEDURE det_3x3_1,det_3x3_2
  END INTERFACE

  INTERFACE invert_matrix
    MODULE PROCEDURE invert_matrix_d,invert_matrix_z
  END INTERFACE

  INTERFACE set_diag
    MODULE PROCEDURE set_diag_scalar_d,set_diag_scalar_z
  END INTERFACE

  INTERFACE swap
    MODULE PROCEDURE swap_scalar,swap_vector
  END INTERFACE

  INTERFACE unit_matrix
    MODULE PROCEDURE unit_matrix_d,unit_matrix_z
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief  Calculation of the angle between the vectors a and b.
!>         The angle is returned in radians.
!> \param a ...
!> \param b ...
!> \retval angle_ab ...
!> \date    14.10.1998
!> \author  MK
!> \version 1.0
! *****************************************************************************
  FUNCTION angle(a,b) RESULT(angle_ab)

    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a, b
    REAL(KIND=dp)                            :: angle_ab

    REAL(KIND=dp)                            :: length_of_a, length_of_b
    REAL(KIND=dp), DIMENSION(SIZE(a, 1))     :: a_norm, b_norm

    length_of_a = SQRT(DOT_PRODUCT(a,a))
    length_of_b = SQRT(DOT_PRODUCT(b,b))

    IF ((length_of_a > eps_geo).AND.(length_of_b > eps_geo)) THEN
       a_norm(:) = a(:)/length_of_a
       b_norm(:) = b(:)/length_of_b
       angle_ab = ACOS(MIN(MAX(DOT_PRODUCT(a_norm,b_norm),-1.0_dp),1.0_dp))
    ELSE
       angle_ab = 0.0_dp
    END IF

  END FUNCTION angle

! *****************************************************************************
!> \brief   The binomial coefficient n over k for 0 <= k <= n is calculated,
!>            otherwise zero is returned.
!> \param n ...
!> \param k ...
!> \retval n_over_k ...
!> \date    08.03.1999
!> \author  MK
!> \version 1.0
! *****************************************************************************
  FUNCTION binomial(n,k) RESULT(n_over_k)

    INTEGER, INTENT(IN)                      :: n, k
    REAL(KIND=dp)                            :: n_over_k

    IF ((k >= 0).AND.(k <= n)) THEN
       n_over_k = fac(n)/(fac(n-k)*fac(k))
    ELSE
       n_over_k = 0.0_dp
    END IF

  END FUNCTION binomial

! *****************************************************************************
!> \brief   The rotation matrix rotmat which rotates a vector about a
!>          rotation axis defined by the vector a is build up.
!>          The rotation angle is phi (radians).
!> \param phi ...
!> \param a ...
!> \param rotmat ...
!> \date    16.10.1998
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE build_rotmat(phi,a,rotmat)
    REAL(KIND=dp), INTENT(IN)                :: phi
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT)                            :: rotmat

    REAL(KIND=dp)                            :: cosp, cost, length_of_a, sinp
    REAL(KIND=dp), DIMENSION(3)              :: d

    length_of_a = SQRT(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
    ! Check the length of the vector a
    IF (length_of_a > eps_geo) THEN

       d(:) = a(:)/length_of_a

       cosp = COS(phi)
       sinp = SIN(phi)
       cost = 1.0_dp - cosp

       rotmat(1,1) = d(1)*d(1)*cost + cosp
       rotmat(1,2) = d(1)*d(2)*cost - d(3)*sinp
       rotmat(1,3) = d(1)*d(3)*cost + d(2)*sinp
       rotmat(2,1) = d(2)*d(1)*cost + d(3)*sinp
       rotmat(2,2) = d(2)*d(2)*cost + cosp
       rotmat(2,3) = d(2)*d(3)*cost - d(1)*sinp
       rotmat(3,1) = d(3)*d(1)*cost - d(2)*sinp
       rotmat(3,2) = d(3)*d(2)*cost + d(1)*sinp
       rotmat(3,3) = d(3)*d(3)*cost + cosp
    ELSE
       CALL unit_matrix(rotmat)
    END IF

  END SUBROUTINE build_rotmat

! *****************************************************************************
!> \brief   Returns the determinante of the 3x3 matrix a.
!> \param a ...
!> \retval det_a ...
!> \date    13.03.2004
!> \author  MK
!> \version 1.0
! *****************************************************************************
  FUNCTION det_3x3_1(a) RESULT(det_a)
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: a
    REAL(KIND=dp)                            :: det_a

    det_a = a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) +&
            a(1,2)*(a(2,3)*a(3,1) - a(2,1)*a(3,3)) +&
            a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))

  END FUNCTION det_3x3_1

! *****************************************************************************
!> \brief   Returns the determinante of the 3x3 matrix a given by its columns.
!> \param a1 ...
!> \param a2 ...
!> \param a3 ...
!> \retval det_a ...
!> \date    13.03.2004
!> \author  MK
!> \version 1.0
! *****************************************************************************
  FUNCTION det_3x3_2(a1,a2,a3) RESULT(det_a)
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a1, a2, a3
    REAL(KIND=dp)                            :: det_a

    det_a = a1(1)*(a2(2)*a3(3) - a3(2)*a2(3)) +&
            a2(1)*(a3(2)*a1(3) - a1(2)*a3(3)) +&
            a3(1)*(a1(2)*a2(3) - a2(2)*a1(3))

  END FUNCTION det_3x3_2

! *****************************************************************************
!> \brief Diagonalize the symmetric n by n matrix a using the LAPACK
!>        library. Only the upper triangle of matrix a is used.
!>        Externals (LAPACK 3.0)
!> \param a ...
!> \param eigval ...
!> \param dac ...
!> \param error ...
!> \date    29.03.1999
!> \par Variables
!>      - a       : Symmetric matrix to be diagonalized (input; upper triangle) ->
!>      -           eigenvectors of the matrix a (output).
!>      - dac     : If true, then the divide-and-conquer algorithm is applied.
!>      - eigval  : Eigenvalues of the matrix a (output).
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE diamat_all(a,eigval,dac,error)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigval
    LOGICAL, INTENT(IN), OPTIONAL            :: dac
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'diamat_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info, liwork, lwork, &
                                                n, nb, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iwork
    INTEGER, EXTERNAL                        :: ilaenv
    LOGICAL                                  :: divide_and_conquer, failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: work

    EXTERNAL dsyev,dsyevd

    failure = .FALSE.
    CALL timeset(routineN,handle)

    ! Get the size of the matrix a
    n = SIZE(a,1)

    ! Check the size of matrix a
    IF (SIZE(a,2) /= n) THEN
       CALL stop_program(moduleN,routineN,__LINE__,&
            "Check the size of matrix a (parameter #1)")
    END IF

    ! Check the size of vector eigval
    IF (SIZE(eigval) /= n) THEN
       CALL stop_program(moduleN,routineN,__LINE__,&
            "The dimension of vector eigval is too small")
    END IF

    ! Check, if the divide-and-conquer algorithm is requested

    IF (PRESENT(dac)) THEN
       divide_and_conquer = dac
    ELSE
       divide_and_conquer = .FALSE.
    END IF

    ! Get the optimal work storage size

    IF (divide_and_conquer) THEN
       lwork = 2*n**2 + 6*n + 1
       liwork = 5*n + 3
    ELSE
       nb = ilaenv(1,"DSYTRD","U",n,-1,-1,-1)
       lwork = (nb + 2)*n
    END IF

    ! Allocate work storage

    ALLOCATE (work(lwork),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (divide_and_conquer) THEN
       ALLOCATE (iwork(liwork),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    ! Diagonalize the matrix a

    IF (divide_and_conquer) THEN
       CALL dsyevd("V","U",n,a,n,eigval,work,lwork,iwork,liwork,info)
    ELSE
       CALL dsyev("V","U",n,a,n,eigval,work,lwork,info)
    END IF

    IF (info /= 0) THEN
       IF (divide_and_conquer) THEN
          CALL stop_program(moduleN,routineN,__LINE__,&
               "The matrix diagonalization with dsyevd failed")
       ELSE
          CALL stop_program(moduleN,routineN,__LINE__,&
               "The matrix diagonalization with dsyev failed")
       END IF
    END IF

    ! Release work storage

    DEALLOCATE (work,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (divide_and_conquer) THEN
       DEALLOCATE (iwork,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE diamat_all

! *****************************************************************************
!> \brief   Returns the dihedral angle, i.e. the angle between the planes
!>          defined by the vectors (-ab,bc) and (cd,-bc).
!>          The dihedral angle is returned in radians.
!> \param ab ...
!> \param bc ...
!> \param cd ...
!> \retval dihedral_angle_abcd ...
!> \date    13.03.2004
!> \author  MK
!> \version 1.0
! *****************************************************************************
  FUNCTION dihedral_angle(ab,bc,cd) RESULT(dihedral_angle_abcd)
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ab, bc, cd
    REAL(KIND=dp)                            :: dihedral_angle_abcd

    REAL(KIND=dp)                            :: det_abcd
    REAL(KIND=dp), DIMENSION(3)              :: abc, bcd

    abc = vector_product(bc,-ab)
    bcd = vector_product(cd,-bc)
    ! Calculate the normal vectors of the planes
    ! defined by the points a,b,c and b,c,d

    det_abcd = det_3x3(abc,bcd,-bc)
    dihedral_angle_abcd = SIGN(1.0_dp,det_abcd)*angle(abc,bcd)

  END FUNCTION dihedral_angle

! *****************************************************************************
!> \brief   Return the diagonal elements of matrix a as a vector.
!> \param a ...
!> \retval a_diag ...
!> \date    20.11.1998
!> \author  MK
!> \version 1.0
! *****************************************************************************
  FUNCTION get_diag(a) RESULT(a_diag)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: a
    REAL(KIND=dp), &
      DIMENSION(MIN(SIZE(a, 1), SIZE(a, 2))) :: a_diag

    INTEGER                                  :: i, n

    n = MIN(SIZE(a,1),SIZE(a,2))

    DO i=1,n
       a_diag(i) = a(i,i)
    END DO

  END FUNCTION get_diag

! *****************************************************************************
!> \brief   Returns the inverse of the 3 x 3 matrix a.
!> \param a ...
!> \retval a_inv ...
!> \date    13.03.2004
!> \author  MK
!> \version 1.0
! *****************************************************************************
  FUNCTION inv_3x3(a) RESULT(a_inv)

    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: a
    REAL(KIND=dp), DIMENSION(3, 3)           :: a_inv

    REAL(KIND=dp)                            :: det_a

    det_a = 1.0_dp/det_3x3(a)

    a_inv(1,1) = (a(2,2)*a(3,3) - a(3,2)*a(2,3))*det_a
    a_inv(2,1) = (a(2,3)*a(3,1) - a(3,3)*a(2,1))*det_a
    a_inv(3,1) = (a(2,1)*a(3,2) - a(3,1)*a(2,2))*det_a

    a_inv(1,2) = (a(1,3)*a(3,2) - a(3,3)*a(1,2))*det_a
    a_inv(2,2) = (a(1,1)*a(3,3) - a(3,1)*a(1,3))*det_a
    a_inv(3,2) = (a(1,2)*a(3,1) - a(3,2)*a(1,1))*det_a

    a_inv(1,3) = (a(1,2)*a(2,3) - a(2,2)*a(1,3))*det_a
    a_inv(2,3) = (a(1,3)*a(2,1) - a(2,3)*a(1,1))*det_a
    a_inv(3,3) = (a(1,1)*a(2,2) - a(2,1)*a(1,2))*det_a

  END FUNCTION inv_3x3

! *****************************************************************************
!> \brief returns inverse of matrix using the lapack routines DGETRF and DGETRI
!> \param a ...
!> \param info ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE invmat(a,info,error)
    REAL(KIND=dp), INTENT(INOUT)             :: a(:,:)
    INTEGER, INTENT(OUT)                     :: info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'invmat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lwork, n
    INTEGER, ALLOCATABLE                     :: ipiv(:)
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE               :: work(:)

    failure = .FALSE.
    n = SIZE(a,1)
    lwork = 20*n
    ALLOCATE (ipiv(n),STAT=info)
    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (work(lwork),STAT=info)
    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)
    ipiv = 0
    work  = 0._dp
    info = 0
    CALL dgetrf(n,n,a,n,ipiv,info)
    IF (info==0) THEN
       CALL dgetri(n,a,n,ipiv,work,lwork,info)
    END IF
    DEALLOCATE (ipiv,work,STAT=info)
    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE invmat

! *****************************************************************************
!> \brief  Compute the inverse of the n by n real matrix a using the LAPACK
!>         library
!> \param a ...
!> \param a_inverse ...
!> \param eval_error ...
!> \param option ...
!> \param error ...
!> \param improve ...
!> \date   23.03.1999
!> \par Variables
!>       - a        : Real matrix to be inverted (input).
!>       - a_inverse: Inverse of the matrix a (output).
!>       - a_lu     : LU factorization of matrix a.
!>       - a_norm   : Norm of matrix a.
!>       - error    : Estimated error of the inversion.
!>       - r_cond   : Reciprocal condition number of the matrix a.
!>       - trans    : "N" => invert a
!>       -            "T" => invert transpose(a)
!> \author MK
!> \version 1.0
!> \note NB add improve argument, used to disable call to dgerfs
! *****************************************************************************
  SUBROUTINE invert_matrix_d(a,a_inverse,eval_error,option,error,improve)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: a
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: a_inverse
    REAL(KIND=dp), INTENT(OUT)               :: eval_error
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL   :: option
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: improve

    CHARACTER(LEN=*), PARAMETER :: routineN = 'invert_matrix_d', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: norm, trans
    CHARACTER(LEN=default_string_length)     :: message
    INTEGER                                  :: info, iter, n, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipiv, iwork
    LOGICAL                                  :: do_improve, failure
    REAL(KIND=dp)                            :: a_norm, old_eval_error, r_cond
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: berr, ferr, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: a_lu, b
    REAL(KIND=dp), EXTERNAL                  :: dlange

    EXTERNAL dgecon,dgerfs,dgetrf,dgetrs

    failure = .FALSE.

    ! Check for optional parameter
    IF (PRESENT(option)) THEN
      trans = option
    ELSE
      trans = "N"
    END IF

    IF (PRESENT(improve)) THEN
      do_improve = improve
    ELSE
      do_improve = .TRUE.
    ENDIF

    ! Get the dimension of matrix a
    n = SIZE(a,1)

    ! Check array dimensions
    IF (n == 0) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Matrix to be inverted of zero size")
    END IF

    IF (n /= SIZE(a,2)) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Check the array bounds of parameter #1")
    END IF

    IF ((n /= SIZE(a_inverse,1)).OR.&
        (n /= SIZE(a_inverse,2))) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Check the array bounds of parameter #2")
    END IF

    ! Allocate work storage
    ALLOCATE (a_lu(n,n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (b(n,n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (berr(n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (ferr(n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (ipiv(n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (iwork(n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (work(4*n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    a_lu(1:n,1:n) = a(1:n,1:n)

    ! Compute the LU factorization of the matrix a
    CALL dgetrf(n,n,a_lu,n,ipiv,info)

    IF (info /= 0) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "The LU factorization in dgetrf failed")
    END IF

    ! Compute the norm of the matrix a

    IF (trans == "N") THEN
      norm = '1'
    ELSE
      norm = 'I'
    END IF

    a_norm = dlange(norm,n,n,a,n,work)

    ! Compute the reciprocal of the condition number of a

    CALL dgecon(norm,n,a_lu,n,a_norm,r_cond,work,iwork,info)

    IF (info /= 0) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "The computation of the condition number in "//&
                        "dgecon failed")
    END IF

    IF (r_cond < EPSILON(0.0_dp)) THEN
      WRITE (message,"(A,ES10.3)") "R_COND =",r_cond
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Bad condition number "//TRIM(message)//" (smaller than the machine "//&
                        "working precision)")
    END IF

    ! Solve a system of linear equations using the LU factorization computed by dgetrf

    CALL unit_matrix(a_inverse)

    CALL dgetrs(trans,n,n,a_lu,n,ipiv,a_inverse,n,info)

    IF (info /= 0) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Solving the system of linear equations in dgetrs "//&
                        "failed")
    END IF

    ! Improve the computed solution iteratively
    CALL unit_matrix(b) ! Initialize right-hand sides

    eval_error = 0.0_dp

    IF (do_improve) THEN
      DO iter=1,10

        CALL dgerfs(trans,n,n,a,n,a_lu,n,ipiv,b,n,a_inverse,n,ferr,berr,&
                    work,iwork,info)

        IF (info /= 0) THEN
          CALL stop_program(moduleN,routineN,__LINE__,&
                            "Improving the computed solution in dgerfs failed")
        END IF

        old_eval_error = eval_error
        eval_error = MAXVAL(ferr)

        IF (ABS(eval_error - old_eval_error) <= EPSILON(1.0_dp)) EXIT

      END DO
    ENDIF

    ! Release work storage
    DEALLOCATE (work,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (iwork,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (ipiv,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (ferr,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (berr,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (b,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (a_lu,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE invert_matrix_d

! *****************************************************************************
!> \brief  Compute the inverse of the n by n complex matrix a using the LAPACK
!>         library
!> \param a ...
!> \param a_inverse ...
!> \param eval_error ...
!> \param option ...
!> \param error ...
!> \date   08.06.2009
!> \par Variables
!>       - a        : Complex matrix to be inverted (input).
!>       - a_inverse: Inverse of the matrix a (output).
!>       - a_lu     : LU factorization of matrix a.
!>       - a_norm   : Norm of matrix a.
!>       - error    : Estimated error of the inversion.
!>       - r_cond   : Reciprocal condition number of the matrix a.
!>       - trans    : "N" => invert a
!>       -            "T" => invert transpose(a)
!> \author MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE invert_matrix_z(a,a_inverse,eval_error,option,error)
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: a
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: a_inverse
    REAL(KIND=dp), INTENT(OUT)               :: eval_error
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL   :: option
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'invert_matrix_z', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: norm, trans
    CHARACTER(LEN=default_string_length)     :: message
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: work
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: a_lu, b
    INTEGER                                  :: info, iter, n, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipiv
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a_norm, old_eval_error, r_cond
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: berr, ferr, rwork
    REAL(KIND=dp), EXTERNAL                  :: zlange

    EXTERNAL zgecon,zgerfs,zgetrf,zgetrs

    failure = .FALSE.

    ! Check for optional parameter
    IF (PRESENT(option)) THEN
      trans = option
    ELSE
      trans = "N"
    END IF

    ! Get the dimension of matrix a
    n = SIZE(a,1)

    ! Check array dimensions
    IF (n == 0) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Matrix to be inverted of zero size")
    END IF

    IF (n /= SIZE(a,2)) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Check the array bounds of parameter #1")
    END IF

    IF ((n /= SIZE(a_inverse,1)).OR.&
        (n /= SIZE(a_inverse,2))) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Check the array bounds of parameter #2")
    END IF

    ! Allocate work storage
    ALLOCATE (a_lu(n,n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (b(n,n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (berr(n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (ferr(n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (ipiv(n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (rwork(2*n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (work(2*n),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    a_lu(1:n,1:n) = a(1:n,1:n)

    ! Compute the LU factorization of the matrix a
    CALL zgetrf(n,n,a_lu,n,ipiv,info)

    IF (info /= 0) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "The LU factorization in dgetrf failed")
    END IF

    ! Compute the norm of the matrix a

    IF (trans == "N") THEN
      norm = '1'
    ELSE
      norm = 'I'
    END IF

    a_norm = zlange(norm,n,n,a,n,work)

    ! Compute the reciprocal of the condition number of a

    CALL zgecon(norm,n,a_lu,n,a_norm,r_cond,work,rwork,info)

    IF (info /= 0) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "The computation of the condition number in "//&
                        "dgecon failed")
    END IF

    IF (r_cond < EPSILON(0.0_dp)) THEN
      WRITE (message,"(A,ES10.3)") "R_COND =",r_cond
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Bad condition number "//TRIM(message)//" (smaller than the machine "//&
                        "working precision)")
    END IF

    ! Solve a system of linear equations using the LU factorization computed by dgetrf

    CALL unit_matrix(a_inverse)

    CALL zgetrs(trans,n,n,a_lu,n,ipiv,a_inverse,n,info)

    IF (info /= 0) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Solving the system of linear equations in dgetrs "//&
                        "failed")
    END IF

    ! Improve the computed solution iteratively
    CALL unit_matrix(b) ! Initialize right-hand sides

    eval_error = 0.0_dp

    DO iter=1,10

      CALL zgerfs(trans,n,n,a,n,a_lu,n,ipiv,b,n,a_inverse,n,ferr,berr,&
                  work,rwork,info)

      IF (info /= 0) THEN
        CALL stop_program(moduleN,routineN,__LINE__,&
                          "Improving the computed solution in dgerfs failed")
      END IF

      old_eval_error = eval_error
      eval_error = MAXVAL(ferr)

      IF (ABS(eval_error - old_eval_error) <= EPSILON(1.0_dp)) EXIT

    END DO

    ! Release work storage
    DEALLOCATE (work,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (rwork,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (ipiv,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (ferr,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (berr,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (b,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (a_lu,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE invert_matrix_z

! *****************************************************************************
!> \brief  Reflection of the vector a through a mirror plane defined by the
!>         normal vector b. The reflected vector a is stored in a_mirror.
!> \param a ...
!> \param b ...
!> \retval a_mirror ...
!> \date    16.10.1998
!> \author  MK
!> \version 1.0
! *****************************************************************************
  FUNCTION reflect_vector(a,b) RESULT(a_mirror)

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a, b
    REAL(KIND=dp), DIMENSION(3)              :: a_mirror

    REAL(KIND=dp)                            :: length_of_b, scapro
    REAL(KIND=dp), DIMENSION(3)              :: d

    length_of_b = SQRT(b(1)*b(1) + b(2)*b(2) + b(3)*b(3))

    IF (length_of_b > eps_geo) THEN

       d(:) = b(:)/length_of_b

       ! Calculate the mirror image a_mirror of the vector a
       scapro = a(1)*d(1) + a(2)*d(2) + a(3)*d(3)

       a_mirror(:) = a(:) - 2.0_dp*scapro*d(:)

    ELSE

       a_mirror(:) = 0.0_dp

    END IF

  END FUNCTION reflect_vector

! *****************************************************************************
!> \brief   Rotation of the vector a about an rotation axis defined by the
!>          vector b. The rotation angle is phi (radians). The rotated vector
!>          a is stored in a_rot.
!> \param a ...
!> \param phi ...
!> \param b ...
!> \retval a_rot ...
!> \date    16.10.1998
!> \author  MK
!> \version 1.0
! *****************************************************************************
  FUNCTION rotate_vector(a,phi,b) RESULT(a_rot)
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a
    REAL(KIND=dp), INTENT(IN)                :: phi
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: b
    REAL(KIND=dp), DIMENSION(3)              :: a_rot

    REAL(KIND=dp)                            :: length_of_b
    REAL(KIND=dp), DIMENSION(3, 3)           :: rotmat

    length_of_b = SQRT(b(1)*b(1) + b(2)*b(2) + b(3)*b(3))
    IF (length_of_b > eps_geo) THEN

       ! Build up the rotation matrix rotmat
       CALL build_rotmat(phi,b,rotmat)

       ! Rotate the vector a by phi about the axis defined by vector b
       a_rot(:) = MATMUL(rotmat,a)

    ELSE

       a_rot(:) = 0.0_dp

    END IF

  END FUNCTION rotate_vector

! *****************************************************************************
!> \brief   Set the diagonal elements of matrix a to b.
!> \param a ...
!> \param b ...
!> \date    20.11.1998
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_diag_scalar_d(a,b)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    REAL(KIND=dp), INTENT(IN)                :: b

    INTEGER                                  :: i, n

    n = MIN(SIZE(a,1),SIZE(a,2))
    DO i=1,n
      a(i,i) = b
    END DO

  END SUBROUTINE set_diag_scalar_d

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
! *****************************************************************************
  SUBROUTINE set_diag_scalar_z(a,b)
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    COMPLEX(KIND=dp), INTENT(IN)             :: b

    INTEGER                                  :: i, n

    n = MIN(SIZE(a,1),SIZE(a,2))
    DO i=1,n
      a(i,i) = b
    END DO

  END SUBROUTINE set_diag_scalar_z

! *****************************************************************************
!> \brief   Symmetrize the matrix a.
!> \param a ...
!> \param option ...
!> \date    16.10.1998
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE symmetrize_matrix(a,option)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    CHARACTER(LEN=*), INTENT(IN)             :: option

    CHARACTER(len=*), PARAMETER :: routineN = 'symmetrize_matrix', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: message
    INTEGER                                  :: i, n

    n = MIN(SIZE(a,1),SIZE(a,2))

    IF (option == "lower_to_upper") THEN
       DO i=1,n-1
          a(i,i+1:n) = a(i+1:n,i)
       END DO
    ELSE IF (option == "upper_to_lower") THEN
       DO i=1,n-1
          a(i+1:n,i) = a(i,i+1:n)
       END DO
    ELSE IF (option == "anti_lower_to_upper") THEN
       DO i=1,n-1
          a(i,i+1:n) = -a(i+1:n,i)
       END DO
    ELSE IF (option == "anti_upper_to_lower") THEN
       DO i=1,n-1
          a(i+1:n,i) = -a(i,i+1:n)
       END DO
    ELSE
       WRITE (message,"(A)")&
            "Invalid option <"//TRIM(option)//"> was specified for parameter #2"
       CALL stop_program(moduleN,routineN,__LINE__,message)
    END IF

  END SUBROUTINE symmetrize_matrix

! *****************************************************************************
!> \brief   Set the matrix a to be a unit matrix.
!> \param a ...
!> \date    16.10.1998
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE unit_matrix_d(a)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a

    a(:,:) = 0.0_dp
    CALL set_diag(a,1.0_dp)

  END SUBROUTINE unit_matrix_d

! *****************************************************************************
!> \brief ...
!> \param a ...
! *****************************************************************************
  SUBROUTINE unit_matrix_z(a)
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a

    a(:,:) = (0.0_dp,0.0_dp)
    CALL set_diag(a,(1.0_dp,0.0_dp))

  END SUBROUTINE unit_matrix_z

! *****************************************************************************
!> \brief   Calculation of the vector product c = a x b.
!> \param a ...
!> \param b ...
!> \retval c ...
!> \date    16.10.1998
!> \author  MK
!> \version 1.0
! *****************************************************************************
  PURE FUNCTION vector_product(a,b) RESULT(c)
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a, b
    REAL(KIND=dp), DIMENSION(3)              :: c

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  END FUNCTION vector_product

! *****************************************************************************
!> \brief computes the greatest common divisor of two number
!> \param a ...
!> \param b ...
!> \author Joost VandeVondele
! *****************************************************************************
  ELEMENTAL FUNCTION gcd(a,b)
    INTEGER, INTENT(IN)                      :: a, b
    INTEGER                                  :: gcd

    INTEGER                                  :: aa, ab, l, rem, s

    aa=ABS(a)
    ab=ABS(b)
    IF (aa<ab) THEN
       s=aa
       l=ab
    ELSE
       s=ab
       l=aa
    ENDIF
    IF (s.NE.0) THEN
       DO
          rem=MOD(l,s)
          IF (rem==0) EXIT
          l=s
          s=rem
       ENDDO
       GCD=s
    ELSE
       GCD=l
    ENDIF
  END FUNCTION gcd

! *****************************************************************************
!> \brief computes the least common multiplier of two numbers
!> \param a ...
!> \param b ...
!> \author Joost VandeVondele
! *****************************************************************************
  ELEMENTAL FUNCTION lcm(a,b)
    INTEGER, INTENT(IN)                      :: a, b
    INTEGER                                  :: lcm

    INTEGER                                  :: tmp

    tmp = gcd(a,b)
    IF (tmp==0) THEN
       lcm = 0
    ELSE
       ! could still overflow if the true lcm is larger than maxint
       lcm = ABS((a/tmp)*b)
    END IF
  END FUNCTION lcm

! *****************************************************************************
!> \brief computes the exponential integral
!>      Ei(x) = Int(exp(-x*t)/t,t=1..infinity)  x>0
!> \param x ...
!> \author JGH (adapted from Numerical recipies)
! *****************************************************************************
  FUNCTION ei(x)
    REAL(dp)                                 :: x, ei

    CHARACTER(len=*), PARAMETER :: routineN = 'ei', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: maxit = 100
    REAL(dp), PARAMETER                      :: eps = EPSILON(0.0_dp), &
                                                fpmin = TINY(0.0_dp)

    INTEGER                                  :: k
    REAL(dp)                                 :: fact, prev, sum1, term

    IF(x <= 0._dp) THEN
       CALL stop_program(moduleN,routineN,__LINE__,"Invalid argument")
    END IF

    IF (x < fpmin) THEN
       ei = LOG(x) + euler
    ELSE IF (x <= -LOG(EPS)) THEN
       sum1 = 0._dp
       fact = 1._dp
       DO k = 1,maxit
          fact = fact*x/REAL(k,dp)
          term = fact/REAL(k,dp)
          sum1 = sum1+term
          IF (term < eps*sum1) EXIT
       END DO
       ei = sum1+LOG(x)+euler
    ELSE
       sum1 = 0._dp
       term = 1._dp
       DO k = 1,maxit
          prev = term
          term = term*REAL(k,dp)/x
          IF (term < eps) EXIT
          IF (term < prev) THEN
             sum1 = sum1+term
          ELSE
             sum1 = sum1-prev
             EXIT
          END IF
       END DO
       ei = EXP(x)*(1._dp+sum1)/x
    END IF

  END FUNCTION ei

! *****************************************************************************
!> \brief ...
!> \param mat1 ...
!> \param mat2 ...
! *****************************************************************************
  FUNCTION matmul_3x3 ( mat1, mat2 )
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: mat1, mat2
    REAL(KIND=dp), DIMENSION(3, 3)           :: matmul_3x3

    matmul_3x3(1,1)=mat1(1,1)*mat2(1,1)+mat1(1,2)*mat2(2,1)+mat1(1,3)*mat2(3,1)
    matmul_3x3(1,2)=mat1(1,1)*mat2(1,2)+mat1(1,2)*mat2(2,2)+mat1(1,3)*mat2(3,2)
    matmul_3x3(1,3)=mat1(1,1)*mat2(1,3)+mat1(1,2)*mat2(2,3)+mat1(1,3)*mat2(3,3)
    matmul_3x3(2,1)=mat1(2,1)*mat2(1,1)+mat1(2,2)*mat2(2,1)+mat1(2,3)*mat2(3,1)
    matmul_3x3(2,2)=mat1(2,1)*mat2(1,2)+mat1(2,2)*mat2(2,2)+mat1(2,3)*mat2(3,2)
    matmul_3x3(2,3)=mat1(2,1)*mat2(1,3)+mat1(2,2)*mat2(2,3)+mat1(2,3)*mat2(3,3)
    matmul_3x3(3,1)=mat1(3,1)*mat2(1,1)+mat1(3,2)*mat2(2,1)+mat1(3,3)*mat2(3,1)
    matmul_3x3(3,2)=mat1(3,1)*mat2(1,2)+mat1(3,2)*mat2(2,2)+mat1(3,3)*mat2(3,2)
    matmul_3x3(3,3)=mat1(3,1)*mat2(1,3)+mat1(3,2)*mat2(2,3)+mat1(3,3)*mat2(3,3)
  END FUNCTION matmul_3x3

! *****************************************************************************
!> \brief ...
!> \param res ...
!> \param mat ...
!> \param vec ...
! *****************************************************************************
  SUBROUTINE matvec_3x3 (res, mat, vec )
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: res
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: mat
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: vec

    res(1)=mat(1,1)*vec(1)+mat(1,2)*vec(2)+mat(1,3)*vec(3)
    res(2)=mat(2,1)*vec(1)+mat(2,2)*vec(2)+mat(2,3)*vec(3)
    res(3)=mat(3,1)*vec(1)+mat(3,2)*vec(2)+mat(3,3)*vec(3)
  END SUBROUTINE matvec_3x3

! *****************************************************************************
!> \brief ...
!> \param vec1 ...
!> \param vec2 ...
! *****************************************************************************
  FUNCTION dotprod_3d ( vec1, vec2 )
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: vec1, vec2
    REAL(KIND=dp)                            :: dotprod_3d

    dotprod_3d = &
         vec1 ( 1 ) * vec2 ( 1 ) &
         + vec1 ( 2 ) * vec2 ( 2 ) &
         + vec1 ( 3 ) * vec2 ( 3 )
  END FUNCTION dotprod_3d

! *****************************************************************************
!> \brief ...
!> \param mat ...
! *****************************************************************************
  FUNCTION transpose_3d ( mat )
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: mat
    REAL(KIND=dp), DIMENSION(3, 3)           :: transpose_3d

    INTEGER                                  :: i, j

    DO i = 1, 3
       DO j = 1, 3
          transpose_3d ( j, i ) = mat ( i, j )
       END DO
    END DO
  END FUNCTION transpose_3d

! *****************************************************************************
!> \brief computes the exponential integral
!>      En(x) = Int(exp(-x*t)/t^n,t=1..infinity)  x>0, n=0,1,..
!>      Note: Ei(-x) = -E1(x)
!> \param n ...
!> \param x ...
!> \par History
!>      05.2007 Created
!> \author Manuel Guidon (adapted from Numerical recipies)
! *****************************************************************************
  FUNCTION expint(n,x)
    INTEGER                                  :: n
    REAL(dp)                                 :: x, expint

    CHARACTER(len=*), PARAMETER :: routineN = 'expint', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: maxit = 100
    REAL(dp), PARAMETER :: eps = 6.e-14_dp, &
      euler = 0.5772156649015328606065120_dp, fpmin = TINY(0.0_dp)

    INTEGER                                  :: i, ii, nm1
    REAL(dp)                                 :: a, b, c, d, del, fact, h, psi

    nm1=n-1

    IF(n.lt.0.OR.x.lt.0.0_dp.OR.(x.eq.0.0_dp.AND.(n.EQ.0.or.n.EQ.1))) THEN
      CALL stop_program(moduleN,routineN,__LINE__,"Invalid argument")
    ELSE IF(n.EQ.0) THEN       !Special case.
      expint=EXP(-x)/x
    ELSE IF(x.EQ.0.0_dp) THEN  !Another special case.
      expint=1.0_dp/nm1
    ELSE IF(x.GT.1.0_dp) THEN  !Lentz's algorithm (§5.2).
      b=x+n
      c=1.0_dp/FPMIN
      d=1.0_dp/b
      h=d
      DO i = 1,MAXIT
        a=-i*(nm1+i)
        b=b+2.0_dp
        d=1.0_dp/(a*d+b)
        c=b+a/c
        del=c*d
        h=h*del
        IF(ABS(del-1.0_dp).LT.EPS) THEN
          expint=h*EXP(-x)
          RETURN
        END IF
      END DO
      CALL stop_program(moduleN,routineN,__LINE__,"continued fraction failed in expint")
    ELSE !Evaluate series.
      IF(nm1.NE.0)THEN  !Set first term.
        expint=1.0_dp/nm1
      ELSE
        expint=-LOG(x)-euler
      END IF
      fact=1.0_dp
      DO i=1,MAXIT
        fact=-fact*x/i
        IF(i.NE.nm1) THEN
          del=-fact/(i-nm1)
        ELSE
          psi=-euler !Compute I(n).
          DO ii=1,nm1
            psi=psi+1.0_dp/ii
          END DO
          del=fact*(-LOG(x)+psi)
        END IF
        expint=expint+del
        IF(ABS(del).LT.ABS(expint)*EPS) RETURN
      END DO
       CALL stop_program(moduleN,routineN,__LINE__,"series failed in expint")
    END IF
    RETURN
  END FUNCTION expint

! *****************************************************************************
!> \brief  Jacobi matrix diagonalization. The eigenvalues are returned in
!>         vector d and the eigenvectors are returned in matrix v in ascending
!>         order.
!>
!> \param a ...
!> \param d ...
!> \param v ...
!> \par History
!>         - Creation (20.11.98, Matthias Krack)
! *****************************************************************************
  SUBROUTINE jacobi(a,d,v)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: d
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: v

    INTEGER                                  :: n

    n = SIZE(d(:))

    ! Diagonalize matrix a
    CALL diag(n,a,d,v)

    ! Sort eigenvalues and eigenvector in ascending order
    CALL eigsrt(n,d,v)

  END SUBROUTINE jacobi

! *****************************************************************************
!> \brief  Diagonalize matrix a. The eigenvalues are returned in vector d
!>         and the eigenvectors are returned in matrix v.
!>
!> \param n ...
!> \param a ...
!> \param d ...
!> \param v ...
!> \par History
!>         - Creation (20.11.98, Matthias Krack)
! *****************************************************************************
  SUBROUTINE diag(n,a,d,v)

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(n, n), &
      INTENT(INOUT)                          :: a
    REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: d
    REAL(KIND=dp), DIMENSION(n, n), &
      INTENT(OUT)                            :: v

    INTEGER                                  :: i, ip, iq
    REAL(KIND=dp)                            :: a_max, c, d_min, g, h, s, t, &
                                                tau, theta, tresh
    REAL(KIND=dp), DIMENSION(n)              :: b, z

    CALL unit_matrix(v(:,:))

    b(:) = get_diag(a(:,:))
    d(:) = b(:)
    z(:) = 0.0_dp

    ! Go for 50 iterations
    DO i=1,50
       d_min = MAX(1.0E-3_dp,MINVAL(ABS(d(:))))
       a_max = 0.0_dp
       DO ip=1,n-1
          a_max = MAX(a_max,MAXVAL(ABS(a(ip,ip+1:n))))
       END DO
       IF (a_max < 1.0E-10_dp*d_min) RETURN
       tresh = MERGE(a_max,0.0_dp,(i < 4))
       DO ip=1,n-1
          DO iq=ip+1,n
             g = 100.0_dp*ABS(a(ip,iq))
             IF ((i > 4).AND.&
                  ((ABS(d(ip)) + g) == ABS(d(ip))).AND.&
                  ((ABS(d(iq)) + g) == ABS(d(iq)))) THEN
                a(ip,iq) = 0.0_dp
             ELSE IF (ABS(a(ip,iq)) > tresh) THEN
                h = d(iq) - d(ip)
                IF ((ABS(h) + g) == ABS(h)) THEN
                   t = a(ip,iq)/h
                ELSE
                   theta = 0.5_dp*h/a(ip,iq)
                   t = 1.0_dp/(ABS(theta) + SQRT(1.0_dp + theta**2))
                   IF (theta < 0.0_dp) t = -t
                END IF
                c = 1.0_dp/SQRT(1.0_dp + t**2)
                s = t*c
                tau = s/(1.0_dp + c)
                h = t*a(ip,iq)
                z(ip) = z(ip) - h
                z(iq) = z(iq) + h
                d(ip) = d(ip) - h
                d(iq) = d(iq) + h
                a(ip,iq) = 0.0_dp
                CALL jrotate(a(1:ip-1,ip),a(1:ip-1,iq),s,tau)
                CALL jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq),s,tau)
                CALL jrotate(a(ip,iq+1:n),a(iq,iq+1:n),s,tau)
                CALL jrotate(v(:,ip),v(:,iq),s,tau)
             END IF
          END DO
       END DO
       b(:) = b(:) + z(:)
       d(:) = b(:)
       z(:) = 0.0_dp
    END DO
    WRITE (*,'(/,T2,A,/)') 'Too many iterations in jacobi'

  END SUBROUTINE diag

! *****************************************************************************
!> \brief  Perform a Jacobi rotation of the vectors a and b.
!>
!> \param a ...
!> \param b ...
!> \param ss ...
!> \param tt ...
!> \par History
!>         - Creation (20.11.98, Matthias Krack)
! *****************************************************************************
  SUBROUTINE jrotate(a,b,ss,tt)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: a, b
    REAL(KIND=dp), INTENT(IN)                :: ss, tt

    REAL(KIND=dp), DIMENSION(SIZE(a))        :: c

    c(:) = a(:)
    a(:) = a(:) - ss*(b(:) + a(:)*tt)
    b(:) = b(:) + ss*(c(:) - b(:)*tt)

  END SUBROUTINE jrotate

! *****************************************************************************
!> \brief Sort the values in vector d in ascending order and swap the
!>        corresponding columns of matrix v.
!>
!> \param n ...
!> \param d ...
!> \param v ...
!> \par History
!>         - Creation (20.11.98, Matthias Krack)
! *****************************************************************************
  SUBROUTINE eigsrt(n,d,v)
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: d
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: v

    INTEGER                                  :: i, j

    DO i=1,n-1
       j = SUM(MINLOC(d(i:n))) + i - 1
       IF (j /= i) THEN
          CALL swap(d(i),d(j))
          CALL swap(v(:,i),v(:,j))
       END IF
    END DO

  END SUBROUTINE eigsrt

! **************************************************************************
!> \brief Swap two scalars
!>
!> \param a ...
!> \param b ...
!> \par History
!>         - Creation (20.11.98, Matthias Krack)
! *****************************************************************************
   SUBROUTINE swap_scalar(a,b)
    REAL(KIND=dp), INTENT(INOUT)             :: a, b

    REAL(KIND=dp)                            :: c

    c = a
    a = b
    b = c

   END SUBROUTINE swap_scalar

! **************************************************************************
!> \brief Swap two vectors
!>
!> \param a ...
!> \param b ...
!> \par History
!>         - Creation (20.11.98, Matthias Krack)
! *****************************************************************************
   SUBROUTINE swap_vector(a,b)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: a, b

    CHARACTER(LEN=*), PARAMETER :: routineN = 'swap_vector', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n
    REAL(KIND=dp)                            :: c

    n = SIZE(a)

    IF (n /= SIZE(b)) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Check the array bounds of the parameters")
    END IF

    DO i=1,n
       c = a(i)
       a(i) = b(i)
       b(i) = c
    END DO

  END SUBROUTINE swap_vector

END MODULE mathlib
