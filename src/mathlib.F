!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/mathlib [1.0] *
!!
!!   NAME
!!     mathlib
!!
!!   FUNCTION
!!     Collection of simple mathematical functions and subroutines
!!
!!   AUTHOR
!!     MK (15.11.1998)
!!
!!   MODIFICATION HISTORY
!!     FUNCTION angle updated and FUNCTION dihedral angle added; cleaned
!!     (13.03.2004,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE mathlib

  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_size
  USE mathconstants,                   ONLY: euler,&
                                             fac
  USE output_utilities,                ONLY: print_message
  USE termination,                     ONLY: get_error_unit,&
                                             stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mathlib'
  REAL(KIND = dp), PARAMETER  :: eps_geo = 1.0E-6_dp

! *** Public subroutines ***

  PUBLIC :: build_rotmat,&
            diamat_all,&
            invert_matrix,&
            power_matrix,&
            set_diag,&
            symmetrize_matrix,&
            unit_matrix

! *** Public functions ***

  PUBLIC :: angle,&
            binomial,&
            det_3x3,&
            dihedral_angle,&
            gcd,&
            inv_3x3,&
            lcm,&
            ei,&
            vector_product

! *****************************************************************************

  INTERFACE det_3x3
    MODULE PROCEDURE det_3x3_1,det_3x3_2
  END INTERFACE

  INTERFACE set_diag
    MODULE PROCEDURE set_diag_scalar,set_diag_vector
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  FUNCTION angle(a,b) RESULT(angle_ab)

!   Purpose: Calculation of the angle between the vectors a and b.
!            The angle is returned in radians.

!   History: - Creation (14.10.1998,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a, b
    REAL(KIND=dp)                            :: angle_ab

    REAL(KIND=dp)                            :: length_of_a, length_of_b
    REAL(KIND=dp), DIMENSION(SIZE(a, 1))     :: a_norm, b_norm

!   ---------------------------------------------------------------------------

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

  FUNCTION binomial(n,k) RESULT(n_over_k)

!   Purpose: The binomial coefficient n over k for 0 <= k <= n is calculated,
!            otherwise zero is returned.

!   History: - Creation (08.03.1999,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: n, k
    REAL(KIND=dp)                            :: n_over_k

!   ---------------------------------------------------------------------------

    IF ((k >= 0).AND.(k <= n)) THEN
      n_over_k = fac(n)/(fac(n-k)*fac(k))
    ELSE
      n_over_k = 0.0_dp
    END IF

  END FUNCTION binomial

! *****************************************************************************

  SUBROUTINE build_rotmat(phi,a,rotmat)

!   Purpose: The rotation matrix rotmat which rotates a vector about a
!            rotation axis defined by the vector a is build up.
!            The rotation angle is phi (radians).

!   History: - Creation (16.10.1998,MK)

!   ***************************************************************************

    REAL(KIND=dp), INTENT(IN)                :: phi
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT)                            :: rotmat

    REAL(KIND=dp)                            :: cosp, cost, length_of_a, sinp
    REAL(KIND=dp), DIMENSION(3)              :: d

!   ---------------------------------------------------------------------------
!   *** Check the length of the vector a ***

    length_of_a = SQRT(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))

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

  FUNCTION det_3x3_1(a) RESULT(det_a)

!   Purpose: Returns the determinante of the 3x3 matrix a.

!   History: - Creation (13.03.2004,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: a
    REAL(KIND=dp)                            :: det_a

!   ---------------------------------------------------------------------------

    det_a = a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) +&
            a(1,2)*(a(2,3)*a(3,1) - a(2,1)*a(3,3)) +&
            a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))

  END FUNCTION det_3x3_1

! *****************************************************************************

  FUNCTION det_3x3_2(a1,a2,a3) RESULT(det_a)

!   Purpose: Returns the determinante of the 3x3 matrix a given by its columns.

!   History: - Creation (13.03.2004,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a1, a2, a3
    REAL(KIND=dp)                            :: det_a

!   ---------------------------------------------------------------------------

    det_a = a1(1)*(a2(2)*a3(3) - a3(2)*a2(3)) +&
            a2(1)*(a3(2)*a1(3) - a1(2)*a3(3)) +&
            a3(1)*(a1(2)*a2(3) - a2(2)*a1(3))

  END FUNCTION det_3x3_2

! *****************************************************************************

  SUBROUTINE diamat_all(a,eigval,dac)

!   Purpose: Diagonalize the symmetric n by n matrix a using the LAPACK
!            library. Only the upper triangle of matrix a is used.

!   History: - Creation (29.03.1999,MK)

!   ***************************************************************************

!   a       : Symmetric matrix to be diagonalized (input; upper triangle) ->
!             eigenvectors of the matrix a (output).
!   dac     : If true, then the divide-and-conquer algorithm is applied.
!   eigval  : Eigenvalues of the matrix a (output).

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigval
    LOGICAL, INTENT(IN), OPTIONAL            :: dac

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'diamat_all'

    INTEGER                                  :: info, istat, liwork, lwork, &
                                                n, nb
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iwork
    INTEGER, EXTERNAL                        :: ilaenv
    LOGICAL                                  :: divide_and_conquer
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: work

!   *** Externals (LAPACK 3.0) ***

    EXTERNAL dsyev,dsyevd

!   ---------------------------------------------------------------------------

!   *** Get the size of the matrix a ***

    n = SIZE(a,1)

!   *** Check the size of matrix a ***

    IF (SIZE(a,2) /= n) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Check the size of matrix a (parameter #1)")
    END IF

!   *** Check the size of vector eigval ***

    IF (SIZE(eigval) /= n) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "The dimension of vector eigval is too small")
    END IF

!   *** Check, if the divide-and-conquer algorithm is requested ***

    IF (PRESENT(dac)) THEN
      divide_and_conquer = dac
    ELSE
      divide_and_conquer = .FALSE.
    END IF

!   *** Get the optimal work storage size ***

    IF (divide_and_conquer) THEN
      lwork = 2*n**2 + 6*n + 1
      liwork = 5*n + 3
    ELSE
      nb = ilaenv(1,"DSYTRD","U",n,-1,-1,-1)
      lwork = (nb + 2)*n
    END IF

!   *** Allocate work storage ***

    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "work",lwork*dp_size)
    IF (divide_and_conquer) THEN
      ALLOCATE (iwork(liwork),STAT=istat)
      IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                       "iwork",liwork*int_size)
    END IF

!   *** Diagonalize the matrix a ***

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

!   *** Release work storage ***

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "work")
    IF (divide_and_conquer) THEN
      DEALLOCATE (iwork,STAT=istat)
      IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                       "iwork")
    END IF

  END SUBROUTINE diamat_all

! *****************************************************************************

  FUNCTION dihedral_angle(ab,bc,cd) RESULT(dihedral_angle_abcd)

!   Purpose: Returns the dihedral angle, i.e. the angle between the planes
!            defined by the vectors (-ab,bc) and (cd,-bc).
!            The dihedral angle is returned in radians.

!   History: - Creation (13.03.2004,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ab, bc, cd
    REAL(KIND=dp)                            :: dihedral_angle_abcd

    REAL(KIND=dp)                            :: det_abcd
    REAL(KIND=dp), DIMENSION(3)              :: abc, bcd

!   ---------------------------------------------------------------------------
!   *** Calculate the normal vectors of the planes ***
!   *** defined by the points a,b,c and b,c,d      ***

    abc = vector_product(bc,-ab)
    bcd = vector_product(cd,-bc)

    det_abcd = det_3x3(abc,bcd,-bc)

    dihedral_angle_abcd = SIGN(1.0_dp,det_abcd)*angle(abc,bcd)

  END FUNCTION dihedral_angle

! *****************************************************************************

  FUNCTION gauss_exponent(l,radius,threshold,prefactor) RESULT(exponent)

!   Purpose: The exponent of a primitive Gaussian function for a given radius
!            and threshold is calculated.

!   History: - Creation (07.03.1999,MK)

!   ***************************************************************************

!   exponent : Exponent of the primitive Gaussian function.
!   l        : Angular momentum quantum number l.
!   prefactor: Prefactor of the Gaussian function (e.g. a contraction
!              coefficient).
!   radius   : Calculated radius of the Gaussian function.
!   threshold: Threshold for radius.

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), INTENT(IN)                :: radius, threshold, prefactor
    REAL(KIND=dp)                            :: exponent

!   ---------------------------------------------------------------------------

    exponent = 0.0_dp

    IF (radius < 1.0E-6_dp) RETURN
    IF (threshold < 1.0E-12_dp) RETURN

    exponent = LOG(ABS(prefactor)*radius**l/threshold)/radius**2

  END FUNCTION gauss_exponent

! *****************************************************************************

  FUNCTION get_diag(a) RESULT(a_diag)

!   Purpose: Return the diagonal elements of matrix a as a vector.

!   History: - Creation (20.11.1998,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: a
    REAL(KIND=dp), &
      DIMENSION(MIN(SIZE(a, 1), SIZE(a, 2))) :: a_diag

    INTEGER                                  :: i, n

!   ---------------------------------------------------------------------------

    n = MIN(SIZE(a,1),SIZE(a,2))

    DO i=1,n
      a_diag(i) = a(i,i)
    END DO

  END FUNCTION get_diag

! *****************************************************************************

  FUNCTION inv_3x3(a) RESULT(a_inv)

!   Purpose: Returns the inverse of the 3 x 3 matrix a.

!   History: - Creation (13.03.2004,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: a
    REAL(KIND=dp), DIMENSION(3, 3)           :: a_inv

    REAL(KIND=dp)                            :: det_a

!   ---------------------------------------------------------------------------

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

  SUBROUTINE invert_matrix(a,a_inverse,error,option)

!   Purpose: Compute the inverse of the n by n matrix a using the LAPACK
!            library.

!   History: - Creation (23.03.1999,MK)

!   ***************************************************************************

!   a        : Matrix to be inverted (input).
!   a_inverse: Inverse of the matrix a (output).
!   a_lu     : LU factorization of matrix a.
!   a_norm   : Norm of matrix a.
!   error    : Estimated error of the inversion.
!   r_cond   : Reciprocal condition number of the matrix a.
!   trans    : "N" => invert a
!              "T" => invert transpose(a)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: a
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: a_inverse
    REAL(KIND=dp), INTENT(OUT)               :: error
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL   :: option

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'invert_matrix'

    CHARACTER(LEN=1)                         :: norm, trans
    CHARACTER(LEN=80)                        :: message
    INTEGER                                  :: info, istat, iter, n, &
                                                output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipiv, iwork
    REAL(KIND=dp)                            :: a_norm, old_error, r_cond
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: berr, ferr, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: a_lu, b
    REAL(KIND=dp), EXTERNAL                  :: dlange

!   *** Externals (LAPACK 3.0) ***

    EXTERNAL dgecon,dgerfs,dgetrf,dgetrs

!   ---------------------------------------------------------------------------

!   *** Check for optional parameter ***

    IF (PRESENT(option)) THEN
      trans = option
    ELSE
      trans = "N"
    END IF

!   *** Get the dimension of matrix a ***

    n = SIZE(a,1)

!   *** Check array dimensions ***

    IF (n == 0) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Matrix to be inverted of zero size")
    ENDIF

    IF (n /= SIZE(a,2)) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Check the array bounds of parameter #1")
    END IF

    IF ((n /= SIZE(a_inverse,1)).OR.&
        (n /= SIZE(a_inverse,2))) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Check the array bounds of parameter #2")
    END IF

!   *** Allocate work storage ***

    ALLOCATE (a_lu(n,n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "a_lu",n*n*dp_size)

    ALLOCATE (b(n,n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "b",n*n*dp_size)

    ALLOCATE (berr(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "berr",n*dp_size)

    ALLOCATE (ferr(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "ferr",n*dp_size)

    ALLOCATE (ipiv(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "ipiv",n*int_size)

    ALLOCATE (iwork(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "iwork",n*int_size)

    ALLOCATE (work(4*n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "work",4*n*int_size)

    a_lu(1:n,1:n) = a(1:n,1:n)

!   *** Compute the LU factorization of the matrix a ***
    CALL dgetrf(n,n,a_lu,n,ipiv,info)

    IF (info /= 0) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "The LU factorization in dgetrf failed")
    END IF

!   *** Compute the norm of the matrix a ***

    IF (trans == "N") THEN
      norm = '1'
    ELSE
      norm = 'I'
    END IF

    a_norm = dlange(norm,n,n,a,n,work)

!   *** Compute the reciprocal of the condition number of a ***

    CALL dgecon(norm,n,a_lu,n,a_norm,r_cond,work,iwork,info)

    IF (info /= 0) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "The computation of the condition number in "//&
                        "dgecon failed")
    END IF

    IF (r_cond < EPSILON(0.0_dp)) THEN
      CALL get_error_unit(output_unit)
      WRITE (message,"(A,ES10.3,A)") "R_COND =",r_cond
      CALL print_message(message,output_unit,1,1,1)
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Bad condition number (smaller than the machine "//&
                        "working precision)")
    END IF

!   *** Solve a system of linear equations using ***
!   *** the LU factorization computed by dgetrf  ***

    CALL unit_matrix(a_inverse)

    CALL dgetrs(trans,n,n,a_lu,n,ipiv,a_inverse,n,info)

    IF (info /= 0) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Solving the system of linear equations in dgetrs "//&
                        "failed")
    END IF

!   *** Improve the computed solution iteratively ***

    CALL unit_matrix(b) ! Initialize right-hand sides

    error = 0.0_dp

    DO iter=1,10

      CALL dgerfs(trans,n,n,a,n,a_lu,n,ipiv,b,n,a_inverse,n,ferr,berr,&
                  work,iwork,info)

      IF (info /= 0) THEN
        CALL stop_program(moduleN,routineN,__LINE__,&
                          "Improving the computed solution in dgerfs failed")
      END IF

      old_error = error
      error = MAXVAL(ferr)

      IF (ABS(error - old_error) <= EPSILON(1.0_dp)) EXIT

    END DO

!   *** Release work storage ***

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "work")
    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "iwork")
    DEALLOCATE (ipiv,STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "ipiv")
    DEALLOCATE (ferr,STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "ferr")
    DEALLOCATE (berr,STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "berr")
    DEALLOCATE (b,STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "b")
    DEALLOCATE (a_lu,STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "a_lu")

  END SUBROUTINE invert_matrix

! *****************************************************************************

  SUBROUTINE power_matrix(a,a_power,exponent,threshold,n_dependent,dac)

!   Purpose: Raise the real symmetric n by n matrix a to the power given by
!            exponent. All eigenvectors with a corresponding eigenvalue lower
!            than threshold are quenched. Only the upper triangle of matrix a
!            is used.

!   History: - Creation (29.03.1999,MK)

!   ***************************************************************************

!   a          : Symmetric matrix to be powered (input; upper triangle) ->
!                Destroyed on exit.
!   a_power    : Power of matrix a => a**exponent (output).
!   dac        : Use the divide-and-conquer algorithm for the diagonalization.
!   exponent   : Matrix exponent (input).
!   n_dependent: Number of the eigenvectors which are linear dependent due to
!                the defined eigval_eps (output).
!   threshold  : Threshold value for eigenvector quenching (input).

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: a_power
    REAL(KIND=dp), INTENT(IN)                :: exponent
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: threshold
    INTEGER, INTENT(OUT), OPTIONAL           :: n_dependent
    LOGICAL, INTENT(IN), OPTIONAL            :: dac

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'power_matrix'

    INTEGER                                  :: i, istat, n, n_dep
    LOGICAL                                  :: divide_and_conquer
    REAL(KIND=dp)                            :: eps_eigval, expa
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigval

    EXTERNAL dsyrk

!   ---------------------------------------------------------------------------

!   *** Define the threshold for the eigenvalue quenching ***

    IF (PRESENT(threshold)) THEN
      eps_eigval = threshold
    ELSE
      eps_eigval = EPSILON(0.0_dp)
    END IF

!   *** Get the dimension of matrix a ***

    n = SIZE(a,1)

!   *** Check array dimensions ***

    IF (n /= SIZE(a,2)) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Check the array bounds of parameter #1")
    END IF

    IF ((n /= SIZE(a_power,1)).OR.&
        (n /= SIZE(a_power,2))) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Check the array bounds of parameter #2")
    END IF

!   *** Allocate work storage ***

    ALLOCATE (eigval(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "eigval",n*dp_size)

!   *** Check, if the divide-and-conquer algorithm is requested ***

    IF (PRESENT(dac)) THEN
      divide_and_conquer = dac
    ELSE
      divide_and_conquer = .FALSE.
    END IF

!   *** Compute the eigenvectors and eigenvalues of the matrix a ***

    CALL diamat_all(a,eigval,divide_and_conquer)

!   *** Build a**exponent with eigenvector quenching ***

    expa = 0.5_dp*exponent

    n_dep = 0

    DO i=1,n
      IF (eigval(i) < eps_eigval) THEN
        a(1:n,i) = 0.0_dp
        n_dep = n_dep + 1
      ELSE
        eigval(i) = eigval(i)**expa
        a(1:n,i) = eigval(i)*a(1:n,i)
      END IF
    END DO

    IF (PRESENT(n_dependent)) THEN
      n_dependent = n_dep
    END IF

!   *** a_power <- a*Transpose(a) ***

    CALL dsyrk("U","N",n,n,1.0_dp,a(1,1),n,0.0_dp,a_power(1,1),n)

!   *** Copy upper triangle of matrix a_power to lower triangle ***

    CALL symmetrize_matrix(a_power,"upper_to_lower")

!   *** Release work storage ***

    DEALLOCATE (eigval,STAT=istat)
    IF (istat /= 0) CALL stop_memory(moduleN,routineN,__LINE__,&
                                     "eigval")

  END SUBROUTINE power_matrix

! *****************************************************************************

  FUNCTION reflect_vector(a,b) RESULT(a_mirror)

!   Purpose: Reflection of the vector a through a mirror plane defined by the
!            normal vector b. The reflected vector a is stored in a_mirror.

!   History: - Creation (16.10.1998,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a, b
    REAL(KIND=dp), DIMENSION(3)              :: a_mirror

    REAL(KIND=dp)                            :: length_of_b, scapro
    REAL(KIND=dp), DIMENSION(3)              :: d

!   ---------------------------------------------------------------------------

    length_of_b = SQRT(b(1)*b(1) + b(2)*b(2) + b(3)*b(3))

    IF (length_of_b > eps_geo) THEN

      d(:) = b(:)/length_of_b

!     *** Calculate the mirror image a_mirror of the vector a ***

      scapro = a(1)*d(1) + a(2)*d(2) + a(3)*d(3)

      a_mirror(:) = a(:) - 2.0_dp*scapro*d(:)

    ELSE

      a_mirror(:) = 0.0_dp

    END IF

  END FUNCTION reflect_vector

! *****************************************************************************

  FUNCTION rotate_vector(a,phi,b) RESULT(a_rot)

!   Purpose: Rotation of the vector a about an rotation axis defined by the
!            vector b. The rotation angle is phi (radians). The rotated vector
!            a is stored in a_rot.

!   History: - Creation (16.10.1998,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a
    REAL(KIND=dp), INTENT(IN)                :: phi
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: b
    REAL(KIND=dp), DIMENSION(3)              :: a_rot

    REAL(KIND=dp)                            :: length_of_b
    REAL(KIND=dp), DIMENSION(3, 3)           :: rotmat

!   ---------------------------------------------------------------------------

    length_of_b = SQRT(b(1)*b(1) + b(2)*b(2) + b(3)*b(3))

    IF (length_of_b > eps_geo) THEN

!     *** Build up the rotation matrix rotmat ***

      CALL build_rotmat(phi,b,rotmat)

!     *** Rotate the vector a by phi about the axis defined by vector b ***

      a_rot(:) = MATMUL(rotmat,a)

    ELSE

      a_rot(:) = 0.0_dp

    END IF

  END FUNCTION rotate_vector

! *****************************************************************************

  SUBROUTINE set_diag_scalar(a,b)

!   Purpose: Set the diagonal elements of matrix a to b.

!   History: - Creation (20.11.1998,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    REAL(KIND=dp), INTENT(IN)                :: b

    INTEGER                                  :: i, n

!   ---------------------------------------------------------------------------

    n = MIN(SIZE(a,1),SIZE(a,2))

    DO i=1,n
      a(i,i) = b
    END DO

  END SUBROUTINE set_diag_scalar

! *****************************************************************************

  SUBROUTINE set_diag_vector(a,b)

!   Purpose: Set the diagonal elements of matrix a to the values of vector b.

!   History: - Creation (20.11.1998,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: b

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'set_diag_vector'

    INTEGER                                  :: i, n

!   ---------------------------------------------------------------------------

    n = MIN(SIZE(a,1),SIZE(a,2))

    IF (SIZE(b) == n) THEN
      DO i=1,n
        a(i,i) = b(i)
      END DO
    ELSE
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "Check the array bounds of the parameters")
    END IF

  END SUBROUTINE set_diag_vector

! *****************************************************************************

  SUBROUTINE symmetrize_matrix(a,option)

!   Purpose: Symmetrize the matrix a.

!   History: - Creation (16.10.1998,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    CHARACTER(LEN=*), INTENT(IN)             :: option

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'symmetrize_matrix'

    CHARACTER(LEN=default_string_length)     :: message
    INTEGER                                  :: i, n

!   ---------------------------------------------------------------------------

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

  SUBROUTINE unit_matrix(a)

!   Purpose: Set the matrix a to be a unit matrix.

!   History: - Creation (16.10.1998,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a

!   ---------------------------------------------------------------------------

    a(:,:) = 0.0_dp

    CALL set_diag(a,1.0_dp)

  END SUBROUTINE unit_matrix

! *****************************************************************************

  FUNCTION vector_product(a,b) RESULT(c)

!   Purpose: Calculation of the vector product c = a x b.

!   History: - Creation (16.10.1998,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a, b
    REAL(KIND=dp), DIMENSION(3)              :: c

!   ---------------------------------------------------------------------------

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  END FUNCTION vector_product

! *****************************************************************************
!!****f* mathlib/gcd [1.0] *
!!
!!   NAME
!!     gcd
!!
!!   SYNOPSIS
!!     FUNCTION gcd(a,b)
!!        INTEGER :: a,b,gcd
!!     END FUNCTION
!!
!!   FUNCTION
!!     computes the greatest common divisor of two number
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  FUNCTION gcd(a,b)
    INTEGER                                  :: a, b, gcd

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
     RETURN
  END FUNCTION gcd

!!****f* mathlib/lcm [1.0] *
!!
!!   NAME
!!     lcm
!!
!!   SYNOPSIS
!!     FUNCTION lcm(a,b)
!!        integer :: lcm,a,b
!!     END FUNCTION
!!
!!   FUNCTION
!!     computes the least common multiplier of two numbers
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  FUNCTION lcm(a,b)
    INTEGER                                  :: a, b, lcm

    INTEGER                                  :: tmp

     tmp = gcd(a,b)
     IF (tmp==0) THEN
        lcm = 0
     ELSE
        ! could overflow
        lcm = ABS(a*b)/tmp
     END IF

     RETURN
  END FUNCTION lcm

! *****************************************************************************
!****f* mathlib/ei [1.0] *
!!
!!   NAME
!!     ei
!!
!!   FUNCTION
!!     computes the exponential integral
!!     Ei(x) = Int(exp(-x*t)/t,t=1..infinity)  x>0
!!
!!   AUTHOR
!!     JGH (adapted from Numerical recipies)
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  FUNCTION ei(x)
    REAL(dp)                                 :: x, ei

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'ei'
    INTEGER, PARAMETER                       :: maxit = 100
    REAL(dp), PARAMETER                      :: eps = 6.e-14_dp, &
                                                fpmin = 1.e-100_dp

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

END MODULE mathlib
