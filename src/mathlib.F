MODULE mathlib

! Purpose: Collection of simple mathematical functions and subroutines.

! History: Creation (15.11.1998, Matthias Krack)

! *****************************************************************************

! Index:

! FUNCTION angle(a,b)
! FUNCTION binomial_coefficient(n,k) RESULT(n_over_k)
! SUBROUTINE build_rotmat(phi,a,rotmat)
! SUBROUTINE diamat(a,eigvec,eigval)
! SUBROUTINE diamat_all(a,eigval,dac)
! FUNCTION gauss_exponent(l,radius,threshold,prefactor) RESULT(exponent)
! FUNCTION get_diag(a) RESULT(a_diag)
! SUBROUTINE invert_matrix(a,a_inverse,error,option)
! SUBROUTINE power_matrix(a,a_power,exponent,threshold,n_dependent,dac)
! FUNCTION reflect_vector(a,b) RESULT(a_mirror)
! FUNCTION rotate_vector(a,phi,b) RESULT(a_rot)
! SUBROUTINE set_diag_scalar(a,b)
! SUBROUTINE set_diag_vector(a,b)
! SUBROUTINE symmetrize_matrix(a,option)
! SUBROUTINE unit_matrix(a)
! FUNCTION vector_product(a,b) RESULT(c)

! *****************************************************************************

  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fac
  USE output_utilities,                ONLY: print_message,&
                                             print_warning
  USE termination,                     ONLY: get_error_unit,&
                                             stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PUBLIC

  REAL(KIND = dp), PARAMETER :: eps_geo = 1.0E-6_dp

! *****************************************************************************

  INTERFACE OPERATOR (.cross.)
    MODULE PROCEDURE vector_product
  END INTERFACE

  INTERFACE OPERATOR (.over.)
    MODULE PROCEDURE binomial_coefficient
  END INTERFACE

  INTERFACE diagonalize_matrix
    MODULE PROCEDURE diamat,diamat_all
  END INTERFACE

  INTERFACE set_diag
    MODULE PROCEDURE set_diag_scalar,set_diag_vector
  END INTERFACE

! *****************************************************************************

CONTAINS

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
    integer                                  :: a, b, lcm

    integer                                  :: tmp

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

  FUNCTION angle(a,b)

!   Purpose: Calculation of the angle between the vectors a and b. The
!            angle is returned in radians.

!   History: - Creation (14.10.1998, Matthias Krack)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a, b
    REAL(KIND=dp)                            :: angle

    REAL(KIND=dp)                            :: length_of_a, length_of_b
    REAL(KIND=dp), DIMENSION(SIZE(a, 1))     :: a_norm, b_norm

!   ---------------------------------------------------------------------------

    length_of_a = SQRT(DOT_PRODUCT(a(:),a(:)))
    length_of_b = SQRT(DOT_PRODUCT(b(:),b(:)))

    IF ((length_of_a > eps_geo).AND.(length_of_b > eps_geo)) THEN
      a_norm(:) = a(:)/length_of_a
      b_norm(:) = b(:)/length_of_b
      angle = ACOS(MIN(MAX(DOT_PRODUCT(a_norm(:),b_norm(:)),-1.0_dp),1.0_dp))
    ELSE
      angle = 0.0_dp
    END IF

  END FUNCTION angle

! *****************************************************************************

  FUNCTION binomial_coefficient(n,k) RESULT(n_over_k)

!   Purpose: The binomial coefficient n over k for 0 <= k <= n is calculated,
!            otherwise zero is returned.

!   History: - Creation (08.03.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: n, k
    REAL(KIND=dp)                            :: n_over_k

!   ---------------------------------------------------------------------------

    IF ((k >= 0).AND.(k <= n)) THEN
      n_over_k = fac(n)/(fac(n-k)*fac(k))
    ELSE
      n_over_k = 0.0_dp
    END IF

  END FUNCTION binomial_coefficient

! *****************************************************************************

  SUBROUTINE build_rotmat(phi,a,rotmat)

!   Purpose: The rotation matrix rotmat which rotates a vector about a
!            rotation axis defined by the vector a is build up.
!            The rotation angle is phi (radians).

!   History: - Creation (16.10.1998, Matthias Krack)

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

      CALL unit_matrix(rotmat(:,:))

    END IF

  END SUBROUTINE build_rotmat

! *****************************************************************************

  SUBROUTINE diamat(a,eigvec,eigval)

!   Purpose: Diagonalize the symmetric n by n matrix a using the LAPACK
!            library. The second dimension of eigvec determines the number
!            of computed eigenvectors and eigenvalues.

!   History: - Creation (29.10.1999, Matthias Krack)

!   ***************************************************************************

!   a       : Symmetric matrix to be diagonalized (input; upper triangle) ->
!             destroyed on exit.
!   eigval  : Eigenvalues of the matrix a (output).
!   eigvec  : Eigenvectors of the matrix a (output).

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: eigvec
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigval

    INTEGER                                  :: info, istat, liwork, lwork, &
                                                m, n, n_eigvec, nb
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ifail, iwork
    INTEGER, EXTERNAL                        :: ilaenv
    REAL(KIND=dp)                            :: abstol
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: work
    REAL(KIND=dp), EXTERNAL                  :: dlamch

!   *** Externals (LAPACK) ***

    EXTERNAL dsyevx

!   ---------------------------------------------------------------------------

!   *** Get the size of the matrix a ***

    n = SIZE(a,1)

!   *** Check the size of matrix a ***

    IF (SIZE(a,2) /= n) THEN
      CALL stop_program("SUBROUTINE diagonalize_matrix (MODULE mathlib)",&
                        "Check the size of matrix a (parameter #1)")
    END IF

!   *** Get the size of the second dimension of matrix eigvec,   ***
!   *** i.e. the requested number of eigenvectors to be computed ***

    n_eigvec = SIZE(eigvec,2)

!   *** Check the size of matrix eigvec ***

    IF (SIZE(eigvec,1) /= n) THEN
      CALL stop_program("SUBROUTINE diagonalize_matrix (MODULE mathlib)",&
                        "Check the size of matrix eigvec (parameter #2)")
    END IF

    IF (n_eigvec < 1) THEN
      CALL stop_program("SUBROUTINE diagonalize_matrix (MODULE mathlib)",&
                        "The number of requested eigenvectors is lower "//&
                        "than the one")
    ELSE IF (n_eigvec > n) THEN
      CALL stop_program("SUBROUTINE diagonalize_matrix (MODULE mathlib)",&
                        "The number of requested eigenvectors is larger "//&
                        "than the order of the matrix a")
    END IF

!   *** Check the size of vector eigval ***

    IF (SIZE(eigval) /= n_eigvec) THEN
      CALL stop_program("SUBROUTINE diagonalize_matrix (MODULE mathlib)",&
                        "The size of vector eigval does not match")
    END IF

!   *** Get the optimal work storage size ***

    nb = MAX(ilaenv(1,"DSYTRD","U",n,-1,-1,-1),&
             ilaenv(1,"DORMTR","U",n,-1,-1,-1))

    lwork = MAX((nb + 3)*n,8*n)+n ! sun bug fix
    liwork = 5*n

!   *** Set the absolute error tolerance for the eigenvalues ***

    abstol = 2.0_dp*dlamch("S")

!   *** Allocate work storage ***

    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE diagonalize_matrix "//&
                                     "(MODULE mathlib)","work",lwork)
    ALLOCATE (iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE diagonalize_matrix "//&
                                     "(MODULE mathlib)","iwork",liwork)
    ALLOCATE (ifail(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE diagonalize_matrix "//&
                                     "(MODULE mathlib)","ifail",n)

!   *** Diagonalize the matrix a ***

    CALL dsyevx("V","I","U",n,a(:,:),n,0.0_dp,1.0_dp,1,n_eigvec,abstol,m,&
                eigval(:),eigvec(:,:),n,work(:),lwork,iwork(:),ifail(:),info)

    IF (info /= 0) THEN
      CALL stop_program("SUBROUTINE diagonalize_matrix (MODULE mathlib)",&
                        "The matrix diagonalization in dsyevx failed")
    END IF

!   *** Release work storage ***

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE diagonalize_matrix "//&
                                     "(MODULE mathlib)","work")
    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE diagonalize_matrix "//&
                                     "(MODULE mathlib)","iwork")
    DEALLOCATE (ifail,STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE diagonalize_matrix "//&
                                     "(MODULE mathlib)","ifail")

  END SUBROUTINE diamat

! *****************************************************************************

  SUBROUTINE diamat_all(a,eigval,dac)

!   Purpose: Diagonalize the symmetric n by n matrix a using the LAPACK
!            library. Only the upper triangle of matrix a is used.

!   History: - Creation (29.03.1999, Matthias Krack)

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

    INTEGER                                  :: info, istat, liwork, lwork, &
                                                n, nb
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iwork
    INTEGER, EXTERNAL                        :: ilaenv
    LOGICAL                                  :: divide_and_conquer
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: work

!   *** Externals (LAPACK) ***
!MK EXTERNAL dsyev,dsyevd

    EXTERNAL dsyev

!   ---------------------------------------------------------------------------

!   *** Get the size of the matrix a ***

    n = SIZE(a,1)

!   *** Check the size of matrix a ***

    IF (SIZE(a,2) /= n) THEN
      CALL stop_program("SUBROUTINE diagonalize_matrix (MODULE mathlib)",&
                        "Check the size of matrix a (parameter #1)")
    END IF

!   *** Check the size of vector eigval ***

    IF (SIZE(eigval) /= n) THEN
      CALL stop_program("SUBROUTINE diagonalize_matrix (MODULE mathlib)",&
                        "The dimension of vector eigval is too small")
    END IF

!   *** Check, if the divide-and-conquer algorithm is requested ***

!MK IF (PRESENT(dac)) THEN
!     divide_and_conquer = dac
!   ELSE
      divide_and_conquer = .FALSE.
!   END IF

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
    IF (istat /= 0) CALL stop_memory("SUBROUTINE diagonalize_matrix "//&
                                     "(MODULE mathlib)","work",lwork)
    IF (divide_and_conquer) THEN
      ALLOCATE (iwork(liwork),STAT=istat)
      IF (istat /= 0) CALL stop_memory("SUBROUTINE diagonalize_matrix "//&
                                       "(MODULE mathlib)","iwork",liwork)
    END IF

!   *** Diagonalize the matrix a ***

!MK IF (divide_and_conquer) THEN
!     CALL dsyevd("V","U",n,a(:,:),n,eigval(:),work(:),lwork,iwork(:),liwork,&
!                 info)
!   ELSE
      CALL dsyev("V","U",n,a(1,1),n,eigval(1),work(1),lwork,info)
!   END IF

    IF (info /= 0) THEN
      IF (divide_and_conquer) THEN
        CALL stop_program("SUBROUTINE diagonalize_matrix (MODULE mathlib)",&
                          "The matrix diagonalization with dsyevd failed")
      ELSE
        CALL stop_program("SUBROUTINE diagonalize_matrix (MODULE mathlib)",&
                          "The matrix diagonalization with dsyev failed")
      END IF
    END IF

!   *** Release work storage ***

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE diagonalize_matrix "//&
                                     "(MODULE mathlib)","work")
    IF (divide_and_conquer) THEN
      DEALLOCATE (iwork,STAT=istat)
      IF (istat /= 0) CALL stop_memory("SUBROUTINE diagonalize_matrix "//&
                                       "(MODULE mathlib)","iwork")
    END IF

  END SUBROUTINE diamat_all

! *****************************************************************************

  FUNCTION gauss_exponent(l,radius,threshold,prefactor) RESULT(exponent)

!   Purpose: The exponent of a primitive Gaussian function for a given radius
!            and threshold is calculated.

!   History: - Creation (07.03.1999, Matthias Krack)

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

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "FUNCTION gauss_exponent (MODULE mathlib)"

    REAL(KIND=dp)                            :: d, r, t

!   ---------------------------------------------------------------------------
!   *** Initialize function value ***

    exponent = 0.0_dp

!   *** Load and check parameter values ***

    IF (l < 0) THEN
      CALL stop_program(routine,"The angular momentum quantum number is "//&
                                "negative")
    END IF

    IF (radius == 0.0_dp) THEN
      CALL stop_program(routine,"The requested radius is zero")
    ELSE
      r = ABS(radius)
    END IF

    IF (threshold == 0.0_dp) THEN
      CALL stop_program(routine,"The requested threshold is zero")
    ELSE
      t = ABS(threshold)
    END IF

    IF (prefactor == 0.0_dp) THEN
      CALL stop_program(routine,"The requested prefactor is zero")
    ELSE
      d = ABS(prefactor)
    END IF

!   *** Calculate exponent ***

    exponent = LOG(d*r**l/t)/r**2

  END FUNCTION gauss_exponent

! *****************************************************************************

  FUNCTION get_diag(a) RESULT(a_diag)

!   Purpose: Return the diagonal elements of matrix a as a vector.

!   History: - Creation (20.11.1998, Matthias Krack)

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

  SUBROUTINE invert_matrix(a,a_inverse,error,option)

!   Purpose: Compute the inverse of the n by n matrix a using the LAPACK
!            library.

!   History: - Creation (23.03.1999, Matthias Krack)

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

    CHARACTER(LEN=1)                         :: norm, trans
    CHARACTER(LEN=80)                        :: message
    INTEGER                                  :: info, istat, iter, n, &
                                                output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipiv, iwork
    REAL(KIND=dp)                            :: a_norm, old_error, r_cond
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: berr, ferr, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: a_lu, b

!   *** Externals (LAPACK) ***
!MK    REAL(KIND = dp), EXTERNAL :: dlange

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

    IF (n /= SIZE(a,2)) THEN
      CALL stop_program("SUBROUTINE invert_matrix (MODULE mathlib)",&
                        "Check the array bounds of parameter #1")
    END IF

    IF ((n /= SIZE(a_inverse,1)).OR.&
        (n /= SIZE(a_inverse,2))) THEN
      CALL stop_program("SUBROUTINE invert_matrix (MODULE mathlib)",&
                        "Check the array bounds of parameter #2")
    END IF

!   *** Allocate work storage ***

    ALLOCATE (a_lu(n,n),STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","a_lu",n*n)

    ALLOCATE (b(n,n),STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","b",n*n)

    ALLOCATE (berr(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","berr",n)

    ALLOCATE (ferr(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","ferr",n)

    ALLOCATE (ipiv(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","ipiv",n)

    ALLOCATE (iwork(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","iwork",n)

    ALLOCATE (work(4*n),STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","work",4*n)

    a_lu(1:n,1:n) = a(1:n,1:n)

!   *** Compute the LU factorization of the matrix a ***

    CALL dgetrf(n,n,a_lu(:,:),n,ipiv(:),info)

    IF (info /= 0) THEN
      CALL stop_program("SUBROUTINE invert_matrix (MODULE mathlib)",&
                        "The LU factorization in dgetrf failed")
    END IF

!   *** Compute the norm of the matrix a ***

    IF (trans == "N") THEN
      norm = '1'
    ELSE
      norm = 'I'
    END IF

!MK    a_norm = dlange(norm,n,n,a(:,:),n,work(:))

!   *** Compute the reciprocal of the condition number of a ***

    CALL dgecon(norm,n,a_lu(:,:),n,a_norm,r_cond,work(:),iwork(:),info)

    IF (info /= 0) THEN
      CALL stop_program("SUBROUTINE invert_matrix (MODULE mathlib)",&
                        "The computation of the condition number in "//&
                        "dgecon failed")
    END IF

    IF (r_cond < EPSILON(0.0_dp)) THEN
      CALL get_error_unit(output_unit)
      WRITE (message,"(A,ES10.3,A)") "R_COND =",r_cond
      CALL print_message(message,output_unit,1,1,1)
      CALL print_warning("SUBROUTINE invert_matrix (MODULE mathlib)",&
                         "Bad condition number (smaller than the machine "//&
                         "working precision)",output_unit)
    END IF

!   *** Solve a system of linear equations using ***
!   *** the LU factorization computed by dgetrf  ***

    CALL unit_matrix(a_inverse(:,:))

    CALL dgetrs(trans,n,n,a_lu(:,:),n,ipiv(:),a_inverse(:,:),n,info)

    IF (info /= 0) THEN
      CALL stop_program("SUBROUTINE invert_matrix (MODULE mathlib)",&
                        "Solving the system of linear equations in dgetrs "//&
                        "failed")
    END IF

!   *** Improve the computed solution iteratively ***

    CALL unit_matrix(b) ! Initialize right-hand sides

    error = 0.0_dp

    DO iter=1,10

      CALL dgerfs(trans,n,n,a(:,:),n,a_lu(:,:),n,ipiv(:),b(:,:),n,&
                  a_inverse(:,:),n,ferr(:),berr(:),work(:),iwork(:),info)

      IF (info /= 0) THEN
        CALL stop_program("SUBROUTINE invert_matrix (MODULE mathlib)",&
                          "Improving the computed solution in dgerfs failed")
      END IF

      old_error = error
      error = MAXVAL(ferr(:))

      IF (ABS(error - old_error) <= EPSILON(0.0_dp)) EXIT

    END DO

!   *** Release work storage ***

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","work")
    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","iwork")
    DEALLOCATE (ipiv,STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","ipiv")
    DEALLOCATE (ferr,STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","ferr")
    DEALLOCATE (berr,STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","berr")
    DEALLOCATE (b,STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","b")
    DEALLOCATE (a_lu,STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE invert_matrix "//&
                                     "(MODULE mathlib)","a_lu")

  END SUBROUTINE invert_matrix

! *****************************************************************************

  SUBROUTINE power_matrix(a,a_power,exponent,threshold,n_dependent,dac)

!   Purpose: Raise the real symmetric n by n matrix a to the power given by
!            exponent. All eigenvectors with a corresponding eigenvalue lower
!            than threshold are quenched. Only the upper triangle of matrix a
!            is used.

!   History: - Creation (29.03.1999, Matthias Krack)

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
      CALL stop_program("SUBROUTINE power_matrix (MODULE mathlib)",&
                        "Check the array bounds of parameter #1")
    END IF

    IF ((n /= SIZE(a_power,1)).OR.&
        (n /= SIZE(a_power,2))) THEN
      CALL stop_program("SUBROUTINE power_matrix (MODULE mathlib)",&
                        "Check the array bounds of parameter #2")
    END IF

!   *** Allocate work storage ***

    ALLOCATE (eigval(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE power_matrix "//&
                                     "(MODULE mathlib)","eigval",n)

!   *** Check, if the divide-and-conquer algorithm is requested ***

    IF (PRESENT(dac)) THEN
!MK   divide_and_conquer = dac
      divide_and_conquer = .FALSE.
    ELSE
      divide_and_conquer = .FALSE.
    END IF

!   *** Compute the eigenvectors and eigenvalues of the matrix a ***

    CALL diagonalize_matrix(a(:,:),eigval(:),divide_and_conquer)

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

    CALL symmetrize_matrix(a_power(:,:),"upper_to_lower")

!   *** Release work storage ***

    DEALLOCATE (eigval,STAT=istat)
    IF (istat /= 0) CALL stop_memory("SUBROUTINE power_matrix "//&
                                     "(MODULE mathlib)","eigval")

  END SUBROUTINE power_matrix

! *****************************************************************************

  FUNCTION reflect_vector(a,b) RESULT(a_mirror)

!   Purpose: Reflection of the vector a through a mirror plane defined by the
!            normal vector b. The reflected vector a is stored in a_mirror.

!   History: - Creation (16.10.1998, Matthias Krack)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a, b
    REAL(KIND=dp), DIMENSION(3)              :: a_mirror

    REAL(KIND=dp)                            :: length_of_b, scapro
    REAL(KIND=dp), DIMENSION(3)              :: d

!   ---------------------------------------------------------------------------
!   *** Check the length of the vector b ***

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

!   History: - Creation (16.10.1998, Matthias Krack)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a
    REAL(KIND=dp), INTENT(IN)                :: phi
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: b
    REAL(KIND=dp), DIMENSION(3)              :: a_rot

    REAL(KIND=dp)                            :: length_of_b
    REAL(KIND=dp), DIMENSION(3, 3)           :: rotmat

!   ---------------------------------------------------------------------------
!   *** Check the length of vector the b ***

    length_of_b = SQRT(b(1)*b(1) + b(2)*b(2) + b(3)*b(3))

    IF (length_of_b > eps_geo) THEN

!     *** Build up the rotation matrix rotmat ***

      CALL build_rotmat(phi,b(:),rotmat(:,:))

!     *** Rotate the vector a by phi about the axis defined by vector b ***

      a_rot(:) = MATMUL(rotmat(:,:),a(:))

    ELSE

      a_rot(:) = 0.0_dp

    END IF

  END FUNCTION rotate_vector

! *****************************************************************************

  SUBROUTINE set_diag_scalar(a,b)

!   Purpose: Set the diagonal elements of matrix a to b.

!   History: - Creation (20.11.1998, Matthias Krack)

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

!   History: - Creation (20.11.1998, Matthias Krack)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: b

    INTEGER                                  :: i, n

!   ---------------------------------------------------------------------------

    n = MIN(SIZE(a,1),SIZE(a,2))

    IF (SIZE(b) == n) THEN
      DO i=1,n
        a(i,i) = b(i)
      END DO
    ELSE
      CALL stop_program("SUBROUTINE set_diag_vector (MODULE mathlib)",&
                        "Check the array bounds of the parameters")
    END IF

  END SUBROUTINE set_diag_vector

! *****************************************************************************

  SUBROUTINE symmetrize_matrix(a,option)

!   Purpose: Symmetrize the matrix a.

!   History: - Creation (16.10.1998, Matthias Krack)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    CHARACTER(LEN=*), INTENT(IN)             :: option

    CHARACTER(LEN=80)                        :: message
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
      CALL stop_program("SUBROUTNE symmetrize_matrix (MODULE mathlib)",message)
    END IF

  END SUBROUTINE symmetrize_matrix

! *****************************************************************************

  SUBROUTINE unit_matrix(a)

!   Purpose: Set the matrix a to be a unit matrix.

!   History: - Creation (16.10.1998, Matthias Krack)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a

!   ---------------------------------------------------------------------------

    a(:,:) = 0.0_dp

    CALL set_diag(a(:,:),1.0_dp)

  END SUBROUTINE unit_matrix

! *****************************************************************************

  FUNCTION vector_product(a,b) RESULT(c)

!   Purpose: Calculation of the vector product c = a x b.

!   History: - Creation (16.10.1998, Matthias Krack)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: a, b
    REAL(KIND=dp), DIMENSION(3)              :: c

!   ---------------------------------------------------------------------------

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  END FUNCTION vector_product

! *****************************************************************************

END MODULE mathlib
