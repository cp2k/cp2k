!!****** cp2k/nrutil [1.0] *
!!
!!   NAME
!!     nrutil
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
    MODULE nrutil
    USE kinds,                           ONLY: dbl
      IMPLICIT NONE
    PRIVATE
    PUBLIC :: swap
      INTEGER, PARAMETER :: i4b = kind(0)
      INTEGER, PARAMETER :: sp = kind(1._dbl)
      INTEGER, PARAMETER :: spc = kind((1._dbl,1._dbl))
      INTEGER, PARAMETER :: lgt = kind(.TRUE.)
      REAL (sp), PARAMETER :: twopi = &
        6.283185307179586476925286766559005768394_sp
      INTEGER (i4b), PARAMETER :: npar_arth = 16, npar2_arth = 8
      INTEGER (i4b), PARAMETER :: npar_geop = 4, npar2_geop = 2
      INTEGER (i4b), PARAMETER :: npar_cumsum = 16
      INTEGER (i4b), PARAMETER :: npar_cumprod = 8
      INTEGER (i4b), PARAMETER :: npar_poly = 8
      INTEGER (i4b), PARAMETER :: npar_polyterm = 8
      INTERFACE array_copy
        MODULE PROCEDURE array_copy_r, array_copy_i
      END INTERFACE
      INTERFACE swap
        MODULE PROCEDURE swap_i, swap_r, swap_rv, swap_c, swap_cv, swap_cm, &
          masked_swap_rs, masked_swap_rv, masked_swap_rm
      END INTERFACE
      INTERFACE reallocate
        MODULE PROCEDURE reallocate_rv, reallocate_rm, reallocate_iv, &
          reallocate_im, reallocate_hv
      END INTERFACE
      INTERFACE imaxloc
        MODULE PROCEDURE imaxloc_r, imaxloc_i
      END INTERFACE
      INTERFACE assert
        MODULE PROCEDURE assert1, assert2, assert3, assert4, assert_v
      END INTERFACE
      INTERFACE assert_eq
        MODULE PROCEDURE assert_eq2, assert_eq3, assert_eq4, assert_eqn
      END INTERFACE
      INTERFACE arth
        MODULE PROCEDURE arth_r, arth_i
      END INTERFACE
      INTERFACE geop
        MODULE PROCEDURE geop_r, geop_i, geop_c
      END INTERFACE
      INTERFACE cumsum
        MODULE PROCEDURE cumsum_r, cumsum_i
      END INTERFACE
      INTERFACE poly
        MODULE PROCEDURE poly_rr, poly_rrv, poly_rc, poly_cc, poly_msk_rrv
      END INTERFACE
      INTERFACE poly_term
        MODULE PROCEDURE poly_term_rr, poly_term_cc
      END INTERFACE
      INTERFACE outerprod
        MODULE PROCEDURE outerprod_r
      END INTERFACE
      INTERFACE outerdiff
        MODULE PROCEDURE outerdiff_r, outerdiff_i
      END INTERFACE
      INTERFACE scatter_add
        MODULE PROCEDURE scatter_add_r
      END INTERFACE
      INTERFACE scatter_max
        MODULE PROCEDURE scatter_max_r
      END INTERFACE
      INTERFACE diagadd
        MODULE PROCEDURE diagadd_rv, diagadd_r
      END INTERFACE
      INTERFACE diagmult
        MODULE PROCEDURE diagmult_rv, diagmult_r
      END INTERFACE
      INTERFACE get_diag
        MODULE PROCEDURE get_diag_rv
      END INTERFACE
      INTERFACE put_diag
        MODULE PROCEDURE put_diag_rv, put_diag_r
      END INTERFACE

!!*****
    CONTAINS
!BL
      SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: src
    REAL(sp), DIMENSION(:), INTENT(OUT)      :: dest
    INTEGER(i4b), INTENT(OUT)                :: n_copied, n_not_copied

        n_copied = min(size(src),size(dest))
        n_not_copied = size(src) - n_copied
        dest(1:n_copied) = src(1:n_copied)
      END SUBROUTINE array_copy_r
!BL
      SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
    INTEGER(i4b), DIMENSION(:), INTENT(IN)   :: src
    INTEGER(i4b), DIMENSION(:), INTENT(OUT)  :: dest
    INTEGER(i4b), INTENT(OUT)                :: n_copied, n_not_copied

        n_copied = min(size(src),size(dest))
        n_not_copied = size(src) - n_copied
        dest(1:n_copied) = src(1:n_copied)
      END SUBROUTINE array_copy_i
!BL
!BL
      SUBROUTINE swap_i(a,b)
    INTEGER(i4b), INTENT(INOUT)              :: a, b

    INTEGER(i4b)                             :: dum

        dum = a
        a = b
        b = dum
      END SUBROUTINE swap_i
!BL
      SUBROUTINE swap_r(a,b)
    REAL(sp), INTENT(INOUT)                  :: a, b

    REAL(sp)                                 :: dum

        dum = a
        a = b
        b = dum
      END SUBROUTINE swap_r
!BL
      SUBROUTINE swap_rv(a,b)
    REAL(sp), DIMENSION(:), INTENT(INOUT)    :: a, b

    REAL(sp), DIMENSION(size(a))             :: dum

        dum = a
        a = b
        b = dum
      END SUBROUTINE swap_rv
!BL
      SUBROUTINE swap_c(a,b)
    COMPLEX(spc), INTENT(INOUT)              :: a, b

    COMPLEX(spc)                             :: dum

        dum = a
        a = b
        b = dum
      END SUBROUTINE swap_c
!BL
      SUBROUTINE swap_cv(a,b)
    COMPLEX(spc), DIMENSION(:), &
      INTENT(INOUT)                          :: a, b

    COMPLEX(spc), DIMENSION(size(a))         :: dum

        dum = a
        a = b
        b = dum
      END SUBROUTINE swap_cv
!BL
      SUBROUTINE swap_cm(a,b)
    COMPLEX(spc), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a, b

    COMPLEX(spc), &
      DIMENSION(size(a, 1), size(a, 2))      :: dum

        dum = a
        a = b
        b = dum
      END SUBROUTINE swap_cm
!BL
      SUBROUTINE masked_swap_rs(a,b,mask)
    REAL(sp), INTENT(INOUT)                  :: a, b
    LOGICAL(lgt), INTENT(IN)                 :: mask

    REAL(sp)                                 :: swp

        IF (mask) THEN
          swp = a
          a = b
          b = swp
        END IF
      END SUBROUTINE masked_swap_rs
!BL
      SUBROUTINE masked_swap_rv(a,b,mask)
    REAL(sp), DIMENSION(:), INTENT(INOUT)    :: a, b
    LOGICAL(lgt), DIMENSION(:), INTENT(IN)   :: mask

    REAL(sp), DIMENSION(size(a))             :: swp

        WHERE (mask)
          swp = a
          a = b
          b = swp
        END WHERE
      END SUBROUTINE masked_swap_rv
!BL
      SUBROUTINE masked_swap_rm(a,b,mask)
    REAL(sp), DIMENSION(:, :), INTENT(INOUT) :: a, b
    LOGICAL(lgt), DIMENSION(:, :), &
      INTENT(IN)                             :: mask

    REAL(sp), &
      DIMENSION(size(a, 1), size(a, 2))      :: swp

        WHERE (mask)
          swp = a
          a = b
          b = swp
        END WHERE
      END SUBROUTINE masked_swap_rm
!BL
!BL
      FUNCTION reallocate_rv(p,n)
    REAL(sp), DIMENSION(:), POINTER          :: p
    INTEGER(i4b), INTENT(IN)                 :: n
    REAL(sp), DIMENSION(:), POINTER          :: reallocate_rv

    INTEGER(i4b)                             :: ierr, nold

        ALLOCATE (reallocate_rv(n),STAT=ierr)
        IF (ierr/=0) CALL nrerror( &
          'reallocate_rv: problem in attempt to allocate memory')
        IF ( .NOT. associated(p)) RETURN
        nold = size(p)
        reallocate_rv(1:min(nold,n)) = p(1:min(nold,n))
        DEALLOCATE (p,STAT=ierr)
        IF (ierr/=0) CALL nrerror( &
          'reallocate_rv: problem in attempt to deallocate memory')
      END FUNCTION reallocate_rv
!BL
      FUNCTION reallocate_iv(p,n)
    INTEGER(i4b), DIMENSION(:), POINTER      :: p
    INTEGER(i4b), INTENT(IN)                 :: n
    INTEGER(i4b), DIMENSION(:), POINTER      :: reallocate_iv

    INTEGER(i4b)                             :: ierr, nold

        ALLOCATE (reallocate_iv(n),STAT=ierr)
        IF (ierr/=0) CALL nrerror( &
          'reallocate_iv: problem in attempt to allocate memory')
        IF ( .NOT. associated(p)) RETURN
        nold = size(p)
        reallocate_iv(1:min(nold,n)) = p(1:min(nold,n))
        DEALLOCATE (p,STAT=ierr)
        IF (ierr/=0) CALL nrerror( &
          'reallocate_iv: problem in attempt to deallocate memory')
      END FUNCTION reallocate_iv
!BL
      FUNCTION reallocate_hv(p,n)
    CHARACTER(1), DIMENSION(:), POINTER      :: p
    INTEGER(i4b), INTENT(IN)                 :: n
    CHARACTER(1), DIMENSION(:), POINTER      :: reallocate_hv

    INTEGER(i4b)                             :: ierr, nold

        ALLOCATE (reallocate_hv(n),STAT=ierr)
        IF (ierr/=0) CALL nrerror( &
          'reallocate_hv: problem in attempt to allocate memory')
        IF ( .NOT. associated(p)) RETURN
        nold = size(p)
        reallocate_hv(1:min(nold,n)) = p(1:min(nold,n))
        DEALLOCATE (p,STAT=ierr)
        IF (ierr/=0) CALL nrerror( &
          'reallocate_hv: problem in attempt to deallocate memory')
      END FUNCTION reallocate_hv
!BL
      FUNCTION reallocate_rm(p,n,m)
    REAL(sp), DIMENSION(:, :), POINTER       :: p
    INTEGER(i4b), INTENT(IN)                 :: n, m
    REAL(sp), DIMENSION(:, :), POINTER       :: reallocate_rm

    INTEGER(i4b)                             :: ierr, mold, nold

        ALLOCATE (reallocate_rm(n,m),STAT=ierr)
        IF (ierr/=0) CALL nrerror( &
          'reallocate_rm: problem in attempt to allocate memory')
        IF ( .NOT. associated(p)) RETURN
        nold = size(p,1)
        mold = size(p,2)
        reallocate_rm(1:min(nold,n),1:min(mold,m)) = p(1:min(nold,n), &
          1:min(mold,m))
        DEALLOCATE (p,STAT=ierr)
        IF (ierr/=0) CALL nrerror( &
          'reallocate_rm: problem in attempt to deallocate memory')
      END FUNCTION reallocate_rm
!BL
      FUNCTION reallocate_im(p,n,m)
    INTEGER(i4b), DIMENSION(:, :), POINTER   :: p
    INTEGER(i4b), INTENT(IN)                 :: n, m
    INTEGER(i4b), DIMENSION(:, :), POINTER   :: reallocate_im

    INTEGER(i4b)                             :: ierr, mold, nold

        ALLOCATE (reallocate_im(n,m),STAT=ierr)
        IF (ierr/=0) CALL nrerror( &
          'reallocate_im: problem in attempt to allocate memory')
        IF ( .NOT. associated(p)) RETURN
        nold = size(p,1)
        mold = size(p,2)
        reallocate_im(1:min(nold,n),1:min(mold,m)) = p(1:min(nold,n), &
          1:min(mold,m))
        DEALLOCATE (p,STAT=ierr)
        IF (ierr/=0) CALL nrerror( &
          'reallocate_im: problem in attempt to deallocate memory')
      END FUNCTION reallocate_im
!BL
      FUNCTION ifirstloc(mask)
    LOGICAL(lgt), DIMENSION(:), INTENT(IN)   :: mask
    INTEGER(i4b)                             :: ifirstloc

    INTEGER(i4b), DIMENSION(1)               :: loc

        loc = maxloc(merge(1,0,mask))
        ifirstloc = loc(1)
        IF ( .NOT. mask(ifirstloc)) ifirstloc = size(mask) + 1
      END FUNCTION ifirstloc
!BL
      FUNCTION imaxloc_r(arr)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: arr
    INTEGER(i4b)                             :: imaxloc_r

    INTEGER(i4b), DIMENSION(1)               :: imax

        imax = maxloc(arr(:))
        imaxloc_r = imax(1)
      END FUNCTION imaxloc_r
!BL
      FUNCTION imaxloc_i(iarr)
    INTEGER(i4b), DIMENSION(:), INTENT(IN)   :: iarr
    INTEGER(i4b)                             :: imaxloc_i

    INTEGER(i4b), DIMENSION(1)               :: imax

        imax = maxloc(iarr(:))
        imaxloc_i = imax(1)
      END FUNCTION imaxloc_i
!BL
      FUNCTION iminloc(arr)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: arr
    INTEGER(i4b)                             :: iminloc

    INTEGER(i4b), DIMENSION(1)               :: imin

        imin = minloc(arr(:))
        iminloc = imin(1)
      END FUNCTION iminloc
!BL
      SUBROUTINE assert1(n1,string)
    LOGICAL, INTENT(IN)                      :: n1
    CHARACTER(len=*), INTENT(IN)             :: string

        IF ( .NOT. n1) THEN
          WRITE (*,*) 'nrerror: an assertion failed with this tag:', string
          STOP 'program terminated by assert1'
        END IF
      END SUBROUTINE assert1
!BL
      SUBROUTINE assert2(n1,n2,string)
    LOGICAL, INTENT(IN)                      :: n1, n2
    CHARACTER(len=*), INTENT(IN)             :: string

        IF ( .NOT. (n1 .AND. n2)) THEN
          WRITE (*,*) 'nrerror: an assertion failed with this tag:', string
          STOP 'program terminated by assert2'
        END IF
      END SUBROUTINE assert2
!BL
      SUBROUTINE assert3(n1,n2,n3,string)
    LOGICAL, INTENT(IN)                      :: n1, n2, n3
    CHARACTER(len=*), INTENT(IN)             :: string

        IF ( .NOT. (n1 .AND. n2 .AND. n3)) THEN
          WRITE (*,*) 'nrerror: an assertion failed with this tag:', string
          STOP 'program terminated by assert3'
        END IF
      END SUBROUTINE assert3
!BL
      SUBROUTINE assert4(n1,n2,n3,n4,string)
    LOGICAL, INTENT(IN)                      :: n1, n2, n3, n4
    CHARACTER(len=*), INTENT(IN)             :: string

        IF ( .NOT. (n1 .AND. n2 .AND. n3 .AND. n4)) THEN
          WRITE (*,*) 'nrerror: an assertion failed with this tag:', string
          STOP 'program terminated by assert4'
        END IF
      END SUBROUTINE assert4
!BL
      SUBROUTINE assert_v(n,string)
    LOGICAL, DIMENSION(:), INTENT(IN)        :: n
    CHARACTER(len=*), INTENT(IN)             :: string

        IF ( .NOT. all(n)) THEN
          WRITE (*,*) 'nrerror: an assertion failed with this tag:', string
          STOP 'program terminated by assert_v'
        END IF
      END SUBROUTINE assert_v
!BL
      FUNCTION assert_eq2(n1,n2,string)
    INTEGER, INTENT(IN)                      :: n1, n2
    CHARACTER(len=*), INTENT(IN)             :: string
    INTEGER                                  :: assert_eq2

        IF (n1==n2) THEN
          assert_eq2 = n1
        ELSE
          WRITE (*,*) 'nrerror: an assert_eq failed with this tag:', string
          STOP 'program terminated by assert_eq2'
        END IF
      END FUNCTION assert_eq2
!BL
      FUNCTION assert_eq3(n1,n2,n3,string)
    INTEGER, INTENT(IN)                      :: n1, n2, n3
    CHARACTER(len=*), INTENT(IN)             :: string
    INTEGER                                  :: assert_eq3

        IF (n1==n2 .AND. n2==n3) THEN
          assert_eq3 = n1
        ELSE
          WRITE (*,*) 'nrerror: an assert_eq failed with this tag:', string
          STOP 'program terminated by assert_eq3'
        END IF
      END FUNCTION assert_eq3
!BL
      FUNCTION assert_eq4(n1,n2,n3,n4,string)
    INTEGER, INTENT(IN)                      :: n1, n2, n3, n4
    CHARACTER(len=*), INTENT(IN)             :: string
    INTEGER                                  :: assert_eq4

        IF (n1==n2 .AND. n2==n3 .AND. n3==n4) THEN
          assert_eq4 = n1
        ELSE
          WRITE (*,*) 'nrerror: an assert_eq failed with this tag:', string
          STOP 'program terminated by assert_eq4'
        END IF
      END FUNCTION assert_eq4
!BL
      FUNCTION assert_eqn(nn,string)
    INTEGER, DIMENSION(:), INTENT(IN)        :: nn
    CHARACTER(len=*), INTENT(IN)             :: string
    INTEGER                                  :: assert_eqn

        IF (all(nn(2:)==nn(1))) THEN
          assert_eqn = nn(1)
        ELSE
          WRITE (*,*) 'nrerror: an assert_eq failed with this tag:', string
          STOP 'program terminated by assert_eqn'
        END IF
      END FUNCTION assert_eqn
!BL
      SUBROUTINE nrerror(string)
    CHARACTER(len=*), INTENT(IN)             :: string

        WRITE (*,*) 'nrerror: ', string
        STOP 'program terminated by nrerror'
      END SUBROUTINE nrerror
!BL
      FUNCTION arth_r(first,increment,n)
    REAL(sp), INTENT(IN)                     :: first, increment
    INTEGER(i4b), INTENT(IN)                 :: n
    REAL(sp), DIMENSION(n)                   :: arth_r

    INTEGER(i4b)                             :: k, k2
    REAL(sp)                                 :: temp

        IF (n>0) arth_r(1) = first
        IF (n<=npar_arth) THEN
          DO k = 2, n
            arth_r(k) = arth_r(k-1) + increment
          END DO
        ELSE
          DO k = 2, npar2_arth
            arth_r(k) = arth_r(k-1) + increment
          END DO
          temp = increment*npar2_arth
          k = npar2_arth
          DO
            IF (k>=n) EXIT
            k2 = k + k
            arth_r(k+1:min(k2,n)) = temp + arth_r(1:min(k,n-k))
            temp = temp + temp
            k = k2
          END DO
        END IF
      END FUNCTION arth_r
!BL
      FUNCTION arth_i(first,increment,n)
    INTEGER(i4b), INTENT(IN)                 :: first, increment, n
    INTEGER(i4b), DIMENSION(n)               :: arth_i

    INTEGER(i4b)                             :: k, k2, temp

        IF (n>0) arth_i(1) = first
        IF (n<=npar_arth) THEN
          DO k = 2, n
            arth_i(k) = arth_i(k-1) + increment
          END DO
        ELSE
          DO k = 2, npar2_arth
            arth_i(k) = arth_i(k-1) + increment
          END DO
          temp = increment*npar2_arth
          k = npar2_arth
          DO
            IF (k>=n) EXIT
            k2 = k + k
            arth_i(k+1:min(k2,n)) = temp + arth_i(1:min(k,n-k))
            temp = temp + temp
            k = k2
          END DO
        END IF
      END FUNCTION arth_i
!BL
!BL
      FUNCTION geop_r(first,factor,n)
    REAL(sp), INTENT(IN)                     :: first, factor
    INTEGER(i4b), INTENT(IN)                 :: n
    REAL(sp), DIMENSION(n)                   :: geop_r

    INTEGER(i4b)                             :: k, k2
    REAL(sp)                                 :: temp

        IF (n>0) geop_r(1) = first
        IF (n<=npar_geop) THEN
          DO k = 2, n
            geop_r(k) = geop_r(k-1)*factor
          END DO
        ELSE
          DO k = 2, npar2_geop
            geop_r(k) = geop_r(k-1)*factor
          END DO
          temp = factor**npar2_geop
          k = npar2_geop
          DO
            IF (k>=n) EXIT
            k2 = k + k
            geop_r(k+1:min(k2,n)) = temp*geop_r(1:min(k,n-k))
            temp = temp*temp
            k = k2
          END DO
        END IF
      END FUNCTION geop_r
!BL
      FUNCTION geop_i(first,factor,n)
    INTEGER(i4b), INTENT(IN)                 :: first, factor, n
    INTEGER(i4b), DIMENSION(n)               :: geop_i

    INTEGER(i4b)                             :: k, k2, temp

        IF (n>0) geop_i(1) = first
        IF (n<=npar_geop) THEN
          DO k = 2, n
            geop_i(k) = geop_i(k-1)*factor
          END DO
        ELSE
          DO k = 2, npar2_geop
            geop_i(k) = geop_i(k-1)*factor
          END DO
          temp = factor**npar2_geop
          k = npar2_geop
          DO
            IF (k>=n) EXIT
            k2 = k + k
            geop_i(k+1:min(k2,n)) = temp*geop_i(1:min(k,n-k))
            temp = temp*temp
            k = k2
          END DO
        END IF
      END FUNCTION geop_i
!BL
      FUNCTION geop_c(first,factor,n)
    COMPLEX(sp), INTENT(IN)                  :: first, factor
    INTEGER(i4b), INTENT(IN)                 :: n
    COMPLEX(sp), DIMENSION(n)                :: geop_c

    COMPLEX(sp)                              :: temp
    INTEGER(i4b)                             :: k, k2

        IF (n>0) geop_c(1) = first
        IF (n<=npar_geop) THEN
          DO k = 2, n
            geop_c(k) = geop_c(k-1)*factor
          END DO
        ELSE
          DO k = 2, npar2_geop
            geop_c(k) = geop_c(k-1)*factor
          END DO
          temp = factor**npar2_geop
          k = npar2_geop
          DO
            IF (k>=n) EXIT
            k2 = k + k
            geop_c(k+1:min(k2,n)) = temp*geop_c(1:min(k,n-k))
            temp = temp*temp
            k = k2
          END DO
        END IF
      END FUNCTION geop_c
!BL
!BL
      RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT (ans)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: arr
    REAL(sp), INTENT(IN), OPTIONAL           :: seed
    REAL(sp), DIMENSION(size(arr))           :: ans

    INTEGER(i4b)                             :: j, n
    REAL(sp)                                 :: sd

        n = size(arr)
        IF (n==0_i4b) RETURN
        sd = 0.0_sp
        IF (present(seed)) sd = seed
        ans(1) = arr(1) + sd
        IF (n<npar_cumsum) THEN
          DO j = 2, n
            ans(j) = ans(j-1) + arr(j)
          END DO
        ELSE
          ans(2:n:2) = cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
          ans(3:n:2) = ans(2:n-1:2) + arr(3:n:2)
        END IF
      END FUNCTION cumsum_r
!BL
      RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT (ans)
    INTEGER(i4b), DIMENSION(:), INTENT(IN)   :: arr
    INTEGER(i4b), INTENT(IN), OPTIONAL       :: seed
    INTEGER(i4b), DIMENSION(size(arr))       :: ans

    INTEGER(i4b)                             :: j, n, sd

        n = size(arr)
        IF (n==0_i4b) RETURN
        sd =  0_i4b
        IF (present(seed)) sd = seed
        ans(1) = arr(1) + sd
        IF (n<npar_cumsum) THEN
          DO j = 2, n
            ans(j) = ans(j-1) + arr(j)
          END DO
        ELSE
          ans(2:n:2) = cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
          ans(3:n:2) = ans(2:n-1:2) + arr(3:n:2)
        END IF
      END FUNCTION cumsum_i
!BL
!BL
      RECURSIVE FUNCTION cumprod(arr,seed) RESULT (ans)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: arr
    REAL(sp), INTENT(IN), OPTIONAL           :: seed
    REAL(sp), DIMENSION(size(arr))           :: ans

    INTEGER(i4b)                             :: j, n
    REAL(sp)                                 :: sd

        n = size(arr)
        IF (n==0_i4b) RETURN
        sd = 1.0_sp
        IF (present(seed)) sd = seed
        ans(1) = arr(1)*sd
        IF (n<npar_cumprod) THEN
          DO j = 2, n
            ans(j) = ans(j-1)*arr(j)
          END DO
        ELSE
          ans(2:n:2) = cumprod(arr(2:n:2)*arr(1:n-1:2),sd)
          ans(3:n:2) = ans(2:n-1:2)*arr(3:n:2)
        END IF
      END FUNCTION cumprod
!BL
!BL
      FUNCTION poly_rr(x,coeffs)
    REAL(sp), INTENT(IN)                     :: x
    REAL(sp), DIMENSION(:), INTENT(IN)       :: coeffs
    REAL(sp)                                 :: poly_rr

    INTEGER(i4b)                             :: i, ierr, n, nn
    REAL(sp)                                 :: pow
    REAL(sp), ALLOCATABLE, DIMENSION(:)      :: vec

        n = size(coeffs)
        IF (n<=0) THEN
          poly_rr = 0.0_sp
        ELSE IF (n<npar_poly) THEN
          poly_rr = coeffs(n)
          DO i = n - 1, 1, -1
            poly_rr = x*poly_rr + coeffs(i)
          END DO
        ELSE
          ALLOCATE (vec(n+1),STAT=ierr)
          IF (ierr/=0) CALL nrerror( &
            'poly_rr: problem in attempt to allocate memory')
          pow = x
          vec(1:n) = coeffs
          DO
            vec(n+1) = 0.0_sp
            nn = ishft(n+1,-1)
            vec(1:nn) = vec(1:n:2) + pow*vec(2:n+1:2)
            IF (nn==1) EXIT
            pow = pow*pow
            n = nn
          END DO
          poly_rr = vec(1)
          DEALLOCATE (vec,STAT=ierr)
          IF (ierr/=0) CALL nrerror( &
            'poly_rr: problem in attempt to deallocate memory')
        END IF
      END FUNCTION poly_rr
!BL
      FUNCTION poly_rc(x,coeffs)
    COMPLEX(spc), INTENT(IN)                 :: x
    REAL(sp), DIMENSION(:), INTENT(IN)       :: coeffs
    COMPLEX(spc)                             :: poly_rc

    COMPLEX(spc)                             :: pow
    COMPLEX(spc), ALLOCATABLE, DIMENSION(:)  :: vec
    INTEGER(i4b)                             :: i, ierr, n, nn

        n = size(coeffs)
        IF (n<=0) THEN
          poly_rc = 0.0_sp
        ELSE IF (n<npar_poly) THEN
          poly_rc = coeffs(n)
          DO i = n - 1, 1, -1
            poly_rc = x*poly_rc + coeffs(i)
          END DO
        ELSE
          ALLOCATE (vec(n+1),STAT=ierr)
          IF (ierr/=0) CALL nrerror( &
            'poly_rc: problem in attempt to allocate memory')
          pow = x
          vec(1:n) = coeffs
          DO
            vec(n+1) = 0.0_sp
            nn = ishft(n+1,-1)
            vec(1:nn) = vec(1:n:2) + pow*vec(2:n+1:2)
            IF (nn==1) EXIT
            pow = pow*pow
            n = nn
          END DO
          poly_rc = vec(1)
          DEALLOCATE (vec,STAT=ierr)
          IF (ierr/=0) CALL nrerror( &
            'poly_rc: problem in attempt to deallocate memory')
        END IF
      END FUNCTION poly_rc
!BL
      FUNCTION poly_cc(x,coeffs)
    COMPLEX(spc), INTENT(IN)                 :: x
    COMPLEX(spc), DIMENSION(:), INTENT(IN)   :: coeffs
    COMPLEX(spc)                             :: poly_cc

    COMPLEX(spc)                             :: pow
    COMPLEX(spc), ALLOCATABLE, DIMENSION(:)  :: vec
    INTEGER(i4b)                             :: i, ierr, n, nn

        n = size(coeffs)
        IF (n<=0) THEN
          poly_cc = 0.0_sp
        ELSE IF (n<npar_poly) THEN
          poly_cc = coeffs(n)
          DO i = n - 1, 1, -1
            poly_cc = x*poly_cc + coeffs(i)
          END DO
        ELSE
          ALLOCATE (vec(n+1),STAT=ierr)
          IF (ierr/=0) CALL nrerror( &
            'poly_cc: problem in attempt to allocate memory')
          pow = x
          vec(1:n) = coeffs
          DO
            vec(n+1) = 0.0_sp
            nn = ishft(n+1,-1)
            vec(1:nn) = vec(1:n:2) + pow*vec(2:n+1:2)
            IF (nn==1) EXIT
            pow = pow*pow
            n = nn
          END DO
          poly_cc = vec(1)
          DEALLOCATE (vec,STAT=ierr)
          IF (ierr/=0) CALL nrerror( &
            'poly_cc: problem in attempt to deallocate memory')
        END IF
      END FUNCTION poly_cc
!BL
      FUNCTION poly_rrv(x,coeffs)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: x, coeffs
    REAL(sp), DIMENSION(size(x))             :: poly_rrv

    INTEGER(i4b)                             :: i, m, n

        m = size(coeffs)
        n = size(x)
        IF (m<=0) THEN
          poly_rrv = 0.0_sp
        ELSE IF (m<n .OR. m<npar_poly) THEN
          poly_rrv = coeffs(m)
          DO i = m - 1, 1, -1
            poly_rrv = x*poly_rrv + coeffs(i)
          END DO
        ELSE
          DO i = 1, n
            poly_rrv(i) = poly_rr(x(i),coeffs)
          END DO
        END IF
      END FUNCTION poly_rrv
!BL
      FUNCTION poly_msk_rrv(x,coeffs,mask)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: x, coeffs
    LOGICAL(lgt), DIMENSION(:), INTENT(IN)   :: mask
    REAL(sp), DIMENSION(size(x))             :: poly_msk_rrv

        poly_msk_rrv = unpack(poly_rrv(pack(x,mask),coeffs),mask,0.0_sp)
      END FUNCTION poly_msk_rrv
!BL
!BL
      RECURSIVE FUNCTION poly_term_rr(a,b) RESULT (u)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: a
    REAL(sp), INTENT(IN)                     :: b
    REAL(sp), DIMENSION(size(a))             :: u

    INTEGER(i4b)                             :: j, n

        n = size(a)
        IF (n<=0) RETURN
        u(1) = a(1)
        IF (n<npar_polyterm) THEN
          DO j = 2, n
            u(j) = a(j) + b*u(j-1)
          END DO
        ELSE
          u(2:n:2) = poly_term_rr(a(2:n:2)+a(1:n-1:2)*b,b*b)
          u(3:n:2) = a(3:n:2) + b*u(2:n-1:2)
        END IF
      END FUNCTION poly_term_rr
!BL
      RECURSIVE FUNCTION poly_term_cc(a,b) RESULT (u)
    COMPLEX(spc), DIMENSION(:), INTENT(IN)   :: a
    COMPLEX(spc), INTENT(IN)                 :: b
    COMPLEX(spc), DIMENSION(size(a))         :: u

    INTEGER(i4b)                             :: j, n

        n = size(a)
        IF (n<=0) RETURN
        u(1) = a(1)
        IF (n<npar_polyterm) THEN
          DO j = 2, n
            u(j) = a(j) + b*u(j-1)
          END DO
        ELSE
          u(2:n:2) = poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b)
          u(3:n:2) = a(3:n:2) + b*u(2:n-1:2)
        END IF
      END FUNCTION poly_term_cc
!BL
!BL
      FUNCTION zroots_unity(n,nn)
    INTEGER(i4b), INTENT(IN)                 :: n, nn
    COMPLEX(spc), DIMENSION(nn)              :: zroots_unity

    INTEGER(i4b)                             :: k
    REAL(sp)                                 :: theta

        zroots_unity(1) = 1.0
        theta = twopi/n
        k = 1
        DO
          IF (k>=nn) EXIT
          zroots_unity(k+1) = cmplx(cos(k*theta),sin(k*theta),spc)
          zroots_unity(k+2:min(2*k,nn)) = zroots_unity(k+1)* &
            zroots_unity(2:min(k,nn-k))
          k = 2*k
        END DO
      END FUNCTION zroots_unity
!BL
      FUNCTION outerprod_r(a,b)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: a, b
    REAL(sp), DIMENSION(size(a), size(b))    :: outerprod_r

        outerprod_r = spread(a,dim=2,ncopies=size(b))* &
          spread(b,dim=1,ncopies=size(a))
      END FUNCTION outerprod_r
!BL
      FUNCTION outerdiv(a,b)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: a, b
    REAL(sp), DIMENSION(size(a), size(b))    :: outerdiv

        outerdiv = spread(a,dim=2,ncopies=size(b))/spread(b,dim=1,ncopies=size &
          (a))
      END FUNCTION outerdiv
!BL
      FUNCTION outersum(a,b)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: a, b
    REAL(sp), DIMENSION(size(a), size(b))    :: outersum

        outersum = spread(a,dim=2,ncopies=size(b)) + &
          spread(b,dim=1,ncopies=size(a))
      END FUNCTION outersum
!BL
      FUNCTION outerdiff_r(a,b)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: a, b
    REAL(sp), DIMENSION(size(a), size(b))    :: outerdiff_r

        outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
          spread(b,dim=1,ncopies=size(a))
      END FUNCTION outerdiff_r
!BL
      FUNCTION outerdiff_i(a,b)
    INTEGER(i4b), DIMENSION(:), INTENT(IN)   :: a, b
    INTEGER(i4b), &
      DIMENSION(size(a), size(b))            :: outerdiff_i

        outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
          spread(b,dim=1,ncopies=size(a))
      END FUNCTION outerdiff_i
!BL
      FUNCTION outerand(a,b)
    LOGICAL(lgt), DIMENSION(:), INTENT(IN)   :: a, b
    LOGICAL(lgt), &
      DIMENSION(size(a), size(b))            :: outerand

        outerand = spread(a,dim=2,ncopies=size(b)) .AND. &
          spread(b,dim=1,ncopies=size(a))
      END FUNCTION outerand
!BL
      SUBROUTINE scatter_add_r(dest,source,dest_index)
    REAL(sp), DIMENSION(:), INTENT(OUT)      :: dest
    REAL(sp), DIMENSION(:), INTENT(IN)       :: source
    INTEGER(i4b), DIMENSION(:), INTENT(IN)   :: dest_index

    INTEGER(i4b)                             :: i, j, m, n

        n = assert_eq2(size(source),size(dest_index),'scatter_add_r')
        m = size(dest)
        DO j = 1, n
          i = dest_index(j)
          IF (i>0 .AND. i<=m) dest(i) = dest(i) + source(j)
        END DO
      END SUBROUTINE scatter_add_r
      SUBROUTINE scatter_max_r(dest,source,dest_index)
    REAL(sp), DIMENSION(:), INTENT(OUT)      :: dest
    REAL(sp), DIMENSION(:), INTENT(IN)       :: source
    INTEGER(i4b), DIMENSION(:), INTENT(IN)   :: dest_index

    INTEGER(i4b)                             :: i, j, m, n

        n = assert_eq2(size(source),size(dest_index),'scatter_max_r')
        m = size(dest)
        DO j = 1, n
          i = dest_index(j)
          IF (i>0 .AND. i<=m) dest(i) = max(dest(i),source(j))
        END DO
      END SUBROUTINE scatter_max_r
!BL
      SUBROUTINE diagadd_rv(mat,diag)
    REAL(sp), DIMENSION(:, :), INTENT(INOUT) :: mat
    REAL(sp), DIMENSION(:), INTENT(IN)       :: diag

    INTEGER(i4b)                             :: j, n

        n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
        DO j = 1, n
          mat(j,j) = mat(j,j) + diag(j)
        END DO
      END SUBROUTINE diagadd_rv
!BL
      SUBROUTINE diagadd_r(mat,diag)
    REAL(sp), DIMENSION(:, :), INTENT(INOUT) :: mat
    REAL(sp), INTENT(IN)                     :: diag

    INTEGER(i4b)                             :: j, n

        n = min(size(mat,1),size(mat,2))
        DO j = 1, n
          mat(j,j) = mat(j,j) + diag
        END DO
      END SUBROUTINE diagadd_r
!BL
      SUBROUTINE diagmult_rv(mat,diag)
    REAL(sp), DIMENSION(:, :), INTENT(INOUT) :: mat
    REAL(sp), DIMENSION(:), INTENT(IN)       :: diag

    INTEGER(i4b)                             :: j, n

        n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
        DO j = 1, n
          mat(j,j) = mat(j,j)*diag(j)
        END DO
      END SUBROUTINE diagmult_rv
!BL
      SUBROUTINE diagmult_r(mat,diag)
    REAL(sp), DIMENSION(:, :), INTENT(INOUT) :: mat
    REAL(sp), INTENT(IN)                     :: diag

    INTEGER(i4b)                             :: j, n

        n = min(size(mat,1),size(mat,2))
        DO j = 1, n
          mat(j,j) = mat(j,j)*diag
        END DO
      END SUBROUTINE diagmult_r
!BL
      FUNCTION get_diag_rv(mat)
    REAL(sp), DIMENSION(:, :), INTENT(IN)    :: mat
    REAL(sp), DIMENSION(size(mat, 1))        :: get_diag_rv

    INTEGER(i4b)                             :: j

        j = assert_eq2(size(mat,1),size(mat,2),'get_diag_rv')
        DO j = 1, size(mat,1)
          get_diag_rv(j) = mat(j,j)
        END DO
      END FUNCTION get_diag_rv
!BL
      SUBROUTINE put_diag_rv(diagv,mat)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: diagv
    REAL(sp), DIMENSION(:, :), INTENT(INOUT) :: mat

    INTEGER(i4b)                             :: j, n

        n = assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
        DO j = 1, n
          mat(j,j) = diagv(j)
        END DO
      END SUBROUTINE put_diag_rv
!BL
      SUBROUTINE put_diag_r(scal,mat)
    REAL(sp), INTENT(IN)                     :: scal
    REAL(sp), DIMENSION(:, :), INTENT(INOUT) :: mat

    INTEGER(i4b)                             :: j, n

        n = min(size(mat,1),size(mat,2))
        DO j = 1, n
          mat(j,j) = scal
        END DO
      END SUBROUTINE put_diag_r
!BL
      SUBROUTINE unit_matrix(mat)
    REAL(sp), DIMENSION(:, :), INTENT(OUT)   :: mat

    INTEGER(i4b)                             :: i, n

        n = min(size(mat,1),size(mat,2))
        mat(:,:) = 0.0_sp
        DO i = 1, n
          mat(i,i) = 1.0_sp
        END DO
      END SUBROUTINE unit_matrix
!BL
      FUNCTION upper_triangle(j,k,extra)
    INTEGER(i4b), INTENT(IN)                 :: j, k
    INTEGER(i4b), INTENT(IN), OPTIONAL       :: extra
    LOGICAL(lgt), DIMENSION(j, k)            :: upper_triangle

    INTEGER(i4b)                             :: n

        n = 0
        IF (present(extra)) n = extra
        upper_triangle = (outerdiff(arth_i(1,1,j),arth_i(1,1,k))<n)
      END FUNCTION upper_triangle
!BL
      FUNCTION lower_triangle(j,k,extra)
    INTEGER(i4b), INTENT(IN)                 :: j, k
    INTEGER(i4b), INTENT(IN), OPTIONAL       :: extra
    LOGICAL(lgt), DIMENSION(j, k)            :: lower_triangle

    INTEGER(i4b)                             :: n

        n = 0
        IF (present(extra)) n = extra
        lower_triangle = (outerdiff(arth_i(1,1,j),arth_i(1,1,k))>-n)
      END FUNCTION lower_triangle
!BL
      FUNCTION vabs(v)
    REAL(sp), DIMENSION(:), INTENT(IN)       :: v
    REAL(sp)                                 :: vabs

        vabs = sqrt(dot_product(v,v))
      END FUNCTION vabs
!BL
    END MODULE nrutil
