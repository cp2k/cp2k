!
!     a simple wrapper program to
!     obtain the C0 coefficients for an NDERIV+1 valued function in the
!     rectangular domain A1<=X1<=B1, A2<=X2<=B2
!
!     input from 'drvrec_mpfr.in'
!     output to 'drvrec_mpfr.out'
!

      SUBROUTINE drvrec(a1,b1,a2,b2,deg,nderiv,c0,esterr_out)
        USE padua2_mpfr
        IMPLICIT NONE

        INTEGER deg,npd,nderiV
        TYPE(mpfr_type) pd1((deg + 1) * (deg + 2) / 2)
        TYPE(mpfr_type) pd2((deg + 1) * (deg + 2) / 2)
        TYPE(mpfr_type) wpd((deg + 1) * (deg + 2) / 2)
        TYPE(mpfr_type) fpd((deg + 1) * (deg + 2) / 2,0:nderiv)
        TYPE(mpfr_type) t1(0:deg,(deg + 1) * (deg + 2) / 2)
        TYPE(mpfr_type) t2(0:deg,(deg + 1) * (deg + 2) / 2)
        TYPE(mpfr_type) c0(0:deg,0:deg,0:nderiv)
        TYPE(mpfr_type) res(0:nderiv),esterr(0:nderiv)
        TYPE(mpfr_type) a1,b1,a2,b2,esterr_out,x1,x2,zero
        INTEGER i,ntg1

        CALL pdpts(deg,pd1,pd2,wpd,npd)
        
        DO i = 1,npd
           x1=((b1 - a1) * pd1(i) + (b1 + a1)) / 2
           x2=((b2 - a2) * pd2(i) + (b2 + a2)) / 2
           CALL f(res,nderiv,x1,x2)
           FPD(I,:) = RES(0:NDERIV)
        ENDDO
        DO i=0,nderiv
           CALL padua2(deg,npd,pd1,pd2,wpd,fpd(:,i),t1,t2,c0(:,:,i),esterr(i))
        ENDDO

        esterr_out=0
        zero=0 
        DO i=0,nderiv
           IF (mpfr_cmp(esterr(i),zero)>0) THEN
              IF (mpfr_cmp(esterr(i),esterr_out)>0) THEN
                   esterr_out=esterr(i)
              ENDIF
           ELSE
              IF (mpfr_cmp(-esterr(i),esterr_out)>0) esterr_out=-esterr(i)
           ENDIF
        ENDDO

      END SUBROUTINE drvrec

      ! evaluate the function, + nderiv other functions in the point X1,X2
      ! write additionally to a file 'history.res' for later checking the quality of the interpolation
      SUBROUTINE f(res,nderiv,x1,x2)
        USE mpfr
        USE mpfr_ops
        USE mpfr_cutoff_gamma
        IMPLICIT NONE
        INTEGER         :: nderiv
        TYPE(mpfr_type) :: x1,x2
        TYPE(mpfr_type) :: res(0:nderiv)

        TYPE(mpfr_type) :: T,r,upper,lower,zero
        TYPE(mpfr_type) :: dummy(0:21)
        INTEGER          :: i

        INTEGER :: functionid,should_output
        COMMON/functiontype/functionid,should_output

        zero=.CONVERT."0"

        ! functionid in this case just selects a different
        ! domain of the same function
        ! 1 = farfield  (R>11)
        ! 2 = nearfield (R<11)
        SELECT CASE(functionid)
        CASE (1)
          ! if R=+Infinity the result is zero
          IF (mpfr_cmp(x2,zero)<=0) THEN
             res=.CONVERT."0"
             return
          ENDIF
          r=11/x2
          upper=r**2 + 11*r + 50
          lower=r**2 - 11*r 
        CASE (2)
          R=x2*11
          upper=r**2 + 11*r + 50
          lower=.CONVERT."0"
        CASE DEFAULT
          STOP "Function ID not implemented"
        END SELECT 

        t=lower+x1*(upper-lower)

        !t is zero, return the limiting expansion
        IF (mpfr_cmp(t,zero)<=0) THEN
           CALL cutoff_gamma_T0(21,R,dummy)
        ELSE
           CALL cutoff_gamma(21,t,r,dummy)
        END IF
        res(0:nderiv)=dummy(0:nderiv)

        IF (should_output>0) THEN
           WRITE(should_output,*) REAL(R),REAL(T),REAL(res)
        ENDIF

      END SUBROUTINE f

      !
      ! this program is a simple interface to obtain the C0 coefficients
      ! with high accuracy (i.e. all digits & correctly rounded),
      ! for a function programmed in subroutine f
      ! this version uses double precision (~16 digits)
      ! input/output but does all calculations 
      ! internally using a user specified number of digits
      !
      PROGRAM drvrec_padua
         USE mpfr
         USE mpfr_ops
         USE padua2_mpfr, ONLY: pd2val_mp=>pd2val
         IMPLICIT NONE
         INTEGER :: deg,nderiv,i,j,digits,deg_ref,k,l
         INTEGER*2 :: precision

         REAL(KIND=dp) :: esterr,A1,B1,A2,B2
         REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: c0

         TYPE(mpfr_type) ::  A1_mp,B1_mp,A2_mp,B2_mp,esterr_mp
         TYPE(mpfr_type), DIMENSION(:,:,:), ALLOCATABLE :: c0_mp

         integer :: functionid,should_output
         COMMON/functiontype/functionid,should_output

         ! keeps a history of evaluated function values. 
         ! useful to investigate the error later on
         should_output=31
         IF (should_output>0) THEN
           OPEN(should_output,FILE="history.dat",ACCESS="APPEND") 
         ENDIF

         ! input definition
         OPEN(17,FILE="drvrec_mpfr.in")
         READ(17,*) digits
         READ(17,*) deg
         READ(17,*) nderiv
         READ(17,*) A1,B1
         READ(17,*) A2,B2
         READ(17,*) functionid
         CLOSE(17)

         precision = log(10.0)/log(2.0)*digits
         CALL mpfr_set_default_precision(precision)

         A1_mp=A1
         B1_mp=B1
         A2_mp=A2
         B2_mp=B2

         ALLOCATE(c0_mp(0:deg,0:deg,0:nderiv))
         CALL drvrec(a1_mp,b1_mp,a2_mp,b2_mp,deg,nderiv,c0_mp,esterr_mp)

         ! convert ot double precision and output
         ALLOCATE(c0(0:deg,0:deg,0:nderiv))
         c0=REAL(c0_mp)
         esterr=REAL(esterr_mp)

         OPEN(17,FILE="drvrec_mpfr.out")
         write(17,*) esterr
         write(17,*) c0
         CLOSE(17)

         IF (should_output>0) THEN
           CLOSE(should_output)
         ENDIF

      END PROGRAM drvrec_padua
