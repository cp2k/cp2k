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
        USE functions, ONLY: evalf=>f
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

!$OMP PARALLEL DO &
!$OMP DEFAULT(NONE) &
!$OMP PRIVATE(i,x1,x2,res) &
!$OMP SHARED(npd,a1,b1,a2,b2,pd1,pd2,fpd,nderiv)
        DO i = 1,npd
           x1=((b1 - a1) * pd1(i) + (b1 + a1)) / 2
           x2=((b2 - a2) * pd2(i) + (b2 + a2)) / 2
           CALL evalf(res,nderiv,x1,x2)
           FPD(I,:) = RES(0:NDERIV)
        ENDDO

!$OMP PARALLEL DO &
!$OMP DEFAULT(NONE) &
!$OMP PRIVATE(i,t1,t2) &
!$OMP SHARED(npd,deg,pd1,pd2,wpd,fpd,c0,esterr)
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
         USE functions, ONLY: functionid,should_output
         IMPLICIT NONE
         INTEGER :: deg,nderiv,i,j,digits,deg_ref,k,l
         INTEGER*2 :: precision

         REAL(KIND=dp) :: esterr,A1,B1,A2,B2
         REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: c0

         TYPE(mpfr_type) ::  A1_mp,B1_mp,A2_mp,B2_mp,esterr_mp
         TYPE(mpfr_type), DIMENSION(:,:,:), ALLOCATABLE :: c0_mp

         ! keeps a history of evaluated function values. 
         ! useful to investigate the error later on
         should_output=31
         IF (should_output>0) THEN
            OPEN(should_output,FILE="function_vals.dat")
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
!$OMP PARALLEL DEFAULT(NONE) SHARED(precision)
         CALL mpfr_set_default_precision(precision)
!$OMP END PARALLEL

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
