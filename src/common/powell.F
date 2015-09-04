!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!******************************************************************************
MODULE powell
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'powell'

  TYPE opt_state_type
    INTEGER            :: state
    INTEGER            :: nvar
    INTEGER            :: iprint
    INTEGER            :: unit
    INTEGER            :: maxfun
    REAL(dp)           :: rhobeg, rhoend
    REAL(dp), DIMENSION(:), POINTER  :: w
    REAL(dp), DIMENSION(:), POINTER  :: xopt
    ! local variables
    INTEGER            :: np, nh, nptm, nftest, idz, itest, nf, nfm, nfmm, &
                          nfsav, knew, kopt, ksave, ktemp
    REAL(dp)           :: rhosq, recip, reciq, fbeg, fopt, diffa, xoptsq, &
                          rho, delta, dsq, dnorm, ratio, temp, tempq, beta, &
                          dx, vquad, diff, diffc, diffb, fsave, detrat, hdiag, &
                          distsq, gisq, gqsq, f, bstep, alpha, dstep
  END TYPE opt_state_type

  PRIVATE
  PUBLIC      :: powell_optimize, opt_state_type

!******************************************************************************

CONTAINS

!******************************************************************************
! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param x ...
!> \param optstate ...
! *****************************************************************************
  SUBROUTINE powell_optimize (n,x,optstate)
    INTEGER                                  :: n
    REAL(dp), DIMENSION(*)                   :: x
    TYPE(opt_state_type)                     :: optstate

    CHARACTER(len=*), PARAMETER :: routineN = 'powell_optimize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npt

    CALL timeset(routineN,handle)

    SELECT CASE (optstate%state)
    CASE (0)
       npt=2*n+1
       ALLOCATE(optstate%w((npt+13)*(npt+n)+3*n*(n+3)/2))
       ALLOCATE(optstate%xopt(n))
       ! Initialize w
       optstate%w = 0.0_dp
       optstate%state = 1
       CALL newuoa (n, x, optstate)
    CASE (1,2)
       CALL newuoa (n, x, optstate)
    CASE (3)
       IF(optstate%unit > 0) THEN
         WRITE(optstate%unit,*) "POWELL| Exceeding maximum number of steps"
       ENDIF
       optstate%state = -1
    CASE (4)
       IF(optstate%unit > 0) THEN
         WRITE(optstate%unit,*) "POWELL| Error in trust region"
       ENDIF
       optstate%state = -1
    CASE (5)
       IF(optstate%unit > 0) THEN
         WRITE(optstate%unit,*) "POWELL| N out of range"
       ENDIF
       optstate%state = -1
    CASE (6,7)
       optstate%state = -1
    CASE (8)
       x(1:n) = optstate%xopt(1:n)
       DEALLOCATE(optstate%w)
       DEALLOCATE(optstate%xopt)
       optstate%state = -1
    CASE DEFAULT
       STOP
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE powell_optimize
!******************************************************************************
! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param x ...
!> \param optstate ...
! *****************************************************************************
  SUBROUTINE newuoa (n,x,optstate)

    INTEGER                                  :: n
    REAL(dp), DIMENSION(*)                   :: x
    TYPE(opt_state_type)                     :: optstate

    INTEGER                                  :: ibmat, id, ifv, igq, ihq, &
                                                ipq, ivl, iw, ixb, ixn, ixo, &
                                                ixp, izmat, maxfun, ndim, np, &
                                                npt, nptm
    REAL(dp)                                 :: rhobeg, rhoend

    maxfun = optstate%maxfun
    rhobeg = optstate%rhobeg
    rhoend = optstate%rhoend

    !
    !     This subroutine seeks the least value of a function of many variab
    !     by a trust region method that forms quadratic models by interpolat
    !     There can be some freedom in the interpolation conditions, which i
    !     taken up by minimizing the Frobenius norm of the change to the sec
    !     derivative of the quadratic model, beginning with a zero matrix. T
    !     arguments of the subroutine are as follows.
    !
    !     N must be set to the number of variables and must be at least two.
    !     NPT is the number of interpolation conditions. Its value must be i
    !       interval [N+2,(N+1)(N+2)/2].
    !     Initial values of the variables must be set in X(1),X(2),...,X(N).
    !       will be changed to the values that give the least calculated F.
    !     RHOBEG and RHOEND must be set to the initial and final values of a
    !       region radius, so both must be positive with RHOEND<=RHOBEG. Typ
    !       RHOBEG should be about one tenth of the greatest expected change
    !       variable, and RHOEND should indicate the accuracy that is requir
    !       the final values of the variables.
    !     The value of IPRINT should be set to 0, 1, 2 or 3, which controls
    !       amount of printing. Specifically, there is no output if IPRINT=0
    !       there is output only at the return if IPRINT=1. Otherwise, each
    !       value of RHO is printed, with the best vector of variables so fa
    !       the corresponding value of the objective function. Further, each
    !       value of F with its variables are output if IPRINT=3.
    !     MAXFUN must be set to an upper bound on the number of calls of CAL
    !     The array W will be used for working space. Its length must be at
    !     (NPT+13)*(NPT+N)+3*N*(N+3)/2.
    !
    !     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must se
    !     the value of the objective function for the variables X(1),X(2),..
    !
    !     Partition the working space array, so that different parts of it c
    !     treated separately by the subroutine that performs the main calcul
    !
    np=n+1
    npt=2*n+1
    nptm=npt-np
    IF (npt < n+2 .OR. npt > ((n+2)*np)/2) THEN
       optstate%state = 5
       RETURN
    END IF
    ndim=npt+n
    ixb=1
    ixo=ixb+n
    ixn=ixo+n
    ixp=ixn+n
    ifv=ixp+n*npt
    igq=ifv+npt
    ihq=igq+n
    ipq=ihq+(n*np)/2
    ibmat=ipq+npt
    izmat=ibmat+ndim*n
    id=izmat+npt*nptm
    ivl=id+n
    iw=ivl+ndim
    !
    !     The above settings provide a partition of W for subroutine NEWUOB.
    !     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements
    !     W plus the space that is needed by the last array of NEWUOB.
    !
    CALL newuob (n,npt,x,rhobeg,rhoend,maxfun,optstate%w(ixb:),optstate%w(ixo:),&
         optstate%w(ixn:),optstate%w(ixp:),optstate%w(ifv:),optstate%w(igq:),optstate%w(ihq:),&
         optstate%w(ipq:),optstate%w(ibmat:),optstate%w(izmat:),ndim,optstate%w(id:),&
         optstate%w(ivl:),optstate%w(iw:),optstate)

    optstate%xopt(1:n) = optstate%w(ixb:ixb+n-1) + optstate%w(ixo:ixo+n-1)


  END SUBROUTINE newuoa

!******************************************************************************
! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param npt ...
!> \param x ...
!> \param rhobeg ...
!> \param rhoend ...
!> \param maxfun ...
!> \param xbase ...
!> \param xopt ...
!> \param xnew ...
!> \param xpt ...
!> \param fval ...
!> \param gq ...
!> \param hq ...
!> \param pq ...
!> \param bmat ...
!> \param zmat ...
!> \param ndim ...
!> \param d ...
!> \param vlag ...
!> \param w ...
!> \param opt ...
! *****************************************************************************
  SUBROUTINE newuob (n,npt,x,rhobeg,rhoend,maxfun,xbase,&
       xopt,xnew,xpt,fval,gq,hq,pq,bmat,zmat,ndim,d,vlag,w,opt)

    INTEGER, INTENT(inout)                   :: n, npt
    REAL(dp), DIMENSION(1:n), INTENT(inout)  :: x
    REAL(dp), INTENT(inout)                  :: rhobeg, rhoend
    INTEGER, INTENT(inout)                   :: maxfun
    REAL(dp), DIMENSION(*), INTENT(inout)    :: xbase, xopt, xnew
    REAL(dp), DIMENSION(npt, *), &
      INTENT(inout)                          :: xpt
    REAL(dp), DIMENSION(*), INTENT(inout)    :: fval, gq, hq, pq
    INTEGER, INTENT(inout)                   :: ndim
    REAL(dp), DIMENSION(npt, *), &
      INTENT(inout)                          :: zmat
    REAL(dp), DIMENSION(ndim, *), &
      INTENT(inout)                          :: bmat
    REAL(dp), DIMENSION(*), INTENT(inout)    :: d, vlag, w
    TYPE(opt_state_type)                     :: opt

    INTEGER                                  :: i, idz, ih, ip, ipt, itemp, &
                                                itest, j, jp, jpt, k, knew, &
                                                kopt, ksave, ktemp, nf, nfm, &
                                                nfmm, nfsav, nftest, nh, np, &
                                                nptm
    REAL(dp) :: alpha, beta, bstep, bsum, crvmin, delta, detrat, diff, diffa, &
      diffb, diffc, distsq, dnorm, dsq, dstep, dx, f, fbeg, fopt, fsave, &
      gisq, gqsq, half, hdiag, one, ratio, recip, reciq, rho, rhosq, sum, &
      suma, sumb, sumz, temp, tempq, tenth, vquad, xipt, xjpt, xoptsq, zero

!
!     The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are ide
!       to the corresponding arguments in SUBROUTINE NEWUOA.
!     XBASE will hold a shift of origin that should reduce the contribut
!       from rounding errors to values of the model and Lagrange functio
!     XOPT will be set to the displacement from XBASE of the vector of
!       variables that provides the least calculated F so far.
!     XNEW will be set to the displacement from XBASE of the vector of
!       variables for the current calculation of F.
!     XPT will contain the interpolation point coordinates relative to X
!     FVAL will hold the values of F at the interpolation points.
!     GQ will hold the gradient of the quadratic model at XBASE.
!     HQ will hold the explicit second derivatives of the quadratic mode
!     PQ will contain the parameters of the implicit second derivatives
!       the quadratic model.
!     BMAT will hold the last N columns of H.
!     ZMAT will hold the factorization of the leading NPT by NPT submatr
!       H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, wh
!       the elements of DZ are plus or minus one, as specified by IDZ.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     D is reserved for trial steps from XOPT.
!     VLAG will contain the values of the Lagrange functions at a new po
!       They are part of a product that requires VLAG to be of length ND
!     The array W will be used for working space. Its length must be at
!       10*NDIM = 10*(NPT+N).

    IF ( opt%state == 1 ) THEN
       ! initialize all variable that will be stored
       np           = 0
       nh           = 0
       nptm         = 0
       nftest       = 0
       idz          = 0
       itest        = 0
       nf           = 0
       nfm          = 0
       nfmm         = 0
       nfsav        = 0
       knew         = 0
       kopt         = 0
       ksave        = 0
       ktemp        = 0
       rhosq        = 0._dp
       recip        = 0._dp
       reciq        = 0._dp
       fbeg         = 0._dp
       fopt         = 0._dp
       diffa        = 0._dp
       xoptsq       = 0._dp
       rho          = 0._dp
       delta        = 0._dp
       dsq          = 0._dp
       dnorm        = 0._dp
       ratio        = 0._dp
       temp         = 0._dp
       tempq        = 0._dp
       beta         = 0._dp
       dx           = 0._dp
       vquad        = 0._dp
       diff         = 0._dp
       diffc        = 0._dp
       diffb        = 0._dp
       fsave        = 0._dp
       detrat       = 0._dp
       hdiag        = 0._dp
       distsq       = 0._dp
       gisq         = 0._dp
       gqsq         = 0._dp
       f            = 0._dp
       bstep        = 0._dp
       alpha        = 0._dp
       dstep        = 0._dp
       !
    END IF

    ipt          = 0
    jpt          = 0
    xipt         = 0._dp
    xjpt         = 0._dp

    half=0.5_dp
    one=1.0_dp
    tenth=0.1_dp
    zero=0.0_dp
    np=n+1
    nh=(n*np)/2
    nptm=npt-np
    nftest=MAX(maxfun,1)

    IF ( opt%state == 2 ) GOTO 1000
    !
    !     Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
    !
    DO j=1,n
       xbase(j)=x(j)
       DO k=1,npt
          xpt(k,j)=zero
       END DO
       DO i=1,ndim
          bmat(i,j)=zero
       END DO
    END DO
    DO ih=1,nh
       hq(ih)=zero
    END DO
    DO k=1,npt
       pq(k)=zero
       DO j=1,nptm
          zmat(k,j)=zero
       END DO
    END DO
    !
    !     Begin the initialization procedure. NF becomes one more than the n
    !     of function values so far. The coordinates of the displacement of
    !     next initial interpolation point from XBASE are set in XPT(NF,.).
    !
    rhosq=rhobeg*rhobeg
    recip=one/rhosq
    reciq=SQRT(half)/rhosq
    nf=0
50  nfm=nf
    nfmm=nf-n
    nf=nf+1
    IF (nfm <= 2*n) THEN
       IF (nfm >= 1 .AND. nfm <= N) THEN
          xpt(nf,nfm)=rhobeg
       ELSE IF (nfm > n) THEN
          xpt(nf,nfmm)=-rhobeg
       END IF
    ELSE
       itemp=(nfmm-1)/n
       jpt=nfm-itemp*n-n
       ipt=jpt+itemp
       IF (ipt > n) THEN
          itemp=jpt
          jpt=ipt-n
          ipt=itemp
       END IF
       xipt=rhobeg
       IF (fval(ipt+np) < fval(ipt+1)) xipt=-xipt
       XJPT=RHOBEG
       IF (fval(jpt+np) < fval(jpt+1)) xjpt=-xjpt
       xpt(nf,ipt)=xipt
       xpt(nf,jpt)=xjpt
    END IF
    !
    !     Calculate the next value of F, label 70 being reached immediately
    !     after this calculation. The least function value so far and its in
    !     are required.
    !
    DO j=1,n
       x(j)=xpt(nf,j)+xbase(j)
    END DO
    GOTO 310
70  fval(nf)=f
    IF (nf == 1) THEN
       fbeg=f
       fopt=f
       kopt=1
    ELSE IF (f < fopt) THEN
       fopt=f
       kopt=nf
    END IF
    !
    !     Set the nonzero initial elements of BMAT and the quadratic model i
    !     the cases when NF is at most 2*N+1.
    !
    IF (NFM <= 2*N) THEN
       IF (nfm >= 1 .AND. nfm <= n) THEN
          gq(nfm)=(f-fbeg)/rhobeg
          IF (npt < nf+n) THEN
             bmat(1,nfm)=-one/rhobeg
             bmat(nf,nfm)=one/rhobeg
             bmat(npt+nfm,nfm)=-half*rhosq
          END IF
       ELSE IF (nfm > n) THEN
          bmat(nf-n,nfmm)=half/rhobeg
          bmat(nf,nfmm)=-half/rhobeg
          zmat(1,nfmm)=-reciq-reciq
          zmat(nf-n,nfmm)=reciq
          zmat(nf,nfmm)=reciq
          ih=(nfmm*(nfmm+1))/2
          temp=(fbeg-f)/rhobeg
          hq(ih)=(gq(nfmm)-temp)/rhobeg
          gq(nfmm)=half*(gq(nfmm)+temp)
       END IF
       !
       !     Set the off-diagonal second derivatives of the Lagrange functions
       !     the initial quadratic model.
       !
    ELSE
       ih=(ipt*(ipt-1))/2+jpt
       IF (xipt < zero) ipt=ipt+n
       IF (xjpt < zero) jpt=jpt+n
       zmat(1,nfmm)=recip
       zmat(nf,nfmm)=recip
       zmat(ipt+1,nfmm)=-recip
       zmat(jpt+1,nfmm)=-recip
       hq(ih)=(fbeg-fval(ipt+1)-fval(jpt+1)+f)/(xipt*xjpt)
    END IF
    IF (nf < npt) GOTO 50
    !
    !     Begin the iterative procedure, because the initial model is comple
    !
    rho=rhobeg
    delta=rho
    idz=1
    diffa=zero
    diffb=zero
    itest=0
    xoptsq=zero
    DO i=1,n
       xopt(i)=xpt(kopt,i)
       xoptsq=xoptsq+xopt(i)**2
    END DO
90  nfsav=nf
    !
    !     Generate the next trust region step and test its length. Set KNEW
    !     to -1 if the purpose of the next F will be to improve the model.
    !
100 knew=0
    CALL trsapp (n,npt,xopt,xpt,gq,hq,pq,delta,d,w,w(np),w(np+n),w(np+2*n),crvmin)
    dsq=zero
    DO i=1,n
       dsq=dsq+d(i)**2
    END DO
    dnorm=MIN(delta,SQRT(dsq))
    IF (dnorm < half*rho) THEN
       knew=-1
       delta=tenth*delta
       ratio=-1.0_dp
       IF (delta <= 1.5_dp*rho) delta=rho
       IF (nf <= nfsav+2) GOTO 460
       temp=0.125_dp*crvmin*rho*rho
       IF (temp <= MAX(diffa,diffb,diffc)) GOTO 460
       GOTO 490
    END IF
    !
    !     Shift XBASE if XOPT may be too far from XBASE. First make the chan
    !     to BMAT that do not depend on ZMAT.
    !
120 IF (dsq <= 1.0e-3_dp*xoptsq) THEN
       tempq=0.25_dp*xoptsq
       DO k=1,npt
          sum=zero
          DO i=1,n
             sum=sum+xpt(k,i)*xopt(i)
          END DO
          temp=pq(k)*sum
          sum=sum-half*xoptsq
          w(npt+k)=sum
          DO i=1,n
             gq(i)=gq(i)+temp*xpt(k,i)
             xpt(k,i)=xpt(k,i)-half*xopt(i)
             vlag(i)=bmat(k,i)
             w(i)=sum*xpt(k,i)+tempq*xopt(i)
             ip=npt+i
             DO j=1,i
                bmat(ip,j)=bmat(ip,j)+vlag(i)*w(j)+w(i)*vlag(j)
             END DO
          END DO
       END DO
       !
       !     Then the revisions of BMAT that depend on ZMAT are calculated.
       !
       DO k=1,nptm
          sumz=zero
          DO i=1,npt
             sumz=sumz+zmat(i,k)
             w(i)=w(npt+i)*zmat(i,k)
          END DO
          DO j=1,n
             sum=tempq*sumz*xopt(j)
             DO i=1,npt
                sum=sum+w(i)*xpt(i,j)
                vlag(j)=sum
                IF (k < idz) sum=-sum
             END DO
             DO i=1,npt
                bmat(i,j)=bmat(i,j)+sum*zmat(i,k)
             END DO
          END DO
          DO i=1,n
             ip=i+npt
             temp=vlag(i)
             IF (k < idz) temp=-temp
             DO j=1,i
                bmat(ip,j)=bmat(ip,j)+temp*vlag(j)
             END DO
          END DO
       END DO
       !
       !     The following instructions complete the shift of XBASE, including
       !     the changes to the parameters of the quadratic model.
       !
       ih=0
       DO j=1,n
          w(j)=zero
          DO k=1,npt
             w(j)=w(j)+pq(k)*xpt(k,j)
             xpt(k,j)=xpt(k,j)-half*xopt(j)
          END DO
          DO i=1,j
             ih=ih+1
             IF (i < j) gq(j)=gq(j)+hq(ih)*xopt(i)
             gq(i)=gq(i)+hq(ih)*xopt(j)
             hq(ih)=hq(ih)+w(i)*xopt(j)+xopt(i)*w(j)
             bmat(npt+i,j)=bmat(npt+j,i)
          END DO
       END DO
       DO j=1,n
          xbase(j)=xbase(j)+xopt(j)
          xopt(j)=zero
       END DO
       xoptsq=zero
    END IF
    !
    !     Pick the model step if KNEW is positive. A different choice of D
    !     may be made later, if the choice of D by BIGLAG causes substantial
    !     cancellation in DENOM.
    !
    IF (knew > 0) THEN
       CALL biglag (n,npt,xopt,xpt,bmat,zmat,idz,ndim,knew,dstep,    &
                   d,alpha,vlag,vlag(npt+1),w,w(np),w(np+n))
    END IF
    !
    !     Calculate VLAG and BETA for the current choice of D. The first NPT
    !     components of W_check will be held in W.
    !
    DO k=1,npt
       suma=zero
       sumb=zero
       sum=zero
       DO j=1,n
          suma=suma+xpt(k,j)*d(j)
          sumb=sumb+xpt(k,j)*xopt(j)
          sum=sum+bmat(k,j)*d(j)
       END DO
       w(k)=suma*(half*suma+sumb)
       vlag(k)=sum
    END DO
    beta=zero
    DO k=1,nptm
       sum=zero
       DO i=1,npt
          sum=sum+zmat(i,k)*w(i)
       END DO
       IF (k < idz) THEN
          beta=beta+sum*sum
          sum=-sum
       ELSE
          beta=beta-sum*sum
       END IF
       DO i=1,npt
          vlag(i)=vlag(i)+sum*zmat(i,k)
       END DO
    END DO
    bsum=zero
    dx=zero
    DO j=1,n
       sum=zero
       DO i=1,npt
          sum=sum+w(i)*bmat(i,j)
       END DO
       bsum=bsum+sum*d(j)
       jp=npt+j
       DO k=1,n
          sum=sum+bmat(jp,k)*d(k)
       END DO
       vlag(jp)=sum
       bsum=bsum+sum*d(j)
       dx=dx+d(j)*xopt(j)
    END DO
    beta=dx*dx+dsq*(xoptsq+dx+dx+half*dsq)+beta-bsum
    vlag(kopt)=vlag(kopt)+one
    !
    !     If KNEW is positive and if the cancellation in DENOM is unacceptab
    !     then BIGDEN calculates an alternative model step, XNEW being used
    !     working space.
    !
    IF (knew > 0) THEN
       temp=one+alpha*beta/vlag(knew)**2
       IF (ABS(temp) <= 0.8_dp) THEN
          CALL bigden (n,npt,xopt,xpt,bmat,zmat,idz,ndim,kopt,      &
                        knew,d,w,vlag,beta,xnew,w(ndim+1),w(6*ndim+1))
       END IF
    END IF
    !
    !     Calculate the next value of the objective function.
    !
290 DO i=1,n
       xnew(i)=xopt(i)+d(i)
       x(i)=xbase(i)+xnew(i)
    END DO
    nf=nf+1
310 IF (nf > nftest) THEN
       !         return to many steps
       nf=nf-1
       opt%state = 3
       CALL get_state
       GOTO 530
    END IF

    CALL get_state

    opt%state = 2

    RETURN

1000 CONTINUE

    CALL set_state

    IF (nf <= npt) GOTO 70
    IF (knew == -1) THEN
       opt%state = 6
       CALL get_state
       GOTO 530
    END IF
    !
    !     Use the quadratic model to predict the change in F due to the step
    !     and set DIFF to the error of this prediction.
    !
    vquad=zero
    ih=0
    DO j=1,n
       vquad=vquad+d(j)*gq(j)
       DO i=1,j
          ih=ih+1
          temp=d(i)*xnew(j)+d(j)*xopt(i)
          IF (i == j) temp=half*temp
          vquad=vquad+temp*hq(ih)
       END DO
    END DO
    DO k=1,npt
       vquad=vquad+pq(k)*w(k)
    END DO
    diff=f-fopt-vquad
    diffc=diffb
    diffb=diffa
    diffa=ABS(diff)
    IF (dnorm > rho) nfsav=nf
    !
    !     Update FOPT and XOPT if the new F is the least value of the object
    !     function so far. The branch when KNEW is positive occurs if D is n
    !     a trust region step.
    !
    fsave=fopt
    IF (f < fopt) THEN
       fopt=f
       xoptsq=zero
       DO i=1,n
          xopt(i)=xnew(i)
          xoptsq=xoptsq+xopt(i)**2
       END DO
    END IF
    ksave=knew
    IF (knew > 0) GOTO 410
    !
    !     Pick the next value of DELTA after a trust region step.
    !
    IF (vquad >= zero) THEN
       ! Return because a trust region step has failed to reduce Q
       opt%state = 4
       CALL get_state
       GOTO 530
    END IF
    ratio=(f-fsave)/vquad
    IF (ratio <= tenth) THEN
       delta=half*dnorm
    ELSE IF (ratio <= 0.7_dp) THEN
       delta=MAX(half*delta,dnorm)
    ELSE
       delta=MAX(half*delta,dnorm+dnorm)
    END IF
    IF (delta <= 1.5_dp*rho) delta=rho
    !
    !     Set KNEW to the index of the next interpolation point to be delete
    !
    rhosq=MAX(tenth*delta,rho)**2
    ktemp=0
    detrat=zero
    IF (f >= fsave) THEN
       ktemp=kopt
       detrat=one
    END IF
    DO k=1,npt
       hdiag=zero
       DO j=1,nptm
          temp=one
          IF (j < idz) temp=-one
          hdiag=hdiag+temp*zmat(k,j)**2
       END DO
       temp=ABS(beta*hdiag+vlag(k)**2)
       distsq=zero
       DO j=1,n
          distsq=distsq+(xpt(k,j)-xopt(j))**2
       END DO
       IF (distsq > rhosq) temp=temp*(distsq/rhosq)**3
       IF (temp > detrat .AND. k /= ktemp) THEN
          detrat=temp
          knew=k
       END IF
    END DO
    IF (knew == 0) GOTO 460
    !
    !     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
    !     can be moved. Begin the updating of the quadratic model, starting
    !     with the explicit second derivative term.
    !
410 CALL update (n,npt,bmat,zmat,idz,ndim,vlag,beta,knew,w)
    fval(knew)=f
    ih=0
    DO i=1,n
       temp=pq(knew)*xpt(knew,i)
       DO j=1,i
          ih=ih+1
          hq(ih)=hq(ih)+temp*xpt(knew,j)
       END DO
    END DO
    pq(knew)=zero
    !
    !     Update the other second derivative parameters, and then the gradie
    !     vector of the model. Also include the new interpolation point.
    !
    DO j=1,nptm
       temp=diff*zmat(knew,j)
       IF (j < idz) temp=-temp
       DO k=1,npt
          pq(k)=pq(k)+temp*zmat(k,j)
       END DO
    END DO
    gqsq=zero
    DO i=1,n
       gq(i)=gq(i)+diff*bmat(knew,i)
       gqsq=gqsq+gq(i)**2
       xpt(knew,i)=xnew(i)
    END DO
    !
    !     If a trust region step makes a small change to the objective funct
    !     then calculate the gradient of the least Frobenius norm interpolan
    !     XBASE, and store it in W, using VLAG for a vector of right hand si
    !
    IF (ksave == 0 .AND. delta == rho) THEN
       IF (ABS(ratio) > 1.0e-2_dp) THEN
          itest=0
       ELSE
          DO k=1,npt
             vlag(k)=fval(k)-fval(kopt)
          END DO
          gisq=zero
          DO i=1,n
             sum=zero
             DO k=1,npt
                sum=sum+bmat(k,i)*vlag(k)
             END DO
             gisq=gisq+sum*sum
             w(i)=sum
          END DO
          !
          !     Test whether to replace the new quadratic model by the least Frobe
          !     norm interpolant, making the replacement if the test is satisfied.
          !
          itest=itest+1
          IF (gqsq < 1.0e2_dp*gisq) itest=0
          IF (itest >= 3) THEN
             DO i=1,n
                gq(i)=w(i)
             END DO
             DO ih=1,nh
                hq(ih)=zero
             END DO
             DO j=1,nptm
                w(j)=zero
                DO k=1,npt
                   w(j)=w(j)+vlag(k)*zmat(k,j)
                END DO
                IF (j < idz) w(j)=-w(j)
             END DO
             DO k=1,npt
                pq(k)=zero
                DO j=1,nptm
                   pq(k)=pq(k)+zmat(k,j)*w(j)
                END DO
             END DO
             itest=0
          END IF
       END IF
    END IF
    IF (f < fsave) kopt=knew
    !
    !     If a trust region step has provided a sufficient decrease in F, th
    !     branch for another trust region calculation. The case KSAVE>0 occu
    !     when the new function value was calculated by a model step.
    !
    IF (f <= fsave+tenth*vquad) GOTO 100
    IF (ksave > 0) GOTO 100
    !
    !     Alternatively, find out if the interpolation points are close enou
    !     to the best point so far.
    !
    knew=0
460 distsq=4.0_dp*delta*delta
    DO k=1,npt
       sum=zero
       DO j=1,n
          sum=sum+(xpt(k,j)-xopt(j))**2
       END DO
       IF (sum > distsq) THEN
          knew=k
          distsq=sum
       END IF
    END DO
    !
    !     If KNEW is positive, then set DSTEP, and branch back for the next
    !     iteration, which will generate a "model step".
    !
    IF (knew > 0) THEN
       dstep=MAX(MIN(tenth*SQRT(distsq),half*delta),rho)
       dsq=dstep*dstep
       GOTO 120
    END IF
    IF (ratio > zero) GOTO 100
    IF (MAX(delta,dnorm) > rho) GOTO 100
    !
    !     The calculations with the current value of RHO are complete. Pick
    !     next values of RHO and DELTA.
    !
490 IF (rho > rhoend) THEN
       delta=half*rho
       ratio=rho/rhoend
       IF (ratio <= 16.0_dp) THEN
          rho=rhoend
       ELSE IF (ratio <= 250.0_dp) THEN
          rho=SQRT(ratio)*rhoend
       ELSE
          rho=tenth*rho
       END IF
       delta=MAX(delta,rho)
       GOTO 90
    END IF
    !
    !     Return from the calculation, after another Newton-Raphson step, if
    !     it is too short to have been tried before.
    !
    IF (knew == -1) GOTO 290
    opt%state = 7
    CALL get_state
530 IF (fopt <= f) THEN
       DO i=1,n
          x(i)=xbase(i)+xopt(i)
       END DO
       f=fopt
    END IF

    CALL get_state

    !******************************************************************************
  CONTAINS
    !******************************************************************************
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE get_state()
      opt%np       = np
      opt%nh       = nh
      opt%nptm     = nptm
      opt%nftest   = nftest
      opt%idz      = idz
      opt%itest    = itest
      opt%nf       = nf
      opt%nfm      = nfm
      opt%nfmm     = nfmm
      opt%nfsav    = nfsav
      opt%knew     = knew
      opt%kopt     = kopt
      opt%ksave    = ksave
      opt%ktemp    = ktemp
      opt%rhosq    = rhosq
      opt%recip    = recip
      opt%reciq    = reciq
      opt%fbeg     = fbeg
      opt%fopt     = fopt
      opt%diffa    = diffa
      opt%xoptsq   = xoptsq
      opt%rho      = rho
      opt%delta    = delta
      opt%dsq      = dsq
      opt%dnorm    = dnorm
      opt%ratio    = ratio
      opt%temp     = temp
      opt%tempq    = tempq
      opt%beta     = beta
      opt%dx       = dx
      opt%vquad    = vquad
      opt%diff     = diff
      opt%diffc    = diffc
      opt%diffb    = diffb
      opt%fsave    = fsave
      opt%detrat   = detrat
      opt%hdiag    = hdiag
      opt%distsq   = distsq
      opt%gisq     = gisq
      opt%gqsq     = gqsq
      opt%f        = f
      opt%bstep    = bstep
      opt%alpha    = alpha
      opt%dstep    = dstep
    END SUBROUTINE get_state

    !******************************************************************************
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE set_state()
      np           = opt%np
      nh           = opt%nh
      nptm         = opt%nptm
      nftest       = opt%nftest
      idz          = opt%idz
      itest        = opt%itest
      nf           = opt%nf
      nfm          = opt%nfm
      nfmm         = opt%nfmm
      nfsav        = opt%nfsav
      knew         = opt%knew
      kopt         = opt%kopt
      ksave        = opt%ksave
      ktemp        = opt%ktemp
      rhosq        = opt%rhosq
      recip        = opt%recip
      reciq        = opt%reciq
      fbeg         = opt%fbeg
      fopt         = opt%fopt
      diffa        = opt%diffa
      xoptsq       = opt%xoptsq
      rho          = opt%rho
      delta        = opt%delta
      dsq          = opt%dsq
      dnorm        = opt%dnorm
      ratio        = opt%ratio
      temp         = opt%temp
      tempq        = opt%tempq
      beta         = opt%beta
      dx           = opt%dx
      vquad        = opt%vquad
      diff         = opt%diff
      diffc        = opt%diffc
      diffb        = opt%diffb
      fsave        = opt%fsave
      detrat       = opt%detrat
      hdiag        = opt%hdiag
      distsq       = opt%distsq
      gisq         = opt%gisq
      gqsq         = opt%gqsq
      f            = opt%f
      bstep        = opt%bstep
      alpha        = opt%alpha
      dstep        = opt%dstep
    END SUBROUTINE set_state

  END SUBROUTINE newuob

!******************************************************************************

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param npt ...
!> \param xopt ...
!> \param xpt ...
!> \param bmat ...
!> \param zmat ...
!> \param idz ...
!> \param ndim ...
!> \param kopt ...
!> \param knew ...
!> \param d ...
!> \param w ...
!> \param vlag ...
!> \param beta ...
!> \param s ...
!> \param wvec ...
!> \param prod ...
! *****************************************************************************
  SUBROUTINE bigden (n,npt,xopt,xpt,bmat,zmat,idz,ndim,kopt,&
       knew,d,w,vlag,beta,s,wvec,prod)

    INTEGER, INTENT(inout)                   :: n, npt
    REAL(dp), DIMENSION(*), INTENT(inout)    :: xopt
    REAL(dp), DIMENSION(npt, *), &
      INTENT(inout)                          :: xpt
    INTEGER, INTENT(inout)                   :: ndim, idz
    REAL(dp), DIMENSION(npt, *), &
      INTENT(inout)                          :: zmat
    REAL(dp), DIMENSION(ndim, *), &
      INTENT(inout)                          :: bmat
    INTEGER, INTENT(inout)                   :: kopt, knew
    REAL(dp), DIMENSION(*), INTENT(inout)    :: d, w, vlag
    REAL(dp), INTENT(inout)                  :: beta
    REAL(dp), DIMENSION(*), INTENT(inout)    :: s
    REAL(dp), DIMENSION(ndim, *), &
      INTENT(inout)                          :: wvec, prod

    REAL(dp), PARAMETER                      :: half = 0.5_dp, one = 1._dp, &
                                                quart = 0.25_dp, two = 2._dp, &
                                                zero = 0._dp

    INTEGER                                  :: i, ip, isave, iterc, iu, j, &
                                                jc, k, ksav, nptm, nw
    REAL(dp) :: alpha, angle, dd, denmax, denold, densav, diff, ds, dstemp, &
      dtest, ss, ssden, sstemp, step, sum, sumold, tau, temp, tempa, tempb, &
      tempc, xoptd, xopts, xoptsq
    REAL(dp), DIMENSION(9)                   :: den, denex, par

!
!     N is the number of variables.
!     NPT is the number of interpolation equations.
!     XOPT is the best interpolation point so far.
!     XPT contains the coordinates of the current interpolation points.
!     BMAT provides the last N columns of H.
!     ZMAT and IDZ give a factorization of the first NPT by NPT submatri
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     KOPT is the index of the optimal interpolation point.
!     KNEW is the index of the interpolation point that is going to be m
!     D will be set to the step from XOPT to the new point, and on entry
!       should be the D that was calculated by the last call of BIGLAG.
!       length of the initial D provides a trust region bound on the fin
!     W will be set to Wcheck for the final choice of D.
!     VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
!     BETA will be set to the value that will occur in the updating form
!       when the KNEW-th interpolation point is moved to its new positio
!     S, WVEC, PROD and the private arrays DEN, DENEX and PAR will be us
!       for working space.
!
!     D is calculated in a way that should provide a denominator with a
!     modulus in the updating formula when the KNEW-th interpolation poi
!     shifted to the new position XOPT+D.
!

    nptm=npt-n-1
    !
    !     Store the first NPT elements of the KNEW-th column of H in W(N+1)
    !     to W(N+NPT).
    !
    DO k=1,npt
       w(n+k)=zero
    END DO
    DO j=1,nptm
       temp=zmat(knew,j)
       IF (j < idz) temp=-temp
       DO k=1,npt
          w(n+k)=w(n+k)+temp*zmat(k,j)
       END DO
    END DO
    alpha=w(n+knew)
    !
    !     The initial search direction D is taken from the last call of BIGL
    !     and the initial S is set below, usually to the direction from X_OP
    !     to X_KNEW, but a different direction to an interpolation point may
    !     be chosen, in order to prevent S from being nearly parallel to D.
    !
    dd=zero
    ds=zero
    ss=zero
    xoptsq=zero
    DO i=1,n
       dd=dd+d(i)**2
       s(i)=xpt(knew,i)-xopt(i)
       ds=ds+d(i)*s(i)
       ss=ss+s(i)**2
       xoptsq=xoptsq+xopt(i)**2
    END DO
    IF (ds*ds > 0.99_dp*dd*ss) THEN
       ksav=knew
       dtest=ds*ds/ss
       DO k=1,npt
          IF (k /= kopt) THEN
             dstemp=zero
             sstemp=zero
             DO i=1,n
                diff=xpt(k,i)-xopt(i)
                dstemp=dstemp+d(i)*diff
                sstemp=sstemp+diff*diff
             END DO
             IF (dstemp*dstemp/sstemp < dtest) THEN
                ksav=k
                dtest=dstemp*dstemp/sstemp
                ds=dstemp
                ss=sstemp
             END IF
          END IF
       END DO
       DO i=1,n
          s(i)=xpt(ksav,i)-xopt(i)
       END DO
    END IF
    ssden=dd*ss-ds*ds
    iterc=0
    densav=zero
    !
    !     Begin the iteration by overwriting S with a vector that has the
    !     required length and direction.
    !
    mainloop : DO
       iterc=iterc+1
       temp=one/SQRT(ssden)
       xoptd=zero
       xopts=zero
       DO i=1,n
          s(i)=temp*(dd*s(i)-ds*d(i))
          xoptd=xoptd+xopt(i)*d(i)
          xopts=xopts+xopt(i)*s(i)
       END DO
       !
       !     Set the coefficients of the first two terms of BETA.
       !
       tempa=half*xoptd*xoptd
       tempb=half*xopts*xopts
       den(1)=dd*(xoptsq+half*dd)+tempa+tempb
       den(2)=two*xoptd*dd
       den(3)=two*xopts*dd
       den(4)=tempa-tempb
       den(5)=xoptd*xopts
       DO i=6,9
          den(i)=zero
       END DO
       !
       !     Put the coefficients of Wcheck in WVEC.
       !
       DO k=1,npt
          tempa=zero
          tempb=zero
          tempc=zero
          DO i=1,n
             tempa=tempa+xpt(k,i)*d(i)
             tempb=tempb+xpt(k,i)*s(i)
             tempc=tempc+xpt(k,i)*xopt(i)
          END DO
          wvec(k,1)=quart*(tempa*tempa+tempb*tempb)
          wvec(k,2)=tempa*tempc
          wvec(k,3)=tempb*tempc
          wvec(k,4)=quart*(tempa*tempa-tempb*tempb)
          wvec(k,5)=half*tempa*tempb
       END DO
       DO i=1,n
          ip=i+npt
          wvec(ip,1)=zero
          wvec(ip,2)=d(i)
          wvec(ip,3)=s(i)
          wvec(ip,4)=zero
          wvec(ip,5)=zero
       END DO
       !
       !     Put the coefficents of THETA*Wcheck in PROD.
       !
       DO jc=1,5
          nw=npt
          IF (jc == 2 .OR. jc == 3) nw=ndim
          DO k=1,npt
             prod(k,jc)=zero
          END DO
          DO j=1,nptm
             sum=zero
             DO k=1,npt
                sum=sum+zmat(k,j)*wvec(k,jc)
             END DO
             IF (j < idz) sum=-sum
             DO k=1,npt
                prod(k,jc)=prod(k,jc)+sum*zmat(k,j)
             END DO
          END DO
          IF (nw == ndim) THEN
             DO k=1,npt
                sum=zero
                DO j=1,n
                   sum=sum+bmat(k,j)*wvec(npt+j,jc)
                END DO
                prod(k,jc)=prod(k,jc)+sum
             END DO
          END IF
          DO j=1,n
             sum=zero
             DO i=1,nw
                sum=sum+bmat(i,j)*wvec(i,jc)
             END DO
             prod(npt+j,jc)=sum
          END DO
       END DO
       !
       !     Include in DEN the part of BETA that depends on THETA.
       !
       DO k=1,ndim
          sum=zero
          DO I=1,5
             par(i)=half*prod(k,i)*wvec(k,i)
             sum=sum+par(i)
          END DO
          den(1)=den(1)-par(1)-sum
          tempa=prod(k,1)*wvec(k,2)+prod(k,2)*wvec(k,1)
          tempb=prod(k,2)*wvec(k,4)+prod(k,4)*wvec(k,2)
          tempc=prod(k,3)*wvec(k,5)+prod(k,5)*wvec(k,3)
          den(2)=den(2)-tempa-half*(tempb+tempc)
          den(6)=den(6)-half*(tempb-tempc)
          tempa=prod(k,1)*wvec(k,3)+prod(k,3)*wvec(k,1)
          tempb=prod(k,2)*wvec(k,5)+prod(k,5)*wvec(k,2)
          tempc=prod(k,3)*wvec(k,4)+prod(k,4)*wvec(k,3)
          den(3)=den(3)-tempa-half*(tempb-tempc)
          den(7)=den(7)-half*(tempb+tempc)
          tempa=prod(k,1)*wvec(k,4)+prod(k,4)*wvec(k,1)
          den(4)=den(4)-tempa-par(2)+par(3)
          tempa=prod(k,1)*wvec(k,5)+prod(k,5)*wvec(k,1)
          tempb=prod(k,2)*wvec(k,3)+prod(k,3)*wvec(k,2)
          den(5)=den(5)-tempa-half*tempb
          den(8)=den(8)-par(4)+par(5)
          tempa=prod(k,4)*wvec(k,5)+prod(k,5)*wvec(k,4)
          den(9)=den(9)-half*tempa
       END DO
       !
       !     Extend DEN so that it holds all the coefficients of DENOM.
       !
       sum=zero
       DO i=1,5
          par(i)=half*prod(knew,i)**2
          sum=sum+par(i)
       END DO
       denex(1)=alpha*den(1)+par(1)+sum
       tempa=two*prod(knew,1)*prod(knew,2)
       tempb=prod(knew,2)*prod(knew,4)
       tempc=prod(knew,3)*prod(knew,5)
       denex(2)=alpha*den(2)+tempa+tempb+tempc
       denex(6)=alpha*den(6)+tempb-tempc
       tempa=two*prod(knew,1)*prod(knew,3)
       tempb=prod(knew,2)*prod(knew,5)
       tempc=prod(knew,3)*prod(knew,4)
       denex(3)=alpha*den(3)+tempa+tempb-tempc
       denex(7)=alpha*den(7)+tempb+tempc
       tempa=two*prod(knew,1)*prod(knew,4)
       denex(4)=alpha*den(4)+tempa+par(2)-par(3)
       tempa=two*prod(knew,1)*prod(knew,5)
       denex(5)=alpha*den(5)+tempa+prod(knew,2)*prod(knew,3)
       denex(8)=alpha*den(8)+par(4)-par(5)
       denex(9)=alpha*den(9)+prod(knew,4)*prod(knew,5)
       !
       !     Seek the value of the angle that maximizes the modulus of DENOM.
       !
       sum=denex(1)+denex(2)+denex(4)+denex(6)+denex(8)
       denold=sum
       denmax=sum
       isave=0
       iu=49
       temp=twopi/REAL(IU+1,dp)
       par(1)=one
       DO i=1,iu
          angle=REAL(i,dp)*temp
          par(2)=COS(angle)
          par(3)=SIN(angle)
          DO j=4,8,2
             par(j)=par(2)*par(j-2)-par(3)*par(j-1)
             par(j+1)=par(2)*par(j-1)+par(3)*par(j-2)
          END DO
          sumold=sum
          sum=zero
          DO j=1,9
             sum=sum+denex(j)*par(j)
          END DO
          IF (ABS(sum) > ABS(denmax)) THEN
             denmax=sum
             isave=i
             tempa=sumold
          ELSE IF (i == isave+1) THEN
             tempb=sum
          END IF
       END DO
       IF (isave == 0) tempa=sum
       IF (isave == iu) tempb=denold
       step=zero
       IF (tempa /= tempb) THEN
          tempa=tempa-denmax
          tempb=tempb-denmax
          step=half*(tempa-tempb)/(tempa+tempb)
       END IF
       angle=temp*(REAL(isave,dp)+step)
       !
       !     Calculate the new parameters of the denominator, the new VLAG vect
       !     and the new D. Then test for convergence.
       !
       par(2)=COS(angle)
       par(3)=SIN(angle)
       DO j=4,8,2
          par(j)=par(2)*par(j-2)-par(3)*par(j-1)
          par(j+1)=par(2)*par(j-1)+par(3)*par(j-2)
       END DO
       beta=zero
       denmax=zero
       DO j=1,9
          beta=beta+den(j)*par(j)
          denmax=denmax+denex(j)*par(j)
       END DO
       DO k=1,ndim
          vlag(k)=zero
          DO j=1,5
             vlag(k)=vlag(k)+prod(k,j)*par(j)
          END DO
       END DO
       tau=vlag(knew)
       dd=zero
       tempa=zero
       tempb=zero
       DO i=1,n
          d(i)=par(2)*d(i)+par(3)*s(i)
          w(i)=xopt(i)+d(i)
          dd=dd+d(i)**2
          tempa=tempa+d(i)*w(i)
          tempb=tempb+w(i)*w(i)
       END DO
       IF (iterc >= n) EXIT mainloop
       IF (iterc >= 1) densav=MAX(densav,denold)
       IF (ABS(denmax) <= 1.1_dp*ABS(densav)) EXIT mainloop
       densav=denmax
       !
       !     Set S to half the gradient of the denominator with respect to D.
       !     Then branch for the next iteration.
       !
       DO i=1,n
          temp=tempa*xopt(i)+tempb*d(i)-vlag(npt+i)
          s(i)=tau*bmat(knew,i)+alpha*temp
       END DO
       DO k=1,npt
          sum=zero
          DO j=1,n
             sum=sum+xpt(k,j)*w(j)
          END DO
          temp=(tau*w(n+k)-alpha*vlag(k))*sum
          DO i=1,n
             s(i)=s(i)+temp*xpt(k,i)
          END DO
       END DO
       ss=zero
       ds=zero
       DO i=1,n
          ss=ss+s(i)**2
          ds=ds+d(i)*s(i)
       END DO
       ssden=dd*ss-ds*ds
       IF (ssden < 1.0e-8_dp*dd*ss) EXIT mainloop
    END DO mainloop
    !
    !     Set the vector W before the RETURN from the subroutine.
    !
    DO k=1,ndim
       w(k)=zero
       DO j=1,5
          w(k)=w(k)+wvec(k,j)*par(j)
       END DO
    END DO
    vlag(kopt)=vlag(kopt)+one

  END SUBROUTINE bigden
!******************************************************************************

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param npt ...
!> \param xopt ...
!> \param xpt ...
!> \param bmat ...
!> \param zmat ...
!> \param idz ...
!> \param ndim ...
!> \param knew ...
!> \param delta ...
!> \param d ...
!> \param alpha ...
!> \param hcol ...
!> \param gc ...
!> \param gd ...
!> \param s ...
!> \param w ...
! *****************************************************************************
  SUBROUTINE biglag (n,npt,xopt,xpt,bmat,zmat,idz,ndim,knew,&
       delta,d,alpha,hcol,gc,gd,s,w)
    INTEGER, INTENT(inout)                   :: n, npt
    REAL(dp), DIMENSION(*), INTENT(inout)    :: xopt
    REAL(dp), DIMENSION(npt, *), &
      INTENT(inout)                          :: xpt
    INTEGER, INTENT(inout)                   :: ndim, idz
    REAL(dp), DIMENSION(npt, *), &
      INTENT(inout)                          :: zmat
    REAL(dp), DIMENSION(ndim, *), &
      INTENT(inout)                          :: bmat
    INTEGER, INTENT(inout)                   :: knew
    REAL(dp), INTENT(inout)                  :: delta
    REAL(dp), DIMENSION(*), INTENT(inout)    :: d
    REAL(dp), INTENT(inout)                  :: alpha
    REAL(dp), DIMENSION(*), INTENT(inout)    :: hcol, gc, gd, s, w

    REAL(dp), PARAMETER                      :: half = 0.5_dp, one = 1._dp, &
                                                zero = 0._dp

    INTEGER                                  :: i, isave, iterc, iu, j, k, &
                                                nptm
    REAL(dp) :: angle, cf1, cf2, cf3, cf4, cf5, cth, dd, delsq, denom, dhd, &
      gg, scale, sp, ss, step, sth, sum, tau, taubeg, taumax, tauold, temp, &
      tempa, tempb

!
!
!     N is the number of variables.
!     NPT is the number of interpolation equations.
!     XOPT is the best interpolation point so far.
!     XPT contains the coordinates of the current interpolation points.
!     BMAT provides the last N columns of H.
!     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     KNEW is the index of the interpolation point that is going to be m
!     DELTA is the current trust region bound.
!     D will be set to the step from XOPT to the new point.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     HCOL, GC, GD, S and W will be used for working space.
!
!     The step D is calculated in a way that attempts to maximize the mo
!     of LFUNC(XOPT+D), subject to the bound ||D|| .LE. DELTA, where LFU
!     the KNEW-th Lagrange function.
!

    delsq=delta*delta
    nptm=npt-n-1
    !
    !     Set the first NPT components of HCOL to the leading elements of th
    !     KNEW-th column of H.
    !
    iterc=0
    DO k=1,npt
       hcol(k)=zero
    END DO
    DO j=1,nptm
       temp=zmat(knew,j)
       IF (j < idz) temp=-temp
       DO k=1,npt
          hcol(k)=hcol(k)+temp*zmat(k,j)
       END DO
    END DO
    alpha=hcol(knew)
    !
    !     Set the unscaled initial direction D. Form the gradient of LFUNC a
    !     XOPT, and multiply D by the second derivative matrix of LFUNC.
    !
    dd=zero
    DO i=1,n
       d(i)=xpt(knew,i)-xopt(i)
       gc(i)=bmat(knew,i)
       gd(i)=zero
       dd=dd+d(i)**2
    END DO
    DO k=1,npt
       temp=zero
       sum=zero
       DO j=1,n
          temp=temp+xpt(k,j)*xopt(j)
          sum=sum+xpt(k,j)*d(j)
       END DO
       temp=hcol(k)*temp
       sum=hcol(k)*sum
       DO i=1,n
          gc(i)=gc(i)+temp*xpt(k,i)
          gd(i)=gd(i)+sum*xpt(k,i)
       END DO
    END DO
    !
    !     Scale D and GD, with a sign change if required. Set S to another
    !     vector in the initial two dimensional subspace.
    !
    gg=zero
    sp=zero
    dhd=zero
    DO i=1,n
       gg=gg+gc(i)**2
       sp=sp+d(i)*gc(i)
       dhd=dhd+d(i)*gd(i)
    END DO
    scale=delta/SQRT(dd)
    IF (sp*dhd < zero) scale=-scale
    temp=zero
    IF (sp*sp > 0.99_dp*dd*gg) temp=one
    tau=scale*(ABS(sp)+half*scale*ABS(dhd))
    IF (gg*delsq < 0.01_dp*tau*tau) temp=one
    DO i=1,n
       d(i)=scale*d(i)
       gd(i)=scale*gd(i)
       s(i)=gc(i)+temp*gd(i)
    END DO
    !
    !     Begin the iteration by overwriting S with a vector that has the
    !     required length and direction, except that termination occurs if
    !     the given D and S are nearly parallel.
    !
    mainloop : DO
       iterc=iterc+1
       dd=zero
       sp=zero
       ss=zero
       DO i=1,n
          dd=dd+d(i)**2
          sp=sp+d(i)*s(i)
          ss=ss+s(i)**2
       END DO
       temp=dd*ss-sp*sp
       IF (temp <= 1.0e-8_dp*dd*ss) EXIT mainloop
       denom=SQRT(temp)
       DO i=1,n
          s(i)=(dd*s(i)-sp*d(i))/denom
          w(i)=zero
       END DO
       !
       !     Calculate the coefficients of the objective function on the circle
       !     beginning with the multiplication of S by the second derivative ma
       !
       DO k=1,npt
          sum=zero
          DO j=1,n
             sum=sum+xpt(k,j)*s(j)
          END DO
          sum=hcol(k)*sum
          DO i=1,n
             w(i)=w(i)+sum*xpt(k,i)
          END DO
       END DO
       cf1=zero
       cf2=zero
       cf3=zero
       cf4=zero
       cf5=zero
       DO i=1,n
          cf1=cf1+s(i)*w(i)
          cf2=cf2+d(i)*gc(i)
          cf3=cf3+s(i)*gc(i)
          cf4=cf4+d(i)*gd(i)
          cf5=cf5+s(i)*gd(i)
       END DO
       cf1=half*cf1
       cf4=half*cf4-cf1
       !
       !     Seek the value of the angle that maximizes the modulus of TAU.
       !
       taubeg=cf1+cf2+cf4
       taumax=taubeg
       tauold=taubeg
       isave=0
       iu=49
       temp=twopi/REAL(iu+1,DP)
       DO i=1,iu
          angle=REAL(i,dp)*temp
          cth=COS(angle)
          sth=SIN(angle)
          tau=cf1+(cf2+cf4*cth)*cth+(cf3+cf5*cth)*sth
          IF (ABS(tau) > ABS(taumax)) THEN
             taumax=tau
             isave=i
             tempa=tauold
          ELSE IF (i == isave+1) THEN
             tempb=taU
          END IF
          tauold=tau
       END DO
       IF (isave == 0) tempa=tau
       IF (isave == iu) tempb=taubeg
       step=zero
       IF (tempa /= tempb) THEN
          tempa=tempa-taumax
          tempb=tempb-taumax
          step=half*(tempa-tempb)/(tempa+tempb)
       END IF
       angle=temp*(REAL(isave,DP)+step)
       !
       !     Calculate the new D and GD. Then test for convergence.
       !
       cth=COS(angle)
       sth=SIN(angle)
       tau=cf1+(cf2+cf4*cth)*cth+(cf3+cf5*cth)*sth
       DO i=1,n
          d(i)=cth*d(i)+sth*s(i)
          gd(i)=cth*gd(i)+sth*w(i)
          s(i)=gc(i)+gd(i)
       END DO
       IF (ABS(tau) <= 1.1_dp*ABS(taubeg)) EXIT mainloop
       IF (iterc >= n) EXIT mainloop
    END DO mainloop

  END SUBROUTINE biglag
!******************************************************************************

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param npt ...
!> \param xopt ...
!> \param xpt ...
!> \param gq ...
!> \param hq ...
!> \param pq ...
!> \param delta ...
!> \param step ...
!> \param d ...
!> \param g ...
!> \param hd ...
!> \param hs ...
!> \param crvmin ...
! *****************************************************************************
  SUBROUTINE trsapp (n,npt,xopt,xpt,gq,hq,pq,delta,step,d,g,hd,hs,crvmin)

    INTEGER, INTENT(INOUT)                   :: n, npt
    REAL(dp), DIMENSION(*), INTENT(INOUT)    :: xopt
    REAL(dp), DIMENSION(npt, *), &
      INTENT(INOUT)                          :: xpt
    REAL(dp), DIMENSION(*), INTENT(INOUT)    :: gq, hq, pq
    REAL(dp), INTENT(INOUT)                  :: delta
    REAL(dp), DIMENSION(*), INTENT(INOUT)    :: step, d, g, hd, hs
    REAL(dp), INTENT(INOUT)                  :: crvmin

    REAL(dp), PARAMETER                      :: half = 0.5_dp, zero = 0.0_dp

    INTEGER                                  :: i, isave, iterc, itermax, &
                                                itersw, iu, j
    LOGICAL                                  :: jump1, jump2
    REAL(dp) :: alpha, angle, angtest, bstep, cf, cth, dd, delsq, dg, dhd, &
      dhs, ds, gg, ggbeg, ggsav, qadd, qbeg, qmin, qnew, qred, qsav, ratio, &
      reduc, sg, sgk, shs, ss, sth, temp, tempa, tempb

!
!   N is the number of variables of a quadratic objective function, Q
!   The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meani
!     in order to define the current quadratic model Q.
!   DELTA is the trust region radius, and has to be positive.
!   STEP will be set to the calculated trial step.
!   The arrays D, G, HD and HS will be used for working space.
!   CRVMIN will be set to the least curvature of H along the conjugate
!     directions that occur, except that it is set to zero if STEP goe
!     all the way to the trust region boundary.
!
!   The calculation of STEP begins with the truncated conjugate gradient
!   method. If the boundary of the trust region is reached, then further
!   changes to STEP may be made, each one being in the 2D space spanned
!   by the current STEP and the corresponding gradient of Q. Thus STEP
!   should provide a substantial reduction to Q within the trust region
!
!   Initialization, which includes setting HD to H times XOPT.
!

    delsq=delta*delta
    iterc=0
    itermax=n
    itersw=itermax
    DO i=1,n
       d(i)=xopt(i)
    END DO
    CALL updatehd
    !
    !   Prepare for the first line search.
    !
    qred=zero
    dd=zero
    DO i=1,n
       step(i)=zero
       hs(i)=zero
       g(i)=gq(i)+hd(i)
       d(i)=-g(i)
       dd=dd+d(i)**2
    END DO
    crvmin=zero
    IF (dd == zero) RETURN
    ds=zero
    ss=zero
    gg=dd
    ggbeg=gg
    !
    !   Calculate the step to the trust region boundary and the product HD
    !
    jump1 = .FALSE.
    jump2 = .FALSE.
    mainloop : DO
       IF ( .NOT. jump2 ) THEN
          IF ( .NOT. jump1 ) THEN
             iterc=iterc+1
             temp=delsq-ss
             bstep=temp/(ds+SQRT(ds*ds+dd*temp))
             CALL updatehd
          END IF
          jump1 = .FALSE.
          IF (iterc <= itersw) THEN
             dhd=zero
             DO j=1,n
                dhd=dhd+d(j)*hd(j)
             END DO
             !
             !     Update CRVMIN and set the step-length ALPHA.
             !
             alpha=bstep
             IF (dhd > zero) THEN
                temp=dhd/dd
                IF (iterc == 1) crvmin=temp
                crvmin=MIN(crvmin,temp)
                alpha=MIN(alpha,gg/dhd)
             END IF
             qadd=alpha*(gg-half*alpha*dhd)
             qred=qred+qadd
             !
             !     Update STEP and HS.
             !
             ggsav=gg
             gg=zero
             DO i=1,n
                step(i)=step(i)+alpha*d(i)
                hs(i)=hs(i)+alpha*hd(i)
                gg=gg+(g(i)+hs(i))**2
             END DO
             !
             !     Begin another conjugate direction iteration if required.
             !
             IF (alpha < bstep) THEN
                IF (qadd <= 0.01_dp*qred) EXIT mainloop
                IF (gg <= 1.0e-4_dp*ggbeg) EXIT mainloop
                IF (iterc == itermax) EXIT mainloop
                temp=gg/ggsav
                dd=zero
                ds=zero
                ss=zero
                DO i=1,n
                   d(i)=temp*d(i)-g(i)-hs(i)
                   dd=dd+d(i)**2
                   ds=ds+d(i)*step(I)
                   ss=ss+step(i)**2
                END DO
                IF (ds <= zero) EXIT mainloop
                IF (ss < delsq) CYCLE mainloop
             END IF
             crvmin=zero
             itersw=iterc
             jump2 = .TRUE.
             IF (gg <= 1.0e-4_dp*ggbeg) EXIT mainloop
          ELSE
             jump2 = .FALSE.
          END IF
       END IF
       !
       !     Test whether an alternative iteration is required.
       !
!!!!  IF (gg <= 1.0e-4_dp*ggbeg) EXIT mainloop
       IF (jump2) THEN
          sg=zero
          shs=zero
          DO i=1,n
             sg=sg+step(i)*g(i)
             shs=shs+step(i)*hs(i)
          END DO
          sgk=sg+shs
          angtest=sgk/SQRT(gg*delsq)
          IF (angtest <= -0.99_dp) EXIT mainloop
          !
          !     Begin the alternative iteration by calculating D and HD and some
          !     scalar products.
          !
          iterc=iterc+1
          temp=SQRT(delsq*gg-sgk*sgk)
          tempa=delsq/temp
          tempb=sgk/temp
          DO i=1,n
             d(i)=tempa*(g(i)+hs(i))-tempb*step(i)
          END DO
          CALL updatehd
          IF (iterc <= itersw) THEN
             jump1 = .TRUE.
             CYCLE mainloop
          END IF
       END IF
       dg=zero
       dhd=zero
       dhs=zero
       DO i=1,n
          dg=dg+d(i)*g(i)
          dhd=dhd+hd(i)*d(i)
          dhs=dhs+hd(i)*step(i)
       END DO
       !
       !     Seek the value of the angle that minimizes Q.
       !
       cf=half*(shs-dhd)
       qbeg=sg+cf
       qsav=qbeg
       qmin=qbeg
       isave=0
       iu=49
       temp=twopi/REAL(iu+1,dp)
       DO i=1,iu
          angle=REAL(i,dp)*temp
          cth=COS(angle)
          sth=SIN(angle)
          qnew=(sg+cf*cth)*cth+(dg+dhs*cth)*sth
          IF (qnew < qmin) THEN
             qmin=qnew
             isave=i
             tempa=qsav
          ELSE IF (i == isave+1) THEN
             tempb=qnew
          END IF
          qsav=qnew
       END DO
       IF (isave == zero) tempa=qnew
       IF (isave == iu) tempb=qbeg
       angle=zero
       IF (tempa /= tempb) THEN
          tempa=tempa-qmin
          tempb=tempb-qmin
          angle=half*(tempa-tempb)/(tempa+tempb)
       END IF
       angle=temp*(REAL(isave,DP)+angle)
       !
       !     Calculate the new STEP and HS. Then test for convergence.
       !
       cth=COS(angle)
       sth=SIN(angle)
       reduc=qbeg-(sg+cf*cth)*cth-(dg+dhs*cth)*sth
       gg=zero
       DO i=1,n
          step(i)=cth*step(i)+sth*d(i)
          hs(i)=cth*hs(i)+sth*hd(i)
          gg=gg+(g(i)+hs(i))**2
       END DO
       qred=qred+reduc
       ratio=reduc/qred
       IF (iterc < itermax .AND. ratio > 0.01_dp) THEN
          jump2 = .TRUE.
       ELSE
          EXIT mainloop
       END IF

       IF (gg <= 1.0e-4_dp*ggbeg) EXIT mainloop

    END DO mainloop

    !*******************************************************************************
  CONTAINS
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE updatehd
    INTEGER                                  :: i, ih, j, k

      DO i=1,n
         hd(i)=zero
      END DO
      DO k=1,npt
         temp=zero
         DO j=1,n
            temp=temp+xpt(k,j)*d(j)
         END DO
         temp=temp*pq(k)
         DO i=1,n
            hd(i)=hd(i)+temp*xpt(k,i)
         END DO
      END DO
      ih=0
      DO j=1,n
         DO i=1,j
            ih=ih+1
            IF (i < j) hd(j)=hd(j)+hq(ih)*d(i)
            hd(i)=hd(i)+hq(ih)*d(j)
         END DO
      END DO
    END SUBROUTINE updatehd

  END SUBROUTINE trsapp
!******************************************************************************

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param npt ...
!> \param bmat ...
!> \param zmat ...
!> \param idz ...
!> \param ndim ...
!> \param vlag ...
!> \param beta ...
!> \param knew ...
!> \param w ...
! *****************************************************************************
  SUBROUTINE update (n,npt,bmat,zmat,idz,ndim,vlag,beta,knew,w)

    INTEGER, INTENT(IN)                      :: n, npt, ndim
    INTEGER, INTENT(INOUT)                   :: idz
    REAL(dp), DIMENSION(npt, *), &
      INTENT(INOUT)                          :: zmat
    REAL(dp), DIMENSION(ndim, *), &
      INTENT(INOUT)                          :: bmat
    REAL(dp), DIMENSION(*), INTENT(INOUT)    :: vlag
    REAL(dp), INTENT(INOUT)                  :: beta
    INTEGER, INTENT(INOUT)                   :: knew
    REAL(dp), DIMENSION(*), INTENT(INOUT)    :: w

    REAL(dp), PARAMETER                      :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: i, iflag, j, ja, jb, jl, jp, &
                                                nptm
    REAL(dp)                                 :: alpha, denom, scala, scalb, &
                                                tau, tausq, temp, tempa, tempb

!   The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
!   interpolation point that has index KNEW. On entry, VLAG contains the
!   components of the vector Theta*Wcheck+e_b of the updating formula
!   (6.11), and BETA holds the value of the parameter that has this na
!   The vector W is used for working space.
!

    nptm=npt-n-1
    !
    !     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
    !
    jl=1
    DO j=2,nptm
       IF (j == idz) THEN
          jl=idz
       ELSE IF (zmat(knew,j) /= zero) THEN
          temp=SQRT(zmat(knew,jl)**2+zmat(knew,j)**2)
          tempa=zmat(knew,jl)/temp
          tempb=zmat(knew,j)/temp
          DO  I=1,NPT
             temp=tempa*zmat(i,jl)+tempb*zmat(i,j)
             zmat(i,j)=tempa*zmat(i,j)-tempb*zmat(i,jl)
             zmat(i,jl)=temp
          END DO
          zmat(knew,j)=zero
       END IF
    END DO
    !
    !   Put the first NPT components of the KNEW-th column of HLAG into W,
    !   and calculate the parameters of the updating formula.
    !
    tempa=zmat(knew,1)
    IF (idz >= 2) tempa=-tempa
    IF (jl > 1) tempb=zmat(knew,jl)
    DO i=1,npt
       w(i)=tempa*zmat(i,1)
       IF (jl > 1) w(i)=w(i)+tempb*zmat(i,jl)
    END DO
    alpha=w(knew)
    tau=vlag(knew)
    tausq=tau*tau
    denom=alpha*beta+tausq
    vlag(knew)=vlag(knew)-one
    !
    !   Complete the updating of ZMAT when there is only one nonzero eleme
    !   in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to
    !   then the first column of ZMAT will be exchanged with another one l
    !
    iflag=0
    IF (JL == 1) THEN
       temp=SQRT(ABS(denom))
       tempb=tempa/temp
       tempa=tau/temp
       DO i=1,npt
          zmat(i,1)=tempa*zmat(i,1)-tempb*vlag(i)
       END DO
       IF (idz == 1 .AND. temp < zero) idz=2
       IF (idz >= 2 .AND. temp >= zero) iflag=1
    ELSE
       !
       !   Complete the updating of ZMAT in the alternative case.
       !
       ja=1
       IF (beta >= zero) ja=jl
       jb=jl+1-ja
       temp=zmat(knew,jb)/denom
       tempa=temp*beta
       tempb=temp*tau
       temp=zmat(knew,ja)
       scala=one/SQRT(ABS(beta)*temp*temp+tausq)
       scalb=scala*SQRT(ABS(denom))
       DO i=1,npt
          zmat(i,ja)=scala*(tau*zmat(i,ja)-temp*vlag(i))
          zmat(i,jb)=scalb*(zmat(i,jb)-tempa*w(i)-tempb*vlag(i))
       END DO
       IF (denom <= zero) THEN
          IF (beta < zero) idz=idz+1
          IF (beta >= zero) iflag=1
       END IF
    END IF
    !
    !   IDZ is reduced in the following case, and usually the first column
    !   of ZMAT is exchanged with a later one.
    !
    IF (iflag == 1) THEN
       idz=idz-1
       DO i=1,npt
          temp=zmat(i,1)
          zmat(i,1)=zmat(i,idz)
          zmat(i,idz)=temp
       END DO
    END IF
    !
    !   Finally, update the matrix BMAT.
    !
    DO j=1,n
       jp=npt+j
       w(jp)=bmat(knew,j)
       tempa=(alpha*vlag(jp)-tau*w(jp))/denom
       tempb=(-beta*w(jp)-tau*vlag(jp))/denom
       DO i=1,jp
          bmat(i,j)=bmat(i,j)+tempa*vlag(i)+tempb*w(i)
          IF (i > npt) bmat(jp,i-npt)=bmat(i,j)
       END DO
    END DO

  END SUBROUTINE update

END MODULE powell

!******************************************************************************
