! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!
#ifndef FMG_INTEGRATE_FULL
    CALL poly_p_eval2b(IF_CHECK(poly_jk,poly_jk(1)),size_jk,REAL(j,dp),&
         IF_CHECK(poly_k,poly_k(0)),size_k,npoly=npoly,grad=grad_val,xi=IF_CHECK(xi,xi(1)))
#endif
    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)

#if defined(FMG_INTEGRATE)||defined(FMG_INTEGRATE_FULL)
    k_vals=0.0_dp
#endif
    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
#if defined(FMG_INTEGRATE)||defined(FMG_INTEGRATE_FULL)
                gval=IF_FLAT(grid(ii+ij+ik+1),grid(ik,ij,ii))*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
#else
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
#if defined(FM_FLAT_GRID)
                grid(ii+ij+ik+1) = grid(ii+ij+ik+1) + p_v*res_k
#else
                IF ( PRESENT ( lgrid ) ) THEN
                  ig = ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig,ithread)=lgrid%r(ig,ithread) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF
#endif
#endif

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
#if defined(FMG_INTEGRATE)||defined(FMG_INTEGRATE_FULL)
                gval=IF_FLAT(grid(ii+ij+ik+1),grid(ik,ij,ii))*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
#else
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
#if defined(FM_FLAT_GRID)
                grid(ii+ij+ik+1) = grid(ii+ij+ik+1) + p_v*res_k
#else
                IF ( PRESENT ( lgrid ) ) THEN
                  ig = ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig,ithread)=lgrid%r(ig,ithread) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF
#endif
#endif
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
#if defined(FMG_INTEGRATE)||defined(FMG_INTEGRATE_FULL)
                gval=IF_FLAT(grid(ii+ij+ik+1),grid(ik,ij,ii))*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
#else
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
#if defined(FM_FLAT_GRID)
                grid(ii+ij+ik+1) = grid(ii+ij+ik+1) + p_v*res_k
#else
                IF ( PRESENT ( lgrid ) ) THEN
                  ig = ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig,ithread)=lgrid%r(ig,ithread) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF
#endif
#endif

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
#if defined(FMG_INTEGRATE)||defined(FMG_INTEGRATE_FULL)
                gval=IF_FLAT(grid(ii+ij+ik+1),grid(ik,ij,ii))*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
#else
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
#if defined(FM_FLAT_GRID)
                grid(ii+ij+ik+1) = grid(ii+ij+ik+1) + p_v*res_k
#else
                IF ( PRESENT ( lgrid ) ) THEN
                  ig = ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig,ithread)=lgrid%r(ig,ithread) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF
#endif
#endif
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
#ifdef FMG_INTEGRATE
    pShift=0
    DO ipoly=1,npoly
        p_v=0.0_dp
            p_v=p_v+poly_k(pShift+0)*k_vals(0)
            p_v=p_v+poly_k(pShift+1)*k_vals(1)
            p_v=p_v+poly_k(pShift+2)*k_vals(2)
        mres(ipoly)=mres(ipoly)+p_v
        pShift=pShift+grad_val+1
    END DO
#elif defined(FMG_INTEGRATE_FULL)
    CALL poly_padd_uneval2b(IF_CHECK(poly_jk,poly_jk(1)),size_jk,REAL(j,dp),IF_CHECK(k_vals,k_vals(0)),&
        size_k,npoly=npoly,grad=grad_val,xi=IF_CHECK(xi,xi(1)))
#endif
