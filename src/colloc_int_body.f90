
    INTEGER, DIMENSION(0:2), PARAMETER       :: permut = (/ 2,1,0 /)
    INTEGER :: grad, i, i0, ii, iiShift, iiShift2, iistart, iistart2, ij, &
      ijShift, iJump, ik, ikShift, ikShift2, ikstart, ikstart2, iend,iend2,&
      imax, imax1, imin, imin1, istart, istart2, j, jend, jJump, jmax, jmax1, jmin, &
      jmin1, jstart, k, kend, kend2, kgrad, kJump, kmax, kmax1, kmin, kmin1, &
      kstart, kstart2, max_j, stat,size_jk,size_k,size_ijk
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: k_bounds
    INTEGER, DIMENSION(0:2)                  :: cellShift, l_shift, l_ub, &
                                                ndim, period, shiftPos
    INTEGER, DIMENSION(2, 0:2)               :: l_bounds
    LOGICAL                                  :: failure, has_overlap, &
                                                k_bounds_alloc, poly_alloc
    REAL(dp) :: cci0, cci1, cci2, ccj0, ccj0_i0, ccj0_i1, ccj0_i2, ccj1, &
      ccj1_i0, ccj1_i1, ccj2, cck0, cck0_0, cck0_0_p, cck0_i, cck0_i2, &
      cck0_ij, cck0_j, cck0_j2, cck0_j_p, cck1, cck1_0, cck1_0_p, cck1_i, &
      cck1_j, cck2, delta_i, delta_j, delta_k, g_scale, i_coeffn_i, icoeff0, &
      ii_coeff0, ii_coeff2, ii_coeff2_jump, ii_coeffn, ii_coeffn_jump, &
      ij_coeff0, ij_coeff0_jump, ij_coeff1, ik_coeff0, ik_coeff0_jump, &
      ik_coeff1, j_coeffn_i, j_coeffn_j, jcoeff0, jj_coeff0, jj_coeff2, &
      jj_coeffn, jk_coeff0, jk_coeff1, k_coeffn_i, k_coeffn_j, k_coeffn_k, &
      kcoeff0, kk_coeff0, kk_coeff2, kk_coeffn, m(0:2,0:2), maxr2, p_kk, &
      r_0
    REAL(dp) :: res_0, res_i, res_j, res_k, scaled_h(0:2,0:2), sqDi, sqDj, &
      sqDk
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: poly_ijk, poly_jk,xi
    REAL(dp), DIMENSION(0:2)                 :: l, normD, p_shift, resPos, &
                                                resPosReal, riPos, rpos, wrPos

#ifdef FMG_INTEGRATE
    INTEGER :: ipoly,pShift
    REAL(dp) :: det, gval,p_v
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: mres, k_vals, poly_k
#elif defined(FMG_INTEGRATE_FULL)
    REAL(dp) :: det, gval
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: k_vals
    INTEGER, PARAMETER :: npoly=1
#else
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: poly_k
    REAL(dp) :: p_v
    INTEGER, PARAMETER :: npoly=1
#endif

    failure=.FALSE.
    k_bounds_alloc=.FALSE.
    poly_alloc=.FALSE.
#ifndef  FMG_INTEGRATE_FULL
    IF (ALL(poly==0.0_dp)) GOTO 21
#endif
    IF (PRESENT(poly_shift)) THEN
        p_shift=poly_shift
    ELSE
        p_shift=0.0_dp
    END IF
    IF (PRESENT(gdim)) THEN
        DO i=0,2
            ndim(permut(i))=gdim(i)
        END DO
    ELSE
        ndim(permut(0))=SIZE(grid,1)
        ndim(permut(1))=SIZE(grid,2)
        ndim(permut(2))=SIZE(grid,3)
    END IF
    g_scale=1.0_dp
    IF (PRESENT(scale)) g_scale=scale
#if defined(FMG_INTEGRATE)||defined(FMG_INTEGRATE_FULL)
    det=(h(0,0)*(h(1,1)*h(2,2)-h(1,2)*h(2,1))&
        -h(1,0)*(h(0,1)*h(2,2)-h(0,2)*h(2,1))&
        +h(2,0)*(h(0,1)*h(1,2)-h(0,2)*h(1,1)))
    g_scale=g_scale*ABS(det)/REAL(ndim(0)*ndim(1)*ndim(2),dp)
#endif
    IF (PRESENT(local_bounds)) THEN
        DO i=0,2
            l_bounds(:,permut(i))=local_bounds(:,i)
        END DO
    ELSE
        l_bounds(1,:)=0
        l_bounds(2,:)=ndim-1
    END IF
    IF (PRESENT(local_shift)) THEN
        DO i=0,2
            l_shift(permut(i))=local_shift(i)
        END DO
    ELSE
        l_shift=0 ! l_bounds(1,:)
    END IF
    l_ub=l_bounds(2,:)-l_bounds(1,:)+l_shift
    DO i=1,3
        CPPrecondition(l_ub(permut(i-1))<SIZE(grid,i),cp_failure_level,routineP,error,failure)
    END DO
    IF (failure) GOTO 21
    
    DO i=0,2
        period(permut(i))=periodic(i)
    END DO
    CPPrecondition(ALL(l_bounds(2,:)<ndim.or.period(:)==1),cp_failure_level,routineP,error,failure)
    CPPrecondition(ALL(l_bounds(1,:)>=0 .or.period(:)==1),cp_failure_level,routineP,error,failure)
    CPPrecondition(ALL(l_bounds(2,:)-l_bounds(1,:)<ndim),cp_failure_level,routineP,error,failure)
    rPos=0.0_dp
    DO j=0,2
        DO i=0,2
            rPos(permut(i))=rPos(permut(i))+h_inv(i,j)*posi(j)
        END DO
    END DO
    cellShift=FLOOR(rPos)*period
    wrPos=rPos-REAL(cellShift,dp)
    riPos=wrPos*ndim
    shiftPos=FLOOR(riPos+0.5_dp)
    resPos=riPos-shiftPos
    normD=1.0_dp/REAL(ndim,dp)
    scaled_h=0.0_dp
    DO j=0,2
        DO i=0,2
            scaled_h(i,permut(j))=h(i,j)*normD(permut(j))
        END DO
    END DO
    resPosReal=0.0_dp
    DO j=0,2
        DO i=0,2
            resPosReal(i)=resPosReal(i)+h(i,j)*(wrPos(permut(j))-normD(permut(j))*REAL(shiftPos(permut(j)),dp))
        END DO
    END DO
        
    maxr2=0.0_dp
    DO j=0,2
        DO i=0,2
            maxr2=maxr2+(scaled_h(i,j))**2 ! guarantee at least the nearest points (this increases the sphere, increase just the box?)
        END DO
    END DO
    maxr2=MAX(max_r2,maxr2)
        
    ! build up quadratic form (ellipsoid)
    m=0.0_dp
    DO j=0,2
        DO i=0,2
            DO k=0,2
                m(i,j)=m(i,j)+scaled_h(k,i)*scaled_h(k,j)
            END DO
        END DO
    END DO
    
    l=0.0_dp
    DO j=0,2
        DO i=0,2
            l(j)=l(j)-2.0*resPos(i)*m(i,j)
        END DO
    END DO
    
    r_0=0.0_dp
    DO i=0,2
        r_0=r_0-0.5*resPos(i)*l(i)
    END DO
    
    ! calc i boundaries
    cci2 = (m(2,2) * m(0,0) * m(1,1) - m(2,2) * m(0,1) ** 2 - m(1,1) * m(0,2) ** 2 &
        + 2.0_dp * m(0,2) * m(0,1) * m(1,2) - m(0,0) * m(1,2) ** 2) / (m(2,2) * m(1,1) - m(1,2) ** 2)
    cci1 = -(-m(2,2) * l(0) * m(1,1) + m(2,2) * m(0,1) * l(1) + l(2) * m(0,2) * m(1,1) &
        + l(0) * m(1,2) ** 2 - l(2) * m(0,1) * m(1,2) - m(1,2) * l(1) * m(0,2)) / (m(2,2) * m(1,1) - m(1,2) ** 2)
    cci0 = -((-4.0_dp * m(2,2) * r_0 * m(1,1) + m(2,2) * l(1) ** 2 + l(2) ** 2 * m(1,1) &
        - 2.0_dp * l(1) * m(1,2) * l(2) + 4.0_dp * r_0 * m(1,2) ** 2) &
        / (m(2,2) * m(1,1) - m(1,2) ** 2)) / 4.0_dp-maxr2
    delta_i=cci1*cci1-4.0_dp*cci2*cci0
    IF (delta_i<=0) GOTO 21
    sqDi=SQRT(delta_i)
    imin=CEILING((-cci1-sqDi)/(2.0_dp*cci2))
    imax=FLOOR((-cci1+sqDi)/(2.0_dp*cci2))
    !! box early return
    
    IF (period(0)==1) THEN
        has_overlap=imax-imin+1>ndim(0).OR.(l_bounds(1,0)==0.and.l_bounds(2,0)==ndim(0)-1)
        IF (.not.has_overlap) THEN
            imin1=MODULO(imin+shiftPos(0),ndim(0))
            imax1=imin1+imax-imin+1
            IF (imin1<l_bounds(1,0)) THEN
                has_overlap=imax1>=l_bounds(1,0)
            ELSE
                has_overlap=imin1<=l_bounds(2,0).OR.(imax1>=ndim(0).and.l_bounds(1,0)<=imax1+ndim(0))
            END IF
            IF (.not.has_overlap) GOTO 21
        END IF
    ELSE
        IF (imax+shiftPos(0)<l_bounds(1,0).or.imin+shiftPos(0)>l_bounds(2,0)) GOTO 21
    END IF
    
    ! j box bounds
    has_overlap=l_bounds(1,1)==0.and.l_bounds(2,1)==ndim(1)-1
    IF (.not.has_overlap) THEN
        ccj2 = (m(0,0) * m(2,2) * m(1,1) - m(0,0) * m(1,2) ** 2 - m(0,1) ** 2 * m(2,2) &
            + 2.0_dp * m(0,1) * m(0,2) * m(1,2) - m(1,1) * m(0,2) ** 2) &
            / (m(0,0) * m(2,2) - m(0,2) ** 2)
        ccj1 = -(-m(0,0) * l(1) * m(2,2) + m(0,0) * l(2) * m(1,2) + l(0) * m(0,1) * m(2,2) &
            - m(0,2) * l(2) * m(0,1) - l(0) * m(0,2) * m(1,2) + l(1) * m(0,2) ** 2) &
            / (m(0,0) * m(2,2) - m(0,2) ** 2)
        ccj0 = (4.0_dp * m(0,0) * m(2,2) * r_0 - m(0,0) * l(2) ** 2 - m(2,2) * l(0) ** 2 &
            + 2.0_dp * m(0,2) * l(2) * l(0) - 4.0_dp * r_0 * m(0,2) ** 2) &
            / (m(0,0) * m(2,2) - m(0,2) ** 2) / 4.0_dp-maxr2
        delta_j=ccj1*ccj1-4.0_dp*ccj2*ccj0
        IF (delta_j<=0) GOTO 21
        sqDj=SQRT(delta_j)
        jmin=CEILING((-ccj1-sqDj)/(2.0_dp*ccj2))
        jmax=FLOOR((-ccj1+sqDj)/(2.0_dp*ccj2))
        IF (period(1)==1) THEN
            IF (jmax-jmin+1<ndim(1)) THEN
                jmin1=MODULO(jmin+shiftPos(1),ndim(1))
                jmax1=jmin1+jmax-jmin+1
                IF (jmin1<l_bounds(1,1)) THEN
                    has_overlap=jmax1>=l_bounds(1,1)
                ELSE
                    has_overlap=jmin1<=l_bounds(2,1).OR.(jmax1>=ndim(1).AND.(l_bounds(1,1)<=jmax1-ndim(1)))
                END IF
                IF (.not.has_overlap) GOTO 21
            END IF
        ELSE
            IF (jmax+shiftPos(1)<l_bounds(1,1).or.jmin+shiftPos(1)>l_bounds(2,1)) GOTO 21
        END IF
    END IF
    
    ! k box bounds
    has_overlap=l_bounds(1,2)==0.and.l_bounds(2,2)==ndim(2)-1
    IF (.not.has_overlap) THEN
        cck2 = (m(0,0) * m(2,2) * m(1,1) - m(0,0) * m(1,2) ** 2 - m(0,1) ** 2 * m(2,2) &
            + 2.0_dp * m(0,1) * m(0,2) * m(1,2) - m(1,1) * m(0,2) ** 2) / (m(0,0) * m(1,1) - m(0,1) ** 2)
        cck1 = (m(0,0) * l(2) * m(1,1) - m(0,0) * m(1,2) * l(1) - m(0,2) * l(0) * m(1,1) &
            + l(0) * m(0,1) * m(1,2) - l(2) * m(0,1) ** 2 + m(0,1) * l(1) * m(0,2)) / (m(0,0) * m(1,1) - m(0,1) ** 2)
        cck0 = (4.0_dp * m(0,0) * m(1,1) * r_0 - m(0,0) * l(1) ** 2 - m(1,1) * l(0) ** 2 &
            + 2.0_dp * m(0,1) * l(1) * l(0) - 4.0_dp * r_0 * m(0,1) ** 2) &
            / (m(0,0) * m(1,1) - m(0,1) ** 2) / 4.0_dp-maxr2
        delta_k=cck1*cck1-4.0_dp*cck2*cck0
        IF (delta_k<=0) GOTO 21
        sqDk=SQRT(delta_k)
        kmin=CEILING((-cck1-sqDk)/(2.0_dp*cck2))
        kmax=FLOOR((-cck1+sqDk)/(2.0_dp*cck2))
    
        IF (period(2)==1) THEN
            IF (kmax-kmin+1<ndim(2)) THEN
                kmin1=MODULO(kmin+shiftPos(2),ndim(2))
                kmax1=kmin1+kmax-kmin+1
                IF (kmin1<l_bounds(1,2)) THEN
                    has_overlap=kmax1>=l_bounds(1,2)
                ELSE
                    has_overlap=kmin1<=l_bounds(2,2).OR.&
                        (kmax1>=ndim(2).AND.(l_bounds(1,2)<=MODULO(kmax1,ndim(2))))
                END IF
                IF (.not.has_overlap) GOTO 21
            END IF
        ELSE
            IF (kmax+shiftPos(2)<l_bounds(1,2).or.kmin+shiftPos(2)>l_bounds(2,2)) GOTO 21
        END IF
    END IF
    
    ! k bounds (cache a la cube_info, or inversely integrate in the collocate loop?)
    ccj2   = (m(2,2) * m(1,1) - m(1,2) ** 2) / m(2,2)
    ccj1_i1=(2 * m(2,2) * m(0,1) - 2 * m(0,2) * m(1,2)) / m(2,2)
    ccj1_i0=(-l(2) * m(1,2) + m(2,2) * l(1)) / m(2,2)
    ccj0_i2=(m(2,2) * m(0,0)-m(0,2) ** 2) / m(2,2)
    ccj0_i1=( m(2,2) * l(0) - m(0,2) * l(2) ) / m(2,2)
    ccj0_i0=(m(2,2) * r_0 - 0.25 * l(2) ** 2) / m(2,2) - maxr2
    cck2   = m(2,2)
    cck1_i = 2 * m(0,2)
    cck1_j = 2 * m(1,2)
    cck1_0 = l(2)
    cck0_i2 = m(0,0)
    cck0_ij = 2 * m(0,1)
    cck0_i = l(0)
    cck0_j2 = m(1,1)
    cck0_j = l(1)
    cck0_0 = r_0 - maxr2
    
    ! find maximum number of j
    max_j=0
    DO i0=0,1
        i=(imin+imax)/2+i0
        ccj1 = ccj1_i1 * i +ccj1_i0
        ccj0 = (ccj0_i2*i+ccj0_i1)*i+ccj0_i0
        delta_j=ccj1*ccj1-4*ccj2*ccj0
        IF (delta_j>=0) THEN
            sqDj=SQRT(delta_j)
            max_j=MAX(max_j,FLOOR((-ccj1+sqDj)/(2.0_dp*ccj2)) &
                        -CEILING((-ccj1-sqDj)/(2.0_dp*ccj2))+1)
        END IF
    END DO
    max_j=max_j+1 ! just to be sure...
    IF (period(1)==0) max_j=MIN(max_j,l_bounds(2,1)-l_bounds(1,1)+1)
    
    IF (period(0)==0) THEN
        imin=MAX(l_bounds(1,0)-shiftPos(0),imin)
        imax=MIN(l_bounds(2,0)-shiftPos(0),imax)
    END IF
    
    ! k bounds (cache a la cube_info?)
    has_overlap=.FALSE.
    ALLOCATE(k_bounds(0:1,0:max_j-1,0:MAX(0,imax-imin+1)),stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
    k_bounds_alloc=.TRUE.
    k_bounds=0
    istart=imin
    iiShift=shiftPos(0)-l_bounds(2,0)+istart
    IF (iiShift>0) iiShift=iiShift+ndim(0)-1
    iiShift=(iiShift/ndim(0))*ndim(0)-shiftPos(0)
    !iiShift=CEILING(REAL(shiftPos(0)+istart-l_bounds(2,0))/REAL(ndim(0)))*ndim(0)-shiftPos(0))
    istart=MAX(iiShift+l_bounds(1,0),istart)
    iend=MIN(iiShift+l_bounds(2,0),imax)
    iJump=ndim(0)-l_bounds(2,0)+l_bounds(1,0)-1
    jJump=ndim(1)-l_bounds(2,1)+l_bounds(1,1)-1
    DO
        DO i=istart,iend
            ccj1 = ccj1_i1 * i +ccj1_i0
            ccj0 = (ccj0_i2*i+ccj0_i1)*i+ccj0_i0
            delta_j=ccj1*ccj1-4*ccj2*ccj0
            IF (delta_j<0) CONTINUE
            sqDj=SQRT(delta_j)
            jmin=CEILING((-ccj1-sqDj)/(2.0_dp*ccj2))
            jmax=FLOOR((-ccj1+sqDj)/(2.0_dp*ccj2))
            cck0_0_p=cck0_0+(cck0_i2*i+cck0_i)*i
            cck0_j_p= cck0_j+cck0_ij*i
            cck1_0_p=cck1_0+cck1_i*i
            IF (period(1)==0) THEN
                jmin=MAX(l_bounds(1,1)-shiftPos(1),jmin)
                jmax=MIN(l_bounds(2,1)-shiftPos(1),jmax)
            END IF
            jstart=jmin
            ijShift=shiftPos(1)+jstart-l_bounds(2,1)
            IF (ijShift>0) ijShift=ijShift+ndim(1)-1
            ijShift=(ijShift/ndim(1))*ndim(1)-shiftPos(1)
            ! ijShift=CEILING(REAL(shiftPos(1)+jstart-l_bounds(2,1))/REAL(ndim(1)))*ndim(1)-shiftPos(1)
            jstart=MAX(ijShift+l_bounds(1,1),jstart)
            jend=MIN(ijShift+l_bounds(2,1),jmax)
            DO
                DO j=jstart,jend
                    cck1=cck1_0_p+cck1_j*j
                    cck0=cck0_0_p+(cck0_j_p+cck0_j2*j)*j
        
                    delta_k=cck1*cck1-4*cck2*cck0
                    IF (delta_k<0) THEN
                        k_bounds(0,j-jmin,i-imin)=0 ! CEILING((-cck1)/(2.0_dp*cck2))
                        k_bounds(1,j-jmin,i-imin)=-1 ! k_bounds(0,j-jmin,i-imin)-1
                    ELSE
                        sqDk=SQRT(delta_k)
                        kmin=CEILING((-cck1-sqDk)/(2.0_dp*cck2))
                        kmax=FLOOR((-cck1+sqDk)/(2.0_dp*cck2))
                        
                        ! ! reduce kmax,kmin
                        ! ! this should be done later if k_bounds are shared by threads with different slices
                        ! ikShift=FLOOR(REAL(shiftPos(2)+kmax-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
                        ! kmax=MIN(kmax,ikShift+l_bounds(2,2))
                        ! ikShift2=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kmin)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
                        ! kmin=MAX(kmin,ikShift2+l_bounds(1,2))
                        
                        k_bounds(0,j-jmin,i-imin)=kmin
                        k_bounds(1,j-jmin,i-imin)=kmax
                        IF (kmax>=kmin) has_overlap=.TRUE.
                    END IF
                END DO
                jstart=jend+jJump+1
                IF (jstart>jmax) EXIT
                jend=MIN(jend+ndim(1),jmax)
            END DO
        END DO
        istart=iend+iJump+1
        IF (istart>imax) EXIT
        iend=MIN(iend+ndim(0),imax)
    END DO
    IF (period(2)==0) THEN
        k_bounds(0,:,:)=MAX(l_bounds(1,2)-shiftPos(2),k_bounds(0,:,:))
        k_bounds(1,:,:)=MIN(l_bounds(2,2)-shiftPos(2),k_bounds(1,:,:))
    END IF
    
    IF (.not.has_overlap) GOTO 21
    
    ! poly x,y,z -> i,j,k
    grad=grad_size3(SIZE(poly)/npoly)
    size_jk=poly_size2(grad)*npoly
    size_k=poly_size1(grad)*npoly
    size_ijk=poly_size3(grad)*npoly
    ALLOCATE(poly_ijk(size_ijk),&
        poly_jk(size_jk),&
        xi(grad+1),stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
#ifdef FMG_INTEGRATE
    ALLOCATE(poly_k(0:size_k-1),mres(npoly),k_vals(0:grad),stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
    mres=0.0_dp
#elif defined(FMG_INTEGRATE_FULL)
    ALLOCATE(k_vals(0:grad),stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
#else
    ALLOCATE(poly_k(0:size_k-1),stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
#endif
poly_alloc=.TRUE.
#ifdef FMG_INTEGRATE_FULL
    CPPreconditionNoFail(SIZE(poly)==poly_size3(grad),cp_failure_level,routineP,error)
    poly_ijk=0.0_dp
#else
    CALL poly_affine_t3(poly,scaled_h,-resPosReal+p_shift,poly_ijk,&
        npoly=npoly,error=error)
#endif
    
    ij_coeff0=EXP(-2.0_dp*alphai*m(0,1))
    ik_coeff0=EXP(-2.0_dp*alphai*m(0,2))
    ii_coeff0=EXP(-alphai*m(0,0))
    jk_coeff0=EXP(-2.0_dp*alphai*m(1,2))
    jj_coeff0=EXP(-alphai*m(1,1))
    kk_coeff0=EXP(-alphai*m(2,2))
    jk_coeff1=1.0_dp/jk_coeff0
    ij_coeff1=1.0_dp/ij_coeff0
    ik_coeff1=1.0_dp/ik_coeff0
    ii_coeff2=ii_coeff0*ii_coeff0
    jj_coeff2=jj_coeff0*jj_coeff0
    kk_coeff2=kk_coeff0*kk_coeff0
    icoeff0=EXP(-alphai*l(0))
    jcoeff0=EXP(-alphai*l(1))
    kcoeff0=EXP(-alphai*l(2))
    res_0=EXP(-alphai*r_0)*g_scale
    
    i_coeffn_i=icoeff0
    j_coeffn_i=jcoeff0
    k_coeffn_i=kcoeff0
    ii_coeffn =i_coeffn_i*ii_coeff0
    res_i=res_0
    
    iJump=ndim(0)-l_bounds(2,0)+l_bounds(1,0)-1
    istart=MAX(0,imin)
    iiShift=shiftPos(0)-l_bounds(2,0)+istart
    IF (iiShift>0) iiShift=iiShift+ndim(0)-1
    iiShift=(iiShift/ndim(0))*ndim(0)-shiftPos(0)
    !iiShift=CEILING(REAL(shiftPos(0)+istart-l_bounds(2,0))/REAL(ndim(0)))*ndim(0)-shiftPos(0)
    istart=MAX(iiShift+l_bounds(1,0),istart)
    iistart=istart-iiShift-l_bounds(1,0)+l_shift(0)
    istart2=MIN(-1,imax)
    iiShift2=shiftPos(0)+istart2-l_bounds(1,0)
    IF (iiShift2<0) iiShift2=iiShift2-ndim(0)+1
    iiShift2=(iiShift2/ndim(0))*ndim(0)-shiftPos(0)
    !iiShift2=FLOOR(REAL(shiftPos(0)+istart2-l_bounds(1,0))/REAL(ndim(0)))*ndim(0)-shiftPos(0)
    istart2=MIN(iiShift2+l_bounds(2,0),istart2)
    iistart2=istart2-iiShift2-l_bounds(1,0)+l_shift(0)

    IF (iJump/=0.and.(iistart+imax-istart>=ndim(0)+l_shift(0) .OR.&
        iistart2+imin-istart2<=l_ub(0)-ndim(0))) THEN
        ! will wrap
        ij_coeff0_jump =ij_coeff0**(iJump)
        ik_coeff0_jump =ik_coeff0**(iJump)
        ii_coeff2_jump =ii_coeff2**(iJump)
        ii_coeffn_jump =ii_coeff0**((iJump)*(iJump-1))
    ELSE
        ij_coeff0_jump = 1.0_dp
        ik_coeff0_jump = 1.0_dp
        ii_coeff2_jump = 1.0_dp
        ii_coeffn_jump = 1.0_dp
    END IF
    
    iend=MIN(iiShift+l_bounds(2,0),imax)
    ii=iistart
    IF (i>0) THEN
        ii_coeffn=i_coeffn_i*ii_coeff0**(2*istart+1)
        j_coeffn_i=jcoeff0*ij_coeff0**istart
        k_coeffn_i=kcoeff0*ik_coeff0**istart
        res_i=res_0*(ii_coeff0**istart*i_coeffn_i)**istart
    END IF
    DO
        DO i=istart,iend
            ! perform j loop
            IF (ABS(res_i)>small) THEN
                CALL j_loop
            END IF
            j_coeffn_i=j_coeffn_i*ij_coeff0
            k_coeffn_i=k_coeffn_i*ik_coeff0
            res_i=res_i*ii_coeffn
            ii_coeffn=ii_coeffn*ii_coeff2
            ii=ii+1
        END DO
        istart=iend+iJump+1
        IF (istart>imax) EXIT
        iend=MIN(iend+ndim(0),imax)
        ii=l_shift(0)
        j_coeffn_i=j_coeffn_i*ij_coeff0_jump
        k_coeffn_i=k_coeffn_i*ik_coeff0_jump
        res_i=res_i*ii_coeffn**(iJump)*ii_coeffn_jump
        ii_coeffn=ii_coeffn*ii_coeff2_jump
    END DO
    
    ! neg i side
    i_coeffn_i=1.0_dp/icoeff0
    j_coeffn_i=jcoeff0
    k_coeffn_i=kcoeff0
    res_i=res_0
    ii_coeffn=i_coeffn_i*ii_coeff0
    
    iend2=MAX(iiShift2+l_bounds(1,0),imin)
    ii=iistart2
    IF (istart2<-1) THEN
        ii_coeffn=i_coeffn_i*ii_coeff0**(-(2*istart2+1))
        j_coeffn_i=jcoeff0*ij_coeff0**(istart2+1)
        k_coeffn_i=kcoeff0*ik_coeff0**(istart2+1)
        res_i=res_0*(ii_coeff0**(-istart2-1)*i_coeffn_i)**(-istart2-1)
    END IF
    DO
        DO i=istart2,iend2,-1
            j_coeffn_i=j_coeffn_i*ij_coeff1
            k_coeffn_i=k_coeffn_i*ik_coeff1
            res_i=res_i*ii_coeffn
            ii_coeffn=ii_coeffn*ii_coeff2

            ! perform j loop
            IF (ABS(res_i)>small) THEN
                CALL j_loop
            END IF
            ii=ii-1
        END DO
        istart2=iend2-iJump-1
        IF (istart2<imin) EXIT
        iend2=MAX(iend2-ndim(0),imin)
        ii=l_ub(0)
        j_coeffn_i=j_coeffn_i/ij_coeff0_jump
        k_coeffn_i=k_coeffn_i/ik_coeff0_jump
        res_i=res_i*ii_coeffn**iJump*ii_coeffn_jump
        ii_coeffn=ii_coeffn*ii_coeff2_jump
    END DO

    ! the final cleanup
21  CONTINUE

    ! dealloc to avoid compiler bug
    IF (k_bounds_alloc) THEN
        DEALLOCATE(k_bounds,stat=stat)
        CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
    END IF        
    IF (poly_alloc) THEN
#ifdef FMG_INTEGRATE
        DEALLOCATE(poly_ijk,poly_jk,poly_k,stat=stat)
        CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
        ! set result
        res=mres
        DEALLOCATE(mres,k_vals,stat=stat)
        CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
    ELSE
        res=0.0_dp
#elif defined(FMG_INTEGRATE_FULL)
        CALL poly_affine_t3t(poly_ijk,scaled_h,-resPosReal+p_shift,poly,&
            npoly=npoly,error=error)
        DEALLOCATE(k_vals,stat=stat)
        CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
        DEALLOCATE(poly_ijk,poly_jk,stat=stat)
        CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
    ELSE
        poly=0.0_dp
#else
        DEALLOCATE(poly_ijk,poly_jk,poly_k,stat=stat)
        CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
#endif
    END IF
    CONTAINS
    
!!!!!!!!!!!!!
    !!!!!!! j loop 
SUBROUTINE j_loop
    ! calculate j bounds
    ccj1 = ccj1_i1 * i +ccj1_i0
    ccj0 = (ccj0_i2*i+ccj0_i1)*i+ccj0_i0
    delta_j=ccj1*ccj1-4*ccj2*ccj0
    IF (delta_j<0) THEN
        RETURN
    END IF
    sqDj=SQRT(delta_j)
    jmin=CEILING((-ccj1-sqDj)/(2.0_dp*ccj2))
    jmax=FLOOR((-ccj1+sqDj)/(2.0_dp*ccj2))
    
#ifdef FMG_INTEGRATE_FULL
    poly_jk=0.0_dp
#else
    CALL poly_p_eval3b(poly_ijk,size_ijk,REAL(i,dp),poly_jk,size_jk,&
            npoly=npoly,grad=grad,xi=xi)
#endif

    IF (period(1)==0) THEN 
        jmin=MAX(l_bounds(1,1)-shiftPos(1),jmin)
        jmax=MIN(l_bounds(2,1)-shiftPos(1),jmax)
    END IF
    
    ! pos j side
    j_coeffn_j=j_coeffn_i
    k_coeffn_j=k_coeffn_i
    jj_coeffn=j_coeffn_j*jj_coeff0
    res_j=res_i

    jJump=ndim(1)-l_bounds(2,1)+l_bounds(1,1)
    jstart=MAX(0,jmin)
    ijShift=shiftPos(1)+jstart-l_bounds(2,1)
    IF (ijShift>0) ijShift=ijShift+ndim(1)-1
    ijShift=(ijShift/ndim(1))*ndim(1)-shiftPos(1)
    !ijShift=CEILING(REAL(shiftPos(1)+jstart-l_bounds(2,1))/REAL(ndim(1)))*ndim(1)-shiftPos(1)
    jstart=MAX(ijShift+l_bounds(1,1),jstart)
    jend=MIN(ijShift+l_bounds(2,1),jmax)
    ij=jstart-ijShift-l_bounds(1,1)+l_shift(1)

    IF (jstart>0) THEN
        k_coeffn_j=k_coeffn_i*jk_coeff0**jstart
        jj_coeffn=j_coeffn_j*jj_coeff0**(2*jstart+1)
        res_j=res_i*(jj_coeff0**jstart*j_coeffn_j)**jstart
    END IF
    DO
        DO j=jstart,jend
            kmin=k_bounds(0,j-jmin,i-imin)
            kmax=k_bounds(1,j-jmin,i-imin)
            ! do k loop
            IF (res_j/=0.0_dp.and.k_coeffn_j/=0.0_dp.and.kmin<=kmax &
                .and.ABS(res_j)>small) THEN
                CALL k_loop
            END IF
            k_coeffn_j=k_coeffn_j*jk_coeff0
            res_j=res_j*jj_coeffn
            jj_coeffn=jj_coeffn*jj_coeff2
            ij=ij+1
        END DO
        jstart=jend+jJump
        IF (jstart>jmax) EXIT
        ij=l_shift(1)
        jend=MIN(jend+ndim(1),jmax)
        IF (jJump/=1) THEN ! remove if?
            k_coeffn_j=k_coeffn_i*jk_coeff0**jstart
            jj_coeffn=j_coeffn_j*jj_coeff0**(2*jstart+1)
            res_j=res_i*(jj_coeff0**jstart*j_coeffn_j)**jstart
        END IF
    END DO

    ! neg j side
    j_coeffn_j=1.0_dp/j_coeffn_i
    k_coeffn_j=k_coeffn_i
    jj_coeffn=j_coeffn_j*jj_coeff0
    res_j=res_i
    
    jstart=MIN(-1,jmax)
    ijShift=shiftPos(1)+jstart-l_bounds(1,1)
    IF (ijShift<0) ijShift=ijShift-ndim(1)+1
    ijShift=(ijShift/ndim(1))*ndim(1)-shiftPos(1)
    !ijShift=FLOOR(REAL(shiftPos(1)+jstart-l_bounds(1,1))/REAL(ndim(1)))*ndim(1)-shiftPos(1))
    jstart=MIN(ijShift+l_bounds(2,1),jstart)
    jend=MAX(ijShift+l_bounds(1,1),jmin)
    ij=jstart-ijShift-l_bounds(1,1)+l_shift(1)
    IF (jstart<-1) THEN
        k_coeffn_j=k_coeffn_i*jk_coeff0**(jstart+1)
        jj_coeffn=j_coeffn_j*jj_coeff0**(-(2*jstart+1))
        res_j=res_i*(jj_coeff0**(-jstart-1)*j_coeffn_j)**(-jstart-1)
    END IF
    DO
        DO j=jstart,jend,-1
            k_coeffn_j=k_coeffn_j*jk_coeff1
            res_j=res_j*jj_coeffn
            jj_coeffn=jj_coeffn*jj_coeff2
        
            kmin=k_bounds(0,j-jmin,i-imin)
            kmax=k_bounds(1,j-jmin,i-imin)
            ! perform k loop
            IF (res_j/=0.0_dp.and.k_coeffn_j/=0.0_dp.and.kmin<=kmax &
                .and.ABS(res_j)>small) THEN
                CALL k_loop
            END IF
            ij=ij-1
        END DO
        jstart=jend-jJump
        IF (jstart<jmin) EXIT
        jend=MAX(jend-ndim(1),jmin)
        ij=l_ub(1)
        IF (jJump/=1) THEN ! remove if?
            k_coeffn_j=k_coeffn_i*jk_coeff0**(jstart+1)
            jj_coeffn=j_coeffn_j*jj_coeff0**(-(2*jstart+1))
            res_j=res_i*(jj_coeff0**(-jstart-1)*j_coeffn_j)**(-jstart-1)
        END IF
    END DO
#ifdef FMG_INTEGRATE_FULL
    CALL poly_padd_uneval3b(poly_ijk,size_ijk,REAL(i,dp),poly_jk,size_jk,&
        npoly=npoly,grad=grad,xi=xi)
    !CALL poly_padd_uneval3(poly_ijk,REAL(i,dp),poly_jk,npoly=npoly,error=error)
#endif
END SUBROUTINE

!!!!!!!!!!!!
    !!!!!!! k loop
SUBROUTINE k_loop()
#ifndef FMG_INTEGRATE_FULL
    CALL poly_p_eval2b(poly_jk,size_jk,REAL(j,dp),poly_k,&
        size_k,npoly=npoly,grad=grad,xi=xi)
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
                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                DO kgrad=1,grad
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(kgrad)=k_vals(kgrad)+p_kk
                END DO
#else
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                DO kgrad=1,grad
                    p_v=p_v+poly_k(kgrad)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                END DO
                grid(ik,ij,ii)=grid(ik,ij,ii)+p_v*res_k
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
                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                DO kgrad=1,grad
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(kgrad)=k_vals(kgrad)+p_kk
                END DO
#else
                p_v=poly_k(0)
                p_kk=k
                DO kgrad=1,grad
                    p_v=p_v+poly_k(kgrad)*p_kk
                    p_kk=p_kk*k
                END DO
                grid(ik,ij,ii)=grid(ik,ij,ii)+p_v*res_k
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
                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                DO kgrad=1,grad
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(kgrad)=k_vals(kgrad)+p_kk
                END DO
#else
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                DO kgrad=1,grad
                    p_v=p_v+poly_k(kgrad)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                END DO
                grid(ik,ij,ii)=grid(ik,ij,ii)+p_v*res_k
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
                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                DO kgrad=1,grad
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(kgrad)=k_vals(kgrad)+p_kk
                END DO
#else
                p_v=poly_k(0)
                p_kk=k
                DO kgrad=1,grad
                    p_v=p_v+poly_k(kgrad)*p_kk
                    p_kk=p_kk*k
                END DO
                grid(ik,ij,ii)=grid(ik,ij,ii)+p_v*res_k
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
        DO kgrad=0,grad
            p_v=p_v+poly_k(pShift+kgrad)*k_vals(kgrad)
        END DO
        mres(ipoly)=mres(ipoly)+p_v
        pShift=pShift+grad+1
    END DO
#elif defined(FMG_INTEGRATE_FULL)
    CALL poly_padd_uneval2b(poly_jk,size_jk,REAL(j,dp),k_vals,&
        size_k,npoly=npoly,grad=grad,xi=xi)
    !CALL poly_padd_uneval2(poly_jk,REAL(j,dp),k_vals,npoly=npoly,error=error)
#endif
END SUBROUTINE

