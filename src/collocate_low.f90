    !
    ! this treats efficiently the orthogonal case
    !
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE collocate_ortho()


!   *** properties of the grid ***

    ! notice we're in the ortho case
    dr(1) = rsgrid%desc%dh(1,1)
    dr(2) = rsgrid%desc%dh(2,2)
    dr(3) = rsgrid%desc%dh(3,3)

!   *** get the sub grid properties for the given radius ***
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)
    cmax=MAXVAL(ub_cube)

!   *** position of the gaussian product
!
!   this is the actual definition of the position on the grid
!   i.e. a point rp(:) gets here grid coordinates
!   MODULO(rp(:)/dr(:),ng(:))+1
!   hence (0.0,0.0,0.0) in real space is rsgrid%lb on the rsgrid ((1,1,1) on grid)
!



    ALLOCATE(map(-cmax:cmax,3))
    CALL compute_cube_center(cubecenter,rsgrid%desc,zeta,zetb,ra,rab)
    roffset(:)    = rp(:) - REAL(cubecenter(:),dp)*dr(:)
!   *** a mapping so that the ig corresponds to the right grid point
    DO i=1,3
      IF ( rsgrid % desc % perd ( i ) == 1 ) THEN
        start=lb_cube(i)
        DO
         offset=MODULO(cubecenter(i)+start,ng(i))+1-start
         length=MIN(ub_cube(i),ng(i)-offset)-start
         DO ig=start,start+length
            map(ig,i) = ig+offset
         END DO
         IF (start+length.GE.ub_cube(i)) EXIT
         start=start+length+1
        END DO
      ELSE
        ! this takes partial grid + border regions into account
        offset=MODULO(cubecenter(i)+lb_cube(i)+rsgrid%desc%lb(i)-rsgrid%lb_local(i),ng(i))+1-lb_cube(i)
        ! check for out of bounds
        IF (ub_cube(i)+offset>UBOUND(grid,i).OR.lb_cube(i)+offset<LBOUND(grid,i)) THEN
           CPASSERT(.FALSE.)
        ENDIF
        DO ig=lb_cube(i),ub_cube(i)
           map(ig,i) = ig+offset
        END DO
      END IF
    ENDDO
    ALLOCATE(pol_z(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_y(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_x(0:lp,-cmax:cmax))

    IF(PRESENT(ir).AND.PRESENT(rsgauge)) CALL collocate_ortho_set_to_0()

#include "prep.f90"

    IF ( PRESENT ( lgrid ) ) THEN
       ig = lgrid%ldim * ithread_l + 1
#include "call_collocate_omp.f90"
    ELSE
#include "call_collocate.f90"
    END IF

    IF(PRESENT(ir).AND.PRESENT(rsgauge)) CALL collocate_gauge_ortho()

    ! deallocation needed to pass around a pgi bug..
    DEALLOCATE(pol_z)
    DEALLOCATE(pol_y)
    DEALLOCATE(pol_x)
    DEALLOCATE(map)

    END SUBROUTINE collocate_ortho

! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE collocate_gauge_ortho()
    INTEGER                                  :: i, igmax, igmin, j, j2, jg, &
                                                jg2, jgmin, k, k2, kg, kg2, &
                                                kgmin, sci
    REAL(KIND=dp)                            :: point(3,4), res(4), x, y, y2, &
                                                z, z2

! notice we're in the ortho case

      dr(1) = rsgrid%desc%dh(1,1)
      dr(2) = rsgrid%desc%dh(2,2)
      dr(3) = rsgrid%desc%dh(3,3)
      !
      sci=1
      kgmin=sphere_bounds(sci)
      sci=sci+1
      DO kg=kgmin,0
         kg2=1-kg
         k=map(kg,3)
         k2=map(kg2,3)
         jgmin=sphere_bounds(sci)
         sci=sci+1
         z  = (REAL( kg,dp)+REAL(cubecenter(3),dp))*dr(3)
         z2 = (REAL(kg2,dp)+REAL(cubecenter(3),dp))*dr(3)
         DO jg=jgmin,0
            jg2=1-jg
            j=map(jg,2)
            j2=map(jg2,2)
            igmin=sphere_bounds(sci)
            sci=sci+1
            igmax=1-igmin
            y  = (REAL( jg,dp)+REAL(cubecenter(2),dp))*dr(2)
            y2 = (REAL(jg2,dp)+REAL(cubecenter(2),dp))*dr(2)
            DO ig=igmin,igmax
               i=map(ig,1)
               x = (REAL(ig,dp)+REAL(cubecenter(1),dp))*dr(1)
               point(1,1) = x;point(2,1) = y ;point(3,1) = z
               point(1,2) = x;point(2,2) = y2;point(3,2) = z
               point(1,3) = x;point(2,3) = y ;point(3,3) = z2
               point(1,4) = x;point(2,4) = y2;point(3,4) = z2
               !
               res(1) = ( point(ir,1) - rb(ir) ) - gauge(i, j, k)
               res(2) = ( point(ir,2) - rb(ir) ) - gauge(i,j2, k)
               res(3) = ( point(ir,3) - rb(ir) ) - gauge(i, j,k2)
               res(4) = ( point(ir,4) - rb(ir) ) - gauge(i,j2,k2)
               !
               grid_tmp(i, j, k) = grid_tmp(i, j, k) + grid(i, j, k) * res(1)
               grid_tmp(i,j2, k) = grid_tmp(i,j2, k) + grid(i,j2, k) * res(2)
               grid_tmp(i, j,k2) = grid_tmp(i, j,k2) + grid(i, j,k2) * res(3)
               grid_tmp(i,j2,k2) = grid_tmp(i,j2,k2) + grid(i,j2,k2) * res(4)
            ENDDO
         ENDDO
      ENDDO
    END SUBROUTINE collocate_gauge_ortho

! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE collocate_ortho_set_to_0()
    INTEGER                                  :: i, igmax, igmin, j, j2, jg, &
                                                jg2, jgmin, k, k2, kg, kg2, &
                                                kgmin, sci

!

      dr(1) = rsgrid%desc%dh(1,1)
      dr(2) = rsgrid%desc%dh(2,2)
      dr(3) = rsgrid%desc%dh(3,3)
      !
      sci=1
      kgmin=sphere_bounds(sci)
      sci=sci+1
      DO kg=kgmin,0
         kg2=1-kg
         k=map(kg,3)
         k2=map(kg2,3)
         jgmin=sphere_bounds(sci)
         sci=sci+1
         DO jg=jgmin,0
            jg2=1-jg
            j=map(jg,2)
            j2=map(jg2,2)
            igmin=sphere_bounds(sci)
            sci=sci+1
            igmax=1-igmin
            DO ig=igmin,igmax
               i=map(ig,1)
               grid(i, j, k) = 0.0_dp
               grid(i,j2, k) = 0.0_dp
               grid(i, j,k2) = 0.0_dp
               grid(i,j2,k2) = 0.0_dp
            ENDDO
         ENDDO
      ENDDO
    END SUBROUTINE collocate_ortho_set_to_0

!
!   this is a general 'optimized' routine to do the collocation
!
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE collocate_general_opt()

    INTEGER :: i, i_index, il, ilx, ily, ilz, index_max(3), index_min(3), &
      ismax, ismin, j, j_index, jl, jlx, jly, jlz, k, k_index, kl, klx, kly, &
      klz, lpx, lpy, lpz, lx, ly, lz, offset(3)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: grid_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: coef_map
    REAL(KIND=dp) :: a, b, c, d, di, dip, dj, djp, dk, dkp, exp0i, exp1i, &
      exp2i, gp(3), hmatgrid(3,3), pointj(3), pointk(3), res, v(3)
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: coef_ijk
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: hmatgridp

!
! transform P_{lxp,lyp,lzp} into a P_{lip,ljp,lkp} such that
! sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-x_p)**lxp (y-y_p)**lyp (z-z_p)**lzp =
! sum_{lip,ljp,lkp} P_{lip,ljp,lkp} (i-i_p)**lip (j-j_p)**ljp (k-k_p)**lkp
!

      ALLOCATE(coef_ijk(((lp+1)*(lp+2)*(lp+3))/6))

      ! aux mapping array to simplify life
      ALLOCATE(coef_map(0:lp,0:lp,0:lp))
      coef_map=HUGE(coef_map)
      lxyz=0
      DO lzp=0,lp
      DO lyp=0,lp-lzp
      DO lxp=0,lp-lzp-lyp
          lxyz=lxyz+1
          coef_ijk(lxyz)=0.0_dp
          coef_map(lxp,lyp,lzp)=lxyz
      ENDDO
      ENDDO
      ENDDO

      ! cell hmat in grid points
      hmatgrid=rsgrid%desc%dh

      ! center in grid coords
      gp=MATMUL(rsgrid%desc%dh_inv,rp)
      cubecenter(:) = FLOOR(gp)

      ! transform using multinomials
      ALLOCATE(hmatgridp(3,3,0:lp))
      hmatgridp(:,:,0)=1.0_dp
      DO k=1,lp
         hmatgridp(:,:,k)=hmatgridp(:,:,k-1)*hmatgrid(:,:)
      ENDDO

      lpx=lp
      DO klx=0,lpx
      DO jlx=0,lpx-klx
      DO ilx=0,lpx-klx-jlx
         lx=ilx+jlx+klx
         lpy=lp-lx
         DO kly=0,lpy
         DO jly=0,lpy-kly
         DO ily=0,lpy-kly-jly
            ly=ily+jly+kly
            lpz=lp-lx-ly
            DO klz=0,lpz
            DO jlz=0,lpz-klz
            DO ilz=0,lpz-klz-jlz
               lz=ilz+jlz+klz

               il=ilx+ily+ilz
               jl=jlx+jly+jlz
               kl=klx+kly+klz
               coef_ijk(coef_map(il,jl,kl))=coef_ijk(coef_map(il,jl,kl))+ coef_xyz(coef_map(lx,ly,lz))* &
                                            hmatgridp(1,1,ilx) * hmatgridp(1,2,jlx) * hmatgridp(1,3,klx) * &
                                            hmatgridp(2,1,ily) * hmatgridp(2,2,jly) * hmatgridp(2,3,kly) * &
                                            hmatgridp(3,1,ilz) * hmatgridp(3,2,jlz) * hmatgridp(3,3,klz) * &
                                            fac(lx)*fac(ly)*fac(lz)/ &
                        (fac(ilx)*fac(ily)*fac(ilz)*fac(jlx)*fac(jly)*fac(jlz)*fac(klx)*fac(kly)*fac(klz))
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL return_cube_nonortho(cube_info,radius,index_min,index_max,rp)

      offset(:)=MODULO(index_min(:)+rsgrid%desc%lb(:)-rsgrid%lb_local(:),ng(:))+1

      ALLOCATE(grid_map(index_min(1):index_max(1)))
      DO i=index_min(1),index_max(1)
        grid_map(i)=MODULO(i,ng(1))+1
        IF (rsgrid % desc % perd ( 1 )==1) THEN
           grid_map(i)=MODULO(i,ng(1))+1
        ELSE
           grid_map(i)=i-index_min(1)+offset(1)
        ENDIF
      ENDDO

      ! go over the grid, but cycle if the point is not within the radius
      DO k=index_min(3),index_max(3)
        dk=k-gp(3)
        pointk=hmatgrid(:,3)*dk

        IF (rsgrid % desc % perd ( 3 )==1) THEN
           k_index=MODULO(k,ng(3))+1
        ELSE
           k_index=k-index_min(3)+offset(3)
        ENDIF

        coef_xyt=0.0_dp
        lxyz = 0
        dkp=1.0_dp
        DO kl=0,lp
           lxy=0
           DO jl=0,lp-kl
              DO il=0,lp-kl-jl
                 lxyz=lxyz+1 ; lxy=lxy+1
                 coef_xyt(lxy)=coef_xyt(lxy)+coef_ijk(lxyz)*dkp
              ENDDO
              lxy=lxy+kl
           ENDDO
           dkp=dkp*dk
        ENDDO

        DO j=index_min(2),index_max(2)
          dj=j-gp(2)
          pointj=pointk+hmatgrid(:,2)*dj
          IF (rsgrid % desc % perd ( 2 )==1) THEN
             j_index=MODULO(j,ng(2))+1
          ELSE
             j_index=j-index_min(2)+offset(2)
          ENDIF

          coef_xtt=0.0_dp
          lxy=0
          djp=1.0_dp
          DO jl=0,lp
            DO il=0,lp-jl
               lxy=lxy+1
               coef_xtt(il)=coef_xtt(il)+coef_xyt(lxy)*djp
            ENDDO
            djp=djp*dj
          ENDDO

          ! find bounds for the inner loop
          ! based on a quadratic equation in i
          ! a*i**2+b*i+c=radius**2
          v=pointj-gp(1)*hmatgrid(:,1)
          a=DOT_PRODUCT(hmatgrid(:,1),hmatgrid(:,1))
          b=2*DOT_PRODUCT(v,hmatgrid(:,1))
          c=DOT_PRODUCT(v,v)
          d=b*b-4*a*(c-radius**2)

          IF (d<0) THEN
              CYCLE
          ELSE
              d=SQRT(d)
              ismin=CEILING((-b-d)/(2*a))
              ismax=FLOOR((-b+d)/(2*a))
          ENDIF
          ! prepare for computing -zetp*rsq
          a=-zetp*a
          b=-zetp*b
          c=-zetp*c
          i=ismin-1

          ! the recursion relation might have to be done
          ! from the center of the gaussian (in both directions)
          ! instead as the current implementation from an edge
          exp2i=EXP((a*i+b)*i+c)
          exp1i=EXP(2*a*i+a+b)
          exp0i=EXP(2*a)

          DO i=ismin,ismax
             di=i-gp(1)

             ! polynomial terms
             res=0.0_dp
             dip=1.0_dp
             DO il=0,lp
                res=res+coef_xtt(il)*dip
                dip=dip*di
             ENDDO

             ! the exponential recursion
             exp2i=exp2i*exp1i
             exp1i=exp1i*exp0i
             res=res*exp2i

             i_index=grid_map(i)
             IF ( PRESENT ( lgrid ) ) THEN
                ig = lgrid%ldim * ithread_l + (k_index-1) * ng(2) * ng(1) + (j_index-1) * ng(1) + (i_index-1) + 1
                lgrid%r(ig)=lgrid%r(ig)+res
             ELSE
                grid(i_index,j_index,k_index)=grid(i_index,j_index,k_index)+res
             ENDIF
          ENDDO
        ENDDO
      ENDDO
      !t2=nanotime_ia32()
      !write(*,*) t2-t1
      ! deallocation needed to pass around a pgi bug..
      DEALLOCATE(coef_ijk)
      DEALLOCATE(coef_map)
      DEALLOCATE(hmatgridp)
      DEALLOCATE(grid_map)

    END SUBROUTINE collocate_general_opt

! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE collocate_general_subpatch()
    INTEGER, DIMENSION(2, 3)                 :: local_b
    INTEGER, DIMENSION(3)                    :: local_s, periodic
    REAL(dp), &
      DIMENSION((lp+1)*(lp+2)*(lp+3)/6)      :: poly_d3

        periodic=1 ! cell%perd
        CALL poly_cp2k2d3(coef_xyz,lp,poly_d3)
        local_b(1,:)=rsgrid%lb_real-rsgrid%desc%lb
        local_b(2,:)=rsgrid%ub_real-rsgrid%desc%lb
        local_s=rsgrid%lb_real-rsgrid%lb_local
        IF (BTEST(subpatch_pattern,0)) local_b(1,1)=local_b(1,1)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,1)) local_b(2,1)=local_b(2,1)+rsgrid%desc%border
        IF (BTEST(subpatch_pattern,2)) local_b(1,2)=local_b(1,2)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,3)) local_b(2,2)=local_b(2,2)+rsgrid%desc%border
        IF (BTEST(subpatch_pattern,4)) local_b(1,3)=local_b(1,3)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,5)) local_b(2,3)=local_b(2,3)+rsgrid%desc%border
        IF (BTEST(subpatch_pattern,0)) local_s(1)=local_s(1)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,2)) local_s(2)=local_s(2)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,4)) local_s(3)=local_s(3)-rsgrid%desc%border
        IF ( PRESENT ( lgrid ) ) THEN
          CALL collocGauss(h=cell%hmat,h_inv=cell%h_inv,&
            grid=grid,poly=poly_d3,alphai=zetp,posi=rp,max_r2=radius*radius,&
            periodic=periodic,gdim=ng,local_bounds=local_b,local_shift=local_s,&
            lgrid=lgrid)
        ELSE
          CALL collocGauss(h=cell%hmat,h_inv=cell%h_inv,&
            grid=grid,poly=poly_d3,alphai=zetp,posi=rp,max_r2=radius*radius,&
            periodic=periodic,gdim=ng,local_bounds=local_b,local_shift=local_s)
        END IF
        ! defaults: local_shift=(/0,0,0/),poly_shift=(/0.0_dp,0.0_dp,0.0_dp/),scale=1.0_dp,

    END SUBROUTINE

! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE collocate_general_wings()
    INTEGER, DIMENSION(2, 3)                 :: local_b
    INTEGER, DIMENSION(3)                    :: periodic
    REAL(dp), &
      DIMENSION((lp+1)*(lp+2)*(lp+3)/6)      :: poly_d3
    REAL(dp), DIMENSION(3)                   :: local_shift, rShifted

        periodic=1 ! cell%perd
        CALL poly_cp2k2d3(coef_xyz,lp,poly_d3)
        local_b(1,:)=0
        local_b(2,:)=MIN(rsgrid%desc%npts-1,rsgrid%ub_local-rsgrid%lb_local)
        local_shift=REAL(rsgrid%desc%lb-rsgrid%lb_local,dp)/REAL(rsgrid%desc%npts,dp)
        rShifted(1)=rp(1)+cell%hmat(1,1)*local_shift(1)&
             +cell%hmat(1,2)*local_shift(2)&
             +cell%hmat(1,3)*local_shift(3)
        rShifted(2)=rp(2)+cell%hmat(2,1)*local_shift(1)&
             +cell%hmat(2,2)*local_shift(2)&
             +cell%hmat(2,3)*local_shift(3)
        rShifted(3)=rp(3)+cell%hmat(3,1)*local_shift(1)&
             +cell%hmat(3,2)*local_shift(2)&
             +cell%hmat(3,3)*local_shift(3)
        IF ( PRESENT ( lgrid ) ) THEN
          CALL collocGauss(h=cell%hmat,h_inv=cell%h_inv,&
            grid=grid,poly=poly_d3,alphai=zetp,posi=rShifted,max_r2=radius*radius,&
            periodic=periodic,gdim=ng,local_bounds=local_b,&
            lgrid=lgrid)
        ELSE
          CALL collocGauss(h=cell%hmat,h_inv=cell%h_inv,&
            grid=grid,poly=poly_d3,alphai=zetp,posi=rShifted,max_r2=radius*radius,&
            periodic=periodic,gdim=ng,local_bounds=local_b)
        END IF
        ! defaults: local_shift=(/0,0,0/),poly_shift=(/0.0_dp,0.0_dp,0.0_dp/),scale=1.0_dp,

    END SUBROUTINE

!
!   this is a general 'reference' routine to do the collocation
!
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE collocate_general()

    INTEGER                                  :: i, index_max(3), &
                                                index_min(3), ipoint(3), j, k
    REAL(KIND=dp)                            :: point(3)

! still hard-wired (see MODULO)

      CPASSERT(ALL(rsgrid % desc % perd==1))

      CALL return_cube_nonortho(cube_info,radius,index_min,index_max,rp)

      ! go over the grid, but cycle if the point is not within the radius
      DO k=index_min(3),index_max(3)
      DO j=index_min(2),index_max(2)
      DO i=index_min(1),index_max(1)
         ! point in real space
         point=MATMUL(cell%hmat,REAL((/i,j,k/),KIND=dp)/ng)
         ! skip if outside of the sphere
         IF (SUM((point-rp)**2)>radius**2) CYCLE
         ! point on the grid (including pbc)
         ipoint=MODULO((/i,j,k/),ng)+1
         ! add to grid
         IF ( PRESENT ( lgrid ) ) THEN
            ig = lgrid%ldim * ithread_l + ipoint(3) * ng(2) * ng(1) + ipoint(2) * ng(1) + ipoint(1) + 1
            lgrid%r(ig)=lgrid%r(ig)+primitive_value(point)
         ELSE
            grid(ipoint(1),ipoint(2),ipoint(3))=grid(ipoint(1),ipoint(2),ipoint(3))+primitive_value(point)
         ENDIF
      ENDDO
      ENDDO
      ENDDO

    END SUBROUTINE collocate_general

! *****************************************************************************
!> \brief ...
!> \param point ...
!> \retval res ...
! *****************************************************************************
    FUNCTION primitive_value(point) RESULT(res)
    REAL(KIND=dp)                            :: point(3), res

    REAL(KIND=dp)                            :: dra(3), drap(3), drb(3), &
                                                drbp(3), myexp

        res=0.0_dp
        myexp=EXP(-zetp*SUM((point-rp)**2))*prefactor
        dra=point-ra
        drb=point-rb
        drap(1)=1.0_dp
        DO lxa=0,la_max_local
        drbp(1)=1.0_dp
        DO lxb=0,lb_max_local
           drap(2)=1.0_dp
           DO lya=0,la_max_local-lxa
           drbp(2)=1.0_dp
           DO lyb=0,lb_max_local-lxb
              drap(3)=1.0_dp
              DO lza=1,MAX(la_min_local-lxa-lya,0)
                 drap(3)=drap(3)*dra(3)
              ENDDO
              DO lza=MAX(la_min_local-lxa-lya,0),la_max_local-lxa-lya
              drbp(3)=1.0_dp
              DO lzb=1,MAX(lb_min_local-lxb-lyb,0)
                 drbp(3)=drbp(3)*drb(3)
              ENDDO
              DO lzb=MAX(lb_min_local-lxb-lyb,0),lb_max_local-lxb-lyb
                ico=coset(lxa,lya,lza)
                jco=coset(lxb,lyb,lzb)
                res=res+pab_local(ico+o1_local,jco+o2_local)*myexp*PRODUCT(drap)*PRODUCT(drbp)
                drbp(3)=drbp(3)*drb(3)
              ENDDO
              drap(3)=drap(3)*dra(3)
              ENDDO
           drbp(2)=drbp(2)*drb(2)
           ENDDO
           drap(2)=drap(2)*dra(2)
           ENDDO
        drbp(1)=drbp(1)*drb(1)
        ENDDO
        drap(1)=drap(1)*dra(1)
        ENDDO

    END FUNCTION primitive_value
