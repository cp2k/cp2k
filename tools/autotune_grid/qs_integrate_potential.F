MODULE qs_integrate_potential
  use kinds
  use util
  use l_utils
  use cube_utils
  use qs_interactions
  use orbital_pointers
  use basic_types
  IMPLICIT NONE

CONTAINS

  SUBROUTINE integrate_all(p,rs_grids,cube_info,l_info,eps_rho_rspace,time)
    TYPE(realspace_grid_type), DIMENSION(:) :: rs_grids
    TYPE(cube_info_type), DIMENSION(:)      :: cube_info
    TYPE(l_info_type)                       :: l_info
    REAL(KIND=dp)                           :: eps_rho_rspace,intsum
    TYPE(primitive_type), DIMENSION(:)      :: p
    REAL                                    :: t1,t2,time
    INTEGER, DIMENSION(:), ALLOCATABLE      :: r

    INTEGER :: i,j,k,l,m

    ALLOCATE(r(SIZE(p)))

    DO i=1,SIZE(rs_grids)
     DO k=LBOUND(rs_grids(i)%r,3),UBOUND(rs_grids(i)%r,3)
     DO l=LBOUND(rs_grids(i)%r,2),UBOUND(rs_grids(i)%r,2)
     DO m=LBOUND(rs_grids(i)%r,1),UBOUND(rs_grids(i)%r,1)
       rs_grids(i)%r(m,l,k)=1.0_dp/sqrt(1.0_dp+(REAL(m+3*k+5*l,dp))**2)
     ENDDO
     ENDDO
     ENDDO
    ENDDO

    CALL CPU_TIME(t1)

    DO i=1,SIZE(p) 
       j=i
       CALL integrate_pgf_product_rspace(p(j)%la_max,p(j)%zeta,p(j)%la_min, &
                                         p(j)%lb_max,p(j)%zetb,p(j)%lb_min, &
                                         p(j)%ra,p(j)%rab,p(j)%rab2, &
                                         rs_grids(p(j)%igrid),cube_info(p(j)%igrid),l_info,eps_rho_rspace,&
                                         intsum)
       p(j)%intsum=intsum
    ENDDO

    CALL CPU_TIME(t2)
    time=t2-t1

  END SUBROUTINE

  SUBROUTINE integrate_pgf_product_rspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2, &
                                          rsgrid,cube_info,l_info,&
                                          eps_rho_rspace,&
                                          intsum)

    INTEGER, INTENT(IN)                      :: la_max
    REAL(KIND=dp), INTENT(IN)                :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(KIND=dp), INTENT(IN)                :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rab
    REAL(KIND=dp), INTENT(IN)                :: rab2
    TYPE(realspace_grid_type)                :: rsgrid
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    REAL(KIND=dp), INTENT(IN)                :: eps_rho_rspace
    REAL(KIND=dp), INTENT(OUT)               :: intsum

    INTEGER :: cmax, coef_max, gridbounds(2,3), i, ico, icoef, ig, ithread_l, &
      jco, k, l, la_max_local, la_min_local, lb_max_local, lb_min_local, &
      length, lx, lx_max, lxa, lxb, lxy, lxy_max, lxyz, lxyz_max, lya, lyb, &
      lza, lzb, offset, start
    INTEGER, DIMENSION(3)                    :: cubecenter, lb_cube, ng, &
                                                ub_cube
    INTEGER, DIMENSION(:), POINTER           :: ly_max, lz_max, sphere_bounds
    INTEGER, DIMENSION(:, :), POINTER        :: map
    LOGICAL                                  :: failure, &
                                                my_map_consistent
    REAL(KIND=dp) :: a, b, binomial_k_lxa, binomial_l_lxb, cutoff, pg, &
      prefactor, radius, rpg, ya, yap, yb, ybp, za, zap, zb, zbp, zetp, prefetch, &
      exp_x0, exp_x1, exp_x2, f, ftza, ftzb
    REAL(KIND=dp), DIMENSION(3)              :: dr, rap, rb, rbp, roffset, rp
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid
    REAL(KIND=dp), POINTER                   :: alpha_old(:,:), dpy(:,:), &
                                                dpz(:,:), polx(:,:), &
                                                poly(:,:), polz(:,:), pzyx(:)

    INTEGER :: lxp,lyp,lzp,lp,iaxis
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: alpha
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: coef_xyz
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: coef_xyt
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: coef_xtt
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: coef_ttz
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: coef_tyz

    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: pol_z
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: pol_y
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: pol_x
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: vab
    REAL(KIND=dp) :: t_exp_1,t_exp_2,t_exp_min_1,t_exp_min_2,t_exp_plus_1,t_exp_plus_2


!   ---------------------------------------------------------------------------

    failure = .FALSE.

    ithread_l=0

    my_map_consistent=.TRUE.


    zetp      = zeta + zetb
    f         = zetb/zetp
    rap(:)    = f*rab(:)
    rbp(:)    = rap(:) - rab(:)
    rp(:)     = ra(:) + rap(:)
    rb(:)     = ra(:)+rab(:)

    cutoff    = 1.0_dp
    prefactor = EXP(-zeta*f*rab2)
    radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,ra=ra,rb=rb,rp=rp,&
                                       zetp=zetp,eps=eps_rho_rspace,&
                                       prefactor=prefactor,cutoff=cutoff)
    prefactor = EXP(-zeta*f*rab2)

    IF (radius .EQ. 0.0_dp ) THEN
      RETURN
    END IF

    la_max_local=la_max
    la_min_local=la_min
    lb_max_local=lb_max
    lb_min_local=lb_min

    coef_max=la_max_local+lb_max_local+1

!   *** properties of the grid ***
    dr(:) = rsgrid%dr(:)
    ng(:) = rsgrid%npts(:)

! WARNING: this resets the lower bounds of grid to 1 (unlike grid => rsgrid%r)
    grid => rsgrid%r(:,:,:)
    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)

    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)
    CALL return_l_info(l_info,la_min_local,la_max_local,lb_min_local,lb_max_local, &
                       ithread_l,lx_max,lxy_max,lxyz_max,ly_max,lz_max,&
                       map,polx,poly,polz,dpy,dpz,alpha_old,pzyx,cmax)

    dr(:) = rsgrid%dr(:)
    ng(:) = rsgrid%npts(:)

    grid => rsgrid%r(:,:,:)
    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)

    cubecenter(:) = FLOOR(rp(:)/dr(:))
    roffset(:) = rp(:) - REAL(cubecenter(:),dp)*dr(:)

    DO i=1,3
      IF ( rsgrid % perd ( i ) == 1 ) THEN
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
        offset=MODULO(cubecenter(i),ng(i))+rsgrid%lb(i)
        offset=offset-rsgrid%lb_local(i)+1
        DO ig=lb_cube(i),ub_cube(i)
           map(ig,i) = ig+offset
        END DO
      END IF
    ENDDO

    INCLUDE 'prep.f90'

    INCLUDE 'call_integrate.f90'

!    CALL integrate_core_default(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax),pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
!                                map(-cmax,1),sphere_bounds(1),lp,cmax,gridbounds(1,1))

    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    ALLOCATE(vab(ncoset(la_max_local),ncoset(lb_max_local)))
    ALLOCATE(coef_ttz(0:la_max_local,0:lb_max_local))
    ALLOCATE(coef_tyz(0:la_max_local,0:lb_max_local,0:la_max_local,0:lb_max_local))


    IF (.TRUE.) THEN

       vab=0.0_dp
       lxyz=0
       DO lzp=0,lp
          coef_tyz=0.0_dp
          DO lyp=0,lp-lzp
             coef_ttz=0.0_dp
             DO lxp=0,lp-lzp-lyp
                lxyz=lxyz+1
                DO lxb=0,lb_max_local
                DO lxa=0,la_max_local
                   coef_ttz(lxa,lxb)=coef_ttz(lxa,lxb)+coef_xyz(lxyz)*alpha(lxp,lxa,lxb,1)
                ENDDO
                ENDDO

             ENDDO

             DO lyb=0,lb_max_local
             DO lya=0,la_max_local
                DO lxb=0,lb_max_local-lyb
                DO lxa=0,la_max_local-lya
                   coef_tyz(lxa,lxb,lya,lyb)=coef_tyz(lxa,lxb,lya,lyb)+coef_ttz(lxa,lxb)*alpha(lyp,lya,lyb,2)
                ENDDO
                ENDDO
             ENDDO
             ENDDO

          ENDDO

          DO lzb=0,lb_max_local
          DO lza=0,la_max_local
             DO lyb=0,lb_max_local-lzb
             DO lya=0,la_max_local-lza
                DO lxb=MAX(lb_min_local-lzb-lyb,0),lb_max_local-lzb-lyb
                jco=coset(lxb,lyb,lzb)
                DO lxa=MAX(la_min_local-lza-lya,0),la_max_local-lza-lya
                   ico=coset(lxa,lya,lza)
                   vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,lya,lyb)*alpha(lzp,lza,lzb,3)
                ENDDO
                ENDDO
             ENDDO
             ENDDO
          ENDDO
          ENDDO

       ENDDO

    ELSE

       DO lzb=0,lb_max_local
       DO lza=0,la_max_local
          DO lyb=0,lb_max_local-lzb
          DO lya=0,la_max_local-lza
             DO lxb=MAX(lb_min_local-lzb-lyb,0),lb_max_local-lzb-lyb
             DO lxa=MAX(la_min_local-lza-lya,0),la_max_local-lza-lya
                ico=coset(lxa,lya,lza)
                jco=coset(lxb,lyb,lzb)
                vab(ico,jco)=0.0_dp
                lxyz=0
                DO lzp=0,lp
                 DO lyp=0,lp-lzp
                  DO lxp=0,lp-lzp-lyp
                   lxyz=lxyz+1
                   vab(ico,jco)=vab(ico,jco)+coef_xyz(lxyz)*alpha(lzp,lza,lzb,3)*alpha(lyp,lya,lyb,2)*alpha(lxp,lxa,lxb,1)
                  ENDDO
                 ENDDO
                ENDDO
             ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
       ENDDO

    ENDIF

    ! compute a kind of forces checksum
    intsum=0.0_dp
    lxyz=0
      DO lxa=0,la_max
        DO lxb=0,lb_max
          DO lya=0,la_max-lxa
            DO lyb=0,lb_max-lxb
              DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
                DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
                  ico=coset(lxa,lya,lza)
                  jco=coset(lxb,lyb,lzb)
                  intsum=intsum+vab(ico,jco)*(coset(lxa,lya,lza)+1000*coset(lxb,lyb,lzb))
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      CALL finish_l_info(polx,poly,polz,pzyx)

  END SUBROUTINE integrate_pgf_product_rspace
END MODULE  qs_integrate_potential
