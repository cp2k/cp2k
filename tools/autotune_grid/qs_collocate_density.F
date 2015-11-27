MODULE qs_collocate_density
  use kinds
  use util
  use l_utils
  use cube_utils
  use qs_interactions
  use orbital_pointers
  use basic_types
  IMPLICIT NONE

CONTAINS

  SUBROUTINE collocate_all_type_1(p,rs_grids,cube_info,l_info,eps_rho_rspace,time)
    TYPE(realspace_grid_type), DIMENSION(:) :: rs_grids
    TYPE(cube_info_type), DIMENSION(:)      :: cube_info
    TYPE(l_info_type)                       :: l_info
    REAL(KIND=dp)                           :: eps_rho_rspace
    TYPE(primitive_type), DIMENSION(:)      :: p
    REAL                                    :: t1,t2,time
    INTEGER, DIMENSION(:), ALLOCATABLE      :: r

    INTEGER :: i,j, lmax_global

    ALLOCATE(r(SIZE(p)))


    DO i=1,SIZE(rs_grids)
       rs_grids(i)%r=0.0_dp
    ENDDO

    CALL CPU_TIME(t1)

    CALL sort_z(rs_grids,p,r)

    ! Find lmax_global over all p values
    lmax_global=0
    DO i=1,SIZE(p) 
       lmax_global=MAX(p(i)%la_max,lmax_global)
    END DO 

    DO i=1,SIZE(p) 
       j=r(i)
       CALL collocate_pgf_product_rspace(p(j)%la_max,p(j)%zeta,p(j)%la_min, &
                                         p(j)%lb_max,p(j)%zetb,p(j)%lb_min, &
                                         p(j)%ra,p(j)%rab,p(j)%rab2, &
                                         rs_grids(p(j)%igrid),cube_info(p(j)%igrid),&
                                         l_info,eps_rho_rspace,lmax_global)
    ENDDO

    CALL CPU_TIME(t2)
    time=t2-t1

  END SUBROUTINE

  SUBROUTINE sort_z(rs_grids,p,r)
    TYPE(realspace_grid_type), DIMENSION(:) :: rs_grids
    TYPE(primitive_type), DIMENSION(:)      :: p
    INTEGER, DIMENSION(:)      :: r
    INTEGER, ALLOCATABLE, DIMENSION(:)      :: z
    INTEGER :: i,n
    REAL(kind=dp) :: zetp,f,rap(3),rp(3)
    INTEGER :: cubecenter(3)
 
    n=SIZE(p)
    ALLOCATE(z(n))
    DO i=1,n
           zetp      = p(i)%zeta + p(i)%zetb
           f         = p(i)%zetb/zetp
           rap(:)    = f*p(i)%rab(:)
           rp(:)     = p(i)%ra(:) + rap(:)
           cubecenter(:) = MODULO(FLOOR(rp(:)/rs_grids(p(i)%igrid)%dr(:)),rs_grids(p(i)%igrid)%npts(:))
           z(i)=10000*p(i)%igrid+cubecenter(3)*10+p(i)%la_max+p(i)%lb_max
    ENDDO
    CALL sort(z,n,r)
  END SUBROUTINE

  SUBROUTINE collocate_pgf_product_rspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2, &
                                          rsgrid,cube_info,l_info,&
                                          eps_rho_rspace,lmax_global)

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
    INTEGER, INTENT(IN)                      :: lmax_global

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
    REAL(KIND=dp) :: a, b, binomial_k_lxa, binomial_l_lxb, cutoff, f, pg, &
      prefactor, radius, rpg, ya, yap, yb, ybp, za, zap, zb, zbp, zetp, prefetch
    REAL(KIND=dp), DIMENSION(3)              :: dr, rap, rb, rbp, roffset, rp
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid
    REAL(KIND=dp), POINTER                   :: alpha_old(:,:), dpy(:,:), &
                                                dpz(:,:), polx(:,:), &
                                                poly(:,:), polz(:,:), pzyx(:)

    INTEGER :: lxp,lyp,lzp,lp,iaxis,lxpm,lypm
    REAL(kind=dp) :: p_ele,ax,ay,az
    REAL(kind=dp), DIMENSION(0:lmax_global*2,0:lmax_global,0:lmax_global,3) :: alpha
    REAL(kind=dp), DIMENSION(((lmax_global*2+1)*(lmax_global*2+2)*(lmax_global*2+3))/6) :: coef_xyz
    REAL(kind=dp), DIMENSION(((lmax_global*2+1)*(lmax_global*2+2))/2) :: coef_xyt
    REAL(kind=dp), DIMENSION(0:lmax_global*2) :: coef_xtt
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: pol_z
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: pol_y
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: pol_x
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

!   *** get the sub grid properties for the given radius ***
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)

!   *** get the l_info logic and arrays ***
    CALL return_l_info(l_info,la_min_local,la_max_local,lb_min_local,lb_max_local,&
                       ithread_l,lx_max, lxy_max,lxyz_max,ly_max,lz_max, &
                       map,polx,poly,polz,dpy,dpz,alpha_old,pzyx,cmax)

!   *** position of the gaussian product
!
!   this is the actual definition of the position on the grid
!   i.e. a point rp(:) gets here grid coordinates
!   MODULO(rp(:)/dr(:),ng(:))+1
!   hence (0.0,0.0,0.0) in real space is rsgrid%lb on the rsgrid ((1,1,1) on grid)
!

    cubecenter(:) = FLOOR(rp(:)/dr(:))
    roffset(:)    = rp(:) - REAL(cubecenter(:),dp)*dr(:)
!   *** a mapping so that the ig corresponds to the right grid point
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

    ! Allocate scratch arrays here
    lp=la_max_local+lb_max_local
    ALLOCATE(pol_z(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_y(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_x(0:lp,-cmax:cmax))

    INCLUDE 'prep.f90'

!
!   compute P_{lxp,lyp,lzp} given P_{lxa,lya,lza,lxb,lyb,lzb} and alpha(ls,lxa,lxb,1)
!   use a three step procedure
!   we don't store zeros, so counting is done using lxyz,lxy in order to have contiguous memory access in collocate_fast.F
!
       lxyz=0
       DO lzp=0,lp
       DO lyp=0,lp-lzp
       DO lxp=0,lp-lzp-lyp
           lxyz=lxyz+1
           coef_xyz(lxyz)=0.0_dp
       ENDDO
       ENDDO
       ENDDO
       DO lzb=0,lb_max_local
       DO lza=0,la_max_local
          lxy=0
          DO lyp=0,lp-lza-lzb
             DO lxp=0,lp-lza-lzb-lyp
                lxy=lxy+1
                coef_xyt(lxy)=0.0_dp
             ENDDO
             lxy=lxy+lza+lzb
          ENDDO
          DO lyb=0,lb_max_local-lzb
          DO lya=0,la_max_local-lza
             lxpm=(lb_max_local-lzb-lyb)+(la_max_local-lza-lya) 
             coef_xtt(0:lxpm)=0.0_dp
             DO lxb=MAX(lb_min_local-lzb-lyb,0),lb_max_local-lzb-lyb
             DO lxa=MAX(la_min_local-lza-lya,0),la_max_local-lza-lya
                ico=coset(lxa,lya,lza)
                jco=coset(lxb,lyb,lzb)
                p_ele=prefactor !pab(ico,jco)
                DO lxp=0,lxa+lxb
                   coef_xtt(lxp)=coef_xtt(lxp)+p_ele*alpha(lxp,lxa,lxb,1)
                ENDDO
             ENDDO
             ENDDO
             lxy=0
             DO lyp=0,lya+lyb
                DO lxp=0,lp-lza-lzb-lya-lyb
                  lxy=lxy+1
                  coef_xyt(lxy)=coef_xyt(lxy)+alpha(lyp,lya,lyb,2)*coef_xtt(lxp)
                ENDDO
                lxy=lxy+lza+lzb+lya+lyb-lyp
             ENDDO
          ENDDO
          ENDDO
          lxyz=0
          DO lzp=0,lza+lzb
             lxy=0
             DO lyp=0,lp-lza-lzb
                DO lxp=0,lp-lza-lzb-lyp
                       lxy=lxy+1 ; lxyz=lxyz+1
                       coef_xyz(lxyz)=coef_xyz(lxyz)+alpha(lzp,lza,lzb,3)*coef_xyt(lxy)
                ENDDO
                lxy=lxy+lza+lzb ; lxyz=lxyz+lza+lzb-lzp
             ENDDO
             DO lyp=lp-lza-lzb+1,lp-lzp
                DO lxp=0,lp-lyp-lzp
                   lxyz=lxyz+1
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       ENDDO

    INCLUDE 'call_collocate.f90'

    CALL finish_l_info(polx,poly,polz,pzyx)

  END SUBROUTINE collocate_pgf_product_rspace
END MODULE 
