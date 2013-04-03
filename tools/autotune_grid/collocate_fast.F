  SUBROUTINE collocate_core_default(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,lp,cmax,gridbounds)

    IMPLICIT NONE

  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND ( 14, 200 )
    integer, INTENT(IN)     :: lp
    integer, INTENT(IN)     :: cmax
    integer, INTENT(IN)     :: gridbounds(2,3)
    integer, INTENT(IN)     :: sphere_bounds(*)
    integer, INTENT(IN)     :: map(-cmax:cmax,1:3)
    real(wp), INTENT(IN)    :: pol_x(0:lp,-cmax:cmax)
    real(wp), INTENT(IN)    :: pol_y(1:2,0:lp,-cmax:0)
    real(wp), INTENT(IN)    :: pol_z(1:2,0:lp,-cmax:0)
    real(wp), INTENT(IN)    :: coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)
    real(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), & 
                                    gridbounds(1,2):gridbounds(2,2), &
                                    gridbounds(1,3):gridbounds(2,3))

    real(wp) ::  coef_xy(2,(lp+1)*(lp+2)/2)
    real(wp) ::  coef_x(4,0:lp)

    integer kg,k,jgmin,jg,j,j2,igmax,ig,i,kgmin,igmin,k2,i2,jg2,kg2
    integer sci,lxp,lyp,lzp
    real(wp) s01,s02,s03,s04,s(4)
    integer lxyz,lxy,lx

    sci=1

    kgmin=sphere_bounds(sci)
    sci=sci+1
    DO kg=kgmin,0
       kg2=1-kg
       k=map(kg,3)
       k2=map(kg2,3)

       coef_xy=0.0_wp
       lxyz = 0
       DO lzp=0,lp
          lxy=0
          DO lyp=0,lp-lzp
             DO lxp=0,lp-lzp-lyp
                lxyz=lxyz+1 ; lxy=lxy+1
                coef_xy(1,lxy)=coef_xy(1,lxy)+coef_xyz(lxyz)*pol_z(1,lzp,kg)
                coef_xy(2,lxy)=coef_xy(2,lxy)+coef_xyz(lxyz)*pol_z(2,lzp,kg)
             ENDDO
             lxy=lxy+lzp
          ENDDO
       ENDDO

       jgmin=sphere_bounds(sci)
       sci=sci+1
       DO jg=jgmin,0
          jg2=1-jg
          j=map(jg,2)
          j2=map(jg2,2)
          igmin=sphere_bounds(sci)
          sci=sci+1
          igmax=1-igmin

          coef_x=0.0_wp
          lxy=0
          DO lyp=0,lp
          DO lxp=0,lp-lyp
             lxy=lxy+1
             coef_x(1,lxp)=coef_x(1,lxp)+coef_xy(1,lxy)*pol_y(1,lyp,jg)
             coef_x(2,lxp)=coef_x(2,lxp)+coef_xy(2,lxy)*pol_y(1,lyp,jg)
             coef_x(3,lxp)=coef_x(3,lxp)+coef_xy(1,lxy)*pol_y(2,lyp,jg)
             coef_x(4,lxp)=coef_x(4,lxp)+coef_xy(2,lxy)*pol_y(2,lyp,jg)
          ENDDO
          ENDDO

          DO ig=igmin,igmax
             i=map(ig,1)
             s01=0.0_wp
             s02=0.0_wp
             s03=0.0_wp
             s04=0.0_wp
             DO lxp=0,lp
                s01=s01+coef_x(1,lxp)*pol_x(lxp,ig)
                s02=s02+coef_x(2,lxp)*pol_x(lxp,ig)
                s03=s03+coef_x(3,lxp)*pol_x(lxp,ig)
                s04=s04+coef_x(4,lxp)*pol_x(lxp,ig)
             ENDDO
             grid(i,j,k) = grid(i,j,k)     + s01 
             grid(i,j2,k) = grid(i,j2,k)   + s03 
             grid(i,j,k2) = grid(i,j,k2)   + s02 
             grid(i,j2,k2) = grid(i,j2,k2) + s04 
          END DO

       END DO
    END DO
  
  END SUBROUTINE collocate_core_default
