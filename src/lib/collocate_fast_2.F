  SUBROUTINE collocate_core_default(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,lp,cmax,gridbounds)
    USE lib_kinds,                        ONLY: wp
    INTEGER, INTENT(IN)                      :: sphere_bounds(*), lp
    REAL(wp), INTENT(IN) :: coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)
    INTEGER, INTENT(IN)                      :: cmax
    REAL(wp), INTENT(IN)                     :: pol_x(0:lp,-cmax:cmax), &
                                                pol_y(1:2,0:lp,-cmax:0), &
                                                pol_z(1:2,0:lp,-cmax:0)
    INTEGER, INTENT(IN)                      :: map(-cmax:cmax,1:3), &
                                                gridbounds(2,3)
    REAL(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), &
      gridbounds(1,2):gridbounds(2,2), gridbounds(1,3):gridbounds(2,3))

    INTEGER                                  :: i, ig, igmax, igmin, j, j2, &
                                                jg, jg2, jgmin, k, k2, kg, &
                                                kg2, kgmin, lxp, lxy, lxyz, &
                                                lyp, lzp, sci
    REAL(wp)                                 :: coef_x(4,0:lp), &
                                                coef_xy(2,(lp+1)*(lp+2)/2), &
                                                s01, s02, s03, s04

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
  SUBROUTINE collocate_core_0(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)
    USE lib_kinds,                        ONLY: wp
    INTEGER, INTENT(IN)                      :: sphere_bounds(*), cmax, &
                                                map(-cmax:cmax,1:3), &
                                                gridbounds(2,3)
    REAL(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), &
      gridbounds(1,2):gridbounds(2,2), gridbounds(1,3):gridbounds(2,3))
    INTEGER, PARAMETER                       :: lp = 0
    REAL(wp), INTENT(IN) :: pol_x(0:lp,-cmax:cmax), pol_y(1:2,0:lp,-cmax:0), &
      pol_z(1:2,0:lp,-cmax:0), coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)

    INTEGER                                  :: i, ig, igmax, igmin, j, j2, &
                                                jg, jg2, jgmin, k, k2, kg, &
                                                kg2, kgmin, lxp, lxy, lxyz, &
                                                lyp, lzp, sci
    REAL(wp)                                 :: coef_x(4,0:lp), &
                                                coef_xy(2,(lp+1)*(lp+2)/2), &
                                                s01, s02, s03, s04

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
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,1)*pol_y(1,0,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,1)*pol_y(2,0,jg)
          DO ig=igmin,igmax
             i=map(ig,1)
             s01=0.0_wp
             s02=0.0_wp
             s03=0.0_wp
             s04=0.0_wp
                s01=s01+coef_x(1,0)*pol_x(0,ig)
                s02=s02+coef_x(2,0)*pol_x(0,ig)
                s03=s03+coef_x(3,0)*pol_x(0,ig)
                s04=s04+coef_x(4,0)*pol_x(0,ig)
             grid(i,j,k) = grid(i,j,k)     + s01
             grid(i,j2,k) = grid(i,j2,k)   + s03
             grid(i,j,k2) = grid(i,j,k2)   + s02
             grid(i,j2,k2) = grid(i,j2,k2) + s04
          END DO
       END DO
    END DO

  END SUBROUTINE collocate_core_0
  SUBROUTINE collocate_core_1(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)
    USE lib_kinds,                        ONLY: wp
    INTEGER, INTENT(IN)                      :: sphere_bounds(*), cmax, &
                                                map(-cmax:cmax,1:3), &
                                                gridbounds(2,3)
    REAL(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), &
      gridbounds(1,2):gridbounds(2,2), gridbounds(1,3):gridbounds(2,3))
    INTEGER, PARAMETER                       :: lp = 1
    REAL(wp), INTENT(IN) :: pol_x(0:lp,-cmax:cmax), pol_y(1:2,0:lp,-cmax:0), &
      pol_z(1:2,0:lp,-cmax:0), coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)

    INTEGER                                  :: i, ig, igmax, igmin, j, j2, &
                                                jg, jg2, jgmin, k, k2, kg, &
                                                kg2, kgmin, lxp, lxy, lxyz, &
                                                lyp, lzp, sci
    REAL(wp)                                 :: coef_x(4,0:lp), &
                                                coef_xy(2,(lp+1)*(lp+2)/2), &
                                                s01, s02, s03, s04

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
                coef_xy(:,lxy)=coef_xy(:,lxy)+coef_xyz(lxyz)*pol_z(:,lzp,kg)
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
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,1)*pol_y(1,0,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,1)*pol_y(2,0,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,2)*pol_y(1,0,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,2)*pol_y(2,0,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,3)*pol_y(1,1,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,3)*pol_y(2,1,jg)
          DO ig=igmin,igmax
             i=map(ig,1)
             s01=0.0_wp
             s02=0.0_wp
             s03=0.0_wp
             s04=0.0_wp
                s01=s01+coef_x(1,0)*pol_x(0,ig)
                s02=s02+coef_x(2,0)*pol_x(0,ig)
                s03=s03+coef_x(3,0)*pol_x(0,ig)
                s04=s04+coef_x(4,0)*pol_x(0,ig)
                s01=s01+coef_x(1,1)*pol_x(1,ig)
                s02=s02+coef_x(2,1)*pol_x(1,ig)
                s03=s03+coef_x(3,1)*pol_x(1,ig)
                s04=s04+coef_x(4,1)*pol_x(1,ig)
             grid(i,j,k) = grid(i,j,k)     + s01
             grid(i,j2,k) = grid(i,j2,k)   + s03
             grid(i,j,k2) = grid(i,j,k2)   + s02
             grid(i,j2,k2) = grid(i,j2,k2) + s04
          END DO
       END DO
    END DO

  END SUBROUTINE collocate_core_1
  SUBROUTINE collocate_core_2(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)
    USE lib_kinds,                        ONLY: wp
    INTEGER, INTENT(IN)                      :: sphere_bounds(*), cmax, &
                                                map(-cmax:cmax,1:3), &
                                                gridbounds(2,3)
    REAL(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), &
      gridbounds(1,2):gridbounds(2,2), gridbounds(1,3):gridbounds(2,3))
    INTEGER, PARAMETER                       :: lp = 2
    REAL(wp), INTENT(IN) :: pol_x(0:lp,-cmax:cmax), pol_y(1:2,0:lp,-cmax:0), &
      pol_z(1:2,0:lp,-cmax:0), coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)

    INTEGER                                  :: i, ig, igmax, igmin, j, j2, &
                                                jg, jg2, jgmin, k, k2, kg, &
                                                kg2, kgmin, sci
    REAL(wp)                                 :: coef_x(4,0:lp), &
                                                coef_xy(2,(lp+1)*(lp+2)/2), &
                                                s01, s02, s03, s04

    sci=1

    kgmin=sphere_bounds(sci)
    sci=sci+1
    DO kg=kgmin,0
       kg2=1-kg
       k=map(kg,3)
       k2=map(kg2,3)
       coef_xy=0.0_wp
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(1)*pol_z(:,0,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(2)*pol_z(:,0,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(3)*pol_z(:,0,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(4)*pol_z(:,0,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(5)*pol_z(:,0,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(6)*pol_z(:,0,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(7)*pol_z(:,1,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(8)*pol_z(:,1,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(9)*pol_z(:,1,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(10)*pol_z(:,2,kg)
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
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,1)*pol_y(1,0,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,1)*pol_y(2,0,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,2)*pol_y(1,0,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,2)*pol_y(2,0,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,3)*pol_y(1,0,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,3)*pol_y(2,0,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,4)*pol_y(1,1,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,4)*pol_y(2,1,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,5)*pol_y(1,1,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,5)*pol_y(2,1,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,6)*pol_y(1,2,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,6)*pol_y(2,2,jg)
          DO ig=igmin,igmax
             i=map(ig,1)
             s01=0.0_wp
             s02=0.0_wp
             s03=0.0_wp
             s04=0.0_wp
                s01=s01+coef_x(1,0)*pol_x(0,ig)
                s02=s02+coef_x(2,0)*pol_x(0,ig)
                s03=s03+coef_x(3,0)*pol_x(0,ig)
                s04=s04+coef_x(4,0)*pol_x(0,ig)
                s01=s01+coef_x(1,1)*pol_x(1,ig)
                s02=s02+coef_x(2,1)*pol_x(1,ig)
                s03=s03+coef_x(3,1)*pol_x(1,ig)
                s04=s04+coef_x(4,1)*pol_x(1,ig)
                s01=s01+coef_x(1,2)*pol_x(2,ig)
                s02=s02+coef_x(2,2)*pol_x(2,ig)
                s03=s03+coef_x(3,2)*pol_x(2,ig)
                s04=s04+coef_x(4,2)*pol_x(2,ig)
             grid(i,j,k) = grid(i,j,k)     + s01
             grid(i,j2,k) = grid(i,j2,k)   + s03
             grid(i,j,k2) = grid(i,j,k2)   + s02
             grid(i,j2,k2) = grid(i,j2,k2) + s04
          END DO
       END DO
    END DO

  END SUBROUTINE collocate_core_2
  SUBROUTINE collocate_core_3(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)
    USE lib_kinds,                        ONLY: wp
    INTEGER, INTENT(IN)                      :: sphere_bounds(*), cmax, &
                                                map(-cmax:cmax,1:3), &
                                                gridbounds(2,3)
    REAL(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), &
      gridbounds(1,2):gridbounds(2,2), gridbounds(1,3):gridbounds(2,3))
    INTEGER, PARAMETER                       :: lp = 3
    REAL(wp), INTENT(IN) :: pol_x(0:lp,-cmax:cmax), pol_y(1:2,0:lp,-cmax:0), &
      pol_z(1:2,0:lp,-cmax:0), coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)

    INTEGER                                  :: i, ig, igmax, igmin, j, j2, &
                                                jg, jg2, jgmin, k, k2, kg, &
                                                kg2, kgmin, sci
    REAL(wp)                                 :: coef_x(4,0:lp), &
                                                coef_xy(2,(lp+1)*(lp+2)/2), &
                                                s01, s02, s03, s04

    sci=1

    kgmin=sphere_bounds(sci)
    sci=sci+1
    DO kg=kgmin,0
       kg2=1-kg
       k=map(kg,3)
       k2=map(kg2,3)
       coef_xy=0.0_wp
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(1)*pol_z(:,0,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(2)*pol_z(:,0,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(3)*pol_z(:,0,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(4)*pol_z(:,0,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(5)*pol_z(:,0,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(6)*pol_z(:,0,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(7)*pol_z(:,0,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(8)*pol_z(:,0,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(9)*pol_z(:,0,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(10)*pol_z(:,0,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(11)*pol_z(:,1,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(12)*pol_z(:,1,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(13)*pol_z(:,1,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(14)*pol_z(:,1,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(15)*pol_z(:,1,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(16)*pol_z(:,1,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(17)*pol_z(:,2,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(18)*pol_z(:,2,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(19)*pol_z(:,2,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(20)*pol_z(:,3,kg)
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
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,1)*pol_y(1,0,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,1)*pol_y(2,0,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,2)*pol_y(1,0,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,2)*pol_y(2,0,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,3)*pol_y(1,0,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,3)*pol_y(2,0,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,4)*pol_y(1,0,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,4)*pol_y(2,0,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,5)*pol_y(1,1,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,5)*pol_y(2,1,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,6)*pol_y(1,1,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,6)*pol_y(2,1,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,7)*pol_y(1,1,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,7)*pol_y(2,1,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,8)*pol_y(1,2,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,8)*pol_y(2,2,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,9)*pol_y(1,2,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,9)*pol_y(2,2,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,10)*pol_y(1,3,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,10)*pol_y(2,3,jg)
          DO ig=igmin,igmax
             i=map(ig,1)
             s01=0.0_wp
             s02=0.0_wp
             s03=0.0_wp
             s04=0.0_wp
                s01=s01+coef_x(1,0)*pol_x(0,ig)
                s02=s02+coef_x(2,0)*pol_x(0,ig)
                s03=s03+coef_x(3,0)*pol_x(0,ig)
                s04=s04+coef_x(4,0)*pol_x(0,ig)
                s01=s01+coef_x(1,1)*pol_x(1,ig)
                s02=s02+coef_x(2,1)*pol_x(1,ig)
                s03=s03+coef_x(3,1)*pol_x(1,ig)
                s04=s04+coef_x(4,1)*pol_x(1,ig)
                s01=s01+coef_x(1,2)*pol_x(2,ig)
                s02=s02+coef_x(2,2)*pol_x(2,ig)
                s03=s03+coef_x(3,2)*pol_x(2,ig)
                s04=s04+coef_x(4,2)*pol_x(2,ig)
                s01=s01+coef_x(1,3)*pol_x(3,ig)
                s02=s02+coef_x(2,3)*pol_x(3,ig)
                s03=s03+coef_x(3,3)*pol_x(3,ig)
                s04=s04+coef_x(4,3)*pol_x(3,ig)
             grid(i,j,k) = grid(i,j,k)     + s01
             grid(i,j2,k) = grid(i,j2,k)   + s03
             grid(i,j,k2) = grid(i,j,k2)   + s02
             grid(i,j2,k2) = grid(i,j2,k2) + s04
          END DO
       END DO
    END DO

  END SUBROUTINE collocate_core_3
  SUBROUTINE collocate_core_4(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)
    USE lib_kinds,                        ONLY: wp
    INTEGER, INTENT(IN)                      :: sphere_bounds(*), cmax, &
                                                map(-cmax:cmax,1:3), &
                                                gridbounds(2,3)
    REAL(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), &
      gridbounds(1,2):gridbounds(2,2), gridbounds(1,3):gridbounds(2,3))
    INTEGER, PARAMETER                       :: lp = 4
    REAL(wp), INTENT(IN) :: pol_x(0:lp,-cmax:cmax), pol_y(1:2,0:lp,-cmax:0), &
      pol_z(1:2,0:lp,-cmax:0), coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)

    INTEGER                                  :: i, ig, igmax, igmin, j, j2, &
                                                jg, jg2, jgmin, k, k2, kg, &
                                                kg2, kgmin, lxp, sci
    REAL(wp)                                 :: coef_x(4,0:lp), &
                                                coef_xy(2,(lp+1)*(lp+2)/2), &
                                                s01, s02, s03, s04

    sci=1

    kgmin=sphere_bounds(sci)
    sci=sci+1
    DO kg=kgmin,0
       kg2=1-kg
       k=map(kg,3)
       k2=map(kg2,3)
       coef_xy=0.0_wp
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(1)*pol_z(:,0,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(2)*pol_z(:,0,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(3)*pol_z(:,0,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(4)*pol_z(:,0,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(5)*pol_z(:,0,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(6)*pol_z(:,0,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(7)*pol_z(:,0,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(8)*pol_z(:,0,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(9)*pol_z(:,0,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(10)*pol_z(:,0,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(11)*pol_z(:,0,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(12)*pol_z(:,0,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(13)*pol_z(:,0,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(14)*pol_z(:,0,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(15)*pol_z(:,0,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(16)*pol_z(:,1,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(17)*pol_z(:,1,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(18)*pol_z(:,1,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(19)*pol_z(:,1,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(20)*pol_z(:,1,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(21)*pol_z(:,1,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(22)*pol_z(:,1,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(23)*pol_z(:,1,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(24)*pol_z(:,1,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(25)*pol_z(:,1,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(26)*pol_z(:,2,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(27)*pol_z(:,2,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(28)*pol_z(:,2,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(29)*pol_z(:,2,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(30)*pol_z(:,2,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(31)*pol_z(:,2,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(32)*pol_z(:,3,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(33)*pol_z(:,3,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(34)*pol_z(:,3,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(35)*pol_z(:,4,kg)
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
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,1)*pol_y(1,0,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,1)*pol_y(2,0,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,2)*pol_y(1,0,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,2)*pol_y(2,0,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,3)*pol_y(1,0,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,3)*pol_y(2,0,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,4)*pol_y(1,0,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,4)*pol_y(2,0,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,5)*pol_y(1,0,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,5)*pol_y(2,0,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,6)*pol_y(1,1,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,6)*pol_y(2,1,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,7)*pol_y(1,1,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,7)*pol_y(2,1,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,8)*pol_y(1,1,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,8)*pol_y(2,1,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,9)*pol_y(1,1,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,9)*pol_y(2,1,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,10)*pol_y(1,2,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,10)*pol_y(2,2,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,11)*pol_y(1,2,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,11)*pol_y(2,2,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,12)*pol_y(1,2,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,12)*pol_y(2,2,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,13)*pol_y(1,3,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,13)*pol_y(2,3,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,14)*pol_y(1,3,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,14)*pol_y(2,3,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,15)*pol_y(1,4,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,15)*pol_y(2,4,jg)
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

  END SUBROUTINE collocate_core_4
  SUBROUTINE collocate_core_5(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)
    USE lib_kinds,                        ONLY: wp
    INTEGER, INTENT(IN)                      :: sphere_bounds(*), cmax, &
                                                map(-cmax:cmax,1:3), &
                                                gridbounds(2,3)
    REAL(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), &
      gridbounds(1,2):gridbounds(2,2), gridbounds(1,3):gridbounds(2,3))
    INTEGER, PARAMETER                       :: lp = 5
    REAL(wp), INTENT(IN) :: pol_x(0:lp,-cmax:cmax), pol_y(1:2,0:lp,-cmax:0), &
      pol_z(1:2,0:lp,-cmax:0), coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)

    INTEGER                                  :: i, ig, igmax, igmin, j, j2, &
                                                jg, jg2, jgmin, k, k2, kg, &
                                                kg2, kgmin, lxp, sci
    REAL(wp)                                 :: coef_x(4,0:lp), &
                                                coef_xy(2,(lp+1)*(lp+2)/2), &
                                                s01, s02, s03, s04

    sci=1

    kgmin=sphere_bounds(sci)
    sci=sci+1
    DO kg=kgmin,0
       kg2=1-kg
       k=map(kg,3)
       k2=map(kg2,3)
       coef_xy=0.0_wp
                coef_xy(1,1)=coef_xy(1,1)+coef_xyz(1)*pol_z(1,0,kg)
                coef_xy(2,1)=coef_xy(2,1)+coef_xyz(1)*pol_z(2,0,kg)
                coef_xy(1,2)=coef_xy(1,2)+coef_xyz(2)*pol_z(1,0,kg)
                coef_xy(2,2)=coef_xy(2,2)+coef_xyz(2)*pol_z(2,0,kg)
                coef_xy(1,3)=coef_xy(1,3)+coef_xyz(3)*pol_z(1,0,kg)
                coef_xy(2,3)=coef_xy(2,3)+coef_xyz(3)*pol_z(2,0,kg)
                coef_xy(1,4)=coef_xy(1,4)+coef_xyz(4)*pol_z(1,0,kg)
                coef_xy(2,4)=coef_xy(2,4)+coef_xyz(4)*pol_z(2,0,kg)
                coef_xy(1,5)=coef_xy(1,5)+coef_xyz(5)*pol_z(1,0,kg)
                coef_xy(2,5)=coef_xy(2,5)+coef_xyz(5)*pol_z(2,0,kg)
                coef_xy(1,6)=coef_xy(1,6)+coef_xyz(6)*pol_z(1,0,kg)
                coef_xy(2,6)=coef_xy(2,6)+coef_xyz(6)*pol_z(2,0,kg)
                coef_xy(1,7)=coef_xy(1,7)+coef_xyz(7)*pol_z(1,0,kg)
                coef_xy(2,7)=coef_xy(2,7)+coef_xyz(7)*pol_z(2,0,kg)
                coef_xy(1,8)=coef_xy(1,8)+coef_xyz(8)*pol_z(1,0,kg)
                coef_xy(2,8)=coef_xy(2,8)+coef_xyz(8)*pol_z(2,0,kg)
                coef_xy(1,9)=coef_xy(1,9)+coef_xyz(9)*pol_z(1,0,kg)
                coef_xy(2,9)=coef_xy(2,9)+coef_xyz(9)*pol_z(2,0,kg)
                coef_xy(1,10)=coef_xy(1,10)+coef_xyz(10)*pol_z(1,0,kg)
                coef_xy(2,10)=coef_xy(2,10)+coef_xyz(10)*pol_z(2,0,kg)
                coef_xy(1,11)=coef_xy(1,11)+coef_xyz(11)*pol_z(1,0,kg)
                coef_xy(2,11)=coef_xy(2,11)+coef_xyz(11)*pol_z(2,0,kg)
                coef_xy(1,12)=coef_xy(1,12)+coef_xyz(12)*pol_z(1,0,kg)
                coef_xy(2,12)=coef_xy(2,12)+coef_xyz(12)*pol_z(2,0,kg)
                coef_xy(1,13)=coef_xy(1,13)+coef_xyz(13)*pol_z(1,0,kg)
                coef_xy(2,13)=coef_xy(2,13)+coef_xyz(13)*pol_z(2,0,kg)
                coef_xy(1,14)=coef_xy(1,14)+coef_xyz(14)*pol_z(1,0,kg)
                coef_xy(2,14)=coef_xy(2,14)+coef_xyz(14)*pol_z(2,0,kg)
                coef_xy(1,15)=coef_xy(1,15)+coef_xyz(15)*pol_z(1,0,kg)
                coef_xy(2,15)=coef_xy(2,15)+coef_xyz(15)*pol_z(2,0,kg)
                coef_xy(1,16)=coef_xy(1,16)+coef_xyz(16)*pol_z(1,0,kg)
                coef_xy(2,16)=coef_xy(2,16)+coef_xyz(16)*pol_z(2,0,kg)
                coef_xy(1,17)=coef_xy(1,17)+coef_xyz(17)*pol_z(1,0,kg)
                coef_xy(2,17)=coef_xy(2,17)+coef_xyz(17)*pol_z(2,0,kg)
                coef_xy(1,18)=coef_xy(1,18)+coef_xyz(18)*pol_z(1,0,kg)
                coef_xy(2,18)=coef_xy(2,18)+coef_xyz(18)*pol_z(2,0,kg)
                coef_xy(1,19)=coef_xy(1,19)+coef_xyz(19)*pol_z(1,0,kg)
                coef_xy(2,19)=coef_xy(2,19)+coef_xyz(19)*pol_z(2,0,kg)
                coef_xy(1,20)=coef_xy(1,20)+coef_xyz(20)*pol_z(1,0,kg)
                coef_xy(2,20)=coef_xy(2,20)+coef_xyz(20)*pol_z(2,0,kg)
                coef_xy(1,21)=coef_xy(1,21)+coef_xyz(21)*pol_z(1,0,kg)
                coef_xy(2,21)=coef_xy(2,21)+coef_xyz(21)*pol_z(2,0,kg)
                coef_xy(1,1)=coef_xy(1,1)+coef_xyz(22)*pol_z(1,1,kg)
                coef_xy(2,1)=coef_xy(2,1)+coef_xyz(22)*pol_z(2,1,kg)
                coef_xy(1,2)=coef_xy(1,2)+coef_xyz(23)*pol_z(1,1,kg)
                coef_xy(2,2)=coef_xy(2,2)+coef_xyz(23)*pol_z(2,1,kg)
                coef_xy(1,3)=coef_xy(1,3)+coef_xyz(24)*pol_z(1,1,kg)
                coef_xy(2,3)=coef_xy(2,3)+coef_xyz(24)*pol_z(2,1,kg)
                coef_xy(1,4)=coef_xy(1,4)+coef_xyz(25)*pol_z(1,1,kg)
                coef_xy(2,4)=coef_xy(2,4)+coef_xyz(25)*pol_z(2,1,kg)
                coef_xy(1,5)=coef_xy(1,5)+coef_xyz(26)*pol_z(1,1,kg)
                coef_xy(2,5)=coef_xy(2,5)+coef_xyz(26)*pol_z(2,1,kg)
                coef_xy(1,7)=coef_xy(1,7)+coef_xyz(27)*pol_z(1,1,kg)
                coef_xy(2,7)=coef_xy(2,7)+coef_xyz(27)*pol_z(2,1,kg)
                coef_xy(1,8)=coef_xy(1,8)+coef_xyz(28)*pol_z(1,1,kg)
                coef_xy(2,8)=coef_xy(2,8)+coef_xyz(28)*pol_z(2,1,kg)
                coef_xy(1,9)=coef_xy(1,9)+coef_xyz(29)*pol_z(1,1,kg)
                coef_xy(2,9)=coef_xy(2,9)+coef_xyz(29)*pol_z(2,1,kg)
                coef_xy(1,10)=coef_xy(1,10)+coef_xyz(30)*pol_z(1,1,kg)
                coef_xy(2,10)=coef_xy(2,10)+coef_xyz(30)*pol_z(2,1,kg)
                coef_xy(1,12)=coef_xy(1,12)+coef_xyz(31)*pol_z(1,1,kg)
                coef_xy(2,12)=coef_xy(2,12)+coef_xyz(31)*pol_z(2,1,kg)
                coef_xy(1,13)=coef_xy(1,13)+coef_xyz(32)*pol_z(1,1,kg)
                coef_xy(2,13)=coef_xy(2,13)+coef_xyz(32)*pol_z(2,1,kg)
                coef_xy(1,14)=coef_xy(1,14)+coef_xyz(33)*pol_z(1,1,kg)
                coef_xy(2,14)=coef_xy(2,14)+coef_xyz(33)*pol_z(2,1,kg)
                coef_xy(1,16)=coef_xy(1,16)+coef_xyz(34)*pol_z(1,1,kg)
                coef_xy(2,16)=coef_xy(2,16)+coef_xyz(34)*pol_z(2,1,kg)
                coef_xy(1,17)=coef_xy(1,17)+coef_xyz(35)*pol_z(1,1,kg)
                coef_xy(2,17)=coef_xy(2,17)+coef_xyz(35)*pol_z(2,1,kg)
                coef_xy(1,19)=coef_xy(1,19)+coef_xyz(36)*pol_z(1,1,kg)
                coef_xy(2,19)=coef_xy(2,19)+coef_xyz(36)*pol_z(2,1,kg)
                coef_xy(1,1)=coef_xy(1,1)+coef_xyz(37)*pol_z(1,2,kg)
                coef_xy(2,1)=coef_xy(2,1)+coef_xyz(37)*pol_z(2,2,kg)
                coef_xy(1,2)=coef_xy(1,2)+coef_xyz(38)*pol_z(1,2,kg)
                coef_xy(2,2)=coef_xy(2,2)+coef_xyz(38)*pol_z(2,2,kg)
                coef_xy(1,3)=coef_xy(1,3)+coef_xyz(39)*pol_z(1,2,kg)
                coef_xy(2,3)=coef_xy(2,3)+coef_xyz(39)*pol_z(2,2,kg)
                coef_xy(1,4)=coef_xy(1,4)+coef_xyz(40)*pol_z(1,2,kg)
                coef_xy(2,4)=coef_xy(2,4)+coef_xyz(40)*pol_z(2,2,kg)
                coef_xy(1,7)=coef_xy(1,7)+coef_xyz(41)*pol_z(1,2,kg)
                coef_xy(2,7)=coef_xy(2,7)+coef_xyz(41)*pol_z(2,2,kg)
                coef_xy(1,8)=coef_xy(1,8)+coef_xyz(42)*pol_z(1,2,kg)
                coef_xy(2,8)=coef_xy(2,8)+coef_xyz(42)*pol_z(2,2,kg)
                coef_xy(1,9)=coef_xy(1,9)+coef_xyz(43)*pol_z(1,2,kg)
                coef_xy(2,9)=coef_xy(2,9)+coef_xyz(43)*pol_z(2,2,kg)
                coef_xy(1,12)=coef_xy(1,12)+coef_xyz(44)*pol_z(1,2,kg)
                coef_xy(2,12)=coef_xy(2,12)+coef_xyz(44)*pol_z(2,2,kg)
                coef_xy(1,13)=coef_xy(1,13)+coef_xyz(45)*pol_z(1,2,kg)
                coef_xy(2,13)=coef_xy(2,13)+coef_xyz(45)*pol_z(2,2,kg)
                coef_xy(1,16)=coef_xy(1,16)+coef_xyz(46)*pol_z(1,2,kg)
                coef_xy(2,16)=coef_xy(2,16)+coef_xyz(46)*pol_z(2,2,kg)
                coef_xy(1,1)=coef_xy(1,1)+coef_xyz(47)*pol_z(1,3,kg)
                coef_xy(2,1)=coef_xy(2,1)+coef_xyz(47)*pol_z(2,3,kg)
                coef_xy(1,2)=coef_xy(1,2)+coef_xyz(48)*pol_z(1,3,kg)
                coef_xy(2,2)=coef_xy(2,2)+coef_xyz(48)*pol_z(2,3,kg)
                coef_xy(1,3)=coef_xy(1,3)+coef_xyz(49)*pol_z(1,3,kg)
                coef_xy(2,3)=coef_xy(2,3)+coef_xyz(49)*pol_z(2,3,kg)
                coef_xy(1,7)=coef_xy(1,7)+coef_xyz(50)*pol_z(1,3,kg)
                coef_xy(2,7)=coef_xy(2,7)+coef_xyz(50)*pol_z(2,3,kg)
                coef_xy(1,8)=coef_xy(1,8)+coef_xyz(51)*pol_z(1,3,kg)
                coef_xy(2,8)=coef_xy(2,8)+coef_xyz(51)*pol_z(2,3,kg)
                coef_xy(1,12)=coef_xy(1,12)+coef_xyz(52)*pol_z(1,3,kg)
                coef_xy(2,12)=coef_xy(2,12)+coef_xyz(52)*pol_z(2,3,kg)
                coef_xy(1,1)=coef_xy(1,1)+coef_xyz(53)*pol_z(1,4,kg)
                coef_xy(2,1)=coef_xy(2,1)+coef_xyz(53)*pol_z(2,4,kg)
                coef_xy(1,2)=coef_xy(1,2)+coef_xyz(54)*pol_z(1,4,kg)
                coef_xy(2,2)=coef_xy(2,2)+coef_xyz(54)*pol_z(2,4,kg)
                coef_xy(1,7)=coef_xy(1,7)+coef_xyz(55)*pol_z(1,4,kg)
                coef_xy(2,7)=coef_xy(2,7)+coef_xyz(55)*pol_z(2,4,kg)
                coef_xy(1,1)=coef_xy(1,1)+coef_xyz(56)*pol_z(1,5,kg)
                coef_xy(2,1)=coef_xy(2,1)+coef_xyz(56)*pol_z(2,5,kg)
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
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,1)*pol_y(1,0,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,1)*pol_y(2,0,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,2)*pol_y(1,0,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,2)*pol_y(2,0,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,3)*pol_y(1,0,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,3)*pol_y(2,0,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,4)*pol_y(1,0,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,4)*pol_y(2,0,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,5)*pol_y(1,0,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,5)*pol_y(2,0,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,6)*pol_y(1,0,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,6)*pol_y(2,0,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,7)*pol_y(1,1,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,7)*pol_y(2,1,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,8)*pol_y(1,1,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,8)*pol_y(2,1,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,9)*pol_y(1,1,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,9)*pol_y(2,1,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,10)*pol_y(1,1,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,10)*pol_y(2,1,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,11)*pol_y(1,1,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,11)*pol_y(2,1,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,12)*pol_y(1,2,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,12)*pol_y(2,2,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,13)*pol_y(1,2,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,13)*pol_y(2,2,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,14)*pol_y(1,2,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,14)*pol_y(2,2,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,15)*pol_y(1,2,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,15)*pol_y(2,2,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,16)*pol_y(1,3,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,16)*pol_y(2,3,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,17)*pol_y(1,3,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,17)*pol_y(2,3,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,18)*pol_y(1,3,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,18)*pol_y(2,3,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,19)*pol_y(1,4,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,19)*pol_y(2,4,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,20)*pol_y(1,4,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,20)*pol_y(2,4,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,21)*pol_y(1,5,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,21)*pol_y(2,5,jg)
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

  END SUBROUTINE collocate_core_5
  SUBROUTINE collocate_core_6(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)
    USE lib_kinds,                        ONLY: wp
    INTEGER, INTENT(IN)                      :: sphere_bounds(*), cmax, &
                                                map(-cmax:cmax,1:3), &
                                                gridbounds(2,3)
    REAL(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), &
      gridbounds(1,2):gridbounds(2,2), gridbounds(1,3):gridbounds(2,3))
    INTEGER, PARAMETER                       :: lp = 6
    REAL(wp), INTENT(IN) :: pol_x(0:lp,-cmax:cmax), pol_y(1:2,0:lp,-cmax:0), &
      pol_z(1:2,0:lp,-cmax:0), coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)

    INTEGER                                  :: i, ig, igmax, igmin, j, j2, &
                                                jg, jg2, jgmin, k, k2, kg, &
                                                kg2, kgmin, lxp, sci
    REAL(wp)                                 :: coef_x(4,0:lp), &
                                                coef_xy(2,(lp+1)*(lp+2)/2), &
                                                s01, s02, s03, s04

    sci=1

    kgmin=sphere_bounds(sci)
    sci=sci+1
    DO kg=kgmin,0
       kg2=1-kg
       k=map(kg,3)
       k2=map(kg2,3)
       coef_xy=0.0_wp
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(1)*pol_z(:,0,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(2)*pol_z(:,0,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(3)*pol_z(:,0,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(4)*pol_z(:,0,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(5)*pol_z(:,0,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(6)*pol_z(:,0,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(7)*pol_z(:,0,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(8)*pol_z(:,0,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(9)*pol_z(:,0,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(10)*pol_z(:,0,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(11)*pol_z(:,0,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(12)*pol_z(:,0,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(13)*pol_z(:,0,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(14)*pol_z(:,0,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(15)*pol_z(:,0,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(16)*pol_z(:,0,kg)
                coef_xy(:,17)=coef_xy(:,17)+coef_xyz(17)*pol_z(:,0,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(18)*pol_z(:,0,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(19)*pol_z(:,0,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(20)*pol_z(:,0,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(21)*pol_z(:,0,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(22)*pol_z(:,0,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(23)*pol_z(:,0,kg)
                coef_xy(:,24)=coef_xy(:,24)+coef_xyz(24)*pol_z(:,0,kg)
                coef_xy(:,25)=coef_xy(:,25)+coef_xyz(25)*pol_z(:,0,kg)
                coef_xy(:,26)=coef_xy(:,26)+coef_xyz(26)*pol_z(:,0,kg)
                coef_xy(:,27)=coef_xy(:,27)+coef_xyz(27)*pol_z(:,0,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(28)*pol_z(:,0,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(29)*pol_z(:,1,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(30)*pol_z(:,1,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(31)*pol_z(:,1,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(32)*pol_z(:,1,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(33)*pol_z(:,1,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(34)*pol_z(:,1,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(35)*pol_z(:,1,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(36)*pol_z(:,1,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(37)*pol_z(:,1,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(38)*pol_z(:,1,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(39)*pol_z(:,1,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(40)*pol_z(:,1,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(41)*pol_z(:,1,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(42)*pol_z(:,1,kg)
                coef_xy(:,17)=coef_xy(:,17)+coef_xyz(43)*pol_z(:,1,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(44)*pol_z(:,1,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(45)*pol_z(:,1,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(46)*pol_z(:,1,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(47)*pol_z(:,1,kg)
                coef_xy(:,24)=coef_xy(:,24)+coef_xyz(48)*pol_z(:,1,kg)
                coef_xy(:,26)=coef_xy(:,26)+coef_xyz(49)*pol_z(:,1,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(50)*pol_z(:,2,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(51)*pol_z(:,2,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(52)*pol_z(:,2,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(53)*pol_z(:,2,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(54)*pol_z(:,2,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(55)*pol_z(:,2,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(56)*pol_z(:,2,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(57)*pol_z(:,2,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(58)*pol_z(:,2,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(59)*pol_z(:,2,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(60)*pol_z(:,2,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(61)*pol_z(:,2,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(62)*pol_z(:,2,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(63)*pol_z(:,2,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(64)*pol_z(:,2,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(65)*pol_z(:,3,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(66)*pol_z(:,3,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(67)*pol_z(:,3,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(68)*pol_z(:,3,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(69)*pol_z(:,3,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(70)*pol_z(:,3,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(71)*pol_z(:,3,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(72)*pol_z(:,3,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(73)*pol_z(:,3,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(74)*pol_z(:,3,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(75)*pol_z(:,4,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(76)*pol_z(:,4,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(77)*pol_z(:,4,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(78)*pol_z(:,4,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(79)*pol_z(:,4,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(80)*pol_z(:,4,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(81)*pol_z(:,5,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(82)*pol_z(:,5,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(83)*pol_z(:,5,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(84)*pol_z(:,6,kg)
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
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,1)*pol_y(1,0,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,1)*pol_y(2,0,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,2)*pol_y(1,0,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,2)*pol_y(2,0,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,3)*pol_y(1,0,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,3)*pol_y(2,0,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,4)*pol_y(1,0,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,4)*pol_y(2,0,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,5)*pol_y(1,0,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,5)*pol_y(2,0,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,6)*pol_y(1,0,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,6)*pol_y(2,0,jg)
             coef_x(1:2,6)=coef_x(1:2,6)+coef_xy(1:2,7)*pol_y(1,0,jg)
             coef_x(3:4,6)=coef_x(3:4,6)+coef_xy(1:2,7)*pol_y(2,0,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,8)*pol_y(1,1,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,8)*pol_y(2,1,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,9)*pol_y(1,1,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,9)*pol_y(2,1,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,10)*pol_y(1,1,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,10)*pol_y(2,1,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,11)*pol_y(1,1,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,11)*pol_y(2,1,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,12)*pol_y(1,1,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,12)*pol_y(2,1,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,13)*pol_y(1,1,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,13)*pol_y(2,1,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,14)*pol_y(1,2,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,14)*pol_y(2,2,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,15)*pol_y(1,2,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,15)*pol_y(2,2,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,16)*pol_y(1,2,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,16)*pol_y(2,2,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,17)*pol_y(1,2,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,17)*pol_y(2,2,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,18)*pol_y(1,2,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,18)*pol_y(2,2,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,19)*pol_y(1,3,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,19)*pol_y(2,3,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,20)*pol_y(1,3,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,20)*pol_y(2,3,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,21)*pol_y(1,3,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,21)*pol_y(2,3,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,22)*pol_y(1,3,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,22)*pol_y(2,3,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,23)*pol_y(1,4,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,23)*pol_y(2,4,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,24)*pol_y(1,4,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,24)*pol_y(2,4,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,25)*pol_y(1,4,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,25)*pol_y(2,4,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,26)*pol_y(1,5,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,26)*pol_y(2,5,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,27)*pol_y(1,5,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,27)*pol_y(2,5,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,28)*pol_y(1,6,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,28)*pol_y(2,6,jg)
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

  END SUBROUTINE collocate_core_6
  SUBROUTINE collocate_core_7(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)
    USE lib_kinds,                        ONLY: wp
    INTEGER, INTENT(IN)                      :: sphere_bounds(*), cmax, &
                                                map(-cmax:cmax,1:3), &
                                                gridbounds(2,3)
    REAL(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), &
      gridbounds(1,2):gridbounds(2,2), gridbounds(1,3):gridbounds(2,3))
    INTEGER, PARAMETER                       :: lp = 7
    REAL(wp), INTENT(IN) :: pol_x(0:lp,-cmax:cmax), pol_y(1:2,0:lp,-cmax:0), &
      pol_z(1:2,0:lp,-cmax:0), coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)

    INTEGER                                  :: i, ig, igmax, igmin, j, j2, &
                                                jg, jg2, jgmin, k, k2, kg, &
                                                kg2, kgmin, sci
    REAL(wp)                                 :: coef_x(4,0:lp), &
                                                coef_xy(2,(lp+1)*(lp+2)/2), &
                                                s01, s02, s03, s04

    sci=1

    kgmin=sphere_bounds(sci)
    sci=sci+1
    DO kg=kgmin,0
       kg2=1-kg
       k=map(kg,3)
       k2=map(kg2,3)
       coef_xy=0.0_wp
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(1)*pol_z(:,0,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(2)*pol_z(:,0,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(3)*pol_z(:,0,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(4)*pol_z(:,0,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(5)*pol_z(:,0,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(6)*pol_z(:,0,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(7)*pol_z(:,0,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(8)*pol_z(:,0,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(9)*pol_z(:,0,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(10)*pol_z(:,0,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(11)*pol_z(:,0,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(12)*pol_z(:,0,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(13)*pol_z(:,0,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(14)*pol_z(:,0,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(15)*pol_z(:,0,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(16)*pol_z(:,0,kg)
                coef_xy(:,17)=coef_xy(:,17)+coef_xyz(17)*pol_z(:,0,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(18)*pol_z(:,0,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(19)*pol_z(:,0,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(20)*pol_z(:,0,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(21)*pol_z(:,0,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(22)*pol_z(:,0,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(23)*pol_z(:,0,kg)
                coef_xy(:,24)=coef_xy(:,24)+coef_xyz(24)*pol_z(:,0,kg)
                coef_xy(:,25)=coef_xy(:,25)+coef_xyz(25)*pol_z(:,0,kg)
                coef_xy(:,26)=coef_xy(:,26)+coef_xyz(26)*pol_z(:,0,kg)
                coef_xy(:,27)=coef_xy(:,27)+coef_xyz(27)*pol_z(:,0,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(28)*pol_z(:,0,kg)
                coef_xy(:,29)=coef_xy(:,29)+coef_xyz(29)*pol_z(:,0,kg)
                coef_xy(:,30)=coef_xy(:,30)+coef_xyz(30)*pol_z(:,0,kg)
                coef_xy(:,31)=coef_xy(:,31)+coef_xyz(31)*pol_z(:,0,kg)
                coef_xy(:,32)=coef_xy(:,32)+coef_xyz(32)*pol_z(:,0,kg)
                coef_xy(:,33)=coef_xy(:,33)+coef_xyz(33)*pol_z(:,0,kg)
                coef_xy(:,34)=coef_xy(:,34)+coef_xyz(34)*pol_z(:,0,kg)
                coef_xy(:,35)=coef_xy(:,35)+coef_xyz(35)*pol_z(:,0,kg)
                coef_xy(:,36)=coef_xy(:,36)+coef_xyz(36)*pol_z(:,0,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(37)*pol_z(:,1,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(38)*pol_z(:,1,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(39)*pol_z(:,1,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(40)*pol_z(:,1,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(41)*pol_z(:,1,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(42)*pol_z(:,1,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(43)*pol_z(:,1,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(44)*pol_z(:,1,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(45)*pol_z(:,1,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(46)*pol_z(:,1,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(47)*pol_z(:,1,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(48)*pol_z(:,1,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(49)*pol_z(:,1,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(50)*pol_z(:,1,kg)
                coef_xy(:,17)=coef_xy(:,17)+coef_xyz(51)*pol_z(:,1,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(52)*pol_z(:,1,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(53)*pol_z(:,1,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(54)*pol_z(:,1,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(55)*pol_z(:,1,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(56)*pol_z(:,1,kg)
                coef_xy(:,24)=coef_xy(:,24)+coef_xyz(57)*pol_z(:,1,kg)
                coef_xy(:,25)=coef_xy(:,25)+coef_xyz(58)*pol_z(:,1,kg)
                coef_xy(:,27)=coef_xy(:,27)+coef_xyz(59)*pol_z(:,1,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(60)*pol_z(:,1,kg)
                coef_xy(:,29)=coef_xy(:,29)+coef_xyz(61)*pol_z(:,1,kg)
                coef_xy(:,31)=coef_xy(:,31)+coef_xyz(62)*pol_z(:,1,kg)
                coef_xy(:,32)=coef_xy(:,32)+coef_xyz(63)*pol_z(:,1,kg)
                coef_xy(:,34)=coef_xy(:,34)+coef_xyz(64)*pol_z(:,1,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(65)*pol_z(:,2,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(66)*pol_z(:,2,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(67)*pol_z(:,2,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(68)*pol_z(:,2,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(69)*pol_z(:,2,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(70)*pol_z(:,2,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(71)*pol_z(:,2,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(72)*pol_z(:,2,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(73)*pol_z(:,2,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(74)*pol_z(:,2,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(75)*pol_z(:,2,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(76)*pol_z(:,2,kg)
                coef_xy(:,17)=coef_xy(:,17)+coef_xyz(77)*pol_z(:,2,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(78)*pol_z(:,2,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(79)*pol_z(:,2,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(80)*pol_z(:,2,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(81)*pol_z(:,2,kg)
                coef_xy(:,24)=coef_xy(:,24)+coef_xyz(82)*pol_z(:,2,kg)
                coef_xy(:,27)=coef_xy(:,27)+coef_xyz(83)*pol_z(:,2,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(84)*pol_z(:,2,kg)
                coef_xy(:,31)=coef_xy(:,31)+coef_xyz(85)*pol_z(:,2,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(86)*pol_z(:,3,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(87)*pol_z(:,3,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(88)*pol_z(:,3,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(89)*pol_z(:,3,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(90)*pol_z(:,3,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(91)*pol_z(:,3,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(92)*pol_z(:,3,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(93)*pol_z(:,3,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(94)*pol_z(:,3,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(95)*pol_z(:,3,kg)
                coef_xy(:,17)=coef_xy(:,17)+coef_xyz(96)*pol_z(:,3,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(97)*pol_z(:,3,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(98)*pol_z(:,3,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(99)*pol_z(:,3,kg)
                coef_xy(:,27)=coef_xy(:,27)+coef_xyz(100)*pol_z(:,3,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(101)*pol_z(:,4,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(102)*pol_z(:,4,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(103)*pol_z(:,4,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(104)*pol_z(:,4,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(105)*pol_z(:,4,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(106)*pol_z(:,4,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(107)*pol_z(:,4,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(108)*pol_z(:,4,kg)
                coef_xy(:,17)=coef_xy(:,17)+coef_xyz(109)*pol_z(:,4,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(110)*pol_z(:,4,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(111)*pol_z(:,5,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(112)*pol_z(:,5,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(113)*pol_z(:,5,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(114)*pol_z(:,5,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(115)*pol_z(:,5,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(116)*pol_z(:,5,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(117)*pol_z(:,6,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(118)*pol_z(:,6,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(119)*pol_z(:,6,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(120)*pol_z(:,7,kg)
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
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,1)*pol_y(1,0,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,1)*pol_y(2,0,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,2)*pol_y(1,0,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,2)*pol_y(2,0,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,3)*pol_y(1,0,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,3)*pol_y(2,0,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,4)*pol_y(1,0,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,4)*pol_y(2,0,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,5)*pol_y(1,0,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,5)*pol_y(2,0,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,6)*pol_y(1,0,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,6)*pol_y(2,0,jg)
             coef_x(1:2,6)=coef_x(1:2,6)+coef_xy(1:2,7)*pol_y(1,0,jg)
             coef_x(3:4,6)=coef_x(3:4,6)+coef_xy(1:2,7)*pol_y(2,0,jg)
             coef_x(1:2,7)=coef_x(1:2,7)+coef_xy(1:2,8)*pol_y(1,0,jg)
             coef_x(3:4,7)=coef_x(3:4,7)+coef_xy(1:2,8)*pol_y(2,0,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,9)*pol_y(1,1,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,9)*pol_y(2,1,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,10)*pol_y(1,1,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,10)*pol_y(2,1,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,11)*pol_y(1,1,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,11)*pol_y(2,1,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,12)*pol_y(1,1,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,12)*pol_y(2,1,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,13)*pol_y(1,1,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,13)*pol_y(2,1,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,14)*pol_y(1,1,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,14)*pol_y(2,1,jg)
             coef_x(1:2,6)=coef_x(1:2,6)+coef_xy(1:2,15)*pol_y(1,1,jg)
             coef_x(3:4,6)=coef_x(3:4,6)+coef_xy(1:2,15)*pol_y(2,1,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,16)*pol_y(1,2,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,16)*pol_y(2,2,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,17)*pol_y(1,2,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,17)*pol_y(2,2,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,18)*pol_y(1,2,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,18)*pol_y(2,2,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,19)*pol_y(1,2,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,19)*pol_y(2,2,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,20)*pol_y(1,2,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,20)*pol_y(2,2,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,21)*pol_y(1,2,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,21)*pol_y(2,2,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,22)*pol_y(1,3,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,22)*pol_y(2,3,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,23)*pol_y(1,3,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,23)*pol_y(2,3,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,24)*pol_y(1,3,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,24)*pol_y(2,3,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,25)*pol_y(1,3,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,25)*pol_y(2,3,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,26)*pol_y(1,3,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,26)*pol_y(2,3,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,27)*pol_y(1,4,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,27)*pol_y(2,4,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,28)*pol_y(1,4,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,28)*pol_y(2,4,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,29)*pol_y(1,4,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,29)*pol_y(2,4,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,30)*pol_y(1,4,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,30)*pol_y(2,4,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,31)*pol_y(1,5,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,31)*pol_y(2,5,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,32)*pol_y(1,5,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,32)*pol_y(2,5,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,33)*pol_y(1,5,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,33)*pol_y(2,5,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,34)*pol_y(1,6,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,34)*pol_y(2,6,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,35)*pol_y(1,6,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,35)*pol_y(2,6,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,36)*pol_y(1,7,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,36)*pol_y(2,7,jg)
          DO ig=igmin,igmax
             i=map(ig,1)
             s01=0.0_wp
             s02=0.0_wp
             s03=0.0_wp
             s04=0.0_wp
                s01=s01+coef_x(1,0)*pol_x(0,ig)
                s02=s02+coef_x(2,0)*pol_x(0,ig)
                s03=s03+coef_x(3,0)*pol_x(0,ig)
                s04=s04+coef_x(4,0)*pol_x(0,ig)
                s01=s01+coef_x(1,1)*pol_x(1,ig)
                s02=s02+coef_x(2,1)*pol_x(1,ig)
                s03=s03+coef_x(3,1)*pol_x(1,ig)
                s04=s04+coef_x(4,1)*pol_x(1,ig)
                s01=s01+coef_x(1,2)*pol_x(2,ig)
                s02=s02+coef_x(2,2)*pol_x(2,ig)
                s03=s03+coef_x(3,2)*pol_x(2,ig)
                s04=s04+coef_x(4,2)*pol_x(2,ig)
                s01=s01+coef_x(1,3)*pol_x(3,ig)
                s02=s02+coef_x(2,3)*pol_x(3,ig)
                s03=s03+coef_x(3,3)*pol_x(3,ig)
                s04=s04+coef_x(4,3)*pol_x(3,ig)
                s01=s01+coef_x(1,4)*pol_x(4,ig)
                s02=s02+coef_x(2,4)*pol_x(4,ig)
                s03=s03+coef_x(3,4)*pol_x(4,ig)
                s04=s04+coef_x(4,4)*pol_x(4,ig)
                s01=s01+coef_x(1,5)*pol_x(5,ig)
                s02=s02+coef_x(2,5)*pol_x(5,ig)
                s03=s03+coef_x(3,5)*pol_x(5,ig)
                s04=s04+coef_x(4,5)*pol_x(5,ig)
                s01=s01+coef_x(1,6)*pol_x(6,ig)
                s02=s02+coef_x(2,6)*pol_x(6,ig)
                s03=s03+coef_x(3,6)*pol_x(6,ig)
                s04=s04+coef_x(4,6)*pol_x(6,ig)
                s01=s01+coef_x(1,7)*pol_x(7,ig)
                s02=s02+coef_x(2,7)*pol_x(7,ig)
                s03=s03+coef_x(3,7)*pol_x(7,ig)
                s04=s04+coef_x(4,7)*pol_x(7,ig)
             grid(i,j,k) = grid(i,j,k)     + s01
             grid(i,j2,k) = grid(i,j2,k)   + s03
             grid(i,j,k2) = grid(i,j,k2)   + s02
             grid(i,j2,k2) = grid(i,j2,k2) + s04
          END DO
       END DO
    END DO

  END SUBROUTINE collocate_core_7
  SUBROUTINE collocate_core_8(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)
    USE lib_kinds,                        ONLY: wp
    INTEGER, INTENT(IN)                      :: sphere_bounds(*), cmax, &
                                                map(-cmax:cmax,1:3), &
                                                gridbounds(2,3)
    REAL(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), &
      gridbounds(1,2):gridbounds(2,2), gridbounds(1,3):gridbounds(2,3))
    INTEGER, PARAMETER                       :: lp = 8
    REAL(wp), INTENT(IN) :: pol_x(0:lp,-cmax:cmax), pol_y(1:2,0:lp,-cmax:0), &
      pol_z(1:2,0:lp,-cmax:0), coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)

    INTEGER                                  :: i, ig, igmax, igmin, j, j2, &
                                                jg, jg2, jgmin, k, k2, kg, &
                                                kg2, kgmin, sci
    REAL(wp)                                 :: coef_x(4,0:lp), &
                                                coef_xy(2,(lp+1)*(lp+2)/2), &
                                                s01, s02, s03, s04

    sci=1

    kgmin=sphere_bounds(sci)
    sci=sci+1
    DO kg=kgmin,0
       kg2=1-kg
       k=map(kg,3)
       k2=map(kg2,3)
       coef_xy=0.0_wp
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(1)*pol_z(:,0,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(2)*pol_z(:,0,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(3)*pol_z(:,0,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(4)*pol_z(:,0,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(5)*pol_z(:,0,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(6)*pol_z(:,0,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(7)*pol_z(:,0,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(8)*pol_z(:,0,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(9)*pol_z(:,0,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(10)*pol_z(:,0,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(11)*pol_z(:,0,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(12)*pol_z(:,0,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(13)*pol_z(:,0,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(14)*pol_z(:,0,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(15)*pol_z(:,0,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(16)*pol_z(:,0,kg)
                coef_xy(:,17)=coef_xy(:,17)+coef_xyz(17)*pol_z(:,0,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(18)*pol_z(:,0,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(19)*pol_z(:,0,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(20)*pol_z(:,0,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(21)*pol_z(:,0,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(22)*pol_z(:,0,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(23)*pol_z(:,0,kg)
                coef_xy(:,24)=coef_xy(:,24)+coef_xyz(24)*pol_z(:,0,kg)
                coef_xy(:,25)=coef_xy(:,25)+coef_xyz(25)*pol_z(:,0,kg)
                coef_xy(:,26)=coef_xy(:,26)+coef_xyz(26)*pol_z(:,0,kg)
                coef_xy(:,27)=coef_xy(:,27)+coef_xyz(27)*pol_z(:,0,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(28)*pol_z(:,0,kg)
                coef_xy(:,29)=coef_xy(:,29)+coef_xyz(29)*pol_z(:,0,kg)
                coef_xy(:,30)=coef_xy(:,30)+coef_xyz(30)*pol_z(:,0,kg)
                coef_xy(:,31)=coef_xy(:,31)+coef_xyz(31)*pol_z(:,0,kg)
                coef_xy(:,32)=coef_xy(:,32)+coef_xyz(32)*pol_z(:,0,kg)
                coef_xy(:,33)=coef_xy(:,33)+coef_xyz(33)*pol_z(:,0,kg)
                coef_xy(:,34)=coef_xy(:,34)+coef_xyz(34)*pol_z(:,0,kg)
                coef_xy(:,35)=coef_xy(:,35)+coef_xyz(35)*pol_z(:,0,kg)
                coef_xy(:,36)=coef_xy(:,36)+coef_xyz(36)*pol_z(:,0,kg)
                coef_xy(:,37)=coef_xy(:,37)+coef_xyz(37)*pol_z(:,0,kg)
                coef_xy(:,38)=coef_xy(:,38)+coef_xyz(38)*pol_z(:,0,kg)
                coef_xy(:,39)=coef_xy(:,39)+coef_xyz(39)*pol_z(:,0,kg)
                coef_xy(:,40)=coef_xy(:,40)+coef_xyz(40)*pol_z(:,0,kg)
                coef_xy(:,41)=coef_xy(:,41)+coef_xyz(41)*pol_z(:,0,kg)
                coef_xy(:,42)=coef_xy(:,42)+coef_xyz(42)*pol_z(:,0,kg)
                coef_xy(:,43)=coef_xy(:,43)+coef_xyz(43)*pol_z(:,0,kg)
                coef_xy(:,44)=coef_xy(:,44)+coef_xyz(44)*pol_z(:,0,kg)
                coef_xy(:,45)=coef_xy(:,45)+coef_xyz(45)*pol_z(:,0,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(46)*pol_z(:,1,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(47)*pol_z(:,1,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(48)*pol_z(:,1,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(49)*pol_z(:,1,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(50)*pol_z(:,1,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(51)*pol_z(:,1,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(52)*pol_z(:,1,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(53)*pol_z(:,1,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(54)*pol_z(:,1,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(55)*pol_z(:,1,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(56)*pol_z(:,1,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(57)*pol_z(:,1,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(58)*pol_z(:,1,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(59)*pol_z(:,1,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(60)*pol_z(:,1,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(61)*pol_z(:,1,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(62)*pol_z(:,1,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(63)*pol_z(:,1,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(64)*pol_z(:,1,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(65)*pol_z(:,1,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(66)*pol_z(:,1,kg)
                coef_xy(:,25)=coef_xy(:,25)+coef_xyz(67)*pol_z(:,1,kg)
                coef_xy(:,26)=coef_xy(:,26)+coef_xyz(68)*pol_z(:,1,kg)
                coef_xy(:,27)=coef_xy(:,27)+coef_xyz(69)*pol_z(:,1,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(70)*pol_z(:,1,kg)
                coef_xy(:,29)=coef_xy(:,29)+coef_xyz(71)*pol_z(:,1,kg)
                coef_xy(:,31)=coef_xy(:,31)+coef_xyz(72)*pol_z(:,1,kg)
                coef_xy(:,32)=coef_xy(:,32)+coef_xyz(73)*pol_z(:,1,kg)
                coef_xy(:,33)=coef_xy(:,33)+coef_xyz(74)*pol_z(:,1,kg)
                coef_xy(:,34)=coef_xy(:,34)+coef_xyz(75)*pol_z(:,1,kg)
                coef_xy(:,36)=coef_xy(:,36)+coef_xyz(76)*pol_z(:,1,kg)
                coef_xy(:,37)=coef_xy(:,37)+coef_xyz(77)*pol_z(:,1,kg)
                coef_xy(:,38)=coef_xy(:,38)+coef_xyz(78)*pol_z(:,1,kg)
                coef_xy(:,40)=coef_xy(:,40)+coef_xyz(79)*pol_z(:,1,kg)
                coef_xy(:,41)=coef_xy(:,41)+coef_xyz(80)*pol_z(:,1,kg)
                coef_xy(:,43)=coef_xy(:,43)+coef_xyz(81)*pol_z(:,1,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(82)*pol_z(:,2,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(83)*pol_z(:,2,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(84)*pol_z(:,2,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(85)*pol_z(:,2,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(86)*pol_z(:,2,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(87)*pol_z(:,2,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(88)*pol_z(:,2,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(89)*pol_z(:,2,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(90)*pol_z(:,2,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(91)*pol_z(:,2,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(92)*pol_z(:,2,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(93)*pol_z(:,2,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(94)*pol_z(:,2,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(95)*pol_z(:,2,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(96)*pol_z(:,2,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(97)*pol_z(:,2,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(98)*pol_z(:,2,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(99)*pol_z(:,2,kg)
                coef_xy(:,25)=coef_xy(:,25)+coef_xyz(100)*pol_z(:,2,kg)
                coef_xy(:,26)=coef_xy(:,26)+coef_xyz(101)*pol_z(:,2,kg)
                coef_xy(:,27)=coef_xy(:,27)+coef_xyz(102)*pol_z(:,2,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(103)*pol_z(:,2,kg)
                coef_xy(:,31)=coef_xy(:,31)+coef_xyz(104)*pol_z(:,2,kg)
                coef_xy(:,32)=coef_xy(:,32)+coef_xyz(105)*pol_z(:,2,kg)
                coef_xy(:,33)=coef_xy(:,33)+coef_xyz(106)*pol_z(:,2,kg)
                coef_xy(:,36)=coef_xy(:,36)+coef_xyz(107)*pol_z(:,2,kg)
                coef_xy(:,37)=coef_xy(:,37)+coef_xyz(108)*pol_z(:,2,kg)
                coef_xy(:,40)=coef_xy(:,40)+coef_xyz(109)*pol_z(:,2,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(110)*pol_z(:,3,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(111)*pol_z(:,3,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(112)*pol_z(:,3,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(113)*pol_z(:,3,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(114)*pol_z(:,3,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(115)*pol_z(:,3,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(116)*pol_z(:,3,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(117)*pol_z(:,3,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(118)*pol_z(:,3,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(119)*pol_z(:,3,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(120)*pol_z(:,3,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(121)*pol_z(:,3,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(122)*pol_z(:,3,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(123)*pol_z(:,3,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(124)*pol_z(:,3,kg)
                coef_xy(:,25)=coef_xy(:,25)+coef_xyz(125)*pol_z(:,3,kg)
                coef_xy(:,26)=coef_xy(:,26)+coef_xyz(126)*pol_z(:,3,kg)
                coef_xy(:,27)=coef_xy(:,27)+coef_xyz(127)*pol_z(:,3,kg)
                coef_xy(:,31)=coef_xy(:,31)+coef_xyz(128)*pol_z(:,3,kg)
                coef_xy(:,32)=coef_xy(:,32)+coef_xyz(129)*pol_z(:,3,kg)
                coef_xy(:,36)=coef_xy(:,36)+coef_xyz(130)*pol_z(:,3,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(131)*pol_z(:,4,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(132)*pol_z(:,4,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(133)*pol_z(:,4,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(134)*pol_z(:,4,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(135)*pol_z(:,4,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(136)*pol_z(:,4,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(137)*pol_z(:,4,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(138)*pol_z(:,4,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(139)*pol_z(:,4,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(140)*pol_z(:,4,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(141)*pol_z(:,4,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(142)*pol_z(:,4,kg)
                coef_xy(:,25)=coef_xy(:,25)+coef_xyz(143)*pol_z(:,4,kg)
                coef_xy(:,26)=coef_xy(:,26)+coef_xyz(144)*pol_z(:,4,kg)
                coef_xy(:,31)=coef_xy(:,31)+coef_xyz(145)*pol_z(:,4,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(146)*pol_z(:,5,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(147)*pol_z(:,5,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(148)*pol_z(:,5,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(149)*pol_z(:,5,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(150)*pol_z(:,5,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(151)*pol_z(:,5,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(152)*pol_z(:,5,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(153)*pol_z(:,5,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(154)*pol_z(:,5,kg)
                coef_xy(:,25)=coef_xy(:,25)+coef_xyz(155)*pol_z(:,5,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(156)*pol_z(:,6,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(157)*pol_z(:,6,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(158)*pol_z(:,6,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(159)*pol_z(:,6,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(160)*pol_z(:,6,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(161)*pol_z(:,6,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(162)*pol_z(:,7,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(163)*pol_z(:,7,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(164)*pol_z(:,7,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(165)*pol_z(:,8,kg)
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
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,1)*pol_y(1,0,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,1)*pol_y(2,0,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,2)*pol_y(1,0,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,2)*pol_y(2,0,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,3)*pol_y(1,0,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,3)*pol_y(2,0,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,4)*pol_y(1,0,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,4)*pol_y(2,0,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,5)*pol_y(1,0,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,5)*pol_y(2,0,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,6)*pol_y(1,0,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,6)*pol_y(2,0,jg)
             coef_x(1:2,6)=coef_x(1:2,6)+coef_xy(1:2,7)*pol_y(1,0,jg)
             coef_x(3:4,6)=coef_x(3:4,6)+coef_xy(1:2,7)*pol_y(2,0,jg)
             coef_x(1:2,7)=coef_x(1:2,7)+coef_xy(1:2,8)*pol_y(1,0,jg)
             coef_x(3:4,7)=coef_x(3:4,7)+coef_xy(1:2,8)*pol_y(2,0,jg)
             coef_x(1:2,8)=coef_x(1:2,8)+coef_xy(1:2,9)*pol_y(1,0,jg)
             coef_x(3:4,8)=coef_x(3:4,8)+coef_xy(1:2,9)*pol_y(2,0,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,10)*pol_y(1,1,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,10)*pol_y(2,1,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,11)*pol_y(1,1,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,11)*pol_y(2,1,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,12)*pol_y(1,1,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,12)*pol_y(2,1,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,13)*pol_y(1,1,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,13)*pol_y(2,1,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,14)*pol_y(1,1,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,14)*pol_y(2,1,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,15)*pol_y(1,1,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,15)*pol_y(2,1,jg)
             coef_x(1:2,6)=coef_x(1:2,6)+coef_xy(1:2,16)*pol_y(1,1,jg)
             coef_x(3:4,6)=coef_x(3:4,6)+coef_xy(1:2,16)*pol_y(2,1,jg)
             coef_x(1:2,7)=coef_x(1:2,7)+coef_xy(1:2,17)*pol_y(1,1,jg)
             coef_x(3:4,7)=coef_x(3:4,7)+coef_xy(1:2,17)*pol_y(2,1,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,18)*pol_y(1,2,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,18)*pol_y(2,2,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,19)*pol_y(1,2,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,19)*pol_y(2,2,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,20)*pol_y(1,2,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,20)*pol_y(2,2,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,21)*pol_y(1,2,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,21)*pol_y(2,2,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,22)*pol_y(1,2,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,22)*pol_y(2,2,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,23)*pol_y(1,2,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,23)*pol_y(2,2,jg)
             coef_x(1:2,6)=coef_x(1:2,6)+coef_xy(1:2,24)*pol_y(1,2,jg)
             coef_x(3:4,6)=coef_x(3:4,6)+coef_xy(1:2,24)*pol_y(2,2,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,25)*pol_y(1,3,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,25)*pol_y(2,3,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,26)*pol_y(1,3,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,26)*pol_y(2,3,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,27)*pol_y(1,3,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,27)*pol_y(2,3,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,28)*pol_y(1,3,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,28)*pol_y(2,3,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,29)*pol_y(1,3,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,29)*pol_y(2,3,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,30)*pol_y(1,3,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,30)*pol_y(2,3,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,31)*pol_y(1,4,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,31)*pol_y(2,4,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,32)*pol_y(1,4,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,32)*pol_y(2,4,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,33)*pol_y(1,4,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,33)*pol_y(2,4,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,34)*pol_y(1,4,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,34)*pol_y(2,4,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,35)*pol_y(1,4,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,35)*pol_y(2,4,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,36)*pol_y(1,5,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,36)*pol_y(2,5,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,37)*pol_y(1,5,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,37)*pol_y(2,5,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,38)*pol_y(1,5,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,38)*pol_y(2,5,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,39)*pol_y(1,5,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,39)*pol_y(2,5,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,40)*pol_y(1,6,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,40)*pol_y(2,6,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,41)*pol_y(1,6,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,41)*pol_y(2,6,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,42)*pol_y(1,6,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,42)*pol_y(2,6,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,43)*pol_y(1,7,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,43)*pol_y(2,7,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,44)*pol_y(1,7,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,44)*pol_y(2,7,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,45)*pol_y(1,8,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,45)*pol_y(2,8,jg)
          DO ig=igmin,igmax
             i=map(ig,1)
             s01=0.0_wp
             s02=0.0_wp
             s03=0.0_wp
             s04=0.0_wp
                s01=s01+coef_x(1,0)*pol_x(0,ig)
                s02=s02+coef_x(2,0)*pol_x(0,ig)
                s03=s03+coef_x(3,0)*pol_x(0,ig)
                s04=s04+coef_x(4,0)*pol_x(0,ig)
                s01=s01+coef_x(1,1)*pol_x(1,ig)
                s02=s02+coef_x(2,1)*pol_x(1,ig)
                s03=s03+coef_x(3,1)*pol_x(1,ig)
                s04=s04+coef_x(4,1)*pol_x(1,ig)
                s01=s01+coef_x(1,2)*pol_x(2,ig)
                s02=s02+coef_x(2,2)*pol_x(2,ig)
                s03=s03+coef_x(3,2)*pol_x(2,ig)
                s04=s04+coef_x(4,2)*pol_x(2,ig)
                s01=s01+coef_x(1,3)*pol_x(3,ig)
                s02=s02+coef_x(2,3)*pol_x(3,ig)
                s03=s03+coef_x(3,3)*pol_x(3,ig)
                s04=s04+coef_x(4,3)*pol_x(3,ig)
                s01=s01+coef_x(1,4)*pol_x(4,ig)
                s02=s02+coef_x(2,4)*pol_x(4,ig)
                s03=s03+coef_x(3,4)*pol_x(4,ig)
                s04=s04+coef_x(4,4)*pol_x(4,ig)
                s01=s01+coef_x(1,5)*pol_x(5,ig)
                s02=s02+coef_x(2,5)*pol_x(5,ig)
                s03=s03+coef_x(3,5)*pol_x(5,ig)
                s04=s04+coef_x(4,5)*pol_x(5,ig)
                s01=s01+coef_x(1,6)*pol_x(6,ig)
                s02=s02+coef_x(2,6)*pol_x(6,ig)
                s03=s03+coef_x(3,6)*pol_x(6,ig)
                s04=s04+coef_x(4,6)*pol_x(6,ig)
                s01=s01+coef_x(1,7)*pol_x(7,ig)
                s02=s02+coef_x(2,7)*pol_x(7,ig)
                s03=s03+coef_x(3,7)*pol_x(7,ig)
                s04=s04+coef_x(4,7)*pol_x(7,ig)
                s01=s01+coef_x(1,8)*pol_x(8,ig)
                s02=s02+coef_x(2,8)*pol_x(8,ig)
                s03=s03+coef_x(3,8)*pol_x(8,ig)
                s04=s04+coef_x(4,8)*pol_x(8,ig)
             grid(i,j,k) = grid(i,j,k)     + s01
             grid(i,j2,k) = grid(i,j2,k)   + s03
             grid(i,j,k2) = grid(i,j,k2)   + s02
             grid(i,j2,k2) = grid(i,j2,k2) + s04
          END DO
       END DO
    END DO

  END SUBROUTINE collocate_core_8
  SUBROUTINE collocate_core_9(grid,coef_xyz,pol_x,pol_y,pol_z,map,sphere_bounds,cmax,gridbounds)
    USE lib_kinds,                        ONLY: wp
    INTEGER, INTENT(IN)                      :: sphere_bounds(*), cmax, &
                                                map(-cmax:cmax,1:3), &
                                                gridbounds(2,3)
    REAL(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1), &
      gridbounds(1,2):gridbounds(2,2), gridbounds(1,3):gridbounds(2,3))
    INTEGER, PARAMETER                       :: lp = 9
    REAL(wp), INTENT(IN) :: pol_x(0:lp,-cmax:cmax), pol_y(1:2,0:lp,-cmax:0), &
      pol_z(1:2,0:lp,-cmax:0), coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)

    INTEGER                                  :: i, ig, igmax, igmin, j, j2, &
                                                jg, jg2, jgmin, k, k2, kg, &
                                                kg2, kgmin, sci
    REAL(wp)                                 :: coef_x(4,0:lp), &
                                                coef_xy(2,(lp+1)*(lp+2)/2), &
                                                s01, s02, s03, s04

    sci=1

    kgmin=sphere_bounds(sci)
    sci=sci+1
    DO kg=kgmin,0
       kg2=1-kg
       k=map(kg,3)
       k2=map(kg2,3)
       coef_xy=0.0_wp
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(1)*pol_z(:,0,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(2)*pol_z(:,0,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(3)*pol_z(:,0,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(4)*pol_z(:,0,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(5)*pol_z(:,0,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(6)*pol_z(:,0,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(7)*pol_z(:,0,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(8)*pol_z(:,0,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(9)*pol_z(:,0,kg)
                coef_xy(:,10)=coef_xy(:,10)+coef_xyz(10)*pol_z(:,0,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(11)*pol_z(:,0,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(12)*pol_z(:,0,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(13)*pol_z(:,0,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(14)*pol_z(:,0,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(15)*pol_z(:,0,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(16)*pol_z(:,0,kg)
                coef_xy(:,17)=coef_xy(:,17)+coef_xyz(17)*pol_z(:,0,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(18)*pol_z(:,0,kg)
                coef_xy(:,19)=coef_xy(:,19)+coef_xyz(19)*pol_z(:,0,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(20)*pol_z(:,0,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(21)*pol_z(:,0,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(22)*pol_z(:,0,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(23)*pol_z(:,0,kg)
                coef_xy(:,24)=coef_xy(:,24)+coef_xyz(24)*pol_z(:,0,kg)
                coef_xy(:,25)=coef_xy(:,25)+coef_xyz(25)*pol_z(:,0,kg)
                coef_xy(:,26)=coef_xy(:,26)+coef_xyz(26)*pol_z(:,0,kg)
                coef_xy(:,27)=coef_xy(:,27)+coef_xyz(27)*pol_z(:,0,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(28)*pol_z(:,0,kg)
                coef_xy(:,29)=coef_xy(:,29)+coef_xyz(29)*pol_z(:,0,kg)
                coef_xy(:,30)=coef_xy(:,30)+coef_xyz(30)*pol_z(:,0,kg)
                coef_xy(:,31)=coef_xy(:,31)+coef_xyz(31)*pol_z(:,0,kg)
                coef_xy(:,32)=coef_xy(:,32)+coef_xyz(32)*pol_z(:,0,kg)
                coef_xy(:,33)=coef_xy(:,33)+coef_xyz(33)*pol_z(:,0,kg)
                coef_xy(:,34)=coef_xy(:,34)+coef_xyz(34)*pol_z(:,0,kg)
                coef_xy(:,35)=coef_xy(:,35)+coef_xyz(35)*pol_z(:,0,kg)
                coef_xy(:,36)=coef_xy(:,36)+coef_xyz(36)*pol_z(:,0,kg)
                coef_xy(:,37)=coef_xy(:,37)+coef_xyz(37)*pol_z(:,0,kg)
                coef_xy(:,38)=coef_xy(:,38)+coef_xyz(38)*pol_z(:,0,kg)
                coef_xy(:,39)=coef_xy(:,39)+coef_xyz(39)*pol_z(:,0,kg)
                coef_xy(:,40)=coef_xy(:,40)+coef_xyz(40)*pol_z(:,0,kg)
                coef_xy(:,41)=coef_xy(:,41)+coef_xyz(41)*pol_z(:,0,kg)
                coef_xy(:,42)=coef_xy(:,42)+coef_xyz(42)*pol_z(:,0,kg)
                coef_xy(:,43)=coef_xy(:,43)+coef_xyz(43)*pol_z(:,0,kg)
                coef_xy(:,44)=coef_xy(:,44)+coef_xyz(44)*pol_z(:,0,kg)
                coef_xy(:,45)=coef_xy(:,45)+coef_xyz(45)*pol_z(:,0,kg)
                coef_xy(:,46)=coef_xy(:,46)+coef_xyz(46)*pol_z(:,0,kg)
                coef_xy(:,47)=coef_xy(:,47)+coef_xyz(47)*pol_z(:,0,kg)
                coef_xy(:,48)=coef_xy(:,48)+coef_xyz(48)*pol_z(:,0,kg)
                coef_xy(:,49)=coef_xy(:,49)+coef_xyz(49)*pol_z(:,0,kg)
                coef_xy(:,50)=coef_xy(:,50)+coef_xyz(50)*pol_z(:,0,kg)
                coef_xy(:,51)=coef_xy(:,51)+coef_xyz(51)*pol_z(:,0,kg)
                coef_xy(:,52)=coef_xy(:,52)+coef_xyz(52)*pol_z(:,0,kg)
                coef_xy(:,53)=coef_xy(:,53)+coef_xyz(53)*pol_z(:,0,kg)
                coef_xy(:,54)=coef_xy(:,54)+coef_xyz(54)*pol_z(:,0,kg)
                coef_xy(:,55)=coef_xy(:,55)+coef_xyz(55)*pol_z(:,0,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(56)*pol_z(:,1,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(57)*pol_z(:,1,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(58)*pol_z(:,1,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(59)*pol_z(:,1,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(60)*pol_z(:,1,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(61)*pol_z(:,1,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(62)*pol_z(:,1,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(63)*pol_z(:,1,kg)
                coef_xy(:,9)=coef_xy(:,9)+coef_xyz(64)*pol_z(:,1,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(65)*pol_z(:,1,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(66)*pol_z(:,1,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(67)*pol_z(:,1,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(68)*pol_z(:,1,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(69)*pol_z(:,1,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(70)*pol_z(:,1,kg)
                coef_xy(:,17)=coef_xy(:,17)+coef_xyz(71)*pol_z(:,1,kg)
                coef_xy(:,18)=coef_xy(:,18)+coef_xyz(72)*pol_z(:,1,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(73)*pol_z(:,1,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(74)*pol_z(:,1,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(75)*pol_z(:,1,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(76)*pol_z(:,1,kg)
                coef_xy(:,24)=coef_xy(:,24)+coef_xyz(77)*pol_z(:,1,kg)
                coef_xy(:,25)=coef_xy(:,25)+coef_xyz(78)*pol_z(:,1,kg)
                coef_xy(:,26)=coef_xy(:,26)+coef_xyz(79)*pol_z(:,1,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(80)*pol_z(:,1,kg)
                coef_xy(:,29)=coef_xy(:,29)+coef_xyz(81)*pol_z(:,1,kg)
                coef_xy(:,30)=coef_xy(:,30)+coef_xyz(82)*pol_z(:,1,kg)
                coef_xy(:,31)=coef_xy(:,31)+coef_xyz(83)*pol_z(:,1,kg)
                coef_xy(:,32)=coef_xy(:,32)+coef_xyz(84)*pol_z(:,1,kg)
                coef_xy(:,33)=coef_xy(:,33)+coef_xyz(85)*pol_z(:,1,kg)
                coef_xy(:,35)=coef_xy(:,35)+coef_xyz(86)*pol_z(:,1,kg)
                coef_xy(:,36)=coef_xy(:,36)+coef_xyz(87)*pol_z(:,1,kg)
                coef_xy(:,37)=coef_xy(:,37)+coef_xyz(88)*pol_z(:,1,kg)
                coef_xy(:,38)=coef_xy(:,38)+coef_xyz(89)*pol_z(:,1,kg)
                coef_xy(:,39)=coef_xy(:,39)+coef_xyz(90)*pol_z(:,1,kg)
                coef_xy(:,41)=coef_xy(:,41)+coef_xyz(91)*pol_z(:,1,kg)
                coef_xy(:,42)=coef_xy(:,42)+coef_xyz(92)*pol_z(:,1,kg)
                coef_xy(:,43)=coef_xy(:,43)+coef_xyz(93)*pol_z(:,1,kg)
                coef_xy(:,44)=coef_xy(:,44)+coef_xyz(94)*pol_z(:,1,kg)
                coef_xy(:,46)=coef_xy(:,46)+coef_xyz(95)*pol_z(:,1,kg)
                coef_xy(:,47)=coef_xy(:,47)+coef_xyz(96)*pol_z(:,1,kg)
                coef_xy(:,48)=coef_xy(:,48)+coef_xyz(97)*pol_z(:,1,kg)
                coef_xy(:,50)=coef_xy(:,50)+coef_xyz(98)*pol_z(:,1,kg)
                coef_xy(:,51)=coef_xy(:,51)+coef_xyz(99)*pol_z(:,1,kg)
                coef_xy(:,53)=coef_xy(:,53)+coef_xyz(100)*pol_z(:,1,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(101)*pol_z(:,2,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(102)*pol_z(:,2,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(103)*pol_z(:,2,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(104)*pol_z(:,2,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(105)*pol_z(:,2,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(106)*pol_z(:,2,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(107)*pol_z(:,2,kg)
                coef_xy(:,8)=coef_xy(:,8)+coef_xyz(108)*pol_z(:,2,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(109)*pol_z(:,2,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(110)*pol_z(:,2,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(111)*pol_z(:,2,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(112)*pol_z(:,2,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(113)*pol_z(:,2,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(114)*pol_z(:,2,kg)
                coef_xy(:,17)=coef_xy(:,17)+coef_xyz(115)*pol_z(:,2,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(116)*pol_z(:,2,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(117)*pol_z(:,2,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(118)*pol_z(:,2,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(119)*pol_z(:,2,kg)
                coef_xy(:,24)=coef_xy(:,24)+coef_xyz(120)*pol_z(:,2,kg)
                coef_xy(:,25)=coef_xy(:,25)+coef_xyz(121)*pol_z(:,2,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(122)*pol_z(:,2,kg)
                coef_xy(:,29)=coef_xy(:,29)+coef_xyz(123)*pol_z(:,2,kg)
                coef_xy(:,30)=coef_xy(:,30)+coef_xyz(124)*pol_z(:,2,kg)
                coef_xy(:,31)=coef_xy(:,31)+coef_xyz(125)*pol_z(:,2,kg)
                coef_xy(:,32)=coef_xy(:,32)+coef_xyz(126)*pol_z(:,2,kg)
                coef_xy(:,35)=coef_xy(:,35)+coef_xyz(127)*pol_z(:,2,kg)
                coef_xy(:,36)=coef_xy(:,36)+coef_xyz(128)*pol_z(:,2,kg)
                coef_xy(:,37)=coef_xy(:,37)+coef_xyz(129)*pol_z(:,2,kg)
                coef_xy(:,38)=coef_xy(:,38)+coef_xyz(130)*pol_z(:,2,kg)
                coef_xy(:,41)=coef_xy(:,41)+coef_xyz(131)*pol_z(:,2,kg)
                coef_xy(:,42)=coef_xy(:,42)+coef_xyz(132)*pol_z(:,2,kg)
                coef_xy(:,43)=coef_xy(:,43)+coef_xyz(133)*pol_z(:,2,kg)
                coef_xy(:,46)=coef_xy(:,46)+coef_xyz(134)*pol_z(:,2,kg)
                coef_xy(:,47)=coef_xy(:,47)+coef_xyz(135)*pol_z(:,2,kg)
                coef_xy(:,50)=coef_xy(:,50)+coef_xyz(136)*pol_z(:,2,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(137)*pol_z(:,3,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(138)*pol_z(:,3,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(139)*pol_z(:,3,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(140)*pol_z(:,3,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(141)*pol_z(:,3,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(142)*pol_z(:,3,kg)
                coef_xy(:,7)=coef_xy(:,7)+coef_xyz(143)*pol_z(:,3,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(144)*pol_z(:,3,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(145)*pol_z(:,3,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(146)*pol_z(:,3,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(147)*pol_z(:,3,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(148)*pol_z(:,3,kg)
                coef_xy(:,16)=coef_xy(:,16)+coef_xyz(149)*pol_z(:,3,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(150)*pol_z(:,3,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(151)*pol_z(:,3,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(152)*pol_z(:,3,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(153)*pol_z(:,3,kg)
                coef_xy(:,24)=coef_xy(:,24)+coef_xyz(154)*pol_z(:,3,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(155)*pol_z(:,3,kg)
                coef_xy(:,29)=coef_xy(:,29)+coef_xyz(156)*pol_z(:,3,kg)
                coef_xy(:,30)=coef_xy(:,30)+coef_xyz(157)*pol_z(:,3,kg)
                coef_xy(:,31)=coef_xy(:,31)+coef_xyz(158)*pol_z(:,3,kg)
                coef_xy(:,35)=coef_xy(:,35)+coef_xyz(159)*pol_z(:,3,kg)
                coef_xy(:,36)=coef_xy(:,36)+coef_xyz(160)*pol_z(:,3,kg)
                coef_xy(:,37)=coef_xy(:,37)+coef_xyz(161)*pol_z(:,3,kg)
                coef_xy(:,41)=coef_xy(:,41)+coef_xyz(162)*pol_z(:,3,kg)
                coef_xy(:,42)=coef_xy(:,42)+coef_xyz(163)*pol_z(:,3,kg)
                coef_xy(:,46)=coef_xy(:,46)+coef_xyz(164)*pol_z(:,3,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(165)*pol_z(:,4,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(166)*pol_z(:,4,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(167)*pol_z(:,4,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(168)*pol_z(:,4,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(169)*pol_z(:,4,kg)
                coef_xy(:,6)=coef_xy(:,6)+coef_xyz(170)*pol_z(:,4,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(171)*pol_z(:,4,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(172)*pol_z(:,4,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(173)*pol_z(:,4,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(174)*pol_z(:,4,kg)
                coef_xy(:,15)=coef_xy(:,15)+coef_xyz(175)*pol_z(:,4,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(176)*pol_z(:,4,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(177)*pol_z(:,4,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(178)*pol_z(:,4,kg)
                coef_xy(:,23)=coef_xy(:,23)+coef_xyz(179)*pol_z(:,4,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(180)*pol_z(:,4,kg)
                coef_xy(:,29)=coef_xy(:,29)+coef_xyz(181)*pol_z(:,4,kg)
                coef_xy(:,30)=coef_xy(:,30)+coef_xyz(182)*pol_z(:,4,kg)
                coef_xy(:,35)=coef_xy(:,35)+coef_xyz(183)*pol_z(:,4,kg)
                coef_xy(:,36)=coef_xy(:,36)+coef_xyz(184)*pol_z(:,4,kg)
                coef_xy(:,41)=coef_xy(:,41)+coef_xyz(185)*pol_z(:,4,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(186)*pol_z(:,5,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(187)*pol_z(:,5,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(188)*pol_z(:,5,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(189)*pol_z(:,5,kg)
                coef_xy(:,5)=coef_xy(:,5)+coef_xyz(190)*pol_z(:,5,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(191)*pol_z(:,5,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(192)*pol_z(:,5,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(193)*pol_z(:,5,kg)
                coef_xy(:,14)=coef_xy(:,14)+coef_xyz(194)*pol_z(:,5,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(195)*pol_z(:,5,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(196)*pol_z(:,5,kg)
                coef_xy(:,22)=coef_xy(:,22)+coef_xyz(197)*pol_z(:,5,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(198)*pol_z(:,5,kg)
                coef_xy(:,29)=coef_xy(:,29)+coef_xyz(199)*pol_z(:,5,kg)
                coef_xy(:,35)=coef_xy(:,35)+coef_xyz(200)*pol_z(:,5,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(201)*pol_z(:,6,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(202)*pol_z(:,6,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(203)*pol_z(:,6,kg)
                coef_xy(:,4)=coef_xy(:,4)+coef_xyz(204)*pol_z(:,6,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(205)*pol_z(:,6,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(206)*pol_z(:,6,kg)
                coef_xy(:,13)=coef_xy(:,13)+coef_xyz(207)*pol_z(:,6,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(208)*pol_z(:,6,kg)
                coef_xy(:,21)=coef_xy(:,21)+coef_xyz(209)*pol_z(:,6,kg)
                coef_xy(:,28)=coef_xy(:,28)+coef_xyz(210)*pol_z(:,6,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(211)*pol_z(:,7,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(212)*pol_z(:,7,kg)
                coef_xy(:,3)=coef_xy(:,3)+coef_xyz(213)*pol_z(:,7,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(214)*pol_z(:,7,kg)
                coef_xy(:,12)=coef_xy(:,12)+coef_xyz(215)*pol_z(:,7,kg)
                coef_xy(:,20)=coef_xy(:,20)+coef_xyz(216)*pol_z(:,7,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(217)*pol_z(:,8,kg)
                coef_xy(:,2)=coef_xy(:,2)+coef_xyz(218)*pol_z(:,8,kg)
                coef_xy(:,11)=coef_xy(:,11)+coef_xyz(219)*pol_z(:,8,kg)
                coef_xy(:,1)=coef_xy(:,1)+coef_xyz(220)*pol_z(:,9,kg)
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
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,1)*pol_y(1,0,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,1)*pol_y(2,0,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,2)*pol_y(1,0,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,2)*pol_y(2,0,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,3)*pol_y(1,0,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,3)*pol_y(2,0,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,4)*pol_y(1,0,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,4)*pol_y(2,0,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,5)*pol_y(1,0,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,5)*pol_y(2,0,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,6)*pol_y(1,0,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,6)*pol_y(2,0,jg)
             coef_x(1:2,6)=coef_x(1:2,6)+coef_xy(1:2,7)*pol_y(1,0,jg)
             coef_x(3:4,6)=coef_x(3:4,6)+coef_xy(1:2,7)*pol_y(2,0,jg)
             coef_x(1:2,7)=coef_x(1:2,7)+coef_xy(1:2,8)*pol_y(1,0,jg)
             coef_x(3:4,7)=coef_x(3:4,7)+coef_xy(1:2,8)*pol_y(2,0,jg)
             coef_x(1:2,8)=coef_x(1:2,8)+coef_xy(1:2,9)*pol_y(1,0,jg)
             coef_x(3:4,8)=coef_x(3:4,8)+coef_xy(1:2,9)*pol_y(2,0,jg)
             coef_x(1:2,9)=coef_x(1:2,9)+coef_xy(1:2,10)*pol_y(1,0,jg)
             coef_x(3:4,9)=coef_x(3:4,9)+coef_xy(1:2,10)*pol_y(2,0,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,11)*pol_y(1,1,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,11)*pol_y(2,1,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,12)*pol_y(1,1,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,12)*pol_y(2,1,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,13)*pol_y(1,1,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,13)*pol_y(2,1,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,14)*pol_y(1,1,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,14)*pol_y(2,1,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,15)*pol_y(1,1,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,15)*pol_y(2,1,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,16)*pol_y(1,1,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,16)*pol_y(2,1,jg)
             coef_x(1:2,6)=coef_x(1:2,6)+coef_xy(1:2,17)*pol_y(1,1,jg)
             coef_x(3:4,6)=coef_x(3:4,6)+coef_xy(1:2,17)*pol_y(2,1,jg)
             coef_x(1:2,7)=coef_x(1:2,7)+coef_xy(1:2,18)*pol_y(1,1,jg)
             coef_x(3:4,7)=coef_x(3:4,7)+coef_xy(1:2,18)*pol_y(2,1,jg)
             coef_x(1:2,8)=coef_x(1:2,8)+coef_xy(1:2,19)*pol_y(1,1,jg)
             coef_x(3:4,8)=coef_x(3:4,8)+coef_xy(1:2,19)*pol_y(2,1,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,20)*pol_y(1,2,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,20)*pol_y(2,2,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,21)*pol_y(1,2,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,21)*pol_y(2,2,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,22)*pol_y(1,2,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,22)*pol_y(2,2,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,23)*pol_y(1,2,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,23)*pol_y(2,2,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,24)*pol_y(1,2,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,24)*pol_y(2,2,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,25)*pol_y(1,2,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,25)*pol_y(2,2,jg)
             coef_x(1:2,6)=coef_x(1:2,6)+coef_xy(1:2,26)*pol_y(1,2,jg)
             coef_x(3:4,6)=coef_x(3:4,6)+coef_xy(1:2,26)*pol_y(2,2,jg)
             coef_x(1:2,7)=coef_x(1:2,7)+coef_xy(1:2,27)*pol_y(1,2,jg)
             coef_x(3:4,7)=coef_x(3:4,7)+coef_xy(1:2,27)*pol_y(2,2,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,28)*pol_y(1,3,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,28)*pol_y(2,3,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,29)*pol_y(1,3,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,29)*pol_y(2,3,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,30)*pol_y(1,3,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,30)*pol_y(2,3,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,31)*pol_y(1,3,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,31)*pol_y(2,3,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,32)*pol_y(1,3,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,32)*pol_y(2,3,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,33)*pol_y(1,3,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,33)*pol_y(2,3,jg)
             coef_x(1:2,6)=coef_x(1:2,6)+coef_xy(1:2,34)*pol_y(1,3,jg)
             coef_x(3:4,6)=coef_x(3:4,6)+coef_xy(1:2,34)*pol_y(2,3,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,35)*pol_y(1,4,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,35)*pol_y(2,4,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,36)*pol_y(1,4,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,36)*pol_y(2,4,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,37)*pol_y(1,4,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,37)*pol_y(2,4,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,38)*pol_y(1,4,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,38)*pol_y(2,4,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,39)*pol_y(1,4,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,39)*pol_y(2,4,jg)
             coef_x(1:2,5)=coef_x(1:2,5)+coef_xy(1:2,40)*pol_y(1,4,jg)
             coef_x(3:4,5)=coef_x(3:4,5)+coef_xy(1:2,40)*pol_y(2,4,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,41)*pol_y(1,5,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,41)*pol_y(2,5,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,42)*pol_y(1,5,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,42)*pol_y(2,5,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,43)*pol_y(1,5,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,43)*pol_y(2,5,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,44)*pol_y(1,5,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,44)*pol_y(2,5,jg)
             coef_x(1:2,4)=coef_x(1:2,4)+coef_xy(1:2,45)*pol_y(1,5,jg)
             coef_x(3:4,4)=coef_x(3:4,4)+coef_xy(1:2,45)*pol_y(2,5,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,46)*pol_y(1,6,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,46)*pol_y(2,6,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,47)*pol_y(1,6,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,47)*pol_y(2,6,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,48)*pol_y(1,6,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,48)*pol_y(2,6,jg)
             coef_x(1:2,3)=coef_x(1:2,3)+coef_xy(1:2,49)*pol_y(1,6,jg)
             coef_x(3:4,3)=coef_x(3:4,3)+coef_xy(1:2,49)*pol_y(2,6,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,50)*pol_y(1,7,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,50)*pol_y(2,7,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,51)*pol_y(1,7,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,51)*pol_y(2,7,jg)
             coef_x(1:2,2)=coef_x(1:2,2)+coef_xy(1:2,52)*pol_y(1,7,jg)
             coef_x(3:4,2)=coef_x(3:4,2)+coef_xy(1:2,52)*pol_y(2,7,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,53)*pol_y(1,8,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,53)*pol_y(2,8,jg)
             coef_x(1:2,1)=coef_x(1:2,1)+coef_xy(1:2,54)*pol_y(1,8,jg)
             coef_x(3:4,1)=coef_x(3:4,1)+coef_xy(1:2,54)*pol_y(2,8,jg)
             coef_x(1:2,0)=coef_x(1:2,0)+coef_xy(1:2,55)*pol_y(1,9,jg)
             coef_x(3:4,0)=coef_x(3:4,0)+coef_xy(1:2,55)*pol_y(2,9,jg)
          DO ig=igmin,igmax
             i=map(ig,1)
             s01=0.0_wp
             s02=0.0_wp
             s03=0.0_wp
             s04=0.0_wp
                s01=s01+coef_x(1,0)*pol_x(0,ig)
                s02=s02+coef_x(2,0)*pol_x(0,ig)
                s03=s03+coef_x(3,0)*pol_x(0,ig)
                s04=s04+coef_x(4,0)*pol_x(0,ig)
                s01=s01+coef_x(1,1)*pol_x(1,ig)
                s02=s02+coef_x(2,1)*pol_x(1,ig)
                s03=s03+coef_x(3,1)*pol_x(1,ig)
                s04=s04+coef_x(4,1)*pol_x(1,ig)
                s01=s01+coef_x(1,2)*pol_x(2,ig)
                s02=s02+coef_x(2,2)*pol_x(2,ig)
                s03=s03+coef_x(3,2)*pol_x(2,ig)
                s04=s04+coef_x(4,2)*pol_x(2,ig)
                s01=s01+coef_x(1,3)*pol_x(3,ig)
                s02=s02+coef_x(2,3)*pol_x(3,ig)
                s03=s03+coef_x(3,3)*pol_x(3,ig)
                s04=s04+coef_x(4,3)*pol_x(3,ig)
                s01=s01+coef_x(1,4)*pol_x(4,ig)
                s02=s02+coef_x(2,4)*pol_x(4,ig)
                s03=s03+coef_x(3,4)*pol_x(4,ig)
                s04=s04+coef_x(4,4)*pol_x(4,ig)
                s01=s01+coef_x(1,5)*pol_x(5,ig)
                s02=s02+coef_x(2,5)*pol_x(5,ig)
                s03=s03+coef_x(3,5)*pol_x(5,ig)
                s04=s04+coef_x(4,5)*pol_x(5,ig)
                s01=s01+coef_x(1,6)*pol_x(6,ig)
                s02=s02+coef_x(2,6)*pol_x(6,ig)
                s03=s03+coef_x(3,6)*pol_x(6,ig)
                s04=s04+coef_x(4,6)*pol_x(6,ig)
                s01=s01+coef_x(1,7)*pol_x(7,ig)
                s02=s02+coef_x(2,7)*pol_x(7,ig)
                s03=s03+coef_x(3,7)*pol_x(7,ig)
                s04=s04+coef_x(4,7)*pol_x(7,ig)
                s01=s01+coef_x(1,8)*pol_x(8,ig)
                s02=s02+coef_x(2,8)*pol_x(8,ig)
                s03=s03+coef_x(3,8)*pol_x(8,ig)
                s04=s04+coef_x(4,8)*pol_x(8,ig)
                s01=s01+coef_x(1,9)*pol_x(9,ig)
                s02=s02+coef_x(2,9)*pol_x(9,ig)
                s03=s03+coef_x(3,9)*pol_x(9,ig)
                s04=s04+coef_x(4,9)*pol_x(9,ig)
             grid(i,j,k) = grid(i,j,k)     + s01
             grid(i,j2,k) = grid(i,j2,k)   + s03
             grid(i,j,k2) = grid(i,j,k2)   + s02
             grid(i,j2,k2) = grid(i,j2,k2) + s04
          END DO
       END DO
    END DO

  END SUBROUTINE collocate_core_9
