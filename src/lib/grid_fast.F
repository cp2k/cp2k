  SUBROUTINE collocate_core(pzyx,polx,   &
                poly,polz,grid,alpha,lx_max,lxy_max,lxyz_max,coef_max, &
                cmax,ly_max,lz_max,gridbounds,&
                map,sphere_bounds)
    IMPLICIT NONE
#if defined(__SGL)
  INTEGER, PARAMETER :: wp = KIND(0.0)
#else
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND ( 14, 200 )
#endif
    integer, INTENT(IN)  :: lx_max,lxy_max,lxyz_max,coef_max,cmax
    integer, INTENT(IN)  :: ly_max(lx_max),lz_max(lxy_max)
    integer, INTENT(IN)  :: map(-cmax:cmax,1:3)
    integer, INTENT(IN)  :: gridbounds(2,3)
    integer, INTENT(IN)  :: sphere_bounds(*)
    real(wp), INTENT(IN) :: pzyx(lxyz_max)
    real(wp), INTENT(IN) :: polz(lxyz_max,-cmax:cmax)
    real(wp), INTENT(IN) :: poly(lxy_max,-cmax:cmax)
    real(wp), INTENT(IN) :: polx(coef_max,-cmax:cmax)
    real(wp), INTENT(IN) :: alpha(coef_max,lx_max)
    real(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1),gridbounds(1,2):gridbounds(2,2),gridbounds(1,3):gridbounds(2,3))

    real(wp) ::  coef(coef_max,4)
    real(wp) ::  pyx(2,lxy_max)
    integer kg,kgmax,k,jgmax,jgmin,jg,j,j2,igmax,ig,i,kgmin,igmin,k2,i2,jg2,kg2
    integer lx,ly,lz,lxy,lxyz,icoef,sci
    real(wp) s1,s2,s3,s4,s5,s6,s7,s8,s9,sa,g1,g2,g3
    real(wp) s1k,s2k,s3k,s4k,s5k,s6k,s7k,s8k,s9k,sak,g1k,g2k
    real(wp) r1,r2,r3,r4,r5,r6,r7,r8
    integer csign


    sci=1
    kgmin=sphere_bounds(sci)
    sci=sci+1
    DO kg=kgmin,0
       kg2=1-kg
       k=map(kg,3)
       k2=map(kg2,3)
       lxyz=0
       DO lxy=1,lxy_max
          pyx(1,lxy)=0.0_wp
          pyx(2,lxy)=0.0_wp
          DO lz=1,lz_max(lxy)
             lxyz=lxyz+1
             pyx(1,lxy)=pyx(1,lxy)+pzyx(lxyz)*polz(lxyz,kg)
             pyx(2,lxy)=pyx(2,lxy)+pzyx(lxyz)*polz(lxyz,kg2)
          ENDDO
       ENDDO
       jgmin=sphere_bounds(sci)
       sci=sci+1
       DO jg=jgmin,0
          jg2=1-jg
          j=map(jg,2)
          j2=map(jg2,2)
          csign=1
          if (j.eq.j2 .or. k.eq.k2) csign=-1 ! skips unrolled routines if there can be a collision
          igmin=sphere_bounds(sci)
          sci=sci+1
          igmax=1-igmin
          SELECT CASE (csign*coef_max)
          CASE(1)
             s1=0.0_wp
             s2=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             lxy=0
             DO ly=1,ly_max(1)
                lxy=lxy+1
                s1=s1+pyx(1,lxy)*poly(lxy,jg)
                s2=s2+pyx(1,lxy)*poly(lxy,jg2)
                s1k=s1k+pyx(2,lxy)*poly(lxy,jg)
                s2k=s2k+pyx(2,lxy)*poly(lxy,jg2)
             ENDDO
             s1=alpha(1,1)*s1
             s2=alpha(1,1)*s2
             s1k=alpha(1,1)*s1k
             s2k=alpha(1,1)*s2k
             DO ig=igmin,igmax,2
                i=map(ig,1)
                i2=map(ig+1,1)
                g1=polx(1,ig)
                g2=polx(1,ig+1)
                r1=grid(i,j,k)
                r2=grid(i,j2,k)
                r3=grid(i,j,k2)
                r4=grid(i,j2,k2)
                r1= r1 + s1*g1
                r2= r2 + s2*g1
                r3= r3 + s1k*g1
                r4= r4 + s2k*g1
                grid(i,j,k)   = r1 
                grid(i,j2,k)  = r2
                grid(i,j,k2)  = r3
                grid(i,j2,k2) = r4 
                r5=grid(i2,j,k)
                r6=grid(i2,j2,k)
                r7=grid(i2,j,k2)
                r8=grid(i2,j2,k2)
                r5 = r5 + s1*g2
                r6 = r6 + s2*g2
                r7 = r7 + s1k*g2
                r8 = r8 + s2k*g2
                grid(i2,j,k)   = r5 
                grid(i2,j2,k)  = r6 
                grid(i2,j,k2)  = r7
                grid(i2,j2,k2) = r8 
             END DO
          CASE(2)
             lxy=0
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             DO lx=1,2
              g1=0.0_wp
              g2=0.0_wp
              g1k=0.0_wp
              g2k=0.0_wp
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 r1=pyx(1,lxy)
                 r2=pyx(2,lxy)
                 r3=poly(lxy,jg)
                 r4=poly(lxy,jg2)
                 g1=g1+r1*r3
                 g2=g2+r1*r4
                 g1k=g1k+r2*r3
                 g2k=g2k+r2*r4
              ENDDO
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              s1=s1+r1*g1
              s2=s2+r2*g1
              s3=s3+r1*g2
              s4=s4+r2*g2
              s1k=s1k+r1*g1k
              s2k=s2k+r2*g1k
              s3k=s3k+r1*g2k
              s4k=s4k+r2*g2k
             ENDDO
             DO ig=igmin,igmax,2
                i=map(ig,1)
                i2=map(ig+1,1)
                g1=polx(1,ig)
                g2=polx(1,ig+1)
                r1=grid(i,j,k)+ s1*g1
                r2=grid(i,j2,k)+ s3*g1
                r3=grid(i,j,k2)+ s1k*g1
                r4=grid(i,j2,k2)+ s3k*g1
                r5=grid(i2,j,k)+ s1*g2
                r6=grid(i2,j2,k)+ s3*g2
                r7=grid(i2,j,k2)+ s1k*g2
                r8=grid(i2,j2,k2)+ s3k*g2
                g1=polx(2,ig)
                g2=polx(2,ig+1)
                grid(i,j,k) = r1+ s2*g1
                grid(i,j2,k) = r2+ s4*g1
                grid(i,j,k2) = r3+ s2k*g1
                grid(i,j2,k2) = r4+ s4k*g1
                grid(i2,j,k) = r5+ s2*g2
                grid(i2,j2,k) = r6+ s4*g2
                grid(i2,j,k2) = r7+ s2k*g2
                grid(i2,j2,k2) = r8+ s4k*g2
             END DO
          CASE(3)
             lxy=0
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s5=0.0_wp
             s6=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             s5k=0.0_wp
             s6k=0.0_wp
             DO lx=1,lx_max
              g1=0.0_wp
              g2=0.0_wp
              g1k=0.0_wp
              g2k=0.0_wp
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 r1=pyx(1,lxy)
                 r2=pyx(2,lxy)
                 r3=poly(lxy,jg)
                 r4=poly(lxy,jg2)
                 g1=g1+r1*r3
                 g2=g2+r1*r4
                 g1k=g1k+r2*r3
                 g2k=g2k+r2*r4
              ENDDO
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              r3=alpha(3,lx)
              s1=s1+r1*g1
              s2=s2+r2*g1
              s3=s3+r3*g1
              s4=s4+r1*g2
              s5=s5+r2*g2
              s6=s6+r3*g2
              s1k=s1k+r1*g1k
              s2k=s2k+r2*g1k
              s3k=s3k+r3*g1k
              s4k=s4k+r1*g2k
              s5k=s5k+r2*g2k
              s6k=s6k+r3*g2k
             ENDDO
             DO ig=igmin,igmax,2
                i=map(ig,1)
                i2=map(ig+1,1)
                g1=polx(1,ig)
                g2=polx(1,ig+1)
                r1=grid(i,j,k)+ s1*g1
                r2=grid(i,j2,k)+ s4*g1
                r3=grid(i,j,k2)+ s1k*g1
                r4=grid(i,j2,k2)+ s4k*g1
                r5=grid(i2,j,k)+ s1*g2
                r6=grid(i2,j2,k)+ s4*g2
                r7=grid(i2,j,k2)+ s1k*g2
                r8=grid(i2,j2,k2)+ s4k*g2
                g1=polx(2,ig)
                g2=polx(2,ig+1)
                r1=r1+s2*g1
                r2=r2+s5*g1
                r3=r3+s2k*g1
                r4=r4+s5k*g1
                r5=r5+s2*g2
                r6=r6+s5*g2
                r7=r7+s2k*g2
                r8=r8+s5k*g2
                g1=polx(3,ig)
                g2=polx(3,ig+1)
                grid(i,j,k) = r1 + s3*g1
                grid(i,j2,k)= r2 + s6*g1
                grid(i,j,k2) = r3 + s3k*g1
                grid(i,j2,k2) = r4 + s6k*g1
                grid(i2,j,k) = r5 + s3*g2
                grid(i2,j2,k)= r6 + s6*g2
                grid(i2,j,k2) = r7 + s3k*g2
                grid(i2,j2,k2) = r8 + s6k*g2
             END DO
          CASE(4)
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s5=0.0_wp
             s6=0.0_wp
             s7=0.0_wp
             s8=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             s5k=0.0_wp
             s6k=0.0_wp
             s7k=0.0_wp
             s8k=0.0_wp
             lxy=0
             DO lx=1,lx_max
              g1=0.0_wp
              g2=0.0_wp
              g1k=0.0_wp
              g2k=0.0_wp
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 r1=pyx(1,lxy)
                 r2=pyx(2,lxy)
                 r3=poly(lxy,jg)
                 r4=poly(lxy,jg2)
                 g1=g1+r1*r3
                 g2=g2+r1*r4
                 g1k=g1k+r2*r3
                 g2k=g2k+r2*r4
              ENDDO
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              r3=alpha(3,lx)
              r4=alpha(4,lx)
              s1=s1+r1*g1
              s2=s2+r2*g1
              s3=s3+r3*g1
              s4=s4+r4*g1
              s5=s5+r1*g2
              s6=s6+r2*g2
              s7=s7+r3*g2
              s8=s8+r4*g2
              s1k=s1k+r1*g1k
              s2k=s2k+r2*g1k
              s3k=s3k+r3*g1k
              s4k=s4k+r4*g1k
              s5k=s5k+r1*g2k
              s6k=s6k+r2*g2k
              s7k=s7k+r3*g2k
              s8k=s8k+r4*g2k
             ENDDO
             DO ig=igmin,igmax,2
                i=map(ig,1)
                i2=map(ig+1,1)
                g1=polx(1,ig)
                g2=polx(1,ig+1)
                r1=grid(i,j,k)+ s1*g1
                r2=grid(i,j2,k)+ s5*g1
                r3=grid(i,j,k2)+ s1k*g1
                r4=grid(i,j2,k2)+ s5k*g1
                r5=grid(i2,j,k)+ s1*g2
                r6=grid(i2,j2,k)+ s5*g2
                r7=grid(i2,j,k2)+ s1k*g2
                r8=grid(i2,j2,k2)+ s5k*g2
                g1=polx(2,ig)
                g2=polx(2,ig+1)
                r1=r1+s2*g1
                r2=r2+s6*g1
                r3=r3+s2k*g1
                r4=r4+s6k*g1
                r5=r5+s2*g2
                r6=r6+s6*g2
                r7=r7+s2k*g2
                r8=r8+s6k*g2
                g1=polx(3,ig)
                g2=polx(3,ig+1)
                r1=r1+s3*g1
                r2=r2+s7*g1
                r3=r3+s3k*g1
                r4=r4+s7k*g1
                r5=r5+s3*g2
                r6=r6+s7*g2
                r7=r7+s3k*g2
                r8=r8+s7k*g2
                g1=polx(4,ig)
                g2=polx(4,ig+1)
                grid(i,j,k) = r1 + s4*g1
                grid(i,j2,k) = r2+ s8*g1
                grid(i,j,k2) = r3+ s4k*g1
                grid(i,j2,k2) = r4+ s8k*g1
                grid(i2,j,k) = r5 + s4*g2
                grid(i2,j2,k) = r6+ s8*g2
                grid(i2,j,k2) = r7+ s4k*g2
                grid(i2,j2,k2) = r8+ s8k*g2
             END DO
          CASE(5)
             lxy=0
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s9=0.0_wp
             s5=0.0_wp
             s6=0.0_wp
             s7=0.0_wp
             s8=0.0_wp
             sa=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             s9k=0.0_wp
             s5k=0.0_wp
             s6k=0.0_wp
             s7k=0.0_wp
             s8k=0.0_wp
             sak=0.0_wp
             DO lx=1,lx_max
              g1=0.0_wp
              g2=0.0_wp
              g1k=0.0_wp
              g2k=0.0_wp
              DO ly=1,ly_max(lx)
                 lxy=lxy+1 
                 r1=pyx(1,lxy)
                 r2=pyx(2,lxy)
                 r3=poly(lxy,jg)
                 r4=poly(lxy,jg2)
                 g1=g1+r1*r3
                 g2=g2+r1*r4
                 g1k=g1k+r2*r3
                 g2k=g2k+r2*r4
              ENDDO
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              r3=alpha(3,lx)
              r4=alpha(4,lx)
              r5=alpha(5,lx)
              s1=s1+r1*g1
              s2=s2+r2*g1
              s3=s3+r3*g1
              s4=s4+r4*g1
              s9=s9+r5*g1
              s5=s5+r1*g2
              s6=s6+r2*g2
              s7=s7+r3*g2
              s8=s8+r4*g2
              sa=sa+r5*g2
              s1k=s1k+r1*g1k
              s2k=s2k+r2*g1k
              s3k=s3k+r3*g1k
              s4k=s4k+r4*g1k
              s9k=s9k+r5*g1k
              s5k=s5k+r1*g2k
              s6k=s6k+r2*g2k
              s7k=s7k+r3*g2k
              s8k=s8k+r4*g2k
              sak=sak+r5*g2k
             ENDDO
             DO ig=igmin,igmax,2
                i=map(ig,1)
                i2=map(ig+1,1)
                g1=polx(1,ig)
                g2=polx(1,ig+1)
                r1= grid(i,j,k) + s1*g1
                r2=grid(i,j2,k) + s5*g1
                r3=grid(i,j,k2) + s1k*g1
                r4=grid(i,j2,k2)+ s5k*g1
                r5= grid(i2,j,k) + s1*g2
                r6=grid(i2,j2,k) + s5*g2
                r7=grid(i2,j,k2) + s1k*g2
                r8=grid(i2,j2,k2)+ s5k*g2
                g1=polx(2,ig)
                g2=polx(2,ig+1)
                r1=r1+s2*g1
                r2=r2+s6*g1
                r3=r3+s2k*g1
                r4=r4+s6k*g1
                r5=r5+s2*g2
                r6=r6+s6*g2
                r7=r7+s2k*g2
                r8=r8+s6k*g2
                g1=polx(3,ig)
                g2=polx(3,ig+1)
                r1=r1+s3*g1
                r2=r2+s7*g1
                r3=r3+s3k*g1
                r4=r4+s7k*g1
                r5=r5+s3*g2
                r6=r6+s7*g2
                r7=r7+s3k*g2
                r8=r8+s7k*g2
                g1=polx(4,ig)
                g2=polx(4,ig+1)
                r1=r1+s4*g1
                r2=r2+s8*g1
                r3=r3+s4k*g1
                r4=r4+s8k*g1
                r5=r5+s4*g2
                r6=r6+s8*g2
                r7=r7+s4k*g2
                r8=r8+s8k*g2
                g1=polx(5,ig)
                g2=polx(5,ig+1)
                grid(i,j,k) = r1 + s9*g1
                grid(i,j2,k)= r2 + sa*g1
                grid(i,j,k2)= r3 + s9k*g1
                grid(i,j2,k2)= r4+ sak*g1
                grid(i2,j,k) = r5 + s9*g2
                grid(i2,j2,k)= r6 + sa*g2
                grid(i2,j,k2)= r7 + s9k*g2
                grid(i2,j2,k2)= r8+ sak*g2
             END DO

          CASE DEFAULT
             coef(:,:)=0.0
             lxy=0
             DO lx=1,lx_max
              g1=0.0_wp
              g2=0.0_wp
              g1k=0.0_wp
              g2k=0.0_wp
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 g1=g1+pyx(1,lxy)*poly(lxy,jg)
                 g2=g2+pyx(1,lxy)*poly(lxy,jg2)
                 g1k=g1k+pyx(2,lxy)*poly(lxy,jg)
                 g2k=g2k+pyx(2,lxy)*poly(lxy,jg2)
              ENDDO
              DO icoef=1,coef_max
                 coef(icoef,1)=coef(icoef,1)+alpha(icoef,lx)*g1
                 coef(icoef,2)=coef(icoef,2)+alpha(icoef,lx)*g2
                 coef(icoef,3)=coef(icoef,3)+alpha(icoef,lx)*g1k
                 coef(icoef,4)=coef(icoef,4)+alpha(icoef,lx)*g2k
              ENDDO
             ENDDO
             DO ig=igmin,igmax
                i=map(ig,1)
                s1=0.0_wp
                s2=0.0_wp
                s1k=0.0_wp
                s2k=0.0_wp
                DO icoef=1,coef_max
                  s1=s1+coef(icoef,1)*polx(icoef,ig)
                  s2=s2+coef(icoef,2)*polx(icoef,ig)
                  s1k=s1k+coef(icoef,3)*polx(icoef,ig)
                  s2k=s2k+coef(icoef,4)*polx(icoef,ig)
                ENDDO
                grid(i,j,k) = grid(i,j,k)     + s1
                grid(i,j2,k) = grid(i,j2,k)   + s2
                grid(i,j,k2) = grid(i,j,k2)   + s1k
                grid(i,j2,k2) = grid(i,j2,k2) + s2k
             END DO
          END SELECT
       END DO
    END DO
  

  END SUBROUTINE

  SUBROUTINE integrate_core(polx,poly,polz,grid,alpha,lx_max,lxy_max,lxyz_max,&
                            coef_max,cmax,ly_max,lz_max,gridbounds,&
                            map,pzyx,la_min,la_max,&
                            lb_min,lb_max,sphere_bounds)
    IMPLICIT NONE
#if defined(__SGL)
  INTEGER, PARAMETER :: wp = KIND(0.0)
#else
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND ( 14, 200 )
#endif


    integer,  INTENT(IN) :: sphere_bounds(*),cmax
    integer,  INTENT(IN) :: lx_max,lxy_max,lxyz_max,coef_max,&
                            la_min,la_max,lb_min,lb_max
    integer,  INTENT(IN) :: ly_max(lx_max),lz_max(lxy_max)
    integer,  INTENT(IN) :: map(-cmax:cmax,1:3)
    integer,  INTENT(IN) :: gridbounds(2,3)
    real(wp), INTENT(IN) :: polz(lxyz_max,-cmax:cmax)
    real(wp), INTENT(IN) :: poly(lxy_max,-cmax:cmax)
    real(wp), INTENT(IN) :: polx(coef_max,-cmax:cmax)
    real(wp) pzyx(lxyz_max)
    real(wp), INTENT(IN) :: alpha(coef_max,lx_max)
    real(wp), INTENT(IN) :: grid(gridbounds(1,1):gridbounds(2,1),&
                            gridbounds(1,2):gridbounds(2,2),&
                            gridbounds(1,3):gridbounds(2,3))
    real(wp) coef(coef_max,4)
    real(wp) pyx(4,lxy_max)

    real(wp) radius2,dz2,dy2
    integer kg,kgmax,k,jgmax,jgmin,jg,j,j2,igmax,ig,i,kgmin,igmin,k2,kg2,jg2
    integer lx,ly,lz,lxy,lxyz,icoef,ico,jco,i2,sci
    real(wp) s1,s2,s3,s4,s5,s6,s7,s8,s9,sa,g1,g2
    real(wp) s1k,s2k,s3k,s4k,s5k,s6k,s7k,s8k,s9k,sak,g1k,g2k
    real(wp) r1,r2,r3,r4,r5,r6

    pzyx(:)=0.0_wp
    sci=1
    kgmin=sphere_bounds(sci)
    sci=sci+1

    DO kg=kgmin,0
       kg2=1-kg
       k=map(kg,3)
       k2=map(kg2,3)
       pyx(:,:)=0.0_wp
       jgmin=sphere_bounds(sci)
       sci=sci+1
       DO jg=jgmin,0
          jg2=1-jg
          j=map(jg,2)
          j2=map(jg2,2)
          igmin=sphere_bounds(sci)
          sci=sci+1
          igmax=1-igmin
          SELECT CASE (coef_max)
          CASE (1)
             s1=0.0_wp
             s2=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             DO ig=igmin,igmax,2
                i =map(ig,1)
                i2 =map(ig+1,1)
                r1=polx(1,ig)
                r2=polx(1,ig+1)
                s1=s1+grid(i,j,k)     *r1
                s2=s2+grid(i,j2,k)    *r1
                s1k=s1k+grid(i,j,k2)  *r1
                s2k=s2k+grid(i,j2,k2) *r1
                s3=s3+grid(i2,j,k)     *r2
                s4=s4+grid(i2,j2,k)    *r2
                s3k=s3k+grid(i2,j,k2)  *r2
                s4k=s4k+grid(i2,j2,k2) *r2
             END DO
             r3=alpha(1,1)
             r1=poly(1,jg)*r3
             r2=poly(1,jg2)*r3
             pyx(1,1)=pyx(1,1)+(s1+s3)*r1
             pyx(2,1)=pyx(2,1)+(s2+s4)*r2
             pyx(3,1)=pyx(3,1)+(s1k+s3k)*r1
             pyx(4,1)=pyx(4,1)+(s2k+s4k)*r2
          CASE (2)
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             DO ig=igmin,igmax
                i=map(ig,1)
                r1=grid(i,j,k)
                r2=grid(i,j2,k)
                r3=grid(i,j,k2)
                r4=grid(i,j2,k2)
                r5=polx(1,ig)
                r6=polx(2,ig)
                s1=s1+r1*r5
                s3=s3+r1*r6
                s2=s2+r2*r5
                s4=s4+r2*r6
                s1k=s1k+r3*r5
                s3k=s3k+r3*r6
                s2k=s2k+r4*r5
                s4k=s4k+r4*r6
             END DO
             lxy=0
             DO lx=1,lx_max
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              g1=r1*s1+r2*s3
              g2=r1*s2+r2*s4
              g1k=r1*s1k+r2*s3k
              g2k=r1*s2k+r2*s4k
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 r1=poly(lxy,jg)
                 r2=poly(lxy,jg2)
                 pyx(1,lxy)=pyx(1,lxy)+g1*r1
                 pyx(2,lxy)=pyx(2,lxy)+g2*r2
                 pyx(3,lxy)=pyx(3,lxy)+g1k*r1
                 pyx(4,lxy)=pyx(4,lxy)+g2k*r2
              ENDDO
             ENDDO
          CASE (3)
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s5=0.0_wp
             s6=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             s5k=0.0_wp
             s6k=0.0_wp
             DO ig=igmin,igmax
                i=map(ig,1)
                r1=grid(i,j,k)
                r2=grid(i,j2,k)
                r5=polx(1,ig)
                s1=s1+r1*r5
                s2=s2+r2*r5
                r3=grid(i,j,k2)
                r4=grid(i,j2,k2)
                s1k=s1k+r3*r5
                s2k=s2k+r4*r5
                r5=polx(2,ig)
                s3=s3+r1*r5
                s4=s4+r2*r5
                s3k=s3k+r3*r5
                s4k=s4k+r4*r5
                r5=polx(3,ig)
                s5=s5+r1*r5
                s6=s6+r2*r5
                s5k=s5k+r3*r5
                s6k=s6k+r4*r5
             END DO
             lxy=0
             DO lx=1,lx_max
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              r3=alpha(3,lx)
              g1=r1*s1+r2*s3+r3*s5
              g2=r1*s2+r2*s4+r3*s6
              g1k=r1*s1k+r2*s3k+r3*s5k
              g2k=r1*s2k+r2*s4k+r3*s6k
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 r1=poly(lxy,jg)
                 r2=poly(lxy,jg2)
                 pyx(1,lxy)=pyx(1,lxy)+g1*r1
                 pyx(2,lxy)=pyx(2,lxy)+g2*r2
                 pyx(3,lxy)=pyx(3,lxy)+g1k*r1
                 pyx(4,lxy)=pyx(4,lxy)+g2k*r2
              ENDDO
             ENDDO
          CASE (4)
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s5=0.0_wp
             s6=0.0_wp
             s7=0.0_wp
             s8=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             s5k=0.0_wp
             s6k=0.0_wp
             s7k=0.0_wp
             s8k=0.0_wp
             DO ig=igmin,igmax
                i=map(ig,1)
                s1=s1+grid(i,j ,k)*polx(1,ig)
                s2=s2+grid(i,j2,k)*polx(1,ig)
                s3=s3+grid(i,j ,k)*polx(2,ig)
                s4=s4+grid(i,j2,k)*polx(2,ig)
                s5=s5+grid(i,j ,k)*polx(3,ig)
                s6=s6+grid(i,j2,k)*polx(3,ig)
                s7=s7+grid(i,j ,k)*polx(4,ig)
                s8=s8+grid(i,j2,k)*polx(4,ig)
                s1k=s1k+grid(i,j ,k2)*polx(1,ig)
                s2k=s2k+grid(i,j2,k2)*polx(1,ig)
                s3k=s3k+grid(i,j ,k2)*polx(2,ig)
                s4k=s4k+grid(i,j2,k2)*polx(2,ig)
                s5k=s5k+grid(i,j ,k2)*polx(3,ig)
                s6k=s6k+grid(i,j2,k2)*polx(3,ig)
                s7k=s7k+grid(i,j ,k2)*polx(4,ig)
                s8k=s8k+grid(i,j2,k2)*polx(4,ig)
             END DO
             lxy=0
             DO lx=1,lx_max
              g1=alpha(1,lx)*s1+alpha(2,lx)*s3+alpha(3,lx)*s5+alpha(4,lx)*s7
              g2=alpha(1,lx)*s2+alpha(2,lx)*s4+alpha(3,lx)*s6+alpha(4,lx)*s8
              g1k=alpha(1,lx)*s1k+alpha(2,lx)*s3k+alpha(3,lx)*s5k+alpha(4,lx)*s7k
              g2k=alpha(1,lx)*s2k+alpha(2,lx)*s4k+alpha(3,lx)*s6k+alpha(4,lx)*s8k
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 pyx(1,lxy)=pyx(1,lxy)+g1*poly(lxy,jg)
                 pyx(2,lxy)=pyx(2,lxy)+g2*poly(lxy,jg2)
                 pyx(3,lxy)=pyx(3,lxy)+g1k*poly(lxy,jg)
                 pyx(4,lxy)=pyx(4,lxy)+g2k*poly(lxy,jg2)
              ENDDO
             ENDDO
           CASE (5)
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s5=0.0_wp
             s6=0.0_wp
             s7=0.0_wp
             s8=0.0_wp
             s9=0.0_wp
             sa=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             s5k=0.0_wp
             s6k=0.0_wp
             s7k=0.0_wp
             s8k=0.0_wp
             s9k=0.0_wp
             sak=0.0_wp
             DO ig=igmin,igmax
                i=map(ig,1)
                g1=grid(i,j,k)
                g2=grid(i,j2,k)
                g1k=grid(i,j,k2)
                g2k=grid(i,j2,k2)
                r1=polx(1,ig)
                r2=polx(2,ig)
                r3=polx(3,ig)
                r4=polx(4,ig)
                r5=polx(5,ig)
                s1=s1+g1*r1
                s2=s2+g2*r1
                s3=s3+g1*r2
                s4=s4+g2*r2
                s5=s5+g1*r3
                s6=s6+g2*r3
                s7=s7+g1*r4
                s8=s8+g2*r4
                s9=s9+g1*r5
                sa=sa+g2*r5
                s1k=s1k+g1k*r1
                s2k=s2k+g2k*r1
                s3k=s3k+g1k*r2
                s4k=s4k+g2k*r2
                s5k=s5k+g1k*r3
                s6k=s6k+g2k*r3
                s7k=s7k+g1k*r4
                s8k=s8k+g2k*r4
                s9k=s9k+g1k*r5
                sak=sak+g2k*r5
             END DO
             lxy=0
             DO lx=1,lx_max
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              r3=alpha(3,lx)
              r4=alpha(4,lx)
              r5=alpha(5,lx)
              g1=r1*s1
              g2=r1*s2
              g1k=r1*s1k
              g2k=r1*s2k
              g1=g1+r2*s3
              g2=g2+r2*s4
              g1k=g1k+r2*s3k
              g2k=g2k+r2*s4k
              g1=g1+r3*s5
              g2=g2+r3*s6
              g1k=g1k+r3*s5k
              g2k=g2k+r3*s6k
              g1=g1+r4*s7
              g2=g2+r4*s8
              g1k=g1k+r4*s7k
              g2k=g2k+r4*s8k
              g1=g1+r5*s9
              g2=g2+r5*sa
              g1k=g1k+r5*s9k
              g2k=g2k+r5*sak
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 r1=poly(lxy,jg)
                 r2=poly(lxy,jg2)
                 pyx(1,lxy)=pyx(1,lxy)+g1*r1
                 pyx(2,lxy)=pyx(2,lxy)+g2*r2
                 pyx(3,lxy)=pyx(3,lxy)+g1k*r1
                 pyx(4,lxy)=pyx(4,lxy)+g2k*r2
              ENDDO
             ENDDO
          CASE DEFAULT
             coef(:,:)=0.0_wp
             DO ig=igmin,igmax
                i=map(ig,1)
                DO icoef=1,coef_max
                  coef(icoef,1)=coef(icoef,1)+grid(i,j,k)*polx(icoef,ig)
                  coef(icoef,2)=coef(icoef,2)+grid(i,j2,k)*polx(icoef,ig)
                  coef(icoef,3)=coef(icoef,3)+grid(i,j,k2)*polx(icoef,ig)
                  coef(icoef,4)=coef(icoef,4)+grid(i,j2,k2)*polx(icoef,ig)
                ENDDO
             END DO
             lxy=0
             DO lx=1,lx_max
              g1=0.0_wp
              g2=0.0_wp
              g1k=0.0_wp
              g2k=0.0_wp
              DO icoef=1,coef_max
                  g1=g1+alpha(icoef,lx)*coef(icoef,1)
                  g2=g2+alpha(icoef,lx)*coef(icoef,2)
                  g1k=g1k+alpha(icoef,lx)*coef(icoef,3)
                  g2k=g2k+alpha(icoef,lx)*coef(icoef,4)
              ENDDO
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 pyx(1,lxy)=pyx(1,lxy)+g1*poly(lxy,jg)
                 pyx(2,lxy)=pyx(2,lxy)+g2*poly(lxy,jg2)
                 pyx(3,lxy)=pyx(3,lxy)+g1k*poly(lxy,jg)
                 pyx(4,lxy)=pyx(4,lxy)+g2k*poly(lxy,jg2)
              ENDDO
             ENDDO
          END SELECT
       END DO ! loop over y, jg
       lxyz=0
       DO lxy=1,lxy_max
          s1=pyx(1,lxy)+pyx(2,lxy)
          s1k=pyx(3,lxy)+pyx(4,lxy)
          DO lz=1,lz_max(lxy)
             lxyz=lxyz+1
             pzyx(lxyz)=pzyx(lxyz)+s1*polz(lxyz,kg)+s1k*polz(lxyz,kg2)
          ENDDO
       ENDDO
    END DO ! loop over z, kg


  END SUBROUTINE


  SUBROUTINE collocate_gf_npbc(grid,xdat,ydat,zdat,bo1,bo2,zlb,zub,ylb,yub,xlb,xub)
    IMPLICIT NONE
    ! Arguments
#if defined(__SGL)
  INTEGER, PARAMETER :: wp = KIND(0.0)
#else
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND ( 14, 200 )
#endif
    INTEGER, INTENT(IN)     :: bo1(2,3), bo2(2,3)
    INTEGER, INTENT(IN)     :: zlb,zub,ylb,yub,xlb,xub
    REAL(wp), INTENT(INOUT) :: grid(bo1(1,1):bo1(2,1),&
                                    bo1(1,2):bo1(2,2),&
                                    bo1(1,3):bo1(2,3))
    REAL(wp), INTENT(IN)    :: xdat(bo2(1,1):bo2(2,1)),&
                               ydat(bo2(1,2):bo2(2,2)),&
                               zdat(bo2(1,3):bo2(2,3))
    ! Local Variables
    INTEGER :: iz, iy, ix
    REAL(wp):: tmp1

    DO iz=zlb,zub
       DO iy=ylb,yub
          tmp1 = zdat(iz) * ydat(iy)
          DO ix=xlb,xub
             grid(ix,iy,iz) = grid(ix,iy,iz) + xdat(ix) * tmp1
          END DO ! Loop on x
       END DO ! Loop on y
    END DO ! Loop on z

  END SUBROUTINE

  SUBROUTINE integrate_gf_npbc(grid,xdat,ydat,zdat,bo1,bo2,zlb,zub,ylb,yub,xlb,xub,force)
    IMPLICIT NONE
    ! Arguments
#if defined(__SGL)
  INTEGER, PARAMETER :: wp = KIND(0.0)
#else
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND ( 14, 200 )
#endif
    INTEGER, INTENT(IN)     :: bo1(2,3), bo2(2,3)
    INTEGER, INTENT(IN)     :: zlb,zub,ylb,yub,xlb,xub
    REAL(wp), INTENT(INOUT) :: grid(bo1(1,1):bo1(2,1),&
                                    bo1(1,2):bo1(2,2),&
                                    bo1(1,3):bo1(2,3))
    REAL(wp), INTENT(IN)    :: xdat(bo2(1,1):bo2(2,1)),&
                               ydat(bo2(1,2):bo2(2,2)),&
                               zdat(bo2(1,3):bo2(2,3))
    REAL(wp), INTENT(INOUT) :: force
    ! Local Variables
    INTEGER :: iz, iy, ix
    REAL(wp):: tmp1

    DO iz=zlb,zub
       DO iy=ylb,yub
          tmp1 = zdat(iz) * ydat(iy)
          DO ix=xlb,xub
             force = force + grid(ix,iy,iz) * xdat(ix) * tmp1
          END DO ! Loop on x
       END DO ! Loop on y
    END DO ! Loop on z

  END SUBROUTINE
