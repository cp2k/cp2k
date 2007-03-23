
!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
! 
!   into 
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
    lp=la_max_local+lb_max_local
    ALLOCATE(coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)) 
    ALLOCATE(alpha(0:lp,0:la_max_local,0:lb_max_local,3))

    ALLOCATE(pol_z(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_y(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_x(0:lp,-cmax:cmax))

!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
    alpha(:,:,:,:)=0.0_dp
    DO iaxis=1,3
    DO lxa=0,la_max_local
    DO lxb=0,lb_max_local
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,iaxis)=alpha(lxa-l+lxb-k,lxa,lxb,iaxis)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(iaxis)-(ra(iaxis)+rab(iaxis)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(iaxis)+rp(iaxis))
       ENDDO
    ENDDO
    ENDDO
    ENDDO

!
!   compute the values of all (x-xp)**lp*exp(..)
!
! still requires the old trick:
!  new trick to avoid to many exps (reuse the result from the previous gridpoint):
!  exp( -a*(x+d)**2)=exp(-a*x**2)*(-2*a*x*d)*exp(-a*d**2)
!  exp(-2*a*(x+d)*d)=exp(-2*a*x*d)*exp(-2*a*d**2)

      iaxis=3
      t_exp_1=EXP(-zetp*dr(iaxis)**2)
      t_exp_2=t_exp_1**2
      t_exp_min_1=EXP(-zetp*(+dr(iaxis)- roffset(iaxis))**2)
      t_exp_min_2=EXP(-2*zetp*(+dr(iaxis)- roffset(iaxis))*(-dr(iaxis)))
      t_exp_plus_1=EXP(-zetp*( - roffset(iaxis))**2)
      t_exp_plus_2=EXP(-2*zetp*( - roffset(iaxis))*(+dr(iaxis)))
      DO ig=0,lb_cube(iaxis),-1
        rpg = REAL(ig,dp)*dr(iaxis) - roffset(iaxis)
        t_exp_min_1=t_exp_min_1*t_exp_min_2*t_exp_1
        t_exp_min_2=t_exp_min_2*t_exp_2
        pg = t_exp_min_1
        ! pg  = EXP(-zetp*rpg**2)
        DO icoef=0,lp
           pol_z(1,icoef,ig)=pg
           pg=pg*(rpg)
        ENDDO

        rpg = REAL(1-ig,dp)*dr(iaxis) - roffset(iaxis)
        t_exp_plus_1=t_exp_plus_1*t_exp_plus_2*t_exp_1
        t_exp_plus_2=t_exp_plus_2*t_exp_2
        pg = t_exp_plus_1
        ! pg  = EXP(-zetp*rpg**2)
        DO icoef=0,lp
           pol_z(2,icoef,ig)=pg
           pg=pg*(rpg)
        ENDDO
      ENDDO

      iaxis=2
      t_exp_1=EXP(-zetp*dr(iaxis)**2)
      t_exp_2=t_exp_1**2
      t_exp_min_1=EXP(-zetp*(+dr(iaxis)- roffset(iaxis))**2)
      t_exp_min_2=EXP(-2*zetp*(+dr(iaxis)- roffset(iaxis))*(-dr(iaxis)))
      t_exp_plus_1=EXP(-zetp*( - roffset(iaxis))**2)
      t_exp_plus_2=EXP(-2*zetp*( - roffset(iaxis))*(+dr(iaxis)))
      DO ig=0,lb_cube(iaxis),-1
        rpg = REAL(ig,dp)*dr(iaxis) - roffset(iaxis)
        t_exp_min_1=t_exp_min_1*t_exp_min_2*t_exp_1
        t_exp_min_2=t_exp_min_2*t_exp_2
        pg = t_exp_min_1
        ! pg  = EXP(-zetp*rpg**2)
        DO icoef=0,lp
           pol_y(1,icoef,ig)=pg
           pg=pg*(rpg)
        ENDDO

        rpg = REAL(1-ig,dp)*dr(iaxis) - roffset(iaxis)
        t_exp_plus_1=t_exp_plus_1*t_exp_plus_2*t_exp_1
        t_exp_plus_2=t_exp_plus_2*t_exp_2
        pg = t_exp_plus_1
        ! pg  = EXP(-zetp*rpg**2)
        DO icoef=0,lp
           pol_y(2,icoef,ig)=pg
           pg=pg*(rpg)
        ENDDO
      ENDDO

      iaxis=1
      t_exp_1=EXP(-zetp*dr(iaxis)**2)
      t_exp_2=t_exp_1**2
      t_exp_min_1=EXP(-zetp*(+dr(iaxis)- roffset(iaxis))**2)
      t_exp_min_2=EXP(-2*zetp*(+dr(iaxis)- roffset(iaxis))*(-dr(iaxis)))
      t_exp_plus_1=EXP(-zetp*( - roffset(iaxis))**2)
      t_exp_plus_2=EXP(-2*zetp*( - roffset(iaxis))*(+dr(iaxis)))
      DO ig=0,lb_cube(1),-1

        rpg = REAL(ig,dp)*dr(1) - roffset(1)
        t_exp_min_1=t_exp_min_1*t_exp_min_2*t_exp_1
        t_exp_min_2=t_exp_min_2*t_exp_2
        pg = t_exp_min_1
        !pg  = EXP(-zetp*rpg**2)
        DO icoef=0,lp
           pol_x(icoef,ig)=pg
           pg=pg*(rpg)
        ENDDO

        rpg = REAL(1-ig,dp)*dr(1) - roffset(1)
        t_exp_plus_1=t_exp_plus_1*t_exp_plus_2*t_exp_1
        t_exp_plus_2=t_exp_plus_2*t_exp_2
        pg = t_exp_plus_1
        ! pg  = EXP(-zetp*rpg**2)
        DO icoef=0,lp
           pol_x(icoef,1-ig)=pg
           pg=pg*(rpg)
        ENDDO
      ENDDO


