!
!   compute the values of all (x-xp)**lp*exp(..)
!
! still requires the old trick:
!  new trick to avoid to many exps (reuse the result from the previous gridpoint):
!  exp( -a*(x+d)**2)=exp(-a*x**2)*(-2*a*x*d)*exp(-a*d**2)
!  exp(-2*a*(x+d)*d)=exp(-2*a*x*d)*exp(-2*a*d**2)

      iaxis = 3
      t_exp_1 = EXP(-zetp*dr(iaxis)**2)
      t_exp_2 = t_exp_1**2
      t_exp_min_1 = EXP(-zetp*(+dr(iaxis) - roffset(iaxis))**2)
      t_exp_min_2 = EXP(-2*zetp*(+dr(iaxis) - roffset(iaxis))*(-dr(iaxis)))
      t_exp_plus_1 = EXP(-zetp*(-roffset(iaxis))**2)
      t_exp_plus_2 = EXP(-2*zetp*(-roffset(iaxis))*(+dr(iaxis)))
      DO ig = 0, lb_cube(iaxis), -1
         rpg = REAL(ig, dp)*dr(iaxis) - roffset(iaxis)
         t_exp_min_1 = t_exp_min_1*t_exp_min_2*t_exp_1
         t_exp_min_2 = t_exp_min_2*t_exp_2
         pg = t_exp_min_1
         ! pg  = EXP(-zetp*rpg**2)
         DO icoef = 0, lp
            pol_z(1, icoef, ig) = pg
            pg = pg*(rpg)
         ENDDO

         rpg = REAL(1 - ig, dp)*dr(iaxis) - roffset(iaxis)
         t_exp_plus_1 = t_exp_plus_1*t_exp_plus_2*t_exp_1
         t_exp_plus_2 = t_exp_plus_2*t_exp_2
         pg = t_exp_plus_1
         ! pg  = EXP(-zetp*rpg**2)
         DO icoef = 0, lp
            pol_z(2, icoef, ig) = pg
            pg = pg*(rpg)
         ENDDO
      ENDDO

      iaxis = 2
      t_exp_1 = EXP(-zetp*dr(iaxis)**2)
      t_exp_2 = t_exp_1**2
      t_exp_min_1 = EXP(-zetp*(+dr(iaxis) - roffset(iaxis))**2)
      t_exp_min_2 = EXP(-2*zetp*(+dr(iaxis) - roffset(iaxis))*(-dr(iaxis)))
      t_exp_plus_1 = EXP(-zetp*(-roffset(iaxis))**2)
      t_exp_plus_2 = EXP(-2*zetp*(-roffset(iaxis))*(+dr(iaxis)))
      DO ig = 0, lb_cube(iaxis), -1
         rpg = REAL(ig, dp)*dr(iaxis) - roffset(iaxis)
         t_exp_min_1 = t_exp_min_1*t_exp_min_2*t_exp_1
         t_exp_min_2 = t_exp_min_2*t_exp_2
         pg = t_exp_min_1
         ! pg  = EXP(-zetp*rpg**2)
         DO icoef = 0, lp
            pol_y(1, icoef, ig) = pg
            pg = pg*(rpg)
         ENDDO

         rpg = REAL(1 - ig, dp)*dr(iaxis) - roffset(iaxis)
         t_exp_plus_1 = t_exp_plus_1*t_exp_plus_2*t_exp_1
         t_exp_plus_2 = t_exp_plus_2*t_exp_2
         pg = t_exp_plus_1
         ! pg  = EXP(-zetp*rpg**2)
         DO icoef = 0, lp
            pol_y(2, icoef, ig) = pg
            pg = pg*(rpg)
         ENDDO
      ENDDO

      iaxis = 1
      t_exp_1 = EXP(-zetp*dr(iaxis)**2)
      t_exp_2 = t_exp_1**2
      t_exp_min_1 = EXP(-zetp*(+dr(iaxis) - roffset(iaxis))**2)
      t_exp_min_2 = EXP(-2*zetp*(+dr(iaxis) - roffset(iaxis))*(-dr(iaxis)))
      t_exp_plus_1 = EXP(-zetp*(-roffset(iaxis))**2)
      t_exp_plus_2 = EXP(-2*zetp*(-roffset(iaxis))*(+dr(iaxis)))
      DO ig = 0, lb_cube(1), -1

         rpg = REAL(ig, dp)*dr(1) - roffset(1)
         t_exp_min_1 = t_exp_min_1*t_exp_min_2*t_exp_1
         t_exp_min_2 = t_exp_min_2*t_exp_2
         pg = t_exp_min_1
         !pg  = EXP(-zetp*rpg**2)
         DO icoef = 0, lp
            pol_x(icoef, ig) = pg
            pg = pg*(rpg)
         ENDDO

         rpg = REAL(1 - ig, dp)*dr(1) - roffset(1)
         t_exp_plus_1 = t_exp_plus_1*t_exp_plus_2*t_exp_1
         t_exp_plus_2 = t_exp_plus_2*t_exp_2
         pg = t_exp_plus_1
         ! pg  = EXP(-zetp*rpg**2)
         DO icoef = 0, lp
            pol_x(icoef, 1 - ig) = pg
            pg = pg*(rpg)
         ENDDO
      ENDDO

