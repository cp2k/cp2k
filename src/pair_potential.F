!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/pair_potential [1.0] *
!!
!!   NAME
!!     pair_potential
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     September 2005 - Introduced the Born-Mayer-Huggins-Fumi-Tosi  Potential
!!                      BMHTF
!!     2006 - Major rewriting of the routines.. Linear scaling setup of splines
!!
!!   SOURCE
!******************************************************************************

MODULE pair_potential

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind, &
                                             get_atomic_kind_set
  USE erf_fn,                          ONLY: erfc
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE fparser,                         ONLY: EvalErrType,&
                                             evalf,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE input_constants,                 ONLY: do_ewald_none
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE pair_potential_types,            ONLY: &
       b4_type, bm_type, compare_pot, ea_type, ft_type, gp_type, gw_type, &
       ip_type, list_pot, lj_charmm_type, lj_type, nn_type, &
       pair_potential_pp_type, pair_potential_single_type, tersoff_type,wl_type
  USE shell_potential_types,           ONLY: get_shell,&
                                             shell_kind_type
  USE splines,                         ONLY: init_spline,&
                                             init_splinexy,&
                                             potential_s,&
                                             spline_env_create,&
                                             spline_environment_type
  USE string_table,                    ONLY: str2id
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  REAL(KIND=dp), PARAMETER :: ifourpi = 1.0_dp / ( 4.0_dp * pi )
  PUBLIC :: potential_f, potential_s
  PUBLIC :: spline_nonbond_control, get_nonbond_storage
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pair_potential'

CONTAINS

  !******************************************************************************
  !!****
  !!****s* pair_potential/ener_pot [1.0] *
  !!
  !!   NAME
  !!     ener_pot
  !!
  !!   FUNCTION
  !!     Evaluates the nonbond potential energy for the implemented FF kinds
  !!
  !!   AUTHOR
  !!
  !!
  !!   MODIFICATION HISTORY
  !!
  !!
  !!   SOURCE
  !******************************************************************************
  FUNCTION ener_pot ( pot, r, vdw_scale14 )
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: vdw_scale14, ener_pot

    INTEGER                                  :: i, index, n
    REAL(KIND=dp)                            :: pp, qq, scale

    IF (pot%type==lj_charmm_type) THEN
       ener_pot = vdw_scale14 * &
            4.0_dp*pot%lj%epsilon*(pot%lj%sigma12*r**(-12)-pot%lj% &
            sigma6*r**(-6)) - pot%energy_cutoff
    ELSE IF (pot%type==lj_type) THEN
       ener_pot = vdw_scale14 * pot%lj%epsilon* &
            (pot%lj%sigma12*r**(-12)-pot%lj%sigma6*r**(-6)) - pot%energy_cutoff
    ELSE IF (pot%type==ip_type) THEN
       ener_pot = 0._dp
       IF ( r > pot % ipbv % rcore ) THEN
          DO i = 2, 15
             ener_pot = ener_pot + pot % ipbv % a(i)/(r**(i-1)*REAL(i-1,dp))
          END DO
       ELSE
          ! use a linear potential
          ener_pot = pot % ipbv % m * r + pot % ipbv % b
       ENDIF
       ener_pot=ener_pot - pot % energy_cutoff
    ELSE IF (pot%type==wl_type) THEN
       ener_pot = pot%willis%a*EXP(-pot%willis%b*r) - pot%willis%c/r**6 &
            - pot%energy_cutoff
    ELSE IF (pot%type==gw_type) THEN
       scale = EXP(pot%goodwin%m*(-(r/pot%goodwin%dc)**pot%goodwin%mc+(pot% &
            goodwin%d/pot%goodwin%dc)**pot%goodwin%mc))
       ener_pot = scale*pot%goodwin%vr0*(pot%goodwin%d/r)**pot%goodwin%m
    ELSE IF (pot%type==ft_type) THEN
       ener_pot = pot%ft%a*EXP(-pot%ft%b*r) - pot%ft%c/r**6 - pot%ft%d/r**8 &
            - pot%energy_cutoff
    ELSE IF (pot%type==ea_type) THEN
       index = INT ( r / pot % eam % drar ) + 1
       IF ( index > 5000 ) THEN
          index = 5000
       ELSEIF ( index < 1 ) THEN
          index = 1
       ENDIF
       qq = r - pot % eam % rval (index)
       pp =  pot % eam % phi ( index ) + &
            qq * pot % eam % phip ( index )
       ener_pot = pp - pot % energy_cutoff
    ELSE IF (pot%type==b4_type) THEN
      IF(r<=pot%buck4r%r1) THEN
        pp = pot%buck4r%a*EXP(-pot%buck4r%b*r)
      ELSEIF(r>pot%buck4r%r1.AND.r<=pot%buck4r%r2) THEN
        pp = 0.0_dp
        DO n = 0,pot%buck4r%npoly1
          pp = pp+pot%buck4r%poly1(n)*r**n
        END DO 
      ELSEIF(r>pot%buck4r%r2.AND.r<=pot%buck4r%r3) THEN
        pp = 0.0_dp
        DO n = 0,pot%buck4r%npoly2
          pp = pp+pot%buck4r%poly2(n)*r**n
        END DO 
      ELSEIF(r>pot%buck4r%r3)THEN
        pp = -pot%buck4r%c/r**6
      END IF
      ener_pot = pp - pot % energy_cutoff
    ELSE IF (pot%type==bm_type) THEN
      ener_pot = pot%buckmo%f0*(pot%buckmo%b1+pot%buckmo%b2)*&
               EXP((pot%buckmo%a1+pot%buckmo%a2-r)/(pot%buckmo%b1+pot%buckmo%b2))&
               - pot%buckmo%c/r**6  &
               + pot%buckmo%d*(EXP(-2._dp*pot%buckmo%beta*(r-pot%buckmo%r0))- &
                               2.0_dp*EXP(-pot%buckmo%beta*(r-pot%buckmo%r0))) &
               - pot%energy_cutoff
    ELSE IF (pot%type==gp_type) THEN
       pot%gp%values(1)=r
       ener_pot = evalf(pot%gp%myid,pot%gp%values) - pot%energy_cutoff
       IF (EvalErrType > 0) &
            CALL stop_program("ener_pot","Error evaluating generic potential energy function")
    ELSE
       ener_pot = 0.0_dp
    END IF
  END FUNCTION ener_pot

  !******************************************************************************
  !!****
  !!****s* pair_potential/ener_coul [1.0] *
  !!
  !!   NAME
  !!     ener_coul
  !!
  !!   FUNCTION
  !!     Evaluates the Coulomb electrostatic energy
  !!
  !!   AUTHOR
  !!
  !!
  !!   MODIFICATION HISTORY
  !!
  !!
  !!   SOURCE
  !******************************************************************************
  FUNCTION ener_coul ( e_cutoff_coul, qi, qj, r, ewald_type,alpha,eps0, ei_scale14, my_do_14 )

!    TYPE(pair_potential_single_type), &
!      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: e_cutoff_coul, qi, qj, r
    INTEGER, INTENT(IN)                      :: ewald_type
    REAL(KIND=dp), INTENT(IN)                :: alpha, eps0, ei_scale14
    LOGICAL, INTENT(IN)                      :: my_do_14
    REAL(KIND=dp)                            :: ener_coul

    ener_coul= 0.0_dp
    IF (.NOT. my_do_14 ) THEN
       IF (ewald_type /= do_ewald_none) THEN
          ener_coul = ei_scale14 * ifourpi / eps0 * qi * qj &
               * erfc ( alpha * r ) / r -  e_cutoff_coul
       ELSE
          ener_coul = ei_scale14 * ifourpi / eps0 * qi * qj / r - e_cutoff_coul
       END IF
    ELSE
       ! to avoid a large numberof points in the spline generation
       ! we evaluate the 1-4 interactions analytically..
       ener_coul = 0.0_dp
    ENDIF

  END FUNCTION ener_coul

  !******************************************************************************
  !!****
  !!****s* pair_potential/spline_nonbond_control [1.0] *
  !!
  !!   NAME
  !!     spline_nonbond_control
  !!
  !!   FUNCTION
  !!     creates the splines for the potentials
  !!
  !!   AUTHOR
  !!
  !!
  !!   MODIFICATION HISTORY
  !!     Teo 2006.05 : Improved speed and accuracy. Linear scaling of the setup
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE spline_nonbond_control ( spline_env, potparm, atomic_kind_set, &
       ewald_env, ei_scale14, vdw_scale14, eps_spline, max_energy, rlow_nb,&
       emax_spline, iw, iw2, iw3, do_14, scale_cutoff, spline_c_coul, &
       spline_cs_coul, spline_sc_coul, error)

    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(ewald_environment_type), OPTIONAL, &
      POINTER                                :: ewald_env
    REAL(KIND=dp), OPTIONAL                  :: ei_scale14, vdw_scale14
    REAL(KIND=dp), INTENT(IN)                :: eps_spline, max_energy, &
                                                rlow_nb, emax_spline
    INTEGER, INTENT(IN)                      :: iw, iw2, iw3
    LOGICAL, INTENT(IN), OPTIONAL            :: do_14
    LOGICAL, INTENT(IN)                      :: scale_cutoff
    TYPE(spline_environment_type), POINTER, OPTIONAL   :: spline_c_coul
    TYPE(spline_environment_type), POINTER, OPTIONAL   :: spline_cs_coul
    TYPE(spline_environment_type), POINTER, OPTIONAL   :: spline_sc_coul
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_nonbond_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ewald_type, handle, i, &
                                                ilevel, j, jx, n, ncount, &
                                                ngp, npoints, ntype
    LOGICAL                                  :: failure, found_locut, is_i_shell, &
                                                is_j_shell, my_do_14
    REAL(KIND=dp) :: alpha, diffmax, dx2, e, e_cutoff_tmp, eps0, hicut, locut, locut_found, &
      my_ei_scale14, my_vdw_scale14, myfac1, myfac2, qci, qcj, qi, qj, qsi, qsj, x, x2, xdum, &
      xdum1, xsav
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(shell_kind_type), POINTER           :: shell_i, shell_j

    failure     = .FALSE.
    n      = 0
    ncount = 0
    ngp    = 0
    my_ei_scale14 = 1.0_dp ; IF (PRESENT(ei_scale14))   my_ei_scale14 = ei_scale14
    my_vdw_scale14= 1.0_dp ; IF (PRESENT(vdw_scale14)) my_vdw_scale14 = vdw_scale14
    my_do_14= .FALSE.      ; IF (PRESENT(do_14))             my_do_14 = do_14
    ntype = SIZE ( atomic_kind_set )
    CALL timeset(routineN,'I','',handle)
    IF (iw3>0) THEN
       WRITE(iw3,'(T2,A,I0,A,I0,A)') 'SPLINE_INFO| Generating: ',&
            (ntype*(ntype+1))/2,' splines due to ',ntype,' different atomic kinds'
    ENDIF
    ! Prescreen for general potential type
    DO i = 1, ntype ! i:  first atom type
       DO j = 1, i  ! j:  second atom type
          IF ((potparm%pot(i,j)%pot%type==gp_type).AND.(potparm%pot(i,j)%pot%undef)) ngp = ngp + 1
       END DO
    END DO
    CALL initf(ngp)
    ngp = 0
    DO i = 1, ntype ! i:  first atom type
       DO j = 1, i  ! j:  second atom type
          IF ((potparm%pot(i,j)%pot%type==gp_type).AND.(potparm%pot(i,j)%pot%undef)) THEN 
             ngp = ngp + 1
             potparm%pot(i,j)%pot%gp%myid = ngp
             CALL parsef(ngp,TRIM(potparm%pot(i,j)%pot%gp%potential),potparm%pot(i,j)%pot%gp%parameters)
          END IF
       END DO
    END DO
    ! real computation of splines
    DO i = 1, ntype ! i:  first atom type
       atomic_kind=> atomic_kind_set(i)
       CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qi,shell_active=is_i_shell,shell=shell_i)
       IF(is_i_shell .AND. .NOT. my_do_14) THEN
         CALL get_shell(shell=shell_i, charge_core=qci, charge_shell=qsi, error=error) 
       ELSE
         qsi = qi
         qci = 0.0_dp
       END IF
       DO j = 1, i  ! j:  second atom type
          IF (iw3>0 .AND. iw<=0) THEN
             IF (MOD(i*(i-1)/2+j,MAX(1,(ntype*(ntype+1))/(2*10)))==0) WRITE(iw3,'(T2,A3,I0)') '...',i*(i-1)/2+j
          ENDIF
          IF (.NOT.potparm%pot(i,j)%pot%undef) CYCLE
          CALL ewald_env_get(ewald_env,ewald_type=ewald_type,alpha=alpha,eps0=eps0,error=error)
          ncount = ncount + 1
          potparm%pot(i,j)%pot%undef = .FALSE.
          atomic_kind=> atomic_kind_set(j)
          CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qj,shell_active=is_j_shell,shell=shell_j)
          IF(is_j_shell .AND. .NOT. my_do_14) THEN
            CALL get_shell(shell=shell_j, charge_core=qcj, charge_shell=qsj, error=error) 
          ELSE
            qsj = qj
            qcj = 0.0_dp
          END IF
          n = spline_env % spltab(i,j)
          locut = rlow_nb
          hicut = SQRT ( potparm %pot ( i, j ) %pot % rcutsq )
          IF (ABS(hicut) <= 1.0E-15_dp) hicut = 1.0E3_dp

          potparm%pot(i,j)%pot %energy_cutoff = 0.0_dp
          potparm%pot(i,j)%pot %e_cutoff_coul = 0.0_dp

          IF (scale_cutoff) THEN
             potparm%pot(i,j)%pot %energy_cutoff = ener_pot(potparm%pot(i,j)%pot,hicut,my_vdw_scale14)
             IF((is_i_shell .OR. is_j_shell) .AND. .NOT. my_do_14) THEN
                e_cutoff_tmp &
                      = ener_coul ( potparm%pot( i, j )%pot%e_cutoff_coul, qsi, qsj, hicut, ewald_type,alpha,eps0,&
                        my_ei_scale14, my_do_14 )+ &
                        ener_coul ( potparm%pot( i, j )%pot%e_cutoff_coul, qsi, qcj, hicut, ewald_type,alpha,eps0,&
                        my_ei_scale14, my_do_14 )+ &
                        ener_coul ( potparm%pot( i, j )%pot%e_cutoff_coul, qci, qsj, hicut, ewald_type,alpha,eps0,&
                        my_ei_scale14, my_do_14 )+ &
                        ener_coul ( potparm%pot( i, j )%pot%e_cutoff_coul, qci, qcj, hicut, ewald_type,alpha,eps0,&
                        my_ei_scale14, my_do_14 )
             ELSE
               e_cutoff_tmp &
                    = ener_coul ( potparm%pot( i, j )%pot%e_cutoff_coul, qi, qj, hicut, ewald_type,alpha,eps0,&
                    my_ei_scale14, my_do_14 )
             END IF
             potparm%pot( i, j )%pot%e_cutoff_coul = e_cutoff_tmp
          END IF

          ! Find the real locut according emax_spline
          ! two pass approach, first with low resolution, refine in a second iteration
          dx2 = (hicut-locut)
          x   = hicut
          locut_found = locut
          found_locut = .FALSE.

          IF((is_i_shell .OR. is_j_shell) .AND. .NOT. my_do_14) THEN
            DO ilevel=1,2
               dx2 = dx2 / 100.0_dp
               DO jx = 1, 100
                  e = ener_pot ( potparm%pot( i, j )%pot, x , my_vdw_scale14 ) &
                    + ener_coul ( potparm%pot( i, j )%pot%e_cutoff_coul, qsi, qsj, x, ewald_type,alpha,eps0,&
                                  my_ei_scale14, my_do_14 )
                  e = e + ener_coul ( 0.0_dp, qsi, qcj, x, ewald_type,alpha,eps0, 1.0_dp, .FALSE. ) + &
                    ener_coul ( 0.0_dp, qci, qsj, x, ewald_type,alpha,eps0, 1.0_dp, .FALSE. ) +&
                    ener_coul ( 0.0_dp, qci, qcj, x, ewald_type,alpha,eps0, 1.0_dp, .FALSE. )  
                  IF (ABS(e) > emax_spline) THEN
                     locut_found = x
                     found_locut = .TRUE.
                     EXIT
                  END IF
                  x = x - dx2
               END DO
               x= x + dx2
            ENDDO
          ELSE
            DO ilevel=1,2
               dx2 = dx2 / 100.0_dp
               DO jx = 1, 100
                  e = ener_pot ( potparm%pot( i, j )%pot, x , my_vdw_scale14 ) &
                    + ener_coul ( potparm%pot( i, j )%pot%e_cutoff_coul, qsi, qsj, &
                      x, ewald_type,alpha,eps0, my_ei_scale14, my_do_14 )
                  IF (ABS(e) > emax_spline) THEN
                     locut_found = x
                     found_locut = .TRUE.
                     EXIT
                  END IF
                  x = x - dx2
               END DO
               x= x + dx2
            ENDDO
          END IF
          locut = locut_found

          npoints = 20
          IF (.NOT. found_locut) npoints = 2  ! Zero spline

          DO WHILE (.TRUE.)
             CALL init_splinexy ( spline_env % spline_data_p ( n ) % spline_data, npoints+1)
             dx2 = (1.0_dp/locut**2-1.0_dp/hicut**2)/REAL(npoints,KIND=dp)
             x2  = 1.0_dp/hicut**2
             DO jx = 1, npoints + 1     ! jx: loop over 1/distance**2
                x = SQRT ( 1.0_dp/x2 )
                e = ener_pot ( potparm%pot( i, j )%pot, x , my_vdw_scale14 ) &
                     + ener_coul ( potparm%pot( i, j )%pot%e_cutoff_coul, qsi, qsj, &
                       x, ewald_type,alpha,eps0, my_ei_scale14, my_do_14 )
                spline_env % spline_data_p ( n ) % spline_data % x(jx) = x2
                spline_env % spline_data_p ( n ) % spline_data % y(jx) = e
                x2 = x2 + dx2
             END DO ! jx: loop over distance
             IF (my_do_14) THEN
                myfac1 = my_ei_scale14 * ifourpi / eps0 * qsi * qsj
                myfac2 = 0.0_dp
                IF (scale_cutoff) THEN
                   myfac2 = my_ei_scale14 * ifourpi / eps0 * qi * qj / hicut
                END IF
                CALL init_spline(spline_env%spline_data_p(n)%spline_data,dx=dx2,&
                     ei_scale14=myfac1, ei_cutoff14=myfac2)
             ELSE
                CALL init_spline(spline_env%spline_data_p(n)%spline_data,dx=dx2)
             END IF

             ! This is the check for required accuracy on spline setup
             dx2 = (hicut - locut)/REAL(5*npoints+1,KIND=dp)
             x2  = locut + dx2
             diffmax = 0.0_dp
             xsav = hicut
             DO jx = 1, 5*npoints
                x = x2
                e = ener_pot ( potparm%pot( i, j )%pot, x , my_vdw_scale14 ) &
                     + ener_coul ( potparm%pot( i, j )%pot%e_cutoff_coul, qsi, qsj, &
                       x, ewald_type,alpha,eps0, my_ei_scale14, my_do_14 )
                IF (e<max_energy) THEN
                   xdum1   = ABS(e-potential_s(spline_env%spline_data_p(n)%spline_data,x*x,xdum))
                   diffmax = MAX(diffmax,xdum1)
                   xsav = MIN(x,xsav)
                END IF
                x2 = x2 + dx2
             END DO

             IF (npoints > 5E4) THEN
                IF (iw>0) THEN
                   WRITE(iw,'(A,I5,A,F15.9)') "SPLINE_INFO| Number of points: ",npoints,&
                        " obtained accuracy ",diffmax
                END IF
                WRITE(*,'(A)')"MM SPLINE: no convergence on required accuracy (adjust EPS_SPLINE and rerun)"
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF
             IF (diffmax > eps_spline ) THEN
                npoints   = CEILING(1.2_dp*REAL(npoints,KIND=dp))
             ELSE
                EXIT
             END IF
          END DO
          IF (iw>0) THEN
             WRITE(iw,*)
             WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Spline number:",i*(i-1)/2+j
             WRITE(iw,'("SPLINE_INFO|",5X,2(A,I7,2X),A,2I7)')"Spline number",ncount,&
                  "Unique number",n," Kinds involved:",i,j
             WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Number of points used in the Splines:",npoints
             WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy requested:",eps_spline
             WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy achieved: ",diffmax
             WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX  used for spline (bohr):               ",&
                  locut,hicut
             WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX used to achieve spline accuracy (bohr):",&
                  xsav,hicut
          END IF

          ! print spline data on file
          IF (iw2>0) THEN
             dx2 = (hicut - locut)/REAL(5000+1,KIND=dp)
             x2  = locut + dx2
             DO jx = 1, 5000
                x = x2
                e = ener_pot ( potparm%pot( i, j )%pot, x , my_vdw_scale14 ) &
                     + ener_coul ( potparm%pot( i, j )%pot%e_cutoff_coul, qsi, qsj, &
                       x, ewald_type,alpha,eps0, my_ei_scale14, my_do_14 )
                WRITE(1000+n,*) x, e, potential_s(spline_env%spline_data_p(n)%spline_data,x*x,xdum)
                x2 = x2 + dx2
             END DO
             WRITE(1000+n,*)
          END IF

!shell-model missing coulomb parts
          IF(.NOT. my_do_14) THEN
            IF(.NOT. is_i_shell ) THEN
              qci=qsi
            END IF
            IF(.NOT. is_j_shell ) THEN
              qcj=qsj
            END IF

            IF((is_i_shell .OR. is_j_shell)) THEN
!  core-core
              CPPostcondition(PRESENT(spline_c_coul),cp_failure_level,routineP,error,failure) 

              npoints = 20
              IF (.NOT. found_locut) npoints = 2  ! Zero spline

              DO WHILE (.TRUE.)
                 CALL init_splinexy ( spline_c_coul % spline_data_p ( n ) % spline_data, npoints+1)
                 dx2 = (1.0_dp/locut**2-1.0_dp/hicut**2)/REAL(npoints,KIND=dp)
                 x2  = 1.0_dp/hicut**2
                 DO jx = 1, npoints + 1     ! jx: loop over 1/distance**2
                    x = SQRT ( 1.0_dp/x2 )
                    e =  ener_coul ( 0.0_dp, qci, qcj, x, ewald_type,alpha,eps0,&
                                     1.0_dp, .FALSE. )
                    spline_c_coul % spline_data_p ( n ) % spline_data % x(jx) = x2
                    spline_c_coul % spline_data_p ( n ) % spline_data % y(jx) = e
                    x2 = x2 + dx2
                 END DO ! jx: loop over distance
                 CALL init_spline(spline_c_coul%spline_data_p(n)%spline_data,dx=dx2)

               ! This is the check for required accuracy on spline setup
                 dx2 = (hicut - locut)/REAL(5*npoints+1,KIND=dp)
                 x2  = locut + dx2
                 diffmax = 0.0_dp
                 xsav = hicut
                 DO jx = 1, 5*npoints
                    x = x2
                    e = ener_coul ( 0.0_dp, qci, qcj, x, ewald_type,alpha,eps0,&
                                    1.0_dp, .FALSE. )
                    IF (e<max_energy) THEN
                       xdum1   = ABS(e-potential_s(spline_c_coul%spline_data_p(n)%spline_data,x*x,xdum))
                       diffmax = MAX(diffmax,xdum1)
                       xsav = MIN(x,xsav)
                    END IF
                    x2 = x2 + dx2
                 END DO

                 IF (npoints > 5E4) THEN
                    IF (iw>0) THEN
                       WRITE(iw,'(A,I5,A,F15.9)') "SPLINE_INFO| Number of points: ",npoints,&
                            " obtained accuracy ",diffmax
                    END IF
                    WRITE(*,'(A)')"MM SPLINE: no convergence on required accuracy (adjust EPS_SPLINE and rerun)"
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                 END IF
                 IF (diffmax > eps_spline ) THEN
                    npoints   = CEILING(1.2_dp*REAL(npoints,KIND=dp))
                 ELSE
                    EXIT
                 END IF
              END DO
              IF (iw>0) THEN
                 WRITE(iw,*)
                 WRITE(iw,'("SPLINE_INFO|",5X,A,I7,5X, A)')"Spline number:",i*(i-1)/2+j, ' Core-Core Coulomb term in Shell-Model'
                 WRITE(iw,'("SPLINE_INFO|",5X,2(A,I7,2X),A,2I7)')"Spline number",ncount,&
                      "Unique number",n," Kinds involved:",i,j
                 WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Number of points used in the Splines:",npoints
                 WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy requested:",eps_spline
                 WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy achieved: ",diffmax
                 WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX  used for spline (bohr):               ",&
                      locut,hicut
                 WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX used to achieve spline accuracy (bohr):",&
                      xsav,hicut
              END IF

            END IF ! if at least one is shell

            IF((is_i_shell .AND. is_j_shell) ) THEN
!  core-shell
              CPPostcondition(PRESENT(spline_cs_coul),cp_failure_level,routineP,error,failure) 
              npoints = 20
              IF (.NOT. found_locut) npoints = 2  ! Zero spline
 
              DO WHILE (.TRUE.)
                 CALL init_splinexy ( spline_cs_coul % spline_data_p ( n ) % spline_data, npoints+1)
                 dx2 = (1.0_dp/locut**2-1.0_dp/hicut**2)/REAL(npoints,KIND=dp)
                 x2  = 1.0_dp/hicut**2
                 DO jx = 1, npoints + 1     ! jx: loop over 1/distance**2
                    x = SQRT ( 1.0_dp/x2 )
                    e =  ener_coul ( 0.0_dp, qci, qsj, x, ewald_type,alpha,eps0,&
                                     1.0_dp, .FALSE. )
                    spline_cs_coul % spline_data_p ( n ) % spline_data % x(jx) = x2
                    spline_cs_coul % spline_data_p ( n ) % spline_data % y(jx) = e
                    x2 = x2 + dx2
                 END DO ! jx: loop over distance
                 CALL init_spline(spline_cs_coul%spline_data_p(n)%spline_data,dx=dx2)
               ! This is the check for required accuracy on spline setup
                 dx2 = (hicut - locut)/REAL(5*npoints+1,KIND=dp)
                 x2  = locut + dx2
                 diffmax = 0.0_dp
                 xsav = hicut
                 DO jx = 1, 5*npoints
                    x = x2
                    e = ener_coul ( 0.0_dp, qci, qsj, x, ewald_type,alpha,eps0,&
                                    1.0_dp, .FALSE. )
                    IF (e<max_energy) THEN
                       xdum1   = ABS(e-potential_s(spline_cs_coul%spline_data_p(n)%spline_data,x*x,xdum))
                       diffmax = MAX(diffmax,xdum1)
                       xsav = MIN(x,xsav)
                    END IF
                    x2 = x2 + dx2
                 END DO
                 IF (npoints > 5E4) THEN
                    IF (iw>0) THEN
                       WRITE(iw,'(A,I5,A,F15.9)') "SPLINE_INFO| Number of points: ",npoints,&
                            " obtained accuracy ",diffmax
                    END IF
                    WRITE(*,'(A)')"MM SPLINE: no convergence on required accuracy (adjust EPS_SPLINE and rerun)"
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                 END IF
                 IF (diffmax > eps_spline ) THEN
                    npoints   = CEILING(1.2_dp*REAL(npoints,KIND=dp))
                 ELSE
                    EXIT
                 END IF
              END DO

              IF (iw>0) THEN
                 WRITE(iw,*)
                 WRITE(iw,'("SPLINE_INFO|",5X,A,I7,5X, A)')"Spline number:",i*(i-1)/2+j, ' Core-Shell Coulomb term in Shell-Model'
                 WRITE(iw,'("SPLINE_INFO|",5X,2(A,I7,2X),A,2I7)')"Spline number",ncount,&
                      "Unique number",n," Kinds involved:",i,j
                 WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Number of points used in the Splines:",npoints
                 WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy requested:",eps_spline
                 WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy achieved: ",diffmax
                 WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX  used for spline (bohr):               ",&
                      locut,hicut
                 WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX used to achieve spline accuracy (bohr):",&
                      xsav,hicut
              END IF
!  shell-core
              CPPostcondition(PRESENT(spline_sc_coul),cp_failure_level,routineP,error,failure) 
              npoints = 20
              IF (.NOT. found_locut) npoints = 2  ! Zero spline

              DO WHILE (.TRUE.)
                 CALL init_splinexy ( spline_sc_coul % spline_data_p ( n ) % spline_data, npoints+1)
                 dx2 = (1.0_dp/locut**2-1.0_dp/hicut**2)/REAL(npoints,KIND=dp)
                 x2  = 1.0_dp/hicut**2
                 DO jx = 1, npoints + 1     ! jx: loop over 1/distance**2
                    x = SQRT ( 1.0_dp/x2 )
                    e =  ener_coul ( 0.0_dp, qsi, qcj, x, ewald_type,alpha,eps0,&
                                     1.0_dp, .FALSE. )
                    spline_sc_coul % spline_data_p ( n ) % spline_data % x(jx) = x2
                    spline_sc_coul % spline_data_p ( n ) % spline_data % y(jx) = e
                    x2 = x2 + dx2
                 END DO ! jx: loop over distance
                 CALL init_spline(spline_sc_coul%spline_data_p(n)%spline_data,dx=dx2)
               ! This is the check for required accuracy on spline setup
                 dx2 = (hicut - locut)/REAL(5*npoints+1,KIND=dp)
                 x2  = locut + dx2
                 diffmax = 0.0_dp
                 xsav = hicut
                 DO jx = 1, 5*npoints
                    x = x2
                    e = ener_coul ( 0.0_dp, qsi, qcj, x, ewald_type,alpha,eps0,&
                                    1.0_dp, .FALSE. )
                    IF (e<max_energy) THEN
                       xdum1   = ABS(e-potential_s(spline_sc_coul%spline_data_p(n)%spline_data,x*x,xdum))
                       diffmax = MAX(diffmax,xdum1)
                       xsav = MIN(x,xsav)
                    END IF
                    x2 = x2 + dx2
                 END DO
                 IF (npoints > 5E4) THEN
                    IF (iw>0) THEN
                       WRITE(iw,'(A,I5,A,F15.9)') "SPLINE_INFO| Number of points: ",npoints,&
                            " obtained accuracy ",diffmax
                    END IF
                    WRITE(*,'(A)')"MM SPLINE: no convergence on required accuracy (adjust EPS_SPLINE and rerun)"
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                 END IF
                 IF (diffmax > eps_spline ) THEN
                    npoints   = CEILING(1.2_dp*REAL(npoints,KIND=dp))
                 ELSE
                    EXIT
                 END IF
              END DO
              IF (iw>0) THEN
                 WRITE(iw,*)
                 WRITE(iw,'("SPLINE_INFO|",5X,A,I7,5X, A)')"Spline number:",i*(i-1)/2+j, ' Shell-Core Coulomb term in Shell-Model'
                 WRITE(iw,'("SPLINE_INFO|",5X,2(A,I7,2X),A,2I7)')"Spline number",ncount,&
                      "Unique number",n," Kinds involved:",i,j
                 WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Number of points used in the Splines:",npoints
                 WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy requested:",eps_spline
                 WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy achieved: ",diffmax
                 WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX  used for spline (bohr):               ",&
                      locut,hicut
                 WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX used to achieve spline accuracy (bohr):",&
                      xsav,hicut
              END IF
            END IF  ! if both are shell

          END IF  ! my_do_14

       END DO  ! jkind
    END DO   !ikind
    CALL finalizef()
    IF (iw3>0) THEN
       WRITE(iw3,'()')
    ENDIF
    IF (iw>0) THEN
       WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Number of pair potential Splines allocated:",&
            MAXVAL(spline_env % spltab)
    END IF
    IF (iw3>0) THEN
       WRITE(iw3,'(T2,A,I7)') 'SPLINE_INFO| Number of unique splines computed:',&
            MAXVAL(spline_env % spltab)
       WRITE(iw3,'(T2,A,I0)') 'SPLINE_INFO| Done'
    ENDIF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE spline_nonbond_control

  !******************************************************************************
  !!****
  !!****s* pair_potential/get_nonbond_storage [1.0] *
  !!
  !!   NAME
  !!     get_nonbond_storage
  !! 
  !!   FUNCTION
  !!     Prescreening of the effective bonds evaluations. linear scaling algorithm
  !!
  !!   AUTHOR
  !!     Teodoro Laino [tlaino] 2006.05
  !!
  !!   MODIFICATION HISTORY
  !!
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE get_nonbond_storage ( spline_env, potparm, atomic_kind_set, &
       spline_c_coul_env, spline_cs_coul_env, spline_sc_coul_env,&
       ei_scale14, vdw_scale14, error)

    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(spline_environment_type), POINTER, OPTIONAL   :: spline_c_coul_env
    TYPE(spline_environment_type), POINTER, OPTIONAL   :: spline_cs_coul_env
    TYPE(spline_environment_type), POINTER, OPTIONAL   :: spline_sc_coul_env
    REAL(KIND=dp), OPTIONAL                  :: ei_scale14, vdw_scale14
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_nonbond_storage', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, idim, iend, istart, j, k, locij, min_val, n, ndim, &
      nk, ntype, nunique, nvar, pot_target, stat, tmpij(2), tmpij0(2)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: Iwork1, Iwork2, my_index
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: tmp_index
    LOGICAL                                  :: at_least_one, check, failure, &
                                                first_val, shell_present
    REAL(KIND=dp)                            :: my_ei_scale14, &
                                                my_vdw_scale14, qi, qj
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Cwork, Rwork, wtmp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pot_par
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    failure     = .FALSE.
    CALL timeset(routineN,'I','',handle)
    ntype = SIZE ( atomic_kind_set )
    my_ei_scale14 = 1.0_dp ; IF (PRESENT(ei_scale14))  my_ei_scale14  = ei_scale14
    my_vdw_scale14= 1.0_dp ; IF (PRESENT(vdw_scale14)) my_vdw_scale14 = vdw_scale14
    DO i = 1, ntype
       DO j = 1, i
          potparm%pot(i,j)%pot%undef = .FALSE.
       END DO
    END DO
    ALLOCATE(tmp_index(ntype,ntype),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !
    nunique = 0
    tmp_index = HUGE(0)
    DO pot_target = MINVAL(list_pot), MAXVAL(list_pot)
       ndim  = 0
       first_val = .FALSE.
       n     = 0
       DO i = 1, ntype
          DO j = 1, i
             n = n + 1
             IF (potparm%pot(i,j)%pot%type==pot_target) THEN
                IF (.NOT.first_val) THEN
                   min_val = n
                   first_val = .TRUE.
                END IF
                tmp_index(i,j) = 1
                tmp_index(j,i) = 1
                ndim = ndim + 1
             END IF
          END DO
       END DO

       IF (ndim==0) CYCLE ! No potential of this kind found
       CALL get_indexes(min_val,ntype,tmpij)
       potparm%pot(tmpij(1),tmpij(2))%pot%undef = .TRUE.
       n     = 1
       ALLOCATE(my_index(ndim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       nvar = 0
       IF (my_ei_scale14 /=0.0_dp) nvar = nvar + 2
       SELECT CASE(pot_target)
       CASE (lj_type,lj_charmm_type)
          IF (my_vdw_scale14/=0.0_dp) nvar = 3 + nvar
       CASE (wl_type)
          nvar = 3 + nvar
       CASE (gw_type)
          nvar = 5 + nvar
       CASE (ea_type)
          nvar = 3 + nvar
       CASE (ft_type)
          nvar = 4 + nvar
       CASE (ip_type)
          nvar = 3 + nvar
       CASE (b4_type)
          nvar = 6 + nvar
       CASE (bm_type)
          nvar = 9 + nvar
       CASE (gp_type)
          nvar = 2 + nvar
       CASE (tersoff_type)
          nvar = 13 + nvar
       CASE (nn_type)
          nvar = nvar
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       IF (nvar /=0) THEN
          ALLOCATE(pot_par(ndim,nvar),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          n  = 0
          nk = 0
          DO i = 1, ntype
             atomic_kind=> atomic_kind_set(i)
             CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qi)
             DO j = 1, i
                atomic_kind=> atomic_kind_set(j)
                CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qj)
                n = n + 1
                IF (potparm%pot(i,j)%pot%type==pot_target) THEN
                   nk = nk + 1
                   my_index(nk) = n
                   SELECT CASE(pot_target)
                   CASE (lj_type,lj_charmm_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%lj%epsilon
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%lj%sigma6
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%lj%sigma12
                   CASE(gp_type)
                      pot_par(nk,1)  = str2id(potparm%pot(i,j)%pot%gp%potential)
                      pot_par(nk,2)  = str2id(potparm%pot(i,j)%pot%gp%variables)
                   CASE (wl_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%willis%a
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%willis%b
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%willis%c
                   CASE (gw_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%goodwin%vr0
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%goodwin%m
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%goodwin%mc
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%goodwin%d
                      pot_par(nk,5)  = potparm%pot(i,j)%pot%goodwin%dc
                   CASE (ea_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%eam%drar
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%eam%drhoar
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%eam%acutal
                   CASE (ft_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%ft%A
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%ft%B
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%ft%C
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%ft%D
                   CASE (ip_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%ipbv%rcore
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%ipbv%m
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%ipbv%b
                   CASE (b4_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%buck4r%a
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%buck4r%b
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%buck4r%c
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%buck4r%r1
                      pot_par(nk,5)  = potparm%pot(i,j)%pot%buck4r%r2
                      pot_par(nk,6)  = potparm%pot(i,j)%pot%buck4r%r3
                   CASE (bm_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%buckmo%f0
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%buckmo%a1
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%buckmo%a2
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%buckmo%b1
                      pot_par(nk,5)  = potparm%pot(i,j)%pot%buckmo%b2
                      pot_par(nk,6)  = potparm%pot(i,j)%pot%buckmo%c
                      pot_par(nk,7)  = potparm%pot(i,j)%pot%buckmo%d
                      pot_par(nk,8)  = potparm%pot(i,j)%pot%buckmo%r0
                      pot_par(nk,9)  = potparm%pot(i,j)%pot%buckmo%beta
                   CASE (tersoff_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%tersoff%A
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%tersoff%B
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%tersoff%lambda1
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%tersoff%lambda2
                      pot_par(nk,5)  = potparm%pot(i,j)%pot%tersoff%alpha
                      pot_par(nk,6)  = potparm%pot(i,j)%pot%tersoff%beta
                      pot_par(nk,7)  = potparm%pot(i,j)%pot%tersoff%n
                      pot_par(nk,8)  = potparm%pot(i,j)%pot%tersoff%c
                      pot_par(nk,9)  = potparm%pot(i,j)%pot%tersoff%d
                      pot_par(nk,10) = potparm%pot(i,j)%pot%tersoff%h
                      pot_par(nk,11) = potparm%pot(i,j)%pot%tersoff%lambda3
                      pot_par(nk,12) = potparm%pot(i,j)%pot%tersoff%bigR
                      pot_par(nk,13) = potparm%pot(i,j)%pot%tersoff%bigD
                   CASE (nn_type)
                      ! no checks
                   CASE DEFAULT
                      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                   END SELECT
                   IF (my_ei_scale14 /=0.0_dp) THEN
                      pot_par(nk,nvar-1)=qi
                      pot_par(nk,nvar  )=qj
                   END IF
                END IF
             END DO
          END DO
          !
          ! MAIN SORTING LOOP
          !
          ALLOCATE(Rwork(ndim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(Iwork1(ndim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(Iwork2(ndim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(wtmp(nvar),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL sort(pot_par(:,1),ndim,Iwork1)
          ! Sort all the other components of the potential
          DO k = 2, nvar
             Rwork = pot_par(:,k)
             DO i = 1, ndim
                pot_par(i,k)=Rwork(Iwork1(i))
             END DO
          END DO
          Iwork2  = my_index
          DO i = 1, ndim
             my_index(i) = Iwork2(Iwork1(i))
          END DO
          ! Iterative sorting
          DO k = 2, nvar
             wtmp(1:k-1) = pot_par(1,1:k-1)
             istart   = 1
             at_least_one = .FALSE.
             DO j = 1, ndim
                Rwork(j)  = pot_par(j,k)
                IF (ALL(pot_par(j,1:k-1)==wtmp(1:k-1))) CYCLE
                iend   = j - 1
                wtmp(1:k-1) = pot_par(j,1:k-1)
                ! If the ordered array has no two same consecutive elements
                ! does not make any sense to proceed ordering the others
                ! related parameters..
                idim = iend-istart+1
                CALL sort(Rwork(istart:iend),idim,Iwork1(istart:iend))
                Iwork1(istart:iend)=Iwork1(istart:iend)-1+istart
                IF (idim/=1) at_least_one = .TRUE.
                istart = j
             END DO
             iend = ndim
             idim = iend-istart+1
             CALL sort(Rwork(istart:iend),idim,Iwork1(istart:iend))
             Iwork1(istart:iend)=Iwork1(istart:iend)-1+istart
             IF (idim/=1) at_least_one = .TRUE.
             pot_par(:,k) = Rwork
             IF (.NOT.at_least_one) EXIT
             ! Sort other components
             DO j = k+1, nvar
                Rwork = pot_par(:,j)
                DO i = 1, ndim
                   pot_par(i,j)=Rwork(Iwork1(i))
                END DO
             END DO
             Iwork2  = my_index
             DO i = 1, ndim
                my_index(i) = Iwork2(Iwork1(i))
             END DO
          END DO
          DEALLOCATE(wtmp,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(Iwork1,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(Iwork2,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(Rwork,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          !
          ! Let's determine the number of unique potentials and tag them
          !
          ALLOCATE(Cwork(nvar),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          Cwork = pot_par(1,:)
          locij = my_index(1)
          CALL get_indexes(locij,ntype,tmpij0)
          istart = 1
          DO j = 1, ndim
             ! Special cases for EAM and IPBV
             locij = my_index(j)
             CALL get_indexes(locij,ntype,tmpij)
             SELECT CASE(pot_target)
             CASE(ea_type,ip_type)
                ! check the array components
                CALL compare_pot(potparm%pot(tmpij(1) ,tmpij(2) )%pot,&
                                 potparm%pot(tmpij0(1),tmpij0(2))%pot,&
                                 check, error)
             CASE(gp_type)
                check = .TRUE.
                IF  (ASSOCIATED(potparm%pot(tmpij(1) ,tmpij(2) )%pot%gp%parameters).AND.&
                     ASSOCIATED(potparm%pot(tmpij0(1),tmpij0(2))%pot%gp%parameters)) THEN
                   IF  (SIZE(potparm%pot(tmpij(1) ,tmpij(2) )%pot%gp%parameters)==&
                        SIZE(potparm%pot(tmpij0(1),tmpij0(2))%pot%gp%parameters)) THEN
                      IF  (ANY(potparm%pot(tmpij(1) ,tmpij(2) )%pot%gp%parameters /= &
                               potparm%pot(tmpij0(1),tmpij0(2))%pot%gp%parameters)) check = .FALSE.
                   END IF
                END IF
                IF  (ASSOCIATED(potparm%pot(tmpij(1) ,tmpij(2) )%pot%gp%values).AND.&
                     ASSOCIATED(potparm%pot(tmpij0(1),tmpij0(2))%pot%gp%values)) THEN
                   IF  (SIZE(potparm%pot(tmpij(1) ,tmpij(2) )%pot%gp%values)==&
                        SIZE(potparm%pot(tmpij0(1),tmpij0(2))%pot%gp%values)) THEN
                      IF (ANY(potparm%pot(tmpij(1) ,tmpij(2) )%pot%gp%values /= &
                              potparm%pot(tmpij0(1),tmpij0(2))%pot%gp%values))     check = .FALSE.
                   END IF
                END IF
             CASE default
                check = .TRUE.
             END SELECT
             IF (ALL(Cwork==pot_par(j,:)).AND.check) CYCLE
             Cwork = pot_par(j,:)
             nunique = nunique + 1
             iend    = j - 1
             min_val = MINVAL(my_index(istart:iend))
             CALL get_indexes(min_val,ntype,tmpij)
             potparm%pot(tmpij(1),tmpij(2))%pot%undef =.TRUE.
             potparm%pot(tmpij(2),tmpij(1))%pot%undef =.TRUE.
             DO i = istart, iend
                locij = my_index(i)
                CALL get_indexes(locij,ntype,tmpij)
                tmp_index(tmpij(1),tmpij(2))=nunique
                tmp_index(tmpij(2),tmpij(1))=nunique
             END DO
             istart = j
             locij = my_index(j)
             CALL get_indexes(locij,ntype,tmpij0)
          END DO
          nunique = nunique + 1
          iend    = ndim
          min_val = MINVAL(my_index(istart:iend))
          CALL get_indexes(min_val,ntype,tmpij)
          potparm%pot(tmpij(1),tmpij(2))%pot%undef =.TRUE.
          potparm%pot(tmpij(2),tmpij(1))%pot%undef =.TRUE.
          DO i = istart, iend
             locij = my_index(i)
             CALL get_indexes(locij,ntype,tmpij)
             tmp_index(tmpij(1),tmpij(2))=nunique
             tmp_index(tmpij(2),tmpij(1))=nunique
          END DO
          DEALLOCATE(Cwork,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(pot_par,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          nunique = nunique + 1
       END IF
       DEALLOCATE(my_index,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO  ! pot_target

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, shell_present=shell_present)
    CALL spline_env_create ( spline_env, ntype, ntab_in=nunique, error=error )
    spline_env%spltab=tmp_index
    IF(shell_present) THEN
      IF(PRESENT(spline_c_coul_env)) THEN
        CALL spline_env_create ( spline_c_coul_env, ntype, ntab_in=nunique, error=error )
        spline_c_coul_env%spltab=tmp_index
      END IF
      IF(PRESENT(spline_cs_coul_env)) THEN
        CALL spline_env_create ( spline_cs_coul_env, ntype, ntab_in=nunique, error=error )
        spline_cs_coul_env%spltab=tmp_index
      END IF
      IF(PRESENT(spline_sc_coul_env)) THEN
        CALL spline_env_create ( spline_sc_coul_env, ntype, ntab_in=nunique, error=error )
        spline_sc_coul_env%spltab=tmp_index
      END IF
    END IF
    DEALLOCATE(tmp_index,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE get_nonbond_storage

  !******************************************************************************
  !!****
  !!****s* pair_potential/get_indexes [1.0] *
  !!
  !!   NAME
  !!     get_indexes
  !!
  !!   FUNCTION
  !!     Gives back the indices of the matrix w.r.t. the collective array index
  !!
  !!   AUTHOR
  !!     Teodoro Laino [tlaino] 2006.05
  !!
  !!   MODIFICATION HISTORY
  !!
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE get_indexes(Inind, ndim, ij)
    INTEGER, INTENT(IN)                      :: Inind, ndim
    INTEGER, DIMENSION(2), INTENT(OUT)       :: ij

    INTEGER                                  :: i, tmp

    tmp = 0
    ij  = HUGE(0)
    DO i = 1, ndim
       tmp = tmp + i
       IF (tmp >= Inind) THEN
          ij(1) = i
          ij(2) = Inind - tmp + i
          EXIT
       END IF
    END DO
  END SUBROUTINE get_indexes

  !******************************************************************************
  !!****
  !!****s* pair_potential/potential_f [1.0] *
  !!
  !!   NAME
  !!     potential_f
  !!
  !!   FUNCTION
  !!    Evaluates the nonbond and electrostaitic potential analytically
  !!
  !!   AUTHOR
  !!     Teodoro Laino [tlaino] 2006.05
  !!
  !!   MODIFICATION HISTORY
  !!
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE potential_f ( rijsq, potparm, qi, qj, iatom, jatom, energy, ewald_env, error )
    REAL(KIND=dp), INTENT(IN)                :: rijsq
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    REAL(KIND=dp), INTENT(IN)                :: qi, qj
    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(KIND=dp), INTENT(OUT)               :: energy
    TYPE(ewald_environment_type), OPTIONAL, &
      POINTER                                :: ewald_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ewald_type
    REAL(KIND=dp)                            :: alpha, eps0, rij

    rij = SQRT ( rijsq )

    ! non-bond energy
    energy = ener_pot ( potparm%pot(iatom,jatom)%pot, rij, 1.0_dp)

    ! real-space ewald non-bond energy
    IF ( PRESENT ( ewald_env ) ) THEN
       CALL ewald_env_get(ewald_env,ewald_type=ewald_type,alpha=alpha,eps0=eps0,error=error)
       energy = energy &
            + ener_coul ( potparm%pot(iatom,jatom )%pot%e_cutoff_coul, qi, qj, rij, ewald_type,alpha,eps0, 1.0_dp,&
            my_do_14=.FALSE.)
    END IF

  END SUBROUTINE potential_f

END MODULE pair_potential

