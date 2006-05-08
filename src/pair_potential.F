!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
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
!!
!!   SOURCE
!******************************************************************************

MODULE pair_potential

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
       get_atomic_kind
  USE erf_fn,                          ONLY: erfc
  USE ewald_environment_types,         ONLY: ewald_env_get,&
       ewald_environment_type
  USE input_constants,                 ONLY: do_ewald_none
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE pair_potential_types,            ONLY: pair_potential_pp_type,&
       pair_potential_single_type
  USE splines,                         ONLY: init_spline,&
       init_splinexy,&
       spline_data_type,&
       spline_environment_type,&
       potential_s,&
       spline_data_release,&
       spline_data_create
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: potential_f, potential_s
  PUBLIC :: spline_nonbond_control
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pair_potential'

CONTAINS


  FUNCTION ener_pot ( pot, qi, qj, r, vdw_scale14 )
    TYPE(pair_potential_single_type), &
         POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: qi, qj, r
    REAL(KIND=dp), OPTIONAL                  :: vdw_scale14
    REAL(KIND=dp)                            :: ener_pot

    INTEGER                                  :: i, index
    REAL(KIND=dp)                            :: pp, qq, scale

    IF (pot%type=='LENNARD-JONES') THEN
       IF(PRESENT(vdw_scale14)) THEN
          ener_pot = vdw_scale14 * &
               4.0_dp*pot%lj%epsilon*(pot%lj%sigma12*r**(-12)-pot%lj% &
               sigma6*r**(-6)) - pot%energy_cutoff
       ELSE
          ener_pot = 4.0_dp*pot%lj%epsilon*(pot%lj%sigma12*r**(-12)-pot%lj% &
               sigma6*r**(-6)) - pot%energy_cutoff
       END IF
    ELSE IF ( pot % type == 'IPBV' ) THEN
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
    ELSE IF (pot%type=='WILLIAMS') THEN
       ener_pot = pot%willis%a*EXP(-pot%willis%b*r) - pot%willis%c/r**6 &
            - pot%energy_cutoff
    ELSE IF (pot%type=='GOODWIN') THEN
       scale = EXP(pot%goodwin%m*(-(r/pot%goodwin%dc)**pot%goodwin%mc+(pot% &
            goodwin%d/pot%goodwin%dc)**pot%goodwin%mc))
       ener_pot = scale*pot%goodwin%vr0*(pot%goodwin%d/r)**pot%goodwin%m
    ELSE IF (pot%type=='BMHFT') THEN
       ener_pot = pot%ft%a*EXP(-pot%ft%b*r) - pot%ft%c/r**6 - pot%ft%d/r**8 &
            - pot%energy_cutoff
    ELSE IF (pot%type=='EAM') THEN
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
    ELSE 
       ener_pot = 0.0_dp
    END IF

  END FUNCTION ener_pot

  !******************************************************************************

  FUNCTION ener_coul ( pot, qi, qj, r, ewald_env, ei_scale14 )

    TYPE(pair_potential_single_type), &
         POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: qi, qj, r
    TYPE(ewald_environment_type), OPTIONAL, &
         POINTER                                :: ewald_env
    REAL(KIND=dp), OPTIONAL                  :: ei_scale14
    REAL(KIND=dp)                            :: ener_coul

    REAL(KIND=dp), PARAMETER :: eps = 1.0E-6_dp, &
         ifourpi = 1.0_dp / ( 4.0_dp * pi )

    INTEGER                                  :: ewald_type
    REAL(KIND=dp)                            :: alpha, eps0, my_ei_scale14

    ener_coul= 0.0_dp
    my_ei_scale14 = 1.0_dp
    IF (PRESENT(ei_scale14))  my_ei_scale14 = ei_scale14
    IF ( PRESENT ( ewald_env ) ) THEN
       CALL ewald_env_get ( ewald_env, alpha = alpha, eps0 = eps0, &
            ewald_type = ewald_type ) 
       IF (ewald_type /= do_ewald_none) THEN
          ener_coul = my_ei_scale14 * ifourpi / eps0 * qi * qj &
               * erfc ( alpha * r ) / r - pot % e_cutoff_coul
       ELSE
          CALL ewald_env_get ( ewald_env, eps0 = eps0)
          ener_coul = my_ei_scale14 * ifourpi / eps0 * qi * qj / r - pot % e_cutoff_coul
       END IF
    ENDIF

  END FUNCTION ener_coul

  !******************************************************************************

  ! creates the splines for the potentials

  SUBROUTINE spline_nonbond_control ( spline_env, potparm, atomic_kind_set, &
       ewald_env, ei_scale14, vdw_scale14, eps_spline, max_energy, rlow_nb,&
       emax_spline, iw, iw2, error)

    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(atomic_kind_type), DIMENSION(:), &
         POINTER                             :: atomic_kind_set
    TYPE(ewald_environment_type), OPTIONAL, &
         POINTER                             :: ewald_env
    REAL(KIND=dp), OPTIONAL                  :: ei_scale14, vdw_scale14
    REAL(KIND=dp), INTENT(IN)                :: eps_spline, max_energy, rlow_nb,&
         emax_spline
    INTEGER, INTENT(IN)                      :: iw, iw2
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_nonbond_control', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: i, j, jx, n, npoints, ntype
    INTEGER                                  :: nsize, nkx, kx, nspline
    REAL(KIND=dp)                            :: dx2, e, &
         hicut, locut, qi, qj, x, x2,&
         xdum,diffmax, xdum1, xsav, ediff
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    LOGICAL                                  :: failure, check_bounds
    
    failure = .FALSE.
    n = 0  
    ntype = SIZE ( atomic_kind_set )
    DO i = 1, ntype ! i:  first atom type

       atomic_kind=> atomic_kind_set(i)
       CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qi)

       DO j = 1, i  ! j:  second atom type

          atomic_kind=> atomic_kind_set(j)
          CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qj)
          n = n + 1
          locut = rlow_nb
          hicut = SQRT ( potparm %pot ( i, j ) %pot % rcutsq )
          IF (ABS(hicut) <= 1.0E-15_dp) hicut = 1.0E3_dp
          potparm%pot(i,j)%pot %energy_cutoff = 0.0_dp
          potparm%pot(i,j)%pot %e_cutoff_coul = 0.0_dp
          IF ( PRESENT ( vdw_scale14 ) ) THEN
             potparm%pot(i,j)%pot %energy_cutoff = ener_pot(potparm%pot(i,j)%pot,qi,qj,hicut,vdw_scale14)
          ELSE
             potparm%pot(i,j)%pot %energy_cutoff = ener_pot(potparm%pot(i,j)%pot,qi,qj,hicut)
          END IF
          IF ( PRESENT ( ei_scale14 ) ) THEN
             potparm%pot(i,j)%pot %e_cutoff_coul &
                  = ener_coul ( potparm%pot( i, j )%pot, qi, qj, hicut, ewald_env, ei_scale14 )
          ELSE
             potparm%pot(i,j)%pot %e_cutoff_coul &
                  = ener_coul ( potparm%pot ( i, j )%pot, qi, qj, hicut, ewald_env )
          ENDIF
          potparm%pot(j,i)%pot %energy_cutoff = potparm%pot(i,j)%pot %energy_cutoff
          potparm%pot(j,i)%pot %e_cutoff_coul = potparm%pot(i,j)%pot %e_cutoff_coul
          ! Find the real locut according emax_spline
          dx2 = (hicut-locut)/10000.0_dp
          x   = hicut
          DO jx = 1, 10000
             IF ( PRESENT ( ei_scale14 ) .AND. PRESENT ( vdw_scale14 ) ) THEN
                e = ener_pot ( potparm%pot( i, j )%pot, qi, qj, x ,vdw_scale14 ) &
                     + ener_coul ( potparm%pot( i, j )%pot, qi, qj, x, ewald_env, ei_scale14 )
             ELSE
                e = ener_pot ( potparm%pot( i, j )%pot, qi, qj, x ) &
                     + ener_coul ( potparm%pot( i, j )%pot, qi, qj, x, ewald_env )
             ENDIF
             IF (ABS(e) > emax_spline) THEN
                locut = x
                EXIT
             END IF
             x = x - dx2
          END DO
          npoints = 20
          IF (jx == 10001) npoints = 2  ! Zero spline
          DO WHILE (.TRUE.)

             CALL init_splinexy ( spline_env % spline_data_p ( n ) % spline_data, npoints+1)

             dx2 = (1.0_dp/locut**2-1.0_dp/hicut**2)/REAL(npoints,KIND=dp)
             x2  = 1.0_dp/hicut**2
             DO jx = 1, npoints + 1 ! jx: loop over 1/distance**2
                x = SQRT ( 1.0_dp/x2 )
                IF ( PRESENT ( ei_scale14 ) .AND. PRESENT ( vdw_scale14 ) ) THEN
                   e = ener_pot ( potparm%pot( i, j )%pot, qi, qj, x ,vdw_scale14 ) &
                        + ener_coul ( potparm%pot( i, j )%pot, qi, qj, x, ewald_env, ei_scale14 )
                ELSE
                   e = ener_pot ( potparm%pot( i, j )%pot, qi, qj, x ) &
                        + ener_coul ( potparm%pot( i, j )%pot, qi, qj, x, ewald_env )
                ENDIF
                spline_env % spline_data_p ( n ) % spline_data % x(jx) = x2
                spline_env % spline_data_p ( n ) % spline_data % y(jx) = e
                x2 = x2 + dx2
             END DO ! jx: loop over distance

             CALL init_spline(spline_env%spline_data_p(n)%spline_data,dx=dx2)

             ! This is the check for required accuracy on spline setup
             dx2 = (hicut - locut)/REAL(5*npoints+1,KIND=dp)
             x2  = locut + dx2
             diffmax = 0.0_dp
             xsav = hicut
             DO jx = 1, 5*npoints
                x = x2
                IF ( PRESENT ( ei_scale14 ) .AND. PRESENT ( vdw_scale14 ) ) THEN
                   e = ener_pot ( potparm%pot( i, j )%pot, qi, qj, x ,vdw_scale14 ) &
                        + ener_coul ( potparm%pot( i, j )%pot, qi, qj, x, ewald_env, ei_scale14 )
                ELSE
                   e = ener_pot ( potparm%pot( i, j )%pot, qi, qj, x ) &
                        + ener_coul ( potparm%pot( i, j )%pot, qi, qj, x, ewald_env )
                ENDIF

                IF (e<max_energy) THEN
                   xdum1   = ABS(e-potential_s(spline_env%spline_data_p(n)%spline_data,x*x,xdum))
                   diffmax = MAX(diffmax,xdum1)
                   xsav = MIN(x,xsav)
                END IF
                x2 = x2 + dx2
             END DO

             IF (npoints > 5E4) THEN
                IF (iw>0) THEN
                   WRITE(iw,'(A,I5,A,F15.9)') "SPLINE_INFO| Number of points: ",npoints," obtained accuracy ",diffmax
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
        
          !check spline existence
          nspline = n
          nsize = SIZE(spline_env%spline_data_p(n)%spline_data%y)
          DO kx = 1, n-1
             nkx = SIZE(spline_env%spline_data_p(kx)%spline_data%y)
             check_bounds = (spline_env%spline_data_p(kx)%spline_data%x1 == &
                             spline_env%spline_data_p(n )%spline_data%x1).AND.&
                            (spline_env%spline_data_p(kx)%spline_data%xn == &
                             spline_env%spline_data_p(n )%spline_data%xn)
             IF (nsize==nkx.AND.check_bounds) THEN
                ediff = MAXVAL(ABS(spline_env%spline_data_p(n)%spline_data%y-&
                     spline_env%spline_data_p(kx)%spline_data%y))
                IF (ediff <= eps_spline) THEN
                   CALL spline_data_release(spline_env%spline_data_p(n)%spline_data,error)
                   CALL spline_data_create ( spline_env%spline_data_p(n)%spline_data)
                   n = n - 1
                   nspline = kx
                   EXIT
                END IF
             END IF
          END DO
          spline_env % spltab ( i, j ) = nspline
          spline_env % spltab ( j, i ) = nspline
          IF (iw>0) THEN
             WRITE(iw,*)
             WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Spline number:",i*(i-1)/2+j
             WRITE(iw,'("SPLINE_INFO|",5X,A,I7,A,2I7)')"Unique Spline number",nspline," Kinds involved:",i,j
             WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Number of points used in the Splines:",npoints
             WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy requested:",eps_spline
             WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy achieved: ",diffmax
             WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX  used for spline (bohr):               ",&
                  locut,hicut
             WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX used to achieve spline accuracy (bohr):",&
                  xsav,hicut
             WRITE(iw,*)
          END IF
          ! print spline data on file
          IF (nspline==n) THEN
             IF (iw2>0) THEN
                dx2 = (hicut - locut)/REAL(5000+1,KIND=dp)
                x2  = locut + dx2
                DO jx = 1, 5000
                   x = x2
                   IF ( PRESENT ( ei_scale14 ) .AND. PRESENT ( vdw_scale14 ) ) THEN
                      e = ener_pot ( potparm%pot( i, j )%pot, qi, qj, x ,vdw_scale14 ) &
                           + ener_coul ( potparm%pot( i, j )%pot, qi, qj, x, ewald_env, ei_scale14 )
                   ELSE
                      e = ener_pot ( potparm%pot( i, j )%pot, qi, qj, x ) &
                           + ener_coul ( potparm%pot( i, j )%pot, qi, qj, x, ewald_env )
                   ENDIF
                   WRITE(1000+n,*) x, e, potential_s(spline_env%spline_data_p(n)%spline_data,x*x,xdum)
                   x2 = x2 + dx2
                END DO
                WRITE(1000+n,*)
             END IF
          END IF
          !
       END DO ! j:  second atom type
    END DO ! i:  first atom type

    IF (iw>0) THEN
       WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Number of pair potential Splines allocated:",&
            n
    END IF

  END SUBROUTINE spline_nonbond_control

  !******************************************************************************

  SUBROUTINE potential_f ( rijsq, potparm, qi, qj, iatom, jatom, energy, ewald_env )
    REAL(KIND=dp), INTENT(IN)                :: rijsq
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    REAL(KIND=dp), INTENT(IN)                :: qi, qj
    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(KIND=dp), INTENT(OUT)               :: energy
    TYPE(ewald_environment_type), OPTIONAL, &
         POINTER                                :: ewald_env

    REAL(KIND=dp)                            :: rij

    rij = SQRT ( rijsq )

    ! non-bond energy
    energy = ener_pot ( potparm%pot(iatom,jatom)%pot, qi, qj, rij )

    ! real-space ewald non-bond energy
    IF ( PRESENT ( ewald_env ) ) THEN
       energy = energy &
            + ener_coul ( potparm%pot(iatom,jatom )%pot, qi, qj, rij, ewald_env )
    END IF

  END SUBROUTINE potential_f

  !******************************************************************************

END MODULE pair_potential

