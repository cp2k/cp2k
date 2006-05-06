!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/manybody_potential [1.0] *
!!
!!   NAME
!!     manybody_potential
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE manybody_potential

  USE kinds,                           ONLY: dp
  USE pair_potential_types,            ONLY: pair_potential_pp_type,&
       pair_potential_single_type
  USE splines,                         ONLY: init_spline,&
       init_splinexy,&
       spline_data_release,&
       spline_data_type,&
       spline_environment_type,&
       splineh,&
       spline_data_release,&
       spline_data_create
  USE termination,                     ONLY: stop_program

#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: rho_f, rho_s
  PUBLIC :: spline_manybody_control
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_potential'

CONTAINS

  FUNCTION manybody_pot ( pot, r )

    TYPE(pair_potential_single_type), &
         INTENT(IN)                             :: pot
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: manybody_pot

    IF (pot%type=='EAM') THEN
       ! No splining necessary for EAM.  Zero splines
       manybody_pot = 0._dp
    ELSE
       WRITE(*,*) "pot%type = ",pot%type
       CALL stop_program("manybody_potential%ener_pot","Potential type unknown"//&
            " pot%type = "//TRIM(pot%type))
    END IF

  END FUNCTION manybody_pot

  !******************************************************************************

  ! creates the splines for the manybody contribution

  SUBROUTINE spline_manybody_control ( spline_env, potparm, eps_spline, max_energy, rlow_nb,&
       emax_spline, iw, iw2, error)

    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    REAL(KIND=dp), INTENT(IN)                :: eps_spline, max_energy, rlow_nb,&
         emax_spline
    INTEGER, INTENT(IN)                      :: iw, iw2
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_manybody_control', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: i, j, jx, n, npoints, ntype
    INTEGER                                  :: nspline, nkx, kx, nsize
    REAL(KIND=dp)                            :: dx2, e, &
         hicut, locut, x, x2,  xdum,diffmax,xdum1,ediff,xsav
    TYPE(spline_data_type), POINTER          :: spline_data
    LOGICAL                                  :: failure, check_bounds

    failure = .FALSE.
    n = 0
    ntype = SIZE ( potparm%pot, 1 )
    n = 0
    DO i = 1, ntype ! i:  first atom type
       DO j = 1, i  ! j:  second atom type
          n = n + 1
          nspline = n
          SELECT CASE ( potparm%pot ( i, j ) % pot%type )
          CASE ( "EAM", "STILLINGER-WEBER", "REBO" )
             spline_data => spline_env % spline_data_p ( n ) % spline_data
             locut = rlow_nb
             hicut = SQRT ( potparm%pot ( i, j )%pot % rcutsq )
             IF (ABS(hicut) <= 1.0E-15_dp) hicut = 1.0E3_dp
             ! Find the real locut according emax_spline
             dx2 = (hicut-locut)/10000.0_dp
             x   = hicut
             DO jx = 1, 10000
                e = manybody_pot ( potparm%pot ( i, j )%pot, x ) 
                IF (ABS(e) > emax_spline) THEN
                   locut = x
                   EXIT
                END IF
                x = x - dx2
             END DO
             npoints = 20
             IF (jx == 10001) npoints = 2  ! Zero spline
             DO WHILE (.TRUE.)

                CALL init_splinexy( spline_data, npoints+1 )

                dx2 = (1.0_dp/locut**2-1.0_dp/hicut**2)/REAL(npoints,KIND=dp)
                x2  = 1.0_dp/hicut**2
                DO jx = 1, npoints + 1 ! jx: loop over 1/distance**2
                   x = SQRT ( 1.0_dp/x2 )
                   e = manybody_pot ( potparm%pot ( i, j )%pot, x ) 
                   spline_data % x(jx) = x2
                   spline_data % y(jx) = e
                   x2 = x2 + dx2
                END DO

                CALL init_spline(spline_data,dx=dx2)

                ! This is the check for required accuracy on spline setup
                dx2 = (hicut - locut)/REAL(5*npoints+1,KIND=dp)
                x2  = locut + dx2
                diffmax = 0.0_dp
                xsav = hicut
                DO jx = 1, 5*npoints
                   x = x2
                   e = manybody_pot ( potparm%pot ( i, j )%pot, x ) 
                   IF (e<max_energy) THEN
                      xdum1   = ABS(e-splineh(spline_env%spline_data_p(n)%spline_data,x*x,xdum))
                      diffmax = MAX(diffmax,xdum1)
                      xsav = MIN(x,xsav)
                   END IF
                   x2 = x2 + dx2
                END DO
                IF (npoints > 5E4) THEN
                   IF (iw>0) THEN
                      WRITE(iw,'(A,I5,A,F15.9)') "SPLINE_INFO| Number of points: ",&
                           npoints," obtained accuracy ",diffmax
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
             nsize = SIZE(spline_env%spline_data_p(n)%spline_data%y)
             DO kx = 1, n-1
                nkx = SIZE(spline_env%spline_data_p(kx)%spline_data%y)
                check_bounds = (spline_env%spline_data_p(kx)%spline_data%x1 == &
                                spline_env%spline_data_p(n )%spline_data%x1).AND.&
                               (spline_env%spline_data_p(kx)%spline_data%xn == &
                                spline_env%spline_data_p(n )%spline_data%xn)                
                IF (nsize==nkx) THEN
                   ediff = MAXVAL(ABS(spline_env%spline_data_p(n)%spline_data%y-&
                        spline_env%spline_data_p(kx)%spline_data%y))
                   IF (ediff <= eps_spline) THEN
                      CALL spline_data_release(spline_env%spline_data_p(n)%spline_data,&
                           error)
                      CALL spline_data_create ( spline_env%spline_data_p(n)%spline_data)
                      n = n - 1
                      nspline = kx
                      EXIT
                   END IF
                END IF
             END DO
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
                      e = manybody_pot ( potparm%pot ( i, j )%pot, x )
                      WRITE(2000+n,*) x, e, splineh(spline_env%spline_data_p(n)%spline_data,x*x,xdum)
                      x2 = x2 + dx2
                   END DO
                   WRITE(2000+n,*)
                END IF
             END IF

          CASE default
             ! Nullify spline_data not involved in interaction
             CALL spline_data_release ( spline_env % spline_data_p ( n ) % spline_data )
             n = n - 1
             nspline = n
          END SELECT
          spline_env % spltab ( i, j ) = nspline
          spline_env % spltab ( j, i ) = nspline          

       END DO ! j:  second atom type
    END DO ! i:  first atom type

    IF (iw>0) THEN
       WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Number of manybody potential Splines allocated:",&
            n
    END IF

  END SUBROUTINE spline_manybody_control

  !******************************************************************************

  !     calculates energy=V(r^2) and fscalar=-V^{prime}(r^2)/r
  !      only kept to have same args as potential()

  SUBROUTINE rho_s ( spline_table, rijsq, energy, fscalar )


    TYPE(spline_data_type), POINTER          :: spline_table
    REAL(KIND=dp), INTENT(IN)                :: rijsq
    REAL(KIND=dp), INTENT(OUT)               :: energy, fscalar

    !------------------------------------------------------------------------------

    energy = splineh ( spline_table, rijsq, fscalar )
    fscalar = -2.0_dp * fscalar

  END SUBROUTINE rho_s

  !******************************************************************************

  SUBROUTINE rho_f ( rijsq, potparm, iatom, jatom, energy )


    REAL(KIND=dp), INTENT(IN)                :: rijsq
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(KIND=dp), INTENT(OUT)               :: energy

    REAL(KIND=dp)                            :: rij

    !------------------------------------------------------------------------------

    rij = SQRT ( rijsq )

    ! non-bond energy
    energy = manybody_pot ( potparm %pot( iatom, jatom )%pot, rij )

  END SUBROUTINE rho_f

  !******************************************************************************

END MODULE manybody_potential

