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

  USE kinds,                           ONLY: dp, default_string_length
  USE mathconstants,                   ONLY: pi
  USE pair_potential_types,            ONLY: pair_potential_single_type,&
                                             pair_potential_pp_type
  USE splines,                         ONLY: init_spline,&
                                             init_splinexy,&
                                             spline_data_type,&
                                             splineh, &
                                             spline_data_release, &
                                             spline_environment_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: rho_f, rho_s
  PUBLIC :: spline_manybody_control

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

FUNCTION manybody_pot ( pot, r )


    TYPE(pair_potential_single_type), INTENT(IN)        :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    REAL(KIND=dp)                                       :: manybody_pot

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

SUBROUTINE spline_manybody_control ( spline_env, potparm )
                                  
    TYPE ( spline_environment_type ), POINTER  :: spline_env 
    TYPE ( pair_potential_pp_type ), POINTER   :: potparm 

    INTEGER                                  :: i, isos, j, jx, n, ntype
    INTEGER                                  :: npoints
    REAL(KIND=dp)                            :: a, cutmax, cutmin, dx2, e, &
                                                hicut, locut, x, x2
    LOGICAL                                  :: have_manybody
    TYPE ( spline_data_type ), POINTER       :: spline_data

!------------------------------------------------------------------------------

  npoints = spline_env % nspline_points

  cutmin = 1.0E+30_dp
  cutmax = 0.0_dp
  n = 0

  ntype = SIZE ( potparm%pot, 1 )

  n = 0
  DO i = 1, ntype ! i:  first atom type
     DO j = i, ntype ! j:  second atom type
        n = n + 1
        spline_env % spltab(i,j) = n
        spline_env % spltab(j,i) = n
! Nullify spline_data not involved in interaction
        SELECT CASE ( potparm%pot ( i, j ) % pot%type )
        CASE ( "EAM", "STILLINGER-WEBER", "REBO" )
        CASE default
          CALL spline_data_release ( spline_env % spline_data_p ( n ) % spline_data )
        END SELECT
        spline_data => spline_env % spline_data_p ( n ) % spline_data
        IF ( .NOT. ASSOCIATED ( spline_data ) ) CYCLE
!ODER: n=j+(i-1)*n-(i-1)*i/2 mit j>=i
        CALL init_splinexy( spline_data, npoints+1 )

        IF (potparm%pot(i,j)%pot%type=='EAM') THEN
           locut =  0.944863_dp ! else put locut to 0.5 Angstrom
        ELSE
           locut = 0.944863_dp  ! else put locut to 0.5 Angstrom
        END IF
        hicut = SQRT ( potparm%pot ( i, j )%pot % rcutsq )
        IF ( hicut > cutmax ) cutmax = hicut
        IF ( locut < cutmin ) cutmin = locut

        dx2 = (hicut**2-locut**2)/real(npoints)
        x2 = locut**2
        DO jx = 1, npoints + 1 ! jx: loop over distance**2
           x = SQRT ( x2 )
           e = manybody_pot ( potparm%pot ( i, j )%pot, x ) 
           spline_data % x(jx) = x2
           spline_data % y(jx) = e
           x2 = x2 + dx2
        END DO ! jx: loop over distance
        CALL init_spline(spline_data,dx=dx2)

     END DO ! j:  second atom type
  END DO ! i:  first atom type

END SUBROUTINE spline_manybody_control

!******************************************************************************

!     calculates energy=V(r^2) and fscalar=-V^{prime}(r^2)/r
!      only kept to have same args as potential()

SUBROUTINE rho_s ( spline_table, rijsq, energy, fscalar )


    REAL(KIND=dp), INTENT(IN)                    :: rijsq
    REAL(KIND=dp), INTENT(OUT)                   :: energy, fscalar
    TYPE ( spline_data_type ), INTENT ( IN )     :: spline_table

!------------------------------------------------------------------------------

  energy = splineh ( spline_table, rijsq, fscalar )
  fscalar = -2.0_dp * fscalar

END SUBROUTINE rho_s

!******************************************************************************

SUBROUTINE rho_f ( rijsq, potparm, iatom, jatom, energy )


    REAL(KIND=dp), INTENT(IN)                    :: rijsq
    TYPE(pair_potential_pp_type), POINTER        :: potparm
    INTEGER, INTENT(IN)                          :: iatom, jatom
    REAL(KIND=dp), INTENT(OUT)                   :: energy

    REAL(KIND=dp)                                :: rij

!------------------------------------------------------------------------------

  rij = SQRT ( rijsq )

! non-bond energy
  energy = manybody_pot ( potparm %pot( iatom, jatom )%pot, rij )

END SUBROUTINE rho_f

!******************************************************************************

END MODULE manybody_potential

