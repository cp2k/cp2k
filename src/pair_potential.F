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
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE pair_potential

  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : zero, pi
  USE particle_types, ONLY : particle_prop_type
  USE splines, ONLY : init_spline, init_splinexy,  &
       spline, spline_1, kill_spline, spline_data, splineh
  USE termination, ONLY : stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_potrange, potential_f, potential_s
  PUBLIC :: potentialparm_type, spline_nonbond_control
  PUBLIC :: lj_pot_type, williams_pot_type, goodwin_pot_type

  TYPE lj_pot_type
    REAL ( dbl ) :: epsilon
    REAL ( dbl ) :: sigma6
    REAL ( dbl ) :: sigma12
  END TYPE lj_pot_type

  TYPE williams_pot_type
    REAL ( dbl ) :: a
    REAL ( dbl ) :: b
    REAL ( dbl ) :: c
  END TYPE williams_pot_type

  TYPE goodwin_pot_type
    REAL ( dbl ) :: vr0
    REAL ( dbl ) :: m, mc
    REAL ( dbl ) :: d, dc
  END TYPE goodwin_pot_type

  TYPE potentialparm_type
    CHARACTER ( LEN = 40 ) :: type
    TYPE (lj_pot_type ) :: lj
    TYPE (williams_pot_type ) :: willis
    TYPE (goodwin_pot_type ) :: goodwin
    REAL ( dbl ) :: rcutsq
    REAL ( dbl ) :: energy_cutoff
    REAL ( dbl ) :: e_cutoff_coul
    REAL ( dbl ) :: lrc
  END TYPE potentialparm_type

! look-up table (i,j)->n
  INTEGER, DIMENSION ( :, : ), ALLOCATABLE, SAVE :: spltab

! splines for interaction n
  TYPE ( spline_data ), DIMENSION ( : ), POINTER, SAVE :: spline_table

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

FUNCTION get_potrange ( potparm ) RESULT ( potrange )

  IMPLICIT NONE

! Result
  REAL ( dbl ) :: potrange

! Arguments
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), INTENT ( IN ) :: potparm

!------------------------------------------------------------------------------

  potrange = MAXVAL( abs ( potparm ( :, : ) % rcutsq ) )
  potrange = SQRT ( potrange )

END FUNCTION get_potrange

!******************************************************************************

FUNCTION ener_pot ( pot, qi, qj, r )

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), INTENT ( IN ) :: qi, qj, r
  TYPE (potentialparm_type ), INTENT ( IN ) :: pot

! Locals
  REAL ( dbl ) :: ener_pot, scale

!------------------------------------------------------------------------------

  IF (pot%type=='LENNARD-JONES') THEN
     ener_pot = 4.0_dbl*pot%lj%epsilon*(pot%lj%sigma12*r**(-12)-pot%lj% &
          sigma6*r**(-6)) - pot%energy_cutoff
  ELSE IF (pot%type=='WILLIAMS') THEN
     ener_pot = pot%willis%a*exp(-pot%willis%b*r) - pot%willis%c/r**6 &
          - pot%energy_cutoff
  ELSE IF (pot%type=='GOODWIN') THEN
     scale = exp(pot%goodwin%m*(-(r/pot%goodwin%dc)**pot%goodwin%mc+(pot% &
          goodwin%d/pot%goodwin%dc)**pot%goodwin%mc))
     ener_pot = scale*pot%goodwin%vr0*(pot%goodwin%d/r)**pot%goodwin%m
  ELSE
     CALL stop_program("pair_potential%ener_pot","Potential type unknown")
  END IF

END FUNCTION ener_pot

!******************************************************************************

FUNCTION ener_coul ( pot, qi, qj, r, ewald_param )

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), INTENT ( IN ) :: qi, qj, r
  TYPE (potentialparm_type ), INTENT ( IN ) :: pot
  TYPE (ewald_parameters_type ), INTENT ( IN ) :: ewald_param

! Locals
  EXTERNAL erfc
  REAL ( dbl ) :: ener_coul, erfc
  REAL ( dbl ), PARAMETER :: ifourpi = 1.0_dbl / ( 4.0_dbl * pi )
  REAL ( dbl ), PARAMETER:: eps = 1.0E-6_dbl

!------------------------------------------------------------------------------

  SELECT CASE( ewald_param % ewald_type )
  CASE( 'EWALD', 'PME', 'SPME' )
     ener_coul = ifourpi / ewald_param % eps0 * qi * qj &
          * erfc ( ewald_param % alpha * r ) / r - pot % e_cutoff_coul
  CASE( 'NONE')
     ener_coul= 0._dbl
  END SELECT

END FUNCTION ener_coul

!******************************************************************************

! creates the splines for the potentials

SUBROUTINE spline_nonbond_control ( potparm, pstat, npoints, ewald_param )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: npoints
  TYPE (potentialparm_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: potparm
  TYPE (particle_prop_type ), INTENT ( IN ) :: pstat ( : )
  TYPE (ewald_parameters_type ), INTENT ( IN ) :: ewald_param

! Locals
  LOGICAL, SAVE :: first_visit = .TRUE.
  INTEGER :: ntype, ntab, n, i, j, jx, isos
  REAL ( dbl ) :: qi, qj, e, ee, x, x2, dx2, locut, hicut, a, cutmax, cutmin

!------------------------------------------------------------------------------

  ntype = SIZE ( potparm, 1 )
  ntab = ( ntype * ntype + ntype ) / 2

  IF ( first_visit ) THEN
     NULLIFY ( spline_table )
     first_visit = .FALSE.
  END IF

  IF ( ASSOCIATED ( spline_table ) ) THEN
     DO i = 1, ntab
        CALL kill_spline ( spline_table(i), 'a' )
     END DO
     DEALLOCATE ( spline_table, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'spline_nonbond_control', &
          'spline_table' )
     NULLIFY ( spline_table )
  END IF

  ALLOCATE ( spline_table ( ntab ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'spline_nonbond_control', &
       'spline_table', ntab )
  DO i = 1, ntab
     NULLIFY ( spline_table(i) % x )
     NULLIFY ( spline_table(i) % y )
     NULLIFY ( spline_table(i) % y2 )
  END DO
  IF ( ALLOCATED ( spltab ) ) THEN
     DEALLOCATE ( spltab, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'spline_nonbond_control', 'spltab' )
  END IF
  ALLOCATE ( spltab ( ntype, ntype ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'spline_nonbond_control', &
       'spltab', ntype ** 2 )

  cutmin = 1.0E+30_dbl
  cutmax = 0.0_dbl
  n = 0
  DO i = 1, ntype ! i:  first atom type
     qi = pstat(i) %charge
     DO j = i, ntype ! j:  second atom type
        qj = pstat(j) %charge
        n = n + 1
        spltab(i,j) = n
        spltab(j,i) = n
!ODER: n=j+(i-1)*n-(i-1)*i/2 mit j>=i
        CALL init_splinexy(spline_table(n),npoints+1)

        IF (potparm(i,j)%type=='LENNARD-JONES') THEN
           locut =  0.5_dbl
        ELSE IF (potparm(i,j)%type=='WILLIAMS') THEN
           locut = 0.5 ! else put locut to 0.5 Angstrom
        ELSE
           locut = 0.5 ! else put locut to 0.5 Angstrom
        END IF
        hicut = sqrt(potparm(i,j)%rcutsq)
        IF ( hicut > cutmax ) cutmax = hicut
        IF ( locut < cutmin ) cutmin = locut

        potparm(i,j) %energy_cutoff = ener_pot(potparm(i,j),qi,qj,hicut)
        potparm(i,j) %e_cutoff_coul &
             = ener_coul ( potparm ( i, j ), qi, qj, hicut, ewald_param )

        dx2 = (hicut**2-locut**2)/real(npoints)
        x2 = locut**2
        DO jx = 1, npoints + 1 ! jx: loop over distance**2
           x = SQRT ( x2 )
           e = ener_pot ( potparm ( i, j ), qi, qj, x ) &
                + ener_coul ( potparm ( i, j ), qi, qj, x, ewald_param )
           spline_table(n) %x(jx) = x2
           spline_table(n) %y(jx) = e
           x2 = x2 + dx2
        END DO ! jx: loop over distance

        CALL init_spline(spline_table(n),dx=dx2)

     END DO ! j:  second atom type
  END DO ! i:  first atom type

END SUBROUTINE spline_nonbond_control

!******************************************************************************

!     calculates energy=V(r^2) and fscalar=-V^{prime}(r^2)/r
!     qi,qj are not used any more. only kept to have same args as potential()

SUBROUTINE potential_s ( rijsq, potparm, iatom, jatom, energy, fscalar )

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), INTENT ( IN ) :: rijsq
  REAL ( dbl ), INTENT ( OUT ) :: energy, fscalar
  INTEGER, INTENT ( IN ) :: iatom, jatom
  TYPE (potentialparm_type ), DIMENSION ( :, : ), INTENT ( IN ) :: potparm

!------------------------------------------------------------------------------

  energy = splineh ( spline_table ( spltab ( iatom, jatom ) ), rijsq, fscalar )
  fscalar = -2.0_dbl * fscalar

END SUBROUTINE potential_s

!******************************************************************************

SUBROUTINE potential_f ( rijsq, potparm, qi, qj, iatom, jatom, energy, &
     ewald_param )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: iatom, jatom
  REAL ( dbl ), INTENT ( IN ) :: rijsq
  REAL ( dbl ), INTENT ( OUT ) :: energy
  REAL ( dbl ), INTENT ( IN ) :: qi, qj
  TYPE (potentialparm_type ), INTENT ( IN ), DIMENSION ( :, : ) :: potparm
  TYPE (ewald_parameters_type ), INTENT ( IN ), OPTIONAL :: ewald_param

! Locals
  REAL ( dbl ) :: rij

!------------------------------------------------------------------------------

  rij = SQRT ( rijsq )

! non-bond energy
  energy = ener_pot ( potparm ( iatom, jatom ), qi, qj, rij )

! real-space ewald non-bond energy
  IF ( PRESENT ( ewald_param ) ) THEN
     energy = energy &
          + ener_coul ( potparm ( iatom, jatom ), qi, qj, rij, ewald_param )
  END IF

END SUBROUTINE potential_f

!******************************************************************************

END MODULE pair_potential

