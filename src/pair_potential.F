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

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE erf_fn,                          ONLY: erfc
  USE ewald_environment_types,         ONLY: ewald_environment_type, &
                                             ewald_env_get
  USE kinds,                           ONLY: dbl, default_string_length
  USE mathconstants,                   ONLY: pi
  USE pair_potential_types,            ONLY: pair_potential_type
  USE splines,                         ONLY: init_spline,&
                                             init_splinexy,&
                                             kill_spline,&
                                             spline_data,&
                                             splineh, &
                                             spline_environment_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: potential_f, potential_s
  PUBLIC :: spline_nonbond_control

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

FUNCTION ener_pot ( pot, qi, qj, r )


    TYPE(pair_potential_type), INTENT(IN)    :: pot
    REAL(dbl), INTENT(IN)                    :: qi, qj, r
    REAL(dbl)                                :: ener_pot

    REAL(dbl)                                :: scale

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
     WRITE(*,*) "pot%type = ",pot%type
     CALL stop_program("pair_potential%ener_pot","Potential type unknown")
  END IF

END FUNCTION ener_pot

!******************************************************************************

FUNCTION ener_coul ( pot, qi, qj, r, ewald_env )


    TYPE(pair_potential_type), INTENT(IN)          :: pot
    REAL(dbl), INTENT(IN)                          :: qi, qj, r
    TYPE(ewald_environment_type), POINTER, OPTIONAL    :: ewald_env
    REAL(dbl)                                :: ener_coul

    REAL(dbl), PARAMETER :: eps = 1.0E-6_dbl, &
      ifourpi = 1.0_dbl / ( 4.0_dbl * pi )
    CHARACTER ( len = default_string_length ) :: ewald_type
    REAL ( dbl ) :: eps0, alpha

!------------------------------------------------------------------------------

  IF ( PRESENT ( ewald_env ) ) THEN
    CALL ewald_env_get ( ewald_env, alpha = alpha, eps0 = eps0, &
                       ewald_type = ewald_type ) 
    SELECT CASE(  ewald_type )
    CASE( 'EWALD', 'PME', 'SPME' )
       ener_coul = ifourpi / eps0 * qi * qj &
            * erfc ( alpha * r ) / r - pot % e_cutoff_coul
    CASE( 'NONE')
       ener_coul= 0._dbl
    END SELECT
  ELSE
    ener_coul= 0._dbl
  ENDIF

END FUNCTION ener_coul

!******************************************************************************

! creates the splines for the potentials

SUBROUTINE spline_nonbond_control ( spline_env, atomic_kind_set, ewald_env )


    TYPE ( spline_environment_type ), INTENT ( INOUT )  :: spline_env 
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                           :: atomic_kind_set
    TYPE(ewald_environment_type), POINTER, OPTIONAL     :: ewald_env

    INTEGER                                  :: i, isos, j, jx, n, ntab, ntype
    INTEGER                                  :: npoints
    REAL(dbl)                                :: a, cutmax, cutmin, dx2, e, &
                                                hicut, locut, qi, qj, x, x2
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
  !  INTEGER, POINTER                         :: spltab ( :, : )
  !  TYPE ( spline_data ), POINTER            :: spline_table ( : )
    TYPE ( pair_potential_type ), POINTER    :: potparm ( :, : )

!------------------------------------------------------------------------------

  potparm => spline_env % potparm
  npoints = spline_env % nspline_points

  ntype = SIZE ( atomic_kind_set )
  ntab = ( ntype * ntype + ntype ) / 2

  IF ( ASSOCIATED ( spline_env % spline_table ) ) THEN
     DO i = 1, ntab
        CALL kill_spline ( spline_env % spline_table(i), 'a' )
     END DO
     DEALLOCATE ( spline_env % spline_table, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'spline_nonbond_control', &
          'spline_table' )
     NULLIFY ( spline_env % spline_table )
  END IF

  ALLOCATE ( spline_env % spline_table ( ntab ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'spline_nonbond_control', &
       'spline_table', ntab )
  DO i = 1, ntab
     NULLIFY ( spline_env % spline_table(i) % x )
     NULLIFY ( spline_env % spline_table(i) % y )
     NULLIFY ( spline_env % spline_table(i) % y2 )
  END DO
  IF ( ASSOCIATED ( spline_env % spltab ) ) THEN
     DEALLOCATE ( spline_env % spltab, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'spline_nonbond_control', 'spltab' )
  END IF
  ALLOCATE ( spline_env % spltab ( ntype, ntype ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'spline_nonbond_control', &
       'spltab', ntype ** 2 )

  cutmin = 1.0E+30_dbl
  cutmax = 0.0_dbl
  n = 0

  DO i = 1, ntype ! i:  first atom type

     atomic_kind=> atomic_kind_set(i)
     CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qi)

     DO j = i, ntype ! j:  second atom type

        atomic_kind=> atomic_kind_set(j)
        CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qj)
        n = n + 1
        spline_env % spltab(i,j) = n
        spline_env % spltab(j,i) = n
!ODER: n=j+(i-1)*n-(i-1)*i/2 mit j>=i
        CALL init_splinexy(spline_env % spline_table(n),npoints+1)

        IF (potparm(i,j)%type=='LENNARD-JONES') THEN
           locut =  0.944863_dbl ! else put locut to 0.5 Angstrom
        ELSE IF (potparm(i,j)%type=='WILLIAMS') THEN
           locut = 0.944863_dbl  ! else put locut to 0.5 Angstrom
        ELSE
           locut = 0.944863_dbl  ! else put locut to 0.5 Angstrom
        END IF
        hicut = SQRT(potparm(i,j)%rcutsq)
        IF ( hicut > cutmax ) cutmax = hicut
        IF ( locut < cutmin ) cutmin = locut
        potparm(i,j) %energy_cutoff = 0.0_dbl
        potparm(i,j) %energy_cutoff = ener_pot(potparm(i,j),qi,qj,hicut)

        potparm(i,j) %e_cutoff_coul = 0.0_dbl
        potparm(i,j) %e_cutoff_coul &
              = ener_coul ( potparm ( i, j ), qi, qj, hicut, ewald_env )

        potparm(j,i) %energy_cutoff = potparm(i,j) %energy_cutoff
        potparm(j,i) %e_cutoff_coul = potparm(i,j) %e_cutoff_coul

        dx2 = (hicut**2-locut**2)/real(npoints)
        x2 = locut**2
        DO jx = 1, npoints + 1 ! jx: loop over distance**2
           x = SQRT ( x2 )
           e = ener_pot ( potparm ( i, j ), qi, qj, x ) &
                + ener_coul ( potparm ( i, j ), qi, qj, x, ewald_env )
           spline_env % spline_table(n) %x(jx) = x2
           spline_env % spline_table(n) %y(jx) = e
           x2 = x2 + dx2
        END DO ! jx: loop over distance

        CALL init_spline(spline_env % spline_table(n),dx=dx2)

     END DO ! j:  second atom type
  END DO ! i:  first atom type

END SUBROUTINE spline_nonbond_control

!******************************************************************************

!     calculates energy=V(r^2) and fscalar=-V^{prime}(r^2)/r
!     qi,qj are not used any more. only kept to have same args as potential()

SUBROUTINE potential_s ( spline_table, spltab, rijsq, iatom, jatom, energy, fscalar )


    REAL(dbl), INTENT(IN)                    :: rijsq
    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(dbl), INTENT(OUT)                   :: energy, fscalar
    INTEGER, INTENT ( IN ), DIMENSION ( :, : ) :: spltab
    TYPE ( spline_data ), INTENT ( IN ), DIMENSION ( : ) :: spline_table

!------------------------------------------------------------------------------

  energy = splineh ( spline_table ( spltab ( iatom, jatom ) ), rijsq, fscalar )
  fscalar = -2.0_dbl * fscalar

END SUBROUTINE potential_s

!******************************************************************************

SUBROUTINE potential_f ( rijsq, potparm, qi, qj, iatom, jatom, energy, &
     ewald_env )


    REAL(dbl), INTENT(IN)                    :: rijsq
    TYPE(pair_potential_type), &
      DIMENSION(:, :), INTENT(IN)            :: potparm
    REAL(dbl), INTENT(IN)                    :: qi, qj
    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(dbl), INTENT(OUT)                   :: energy
    TYPE(ewald_environment_type), &
      POINTER, OPTIONAL                      :: ewald_env

    REAL(dbl)                                :: rij

!------------------------------------------------------------------------------

  rij = SQRT ( rijsq )

! non-bond energy
  energy = ener_pot ( potparm ( iatom, jatom ), qi, qj, rij )

! real-space ewald non-bond energy
  IF ( PRESENT ( ewald_env ) ) THEN
     energy = energy &
          + ener_coul ( potparm ( iatom, jatom ), qi, qj, rij, ewald_env )
  END IF

END SUBROUTINE potential_f

!******************************************************************************

END MODULE pair_potential

