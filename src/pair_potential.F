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
                                             splineh

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: potential_f, potential_s
  PUBLIC :: spline_nonbond_control

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

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
       ener_pot=ener_pot - pot % energy_cutoff
     ELSE
! use a linear potential  
       ener_pot = pot % ipbv % m * r + pot % ipbv % b
     ENDIF
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
    REAL(KIND=dp)                            :: alpha, eps0

!------------------------------------------------------------------------------

  IF ( PRESENT ( ewald_env ) .AND. .NOT.PRESENT(ei_scale14) ) THEN
    CALL ewald_env_get ( ewald_env, alpha = alpha, eps0 = eps0, &
                       ewald_type = ewald_type ) 
    ener_coul= 0.0_dp
    IF (ewald_type /= do_ewald_none) &
         ener_coul = ifourpi / eps0 * qi * qj &
         * erfc ( alpha * r ) / r - pot % e_cutoff_coul
  ELSE
    CALL ewald_env_get ( ewald_env, alpha = alpha, eps0 = eps0, &
                       ewald_type = ewald_type ) 
         ener_coul = ei_scale14 * ifourpi / eps0 * qi * qj / r - pot % e_cutoff_coul
  ENDIF

END FUNCTION ener_coul

!******************************************************************************

! creates the splines for the potentials

SUBROUTINE spline_nonbond_control ( spline_env, potparm, atomic_kind_set, &
                                    ewald_env, ei_scale14, vdw_scale14 )


    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(ewald_environment_type), OPTIONAL, &
      POINTER                                :: ewald_env
    REAL(KIND=dp), OPTIONAL                  :: ei_scale14, vdw_scale14

    INTEGER                                  :: i, j, jx, n, npoints, ntype
    REAL(KIND=dp)                            :: cutmax, cutmin, dx2, e, &
                                                hicut, locut, qi, qj, x, x2
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!------------------------------------------------------------------------------

  npoints = spline_env % nspline_points

  cutmin = 1.0E+30_dp
  cutmax = 0.0_dp
  n = 0

  ntype = SIZE ( atomic_kind_set )
  DO i = 1, ntype ! i:  first atom type

     atomic_kind=> atomic_kind_set(i)
     CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qi)

     DO j = i, ntype ! j:  second atom type

        atomic_kind=> atomic_kind_set(j)
        CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qj)
        n = n + 1
        spline_env % spltab ( i, j ) = n
        spline_env % spltab ( j, i ) = n
!ODER: n=j+(i-1)*n-(i-1)*i/2 mit j>=i
        CALL init_splinexy ( spline_env % spline_data_p ( n ) % spline_data, npoints+1)

        IF (potparm%pot(i,j)%pot%type=='LENNARD-JONES') THEN
           locut =  0.944863_dp ! else put locut to 0.5 Angstrom
        ELSE IF (potparm%pot(i,j)%pot%type=='WILLIAMS') THEN
           locut = 0.944863_dp  ! else put locut to 0.5 Angstrom
        ELSE IF (potparm%pot(i,j)%pot%type=='IPBV') THEN
           locut = 0.944863_dp  ! 0.5 Angstrom
        ELSE
           locut = 0.944863_dp  ! 0.5 Angstrom
        END IF
        hicut = SQRT ( potparm %pot ( i, j ) %pot % rcutsq )
        IF ( hicut > cutmax ) cutmax = hicut
        IF ( locut < cutmin ) cutmin = locut
        potparm%pot(i,j)%pot %energy_cutoff = 0.0_dp
        IF ( PRESENT ( vdw_scale14 ) ) THEN
          potparm%pot(i,j)%pot %energy_cutoff = ener_pot(potparm%pot(i,j)%pot,qi,qj,hicut,vdw_scale14)
        ELSE
          potparm%pot(i,j)%pot %energy_cutoff = ener_pot(potparm%pot(i,j)%pot,qi,qj,hicut)
        END IF
        potparm%pot(i,j)%pot %e_cutoff_coul = 0.0_dp

        IF ( PRESENT ( ei_scale14 ) ) THEN
          potparm%pot(i,j)%pot %e_cutoff_coul &
          = ener_coul ( potparm%pot( i, j )%pot, qi, qj, hicut, ewald_env, ei_scale14 )
        ELSE
          potparm%pot(i,j)%pot %e_cutoff_coul &
          = ener_coul ( potparm%pot ( i, j )%pot, qi, qj, hicut, ewald_env )
        ENDIF
        potparm%pot(j,i)%pot %energy_cutoff = potparm%pot(i,j)%pot %energy_cutoff
        potparm%pot(j,i)%pot %e_cutoff_coul = potparm%pot(i,j)%pot %e_cutoff_coul

        dx2 = (hicut**2-locut**2)/REAL(npoints)
        x2 = locut**2
        DO jx = 1, npoints + 1 ! jx: loop over distance**2
           x = SQRT ( x2 )
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

        CALL init_spline(spline_env % spline_data_p ( n ) % spline_data,dx=dx2)

     END DO ! j:  second atom type
  END DO ! i:  first atom type

END SUBROUTINE spline_nonbond_control

!******************************************************************************

!     calculates energy=V(r^2) and fscalar=-V^{prime}(r^2)/r
!     qi,qj are not used any more. only kept to have same args as potential()

SUBROUTINE potential_s ( spline_table, rijsq, energy, fscalar )


    TYPE(spline_data_type), INTENT(IN)       :: spline_table
    REAL(KIND=dp), INTENT(IN)                :: rijsq
    REAL(KIND=dp), INTENT(OUT)               :: energy, fscalar

!------------------------------------------------------------------------------

  energy = splineh ( spline_table, rijsq, fscalar )
  fscalar = -2.0_dp * fscalar

END SUBROUTINE potential_s

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

