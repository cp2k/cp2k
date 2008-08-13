!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      September 2005 - Introduced the Born-Mayer-Huggins-Fumi-Tosi  Potential (BMHTF)
!>      2006 - Major rewriting of the routines.. Linear scaling setup of splines
!> \author CJM
! *****************************************************************************
MODULE pair_potential_util

  USE erf_fn,                          ONLY: erfc
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE f77_blas
  USE fparser,                         ONLY: EvalErrType,&
                                             evalf,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE input_constants,                 ONLY: do_ewald_none
  USE kinds,                           ONLY: dp
  USE pair_potential_types,            ONLY: &
       b4_type, bm_type, ea_type, ft_type, gp_type, gw_type, ip_type, &
       lj_charmm_type, lj_type, not_initialized, pair_potential_pp_type, &
       pair_potential_single_type, wl_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pair_potential_util'
  REAL(KIND=dp), PARAMETER, PRIVATE    :: MIN_HICUT_VALUE=1.0E-15_dp,&
                                          DEFAULT_HICUT_VALUE=1.0E3_dp

  PUBLIC :: potential_f
  PUBLIC :: ener_coul, ener_pot

CONTAINS

! *****************************************************************************
!> \brief Evaluates the nonbond potential energy for the implemented FF kinds
! *****************************************************************************
  FUNCTION ener_pot ( pot, r, energy_cutoff, vdw_scale14 ) RESULT(value)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: r, energy_cutoff, vdw_scale14
    REAL(KIND=dp)                            :: value

    INTEGER                                  :: i, index, j, n
    REAL(KIND=dp)                            :: lvalue, pp, qq, scale

    value = 0.0_dp
    DO j = 1, SIZE(pot%type)
       ! A lower boundary for the potential definition was defined
       IF ((pot%set(j)%rmin /= not_initialized).AND.( r <  pot%set(j)%rmin )) CYCLE
       ! An upper boundary for the potential definition was defined
       IF ((pot%set(j)%rmax /= not_initialized).AND.( r >= pot%set(j)%rmax )) CYCLE
       ! If within limits let's compute the potential...
       IF (pot%type(j)==lj_charmm_type) THEN
          lvalue = vdw_scale14 * &
               4.0_dp*pot%set(j)%lj%epsilon*(pot%set(j)%lj%sigma12*r**(-12)-pot%set(j)%lj% &
               sigma6*r**(-6))
       ELSE IF (pot%type(j)==lj_type) THEN
          lvalue = vdw_scale14 * pot%set(j)%lj%epsilon* &
               (pot%set(j)%lj%sigma12*r**(-12)-pot%set(j)%lj%sigma6*r**(-6))
       ELSE IF (pot%type(j)==ip_type) THEN
          lvalue = 0._dp
          IF ( r > pot%set(j)%ipbv%rcore ) THEN
             DO i = 2, 15
                lvalue = lvalue + pot%set(j)%ipbv%a(i)/(r**(i-1)*REAL(i-1,dp))
             END DO
          ELSE
             ! use a linear potential
             lvalue = pot%set(j)%ipbv%m * r + pot%set(j)%ipbv%b
          ENDIF
          lvalue=lvalue
       ELSE IF (pot%type(j)==wl_type) THEN
          lvalue = pot%set(j)%willis%a*EXP(-pot%set(j)%willis%b*r) - pot%set(j)%willis%c/r**6 
       ELSE IF (pot%type(j)==gw_type) THEN
          scale = EXP(pot%set(j)%goodwin%m*(-(r/pot%set(j)%goodwin%dc)**pot%set(j)%goodwin%mc+&
               (pot%set(j)%goodwin%d/pot%set(j)%goodwin%dc)**pot%set(j)%goodwin%mc))
          lvalue = scale*pot%set(j)%goodwin%vr0*(pot%set(j)%goodwin%d/r)**pot%set(j)%goodwin%m
       ELSE IF (pot%type(j)==ft_type) THEN
          lvalue = pot%set(j)%ft%a*EXP(-pot%set(j)%ft%b*r) - pot%set(j)%ft%c/r**6 - pot%set(j)%ft%d/r**8
       ELSE IF (pot%type(j)==ea_type) THEN
          index = INT ( r / pot % set(j) % eam % drar ) + 1
          IF ( index > 5000 ) THEN
             index = 5000
          ELSEIF ( index < 1 ) THEN
             index = 1
          ENDIF
          qq = r - pot % set(j) % eam % rval (index)
          pp =  pot % set(j) % eam % phi ( index ) + &
               qq * pot % set(j) % eam % phip ( index )
          lvalue = pp
       ELSE IF (pot%type(j)==b4_type) THEN
         IF(r<=pot%set(j)%buck4r%r1) THEN
           pp = pot%set(j)%buck4r%a*EXP(-pot%set(j)%buck4r%b*r)
         ELSEIF(r>pot%set(j)%buck4r%r1.AND.r<=pot%set(j)%buck4r%r2) THEN
           pp = 0.0_dp
           DO n = 0,pot%set(j)%buck4r%npoly1
             pp = pp+pot%set(j)%buck4r%poly1(n)*r**n
           END DO 
         ELSEIF(r>pot%set(j)%buck4r%r2.AND.r<=pot%set(j)%buck4r%r3) THEN
           pp = 0.0_dp
           DO n = 0,pot%set(j)%buck4r%npoly2
             pp = pp+pot%set(j)%buck4r%poly2(n)*r**n
           END DO 
         ELSEIF(r>pot%set(j)%buck4r%r3)THEN
           pp = -pot%set(j)%buck4r%c/r**6
         END IF
         lvalue = pp
       ELSE IF (pot%type(j)==bm_type) THEN
         lvalue = pot%set(j)%buckmo%f0*(pot%set(j)%buckmo%b1+pot%set(j)%buckmo%b2)*&
                  EXP((pot%set(j)%buckmo%a1+pot%set(j)%buckmo%a2-r)/(pot%set(j)%buckmo%b1+pot%set(j)%buckmo%b2))&
                  - pot%set(j)%buckmo%c/r**6  &
                  + pot%set(j)%buckmo%d*(EXP(-2._dp*pot%set(j)%buckmo%beta*(r-pot%set(j)%buckmo%r0))- &
                                  2.0_dp*EXP(-pot%set(j)%buckmo%beta*(r-pot%set(j)%buckmo%r0)))
       ELSE IF (pot%type(j)==gp_type) THEN
          pot%set(j)%gp%values(1)=r
          lvalue = evalf(pot%set(j)%gp%myid,pot%set(j)%gp%values)
          IF (EvalErrType > 0) &
               CALL stop_program("ener_pot","Error evaluating generic potential energy function")
       ELSE
          lvalue = 0.0_dp
       END IF
       value = value + lvalue
    END DO
    value = value - energy_cutoff
  END FUNCTION ener_pot

! *****************************************************************************
!> \brief Evaluates the Coulomb electrostatic energy
! *****************************************************************************
  FUNCTION ener_coul ( e_cutoff_coul, fac, r, ewald_type, alpha, my_do_14 )

    REAL(KIND=dp), INTENT(IN)                :: e_cutoff_coul, fac, r
    INTEGER, INTENT(IN)                      :: ewald_type
    REAL(KIND=dp), INTENT(IN)                :: alpha
    LOGICAL, INTENT(IN)                      :: my_do_14
    REAL(KIND=dp)                            :: ener_coul

    ener_coul= 0.0_dp
    IF (.NOT. my_do_14 ) THEN
       IF (ewald_type /= do_ewald_none) THEN
          ener_coul = fac * erfc ( alpha * r ) / r -  e_cutoff_coul
       ELSE
          ener_coul = fac / r - e_cutoff_coul
       END IF
    ELSE
       ! to avoid a large number of points in the spline generation
       ! we evaluate the 1-4 interactions analytically..
       ener_coul = 0.0_dp
    ENDIF

  END FUNCTION ener_coul

! *****************************************************************************
!> \brief Evaluates the nonbond and electrostaitic potential analytically
!> \author Teodoro Laino [tlaino] 2006.05
! *****************************************************************************
  SUBROUTINE potential_f ( rijsq, potparm, qi, qj, iatom, jatom, energy, ewald_env, error )
    REAL(KIND=dp), INTENT(IN)                :: rijsq
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    REAL(KIND=dp), INTENT(IN)                :: qi, qj
    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(KIND=dp), INTENT(OUT)               :: energy
    TYPE(ewald_environment_type), OPTIONAL, &
      POINTER                                :: ewald_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ewald_type, i, max_id
    LOGICAL                                  :: do_multipoles
    REAL(KIND=dp)                            :: alpha, e_cutoff_coul, &
                                                energy_cutoff, f, hicut, rij

    rij = SQRT (rijsq)
    hicut = SQRT (potparm%pot(iatom, jatom)%pot%rcutsq )
    IF (ABS(hicut) <= MIN_HICUT_VALUE) hicut = DEFAULT_HICUT_VALUE
    ! If genpot we need to get ready for the evaluation of the potential
    IF (ANY(potparm%pot(iatom,jatom)%pot%type==gp_type)) THEN
       max_id = -HUGE(0)
       DO i = 1, SIZE(potparm%pot(iatom,jatom)%pot%type)
          IF (potparm%pot(iatom,jatom)%pot%type(i)/=gp_type) CYCLE
          max_id = MAX(max_id, potparm%pot(iatom,jatom)%pot%set(i)%gp%myid)
       END DO
       CALL initf(max_id)
       DO i = 1, SIZE(potparm%pot(iatom,jatom)%pot%type)
          IF (potparm%pot(iatom,jatom)%pot%type(i)/=gp_type) CYCLE
          CALL parsef(potparm%pot(iatom,jatom)%pot%set(i)%gp%myid,&
                      TRIM(potparm%pot(iatom,jatom)%pot%set(i)%gp%potential),&
                      potparm%pot(iatom,jatom)%pot%set(i)%gp%parameters)
       END DO
    END IF
    energy_cutoff = ener_pot ( potparm%pot(iatom,jatom)%pot, hicut, 0.0_dp, 1.0_dp)

    ! non-bond energy
    energy = ener_pot ( potparm%pot(iatom,jatom)%pot, rij, energy_cutoff, 1.0_dp)
    IF (ANY(potparm%pot(iatom,jatom)%pot%type==gp_type)) THEN
       CALL finalizef()
    ENDIF
    ! real-space ewald non-bond energy
    IF ( PRESENT ( ewald_env ) ) THEN
       CALL ewald_env_get(ewald_env,ewald_type=ewald_type,alpha=alpha,&
            do_multipoles=do_multipoles,error=error)
       IF (.NOT.do_multipoles) THEN
          f = qi * qj
          e_cutoff_coul = ener_coul ( 0.0_dp, f, hicut, ewald_type, alpha, my_do_14=.FALSE.)
          energy = energy + ener_coul (e_cutoff_coul, f, rij, ewald_type, alpha, my_do_14=.FALSE.)
       END IF
    END IF

  END SUBROUTINE potential_f

END MODULE pair_potential_util

