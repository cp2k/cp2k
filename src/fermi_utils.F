!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief deal with the Fermi distribution, compute it, fix mu, get derivs
!> \author Joost VandeVondele
!> \date 09.2008
! *****************************************************************************
MODULE fermi_utils
  
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Fermi, FermiFixed, FermiFixedDeriv, Fermikp

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fermi_utils'

CONTAINS
! *****************************************************************************
!> \brief   returns occupations according to Fermi-Dirac statistics
!>          for a given set of energies and fermi level.
!>          Note that singly occupied orbitals are assumed
!> \param   f occupations
!> \param   N total number of electrons (output)
!> \param kTS ...
!> \param   e eigenvalues
!> \param   mu Fermi level (input)
!> \param   T  electronic temperature
!> \param   maxocc maximum occupation of an orbital
!> \param   estate excited state in core level spectroscopy
!> \param   festate occupation of the excited state in core level spectroscopy
!> \date    09.2008
!> \par History
!>          - Made estate and festate optional (LT, 2014/02/26)
!> \author  Joost VandeVondele
! *****************************************************************************
 SUBROUTINE Fermi(f,N,kTS,e,mu,T,maxocc,estate,festate)

    REAL(KIND=dp), INTENT(out)               :: f(:), N, kTS
    REAL(KIND=dp), INTENT(IN)                :: e(:), mu, T, maxocc
    INTEGER, INTENT(IN), OPTIONAL            :: estate
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: festate

    CHARACTER(len=*), PARAMETER :: routineN = 'Fermi', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: I, Nstate
    REAL(KIND=dp)                            :: arg, occupation, term1, &
                                                term2, tmp, tmp2, tmp3, tmp4, &
                                                tmplog

   Nstate=SIZE(e)
   kTS = 0.0_dp
   ! kTS is the entropic contribution to the energy i.e. -TS
   ! kTS= kT*[f ln f + (1-f) ln (1-f)]

   DO i = 1, Nstate
      IF (PRESENT(estate) .AND. PRESENT(festate)) THEN
         IF (i == estate) THEN
            occupation = festate
         ELSE
            occupation = maxocc
         END IF
      ELSE
         occupation = maxocc
      END IF
      ! have the result of exp go to zero instead of overflowing
      IF (e(i) > mu) THEN
         arg  = -(e(i)-mu)/T
         ! tmp is smaller than 1
         tmp  = EXP(arg)
         tmp4 = tmp + 1.0_dp
         tmp2 = tmp/tmp4
         tmp3 = 1.0_dp/tmp4
         f(i) = occupation*tmp2
         ! log(1+eps), might need to be written more accurately
         tmplog = -LOG(tmp4)
         term1= tmp2 * ( arg + tmplog )
         term2= tmp3 * tmplog
         kTS = kTS + T*occupation*(term1 + term2)
      ELSE
         arg  = (e(i)-mu)/T
         ! tmp is smaller than 1
         tmp  = EXP(arg)
         tmp4 = tmp + 1.0_dp
         tmp2 = 1.0_dp/tmp4
         tmp3 = tmp/tmp4
         f(i) = occupation*tmp2
         tmplog = -LOG(tmp4)
         term1= tmp2 * tmplog
         term2= tmp3 * (arg + tmplog )
         kTS = kTS + T*occupation*(term1 + term2)
      END IF
   END DO

   N = accurate_sum(f)

 END SUBROUTINE Fermi

! *****************************************************************************
!> \brief Fermi function for KP cases
!> \param f   Occupation numbers
!> \param nel Number of electrons (total)
!> \param kTS Entropic energy contribution
!> \param e   orbital (band) energies
!> \param mu  chemical potential
!> \param wk  kpoint weights
!> \param t   Temperature
!> \param maxocc Maximum occupation
! *****************************************************************************
 SUBROUTINE Fermi2(f,nel,kTS,e,mu,wk,t,maxocc)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: f
    REAL(KIND=dp), INTENT(OUT)               :: nel, kTS
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: e
    REAL(KIND=dp), INTENT(IN)                :: mu
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: wk
    REAL(KIND=dp), INTENT(IN)                :: t, maxocc

    CHARACTER(len=*), PARAMETER :: routineN = 'Fermi2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, is, nkp, nmo
    REAL(KIND=dp)                            :: arg, beta, term1, term2, tmp, &
                                                tmp2, tmp3, tmp4, tmplog

   nmo = SIZE(e,1)
   nkp = SIZE(e,2)
   kTS = 0.0_dp
   ! kTS is the entropic contribution to the energy i.e. -TS
   ! kTS= kT*[f ln f + (1-f) ln (1-f)]
   IF(t > 1.0e-14_dp) THEN
      beta = 1.0_dp/t
      DO ik=1,nkp
         DO is=1,nmo
            IF (e(is,ik) > mu) THEN
               arg  = -(e(is,ik)-mu)*beta
               tmp  = EXP(arg)
               tmp4 = tmp + 1.0_dp
               tmp2 = tmp/tmp4
               tmp3 = 1.0_dp/tmp4
               f(is,ik) = maxocc*tmp2
               tmplog = -LOG(tmp4)
               term1 = tmp2 * ( arg + tmplog )
               term2 = tmp3 * tmplog
               kTS = kTS + t*maxocc*(term1 + term2)
            ELSE
               arg  = (e(is,ik)-mu)*beta
               tmp  = EXP(arg)
               tmp4 = tmp + 1.0_dp
               tmp2 = 1.0_dp/tmp4
               tmp3 = tmp/tmp4
               f(is,ik) = maxocc*tmp2
               tmplog = -LOG(tmp4)
               term1= tmp2 * tmplog
               term2= tmp3 * (arg + tmplog )
               kTS = kTS + T*maxocc*(term1 + term2)
            END IF
         END DO
      END DO
   ELSE
      beta = HUGE(1.0_dp)
      DO ik=1,nkp
         DO is=1,nmo
            IF (e(is,ik) <= mu) THEN
               f(is,ik) = maxocc
            ELSE
               f(is,ik) = 0.0_dp
            END IF
         END DO
      END DO
   END IF

   nel = 0.0_dp
   DO ik=1,nkp
      nel = nel + accurate_sum(f(1:nmo,ik)) * wk(ik)
   END DO

 END SUBROUTINE Fermi2

! *****************************************************************************
!> \brief   returns occupations according to Fermi-Dirac statistics
!>          for a given set of energies and number of electrons.
!>          Note that singly occupied orbitals are assumed.
!>          could fail if the fermi level lies out of the range of eigenvalues
!>          (to be fixed)
!> \param   f occupations
!> \param   mu Fermi level (output)
!> \param kTS ...
!> \param   e eigenvalues
!> \param   N total number of electrons (input)
!> \param   T  electronic temperature
!> \param   maxocc maximum occupation of an orbital
!> \param   estate excited state in core level spectroscopy
!> \param   festate occupation of the excited state in core level spectroscopy
!> \date    09.2008
!> \par History
!>          - Made estate and festate optional (LT, 2014/02/26)
!> \author  Joost VandeVondele
! *****************************************************************************
 SUBROUTINE FermiFixed(f,mu,kTS,e,N,T,maxocc,estate,festate)
    REAL(KIND=dp), INTENT(OUT)               :: f(:), mu, kTS
    REAL(KIND=dp), INTENT(IN)                :: e(:), N, T, maxocc
    INTEGER, INTENT(IN), OPTIONAL            :: estate
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: festate

    CHARACTER(len=*), PARAMETER :: routineN = 'FermiFixed', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iter, my_estate
    REAL(KIND=dp)                            :: mu_max, mu_min, mu_now, &
                                                my_festate, N_max, N_min, &
                                                N_now

    IF (PRESENT(estate) .AND. PRESENT(festate)) THEN
       my_estate = estate
       my_festate = festate
    ELSE
       my_estate = NINT(maxocc)
       my_festate = my_estate
    END IF

! bisection search to find N
! first bracket

    mu_min=MINVAL(e)
    iter=0
    DO
      iter=iter+1
      CALL Fermi(f, N_min, kTS, e, mu_min, T, maxocc, my_estate, my_festate)
      IF (N_min>N .OR. iter>20) THEN
         mu_min=mu_min-T
      ELSE
         EXIT
      ENDIF
    ENDDO

    mu_max=MAXVAL(e)
    iter=0
    DO
      iter=iter+1
      CALL Fermi(f, N_max, kTS, e, mu_max, T, maxocc, my_estate, my_festate)
      IF (N_max<N .OR. iter>20) THEN
         mu_max=mu_max+T
      ELSE
         EXIT
      ENDIF
    ENDDO

    ! now bisect
    iter=0
    DO WHILE( mu_max-mu_min > EPSILON(mu)*MAX(1.0_dp,ABS(mu_max),ABS(mu_min)) )
       iter=iter+1
       mu_now = (mu_max+mu_min)/2.0_dp
       CALL Fermi(f, N_now, kTS, e, mu_now, T, maxocc, my_estate, my_festate)
       iter = iter + 1
       IF (N_now<=N) THEN
          mu_min = mu_now
          N_min  = N_now
       ELSE
          mu_max = mu_now
          N_max  = N_now
       ENDIF
       IF (iter>200) EXIT
    ENDDO

    mu=(mu_max+mu_min)/2.0_dp
    CALL Fermi(f, N_now, kTS, e, mu, T, maxocc, my_estate, my_festate)

 END SUBROUTINE FermiFixed

! *****************************************************************************
!> \brief Bisection search to find mu for a given nel (kpoint case)
!> \param f   Occupation numbers
!> \param mu  chemical potential
!> \param kTS Entropic energy contribution
!> \param e   orbital (band) energies
!> \param nel Number of electrons (total)
!> \param wk  kpoint weights
!> \param t   Temperature
!> \param maxocc Maximum occupation
! *****************************************************************************
 SUBROUTINE Fermikp(f,mu,kTS,e,nel,wk,t,maxocc)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: f
    REAL(KIND=dp), INTENT(OUT)               :: mu, kTS
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: e
    REAL(KIND=dp), INTENT(IN)                :: nel
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: wk
    REAL(KIND=dp), INTENT(IN)                :: t, maxocc

    CHARACTER(len=*), PARAMETER :: routineN = 'Fermikp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iter
    REAL(KIND=dp)                            :: mu_max, mu_min, mu_now, &
                                                N_max, N_min, N_now

! bisection search to find mu for a given nel

    mu_min=MINVAL(e) - 0.5_dp
    mu_max=MAXVAL(e) + 0.5_dp
    iter=0
    DO WHILE( mu_max-mu_min > EPSILON(mu)*MAX(1.0_dp,ABS(mu_max),ABS(mu_min)) )
       iter=iter+1
       mu_now = (mu_max+mu_min)/2.0_dp
       CALL Fermi2(f,N_now,kTS,e,mu_now,wk,T,maxocc)
       iter = iter + 1
       IF (N_now<=nel) THEN
          mu_min = mu_now
          N_min  = N_now
       ELSE
          mu_max = mu_now
          N_max  = N_now
       ENDIF
       IF (iter>200) EXIT
    ENDDO

    mu=(mu_max+mu_min)/2.0_dp
    CALL Fermi2(f,N_now,kTS,e,mu,wk,T,maxocc)

 END SUBROUTINE Fermikp

! *****************************************************************************
!> \brief   returns f and dfde for a given set of energies and number of electrons
!>          it is a numerical derivative, trying to use a reasonable step length
!>          it ought to yield an accuracy of approximately EPSILON()^(2/3) (~10^-11)
!>          l ~ 10*T yields best accuracy
!>          Note that singly occupied orbitals are assumed.
!>          To be fixed: this could be parallellized for better efficiency
!> \param   dfde derivatives of the occupation numbers with respect to the eigenvalues
!>               the ith column is the derivative of f wrt to e_i
!> \param   f occupations
!> \param   mu Fermi level (input)
!> \param kTS ...
!> \param   e eigenvalues
!> \param   N total number of electrons (output)
!> \param   T  electronic temperature
!> \param maxocc ...
!> \param   l  typical lenght scale (~ 10 * T)
!> \param estate ...
!> \param festate ...
!> \date    09.2008
!> \par   History
!>          - Made estate and festate optional (LT, 2014/02/26)
!>          - Changed order of input, so l is before the two optional varaibles
!>            (LT, 2014/02/26)
!> \author  Joost VandeVondele
! *****************************************************************************
 SUBROUTINE FermiFixedDeriv(dfde,f,mu,kTS,e,N,T,maxocc,l,estate,festate)
    REAL(KIND=dp), INTENT(OUT)               :: dfde(:,:), f(:), mu, kTS
    REAL(KIND=dp), INTENT(IN)                :: e(:), N, T, maxocc, l
    INTEGER, INTENT(IN), OPTIONAL            :: estate
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: festate

    CHARACTER(len=*), PARAMETER :: routineN = 'FermiFixedDeriv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, I, my_estate, Nstate
    REAL(KIND=dp)                            :: h, mux, my_festate
    REAL(KIND=dp), ALLOCATABLE               :: ex(:), fx(:)

    CALL timeset(routineN, handle)

    Nstate=SIZE(e)
    ALLOCATE(ex(Nstate), fx(Nstate))

    IF (PRESENT(estate) .AND. PRESENT(festate)) THEN
       my_estate = estate
       my_festate = festate
    ELSE
       my_estate = NINT(maxocc)
       my_festate = my_estate
    END IF

    DO I=1,Nstate
       ! NR 5.7.8
       ! the problem here is that each f_i 'seems to have' a different lenght scale
       ! and it would be to expensive to compute each single df_i/de_i using a finite difference
       h=(EPSILON(h)**(1.0_dp/3.0_dp))*l
       ! get an exact machine representable number close to this h
       h=2.0_dp**EXPONENT(h)
       ! this should write three times the same number
       ! write(*,*) h,(e(i)+h)-e(i),(e(i)-h)-e(i)
       ! and the symmetric finite difference
       ex(:)=e
       ex(i)=e(i)+h
       CALL FermiFixed(fx, mux, kTS, ex, N, T, maxocc, my_estate, my_festate)
       dfde(:,I)=fx
       ex(i)=e(i)-h
       CALL FermiFixed(fx, mux, kTS, ex, N, T, maxocc, my_estate, my_festate)
       dfde(:,I)=(dfde(:,I)-fx)/(2.0_dp*h)
    ENDDO
    DEALLOCATE(ex, fx)

    CALL FermiFixed(f, mu, kTS, e, N, T, maxocc, my_estate, my_festate)

    CALL timestop(handle)

 END SUBROUTINE FermiFixedDeriv

END MODULE fermi_utils
