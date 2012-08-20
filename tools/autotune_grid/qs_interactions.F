MODULE qs_interactions
 USE kinds
 USE orbital_pointers
 IMPLICIT NONE
CONTAINS
    FUNCTION exp_radius(l,alpha,threshold,prefactor,epsin) RESULT(radius)
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), INTENT(IN)                :: alpha, threshold, prefactor
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: epsin
    REAL(KIND=dp)                            :: radius

    INTEGER, PARAMETER                       :: maxiter = 5000

    INTEGER                                  :: iter
    REAL(KIND=dp)                            :: a, ar2, d, epsiter, g, r, &
                                                rhigh, rlow, rmid, t

!   ---------------------------------------------------------------------------

    IF (PRESENT(epsin)) THEN
       epsiter=epsin
    ELSE
       epsiter=EPSILON(epsiter)*(1.0E-12_dp / 2.22044604925031E-16_dp)
    ENDIF

!   *** Initialize function value ***

    radius = 0.0_dp

!   *** Load and check parameter values ***

    IF (l < 0) THEN
    END IF

    IF (alpha == 0.0_dp) THEN
    ELSE
      a = ABS(alpha)
    END IF

    IF (threshold == 0.0_dp) THEN
    ELSE
      t = ABS(threshold)
    END IF

    IF (prefactor == 0.0_dp) THEN
      RETURN
    ELSE
      d = ABS(prefactor)
    END IF

!   *** Calculate the Maximum g(r) ***
    r = SQRT(0.5_dp*REAL(l,dp)/a)
    ar2 = a*r*r

    IF (l == 0) THEN
      g = d
    ELSE
      g = d*r**l*EXP(-ar2)
    END IF

    IF (t > g) THEN
      RETURN
    END IF

    rlow  = r
    rhigh = 2.0_dp*rlow+1.0_dp
    iter=0
    DO 
       iter=iter+1
       g = d*rhigh**l*EXP(-a*rhigh**2)
       IF (g < t) EXIT
       rlow = rhigh
       rhigh =  2.0_dp*rlow+1.0_dp 
    ENDDO

    DO iter=1,maxiter
       rmid= (rlow+rhigh)*0.5_dp
       ar2 = a*rmid*rmid
       g   = d*rmid**l*EXP(-ar2)
       IF  (g.lt.t) THEN
           rhigh=rmid
       ELSE
           rlow=rmid
       ENDIF
       IF (ABS(rhigh-rlow).lt.epsiter) THEN
          radius=rhigh
          RETURN
       ENDIF
    ENDDO 
     
  END FUNCTION exp_radius


 FUNCTION exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,o1,o2,ra,rb,rp,&
                          zetp,eps,prefactor,cutoff) RESULT(radius)

    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, lb_max
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: pab
    INTEGER, OPTIONAL                        :: o1, o2
    REAL(KIND=dp), INTENT(IN)                :: ra(3), rb(3), rp(3), zetp, &
                                                eps, prefactor, cutoff
    REAL(KIND=dp)                            :: radius,rad_a,rad_b

    INTEGER                                  :: d, i, ico, j, jco, k, la(3), &
                                                lb(3), lxa, lxb, lya, lyb, &
                                                lza, lzb
    REAL(KIND=dp)                            :: bini, binj, coef(0:20), &
                                                polycoef(0:60), &
                                                prefactor_local, s1, s2


    ! get the local prefactor, we'll now use the largest density matrix element of the block to screen
    IF (PRESENT(pab)) THEN
       prefactor_local=cutoff
       DO lxa=0,la_max
       DO lxb=0,lb_max
          DO lya=0,la_max-lxa
          DO lyb=0,lb_max-lxb
             DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
             DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
                la=(/lxa,lya,lza/)
                lb=(/lxb,lyb,lzb/)
                ico=coset(lxa,lya,lza)
                jco=coset(lxb,lyb,lzb)
                prefactor_local=MAX(ABS(pab(o1+ico,o2+jco)),prefactor_local)
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       prefactor_local=prefactor*prefactor_local
    ELSE
       prefactor_local=prefactor*MAX(1.0_dp,cutoff)
    ENDIF

    rad_a=SQRT(SUM((ra-rp)**2))
    rad_b=SQRT(SUM((rb-rp)**2))

    polycoef(0:la_max+lb_max)=0.0_dp
    DO lxa=0,la_max
    DO lxb=0,lb_max
       coef(0:la_max+lb_max)=0.0_dp
       bini=1.0_dp
       s1=1.0_dp
       DO i=0,lxa
          binj=1.0_dp
          s2=1.0_dp
          DO j=0,lxb
             coef(lxa+lxb-i-j)=coef(lxa+lxb-i-j) + bini*binj*s1*s2
             binj=(binj*(lxb-j))/(j+1)
             s2=s2*(rad_b)
          ENDDO
          bini=(bini*(lxa-i))/(i+1)
          s1=s1*(rad_a)
       ENDDO
       DO i=0,lxa+lxb
          polycoef(i)=MAX(polycoef(i),coef(i))
       ENDDO
    ENDDO
    ENDDO

    polycoef(0:la_max+lb_max)=polycoef(0:la_max+lb_max)*prefactor_local
    radius=0.0_dp
    DO i=0,la_max+lb_max
          radius=MAX(radius,exp_radius(i,zetp,eps,polycoef(i),1.0E-2_dp) )
    ENDDO

  END FUNCTION exp_radius_very_extended

END MODULE qs_interactions
