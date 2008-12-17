!
! wrap all possible functions together, they need unique ids 
!
MODULE functions
    USE mpfr
    USE mpfr_ops
    USE mpfr_cutoff_gamma
    USE mpfr_yukawa
    USE function_types
    IMPLICIT NONE

    INTEGER, SAVE :: functionid=-1
    INTEGER, SAVE :: should_output=-1

CONTAINS

      ! evaluate the function, + nderiv other functions in the point X1,X2 in a reduced domain 
      ! as a first step, this transforms X1,X2 into the natural variables of the original domain (e.g. R,T)
      ! write additionally to an open unit (should_output>0) connected to a file, for later checking the quality of the interpolation
      SUBROUTINE f(res,nderiv,x1,x2)
        INTEGER         :: nderiv
        TYPE(mpfr_type) :: x1,x2
        TYPE(mpfr_type) :: res(0:nderiv)

        TYPE(mpfr_type) :: T,r,upper,lower,zero
        TYPE(mpfr_type) :: dummy(0:21)
        INTEGER          :: i

        zero=.CONVERT."0"

        SELECT CASE(functionid)
        CASE(fun_trunc_coulomb_farfield,fun_trunc_coulomb_nearfield)
            ! 1 = farfield  (R>11)
            ! 2 = nearfield (R<11)
            SELECT CASE(functionid)
            CASE (fun_trunc_coulomb_farfield)
              ! if R=+Infinity the result is zero
              IF (mpfr_cmp(x2,zero)<=0) THEN
                 res=.CONVERT."0"
                 return
              ENDIF
              r=11/x2
              upper=r**2 + 11*r + 50
              lower=r**2 - 11*r
            CASE (fun_trunc_coulomb_nearfield)
              R=x2*11
              upper=r**2 + 11*r + 50
              lower=.CONVERT."0"
            CASE DEFAULT
              STOP "Function ID not implemented"
            END SELECT

            t=lower+x1*(upper-lower)

            !t is zero, return the limiting expansion
            IF (mpfr_cmp(t,zero)<=0) THEN
               CALL cutoff_gamma_T0(21,R,dummy)
            ELSE
               CALL cutoff_gamma(21,t,r,dummy)
            END IF
            res(0:nderiv)=dummy(0:nderiv)
        CASE(fun_yukawa)
           ! deal with infinite T,R locally
           IF((X1==zero) .OR. (X2==zero)) THEN
              res=zero
              return
           ELSE
              T=(1/X1-1)**2
              R=(1/X2-1)**2
              CALL yukawa_gn_all(nderiv,T,R,dummy)
           ENDIF
           res(0:nderiv)=dummy(0:nderiv)
        CASE DEFAULT
           STOP "Function ID not implemented"
        END SELECT

        IF (should_output>0) THEN
!$OMP CRITICAL
           WRITE(should_output,*) REAL(R),REAL(T),REAL(res(0:nderiv))
!$OMP END CRITICAL
        ENDIF

      END SUBROUTINE f

END MODULE functions
