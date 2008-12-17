!
! a simple interface to the various functions, useful for testing 
! reads from 'evalf.in' writes to 'evalf.out'
!
PROGRAM main
  USE mpfr_cutoff_gamma
  IMPLICIT NONE

  REAL(dp)       :: T,R,dummy(0:21)
  INTEGER*2      :: precision
  TYPE(mpfr_type):: Tval_mp, Rval_mp, dummy_mp(0:21)
  INTEGER        :: i

  OPEN(UNIT=77,FILE="evalf.in")
  READ(77,*) R,T
  CLOSE(77)
  precision = log(10.0)/log(2.0)*4096
  CALL mpfr_set_default_precision(precision)
  Tval_mp = T
  Rval_mp = R
  IF(T==0.0_dp) THEN
    CALL cutoff_gamma_T0(21,Rval_mp,dummy_mp)
  ELSE
    CALL cutoff_gamma(21,Tval_mp,Rval_mp,dummy_mp)   
  END IF
  DO i=0,21
    dummy(i)=mp_to_real(dummy_mp(i))
  END DO
  OPEN(UNIT=78,FILE="evalf.out")
  write(78,*) R,T,dummy(:)
  CLOSE(78)
END PROGRAM main
