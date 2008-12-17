#if defined(__T_C_G0)
  USE T_C_G0
#endif
#if defined(__YUKAWA)
  USE yukawa
#endif
  IMPLICIT NONE
  INTEGER :: Nder
  REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: ref
  REAL(KIND=dp) :: t1,t2,R,T
  REAL(KIND=dp), ALLOCATABLE :: RES(:),MAXERR(:)
  INTEGER :: I,J,Nline
  LOGICAL :: use_gamma


#if defined(__T_C_G0)
  OPEN(12,FILE="T_C_G.dat")
#endif
#if defined(__YUKAWA)
  OPEN(12,FILE="yukawa.dat")
#endif

  Nder=MIN(9,nderiv_max)
  CALL INIT(Nder,12)
  CLOSE(12)
  ALLOCATE(RES(0:Nder),MAXERR(0:Nder))

  Nline=0
  OPEN(13,FILE="history.dat")
  DO
     READ(13,*,END=999) R,T 
     Nline=Nline+1
  ENDDO
999 CONTINUE
  REWIND(13)
  ALLOCATE(ref(-2:Nder,Nline))
  DO I=1,Nline
     READ(13,*) ref(:,I)
  ENDDO
  CLOSE(13)

  maxerr=0
  write(6,'(A)')       "List of points with largest errors so far"
  write(6,'(A)')       "-----------------------------------------"
  CALL CPU_TIME(T1)
  DO I=1,Nline
   R=ref(-2,I)
   T=ref(-1,I)
#if defined(__T_C_G0)
   CALL T_C_G0_n(RES,use_gamma,R,T,Nder)
   ! reference points can be outside the expected domain due to 
   ! the fact that they are rounded from arbitrary precision to 
   ! double precision
   IF (.NOT.use_gamma) THEN
     IF (ANY(RES.NE.0)) THEN
       IF (ANY((abs(RES-REF(0:Nder,I)))>(maxerr))) THEN
         maxerr=MAX(maxerr,abs(RES-REF(0:Nder,I)))
         write(6,*) R,T,RES(0:Nder),REF(0:Nder,I)
       ENDIF
     ENDIF
   ELSE
     ! write(6,*) "oops",R,T,R**2-11*R-T
   ENDIF
#endif
#if defined(__YUKAWA)
   CALL yukawa_n(RES,R,T,Nder)
   IF (ANY((abs(RES-REF(0:Nder,I)))>(maxerr))) THEN
     maxerr=MAX(maxerr,abs(RES-REF(0:Nder,I)))
     write(6,*) R,T,RES(0:Nder),REF(0:Nder,I)
   ENDIF
#endif
  ENDDO
  CALL CPU_TIME(T2)
  write(6,'(A)')       "Summary of generation"
  write(6,'(A)')       "---------------------"
  write(6,'(A,I8)')    "Degree of the interpolation:                 ",degree
  write(6,'(A,E8.3)')  "Target absolute error:                       ",target_error
  write(6,'(A,I8)')    "Number of domains after bisection:           ",patches
  write(6,'(A,I8)')    "Table size for polynomial coefficients [Kb]: ",(8*patches*(Nder+1)*(degree+1)*(degree+2)/2)/(1024)
  write(6,'(A,I8)')    "Number of reference points:                  ",Nline
  write(6,'(A,I8)')    "Tested up to derivative:                     ",Nder
  write(6,'(A,F8.3)')  "Time for evaluation [s]:                     ",t2-t1
  write(6,'(A,E8.3)')  "Maximum error over all reference points:     ",MAXVAL(maxerr)
  write(6,'(A)')       "Maximum error per derivative:"
  write(6,'(5E18.4)')  maxerr

END

