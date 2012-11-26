PROGRAM tiny_gen
USE mults
IMPLICIT NONE

INTEGER :: M,N,K
CHARACTER(LEN=80) :: arg
INTEGER :: mi,mf,ni,nf,ki,kf,iloop,ku,nu,mu,transpose_flavor, data_type

CALL GET_COMMAND_ARGUMENT(1,arg)
READ(arg,*) M
CALL GET_COMMAND_ARGUMENT(2,arg)
READ(arg,*) N
CALL GET_COMMAND_ARGUMENT(3,arg)
READ(arg,*) K
CALL GET_COMMAND_ARGUMENT(4,arg)
READ(arg,*) transpose_flavor
CALL GET_COMMAND_ARGUMENT(5,arg)
READ(arg,*) data_type

   CALL loop_variants(1)

   CALL write_test_fun(M,N,K,transpose_flavor,data_type)

  write(6,*)                    " PROGRAM tiny_find"
  write(6,*)                    "    IMPLICIT NONE"
  write(6,'(A,I0,A,I0,A,I0,A)') "    INTEGER, PARAMETER :: M=",M,",N=",N,",K=",K,",Nmin=2"
  CALL write_matrix_defs(M,N,K,transpose_flavor,data_type)
  write(6,*)                    "    REAL         :: timing(6,M,N,K), best_time, test"
  write(6,*)                    "    REAL(KIND=KIND(0.D0)) :: flops,gflop"
  write(6,*)                    "    INTEGER      :: imin,Niter,iloop,i,j,l, best_loop,best_mu, best_nu, best_ku"
  write(6,*)                    "    INTERFACE"
  write(6,*)                    "      SUBROUTINE X(A,B,C)"
  CALL write_matrix_defs(M,N,K,transpose_flavor,data_type)
  write(6,*)                    "      END SUBROUTINE"
  write(6,*)                    "    END INTERFACE"
  CALL loop_variants(2)
  write(6,*)                    ""
  write(6,*)                    "    flops=2*REAL(M,KIND=KIND(0.D0))*N*K"
  write(6,*)                    "    gflop=1000.0D0*1000.0D0*1000.0D0"
  write(6,*)                    "    ! assume we would like to do 1 Gflop for testing a subroutine"
  write(6,*)                    "    Niter=MAX(1,CEILING(MIN(100000000.0D0,1*gflop/flops)))"
  write(6,*)                    ""
  write(6,*)                    "    best_time=HUGE(best_time)"
  write(6,*)                    "    timing=best_time"
  write(6,*)                    "    C=0 ; A=0 ; B=0  "
  write(6,*)                    ""
  write(6,*)                    "    DO imin=1,Nmin"
  CALL loop_variants(3)
  write(6,*)                    "    ENDDO"
  write(6,*)                    ""
  write(6,*)                    "    DO iloop=1,6"
  write(6,*)                    "    DO i=1,M"
  write(6,*)                    "    DO j=1,N"
  write(6,*)                    "    DO l=1,K"
  write(6,*)                    "       IF (timing(iloop,i,j,l)< best_time) THEN"
  write(6,*)                    "          best_time=timing(iloop,i,j,l)"
  write(6,*)                    "          best_loop=iloop"
  write(6,*)                    "          best_mu=i"
  write(6,*)                    "          best_nu=j"
  write(6,*)                    "          best_ku=l"
  write(6,*)                    "       ENDIF"
  write(6,*)                    "    ENDDO"
  write(6,*)                    "    ENDDO"
  write(6,*)                    "    ENDDO"
  write(6,*)                    "    ENDDO"
  write(6,'(A112)')             '    write(6,''(4I4,F12.6,F12.3)'') '//&
                  'best_loop,best_mu, best_nu, best_ku, best_time, (flops*Niter/best_time)/gflop'
  write(6,*)                    "END PROGRAM tiny_find"

CONTAINS

 SUBROUTINE loop_variants(itype)
    INTEGER :: itype

    mi=1 ; mf=M ; ni=1 ; nf=N ; ki=1 ; kf=K
    
    DO iloop=1,6
    DO mu=1,M
    DO nu=1,N
    DO ku=1,K
    ! unroll only if no cleanup is needed
    IF (MOD(M,mu).NE.0) CYCLE
    IF (MOD(N,nu).NE.0) CYCLE
    IF (MOD(K,ku).NE.0) CYCLE
    ! do not unroll with more than 32 mults in the inner loop (some more unrolling can be faster, but we reduce size and speedup lib generation)
    IF (mu*nu*ku>32)   CYCLE

       SELECT CASE(itype)
       CASE(1)
         write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
              "   SUBROUTINE smm_"//trstr(transpose_flavor,data_type)//"_",M,"_",N,"_",K,"_",iloop,"_",mu,"_",nu,"_",ku,"(A,B,C)"
         CALL write_matrix_defs(M,N,K,transpose_flavor,data_type)
         write(6,'(A)')                "      INTEGER ::i,j,l"
         CALL smm_inner(mi,mf,ni,nf,ki,kf,iloop,mu,nu,ku,transpose_flavor,data_type)
         write(6,*) "   END SUBROUTINE"
       CASE(2)
         write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
              "   PROCEDURE(X) :: smm_"//trstr(transpose_flavor,data_type)//"_",M,"_",N,"_",K,"_",iloop,"_",mu,"_",nu,"_",ku
       CASE(3)
         write(6,'(A,I0,A,I0,A,I0,A,I0,A)') "       timing(",iloop,",",mu,",",nu,",",ku,")= &"
         write(6,'(A,I0,A,I0,A,I0,A,I0,A)') "       MIN(timing(",iloop,",",mu,",",nu,",",ku,"), &"
         write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
              "   TEST(smm_"//trstr(transpose_flavor,data_type)//"_",M,"_",N,"_",K,"_",iloop,"_",mu,"_",nu,"_",ku,",A,B,C,Niter))"
         write(6,'(A,I0,A,I0,A,I0,A,I0,A)') 'write(6,''(4I4,F12.6,F12.3)'') ',iloop,",",mu,",",nu,",",ku,",& "
         write(6,'(A,I0,A,I0,A,I0,A,I0,A)') "timing(",iloop,",",mu,",",nu,",",ku,"),&"
         write(6,'(A,I0,A,I0,A,I0,A,I0,A)') "flops*Niter/gflop/timing(",iloop,",",mu,",",nu,",",ku,")"
                    
       END SELECT
    ENDDO
    ENDDO
    ENDDO
    ENDDO

 END SUBROUTINE

END PROGRAM tiny_gen
