!
! generate a benchmark for the following cases
!
! 1) tiny_gemm
! 2) matmul
! 3) dgemm
! 4) multrec 1
! 5) multrec 2
! 6) multrec 3
! 7) multrec 4
!
PROGRAM small_gen
   USE mults
   USE multrec_gen
   IMPLICIT NONE
   
   INTEGER :: M,N,K,nline,iline,max_dim, best_square(4),i,isquare, opts(4), transpose_flavor, data_type
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: tiny_opts
   REAL, DIMENSION(:), ALLOCATABLE :: tiny_perf,square_perf
   REAL :: tmp
   CHARACTER(LEN=1024) :: arg,filename,line,label
   
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
   
   ! generation of the tiny version
   write(label,'(A,I0)') "_",1
   CALL mult_versions(M,N,K,1,label,transpose_flavor,data_type)

   ! generation of the matmul version
   write(label,'(A,I0)') "_",2
   CALL mult_versions(M,N,K,2,label,transpose_flavor,data_type)

   ! generation of the dgemm version
   write(label,'(A,I0)') "_",3
   CALL mult_versions(M,N,K,3,label,transpose_flavor,data_type)

   ! generation of the multrec versions (4)
   DO isquare=1,4
      write(label,'(A,I0)') "_",3+isquare
      CALL mult_versions(M,N,K,3+isquare,label,transpose_flavor,data_type)
   ENDDO

   ! test function
   CALL write_test_fun(M,N,K,transpose_flavor,data_type)

   ! test program
   write(6,*)                    " PROGRAM small_find"
   write(6,*)                    "    IMPLICIT NONE"
   write(6,'(A,I0,A,I0,A,I0,A,I0)') "    INTEGER, PARAMETER :: M=",&
         M,",N=",N,",K=",K,",Nmin=5,versions=",3+SIZE(best_square)
   CALL write_matrix_defs(M,N,K,transpose_flavor,data_type)
   write(6,*)                    "    REAL         :: timing(versions), best_time, test"
   write(6,*)                    "    REAL(KIND=KIND(0.D0)) :: flops,gflop"
   write(6,*)                    "    INTEGER      :: imin,Niter,iloop,i,j,l, best_loop,best_mu, best_nu, best_ku"
   write(6,*)                    "    INTERFACE"
   write(6,*)                    "      SUBROUTINE X(A,B,C)"
   CALL write_matrix_defs(M,N,K,transpose_flavor,data_type)
   write(6,*)                    "      END SUBROUTINE"
   write(6,*)                    "    END INTERFACE"
   DO i=1,3+SIZE(best_square)
      write(6,'(A,I0,A,I0,A,I0,A,I0)') "PROCEDURE(X) :: smm_"//trstr(transpose_flavor,data_type)//"_",&
          M,"_",N,"_",K,"_",i
   ENDDO
   write(6,*)                    ""
   write(6,*)                    "    flops=2*REAL(M,KIND=KIND(0.D0))*N*K"
   write(6,*)                    "    gflop=1000.0D0*1000.0D0*1000.0D0"
   write(6,*)                    "    ! assume we would like to do 10 Gflop for testing a subroutine"
   write(6,*)                    "    Niter=MAX(1,CEILING(MIN(100000000.0D0,5*gflop/flops)))"
   write(6,*)                    ""
   write(6,*)                    "    best_time=HUGE(best_time)"
   write(6,*)                    "    timing=best_time"
   write(6,*)                    "    C=0 ; A=0 ; B=0  "
   write(6,*)                    ""
   write(6,*)                    "    DO imin=1,Nmin"
   DO i=1,3+SIZE(best_square)
          write(6,*) "       timing(",i,")= &"
          write(6,*) "       MIN(timing(",i,"), &"
          write(6,'(A,I0,A,I0,A,I0,A,I0,A)') "   TEST(smm_"//trstr(transpose_flavor,data_type)//"_",&
     M,"_",N,"_",K,"_",i,",A,B,C,Niter))"
          write(6,*) 'write(6,''(1I4,F12.6,F12.3)'') ',i,',&'
          write(6,*) "timing(",i,"),&"
          write(6,*) "flops*Niter/gflop/timing(",i,")"
   ENDDO
   write(6,*)                    "    ENDDO"
   write(6,*)                    ""
   write(6,*)                    "    DO iloop=1,versions"
   write(6,*)                    "       IF (timing(iloop)< best_time) THEN"
   write(6,*)                    "          best_time=timing(iloop)"
   write(6,*)                    "          best_loop=iloop"
   write(6,*)                    "       ENDIF"
   write(6,*)                    "    ENDDO"
   write(6,*)                    '    write(6,''(1I4,F12.6,F12.3)'') '//&
                   'best_loop, best_time, (flops*Niter/best_time)/gflop'
   write(6,*)                    "END PROGRAM small_find"

END PROGRAM small_gen
