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
! 8) multvec
!
PROGRAM small_gen
   USE mults
   USE multrec_gen
   IMPLICIT NONE

   INTEGER :: M,N,K,nline,iline,max_dim, best_square(4),i,isquare, opts(4), transpose_flavor, data_type, SIMD_size
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: tiny_opts
   REAL, DIMENSION(:), ALLOCATABLE :: tiny_perf,square_perf
   CHARACTER(LEN=1024) :: arg,filename,line,label
   INTEGER :: ibest_square=3

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
   CALL GET_COMMAND_ARGUMENT(6,arg)
   READ(arg,*) SIMD_size
   CALL GET_COMMAND_ARGUMENT(7,filename)

   ! generation of the tiny version
   write(label,'(A,I0)') "_",1
   CALL mult_versions(M,N,K,1,label,transpose_flavor,data_type,SIMD_size,filename)

   ! generation of the matmul version
   write(label,'(A,I0)') "_",2
   CALL mult_versions(M,N,K,2,label,transpose_flavor,data_type,SIMD_size,filename)

   ! generation of the dgemm version
   write(label,'(A,I0)') "_",3
   CALL mult_versions(M,N,K,3,label,transpose_flavor,data_type,SIMD_size,filename)

   ! generation of the multrec versions (4)
   DO isquare=1,SIZE(best_square)
      write(label,'(A,I0)') "_",ibest_square+isquare
      CALL mult_versions(M,N,K,ibest_square+isquare,label,transpose_flavor,data_type,SIMD_size,filename)
   ENDDO
   
   ! generation of the vector version, 
   ! only in the case of SIMD_size=32(i.e. AVX) and SIMD_size=64(i.e. MIC)
   IF ((SIMD_size==32 .OR. SIMD_size==64) .AND. transpose_flavor==1 .AND. data_type<=2) THEN
      ibest_square=ibest_square+1
      write(label,'(A,I0)') "_",ibest_square+SIZE(best_square)
      CALL mult_versions(M,N,K,ibest_square+SIZE(best_square),label,&
           transpose_flavor,data_type,SIMD_size,filename,stack_size_label="")
   ENDIF

   ! test function
   CALL write_test_fun(M,N,K,transpose_flavor,data_type)

   ! test program
   write(6,'(A)')                    " PROGRAM small_find"
   write(6,'(A)')                    "    IMPLICIT NONE"
   write(6,'(A,I0,A,I0,A,I0,A,I0)') "     INTEGER, PARAMETER :: M=",&
         M,",N=",N,",K=",K,",Nmin=5,versions=",ibest_square+SIZE(best_square)

   CALL write_matrix_defs(M,N,K,transpose_flavor,data_type,.FALSE.,padding=.TRUE.)

   write(6,'(A)')                    "    REAL         :: timing(versions), best_time, test"
   write(6,'(A)')                    "    REAL(KIND=KIND(0.D0)) :: flops,gflop"
   write(6,'(A)')                    "    INTEGER      :: imin,Niter,iloop,best_loop"
   write(6,'(A)')                    "    INTERFACE"
   write(6,'(A)')                    "      SUBROUTINE X("//TRIM(trparam())//")"
   CALL write_matrix_defs(data_type=data_type,write_intent=.TRUE.)
   write(6,'(A)')                    "      END SUBROUTINE"
   write(6,'(A)')                    "    END INTERFACE"
   DO i=1,ibest_square+SIZE(best_square)
      write(6,'(A,I0,A,I0,A,I0,A,I0)') "PROCEDURE(X) :: smm_"//trstr(transpose_flavor,data_type)//"_",&
          M,"_",N,"_",K,"_",i
   ENDDO
   write(6,'(A)')                    ""
   write(6,'(A)')                    "    flops=2*REAL(M,KIND=KIND(0.D0))*N*K"
   write(6,'(A)')                    "    gflop=1000.0D0*1000.0D0*1000.0D0"
   write(6,'(A)')                    "    ! assume we would like to do 1 Gflop for testing a subroutine"
   write(6,'(A)')                    "    Niter=MAX(1,CEILING(MIN(100000000.0D0,1*gflop/flops)))"
   write(6,'(A)')                    ""
   write(6,'(A)')                    "    best_time=HUGE(best_time)"
   write(6,'(A)')                    "    timing=best_time"
   write(6,'(A)')                    "    C=0 ; A=0 ; B=0  "
   write(6,'(A)')                    ""
   write(6,'(A)')                    "    DO imin=1,Nmin"
   DO i=1,ibest_square+SIZE(best_square)
          write(6,*) "       timing(",i,")= &"
          write(6,*) "       MIN(timing(",i,"), &"
          write(6,'(A,I0,A,I0,A,I0,A,I0,A)') "   TEST(smm_"//trstr(transpose_flavor,data_type)//"_",&
     M,"_",N,"_",K,"_",i,",A,B,C,Niter))"
          write(6,*) 'write(6,''(1I4,F12.6,F12.3)'') ',i,',&'
          write(6,*) "timing(",i,"),&"
          write(6,*) "flops*Niter/gflop/timing(",i,")"
   ENDDO
   write(6,'(A)')                    "    ENDDO"
   write(6,'(A)')                    ""
   write(6,'(A)')                    "    DO iloop=1,versions"
   write(6,'(A)')                    "       IF (timing(iloop)< best_time) THEN"
   write(6,'(A)')                    "          best_time=timing(iloop)"
   write(6,'(A)')                    "          best_loop=iloop"
   write(6,'(A)')                    "       ENDIF"
   write(6,'(A)')                    "    ENDDO"
   write(6,'(A84)')             '    write(6,''(1I4,F12.6,F12.3)'') '//&
                   'best_loop, best_time, (flops*Niter/best_time)/gflop'
   write(6,'(A)')                    "END PROGRAM small_find"

END PROGRAM small_gen
