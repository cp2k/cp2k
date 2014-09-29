PROGRAM tiny_gen
  USE mults
  IMPLICIT NONE

  INTEGER :: M,N,K,transpose_flavor,data_type
  INTEGER, PARAMETER :: Nloop=6
  CHARACTER(LEN=80) :: arg
  INTEGER :: Nk,mi,mf,ni,nf,ki,kf,iloop,ku,nu,mu

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
  write(6,'(A)')                    "#ifdef __INTEL_OFFLOAD"
  write(6,'(A)')                    "!dir$ attributes offload:mic :: run_kernels"
  write(6,'(A)')                    "#endif"
  write(6,'(A,I0,A,I0,A,I0,A)')     "SUBROUTINE tiny_find_",M,"_",N,"_",K,"(unit)"
  write(6,'(A)')                    "  IMPLICIT NONE"
  write(6,'(A)')                    "  INTEGER :: unit ! Output unit"
  write(6,'(A,I0,A,I0,A,I0)')       "  INTEGER, PARAMETER :: M=",M,",N=",N,",K=",K
  write(6,'(A)')                    "  CHARACTER(len=64) :: filename"
  CALL write_matrix_defs(M,N,K,transpose_flavor,data_type,.FALSE.)
  write(6,'(A)')                    "  INTERFACE"
  write(6,'(A)')                    "    SUBROUTINE X(A,B,C)"
  CALL write_matrix_defs(data_type=data_type,write_intent=.FALSE.)
  write(6,'(A)')                    "    END SUBROUTINE"
  write(6,'(A)')                    "  END INTERFACE"
  CALL loop_variants(2)
  write(6,'(A,I0,A,I0)')            "  INTEGER, PARAMETER :: Nmin=2, Nk=",Nk,", Nloop=",Nloop
  write(6,'(A)')                    "  TYPE t_kernels"
  write(6,'(A)')                    "    PROCEDURE(X), POINTER, NOPASS :: ptr"
  write(6,'(A)')                    "  END TYPE t_kernels"
  write(6,'(A)')                    "  TYPE(t_kernels) :: kernels(Nk,Nloop)"
  write(6,'(A)')                    "  INTEGER :: mnk(3,Nk) ! mu, nu, ku"
  CALL loop_variants(3)
  write(6,'(A,I0,A,I0,A,I0,A)')     "  filename='tiny_find_",M,"_",N,"_",K,".out'"
  write(6,'(A)')                    "  C = 0 ; A = 0 ; B = 0"
  write(6,'(A)')                    "  CALL run_kernels(filename,unit,M,N,K,A,B,C,Nmin,Nk,Nloop,kernels,mnk)"
  write(6,'(A,I0,A,I0,A,I0)')       "END SUBROUTINE tiny_find_",M,"_",N,"_",K

CONTAINS

  SUBROUTINE loop_variants(itype)

    INTEGER :: itype

    mi=1 ; mf=M ; ni=1 ; nf=N ; ki=1 ; kf=K

    DO iloop=1,Nloop

       IF (itype.ne.1) Nk=0

       DO mu=1,M
          DO nu=1,N
             DO ku=1,K
                ! unroll only if no cleanup is needed
                IF (MOD(M,mu).NE.0) CYCLE
                IF (MOD(N,nu).NE.0) CYCLE
                IF (MOD(K,ku).NE.0) CYCLE
                ! do not unroll with more than 32 mults in the inner loop 
                ! (some more unrolling can be faster, but we reduce size and speedup lib generation)
                IF (mu*nu*ku>32)   CYCLE

                SELECT CASE(itype)
                CASE(1)
                   write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
                        "   SUBROUTINE smm_"//trstr(transpose_flavor,data_type)//"_",M,"_",N,"_",K,"_", &
                        iloop,"_",mu,"_",nu,"_",ku,"(A,B,C)"
                   CALL write_matrix_defs(M,N,K,transpose_flavor,data_type,.FALSE.)
                   write(6,'(A)')                "      INTEGER ::i,j,l"
                   CALL smm_inner(mi,mf,ni,nf,ki,kf,iloop,mu,nu,ku,transpose_flavor,data_type)
                   write(6,'(A)') "   END SUBROUTINE"
                CASE(2)
                   Nk=Nk+1
                   write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
                        "  PROCEDURE(X) :: smm_"//trstr(transpose_flavor,data_type)//"_",M,"_",N,"_",K,"_",iloop, &
                        "_",mu,"_",nu,"_",ku
                CASE(3)
                   Nk=Nk+1
                   write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)')  "  kernels(",Nk,",",iloop, & 
                        ")%ptr => smm_"//trstr(transpose_flavor,data_type)//"_",M,"_",N,"_",K,"_", &
                        iloop,"_",mu,"_",nu,"_",ku
                   IF (iloop.eq.1) THEN
                      write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)')              "  mnk(1,",Nk,")=",mu, & 
                                                                              " ; mnk(2,",Nk,")=",nu, &
                                                                              " ; mnk(3,",Nk,")=",ku
                   END IF
                END SELECT
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
  END SUBROUTINE loop_variants

END PROGRAM tiny_gen
