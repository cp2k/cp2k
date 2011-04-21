MODULE mults

  IMPLICIT NONE

CONTAINS
  FUNCTION trdat(data_type)
    INTEGER :: data_type
    CHARACTER(LEN=25) :: trdat
    SELECT CASE(data_type)
    CASE(1)
      trdat="REAL(KIND=KIND(0.0D0))"
    CASE(2)
      trdat="REAL(KIND=KIND(0.0))"
    CASE(3)
      trdat="COMPLEX(KIND=KIND(0.0D0))"
    CASE(4)
      trdat="COMPLEX(KIND=KIND(0.0))"
    END SELECT
  END FUNCTION
  FUNCTION trgemm(data_type)
    INTEGER :: data_type
    CHARACTER(LEN=5) :: trgemm
    SELECT CASE(data_type)
    CASE(1)
      trgemm="DGEMM"
    CASE(2)
      trgemm="SGEMM"
    CASE(3)
      trgemm="ZGEMM"
    CASE(4)
      trgemm="CGEMM"
    END SELECT
  END FUNCTION
  FUNCTION trstr(tranpose_flavor,data_type)
    INTEGER :: tranpose_flavor, data_type
    CHARACTER(LEN=3) :: trstr
    CHARACTER(LEN=1) :: dstr
    SELECT CASE(data_type)
    CASE(1)
     dstr="d"
    CASE(2)
     dstr="s"
    CASE(3)
     dstr="z"
    CASE(4)
     dstr="c"
    END SELECT
    SELECT CASE(tranpose_flavor)
    CASE(1)
     trstr=dstr//"nn"
    CASE(2)
     trstr=dstr//"tn"
    CASE(3)
     trstr=dstr//"nt"
    CASE(4)
     trstr=dstr//"tt"
    END SELECT
  END FUNCTION trstr

  SUBROUTINE write_matrix_defs(M,N,K,tranpose_flavor,data_type)
   INTEGER M,N,K,tranpose_flavor,data_type
   SELECT CASE(tranpose_flavor)
   CASE(1)
     write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
        "      "//trdat(data_type)//":: C(",M,",",N,"), B(",K,",",N,"), A(",M,",",K,")"
   CASE(2)
     write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
        "      "//trdat(data_type)//":: C(",M,",",N,"), B(",K,",",N,"), A(",K,",",M,")"
   CASE(3)
     write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
        "      "//trdat(data_type)//":: C(",M,",",N,"), B(",N,",",K,"), A(",M,",",K,")"
   CASE(4)
     write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') & 
        "      "//trdat(data_type)//":: C(",M,",",N,"), B(",N,",",K,"), A(",K,",",M,")"
   END SELECT
  END SUBROUTINE write_matrix_defs

  SUBROUTINE write_test_fun(M,N,K,transpose_flavor,data_type)
   INTEGER :: M,N,K,transpose_flavor,data_type
   write(6,*)                    "FUNCTION TEST(X,A,B,C,N) RESULT(res)"
   write(6,*)                    "      "//trdat(data_type)//" ::C(*), A(*), B(*)"
   write(6,*)                    "   INTEGER :: N"
   write(6,*)                    "   REAL :: t2,t1,res"
   write(6,*)                    "   INTERFACE"
   write(6,*)                    "     SUBROUTINE X(A,B,C)"
   CALL write_matrix_defs(M,N,K,transpose_flavor,data_type)
   write(6,*)                    "     END SUBROUTINE"
   write(6,*)                    "   END INTERFACE"
   write(6,*)                    "   INTEGER :: i"

   write(6,*)                    "   CALL CPU_TIME(t1)"
   write(6,*)                    "   DO i=1,N"
   write(6,*)                    "      CALL X(A,B,C)"
   write(6,*)                    "   ENDDO"
   write(6,*)                    "   CALL CPU_TIME(t2)"
   write(6,*)                    "   res=REAL(t2-t1,KIND=KIND(res))"
   write(6,*)                    "END FUNCTION"
  END SUBROUTINE

  SUBROUTINE smm_inner(mi,mf,ni,nf,ki,kf,iloop,mu,nu,ku,transpose_flavor,data_type)
     INTEGER :: mi,mf,ni,nf,ki,kf,iloop,mu,nu,ku,transpose_flavor,data_type
     INTEGER :: im,in,ik,ido
     INTEGER :: loop_order(3,6),have_loops

     loop_order(:,1)=(/1,2,3/)
     loop_order(:,2)=(/2,1,3/)
     loop_order(:,3)=(/2,3,1/)
     loop_order(:,4)=(/1,3,2/)
     loop_order(:,5)=(/3,1,2/)
     loop_order(:,6)=(/3,2,1/)
     have_loops=0
     CALL out_loop(mi,mf,ni,nf,ki,kf,mu,nu,ku,loop_order(1,iloop),have_loops)
     CALL out_loop(mi,mf,ni,nf,ki,kf,mu,nu,ku,loop_order(2,iloop),have_loops)
     CALL out_loop(mi,mf,ni,nf,ki,kf,mu,nu,ku,loop_order(3,iloop),have_loops) 
     ! what is the fastest order for these loops ? Does it matter ?
     DO im=0,mu-1
     DO in=0,nu-1
     DO ik=0,ku-1
        SELECT CASE(transpose_flavor) 
        CASE(1)
          write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
           "        C(i+",im,",j+",in,")=C(i+",im,",j+",in,")+A(i+",im,",l+",ik,")*B(l+",ik,",j+",in,")"
        CASE(2)
          write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
           "        C(i+",im,",j+",in,")=C(i+",im,",j+",in,")+A(l+",ik,",i+",im,")*B(l+",ik,",j+",in,")"
        CASE(3)
          write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
           "        C(i+",im,",j+",in,")=C(i+",im,",j+",in,")+A(i+",im,",l+",ik,")*B(j+",in,",l+",ik,")"
        CASE(4)
          write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
           "        C(i+",im,",j+",in,")=C(i+",im,",j+",in,")+A(l+",ik,",i+",im,")*B(j+",in,",l+",ik,")"
        END SELECT
     ENDDO
     ENDDO
     ENDDO
     DO ido=1,have_loops
     write(6,*) "     ENDDO "
     ENDDO
  END SUBROUTINE smm_inner

  SUBROUTINE out_loop(mi,mf,ni,nf,ki,kf,mu,nu,ku,ichoice,have_loops)
     INTEGER :: mi,mf,ni,nf,ki,kf,ichoice,mu,nu,ku,have_loops
     IF (ichoice==1) THEN
        IF (nf-ni+1>nu) THEN
           write(6,*) "     DO j=",ni,",",nf,",",nu
           have_loops=have_loops+1
        ELSE
           write(6,*) "     j=",ni 
        ENDIF
     ENDIF
     IF (ichoice==2) THEN
        IF (mf-mi+1>mu) THEN
           write(6,*) "     DO i=",mi,",",mf,",",mu
           have_loops=have_loops+1
        ELSE
           write(6,*) "     i=",mi 
        ENDIF
     ENDIF
     IF (ichoice==3) THEN
        IF (kf-ki+1>ku) THEN
           write(6,*) "     DO l=",ki,",",kf,",",ku
           have_loops=have_loops+1
        ELSE
           write(6,*) "     l=",ki 
        ENDIF
     ENDIF
  END SUBROUTINE

END MODULE mults
