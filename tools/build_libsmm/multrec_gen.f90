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
! 8) Vector version
!
MODULE multrec_gen
  USE mults
  IMPLICIT NONE
CONTAINS
  SUBROUTINE find_tiny_opts(opts,tiny_opts,m,n,k)
     INTEGER, INTENT(OUT) :: opts(4)
     INTEGER, DIMENSION(:,:) :: tiny_opts
     INTEGER :: m,n,k
     INTEGER :: i
     opts=(/5,1,1,1/)
     DO i=1,SIZE(tiny_opts,2)
        IF (ALL(tiny_opts(1:3,i)==(/m,n,k/))) opts=tiny_opts(4:7,i)
     ENDDO
  END SUBROUTINE find_tiny_opts

  RECURSIVE SUBROUTINE MULTREC(mi,mf,ni,nf,ki,kf,block_size,tiny_opts,transpose_flavor,data_type)
    INTEGER :: mi,mf,ni,nf,ki,kf,block_size, tiny_opts(:,:), transpose_flavor,data_type
    INTEGER :: M,N,K,opts(4)
    INTEGER :: cut,s1

    M=mf-mi+1
    N=nf-ni+1
    K=kf-ki+1

    ! small sizes are done directly, otherwise we recurse
    IF (M<=block_size .AND. N<=block_size .AND. K<=block_size) THEN
       CALL find_tiny_opts(opts,tiny_opts,m,n,k)
       CALL smm_inner(mi,mf,ni,nf,ki,kf,opts(1),opts(2),opts(3),opts(4),transpose_flavor,data_type)
    ELSE
       ! a three case recursion
       IF (M>=MAX(N,K)) cut=1
       IF (K>=MAX(N,M)) cut=2
       IF (N>=MAX(M,K)) cut=3
       SELECT CASE(cut)
       CASE(1)
          s1=((M/2+block_size-1)/block_size)*block_size
          CALL MULTREC(mi,mi+s1-1,ni,nf,ki,kf,block_size,tiny_opts,transpose_flavor,data_type)
          CALL MULTREC(mi+s1,mf,ni,nf,ki,kf,block_size,tiny_opts,transpose_flavor,data_type)
       CASE(2)
          s1=((K/2+block_size-1)/block_size)*block_size
          CALL MULTREC(mi,mf,ni,nf,ki,ki+s1-1,block_size,tiny_opts,transpose_flavor,data_type)
          CALL MULTREC(mi,mf,ni,nf,ki+s1,kf,block_size,tiny_opts,transpose_flavor,data_type)
       CASE(3)
          s1=((N/2+block_size-1)/block_size)*block_size
          CALL MULTREC(mi,mf,ni,ni+s1-1,ki,kf,block_size,tiny_opts,transpose_flavor,data_type)
          CALL MULTREC(mi,mf,ni+s1,nf,ki,kf,block_size,tiny_opts,transpose_flavor,data_type)
       END SELECT
    ENDIF
  END SUBROUTINE MULTREC

  FUNCTION trsum(last)
    LOGICAL :: last
    CHARACTER(LEN=25) :: trsum
    IF (last) THEN
       trsum=""
    ELSE
       trsum="+ &"
    ENDIF
  END FUNCTION trsum
  
  SUBROUTINE write_subroutine_stack(label,M,N,K,transpose_flavor,data_type,stack_size_label,Cbuffer_row,Cbuffer_col)
    CHARACTER(LEN=*) :: label
    INTEGER :: M,N,K,transpose_flavor,data_type
    CHARACTER(LEN=*), OPTIONAL :: stack_size_label
    INTEGER, OPTIONAL :: Cbuffer_row, Cbuffer_col

    IF (PRESENT(stack_size_label).AND.stack_size_label/="") THEN
       write(6,'(A,I0,A,I0,A,I0,A)') "   SUBROUTINE smm_"//trstr(transpose_flavor,data_type)//"_",&
            M,"_",N,"_",K,"_stack"//TRIM(label)//"("//TRIM(trparam(stack_size_label))//")"
       CALL write_stack_params(data_type,stack_size_label)
       write(6,'(A)')                    "      INTEGER            :: sp"
       IF (PRESENT(Cbuffer_row).AND.PRESENT(Cbuffer_col)) THEN
          write(6,'(A,I0,A,I0,A)')   "       "//trdat(data_type)//":: Cbuffer(",Cbuffer_row,",",Cbuffer_col,")"
       ENDIF
       write(6,'(A)')                    "      DO sp = 1, "//TRIM(stack_size_label)
       IF (PRESENT(Cbuffer_row).AND.PRESENT(Cbuffer_col)) THEN
          write(6,'(A,I0,A,I0,A,I0,A)') "         CALL smm_"//trstr(transpose_flavor,data_type)//"_",&
               M,"_",N,"_",K,TRIM(label)//"_buffer(A(params(p_a_first,sp)),B(params(p_b_first,sp)),C(params(p_c_first,sp)),Cbuffer)"
       ELSE
       write(6,'(A,I0,A,I0,A,I0,A)') "         CALL smm_"//trstr(transpose_flavor,data_type)//"_",&
            M,"_",N,"_",K,TRIM(label)//"(A(params(p_a_first,sp)),B(params(p_b_first,sp)),C(params(p_c_first,sp)))"
       ENDIF
       write(6,'(A)')                  "      ENDDO"
       write(6,'(A)') "   END SUBROUTINE"
    ENDIF
  END SUBROUTINE write_subroutine_stack

  SUBROUTINE MULTVECTOR(label,M,N,K,transpose_flavor,data_type,nSIMD,stride,stack_size_label)
    INTEGER :: M,N,K,sj,je,ji,sl,le,li
    INTEGER :: transpose_flavor,data_type,nSIMD,stride
    INTEGER :: multElements,modElements
    CHARACTER(LEN=*) :: label
    CHARACTER(LEN=*), OPTIONAL :: stack_size_label

    multElements=(M/nSIMD)*nSIMD
    modElements=MOD(M,nSIMD)

    IF (modElements>0) THEN
       if (PRESENT(stack_size_label).AND.stack_size_label/="") THEN
          write(6,'(A,I0,A,I0,A,I0,A)')    "   PURE SUBROUTINE smm_"//trstr(transpose_flavor,data_type)//"_",&
               M,"_",N,"_",K,TRIM(label)//"_buffer(A,B,C,Cbuffer)"
          write(6,'(A,I0,A,I0,A)') "      "//trdat(data_type,"INOUT")//" :: Cbuffer(",nSIMD,",",MIN(stride,N),")"
       ELSE
          write(6,'(A,I0,A,I0,A,I0,A)')    "   PURE SUBROUTINE smm_"//trstr(transpose_flavor,data_type)//"_",&
               M,"_",N,"_",K,TRIM(label)//"(A,B,C)"
          write(6,'(A,I0,A,I0,A)') "      "//trdat(data_type)//" :: Cbuffer(",nSIMD,",",MIN(stride,N),")"
       ENDIF
    ELSE
       if (PRESENT(stack_size_label)) THEN
          write(6,'(A,I0,A,I0,A,I0,A)')    "   PURE SUBROUTINE smm_"//trstr(transpose_flavor,data_type)//"_",&
               M,"_",N,"_",K,TRIM(label)//"(A,B,C)"
       ELSE
          RETURN
       ENDIF
    ENDIF

    CALL write_matrix_defs(M,N,K,transpose_flavor,data_type,.TRUE.) 
    
    write(6,'(A)') "      INTEGER :: i"

    sj=stride ! blocking dimension in N
    sl=stride ! blocking dimension in K

    DO je=1,N,sj

       DO le=1,K,sl

          IF (multElements>0) THEN
             write(6,'(A,I0,A,I0,A,I0)') "     DO i=",1,",",multElements,",",1
             DO ji=je,MIN(je+sj-1,N),1
                write(6,'(A,I0,A,I0,A)') "       C(i,",ji,")=C(i,",ji,")+ &"
                DO li=le,MIN(le+sl-1,K),1
                   write (6,'(A,I0,A,I0,A,I0,A)') "         A(i,",&
                     li,")*B(",li,",",ji,")"//trsum(li==MIN(le+sl-1,K))
                ENDDO
             ENDDO
             write(6,'(A)') "     ENDDO "
          ENDIF

          ! consider remaining elements
          IF (modElements>0) THEN
             write(6,'(A,I0,A,I0,A,I0)') "     DO i=",1,",",nSIMD,",",1
             DO ji=je,MIN(je+sj-1,N),1
                IF (le>1) THEN
                   write(6,'(A,I0,A,I0,A)') "       Cbuffer(i,",&
                     MOD(ji-1,sj)+1,")=Cbuffer(i,",MOD(ji-1,sj)+1,")+ &"
                ELSE
                   write(6,'(A,I0,A,I0,A,I0,A)') "       Cbuffer(i,",&
                     MOD(ji-1,sj)+1,")=C(i+",multElements,",",ji,")+ &"
                ENDIF

                DO li=le,MIN(le+sl-1,K),1
                   write (6,'(A,I0,A,I0,A,I0,A,I0,A)') "         A(i+",&
                     multElements,",",li,")*B(",li,",",ji,")"//trsum(li==MIN(le+sl-1,K))
                ENDDO
             ENDDO
             write(6,'(A)') "     ENDDO "
          
          ENDIF
       ENDDO

       ! copy the remaining elements
       IF (modElements>0) THEN
          write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') "      C(",&
            1+multElements,":",M,",",je,":",MIN(je+sj-1,N),")=Cbuffer(", &
            1,":",modElements,",1:",MIN(je+sj,N+1)-je,")"
          
       ENDIF

    ENDDO

    write(6,'(A)') "   END SUBROUTINE"

    if (PRESENT(stack_size_label)) THEN

       IF (modElements>0.AND.nSIMD>0.AND.stride>0) THEN
          CALL write_subroutine_stack(label,M,N,K,transpose_flavor,data_type,stack_size_label,nSIMD,MIN(stride,N))
       ELSE
          CALL write_subroutine_stack(label,M,N,K,transpose_flavor,data_type,stack_size_label)
       ENDIF
    ENDIF

  END SUBROUTINE MULTVECTOR

  SUBROUTINE mult_versions(M,N,K,version,label,transpose_flavor,data_type,SIMD_size,filename,&
                           stack_size_label,write_buffer_interface)
     INTEGER :: M,N,K,version,transpose_flavor,data_type,SIMD_size
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: tiny_opts
     INTEGER :: best_square(4)
     REAL, ALLOCATABLE, DIMENSION(:)      :: tiny_perf,square_perf
     CHARACTER(LEN=1024) :: filename,line
     CHARACTER(LEN=*), OPTIONAL :: stack_size_label
     LOGICAL, OPTIONAL :: write_buffer_interface
     CHARACTER(LEN=*) :: label
     INTEGER :: opts(4),blocksize,i,iline,nline,max_dim,isquare
     REAL :: tmp
     INTEGER :: size_type, nSIMD
     INTEGER, PARAMETER :: stride=8 ! used for the unrolling
     size_type=0; nSIMD=0

     ! only in the case of SIMD_size=32(i.e. AVX) and SIMD_size=64(i.e. MIC)
     IF ((SIMD_size==32 .OR. SIMD_size==64) .AND. transpose_flavor==1 .AND. data_type<=2 .AND. &
          (LABEL=="" .OR. version==8)) THEN

        SELECT CASE(data_type)
        CASE(1)
           size_type=8 !double precision bytes
        CASE(2)
           size_type=4 !single precision bytes
        END SELECT
     
        nSIMD=SIMD_size/size_type
     ENDIF

     !
     ! filename is the result of tiny optimization (cat tiny_gen_optimal.out)
     ! 1 1 1    5   1   1   1    0.376023       0.532
     !
     OPEN(UNIT=10,FILE=filename)
     REWIND(10)
     nline=0
     DO
       READ(10,*,END=999) line
       nline=nline+1
     ENDDO
999  CONTINUE
     ALLOCATE(tiny_opts(7,nline))
     ALLOCATE(tiny_perf(nline))
     REWIND(10)
     DO iline=1,nline
       READ(10,'(A1024)',END=999) line
       READ(line,*) tiny_opts(:,iline),tmp,tiny_perf(iline)
     ENDDO
     CLOSE(10)

     ! find square sizes that give good performance with tiny opts
     max_dim=MAXVAL(tiny_opts(1:3,:))
     ALLOCATE(square_perf(max_dim))
     square_perf=-1
     DO iline=1,nline
        IF (tiny_opts(1,iline)==tiny_opts(2,iline) .AND. tiny_opts(1,iline)==tiny_opts(3,iline)) THEN
           square_perf(tiny_opts(1,iline))=tiny_perf(iline)
        ENDIF
     ENDDO
     best_square=-1
     DO isquare=1,SIZE(best_square)
        tmp=-HUGE(tmp)
        DO i=1,max_dim
           IF (square_perf(i)>tmp .AND. .NOT. ANY(best_square.EQ.i)) THEN
              tmp=square_perf(i)
              best_square(isquare)=i
           ENDIF
        ENDDO
     ENDDO
     IF (ANY(best_square<1)) ERROR STOP "tiny opts file needs sufficiently many square sizes"

     IF (version.ge.1.and.version.le.7) THEN
        IF (version.ne.3) THEN
           write(6,'(A,I0,A,I0,A,I0,A)')    "   PURE SUBROUTINE smm_"//trstr(transpose_flavor,data_type)//"_",&
                M,"_",N,"_",K,TRIM(label)//"(A,B,C)"
           ELSE
              write(6,'(A,I0,A,I0,A,I0,A)')    "   SUBROUTINE smm_"//trstr(transpose_flavor,data_type)//"_",&
                   M,"_",N,"_",K,TRIM(label)//"(A,B,C)"
           ENDIF
        CALL write_matrix_defs(M,N,K,transpose_flavor,data_type,.TRUE.) 
     ENDIF

     SELECT CASE(version)
     CASE(1) 
       ! generation of the tiny version
       write(6,'(A)')                     "      INTEGER ::i,j,l"
       CALL find_tiny_opts(opts,tiny_opts,m,n,k)
       CALL smm_inner(1,M,1,N,1,K,opts(1),opts(2),opts(3),opts(4),transpose_flavor,data_type)
     CASE(2)
       ! generation of the matmul version
       SELECT CASE(transpose_flavor)
       CASE(1)
          write(6,'(A)')                  "      C = C + MATMUL(A,B) ! so easy"
       CASE(2)
          write(6,'(A)')                  "      C = C + MATMUL(TRANSPOSE(A),B) ! so easy"
       CASE(3)
          write(6,'(A)')                  "      C = C + MATMUL(A,TRANSPOSE(B)) ! so easy"
       CASE(4)          
          write(6,'(A)')                  "      C = C + MATMUL(TRANSPOSE(A),TRANSPOSE(B)) ! so easy"
       END SELECT
     CASE(3)
       ! generation of the gemm version
       WRITE(6,'(A)') "      "//trdat(data_type)//", PARAMETER :: one=1"
       SELECT CASE(transpose_flavor)
       CASE(1)
         write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
               "      CALL "//trgemm(data_type)//"('N','N',",M,",",N,",",K,",one,A,",M,",B,",K,",one,C,",M,")"
       CASE(2)
         write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
               "      CALL "//trgemm(data_type)//"('T','N',",M,",",N,",",K,",one,A,",K,",B,",K,",one,C,",M,")"
       CASE(3)
         write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
               "      CALL "//trgemm(data_type)//"('N','T',",M,",",N,",",K,",one,A,",M,",B,",N,",one,C,",M,")"
       CASE(4)
         write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
               "      CALL "//trgemm(data_type)//"('T','T',",M,",",N,",",K,",one,A,",K,",B,",N,",one,C,",M,")"
       END SELECT
     CASE(4,5,6,7)
       isquare=version-3
       ! generation of the multrec versions
       write(6,'(A)')                     "      INTEGER ::i,j,l"
       blocksize=best_square(isquare)
       CALL MULTREC(1,M,1,N,1,K,blocksize,tiny_opts,transpose_flavor,data_type)
     CASE(8)
        ! generation of the vector version
        IF (nSIMD>0) THEN
           CALL MULTVECTOR(label,M,N,K,transpose_flavor,data_type,nSIMD,stride,stack_size_label)
           IF (PRESENT(write_buffer_interface).AND.write_buffer_interface) THEN
              CALL MULTVECTOR(label,M,N,K,transpose_flavor,data_type,nSIMD,stride)
           ENDIF

        ENDIF
     CASE DEFAULT
       ERROR STOP "MISSING CASE mult_versions"
     END SELECT

     IF (version.ge.1.and.version.le.7) THEN
        write(6,'(A)') "   END SUBROUTINE"
        CALL write_subroutine_stack(label,M,N,K,transpose_flavor,data_type,stack_size_label)
     ENDIF

  END SUBROUTINE mult_versions

END MODULE multrec_gen
