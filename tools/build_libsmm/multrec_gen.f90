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

  RECURSIVE SUBROUTINE MULTREC(mi,mf,ni,nf,ki,kf,block_size,tiny_opts,transpose_flavor)
    INTEGER :: mi,mf,ni,nf,ki,kf,block_size, tiny_opts(:,:), transpose_flavor
    INTEGER :: M,N,K,opts(4)
    INTEGER :: cut,s1

    M=mf-mi+1
    N=nf-ni+1
    K=kf-ki+1

    ! small sizes are done directly, otherwise we recurse
    IF (M<=block_size .AND. N<=block_size .AND. K<=block_size) THEN
       CALL find_tiny_opts(opts,tiny_opts,m,n,k)
       CALL smm_inner(mi,mf,ni,nf,ki,kf,opts(1),opts(2),opts(3),opts(4),transpose_flavor)
    ELSE
       ! a three case recursion
       IF (M>=MAX(N,K)) cut=1
       IF (K>=MAX(N,M)) cut=2
       IF (N>=MAX(M,K)) cut=3
       SELECT CASE(cut)
       CASE(1)
          s1=((M/2+block_size-1)/block_size)*block_size
          CALL MULTREC(mi,mi+s1-1,ni,nf,ki,kf,block_size,tiny_opts,transpose_flavor)
          CALL MULTREC(mi+s1,mf,ni,nf,ki,kf,block_size,tiny_opts,transpose_flavor)
       CASE(2)
          s1=((K/2+block_size-1)/block_size)*block_size
          CALL MULTREC(mi,mf,ni,nf,ki,ki+s1-1,block_size,tiny_opts,transpose_flavor)
          CALL MULTREC(mi,mf,ni,nf,ki+s1,kf,block_size,tiny_opts,transpose_flavor)
       CASE(3)
          s1=((N/2+block_size-1)/block_size)*block_size
          CALL MULTREC(mi,mf,ni,ni+s1-1,ki,kf,block_size,tiny_opts,transpose_flavor)
          CALL MULTREC(mi,mf,ni+s1,nf,ki,kf,block_size,tiny_opts,transpose_flavor)
       END SELECT
    ENDIF
  END SUBROUTINE MULTREC

  SUBROUTINE mult_versions(M,N,K,version,label,transpose_flavor)
     INTEGER :: m,n,k,version,transpose_flavor
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: tiny_opts
     INTEGER :: best_square(4)
     REAL, ALLOCATABLE, DIMENSION(:)      :: tiny_perf,square_perf
     CHARACTER(LEN=1024) :: filename,line
     CHARACTER(LEN=*) :: label
     INTEGER :: opts(4),blocksize,i,iline,nline,max_dim,isquare
     REAL :: tmp

     !
     ! filename is the result of tiny optimization (cat tiny_gen_optimal.out)
     ! 1 1 1    5   1   1   1    0.376023       0.532
     !
     filename="tiny_gen_optimal.out"
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
     IF (ANY(best_square<1)) STOP "tiny opts file needs sufficiently many square sizes"

     SELECT CASE(version)
     CASE(1) 
       ! generation of the tiny version
       write(6,'(A,I0,A,I0,A,I0,A,A)')    "   SUBROUTINE smm_d"//trstr(transpose_flavor)//"_",M,"_",N,"_",K,TRIM(label),"(A,B,C)"
       CALL write_matrix_defs(M,N,K,transpose_flavor)
       write(6,'(A)')                     "      INTEGER ::i,j,l"
       CALL find_tiny_opts(opts,tiny_opts,m,n,k)
       CALL smm_inner(1,M,1,N,1,K,opts(1),opts(2),opts(3),opts(4),transpose_flavor)
       write(6,'(A)') "   END SUBROUTINE"
     CASE(2)
       ! generation of the matmul version
       write(6,'(A,I0,A,I0,A,I0,A,A)')    "   SUBROUTINE smm_d"//trstr(transpose_flavor)//"_",M,"_",N,"_",K,TRIM(label),"(A,B,C)"
       CALL write_matrix_defs(M,N,K,transpose_flavor)
       SELECT CASE(transpose_flavor)
       CASE(1)
          write(6,'(A)')                         "      C = C + MATMUL(A,B) ! so easy"
       CASE(2)
          write(6,'(A)')                         "      C = C + MATMUL(TRANSPOSE(A),B) ! so easy"
       CASE(3)
          write(6,'(A)')                         "      C = C + MATMUL(A,TRANSPOSE(B)) ! so easy"
       CASE(4)
          write(6,'(A)')                         "      C = C + MATMUL(TRANSPOSE(A),TRANSPOSE(B)) ! so easy"
       END SELECT
       write(6,'(A)') "   END SUBROUTINE"
     CASE(3)
       ! generation of the dgemm version
       write(6,'(A,I0,A,I0,A,I0,A,A)')    "   SUBROUTINE smm_d"//trstr(transpose_flavor)//"_",M,"_",N,"_",K,TRIM(label),"(A,B,C)"
       CALL write_matrix_defs(M,N,K,transpose_flavor)
       SELECT CASE(transpose_flavor)
       CASE(1)
         write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
               "      CALL DGEMM('N','N',",M,",",N,",",K,",1.0D0,A,",M,",B,",K,",1.0D0,C,",M,")"
       CASE(2)
         write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
               "      CALL DGEMM('T','N',",M,",",N,",",K,",1.0D0,A,",K,",B,",K,",1.0D0,C,",M,")"
       CASE(3)
         write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
               "      CALL DGEMM('N','T',",M,",",N,",",K,",1.0D0,A,",M,",B,",N,",1.0D0,C,",M,")"
       CASE(4)
         write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
               "      CALL DGEMM('T','T',",M,",",N,",",K,",1.0D0,A,",K,",B,",N,",1.0D0,C,",M,")"
       END SELECT
       write(6,'(A)') "   END SUBROUTINE"
     CASE(4,5,6,7)
       isquare=version-3
       ! generation of the multrec versions
       write(6,'(A,I0,A,I0,A,I0,A,A)')    "   SUBROUTINE smm_d"//trstr(transpose_flavor)//"_",M,"_",N,"_",K,TRIM(label),"(A,B,C)"
       CALL write_matrix_defs(M,N,K,transpose_flavor)
       write(6,'(A)')                     "      INTEGER ::i,j,l"
       blocksize=best_square(isquare)
       CALL MULTREC(1,M,1,N,1,K,blocksize,tiny_opts,transpose_flavor)
       write(6,'(A)') "   END SUBROUTINE"
     CASE DEFAULT
       STOP "MISSING CASE mult_versions"
     END SELECT
  END SUBROUTINE mult_versions

END MODULE multrec_gen
