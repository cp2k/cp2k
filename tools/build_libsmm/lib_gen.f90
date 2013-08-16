PROGRAM lib_gen
   USE mults
   USE multrec_gen
   IMPLICIT NONE
   
   INTEGER :: M,N,K
   CHARACTER(LEN=1024) :: arg,filename_small,filename_tiny,line,label
   REAL :: tmp
   INTEGER, DIMENSION(:,:), ALLOCATABLE ::  small_opts
   REAL, DIMENSION(:), ALLOCATABLE :: small_perf
   INTEGER :: opt,iline,nline,transpose_flavor,data_type, SIMD_size
   CHARACTER(LEN=10), PARAMETER :: stack_size_label = "stack_size"
   
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
   CALL GET_COMMAND_ARGUMENT(7,filename_small)
   CALL GET_COMMAND_ARGUMENT(8,filename_tiny)

   !
   ! filename is the result of small optimization (cat small_gen_optimal.out )
   ! 6 13 6    4    0.756046       6.613
   !
   OPEN(UNIT=10,FILE=filename_small)
   REWIND(10)
   nline=0
   DO
     READ(10,*,END=999) line
     nline=nline+1
   ENDDO
999  CONTINUE
   ALLOCATE(small_opts(4,nline))
   ALLOCATE(small_perf(nline))
   REWIND(10)
   DO iline=1,nline
     READ(10,'(A1024)',END=999) line
     READ(line,*) small_opts(:,iline),tmp,small_perf(iline)
   ENDDO
   CLOSE(10)

   CALL find_small_opt(opt,small_opts,m,n,k)

   label=""
   CALL mult_versions(M,N,K,opt,label,transpose_flavor,data_type,SIMD_size,filename_tiny,stack_size_label,.TRUE.)

CONTAINS

  SUBROUTINE find_small_opt(opt,small_opts,m,n,k)
     INTEGER, INTENT(OUT) :: opt
     INTEGER, DIMENSION(:,:) :: small_opts
     INTEGER :: m,n,k
     INTEGER :: i
     ! by default we call dgemm (but this should never happen)
     opt=3
     DO i=1,SIZE(small_opts,2)
        IF (ALL(small_opts(1:3,i)==(/m,n,k/))) opt=small_opts(4,i)
     ENDDO
  END SUBROUTINE find_small_opt
   
END PROGRAM lib_gen
