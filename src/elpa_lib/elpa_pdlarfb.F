MODULE elpa_pdlarfb

  USE elpa1
  USE tum_utils

    IMPLICIT NONE
#if defined (__parallel) 
    INCLUDE 'mpif.h'

    PRIVATE

    PUBLIC :: tum_pdlarfb_1dcomm
    PUBLIC :: tum_pdlarft_pdlarfb_1dcomm
    PUBLIC :: tum_pdlarft_set_merge_1dcomm
    PUBLIC :: tum_pdlarft_tree_merge_1dcomm
    PUBLIC :: tum_pdlarfl_1dcomm
    PUBLIC :: tum_pdlarfl2_tmatrix_1dcomm
    PUBLIC :: tum_tmerge_pdlarfb_1dcomm


CONTAINS

SUBROUTINE tum_pdlarfb_1dcomm(m,mb,n,k,a,lda,v,ldv,tau,t,ldt,baseidx,idx,rev,mpicomm,work,lwork)
    
    USE tum_utils
    INTEGER                                  :: lda, ldv, ldt, lwork

! input variables (local)

    INTEGER m,mb,n,k,baseidx,idx,rev,mpicomm
    DOUBLE PRECISION a(lda,*),v(ldv,*),tau(*),t(ldt,*),work(k,*)

    ! input variables (global)
 
    ! output variables (global)

    ! derived input variables from TUM_PQRPARAM

    ! local scalars
    INTEGER localsize,offset,baseoffset
    INTEGER mpirank,mpiprocs,mpierr

        IF (idx .LE. 1) RETURN

    IF (n .LE. 0) RETURN ! nothing to do

    IF (k .EQ. 1) THEN
        CALL tum_pdlarfl_1dcomm(v,1,baseidx,a,lda,tau(1), &
                                work,lwork,m,n,idx,mb,rev,mpicomm)
        RETURN
    ELSE IF (k .EQ. 2) THEN
        CALL tum_pdlarfl2_tmatrix_1dcomm(v,ldv,baseidx,a,lda,t,ldt, &
                                 work,lwork,m,n,idx,mb,rev,mpicomm)
        RETURN
    END IF

    IF (lwork .EQ. -1) THEN
        work(1,1) = DBLE(2*k*n)
        RETURN
    END IF
 
    !print *,'updating trailing matrix with k=',k

    CALL MPI_Comm_rank(mpicomm,mpirank,mpierr)
    CALL MPI_Comm_size(mpicomm,mpiprocs,mpierr)

    ! use baseidx as idx here, otherwise the upper triangle part will be lost
    ! during the calculation, especially in the reversed case
    CALL local_size_offset_1d(m,mb,baseidx,baseidx,rev,mpirank,mpiprocs, &
                                localsize,baseoffset,offset)

    ! Z' = Y' * A
    IF (localsize .GT. 0) THEN
        CALL dgemm("Trans","Notrans",k,n,localsize,1.0d0,v(baseoffset,1),ldv,a(offset,1),lda,0.0d0,work(1,1),k)
    ELSE
        work(1:k,1:n) = 0.0d0
    END IF

    ! data exchange
    CALL mpi_allreduce(work(1,1),work(1,n+1),k*n,mpi_real8,mpi_sum,mpicomm,mpierr)
    
    CALL tum_pdlarfb_kernel_local(localsize,n,k,a(offset,1),lda,v(baseoffset,1),ldv,t,ldt,work(1,n+1),k)
END SUBROUTINE tum_pdlarfb_1dcomm 

! generalized pdlarfl2 version
! TODO: include T merge here (seperate by "old" and "new" index)
SUBROUTINE tum_pdlarft_pdlarfb_1dcomm(m,mb,n,oldk,k,v,ldv,tau,t,ldt,a,lda,baseidx,rev,mpicomm,work,lwork)
    USE tum_utils
    INTEGER                                  :: ldv, ldt, lda, lwork

! input variables (local)

    INTEGER m,mb,n,k,oldk,baseidx,rev,mpicomm
    DOUBLE PRECISION v(ldv,*),tau(*),t(ldt,*),work(k,*),a(lda,*)

    ! input variables (global)
 
    ! output variables (global)

    ! derived input variables from TUM_PQRPARAM

    ! local scalars
    INTEGER localsize,offset,baseoffset
    INTEGER mpirank,mpiprocs,mpierr
    INTEGER icol

    INTEGER sendoffset,recvoffset,sendsize

    sendoffset = 1
    sendsize = k*(k+n+oldk)
    recvoffset = sendoffset+(k+n+oldk)

    IF (lwork .EQ. -1) THEN
        work(1,1) = DBLE(2*(k*k+k*n+oldk))
        RETURN
    END IF

    CALL MPI_Comm_rank(mpicomm,mpirank,mpierr)
    CALL MPI_Comm_size(mpicomm,mpiprocs,mpierr)

    CALL local_size_offset_1d(m,mb,baseidx,baseidx,rev,mpirank,mpiprocs, &
                                localsize,baseoffset,offset)

    IF (localsize .GT. 0) THEN
            ! calculate inner product of householdervectors
            CALL dsyrk("Upper","Trans",k,localsize,1.0d0,v(baseoffset,1),ldv,0.0d0,work(1,1),k)

            ! calculate matrix matrix product of householder vectors and target matrix 
            ! Z' = Y' * A
            CALL dgemm("Trans","Notrans",k,n,localsize,1.0d0,v(baseoffset,1),ldv,a(offset,1),lda,0.0d0,work(1,k+1),k)

            ! TODO: reserved for T merge parts
            work(1:k,n+k+1:n+k+oldk) = 0.0d0
    ELSE
        work(1:k,1:(n+k+oldk)) = 0.0d0
    END IF

    ! exchange data
    CALL mpi_allreduce(work(1,sendoffset),work(1,recvoffset),sendsize,mpi_real8,mpi_sum,mpicomm,mpierr)

        ! generate T matrix (pdlarft)
        t(1:k,1:k) = 0.0d0 ! DEBUG: clear buffer first

        ! T1 = tau1
        ! | tauk  Tk-1' * (-tauk * Y(:,1,k+1:n) * Y(:,k))' |
        ! | 0           Tk-1                           |
        t(k,k) = tau(k)
        DO icol=k-1,1,-1
            t(icol,icol+1:k) = -tau(icol)*work(icol,recvoffset+icol:recvoffset+k-1)
            CALL dtrmv("Upper","Trans","Nonunit",k-icol,t(icol+1,icol+1),ldt,t(icol,icol+1),ldt)
            t(icol,icol) = tau(icol)
        END DO

        ! TODO: elmroth and gustavson
 
        ! update matrix (pdlarfb)
        ! Z' = T * Z'
        CALL dtrmm("Left","Upper","Notrans","Nonunit",k,n,1.0d0,t,ldt,work(1,recvoffset+k),k)

        ! A = A - Y * V'
        CALL dgemm("Notrans","Notrans",localsize,n,k,-1.0d0,v(baseoffset,1),ldv,work(1,recvoffset+k),k,1.0d0,a(offset,1),lda)

END SUBROUTINE tum_pdlarft_pdlarfb_1dcomm

SUBROUTINE tum_pdlarft_set_merge_1dcomm(m,mb,n,blocksize,v,ldv,t,ldt,baseidx,rev,mpicomm,work,lwork)
    USE tum_utils
    INTEGER                                  :: ldv, ldt, lwork

! input variables (local)

    INTEGER m,mb,n,blocksize,baseidx,rev,mpicomm
    DOUBLE PRECISION v(ldv,*),t(ldt,*),work(n,*)

    ! input variables (global)
 
    ! output variables (global)

    ! derived input variables from TUM_PQRPARAM

    ! local scalars
    INTEGER localsize,offset,baseoffset
    INTEGER mpirank,mpiprocs,mpierr

    IF (lwork .EQ. -1) THEN
        work(1,1) = DBLE(2*n*n)
        RETURN
    END IF
 
    CALL MPI_Comm_rank(mpicomm,mpirank,mpierr)
    CALL MPI_Comm_size(mpicomm,mpiprocs,mpierr)

    CALL local_size_offset_1d(m,mb,baseidx,baseidx,rev,mpirank,mpiprocs, &
                                localsize,baseoffset,offset)

    IF (localsize .GT. 0) THEN
        CALL dsyrk("Upper","Trans",n,localsize,1.0d0,v(baseoffset,1),ldv,0.0d0,work(1,1),n)
    ELSE
        work(1:n,1:n) = 0.0d0
    END IF
 
    CALL mpi_allreduce(work(1,1),work(1,n+1),n*n,mpi_real8,mpi_sum,mpicomm,mpierr)

        ! skip Y4'*Y4 part
        offset = MOD(n,blocksize)
        IF (offset .EQ. 0) offset=blocksize
        CALL tum_tmerge_set_kernel(n,blocksize,t,ldt,work(1,n+1+offset),n)

END SUBROUTINE tum_pdlarft_set_merge_1dcomm

SUBROUTINE tum_pdlarft_tree_merge_1dcomm(m,mb,n,blocksize,treeorder,v,ldv,t,ldt,baseidx,rev,mpicomm,work,lwork)
    USE tum_utils
    INTEGER                                  :: ldv, ldt, lwork

! input variables (local)

    INTEGER m,mb,n,blocksize,treeorder,baseidx,rev,mpicomm
    DOUBLE PRECISION v(ldv,*),t(ldt,*),work(n,*)

    ! input variables (global)
 
    ! output variables (global)

    ! derived input variables from TUM_PQRPARAM

    ! local scalars
    INTEGER localsize,offset,baseoffset
    INTEGER mpirank,mpiprocs,mpierr

    IF (lwork .EQ. -1) THEN
        work(1,1) = DBLE(2*n*n)
        RETURN
    END IF

    IF (n .LE. blocksize) RETURN ! nothing to do
 
    CALL MPI_Comm_rank(mpicomm,mpirank,mpierr)
    CALL MPI_Comm_size(mpicomm,mpiprocs,mpierr)

    CALL local_size_offset_1d(m,mb,baseidx,baseidx,rev,mpirank,mpiprocs, &
                                localsize,baseoffset,offset)

    IF (localsize .GT. 0) THEN
        CALL dsyrk("Upper","Trans",n,localsize,1.0d0,v(baseoffset,1),ldv,0.0d0,work(1,1),n)
    ELSE
        work(1:n,1:n) = 0.0d0
    END IF
 
    CALL mpi_allreduce(work(1,1),work(1,n+1),n*n,mpi_real8,mpi_sum,mpicomm,mpierr)

        ! skip Y4'*Y4 part
        offset = MOD(n,blocksize)
        IF (offset .EQ. 0) offset=blocksize
        CALL tum_tmerge_tree_kernel(n,blocksize,treeorder,t,ldt,work(1,n+1+offset),n)

END SUBROUTINE tum_pdlarft_tree_merge_1dcomm

! apply householder vector to the left 
! - assume unitary matrix
! - assume right positions for v
SUBROUTINE tum_pdlarfl_1dcomm(v,incv,baseidx,a,lda,tau,work,lwork,m,n,idx,mb,rev,mpicomm)
    USE ELPA1
    USE tum_utils
    INTEGER                                  :: incv, baseidx, lda, lwork

! input variables (local)

    DOUBLE PRECISION v(*),a(lda,*),work(*)

    ! input variables (global)
    INTEGER m,n,mb,rev,idx,mpicomm
    DOUBLE PRECISION tau
 
    ! output variables (global)
 
    ! local scalars
    INTEGER mpierr,mpirank,mpiprocs
    INTEGER sendsize,recvsize,icol
    INTEGER local_size,local_offset
    INTEGER v_local_offset

    ! external functions
    DOUBLE PRECISION ddot
    EXTERNAL dgemv,dger,ddot
  
    CALL MPI_Comm_rank(mpicomm, mpirank, mpierr)
    CALL MPI_Comm_size(mpicomm, mpiprocs, mpierr)

    sendsize = n
    recvsize = sendsize

    IF (lwork .EQ. -1) THEN
        work(1) = DBLE(sendsize + recvsize)
        RETURN
    END IF
 
    IF (n .LE. 0) RETURN

        IF (idx .LE. 1) RETURN
 
    CALL local_size_offset_1d(m,mb,baseidx,idx,rev,mpirank,mpiprocs, &
                              local_size,v_local_offset,local_offset)
 
    !print *,'hl ref',local_size,n

    v_local_offset = v_local_offset * incv

    IF (local_size > 0) THEN
        
        DO icol=1,n
            work(icol) = DOT_PRODUCT(v(v_local_offset:v_local_offset+local_size-1),a(local_offset:local_offset+local_size-1,icol))

        END DO
    ELSE
        work(1:n) = 0.0d0
    END IF
 
    CALL mpi_allreduce(work, work(sendsize+1), sendsize, mpi_real8, mpi_sum, mpicomm, mpierr)

    IF (local_size > 0) THEN

         DO icol=1,n
             a(local_offset:local_offset+local_size-1,icol) = &
             a(local_offset:local_offset+local_size-1,icol) &
             - tau*work(sendsize+icol)*v(v_local_offset:v_local_offset &
             +local_size-1)
         ENDDO
    END IF

END SUBROUTINE tum_pdlarfl_1dcomm

SUBROUTINE tum_pdlarfl2_tmatrix_1dcomm(v,ldv,baseidx,a,lda,t,ldt,work,lwork,m,n,idx,mb,rev,mpicomm)
    USE ELPA1
    USE tum_utils
    INTEGER                                  :: ldv, baseidx, lda, ldt, lwork

! input variables (local)

    DOUBLE PRECISION v(ldv,*),a(lda,*),work(*),t(ldt,*)

    ! input variables (global)
    INTEGER m,n,mb,rev,idx,mpicomm
 
    ! output variables (global)
 
    ! local scalars
    INTEGER mpierr,mpirank,mpiprocs,mpirank_top1,mpirank_top2
    INTEGER dgemv1_offset,dgemv2_offset
    INTEGER sendsize, recvsize
    INTEGER local_size1,local_offset1
    INTEGER local_size2,local_offset2
    INTEGER local_size_dger,local_offset_dger
    INTEGER v1_local_offset,v2_local_offset
    INTEGER v_local_offset_dger
    DOUBLE PRECISION hvdot
    INTEGER irow,icol,v1col,v2col

    ! external functions
    DOUBLE PRECISION ddot
    EXTERNAL dgemv,dger,ddot,daxpy

    CALL MPI_Comm_rank(mpicomm, mpirank, mpierr)
    CALL MPI_Comm_size(mpicomm, mpiprocs, mpierr)
 
    sendsize = 2*n
    recvsize = sendsize

    IF (lwork .EQ. -1) THEN
        work(1) = sendsize + recvsize
        RETURN
    END IF
 
    dgemv1_offset = 1
    dgemv2_offset = dgemv1_offset + n

        ! in 2x2 matrix case only one householder vector was generated
        IF (idx .LE. 2) THEN
            CALL tum_pdlarfl_1dcomm(v(1,2),1,baseidx,a,lda,t(2,2), &
                                    work,lwork,m,n,idx,mb,rev,mpicomm)
            RETURN
        END IF

        CALL local_size_offset_1d(m,mb,baseidx,idx,rev,mpirank,mpiprocs, &
                                  local_size1,v1_local_offset,local_offset1)
        CALL local_size_offset_1d(m,mb,baseidx,idx-1,rev,mpirank,mpiprocs, &
                                  local_size2,v2_local_offset,local_offset2)

        v1_local_offset = v1_local_offset * 1
        v2_local_offset = v2_local_offset * 1

        v1col = 2
        v2col = 1

        ! keep buffers clean in case that local_size1/local_size2 are zero
        work(1:sendsize) = 0.0d0

        CALL dgemv("Trans",local_size1,n,1.0d0,a(local_offset1,1),lda,&
        v(v1_local_offset,v1col),1,0.0d0,work(dgemv1_offset),1)
        CALL dgemv("Trans",local_size2,n,t(v2col,v2col),a(local_offset2,1),lda,&
        v(v2_local_offset,v2col),1,0.0d0,work(dgemv2_offset),1)

        CALL mpi_allreduce(work, work(sendsize+1), sendsize, mpi_real8, &
        mpi_sum, mpicomm, mpierr)
  
        ! update second vector
        CALL daxpy(n,t(1,2),work(sendsize+dgemv1_offset),1, &
        work(sendsize+dgemv2_offset),1)

        CALL local_size_offset_1d(m,mb,baseidx,idx-2,rev,mpirank,mpiprocs, &
             local_size_dger,v_local_offset_dger,local_offset_dger)

        ! get ranks of processes with topelements
        mpirank_top1 = MOD((idx-1)/mb,mpiprocs)
        mpirank_top2 = MOD((idx-2)/mb,mpiprocs)

        IF (mpirank_top1 .EQ. mpirank) local_offset1 = local_size1
        IF (mpirank_top2 .EQ. mpirank) THEN
            local_offset2 = local_size2
            v2_local_offset = local_size2
        END IF

    ! use hvdot as temporary variable
    hvdot = t(v1col,v1col)
    DO icol=1,n
        ! make use of "1" entries in householder vectors
        IF (mpirank_top1 .EQ. mpirank) THEN
            a(local_offset1,icol) = a(local_offset1,icol) &
                                    - work(sendsize+dgemv1_offset+icol-1)*hvdot
        END IF

        IF (mpirank_top2 .EQ. mpirank) THEN
            a(local_offset2,icol) = a(local_offset2,icol)- & 
            v(v2_local_offset,v1col) * work(sendsize+dgemv1_offset+icol-1) &
             *hvdot - work(sendsize+dgemv2_offset+icol-1)
        END IF

        DO irow=1,local_size_dger
          a(local_offset_dger+irow-1,icol) = a(local_offset_dger+irow-1,icol) &
          - work(sendsize+dgemv1_offset+icol-1)*v(v_local_offset_dger+irow-1, &
            v1col)*hvdot - work(sendsize+dgemv2_offset+icol-1)* & 
            v(v_local_offset_dger+irow-1,v2col) 
        END DO
    END DO

END SUBROUTINE tum_pdlarfl2_tmatrix_1dcomm

! generalized pdlarfl2 version
! TODO: include T merge here (seperate by "old" and "new" index)
SUBROUTINE tum_tmerge_pdlarfb_1dcomm(m,mb,n,oldk,k,v,ldv,t,ldt,a,lda,baseidx,rev,updatemode,mpicomm,work,lwork)
    USE tum_utils
    INTEGER                                  :: ldv, ldt, lda, lwork

! input variables (local)

    DOUBLE PRECISION v(ldv,*),t(ldt,*),work(*),a(lda,*)

    ! input variables (global)
    INTEGER m,mb,n,k,oldk,baseidx,rev,updatemode,mpicomm
 
    ! output variables (global)

    ! derived input variables from TUM_PQRPARAM

    ! local scalars
    INTEGER localsize,offset,baseoffset
    INTEGER mpirank,mpiprocs,mpierr

    INTEGER sendoffset,recvoffset,sendsize
    INTEGER updateoffset,updatelda,updatesize
    INTEGER mergeoffset,mergelda,mergesize
    INTEGER tgenoffset,tgenlda,tgensize

        IF (updatemode .EQ. ICHAR('I')) THEN
            updatelda = oldk+k
        ELSE
            updatelda = k
        END IF

        updatesize = updatelda*n

        mergelda = k
        mergesize = mergelda*oldk

        tgenlda = 0
        tgensize = 0

        sendsize = updatesize + mergesize + tgensize

    IF (lwork .EQ. -1) THEN
        work(1) = DBLE(2*sendsize)
        RETURN
    END IF

    CALL MPI_Comm_rank(mpicomm,mpirank,mpierr)
    CALL MPI_Comm_size(mpicomm,mpiprocs,mpierr)
 
    ! use baseidx as idx here, otherwise the upper triangle part will be lost
    ! during the calculation, especially in the reversed case
    CALL local_size_offset_1d(m,mb,baseidx,baseidx,rev,mpirank,mpiprocs, &
                                localsize,baseoffset,offset)

    sendoffset = 1

        IF (oldk .GT. 0) THEN
            updateoffset = 0
            mergeoffset = updateoffset + updatesize
            tgenoffset = mergeoffset + mergesize
            
            sendsize = updatesize + mergesize + tgensize

            !print *,'sendsize',sendsize,updatesize,mergesize,tgensize
            !print *,'merging nr of rotations', oldk+k
 
            IF (localsize .GT. 0) THEN
                ! calculate matrix matrix product of householder vectors and target matrix 

                IF (updatemode .EQ. ICHAR('I')) THEN
                    ! Z' = (Y1,Y2)' * A
                    CALL dgemm("Trans","Notrans",k+oldk,n,localsize,1.0d0, &
                    v(baseoffset,1),ldv,a(offset,1),lda,0.0d0, &
                    work(sendoffset+updateoffset),updatelda)
                ELSE
                    ! Z' = Y1' * A
                    CALL dgemm("Trans","Notrans",k,n,localsize,1.0d0, &
                    v(baseoffset,1),ldv,a(offset,1),lda,0.0d0, &
                    work(sendoffset+updateoffset),updatelda) 
                END IF

                ! calculate parts needed for T merge
                CALL dgemm("Trans","Notrans",k,oldk,localsize,1.0d0, &
                v(baseoffset,1),ldv,v(baseoffset,k+1),ldv,0.0d0, &
                work(sendoffset+mergeoffset),mergelda) 

            ELSE
                ! cleanup buffer
                work(sendoffset:sendoffset+sendsize-1) = 0.0d0
            END IF
        ELSE
            ! do not calculate parts for T merge as there is nothing to merge

            updateoffset = 0
            
            tgenoffset = updateoffset + updatesize
            
            sendsize = updatesize + tgensize
 
            IF (localsize .GT. 0) THEN
                ! calculate matrix matrix product of householder vectors and target matrix 
                ! Z' = (Y1)' * A
                CALL dgemm("Trans","Notrans",k,n,localsize,1.0d0, &
                v(baseoffset,1),ldv,a(offset,1),lda,0.0d0, &
                work(sendoffset+updateoffset),updatelda) 

            ELSE
                ! cleanup buffer
                work(sendoffset:sendoffset+sendsize-1) = 0.0d0
            END IF

        END IF

    recvoffset = sendoffset + sendsize

    IF (sendsize .LE. 0) RETURN ! nothing to do

    ! exchange data
    CALL mpi_allreduce(work(sendoffset),work(recvoffset),sendsize, &
    mpi_real8,mpi_sum,mpicomm,mpierr)
 
    updateoffset = recvoffset+updateoffset
    mergeoffset = recvoffset+mergeoffset
    tgenoffset = recvoffset+tgenoffset

        IF (oldk .GT. 0) THEN
            CALL tum_pdlarft_merge_kernel_local(oldk,k,t,ldt,work(mergeoffset),mergelda)

            IF (localsize .GT. 0) THEN
                IF (updatemode .EQ. ICHAR('I')) THEN

                    ! update matrix (pdlarfb) with complete T
                    CALL tum_pdlarfb_kernel_local(localsize,n,k+oldk, &
                    a(offset,1),lda,v(baseoffset,1),ldv,t(1,1),ldt, &
                    work(updateoffset),updatelda) 
                ELSE
                    ! update matrix (pdlarfb) with small T (same as update with no old T TODO)
                    CALL tum_pdlarfb_kernel_local(localsize,n,k,a(offset,1), &
                    lda,v(baseoffset,1),ldv,t(1,1),ldt,work(updateoffset), &
                    updatelda) 
                END IF
            END IF
        ELSE
            IF (localsize .GT. 0) THEN
                ! update matrix (pdlarfb) with small T
                CALL tum_pdlarfb_kernel_local(localsize,n,k,a(offset,1), &
                lda,v(baseoffset,1),ldv,t(1,1),ldt,work(updateoffset), &
                updatelda)
            END IF
        END IF

END SUBROUTINE tum_tmerge_pdlarfb_1dcomm

#endif    
END MODULE elpa_pdlarfb
