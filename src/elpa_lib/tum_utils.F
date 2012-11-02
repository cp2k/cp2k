MODULE tum_utils
 
    IMPLICIT NONE

    PRIVATE
#if defined (__parallel)
    INCLUDE "mpif.h"
 
    PUBLIC :: local_size_offset_1d
    PUBLIC :: reverse_vector_local
    PUBLIC :: reverse_matrix_local
    PUBLIC :: reverse_matrix_1dcomm
    PUBLIC :: reverse_matrix_2dcomm_ref

CONTAINS

! rev parameter is critical, even in rev only mode!
! pdgeqrf_2dcomm uses rev=0 version to determine the process columns 
! involved in the qr decomposition
SUBROUTINE local_size_offset_1d(n,nb,baseidx,idx,rev,rank,nprocs, &
                                lsize,baseoffset,offset)
    USE elpa1
    INTEGER                                  :: n, nb, baseidx, idx, rev, &
                                                rank, nprocs, lsize, &
                                                baseoffset, offset

    INTEGER                                  :: rank_idx

! input
! output
! local scalars

    rank_idx = MOD((idx-1)/nb,nprocs)

    ! calculate local size and offsets
    IF (rev .EQ. 1) THEN
        IF (idx > 0) THEN
            lsize = local_index(idx,rank,nprocs,nb,-1)
        ELSE 
            lsize = 0
        END IF

        baseoffset = 1
        offset = 1
    ELSE
        offset = local_index(idx,rank,nprocs,nb,1)
        baseoffset = local_index(baseidx,rank,nprocs,nb,1)

        lsize = local_index(n,rank,nprocs,nb,-1)
        !print *,'baseidx,idx',baseidx,idx,lsize,n

        lsize = lsize - offset + 1
 
        baseoffset = offset - baseoffset + 1
    END IF

END SUBROUTINE local_size_offset_1d


SUBROUTINE reverse_vector_local(n,x,incx,work,lwork)

    ! input
    INTEGER                                  :: n, incx, lwork

    DOUBLE PRECISION x(*),work(*)

    ! local scalars
    DOUBLE PRECISION temp
    INTEGER srcoffset,destoffset,ientry

    IF (lwork .EQ. -1) THEN
        work(1) = 0.0d0
        RETURN
    END IF
    
    DO ientry=1,n/2
        srcoffset=1+(ientry-1)*incx
        destoffset=1+(n-ientry)*incx
        
        temp = x(srcoffset)
        x(srcoffset) = x(destoffset)
        x(destoffset) = temp
    END DO

END SUBROUTINE reverse_vector_local

SUBROUTINE reverse_matrix_local(trans,m,n,a,lda,work,lwork)

    ! input
    INTEGER                                  :: trans, m, n, lda, lwork

    DOUBLE PRECISION a(lda,*),work(*)

    ! local scalars
    DOUBLE PRECISION dworksize(1)
    INTEGER incx
    INTEGER dimsize
    INTEGER i

    IF (trans .EQ. 1) THEN
        incx = lda
        dimsize = n
    ELSE
        incx = 1
        dimsize = m
    END IF

    IF (lwork .EQ. -1) THEN
        CALL reverse_vector_local(dimsize,a,incx,dworksize,-1)
        work(1) = dworksize(1)
        RETURN
    END IF

    IF (trans .EQ. 1) THEN
        DO i=1,m
            CALL reverse_vector_local(dimsize,a(i,1),incx,work,lwork)
        END DO
    ELSE
        DO i=1,n
            CALL reverse_vector_local(dimsize,a(1,i),incx,work,lwork)
        END DO
    END IF
    
END SUBROUTINE reverse_matrix_local

SUBROUTINE reverse_matrix_2dcomm_ref(m,n,mb,nb,a,lda,work,lwork,mpicomm_cols,mpicomm_rows)

    ! input
    INTEGER                                  :: m, n, mb, nb, lda, lwork, &
                                                mpicomm_cols, mpicomm_rows

    DOUBLE PRECISION a(lda,*),work(*)

    ! local scalars
    DOUBLE PRECISION reverse_column_size(1)
    DOUBLE PRECISION reverse_row_size(1)

    INTEGER mpirank_cols,mpirank_rows
    INTEGER mpiprocs_cols,mpiprocs_rows
    INTEGER mpierr
    INTEGER lrows,lcols,offset,baseoffset

    CALL MPI_Comm_rank(mpicomm_cols,mpirank_cols,mpierr)
    CALL MPI_Comm_rank(mpicomm_rows,mpirank_rows,mpierr)
    CALL MPI_Comm_size(mpicomm_cols,mpiprocs_cols,mpierr)
    CALL MPI_Comm_size(mpicomm_rows,mpiprocs_rows,mpierr)
 
    CALL local_size_offset_1d(m,mb,1,1,0,mpirank_cols,mpiprocs_cols, &
                                  lrows,baseoffset,offset)
  
    CALL local_size_offset_1d(n,nb,1,1,0,mpirank_rows,mpiprocs_rows, &
                                  lcols,baseoffset,offset)
 
    IF (lwork .EQ. -1) THEN
        CALL reverse_matrix_1dcomm(0,m,lcols,mb,a,lda,reverse_column_size,-1,mpicomm_cols)
        CALL reverse_matrix_1dcomm(1,lrows,n,nb,a,lda,reverse_row_size,-1,mpicomm_rows)
        work(1) = MAX(reverse_column_size(1),reverse_row_size(1))
        RETURN
    END IF
 
    CALL reverse_matrix_1dcomm(0,m,lcols,mb,a,lda,work,lwork,mpicomm_cols)
    CALL reverse_matrix_1dcomm(1,lrows,n,nb,a,lda,work,lwork,mpicomm_rows)
END SUBROUTINE reverse_matrix_2dcomm_ref

! b: if trans = 'N': b is size of block distribution between rows
! b: if trans = 'T': b is size of block distribution between columns
SUBROUTINE reverse_matrix_1dcomm(trans,m,n,b,a,lda,work,lwork,mpicomm)
    
   
    INTEGER                                  :: trans, m, n, b, lda, lwork, &
                                                mpicomm

! input

    DOUBLE PRECISION a(lda,*),work(*)
 
    ! local scalars
    INTEGER mpirank,mpiprocs,mpierr,mpistatus(MPI_STATUS_SIZE)
    INTEGER nr_blocks,dest_process,src_process,step
    INTEGER lsize,baseoffset,offset
    INTEGER current_index,destblk,srcblk,icol,next_index
    INTEGER sendcount,recvcount
    INTEGER sendoffset,recvoffset
    INTEGER newmatrix_offset,work_offset
    INTEGER lcols,lrows,lroffset,lcoffset,dimsize,fixedsize
    DOUBLE PRECISION dworksize(1)

    CALL MPI_Comm_rank(mpicomm, mpirank, mpierr)
    CALL MPI_Comm_size(mpicomm, mpiprocs, mpierr)
  
    IF (trans .EQ. 1) THEN
        CALL local_size_offset_1d(n,b,1,1,0,mpirank,mpiprocs, &
                                  lcols,baseoffset,lcoffset)
        lrows = m
    ELSE
        CALL local_size_offset_1d(m,b,1,1,0,mpirank,mpiprocs, &
                                  lrows,baseoffset,lroffset)
        lcols = n
    END IF
                          
    IF (lwork .EQ. -1) THEN
        CALL reverse_matrix_local(trans,lrows,lcols,a,MAX(lrows,lcols),dworksize,-1)
        work(1) = DBLE(3*lrows*lcols) + dworksize(1)
        RETURN
    END IF

    sendoffset = 1
    recvoffset = sendoffset + lrows*lcols
    newmatrix_offset = recvoffset + lrows*lcols
    work_offset = newmatrix_offset + lrows*lcols

    IF (trans .EQ. 1) THEN
        dimsize = n
        fixedsize = m
    ELSE
        dimsize = m
        fixedsize = n
    END IF

    IF (dimsize .LE. 1) THEN
        RETURN ! nothing to do
    END IF
 
    ! 1. adjust step size to remainder size
    nr_blocks = dimsize / b
    nr_blocks = nr_blocks * b
    step = dimsize - nr_blocks
    IF (step .EQ. 0) step = b

    ! 2. iterate over destination blocks starting with process 0
    current_index = 1
    DO WHILE (current_index .LE. dimsize)
        destblk = (current_index-1) / b
        dest_process = MOD(destblk,mpiprocs)
        srcblk = (dimsize-current_index) / b
        src_process = MOD(srcblk,mpiprocs)

        next_index = current_index+step

        ! block for dest_process is located on mpirank if lsize > 0
        CALL local_size_offset_1d(dimsize-current_index+1,b,dimsize-next_index+2,dimsize-next_index+2,0, &
                                  src_process,mpiprocs,lsize,baseoffset,offset)

        sendcount = lsize*fixedsize
        recvcount = sendcount

        ! TODO: this send/recv stuff seems to blow up on BlueGene/P 
        ! TODO: is there actually room for the requested matrix part? the target
        ! process might not have any parts at all (thus no room)
        IF ((src_process .EQ. mpirank) .AND. (dest_process .EQ. src_process)) THEN
                ! 5. pack data
                IF (trans .EQ. 1) THEN
                    DO icol=offset,offset+lsize-1
                        work(sendoffset+(icol-offset)*lrows:sendoffset+(icol-offset+1)*lrows-1) = &
                            a(1:lrows,icol)
                    END DO
                ELSE
                    DO icol=1,lcols
                        work(sendoffset+(icol-1)*lsize:sendoffset+icol*lsize-1) = &
                            a(offset:offset+lsize-1,icol)
                    END DO
                END IF
 
                ! 7. reverse data
                IF (trans .EQ. 1) THEN
                    CALL reverse_matrix_local(1,lrows,lsize,work(sendoffset),lrows,work(work_offset),lwork)
                ELSE
                    CALL reverse_matrix_local(0,lsize,lcols,work(sendoffset),lsize,work(work_offset),lwork)
                END IF

                ! 8. store in temp matrix
                IF (trans .EQ. 1) THEN
                    DO icol=1,lsize
                        work(newmatrix_offset+(icol-1)*lrows:newmatrix_offset+icol*lrows-1) = &
                            work(sendoffset+(icol-1)*lrows:sendoffset+icol*lrows-1)
                    END DO

                    newmatrix_offset = newmatrix_offset + lsize*lrows
                ELSE
                    DO icol=1,lcols
                        work(newmatrix_offset+(icol-1)*lrows:newmatrix_offset+(icol-1)*lrows+lsize-1) = &
                            work(sendoffset+(icol-1)*lsize:sendoffset+icol*lsize-1)
                    END DO

                    newmatrix_offset = newmatrix_offset + lsize
                END IF
        ELSE

            IF (dest_process .EQ. mpirank) THEN
                ! 6b. call MPI_Recv
                CALL MPI_Recv(work(recvoffset), recvcount, mpi_real8, &
                              src_process, current_index, mpicomm, mpistatus, mpierr)

                ! 7. reverse data
                IF (trans .EQ. 1) THEN
                    CALL reverse_matrix_local(1,lrows,lsize,work(recvoffset),lrows,work(work_offset),lwork)
                ELSE
                    CALL reverse_matrix_local(0,lsize,lcols,work(recvoffset),lsize,work(work_offset),lwork)
                END IF

                ! 8. store in temp matrix
                IF (trans .EQ. 1) THEN
                    DO icol=1,lsize
                        work(newmatrix_offset+(icol-1)*lrows:newmatrix_offset+icol*lrows-1) = &
                            work(recvoffset+(icol-1)*lrows:recvoffset+icol*lrows-1)
                    END DO

                    newmatrix_offset = newmatrix_offset + lsize*lrows
                ELSE
                    DO icol=1,lcols
                        work(newmatrix_offset+(icol-1)*lrows:newmatrix_offset+(icol-1)*lrows+lsize-1) = &
                            work(recvoffset+(icol-1)*lsize:recvoffset+icol*lsize-1)
                    END DO

                    newmatrix_offset = newmatrix_offset + lsize
                END IF
            END IF

            IF (src_process .EQ. mpirank) THEN
                ! 5. pack data
                IF (trans .EQ. 1) THEN
                    DO icol=offset,offset+lsize-1
                        work(sendoffset+(icol-offset)*lrows:sendoffset+(icol-offset+1)*lrows-1) = &
                            a(1:lrows,icol)
                    END DO
                ELSE
                    DO icol=1,lcols
                        work(sendoffset+(icol-1)*lsize:sendoffset+icol*lsize-1) = &
                            a(offset:offset+lsize-1,icol)
                    END DO
                END IF

                ! 6a. call MPI_Send
                CALL MPI_Send(work(sendoffset), sendcount, mpi_real8, &
                                  dest_process, current_index, mpicomm, mpierr)
            END IF
        END IF

        current_index = next_index
    END DO

   ! 9. copy temp matrix to real matrix
   newmatrix_offset = recvoffset + lrows*lcols
   DO icol=1,lcols
        a(1:lrows,icol) = &
            work(newmatrix_offset+(icol-1)*lrows:newmatrix_offset+icol*lrows-1)
   END DO

END SUBROUTINE reverse_matrix_1dcomm

#endif  
END MODULE
