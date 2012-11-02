MODULE elpa_pdgeqrf

  USE elpa1
  USE elpa_pdlarfb
  USE tum_utils

    IMPLICIT NONE

    PRIVATE
#if defined (__parallel)
    INCLUDE "mpif.h"

    PUBLIC :: tum_pdgeqrf_2dcomm
    PUBLIC :: tum_pqrparam_init
    PUBLIC :: tum_pdlarfg2_1dcomm_check

CONTAINS

SUBROUTINE tum_pdgeqrf_2dcomm(a,lda,v,ldv,tau,t,ldt,work,lwork,m,n,mb,nb, &
            rowidx,colidx,rev,trans,PQRPARAM,mpicomm_rows,mpicomm_cols, &
            blockheuristic)
    INTEGER                                  :: eps_, gmode_, rank_

! parameter setup

    PARAMETER   (gmode_ = 1,rank_ = 2,eps_=3)

    ! input variables (local)
    INTEGER lda,lwork,ldv,ldt
    DOUBLE PRECISION a(lda,*),v(ldv,*),tau(*),work(*),t(ldt,*)

    ! input variables (global)
    INTEGER m,n,mb,nb,rowidx,colidx,rev,trans,mpicomm_cols,mpicomm_rows
    INTEGER PQRPARAM(*)
 
    ! output variables (global)
    DOUBLE PRECISION blockheuristic(*)

    ! input variables derived from PQRPARAM
    INTEGER updatemode,tmerge,size2d
    ! local scalars
    INTEGER mpierr,mpirank_cols,broadcast_size,mpirank_rows
    INTEGER mpirank_cols_qr,mpiprocs_cols
    INTEGER lcols_temp,lcols,icol,lastcol
    INTEGER baseoffset,offset,idx,voffset
    INTEGER update_voffset,update_tauoffset
    INTEGER update_lcols
    INTEGER work_offset

    DOUBLE PRECISION dbroadcast_size(1),dtmat_bcast_size(1)
    DOUBLE PRECISION pdgeqrf_size(1),pdlarft_size(1),pdlarfb_size(1), &
    tmerge_pdlarfb_size(1)
    INTEGER temptau_offset,temptau_size,broadcast_offset,tmat_bcast_size
    INTEGER remaining_cols
    INTEGER total_cols
    INTEGER incremental_update_size ! needed for incremental update mode
 
    size2d = PQRPARAM(1)
    updatemode = PQRPARAM(2)
    tmerge = PQRPARAM(3)

    ! copy value before we are going to filter it
    total_cols = n

    CALL mpi_comm_rank(mpicomm_cols,mpirank_cols,mpierr)
    CALL mpi_comm_rank(mpicomm_rows,mpirank_rows,mpierr)
    CALL mpi_comm_size(mpicomm_cols,mpiprocs_cols,mpierr)
  
   
    CALL tum_pdgeqrf_1dcomm(a,lda,v,ldv,tau,t,ldt,pdgeqrf_size(1),-1,m, & 
    total_cols,mb,rowidx,rowidx,rev,trans,PQRPARAM(4),mpicomm_rows, &
    blockheuristic)
    CALL tum_pdgeqrf_pack_unpack(v,ldv,dbroadcast_size(1),-1,m,total_cols, &
    mb,rowidx,rowidx,rev,0,mpicomm_rows)
    CALL tum_pdgeqrf_pack_unpack_tmatrix(tau,t,ldt,dtmat_bcast_size(1),-1, &
    total_cols,0)
    pdlarft_size(1) = 0.0d0
    CALL tum_pdlarfb_1dcomm(m,mb,total_cols,total_cols,a,lda,v,ldv,tau,t,&
    ldt,rowidx,rowidx,rev,mpicomm_rows,pdlarfb_size(1),-1)
    CALL tum_tmerge_pdlarfb_1dcomm(m,mb,total_cols,total_cols,total_cols,v, &
    ldv,t,ldt,a,lda,rowidx,rev,updatemode,mpicomm_rows, &
    tmerge_pdlarfb_size(1),-1)


    temptau_offset = 1
    temptau_size = total_cols
    broadcast_offset = temptau_offset + temptau_size 
    broadcast_size = dbroadcast_size(1) + dtmat_bcast_size(1)
    work_offset = broadcast_offset + broadcast_size

    IF (lwork .EQ. -1) THEN
        work(1) = (DBLE(temptau_size) + DBLE(broadcast_size) + &
        MAX(pdgeqrf_size(1),pdlarft_size(1),pdlarfb_size(1), &
        tmerge_pdlarfb_size(1)))
        RETURN
    END IF

    lastcol = colidx-total_cols+1
    voffset = total_cols
  
    incremental_update_size = 0
 
    ! clear v buffer: just ensure that there is no junk in the upper triangle
    ! part, otherwise pdlarfb gets some problems
    ! pdlarfl(2) do not have these problems as they are working more on a vector
    ! basis
    v(1:ldv,1:total_cols) = 0.0d0
 
    icol = colidx

    remaining_cols = total_cols

    !print *,'start decomposition',m,rowidx,colidx
 
    DO WHILE (remaining_cols .GT. 0)

        ! determine rank of process column with next qr block
        mpirank_cols_qr = MOD((icol-1)/nb,mpiprocs_cols)

        ! lcols can't be larger than than nb
        ! exception: there is only one process column
 
        ! however, we might not start at the first local column.
        ! therefore assume a matrix of size (1xlcols) starting at (1,icol)
        ! determine the real amount of local columns
        lcols_temp = MIN(nb,(icol-lastcol+1))

        ! blocking parameter
        lcols_temp = MAX(MIN(lcols_temp,size2d),1)

        ! determine size from last decomposition column
        !  to first decomposition column
        CALL local_size_offset_1d(icol,nb,icol-lcols_temp+1,icol-lcols_temp+1, &
             0, mpirank_cols_qr,mpiprocs_cols, lcols,baseoffset,offset)
 
        voffset = remaining_cols - lcols + 1

        idx = rowidx - colidx + icol

        IF (mpirank_cols .EQ. mpirank_cols_qr) THEN 
            ! qr decomposition part

            tau(offset:offset+lcols-1) = 0.0d0

            CALL tum_pdgeqrf_1dcomm(a(1,offset),lda,v(1,voffset),ldv, &
             tau(offset),t(voffset,voffset),ldt,work(work_offset),lwork, &
             m,lcols,mb,rowidx,idx,rev,trans,PQRPARAM(4),mpicomm_rows, &
             blockheuristic)
            !print *,'offset voffset',offset,voffset,idx
    
            ! pack broadcast buffer (v + tau)
            CALL tum_pdgeqrf_pack_unpack(v(1,voffset),ldv, &
            work(broadcast_offset),lwork,m,lcols,mb,rowidx,&
            idx,rev,0,mpicomm_rows)
     
            ! determine broadcast size
            CALL tum_pdgeqrf_pack_unpack(v(1,voffset),ldv, &
            dbroadcast_size(1),-1,m,lcols,mb,rowidx,idx,rev,0,&
            mpicomm_rows)
            broadcast_size = dbroadcast_size(1)
  
            !if (mpirank_rows .eq. 0) then
            ! pack tmatrix into broadcast buffer and calculate new size
            CALL tum_pdgeqrf_pack_unpack_tmatrix(tau(offset), &
            t(voffset,voffset),ldt,work(broadcast_offset+ &
            broadcast_size),lwork,lcols,0)
            CALL tum_pdgeqrf_pack_unpack_tmatrix(tau(offset),&
            t(voffset,voffset),ldt,dtmat_bcast_size(1),-1,lcols,0)
            broadcast_size = broadcast_size + dtmat_bcast_size(1)
            !end if
 
            ! initiate broadcast (send part)
            CALL MPI_Bcast(work(broadcast_offset),broadcast_size,mpi_real8, &
                           mpirank_cols_qr,mpicomm_cols,mpierr)

            ! copy tau parts into temporary tau buffer
            work(temptau_offset+voffset-1:temptau_offset+(voffset-1)+lcols-1)= &
            tau(offset:offset+lcols-1)

            !print *,'generated tau:', tau(offset)
        ELSE
            ! vector exchange part

            ! determine broadcast size
            CALL tum_pdgeqrf_pack_unpack(v(1,voffset),ldv,dbroadcast_size(1), &
            -1,m,lcols,mb,rowidx,idx,rev,1,mpicomm_rows)
            broadcast_size = dbroadcast_size(1)
 
            CALL tum_pdgeqrf_pack_unpack_tmatrix(work(temptau_offset+ & 
            voffset-1),t(voffset,voffset),ldt,dtmat_bcast_size(1),-1,lcols,0) 
            tmat_bcast_size = dtmat_bcast_size(1)

            !print *,'broadcast_size (nonqr)',broadcast_size
            broadcast_size = dbroadcast_size(1) + dtmat_bcast_size(1)
 
            ! initiate broadcast (recv part)
            CALL MPI_Bcast(work(broadcast_offset),broadcast_size,mpi_real8, &
                           mpirank_cols_qr,mpicomm_cols,mpierr)
            
            ! last n*n elements in buffer are (still empty) T matrix elements
            ! fetch from first process in each column
 
            ! unpack broadcast buffer (v + tau)
            CALL tum_pdgeqrf_pack_unpack(v(1,voffset),ldv, &
            work(broadcast_offset),lwork,m,lcols,mb,rowidx,idx, &
            rev,1,mpicomm_rows)
 
            ! now send t matrix to other processes in our process column
            broadcast_size = dbroadcast_size(1)
            tmat_bcast_size = dtmat_bcast_size(1)

            ! t matrix should now be available on all processes => unpack
            CALL tum_pdgeqrf_pack_unpack_tmatrix(work(temptau_offset+ &
            voffset-1),t(voffset,voffset),ldt,work(broadcast_offset+ &
            broadcast_size),lwork,lcols,1)
        END IF

        remaining_cols = remaining_cols - lcols
 
        ! apply householder vectors to whole trailing matrix parts (if any)

        update_voffset = voffset
        update_tauoffset = icol
        update_lcols = lcols
        incremental_update_size = incremental_update_size + lcols
 
        icol = icol - lcols
        ! count colums from first column of global block to current index
        CALL local_size_offset_1d(icol,nb,colidx-n+1,colidx-n+1,0, &
                                      mpirank_cols,mpiprocs_cols, &
                                      lcols,baseoffset,offset)

        IF (lcols .GT. 0) THEN

            !print *,'updating trailing matrix'

                        IF (updatemode .EQ. ICHAR('I')) THEN
                                PRINT *,'pdgeqrf_2dcomm: incremental update not yet implemented! rev=1'
                        ELSE IF (updatemode .EQ. ICHAR('F')) THEN
                                ! full update no merging
                                CALL tum_pdlarfb_1dcomm(m,mb,lcols,& 
                                update_lcols,a(1,offset),lda, &
                                v(1,update_voffset),ldv, &
                                work(temptau_offset+update_voffset-1), & 
                                t(update_voffset,update_voffset),ldt, &
                                rowidx,idx,1,mpicomm_rows,work(work_offset), & 
                                lwork)
                        ELSE 
                                ! full update + merging default
                                CALL tum_tmerge_pdlarfb_1dcomm(m,mb,lcols,n- & 
                                (update_voffset+update_lcols-1),update_lcols,&
                                v(1,update_voffset),ldv, &
                                t(update_voffset,update_voffset),ldt, &
                                a(1,offset),lda,rowidx,1,updatemode, &
                                mpicomm_rows,work(work_offset),lwork)
                        END IF
        ELSE
                        IF (updatemode .EQ. ICHAR('I')) THEN
                        PRINT *,'sole merging of (incremental) T matrix', &
                        mpirank_cols, n- & 
                        (update_voffset+incremental_update_size-1)
                        CALL tum_tmerge_pdlarfb_1dcomm(m,mb,0,n- &
                        (update_voffset+incremental_update_size-1), &
                        incremental_update_size,v(1,update_voffset),ldv, &
                        t(update_voffset,update_voffset),ldt, &
                        a,lda,rowidx,1,updatemode,mpicomm_rows,&
                        work(work_offset),lwork)

                                ! reset for upcoming incremental updates
                                incremental_update_size = 0
                        ELSE IF (updatemode .EQ. ICHAR('M')) THEN
                                ! final merge
                         CALL tum_tmerge_pdlarfb_1dcomm(m,mb,0,n- & 
                         (update_voffset+update_lcols-1),update_lcols, &
                         v(1,update_voffset),ldv, &
                         t(update_voffset,update_voffset),ldt, & 
                         a,lda,rowidx,1,updatemode,mpicomm_rows, &
                         work(work_offset),lwork)
                        ELSE 
                                ! full updatemode - nothing to update
                        END IF

                        ! reset for upcoming incremental updates
                        incremental_update_size = 0
        END IF
    END DO
 
    IF ((tmerge .GT. 0) .AND. (updatemode .EQ. ICHAR('F'))) THEN
        ! finally merge all small T parts
        CALL tum_pdlarft_tree_merge_1dcomm(m,mb,n,size2d,tmerge,v,& 
        ldv,t,ldt,rowidx,rev,mpicomm_rows,work,lwork)
    END IF

    !print *,'stop decomposition',rowidx,colidx
END SUBROUTINE tum_pdgeqrf_2dcomm

SUBROUTINE tum_pdgeqrf_1dcomm(a,lda,v,ldv,tau,t,ldt,work,lwork, &
m,n,mb,baseidx,rowidx,rev,trans,PQRPARAM,mpicomm,blockheuristic)
!    use ELPA1
  
 
    ! parameter setup
    INTEGER                                  :: eps_, gmode_, rank_

    PARAMETER   (gmode_ = 1,rank_ = 2,eps_=3)

    ! input variables (local)
    INTEGER lda,lwork,ldv,ldt
    DOUBLE PRECISION a(lda,*),v(ldv,*),tau(*),t(ldt,*),work(*)

    ! input variables (global)
    INTEGER m,n,mb,baseidx,rowidx,rev,trans,mpicomm
    INTEGER PQRPARAM(*)

    ! derived input variables
  
    ! derived further input variables from TUM_PQRPARAM
    INTEGER size1d,updatemode,tmerge

    ! output variables (global)
    DOUBLE PRECISION blockheuristic(*)
    ! local scalars
    INTEGER nr_blocks,remainder,current_block,aoffset,idx,updatesize
    DOUBLE PRECISION pdgeqr2_size(1),pdlarfb_size(1),tmerge_tree_size(1)

    size1d = MAX(MIN(PQRPARAM(1),n),1)
    updatemode = PQRPARAM(2)
    tmerge = PQRPARAM(3)

    IF (lwork .EQ. -1) THEN
        CALL tum_pdgeqr2_1dcomm(a,lda,v,ldv,tau,t,ldt,pdgeqr2_size,-1, & 
                                m,size1d,mb,baseidx,baseidx,rev,trans,PQRPARAM(4),mpicomm,blockheuristic)

        ! reserve more space for incremental mode
        CALL tum_tmerge_pdlarfb_1dcomm(m,mb,n,n,n,v,ldv,t,ldt, &
                                       a,lda,baseidx,rev,updatemode,mpicomm,pdlarfb_size,-1)
 
        CALL tum_pdlarft_tree_merge_1dcomm(m,mb,n,size1d,tmerge,v,ldv,t,ldt,baseidx,rev,mpicomm,tmerge_tree_size,-1)

        work(1) = MAX(pdlarfb_size(1),pdgeqr2_size(1),tmerge_tree_size(1))
        RETURN
    END IF

        nr_blocks = n / size1d
        remainder = n - nr_blocks*size1d
  
        current_block = 0
        DO WHILE (current_block .LT. nr_blocks)
            idx = rowidx-current_block*size1d
            updatesize = n-(current_block+1)*size1d
            aoffset = 1+updatesize

            CALL tum_pdgeqr2_1dcomm(a(1,aoffset),lda,v(1,aoffset),ldv,tau(aoffset),t(aoffset,aoffset),ldt,work,lwork, & 
                                    m,size1d,mb,baseidx,idx,1,trans,PQRPARAM(4),mpicomm,blockheuristic)

            IF (updatemode .EQ. ICHAR('M')) THEN
                ! full update + merging
                CALL tum_tmerge_pdlarfb_1dcomm(m,mb,updatesize,current_block*size1d,size1d, & 
                                               v(1,aoffset),ldv,t(aoffset,aoffset),ldt, &
                                               a,lda,baseidx,1,ICHAR('F'),mpicomm,work,lwork)
            ELSE IF (updatemode .EQ. ICHAR('I')) THEN
                IF (updatesize .GE. size1d) THEN
                    ! incremental update + merging
                    CALL tum_tmerge_pdlarfb_1dcomm(m,mb,size1d,current_block*size1d,size1d, & 
                                                   v(1,aoffset),ldv,t(aoffset,aoffset),ldt, &
                                                   a(1,aoffset-size1d),lda,baseidx,1,updatemode,mpicomm,work,lwork)

                ELSE ! only remainder left
                    ! incremental update + merging
                    CALL tum_tmerge_pdlarfb_1dcomm(m,mb,remainder,current_block*size1d,size1d, & 
                                                   v(1,aoffset),ldv,t(aoffset,aoffset),ldt, &
                                                   a(1,1),lda,baseidx,1,updatemode,mpicomm,work,lwork)
                END IF
            ELSE ! full update no merging is default
                ! full update no merging
                CALL tum_pdlarfb_1dcomm(m,mb,updatesize,size1d,a,lda,v(1,aoffset),ldv, &
                                        tau(aoffset),t(aoffset,aoffset),ldt,baseidx,idx,1,mpicomm,work,lwork)
            END IF

            ! move on to next block
            current_block = current_block+1
        END DO

        IF (remainder .GT. 0) THEN
            aoffset = 1
            idx = rowidx-size1d*nr_blocks
            CALL tum_pdgeqr2_1dcomm(a(1,aoffset),lda,v,ldv,tau,t,ldt,work,lwork, &
                                    m,remainder,mb,baseidx,idx,1,trans,PQRPARAM(4),mpicomm,blockheuristic)

            IF ((updatemode .EQ. ICHAR('I')) .OR. (updatemode .EQ. ICHAR('M'))) THEN
                ! final merging
                CALL tum_tmerge_pdlarfb_1dcomm(m,mb,0,size1d*nr_blocks,remainder, & 
                                               v,ldv,t,ldt, &
                                               a,lda,baseidx,1,updatemode,mpicomm,work,lwork) ! updatemode argument does not matter
            END IF
        END IF
 
    IF ((tmerge .GT. 0) .AND. (updatemode .EQ. ICHAR('F'))) THEN
        ! finally merge all small T parts
        CALL tum_pdlarft_tree_merge_1dcomm(m,mb,n,size1d,tmerge,v,ldv,t,ldt,baseidx,rev,mpicomm,work,lwork)
    END IF
END SUBROUTINE tum_pdgeqrf_1dcomm

! local a and tau are assumed to be positioned at the right column from a local
! perspective
! TODO: if local amount of data turns to zero the algorithm might produce wrong
! results (probably due to old buffer contents)
SUBROUTINE tum_pdgeqr2_1dcomm(a,lda,v,ldv,tau,t,ldt,work,lwork,m,n,mb,baseidx, &
                              rowidx,rev,trans,PQRPARAM,mpicomm,blockheuristic)
    
    INTEGER                                  :: eps_, gmode_, rank_, upmode1_

! parameter setup

    PARAMETER   (gmode_ = 1,rank_ = 2,eps_=3, upmode1_=4)

    ! input variables (local)
    INTEGER lda,lwork,ldv,ldt
    DOUBLE PRECISION a(lda,*),v(ldv,*),tau(*),t(ldt,*),work(*)

    ! input variables (global)
    INTEGER m,n,mb,baseidx,rowidx,rev,trans,mpicomm
    INTEGER PQRPARAM(*)
 
    ! output variables (global)
    DOUBLE PRECISION blockheuristic(*)
 
    ! derived further input variables from TUM_PQRPARAM
    INTEGER maxrank,hgmode,updatemode

    ! local scalars
    INTEGER icol,incx,idx
    DOUBLE PRECISION pdlarfg_size(1),pdlarf_size(1),total_size
    DOUBLE PRECISION pdlarfg2_size(1),pdlarfgk_size(1),pdlarfl2_size(1)
    DOUBLE PRECISION pdlarft_size(1),pdlarfb_size(1),pdlarft_pdlarfb_size(1),tmerge_pdlarfb_size(1)
    INTEGER mpirank,mpiprocs,mpierr
    INTEGER rank,lastcol,actualrank,nextrank
    INTEGER update_cols,decomposition_cols
    INTEGER current_column
 
    maxrank = MIN(PQRPARAM(1),n)
    updatemode = PQRPARAM(2)
    hgmode = PQRPARAM(4)

    CALL MPI_Comm_rank(mpicomm, mpirank, mpierr)
    CALL MPI_Comm_size(mpicomm, mpiprocs, mpierr)

    IF (trans .EQ. 1) THEN
        incx = lda
    ELSE
        incx = 1
    END IF
    
    IF (lwork .EQ. -1) THEN
        CALL tum_pdlarfg_1dcomm(a,incx,tau(1),pdlarfg_size(1),-1,n,rowidx,mb,hgmode,rev,mpicomm)
        CALL tum_pdlarfl_1dcomm(v,1,baseidx,a,lda,tau(1),pdlarf_size(1),-1,m,n,rowidx,mb,rev,mpicomm)
        CALL tum_pdlarfg2_1dcomm_ref(a,lda,tau,t,ldt,v,ldv,baseidx,pdlarfg2_size(1),-1,m,rowidx,mb,PQRPARAM,rev,mpicomm,actualrank)
        CALL tum_pdlarfgk_1dcomm(a,lda,tau,t,ldt,v,ldv,baseidx,pdlarfgk_size(1),-1,m,n,rowidx,mb,PQRPARAM,rev,mpicomm,actualrank)
        CALL tum_pdlarfl2_tmatrix_1dcomm(v,ldv,baseidx,a,lda,t,ldt,pdlarfl2_size(1),-1,m,n,rowidx,mb,rev,mpicomm)
        pdlarft_size(1) = 0.0d0
        CALL tum_pdlarfb_1dcomm(m,mb,n,n,a,lda,v,ldv,tau,t,ldt,baseidx,rowidx,1,mpicomm,pdlarfb_size(1),-1)
        pdlarft_pdlarfb_size(1) = 0.0d0
        CALL tum_tmerge_pdlarfb_1dcomm(m,mb,n,n,n,v,ldv,t,ldt,a,lda,rowidx,rev,updatemode,mpicomm,tmerge_pdlarfb_size(1),-1)

        total_size = MAX(pdlarfg_size(1),pdlarf_size(1),pdlarfg2_size(1), &
        pdlarfgk_size(1),pdlarfl2_size(1),pdlarft_size(1),pdlarfb_size(1), &
        pdlarft_pdlarfb_size(1),tmerge_pdlarfb_size(1))

        work(1) = total_size
        RETURN
    END IF

        icol = 1
        lastcol = MIN(rowidx,n)
        decomposition_cols = lastcol
        update_cols = n
        DO WHILE (decomposition_cols .GT. 0) ! local qr block
            icol = lastcol-decomposition_cols+1
            idx = rowidx-icol+1

            ! get possible rank size
            ! limited by number of columns and remaining rows
            rank = MIN(n-icol+1,maxrank,idx)

            current_column = n-icol+1-rank+1

            IF (rank .EQ. 1) THEN

                CALL tum_pdlarfg_1dcomm(a(1,current_column),incx, &
                                        tau(current_column),work,lwork, &
                                        m,idx,mb,hgmode,1,mpicomm)

                v(1:ldv,current_column) = 0.0d0
                CALL tum_pdlarfg_copy_1dcomm(a(1,current_column),incx, &
                                             v(1,current_column),1, &
                                             m,baseidx,idx,mb,1,mpicomm)

                ! initialize t matrix part
                t(current_column,current_column) = tau(current_column)

                actualrank = 1

            ELSE IF (rank .EQ. 2) THEN 
                CALL tum_pdlarfg2_1dcomm_ref(a(1,current_column),lda,tau(current_column), &
                                             t(current_column,current_column),ldt,v(1,current_column),ldv, &
                                            baseidx,work,lwork,m,idx,mb,PQRPARAM,1,mpicomm,actualrank)
            
            ELSE 
                CALL tum_pdlarfgk_1dcomm(a(1,current_column),lda,tau(current_column), &
                                         t(current_column,current_column),ldt,v(1,current_column),ldv, &
                                         baseidx,work,lwork,m,rank,idx,mb,PQRPARAM,1,mpicomm,actualrank)
 
            END IF
  
            blockheuristic(actualrank) = blockheuristic(actualrank) + 1

            ! the blocked decomposition versions already updated their non
            ! decomposed parts using their information after communication
            update_cols = decomposition_cols - rank
            decomposition_cols = decomposition_cols - actualrank
 
            ! needed for incremental update
            nextrank = MIN(n-(lastcol-decomposition_cols+1)+1,maxrank, &
            rowidx-(lastcol-decomposition_cols+1)+1)

            IF (current_column .GT. 1) THEN
                idx = rowidx-icol+1
 
                IF (updatemode .EQ. ICHAR('I')) THEN
                    ! incremental update + merging
                    CALL tum_tmerge_pdlarfb_1dcomm(m,mb,nextrank- &
                    (rank-actualrank),n-(current_column+rank-1), &
                    actualrank,v(1,current_column+(rank-actualrank)), &
                    ldv, t(current_column+(rank-actualrank), &
                    current_column+(rank-actualrank)),ldt, &
                    a(1,current_column-nextrank+(rank-actualrank)), &
                    lda,baseidx,rev,updatemode,mpicomm,work,lwork)
                ELSE
                    ! full update + merging
                    CALL tum_tmerge_pdlarfb_1dcomm(m,mb,update_cols,n- &
                    (current_column+rank-1),actualrank,v(1,current_column+ &
                    (rank-actualrank)),ldv, &
                    t(current_column+(rank-actualrank), &
                    current_column+(rank-actualrank)),ldt, &
                    a(1,1),lda,baseidx,rev,updatemode,mpicomm,work,lwork)
                END IF
            ELSE
                CALL tum_tmerge_pdlarfb_1dcomm(m,mb,0,n- &
                (current_column+rank-1),actualrank, &
                v(1,current_column+(rank-actualrank)),ldv, &
                t(current_column+(rank-actualrank),current_column+ &
                (rank-actualrank)),ldt, &
                a,lda,baseidx,rev,updatemode,mpicomm,work,lwork)
            END IF

        END DO
END SUBROUTINE tum_pdgeqr2_1dcomm

! incx == 1: column major
! incx != 1: row major
SUBROUTINE tum_pdlarfg_1dcomm(x,incx,tau,work,lwork,n,idx,nb,hgmode,rev,mpi_comm)
    
    
    INTEGER                                  :: eps_, gmode_, rank_

! parameter setup

    PARAMETER   (gmode_ = 1,rank_ = 2,eps_=3)

    ! input variables (local)
    INTEGER incx,lwork,hgmode
    DOUBLE PRECISION x(*),work(*)

    ! input variables (global)
    INTEGER mpi_comm,nb,idx,n,rev

    ! output variables (global)
    DOUBLE PRECISION tau
    ! local scalars
    INTEGER mpierr,mpirank,mpiprocs,mpirank_top
    INTEGER sendsize,recvsize
    INTEGER local_size,local_offset,baseoffset
    INTEGER topidx,top,iproc
    DOUBLE PRECISION alpha,xnorm,dot,xf

    ! external functions
    DOUBLE PRECISION ddot,dlapy2,dnrm2
    EXTERNAL ddot,dscal,dlapy2,dnrm2

    ! intrinsic
    INTRINSIC sign

        IF (idx .LE. 1) THEN
                tau = 0.0d0
                RETURN
        END IF

    CALL MPI_Comm_rank(mpi_comm, mpirank, mpierr)
    CALL MPI_Comm_size(mpi_comm, mpiprocs, mpierr)

    ! calculate expected work size and store in work(1)
    IF (hgmode .EQ. ICHAR('s')) THEN
        ! allreduce (MPI_SUM)
        sendsize = 2
        recvsize = sendsize
    ELSE IF (hgmode .EQ. ICHAR('x')) THEN
        ! alltoall
        sendsize = mpiprocs*2
        recvsize = sendsize
    ELSE IF (hgmode .EQ. ICHAR('g')) THEN
        ! allgather
        sendsize = 2
        recvsize = mpiprocs*sendsize
    ELSE
        ! no exchange at all (benchmarking)
        sendsize = 2
        recvsize = sendsize
    END IF

    IF (lwork .EQ. -1) THEN
        work(1) = DBLE(sendsize + recvsize)
        RETURN
    END IF
 
    ! Processor id for global index of top element
    mpirank_top = MOD((idx-1)/nb,mpiprocs)
    IF (mpirank .EQ. mpirank_top) THEN
        topidx = local_index(idx,mpirank_top,mpiprocs,nb,0)
        top = 1+(topidx-1)*incx
    END IF
 
        CALL local_size_offset_1d(n,nb,idx,idx-1,rev,mpirank,mpiprocs, &
                                                          local_size,baseoffset,local_offset)

    local_offset = local_offset * incx

    ! calculate and exchange information
    IF (hgmode .EQ. ICHAR('s')) THEN
        IF (mpirank .EQ. mpirank_top) THEN
            alpha = x(top)
        ELSE 
            alpha = 0.0d0
        END IF

        dot = ddot(local_size, &
                   x(local_offset), incx, &
                   x(local_offset), incx)

        work(1) = alpha
        work(2) = dot
        
        CALL mpi_allreduce(work(1),work(sendsize+1), &
                           sendsize,mpi_real8,mpi_sum, &
                           mpi_comm,mpierr)

        alpha = work(sendsize+1)
        xnorm = SQRT(work(sendsize+2))
    ELSE IF (hgmode .EQ. ICHAR('x')) THEN
        IF (mpirank .EQ. mpirank_top) THEN
            alpha = x(top)
        ELSE 
            alpha = 0.0d0
        END IF

        xnorm = dnrm2(local_size, x(local_offset), incx)

        DO iproc=0,mpiprocs-1
            work(2*iproc+1) = alpha
            work(2*iproc+2) = xnorm
        END DO

        CALL mpi_alltoall(work(1),2,mpi_real8, &
                          work(sendsize+1),2,mpi_real8, &
                          mpi_comm,mpierr)

        ! extract alpha value
        alpha = work(sendsize+1+mpirank_top*2)

        ! copy norm parts of buffer to beginning
        DO iproc=0,mpiprocs-1
            work(iproc+1) = work(sendsize+1+2*iproc+1)
        END DO

        xnorm = dnrm2(mpiprocs, work(1), 1)
    ELSE IF (hgmode .EQ. ICHAR('g')) THEN
        IF (mpirank .EQ. mpirank_top) THEN
            alpha = x(top)
        ELSE 
            alpha = 0.0d0
        END IF

        xnorm = dnrm2(local_size, x(local_offset), incx)
        work(1) = alpha
        work(2) = xnorm

        ! allgather
        CALL mpi_allgather(work(1),sendsize,mpi_real8, &
                          work(sendsize+1),sendsize,mpi_real8, &
                          mpi_comm,mpierr)

        ! extract alpha value
        alpha = work(sendsize+1+mpirank_top*2)
 
        ! copy norm parts of buffer to beginning
        DO iproc=0,mpiprocs-1
            work(iproc+1) = work(sendsize+1+2*iproc+1)
        END DO

        xnorm = dnrm2(mpiprocs, work(1), 1)
    ELSE
        ! dnrm2
        xnorm = dnrm2(local_size, x(local_offset), incx)

        IF (mpirank .EQ. mpirank_top) THEN
            alpha = x(top)
        ELSE 
            alpha = 0.0d0
        END IF

        ! no exchange at all (benchmarking)
 
        xnorm = 0.0d0
    END IF

    !print *,'ref hg:', idx,xnorm,alpha
    !print *,x(1:n)

    ! calculate householder information
    IF (xnorm .EQ. 0.0d0) THEN
        ! H = I

        tau = 0.0d0
    ELSE
        ! General case

        CALL hh_transform_real(alpha,xnorm**2,xf,tau)
        IF (mpirank .EQ. mpirank_top) THEN
            x(top) = alpha
        END IF

        CALL dscal(local_size, xf, &
                   x(local_offset), incx)
 
        ! TODO: reimplement norm rescale method of 
        ! original PDLARFG using mpi?

    END IF

    ! useful for debugging
    !print *,'hg:mpirank,idx,beta,alpha:',mpirank,idx,beta,alpha,1.0d0/(beta+alpha),tau
    !print *,x(1:n)
END SUBROUTINE tum_pdlarfg_1dcomm

SUBROUTINE tum_pdlarfg2_1dcomm_ref(a,lda,tau,t,ldt,v,ldv,baseidx,work,lwork,m,idx,mb,PQRPARAM,rev,mpicomm,actualk)
 
    ! parameter setup
    INTEGER                                  :: eps_, gmode_, rank_, upmode1_

    PARAMETER   (gmode_ = 1,rank_ = 2,eps_=3, upmode1_=4)

    ! input variables (local)
    INTEGER lda,lwork,ldv,ldt
    DOUBLE PRECISION a(lda,*),v(ldv,*),tau(*),work(*),t(ldt,*)

    ! input variables (global)
    INTEGER m,idx,baseidx,mb,rev,mpicomm
    INTEGER PQRPARAM(*)
 
    ! output variables (global)
    INTEGER actualk
 
    ! derived input variables from TUM_PQRPARAM
    INTEGER eps
    ! local scalars
    DOUBLE PRECISION dseedwork_size(1)
    INTEGER seedwork_size,seed_size
    INTEGER seedwork_offset,seed_offset
    LOGICAL accurate

    CALL tum_pdlarfg2_1dcomm_seed(a,lda,dseedwork_size(1),-1,work,m,mb,idx,rev,mpicomm)
    seedwork_size = dseedwork_size(1)
    seed_size = seedwork_size
 
    IF (lwork .EQ. -1) THEN
        work(1) = seedwork_size + seed_size
        RETURN
    END IF

    seedwork_offset = 1
    seed_offset = seedwork_offset + seedwork_size

    eps = PQRPARAM(3)

    ! check for border cases (only a 2x2 matrix left)
        IF (idx .LE. 1) THEN
                tau(1:2) = 0.0d0
                t(1:2,1:2) = 0.0d0
                RETURN
        END IF

    CALL tum_pdlarfg2_1dcomm_seed(a,lda,work(seedwork_offset),lwork,work(seed_offset),m,mb,idx,rev,mpicomm)

        IF (eps .GT. 0) THEN
            accurate = tum_pdlarfg2_1dcomm_check(work(seed_offset),eps)
        ELSE
            accurate = .TRUE.
        END IF

        CALL tum_pdlarfg2_1dcomm_vector(a(1,2),1,tau(2),work(seed_offset), &
                                        m,mb,idx,0,1,mpicomm)

        CALL tum_pdlarfg_copy_1dcomm(a(1,2),1, &
                                     v(1,2),1, &
                                     m,baseidx,idx,mb,1,mpicomm)

        CALL tum_pdlarfg2_1dcomm_update(v(1,2),1,baseidx,a(1,1),lda,work(seed_offset),m,idx,mb,rev,mpicomm)

        ! check for 2x2 matrix case => only one householder vector will be
        ! generated
        IF (idx .GT. 2) THEN
            IF (accurate .EQV. .TRUE.) THEN
                CALL tum_pdlarfg2_1dcomm_vector(a(1,1),1,tau(1),work(seed_offset), &
                                                m,mb,idx-1,1,1,mpicomm)

                CALL tum_pdlarfg_copy_1dcomm(a(1,1),1, &
                                             v(1,1),1, &
                                             m,baseidx,idx-1,mb,1,mpicomm)

                ! generate fuse element
                CALL tum_pdlarfg2_1dcomm_finalize_tmatrix(work(seed_offset),tau,t,ldt)

                actualk = 2
            ELSE
                t(1,1) = 0.0d0
                t(1,2) = 0.0d0
                t(2,2) = tau(2)

                actualk = 1
            END IF
        ELSE
            t(1,1) = 0.0d0
            t(1,2) = 0.0d0
            t(2,2) = tau(2)

            ! no more vectors to create

            tau(1) = 0.0d0

            actualk = 2

            !print *,'rank2: no more data'
        END IF
END SUBROUTINE tum_pdlarfg2_1dcomm_ref

SUBROUTINE tum_pdlarfg2_1dcomm_seed(a,lda,work,lwork,seed,n,nb,idx,rev,mpicomm)
    
    INTEGER                                  :: lda, lwork

! input variables (local)

    DOUBLE PRECISION a(lda,*),work(*),seed(*)

    ! input variables (global)
    INTEGER n,nb,idx,rev,mpicomm
 
    ! output variables (global)

    ! external functions
    DOUBLE PRECISION ddot
    EXTERNAL ddot


    ! local scalars
    DOUBLE PRECISION top11,top21,top12,top22
    DOUBLE PRECISION dot11,dot12,dot22
    INTEGER mpirank,mpiprocs,mpierr
    INTEGER mpirank_top11,mpirank_top21
    INTEGER top11_offset,top21_offset
    INTEGER baseoffset
    INTEGER local_offset1,local_size1
    INTEGER local_offset2,local_size2

    IF (lwork .EQ. -1) THEN
        work(1) = DBLE(8)
        RETURN
    END IF
  
    CALL MPI_Comm_rank(mpicomm, mpirank, mpierr)
    CALL MPI_Comm_size(mpicomm, mpiprocs, mpierr)

        CALL local_size_offset_1d(n,nb,idx,idx-1,rev,mpirank,mpiprocs, &
                              local_size1,baseoffset,local_offset1)

        CALL local_size_offset_1d(n,nb,idx,idx-2,rev,mpirank,mpiprocs, &
                              local_size2,baseoffset,local_offset2)

        mpirank_top11 = MOD((idx-1)/nb,mpiprocs)
        mpirank_top21 = MOD((idx-2)/nb,mpiprocs)

        top11_offset = local_index(idx,mpirank_top11,mpiprocs,nb,0)
        top21_offset = local_index(idx-1,mpirank_top21,mpiprocs,nb,0)

        IF (mpirank_top11 .EQ. mpirank) THEN
            top11 = a(top11_offset,2)
            top12 = a(top11_offset,1)
        ELSE
            top11 = 0.0d0
            top12 = 0.0d0
        END IF

        IF (mpirank_top21 .EQ. mpirank) THEN
            top21 = a(top21_offset,2)
            top22 = a(top21_offset,1)
        ELSE
            top21 = 0.0d0
            top22 = 0.0d0
        END IF

        ! calculate 3 dot products
        dot11 = ddot(local_size1,a(local_offset1,2),1,a(local_offset1,2),1)
        dot12 = ddot(local_size1,a(local_offset1,2),1,a(local_offset1,1),1)
        dot22 = ddot(local_size2,a(local_offset2,1),1,a(local_offset2,1),1)

    ! store results in work buffer
    work(1) = top11
    work(2) = dot11
    work(3) = top12
    work(4) = dot12
    work(5) = top21
    work(6) = top22
    work(7) = dot22
    work(8) = 0.0d0 ! fill up buffer

    ! exchange partial results
    CALL mpi_allreduce(work, seed, 8, mpi_real8, mpi_sum, &
                       mpicomm, mpierr)
END SUBROUTINE tum_pdlarfg2_1dcomm_seed

LOGICAL FUNCTION tum_pdlarfg2_1dcomm_check(seed,eps)
    IMPLICIT NONE

    ! input variables
    DOUBLE PRECISION seed(*)
    INTEGER eps

    ! local scalars
    DOUBLE PRECISION epsd,first,second,first_second,estimate
    LOGICAL accurate
    DOUBLE PRECISION dot11,dot12,dot22
    DOUBLE PRECISION top11,top12,top21,top22
 
    EPSD = EPS
  
    top11 = seed(1)
    dot11 = seed(2)
    top12 = seed(3)
    dot12 = seed(4)
        
    top21 = seed(5)
    top22 = seed(6)
    dot22 = seed(7)

    ! reconstruct the whole inner products 
    ! (including squares of the top elements)
    first = dot11 + top11*top11
    second = dot22 + top22*top22 + top12*top12
    first_second = dot12 + top11*top12
 
    ! zero Householder vector (zero norm) case
    IF (first*second .EQ. 0.0d0) THEN
       tum_pdlarfg2_1dcomm_check = .FALSE.
       RETURN
    END IF

    estimate = ABS((first_second*first_second)/(first*second)) 

    !print *,'estimate:',estimate
    
    ! if accurate the following check holds
    accurate = (estimate .LE. (epsd/(1.0d0+epsd)))
  
    tum_pdlarfg2_1dcomm_check = accurate
END FUNCTION tum_pdlarfg2_1dcomm_check

! id=0: first vector
! id=1: second vector
SUBROUTINE tum_pdlarfg2_1dcomm_vector(x,incx,tau,seed,n,nb,idx,id,rev,mpicomm)
    USE ELPA1
    USE tum_utils
    INTEGER                                  :: incx

! input variables (local)

    DOUBLE PRECISION x(*),seed(*),tau

    ! input variables (global)
    INTEGER n,nb,idx,id,rev,mpicomm
 
    ! output variables (global)

    ! external functions
    DOUBLE PRECISION dlapy2
    EXTERNAL dlapy2,dscal

    ! local scalars
    INTEGER mpirank,mpirank_top,mpiprocs,mpierr
    DOUBLE PRECISION alpha,dot,beta,xnorm
    INTEGER local_size,baseoffset,local_offset,top,topidx

    CALL MPI_Comm_rank(mpicomm, mpirank, mpierr)
    CALL MPI_Comm_size(mpicomm, mpiprocs, mpierr)

    CALL local_size_offset_1d(n,nb,idx,idx-1,rev,mpirank,mpiprocs, &
                                  local_size,baseoffset,local_offset)

    local_offset = local_offset * incx

    ! Processor id for global index of top element
    mpirank_top = MOD((idx-1)/nb,mpiprocs)
    IF (mpirank .EQ. mpirank_top) THEN
        topidx = local_index(idx,mpirank_top,mpiprocs,nb,0)
        top = 1+(topidx-1)*incx
    END IF

    alpha = seed(id*5+1)
    dot = seed(id*5+2)
    
    xnorm = SQRT(dot)

    IF (xnorm .EQ. 0.0d0) THEN
        ! H = I

        tau = 0.0d0
    ELSE
        ! General case

        beta = SIGN(dlapy2(alpha, xnorm), alpha)
        tau = (beta+alpha) / beta

        !print *,'hg2',tau,xnorm,alpha
        
        CALL dscal(local_size, 1.0d0/(beta+alpha), &
                   x(local_offset), incx)
 
        ! TODO: reimplement norm rescale method of 
        ! original PDLARFG using mpi?

        IF (mpirank .EQ. mpirank_top) THEN
            x(top) = -beta
        END IF

        seed(8) = beta
    END IF
END SUBROUTINE tum_pdlarfg2_1dcomm_vector

SUBROUTINE tum_pdlarfg2_1dcomm_update(v,incv,baseidx,a,lda,seed,n,idx,nb,rev,mpicomm)
    USE ELPA1
    USE tum_utils
    INTEGER                                  :: incv, lda

! input variables (local)

    DOUBLE PRECISION v(*),a(lda,*),seed(*)

    ! input variables (global)
    INTEGER n,baseidx,idx,nb,rev,mpicomm
 
    ! output variables (global)

    ! external functions
    EXTERNAL daxpy

    ! local scalars
    INTEGER mpirank,mpiprocs,mpierr
    INTEGER local_size,local_offset,baseoffset
    DOUBLE PRECISION z,coeff,beta
    DOUBLE PRECISION dot11,dot12,dot22
    DOUBLE PRECISION top11,top12,top21,top22
 
    CALL MPI_Comm_rank(mpicomm, mpirank, mpierr)
    CALL MPI_Comm_size(mpicomm, mpiprocs, mpierr)


    ! seed should be updated by previous householder generation
    ! Update inner product of this column and next column vector
    top11 = seed(1)
    dot11 = seed(2)
    top12 = seed(3)
    dot12 = seed(4)
        
    top21 = seed(5)
    top22 = seed(6)
    dot22 = seed(7)
    beta = seed(8)
    
    CALL local_size_offset_1d(n,nb,baseidx,idx,rev,mpirank,mpiprocs, &
                              local_size,baseoffset,local_offset)
    baseoffset = baseoffset * incv

    ! zero Householder vector (zero norm) case
    IF (beta .EQ. 0.0d0) THEN
        RETURN
    END IF
    z = (dot12 + top11 * top12) / beta + top12

    !print *,'hg2 update:',baseidx,idx,mpirank,local_size

    CALL daxpy(local_size, -z, v(baseoffset),1, a(local_offset,1),1)
    
    ! prepare a full dot22 for update
    dot22 = dot22 + top22*top22

    ! calculate coefficient
    COEFF = z / (top11 + beta)
  
    ! update inner product of next vector
    dot22 = dot22 - coeff * (2*dot12 - coeff*dot11)
 
    ! update dot12 value to represent update with first vector 
    ! (needed for T matrix)
    dot12 = dot12 - COEFF * dot11 
    
    ! update top element of next vector
    top22 = top22 - coeff * top21
    seed(6) = top22

    ! restore separated dot22 for vector generation
    seed(7) = dot22  - top22*top22 

    !------------------------------------------------------
    ! prepare elements for T matrix
    seed(4) = dot12

    ! prepare dot matrix for fuse element of T matrix
    ! replace top11 value with -beta1
    seed(1) = beta
END SUBROUTINE tum_pdlarfg2_1dcomm_update

! run this function after second vector
SUBROUTINE tum_pdlarfg2_1dcomm_finalize_tmatrix(seed,tau,t,ldt)

    INTEGER                                  :: ldt

    DOUBLE PRECISION seed(*),t(ldt,*),tau(*)
    DOUBLE PRECISION dot12,beta1,top21,beta2
 
    beta1 = seed(1)
    dot12 = seed(4)
    top21 = seed(5)
    beta2 = seed(8)
 
    !print *,'beta1 beta2',beta1,beta2

    dot12 = dot12 / beta2 + top21
    dot12 = -(dot12 / beta1)

    t(1,1) = tau(1)
    t(1,2) = dot12
    t(2,2) = tau(2)
END SUBROUTINE tum_pdlarfg2_1dcomm_finalize_tmatrix

SUBROUTINE tum_pdlarfgk_1dcomm(a,lda,tau,t,ldt,v,ldv,baseidx,work,lwork,m,k,idx,mb,PQRPARAM,rev,mpicomm,actualk)

 
    ! parameter setup

    ! input variables (local)
    INTEGER                                  :: lda, ldt, ldv, lwork

    DOUBLE PRECISION a(lda,*),v(ldv,*),tau(*),work(*),t(ldt,*)

    ! input variables (global)
    INTEGER m,k,idx,baseidx,mb,rev,mpicomm
    INTEGER PQRPARAM(*)
 
    ! output variables (global)
    INTEGER actualk

    ! local scalars
    INTEGER ivector
    DOUBLE PRECISION pdlarfg_size(1),pdlarf_size(1)
    DOUBLE PRECISION pdlarfgk_1dcomm_seed_size(1),pdlarfgk_1dcomm_check_size(1)
    DOUBLE PRECISION pdlarfgk_1dcomm_update_size(1)
    INTEGER seedC_size,seedC_offset
    INTEGER seedD_size,seedD_offset
    INTEGER work_offset

    seedC_size = k*k
    seedC_offset = 1
    seedD_size = k*k
    seedD_offset = seedC_offset + seedC_size
    work_offset = seedD_offset + seedD_size

    IF (lwork .EQ. -1) THEN
        CALL tum_pdlarfg_1dcomm(a,1,tau(1),pdlarfg_size(1),-1,m,baseidx,mb,PQRPARAM(4),rev,mpicomm)
        CALL tum_pdlarfl_1dcomm(v,1,baseidx,a,lda,tau(1),pdlarf_size(1),-1,m,k,baseidx,mb,rev,mpicomm)
        CALL tum_pdlarfgk_1dcomm_seed(a,lda,baseidx,pdlarfgk_1dcomm_seed_size(1),-1,work,work,m,k,mb,mpicomm)
        !call tum_pdlarfgk_1dcomm_check(work,work,k,PQRPARAM,pdlarfgk_1dcomm_check_size(1),-1,actualk)
        CALL tum_pdlarfgk_1dcomm_check_improved(work,work,k,PQRPARAM,pdlarfgk_1dcomm_check_size(1),-1,actualk)
        CALL tum_pdlarfgk_1dcomm_update(a,lda,baseidx,pdlarfgk_1dcomm_update_size(1),-1,work,work,k,k,1,work,m,mb,rev,mpicomm)
        work(1) = MAX(pdlarfg_size(1),pdlarf_size(1), &
        pdlarfgk_1dcomm_seed_size(1),pdlarfgk_1dcomm_check_size(1), &
        pdlarfgk_1dcomm_update_size(1)) + DBLE(seedC_size + seedD_size);
        RETURN
    END IF

        CALL tum_pdlarfgk_1dcomm_seed(a(1,1),lda,idx,work(work_offset), & 
        lwork,work(seedC_offset),work(seedD_offset),m,k,mb,mpicomm)
        !call tum_pdlarfgk_1dcomm_check(work(seedC_offset),work(seedD_offset),k,PQRPARAM,work(work_offset),lwork,actualk)
        CALL tum_pdlarfgk_1dcomm_check_improved(work(seedC_offset),&
        work(seedD_offset),k,PQRPARAM,work(work_offset),lwork,actualk)
 
        !print *,'possible rank:', actualk

        ! override useful for debugging
        !actualk = 1
        !actualk = k
        !actualk= min(actualk,2)
        DO ivector=1,actualk
           CALL tum_pdlarfgk_1dcomm_vector(a(1,k-ivector+1),1,idx, &
           tau(k-ivector+1), work(seedC_offset),work(seedD_offset),k, &
           ivector,m,mb,rev,mpicomm)

           CALL tum_pdlarfgk_1dcomm_update(a(1,1),lda,idx,work(work_offset), &
           lwork,work(seedC_offset), work(seedD_offset),k,actualk,ivector,tau,m,mb,rev,mpicomm)

            CALL tum_pdlarfg_copy_1dcomm(a(1,k-ivector+1),1, &
                                         v(1,k-ivector+1),1, &
                                         m,baseidx,idx-ivector+1,mb,1,mpicomm)
        END DO

        ! generate final T matrix and convert preliminary tau values into real ones
        CALL tum_pdlarfgk_1dcomm_generateT(work(seedC_offset), &
        work(seedD_offset),k,actualk,tau,t,ldt)

END SUBROUTINE tum_pdlarfgk_1dcomm

SUBROUTINE tum_pdlarfgk_1dcomm_seed(a,lda,baseidx,work,lwork,seedC,seedD,m,k,mb,mpicomm)
    USE ELPA1
    USE tum_utils
    INTEGER                                  :: lda, lwork

! parameter setup
! input variables (local)

    DOUBLE PRECISION a(lda,*), work(*)

    ! input variables (global)
    INTEGER m,k,baseidx,mb,mpicomm
    DOUBLE PRECISION seedC(k,*),seedD(k,*)
 
    ! output variables (global)

    ! derived input variables from TUM_PQRPARAM

    ! local scalars
    INTEGER mpierr,mpirank,mpiprocs,mpirank_top
    INTEGER icol,irow,lidx,remsize
    INTEGER remaining_rank

    INTEGER C_size,D_size,sendoffset,recvoffset,sendrecv_size
    INTEGER localoffset,localsize,baseoffset
 
    CALL MPI_Comm_rank(mpicomm, mpirank, mpierr)
    CALL MPI_Comm_size(mpicomm, mpiprocs, mpierr)

    C_size = k*k
    D_size = k*k
    sendoffset = 1
    sendrecv_size = C_size+D_size
    recvoffset = sendoffset + sendrecv_size

    IF (lwork .EQ. -1) THEN
        work(1) = DBLE(2*sendrecv_size)
        RETURN
    END IF
  
    ! clear buffer
    work(sendoffset:sendoffset+sendrecv_size-1)=0.0d0

    ! collect C part
    DO icol=1,k

        remaining_rank = k
        DO WHILE (remaining_rank .GT. 0)
            irow = k - remaining_rank + 1
            lidx = baseidx - remaining_rank + 1

            ! determine chunk where the current top element is located
            mpirank_top = MOD((lidx-1)/mb,mpiprocs) 

            ! limit max number of remaining elements of this chunk to the block
            ! distribution parameter
            remsize = MIN(remaining_rank,mb)

            ! determine the number of needed elements in this chunk 
            CALL local_size_offset_1d(lidx+remsize-1,mb, &
                                      lidx,lidx,0, &
                                      mpirank_top,mpiprocs, &
                                      localsize,baseoffset,localoffset)

            !print *,'local rank',localsize,localoffset

            IF (mpirank .EQ. mpirank_top) THEN
                ! copy elements to buffer
                work(sendoffset+(icol-1)*k+irow-1:sendoffset+(icol-1)*k+irow-1+localsize-1) &
                            = a(localoffset:localoffset+remsize-1,icol)
            END IF

            ! jump to next chunk
            remaining_rank = remaining_rank - localsize
        END DO
    END DO

    ! collect D part
        CALL local_size_offset_1d(m,mb,baseidx-k,baseidx-k,1, &
                                                          mpirank,mpiprocs, &
                                                          localsize,baseoffset,localoffset)
 
    !print *,'localsize',localsize,localoffset
    IF (localsize > 0) THEN
        CALL dsyrk("Upper", "Trans", k, localsize, &
                   1.0d0, a(localoffset,1), lda, &
                   0.0d0, work(sendoffset+C_size), k)
    ELSE
        work(sendoffset+C_size:sendoffset+C_size+k*k-1) = 0.0d0
    END IF
 
    ! TODO: store symmetric part more efficiently

    ! allreduce operation on results
    CALL mpi_allreduce(work(sendoffset),work(recvoffset),sendrecv_size, &
                       mpi_real8,mpi_sum,mpicomm,mpierr)

    ! unpack result from buffer into seedC and seedD
    seedC(1:k,1:k) = 0.0d0
    DO icol=1,k
        seedC(1:k,icol) = work(recvoffset+(icol-1)*k:recvoffset+icol*k-1)
    END DO
 
    seedD(1:k,1:k) = 0.0d0
    DO icol=1,k
        seedD(1:k,icol) = work(recvoffset+C_size+(icol-1)*k:recvoffset+C_size+icol*k-1)
    END DO
END SUBROUTINE tum_pdlarfgk_1dcomm_seed

! k is assumed to be larger than two
SUBROUTINE tum_pdlarfgk_1dcomm_check_improved(seedC,seedD,k,PQRPARAM,work,lwork,possiblerank)

    ! input variables (global)
    INTEGER                                  :: k, PQRPARAM(*), lwork

    DOUBLE PRECISION seedC(k,*),seedD(k,*),work(k,*)
 
    ! output variables (global)
    INTEGER possiblerank
 
    ! derived input variables from TUM_PQRPARAM
    INTEGER eps

    ! local variables
    INTEGER i,j,l
    DOUBLE PRECISION sum_squares,diagonal_square,relative_error,epsd,diagonal_root
    DOUBLE PRECISION dreverse_matrix_work(1)

    ! external functions
    DOUBLE PRECISION ddot,dlapy2,dnrm2
    EXTERNAL ddot,dscal,dlapy2,dnrm2

    IF (lwork .EQ. -1) THEN
        CALL reverse_matrix_local(1,k,k,work,k,dreverse_matrix_work,-1)
        work(1,1) = DBLE(k*k) + dreverse_matrix_work(1)
        RETURN
    END IF

    eps = PQRPARAM(3)

    IF (eps .EQ. 0) THEN 
        possiblerank = k
        RETURN
    END IF

    epsd = DBLE(eps)
 
    ! build complete inner product from seedC and seedD
    ! copy seedD to work
    work(:,1:k) = seedD(:,1:k)

    ! add inner products of seedC to work
    CALL dsyrk("Upper", "Trans", k, k, &
               1.0d0, seedC(1,1), k, &
               1.0d0, work, k)

        ! TODO: optimize this part!
        CALL reverse_matrix_local(0,k,k,work(1,1),k,work(1,k+1),lwork-2*k)
        CALL reverse_matrix_local(1,k,k,work(1,1),k,work(1,k+1),lwork-2*k)

    ! transpose matrix
        DO i=1,k
          DO j=i+1,k
            work(i,j) = work(j,i)
          END DO
        END DO


    ! do cholesky decomposition
    i = 0
    DO WHILE ((i .LT. k))
      i = i + 1
      
      diagonal_square = ABS(work(i,i))
      diagonal_root  = SQRT(diagonal_square)
      
      ! zero Householder vector (zero norm) case
      IF ((ABS(diagonal_square) .EQ. 0.0d0) .OR. (ABS(diagonal_root) .EQ. 0.0d0)) THEN
        possiblerank = MAX(i-1,1)
        RETURN
      END IF
    
      ! check if relative error is bounded for each Householder vector
      ! Householder i is stable iff Househoulder i-1 is "stable" and the accuracy criterion
      ! holds.
      ! first Householder vector is considered as "stable".
          
      DO j=i+1,k
          work(i,j) = work(i,j) / diagonal_root
          DO l=i+1,j
              work(l,j) = work(l,j) - work(i,j) * work(i,l)
          END DO
      END DO
      !print *,'cholesky step done'
    
      ! build sum of squares
      IF(i .EQ. 1) THEN
        sum_squares = 0.0d0
      ELSE
        sum_squares = ddot(i-1,work(1,i),1,work(1,i),1)
      END IF
      !relative_error = sum_squares / diagonal_square
      !print *,'error ',i,sum_squares,diagonal_square,relative_error
      
      IF (sum_squares .GE. (epsd * diagonal_square)) THEN
        possiblerank = MAX(i-1,1)
        RETURN
      END IF
    END DO

    possiblerank = i
    !print *,'possible rank', possiblerank
END SUBROUTINE tum_pdlarfgk_1dcomm_check_improved

! TODO: zero Householder vector (zero norm) case
! - check alpha values as well (from seedC)
SUBROUTINE tum_pdlarfgk_1dcomm_check(seedC,seedD,k,PQRPARAM,work,lwork,possiblerank)
    USE tum_utils
    INTEGER                                  :: k, PQRPARAM(*), lwork

! parameter setup
! input variables (local)
! input variables (global)

    DOUBLE PRECISION seedC(k,*),seedD(k,*),work(k,*)
 
    ! output variables (global)
    INTEGER possiblerank

    ! derived input variables from TUM_PQRPARAM
    INTEGER eps

    ! local scalars
    INTEGER icol,isqr,iprod
    DOUBLE PRECISION epsd,sum_sqr,sum_products,diff,temp,ortho,ortho_sum
    DOUBLE PRECISION dreverse_matrix_work(1)

    IF (lwork .EQ. -1) THEN
        CALL reverse_matrix_local(1,k,k,work,k,dreverse_matrix_work,-1)
        work(1,1) = DBLE(k*k) + dreverse_matrix_work(1)
        RETURN
    END IF

    eps = PQRPARAM(3)

    IF (eps .EQ. 0) THEN 
        possiblerank = k
        RETURN
    END IF

    epsd = DBLE(eps)


    ! copy seedD to work
    work(:,1:k) = seedD(:,1:k)

    ! add inner products of seedC to work
    CALL dsyrk("Upper", "Trans", k, k, &
               1.0d0, seedC(1,1), k, &
               1.0d0, work, k)

        ! TODO: optimize this part!
        CALL reverse_matrix_local(0,k,k,work(1,1),k,work(1,k+1),lwork-2*k)
        CALL reverse_matrix_local(1,k,k,work(1,1),k,work(1,k+1),lwork-2*k)

        ! transpose matrix
        DO icol=1,k
                DO isqr=icol+1,k
                        work(icol,isqr) = work(isqr,icol)
                END DO
        END DO

    ! work contains now the full inner product of the global (sub-)matrix
    DO icol=1,k
        ! zero Householder vector (zero norm) case
        IF (ABS(work(icol,icol)) .EQ. 0.0d0) THEN
            !print *,'too small ', icol, work(icol,icol)
            possiblerank = MAX(icol,1)
            RETURN
        END IF

        sum_sqr = 0.0d0
        DO isqr=1,icol-1
            sum_products = 0.0d0
            DO iprod=1,isqr-1
                sum_products = sum_products + work(iprod,isqr)*work(iprod,icol)
            END DO

            !print *,'divisor',icol,isqr,work(isqr,isqr)
            temp = (work(isqr,icol) - sum_products)/work(isqr,isqr)
            work(isqr,icol) = temp
            sum_sqr = sum_sqr + temp*temp
        END DO

        ! calculate diagonal value
        diff = work(icol,icol) - sum_sqr
        IF (diff .LT. 0.0d0) THEN
            ! we definitely have a problem now
            possiblerank = icol-1 ! only decompose to previous column (including)
            RETURN
        END IF
        work(icol,icol) = SQRT(diff)

        ! calculate orthogonality
        ortho = 0.0d0
        DO isqr=1,icol-1
            ortho_sum = 0.0d0
            DO iprod=isqr,icol-1
                temp = work(isqr,iprod)*work(isqr,iprod)
                !print *,'ortho ', work(iprod,iprod)
                temp = temp / (work(iprod,iprod)*work(iprod,iprod))
                ortho_sum = ortho_sum + temp
            END DO
            ortho = ortho + ortho_sum * (work(isqr,icol)*work(isqr,icol))
        END DO
 
        ! ---------------- with division by zero ----------------------- !

        !ortho = ortho / diff;

        ! if current estimate is not accurate enough, the following check holds
        !if (ortho .gt. epsd) then
        !    possiblerank = icol-1 ! only decompose to previous column (including)
        !    return
        !end if
 
        ! ---------------- without division by zero ----------------------- !

        ! if current estimate is not accurate enough, the following check holds
        IF (ortho .GT. epsd * diff) THEN
            possiblerank = icol-1 ! only decompose to previous column (including)
            RETURN
        END IF
    END DO

    ! if we get to this point, the accuracy condition holds for the whole block
    possiblerank = k
END SUBROUTINE tum_pdlarfgk_1dcomm_check

!sidx: seed idx
!k: max rank used during seed phase
!rank: actual rank (k >= rank)
SUBROUTINE tum_pdlarfgk_1dcomm_vector(x,incx,baseidx,tau,seedC,seedD,k,sidx,n,nb,rev,mpicomm)
    USE ELPA1
    USE tum_utils
    INTEGER                                  :: incx

! input variables (local)

    DOUBLE PRECISION x(*),tau

    ! input variables (global)
    INTEGER n,nb,baseidx,rev,mpicomm,k,sidx
    DOUBLE PRECISION seedC(k,*),seedD(k,*)
 
    ! output variables (global)

    ! external functions
    DOUBLE PRECISION dlapy2,dnrm2
    EXTERNAL dlapy2,dscal,dnrm2

    ! local scalars
    INTEGER mpirank,mpirank_top,mpiprocs,mpierr
    DOUBLE PRECISION alpha,dot,beta,xnorm
    INTEGER local_size,baseoffset,local_offset,top,topidx
    INTEGER lidx

    CALL MPI_Comm_rank(mpicomm, mpirank, mpierr)
    CALL MPI_Comm_size(mpicomm, mpiprocs, mpierr)

        lidx = baseidx-sidx+1
        CALL local_size_offset_1d(n,nb,baseidx,lidx-1,rev,mpirank,mpiprocs, &
                                                          local_size,baseoffset,local_offset)
 
    local_offset = local_offset * incx

    ! Processor id for global index of top element
    mpirank_top = MOD((lidx-1)/nb,mpiprocs)
    IF (mpirank .EQ. mpirank_top) THEN
        topidx = local_index((lidx),mpirank_top,mpiprocs,nb,0)
        top = 1+(topidx-1)*incx
    END IF

        alpha = seedC(k-sidx+1,k-sidx+1)
        dot = seedD(k-sidx+1,k-sidx+1)
        ! assemble actual norm from both seed parts
        xnorm = dlapy2(SQRT(dot), dnrm2(k-sidx,seedC(1,k-sidx+1),1))
    
    IF (xnorm .EQ. 0.0d0) THEN
        tau = 0.0d0
    ELSE
        ! General case

        beta = SIGN(dlapy2(alpha, xnorm), alpha)
        ! store a preliminary version of beta in tau
        tau = beta
        
        ! update global part
        CALL dscal(local_size, 1.0d0/(beta+alpha), &
                   x(local_offset), incx)

        ! do not update local part here due to
        ! dependency of c vector during update process

        ! TODO: reimplement norm rescale method of 
        ! original PDLARFG using mpi?

        IF (mpirank .EQ. mpirank_top) THEN
            x(top) = -beta
        END IF
    END IF
  
END SUBROUTINE tum_pdlarfgk_1dcomm_vector

!k: original max rank used during seed function
!rank: possible rank as from check function
! TODO: if rank is less than k, reduce buffersize in such a way
! that only the required entries for the next pdlarfg steps are
! computed
SUBROUTINE tum_pdlarfgk_1dcomm_update(a,lda,baseidx,work,lwork,seedC,seedD,k,rank,sidx,tau,n,nb,rev,mpicomm)
    USE ELPA1
    USE tum_utils
    INTEGER                                  :: eps_, gmode_, rank_, upmode1_

! parameter setup

    PARAMETER   (gmode_ = 1,rank_ = 2,eps_=3, upmode1_=4)

    ! input variables (local)
    INTEGER lda,lwork
    DOUBLE PRECISION a(lda,*),work(*)

    ! input variables (global)
    INTEGER k,rank,sidx,n,baseidx,nb,rev,mpicomm
    DOUBLE PRECISION beta
 
    ! output variables (global)
    DOUBLE PRECISION seedC(k,*),seedD(k,*),tau(*)

    ! derived input variables from TUM_PQRPARAM

    ! local scalars
    DOUBLE PRECISION alpha
    INTEGER coffset,zoffset,yoffset,voffset,buffersize
    INTEGER mpirank,mpierr,mpiprocs,mpirank_top
    INTEGER localsize,baseoffset,localoffset,topidx
    INTEGER lidx

    IF (lwork .EQ. -1) THEN
        ! buffer for c,z,y,v
        work(1) = 4*k
        RETURN
    END IF

    ! nothing to update anymore
    IF (sidx .GT. rank) RETURN

    CALL MPI_Comm_rank(mpicomm, mpirank, mpierr)
    CALL MPI_Comm_size(mpicomm, mpiprocs, mpierr)

        lidx = baseidx-sidx
        IF (lidx .LT. 1) RETURN

    CALL local_size_offset_1d(n,nb,baseidx,lidx,rev,mpirank,mpiprocs, &
                              localsize,baseoffset,localoffset)

    coffset = 1
    zoffset = coffset + k
    yoffset = zoffset + k
    voffset = yoffset + k
    buffersize = k - sidx

    ! finalize tau values
        alpha = seedC(k-sidx+1,k-sidx+1)
        beta = tau(k-sidx+1)

    ! zero Householder vector (zero norm) case
    !print *,'k update: alpha,beta',alpha,beta
    IF ((beta .EQ. 0.0d0) .OR. (alpha .EQ. 0.0d0))  THEN
        tau(k-sidx+1) = 0.0d0
        seedC(k,k-sidx+1) = 0.0d0
        RETURN
    END IF
            
    tau(k-sidx+1) = (beta+alpha) / beta

    ! ---------------------------------------
        ! calculate c vector (extra vector or encode in seedC/seedD?
        work(coffset:coffset+buffersize-1) = seedD(1:buffersize,k-sidx+1)
        CALL dgemv("Trans", buffersize+1, buffersize, &
               1.0d0,seedC(1,1),k,seedC(1,k-sidx+1),1, &
               1.0d0,work(coffset),1)

        ! calculate z using tau,seedD,seedC and c vector
        work(zoffset:zoffset+buffersize-1) = seedC(k-sidx+1,1:buffersize)
        CALL daxpy(buffersize, 1.0d0/beta, work(coffset), 1, work(zoffset), 1)
 
        ! update A1(local copy) and generate part of householder vectors for use
        CALL daxpy(buffersize, -1.0d0, work(zoffset),1,seedC(k-sidx+1,1),k)
        CALL dscal(buffersize, 1.0d0/(alpha+beta), seedC(1,k-sidx+1),1)
        CALL dger(buffersize, buffersize, -1.0d0, seedC(1,k-sidx+1),1, work(zoffset), 1, seedC(1,1), k)
 
        ! update A global (householder vector already generated by pdlarfgk)
        mpirank_top = MOD(lidx/nb,mpiprocs)
        IF (mpirank .EQ. mpirank_top) THEN
            ! handle first row separately
            topidx = local_index(lidx+1,mpirank_top,mpiprocs,nb,0)
            CALL daxpy(buffersize,-1.0d0,work(zoffset),1,a(topidx,1),lda)
        END IF

        CALL dger(localsize, buffersize,-1.0d0, &
              a(localoffset,k-sidx+1),1,work(zoffset),1, &
              a(localoffset,1),lda)
        
        ! update D (symmetric) => two buffer vectors of size rank
        ! generate y vector
        work(yoffset:yoffset+buffersize-1) = 0.d0
        CALL daxpy(buffersize,1.0d0/(alpha+beta),work(zoffset),1,work(yoffset),1)

        ! generate v vector
        work(voffset:voffset+buffersize-1) = seedD(1:buffersize,k-sidx+1)
        CALL daxpy(buffersize, -0.5d0*seedD(k-sidx+1,k-sidx+1), work(yoffset), 1, work(voffset),1)

        ! symmetric update of D using y and v
        CALL dsyr2("Upper", buffersize,-1.0d0, &
                   work(yoffset),1,work(voffset),1, &
                   seedD(1,1), k)
  
    ! prepare T matrix inner products
    ! D_k(1:k,k+1:n) = D_(k-1)(1:k,k+1:n) - D_(k-1)(1:k,k) * y'
    ! store coefficient 1.0d0/(alpha+beta) in C diagonal elements
        CALL dger(k-sidx,sidx,-1.0d0,work(yoffset),1,seedD(k-sidx+1,k-sidx+1),k,seedD(1,k-sidx+1),k)
        seedC(k,k-sidx+1) = 1.0d0/(alpha+beta)

END SUBROUTINE tum_pdlarfgk_1dcomm_update


SUBROUTINE tum_pdlarfgk_1dcomm_generateT(seedC,seedD,k,actualk,tau,t,ldt)

    INTEGER                                  :: k, actualk, ldt

    DOUBLE PRECISION seedC(k,*),seedD(k,*),tau(*),t(ldt,*)

    INTEGER irow,icol
    DOUBLE PRECISION column_coefficient

        !print *,'reversed on the fly T generation NYI'

        DO icol=1,actualk-1
            ! calculate inner product of householder vector parts in seedC
            ! (actually calculating more than necessary, if actualk < k)
            ! => a lot of junk from row 1 to row k-actualk
            CALL dtrmv('Upper','Trans','Unit',k-icol,seedC(1,1),k,seedC(1,k-icol+1),1)
        
            ! add scaled D parts to current column of C (will become later T rows)
            column_coefficient = seedC(k,k-icol+1)
            DO irow=k-actualk+1,k-1
                seedC(irow,k-icol+1) = ( seedC(irow,k-icol+1) ) +  ( seedD(irow,k-icol+1) * column_coefficient * seedC(k,irow) )
            END DO
        END DO
 
        CALL tum_dlarft_kernel(actualk,tau(k-actualk+1),seedC(k-actualk+1,k-actualk+2),k,t(k-actualk+1,k-actualk+1),ldt)

END SUBROUTINE tum_pdlarfgk_1dcomm_generateT

!direction=0: pack into work buffer
!direction=1: unpack from work buffer
SUBROUTINE tum_pdgeqrf_pack_unpack(v,ldv,work,lwork,m,n,mb,baseidx,rowidx,rev,direction,mpicomm)
    USE ELPA1
    USE tum_utils
    INTEGER                                  :: ldv, lwork

! input variables (local)

    DOUBLE PRECISION v(ldv,*), work(*)

    ! input variables (global)
    INTEGER m,n,mb,baseidx,rowidx,rev,direction,mpicomm
 
    ! output variables (global)
 
    ! local scalars
    INTEGER mpierr,mpirank,mpiprocs
    INTEGER buffersize,icol
    INTEGER local_size,baseoffset,offset

    ! external functions
 
    CALL mpi_comm_rank(mpicomm,mpirank,mpierr)
    CALL mpi_comm_size(mpicomm,mpiprocs,mpierr)
  
    CALL local_size_offset_1d(m,mb,baseidx,rowidx,rev,mpirank,mpiprocs, &
                                  local_size,baseoffset,offset)

    !print *,'pack/unpack',local_size,baseoffset,offset

    ! rough approximate for buffer size
    IF (lwork .EQ. -1) THEN
        buffersize = local_size * n ! vector elements
        work(1) = DBLE(buffersize)
        RETURN
    END IF

    IF (direction .EQ. 0) THEN
        ! copy v part to buffer (including zeros)
        DO icol=1,n
            work(1+local_size*(icol-1):local_size*icol) = v(baseoffset:baseoffset+local_size-1,icol)
        END DO
    ELSE
        ! copy v part from buffer (including zeros)
        DO icol=1,n
            v(baseoffset:baseoffset+local_size-1,icol) = work(1+local_size*(icol-1):local_size*icol)
        END DO
    END IF

    RETURN

END SUBROUTINE tum_pdgeqrf_pack_unpack

!direction=0: pack into work buffer
!direction=1: unpack from work buffer
SUBROUTINE tum_pdgeqrf_pack_unpack_tmatrix(tau,t,ldt,work,lwork,n,direction)
    USE ELPA1
    USE tum_utils
    INTEGER                                  :: ldt, lwork

! input variables (local)

    DOUBLE PRECISION work(*), t(ldt,*),tau(*)

    ! input variables (global)
    INTEGER n,direction
 
    ! output variables (global)
 
    ! local scalars
    INTEGER icol

    ! external functions
 
    IF (lwork .EQ. -1) THEN
        work(1) = DBLE(n*n)
        RETURN
    END IF

    IF (direction .EQ. 0) THEN
        ! append t matrix to buffer (including zeros)
        DO icol=1,n
            work(1+(icol-1)*n:icol*n) = t(1:n,icol)
        END DO
    ELSE
        ! append t matrix from buffer (including zeros)
        DO icol=1,n
            t(1:n,icol) = work(1+(icol-1)*n:icol*n)
            tau(icol) = t(icol,icol)
        END DO
    END IF

END SUBROUTINE tum_pdgeqrf_pack_unpack_tmatrix


! TODO: encode following functionality
!   - Direction? BOTTOM UP or TOP DOWN ("Up", "Down")
!        => influences all related kernels (including DLARFT / DLARFB)
!   - rank-k parameter (k=1,2,...,b)
!        => influences possible update strategies
!        => parameterize the function itself? (FUNCPTR, FUNCARG)
!   - Norm mode? Allreduce, Allgather, AlltoAll, "AllHouse", (ALLNULL = benchmarking local kernels)
!   - subblocking 
!         (maximum block size bounded by data distribution along rows)
!   - blocking method (householder vectors only or compact WY?)
!   - update strategy of trailing parts (incremental, complete) 
!        - difference for subblocks and normal blocks? (UPDATE and UPDATESUB)
!        o "Incremental"
!        o "Full"
!   - final T generation (recursive: subblock wise, block wise, end) (TMERGE)
!        ' (implicitly given by / influences update strategies?)
!        => alternative: during update: iterate over sub t parts
!           => advantage: smaller (cache aware T parts)
!           => disadvantage: more memory write backs 
!                (number of T parts * matrix elements)
!   - partial/sub T generation (TGEN)
!        o add vectors right after creation (Vector)
!        o add set of vectors (Set)
!   - bcast strategy of householder vectors to other process columns 
!        (influences T matrix generation and trailing update 
!         in other process columns)
!        o no broadcast (NONE = benchmarking?, 
!            or not needed due to 1D process grid)
!        o after every housegen (VECTOR)
!        o after every subblk   (SUBBLOCK)
!        o after full local column block decomposition (BLOCK)
!  LOOP Housegen -> BCAST -> GENT/EXTENDT -> LOOP HouseLeft

!subroutine tum_pqrparam_init(PQRPARAM, DIRECTION, RANK, NORMMODE, &
!                             SUBBLK, UPDATE, TGEN, BCAST)

! gmode: control communication pattern of dlarfg
! maxrank: control max number of householder vectors per communication
! eps: error threshold (integer)
! update*: control update pattern in pdgeqr2_1dcomm ('incremental','full','merge')
!               merging = full update with tmatrix merging
! tmerge*: 0: do not merge, 1: incremental merge, >1: recursive merge
!               only matters if update* == full
SUBROUTINE tum_pqrparam_init(pqrparam,size2d,update2d,tmerge2d,size1d,update1d,tmerge1d,maxrank,update,eps,hgmode)


    ! input
    INTEGER                                  :: PQRPARAM(*), size2d
    CHARACTER                                :: update2d
    INTEGER                                  :: tmerge2d, size1d
    CHARACTER                                :: update1d
    INTEGER                                  :: tmerge1d, maxrank
    CHARACTER                                :: update
    INTEGER                                  :: eps
    CHARACTER                                :: hgmode

! output

    PQRPARAM(1) = size2d
    PQRPARAM(2) = ICHAR(update2d)
    PQRPARAM(3) = tmerge2d
    ! TODO: broadcast T yes/no

    PQRPARAM(4) = size1d
    PQRPARAM(5) = ICHAR(update1d)
    PQRPARAM(6) = tmerge1d

    PQRPARAM(7) = maxrank
    PQRPARAM(8) = ICHAR(update)
    PQRPARAM(9) = eps
    PQRPARAM(10) = ICHAR(hgmode)

END SUBROUTINE tum_pqrparam_init


SUBROUTINE tum_pdlarfg_copy_1dcomm(x,incx,v,incv,n,baseidx,idx,nb,rev,mpicomm)
    USE ELPA1
    USE tum_utils
    INTEGER                                  :: incx, incv

! input variables (local)

    DOUBLE PRECISION x(*), v(*)

    ! input variables (global)
    INTEGER baseidx,idx,rev,nb,n
    INTEGER mpicomm
 
    ! output variables (global)
 
    ! local scalars
    INTEGER mpierr,mpiprocs
    INTEGER mpirank,mpirank_top
    INTEGER irow,x_offset
    INTEGER v_offset,local_size


    CALL MPI_Comm_rank(mpicomm, mpirank, mpierr)
    CALL MPI_Comm_size(mpicomm, mpiprocs, mpierr)
 
    CALL local_size_offset_1d(n,nb,baseidx,idx,rev,mpirank,mpiprocs, &
                              local_size,v_offset,x_offset)
    v_offset = v_offset * incv

    !print *,'copy:',mpirank,baseidx,v_offset,x_offset,local_size

    ! copy elements
    DO irow=1,local_size
        v((irow-1)*incv+v_offset) = x((irow-1)*incx+x_offset)
    END DO

    ! replace top element to build an unitary vector
    mpirank_top = MOD((idx-1)/nb,mpiprocs)
    IF (mpirank .EQ. mpirank_top) THEN
        v(local_size*incv) = 1.0d0
    END IF

END SUBROUTINE tum_pdlarfg_copy_1dcomm


#endif    
END MODULE elpa_pdgeqrf
