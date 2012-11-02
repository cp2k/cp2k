! ELPA2 -- 2-stage solver for ELPA
! 
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".

MODULE elpa2

! Version 1.1.2, 2011-02-21

  USE elpa1
  USE elpa_pdgeqrf

  IMPLICIT NONE

  PRIVATE ! By default, all routines contained are private

#if defined(__parallel)
  INCLUDE "mpif.h"

  ! The following routines are public:

  PUBLIC :: solve_evp_real_2stage
  PUBLIC :: solve_evp_complex_2stage

  PUBLIC :: bandred_real
  PUBLIC :: tridiag_band_real
  PUBLIC :: trans_ev_tridi_to_band_real
  PUBLIC :: trans_ev_band_to_full_real

  PUBLIC :: bandred_complex
  PUBLIC :: tridiag_band_complex
  PUBLIC :: trans_ev_tridi_to_band_complex
  PUBLIC :: trans_ev_band_to_full_complex
  
  PUBLIC :: band_band_real
  PUBLIC :: divide_band
  
!-------------------------------------------------------------------------------  

  INTEGER, PUBLIC :: which_qr_decomposition = 1     ! defines, which QR-decomposition algorithm will be used
                                                    ! 0 for unblocked
                                                    ! 1 for blocked (maxrank: nblk)

!-------------------------------------------------------------------------------

  ! The following array contains the Householder vectors of the
  ! transformation band -> tridiagonal.
  ! It is allocated and set in tridiag_band_real and used in
  ! trans_ev_tridi_to_band_real.
  ! It must be deallocated by the user after trans_ev_tridi_to_band_real!

  REAL*8, ALLOCATABLE :: hh_trans_real(:,:)
  COMPLEX*16, ALLOCATABLE :: hh_trans_complex(:,:)

!-------------------------------------------------------------------------------

!******
CONTAINS

SUBROUTINE solve_evp_real_2stage(na, nev, a, lda, ev, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)

!-------------------------------------------------------------------------------
!  solve_evp_real_2stage: Solves the real eigenvalue problem with a 2 stage approach
!
!
!  na          Order of matrix a
!
!  nev         Number of eigenvalues needed
!
!  a(lda,*)    Distributed matrix for which eigenvalues are to be computed.
!              Distribution is like in Scalapack.
!              The full matrix must be set (not only one half like in scalapack).
!              Destroyed on exit (upper and lower half).
!
!  lda         Leading dimension of a
!
!  ev(na)      On output: eigenvalues of a, every processor gets the complete set
!
!  q(ldq,*)    On output: Eigenvectors of a
!              Distribution is like in Scalapack.
!              Must be always dimensioned to the full size (corresponding to (na,na))
!              even if only a part of the eigenvalues is needed.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm_all
!              MPI-Communicator for the total processor set
!
!-------------------------------------------------------------------------------


    INTEGER, INTENT(in)                      :: na, nev, lda
    REAL*8, INTENT(inout)                    :: a(lda,*), ev(na)
    INTEGER, INTENT(in)                      :: ldq
    REAL*8, INTENT(inout)                    :: q(ldq,*)
    INTEGER, INTENT(in)                      :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols, mpi_comm_all

    INTEGER                                  :: mpierr, my_pcol, my_pe, &
                                                my_prow, n_pes, nbw, np_cols, &
                                                np_rows, num_blocks
    REAL*8                                   :: ttt0, ttt1, ttts
    REAL*8, ALLOCATABLE                      :: e(:), tmat(:,:,:)

   CALL mpi_comm_rank(mpi_comm_all,my_pe,mpierr)
   CALL mpi_comm_size(mpi_comm_all,n_pes,mpierr)

   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Choose bandwidth, must be a multiple of nblk, set to a value >= 32

   nbw = (31/nblk+1)*nblk

   num_blocks = (na-1)/nbw + 1

   ALLOCATE(tmat(nbw,nbw,num_blocks))

   ! Reduction full -> band

   ttt0 = MPI_Wtime()
   ttts = ttt0
   CALL bandred_real(na, a, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols, tmat)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) &
      PRINT 1,'Time bandred_real               :',ttt1-ttt0

   ! Reduction band -> tridiagonal

   ALLOCATE(e(na))

   ttt0 = MPI_Wtime()
   CALL tridiag_band_real(na, nbw, nblk, a, lda, ev, e, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) &
      PRINT 1,'Time tridiag_band_real          :',ttt1-ttt0

   CALL mpi_bcast(ev,na,MPI_REAL8,0,mpi_comm_all,mpierr)
   CALL mpi_bcast(e,na,MPI_REAL8,0,mpi_comm_all,mpierr)

   ttt1 = MPI_Wtime()
   time_evp_fwd = ttt1-ttts

   ! Solve tridiagonal system

   ttt0 = MPI_Wtime()
   CALL solve_tridi(na, nev, ev, e, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) &
      PRINT 1,'Time solve_tridi                :',ttt1-ttt0
   time_evp_solve = ttt1-ttt0
   ttts = ttt1

   DEALLOCATE(e)

   ! Backtransform stage 1

   ttt0 = MPI_Wtime()
   CALL trans_ev_tridi_to_band_real(na, nev, nblk, nbw, q, ldq, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) &
      PRINT 1,'Time trans_ev_tridi_to_band_real:',ttt1-ttt0

   ! We can now deallocate the stored householder vectors
   DEALLOCATE(hh_trans_real)

   ! Backtransform stage 2

   ttt0 = MPI_Wtime()
   CALL trans_ev_band_to_full_real(na, nev, nblk, nbw, a, lda, tmat, q, ldq, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) &
      PRINT 1,'Time trans_ev_band_to_full_real :',ttt1-ttt0
   time_evp_back = ttt1-ttts

   DEALLOCATE(tmat)

1  FORMAT(a,f10.3)

END SUBROUTINE solve_evp_real_2stage

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------

SUBROUTINE solve_evp_complex_2stage(na, nev, a, lda, ev, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)

!-------------------------------------------------------------------------------
!  solve_evp_complex_2stage: Solves the complex eigenvalue problem with a 2 stage approach
!
!
!  na          Order of matrix a
!
!  nev         Number of eigenvalues needed
!
!  a(lda,*)    Distributed matrix for which eigenvalues are to be computed.
!              Distribution is like in Scalapack.
!              The full matrix must be set (not only one half like in scalapack).
!              Destroyed on exit (upper and lower half).
!
!  lda         Leading dimension of a
!
!  ev(na)      On output: eigenvalues of a, every processor gets the complete set
!
!  q(ldq,*)    On output: Eigenvectors of a
!              Distribution is like in Scalapack.
!              Must be always dimensioned to the full size (corresponding to (na,na))
!              even if only a part of the eigenvalues is needed.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm_all
!              MPI-Communicator for the total processor set
!
!-------------------------------------------------------------------------------


    INTEGER, INTENT(in)                      :: na, nev, lda
    COMPLEX*16, INTENT(inout)                :: a(lda,*)
    REAL*8, INTENT(inout)                    :: ev(na)
    INTEGER, INTENT(in)                      :: ldq
    COMPLEX*16, INTENT(inout)                :: q(ldq,*)
    INTEGER, INTENT(in)                      :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols, mpi_comm_all

    COMPLEX*16, ALLOCATABLE                  :: tmat(:,:,:)
    INTEGER                                  :: l_cols, l_cols_nev, l_rows, &
                                                mpierr, my_pcol, my_prow, &
                                                nbw, np_cols, np_rows, &
                                                num_blocks
    REAL*8                                   :: ttt0, ttt1, ttts
    REAL*8, ALLOCATABLE                      :: e(:), q_real(:,:)

   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Choose bandwidth, must be a multiple of nblk, set to a value >= 32

   nbw = (31/nblk+1)*nblk

   num_blocks = (na-1)/nbw + 1

   ALLOCATE(tmat(nbw,nbw,num_blocks))

   ! Reduction full -> band

   ttt0 = MPI_Wtime()
   ttts = ttt0
   CALL bandred_complex(na, a, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols, tmat)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) &
      PRINT 1,'Time bandred_complex               :',ttt1-ttt0

   ! Reduction band -> tridiagonal

   ALLOCATE(e(na))

   ttt0 = MPI_Wtime()
   CALL tridiag_band_complex(na, nbw, nblk, a, lda, ev, e, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) &
      PRINT 1,'Time tridiag_band_complex          :',ttt1-ttt0

   CALL mpi_bcast(ev,na,MPI_REAL8,0,mpi_comm_all,mpierr)
   CALL mpi_bcast(e,na,MPI_REAL8,0,mpi_comm_all,mpierr)

   ttt1 = MPI_Wtime()
   time_evp_fwd = ttt1-ttts

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q
   l_cols_nev = local_index(nev, my_pcol, np_cols, nblk, -1) ! Local columns corresponding to nev

   ALLOCATE(q_real(l_rows,l_cols))

   ! Solve tridiagonal system

   ttt0 = MPI_Wtime()
   CALL solve_tridi(na, nev, ev, e, q_real, UBOUND(q_real,1), nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times)  &
      PRINT 1,'Time solve_tridi                   :',ttt1-ttt0
   time_evp_solve = ttt1-ttt0
   ttts = ttt1

   q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)

   DEALLOCATE(e, q_real)

   ! Backtransform stage 1

   ttt0 = MPI_Wtime()
   CALL trans_ev_tridi_to_band_complex(na, nev, nblk, nbw, q, ldq, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) &
      PRINT 1,'Time trans_ev_tridi_to_band_complex:',ttt1-ttt0

   ! We can now deallocate the stored householder vectors
   DEALLOCATE(hh_trans_complex)

   ! Backtransform stage 2

   ttt0 = MPI_Wtime()
   CALL trans_ev_band_to_full_complex(na, nev, nblk, nbw, a, lda, tmat, q, ldq, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) &
      PRINT 1,'Time trans_ev_band_to_full_complex :',ttt1-ttt0
   time_evp_back = ttt1-ttts

   DEALLOCATE(tmat)

1  FORMAT(a,f10.3)

END SUBROUTINE solve_evp_complex_2stage

!-------------------------------------------------------------------------------

SUBROUTINE bandred_real(na, a, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols, tmat)

!-------------------------------------------------------------------------------
!  bandred_real: Reduces a distributed symmetric matrix to band form
!
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to Scalapack, a(:,:) must be set completely (upper and lower half)
!              a(:,:) is overwritten on exit with the band and the Householder vectors
!              in the upper half.
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith of output matrix
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  tmat(nbw,nbw,num_blocks)    where num_blocks = (na-1)/nbw + 1
!              Factors for the Householder vectors (returned), needed for back transformation
!
!-------------------------------------------------------------------------------


    INTEGER                                  :: na, lda
    REAL*8                                   :: a(lda,*)
    INTEGER                                  :: nblk, nbw, mpi_comm_rows, &
                                                mpi_comm_cols
    REAL*8                                   :: tmat(nbw,nbw,*)

    INTEGER :: cur_pcol, i, istep, j, l_cols, l_cols_tile, l_rows, &
      l_rows_tile, lc, lce, lch, lcs, lcx, lr, lre, mpierr, my_pcol, my_prow, &
      n_cols, ncol, nlc, np_cols, np_rows, nrow, pcol, PQRPARAM(11), prow, &
      tile_size, work_size
    REAL*8                                   :: aux1(nbw), aux2(nbw), &
                                                dwork_size(1), tau, &
                                                vav(nbw,nbw), vnorm2, vrl, xf
    REAL*8, ALLOCATABLE                      :: blockheuristic(:), &
                                                tauvector(:), tmp(:,:), &
                                                umc(:,:), vmr(:,:), vr(:), &
                                                work_blocked(:)

! needed for blocked QR decomposition

   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Semibandwith nbw must be a multiple of blocksize nblk

   IF(MOD(nbw,nblk)/=0) THEN
      IF(my_prow==0 .AND. my_pcol==0) THEN
         PRINT *,'ERROR: nbw=',nbw,', nblk=',nblk
         PRINT *,'ELPA2 works only for nbw==n*nblk'
         CALL mpi_abort(mpi_comm_world,0,mpierr)
      ENDIF
   ENDIF

   ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

   tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
   tile_size = ((128*MAX(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

   l_rows_tile = tile_size/np_rows ! local rows of a tile
   l_cols_tile = tile_size/np_cols ! local cols of a tile
 
   IF (which_qr_decomposition == 1) THEN
      CALL tum_pqrparam_init(pqrparam,    nblk,'M',0,   nblk,'M',0,   nblk,'M',1,'s')
      
      ALLOCATE(tauvector(na))
      ALLOCATE(blockheuristic(nblk))
      l_rows = local_index(na, my_prow, np_rows, nblk, -1)
      CALL tum_pdgeqrf_2dcomm(a, lda, vmr, MAX(l_rows,1), tauvector(1), &
      tmat(1,1,1), nbw, dwork_size(1), -1, na, nbw, nblk, nblk, na, na, &
      1, 0, PQRPARAM, mpi_comm_rows, mpi_comm_cols, blockheuristic)
      work_size = dwork_size(1)
      ALLOCATE(work_blocked(work_size))
      
      work_blocked = 0.0d0
   ENDIF

   DO istep = (na-1)/nbw, 1, -1

      n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

      ! Number of local columns/rows of remaining matrix
      l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
      l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)

      ! Allocate vmr and umc to their exact sizes so that they can be used in bcasts and reduces

      ALLOCATE(vmr(MAX(l_rows,1),2*n_cols))
      ALLOCATE(umc(MAX(l_cols,1),2*n_cols))

      ALLOCATE(vr(l_rows+1))

      vmr(1:l_rows,1:n_cols) = 0.
      vr(:) = 0
      tmat(:,:,istep) = 0

      ! Reduce current block to lower triangular form
      IF (which_qr_decomposition == 1) THEN
         CALL tum_pdgeqrf_2dcomm(a, lda, vmr, MAX(l_rows,1), tauvector(1), &
         tmat(1,1,istep), nbw, work_blocked, work_size, na, n_cols, nblk, &
         nblk, istep*nbw+n_cols-nbw, istep*nbw+n_cols, 1, 0, PQRPARAM, &
         mpi_comm_rows, mpi_comm_cols, blockheuristic)
      ELSE

      DO lc = n_cols, 1, -1

         ncol = istep*nbw + lc ! absolute column number of householder vector
         nrow = ncol - nbw ! Absolute number of pivot row

         lr  = local_index(nrow, my_prow, np_rows, nblk, -1) ! current row length
         lch = local_index(ncol, my_pcol, np_cols, nblk, -1) ! HV local column number

         tau = 0

         IF(nrow == 1) EXIT ! Nothing to do

         cur_pcol = pcol(ncol) ! Processor column owning current block

         IF(my_pcol==cur_pcol) THEN

            ! Get vector to be transformed; distribute last element and norm of
            ! remaining elements to all procs in current column

            vr(1:lr) = a(1:lr,lch) ! vector to be transformed

            IF(my_prow==prow(nrow)) THEN
               aux1(1) = DOT_PRODUCT(vr(1:lr-1),vr(1:lr-1))
               aux1(2) = vr(lr)
            ELSE
               aux1(1) = DOT_PRODUCT(vr(1:lr),vr(1:lr))
               aux1(2) = 0.
            ENDIF

            CALL mpi_allreduce(aux1,aux2,2,MPI_REAL8,MPI_SUM,mpi_comm_rows, &
            mpierr)

            vnorm2 = aux2(1)
            vrl    = aux2(2)

            ! Householder transformation

            CALL hh_transform_real(vrl, vnorm2, xf, tau)

            ! Scale vr and store Householder vector for back transformation

            vr(1:lr) = vr(1:lr) * xf
            IF(my_prow==prow(nrow)) THEN
               a(1:lr-1,lch) = vr(1:lr-1)
               a(lr,lch) = vrl
               vr(lr) = 1.
            ELSE
               a(1:lr,lch) = vr(1:lr)
            ENDIF

         ENDIF

         ! Broadcast Householder vector and tau along columns

         vr(lr+1) = tau
         CALL MPI_Bcast(vr,lr+1,MPI_REAL8,cur_pcol,mpi_comm_cols,mpierr)
         vmr(1:lr,lc) = vr(1:lr)
         tau = vr(lr+1)
         tmat(lc,lc,istep) = tau ! Store tau in diagonal of tmat

         ! Transform remaining columns in current block with Householder vector

         ! Local dot product

         aux1 = 0

         nlc = 0 ! number of local columns
         DO j=1,lc-1
            lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
            IF(lcx>0) THEN
               nlc = nlc+1
               IF(lr>0) aux1(nlc) = DOT_PRODUCT(vr(1:lr),a(1:lr,lcx))
            ENDIF
         ENDDO

         ! Get global dot products
         IF(nlc>0) CALL mpi_allreduce(aux1,aux2,nlc,MPI_REAL8,MPI_SUM, &
         mpi_comm_rows,mpierr)

         ! Transform

         nlc = 0
         DO j=1,lc-1
            lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
            IF(lcx>0) THEN
               nlc = nlc+1
               a(1:lr,lcx) = a(1:lr,lcx) - tau*aux2(nlc)*vr(1:lr)
            ENDIF
         ENDDO

      ENDDO

      ! Calculate scalar products of stored Householder vectors.
      ! This can be done in different ways, we use dsyrk

      vav = 0
      IF(l_rows>0) &
         CALL dsyrk('U','T',n_cols,l_rows,1.d0,vmr,UBOUND(vmr,1),0.d0, &
         vav,UBOUND(vav,1))
      CALL symm_matrix_allreduce(n_cols,vav,UBOUND(vav,1),mpi_comm_rows)

      ! Calculate triangular matrix T for block Householder Transformation

      DO lc=n_cols,1,-1
         tau = tmat(lc,lc,istep)
         IF(lc<n_cols) THEN
            CALL dtrmv('U','T','N',n_cols-lc,tmat(lc+1,lc+1,istep), &
            UBOUND(tmat,1),vav(lc+1,lc),1)
            tmat(lc,lc+1:n_cols,istep) = -tau * vav(lc+1:n_cols,lc)
         ENDIF
      ENDDO
      
      ENDIF      

      ! Transpose vmr -> vmc (stored in umc, second half)

      CALL elpa_transpose_vectors  (vmr, UBOUND(vmr,1), mpi_comm_rows, &
      umc(1,n_cols+1), UBOUND(umc,1), mpi_comm_cols, & 
      1, istep*nbw, n_cols, nblk)

      ! Calculate umc = A**T * vmr
      ! Note that the distributed A has to be transposed
      ! Opposed to direct tridiagonalization there is no need to use the cache locality
      ! of the tiles, so we can use strips of the matrix

      umc(1:l_cols,1:n_cols) = 0.d0
      vmr(1:l_rows,n_cols+1:2*n_cols) = 0
      IF(l_cols>0 .AND. l_rows>0) THEN
         DO i=0,(istep*nbw-1)/tile_size

            lcs = i*l_cols_tile+1
            lce = MIN(l_cols,(i+1)*l_cols_tile)
            IF(lce<lcs) CYCLE

            lre = MIN(l_rows,(i+1)*l_rows_tile)
            CALL DGEMM('T','N',lce-lcs+1,n_cols,lre,1.d0,a(1,lcs),UBOUND(a,1), &
                       vmr,UBOUND(vmr,1),1.d0,umc(lcs,1),UBOUND(umc,1))

            IF(i==0) CYCLE
            lre = MIN(l_rows,i*l_rows_tile)
            CALL DGEMM('N','N',lre,n_cols,lce-lcs+1,1.d0,a(1,lcs),lda, &
            umc(lcs,n_cols+1),UBOUND(umc,1),1.d0,vmr(1,n_cols+1),UBOUND(vmr,1))
         ENDDO
      ENDIF

      ! Sum up all ur(:) parts along rows and add them to the uc(:) parts
      ! on the processors containing the diagonal
      ! This is only necessary if ur has been calculated, i.e. if the
      ! global tile size is smaller than the global remaining matrix

      IF(tile_size < istep*nbw) THEN
         CALL elpa_reduce_add_vectors  (vmr(1,n_cols+1),UBOUND(vmr,1), &
         mpi_comm_rows, umc, UBOUND(umc,1), mpi_comm_cols, &
         istep*nbw, n_cols, nblk)
      ENDIF

      IF(l_cols>0) THEN
         ALLOCATE(tmp(l_cols,n_cols))
         CALL mpi_allreduce(umc,tmp,l_cols*n_cols,MPI_REAL8,MPI_SUM, &
         mpi_comm_rows,mpierr)
         umc(1:l_cols,1:n_cols) = tmp(1:l_cols,1:n_cols)
         DEALLOCATE(tmp)
      ENDIF

      ! U = U * Tmat**T

      CALL dtrmm('Right','Upper','Trans','Nonunit',l_cols,n_cols,1.d0, &
      tmat(1,1,istep),UBOUND(tmat,1),umc,UBOUND(umc,1))

      ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

      CALL dgemm('T','N',n_cols,n_cols,l_cols,1.d0,umc,UBOUND(umc,1), & 
      umc(1,n_cols+1),UBOUND(umc,1),0.d0,vav,UBOUND(vav,1))
      CALL dtrmm('Right','Upper','Trans','Nonunit',n_cols,n_cols,1.d0, &
      tmat(1,1,istep),UBOUND(tmat,1),vav,UBOUND(vav,1))

      CALL symm_matrix_allreduce(n_cols,vav,UBOUND(vav,1),mpi_comm_cols)

      ! U = U - 0.5 * V * VAV
      CALL dgemm('N','N',l_cols,n_cols,n_cols,-0.5d0,umc(1,n_cols+1), &
      UBOUND(umc,1),vav,UBOUND(vav,1),1.d0,umc,UBOUND(umc,1))

      ! Transpose umc -> umr (stored in vmr, second half)

       CALL elpa_transpose_vectors  (umc, UBOUND(umc,1), mpi_comm_cols, &
      vmr(1,n_cols+1), UBOUND(vmr,1), mpi_comm_rows, &
      1, istep*nbw, n_cols, nblk)

      ! A = A - V*U**T - U*V**T

      DO i=0,(istep*nbw-1)/tile_size
         lcs = i*l_cols_tile+1
         lce = MIN(l_cols,(i+1)*l_cols_tile)
         lre = MIN(l_rows,(i+1)*l_rows_tile)
         IF(lce<lcs .OR. lre<1) CYCLE
         CALL dgemm('N','T',lre,lce-lcs+1,2*n_cols,-1.d0, &
                    vmr,UBOUND(vmr,1),umc(lcs,1),UBOUND(umc,1), &
                    1.d0,a(1,lcs),lda)
      ENDDO

      DEALLOCATE(vmr, umc, vr)

   ENDDO
 
   IF (which_qr_decomposition == 1) THEN
      DEALLOCATE(work_blocked)
      DEALLOCATE(tauvector)
   ENDIF

END SUBROUTINE bandred_real

!-------------------------------------------------------------------------------

SUBROUTINE symm_matrix_allreduce(n,a,lda,comm)

!-------------------------------------------------------------------------------
!  symm_matrix_allreduce: Does an mpi_allreduce for a symmetric matrix A.
!  On entry, only the upper half of A needs to be set
!  On exit, the complete matrix is set
!-------------------------------------------------------------------------------

    INTEGER                                  :: n, lda
    REAL*8                                   :: a(lda,*)
    INTEGER                                  :: comm

    INTEGER                                  :: i, mpierr, nc
    REAL*8                                   :: h1(n*n), h2(n*n)

   nc = 0
   DO i=1,n
      h1(nc+1:nc+i) = a(1:i,i)
      nc = nc+i
   ENDDO

   CALL mpi_allreduce(h1,h2,nc,MPI_REAL8,MPI_SUM,comm,mpierr)

   nc = 0
   DO i=1,n
      a(1:i,i) = h2(nc+1:nc+i)
      a(i,1:i-1) = a(1:i-1,i)
      nc = nc+i
   ENDDO

END SUBROUTINE symm_matrix_allreduce

!-------------------------------------------------------------------------------

SUBROUTINE trans_ev_band_to_full_real(na, nqc, nblk, nbw, a, lda, tmat, q, &
ldq, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  trans_ev_band_to_full_real:
!  Transforms the eigenvectors of a band matrix back to the eigenvectors of the original matrix
!
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nqc         Number of columns of matrix q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith
!
!  a(lda,*)    Matrix containing the Householder vectors (i.e. matrix a after bandred_real)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a
!
!  tmat(nbw,nbw,.) Factors returned by bandred_real
!
!  q           On input: Eigenvectors of band matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------


    INTEGER                                  :: na, nqc, nblk, nbw, lda
    REAL*8                                   :: a(lda,*), tmat(nbw, nbw, *)
    INTEGER                                  :: ldq
    REAL*8                                   :: q(ldq,*)
    INTEGER                                  :: mpi_comm_rows, mpi_comm_cols

    INTEGER :: i, istep, l_colh, l_cols, l_rows, lc, max_blocks_col, &
      max_blocks_row, max_local_cols, max_local_rows, mpierr, my_pcol, &
      my_prow, n_cols, nb, ncol, np_cols, np_rows, nrow, ns, pcol, prow
    REAL*8, ALLOCATABLE                      :: hvb(:), hvm(:,:), tmp1(:), &
                                                tmp2(:)

   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)


   max_blocks_row = ((na -1)/nblk)/np_rows + 1  ! Rows of A
   max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q!

   max_local_rows = max_blocks_row*nblk
   max_local_cols = max_blocks_col*nblk

   ALLOCATE(tmp1(max_local_cols*nbw))
   ALLOCATE(tmp2(max_local_cols*nbw))
   ALLOCATE(hvb(max_local_rows*nbw))
   ALLOCATE(hvm(max_local_rows,nbw))

   hvm = 0   ! Must be set to 0 !!!
   hvb = 0   ! Safety only

   l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q

   DO istep=1,(na-1)/nbw

      n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

      ! Broadcast all Householder vectors for current step compressed in hvb

      nb = 0
      ns = 0

      DO lc = 1, n_cols
         ncol = istep*nbw + lc ! absolute column number of householder vector
         nrow = ncol - nbw ! absolute number of pivot row

         l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast
         l_colh = local_index(ncol  , my_pcol, np_cols, nblk, -1) ! HV local column number

         IF(my_pcol==pcol(ncol)) hvb(nb+1:nb+l_rows) = a(1:l_rows,l_colh)

         nb = nb+l_rows

         IF(lc==n_cols .OR. MOD(ncol,nblk)==0) THEN
            CALL MPI_Bcast(hvb(ns+1),nb-ns,MPI_REAL8,pcol(ncol),mpi_comm_cols,mpierr)
            ns = nb
         ENDIF
      ENDDO

      ! Expand compressed Householder vectors into matrix hvm

      nb = 0
      DO lc = 1, n_cols
         nrow = (istep-1)*nbw+lc ! absolute number of pivot row
         l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast

         hvm(1:l_rows,lc) = hvb(nb+1:nb+l_rows)
         IF(my_prow==prow(nrow)) hvm(l_rows+1,lc) = 1.

         nb = nb+l_rows
      ENDDO

      l_rows = local_index(MIN(na,(istep+1)*nbw), my_prow, np_rows, nblk, -1)

      ! Q = Q - V * T**T * V**T * Q

      IF(l_rows>0) THEN
         CALL dgemm('T','N',n_cols,l_cols,l_rows,1.d0,hvm,UBOUND(hvm,1), &
                    q,ldq,0.d0,tmp1,n_cols)
      ELSE
         tmp1(1:l_cols*n_cols) = 0
      ENDIF
      CALL mpi_allreduce(tmp1,tmp2,n_cols*l_cols,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
      IF(l_rows>0) THEN
         CALL dtrmm('L','U','T','N',n_cols,l_cols,1.0d0,tmat(1,1,istep),UBOUND(tmat,1),tmp2,n_cols)
         CALL dgemm('N','N',l_rows,l_cols,n_cols,-1.d0,hvm,UBOUND(hvm,1), &
                    tmp2,n_cols,1.d0,q,ldq)
      ENDIF

   ENDDO

   DEALLOCATE(tmp1, tmp2, hvb, hvm)


END SUBROUTINE trans_ev_band_to_full_real

! --------------------------------------------------------------------------------------------------

SUBROUTINE tridiag_band_real(na, nb, nblk, a, lda, d, e, mpi_comm_rows, mpi_comm_cols, mpi_comm)

!-------------------------------------------------------------------------------
! tridiag_band_real:
! Reduces a real symmetric band matrix to tridiagonal form
!
!  na          Order of matrix a
!
!  nb          Semi bandwith
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  a(lda,*)    Distributed system matrix reduced to banded form in the upper diagonal
!
!  lda         Leading dimension of a
!
!  d(na)       Diagonal of tridiagonal matrix, set only on PE 0 (output)
!
!  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0 (output)
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm
!              MPI-Communicator for the total processor set
!-------------------------------------------------------------------------------


    INTEGER, INTENT(in)                      :: na, nb, nblk, lda
    REAL*8, INTENT(in)                       :: a(lda,*)
    REAL*8, INTENT(out)                      :: d(na), e(na)
    INTEGER, INTENT(in)                      :: mpi_comm_rows, mpi_comm_cols, &
                                                mpi_comm

    COMPLEX*16                               :: c_a(1,1), c_ab(1,1)
    INTEGER :: i, iblk, ireq_ab, ireq_hv, istep, iter, local_size, &
      max_blk_size, max_threads, mpi_status(MPI_STATUS_SIZE), mpierr, &
      my_block_e, my_block_s, my_pcol, my_pe, my_prow, my_thread, n, n_off, &
      n_pes, na_s, nblocks, nblocks_total, nc, ne, np_cols, np_rows, nr, ns, &
      nt, num_chunks, num_hh_vecs, nx
    INTEGER, ALLOCATABLE :: block_limits(:), global_id(:,:), &
      global_id_tmp(:,:), hh_cnt(:), hh_dst(:), ireq_hhr(:), ireq_hhs(:), &
      limits(:), mpi_statuses(:,:), omp_block_limits(:), snd_limits(:,:)
    REAL*8                                   :: ab_s(1+nb), h(nb), hd(nb), &
                                                hf, hs(nb), hv(nb), &
                                                hv_new(nb), hv_s(nb), tau, &
                                                tau_new, vnorm2, x
    REAL*8, ALLOCATABLE                      :: ab(:,:), hh_gath(:,:,:), &
                                                hh_send(:,:,:), hv_t(:,:), &
                                                tau_t(:)

! set only on PE 0
! dummies for calling redist_band

!$  integer :: omp_get_max_threads


    CALL mpi_comm_rank(mpi_comm,my_pe,mpierr)
    CALL mpi_comm_size(mpi_comm,n_pes,mpierr)

    CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
    CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
    CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
    CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

    ! Get global_id mapping 2D procssor coordinates to global id

    ALLOCATE(global_id(0:np_rows-1,0:np_cols-1))
    ALLOCATE(global_id_tmp(0:np_rows-1,0:np_cols-1))
    global_id(:,:) = 0
    global_id(my_prow, my_pcol) = my_pe

    global_id_tmp(:,:) = global_id(:,:)
    CALL mpi_allreduce(global_id_tmp, global_id, np_rows*np_cols, mpi_integer, mpi_sum, mpi_comm, mpierr)
    DEALLOCATE(global_id_tmp)


    ! Total number of blocks in the band:

    nblocks_total = (na-1)/nb + 1

    ! Set work distribution

    ALLOCATE(block_limits(0:n_pes))
    CALL divide_band(nblocks_total, n_pes, block_limits)

    ! nblocks: the number of blocks for my task
    nblocks = block_limits(my_pe+1) - block_limits(my_pe)

    ! allocate the part of the band matrix which is needed by this PE
    ! The size is 1 block larger than needed to avoid extensive shifts
    ALLOCATE(ab(2*nb,(nblocks+1)*nb))
    ab = 0 ! needed for lower half, the extra block should also be set to 0 for safety

    ! n_off: Offset of ab within band
    n_off = block_limits(my_pe)*nb

    ! Redistribute band in a to ab
    CALL redist_band(.TRUE., a, c_a, lda, na, nblk, nb, mpi_comm_rows, mpi_comm_cols, mpi_comm, ab, c_ab)

    ! Calculate the workload for each sweep in the back transformation
    ! and the space requirements to hold the HH vectors

    ALLOCATE(limits(0:np_rows))
    CALL determine_workload(na, nb, np_rows, limits)
    max_blk_size = MAXVAL(limits(1:np_rows) - limits(0:np_rows-1))

    num_hh_vecs = 0
    num_chunks  = 0
    nx = na
    DO n = 1, nblocks_total
      CALL determine_workload(nx, nb, np_rows, limits)
      local_size = limits(my_prow+1) - limits(my_prow)
      ! add to number of householder vectors
      ! please note: for nx==1 the one and only HH vector is 0 and is neither calculated nor send below!
      IF(MOD(n-1,np_cols) == my_pcol .AND. local_size>0 .AND. nx>1) THEN
        num_hh_vecs = num_hh_vecs + local_size
        num_chunks  = num_chunks+1
      ENDIF
      nx = nx - nb
    ENDDO

    ! Allocate space for HH vectors

    ALLOCATE(hh_trans_real(nb,num_hh_vecs))

    ! Allocate and init MPI requests

    ALLOCATE(ireq_hhr(num_chunks)) ! Recv requests
    ALLOCATE(ireq_hhs(nblocks))    ! Send requests

    num_hh_vecs = 0
    num_chunks  = 0
    nx = na
    nt = 0
    DO n = 1, nblocks_total
      CALL determine_workload(nx, nb, np_rows, limits)
      local_size = limits(my_prow+1) - limits(my_prow)
      IF(MOD(n-1,np_cols) == my_pcol .AND. local_size>0 .AND. nx>1) THEN
        num_chunks  = num_chunks+1
        CALL mpi_irecv(hh_trans_real(1,num_hh_vecs+1), nb*local_size, mpi_real8, nt, &
                       10+n-block_limits(nt), mpi_comm, ireq_hhr(num_chunks), mpierr)
        num_hh_vecs = num_hh_vecs + local_size
      ENDIF
      nx = nx - nb
      IF(n == block_limits(nt+1)) THEN
        nt = nt + 1
      ENDIF
    ENDDO

    ireq_hhs(:) = MPI_REQUEST_NULL

    ! Buffers for gathering/sending the HH vectors

    ALLOCATE(hh_gath(nb,max_blk_size,nblocks)) ! gathers HH vectors
    ALLOCATE(hh_send(nb,max_blk_size,nblocks)) ! send buffer for HH vectors
    hh_gath(:,:,:) = 0
    hh_send(:,:,:) = 0

    ! Some counters

    ALLOCATE(hh_cnt(nblocks))
    ALLOCATE(hh_dst(nblocks))

    hh_cnt(:) = 1 ! The first transfomation vector is always 0 and not calculated at all
    hh_dst(:) = 0 ! PE number for receive

    ireq_ab = MPI_REQUEST_NULL
    ireq_hv = MPI_REQUEST_NULL

    ! Limits for sending

    ALLOCATE(snd_limits(0:np_rows,nblocks))

    DO iblk=1,nblocks
      CALL determine_workload(na-(iblk+block_limits(my_pe)-1)*nb, nb, np_rows, snd_limits(:,iblk))
    ENDDO

    ! OpenMP work distribution:

    max_threads = 1
!$ max_threads = omp_get_max_threads()

    ! For OpenMP we need at least 2 blocks for every thread
    max_threads = MIN(max_threads, nblocks/2)
    IF(max_threads==0) max_threads = 1

    ALLOCATE(omp_block_limits(0:max_threads))

    ! Get the OpenMP block limits
    CALL divide_band(nblocks, max_threads, omp_block_limits)

    ALLOCATE(hv_t(nb,max_threads), tau_t(max_threads))
    hv_t = 0
    tau_t = 0

    ! ---------------------------------------------------------------------------
    ! Start of calculations

    na_s = block_limits(my_pe)*nb + 1

    IF(my_pe>0 .AND. na_s<=na) THEN
      ! send first column to previous PE
      ! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
      ab_s(1:nb+1) = ab(1:nb+1,na_s-n_off)
      CALL mpi_isend(ab_s,nb+1,mpi_real8,my_pe-1,1,mpi_comm,ireq_ab,mpierr)
    ENDIF

    DO istep=1,na-1-block_limits(my_pe)*nb

      IF(my_pe==0) THEN
        n = MIN(na-na_s,nb) ! number of rows to be reduced
        hv(:) = 0
        tau = 0
        ! The last step (istep=na-1) is only needed for sending the last HH vectors.
        ! We don't want the sign of the last element flipped (analogous to the other sweeps)
        IF(istep < na-1) THEN
          ! Transform first column of remaining matrix
          vnorm2 = SUM(ab(3:n+1,na_s-n_off)**2)
          CALL hh_transform_real(ab(2,na_s-n_off),vnorm2,hf,tau)
          hv(1) = 1
          hv(2:n) = ab(3:n+1,na_s-n_off)*hf
        ENDIF
        d(istep) = ab(1,na_s-n_off)
        e(istep) = ab(2,na_s-n_off)
        IF(istep == na-1) THEN
          d(na) = ab(1,na_s+1-n_off)
          e(na) = 0
        ENDIF
      ELSE
        IF(na>na_s) THEN
          ! Receive Householder vector from previous task, from PE owning subdiagonal
          CALL mpi_recv(hv,nb,mpi_real8,my_pe-1,2,mpi_comm,mpi_status,mpierr)
          tau = hv(1)
          hv(1) = 1.
        ENDIF
      ENDIF

      na_s = na_s+1
      IF(na_s-n_off > nb) THEN
        ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
        ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0
        n_off = n_off + nb
      ENDIF

      IF(max_threads > 1) THEN

        ! Codepath for OpenMP

        ! Please note that in this case it is absolutely necessary to have at least 2 blocks per thread!
        ! Every thread is one reduction cycle behind its predecessor and thus starts one step later.
        ! This simulates the behaviour of the MPI tasks which also work after each other.
        ! The code would be considerably easier, if the MPI communication would be made within
        ! the parallel region - this is avoided here since this would require 
        ! MPI_Init_thread(MPI_THREAD_MULTIPLE) at the start of the program.

        hv_t(:,1) = hv
        tau_t(1) = tau

        DO iter = 1, 2

          ! iter=1 : work on first block
          ! iter=2 : work on remaining blocks
          ! This is done in 2 iterations so that we have a barrier in between:
          ! After the first iteration, it is guaranteed that the last row of the last block
          ! is completed by the next thread.
          ! After the first iteration it is also the place to exchange the last row
          ! with MPI calls

!$omp parallel do private(my_thread, my_block_s, my_block_e, iblk, ns, ne, hv, tau, &
!$omp nc, nr, hs, hd, vnorm2, hf, x, h, i), schedule(static,1), num_threads(max_threads)
          DO my_thread = 1, max_threads

            IF(iter == 1) THEN
              my_block_s = omp_block_limits(my_thread-1) + 1
              my_block_e = my_block_s
            ELSE
              my_block_s = omp_block_limits(my_thread-1) + 2
              my_block_e = omp_block_limits(my_thread)
            ENDIF

            DO iblk = my_block_s, my_block_e

              ns = na_s + (iblk-1)*nb - n_off - my_thread + 1 ! first column in block
              ne = ns+nb-1                    ! last column in block

              IF(istep<my_thread .OR. ns+n_off>na) EXIT

              hv = hv_t(:,my_thread)
              tau = tau_t(my_thread)

              ! Store Householder vector for back transformation

              hh_cnt(iblk) = hh_cnt(iblk) + 1

              hh_gath(1   ,hh_cnt(iblk),iblk) = tau
              hh_gath(2:nb,hh_cnt(iblk),iblk) = hv(2:nb)

              nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
              nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                            ! Note that nr>=0 implies that diagonal block is full (nc==nb)!

              ! Transform diagonal block

              CALL DSYMV('L',nc,tau,ab(1,ns),2*nb-1,hv,1,0.d0,hd,1)

              x = DOT_PRODUCT(hv(1:nc),hd(1:nc))*tau
              hd(1:nc) = hd(1:nc) - 0.5*x*hv(1:nc)

              CALL DSYR2('L',nc,-1.d0,hd,1,hv,1,ab(1,ns),2*nb-1)

              hv_t(:,my_thread) = 0
              tau_t(my_thread)  = 0

              IF(nr<=0) CYCLE ! No subdiagonal block present any more

              ! Transform subdiagonal block

              CALL DGEMV('N',nr,nb,tau,ab(nb+1,ns),2*nb-1,hv,1,0.d0,hs,1)

              IF(nr>1) THEN

                ! complete (old) Householder transformation for first column

                ab(nb+1:nb+nr,ns) = ab(nb+1:nb+nr,ns) - hs(1:nr) ! Note: hv(1) == 1

                ! calculate new Householder transformation for first column
                ! (stored in hv_t(:,my_thread) and tau_t(my_thread))

                vnorm2 = SUM(ab(nb+2:nb+nr,ns)**2)
                CALL hh_transform_real(ab(nb+1,ns),vnorm2,hf,tau_t(my_thread))
                hv_t(1   ,my_thread) = 1.
                hv_t(2:nr,my_thread) = ab(nb+2:nb+nr,ns)*hf
                ab(nb+2:,ns) = 0

                ! update subdiagonal block for old and new Householder transformation
                ! This way we can use a nonsymmetric rank 2 update which is (hopefully) faster

                CALL DGEMV('T',nr,nb-1,tau_t(my_thread),ab(nb,ns+1),2*nb-1,hv_t(1,my_thread),1,0.d0,h(2),1)
                x = DOT_PRODUCT(hs(1:nr),hv_t(1:nr,my_thread))*tau_t(my_thread)
                h(2:nb) = h(2:nb) - x*hv(2:nb)
                ! Unfortunately there is no BLAS routine like DSYR2 for a nonsymmetric rank 2 update ("DGER2")
                DO i=2,nb
                  ab(2+nb-i:1+nb+nr-i,i+ns-1) = ab(2+nb-i:1+nb+nr-i,i+ns-1) - hv_t(1:nr,my_thread)*h(i) - hs(1:nr)*hv(i)
                ENDDO

              ELSE

                ! No new Householder transformation for nr=1, just complete the old one
                ab(nb+1,ns) = ab(nb+1,ns) - hs(1) ! Note: hv(1) == 1
                DO i=2,nb
                  ab(2+nb-i,i+ns-1) = ab(2+nb-i,i+ns-1) - hs(1)*hv(i)
                ENDDO
                ! For safety: there is one remaining dummy transformation (but tau is 0 anyways)
                hv_t(1,my_thread) = 1.

              ENDIF

            ENDDO

          ENDDO ! my_thread
!$omp end parallel do

          IF (iter==1) THEN
            ! We are at the end of the first block

            ! Send our first column to previous PE
            IF(my_pe>0 .AND. na_s <= na) THEN
              CALL mpi_wait(ireq_ab,mpi_status,mpierr)
              ab_s(1:nb+1) = ab(1:nb+1,na_s-n_off)
              CALL mpi_isend(ab_s,nb+1,mpi_real8,my_pe-1,1,mpi_comm,ireq_ab,mpierr)
            ENDIF

            ! Request last column from next PE
            ne = na_s + nblocks*nb - (max_threads-1) - 1
            IF(istep>=max_threads .AND. ne <= na) THEN
              CALL mpi_recv(ab(1,ne-n_off),nb+1,mpi_real8,my_pe+1,1,mpi_comm,mpi_status,mpierr)
            ENDIF

          ELSE
            ! We are at the end of all blocks

            ! Send last HH vector and TAU to next PE if it has been calculated above
            ne = na_s + nblocks*nb - (max_threads-1) - 1
            IF(istep>=max_threads .AND. ne < na) THEN
              CALL mpi_wait(ireq_hv,mpi_status,mpierr)
              hv_s(1) = tau_t(max_threads)
              hv_s(2:) = hv_t(2:,max_threads)
              CALL mpi_isend(hv_s,nb,mpi_real8,my_pe+1,2,mpi_comm,ireq_hv,mpierr)
            ENDIF

            ! "Send" HH vector and TAU to next OpenMP thread
            DO my_thread = max_threads, 2, -1
              hv_t(:,my_thread) = hv_t(:,my_thread-1)
              tau_t(my_thread)  = tau_t(my_thread-1)
            ENDDO

          ENDIF
        ENDDO ! iter

      ELSE

        ! Codepath for 1 thread without OpenMP

        ! The following code is structured in a way to keep waiting times for
        ! other PEs at a minimum, especially if there is only one block.
        ! For this reason, it requests the last column as late as possible
        ! and sends the Householder vector and the first column as early
        ! as possible.

        DO iblk=1,nblocks

          ns = na_s + (iblk-1)*nb - n_off ! first column in block
          ne = ns+nb-1                    ! last column in block

          IF(ns+n_off>na) EXIT

          ! Store Householder vector for back transformation

          hh_cnt(iblk) = hh_cnt(iblk) + 1

          hh_gath(1   ,hh_cnt(iblk),iblk) = tau
          hh_gath(2:nb,hh_cnt(iblk),iblk) = hv(2:nb)

          nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
          nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                        ! Note that nr>=0 implies that diagonal block is full (nc==nb)!

          ! Multiply diagonal block and subdiagonal block with Householder vector

          IF(iblk==nblocks .AND. nc==nb) THEN

            ! We need the last column from the next PE.
            ! First do the matrix multiplications without last column ...

            ! Diagonal block, the contribution of the last element is added below!
            ab(1,ne) = 0
            CALL DSYMV('L',nc,tau,ab(1,ns),2*nb-1,hv,1,0.d0,hd,1)

            ! Subdiagonal block
            IF(nr>0) CALL DGEMV('N',nr,nb-1,tau,ab(nb+1,ns),2*nb-1,hv,1,0.d0,hs,1)

            ! ... then request last column ...
            CALL mpi_recv(ab(1,ne),nb+1,mpi_real8,my_pe+1,1,mpi_comm,mpi_status,mpierr)

            ! ... and complete the result
            hs(1:nr) = hs(1:nr) + ab(2:nr+1,ne)*tau*hv(nb)
            hd(nb) = hd(nb) + ab(1,ne)*hv(nb)*tau

          ELSE

            ! Normal matrix multiply
            CALL DSYMV('L',nc,tau,ab(1,ns),2*nb-1,hv,1,0.d0,hd,1)
            IF(nr>0) CALL DGEMV('N',nr,nb,tau,ab(nb+1,ns),2*nb-1,hv,1,0.d0,hs,1)

          ENDIF

          ! Calculate first column of subdiagonal block and calculate new
          ! Householder transformation for this column

          hv_new(:) = 0 ! Needed, last rows must be 0 for nr < nb
          tau_new = 0

          IF(nr>0) THEN

            ! complete (old) Householder transformation for first column

            ab(nb+1:nb+nr,ns) = ab(nb+1:nb+nr,ns) - hs(1:nr) ! Note: hv(1) == 1

            ! calculate new Householder transformation ...
            IF(nr>1) THEN
              vnorm2 = SUM(ab(nb+2:nb+nr,ns)**2)
              CALL hh_transform_real(ab(nb+1,ns),vnorm2,hf,tau_new)
              hv_new(1) = 1.
              hv_new(2:nr) = ab(nb+2:nb+nr,ns)*hf
              ab(nb+2:,ns) = 0
            ENDIF

            ! ... and send it away immediatly if this is the last block

            IF(iblk==nblocks) THEN
              CALL mpi_wait(ireq_hv,mpi_status,mpierr)
              hv_s(1) = tau_new
              hv_s(2:) = hv_new(2:)
              CALL mpi_isend(hv_s,nb,mpi_real8,my_pe+1,2,mpi_comm,ireq_hv,mpierr)
            ENDIF

          ENDIF


          ! Transform diagonal block
          x = DOT_PRODUCT(hv(1:nc),hd(1:nc))*tau
          hd(1:nc) = hd(1:nc) - 0.5*x*hv(1:nc)

          IF(my_pe>0 .AND. iblk==1) THEN

            ! The first column of the diagonal block has to be send to the previous PE
            ! Calculate first column only ...

            ab(1:nc,ns) = ab(1:nc,ns) - hd(1:nc)*hv(1) - hv(1:nc)*hd(1)

            ! ... send it away ...

            CALL mpi_wait(ireq_ab,mpi_status,mpierr)
            ab_s(1:nb+1) = ab(1:nb+1,ns)
            CALL mpi_isend(ab_s,nb+1,mpi_real8,my_pe-1,1,mpi_comm,ireq_ab,mpierr)

            ! ... and calculate remaining columns with rank-2 update
            IF(nc>1) CALL DSYR2('L',nc-1,-1.d0,hd(2),1,hv(2),1,ab(1,ns+1),2*nb-1)
          ELSE
            ! No need to  send, just a rank-2 update
            CALL DSYR2('L',nc,-1.d0,hd,1,hv,1,ab(1,ns),2*nb-1)
          ENDIF

          ! Do the remaining double Householder transformation on the subdiagonal block cols 2 ... nb

          IF(nr>0) THEN
            IF(nr>1) THEN
              CALL DGEMV('T',nr,nb-1,tau_new,ab(nb,ns+1),2*nb-1,hv_new,1,0.d0,h(2),1)
              x = DOT_PRODUCT(hs(1:nr),hv_new(1:nr))*tau_new
              h(2:nb) = h(2:nb) - x*hv(2:nb)
              ! Unfortunately there is no BLAS routine like DSYR2 for a nonsymmetric rank 2 update ("DGER2")
              DO i=2,nb
                ab(2+nb-i:1+nb+nr-i,i+ns-1) = ab(2+nb-i:1+nb+nr-i,i+ns-1) - hv_new(1:nr)*h(i) - hs(1:nr)*hv(i)
              ENDDO
            ELSE
              ! No double Householder transformation for nr=1, just complete the row
              DO i=2,nb
                ab(2+nb-i,i+ns-1) = ab(2+nb-i,i+ns-1) - hs(1)*hv(i)
              ENDDO
            ENDIF
          ENDIF

          ! Use new HH vector for the next block
          hv(:) = hv_new(:)
          tau = tau_new

        ENDDO

      ENDIF

      DO iblk = 1, nblocks

        IF(hh_dst(iblk) >= np_rows) EXIT
        IF(snd_limits(hh_dst(iblk)+1,iblk) == snd_limits(hh_dst(iblk),iblk)) EXIT

        IF(hh_cnt(iblk) == snd_limits(hh_dst(iblk)+1,iblk)-snd_limits(hh_dst(iblk),iblk)) THEN
          ! Wait for last transfer to finish
          CALL mpi_wait(ireq_hhs(iblk), mpi_status, mpierr)
          ! Copy vectors into send buffer
          hh_send(:,1:hh_cnt(iblk),iblk) = hh_gath(:,1:hh_cnt(iblk),iblk)
          ! Send to destination
          CALL mpi_isend(hh_send(1,1,iblk), nb*hh_cnt(iblk), mpi_real8, &
                         global_id(hh_dst(iblk),MOD(iblk+block_limits(my_pe)-1,np_cols)), &
                         10+iblk, mpi_comm, ireq_hhs(iblk), mpierr)
          ! Reset counter and increase destination row
          hh_cnt(iblk) = 0
          hh_dst(iblk) = hh_dst(iblk)+1
        ENDIF

      ENDDO

    ENDDO

    ! Finish the last outstanding requests
    CALL mpi_wait(ireq_ab,mpi_status,mpierr)
    CALL mpi_wait(ireq_hv,mpi_status,mpierr)

    ALLOCATE(mpi_statuses(MPI_STATUS_SIZE,MAX(nblocks,num_chunks)))
    CALL mpi_waitall(nblocks, ireq_hhs, mpi_statuses, mpierr)
    CALL mpi_waitall(num_chunks, ireq_hhr, mpi_statuses, mpierr)
    DEALLOCATE(mpi_statuses)

    CALL mpi_barrier(mpi_comm,mpierr)

    DEALLOCATE(ab)
    DEALLOCATE(ireq_hhr, ireq_hhs)
    DEALLOCATE(hh_cnt, hh_dst)
    DEALLOCATE(hh_gath, hh_send)
    DEALLOCATE(limits, snd_limits)
    DEALLOCATE(block_limits)
    DEALLOCATE(global_id)

END SUBROUTINE tridiag_band_real

! --------------------------------------------------------------------------------------------------

SUBROUTINE trans_ev_tridi_to_band_real(na, nev, nblk, nbw, q, ldq, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  trans_ev_tridi_to_band_real:
!  Transforms the eigenvectors of a tridiagonal matrix back to the eigenvectors of the band matrix
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nev         Number eigenvectors to compute (= columns of matrix q)
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nb          semi bandwith
!
!  q           On input: Eigenvectors of tridiagonal matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns/both
!
!-------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(in) :: na, nev, nblk, nbw, ldq, mpi_comm_rows, mpi_comm_cols
    REAL*8 q(ldq,*)

    INTEGER np_rows, my_prow, np_cols, my_pcol

    INTEGER i, j, ip, sweep, nbuf, l_nev, a_dim2
    INTEGER current_n, current_local_n, current_n_start, current_n_end
    INTEGER next_n, next_local_n, next_n_start, next_n_end
    INTEGER bottom_msg_length, top_msg_length, next_top_msg_length
    INTEGER thread_width, stripe_width, stripe_count, csw
    INTEGER num_result_blocks, num_result_buffers, num_bufs_recvd
    INTEGER a_off, current_tv_off, max_blk_size, b_off, b_len
    INTEGER mpierr, src, src_offset, dst, offset, nfact, num_blk
    INTEGER mpi_status(MPI_STATUS_SIZE)
    LOGICAL flag

    REAL*8, ALLOCATABLE :: a(:,:,:,:), row(:)
    REAL*8, ALLOCATABLE :: top_border_send_buffer(:,:), top_border_recv_buffer(:,:)
    REAL*8, ALLOCATABLE :: bottom_border_send_buffer(:,:), bottom_border_recv_buffer(:,:)
    REAL*8, ALLOCATABLE :: result_buffer(:,:,:)
    REAL*8, ALLOCATABLE :: bcast_buffer(:,:)

    INTEGER n_off
    INTEGER, ALLOCATABLE :: result_send_request(:), result_recv_request(:), limits(:)
    INTEGER, ALLOCATABLE :: top_send_request(:), bottom_send_request(:)
    INTEGER, ALLOCATABLE :: top_recv_request(:), bottom_recv_request(:)
    INTEGER, ALLOCATABLE :: mpi_statuses(:,:)

    ! MPI send/recv tags, arbitrary

    INTEGER, PARAMETER :: bottom_recv_tag = 111
    INTEGER, PARAMETER :: top_recv_tag    = 222
    INTEGER, PARAMETER :: result_recv_tag = 333

    INTEGER :: max_threads, my_thread
!$  integer :: omp_get_max_threads

    ! Just for measuring the kernel performance
    REAL*8 kernel_time
    INTEGER*8 kernel_flops


    kernel_time = 1.d-100
    kernel_flops = 0

    max_threads = 1
!$  max_threads = omp_get_max_threads()

    CALL MPI_Comm_rank(mpi_comm_rows, my_prow, mpierr)
    CALL MPI_Comm_size(mpi_comm_rows, np_rows, mpierr)
    CALL MPI_Comm_rank(mpi_comm_cols, my_pcol, mpierr)
    CALL MPI_Comm_size(mpi_comm_cols, np_cols, mpierr)

    IF(MOD(nbw,nblk)/=0) THEN
      IF(my_prow==0 .AND. my_pcol==0) THEN
         PRINT *,'ERROR: nbw=',nbw,', nblk=',nblk
         PRINT *,'band backtransform works only for nbw==n*nblk'
         CALL mpi_abort(mpi_comm_world,0,mpierr)
      ENDIF
    ENDIF

    nfact = nbw / nblk


    ! local number of eigenvectors
    l_nev = local_index(nev, my_pcol, np_cols, nblk, -1)

    IF(l_nev==0) THEN
        thread_width = 0
        stripe_width = 0
        stripe_count = 0
    ELSE
        ! Suggested stripe width is 48 since 48*64 real*8 numbers should fit into
        ! every primary cache
        thread_width = (l_nev-1)/max_threads + 1 ! number of eigenvectors per OMP thread
        stripe_width = 48 ! Must be a multiple of 4
        stripe_count = (thread_width-1)/stripe_width + 1
        ! Adapt stripe width so that last one doesn't get too small
        stripe_width = (thread_width-1)/stripe_count + 1
        stripe_width = ((stripe_width+3)/4)*4 ! Must be a multiple of 4 !!!
    ENDIF

    ! Determine the matrix distribution at the beginning

    ALLOCATE(limits(0:np_rows))

    CALL determine_workload(na, nbw, np_rows, limits)

    max_blk_size = MAXVAL(limits(1:np_rows) - limits(0:np_rows-1))

    a_dim2 = max_blk_size + nbw

!DEC$ ATTRIBUTES ALIGN: 64:: a
    ALLOCATE(a(stripe_width,a_dim2,stripe_count,max_threads))
    ! a(:,:,:,:) should be set to 0 in a parallel region, not here!

    ALLOCATE(row(l_nev))
    row(:) = 0

    ! Copy q from a block cyclic distribution into a distribution with contiguous rows,
    ! and transpose the matrix using stripes of given stripe_width for cache blocking.

    ! The peculiar way it is done below is due to the fact that the last row should be
    ! ready first since it is the first one to start below

    ! Please note about the OMP usage below:
    ! This is not for speed, but because we want the matrix a in the memory and
    ! in the cache of the correct thread (if possible)

!$omp parallel do private(my_thread), schedule(static, 1)
    DO my_thread = 1, max_threads
        a(:,:,:,my_thread) = 0 ! if possible, do first touch allocation!
    ENDDO

    DO ip = np_rows-1, 0, -1
        IF(my_prow == ip) THEN
            ! Receive my rows which have not yet been received
            src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
            DO i=limits(ip)+1,limits(ip+1)
                src = MOD((i-1)/nblk, np_rows)
                IF(src < my_prow) THEN
                    CALL MPI_Recv(row, l_nev, MPI_REAL8, src, 0, mpi_comm_rows, mpi_status, mpierr)
!$omp parallel do private(my_thread), schedule(static, 1)
                    DO my_thread = 1, max_threads
                        CALL unpack_row(row,i-limits(ip),my_thread)
                    ENDDO
                ELSEIF(src==my_prow) THEN
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
!$omp parallel do private(my_thread), schedule(static, 1)
                    DO my_thread = 1, max_threads
                        CALL unpack_row(row,i-limits(ip),my_thread)
                    ENDDO
                ENDIF
            ENDDO
            ! Send all rows which have not yet been send
            src_offset = 0
            DO dst = 0, ip-1
              DO i=limits(dst)+1,limits(dst+1)
                IF(MOD((i-1)/nblk, np_rows) == my_prow) THEN
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
                    CALL MPI_Send(row, l_nev, MPI_REAL8, dst, 0, mpi_comm_rows, mpierr)
                ENDIF
              ENDDO
            ENDDO
        ELSE IF(my_prow < ip) THEN
            ! Send all rows going to PE ip
            src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
            DO i=limits(ip)+1,limits(ip+1)
                src = MOD((i-1)/nblk, np_rows)
                IF(src == my_prow) THEN
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
                    CALL MPI_Send(row, l_nev, MPI_REAL8, ip, 0, mpi_comm_rows, mpierr)
                ENDIF
            ENDDO
            ! Receive all rows from PE ip
            DO i=limits(my_prow)+1,limits(my_prow+1)
                src = MOD((i-1)/nblk, np_rows)
                IF(src == ip) THEN
                    CALL MPI_Recv(row, l_nev, MPI_REAL8, src, 0, mpi_comm_rows, mpi_status, mpierr)
!$omp parallel do private(my_thread), schedule(static, 1)
                    DO my_thread = 1, max_threads
                        CALL unpack_row(row,i-limits(my_prow),my_thread)
                    ENDDO
                ENDIF
            ENDDO
        ENDIF
    ENDDO


    ! Set up result buffer queue

    num_result_blocks = ((na-1)/nblk + np_rows - my_prow) / np_rows

    num_result_buffers = 4*nfact
    ALLOCATE(result_buffer(l_nev,nblk,num_result_buffers))

    ALLOCATE(result_send_request(num_result_buffers))
    ALLOCATE(result_recv_request(num_result_buffers))
    result_send_request(:) = MPI_REQUEST_NULL
    result_recv_request(:) = MPI_REQUEST_NULL

    ! Queue up buffers

    IF(my_prow > 0 .AND. l_nev>0) THEN ! note: row 0 always sends
        DO j = 1, MIN(num_result_buffers, num_result_blocks)
            CALL MPI_Irecv(result_buffer(1,1,j), l_nev*nblk, MPI_REAL8, 0, result_recv_tag, &
                           mpi_comm_rows, result_recv_request(j), mpierr)
        ENDDO
    ENDIF

    num_bufs_recvd = 0 ! No buffers received yet

    ! Initialize top/bottom requests

    ALLOCATE(top_send_request(stripe_count))
    ALLOCATE(top_recv_request(stripe_count))
    ALLOCATE(bottom_send_request(stripe_count))
    ALLOCATE(bottom_recv_request(stripe_count))

    top_send_request(:) = MPI_REQUEST_NULL
    top_recv_request(:) = MPI_REQUEST_NULL
    bottom_send_request(:) = MPI_REQUEST_NULL
    bottom_recv_request(:) = MPI_REQUEST_NULL

    ALLOCATE(top_border_send_buffer(stripe_width*nbw*max_threads, stripe_count))
    ALLOCATE(top_border_recv_buffer(stripe_width*nbw*max_threads, stripe_count))
    ALLOCATE(bottom_border_send_buffer(stripe_width*nbw*max_threads, stripe_count))
    ALLOCATE(bottom_border_recv_buffer(stripe_width*nbw*max_threads, stripe_count))

    top_border_send_buffer(:,:) = 0
    top_border_recv_buffer(:,:) = 0
    bottom_border_send_buffer(:,:) = 0
    bottom_border_recv_buffer(:,:) = 0

    ! Initialize broadcast buffer

    ALLOCATE(bcast_buffer(nbw, max_blk_size))
    bcast_buffer = 0

    current_tv_off = 0 ! Offset of next row to be broadcast


    ! ------------------- start of work loop -------------------

    a_off = 0 ! offset in A (to avoid unnecessary shifts)

    top_msg_length = 0
    bottom_msg_length = 0

    DO sweep = 0, (na-1)/nbw

        current_n = na - sweep*nbw
        CALL determine_workload(current_n, nbw, np_rows, limits)
        current_n_start = limits(my_prow)
        current_n_end   = limits(my_prow+1)
        current_local_n = current_n_end - current_n_start

        next_n = MAX(current_n - nbw, 0)
        CALL determine_workload(next_n, nbw, np_rows, limits)
        next_n_start = limits(my_prow)
        next_n_end   = limits(my_prow+1)
        next_local_n = next_n_end - next_n_start

        IF(next_n_end < next_n) THEN
            bottom_msg_length = current_n_end - next_n_end
        ELSE
            bottom_msg_length = 0
        ENDIF

        IF(next_local_n > 0) THEN
            next_top_msg_length = current_n_start - next_n_start
        ELSE
            next_top_msg_length = 0
        ENDIF

        IF(sweep==0 .AND. current_n_end < current_n .AND. l_nev > 0) THEN
            DO i = 1, stripe_count
                csw = MIN(stripe_width, thread_width-(i-1)*stripe_width) ! "current_stripe_width"
                b_len = csw*nbw*max_threads
                CALL MPI_Irecv(bottom_border_recv_buffer(1,i), b_len, MPI_REAL8, my_prow+1, bottom_recv_tag, &
                           mpi_comm_rows, bottom_recv_request(i), mpierr)
            ENDDO
        ENDIF

        IF(current_local_n > 1) THEN
            IF(my_pcol == MOD(sweep,np_cols)) THEN
                bcast_buffer(:,1:current_local_n) = hh_trans_real(:,current_tv_off+1:current_tv_off+current_local_n)
                current_tv_off = current_tv_off + current_local_n
            ENDIF
            CALL mpi_bcast(bcast_buffer, nbw*current_local_n, MPI_REAL8, MOD(sweep,np_cols), mpi_comm_cols, mpierr)
        ELSE
            ! for current_local_n == 1 the one and only HH vector is 0 and not stored in hh_trans_real
            bcast_buffer(:,1) = 0
        ENDIF

        IF(l_nev == 0) CYCLE

        IF(current_local_n > 0) THEN

          DO i = 1, stripe_count

            ! Get real stripe width for strip i;
            ! The last OpenMP tasks may have an even smaller stripe with,
            ! but we don't care about this, i.e. we send/recv a bit too much in this case.
            ! csw: current_stripe_width

            csw = MIN(stripe_width, thread_width-(i-1)*stripe_width)

            !wait_b
            IF(current_n_end < current_n) THEN
                CALL MPI_Wait(bottom_recv_request(i), mpi_status, mpierr)
!$omp parallel do private(my_thread, n_off, b_len, b_off), schedule(static, 1)
                DO my_thread = 1, max_threads
                    n_off = current_local_n+a_off
                    b_len = csw*nbw
                    b_off = (my_thread-1)*b_len
                    a(1:csw,n_off+1:n_off+nbw,i,my_thread) = &
                      RESHAPE(bottom_border_recv_buffer(b_off+1:b_off+b_len,i), (/ csw, nbw /))
                ENDDO
                IF(next_n_end < next_n) THEN
                    CALL MPI_Irecv(bottom_border_recv_buffer(1,i), csw*nbw*max_threads, &
                                   MPI_REAL8, my_prow+1, bottom_recv_tag, &
                                   mpi_comm_rows, bottom_recv_request(i), mpierr)
                ENDIF
            ENDIF

            IF(current_local_n <= bottom_msg_length + top_msg_length) THEN

                !wait_t
                IF(top_msg_length>0) THEN
                    CALL MPI_Wait(top_recv_request(i), mpi_status, mpierr)
                ENDIF

                !compute
!$omp parallel do private(my_thread, n_off, b_len, b_off), schedule(static, 1)
                DO my_thread = 1, max_threads
                    IF(top_msg_length>0) THEN
                        b_len = csw*top_msg_length
                        b_off = (my_thread-1)*b_len
                        a(1:csw,a_off+1:a_off+top_msg_length,i,my_thread) = &
                          RESHAPE(top_border_recv_buffer(b_off+1:b_off+b_len,i), (/ csw, top_msg_length /))
                    ENDIF
                    CALL compute_hh_trafo(0, current_local_n, i, my_thread)
                ENDDO

                !send_b
                CALL MPI_Wait(bottom_send_request(i), mpi_status, mpierr)
                IF(bottom_msg_length>0) THEN
                    n_off = current_local_n+nbw-bottom_msg_length+a_off
                    b_len = csw*bottom_msg_length*max_threads
                    bottom_border_send_buffer(1:b_len,i) = &
                        RESHAPE(a(1:csw,n_off+1:n_off+bottom_msg_length,i,:), (/ b_len /))
                    CALL MPI_Isend(bottom_border_send_buffer(1,i), b_len, MPI_REAL8, my_prow+1, &
                                   top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)
                ENDIF

            ELSE

                !compute
!$omp parallel do private(my_thread, b_len, b_off), schedule(static, 1)
                DO my_thread = 1, max_threads
                    CALL compute_hh_trafo(current_local_n - bottom_msg_length, bottom_msg_length, i, my_thread)
                ENDDO

                !send_b
                CALL MPI_Wait(bottom_send_request(i), mpi_status, mpierr)
                IF(bottom_msg_length > 0) THEN
                    n_off = current_local_n+nbw-bottom_msg_length+a_off
                    b_len = csw*bottom_msg_length*max_threads
                    bottom_border_send_buffer(1:b_len,i) = &
                      RESHAPE(a(1:csw,n_off+1:n_off+bottom_msg_length,i,:), (/ b_len /))
                    CALL MPI_Isend(bottom_border_send_buffer(1,i), b_len, MPI_REAL8, my_prow+1, &
                                   top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)
                ENDIF

                !compute
!$omp parallel do private(my_thread), schedule(static, 1)
                DO my_thread = 1, max_threads
                    CALL compute_hh_trafo(top_msg_length, current_local_n-top_msg_length-bottom_msg_length, i, my_thread)
                ENDDO

                !wait_t
                IF(top_msg_length>0) THEN
                    CALL MPI_Wait(top_recv_request(i), mpi_status, mpierr)
                ENDIF

                !compute
!$omp parallel do private(my_thread, b_len, b_off), schedule(static, 1)
                DO my_thread = 1, max_threads
                    IF(top_msg_length>0) THEN
                        b_len = csw*top_msg_length
                        b_off = (my_thread-1)*b_len
                        a(1:csw,a_off+1:a_off+top_msg_length,i,my_thread) = &
                          RESHAPE(top_border_recv_buffer(b_off+1:b_off+b_len,i), (/ csw, top_msg_length /))
                    ENDIF
                    CALL compute_hh_trafo(0, top_msg_length, i, my_thread)
                ENDDO
            ENDIF

            IF(next_top_msg_length > 0) THEN
                !request top_border data
                b_len = csw*next_top_msg_length*max_threads
                CALL MPI_Irecv(top_border_recv_buffer(1,i), b_len, MPI_REAL8, my_prow-1, &
                               top_recv_tag, mpi_comm_rows, top_recv_request(i), mpierr)
            ENDIF

            !send_t
            IF(my_prow > 0) THEN
                CALL MPI_Wait(top_send_request(i), mpi_status, mpierr)
                b_len = csw*nbw*max_threads
                top_border_send_buffer(1:b_len,i) = RESHAPE(a(1:csw,a_off+1:a_off+nbw,i,:), (/ b_len /))
                CALL MPI_Isend(top_border_send_buffer(1,i), b_len, MPI_REAL8, &
                               my_prow-1, bottom_recv_tag, &
                               mpi_comm_rows, top_send_request(i), mpierr)
            ENDIF

            ! Care that there are not too many outstanding top_recv_request's
            IF(stripe_count > 1) THEN
                IF(i>1) THEN
                    CALL MPI_Wait(top_recv_request(i-1), mpi_status, mpierr)
                ELSE
                    CALL MPI_Wait(top_recv_request(stripe_count), mpi_status, mpierr)
                ENDIF
            ENDIF

          ENDDO

          top_msg_length = next_top_msg_length

        ELSE
            ! wait for last top_send_request
          DO i = 1, stripe_count
            CALL MPI_Wait(top_send_request(i), mpi_status, mpierr)
          ENDDO
        ENDIF

        ! Care about the result

        IF(my_prow == 0) THEN

            ! topmost process sends nbw rows to destination processes

            DO j=0,nfact-1

                num_blk = sweep*nfact+j ! global number of destination block, 0 based
                IF(num_blk*nblk >= na) EXIT

                nbuf = MOD(num_blk, num_result_buffers) + 1 ! buffer number to get this block

                CALL MPI_Wait(result_send_request(nbuf), mpi_status, mpierr)

                dst = MOD(num_blk, np_rows)

                IF(dst == 0) THEN
                    DO i = 1, MIN(na - num_blk*nblk, nblk)
                        CALL pack_row(row, j*nblk+i+a_off)
                        q((num_blk/np_rows)*nblk+i,1:l_nev) = row(:)
                    ENDDO
                ELSE
                    DO i = 1, nblk
                        CALL pack_row(result_buffer(:,i,nbuf),j*nblk+i+a_off)
                    ENDDO
                    CALL MPI_Isend(result_buffer(1,1,nbuf), l_nev*nblk, MPI_REAL8, dst, &
                                   result_recv_tag, mpi_comm_rows, result_send_request(nbuf), mpierr)
                ENDIF
            ENDDO

        ELSE

           ! receive and store final result

            DO j = num_bufs_recvd, num_result_blocks-1

                nbuf = MOD(j, num_result_buffers) + 1 ! buffer number to get this block

                ! If there is still work to do, just test for the next result request
                ! and leave the loop if it is not ready, otherwise wait for all
                ! outstanding requests

                IF(next_local_n > 0) THEN
                    CALL MPI_Test(result_recv_request(nbuf), flag, mpi_status, mpierr)
                    IF(.not.flag) EXIT
                ELSE
                    CALL MPI_Wait(result_recv_request(nbuf), mpi_status, mpierr)
                ENDIF

                ! Fill result buffer into q
                num_blk = j*np_rows + my_prow ! global number of current block, 0 based
                DO i = 1, MIN(na - num_blk*nblk, nblk)
                    q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
                ENDDO

                ! Queue result buffer again if there are outstanding blocks left
                IF(j+num_result_buffers < num_result_blocks) &
                    CALL MPI_Irecv(result_buffer(1,1,nbuf), l_nev*nblk, MPI_REAL8, 0, result_recv_tag, &
                                   mpi_comm_rows, result_recv_request(nbuf), mpierr)

            ENDDO
            num_bufs_recvd = j

        ENDIF

        ! Shift the remaining rows to the front of A (if necessary)

        offset = nbw - top_msg_length
        IF(offset<0) THEN
            PRINT *,'internal error, offset for shifting = ',offset
            CALL MPI_Abort(MPI_COMM_WORLD, 1, mpierr)
        ENDIF
        a_off = a_off + offset
        IF(a_off + next_local_n + nbw > a_dim2) THEN
!$omp parallel do private(my_thread, i, j), schedule(static, 1)
            DO my_thread = 1, max_threads
                DO i = 1, stripe_count
                    DO j = top_msg_length+1, top_msg_length+next_local_n
                       A(:,j,i,my_thread) = A(:,j+a_off,i,my_thread)
                    ENDDO
                ENDDO
            ENDDO
            a_off = 0
        ENDIF

    ENDDO

    ! Just for safety:
    IF(ANY(top_send_request    /= MPI_REQUEST_NULL)) PRINT *,'*** ERROR top_send_request ***',my_prow,my_pcol
    IF(ANY(bottom_send_request /= MPI_REQUEST_NULL)) PRINT *,'*** ERROR bottom_send_request ***',my_prow,my_pcol
    IF(ANY(top_recv_request    /= MPI_REQUEST_NULL)) PRINT *,'*** ERROR top_recv_request ***',my_prow,my_pcol
    IF(ANY(bottom_recv_request /= MPI_REQUEST_NULL)) PRINT *,'*** ERROR bottom_recv_request ***',my_prow,my_pcol

    IF(my_prow == 0) THEN
        ALLOCATE(mpi_statuses(MPI_STATUS_SIZE,num_result_buffers))
        CALL MPI_Waitall(num_result_buffers, result_send_request, mpi_statuses, mpierr)
        DEALLOCATE(mpi_statuses)
    ENDIF

    IF(ANY(result_send_request /= MPI_REQUEST_NULL)) PRINT *,'*** ERROR result_send_request ***',my_prow,my_pcol
    IF(ANY(result_recv_request /= MPI_REQUEST_NULL)) PRINT *,'*** ERROR result_recv_request ***',my_prow,my_pcol

    IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) &
        PRINT '(" Kernel time:",f10.3," MFlops: ",f10.3)', kernel_time, kernel_flops/kernel_time*1.d-6

    ! deallocate all working space

    DEALLOCATE(a)
    DEALLOCATE(row)
    DEALLOCATE(limits)
    DEALLOCATE(result_send_request)
    DEALLOCATE(result_recv_request)
    DEALLOCATE(top_border_send_buffer)
    DEALLOCATE(top_border_recv_buffer)
    DEALLOCATE(bottom_border_send_buffer)
    DEALLOCATE(bottom_border_recv_buffer)
    DEALLOCATE(result_buffer)
    DEALLOCATE(bcast_buffer)
    DEALLOCATE(top_send_request)
    DEALLOCATE(top_recv_request)
    DEALLOCATE(bottom_send_request)
    DEALLOCATE(bottom_recv_request)

CONTAINS

    SUBROUTINE pack_row(row, n)
    REAL*8                                   :: row(:)
    INTEGER                                  :: n

    INTEGER                                  :: i, nl, noff, nt

        DO nt = 1, max_threads
            DO i = 1, stripe_count
                noff = (nt-1)*thread_width + (i-1)*stripe_width
                nl   = MIN(stripe_width, nt*thread_width-noff, l_nev-noff)
                IF(nl<=0) EXIT
                row(noff+1:noff+nl) = a(1:nl,n,i,nt)
            ENDDO
        ENDDO

    END SUBROUTINE

    SUBROUTINE unpack_row(row, n, my_thread)

        ! Private variables in OMP regions (my_thread) should better be in the argument list!
    REAL*8, INTENT(in)                       :: row(:)
    INTEGER, INTENT(in)                      :: n, my_thread

    INTEGER                                  :: i, nl, noff

        DO i=1,stripe_count
            noff = (my_thread-1)*thread_width + (i-1)*stripe_width
            nl   = MIN(stripe_width, my_thread*thread_width-noff, l_nev-noff)
            IF(nl<=0) EXIT
            a(1:nl,n,i,my_thread) = row(noff+1:noff+nl)
        ENDDO

    END SUBROUTINE

    SUBROUTINE compute_hh_trafo(off, ncols, istripe, my_thread)

        ! Private variables in OMP regions (my_thread) should better be in the argument list!
    INTEGER, INTENT(in)                      :: off, ncols, istripe, my_thread

    INTEGER                                  :: j, nl, noff
    REAL*8                                   :: ttt, w(nbw,6)

        ttt = mpi_wtime()
        IF(istripe<stripe_count) THEN
          nl = stripe_width
        ELSE
          noff = (my_thread-1)*thread_width + (istripe-1)*stripe_width
          nl = MIN(my_thread*thread_width-noff, l_nev-noff)
          IF(nl<=0) RETURN
        ENDIF
        
        !FORTRAN CODE / X86 INRINISIC CODE / BG ASSEMBLER USING 2 HOUSEHOLDER VECTORS
        DO j = ncols, 2, -2
            w(:,1) = bcast_buffer(1:nbw,j+off)
            w(:,2) = bcast_buffer(1:nbw,j+off-1)
            CALL double_hh_trafo_2hv(a(1,j+off+a_off-1,istripe,my_thread), w, nbw, nl, stripe_width, nbw)
        ENDDO
        IF(j==1) CALL single_hh_trafo(a(1,1+off+a_off,istripe,my_thread),bcast_buffer(1,off+1), nbw, nl, stripe_width)
                
        ! X86 INTRINSIC CODE, USING 4 HOUSEHOLDER VECTORS
        !do j = ncols, 4, -4
        !    w(:,1) = bcast_buffer(1:nbw,j+off)
        !    w(:,2) = bcast_buffer(1:nbw,j+off-1)
        !    w(:,3) = bcast_buffer(1:nbw,j+off-2)
        !    w(:,4) = bcast_buffer(1:nbw,j+off-3)
        !    call double_hh_trafo_4hv(a(1,j+off+a_off-3,istripe,my_thread), w, nbw, nl, stripe_width, nbw)
        !enddo
        !do jj = j, 2, -2
        !    w(:,1) = bcast_buffer(1:nbw,jj+off)
        !    w(:,2) = bcast_buffer(1:nbw,jj+off-1)
        !    call double_hh_trafo_2hv(a(1,jj+off+a_off-1,istripe,my_thread), w, nbw, nl, stripe_width, nbw)
        !enddo
        !if(jj==1) call single_hh_trafo(a(1,1+off+a_off,istripe,my_thread),bcast_buffer(1,off+1), nbw, nl, stripe_width)

        ! X86 INTRINSIC CODE, USING 6 HOUSEHOLDER VECTORS
        !do j = ncols, 6, -6
        !    w(:,1) = bcast_buffer(1:nbw,j+off)
        !    w(:,2) = bcast_buffer(1:nbw,j+off-1)
        !    w(:,3) = bcast_buffer(1:nbw,j+off-2)
        !    w(:,4) = bcast_buffer(1:nbw,j+off-3)
        !    w(:,5) = bcast_buffer(1:nbw,j+off-4)
        !    w(:,6) = bcast_buffer(1:nbw,j+off-5)
        !    call double_hh_trafo_6hv(a(1,j+off+a_off-5,istripe,my_thread), w, nbw, nl, stripe_width, nbw)
        !enddo
        !do jj = j, 4, -4
        !    w(:,1) = bcast_buffer(1:nbw,jj+off)
        !    w(:,2) = bcast_buffer(1:nbw,jj+off-1)
        !    w(:,3) = bcast_buffer(1:nbw,jj+off-2)
        !    w(:,4) = bcast_buffer(1:nbw,jj+off-3)
        !    call double_hh_trafo_4hv(a(1,jj+off+a_off-3,istripe,my_thread), w, nbw, nl, stripe_width, nbw)
        !enddo
        !do jjj = jj, 2, -2
        !    w(:,1) = bcast_buffer(1:nbw,jjj+off)
        !    w(:,2) = bcast_buffer(1:nbw,jjj+off-1)
        !    call double_hh_trafo_2hv(a(1,jjj+off+a_off-1,istripe,my_thread), w, nbw, nl, stripe_width, nbw)
        !enddo
        !if(jjj==1) call single_hh_trafo(a(1,1+off+a_off,istripe,my_thread),bcast_buffer(1,off+1), nbw, nl, stripe_width)
        
        IF(my_thread==1) THEN
            kernel_flops = kernel_flops + 4*INT(nl,8)*INT(ncols,8)*INT(nbw,8)
            kernel_time  = kernel_time + mpi_wtime()-ttt
        ENDIF

    END SUBROUTINE

END SUBROUTINE

!-------------------------------------------------------------------------------

SUBROUTINE single_hh_trafo(q, hh, nb, nq, ldq)

    ! Perform single real Householder transformation.
    ! This routine is not performance critical and thus it is coded here in Fortran

    REAL*8                                   :: hh(*)
    INTEGER                                  :: nb, nq, ldq
    REAL*8                                   :: q(ldq, *)

    INTEGER                                  :: i
    REAL*8                                   :: v(nq)

! v = q * hh

    v(:) = q(1:nq,1)
    DO i=2,nb
        v(:) = v(:) + q(1:nq,i) * hh(i)
    ENDDO

    ! v = v * tau
    v(:) = v(:) * hh(1)

    ! q = q - v * hh**T
    q(1:nq,1) = q(1:nq,1) - v(:)
    DO i=2,nb
        q(1:nq,i) = q(1:nq,i) - v(:) * hh(i)
    ENDDO

END SUBROUTINE

!-------------------------------------------------------------------------------

SUBROUTINE determine_workload(na, nb, nprocs, limits)

    INTEGER, INTENT(in)                      :: na, nb, nprocs
    INTEGER, INTENT(out)                     :: limits(0:nprocs)

    INTEGER                                  :: i

    IF(na <= 0) THEN
        limits(:) = 0
        RETURN
    ENDIF

    IF(nb*nprocs > na) THEN
        ! there is not enough work for all
        DO i = 0, nprocs
            limits(i) = MIN(na, i*nb)
        ENDDO
    ELSE
        DO i = 0, nprocs
            limits(i) = (i*na)/nprocs
        ENDDO
    ENDIF

END SUBROUTINE

!-------------------------------------------------------------------------------

SUBROUTINE bandred_complex(na, a, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols, tmat)

!-------------------------------------------------------------------------------
!  bandred_complex: Reduces a distributed hermitian matrix to band form
!
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to Scalapack, a(:,:) must be set completely (upper and lower half)
!              a(:,:) is overwritten on exit with the band and the Householder vectors
!              in the upper half.
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith of output matrix
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  tmat(nbw,nbw,num_blocks)    where num_blocks = (na-1)/nbw + 1
!              Factors for the Householder vectors (returned), needed for back transformation
!
!-------------------------------------------------------------------------------


    INTEGER                                  :: na, lda
    COMPLEX*16                               :: a(lda,*)
    INTEGER                                  :: nblk, nbw, mpi_comm_rows, &
                                                mpi_comm_cols
    COMPLEX*16                               :: tmat(nbw,nbw,*)

    COMPLEX*16, PARAMETER                    :: CONE = (1.d0,0.d0), &
                                                CZERO = (0.d0,0.d0)

    COMPLEX*16                               :: tau, vrl, xf
    COMPLEX*16, ALLOCATABLE                  :: tmp(:,:), umc(:,:), vmr(:,:), &
                                                vr(:)
    INTEGER :: cur_pcol, i, istep, j, l_cols, l_cols_tile, l_rows, &
      l_rows_tile, lc, lce, lch, lcs, lcx, lr, lre, mpierr, my_pcol, my_prow, &
      n_cols, ncol, nlc, np_cols, np_rows, nrow, pcol, prow, tile_size
    COMPLEX*16                               :: vav(nbw,nbw), aux2(nbw), &
                                                aux1(nbw)
    REAL*8                                   :: vnorm2

   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Semibandwith nbw must be a multiple of blocksize nblk

   IF(MOD(nbw,nblk)/=0) THEN
      IF(my_prow==0 .AND. my_pcol==0) THEN
         PRINT *,'ERROR: nbw=',nbw,', nblk=',nblk
         PRINT *,'ELPA2 works only for nbw==n*nblk'
         CALL mpi_abort(mpi_comm_world,0,mpierr)
      ENDIF
   ENDIF

   ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

   tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
   tile_size = ((128*MAX(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

   l_rows_tile = tile_size/np_rows ! local rows of a tile
   l_cols_tile = tile_size/np_cols ! local cols of a tile

   DO istep = (na-1)/nbw, 1, -1

      n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

      ! Number of local columns/rows of remaining matrix
      l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
      l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)

      ! Allocate vmr and umc to their exact sizes so that they can be used in bcasts and reduces

      ALLOCATE(vmr(MAX(l_rows,1),2*n_cols))
      ALLOCATE(umc(MAX(l_cols,1),2*n_cols))

      ALLOCATE(vr(l_rows+1))

      vmr(1:l_rows,1:n_cols) = 0.
      vr(:) = 0
      tmat(:,:,istep) = 0

      ! Reduce current block to lower triangular form

      DO lc = n_cols, 1, -1

         ncol = istep*nbw + lc ! absolute column number of householder vector
         nrow = ncol - nbw ! Absolute number of pivot row

         lr  = local_index(nrow, my_prow, np_rows, nblk, -1) ! current row length
         lch = local_index(ncol, my_pcol, np_cols, nblk, -1) ! HV local column number

         tau = 0

         IF(nrow == 1) EXIT ! Nothing to do

         cur_pcol = pcol(ncol) ! Processor column owning current block

         IF(my_pcol==cur_pcol) THEN

            ! Get vector to be transformed; distribute last element and norm of
            ! remaining elements to all procs in current column

            vr(1:lr) = a(1:lr,lch) ! vector to be transformed

            IF(my_prow==prow(nrow)) THEN
               aux1(1) = DOT_PRODUCT(vr(1:lr-1),vr(1:lr-1))
               aux1(2) = vr(lr)
            ELSE
               aux1(1) = DOT_PRODUCT(vr(1:lr),vr(1:lr))
               aux1(2) = 0.
            ENDIF

            CALL mpi_allreduce(aux1,aux2,2,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)

            vnorm2 = aux2(1)
            vrl    = aux2(2)

            ! Householder transformation

            CALL hh_transform_complex(vrl, vnorm2, xf, tau)

            ! Scale vr and store Householder vector for back transformation

            vr(1:lr) = vr(1:lr) * xf
            IF(my_prow==prow(nrow)) THEN
               a(1:lr-1,lch) = vr(1:lr-1)
               a(lr,lch) = vrl
               vr(lr) = 1.
            ELSE
               a(1:lr,lch) = vr(1:lr)
            ENDIF

         ENDIF

         ! Broadcast Householder vector and tau along columns

         vr(lr+1) = tau
         CALL MPI_Bcast(vr,lr+1,MPI_DOUBLE_COMPLEX,cur_pcol,mpi_comm_cols,mpierr)
         vmr(1:lr,lc) = vr(1:lr)
         tau = vr(lr+1)
         tmat(lc,lc,istep) = CONJG(tau) ! Store tau in diagonal of tmat

         ! Transform remaining columns in current block with Householder vector

         ! Local dot product

         aux1 = 0

         nlc = 0 ! number of local columns
         DO j=1,lc-1
            lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
            IF(lcx>0) THEN
               nlc = nlc+1
               aux1(nlc) = DOT_PRODUCT(vr(1:lr),a(1:lr,lcx))
            ENDIF
         ENDDO

         ! Get global dot products
         IF(nlc>0) CALL mpi_allreduce(aux1,aux2,nlc,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)

         ! Transform

         nlc = 0
         DO j=1,lc-1
            lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
            IF(lcx>0) THEN
               nlc = nlc+1
               a(1:lr,lcx) = a(1:lr,lcx) - CONJG(tau)*aux2(nlc)*vr(1:lr)
            ENDIF
         ENDDO

      ENDDO

      ! Calculate scalar products of stored Householder vectors.
      ! This can be done in different ways, we use zherk

      vav = 0
      IF(l_rows>0) &
         CALL zherk('U','C',n_cols,l_rows,CONE,vmr,UBOUND(vmr,1),CZERO,vav,UBOUND(vav,1))
      CALL herm_matrix_allreduce(n_cols,vav,UBOUND(vav,1),mpi_comm_rows)

      ! Calculate triangular matrix T for block Householder Transformation

      DO lc=n_cols,1,-1
         tau = tmat(lc,lc,istep)
         IF(lc<n_cols) THEN
            CALL ztrmv('U','C','N',n_cols-lc,tmat(lc+1,lc+1,istep),UBOUND(tmat,1),vav(lc+1,lc),1)
            tmat(lc,lc+1:n_cols,istep) = -tau * CONJG(vav(lc+1:n_cols,lc))
         ENDIF
      ENDDO

      ! Transpose vmr -> vmc (stored in umc, second half)

      CALL elpa_transpose_vectors  (vmr, 2*UBOUND(vmr,1), mpi_comm_rows, &
                                    umc(1,n_cols+1), 2*UBOUND(umc,1), mpi_comm_cols, &
                                    1, 2*istep*nbw, n_cols, 2*nblk)

      ! Calculate umc = A**T * vmr
      ! Note that the distributed A has to be transposed
      ! Opposed to direct tridiagonalization there is no need to use the cache locality
      ! of the tiles, so we can use strips of the matrix

      umc(1:l_cols,1:n_cols) = 0.d0
      vmr(1:l_rows,n_cols+1:2*n_cols) = 0
      IF(l_cols>0 .AND. l_rows>0) THEN
         DO i=0,(istep*nbw-1)/tile_size

            lcs = i*l_cols_tile+1
            lce = MIN(l_cols,(i+1)*l_cols_tile)
            IF(lce<lcs) CYCLE

            lre = MIN(l_rows,(i+1)*l_rows_tile)
            CALL ZGEMM('C','N',lce-lcs+1,n_cols,lre,CONE,a(1,lcs),UBOUND(a,1), &
                       vmr,UBOUND(vmr,1),CONE,umc(lcs,1),UBOUND(umc,1))

            IF(i==0) CYCLE
            lre = MIN(l_rows,i*l_rows_tile)
            CALL ZGEMM('N','N',lre,n_cols,lce-lcs+1,CONE,a(1,lcs),lda, &
                       umc(lcs,n_cols+1),UBOUND(umc,1),CONE,vmr(1,n_cols+1),UBOUND(vmr,1))
         ENDDO
      ENDIF

      ! Sum up all ur(:) parts along rows and add them to the uc(:) parts
      ! on the processors containing the diagonal
      ! This is only necessary if ur has been calculated, i.e. if the
      ! global tile size is smaller than the global remaining matrix

      IF(tile_size < istep*nbw) THEN
         CALL elpa_reduce_add_vectors  (vmr(1,n_cols+1),2*UBOUND(vmr,1),mpi_comm_rows, &
                                        umc, 2*UBOUND(umc,1), mpi_comm_cols, &
                                        2*istep*nbw, n_cols, 2*nblk)
      ENDIF

      IF(l_cols>0) THEN
         ALLOCATE(tmp(l_cols,n_cols))
         CALL mpi_allreduce(umc,tmp,l_cols*n_cols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)
         umc(1:l_cols,1:n_cols) = tmp(1:l_cols,1:n_cols)
         DEALLOCATE(tmp)
      ENDIF

      ! U = U * Tmat**T

      CALL ztrmm('Right','Upper','C','Nonunit',l_cols,n_cols,CONE,tmat(1,1,istep),UBOUND(tmat,1),umc,UBOUND(umc,1))

      ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

      CALL zgemm('C','N',n_cols,n_cols,l_cols,CONE,umc,UBOUND(umc,1),umc(1,n_cols+1),UBOUND(umc,1),CZERO,vav,UBOUND(vav,1))
      CALL ztrmm('Right','Upper','C','Nonunit',n_cols,n_cols,CONE,tmat(1,1,istep),UBOUND(tmat,1),vav,UBOUND(vav,1))

      CALL herm_matrix_allreduce(n_cols,vav,UBOUND(vav,1),mpi_comm_cols)

      ! U = U - 0.5 * V * VAV
      CALL zgemm('N','N',l_cols,n_cols,n_cols,(-0.5d0,0.d0),umc(1,n_cols+1),UBOUND(umc,1),vav,UBOUND(vav,1),CONE,umc,UBOUND(umc,1))

      ! Transpose umc -> umr (stored in vmr, second half)

       CALL elpa_transpose_vectors  (umc, 2*UBOUND(umc,1), mpi_comm_cols, &
                                     vmr(1,n_cols+1), 2*UBOUND(vmr,1), mpi_comm_rows, &
                                     1, 2*istep*nbw, n_cols, 2*nblk)

      ! A = A - V*U**T - U*V**T

      DO i=0,(istep*nbw-1)/tile_size
         lcs = i*l_cols_tile+1
         lce = MIN(l_cols,(i+1)*l_cols_tile)
         lre = MIN(l_rows,(i+1)*l_rows_tile)
         IF(lce<lcs .OR. lre<1) CYCLE
         CALL zgemm('N','C',lre,lce-lcs+1,2*n_cols,-CONE, &
                    vmr,UBOUND(vmr,1),umc(lcs,1),UBOUND(umc,1), &
                    CONE,a(1,lcs),lda)
      ENDDO

      DEALLOCATE(vmr, umc, vr)

   ENDDO

END SUBROUTINE bandred_complex

!-------------------------------------------------------------------------------

SUBROUTINE herm_matrix_allreduce(n,a,lda,comm)

!-------------------------------------------------------------------------------
!  herm_matrix_allreduce: Does an mpi_allreduce for a hermitian matrix A.
!  On entry, only the upper half of A needs to be set
!  On exit, the complete matrix is set
!-------------------------------------------------------------------------------

    INTEGER                                  :: n, lda
    COMPLEX*16                               :: a(lda,*)
    INTEGER                                  :: comm

    INTEGER                                  :: i, mpierr, nc
    COMPLEX*16                               :: h2(n*n), h1(n*n)

   nc = 0
   DO i=1,n
      h1(nc+1:nc+i) = a(1:i,i)
      nc = nc+i
   ENDDO

   CALL mpi_allreduce(h1,h2,nc,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,mpierr)

   nc = 0
   DO i=1,n
      a(1:i,i) = h2(nc+1:nc+i)
      a(i,1:i-1) = CONJG(a(1:i-1,i))
      nc = nc+i
   ENDDO

END SUBROUTINE herm_matrix_allreduce

!-------------------------------------------------------------------------------

SUBROUTINE trans_ev_band_to_full_complex(na, nqc, nblk, nbw, a, lda, tmat, q, ldq, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  trans_ev_band_to_full_complex:
!  Transforms the eigenvectors of a band matrix back to the eigenvectors of the original matrix
!
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nqc         Number of columns of matrix q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith
!
!  a(lda,*)    Matrix containing the Householder vectors (i.e. matrix a after bandred_complex)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a
!
!  tmat(nbw,nbw,.) Factors returned by bandred_complex
!
!  q           On input: Eigenvectors of band matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------


    INTEGER                                  :: na, nqc, nblk, nbw, lda
    COMPLEX*16                               :: a(lda,*), tmat(nbw, nbw, *)
    INTEGER                                  :: ldq
    COMPLEX*16                               :: q(ldq,*)
    INTEGER                                  :: mpi_comm_rows, mpi_comm_cols

    COMPLEX*16, PARAMETER                    :: CONE = (1.d0,0.d0), &
                                                CZERO = (0.d0,0.d0)

    COMPLEX*16, ALLOCATABLE                  :: hvb(:), hvm(:,:), tmp1(:), &
                                                tmp2(:)
    INTEGER :: i, istep, l_colh, l_cols, l_rows, lc, max_blocks_col, &
      max_blocks_row, max_local_cols, max_local_rows, mpierr, my_pcol, &
      my_prow, n_cols, nb, ncol, np_cols, np_rows, nrow, ns, pcol, prow

   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)


   max_blocks_row = ((na -1)/nblk)/np_rows + 1  ! Rows of A
   max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q!

   max_local_rows = max_blocks_row*nblk
   max_local_cols = max_blocks_col*nblk

   ALLOCATE(tmp1(max_local_cols*nbw))
   ALLOCATE(tmp2(max_local_cols*nbw))
   ALLOCATE(hvb(max_local_rows*nbw))
   ALLOCATE(hvm(max_local_rows,nbw))

   hvm = 0   ! Must be set to 0 !!!
   hvb = 0   ! Safety only

   l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q

   DO istep=1,(na-1)/nbw

      n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

      ! Broadcast all Householder vectors for current step compressed in hvb

      nb = 0
      ns = 0

      DO lc = 1, n_cols
         ncol = istep*nbw + lc ! absolute column number of householder vector
         nrow = ncol - nbw ! absolute number of pivot row

         l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast
         l_colh = local_index(ncol  , my_pcol, np_cols, nblk, -1) ! HV local column number

         IF(my_pcol==pcol(ncol)) hvb(nb+1:nb+l_rows) = a(1:l_rows,l_colh)

         nb = nb+l_rows

         IF(lc==n_cols .OR. MOD(ncol,nblk)==0) THEN
            CALL MPI_Bcast(hvb(ns+1),nb-ns,MPI_DOUBLE_COMPLEX,pcol(ncol),mpi_comm_cols,mpierr)
            ns = nb
         ENDIF
      ENDDO

      ! Expand compressed Householder vectors into matrix hvm

      nb = 0
      DO lc = 1, n_cols
         nrow = (istep-1)*nbw+lc ! absolute number of pivot row
         l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast

         hvm(1:l_rows,lc) = hvb(nb+1:nb+l_rows)
         IF(my_prow==prow(nrow)) hvm(l_rows+1,lc) = 1.

         nb = nb+l_rows
      ENDDO

      l_rows = local_index(MIN(na,(istep+1)*nbw), my_prow, np_rows, nblk, -1)

      ! Q = Q - V * T**T * V**T * Q

      IF(l_rows>0) THEN
         CALL zgemm('C','N',n_cols,l_cols,l_rows,CONE,hvm,UBOUND(hvm,1), &
                    q,ldq,CZERO,tmp1,n_cols)
      ELSE
         tmp1(1:l_cols*n_cols) = 0
      ENDIF
      CALL mpi_allreduce(tmp1,tmp2,n_cols*l_cols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)
      IF(l_rows>0) THEN
         CALL ztrmm('L','U','C','N',n_cols,l_cols,CONE,tmat(1,1,istep),UBOUND(tmat,1),tmp2,n_cols)
         CALL zgemm('N','N',l_rows,l_cols,n_cols,-CONE,hvm,UBOUND(hvm,1), &
                    tmp2,n_cols,CONE,q,ldq)
      ENDIF

   ENDDO

   DEALLOCATE(tmp1, tmp2, hvb, hvm)


END SUBROUTINE trans_ev_band_to_full_complex

!---------------------------------------------------------------------------------------------------

SUBROUTINE tridiag_band_complex(na, nb, nblk, a, lda, d, e, mpi_comm_rows, mpi_comm_cols, mpi_comm)

!-------------------------------------------------------------------------------
! tridiag_band_complex:
! Reduces a complex hermitian band matrix to tridiagonal form
!
!  na          Order of matrix a
!
!  nb          Semi bandwith
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  a(lda,*)    Distributed system matrix reduced to banded form in the upper diagonal
!
!  lda         Leading dimension of a
!
!  d(na)       Diagonal of tridiagonal matrix, set only on PE 0 (output)
!
!  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0 (output)
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm
!              MPI-Communicator for the total processor set
!-------------------------------------------------------------------------------


    INTEGER, INTENT(in)                      :: na, nb, nblk, lda
    COMPLEX*16, INTENT(in)                   :: a(lda,*)
    REAL*8, INTENT(out)                      :: d(na), e(na)
    INTEGER, INTENT(in)                      :: mpi_comm_rows, mpi_comm_cols, &
                                                mpi_comm

    COMPLEX*16                               :: hf, tau, tau_new, x
    COMPLEX*16, ALLOCATABLE                  :: ab(:,:), hh_gath(:,:,:), &
                                                hh_send(:,:,:), hv_t(:,:), &
                                                tau_t(:)
    INTEGER :: i, iblk, ireq_ab, ireq_hv, istep, iter, local_size, &
      max_blk_size, max_threads, mpi_status(MPI_STATUS_SIZE), mpierr, &
      my_block_e, my_block_s, my_pcol, my_pe, my_prow, my_thread, n, n_off, &
      n_pes, na_s, nblocks, nblocks_total, nc, ne, np_cols, np_rows, nr, ns, &
      nt, num_chunks, num_hh_vecs, nx
    INTEGER, ALLOCATABLE :: block_limits(:), global_id(:,:), &
      global_id_tmp(:,:), hh_cnt(:), hh_dst(:), ireq_hhr(:), ireq_hhs(:), &
      limits(:), mpi_statuses(:,:), omp_block_limits(:), snd_limits(:,:)
    COMPLEX*16                               :: hv_s(nb), hv_new(nb), hv(nb), &
                                                hs(nb), hd(nb), h(nb), &
                                                ab_s(1+nb)
    REAL*8                                   :: r_a(1,1), r_ab(1,1), vnorm2

! set only on PE 0
! dummies for calling redist_band

!$  integer :: omp_get_max_threads


    CALL mpi_comm_rank(mpi_comm,my_pe,mpierr)
    CALL mpi_comm_size(mpi_comm,n_pes,mpierr)

    CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
    CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
    CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
    CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

    ! Get global_id mapping 2D procssor coordinates to global id

    ALLOCATE(global_id(0:np_rows-1,0:np_cols-1))
    ALLOCATE(global_id_tmp(0:np_rows-1,0:np_cols-1))
    global_id(:,:) = 0
    global_id(my_prow, my_pcol) = my_pe

    global_id_tmp(:,:) = global_id(:,:)
    CALL mpi_allreduce(global_id_tmp, global_id, np_rows*np_cols, mpi_integer, mpi_sum, mpi_comm, mpierr)
    DEALLOCATE(global_id_tmp)


    ! Total number of blocks in the band:

    nblocks_total = (na-1)/nb + 1

    ! Set work distribution

    ALLOCATE(block_limits(0:n_pes))
    CALL divide_band(nblocks_total, n_pes, block_limits)

    ! nblocks: the number of blocks for my task
    nblocks = block_limits(my_pe+1) - block_limits(my_pe)

    ! allocate the part of the band matrix which is needed by this PE
    ! The size is 1 block larger than needed to avoid extensive shifts
    ALLOCATE(ab(2*nb,(nblocks+1)*nb))
    ab = 0 ! needed for lower half, the extra block should also be set to 0 for safety

    ! n_off: Offset of ab within band
    n_off = block_limits(my_pe)*nb

    ! Redistribute band in a to ab
    CALL redist_band(.FALSE., r_a, a, lda, na, nblk, nb, mpi_comm_rows, mpi_comm_cols, mpi_comm, r_ab, ab)

    ! Calculate the workload for each sweep in the back transformation
    ! and the space requirements to hold the HH vectors

    ALLOCATE(limits(0:np_rows))
    CALL determine_workload(na, nb, np_rows, limits)
    max_blk_size = MAXVAL(limits(1:np_rows) - limits(0:np_rows-1))

    num_hh_vecs = 0
    num_chunks  = 0
    nx = na
    DO n = 1, nblocks_total
      CALL determine_workload(nx, nb, np_rows, limits)
      local_size = limits(my_prow+1) - limits(my_prow)
      ! add to number of householder vectors
      ! please note: for nx==1 the one and only HH vector is 0 and is neither calculated nor send below!
      IF(MOD(n-1,np_cols) == my_pcol .AND. local_size>0 .AND. nx>1) THEN
        num_hh_vecs = num_hh_vecs + local_size
        num_chunks  = num_chunks+1
      ENDIF
      nx = nx - nb
    ENDDO

    ! Allocate space for HH vectors

    ALLOCATE(hh_trans_complex(nb,num_hh_vecs))

    ! Allocate and init MPI requests

    ALLOCATE(ireq_hhr(num_chunks)) ! Recv requests
    ALLOCATE(ireq_hhs(nblocks))    ! Send requests

    num_hh_vecs = 0
    num_chunks  = 0
    nx = na
    nt = 0
    DO n = 1, nblocks_total
      CALL determine_workload(nx, nb, np_rows, limits)
      local_size = limits(my_prow+1) - limits(my_prow)
      IF(MOD(n-1,np_cols) == my_pcol .AND. local_size>0 .AND. nx>1) THEN
        num_chunks  = num_chunks+1
        CALL mpi_irecv(hh_trans_complex(1,num_hh_vecs+1), nb*local_size, MPI_COMPLEX16, nt, &
                       10+n-block_limits(nt), mpi_comm, ireq_hhr(num_chunks), mpierr)
        num_hh_vecs = num_hh_vecs + local_size
      ENDIF
      nx = nx - nb
      IF(n == block_limits(nt+1)) THEN
        nt = nt + 1
      ENDIF
    ENDDO

    ireq_hhs(:) = MPI_REQUEST_NULL

    ! Buffers for gathering/sending the HH vectors

    ALLOCATE(hh_gath(nb,max_blk_size,nblocks)) ! gathers HH vectors
    ALLOCATE(hh_send(nb,max_blk_size,nblocks)) ! send buffer for HH vectors
    hh_gath(:,:,:) = 0
    hh_send(:,:,:) = 0

    ! Some counters

    ALLOCATE(hh_cnt(nblocks))
    ALLOCATE(hh_dst(nblocks))

    hh_cnt(:) = 1 ! The first transfomation vector is always 0 and not calculated at all
    hh_dst(:) = 0 ! PE number for receive

    ireq_ab = MPI_REQUEST_NULL
    ireq_hv = MPI_REQUEST_NULL

    ! Limits for sending

    ALLOCATE(snd_limits(0:np_rows,nblocks))

    DO iblk=1,nblocks
      CALL determine_workload(na-(iblk+block_limits(my_pe)-1)*nb, nb, np_rows, snd_limits(:,iblk))
    ENDDO

    ! OpenMP work distribution:

    max_threads = 1
!$ max_threads = omp_get_max_threads()

    ! For OpenMP we need at least 2 blocks for every thread
    max_threads = MIN(max_threads, nblocks/2)
    IF(max_threads==0) max_threads = 1

    ALLOCATE(omp_block_limits(0:max_threads))

    ! Get the OpenMP block limits
    CALL divide_band(nblocks, max_threads, omp_block_limits)

    ALLOCATE(hv_t(nb,max_threads), tau_t(max_threads))
    hv_t = 0
    tau_t = 0

    ! ---------------------------------------------------------------------------
    ! Start of calculations

    na_s = block_limits(my_pe)*nb + 1

    IF(my_pe>0 .AND. na_s<=na) THEN
      ! send first column to previous PE
      ! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
      ab_s(1:nb+1) = ab(1:nb+1,na_s-n_off)
      CALL mpi_isend(ab_s,nb+1,MPI_COMPLEX16,my_pe-1,1,mpi_comm,ireq_ab,mpierr)
    ENDIF

    DO istep=1,na-1-block_limits(my_pe)*nb

      IF(my_pe==0) THEN
        n = MIN(na-na_s,nb) ! number of rows to be reduced
        hv(:) = 0
        tau = 0
        ! Transform first column of remaining matrix
        ! Opposed to the real case, the last step (istep=na-1) is needed here for making
        ! the last subdiagonal element a real number
        vnorm2 = SUM(DBLE(ab(3:n+1,na_s-n_off))**2+dimag(ab(3:n+1,na_s-n_off))**2)
        IF(n<2) vnorm2 = 0. ! Safety only
        CALL hh_transform_complex(ab(2,na_s-n_off),vnorm2,hf,tau)

        hv(1) = 1
        hv(2:n) = ab(3:n+1,na_s-n_off)*hf

        d(istep) = ab(1,na_s-n_off)
        e(istep) = ab(2,na_s-n_off)
        IF(istep == na-1) THEN
          d(na) = ab(1,na_s+1-n_off)
          e(na) = 0
        ENDIF
      ELSE
        IF(na>na_s) THEN
          ! Receive Householder vector from previous task, from PE owning subdiagonal
          CALL mpi_recv(hv,nb,MPI_COMPLEX16,my_pe-1,2,mpi_comm,mpi_status,mpierr)
          tau = hv(1)
          hv(1) = 1.
        ENDIF
      ENDIF

      na_s = na_s+1
      IF(na_s-n_off > nb) THEN
        ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
        ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0
        n_off = n_off + nb
      ENDIF

      IF(max_threads > 1) THEN

        ! Codepath for OpenMP

        ! Please note that in this case it is absolutely necessary to have at least 2 blocks per thread!
        ! Every thread is one reduction cycle behind its predecessor and thus starts one step later.
        ! This simulates the behaviour of the MPI tasks which also work after each other.
        ! The code would be considerably easier, if the MPI communication would be made within
        ! the parallel region - this is avoided here since this would require 
        ! MPI_Init_thread(MPI_THREAD_MULTIPLE) at the start of the program.

        hv_t(:,1) = hv
        tau_t(1) = tau

        DO iter = 1, 2

          ! iter=1 : work on first block
          ! iter=2 : work on remaining blocks
          ! This is done in 2 iterations so that we have a barrier in between:
          ! After the first iteration, it is guaranteed that the last row of the last block
          ! is completed by the next thread.
          ! After the first iteration it is also the place to exchange the last row
          ! with MPI calls

!$omp parallel do private(my_thread, my_block_s, my_block_e, iblk, ns, ne, hv, tau, &
!$omp nc, nr, hs, hd, vnorm2, hf, x, h, i), schedule(static,1), num_threads(max_threads)
          DO my_thread = 1, max_threads

            IF(iter == 1) THEN
              my_block_s = omp_block_limits(my_thread-1) + 1
              my_block_e = my_block_s
            ELSE
              my_block_s = omp_block_limits(my_thread-1) + 2
              my_block_e = omp_block_limits(my_thread)
            ENDIF

            DO iblk = my_block_s, my_block_e

              ns = na_s + (iblk-1)*nb - n_off - my_thread + 1 ! first column in block
              ne = ns+nb-1                    ! last column in block

              IF(istep<my_thread .OR. ns+n_off>na) EXIT

              hv = hv_t(:,my_thread)
              tau = tau_t(my_thread)

              ! Store Householder vector for back transformation

              hh_cnt(iblk) = hh_cnt(iblk) + 1

              hh_gath(1   ,hh_cnt(iblk),iblk) = tau
              hh_gath(2:nb,hh_cnt(iblk),iblk) = hv(2:nb)

              nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
              nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                            ! Note that nr>=0 implies that diagonal block is full (nc==nb)!

              ! Transform diagonal block

              CALL ZHEMV('L',nc,tau,ab(1,ns),2*nb-1,hv,1,(0.d0,0.d0),hd,1)

              x = DOT_PRODUCT(hv(1:nc),hd(1:nc))*CONJG(tau)
              hd(1:nc) = hd(1:nc) - 0.5*x*hv(1:nc)

              CALL ZHER2('L',nc,(-1.d0,0.d0),hd,1,hv,1,ab(1,ns),2*nb-1)

              hv_t(:,my_thread) = 0
              tau_t(my_thread)  = 0

              IF(nr<=0) CYCLE ! No subdiagonal block present any more

              ! Transform subdiagonal block

              CALL ZGEMV('N',nr,nb,tau,ab(nb+1,ns),2*nb-1,hv,1,(0.d0,0.d0),hs,1)

              IF(nr>1) THEN

                ! complete (old) Householder transformation for first column

                ab(nb+1:nb+nr,ns) = ab(nb+1:nb+nr,ns) - hs(1:nr) ! Note: hv(1) == 1

                ! calculate new Householder transformation for first column
                ! (stored in hv_t(:,my_thread) and tau_t(my_thread))

                vnorm2 = SUM(DBLE(ab(nb+2:nb+nr,ns))**2+dimag(ab(nb+2:nb+nr,ns))**2)
                CALL hh_transform_complex(ab(nb+1,ns),vnorm2,hf,tau_t(my_thread))
                hv_t(1   ,my_thread) = 1.
                hv_t(2:nr,my_thread) = ab(nb+2:nb+nr,ns)*hf
                ab(nb+2:,ns) = 0

                ! update subdiagonal block for old and new Householder transformation
                ! This way we can use a nonsymmetric rank 2 update which is (hopefully) faster

                CALL ZGEMV('C',nr,nb-1,tau_t(my_thread),ab(nb,ns+1),2*nb-1, &
                hv_t(1,my_thread),1,(0.d0,0.d0),h(2),1)
                x = DOT_PRODUCT(hs(1:nr),hv_t(1:nr,my_thread))* &
                tau_t(my_thread)
                h(2:nb) = h(2:nb) - x*hv(2:nb)
                ! Unfortunately there is no BLAS routine like DSYR2 for a nonsymmetric rank 2 update ("DGER2")
                DO i=2,nb
                  ab(2+nb-i:1+nb+nr-i,i+ns-1) = ab(2+nb-i:1+nb+nr-i,i+ns-1) &
                  - hv_t(1:nr,my_thread)*CONJG(h(i)) - hs(1:nr)*CONJG(hv(i))
                ENDDO

              ELSE

                ! No new Householder transformation for nr=1, just complete the old one
                ab(nb+1,ns) = ab(nb+1,ns) - hs(1) ! Note: hv(1) == 1
                DO i=2,nb
                  ab(2+nb-i,i+ns-1) = ab(2+nb-i,i+ns-1) - hs(1)*CONJG(hv(i))
                ENDDO
                ! For safety: there is one remaining dummy transformation (but tau is 0 anyways)
                hv_t(1,my_thread) = 1.

              ENDIF

            ENDDO

          ENDDO ! my_thread
!$omp end parallel do

          IF (iter==1) THEN
            ! We are at the end of the first block

            ! Send our first column to previous PE
            IF(my_pe>0 .AND. na_s <= na) THEN
              CALL mpi_wait(ireq_ab,mpi_status,mpierr)
              ab_s(1:nb+1) = ab(1:nb+1,na_s-n_off)
              CALL mpi_isend(ab_s,nb+1,MPI_COMPLEX16,my_pe-1,1,mpi_comm,ireq_ab,mpierr)
            ENDIF

            ! Request last column from next PE
            ne = na_s + nblocks*nb - (max_threads-1) - 1
            IF(istep>=max_threads .AND. ne <= na) THEN
              CALL mpi_recv(ab(1,ne-n_off),nb+1,MPI_COMPLEX16, &
              my_pe+1,1,mpi_comm,mpi_status,mpierr)
            ENDIF

          ELSE
            ! We are at the end of all blocks

            ! Send last HH vector and TAU to next PE if it has been calculated above
            ne = na_s + nblocks*nb - (max_threads-1) - 1
            IF(istep>=max_threads .AND. ne < na) THEN
              CALL mpi_wait(ireq_hv,mpi_status,mpierr)
              hv_s(1) = tau_t(max_threads)
              hv_s(2:) = hv_t(2:,max_threads)
              CALL mpi_isend(hv_s,nb,MPI_COMPLEX16,my_pe+1,2,mpi_comm,ireq_hv,mpierr)
            ENDIF

            ! "Send" HH vector and TAU to next OpenMP thread
            DO my_thread = max_threads, 2, -1
              hv_t(:,my_thread) = hv_t(:,my_thread-1)
              tau_t(my_thread)  = tau_t(my_thread-1)
            ENDDO

          ENDIF
        ENDDO ! iter

      ELSE

        ! Codepath for 1 thread without OpenMP

        ! The following code is structured in a way to keep waiting times for
        ! other PEs at a minimum, especially if there is only one block.
        ! For this reason, it requests the last column as late as possible
        ! and sends the Householder vector and the first column as early
        ! as possible.

        DO iblk=1,nblocks

          ns = na_s + (iblk-1)*nb - n_off ! first column in block
          ne = ns+nb-1                    ! last column in block

          IF(ns+n_off>na) EXIT

          ! Store Householder vector for back transformation

          hh_cnt(iblk) = hh_cnt(iblk) + 1

          hh_gath(1   ,hh_cnt(iblk),iblk) = tau
          hh_gath(2:nb,hh_cnt(iblk),iblk) = hv(2:nb)

          nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
          nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                        ! Note that nr>=0 implies that diagonal block is full (nc==nb)!

          ! Multiply diagonal block and subdiagonal block with Householder vector

          IF(iblk==nblocks .AND. nc==nb) THEN

            ! We need the last column from the next PE.
            ! First do the matrix multiplications without last column ...

            ! Diagonal block, the contribution of the last element is added below!
            ab(1,ne) = 0
            CALL ZHEMV('L',nc,tau,ab(1,ns),2*nb-1,hv,1,(0.d0,0.d0),hd,1)

            ! Subdiagonal block
            IF(nr>0) CALL ZGEMV('N',nr,nb-1,tau,ab(nb+1,ns),2*nb-1,hv,1,(0.d0,0.d0),hs,1)

            ! ... then request last column ...
            CALL mpi_recv(ab(1,ne),nb+1,MPI_COMPLEX16,my_pe+1,1,mpi_comm,mpi_status,mpierr)

            ! ... and complete the result
            hs(1:nr) = hs(1:nr) + ab(2:nr+1,ne)*tau*hv(nb)
            hd(nb) = hd(nb) + ab(1,ne)*hv(nb)*tau

          ELSE

            ! Normal matrix multiply
            CALL ZHEMV('L',nc,tau,ab(1,ns),2*nb-1,hv,1,(0.d0,0.d0),hd,1)
            IF(nr>0) CALL ZGEMV('N',nr,nb,tau,ab(nb+1,ns),2*nb-1,hv,1,(0.d0,0.d0),hs,1)

          ENDIF

          ! Calculate first column of subdiagonal block and calculate new
          ! Householder transformation for this column

          hv_new(:) = 0 ! Needed, last rows must be 0 for nr < nb
          tau_new = 0

          IF(nr>0) THEN

            ! complete (old) Householder transformation for first column

            ab(nb+1:nb+nr,ns) = ab(nb+1:nb+nr,ns) - hs(1:nr) ! Note: hv(1) == 1

            ! calculate new Householder transformation ...
            IF(nr>1) THEN
              vnorm2 = SUM(DBLE(ab(nb+2:nb+nr,ns))**2+dimag(ab(nb+2:nb+nr,ns))**2)
              CALL hh_transform_complex(ab(nb+1,ns),vnorm2,hf,tau_new)
              hv_new(1) = 1.
              hv_new(2:nr) = ab(nb+2:nb+nr,ns)*hf
              ab(nb+2:,ns) = 0
            ENDIF

            ! ... and send it away immediatly if this is the last block

            IF(iblk==nblocks) THEN
              CALL mpi_wait(ireq_hv,mpi_status,mpierr)
              hv_s(1) = tau_new
              hv_s(2:) = hv_new(2:)
              CALL mpi_isend(hv_s,nb,MPI_COMPLEX16,my_pe+1,2,mpi_comm,ireq_hv,mpierr)
            ENDIF

          ENDIF


          ! Transform diagonal block
          x = DOT_PRODUCT(hv(1:nc),hd(1:nc))*CONJG(tau)
          hd(1:nc) = hd(1:nc) - 0.5*x*hv(1:nc)

          IF(my_pe>0 .AND. iblk==1) THEN

            ! The first column of the diagonal block has to be send to the previous PE
            ! Calculate first column only ...

            ab(1:nc,ns) = ab(1:nc,ns) - hd(1:nc)*CONJG(hv(1)) - hv(1:nc)*CONJG(hd(1))

            ! ... send it away ...

            CALL mpi_wait(ireq_ab,mpi_status,mpierr)
            ab_s(1:nb+1) = ab(1:nb+1,ns)
            CALL mpi_isend(ab_s,nb+1,MPI_COMPLEX16,my_pe-1,1,mpi_comm,ireq_ab,mpierr)

            ! ... and calculate remaining columns with rank-2 update
            IF(nc>1) CALL ZHER2('L',nc-1,(-1.d0,0.d0),hd(2),1,hv(2),1,ab(1,ns+1),2*nb-1)
          ELSE
            ! No need to  send, just a rank-2 update
            CALL ZHER2('L',nc,(-1.d0,0.d0),hd,1,hv,1,ab(1,ns),2*nb-1)
          ENDIF

          ! Do the remaining double Householder transformation on the subdiagonal block cols 2 ... nb

          IF(nr>0) THEN
            IF(nr>1) THEN
              CALL ZGEMV('C',nr,nb-1,tau_new,ab(nb,ns+1),2*nb-1,hv_new,1,(0.d0,0.d0),h(2),1)
              x = DOT_PRODUCT(hs(1:nr),hv_new(1:nr))*tau_new
              h(2:nb) = h(2:nb) - x*hv(2:nb)
              ! Unfortunately there is no BLAS routine like DSYR2 for a nonsymmetric rank 2 update ("DGER2")
              DO i=2,nb
                ab(2+nb-i:1+nb+nr-i,i+ns-1) = ab(2+nb-i:1+nb+nr-i,i+ns-1) - hv_new(1:nr)*CONJG(h(i)) - hs(1:nr)*CONJG(hv(i))
              ENDDO
            ELSE
              ! No double Householder transformation for nr=1, just complete the row
              DO i=2,nb
                ab(2+nb-i,i+ns-1) = ab(2+nb-i,i+ns-1) - hs(1)*CONJG(hv(i))
              ENDDO
            ENDIF
          ENDIF

          ! Use new HH vector for the next block
          hv(:) = hv_new(:)
          tau = tau_new

        ENDDO

      ENDIF

      DO iblk = 1, nblocks

        IF(hh_dst(iblk) >= np_rows) EXIT
        IF(snd_limits(hh_dst(iblk)+1,iblk) == snd_limits(hh_dst(iblk),iblk)) EXIT

        IF(hh_cnt(iblk) == snd_limits(hh_dst(iblk)+1,iblk)-snd_limits(hh_dst(iblk),iblk)) THEN
          ! Wait for last transfer to finish
          CALL mpi_wait(ireq_hhs(iblk), mpi_status, mpierr)
          ! Copy vectors into send buffer
          hh_send(:,1:hh_cnt(iblk),iblk) = hh_gath(:,1:hh_cnt(iblk),iblk)
          ! Send to destination
          CALL mpi_isend(hh_send(1,1,iblk), nb*hh_cnt(iblk), MPI_COMPLEX16, &
                         global_id(hh_dst(iblk),MOD(iblk+block_limits(my_pe)-1,np_cols)), &
                         10+iblk, mpi_comm, ireq_hhs(iblk), mpierr)
          ! Reset counter and increase destination row
          hh_cnt(iblk) = 0
          hh_dst(iblk) = hh_dst(iblk)+1
        ENDIF

      ENDDO

    ENDDO

    ! Finish the last outstanding requests
    CALL mpi_wait(ireq_ab,mpi_status,mpierr)
    CALL mpi_wait(ireq_hv,mpi_status,mpierr)

    ALLOCATE(mpi_statuses(MPI_STATUS_SIZE,MAX(nblocks,num_chunks)))
    CALL mpi_waitall(nblocks, ireq_hhs, mpi_statuses, mpierr)
    CALL mpi_waitall(num_chunks, ireq_hhr, mpi_statuses, mpierr)
    DEALLOCATE(mpi_statuses)

    CALL mpi_barrier(mpi_comm,mpierr)

    DEALLOCATE(ab)
    DEALLOCATE(ireq_hhr, ireq_hhs)
    DEALLOCATE(hh_cnt, hh_dst)
    DEALLOCATE(hh_gath, hh_send)
    DEALLOCATE(limits, snd_limits)
    DEALLOCATE(block_limits)
    DEALLOCATE(global_id)

END SUBROUTINE tridiag_band_complex

!---------------------------------------------------------------------------------------------------

SUBROUTINE trans_ev_tridi_to_band_complex(na, nev, nblk, nbw, q, ldq, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  trans_ev_tridi_to_band_complex:
!  Transforms the eigenvectors of a tridiagonal matrix back to the eigenvectors of the band matrix
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nev         Number eigenvectors to compute (= columns of matrix q)
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nb          semi bandwith
!
!  q           On input: Eigenvectors of tridiagonal matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns/both
!
!-------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(in) :: na, nev, nblk, nbw, ldq, mpi_comm_rows, mpi_comm_cols
    COMPLEX*16 q(ldq,*)

    INTEGER np_rows, my_prow, np_cols, my_pcol

    INTEGER i, j, ip, sweep, nbuf, l_nev, a_dim2
    INTEGER current_n, current_local_n, current_n_start, current_n_end
    INTEGER next_n, next_local_n, next_n_start, next_n_end
    INTEGER bottom_msg_length, top_msg_length, next_top_msg_length
    INTEGER thread_width, stripe_width, stripe_count, csw
    INTEGER num_result_blocks, num_result_buffers, num_bufs_recvd
    INTEGER a_off, current_tv_off, max_blk_size, b_off, b_len
    INTEGER mpierr, src, src_offset, dst, offset, nfact, num_blk
    INTEGER mpi_status(MPI_STATUS_SIZE)
    LOGICAL flag

    COMPLEX*16, ALLOCATABLE :: a(:,:,:,:), row(:)
    COMPLEX*16, ALLOCATABLE :: top_border_send_buffer(:,:), top_border_recv_buffer(:,:)
    COMPLEX*16, ALLOCATABLE :: bottom_border_send_buffer(:,:), bottom_border_recv_buffer(:,:)
    COMPLEX*16, ALLOCATABLE :: result_buffer(:,:,:)
    COMPLEX*16, ALLOCATABLE :: bcast_buffer(:,:)

    INTEGER n_off
    INTEGER, ALLOCATABLE :: result_send_request(:), result_recv_request(:), limits(:)
    INTEGER, ALLOCATABLE :: top_send_request(:), bottom_send_request(:)
    INTEGER, ALLOCATABLE :: top_recv_request(:), bottom_recv_request(:)
    INTEGER, ALLOCATABLE :: mpi_statuses(:,:)

    ! MPI send/recv tags, arbitrary

    INTEGER, PARAMETER :: bottom_recv_tag = 111
    INTEGER, PARAMETER :: top_recv_tag    = 222
    INTEGER, PARAMETER :: result_recv_tag = 333

    INTEGER :: max_threads, my_thread
!$  integer :: omp_get_max_threads

    ! Just for measuring the kernel performance
    REAL*8 kernel_time
    INTEGER*8 kernel_flops


    kernel_time = 1.d-100
    kernel_flops = 0

    max_threads = 1
!$  max_threads = omp_get_max_threads()

    CALL MPI_Comm_rank(mpi_comm_rows, my_prow, mpierr)
    CALL MPI_Comm_size(mpi_comm_rows, np_rows, mpierr)
    CALL MPI_Comm_rank(mpi_comm_cols, my_pcol, mpierr)
    CALL MPI_Comm_size(mpi_comm_cols, np_cols, mpierr)

    IF(MOD(nbw,nblk)/=0) THEN
      IF(my_prow==0 .AND. my_pcol==0) THEN
         PRINT *,'ERROR: nbw=',nbw,', nblk=',nblk
         PRINT *,'band backtransform works only for nbw==n*nblk'
         CALL mpi_abort(mpi_comm_world,0,mpierr)
      ENDIF
    ENDIF

    nfact = nbw / nblk


    ! local number of eigenvectors
    l_nev = local_index(nev, my_pcol, np_cols, nblk, -1)

    IF(l_nev==0) THEN
        thread_width = 0
        stripe_width = 0
        stripe_count = 0
    ELSE
        ! Suggested stripe width is 48 - should this be reduced for the complex case ???
        thread_width = (l_nev-1)/max_threads + 1 ! number of eigenvectors per OMP thread
        stripe_width = 48 ! Must be a multiple of 4
        stripe_count = (thread_width-1)/stripe_width + 1
        ! Adapt stripe width so that last one doesn't get too small
        stripe_width = (thread_width-1)/stripe_count + 1
        stripe_width = ((stripe_width+3)/4)*4 ! Must be a multiple of 4 !!!
    ENDIF

    ! Determine the matrix distribution at the beginning

    ALLOCATE(limits(0:np_rows))

    CALL determine_workload(na, nbw, np_rows, limits)

    max_blk_size = MAXVAL(limits(1:np_rows) - limits(0:np_rows-1))

    a_dim2 = max_blk_size + nbw

    ALLOCATE(a(stripe_width,a_dim2,stripe_count,max_threads))
    ! a(:,:,:,:) should be set to 0 in a parallel region, not here!

    ALLOCATE(row(l_nev))
    row(:) = 0

    ! Copy q from a block cyclic distribution into a distribution with contiguous rows,
    ! and transpose the matrix using stripes of given stripe_width for cache blocking.

    ! The peculiar way it is done below is due to the fact that the last row should be
    ! ready first since it is the first one to start below

    ! Please note about the OMP usage below:
    ! This is not for speed, but because we want the matrix a in the memory and
    ! in the cache of the correct thread (if possible)

!$omp parallel do private(my_thread), schedule(static, 1)
    DO my_thread = 1, max_threads
        a(:,:,:,my_thread) = 0 ! if possible, do first touch allocation!
    ENDDO

    DO ip = np_rows-1, 0, -1
        IF(my_prow == ip) THEN
            ! Receive my rows which have not yet been received
            src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
            DO i=limits(ip)+1,limits(ip+1)
                src = MOD((i-1)/nblk, np_rows)
                IF(src < my_prow) THEN
                    CALL MPI_Recv(row, l_nev, MPI_COMPLEX16, src, 0, mpi_comm_rows, mpi_status, mpierr)
!$omp parallel do private(my_thread), schedule(static, 1)
                    DO my_thread = 1, max_threads
                        CALL unpack_row(row,i-limits(ip),my_thread)
                    ENDDO
                ELSEIF(src==my_prow) THEN
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
!$omp parallel do private(my_thread), schedule(static, 1)
                    DO my_thread = 1, max_threads
                        CALL unpack_row(row,i-limits(ip),my_thread)
                    ENDDO
                ENDIF
            ENDDO
            ! Send all rows which have not yet been send
            src_offset = 0
            DO dst = 0, ip-1
              DO i=limits(dst)+1,limits(dst+1)
                IF(MOD((i-1)/nblk, np_rows) == my_prow) THEN
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
                    CALL MPI_Send(row, l_nev, MPI_COMPLEX16, dst, 0, mpi_comm_rows, mpierr)
                ENDIF
              ENDDO
            ENDDO
        ELSE IF(my_prow < ip) THEN
            ! Send all rows going to PE ip
            src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
            DO i=limits(ip)+1,limits(ip+1)
                src = MOD((i-1)/nblk, np_rows)
                IF(src == my_prow) THEN
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
                    CALL MPI_Send(row, l_nev, MPI_COMPLEX16, ip, 0, mpi_comm_rows, mpierr)
                ENDIF
            ENDDO
            ! Receive all rows from PE ip
            DO i=limits(my_prow)+1,limits(my_prow+1)
                src = MOD((i-1)/nblk, np_rows)
                IF(src == ip) THEN
                    CALL MPI_Recv(row, l_nev, MPI_COMPLEX16, src, 0, mpi_comm_rows, mpi_status, mpierr)
!$omp parallel do private(my_thread), schedule(static, 1)
                    DO my_thread = 1, max_threads
                        CALL unpack_row(row,i-limits(my_prow),my_thread)
                    ENDDO
                ENDIF
            ENDDO
        ENDIF
    ENDDO


    ! Set up result buffer queue

    num_result_blocks = ((na-1)/nblk + np_rows - my_prow) / np_rows

    num_result_buffers = 4*nfact
    ALLOCATE(result_buffer(l_nev,nblk,num_result_buffers))

    ALLOCATE(result_send_request(num_result_buffers))
    ALLOCATE(result_recv_request(num_result_buffers))
    result_send_request(:) = MPI_REQUEST_NULL
    result_recv_request(:) = MPI_REQUEST_NULL

    ! Queue up buffers

    IF(my_prow > 0 .AND. l_nev>0) THEN ! note: row 0 always sends
        DO j = 1, MIN(num_result_buffers, num_result_blocks)
            CALL MPI_Irecv(result_buffer(1,1,j), l_nev*nblk, MPI_COMPLEX16, 0, result_recv_tag, &
                           mpi_comm_rows, result_recv_request(j), mpierr)
        ENDDO
    ENDIF

    num_bufs_recvd = 0 ! No buffers received yet

    ! Initialize top/bottom requests

    ALLOCATE(top_send_request(stripe_count))
    ALLOCATE(top_recv_request(stripe_count))
    ALLOCATE(bottom_send_request(stripe_count))
    ALLOCATE(bottom_recv_request(stripe_count))

    top_send_request(:) = MPI_REQUEST_NULL
    top_recv_request(:) = MPI_REQUEST_NULL
    bottom_send_request(:) = MPI_REQUEST_NULL
    bottom_recv_request(:) = MPI_REQUEST_NULL

    ALLOCATE(top_border_send_buffer(stripe_width*nbw*max_threads, stripe_count))
    ALLOCATE(top_border_recv_buffer(stripe_width*nbw*max_threads, stripe_count))
    ALLOCATE(bottom_border_send_buffer(stripe_width*nbw*max_threads, stripe_count))
    ALLOCATE(bottom_border_recv_buffer(stripe_width*nbw*max_threads, stripe_count))

    top_border_send_buffer(:,:) = 0
    top_border_recv_buffer(:,:) = 0
    bottom_border_send_buffer(:,:) = 0
    bottom_border_recv_buffer(:,:) = 0

    ! Initialize broadcast buffer

    ALLOCATE(bcast_buffer(nbw, max_blk_size))
    bcast_buffer = 0

    current_tv_off = 0 ! Offset of next row to be broadcast


    ! ------------------- start of work loop -------------------

    a_off = 0 ! offset in A (to avoid unnecessary shifts)

    top_msg_length = 0
    bottom_msg_length = 0

    DO sweep = 0, (na-1)/nbw

        current_n = na - sweep*nbw
        CALL determine_workload(current_n, nbw, np_rows, limits)
        current_n_start = limits(my_prow)
        current_n_end   = limits(my_prow+1)
        current_local_n = current_n_end - current_n_start

        next_n = MAX(current_n - nbw, 0)
        CALL determine_workload(next_n, nbw, np_rows, limits)
        next_n_start = limits(my_prow)
        next_n_end   = limits(my_prow+1)
        next_local_n = next_n_end - next_n_start

        IF(next_n_end < next_n) THEN
            bottom_msg_length = current_n_end - next_n_end
        ELSE
            bottom_msg_length = 0
        ENDIF

        IF(next_local_n > 0) THEN
            next_top_msg_length = current_n_start - next_n_start
        ELSE
            next_top_msg_length = 0
        ENDIF

        IF(sweep==0 .AND. current_n_end < current_n .AND. l_nev > 0) THEN
            DO i = 1, stripe_count
                csw = MIN(stripe_width, thread_width-(i-1)*stripe_width) ! "current_stripe_width"
                b_len = csw*nbw*max_threads
                CALL MPI_Irecv(bottom_border_recv_buffer(1,i), b_len, MPI_COMPLEX16, my_prow+1, bottom_recv_tag, &
                           mpi_comm_rows, bottom_recv_request(i), mpierr)
            ENDDO
        ENDIF

        IF(current_local_n > 1) THEN
            IF(my_pcol == MOD(sweep,np_cols)) THEN
                bcast_buffer(:,1:current_local_n) = hh_trans_complex(:,current_tv_off+1:current_tv_off+current_local_n)
                current_tv_off = current_tv_off + current_local_n
            ENDIF
            CALL mpi_bcast(bcast_buffer, nbw*current_local_n, MPI_COMPLEX16, MOD(sweep,np_cols), mpi_comm_cols, mpierr)
        ELSE
            ! for current_local_n == 1 the one and only HH vector is 0 and not stored in hh_trans_complex
            bcast_buffer(:,1) = 0
        ENDIF

        IF(l_nev == 0) CYCLE

        IF(current_local_n > 0) THEN

          DO i = 1, stripe_count

            ! Get real stripe width for strip i;
            ! The last OpenMP tasks may have an even smaller stripe with,
            ! but we don't care about this, i.e. we send/recv a bit too much in this case.
            ! csw: current_stripe_width

            csw = MIN(stripe_width, thread_width-(i-1)*stripe_width)

            !wait_b
            IF(current_n_end < current_n) THEN
                CALL MPI_Wait(bottom_recv_request(i), mpi_status, mpierr)
!$omp parallel do private(my_thread, n_off, b_len, b_off), schedule(static, 1)
                DO my_thread = 1, max_threads
                    n_off = current_local_n+a_off
                    b_len = csw*nbw
                    b_off = (my_thread-1)*b_len
                    a(1:csw,n_off+1:n_off+nbw,i,my_thread) = &
                      RESHAPE(bottom_border_recv_buffer(b_off+1:b_off+b_len,i), (/ csw, nbw /))
                ENDDO
                IF(next_n_end < next_n) THEN
                    CALL MPI_Irecv(bottom_border_recv_buffer(1,i), csw*nbw*max_threads, &
                                   MPI_COMPLEX16, my_prow+1, bottom_recv_tag, &
                                   mpi_comm_rows, bottom_recv_request(i), mpierr)
                ENDIF
            ENDIF

            IF(current_local_n <= bottom_msg_length + top_msg_length) THEN

                !wait_t
                IF(top_msg_length>0) THEN
                    CALL MPI_Wait(top_recv_request(i), mpi_status, mpierr)
                ENDIF

                !compute
!$omp parallel do private(my_thread, n_off, b_len, b_off), schedule(static, 1)
                DO my_thread = 1, max_threads
                    IF(top_msg_length>0) THEN
                        b_len = csw*top_msg_length
                        b_off = (my_thread-1)*b_len
                        a(1:csw,a_off+1:a_off+top_msg_length,i,my_thread) = &
                          RESHAPE(top_border_recv_buffer(b_off+1:b_off+b_len,i), (/ csw, top_msg_length /))
                    ENDIF
                    CALL compute_hh_trafo(0, current_local_n, i, my_thread)
                ENDDO

                !send_b
                CALL MPI_Wait(bottom_send_request(i), mpi_status, mpierr)
                IF(bottom_msg_length>0) THEN
                    n_off = current_local_n+nbw-bottom_msg_length+a_off
                    b_len = csw*bottom_msg_length*max_threads
                    bottom_border_send_buffer(1:b_len,i) = &
                        RESHAPE(a(1:csw,n_off+1:n_off+bottom_msg_length,i,:), (/ b_len /))
                    CALL MPI_Isend(bottom_border_send_buffer(1,i), b_len, MPI_COMPLEX16, my_prow+1, &
                                   top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)
                ENDIF

            ELSE

                !compute
!$omp parallel do private(my_thread, b_len, b_off), schedule(static, 1)
                DO my_thread = 1, max_threads
                    CALL compute_hh_trafo(current_local_n - bottom_msg_length, bottom_msg_length, i, my_thread)
                ENDDO

                !send_b
                CALL MPI_Wait(bottom_send_request(i), mpi_status, mpierr)
                IF(bottom_msg_length > 0) THEN
                    n_off = current_local_n+nbw-bottom_msg_length+a_off
                    b_len = csw*bottom_msg_length*max_threads
                    bottom_border_send_buffer(1:b_len,i) = &
                      RESHAPE(a(1:csw,n_off+1:n_off+bottom_msg_length,i,:), (/ b_len /))
                    CALL MPI_Isend(bottom_border_send_buffer(1,i), b_len, MPI_COMPLEX16, my_prow+1, &
                                   top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)
                ENDIF

                !compute
!$omp parallel do private(my_thread), schedule(static, 1)
                DO my_thread = 1, max_threads
                    CALL compute_hh_trafo(top_msg_length, current_local_n-top_msg_length-bottom_msg_length, i, my_thread)
                ENDDO

                !wait_t
                IF(top_msg_length>0) THEN
                    CALL MPI_Wait(top_recv_request(i), mpi_status, mpierr)
                ENDIF

                !compute
!$omp parallel do private(my_thread, b_len, b_off), schedule(static, 1)
                DO my_thread = 1, max_threads
                    IF(top_msg_length>0) THEN
                        b_len = csw*top_msg_length
                        b_off = (my_thread-1)*b_len
                        a(1:csw,a_off+1:a_off+top_msg_length,i,my_thread) = &
                          RESHAPE(top_border_recv_buffer(b_off+1:b_off+b_len,i), (/ csw, top_msg_length /))
                    ENDIF
                    CALL compute_hh_trafo(0, top_msg_length, i, my_thread)
                ENDDO
            ENDIF

            IF(next_top_msg_length > 0) THEN
                !request top_border data
                b_len = csw*next_top_msg_length*max_threads
                CALL MPI_Irecv(top_border_recv_buffer(1,i), b_len, MPI_COMPLEX16, my_prow-1, &
                               top_recv_tag, mpi_comm_rows, top_recv_request(i), mpierr)
            ENDIF

            !send_t
            IF(my_prow > 0) THEN
                CALL MPI_Wait(top_send_request(i), mpi_status, mpierr)
                b_len = csw*nbw*max_threads
                top_border_send_buffer(1:b_len,i) = RESHAPE(a(1:csw,a_off+1:a_off+nbw,i,:), (/ b_len /))
                CALL MPI_Isend(top_border_send_buffer(1,i), b_len, MPI_COMPLEX16, &
                               my_prow-1, bottom_recv_tag, &
                               mpi_comm_rows, top_send_request(i), mpierr)
            ENDIF

            ! Care that there are not too many outstanding top_recv_request's
            IF(stripe_count > 1) THEN
                IF(i>1) THEN
                    CALL MPI_Wait(top_recv_request(i-1), mpi_status, mpierr)
                ELSE
                    CALL MPI_Wait(top_recv_request(stripe_count), mpi_status, mpierr)
                ENDIF
            ENDIF

          ENDDO

          top_msg_length = next_top_msg_length

        ELSE
            ! wait for last top_send_request
          DO i = 1, stripe_count
            CALL MPI_Wait(top_send_request(i), mpi_status, mpierr)
          ENDDO
        ENDIF

        ! Care about the result

        IF(my_prow == 0) THEN

            ! topmost process sends nbw rows to destination processes

            DO j=0,nfact-1

                num_blk = sweep*nfact+j ! global number of destination block, 0 based
                IF(num_blk*nblk >= na) EXIT

                nbuf = MOD(num_blk, num_result_buffers) + 1 ! buffer number to get this block

                CALL MPI_Wait(result_send_request(nbuf), mpi_status, mpierr)

                dst = MOD(num_blk, np_rows)

                IF(dst == 0) THEN
                    DO i = 1, MIN(na - num_blk*nblk, nblk)
                        CALL pack_row(row, j*nblk+i+a_off)
                        q((num_blk/np_rows)*nblk+i,1:l_nev) = row(:)
                    ENDDO
                ELSE
                    DO i = 1, nblk
                        CALL pack_row(result_buffer(:,i,nbuf),j*nblk+i+a_off)
                    ENDDO
                    CALL MPI_Isend(result_buffer(1,1,nbuf), l_nev*nblk, MPI_COMPLEX16, dst, &
                                   result_recv_tag, mpi_comm_rows, result_send_request(nbuf), mpierr)
                ENDIF
            ENDDO

        ELSE

           ! receive and store final result

            DO j = num_bufs_recvd, num_result_blocks-1

                nbuf = MOD(j, num_result_buffers) + 1 ! buffer number to get this block

                ! If there is still work to do, just test for the next result request
                ! and leave the loop if it is not ready, otherwise wait for all
                ! outstanding requests

                IF(next_local_n > 0) THEN
                    CALL MPI_Test(result_recv_request(nbuf), flag, mpi_status, mpierr)
                    IF(.not.flag) EXIT
                ELSE
                    CALL MPI_Wait(result_recv_request(nbuf), mpi_status, mpierr)
                ENDIF

                ! Fill result buffer into q
                num_blk = j*np_rows + my_prow ! global number of current block, 0 based
                DO i = 1, MIN(na - num_blk*nblk, nblk)
                    q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
                ENDDO

                ! Queue result buffer again if there are outstanding blocks left
                IF(j+num_result_buffers < num_result_blocks) &
                    CALL MPI_Irecv(result_buffer(1,1,nbuf), l_nev*nblk, MPI_COMPLEX16, 0, result_recv_tag, &
                                   mpi_comm_rows, result_recv_request(nbuf), mpierr)

            ENDDO
            num_bufs_recvd = j

        ENDIF

        ! Shift the remaining rows to the front of A (if necessary)

        offset = nbw - top_msg_length
        IF(offset<0) THEN
            PRINT *,'internal error, offset for shifting = ',offset
            CALL MPI_Abort(MPI_COMM_WORLD, 1, mpierr)
        ENDIF
        a_off = a_off + offset
        IF(a_off + next_local_n + nbw > a_dim2) THEN
!$omp parallel do private(my_thread, i, j), schedule(static, 1)
            DO my_thread = 1, max_threads
                DO i = 1, stripe_count
                    DO j = top_msg_length+1, top_msg_length+next_local_n
                       A(:,j,i,my_thread) = A(:,j+a_off,i,my_thread)
                    ENDDO
                ENDDO
            ENDDO
            a_off = 0
        ENDIF

    ENDDO

    ! Just for safety:
    IF(ANY(top_send_request    /= MPI_REQUEST_NULL)) PRINT *,'*** ERROR top_send_request ***',my_prow,my_pcol
    IF(ANY(bottom_send_request /= MPI_REQUEST_NULL)) PRINT *,'*** ERROR bottom_send_request ***',my_prow,my_pcol
    IF(ANY(top_recv_request    /= MPI_REQUEST_NULL)) PRINT *,'*** ERROR top_recv_request ***',my_prow,my_pcol
    IF(ANY(bottom_recv_request /= MPI_REQUEST_NULL)) PRINT *,'*** ERROR bottom_recv_request ***',my_prow,my_pcol

    IF(my_prow == 0) THEN
        ALLOCATE(mpi_statuses(MPI_STATUS_SIZE,num_result_buffers))
        CALL MPI_Waitall(num_result_buffers, result_send_request, mpi_statuses, mpierr)
        DEALLOCATE(mpi_statuses)
    ENDIF

    IF(ANY(result_send_request /= MPI_REQUEST_NULL)) PRINT *,'*** ERROR result_send_request ***',my_prow,my_pcol
    IF(ANY(result_recv_request /= MPI_REQUEST_NULL)) PRINT *,'*** ERROR result_recv_request ***',my_prow,my_pcol

    IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) &
        PRINT '(" Kernel time:",f10.3," MFlops: ",f10.3)', kernel_time, kernel_flops/kernel_time*1.d-6

    ! deallocate all working space

    DEALLOCATE(a)
    DEALLOCATE(row)
    DEALLOCATE(limits)
    DEALLOCATE(result_send_request)
    DEALLOCATE(result_recv_request)
    DEALLOCATE(top_border_send_buffer)
    DEALLOCATE(top_border_recv_buffer)
    DEALLOCATE(bottom_border_send_buffer)
    DEALLOCATE(bottom_border_recv_buffer)
    DEALLOCATE(result_buffer)
    DEALLOCATE(bcast_buffer)
    DEALLOCATE(top_send_request)
    DEALLOCATE(top_recv_request)
    DEALLOCATE(bottom_send_request)
    DEALLOCATE(bottom_recv_request)

CONTAINS

    SUBROUTINE pack_row(row, n)
    COMPLEX*16                               :: row(:)
    INTEGER                                  :: n

    INTEGER                                  :: i, nl, noff, nt

        DO nt = 1, max_threads
            DO i = 1, stripe_count
                noff = (nt-1)*thread_width + (i-1)*stripe_width
                nl   = MIN(stripe_width, nt*thread_width-noff, l_nev-noff)
                IF(nl<=0) EXIT
                row(noff+1:noff+nl) = a(1:nl,n,i,nt)
            ENDDO
        ENDDO

    END SUBROUTINE

    SUBROUTINE unpack_row(row, n, my_thread)

        ! Private variables in OMP regions (my_thread) should better be in the argument list!
    COMPLEX*16, INTENT(in)                   :: row(:)
    INTEGER, INTENT(in)                      :: n, my_thread

    INTEGER                                  :: i, nl, noff

        DO i=1,stripe_count
            noff = (my_thread-1)*thread_width + (i-1)*stripe_width
            nl   = MIN(stripe_width, my_thread*thread_width-noff, l_nev-noff)
            IF(nl<=0) EXIT
            a(1:nl,n,i,my_thread) = row(noff+1:noff+nl)
        ENDDO

    END SUBROUTINE

    SUBROUTINE compute_hh_trafo(off, ncols, istripe, my_thread)

        ! Private variables in OMP regions (my_thread) should better be in the argument list!
    INTEGER, INTENT(in)                      :: off, ncols, istripe, my_thread

    INTEGER                                  :: j, nl, noff
    REAL*8                                   :: ttt

        ttt = mpi_wtime()
        IF(istripe<stripe_count) THEN
          nl = stripe_width
        ELSE
          noff = (my_thread-1)*thread_width + (istripe-1)*stripe_width
          nl = MIN(my_thread*thread_width-noff, l_nev-noff)
          IF(nl<=0) RETURN
        ENDIF
        DO j = ncols, 1, -1
          CALL single_hh_trafo_complex(a(1,j+off+a_off,istripe,my_thread),bcast_buffer(1,j+off),nbw,nl,stripe_width)
        ENDDO
        IF(my_thread==1) THEN
          kernel_flops = kernel_flops + 4*4*INT(nl,8)*INT(ncols,8)*INT(nbw,8)
          kernel_time  = kernel_time + mpi_wtime()-ttt
        ENDIF

    END SUBROUTINE

END SUBROUTINE

! --------------------------------------------------------------------------------------------------
! redist_band: redistributes band from 2D block cyclic form to 1D band

SUBROUTINE redist_band(l_real, r_a, c_a, lda, na, nblk, nbw, mpi_comm_rows, mpi_comm_cols, mpi_comm, r_ab, c_ab)

    LOGICAL, INTENT(in)                      :: l_real
    INTEGER, INTENT(in)                      :: lda
    COMPLEX*16, INTENT(in)                   :: c_a(lda, *)
    REAL*8, INTENT(in)                       :: r_a(lda, *)
    INTEGER, INTENT(in)                      :: na, nblk, nbw, mpi_comm_rows, &
                                                mpi_comm_cols, mpi_comm
    REAL*8, INTENT(out)                      :: r_ab(:,:)
    COMPLEX*16, INTENT(out)                  :: c_ab(:,:)

    COMPLEX*16, ALLOCATABLE                  :: c_buf(:,:), c_rbuf(:,:,:), &
                                                c_sbuf(:,:,:)
    INTEGER :: i, il, is, j, jl, js, l_cols, l_rows, mpierr, my_pcol, my_pe, &
      my_prow, n_off, n_pes, nblocks_total, nfact, np, np_cols, np_rows, npc, &
      npr
    INTEGER, ALLOCATABLE :: block_limits(:), global_id(:,:), &
      global_id_tmp(:,:), ncnt_r(:), ncnt_s(:), nstart_r(:), nstart_s(:)
    REAL*8, ALLOCATABLE                      :: r_buf(:,:), r_rbuf(:,:,:), &
                                                r_sbuf(:,:,:)

   CALL mpi_comm_rank(mpi_comm,my_pe,mpierr)
   CALL mpi_comm_size(mpi_comm,n_pes,mpierr)

   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Get global_id mapping 2D procssor coordinates to global id

   ALLOCATE(global_id(0:np_rows-1,0:np_cols-1))
   ALLOCATE(global_id_tmp(0:np_rows-1,0:np_cols-1))
   global_id(:,:) = 0
   global_id(my_prow, my_pcol) = my_pe

   global_id_tmp(:,:) = global_id(:,:)
   CALL mpi_allreduce(global_id_tmp, global_id, np_rows*np_cols, mpi_integer, mpi_sum, mpi_comm, mpierr)
   DEALLOCATE(global_id_tmp)


   ! Set work distribution

   nblocks_total = (na-1)/nbw + 1

   ALLOCATE(block_limits(0:n_pes))
   CALL divide_band(nblocks_total, n_pes, block_limits)


   ALLOCATE(ncnt_s(0:n_pes-1))
   ALLOCATE(nstart_s(0:n_pes-1))
   ALLOCATE(ncnt_r(0:n_pes-1))
   ALLOCATE(nstart_r(0:n_pes-1))


   nfact = nbw/nblk

   ! Count how many blocks go to which PE

   ncnt_s(:) = 0
   np = 0 ! receiver PE number
   DO j=0,(na-1)/nblk ! loop over rows of blocks
      IF(j/nfact==block_limits(np+1)) np = np+1
      IF(MOD(j,np_rows) == my_prow) THEN
         DO i=0,nfact
            IF(MOD(i+j,np_cols) == my_pcol) THEN
               ncnt_s(np) = ncnt_s(np) + 1
            ENDIF
         ENDDO
      ENDIF
   ENDDO

   ! Allocate send buffer

   IF(l_real) THEN
      ALLOCATE(r_sbuf(nblk,nblk,SUM(ncnt_s)))
      r_sbuf(:,:,:) = 0.
   ELSE
      ALLOCATE(c_sbuf(nblk,nblk,SUM(ncnt_s)))
      c_sbuf(:,:,:) = 0.
   ENDIF

   ! Determine start offsets in send buffer

   nstart_s(0) = 0
   DO i=1,n_pes-1
      nstart_s(i) = nstart_s(i-1) + ncnt_s(i-1)
   ENDDO

   ! Fill send buffer

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of a

   np = 0
   DO j=0,(na-1)/nblk ! loop over rows of blocks
      IF(j/nfact==block_limits(np+1)) np = np+1
      IF(MOD(j,np_rows) == my_prow) THEN
         DO i=0,nfact
            IF(MOD(i+j,np_cols) == my_pcol) THEN
               nstart_s(np) = nstart_s(np) + 1
               js = (j/np_rows)*nblk
               is = ((i+j)/np_cols)*nblk
               jl = MIN(nblk,l_rows-js)
               il = MIN(nblk,l_cols-is)
               IF(l_real) THEN
                  r_sbuf(1:jl,1:il,nstart_s(np)) = r_a(js+1:js+jl,is+1:is+il)
               ELSE
                  c_sbuf(1:jl,1:il,nstart_s(np)) = c_a(js+1:js+jl,is+1:is+il)
               ENDIF
            ENDIF
         ENDDO
      ENDIF
   ENDDO

   ! Count how many blocks we get from which PE

   ncnt_r(:) = 0
   DO j=block_limits(my_pe)*nfact,MIN(block_limits(my_pe+1)*nfact-1,(na-1)/nblk)
      npr = MOD(j,np_rows)
      DO i=0,nfact
         npc = MOD(i+j,np_cols)
         np = global_id(npr,npc)
         ncnt_r(np) = ncnt_r(np) + 1
      ENDDO
   ENDDO

   ! Allocate receive buffer

   IF(l_real) THEN
      ALLOCATE(r_rbuf(nblk,nblk,SUM(ncnt_r)))
   ELSE
      ALLOCATE(c_rbuf(nblk,nblk,SUM(ncnt_r)))
   ENDIF

   ! Set send counts/send offsets, receive counts/receive offsets
   ! now actually in variables, not in blocks

   ncnt_s(:) = ncnt_s(:)*nblk*nblk

   nstart_s(0) = 0
   DO i=1,n_pes-1
      nstart_s(i) = nstart_s(i-1) + ncnt_s(i-1)
   ENDDO

   ncnt_r(:) = ncnt_r(:)*nblk*nblk

   nstart_r(0) = 0
   DO i=1,n_pes-1
      nstart_r(i) = nstart_r(i-1) + ncnt_r(i-1)
   ENDDO

   ! Exchange all data with MPI_Alltoallv

   IF(l_real) THEN
      CALL MPI_Alltoallv(r_sbuf,ncnt_s,nstart_s,MPI_REAL8,r_rbuf,ncnt_r,nstart_r,MPI_REAL8,mpi_comm,mpierr)
   ELSE
      CALL MPI_Alltoallv(c_sbuf,ncnt_s,nstart_s,MPI_COMPLEX16,c_rbuf,ncnt_r,nstart_r,MPI_COMPLEX16,mpi_comm,mpierr)
   ENDIF

   ! set band from receive buffer

   ncnt_r(:) = ncnt_r(:)/(nblk*nblk)

   nstart_r(0) = 0
   DO i=1,n_pes-1
      nstart_r(i) = nstart_r(i-1) + ncnt_r(i-1)
   ENDDO

   IF(l_real) THEN
      ALLOCATE(r_buf((nfact+1)*nblk,nblk))
   ELSE
      ALLOCATE(c_buf((nfact+1)*nblk,nblk))
   ENDIF

   ! n_off: Offset of ab within band
   n_off = block_limits(my_pe)*nbw

   DO j=block_limits(my_pe)*nfact,MIN(block_limits(my_pe+1)*nfact-1,(na-1)/nblk)
      npr = MOD(j,np_rows)
      DO i=0,nfact
         npc = MOD(i+j,np_cols)
         np = global_id(npr,npc)
         nstart_r(np) = nstart_r(np) + 1
         IF(l_real) THEN
            r_buf(i*nblk+1:i*nblk+nblk,:) = TRANSPOSE(r_rbuf(:,:,nstart_r(np)))
         ELSE
            c_buf(i*nblk+1:i*nblk+nblk,:) = CONJG(TRANSPOSE(c_rbuf(:,:,nstart_r(np))))
         ENDIF
      ENDDO
      DO i=1,MIN(nblk,na-j*nblk)
         IF(l_real) THEN
            r_ab(1:nbw+1,i+j*nblk-n_off) = r_buf(i:i+nbw,i)
         ELSE
            c_ab(1:nbw+1,i+j*nblk-n_off) = c_buf(i:i+nbw,i)
         ENDIF
      ENDDO
   ENDDO

   DEALLOCATE(ncnt_s, nstart_s)
   DEALLOCATE(ncnt_r, nstart_r)
   DEALLOCATE(global_id)
   DEALLOCATE(block_limits)
   IF(l_real) THEN
      DEALLOCATE(r_sbuf, r_rbuf, r_buf)
   ELSE
      DEALLOCATE(c_sbuf, c_rbuf, c_buf)
   ENDIF

END SUBROUTINE

!---------------------------------------------------------------------------------------------------
! divide_band: sets the work distribution in band
! Proc n works on blocks block_limits(n)+1 .. block_limits(n+1)

SUBROUTINE divide_band(nblocks_total, n_pes, block_limits)

    INTEGER, INTENT(in)                      :: nblocks_total, n_pes
    INTEGER, INTENT(out)                     :: block_limits(0:n_pes)

    INTEGER                                  :: n, nblocks, nblocks_left

! total number of blocks in band
! number of PEs for division

   block_limits(0) = 0
   IF(nblocks_total < n_pes) THEN
      ! Not enough work for all: The first tasks get exactly 1 block
      DO n=1,n_pes
         block_limits(n) = MIN(nblocks_total,n)
      ENDDO
   ELSE
      ! Enough work for all. If there is no exact loadbalance,
      ! the LAST tasks get more work since they are finishing earlier!
      nblocks = nblocks_total/n_pes
      nblocks_left = nblocks_total - n_pes*nblocks
      DO n=1,n_pes
         IF(n<=n_pes-nblocks_left) THEN
            block_limits(n) = block_limits(n-1) + nblocks
         ELSE
            block_limits(n) = block_limits(n-1) + nblocks + 1
         ENDIF
      ENDDO
   ENDIF

END SUBROUTINE

!-------------------------------------------------------------------------------

SUBROUTINE band_band_real(na, nb, nb2, ab, ab2, d, e, mpi_comm)

!-------------------------------------------------------------------------------
! band_band_real:
! Reduces a real symmetric banded matrix to a real symmetric matrix with smaller bandwidth. Householder transformations are not stored.
! Matrix size na and original bandwidth nb have to be a multiple of the target bandwidth nb2. (Hint: expand your matrix with zero entries, if this 
! requirement doesn't hold)
!
!  na          Order of matrix
!
!  nb          Semi bandwidth of original matrix
!
!  nb2         Semi bandwidth of target matrix
!
!  ab          Input matrix with bandwidth nb. The leading dimension of the banded matrix has to be 2*nb. The parallel data layout 
!              has to be accordant to divide_band(), i.e. the matrix columns block_limits(n)*nb+1 to min(na, block_limits(n+1)*nb) 
!              are located on rank n.
!
!  ab2         Output matrix with bandwidth nb2. The leading dimension of the banded matrix is 2*nb2. The parallel data layout is
!              accordant to divide_band(), i.e. the matrix columns block_limits(n)*nb2+1 to min(na, block_limits(n+1)*nb2) are located
!              on rank n.
!
!  d(na)       Diagonal of tridiagonal matrix, set only on PE 0, set only if ab2 = 1 (output)
!
!  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0, set only if ab2 = 1 (output)
!
!  mpi_comm
!              MPI-Communicator for the total processor set
!-------------------------------------------------------------------------------


    INTEGER, INTENT(in)                      :: na, nb, nb2
    REAL*8, INTENT(inout)                    :: ab(2*nb,*), ab2(2*nb2,*)
    REAL*8, INTENT(out)                      :: d(na), e(na)
    INTEGER, INTENT(in)                      :: mpi_comm

    INTEGER :: dest, i, iblk, info, ireq_ab, ireq_hv, istep, lwork, &
      mpi_status(MPI_STATUS_SIZE), mpierr, my_pe, n, n_off, n_pes, na_s, &
      nblocks, nblocks2, nblocks_total, nblocks_total2, nc, ne, nr, ns
    INTEGER, ALLOCATABLE                     :: block_limits(:), &
                                                block_limits2(:), &
                                                ireq_ab2(:), mpi_statuses(:,:)
    REAL*8 :: ab_r(1+nb,nb2), ab_s(1+nb,nb2), ab_s2(2*nb2,nb2), hv(nb,nb2), &
      hv_new(nb,nb2), hv_s(nb,nb2), tau(nb2), tau_new(nb2), w(nb,nb2), &
      w_new(nb,nb2), work(nb*nb2), work2(nb2*nb2)

! set only on PE 0
!----------------
!----------------

   CALL mpi_comm_rank(mpi_comm,my_pe,mpierr)
   CALL mpi_comm_size(mpi_comm,n_pes,mpierr)

   ! Total number of blocks in the band:
   nblocks_total = (na-1)/nb + 1
   nblocks_total2 = (na-1)/nb2 + 1
   
   ! Set work distribution
   ALLOCATE(block_limits(0:n_pes))
   CALL divide_band(nblocks_total, n_pes, block_limits)
   
   ALLOCATE(block_limits2(0:n_pes))
   CALL divide_band(nblocks_total2, n_pes, block_limits2)

   ! nblocks: the number of blocks for my task
   nblocks = block_limits(my_pe+1) - block_limits(my_pe)
   nblocks2 = block_limits2(my_pe+1) - block_limits2(my_pe)
   
   ALLOCATE(ireq_ab2(1:nblocks2))
   ireq_ab2 = MPI_REQUEST_NULL
   IF(nb2>1) THEN
       DO i=0,nblocks2-1
           CALL mpi_irecv(ab2(1,i*nb2+1),2*nb2*nb2,mpi_real8,0,3,mpi_comm,ireq_ab2(i+1),mpierr)
       ENDDO
   ENDIF

   ! n_off: Offset of ab within band
   n_off = block_limits(my_pe)*nb
   lwork = nb*nb2
   dest = 0

   ireq_ab = MPI_REQUEST_NULL
   ireq_hv = MPI_REQUEST_NULL
   
   ! ---------------------------------------------------------------------------
   ! Start of calculations

   na_s = block_limits(my_pe)*nb + 1

   IF(my_pe>0 .AND. na_s<=na) THEN
      ! send first nb2 columns to previous PE
      ! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
      DO i=1,nb2
              ab_s(1:nb+1,i) = ab(1:nb+1,na_s-n_off+i-1)
      ENDDO
      CALL mpi_isend(ab_s,(nb+1)*nb2,mpi_real8,my_pe-1,1,mpi_comm,ireq_ab,mpierr)
   ENDIF
   
   DO istep=1,na/nb2
   
      IF(my_pe==0) THEN
      
         n = MIN(na-na_s-nb2+1,nb) ! number of rows to be reduced
         hv(:,:) = 0
         tau(:) = 0
         
         ! The last step (istep=na-1) is only needed for sending the last HH vectors.
         ! We don't want the sign of the last element flipped (analogous to the other sweeps)
         IF(istep < na/nb2) THEN
            
            ! Transform first block column of remaining matrix
            CALL dgeqrf(n, nb2, ab(1+nb2,na_s-n_off), 2*nb-1, tau, work, lwork, info);
                        
            DO i=1,nb2
                    hv(i,i) = 1.0
                    hv(i+1:n,i) = ab(1+nb2+1:1+nb2+n-i,na_s-n_off+i-1)
                    ab(1+nb2+1:2*nb,na_s-n_off+i-1) = 0
            ENDDO
            
         ENDIF
         
         IF(nb2==1) THEN
            d(istep) = ab(1,na_s-n_off)
            e(istep) = ab(2,na_s-n_off)
            IF(istep == na) THEN
                    e(na) = 0
            ENDIF
         ELSE
            ab_s2 = 0
            ab_s2(:,:) = ab(1:nb2+1,na_s-n_off:na_s-n_off+nb2-1)
            IF(block_limits2(dest+1)<istep) THEN
                    dest = dest+1
            ENDIF
            CALL mpi_send(ab_s2,2*nb2*nb2,mpi_real8,dest,3,mpi_comm,mpierr)
         ENDIF
         
      ELSE
         IF(na>na_s+nb2-1) THEN
            ! Receive Householder vectors from previous task, from PE owning subdiagonal
            CALL mpi_recv(hv,nb*nb2,mpi_real8,my_pe-1,2,mpi_comm,mpi_status,mpierr)
            DO i=1,nb2
                       tau(i) = hv(i,i)
                       hv(i,i) = 1.
            ENDDO
         ENDIF
      ENDIF
      
      na_s = na_s+nb2
      IF(na_s-n_off > nb) THEN
         ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
         ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0
         n_off = n_off + nb
      ENDIF
      
      DO iblk=1,nblocks
         ns = na_s + (iblk-1)*nb - n_off ! first column in block
         ne = ns+nb-nb2                    ! last column in block

         IF(ns+n_off>na) EXIT
         
         nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
         nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                       ! Note that nr>=0 implies that diagonal block is full (nc==nb)!
                                       
         CALL wy_gen(nc,nb2,w,hv,tau,work,nb)

         IF(iblk==nblocks .AND. nc==nb) THEN
             !request last nb2 columns
             CALL mpi_recv(ab_r,(nb+1)*nb2,mpi_real8,my_pe+1,1,mpi_comm,mpi_status,mpierr)
             DO i=1,nb2
                 ab(1:nb+1,ne+i-1) = ab_r(:,i)
             ENDDO
         ENDIF
         
         hv_new(:,:) = 0 ! Needed, last rows must be 0 for nr < nb
         tau_new(:) = 0
         
         IF(nr>0) THEN
             CALL wy_right(nr,nb,nb2,ab(nb+1,ns),2*nb-1,w,hv,work,nb)
             
             CALL dgeqrf(nr,nb2,ab(nb+1,ns),2*nb-1,tau_new,work,lwork,info);
                          
             DO i=1,nb2
                      hv_new(i,i) = 1.0
                      hv_new(i+1:,i) = ab(nb+2:2*nb-i+1,ns+i-1)
                      ab(nb+2:,ns+i-1) = 0
             ENDDO
             
             !send hh-vector
             IF(iblk==nblocks) THEN
                 CALL mpi_wait(ireq_hv,mpi_status,mpierr)
                 hv_s = hv_new
                 DO i=1,nb2
                     hv_s(i,i) = tau_new(i)
                 ENDDO
                 CALL mpi_isend(hv_s,nb*nb2,mpi_real8,my_pe+1,2,mpi_comm,ireq_hv,mpierr)
             ENDIF
             
         ENDIF
         
         CALL wy_symm(nc,nb2,ab(1,ns),2*nb-1,w,hv,work,work2,nb)
         
         IF(my_pe>0 .AND. iblk==1) THEN
             !send first nb2 columns to previous PE
             CALL mpi_wait(ireq_ab,mpi_status,mpierr)
             DO i=1,nb2
                 ab_s(1:nb+1,i) = ab(1:nb+1,ns+i-1)
             ENDDO
             CALL mpi_isend(ab_s,(nb+1)*nb2,mpi_real8,my_pe-1,1,mpi_comm,ireq_ab,mpierr)
         ENDIF
         
         IF(nr>0) THEN
             CALL wy_gen(nr,nb2,w_new,hv_new,tau_new,work,nb)
             CALL wy_left(nb-nb2,nr,nb2,ab(nb+1-nb2,ns+nb2),2*nb-1,w_new,hv_new,work,nb)
         ENDIF
         
         ! Use new HH vector for the next block
         hv(:,:) = hv_new(:,:)
         tau = tau_new
         
     ENDDO

   ENDDO

   ! Finish the last outstanding requests
   CALL mpi_wait(ireq_ab,mpi_status,mpierr)
   CALL mpi_wait(ireq_hv,mpi_status,mpierr)
   ALLOCATE(mpi_statuses(MPI_STATUS_SIZE,nblocks2))
   CALL mpi_waitall(nblocks2,ireq_ab2,mpi_statuses,mpierr)
   DEALLOCATE(mpi_statuses)

   CALL mpi_barrier(mpi_comm,mpierr)

   DEALLOCATE(block_limits)
   DEALLOCATE(block_limits2)
   DEALLOCATE(ireq_ab2)

END SUBROUTINE

! --------------------------------------------------------------------------------------------------

SUBROUTINE wy_gen(n, nb, W, Y, tau, mem, lda)
    
    INTEGER, INTENT(in)                      :: n, nb
    REAL*8, INTENT(in)                       :: tau(nb), mem(nb)
    INTEGER, INTENT(in)                      :: lda
    REAL*8, INTENT(in)                       :: Y(lda,nb)
    REAL*8, INTENT(out)                      :: W(lda,nb)

    INTEGER                                  :: i

!length of householder-vectors
!number of householder-vectors
!leading dimension of Y and W
!matrix containing nb householder-vectors of length b
!tau values
!output matrix W
!memory for a temporary matrix of size nb

    W(1:n,1) = tau(1)*Y(1:n,1)
    DO i=2,nb
        W(1:n,i) = tau(i)*Y(1:n,i)
        CALL DGEMV('T',n,i-1,1.d0,Y,lda,W(1,i),1,0.d0,mem,1)
        CALL DGEMV('N',n,i-1,-1.d0,W,lda,mem,1,1.d0,W(1,i),1)
    ENDDO
    
END SUBROUTINE

! --------------------------------------------------------------------------------------------------

SUBROUTINE wy_left(n, m, nb, A, lda, W, Y, mem, lda2)

    INTEGER, INTENT(in)                      :: n, m, nb, lda
    REAL*8, INTENT(inout)                    :: A(lda,*)
    REAL*8, INTENT(in)                       :: W(m,nb), Y(m,nb)
    REAL*8, INTENT(inout)                    :: mem(n,nb)
    INTEGER, INTENT(in)                      :: lda2

!width of the matrix A
!length of matrix W and Y
!width of matrix W and Y
!leading dimension of A
!leading dimension of W and Y
!matrix to be transformed
!blocked transformation matrix W
!blocked transformation matrix Y
!memory for a temporary matrix of size n x nb

    CALL DGEMM('T', 'N', nb, n, m, 1.d0, W, lda2, A, lda, 0.d0, mem, nb)
    CALL DGEMM('N', 'N', m, n, nb, -1.d0, Y, lda2, mem, nb, 1.d0, A, lda)

END SUBROUTINE

! --------------------------------------------------------------------------------------------------

SUBROUTINE wy_right(n, m, nb, A, lda, W, Y, mem, lda2)

    INTEGER, INTENT(in)                      :: n, m, nb, lda
    REAL*8, INTENT(inout)                    :: A(lda,*)
    REAL*8, INTENT(in)                       :: W(m,nb), Y(m,nb)
    REAL*8, INTENT(inout)                    :: mem(n,nb)
    INTEGER, INTENT(in)                      :: lda2

!height of the matrix A
!length of matrix W and Y
!width of matrix W and Y
!leading dimension of A
!leading dimension of W and Y
!matrix to be transformed
!blocked transformation matrix W
!blocked transformation matrix Y
!memory for a temporary matrix of size n x nb

    CALL DGEMM('N', 'N', n, nb, m, 1.d0, A, lda, W, lda2, 0.d0, mem, n)
    CALL DGEMM('N', 'T', n, m, nb, -1.d0, mem, n, Y, lda2, 1.d0, A, lda)

END SUBROUTINE

! --------------------------------------------------------------------------------------------------

SUBROUTINE wy_symm(n, nb, A, lda, W, Y, mem, mem2, lda2)

    INTEGER, INTENT(in)                      :: n, nb, lda
    REAL*8, INTENT(inout)                    :: A(lda,*)
    REAL*8, INTENT(in)                       :: W(n,nb), Y(n,nb)
    REAL*8                                   :: mem(n,nb), mem2(nb,nb)
    INTEGER, INTENT(in)                      :: lda2

!width/heigth of the matrix A; length of matrix W and Y
!width of matrix W and Y
!leading dimension of A
!leading dimension of W and Y
!matrix to be transformed
!blocked transformation matrix W
!blocked transformation matrix Y
!memory for a temporary matrix of size n x nb
!memory for a temporary matrix of size nb x nb

    CALL DSYMM('L', 'L', n, nb, 1.d0, A, lda, W, lda2, 0.d0, mem, n)
    CALL DGEMM('T', 'N', nb, nb, n, 1.d0, mem, n, W, lda2, 0.d0, mem2, nb)
    CALL DGEMM('N', 'N', n, nb, nb, -0.5d0, Y, lda2, mem2, nb, 1.d0, mem, n)
    CALL DSYR2K('L', 'N', n, nb, -1.d0, Y, lda2, mem, n, 1.d0, A, lda)

END SUBROUTINE

! --------------------------------------------------------------------------------------------------

#endif 
END MODULE elpa2
