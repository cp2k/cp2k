! ELPA1 -- Faster replacements for ScaLAPACK symmetric eigenvalue routines
! 
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".

MODULE elpa1

! Version 1.1.2, 2011-02-21

  IMPLICIT NONE

  PRIVATE ! By default, all routines contained are private

#if defined(__parallel) 
  INCLUDE "mpif.h"


  ! The following routines are public:

  PUBLIC :: get_elpa_row_col_comms     ! Sets MPI row/col communicators

  PUBLIC :: solve_evp_real             ! Driver routine for real eigenvalue problem
  PUBLIC :: solve_evp_complex          ! Driver routine for complex eigenvalue problem

  PUBLIC :: tridiag_real               ! Transform real symmetric matrix to tridiagonal form
  PUBLIC :: trans_ev_real              ! Transform eigenvectors of a tridiagonal matrix back
  PUBLIC :: mult_at_b_real             ! Multiply real matrices A**T * B

  PUBLIC :: tridiag_complex            ! Transform complex hermitian matrix to tridiagonal form
  PUBLIC :: trans_ev_complex           ! Transform eigenvectors of a tridiagonal matrix back
  PUBLIC :: mult_ah_b_complex          ! Multiply complex matrices A**H * B

  PUBLIC :: solve_tridi                ! Solve tridiagonal eigensystem with divide and conquer method

  PUBLIC :: cholesky_real              ! Cholesky factorization of a real matrix
  PUBLIC :: invert_trm_real            ! Invert real triangular matrix

  PUBLIC :: cholesky_complex           ! Cholesky factorization of a complex matrix
  PUBLIC :: invert_trm_complex         ! Invert complex triangular matrix

  PUBLIC :: local_index                ! Get local index of a block cyclic distributed matrix
  PUBLIC :: least_common_multiple      ! Get least common multiple

  PUBLIC :: hh_transform_real
  PUBLIC :: hh_transform_complex

!-------------------------------------------------------------------------------

  ! Timing results, set by every call to solve_evp_xxx

  REAL*8, PUBLIC :: time_evp_fwd    ! forward transformations (to tridiagonal form)
  REAL*8, PUBLIC :: time_evp_solve  ! time for solving the tridiagonal system
  REAL*8, PUBLIC :: time_evp_back   ! time for back transformations of eigenvectors

  ! Set elpa_print_times to .true. for explicit timing outputs

  LOGICAL, PUBLIC :: elpa_print_times = .FALSE.

!-------------------------------------------------------------------------------
CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE get_elpa_row_col_comms(mpi_comm_global, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
! get_elpa_row_col_comms:
! All ELPA routines need MPI communicators for communicating within
! rows or columns of processes, these are set here.
! mpi_comm_rows/mpi_comm_cols can be free'd with MPI_Comm_free if not used any more.
!
!
!  mpi_comm_global   Global communicator for the calculations (in)
!
!  my_prow           Row coordinate of the calling process in the process grid (in)
!
!  my_pcol           Column coordinate of the calling process in the process grid (in)
!
!  mpi_comm_rows     Communicator for communicating within rows of processes (out)
!
!  mpi_comm_cols     Communicator for communicating within columns of processes (out)
!
!-------------------------------------------------------------------------------


    INTEGER, INTENT(in)                      :: mpi_comm_global, my_prow, &
                                                my_pcol
    INTEGER, INTENT(out)                     :: mpi_comm_rows, mpi_comm_cols

    INTEGER                                  :: mpierr

! mpi_comm_rows is used for communicating WITHIN rows, i.e. all processes
! having the same column coordinate share one mpi_comm_rows.
! So the "color" for splitting is my_pcol and the "key" is my row coordinate.
! Analogous for mpi_comm_cols

   CALL mpi_comm_split(mpi_comm_global,my_pcol,my_prow,mpi_comm_rows,mpierr)
   CALL mpi_comm_split(mpi_comm_global,my_prow,my_pcol,mpi_comm_cols,mpierr)

END SUBROUTINE get_elpa_row_col_comms

!-------------------------------------------------------------------------------

SUBROUTINE solve_evp_real(na, nev, a, lda, ev, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  solve_evp_real: Solves the real eigenvalue problem
!
!
!  na          Order of matrix a
!
!  nev         Number of eigenvalues needed.
!              The smallest nev eigenvalues/eigenvectors are calculated.
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
!
!-------------------------------------------------------------------------------


    INTEGER, INTENT(in)                      :: na, nev, lda
    REAL*8                                   :: a(lda,*), ev(na)
    INTEGER, INTENT(in)                      :: ldq
    REAL*8                                   :: q(ldq,*)
    INTEGER, INTENT(in)                      :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols

    INTEGER                                  :: mpierr, my_pcol, my_prow
    REAL*8                                   :: ttt0, ttt1
    REAL*8, ALLOCATABLE                      :: e(:), tau(:)

   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)

   ALLOCATE(e(na), tau(na))

   ttt0 = MPI_Wtime()
   CALL tridiag_real(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols, ev, e, tau)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) PRINT *,'Time tridiag_real :',ttt1-ttt0
   time_evp_fwd = ttt1-ttt0

   ttt0 = MPI_Wtime()
   CALL solve_tridi(na, nev, ev, e, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) PRINT *,'Time solve_tridi  :',ttt1-ttt0
   time_evp_solve = ttt1-ttt0

   ttt0 = MPI_Wtime()
   CALL trans_ev_real(na, nev, a, lda, tau, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) PRINT *,'Time trans_ev_real:',ttt1-ttt0
   time_evp_back = ttt1-ttt0

   DEALLOCATE(e, tau)

END SUBROUTINE solve_evp_real

!-------------------------------------------------------------------------------


SUBROUTINE solve_evp_complex(na, nev, a, lda, ev, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  solve_evp_complex: Solves the complex eigenvalue problem
!
!
!  na          Order of matrix a
!
!  nev         Number of eigenvalues needed
!              The smallest nev eigenvalues/eigenvectors are calculated.
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
!
!-------------------------------------------------------------------------------


    INTEGER, INTENT(in)                      :: na, nev, lda
    COMPLEX*16                               :: a(lda,*)
    REAL*8                                   :: ev(na)
    INTEGER, INTENT(in)                      :: ldq
    COMPLEX*16                               :: q(ldq,*)
    INTEGER, INTENT(in)                      :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols

    COMPLEX*16, ALLOCATABLE                  :: tau(:)
    INTEGER                                  :: l_cols, l_cols_nev, l_rows, &
                                                mpierr, my_pcol, my_prow, &
                                                np_cols, np_rows
    REAL*8                                   :: ttt0, ttt1
    REAL*8, ALLOCATABLE                      :: e(:), q_real(:,:)

   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q

   l_cols_nev = local_index(nev, my_pcol, np_cols, nblk, -1) ! Local columns corresponding to nev

   ALLOCATE(e(na), tau(na))
   ALLOCATE(q_real(l_rows,l_cols))

   ttt0 = MPI_Wtime()
   CALL tridiag_complex(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols, ev, e, tau)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) PRINT *,'Time tridiag_complex :',ttt1-ttt0
   time_evp_fwd = ttt1-ttt0

   ttt0 = MPI_Wtime()
   CALL solve_tridi(na, nev, ev, e, q_real, l_rows, nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) PRINT *,'Time solve_tridi     :',ttt1-ttt0
   time_evp_solve = ttt1-ttt0

   ttt0 = MPI_Wtime()
   q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)

   CALL trans_ev_complex(na, nev, a, lda, tau, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   IF(my_prow==0 .AND. my_pcol==0 .AND. elpa_print_times) PRINT *,'Time trans_ev_complex:',ttt1-ttt0
   time_evp_back = ttt1-ttt0

   DEALLOCATE(q_real)
   DEALLOCATE(e, tau)

END SUBROUTINE solve_evp_complex

!-------------------------------------------------------------------------------

SUBROUTINE tridiag_real(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols, d, e, tau)

!-------------------------------------------------------------------------------
!  tridiag_real: Reduces a distributed symmetric matrix to tridiagonal form
!                (like Scalapack Routine PDSYTRD)
!
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to PDSYTRD, a(:,:) must be set completely (upper and lower half)
!              a(:,:) is overwritten on exit with the Householder vectors
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  d(na)       Diagonal elements (returned), identical on all processors
!
!  e(na)       Off-Diagonal elements (returned), identical on all processors
!
!  tau(na)     Factors for the Householder vectors (returned), needed for back transformation
!
!-------------------------------------------------------------------------------


    INTEGER                                  :: na, lda
    REAL*8                                   :: a(lda,*)
    INTEGER                                  :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols
    REAL*8                                   :: d(na), e(na), tau(na)

    INTEGER, PARAMETER                       :: max_stored_rows = 32

    INTEGER :: i, istep, j, l_cols, l_cols_tile, l_rows, l_rows_tile, lce, &
      lcs, lre, lrs, max_blocks_col, max_blocks_row, max_local_cols, &
      max_local_rows, max_threads, mpierr, my_pcol, my_prow, my_thread, &
      n_iter, n_threads, np_cols, np_rows, nstor, tile_size, totalblocks

!$ integer omp_get_thread_num, omp_get_num_threads, omp_get_max_threads

   REAL*8 vav, vnorm2, x, aux(2*max_stored_rows), aux1(2), aux2(2), vrl, xf

   REAL*8, ALLOCATABLE:: tmp(:), vr(:), vc(:), ur(:), uc(:), vur(:,:), uvc(:,:)
   REAL*8, ALLOCATABLE:: ur_p(:,:), uc_p(:,:)

   INTEGER pcol, prow
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

   tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
   tile_size = ((128*MAX(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

   l_rows_tile = tile_size/np_rows ! local rows of a tile
   l_cols_tile = tile_size/np_cols ! local cols of a tile


   totalblocks = (na-1)/nblk + 1
   max_blocks_row = (totalblocks-1)/np_rows + 1
   max_blocks_col = (totalblocks-1)/np_cols + 1

   max_local_rows = max_blocks_row*nblk
   max_local_cols = max_blocks_col*nblk

   ALLOCATE(tmp(MAX(max_local_rows,max_local_cols)))
   ALLOCATE(vr(max_local_rows+1))
   ALLOCATE(ur(max_local_rows))
   ALLOCATE(vc(max_local_cols))
   ALLOCATE(uc(max_local_cols))

   max_threads = 1
!$ max_threads = omp_get_max_threads()
   ALLOCATE(ur_p(max_local_rows,0:max_threads-1))
   ALLOCATE(uc_p(max_local_cols,0:max_threads-1))

   tmp = 0
   vr = 0
   ur = 0
   vc = 0
   uc = 0

   ALLOCATE(vur(max_local_rows,2*max_stored_rows))
   ALLOCATE(uvc(max_local_cols,2*max_stored_rows))

   d(:) = 0
   e(:) = 0
   tau(:) = 0

   nstor = 0

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a
   IF(my_prow==prow(na) .AND. my_pcol==pcol(na)) d(na) = a(l_rows,l_cols)

   DO istep=na,3,-1

      ! Calculate number of local rows and columns of the still remaining matrix
      ! on the local processor

      l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
      l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

      ! Calculate vector for Householder transformation on all procs
      ! owning column istep

      IF(my_pcol==pcol(istep)) THEN

         ! Get vector to be transformed; distribute last element and norm of
         ! remaining elements to all procs in current column

         vr(1:l_rows) = a(1:l_rows,l_cols+1)
         IF(nstor>0 .AND. l_rows>0) THEN
            CALL DGEMV('N',l_rows,2*nstor,1.d0,vur,UBOUND(vur,1), &
                       uvc(l_cols+1,1),UBOUND(uvc,1),1.d0,vr,1)
         ENDIF

         IF(my_prow==prow(istep-1)) THEN
            aux1(1) = DOT_PRODUCT(vr(1:l_rows-1),vr(1:l_rows-1))
            aux1(2) = vr(l_rows)
         ELSE
            aux1(1) = DOT_PRODUCT(vr(1:l_rows),vr(1:l_rows))
            aux1(2) = 0.
         ENDIF

         CALL mpi_allreduce(aux1,aux2,2,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)

         vnorm2 = aux2(1)
         vrl    = aux2(2)

         ! Householder transformation

         CALL hh_transform_real(vrl, vnorm2, xf, tau(istep))

         ! Scale vr and store Householder vector for back transformation

         vr(1:l_rows) = vr(1:l_rows) * xf
         IF(my_prow==prow(istep-1)) THEN
            vr(l_rows) = 1.
            e(istep-1) = vrl
         ENDIF
         a(1:l_rows,l_cols+1) = vr(1:l_rows) ! store Householder vector for back transformation

      ENDIF

      ! Broadcast the Householder vector (and tau) along columns

      IF(my_pcol==pcol(istep)) vr(l_rows+1) = tau(istep)
      CALL MPI_Bcast(vr,l_rows+1,MPI_REAL8,pcol(istep),mpi_comm_cols,mpierr)
      tau(istep) =  vr(l_rows+1)

      ! Transpose Householder vector vr -> vc

      CALL elpa_transpose_vectors  (vr, UBOUND(vr,1), mpi_comm_rows, &
                                    vc, UBOUND(vc,1), mpi_comm_cols, &
                                    1, istep-1, 1, nblk)


      ! Calculate u = (A + VU**T + UV**T)*v

      ! For cache efficiency, we use only the upper half of the matrix tiles for this,
      ! thus the result is partly in uc(:) and partly in ur(:)

      uc(1:l_cols) = 0
      ur(1:l_rows) = 0
      IF(l_rows>0 .AND. l_cols>0) THEN

!$OMP PARALLEL PRIVATE(my_thread,n_threads,n_iter,i,lcs,lce,j,lrs,lre)
         my_thread = 0
         n_threads = 1
!$       my_thread = omp_get_thread_num()
!$       n_threads = omp_get_num_threads()
         n_iter = 0

         uc_p(1:l_cols,my_thread) = 0.
         ur_p(1:l_rows,my_thread) = 0.

         DO i=0,(istep-2)/tile_size
            lcs = i*l_cols_tile+1
            lce = MIN(l_cols,(i+1)*l_cols_tile)
            IF(lce<lcs) CYCLE
            DO j=0,i
               lrs = j*l_rows_tile+1
               lre = MIN(l_rows,(j+1)*l_rows_tile)
               IF(lre<lrs) CYCLE
               IF(MOD(n_iter,n_threads) == my_thread) THEN
                 CALL DGEMV('T',lre-lrs+1,lce-lcs+1,1.d0,a(lrs,lcs),lda,vr(lrs),1,1.d0,uc_p(lcs,my_thread),1)
                 IF(i/=j) CALL DGEMV('N',lre-lrs+1,lce-lcs+1,1.d0,a(lrs,lcs),lda,vc(lcs),1,1.d0,ur_p(lrs,my_thread),1)
               ENDIF
               n_iter = n_iter+1
            ENDDO
         ENDDO
!$OMP END PARALLEL

         DO i=0,max_threads-1
            uc(1:l_cols) = uc(1:l_cols) + uc_p(1:l_cols,i)
            ur(1:l_rows) = ur(1:l_rows) + ur_p(1:l_rows,i)
         ENDDO

         IF(nstor>0) THEN
            CALL DGEMV('T',l_rows,2*nstor,1.d0,vur,UBOUND(vur,1),vr,1,0.d0,aux,1)
            CALL DGEMV('N',l_cols,2*nstor,1.d0,uvc,UBOUND(uvc,1),aux,1,1.d0,uc,1)
         ENDIF

      ENDIF

      ! Sum up all ur(:) parts along rows and add them to the uc(:) parts
      ! on the processors containing the diagonal
      ! This is only necessary if ur has been calculated, i.e. if the
      ! global tile size is smaller than the global remaining matrix

      IF(tile_size < istep-1) THEN
         CALL elpa_reduce_add_vectors  (ur, UBOUND(ur,1), mpi_comm_rows, &
                                        uc, UBOUND(uc,1), mpi_comm_cols, &
                                        istep-1, 1, nblk)
      ENDIF

      ! Sum up all the uc(:) parts, transpose uc -> ur

      IF(l_cols>0) THEN
         tmp(1:l_cols) = uc(1:l_cols)
         CALL mpi_allreduce(tmp,uc,l_cols,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
      ENDIF

      CALL elpa_transpose_vectors  (uc, UBOUND(uc,1), mpi_comm_cols, &
                                    ur, UBOUND(ur,1), mpi_comm_rows, &
                                    1, istep-1, 1, nblk)

      ! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )

      x = 0
      IF(l_cols>0) x = DOT_PRODUCT(vc(1:l_cols),uc(1:l_cols))
      CALL mpi_allreduce(x,vav,1,MPI_REAL8,MPI_SUM,mpi_comm_cols,mpierr)

      ! store u and v in the matrices U and V
      ! these matrices are stored combined in one here

      DO j=1,l_rows
         vur(j,2*nstor+1) = tau(istep)*vr(j)
         vur(j,2*nstor+2) = 0.5*tau(istep)*vav*vr(j) - ur(j)
      ENDDO
      DO j=1,l_cols
         uvc(j,2*nstor+1) = 0.5*tau(istep)*vav*vc(j) - uc(j)
         uvc(j,2*nstor+2) = tau(istep)*vc(j)
      ENDDO

      nstor = nstor+1

      ! If the limit of max_stored_rows is reached, calculate A + VU**T + UV**T

      IF(nstor==max_stored_rows .OR. istep==3) THEN

         DO i=0,(istep-2)/tile_size
            lcs = i*l_cols_tile+1
            lce = MIN(l_cols,(i+1)*l_cols_tile)
            lrs = 1
            lre = MIN(l_rows,(i+1)*l_rows_tile)
            IF(lce<lcs .OR. lre<lrs) CYCLE
            CALL dgemm('N','T',lre-lrs+1,lce-lcs+1,2*nstor,1.d0, &
                       vur(lrs,1),UBOUND(vur,1),uvc(lcs,1),UBOUND(uvc,1), &
                       1.d0,a(lrs,lcs),lda)
         ENDDO

         nstor = 0

      ENDIF

      IF(my_prow==prow(istep-1) .AND. my_pcol==pcol(istep-1)) THEN
         IF(nstor>0) a(l_rows,l_cols) = a(l_rows,l_cols) &
                        + DOT_PRODUCT(vur(l_rows,1:2*nstor),uvc(l_cols,1:2*nstor))
         d(istep-1) = a(l_rows,l_cols)
      ENDIF

   ENDDO

   ! Store e(1) and d(1)

   IF(my_prow==prow(1) .AND. my_pcol==pcol(2)) e(1) = a(1,l_cols) ! use last l_cols value of loop above
   IF(my_prow==prow(1) .AND. my_pcol==pcol(1)) d(1) = a(1,1)

   DEALLOCATE(tmp, vr, ur, vc, uc, vur, uvc)

   ! distribute the arrays d and e to all processors

   ALLOCATE(tmp(na))
   tmp = d
   CALL mpi_allreduce(tmp,d,na,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
   tmp = d
   CALL mpi_allreduce(tmp,d,na,MPI_REAL8,MPI_SUM,mpi_comm_cols,mpierr)
   tmp = e
   CALL mpi_allreduce(tmp,e,na,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
   tmp = e
   CALL mpi_allreduce(tmp,e,na,MPI_REAL8,MPI_SUM,mpi_comm_cols,mpierr)
   DEALLOCATE(tmp)

END SUBROUTINE tridiag_real

!-------------------------------------------------------------------------------

SUBROUTINE trans_ev_real(na, nqc, a, lda, tau, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  trans_ev_real: Transforms the eigenvectors of a tridiagonal matrix back
!                 to the eigenvectors of the original matrix
!                 (like Scalapack Routine PDORMTR)
!
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nqc         Number of columns of matrix q
!
!  a(lda,*)    Matrix containing the Householder vectors (i.e. matrix a after tridiag_real)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a
!
!  tau(na)     Factors of the Householder vectors
!
!  q           On input: Eigenvectors of tridiagonal matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------


    INTEGER                                  :: na, nqc, lda
    REAL*8                                   :: a(lda,*), tau(na)
    INTEGER                                  :: ldq
    REAL*8                                   :: q(ldq,*)
    INTEGER                                  :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols

    INTEGER :: cur_pcol, i, ic, ice, ics, istep, l_colh, l_cols, l_rows, &
      max_blocks_col, max_blocks_row, max_local_cols, max_local_rows, &
      max_stored_rows, mpierr, my_pcol, my_prow, n, nb, nc, np_cols, np_rows, &
      nstor, pcol, prow, totalblocks
    REAL*8, ALLOCATABLE                      :: h1(:), h2(:), hvb(:), &
                                                hvm(:,:), tmat(:,:), tmp1(:), &
                                                tmp2(:)

   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)


   totalblocks = (na-1)/nblk + 1
   max_blocks_row = (totalblocks-1)/np_rows + 1
   max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q!

   max_local_rows = max_blocks_row*nblk
   max_local_cols = max_blocks_col*nblk


   max_stored_rows = (63/nblk+1)*nblk

   ALLOCATE(tmat(max_stored_rows,max_stored_rows))
   ALLOCATE(h1(max_stored_rows*max_stored_rows))
   ALLOCATE(h2(max_stored_rows*max_stored_rows))
   ALLOCATE(tmp1(max_local_cols*max_stored_rows))
   ALLOCATE(tmp2(max_local_cols*max_stored_rows))
   ALLOCATE(hvb(max_local_rows*nblk))
   ALLOCATE(hvm(max_local_rows,max_stored_rows))

   hvm = 0   ! Must be set to 0 !!!
   hvb = 0   ! Safety only

   l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q

   nstor = 0

   DO istep=1,na,nblk

      ics = MAX(istep,3)
      ice = MIN(istep+nblk-1,na)
      IF(ice<ics) CYCLE

      cur_pcol = pcol(istep)

      nb = 0
      DO ic=ics,ice

         l_colh = local_index(ic  , my_pcol, np_cols, nblk, -1) ! Column of Householder vector
         l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder vector


         IF(my_pcol==cur_pcol) THEN
            hvb(nb+1:nb+l_rows) = a(1:l_rows,l_colh)
            IF(my_prow==prow(ic-1)) THEN
               hvb(nb+l_rows) = 1.
            ENDIF
         ENDIF

         nb = nb+l_rows
      ENDDO

      IF(nb>0) &
         CALL MPI_Bcast(hvb,nb,MPI_REAL8,cur_pcol,mpi_comm_cols,mpierr)

      nb = 0
      DO ic=ics,ice
         l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder vector
         hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
         nstor = nstor+1
         nb = nb+l_rows
      ENDDO

      ! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
      IF(nstor+nblk>max_stored_rows .OR. istep+nblk>na .OR. (na/np_rows<=256 .AND. nstor>=32)) THEN

         ! Calculate scalar products of stored vectors.
         ! This can be done in different ways, we use dsyrk

         tmat = 0
         IF(l_rows>0) &
            CALL dsyrk('U','T',nstor,l_rows,1.d0,hvm,UBOUND(hvm,1),0.d0,tmat,max_stored_rows)

         nc = 0
         DO n=1,nstor-1
            h1(nc+1:nc+n) = tmat(1:n,n+1)
            nc = nc+n
         ENDDO

         IF(nc>0) CALL mpi_allreduce(h1,h2,nc,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)

         ! Calculate triangular matrix T

         nc = 0
         tmat(1,1) = tau(ice-nstor+1)
         DO n=1,nstor-1
            CALL dtrmv('L','T','N',n,tmat,max_stored_rows,h2(nc+1),1)
            tmat(n+1,1:n) = -h2(nc+1:nc+n)*tau(ice-nstor+n+1)
            tmat(n+1,n+1) = tau(ice-nstor+n+1)
            nc = nc+n
         ENDDO

         ! Q = Q - V * T * V**T * Q

         IF(l_rows>0) THEN
            CALL dgemm('T','N',nstor,l_cols,l_rows,1.d0,hvm,UBOUND(hvm,1), &
                       q,ldq,0.d0,tmp1,nstor)
         ELSE
            tmp1(1:l_cols*nstor) = 0
         ENDIF
         CALL mpi_allreduce(tmp1,tmp2,nstor*l_cols,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
         IF(l_rows>0) THEN
            CALL dtrmm('L','L','N','N',nstor,l_cols,1.0d0,tmat,max_stored_rows,tmp2,nstor)
            CALL dgemm('N','N',l_rows,l_cols,nstor,-1.d0,hvm,UBOUND(hvm,1), &
                       tmp2,nstor,1.d0,q,ldq)
         ENDIF
         nstor = 0
      ENDIF

   ENDDO

   DEALLOCATE(tmat, h1, h2, tmp1, tmp2, hvb, hvm)


END SUBROUTINE trans_ev_real

!-------------------------------------------------------------------------------

SUBROUTINE mult_at_b_real(uplo_a, uplo_c, na, ncb, a, lda, b, ldb, nblk, mpi_comm_rows, mpi_comm_cols, c, ldc)

!-------------------------------------------------------------------------------
!  mult_at_b_real:  Performs C := A**T * B
!
!      where:  A is a square matrix (na,na) which is optionally upper or lower triangular
!              B is a (na,ncb) matrix
!              C is a (na,ncb) matrix where optionally only the upper or lower
!              triangle may be computed
!
!
!  uplo_a      'U' if A is upper triangular
!              'L' if A is lower triangular
!              anything else if A is a full matrix
!              Please note: This pertains to the original A (as set in the calling program)
!              whereas the transpose of A is used for calculations
!              If uplo_a is 'U' or 'L', the other triangle is not used at all,
!              i.e. it may contain arbitrary numbers
!
!  uplo_c      'U' if only the upper diagonal part of C is needed
!              'L' if only the upper diagonal part of C is needed
!              anything else if the full matrix C is needed
!              Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
!              written to a certain extent, i.e. one shouldn't rely on the content there!
!
!  na          Number of rows/columns of A, number of rows of B and C
!
!  ncb         Number of columns  of B and C
!
!  a           Matrix A
!
!  lda         Leading dimension of a
!
!  b           Matrix B
!
!  ldb         Leading dimension of b
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  c           Matrix C
!
!  ldc         Leading dimension of c
!
!-------------------------------------------------------------------------------


    CHARACTER(len=1)                         :: uplo_a, uplo_c
    INTEGER                                  :: na, ncb, lda
    REAL*8                                   :: a(lda,*)
    INTEGER                                  :: ldb
    REAL*8                                   :: b(ldb,*)
    INTEGER                                  :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols, ldc
    REAL*8                                   :: c(ldc,*)

    INTEGER :: gcol, gcol_min, goff, l_cols, l_rows, l_rows_np, lce, lcs, &
      lre, lrs, mpierr, my_pcol, my_prow, n, n_aux_bc, nb, nblk_mult, noff, &
      np, np_bc, np_cols, np_rows, nr_done, nstor, nvals
    INTEGER, ALLOCATABLE                     :: lre_save(:), lrs_save(:)
    LOGICAL                                  :: a_lower, a_upper, c_lower, &
                                                c_upper
    REAL*8, ALLOCATABLE                      :: aux_bc(:), aux_mat(:,:), &
                                                tmp1(:,:), tmp2(:,:)

   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   l_rows = local_index(na,  my_prow, np_rows, nblk, -1) ! Local rows of a and b
   l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b

   ! Block factor for matrix multiplications, must be a multiple of nblk

   IF(na/np_rows<=256) THEN
      nblk_mult = (31/nblk+1)*nblk
   ELSE
      nblk_mult = (63/nblk+1)*nblk
   ENDIF

   ALLOCATE(aux_mat(l_rows,nblk_mult))
   ALLOCATE(aux_bc(l_rows*nblk))
   ALLOCATE(lrs_save(nblk))
   ALLOCATE(lre_save(nblk))

   a_lower = .FALSE.
   a_upper = .FALSE.
   c_lower = .FALSE.
   c_upper = .FALSE.

   IF(uplo_a=='u' .OR. uplo_a=='U') a_upper = .TRUE.
   IF(uplo_a=='l' .OR. uplo_a=='L') a_lower = .TRUE.
   IF(uplo_c=='u' .OR. uplo_c=='U') c_upper = .TRUE.
   IF(uplo_c=='l' .OR. uplo_c=='L') c_lower = .TRUE.

   ! Build up the result matrix by processor rows

   DO np = 0, np_rows-1

      ! In this turn, procs of row np assemble the result

      l_rows_np = local_index(na, np, np_rows, nblk, -1) ! local rows on receiving processors

      nr_done = 0 ! Number of rows done
      aux_mat = 0
      nstor = 0   ! Number of columns stored in aux_mat

      ! Loop over the blocks on row np

      DO nb=0,(l_rows_np-1)/nblk

         goff  = nb*np_rows + np ! Global offset in blocks corresponding to nb

         ! Get the processor column which owns this block (A is transposed, so we need the column)
         ! and the offset in blocks within this column.
         ! The corresponding block column in A is then broadcast to all for multiplication with B

         np_bc = MOD(goff,np_cols)
         noff = goff/np_cols
         n_aux_bc = 0

         ! Gather up the complete block column of A on the owner

         DO n = 1, MIN(l_rows_np-nb*nblk,nblk) ! Loop over columns to be broadcast

            gcol = goff*nblk + n ! global column corresponding to n
            IF(nstor==0 .AND. n==1) gcol_min = gcol

            lrs = 1       ! 1st local row number for broadcast
            lre = l_rows  ! last local row number for broadcast
            IF(a_lower) lrs = local_index(gcol, my_prow, np_rows, nblk, +1)
            IF(a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

            IF(lrs<=lre) THEN
               nvals = lre-lrs+1
               IF(my_pcol == np_bc) aux_bc(n_aux_bc+1:n_aux_bc+nvals) = a(lrs:lre,noff*nblk+n)
               n_aux_bc = n_aux_bc + nvals
            ENDIF

            lrs_save(n) = lrs
            lre_save(n) = lre

         ENDDO

         ! Broadcast block column

         CALL MPI_Bcast(aux_bc,n_aux_bc,MPI_REAL8,np_bc,mpi_comm_cols,mpierr)

         ! Insert what we got in aux_mat

         n_aux_bc = 0
         DO n = 1, MIN(l_rows_np-nb*nblk,nblk)
            nstor = nstor+1
            lrs = lrs_save(n)
            lre = lre_save(n)
            IF(lrs<=lre) THEN
               nvals = lre-lrs+1
               aux_mat(lrs:lre,nstor) = aux_bc(n_aux_bc+1:n_aux_bc+nvals)
               n_aux_bc = n_aux_bc + nvals
            ENDIF
         ENDDO

         ! If we got nblk_mult columns in aux_mat or this is the last block
         ! do the matrix multiplication

         IF(nstor==nblk_mult .OR. nb*nblk+nblk >= l_rows_np) THEN

            lrs = 1       ! 1st local row number for multiply
            lre = l_rows  ! last local row number for multiply
            IF(a_lower) lrs = local_index(gcol_min, my_prow, np_rows, nblk, +1)
            IF(a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

            lcs = 1       ! 1st local col number for multiply
            lce = l_cols  ! last local col number for multiply
            IF(c_upper) lcs = local_index(gcol_min, my_pcol, np_cols, nblk, +1)
            IF(c_lower) lce = MIN(local_index(gcol, my_pcol, np_cols, nblk, -1),l_cols)

            IF(lcs<=lce) THEN
               ALLOCATE(tmp1(nstor,lcs:lce),tmp2(nstor,lcs:lce))
               IF(lrs<=lre) THEN
                  CALL dgemm('T','N',nstor,lce-lcs+1,lre-lrs+1,1.d0,aux_mat(lrs,1),UBOUND(aux_mat,1), &
                             b(lrs,lcs),ldb,0.d0,tmp1,nstor)
               ELSE
                  tmp1 = 0
               ENDIF

               ! Sum up the results and send to processor row np
               CALL mpi_reduce(tmp1,tmp2,nstor*(lce-lcs+1),MPI_REAL8,MPI_SUM,np,mpi_comm_rows,mpierr)

               ! Put the result into C
               IF(my_prow==np) c(nr_done+1:nr_done+nstor,lcs:lce) = tmp2(1:nstor,lcs:lce)

               DEALLOCATE(tmp1,tmp2)
            ENDIF

            nr_done = nr_done+nstor
            nstor=0
            aux_mat(:,:)=0
         ENDIF
      ENDDO
   ENDDO

   DEALLOCATE(aux_mat, aux_bc, lrs_save, lre_save)

END SUBROUTINE mult_at_b_real

!-------------------------------------------------------------------------------

SUBROUTINE tridiag_complex(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols, d, e, tau)

!-------------------------------------------------------------------------------
!  tridiag_complex: Reduces a distributed hermitian matrix to tridiagonal form
!                   (like Scalapack Routine PZHETRD)
!
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to PZHETRD, a(:,:) must be set completely (upper and lower half)
!              a(:,:) is overwritten on exit with the Householder vectors
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  d(na)       Diagonal elements (returned), identical on all processors
!
!  e(na)       Off-Diagonal elements (returned), identical on all processors
!
!  tau(na)     Factors for the Householder vectors (returned), needed for back transformation
!
!-------------------------------------------------------------------------------


    INTEGER                                  :: na, lda
    COMPLEX*16                               :: a(lda,*)
    INTEGER                                  :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols
    REAL*8                                   :: d(na), e(na)
    COMPLEX*16                               :: tau(na)

    COMPLEX*16, PARAMETER                    :: CONE = (1.d0,0.d0), &
                                                CZERO = (0.d0,0.d0)
    INTEGER, PARAMETER                       :: max_stored_rows = 32

    INTEGER :: i, istep, j, l_cols, l_cols_tile, l_rows, l_rows_tile, lce, &
      lcs, lre, lrs, max_blocks_col, max_blocks_row, max_local_cols, &
      max_local_rows, max_threads, mpierr, my_pcol, my_prow, my_thread, &
      n_iter, n_threads, np_cols, np_rows, nstor, tile_size, totalblocks

!$ integer omp_get_thread_num, omp_get_num_threads, omp_get_max_threads

   REAL*8 vnorm2
   COMPLEX*16 vav, xc, aux(2*max_stored_rows),  aux1(2), aux2(2), vrl, xf

   COMPLEX*16, ALLOCATABLE:: tmp(:), vr(:), vc(:), ur(:), uc(:), vur(:,:), uvc(:,:)
   COMPLEX*16, ALLOCATABLE:: ur_p(:,:), uc_p(:,:)
   REAL*8, ALLOCATABLE:: tmpr(:)

   INTEGER pcol, prow
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

   tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
   tile_size = ((128*MAX(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

   l_rows_tile = tile_size/np_rows ! local rows of a tile
   l_cols_tile = tile_size/np_cols ! local cols of a tile


   totalblocks = (na-1)/nblk + 1
   max_blocks_row = (totalblocks-1)/np_rows + 1
   max_blocks_col = (totalblocks-1)/np_cols + 1

   max_local_rows = max_blocks_row*nblk
   max_local_cols = max_blocks_col*nblk

   ALLOCATE(tmp(MAX(max_local_rows,max_local_cols)))
   ALLOCATE(vr(max_local_rows+1))
   ALLOCATE(ur(max_local_rows))
   ALLOCATE(vc(max_local_cols))
   ALLOCATE(uc(max_local_cols))

   max_threads = 1
!$ max_threads = omp_get_max_threads()
   ALLOCATE(ur_p(max_local_rows,0:max_threads-1))
   ALLOCATE(uc_p(max_local_cols,0:max_threads-1))

   tmp = 0
   vr = 0
   ur = 0
   vc = 0
   uc = 0

   ALLOCATE(vur(max_local_rows,2*max_stored_rows))
   ALLOCATE(uvc(max_local_cols,2*max_stored_rows))

   d(:) = 0
   e(:) = 0
   tau(:) = 0

   nstor = 0

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a
   IF(my_prow==prow(na) .AND. my_pcol==pcol(na)) d(na) = a(l_rows,l_cols)

   DO istep=na,3,-1

      ! Calculate number of local rows and columns of the still remaining matrix
      ! on the local processor

      l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
      l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

      ! Calculate vector for Householder transformation on all procs
      ! owning column istep

      IF(my_pcol==pcol(istep)) THEN

         ! Get vector to be transformed; distribute last element and norm of
         ! remaining elements to all procs in current column

         vr(1:l_rows) = a(1:l_rows,l_cols+1)
         IF(nstor>0 .AND. l_rows>0) THEN
            aux(1:2*nstor) = CONJG(uvc(l_cols+1,1:2*nstor))
            CALL ZGEMV('N',l_rows,2*nstor,CONE,vur,UBOUND(vur,1), &
                       aux,1,CONE,vr,1)
         ENDIF

         IF(my_prow==prow(istep-1)) THEN
            aux1(1) = DOT_PRODUCT(vr(1:l_rows-1),vr(1:l_rows-1))
            aux1(2) = vr(l_rows)
         ELSE
            aux1(1) = DOT_PRODUCT(vr(1:l_rows),vr(1:l_rows))
            aux1(2) = 0.
         ENDIF

         CALL mpi_allreduce(aux1,aux2,2,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)

         vnorm2 = aux2(1)
         vrl    = aux2(2)

         ! Householder transformation

         CALL hh_transform_complex(vrl, vnorm2, xf, tau(istep))

         ! Scale vr and store Householder vector for back transformation

         vr(1:l_rows) = vr(1:l_rows) * xf
         IF(my_prow==prow(istep-1)) THEN
            vr(l_rows) = 1.
            e(istep-1) = vrl
         ENDIF
         a(1:l_rows,l_cols+1) = vr(1:l_rows) ! store Householder vector for back transformation

      ENDIF

      ! Broadcast the Householder vector (and tau) along columns

      IF(my_pcol==pcol(istep)) vr(l_rows+1) = tau(istep)
      CALL MPI_Bcast(vr,l_rows+1,MPI_DOUBLE_COMPLEX,pcol(istep),mpi_comm_cols,mpierr)
      tau(istep) =  vr(l_rows+1)

      ! Transpose Householder vector vr -> vc

      CALL elpa_transpose_vectors  (vr, 2*UBOUND(vr,1), mpi_comm_rows, &
                                    vc, 2*UBOUND(vc,1), mpi_comm_cols, &
                                    1, 2*(istep-1), 1, 2*nblk)

      ! Calculate u = (A + VU**T + UV**T)*v

      ! For cache efficiency, we use only the upper half of the matrix tiles for this,
      ! thus the result is partly in uc(:) and partly in ur(:)

      uc(1:l_cols) = 0
      ur(1:l_rows) = 0
      IF(l_rows>0 .AND. l_cols>0) THEN

!$OMP PARALLEL PRIVATE(my_thread,n_threads,n_iter,i,lcs,lce,j,lrs,lre)
         my_thread = 0
         n_threads = 1
!$       my_thread = omp_get_thread_num()
!$       n_threads = omp_get_num_threads()
         n_iter = 0

         uc_p(1:l_cols,my_thread) = 0.
         ur_p(1:l_rows,my_thread) = 0.

         DO i=0,(istep-2)/tile_size
            lcs = i*l_cols_tile+1
            lce = MIN(l_cols,(i+1)*l_cols_tile)
            IF(lce<lcs) CYCLE
            DO j=0,i
               lrs = j*l_rows_tile+1
               lre = MIN(l_rows,(j+1)*l_rows_tile)
               IF(lre<lrs) CYCLE
               IF(MOD(n_iter,n_threads) == my_thread) THEN
                  CALL ZGEMV('C',lre-lrs+1,lce-lcs+1,CONE,a(lrs,lcs),lda,vr(lrs),1,CONE,uc_p(lcs,my_thread),1)
                  IF(i/=j) CALL ZGEMV('N',lre-lrs+1,lce-lcs+1,CONE,a(lrs,lcs),lda,vc(lcs),1,CONE,ur_p(lrs,my_thread),1)
               ENDIF
               n_iter = n_iter+1
            ENDDO
         ENDDO
!$OMP END PARALLEL

         DO i=0,max_threads-1
            uc(1:l_cols) = uc(1:l_cols) + uc_p(1:l_cols,i)
            ur(1:l_rows) = ur(1:l_rows) + ur_p(1:l_rows,i)
         ENDDO

         IF(nstor>0) THEN
            CALL ZGEMV('C',l_rows,2*nstor,CONE,vur,UBOUND(vur,1),vr,1,CZERO,aux,1)
            CALL ZGEMV('N',l_cols,2*nstor,CONE,uvc,UBOUND(uvc,1),aux,1,CONE,uc,1)
         ENDIF

      ENDIF

      ! Sum up all ur(:) parts along rows and add them to the uc(:) parts
      ! on the processors containing the diagonal
      ! This is only necessary if ur has been calculated, i.e. if the
      ! global tile size is smaller than the global remaining matrix

      IF(tile_size < istep-1) THEN
         CALL elpa_reduce_add_vectors  (ur, 2*UBOUND(ur,1), mpi_comm_rows, &
                                        uc, 2*UBOUND(uc,1), mpi_comm_cols, &
                                        2*(istep-1), 1, 2*nblk)
      ENDIF

      ! Sum up all the uc(:) parts, transpose uc -> ur

      IF(l_cols>0) THEN
         tmp(1:l_cols) = uc(1:l_cols)
         CALL mpi_allreduce(tmp,uc,l_cols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)
      ENDIF

      CALL elpa_transpose_vectors  (uc, 2*UBOUND(uc,1), mpi_comm_cols, &
                                    ur, 2*UBOUND(ur,1), mpi_comm_rows, &
                                    1, 2*(istep-1), 1, 2*nblk)

      ! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )

      xc = 0
      IF(l_cols>0) xc = DOT_PRODUCT(vc(1:l_cols),uc(1:l_cols))
      CALL mpi_allreduce(xc,vav,1,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_cols,mpierr)

      ! store u and v in the matrices U and V
      ! these matrices are stored combined in one here

      DO j=1,l_rows
         vur(j,2*nstor+1) = CONJG(tau(istep))*vr(j)
         vur(j,2*nstor+2) = 0.5*CONJG(tau(istep))*vav*vr(j) - ur(j)
      ENDDO
      DO j=1,l_cols
         uvc(j,2*nstor+1) = 0.5*CONJG(tau(istep))*vav*vc(j) - uc(j)
         uvc(j,2*nstor+2) = CONJG(tau(istep))*vc(j)
      ENDDO

      nstor = nstor+1

      ! If the limit of max_stored_rows is reached, calculate A + VU**T + UV**T

      IF(nstor==max_stored_rows .OR. istep==3) THEN

         DO i=0,(istep-2)/tile_size
            lcs = i*l_cols_tile+1
            lce = MIN(l_cols,(i+1)*l_cols_tile)
            lrs = 1
            lre = MIN(l_rows,(i+1)*l_rows_tile)
            IF(lce<lcs .OR. lre<lrs) CYCLE
            CALL ZGEMM('N','C',lre-lrs+1,lce-lcs+1,2*nstor,CONE, &
                       vur(lrs,1),UBOUND(vur,1),uvc(lcs,1),UBOUND(uvc,1), &
                       CONE,a(lrs,lcs),lda)
         ENDDO

         nstor = 0

      ENDIF

      IF(my_prow==prow(istep-1) .AND. my_pcol==pcol(istep-1)) THEN
         IF(nstor>0) a(l_rows,l_cols) = a(l_rows,l_cols) &
                        + DOT_PRODUCT(vur(l_rows,1:2*nstor),uvc(l_cols,1:2*nstor))
         d(istep-1) = a(l_rows,l_cols)
      ENDIF

   ENDDO

   ! Store e(1) and d(1)

   IF(my_pcol==pcol(2)) THEN
      IF(my_prow==prow(1)) THEN
         ! We use last l_cols value of loop above
         vrl = a(1,l_cols)
         CALL hh_transform_complex(vrl, 0.d0, xf, tau(2))
         e(1) = vrl
         a(1,l_cols) = 1. ! for consistency only
      ENDIF
      CALL mpi_bcast(tau(2),1,MPI_DOUBLE_COMPLEX,prow(1),mpi_comm_rows,mpierr)
   ENDIF
   CALL mpi_bcast(tau(2),1,MPI_DOUBLE_COMPLEX,pcol(2),mpi_comm_cols,mpierr)

   IF(my_prow==prow(1) .AND. my_pcol==pcol(1)) d(1) = a(1,1)

   DEALLOCATE(tmp, vr, ur, vc, uc, vur, uvc)

   ! distribute the arrays d and e to all processors

   ALLOCATE(tmpr(na))
   tmpr = d
   CALL mpi_allreduce(tmpr,d,na,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
   tmpr = d
   CALL mpi_allreduce(tmpr,d,na,MPI_REAL8,MPI_SUM,mpi_comm_cols,mpierr)
   tmpr = e
   CALL mpi_allreduce(tmpr,e,na,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
   tmpr = e
   CALL mpi_allreduce(tmpr,e,na,MPI_REAL8,MPI_SUM,mpi_comm_cols,mpierr)
   DEALLOCATE(tmpr)

END SUBROUTINE tridiag_complex

!-------------------------------------------------------------------------------

SUBROUTINE trans_ev_complex(na, nqc, a, lda, tau, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  trans_ev_complex: Transforms the eigenvectors of a tridiagonal matrix back
!                    to the eigenvectors of the original matrix
!                    (like Scalapack Routine PZUNMTR)
!
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nqc         Number of columns of matrix q
!
!  a(lda,*)    Matrix containing the Householder vectors (i.e. matrix a after tridiag_complex)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a
!
!  tau(na)     Factors of the Householder vectors
!
!  q           On input: Eigenvectors of tridiagonal matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------


    INTEGER                                  :: na, nqc, lda
    COMPLEX*16                               :: a(lda,*), tau(na)
    INTEGER                                  :: ldq
    COMPLEX*16                               :: q(ldq,*)
    INTEGER                                  :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols

    COMPLEX*16, PARAMETER                    :: CONE = (1.d0,0.d0), &
                                                CZERO = (0.d0,0.d0)

    COMPLEX*16, ALLOCATABLE                  :: h1(:), h2(:), hvb(:), &
                                                hvm(:,:), tmat(:,:), tmp1(:), &
                                                tmp2(:)
    INTEGER :: cur_pcol, i, ic, ice, ics, istep, l_colh, l_cols, l_rows, &
      max_blocks_col, max_blocks_row, max_local_cols, max_local_rows, &
      max_stored_rows, mpierr, my_pcol, my_prow, n, nb, nc, np_cols, np_rows, &
      nstor, pcol, prow, totalblocks

   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)


   totalblocks = (na-1)/nblk + 1
   max_blocks_row = (totalblocks-1)/np_rows + 1
   max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q!

   max_local_rows = max_blocks_row*nblk
   max_local_cols = max_blocks_col*nblk


   max_stored_rows = (63/nblk+1)*nblk

   ALLOCATE(tmat(max_stored_rows,max_stored_rows))
   ALLOCATE(h1(max_stored_rows*max_stored_rows))
   ALLOCATE(h2(max_stored_rows*max_stored_rows))
   ALLOCATE(tmp1(max_local_cols*max_stored_rows))
   ALLOCATE(tmp2(max_local_cols*max_stored_rows))
   ALLOCATE(hvb(max_local_rows*nblk))
   ALLOCATE(hvm(max_local_rows,max_stored_rows))

   hvm = 0   ! Must be set to 0 !!!
   hvb = 0   ! Safety only

   l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q

   nstor = 0

   ! In the complex case tau(2) /= 0
   IF(my_prow == prow(1)) THEN
      q(1,1:l_cols) = q(1,1:l_cols)*((1.d0,0.d0)-tau(2))
   ENDIF

   DO istep=1,na,nblk

      ics = MAX(istep,3)
      ice = MIN(istep+nblk-1,na)
      IF(ice<ics) CYCLE

      cur_pcol = pcol(istep)

      nb = 0
      DO ic=ics,ice

         l_colh = local_index(ic  , my_pcol, np_cols, nblk, -1) ! Column of Householder vector
         l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder vector


         IF(my_pcol==cur_pcol) THEN
            hvb(nb+1:nb+l_rows) = a(1:l_rows,l_colh)
            IF(my_prow==prow(ic-1)) THEN
               hvb(nb+l_rows) = 1.
            ENDIF
         ENDIF

         nb = nb+l_rows
      ENDDO

      IF(nb>0) &
         CALL MPI_Bcast(hvb,nb,MPI_DOUBLE_COMPLEX,cur_pcol,mpi_comm_cols,mpierr)

      nb = 0
      DO ic=ics,ice
         l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder vector
         hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
         nstor = nstor+1
         nb = nb+l_rows
      ENDDO

      ! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
      IF(nstor+nblk>max_stored_rows .OR. istep+nblk>na .OR. (na/np_rows<=256 .AND. nstor>=32)) THEN

         ! Calculate scalar products of stored vectors.
         ! This can be done in different ways, we use zherk

         tmat = 0
         IF(l_rows>0) &
            CALL zherk('U','C',nstor,l_rows,CONE,hvm,UBOUND(hvm,1),CZERO,tmat,max_stored_rows)

         nc = 0
         DO n=1,nstor-1
            h1(nc+1:nc+n) = tmat(1:n,n+1)
            nc = nc+n
         ENDDO

         IF(nc>0) CALL mpi_allreduce(h1,h2,nc,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)

         ! Calculate triangular matrix T

         nc = 0
         tmat(1,1) = tau(ice-nstor+1)
         DO n=1,nstor-1
            CALL ztrmv('L','C','N',n,tmat,max_stored_rows,h2(nc+1),1)
            tmat(n+1,1:n) = -CONJG(h2(nc+1:nc+n))*tau(ice-nstor+n+1)
            tmat(n+1,n+1) = tau(ice-nstor+n+1)
            nc = nc+n
         ENDDO

         ! Q = Q - V * T * V**T * Q

         IF(l_rows>0) THEN
            CALL zgemm('C','N',nstor,l_cols,l_rows,CONE,hvm,UBOUND(hvm,1), &
                       q,ldq,CZERO,tmp1,nstor)
         ELSE
            tmp1(1:l_cols*nstor) = 0
         ENDIF
         CALL mpi_allreduce(tmp1,tmp2,nstor*l_cols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)
         IF(l_rows>0) THEN
            CALL ztrmm('L','L','N','N',nstor,l_cols,CONE,tmat,max_stored_rows,tmp2,nstor)
            CALL zgemm('N','N',l_rows,l_cols,nstor,-CONE,hvm,UBOUND(hvm,1), &
                       tmp2,nstor,CONE,q,ldq)
         ENDIF
         nstor = 0
      ENDIF

   ENDDO

   DEALLOCATE(tmat, h1, h2, tmp1, tmp2, hvb, hvm)


END SUBROUTINE trans_ev_complex

!-------------------------------------------------------------------------------

SUBROUTINE mult_ah_b_complex(uplo_a, uplo_c, na, ncb, a, lda, b, ldb, nblk, mpi_comm_rows, mpi_comm_cols, c, ldc)

!-------------------------------------------------------------------------------
!  mult_ah_b_complex:  Performs C := A**H * B
!
!      where:  A is a square matrix (na,na) which is optionally upper or lower triangular
!              B is a (na,ncb) matrix
!              C is a (na,ncb) matrix where optionally only the upper or lower
!              triangle may be computed
!
!
!  uplo_a      'U' if A is upper triangular
!              'L' if A is lower triangular
!              anything else if A is a full matrix
!              Please note: This pertains to the original A (as set in the calling program)
!              whereas the transpose of A is used for calculations
!              If uplo_a is 'U' or 'L', the other triangle is not used at all,
!              i.e. it may contain arbitrary numbers
!
!  uplo_c      'U' if only the upper diagonal part of C is needed
!              'L' if only the upper diagonal part of C is needed
!              anything else if the full matrix C is needed
!              Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
!              written to a certain extent, i.e. one shouldn't rely on the content there!
!
!  na          Number of rows/columns of A, number of rows of B and C
!
!  ncb         Number of columns  of B and C
!
!  a           Matrix A
!
!  lda         Leading dimension of a
!
!  b           Matrix B
!
!  ldb         Leading dimension of b
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  c           Matrix C
!
!  ldc         Leading dimension of c
!
!-------------------------------------------------------------------------------


    CHARACTER(len=1)                         :: uplo_a, uplo_c
    INTEGER                                  :: na, ncb, lda
    COMPLEX*16                               :: a(lda,*)
    INTEGER                                  :: ldb
    COMPLEX*16                               :: b(ldb,*)
    INTEGER                                  :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols, ldc
    COMPLEX*16                               :: c(ldc,*)

    COMPLEX*16, ALLOCATABLE                  :: aux_bc(:), aux_mat(:,:), &
                                                tmp1(:,:), tmp2(:,:)
    INTEGER :: gcol, gcol_min, goff, l_cols, l_rows, l_rows_np, lce, lcs, &
      lre, lrs, mpierr, my_pcol, my_prow, n, n_aux_bc, nb, nblk_mult, noff, &
      np, np_bc, np_cols, np_rows, nr_done, nstor, nvals
    INTEGER, ALLOCATABLE                     :: lre_save(:), lrs_save(:)
    LOGICAL                                  :: a_lower, a_upper, c_lower, &
                                                c_upper

   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   l_rows = local_index(na,  my_prow, np_rows, nblk, -1) ! Local rows of a and b
   l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b

   ! Block factor for matrix multiplications, must be a multiple of nblk

   IF(na/np_rows<=256) THEN
      nblk_mult = (31/nblk+1)*nblk
   ELSE
      nblk_mult = (63/nblk+1)*nblk
   ENDIF

   ALLOCATE(aux_mat(l_rows,nblk_mult))
   ALLOCATE(aux_bc(l_rows*nblk))
   ALLOCATE(lrs_save(nblk))
   ALLOCATE(lre_save(nblk))

   a_lower = .FALSE.
   a_upper = .FALSE.
   c_lower = .FALSE.
   c_upper = .FALSE.

   IF(uplo_a=='u' .OR. uplo_a=='U') a_upper = .TRUE.
   IF(uplo_a=='l' .OR. uplo_a=='L') a_lower = .TRUE.
   IF(uplo_c=='u' .OR. uplo_c=='U') c_upper = .TRUE.
   IF(uplo_c=='l' .OR. uplo_c=='L') c_lower = .TRUE.

   ! Build up the result matrix by processor rows

   DO np = 0, np_rows-1

      ! In this turn, procs of row np assemble the result

      l_rows_np = local_index(na, np, np_rows, nblk, -1) ! local rows on receiving processors

      nr_done = 0 ! Number of rows done
      aux_mat = 0
      nstor = 0   ! Number of columns stored in aux_mat

      ! Loop over the blocks on row np

      DO nb=0,(l_rows_np-1)/nblk

         goff  = nb*np_rows + np ! Global offset in blocks corresponding to nb

         ! Get the processor column which owns this block (A is transposed, so we need the column)
         ! and the offset in blocks within this column.
         ! The corresponding block column in A is then broadcast to all for multiplication with B

         np_bc = MOD(goff,np_cols)
         noff = goff/np_cols
         n_aux_bc = 0

         ! Gather up the complete block column of A on the owner

         DO n = 1, MIN(l_rows_np-nb*nblk,nblk) ! Loop over columns to be broadcast

            gcol = goff*nblk + n ! global column corresponding to n
            IF(nstor==0 .AND. n==1) gcol_min = gcol

            lrs = 1       ! 1st local row number for broadcast
            lre = l_rows  ! last local row number for broadcast
            IF(a_lower) lrs = local_index(gcol, my_prow, np_rows, nblk, +1)
            IF(a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

            IF(lrs<=lre) THEN
               nvals = lre-lrs+1
               IF(my_pcol == np_bc) aux_bc(n_aux_bc+1:n_aux_bc+nvals) = a(lrs:lre,noff*nblk+n)
               n_aux_bc = n_aux_bc + nvals
            ENDIF

            lrs_save(n) = lrs
            lre_save(n) = lre

         ENDDO

         ! Broadcast block column

         CALL MPI_Bcast(aux_bc,n_aux_bc,MPI_DOUBLE_COMPLEX,np_bc,mpi_comm_cols,mpierr)

         ! Insert what we got in aux_mat

         n_aux_bc = 0
         DO n = 1, MIN(l_rows_np-nb*nblk,nblk)
            nstor = nstor+1
            lrs = lrs_save(n)
            lre = lre_save(n)
            IF(lrs<=lre) THEN
               nvals = lre-lrs+1
               aux_mat(lrs:lre,nstor) = aux_bc(n_aux_bc+1:n_aux_bc+nvals)
               n_aux_bc = n_aux_bc + nvals
            ENDIF
         ENDDO

         ! If we got nblk_mult columns in aux_mat or this is the last block
         ! do the matrix multiplication

         IF(nstor==nblk_mult .OR. nb*nblk+nblk >= l_rows_np) THEN

            lrs = 1       ! 1st local row number for multiply
            lre = l_rows  ! last local row number for multiply
            IF(a_lower) lrs = local_index(gcol_min, my_prow, np_rows, nblk, +1)
            IF(a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

            lcs = 1       ! 1st local col number for multiply
            lce = l_cols  ! last local col number for multiply
            IF(c_upper) lcs = local_index(gcol_min, my_pcol, np_cols, nblk, +1)
            IF(c_lower) lce = MIN(local_index(gcol, my_pcol, np_cols, nblk, -1),l_cols)

            IF(lcs<=lce) THEN
               ALLOCATE(tmp1(nstor,lcs:lce),tmp2(nstor,lcs:lce))
               IF(lrs<=lre) THEN
                  CALL zgemm('C','N',nstor,lce-lcs+1,lre-lrs+1,(1.d0,0.d0),aux_mat(lrs,1),UBOUND(aux_mat,1), &
                             b(lrs,lcs),ldb,(0.d0,0.d0),tmp1,nstor)
               ELSE
                  tmp1 = 0
               ENDIF

               ! Sum up the results and send to processor row np
               CALL mpi_reduce(tmp1,tmp2,nstor*(lce-lcs+1),MPI_DOUBLE_COMPLEX,MPI_SUM,np,mpi_comm_rows,mpierr)

               ! Put the result into C
               IF(my_prow==np) c(nr_done+1:nr_done+nstor,lcs:lce) = tmp2(1:nstor,lcs:lce)

               DEALLOCATE(tmp1,tmp2)
            ENDIF

            nr_done = nr_done+nstor
            nstor=0
            aux_mat(:,:)=0
         ENDIF
      ENDDO
   ENDDO

   DEALLOCATE(aux_mat, aux_bc, lrs_save, lre_save)

END SUBROUTINE mult_ah_b_complex

!-------------------------------------------------------------------------------

SUBROUTINE solve_tridi( na, nev, d, e, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols )

   IMPLICIT NONE

   INTEGER  na, nev, ldq, nblk, mpi_comm_rows, mpi_comm_cols
   REAL*8 d(na), e(na), q(ldq,*)

   INTEGER i, j, n, np, nc, nev1, l_cols, l_rows
   INTEGER my_prow, my_pcol, np_rows, np_cols, mpierr

   INTEGER, ALLOCATABLE :: limits(:), l_col(:), p_col(:), l_col_bc(:), p_col_bc(:)


   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q

   ! Set Q to 0

   q(1:l_rows, 1:l_cols) = 0.

   ! Get the limits of the subdivisons, each subdivison has as many cols
   ! as fit on the respective processor column

   ALLOCATE(limits(0:np_cols))

   limits(0) = 0
   DO np=0,np_cols-1
      nc = local_index(na, np, np_cols, nblk, -1) ! number of columns on proc column np

      ! Check for the case that a column has have zero width.
      ! This is not supported!
      ! Scalapack supports it but delivers no results for these columns,
      ! which is rather annoying
      IF(nc==0) THEN
         PRINT *,'ERROR: Problem contains processor column with zero width'
         CALL MPI_Abort(MPI_COMM_WORLD,1,mpierr)
      ENDIF

      limits(np+1) = limits(np) + nc
   ENDDO

   ! Subdivide matrix by subtracting rank 1 modifications

   DO i=1,np_cols-1
      n = limits(i)
      d(n) = d(n)-ABS(e(n))
      d(n+1) = d(n+1)-ABS(e(n))
   ENDDO

   ! Solve sub problems on processsor columns

   nc = limits(my_pcol) ! column after which my problem starts

   IF(np_cols>1) THEN
      nev1 = l_cols ! all eigenvectors are needed
   ELSE
      nev1 = MIN(nev,l_cols)
   ENDIF
   CALL solve_tridi_col(l_cols, nev1, nc, d(nc+1), e(nc+1), q, ldq, nblk, mpi_comm_rows)


   ! If there is only 1 processor column, we are done

   IF(np_cols==1) THEN
      DEALLOCATE(limits)
      RETURN
   ENDIF

   ! Set index arrays for Q columns

   ! Dense distribution scheme:

   ALLOCATE(l_col(na))
   ALLOCATE(p_col(na))

   n = 0
   DO np=0,np_cols-1
      nc = local_index(na, np, np_cols, nblk, -1)
      DO i=1,nc
         n = n+1
         l_col(n) = i
         p_col(n) = np
      ENDDO
   ENDDO

   ! Block cyclic distribution scheme, only nev columns are set:

   ALLOCATE(l_col_bc(na))
   ALLOCATE(p_col_bc(na))
   p_col_bc(:) = -1
   l_col_bc(:) = -1

   DO i = 0, na-1, nblk*np_cols
      DO j = 0, np_cols-1
         DO n = 1, nblk
            IF(i+j*nblk+n <= MIN(nev,na)) THEN
               p_col_bc(i+j*nblk+n) = j
               l_col_bc(i+j*nblk+n) = i/np_cols + n
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   ! Recursively merge sub problems

   CALL merge_recursive(0, np_cols)

   DEALLOCATE(limits,l_col,p_col,l_col_bc,p_col_bc)

CONTAINS
RECURSIVE SUBROUTINE merge_recursive(np_off, nprocs)


   ! noff is always a multiple of nblk_ev
   ! nlen-noff is always > nblk_ev

    INTEGER                                  :: np_off, nprocs

    INTEGER                                  :: mpi_status(mpi_status_size), &
                                                n, nlen, nmid, noff, np1, np2

   IF(nprocs<=1) THEN
      ! Safety check only
      PRINT *,"INTERNAL error merge_recursive: nprocs=",nprocs
      CALL mpi_abort(MPI_COMM_WORLD,1,mpierr)
   ENDIF

   ! Split problem into 2 subproblems of size np1 / np2

   np1 = nprocs/2
   np2 = nprocs-np1

   IF(np1 > 1) CALL merge_recursive(np_off, np1)
   IF(np2 > 1) CALL merge_recursive(np_off+np1, np2)


   noff = limits(np_off)
   nmid = limits(np_off+np1) - noff
   nlen = limits(np_off+nprocs) - noff

   IF(my_pcol==np_off) THEN
      DO n=np_off+np1,np_off+nprocs-1
         CALL mpi_send(d(noff+1),nmid,MPI_REAL8,n,1,mpi_comm_cols,mpierr)
      ENDDO
   ENDIF
   IF(my_pcol>=np_off+np1 .AND. my_pcol<np_off+nprocs) THEN
      CALL mpi_recv(d(noff+1),nmid,MPI_REAL8,np_off,1,mpi_comm_cols,mpi_status,mpierr)
   ENDIF

   IF(my_pcol==np_off+np1) THEN
      DO n=np_off,np_off+np1-1
         CALL mpi_send(d(noff+nmid+1),nlen-nmid,MPI_REAL8,n,1,mpi_comm_cols,mpierr)
      ENDDO
   ENDIF
   IF(my_pcol>=np_off .AND. my_pcol<np_off+np1) THEN
      CALL mpi_recv(d(noff+nmid+1),nlen-nmid,MPI_REAL8,np_off+np1,1,mpi_comm_cols,mpi_status,mpierr)
   ENDIF

   IF(nprocs == np_cols) THEN

      ! Last merge, result distribution must be block cyclic, noff==0,
      ! p_col_bc is set so that only nev eigenvalues are calculated

      CALL merge_systems(nlen, nmid, d(noff+1), e(noff+nmid), q, ldq, noff, &
                         nblk, mpi_comm_rows, mpi_comm_cols, l_col, p_col, &
                         l_col_bc, p_col_bc, np_off, nprocs )

   ELSE

      ! Not last merge, leave dense column distribution

      CALL merge_systems(nlen, nmid, d(noff+1), e(noff+nmid), q, ldq, noff, &
                         nblk, mpi_comm_rows, mpi_comm_cols, l_col(noff+1), p_col(noff+1), &
                         l_col(noff+1), p_col(noff+1), np_off, nprocs )
   ENDIF

END SUBROUTINE merge_recursive

END SUBROUTINE solve_tridi

!-------------------------------------------------------------------------------

SUBROUTINE solve_tridi_col( na, nev, nqoff, d, e, q, ldq, nblk, mpi_comm_rows )

   ! Solves the symmetric, tridiagonal eigenvalue problem on one processor column
   ! with the divide and conquer method.
   ! Works best if the number of processor rows is a power of 2!


    INTEGER                                  :: na, nev, nqoff
    REAL*8                                   :: d(na), e(na)
    INTEGER                                  :: ldq
    REAL*8                                   :: q(ldq,*)
    INTEGER                                  :: nblk, mpi_comm_rows

    INTEGER, PARAMETER                       :: min_submatrix_size = 16 

    INTEGER                                  :: i, max_size, mpierr, my_prow, &
                                                n, ndiv, nlen, nmid, noff, &
                                                np, np_rows
    INTEGER, ALLOCATABLE                     :: l_col(:), limits(:), &
                                                p_col_i(:), p_col_o(:)
    REAL*8, ALLOCATABLE                      :: qmat1(:,:), qmat2(:,:)

! Minimum size of the submatrices to be used

   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)

   ! Calculate the number of subdivisions needed.

   n = na
   ndiv = 1
   DO WHILE(2*ndiv<=np_rows .AND. n>2*min_submatrix_size)
      n = ((n+3)/4)*2 ! the bigger one of the two halves, we want EVEN boundaries
      ndiv = ndiv*2
   ENDDO

   ! If there is only 1 processor row and not all eigenvectors are needed
   ! and the matrix size is big enough, then use 2 subdivisions
   ! so that merge_systems is called once and only the needed
   ! eigenvectors are calculated for the final problem.

   IF(np_rows==1 .AND. nev<na .AND. na>2*min_submatrix_size) ndiv = 2

   ALLOCATE(limits(0:ndiv))

   limits(0) = 0
   limits(ndiv) = na

   n = ndiv
   DO WHILE(n>1)
      n = n/2 ! n is always a power of 2
      DO i=0,ndiv-1,2*n
         ! We want to have even boundaries (for cache line alignments)
         limits(i+n) = limits(i) + ((limits(i+2*n)-limits(i)+3)/4)*2
      ENDDO
   ENDDO

   ! Calculate the maximum size of a subproblem

   max_size = 0
   DO i=1,ndiv
      max_size = MAX(max_size,limits(i)-limits(i-1))
   ENDDO

   ! Subdivide matrix by subtracting rank 1 modifications

   DO i=1,ndiv-1
      n = limits(i)
      d(n) = d(n)-ABS(e(n))
      d(n+1) = d(n+1)-ABS(e(n))
   ENDDO

   IF(np_rows==1)    THEN

      ! For 1 processor row there may be 1 or 2 subdivisions

      DO n=0,ndiv-1
         noff = limits(n)        ! Start of subproblem
         nlen = limits(n+1)-noff ! Size of subproblem

         CALL solve_tridi_single(nlen,d(noff+1),e(noff+1),q(nqoff+noff+1,noff+1),UBOUND(q,1))
      ENDDO

   ELSE

      ! Solve sub problems in parallel with solve_tridi_single
      ! There is at maximum 1 subproblem per processor

      ALLOCATE(qmat1(max_size,max_size))
      ALLOCATE(qmat2(max_size,max_size))

      qmat1 = 0 ! Make sure that all elements are defined

      IF(my_prow < ndiv) THEN

         noff = limits(my_prow)        ! Start of subproblem
         nlen = limits(my_prow+1)-noff ! Size of subproblem

         CALL solve_tridi_single(nlen,d(noff+1),e(noff+1),qmat1,UBOUND(qmat1,1))

      ENDIF

      ! Fill eigenvectors in qmat1 into global matrix q

      DO np = 0, ndiv-1

         noff = limits(np)
         nlen = limits(np+1)-noff

         CALL MPI_Bcast(d(noff+1),nlen,MPI_REAL8,np,mpi_comm_rows,mpierr)
         qmat2 = qmat1
         CALL MPI_Bcast(qmat2,max_size*max_size,MPI_REAL8,np,mpi_comm_rows,mpierr)

         DO i=1,nlen
            CALL distribute_global_column(qmat2(1,i), q(1,noff+i), nqoff+noff, nlen, my_prow, np_rows, nblk)
         ENDDO

      ENDDO

      DEALLOCATE(qmat1, qmat2)

   ENDIF

   ! Allocate and set index arrays l_col and p_col

   ALLOCATE(l_col(na), p_col_i(na),  p_col_o(na))

   DO i=1,na
      l_col(i) = i
      p_col_i(i) = 0
      p_col_o(i) = 0
   ENDDO

   ! Merge subproblems

   n = 1
   DO WHILE(n<ndiv) ! if ndiv==1, the problem was solved by single call to solve_tridi_single

      DO i=0,ndiv-1,2*n

         noff = limits(i)
         nmid = limits(i+n) - noff
         nlen = limits(i+2*n) - noff

        IF(nlen == na) THEN
           ! Last merge, set p_col_o=-1 for unneeded (output) eigenvectors
           p_col_o(nev+1:na) = -1
        ENDIF

         CALL merge_systems(nlen, nmid, d(noff+1), e(noff+nmid), q, ldq, nqoff+noff, nblk, &
                            mpi_comm_rows, mpi_comm_self, l_col(noff+1), p_col_i(noff+1), &
                            l_col(noff+1), p_col_o(noff+1), 0, 1)

      ENDDO

      n = 2*n

   ENDDO

   DEALLOCATE(limits, l_col, p_col_i, p_col_o)

END SUBROUTINE solve_tridi_col

!-------------------------------------------------------------------------------

SUBROUTINE solve_tridi_single(nlen, d, e, q, ldq)

   ! Solves the symmetric, tridiagonal eigenvalue problem on a single processor.
   ! Takes precautions if DSTEDC fails or if the eigenvalues are not ordered correctly.


    INTEGER                                  :: nlen
    REAL*8                                   :: d(nlen), e(nlen)
    INTEGER                                  :: ldq
    REAL*8                                   :: q(ldq,nlen)

    INTEGER                                  :: i, info, j, liwork, lwork, &
                                                mpierr
    INTEGER, ALLOCATABLE                     :: iwork(:)
    REAL*8                                   :: dtmp
    REAL*8, ALLOCATABLE                      :: ds(:), es(:), qtmp(:), work(:)

   ALLOCATE(ds(nlen), es(nlen))

   ! Save d and e for the case that dstedc fails

   ds(:) = d(:)
   es(:) = e(:)

   ! First try dstedc, this is normally faster but it may fail sometimes (why???)

   lwork = 1 + 4*nlen + nlen**2
   liwork =  3 + 5*nlen
   ALLOCATE(work(lwork), iwork(liwork))
   CALL dstedc('I',nlen,d,e,q,ldq,work,lwork,iwork,liwork,info)

   IF(info /= 0) THEN

      ! DSTEDC failed, try DSTEQR. The workspace is enough for DSTEQR.

      PRINT '(a,i8,a)','Warning: Lapack routine DSTEDC failed, info= ',info,', Trying DSTEQR!'

      d(:) = ds(:)
      e(:) = es(:)
      CALL dsteqr('I',nlen,d,e,q,ldq,work,info)

      ! If DSTEQR fails also, we don't know what to do further ...
      IF(info /= 0) THEN
         PRINT '(a,i8,a)','ERROR: Lapack routine DSTEQR failed, info= ',info,', Aborting!'
         CALL mpi_abort(mpi_comm_world,0,mpierr)
      ENDIF

   END IF

   DEALLOCATE(work,iwork,ds,es)

   ! Check if eigenvalues are monotonically increasing
   ! This seems to be not always the case  (in the IBM implementation of dstedc ???)

   DO i=1,nlen-1
      IF(d(i+1)<d(i)) THEN
         IF (ABS(d(i+1) - d(i)) / ABS(d(i+1) + d(i)) > 1d-14) THEN
            PRINT '(a,i8,2g25.16)','***WARNING: Monotony error dste**:',i+1,d(i),d(i+1)
         ELSE
            PRINT '(a,i8,2g25.16)','Info: Monotony error dste{dc,qr}:',i+1,d(i),d(i+1)
            PRINT '(a)', 'The eigenvalues from a lapack call are not sorted to machine precision.'
            PRINT '(a)', 'In this extent, this is completely harmless.'
            PRINT '(a)', 'Still, we keep this info message just in case.'
         END IF
         ALLOCATE(qtmp(nlen))
         dtmp = d(i+1)
         qtmp(1:nlen) = q(1:nlen,i+1)
         DO j=i,1,-1
            IF(dtmp<d(j)) THEN
               d(j+1)        = d(j)
               q(1:nlen,j+1) = q(1:nlen,j)
            ELSE
               EXIT ! Loop
            ENDIF
         ENDDO
         d(j+1)        = dtmp
         q(1:nlen,j+1) = qtmp(1:nlen)
         DEALLOCATE(qtmp)
      ENDIF
   ENDDO

END SUBROUTINE solve_tridi_single

!-------------------------------------------------------------------------------

SUBROUTINE merge_systems( na, nm, d, e, q, ldq, nqoff, nblk, mpi_comm_rows, mpi_comm_cols, &
                          l_col, p_col, l_col_out, p_col_out, npc_0, npc_n)

   IMPLICIT NONE

   INTEGER  na, nm, ldq, nqoff, nblk, mpi_comm_rows, mpi_comm_cols, npc_0, npc_n
   INTEGER  l_col(na), p_col(na), l_col_out(na), p_col_out(na)
   REAL*8 d(na), e, q(ldq,*)

   INTEGER, PARAMETER :: max_strip=128

   REAL*8 beta, sig, s, c, t, tau, rho, eps, tol, dlamch, dlapy2, qtrans(2,2), dmax, zmax, d1new, d2new
   REAL*8 z(na), d1(na), d2(na), z1(na), delta(na), dbase(na), ddiff(na), ev_scale(na), tmp(na)
   REAL*8 d1u(na), zu(na), d1l(na), zl(na)
   REAL*8, ALLOCATABLE :: qtmp1(:,:), qtmp2(:,:), ev(:,:)
   REAL*8, ALLOCATABLE :: z_p(:,:)

   INTEGER i, j, na1, na2, l_rows, l_cols, l_rqs, l_rqe, l_rqm, ns, info
   INTEGER l_rnm, nnzu, nnzl, ndef, ncnt, max_local_cols, l_cols_qreorg, np, l_idx, nqcols1, nqcols2
   INTEGER my_proc, n_procs, my_prow, my_pcol, np_rows, np_cols, mpierr, mpi_status(mpi_status_size)
   INTEGER np_next, np_prev, np_rem
   INTEGER idx(na), idx1(na), idx2(na)
   INTEGER coltyp(na), idxq1(na), idxq2(na)

   INTEGER max_threads, my_thread
!$ integer omp_get_max_threads, omp_get_thread_num

   max_threads = 1
!$ max_threads = omp_get_max_threads()
   ALLOCATE(z_p(na,0:max_threads-1))

   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! If my processor column isn't in the requested set, do nothing

   IF(my_pcol<npc_0 .OR. my_pcol>=npc_0+npc_n) RETURN

   ! Determine number of "next" and "prev" column for ring sends

   IF(my_pcol == npc_0+npc_n-1) THEN
      np_next = npc_0
   ELSE
      np_next = my_pcol + 1
   ENDIF

   IF(my_pcol == npc_0) THEN
      np_prev = npc_0+npc_n-1
   ELSE
      np_prev = my_pcol - 1
   ENDIF

   CALL check_monotony(nm,d,'Input1')
   CALL check_monotony(na-nm,d(nm+1),'Input2')

   ! Get global number of processors and my processor number.
   ! Please note that my_proc does not need to match any real processor number,
   ! it is just used for load balancing some loops.

   n_procs = np_rows*npc_n
   my_proc = my_prow*npc_n + (my_pcol-npc_0) ! Row major


   ! Local limits of the rows of Q

   l_rqs = local_index(nqoff+1 , my_prow, np_rows, nblk, +1) ! First row of Q
   l_rqm = local_index(nqoff+nm, my_prow, np_rows, nblk, -1) ! Last row <= nm
   l_rqe = local_index(nqoff+na, my_prow, np_rows, nblk, -1) ! Last row of Q

   l_rnm  = l_rqm-l_rqs+1 ! Number of local rows <= nm
   l_rows = l_rqe-l_rqs+1 ! Total number of local rows


   ! My number of local columns

   l_cols = COUNT(p_col(1:na)==my_pcol)

   ! Get max number of local columns

   max_local_cols = 0
   DO np = npc_0, npc_0+npc_n-1
      max_local_cols = MAX(max_local_cols,COUNT(p_col(1:na)==np))
   ENDDO



   ! Calculations start here

   beta = ABS(e)
   sig  = SIGN(1.d0,e)

   ! Calculate rank-1 modifier z

   z(:) = 0

   IF(MOD((nqoff+nm-1)/nblk,np_rows)==my_prow) THEN
      ! nm is local on my row
      DO i = 1, na
         IF(p_col(i)==my_pcol) z(i) = q(l_rqm,l_col(i))
      ENDDO
   ENDIF

   IF(MOD((nqoff+nm)/nblk,np_rows)==my_prow) THEN
      ! nm+1 is local on my row
      DO i = 1, na
         IF(p_col(i)==my_pcol) z(i) = z(i) + sig*q(l_rqm+1,l_col(i))
      ENDDO
   ENDIF

   CALL global_gather(z, na)

   ! Normalize z so that norm(z) = 1.  Since z is the concatenation of
   ! two normalized vectors, norm2(z) = sqrt(2).

   z = z/SQRT(2.0d0)
   rho = 2.*beta

   ! Calculate index for merging both systems by ascending eigenvalues

   CALL DLAMRG( nm, na-nm, d, 1, 1, idx )

   ! Calculate the allowable deflation tolerance

   zmax = MAXVAL(ABS(z))
   dmax = MAXVAL(ABS(d))
   EPS = DLAMCH( 'Epsilon' )
   TOL = 8.*EPS*MAX(dmax,zmax)

   ! If the rank-1 modifier is small enough, no more needs to be done
   ! except to reorganize D and Q

   IF( RHO*zmax <= TOL ) THEN

      ! Rearrange eigenvalues

      tmp = d
      DO i=1,na
         d(i) = tmp(idx(i))
      ENDDO

      ! Rearrange eigenvectors

      CALL resort_ev(idx)

      RETURN
   ENDIF

   ! Merge and deflate system

   na1 = 0
   na2 = 0

   ! COLTYP:
   ! 1 : non-zero in the upper half only;
   ! 2 : dense;
   ! 3 : non-zero in the lower half only;
   ! 4 : deflated.

   coltyp(1:nm) = 1
   coltyp(nm+1:na) = 3

   DO i=1,na

      IF(rho*ABS(z(idx(i))) <= tol) THEN

         ! Deflate due to small z component.

         na2 = na2+1
         d2(na2)   = d(idx(i))
         idx2(na2) = idx(i)
         coltyp(idx(i)) = 4

      ELSE IF(na1>0) THEN

         ! Check if eigenvalues are close enough to allow deflation.

         S = Z(idx(i))
         C = Z1(na1)

         ! Find sqrt(a**2+b**2) without overflow or
         ! destructive underflow.

         TAU = DLAPY2( C, S )
         T = D1(na1) - D(idx(i))
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ) <= TOL ) THEN

            ! Deflation is possible.

            na2 = na2+1

            Z1(na1) = TAU

            d2new = D(idx(i))*C**2 + D1(na1)*S**2
            d1new = D(idx(i))*S**2 + D1(na1)*C**2

            ! D(idx(i)) >= D1(na1) and C**2 + S**2 == 1.0
            ! This means that after the above transformation it must be
            !    D1(na1) <= d1new <= D(idx(i))
            !    D1(na1) <= d2new <= D(idx(i))
            !
            ! D1(na1) may get bigger but it is still smaller than the next D(idx(i+1))
            ! so there is no problem with sorting here.
            ! d2new <= D(idx(i)) which means that it might be smaller than D2(na2-1)
            ! which makes a check (and possibly a resort) necessary.
            !
            ! The above relations may not hold exactly due to numeric differences
            ! so they have to be enforced in order not to get troubles with sorting.


            IF(d1new<D1(na1)  ) d1new = D1(na1)
            IF(d1new>D(idx(i))) d1new = D(idx(i))

            IF(d2new<D1(na1)  ) d2new = D1(na1)
            IF(d2new>D(idx(i))) d2new = D(idx(i))

            D1(na1) = d1new

            DO j=na2-1,1,-1
               IF(d2new<d2(j)) THEN
                  d2(j+1)   = d2(j)
                  idx2(j+1) = idx2(j)
               ELSE
                  EXIT ! Loop
               ENDIF
            ENDDO

            d2(j+1)   = d2new
            idx2(j+1) = idx(i)

            qtrans(1,1) = C; qtrans(1,2) =-S
            qtrans(2,1) = S; qtrans(2,2) = C

            CALL transform_columns(idx(i), idx1(na1))

            IF(coltyp(idx(i))==1 .AND. coltyp(idx1(na1))/=1) coltyp(idx1(na1)) = 2
            IF(coltyp(idx(i))==3 .AND. coltyp(idx1(na1))/=3) coltyp(idx1(na1)) = 2

            coltyp(idx(i)) = 4

         ELSE
            na1 = na1+1
            d1(na1) = d(idx(i))
            z1(na1) = z(idx(i))
            idx1(na1) = idx(i)
         ENDIF
      ELSE
         na1 = na1+1
         d1(na1) = d(idx(i))
         z1(na1) = z(idx(i))
         idx1(na1) = idx(i)
      ENDIF

   ENDDO
   CALL check_monotony(na1,d1,'Sorted1')
   CALL check_monotony(na2,d2,'Sorted2')

   IF(na1==1 .OR. na1==2) THEN
      ! if(my_proc==0) print *,'--- Remark solve_tridi: na1==',na1,' proc==',myid

      IF(na1==1) THEN
         d(1) = d1(1) + rho*z1(1)**2 ! solve secular equation
      ELSE ! na1==2
         CALL DLAED5(1, d1, z1, qtrans(1,1), rho, d(1))
         CALL DLAED5(2, d1, z1, qtrans(1,2), rho, d(2))

         CALL transform_columns(idx1(1), idx1(2))
      ENDIF

      ! Add the deflated eigenvalues
      d(na1+1:na) = d2(1:na2)

      ! Calculate arrangement of all eigenvalues  in output

      CALL DLAMRG( na1, na-na1, d, 1, 1, idx )

      ! Rearrange eigenvalues

      tmp = d
      DO i=1,na
         d(i) = tmp(idx(i))
      ENDDO

      ! Rearrange eigenvectors

      DO i=1,na
         IF(idx(i)<=na1) THEN
            idxq1(i) = idx1(idx(i))
         ELSE
            idxq1(i) = idx2(idx(i)-na1)
         ENDIF
      ENDDO

      CALL resort_ev(idxq1)

   ELSE IF(na1>2) THEN

      ! Solve secular equation

      z(1:na1) = 1
      z_p(1:na1,:) = 1
      dbase(1:na1) = 0
      ddiff(1:na1) = 0

      info = 0

!$OMP PARALLEL PRIVATE(i,my_thread,delta,s,info,j)
      my_thread = 0
!$    my_thread = omp_get_thread_num()
!$OMP DO
      DO i = my_proc+1, na1, n_procs ! work distributed over all processors

         CALL DLAED4(na1, i, d1, z1, delta, rho, s, info) ! s is not used!

         IF(info/=0) THEN
            ! If DLAED4 fails (may happen especially for LAPACK versions before 3.2)
            ! use the more stable bisection algorithm in solve_secular_equation
            ! print *,'ERROR DLAED4 n=',na1,'i=',i,' Using Bisection'
            CALL solve_secular_equation(na1, i, d1, z1, delta, rho, s)
         ENDIF

         ! Compute updated z

         DO j=1,na1
            IF(i/=j)  z_p(j,my_thread) = z_p(j,my_thread)*( delta(j) / (d1(j)-d1(i)) )
         ENDDO
         z_p(i,my_thread) = z_p(i,my_thread)*delta(i)

         ! store dbase/ddiff

         IF(i<na1) THEN
            IF(ABS(delta(i+1)) < ABS(delta(i))) THEN
               dbase(i) = d1(i+1)
               ddiff(i) = delta(i+1)
            ELSE
               dbase(i) = d1(i)
               ddiff(i) = delta(i)
            ENDIF
         ELSE
            dbase(i) = d1(i)
            ddiff(i) = delta(i)
         ENDIF
      ENDDO
!$OMP END PARALLEL
      DO i = 0, max_threads-1
         z(1:na1) = z(1:na1)*z_p(1:na1,i)
      ENDDO

      CALL global_product(z, na1)
      z(1:na1) = SIGN( SQRT( -z(1:na1) ), z1(1:na1) )

      CALL global_gather(dbase, na1)
      CALL global_gather(ddiff, na1)
      d(1:na1) = dbase(1:na1) - ddiff(1:na1)

      ! Calculate scale factors for eigenvectors

      ev_scale(:) = 0

!$OMP PARALLEL DO PRIVATE(i,tmp)
      DO i = my_proc+1, na1, n_procs ! work distributed over all processors

         ! tmp(1:na1) = z(1:na1) / delta(1:na1,i)  ! original code
         ! tmp(1:na1) = z(1:na1) / (d1(1:na1)-d(i))! bad results

         ! All we want to calculate is tmp = (d1(1:na1)-dbase(i))+ddiff(i)
         ! in exactly this order, but we want to prevent compiler optimization

         tmp(1:na1) = d1(1:na1)-dbase(i)
         CALL v_add_s(tmp,na1,ddiff(i))
         tmp(1:na1) = z(1:na1) / tmp(1:na1)
         ev_scale(i) = 1.0/SQRT(DOT_PRODUCT(tmp(1:na1),tmp(1:na1)))
      ENDDO
!$OMP END PARALLEL DO
      CALL global_gather(ev_scale, na1)

      ! Add the deflated eigenvalues
      d(na1+1:na) = d2(1:na2)

      ! Calculate arrangement of all eigenvalues  in output

      CALL DLAMRG( na1, na-na1, d, 1, 1, idx )

      ! Rearrange eigenvalues

      tmp = d
      DO i=1,na
         d(i) = tmp(idx(i))
      ENDDO
      CALL check_monotony(na,d,'Output')

      ! Eigenvector calculations


      ! Calculate the number of columns in the new local matrix Q
      ! which are updated from non-deflated/deflated eigenvectors.
      ! idxq1/2 stores the global column numbers.

      nqcols1 = 0 ! number of non-deflated eigenvectors
      nqcols2 = 0 ! number of deflated eigenvectors
      DO i = 1, na
         IF(p_col_out(i)==my_pcol) THEN
            IF(idx(i)<=na1) THEN
               nqcols1 = nqcols1+1
               idxq1(nqcols1) = i
            ELSE
               nqcols2 = nqcols2+1
               idxq2(nqcols2) = i
            ENDIF
         ENDIF
      ENDDO

      ALLOCATE(ev(max_local_cols,MIN(max_strip,MAX(1,nqcols1))))
      ALLOCATE(qtmp1(MAX(1,l_rows),max_local_cols))
      ALLOCATE(qtmp2(MAX(1,l_rows),MIN(max_strip,MAX(1,nqcols1))))

      ! Gather nonzero upper/lower components of old matrix Q
      ! which are needed for multiplication with new eigenvectors

      qtmp1 = 0 ! May contain empty (unset) parts
      qtmp2 = 0 ! Not really needed

      nnzu = 0
      nnzl = 0
      DO i = 1, na1
         l_idx = l_col(idx1(i))
         IF(p_col(idx1(i))==my_pcol) THEN
            IF(coltyp(idx1(i))==1 .OR. coltyp(idx1(i))==2) THEN
               nnzu = nnzu+1
               qtmp1(1:l_rnm,nnzu) = q(l_rqs:l_rqm,l_idx)
            ENDIF
            IF(coltyp(idx1(i))==3 .OR. coltyp(idx1(i))==2) THEN
               nnzl = nnzl+1
               qtmp1(l_rnm+1:l_rows,nnzl) = q(l_rqm+1:l_rqe,l_idx)
            ENDIF
         ENDIF
      ENDDO

      ! Gather deflated eigenvalues behind nonzero components

      ndef = MAX(nnzu,nnzl)
      DO i = 1, na2
         l_idx = l_col(idx2(i))
         IF(p_col(idx2(i))==my_pcol) THEN
            ndef = ndef+1
            qtmp1(1:l_rows,ndef) = q(l_rqs:l_rqe,l_idx)
         ENDIF
      ENDDO

      l_cols_qreorg = ndef ! Number of columns in reorganized matrix

      ! Set (output) Q to 0, it will sum up new Q

      DO i = 1, na
         IF(p_col_out(i)==my_pcol) q(l_rqs:l_rqe,l_col_out(i)) = 0
      ENDDO


      np_rem = my_pcol

      DO np = 1, npc_n

         ! Do a ring send of qtmp1

         IF(np>1) THEN

            IF(np_rem==npc_0) THEN
               np_rem = npc_0+npc_n-1
            ELSE
               np_rem = np_rem-1
            ENDIF

            CALL MPI_Sendrecv_replace(qtmp1, l_rows*max_local_cols, MPI_REAL8, &
                                      np_next, 1111, np_prev, 1111, &
                                      mpi_comm_cols, mpi_status, mpierr)
         ENDIF

         ! Gather the parts in d1 and z which are fitting to qtmp1.
         ! This also delivers nnzu/nnzl for proc np_rem

         nnzu = 0
         nnzl = 0
         DO i=1,na1
            IF(p_col(idx1(i))==np_rem) THEN
               IF(coltyp(idx1(i))==1 .OR. coltyp(idx1(i))==2) THEN
                  nnzu = nnzu+1
                  d1u(nnzu) = d1(i)
                  zu (nnzu) = z (i)
               ENDIF
               IF(coltyp(idx1(i))==3 .OR. coltyp(idx1(i))==2) THEN
                  nnzl = nnzl+1
                  d1l(nnzl) = d1(i)
                  zl (nnzl) = z (i)
               ENDIF
            ENDIF
         ENDDO

         ! Set the deflated eigenvectors in Q (comming from proc np_rem)

         ndef = MAX(nnzu,nnzl) ! Remote counter in input matrix
         DO i = 1, na
            j = idx(i)
            IF(j>na1) THEN
               IF(p_col(idx2(j-na1))==np_rem) THEN
                  ndef = ndef+1
                  IF(p_col_out(i)==my_pcol) &
                     q(l_rqs:l_rqe,l_col_out(i)) = qtmp1(1:l_rows,ndef)
               ENDIF
            ENDIF
         ENDDO

         DO ns = 0, nqcols1-1, max_strip ! strimining loop

            ncnt = MIN(max_strip,nqcols1-ns) ! number of columns in this strip

            ! Get partial result from (output) Q

            DO i = 1, ncnt
               qtmp2(1:l_rows,i) = q(l_rqs:l_rqe,l_col_out(idxq1(i+ns)))
            ENDDO

            ! Compute eigenvectors of the rank-1 modified matrix.
            ! Parts for multiplying with upper half of Q:

            DO i = 1, ncnt
               j = idx(idxq1(i+ns))
               ! Calculate the j-th eigenvector of the deflated system
               ! See above why we are doing it this way!
               tmp(1:nnzu) = d1u(1:nnzu)-dbase(j)
               CALL v_add_s(tmp,nnzu,ddiff(j))
               ev(1:nnzu,i) = zu(1:nnzu) / tmp(1:nnzu) * ev_scale(j)
            ENDDO

            ! Multiply old Q with eigenvectors (upper half)

            IF(l_rnm>0 .AND. ncnt>0 .AND. nnzu>0) &
               CALL dgemm('N','N',l_rnm,ncnt,nnzu,1.d0,qtmp1,UBOUND(qtmp1,1),ev,UBOUND(ev,1), &
                          1.d0,qtmp2(1,1),UBOUND(qtmp2,1))

            ! Compute eigenvectors of the rank-1 modified matrix.
            ! Parts for multiplying with lower half of Q:

            DO i = 1, ncnt
               j = idx(idxq1(i+ns))
               ! Calculate the j-th eigenvector of the deflated system
               ! See above why we are doing it this way!
               tmp(1:nnzl) = d1l(1:nnzl)-dbase(j)
               CALL v_add_s(tmp,nnzl,ddiff(j))
               ev(1:nnzl,i) = zl(1:nnzl) / tmp(1:nnzl) * ev_scale(j)
            ENDDO

            ! Multiply old Q with eigenvectors (lower half)

            IF(l_rows-l_rnm>0 .AND. ncnt>0 .AND. nnzl>0) &
               CALL dgemm('N','N',l_rows-l_rnm,ncnt,nnzl,1.d0,qtmp1(l_rnm+1,1),UBOUND(qtmp1,1),ev,UBOUND(ev,1), &
                          1.d0,qtmp2(l_rnm+1,1),UBOUND(qtmp2,1))

            ! Put partial result into (output) Q

            DO i = 1, ncnt
               q(l_rqs:l_rqe,l_col_out(idxq1(i+ns))) = qtmp2(1:l_rows,i)
            ENDDO

         ENDDO
      ENDDO

      DEALLOCATE(ev, qtmp1, qtmp2)

   ENDIF

!-------------------------------------------------------------------------------

CONTAINS
SUBROUTINE resort_ev(idx_ev)


    INTEGER                                  :: idx_ev(*)

    INTEGER                                  :: i, l_cols_out, lc1, lc2, nc, &
                                                pc1, pc2
    REAL*8, ALLOCATABLE                      :: qtmp(:,:)

   IF(l_rows==0) RETURN ! My processor column has no work to do

   ! Resorts eigenvectors so that q_new(:,i) = q_old(:,idx_ev(i))

   l_cols_out = COUNT(p_col_out(1:na)==my_pcol)
   ALLOCATE(qtmp(l_rows,l_cols_out))


   nc = 0

   DO i=1,na

      pc1 = p_col(idx_ev(i))
      lc1 = l_col(idx_ev(i))
      pc2 = p_col_out(i)

      IF(pc2<0) CYCLE ! This column is not needed in output

      IF(pc2==my_pcol) nc = nc+1 ! Counter for output columns

      IF(pc1==my_pcol) THEN
         IF(pc2==my_pcol) THEN
            ! send and recieve column are local
            qtmp(1:l_rows,nc) = q(l_rqs:l_rqe,lc1)
         ELSE
            CALL mpi_send(q(l_rqs,lc1),l_rows,MPI_REAL8,pc2,MOD(i,4096),mpi_comm_cols,mpierr)
         ENDIF
      ELSE IF(pc2==my_pcol) THEN
         CALL mpi_recv(qtmp(1,nc),l_rows,MPI_REAL8,pc1,MOD(i,4096),mpi_comm_cols,mpi_status,mpierr)
      ENDIF
   ENDDO

   ! Insert qtmp into (output) q

   nc = 0

   DO i=1,na

      pc2 = p_col_out(i)
      lc2 = l_col_out(i)

      IF(pc2==my_pcol) THEN
         nc = nc+1
         q(l_rqs:l_rqe,lc2) = qtmp(1:l_rows,nc)
      ENDIF
   ENDDO

   DEALLOCATE(qtmp)

END SUBROUTINE resort_ev

SUBROUTINE transform_columns(col1, col2)


    INTEGER                                  :: col1, col2

    INTEGER                                  :: lc1, lc2, pc1, pc2

   IF(l_rows==0) RETURN ! My processor column has no work to do

   pc1 = p_col(col1)
   lc1 = l_col(col1)
   pc2 = p_col(col2)
   lc2 = l_col(col2)

   IF(pc1==my_pcol) THEN
      IF(pc2==my_pcol) THEN
         ! both columns are local
         tmp(1:l_rows)      = q(l_rqs:l_rqe,lc1)*qtrans(1,1) + q(l_rqs:l_rqe,lc2)*qtrans(2,1)
         q(l_rqs:l_rqe,lc2) = q(l_rqs:l_rqe,lc1)*qtrans(1,2) + q(l_rqs:l_rqe,lc2)*qtrans(2,2)
         q(l_rqs:l_rqe,lc1) = tmp(1:l_rows)
      ELSE
         CALL mpi_sendrecv(q(l_rqs,lc1),l_rows,MPI_REAL8,pc2,1, &
                           tmp,l_rows,MPI_REAL8,pc2,1, &
                           mpi_comm_cols,mpi_status,mpierr)
         q(l_rqs:l_rqe,lc1) = q(l_rqs:l_rqe,lc1)*qtrans(1,1) + tmp(1:l_rows)*qtrans(2,1)
      ENDIF
   ELSE IF(pc2==my_pcol) THEN
      CALL mpi_sendrecv(q(l_rqs,lc2),l_rows,MPI_REAL8,pc1,1, &
                        tmp,l_rows,MPI_REAL8,pc1,1, &
                        mpi_comm_cols,mpi_status,mpierr)
      q(l_rqs:l_rqe,lc2) = tmp(1:l_rows)*qtrans(1,2) + q(l_rqs:l_rqe,lc2)*qtrans(2,2)
   ENDIF

END SUBROUTINE transform_columns

SUBROUTINE global_gather(z, n)

   ! This routine sums up z over all processors.
   ! It should only be used for gathering distributed results,
   ! i.e. z(i) should be nonzero on exactly 1 processor column,
   ! otherways the results may be numerically different on different columns


    INTEGER                                  :: n
    REAL*8                                   :: z(n)

    REAL*8                                   :: tmp(n)

   IF(npc_n==1 .AND. np_rows==1) RETURN ! nothing to do

   ! Do an mpi_allreduce over processor rows

   CALL mpi_allreduce(z, tmp, n, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)

   ! If only 1 processor column, we are done
   IF(npc_n==1) THEN
      z(:) = tmp(:)
      RETURN
   ENDIF

   ! If all processor columns are involved, we can use mpi_allreduce
   IF(npc_n==np_cols) THEN
      CALL mpi_allreduce(tmp, z, n, MPI_REAL8, MPI_SUM, mpi_comm_cols, mpierr)
      RETURN
   ENDIF

   ! Do a ring send over processor columns
   z(:) = 0
   DO np = 1, npc_n
      z(:) = z(:) + tmp(:)
      CALL MPI_Sendrecv_replace(z, n, MPI_REAL8, np_next, 1111, np_prev, 1111, &
                                mpi_comm_cols, mpi_status, mpierr)
   ENDDO

END SUBROUTINE global_gather

SUBROUTINE global_product(z, n)

   ! This routine calculates the global product of z.


    INTEGER                                  :: n
    REAL*8                                   :: z(n)

    REAL*8                                   :: tmp(n)

   IF(npc_n==1 .AND. np_rows==1) RETURN ! nothing to do

   ! Do an mpi_allreduce over processor rows

   CALL mpi_allreduce(z, tmp, n, MPI_REAL8, MPI_PROD, mpi_comm_rows, mpierr)

   ! If only 1 processor column, we are done
   IF(npc_n==1) THEN
      z(:) = tmp(:)
      RETURN
   ENDIF

   ! If all processor columns are involved, we can use mpi_allreduce
   IF(npc_n==np_cols) THEN
      CALL mpi_allreduce(tmp, z, n, MPI_REAL8, MPI_PROD, mpi_comm_cols, mpierr)
      RETURN
   ENDIF

   ! We send all vectors to the first proc, do the product there
   ! and redistribute the result.

   IF(my_pcol == npc_0) THEN
      z(1:n) = tmp(1:n)
      DO np = npc_0+1, npc_0+npc_n-1
         CALL mpi_recv(tmp,n,MPI_REAL8,np,1111,mpi_comm_cols,mpi_status,mpierr)
         z(1:n) = z(1:n)*tmp(1:n)
      ENDDO
      DO np = npc_0+1, npc_0+npc_n-1
         CALL mpi_send(z,n,MPI_REAL8,np,1111,mpi_comm_cols,mpierr)
      ENDDO
   ELSE
      CALL mpi_send(tmp,n,MPI_REAL8,npc_0,1111,mpi_comm_cols,mpierr)
      CALL mpi_recv(z  ,n,MPI_REAL8,npc_0,1111,mpi_comm_cols,mpi_status,mpierr)
   ENDIF

END SUBROUTINE global_product

SUBROUTINE check_monotony(n,d,text)

! This is a test routine for checking if the eigenvalues are monotonically increasing.
! It is for debug purposes only, an error should never be triggered!


    INTEGER                                  :: n
    REAL*8                                   :: d(n)
    CHARACTER(len=*)                         :: text

    INTEGER                                  :: i

   DO i=1,n-1
      IF(d(i+1)<d(i)) THEN
         PRINT '(a,a,i8,2g25.17)','Monotony error on ',text,i,d(i),d(i+1)
         CALL mpi_abort(mpi_comm_world,0,mpierr)
      ENDIF
   ENDDO

END SUBROUTINE check_monotony

END SUBROUTINE merge_systems

!-------------------------------------------------------------------------------

SUBROUTINE v_add_s(v,n,s)
    INTEGER                                  :: n
    REAL*8                                   :: v(n), s

   v(:) = v(:) + s
END SUBROUTINE v_add_s

!-------------------------------------------------------------------------------

SUBROUTINE distribute_global_column(g_col, l_col, noff, nlen, my_prow, np_rows, nblk)


    REAL*8                                   :: l_col(*)
    INTEGER                                  :: noff, nlen
    REAL*8                                   :: g_col(nlen)
    INTEGER                                  :: my_prow, np_rows, nblk

    INTEGER                                  :: g_off, jb, je, js, l_off, &
                                                nbe, nbs

   nbs = noff/(nblk*np_rows)
   nbe = (noff+nlen-1)/(nblk*np_rows)

   DO jb = nbs, nbe

      g_off = jb*nblk*np_rows + nblk*my_prow
      l_off = jb*nblk

      js = MAX(noff+1-g_off,1)
      je = MIN(noff+nlen-g_off,nblk)

      IF(je<js) CYCLE

      l_col(l_off+js:l_off+je) = g_col(g_off+js-noff:g_off+je-noff)

  ENDDO

END SUBROUTINE distribute_global_column

!-------------------------------------------------------------------------------

SUBROUTINE solve_secular_equation(n, i, d, z, delta, rho, dlam)

!-------------------------------------------------------------------------------
! This routine solves the secular equation of a symmetric rank 1 modified
! diagonal matrix:
!
!    1. + rho*SUM(z(:)**2/(d(:)-x)) = 0
!
! It does the same as the LAPACK routine DLAED4 but it uses a bisection technique
! which is more robust (it always yields a solution) but also slower
! than the algorithm used in DLAED4.
!
! The same restictions than in DLAED4 hold, namely:
!
!   rho > 0   and   d(i+1) > d(i)
!
! but this routine will not terminate with error if these are not satisfied
! (it will normally converge to a pole in this case).
!
! The output in DELTA(j) is always (D(j) - lambda_I), even for the cases
! N=1 and N=2 which is not compatible with DLAED4.
! Thus this routine shouldn't be used for these cases as a simple replacement
! of DLAED4.
!
! The arguments are the same as in DLAED4 (with the exception of the INFO argument):
!
!
!  N      (input) INTEGER
!         The length of all arrays.
!
!  I      (input) INTEGER
!         The index of the eigenvalue to be computed.  1 <= I <= N.
!
!  D      (input) DOUBLE PRECISION array, dimension (N)
!         The original eigenvalues.  It is assumed that they are in
!         order, D(I) < D(J)  for I < J.
!
!  Z      (input) DOUBLE PRECISION array, dimension (N)
!         The components of the updating vector.
!
!  DELTA  (output) DOUBLE PRECISION array, dimension (N)
!         DELTA contains (D(j) - lambda_I) in its  j-th component.
!         See remark above about DLAED4 compatibility!
!
!  RHO    (input) DOUBLE PRECISION
!         The scalar in the symmetric updating formula.
!
!  DLAM   (output) DOUBLE PRECISION
!         The computed lambda_I, the I-th updated eigenvalue.
!-------------------------------------------------------------------------------



    INTEGER                                  :: n, i
    REAL*8                                   :: d(n), z(n), delta(n), rho, &
                                                dlam

    INTEGER                                  :: iter
    REAL*8                                   :: a, b, dshift, x, y

! In order to obtain sufficient numerical accuracy we have to shift the problem
! either by d(i) or d(i+1), whichever is closer to the solution
! Upper and lower bound of the shifted solution interval are a and b

   IF(i==n) THEN

      ! Special case: Last eigenvalue
      ! We shift always by d(n), lower bound is d(n),
      ! upper bound is determined by a guess:

      dshift = d(n)
      delta(:) = d(:) - dshift

      a = 0. ! delta(n)
      b = rho*SUM(z(:)**2) + 1. ! rho*SUM(z(:)**2) is the lower bound for the guess

   ELSE

      ! Other eigenvalues: lower bound is d(i), upper bound is d(i+1)
      ! We check the sign of the function in the midpoint of the interval
      ! in order to determine if eigenvalue is more close to d(i) or d(i+1)

      x = 0.5*(d(i)+d(i+1))
      y = 1. + rho*SUM(z(:)**2/(d(:)-x))

      IF(y>0) THEN
         ! solution is next to d(i)
         dshift = d(i)
      ELSE
         ! solution is next to d(i+1)
         dshift = d(i+1)
      ENDIF

      delta(:) = d(:) - dshift
      a = delta(i)
      b = delta(i+1)

   ENDIF

   ! Bisection:

   DO iter=1,200

      ! Interval subdivision

      x = 0.5*(a+b)

      IF(x==a .OR. x==b) EXIT   ! No further interval subdivisions possible
      IF(ABS(x) < 1.d-200) EXIT ! x next to pole

      ! evaluate value at x

      y = 1. + rho*SUM(z(:)**2/(delta(:)-x))

      IF(y==0) THEN
         ! found exact solution
         EXIT
      ELSEIF(y>0) THEN
         b = x
      ELSE
         a = x
      ENDIF

   ENDDO

   ! Solution:

   dlam = x + dshift
   delta(:) = delta(:) - x

END SUBROUTINE solve_secular_equation

!-------------------------------------------------------------------------------

INTEGER FUNCTION local_index(idx, my_proc, num_procs, nblk, iflag)

!-------------------------------------------------------------------------------
!  local_index: returns the local index for a given global index
!               If the global index has no local index on the
!               processor my_proc behaviour is defined by iflag
!
!  Parameters
!
!  idx         Global index
!
!  my_proc     Processor row/column for which to calculate the local index
!
!  num_procs   Total number of processors along row/column
!
!  nblk        Blocksize
!
!  iflag       Controls the behaviour if idx is not on local processor
!              iflag< 0 : Return last local index before that row/col
!              iflag==0 : Return 0
!              iflag> 0 : Return next local index after that row/col
!-------------------------------------------------------------------------------

   IMPLICIT NONE

   INTEGER idx, my_proc, num_procs, nblk, iflag

   INTEGER iblk

   iblk = (idx-1)/nblk  ! global block number, 0 based

   IF(MOD(iblk,num_procs) == my_proc) THEN

      ! block is local, always return local row/col number

      local_index = (iblk/num_procs)*nblk + MOD(idx-1,nblk) + 1

   ELSE

      ! non local block

      IF(iflag == 0) THEN

         local_index = 0

      ELSE

         local_index = (iblk/num_procs)*nblk

         IF(MOD(iblk,num_procs) > my_proc) local_index = local_index + nblk

         IF(iflag>0) local_index = local_index + 1
      ENDIF
   ENDIF

END FUNCTION local_index

!-------------------------------------------------------------------------------

SUBROUTINE cholesky_real(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  cholesky_real: Cholesky factorization of a real symmetric matrix
!
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be factorized.
!              Distribution is like in Scalapack.
!              Only upper triangle is needs to be set.
!              On return, the upper triangle contains the Cholesky factor
!              and the lower triangle is set to 0.
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------


    INTEGER                                  :: na, lda
    REAL*8                                   :: a(lda,*)
    INTEGER                                  :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols

    INTEGER :: i, info, l_col1, l_cols, l_cols_tile, l_colx, l_row1, l_rows, &
      l_rows_tile, l_rowx, lce, lcs, lre, lrs, mpierr, my_pcol, my_prow, n, &
      nc, np_cols, np_rows, pcol, prow, tile_size
    REAL*8, ALLOCATABLE                      :: tmatc(:,:), tmatr(:,:), &
                                                tmp1(:), tmp2(:,:)

   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

   tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
   tile_size = ((128*MAX(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

   l_rows_tile = tile_size/np_rows ! local rows of a tile
   l_cols_tile = tile_size/np_cols ! local cols of a tile


   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

   ALLOCATE(tmp1(nblk*nblk))
   ALLOCATE(tmp2(nblk,nblk))
   tmp1 = 0
   tmp2 = 0

   ALLOCATE(tmatr(l_rows,nblk))
   ALLOCATE(tmatc(l_cols,nblk))
   tmatr = 0
   tmatc = 0


   DO n = 1, na, nblk

      ! Calculate first local row and column of the still remaining matrix
      ! on the local processor

      l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
      l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

      l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
      l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

      IF(n+nblk > na) THEN

         ! This is the last step, just do a Cholesky-Factorization
         ! of the remaining block

         IF(my_prow==prow(n) .AND. my_pcol==pcol(n)) THEN

            CALL dpotrf('U',na-n+1,a(l_row1,l_col1),lda,info)
            IF(info/=0) THEN
               PRINT *,"Error in dpotrf"
               CALL MPI_Abort(MPI_COMM_WORLD,1,mpierr)
            ENDIF

         ENDIF

         EXIT ! Loop

      ENDIF


      IF(my_prow==prow(n)) THEN

         IF(my_pcol==pcol(n)) THEN

            ! The process owning the upper left remaining block does the
            ! Cholesky-Factorization of this block

            CALL dpotrf('U',nblk,a(l_row1,l_col1),lda,info)
            IF(info/=0) THEN
               PRINT *,"Error in dpotrf"
               CALL MPI_Abort(MPI_COMM_WORLD,1,mpierr)
            ENDIF

            nc = 0
            DO i=1,nblk
               tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
               nc = nc+i
            ENDDO
         ENDIF

         CALL MPI_Bcast(tmp1,nblk*(nblk+1)/2,MPI_REAL8,pcol(n),mpi_comm_cols,mpierr)

         nc = 0
         DO i=1,nblk
            tmp2(1:i,i) = tmp1(nc+1:nc+i)
            nc = nc+i
         ENDDO

         IF(l_cols-l_colx+1>0) &
            CALL dtrsm('L','U','T','N',nblk,l_cols-l_colx+1,1.d0,tmp2,UBOUND(tmp2,1),a(l_row1,l_colx),lda)

      ENDIF

      DO i=1,nblk

         IF(my_prow==prow(n)) tmatc(l_colx:l_cols,i) = a(l_row1+i-1,l_colx:l_cols)
         IF(l_cols-l_colx+1>0) &
            CALL MPI_Bcast(tmatc(l_colx,i),l_cols-l_colx+1,MPI_REAL8,prow(n),mpi_comm_rows,mpierr)

      ENDDO

      CALL elpa_transpose_vectors  (tmatc, UBOUND(tmatc,1), mpi_comm_cols, &
                                    tmatr, UBOUND(tmatr,1), mpi_comm_rows, &
                                    n, na, nblk, nblk)

      DO i=0,(na-1)/tile_size
         lcs = MAX(l_colx,i*l_cols_tile+1)
         lce = MIN(l_cols,(i+1)*l_cols_tile)
         lrs = l_rowx
         lre = MIN(l_rows,(i+1)*l_rows_tile)
         IF(lce<lcs .OR. lre<lrs) CYCLE
         CALL DGEMM('N','T',lre-lrs+1,lce-lcs+1,nblk,-1.d0, &
                    tmatr(lrs,1),UBOUND(tmatr,1),tmatc(lcs,1),UBOUND(tmatc,1), &
                    1.d0,a(lrs,lcs),lda)
      ENDDO

   ENDDO

   DEALLOCATE(tmp1, tmp2, tmatr, tmatc)

   ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

   DO i=1,na
      IF(my_pcol==pcol(i)) THEN
         ! column i is on local processor
         l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
         l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
         a(l_row1:l_rows,l_col1) = 0
      ENDIF
   ENDDO

END SUBROUTINE cholesky_real

!-------------------------------------------------------------------------------

SUBROUTINE invert_trm_real(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  invert_trm_real: Inverts a upper triangular matrix
!
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be inverted.
!              Distribution is like in Scalapack.
!              Only upper triangle is needs to be set.
!              The lower triangle is not referenced.
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------


    INTEGER                                  :: na, lda
    REAL*8                                   :: a(lda,*)
    INTEGER                                  :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols

    INTEGER :: i, info, l_col1, l_cols, l_colx, l_row1, l_rows, l_rowx, &
      mpierr, my_pcol, my_prow, n, nb, nc, np_cols, np_rows, ns, pcol, prow
    REAL*8, ALLOCATABLE                      :: tmat1(:,:), tmat2(:,:), &
                                                tmp1(:), tmp2(:,:)

   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

   ALLOCATE(tmp1(nblk*nblk))
   ALLOCATE(tmp2(nblk,nblk))
   tmp1 = 0
   tmp2 = 0

   ALLOCATE(tmat1(l_rows,nblk))
   ALLOCATE(tmat2(nblk,l_cols))
   tmat1 = 0
   tmat2 = 0


   ns = ((na-1)/nblk)*nblk + 1

   DO n = ns,1,-nblk

      l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
      l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

      nb = nblk
      IF(na-n+1 < nblk) nb = na-n+1

      l_rowx = local_index(n+nb, my_prow, np_rows, nblk, +1)
      l_colx = local_index(n+nb, my_pcol, np_cols, nblk, +1)


      IF(my_prow==prow(n)) THEN

         IF(my_pcol==pcol(n)) THEN

            CALL DTRTRI('U','N',nb,a(l_row1,l_col1),lda,info)
            IF(info/=0) THEN
               PRINT *,"Error in DTRTRI"
               CALL MPI_Abort(MPI_COMM_WORLD,1,mpierr)
            ENDIF

            nc = 0
            DO i=1,nb
               tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
               nc = nc+i
            ENDDO
         ENDIF

         CALL MPI_Bcast(tmp1,nb*(nb+1)/2,MPI_REAL8,pcol(n),mpi_comm_cols,mpierr)

         nc = 0
         DO i=1,nb
            tmp2(1:i,i) = tmp1(nc+1:nc+i)
            nc = nc+i
         ENDDO

         IF(l_cols-l_colx+1>0) &
            CALL DTRMM('L','U','N','N',nb,l_cols-l_colx+1,1.d0,tmp2,UBOUND(tmp2,1),a(l_row1,l_colx),lda)

         IF(l_colx<=l_cols)   tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
         IF(my_pcol==pcol(n)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0

      ENDIF

      IF(l_row1>1) THEN
         IF(my_pcol==pcol(n)) THEN
            tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
            a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
         ENDIF

         DO i=1,nb
            CALL MPI_Bcast(tmat1(1,i),l_row1-1,MPI_REAL8,pcol(n),mpi_comm_cols,mpierr)
         ENDDO
      ENDIF

      IF(l_cols-l_col1+1>0) &
         CALL MPI_Bcast(tmat2(1,l_col1),(l_cols-l_col1+1)*nblk,MPI_REAL8,prow(n),mpi_comm_rows,mpierr)

      IF(l_row1>1 .AND. l_cols-l_col1+1>0) &
         CALL dgemm('N','N',l_row1-1,l_cols-l_col1+1,nb, -1.d0, &
                    tmat1,UBOUND(tmat1,1),tmat2(1,l_col1),UBOUND(tmat2,1), &
                    1.d0, a(1,l_col1),lda)

   ENDDO

   DEALLOCATE(tmp1, tmp2, tmat1, tmat2)

END SUBROUTINE invert_trm_real

!-------------------------------------------------------------------------------

SUBROUTINE cholesky_complex(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  cholesky_complex: Cholesky factorization of a complex hermitian matrix
!
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be factorized.
!              Distribution is like in Scalapack.
!              Only upper triangle is needs to be set.
!              On return, the upper triangle contains the Cholesky factor
!              and the lower triangle is set to 0.
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------


    INTEGER                                  :: na, lda
    COMPLEX*16                               :: a(lda,*)
    INTEGER                                  :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols

    COMPLEX*16, ALLOCATABLE                  :: tmatc(:,:), tmatr(:,:), &
                                                tmp1(:), tmp2(:,:)
    INTEGER :: i, info, l_col1, l_cols, l_cols_tile, l_colx, l_row1, l_rows, &
      l_rows_tile, l_rowx, lce, lcs, lre, lrs, mpierr, my_pcol, my_prow, n, &
      nc, np_cols, np_rows, pcol, prow, tile_size

   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

   tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
   tile_size = ((128*MAX(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

   l_rows_tile = tile_size/np_rows ! local rows of a tile
   l_cols_tile = tile_size/np_cols ! local cols of a tile


   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

   ALLOCATE(tmp1(nblk*nblk))
   ALLOCATE(tmp2(nblk,nblk))
   tmp1 = 0
   tmp2 = 0

   ALLOCATE(tmatr(l_rows,nblk))
   ALLOCATE(tmatc(l_cols,nblk))
   tmatr = 0
   tmatc = 0


   DO n = 1, na, nblk

      ! Calculate first local row and column of the still remaining matrix
      ! on the local processor

      l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
      l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

      l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
      l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

      IF(n+nblk > na) THEN

         ! This is the last step, just do a Cholesky-Factorization
         ! of the remaining block

         IF(my_prow==prow(n) .AND. my_pcol==pcol(n)) THEN

            CALL zpotrf('U',na-n+1,a(l_row1,l_col1),lda,info)
            IF(info/=0) THEN
               PRINT *,"Error in zpotrf"
               CALL MPI_Abort(MPI_COMM_WORLD,1,mpierr)
            ENDIF

         ENDIF

         EXIT ! Loop

      ENDIF


      IF(my_prow==prow(n)) THEN

         IF(my_pcol==pcol(n)) THEN

            ! The process owning the upper left remaining block does the
            ! Cholesky-Factorization of this block

            CALL zpotrf('U',nblk,a(l_row1,l_col1),lda,info)
            IF(info/=0) THEN
               PRINT *,"Error in zpotrf"
               CALL MPI_Abort(MPI_COMM_WORLD,1,mpierr)
            ENDIF

            nc = 0
            DO i=1,nblk
               tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
               nc = nc+i
            ENDDO
         ENDIF

         CALL MPI_Bcast(tmp1,nblk*(nblk+1)/2,MPI_DOUBLE_COMPLEX,pcol(n),mpi_comm_cols,mpierr)

         nc = 0
         DO i=1,nblk
            tmp2(1:i,i) = tmp1(nc+1:nc+i)
            nc = nc+i
         ENDDO

         IF(l_cols-l_colx+1>0) &
            CALL ztrsm('L','U','C','N',nblk,l_cols-l_colx+1,(1.d0,0.d0),tmp2,UBOUND(tmp2,1),a(l_row1,l_colx),lda)

      ENDIF

      DO i=1,nblk

         IF(my_prow==prow(n)) tmatc(l_colx:l_cols,i) = CONJG(a(l_row1+i-1,l_colx:l_cols))
         IF(l_cols-l_colx+1>0) &
            CALL MPI_Bcast(tmatc(l_colx,i),l_cols-l_colx+1,MPI_DOUBLE_COMPLEX,prow(n),mpi_comm_rows,mpierr)

      ENDDO

      CALL elpa_transpose_vectors  (tmatc, 2*UBOUND(tmatc,1), mpi_comm_cols, &
                                    tmatr, 2*UBOUND(tmatr,1), mpi_comm_rows, &
                                    2*n-1, 2*na, nblk, 2*nblk)


      DO i=0,(na-1)/tile_size
         lcs = MAX(l_colx,i*l_cols_tile+1)
         lce = MIN(l_cols,(i+1)*l_cols_tile)
         lrs = l_rowx
         lre = MIN(l_rows,(i+1)*l_rows_tile)
         IF(lce<lcs .OR. lre<lrs) CYCLE
         CALL ZGEMM('N','C',lre-lrs+1,lce-lcs+1,nblk,(-1.d0,0.d0), &
                    tmatr(lrs,1),UBOUND(tmatr,1),tmatc(lcs,1),UBOUND(tmatc,1), &
                    (1.d0,0.d0),a(lrs,lcs),lda)
      ENDDO

   ENDDO

   DEALLOCATE(tmp1, tmp2, tmatr, tmatc)

   ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

   DO i=1,na
      IF(my_pcol==pcol(i)) THEN
         ! column i is on local processor
         l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
         l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
         a(l_row1:l_rows,l_col1) = 0
      ENDIF
   ENDDO

END SUBROUTINE cholesky_complex

!-------------------------------------------------------------------------------

SUBROUTINE invert_trm_complex(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  invert_trm_complex: Inverts a upper triangular matrix
!
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be inverted.
!              Distribution is like in Scalapack.
!              Only upper triangle is needs to be set.
!              The lower triangle is not referenced.
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------


    INTEGER                                  :: na, lda
    COMPLEX*16                               :: a(lda,*)
    INTEGER                                  :: nblk, mpi_comm_rows, &
                                                mpi_comm_cols

    COMPLEX*16, ALLOCATABLE                  :: tmat1(:,:), tmat2(:,:), &
                                                tmp1(:), tmp2(:,:)
    INTEGER :: i, info, l_col1, l_cols, l_colx, l_row1, l_rows, l_rowx, &
      mpierr, my_pcol, my_prow, n, nb, nc, np_cols, np_rows, ns, pcol, prow

   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   CALL mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   CALL mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   CALL mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   CALL mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

   ALLOCATE(tmp1(nblk*nblk))
   ALLOCATE(tmp2(nblk,nblk))
   tmp1 = 0
   tmp2 = 0

   ALLOCATE(tmat1(l_rows,nblk))
   ALLOCATE(tmat2(nblk,l_cols))
   tmat1 = 0
   tmat2 = 0


   ns = ((na-1)/nblk)*nblk + 1

   DO n = ns,1,-nblk

      l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
      l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

      nb = nblk
      IF(na-n+1 < nblk) nb = na-n+1

      l_rowx = local_index(n+nb, my_prow, np_rows, nblk, +1)
      l_colx = local_index(n+nb, my_pcol, np_cols, nblk, +1)


      IF(my_prow==prow(n)) THEN

         IF(my_pcol==pcol(n)) THEN

            CALL ZTRTRI('U','N',nb,a(l_row1,l_col1),lda,info)
            IF(info/=0) THEN
               PRINT *,"Error in ZTRTRI"
               CALL MPI_Abort(MPI_COMM_WORLD,1,mpierr)
            ENDIF

            nc = 0
            DO i=1,nb
               tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
               nc = nc+i
            ENDDO
         ENDIF

         CALL MPI_Bcast(tmp1,nb*(nb+1)/2,MPI_DOUBLE_COMPLEX,pcol(n),mpi_comm_cols,mpierr)

         nc = 0
         DO i=1,nb
            tmp2(1:i,i) = tmp1(nc+1:nc+i)
            nc = nc+i
         ENDDO

         IF(l_cols-l_colx+1>0) &
            CALL ZTRMM('L','U','N','N',nb,l_cols-l_colx+1,(1.d0,0.d0),tmp2,UBOUND(tmp2,1),a(l_row1,l_colx),lda)

         IF(l_colx<=l_cols)   tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
         IF(my_pcol==pcol(n)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0

      ENDIF

      IF(l_row1>1) THEN
         IF(my_pcol==pcol(n)) THEN
            tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
            a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
         ENDIF

         DO i=1,nb
            CALL MPI_Bcast(tmat1(1,i),l_row1-1,MPI_DOUBLE_COMPLEX,pcol(n),mpi_comm_cols,mpierr)
         ENDDO
      ENDIF

      IF(l_cols-l_col1+1>0) &
         CALL MPI_Bcast(tmat2(1,l_col1),(l_cols-l_col1+1)*nblk,MPI_DOUBLE_COMPLEX,prow(n),mpi_comm_rows,mpierr)

      IF(l_row1>1 .AND. l_cols-l_col1+1>0) &
         CALL ZGEMM('N','N',l_row1-1,l_cols-l_col1+1,nb, (-1.d0,0.d0), &
                    tmat1,UBOUND(tmat1,1),tmat2(1,l_col1),UBOUND(tmat2,1), &
                    (1.d0,0.d0), a(1,l_col1),lda)

   ENDDO

   DEALLOCATE(tmp1, tmp2, tmat1, tmat2)

END SUBROUTINE invert_trm_complex

! --------------------------------------------------------------------------------------------------

INTEGER FUNCTION least_common_multiple(a, b)

   ! Returns the least common multiple of a and b
   ! There may be more efficient ways to do this, we use the most simple approach

   IMPLICIT NONE
   INTEGER, INTENT(in) :: a, b

   DO least_common_multiple = a, a*(b-1), a
      IF(MOD(least_common_multiple,b)==0) EXIT
   ENDDO
   ! if the loop is left regularly, least_common_multiple = a*b

END FUNCTION

! --------------------------------------------------------------------------------------------------

SUBROUTINE hh_transform_real(alpha, xnorm_sq, xf, tau)

   ! Similar to LAPACK routine DLARFP, but uses ||x||**2 instead of x(:)
   ! and returns the factor xf by which x has to be scaled.
   ! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
   ! since this would be expensive for the parallel implementation.

    REAL*8, INTENT(inout)                    :: alpha
    REAL*8, INTENT(in)                       :: xnorm_sq
    REAL*8, INTENT(out)                      :: xf, tau

    REAL*8                                   :: BETA

   IF( XNORM_SQ==0. ) THEN

      IF( ALPHA>=0. ) THEN
         TAU = 0.
      ELSE
         TAU = 2.
         ALPHA = -ALPHA
      ENDIF
      XF = 0.

   ELSE

      BETA = SIGN( SQRT( ALPHA**2 + XNORM_SQ ), ALPHA )
      ALPHA = ALPHA + BETA
      IF( BETA<0 ) THEN
         BETA = -BETA
         TAU = -ALPHA / BETA
      ELSE
         ALPHA = XNORM_SQ / ALPHA
         TAU = ALPHA / BETA
         ALPHA = -ALPHA
      END IF
      XF = 1./ALPHA
      ALPHA = BETA

   ENDIF

END SUBROUTINE


! --------------------------------------------------------------------------------------------------

SUBROUTINE hh_transform_complex(alpha, xnorm_sq, xf, tau)

   ! Similar to LAPACK routine ZLARFP, but uses ||x||**2 instead of x(:)
   ! and returns the factor xf by which x has to be scaled.
   ! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
   ! since this would be expensive for the parallel implementation.

    COMPLEX*16, INTENT(inout)                :: alpha
    REAL*8, INTENT(in)                       :: xnorm_sq
    COMPLEX*16, INTENT(out)                  :: xf, tau

    REAL*8                                   :: ALPHI, ALPHR, BETA

   ALPHR = DBLE( ALPHA )
   ALPHI = DIMAG( ALPHA )

   IF( XNORM_SQ==0. .AND. ALPHI==0. ) THEN

      IF( ALPHR>=0. ) THEN
         TAU = 0.
      ELSE
         TAU = 2.
         ALPHA = -ALPHA
      ENDIF
      XF = 0.

   ELSE

      BETA = SIGN( SQRT( ALPHR**2 + ALPHI**2 + XNORM_SQ ), ALPHR )
      ALPHA = ALPHA + BETA
      IF( BETA<0 ) THEN
         BETA = -BETA
         TAU = -ALPHA / BETA
      ELSE
         ALPHR = ALPHI * (ALPHI/DBLE( ALPHA ))
         ALPHR = ALPHR + XNORM_SQ/DBLE( ALPHA )
         TAU = DCMPLX( ALPHR/BETA, -ALPHI/BETA )
         ALPHA = DCMPLX( -ALPHR, ALPHI )
      END IF
      XF = 1./ALPHA
      ALPHA = BETA

   ENDIF

END SUBROUTINE

! --------------------------------------------------------------------------------------------------

#endif 
END MODULE elpa1

! --------------------------------------------------------------------------------------------------
! Please note that the following routines are outside of the module ELPA1
! so that they can be used with real or complex data
! --------------------------------------------------------------------------------------------------
#if defined(__parallel)

SUBROUTINE elpa_transpose_vectors(vmat_s,ld_s,comm_s,vmat_t,ld_t,comm_t,nvs,nvr,nvc,nblk)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------

   USE elpa1 ! for least_common_multiple
    INTEGER, INTENT(in)                      :: ld_s, comm_s, ld_t, comm_t, &
                                                nvs, nvr, nvc
    REAL*8, INTENT(inout)                    :: vmat_t(ld_t,nvc)
    REAL*8, INTENT(in)                       :: vmat_s(ld_s,nvc)
    INTEGER, INTENT(in)                      :: nblk

    INTEGER                                  :: i, ips, ipt, k, lc, lcm_s_t, &
                                                mpierr, myps, mypt, n, &
                                                nblks_comm, nblks_skip, &
                                                nblks_tot, nl, nps, npt, ns
    REAL*8, ALLOCATABLE                      :: aux(:)

! for least_common_multiple
! for least_common_multiple
! for least_common_multiple
! for least_common_multiple
! for least_common_multiple

   CALL mpi_comm_rank(comm_s,myps,mpierr)
   CALL mpi_comm_size(comm_s,nps ,mpierr)
   CALL mpi_comm_rank(comm_t,mypt,mpierr)
   CALL mpi_comm_size(comm_t,npt ,mpierr)

   ! The basic idea of this routine is that for every block (in the block cyclic
   ! distribution), the processor within comm_t which owns the diagonal
   ! broadcasts its values of vmat_s to all processors within comm_t.
   ! Of course this has not to be done for every block separately, since
   ! the communictation pattern repeats in the global matrix after
   ! the least common multiple of (nps,npt) blocks

   lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

   nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

   ! Get the number of blocks to be skipped at the begin.
   ! This must be a multiple of lcm_s_t (else it is getting complicated),
   ! thus some elements before nvs will be accessed/set.

   nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

   ALLOCATE(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ))

   DO n = 0, lcm_s_t-1

      ips = MOD(n,nps)
      ipt = MOD(n,npt)

      IF(mypt == ipt) THEN

         nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
         IF(nblks_comm==0) CYCLE

         IF(myps == ips) THEN
            k = 0
            DO lc=1,nvc
               DO i = nblks_skip+n, nblks_tot-1, lcm_s_t
                  ns = (i/nps)*nblk ! local start of block i
                  nl = MIN(nvr-i*nblk,nblk) ! length
                  aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
                  k = k+nblk
               ENDDO
            ENDDO
         ENDIF

         CALL MPI_Bcast(aux,nblks_comm*nblk*nvc,MPI_REAL8,ips,comm_s,mpierr)

         k = 0
         DO lc=1,nvc
            DO i = nblks_skip+n, nblks_tot-1, lcm_s_t
               ns = (i/npt)*nblk ! local start of block i
               nl = MIN(nvr-i*nblk,nblk) ! length
               vmat_t(ns+1:ns+nl,lc) = aux(k+1:k+nl)
               k = k+nblk
            ENDDO
         ENDDO

      ENDIF

   ENDDO

   DEALLOCATE(aux)

END SUBROUTINE

!-------------------------------------------------------------------------------

SUBROUTINE elpa_reduce_add_vectors(vmat_s,ld_s,comm_s,vmat_t,ld_t,comm_t,nvr,nvc,nblk)

!-------------------------------------------------------------------------------
! This routine does a reduce of all vectors in vmat_s over the communicator comm_t.
! The result of the reduce is gathered on the processors owning the diagonal
! and added to the array of vectors vmat_t (which is distributed over comm_t).
!
! Opposed to elpa_transpose_vectors, there is NO identical copy of vmat_s
! in the different members within vmat_t (else a reduce wouldn't be necessary).
! After this routine, an allreduce of vmat_t has to be done.
!
! vmat_s    array of vectors to be reduced and added
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors to which vmat_s is added
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------

   USE elpa1 ! for least_common_multiple
    INTEGER, INTENT(in)                      :: ld_s, comm_s, ld_t, comm_t, &
                                                nvr, nvc
    REAL*8, INTENT(inout)                    :: vmat_t(ld_t,nvc)
    REAL*8, INTENT(in)                       :: vmat_s(ld_s,nvc)
    INTEGER, INTENT(in)                      :: nblk

    INTEGER                                  :: i, ips, ipt, k, lc, lcm_s_t, &
                                                mpierr, myps, mypt, n, &
                                                nblks_tot, nl, nps, npt, ns
    REAL*8, ALLOCATABLE                      :: aux1(:), aux2(:)

! for least_common_multiple
! for least_common_multiple
! for least_common_multiple
! for least_common_multiple
! for least_common_multiple

   CALL mpi_comm_rank(comm_s,myps,mpierr)
   CALL mpi_comm_size(comm_s,nps ,mpierr)
   CALL mpi_comm_rank(comm_t,mypt,mpierr)
   CALL mpi_comm_size(comm_t,npt ,mpierr)

   ! Look to elpa_transpose_vectors for the basic idea!

   ! The communictation pattern repeats in the global matrix after
   ! the least common multiple of (nps,npt) blocks

   lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

   nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

   ALLOCATE(aux1( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ))
   ALLOCATE(aux2( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ))
   aux1(:) = 0
   aux2(:) = 0

   DO n = 0, lcm_s_t-1

      ips = MOD(n,nps)
      ipt = MOD(n,npt)

      IF(myps == ips) THEN

         k = 0
         DO lc=1,nvc
            DO i = n, nblks_tot-1, lcm_s_t
               ns = (i/nps)*nblk ! local start of block i
               nl = MIN(nvr-i*nblk,nblk) ! length
               aux1(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
               k = k+nblk
            ENDDO
         ENDDO

         IF(k>0) CALL mpi_reduce(aux1,aux2,k,MPI_REAL8,MPI_SUM,ipt,comm_t,mpierr)

         IF(mypt == ipt) THEN
            k = 0
            DO lc=1,nvc
               DO i = n, nblks_tot-1, lcm_s_t
                  ns = (i/npt)*nblk ! local start of block i
                  nl = MIN(nvr-i*nblk,nblk) ! length
                  vmat_t(ns+1:ns+nl,lc) = vmat_t(ns+1:ns+nl,lc) + aux2(k+1:k+nl)
                  k = k+nblk
               ENDDO
            ENDDO
         ENDIF

      ENDIF

   ENDDO

   DEALLOCATE(aux1)
   DEALLOCATE(aux2)
END SUBROUTINE
#endif
!-------------------------------------------------------------------------------
