!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006 CP2K developers group                           !
!-----------------------------------------------------------------------------!
!!****h* cp2k/cp_fm_diag [1.0] *
!!
!!   NAME
!!     cp_fm_diag
!!
!!   FUNCTION
!!     used for collecting some of the diagonalization shemes available for
!!     cp_fm_type. cp_fm_power also moved here as it is very related
!!
!!   NOTES
!!     first version : most routines imported
!!
!!   AUTHOR
!!     Joost VandeVondele (2003-08)
!!
!!   MODIFICATION HISTORY
!!     - unused Jacobi routines removed, cosmetics (05.04.06,MK)
!!
!!******************************************************************************

MODULE cp_fm_diag

  USE cp_blacs_calls,                  ONLY: cp_blacs_gridexit,&
                                             cp_blacs_gridinit
  USE cp_fm_basic_linalg,              ONLY: cp_fm_syrk,&
                                             cp_fm_upper_to_full
  USE cp_fm_types,                     ONLY: cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_bcast
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: moduleN = "cp_fm_diag"

  ! Public subroutines

  PUBLIC :: cp_fm_block_jacobi,&
            cp_fm_power, &
            cp_fm_syevd, &
            cp_fm_syevx

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE cp_fm_syevd(matrix,eigenvectors,eigenvalues,error)

    ! Computes all eigenvalues and vectors of a real symmetric matrix
    ! should be quite a bit faster than syevx for that case
    ! especially in parallel with thightly clustered evals
    ! needs more workspace in the worst case, but much better distributed

    TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "cp_fm_syevd"

    INTEGER                                  :: handle, info, istat, liwork, &
                                                lwork, mypcol, myprow, n
    INTEGER, DIMENSION(9)                    :: descm, descv
    INTEGER, DIMENSION(:), POINTER           :: iwork
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m, v
    TYPE(cp_blacs_env_type), POINTER         :: context

    ! -------------------------------------------------------------------------

    CALL timeset("cp_fm_syevd","I","",handle)

    failure = .FALSE.

    context => matrix%matrix_struct%context
    myprow = context%mepos(1)
    mypcol = context%mepos(2)
    n = matrix%matrix_struct%nrow_global
    IF (SIZE(eigenvalues,1) < n) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Incorrect size for eigenvalue array")
    END IF

    m => matrix%local_data

#if defined(__SCALAPACK)
    descm(:) = matrix%matrix_struct%descriptor(:)
    v => eigenvectors%local_data
    descv(:) = eigenvectors%matrix_struct%descriptor(:)

    liwork = 7*n + 8*context%num_pe(2) + 2
    ALLOCATE(iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "iwork",liwork*int_size)
    ! work space query

    lwork = -1
    ALLOCATE(work(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work",dp_size)

    CALL pdsyevd('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv, &
                 work(1),lwork,iwork(1),liwork,info)

    ! look here for a PDORMTR warning :-)
    ! this routine seems to need more workspace than reported
    ! (? lapack bug ...)
    ! arbitrary additional memory  ... we give 100000 more words
    ! (it seems to depend on the block size used)

    lwork = NINT(work(1)+100000)
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work")
    ALLOCATE(work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "lwork",lwork*dp_size)

    CALL pdsyevd('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv, &
                 work(1),lwork,iwork(1),liwork,info)

#else

#if defined(__DEC)
    !if using CXML, LWORK must be at least 1 + 5*N + 2*N*lg N + 3*N**2,
    ! where lg( N ) = smallest integer k such that 2**k >= N
    lwork=1+5*n+2*n*20+3*n**2
    liwork = 5*n + 2
#else
    ! for n=4 with my lapack (linux) I need to add 50 to the value of
    ! lwork. 5000 seems to work in IRIX
    ! unfortunately query is not supported. I did not test n dependence [fawzi]
    lwork = 1 + 6*n + 2*n**2 + 5000
    liwork = 5*n + 3
#endif
    ALLOCATE(work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work",lwork*dp_size)
    ALLOCATE(iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "iwork",liwork*int_size)
    CALL dsyevd('V','U', n, m(1,1), n, eigenvalues(1), work(1), lwork, iwork(1), liwork, info)
    CALL cp_fm_to_fm(matrix,eigenvectors)
#endif

    IF (info /= 0) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Matrix diagonalization failed")
    END IF

    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "iwork")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_syevd

  ! ***************************************************************************

  SUBROUTINE cp_fm_syevx(matrix,eigenvectors,eigenvalues,neig,work_syevx,error)

    ! Diagonalise the symmetric n by n matrix using the LAPACK library.
    ! if bmatrix present, use general solver, will overwrite bmatrix
    ! with its cholesky decom.

    ! - Creation (06.06.2001, Matthias Krack)
    ! - Mod      (05.2001, Joost VandeVondele)

    ! *************************************************************************

    TYPE(cp_fm_type), POINTER                    :: eigenvectors,matrix
    REAL(KIND = dp), INTENT(IN)                  :: work_syevx
    INTEGER, INTENT(IN)                          :: neig
    REAL(KIND = dp), DIMENSION(:), INTENT(OUT)   :: eigenvalues
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "cp_fm_syevx",&
                                   routineP = moduleN//":"//routineN
    REAL(KIND = dp), PARAMETER  :: orfac = -1.0_dp,&
                                   vl = 0.0_dp,&
                                   vu = 0.0_dp

    REAL(KIND = dp) :: abstol
    INTEGER  :: handle,info,istat,liwork,lwork,m,mypcol,myprow,n,nb,&
                nn,np0,npcol,npe,nprow,nq0,nz,output_unit,itype
    LOGICAL  :: ionode

    INTEGER, DIMENSION(9) :: desca,descz,descb

    TYPE(cp_blacs_env_type), POINTER           :: context
    TYPE(cp_logger_type), POINTER              :: logger
    REAL(KIND = dp), DIMENSION(:), ALLOCATABLE :: gap,w,work
    INTEGER, DIMENSION(:), ALLOCATABLE         :: iclustr,ifail,iwork
    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: a,z,b

    REAL(KIND = dp), EXTERNAL :: dlamch

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL  :: iceil,numroc
#else
    INTEGER, EXTERNAL  :: ilaenv
#endif
    ! -------------------------------------------------------------------------

    IF (neig == 0) RETURN

    CALL timeset("cp_fm_syevx","I","",handle)

    abstol = 2.0_dp*dlamch("S")

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source

    context =>  matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

    n = matrix%matrix_struct%nrow_global

    ALLOCATE (w(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "w",n*dp_size)
    w(:) = 0.0_dp

    eigenvalues(:) = 0.0_dp
#if defined(__SCALAPACK)

    IF (matrix%matrix_struct%nrow_block /= matrix%matrix_struct%ncol_block) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Invalid blocksize (no square blocks)")
    END IF

    a => matrix%local_data
    desca(:) = matrix%matrix_struct%descriptor(:)
    z => eigenvectors%local_data
    descz(:) = eigenvectors%matrix_struct%descriptor(:)

    ! Get the optimal work storage size

    npe = nprow*npcol
    nb = matrix%matrix_struct%nrow_block
    nn = MAX(n,nb,2)
    np0 = numroc(nn,nb,0,0,nprow)
    nq0 = MAX(numroc(nn,nb,0,0,npcol),nb)

    lwork = 5*n + MAX(5*nn,np0*nq0) + iceil(neig,npe)*nn + 2*nb*nb +&
            INT(work_syevx*REAL((neig - 1)*n,dp)) !!!! allocates a full matrix on every CPU !!!!!
    liwork = 6*MAX(N,npe+1,4)

    ALLOCATE (gap(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "gap",npe*dp_size)
    gap = 0.0_dp
    ALLOCATE (iclustr(2*npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "iclustr",2*npe*int_size)
    iclustr = 0
    ALLOCATE (ifail(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "ifail",n*int_size)
    ifail = 0
    ALLOCATE (iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "iwork",liwork*int_size)
    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work",lwork*dp_size)

    CALL pdsyevx("V","I","U",n,a(1,1),1,1,desca,vl,vu,1,neig,abstol,m,nz,w(1),orfac,&
                 z(1,1),1,1,descz,work(1),lwork,iwork(1),liwork,ifail(1),iclustr(1),gap,info)


    ! Error handling

    IF (info /= 0) THEN
      IF (ionode) THEN
        output_unit = cp_logger_get_unit_nr(logger,&
                                            cp_failure_level,&
                                            routineP,local=.FALSE.)
        WRITE (unit=output_unit,FMT="(/,(T3,A,T12,1X,I10))")&
          "info    = ",info,&
          "lwork   = ",lwork,&
          "liwork  = ",liwork,&
          "nz      = ",nz
        IF (info > 0) THEN
          WRITE (unit=output_unit,FMT="(/,T3,A,(T12,6(1X,I10)))")&
            "ifail   = ",ifail
          WRITE (unit=output_unit,FMT="(/,T3,A,(T12,6(1X,I10)))")&
            "iclustr = ",iclustr
          WRITE (unit=output_unit,FMT="(/,T3,A,(T12,6(1X,E10.3)))")&
            "gap     = ",gap
        END IF
      END IF
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Error in pdsyevx (ScaLAPACK)")
    END IF

    ! Release work storage

    DEALLOCATE (gap,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "gap")
    DEALLOCATE (iclustr,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "iclustr")
    DEALLOCATE (ifail,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "ifail")
    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "iwork")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work")

#else

    a => matrix%local_data
    z => eigenvectors%local_data

    ! Get the optimal work storage size

    nb = MAX(ilaenv(1,"DSYTRD","U",n,-1,-1,-1),&
             ilaenv(1,"DORMTR","U",n,-1,-1,-1))

    lwork = MAX((nb + 3)*n,8*n)+n ! sun bug fix
    liwork = 5*n

    ALLOCATE (ifail(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "ifail",n*int_size)
    ifail = 0
    ALLOCATE (iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "iwork",liwork*int_size)
    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work",lwork*dp_size)
    info = 0

    CALL dsyevx("V","I","U",n,a(1,1),n,vl,vu,1,neig,abstol,m,w,z(1,1),n,work(1),lwork,&
                iwork(1),ifail(1),info)

    ! Error handling

    IF (info /= 0) THEN
      output_unit = cp_logger_get_unit_nr(logger,&
                    cp_failure_level,routineP,local=.FALSE.)
      WRITE (unit=output_unit,FMT="(/,(T3,A,T12,1X,I10))")&
        "info    = ",info
      IF (info > 0) THEN
        WRITE (unit=output_unit,FMT="(/,T3,A,(T12,6(1X,I10)))")&
          "ifail   = ",ifail
      END IF
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Error in dsyevx")
    END IF

    ! Release work storage

    DEALLOCATE (ifail,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "ifail")
    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "iwork")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work")

#endif
    eigenvalues(1:neig) = w(1:neig)
    DEALLOCATE (w,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "w")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_syevx

  ! ***************************************************************************

  SUBROUTINE cp_fm_power(matrix,work,exponent,threshold,n_dependent)

   ! Raise the real symmetric n by n matrix to the power given by
   ! exponent. All eigenvectors with a corresponding eigenvalue lower
   ! than threshold are quenched. result in matrix

   ! - Creation (29.03.1999, Matthias Krack)
   ! - Parallelised using BLACS and ScaLAPACK (06.06.2001,MK)

   ! **************************************************************************

   TYPE(cp_fm_type), POINTER                 :: matrix,work
   REAL(KIND = dp), INTENT(IN)               :: exponent,threshold
   INTEGER, INTENT(OUT)                      :: n_dependent

   CHARACTER(LEN=*), PARAMETER :: routineN = "cp_fm_power"

   REAL(KIND = dp) :: condition_number,f,p
   INTEGER         :: handle,icol_global,icol_local,ipcol,iprow,irow_global,&
                      irow_local,istat,mypcol,myprow,ncol_block,ncol_global,npcol,&
                      nprow,nrow_block,nrow_global

   REAL(KIND = dp), DIMENSION(:), ALLOCATABLE :: eigenvalues
   REAL(KIND = dp), DIMENSION(:,:), POINTER   :: eigenvectors
   TYPE(cp_blacs_env_type), POINTER           :: context

#if defined(__SCALAPACK)
   INTEGER, EXTERNAL :: indxg2l,indxg2p
#endif
    ! -------------------------------------------------------------------------

    CALL timeset("cp_fm_power","I","",handle)

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    n_dependent = 0
    p = 0.5_dp*exponent

    nrow_global = matrix%matrix_struct%nrow_global
    ncol_global = matrix%matrix_struct%ncol_global

    ALLOCATE (eigenvalues(ncol_global),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
                       "eigenvalues",ncol_global*dp_size)
    END IF
    eigenvalues(:) = 0.0_dp

    ! Compute the eigenvectors and eigenvalues

    CALL cp_fm_syevd(matrix,work,eigenvalues)

#if defined(__SCALAPACK)
    nrow_block = work%matrix_struct%nrow_block
    ncol_block = work%matrix_struct%ncol_block

    eigenvectors => work%local_data

    ! Build matrix**exponent with eigenvector quenching

    p = 0.5_dp*exponent

    n_dependent = 0

    DO icol_global=1,ncol_global

      IF (eigenvalues(icol_global) < threshold) THEN

        n_dependent = n_dependent + 1

        ipcol = indxg2p(icol_global,ncol_block,mypcol,&
                        work%matrix_struct%first_p_pos(2),npcol)

        IF (mypcol == ipcol) THEN
          icol_local = indxg2l(icol_global,ncol_block,mypcol,&
                               work%matrix_struct%first_p_pos(2),npcol)
          DO irow_global=1,nrow_global
            iprow = indxg2p(irow_global,nrow_block,myprow,&
                            work%matrix_struct%first_p_pos(1),nprow)
            IF (myprow == iprow) THEN
              irow_local = indxg2l(irow_global,nrow_block,myprow,&
                                   work%matrix_struct%first_p_pos(1),nprow)
              eigenvectors(irow_local,icol_local) = 0.0_dp
            END IF
          END DO
        END IF

      ELSE

        f = eigenvalues(icol_global)**p

        ipcol = indxg2p(icol_global,ncol_block,mypcol,&
                        work%matrix_struct%first_p_pos(2),npcol)

        IF (mypcol == ipcol) THEN
          icol_local = indxg2l(icol_global,ncol_block,mypcol,&
                               work%matrix_struct%first_p_pos(2),npcol)
          DO irow_global=1,nrow_global
            iprow = indxg2p(irow_global,nrow_block,myprow,&
                            work%matrix_struct%first_p_pos(1),nprow)
            IF (myprow == iprow) THEN
              irow_local = indxg2l(irow_global,nrow_block,myprow,&
                                   work%matrix_struct%first_p_pos(1),nprow)
              eigenvectors(irow_local,icol_local) =&
                f*eigenvectors(irow_local,icol_local)
            END IF
          END DO
        END IF

      END IF

    END DO

#else

    eigenvectors => work%local_data

    ! Build matrix**exponent with eigenvector quenching

    DO icol_global=1,ncol_global

      IF (eigenvalues(icol_global) < threshold) THEN

        n_dependent = n_dependent + 1
        eigenvectors(1:nrow_global,icol_global) = 0.0_dp

      ELSE

        f = eigenvalues(icol_global)**p
        eigenvectors(1:nrow_global,icol_global) =&
          f*eigenvectors(1:nrow_global,icol_global)

      END IF

    END DO

#endif
    CALL cp_fm_syrk("U","N",ncol_global,1.0_dp,work,1,1,0.0_dp,matrix)
    CALL cp_fm_upper_to_full(matrix,work)

    ! Print some warnings/notes

    IF (matrix%matrix_struct%para_env%mepos == 0) THEN
      condition_number = ABS(eigenvalues(ncol_global)/eigenvalues(1))
      WRITE (UNIT=cp_logger_get_default_unit_nr(),FMT="(/,(T2,A,ES15.6))")&
        "CP_FM_POWER: smallest eigenvalue:",eigenvalues(1),&
        "CP_FM_POWER: largest eigenvalue: ",eigenvalues(ncol_global),&
        "CP_FM_POWER: condition number:   ",condition_number
      IF (eigenvalues(1) <= 0.0_dp) THEN
        WRITE (UNIT=cp_logger_get_default_unit_nr(),FMT="(/,T2,A)")&
          "WARNING: matrix has a negative eigenvalue, tighten EPS_DEFAULT"
      END IF
      IF (condition_number > 1.0E12_dp) THEN
        WRITE (UNIT=cp_logger_get_default_unit_nr(),FMT="(/,T2,A)")&
          "WARNING: high condition number => possibly ill-conditioned matrix"
      END IF
    END IF

    DEALLOCATE (eigenvalues,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"eigenvalues")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_power

  ! ***************************************************************************

  SUBROUTINE cp_fm_block_jacobi(matrix,eigenvectors,eigval,thresh,&
                                start_sec_block)

    ! Calculates block diagonalizazion from full symmetric matrix
    ! It has its origin in cp_fm_syevx. This routine rotates only elements
    ! which are larger than a threshold thresh.
    ! start_sec_block is the start of the second block.
    ! IT DOES ONLY ONE SWEEP!

    ! - Creation (07.10.2002, Martin Fengler)
    ! - Cosmetics (05.04.06,MK)

    ! *************************************************************************

    TYPE(cp_fm_type), POINTER                 :: eigenvectors,matrix
    REAL(KIND = dp), DIMENSION(:), INTENT(IN) :: eigval
    INTEGER, INTENT(IN)                       :: start_sec_block
    REAL(KIND = dp), INTENT(IN)               :: thresh

    INTEGER :: handle
    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: a,ev

    REAL(KIND = dp) :: tan_theta,tau,c,s
    INTEGER  :: I,J,q,p,bs,N,q_loc
    REAL(KIND = dp), DIMENSION (:),ALLOCATABLE :: c_ip

#if defined(__SCALAPACK)
    TYPE(cp_blacs_env_type), POINTER :: context

    INTEGER :: myprow,mypcol,nprow,npcol,ictxt_loc,block_dim_row,block_dim_col,&
               info,ev_row_block_size,iam,counter,status,ierr
    INTEGER :: source,allgrp,mynumrows,ictxt,temp_row,temp_col,temp_nprow,&
               temp_npcol,mype,npe,a_row,a_col,a_nb,a_mb,iglob,jglob,&
               i_loc,j_loc

    REAL(KIND = dp), DIMENSION(:,:), ALLOCATABLE :: a_loc, ev_loc
    REAL(KIND = dp), DIMENSION(:), ALLOCATABLE   :: diagonale
    INTEGER, DIMENSION(9)                        :: desca,descz,desc_a_block,&
                                                    desc_ev_loc

    INTEGER, EXTERNAL :: numroc,indxl2g,indxg2l
#endif

    ! -------------------------------------------------------------------------

    CALL timeset("cp_fm_block_jacobi","I","",handle)

#if defined(__SCALAPACK)

    ! Umgebung fuer die Kohn-Sham-Matrix und die Eigenvektoren uebernehmen

    context =>  matrix%matrix_struct%context
    source =    matrix%matrix_struct%para_env%source ! globenv%ionode
    allgrp =    matrix%matrix_struct%para_env%group ! globenv%ionode

    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

    N = matrix%matrix_struct%nrow_global

    A => matrix%local_data
    desca(:) = matrix%matrix_struct%descriptor(:)
    EV => eigenvectors%local_data
    descz(:) = eigenvectors%matrix_struct%descriptor(:)

! Kopiere den Block, der wegrotiert werden soll zunaechst auf den Masterprozessor, und anschliessend
! per Broadcast an alle Prozis
! ACHTUNG start_sec_block sagt aus WO der ZWEITE Block STARTET!!!
! Der Block wird mitsamt dem OO-Block bearbeitet

    block_dim_row = start_sec_block-1
    block_dim_col = N - block_dim_row
    ALLOCATE(A_loc(block_dim_row,block_dim_col))

    mype=matrix%matrix_struct%para_env%mepos 
    npe=matrix%matrix_struct%para_env%num_pe 
    ictxt    =matrix%matrix_struct%context%group
    ! get a new context
    ictxt_loc=matrix%matrix_struct%para_env%group
    CALL cp_blacs_gridinit(ictxt_loc,'Row-major',NPROW*NPCOL,1)

    CALL descinit(desc_a_block,block_dim_row,block_dim_col,block_dim_row,&
                  block_dim_col,0,0,ictxt_loc,block_dim_row,info)

    CALL pdgemr2d(block_dim_row,block_dim_col,A,1,start_sec_block,desca,&
                  A_loc,1,1,desc_a_block,ictxt)
    ! Jetzt sind eigentlich nur im Master-Prozess Daten reingekommen
    CALL mp_bcast(A_loc,0,allgrp)

    ! Da nun jeder ueber den oberen Block verfuegt, koennen wir jetzt die Eigenvektoren so umsortieren, dass
    ! man ein NN*1 Prozessgrid hat, und somit jeder Prozessor ueber ein Buendel von Zeilel verfuegt,
    ! die selbst unabhaengig modifizieren darf.

    ! Aufsetzen der Eigenvektorverteilung
    iam = mype
    ev_row_block_size = n/(nprow*npcol)
    mynumrows = NUMROC(N,ev_row_block_size,iam,0,NPROW*NPCOL)

    ALLOCATE(EV_loc(mynumrows,N),c_ip(mynumrows))

    CALL descinit(desc_ev_loc,N,N,ev_row_block_size,N,0,0,ictxt_loc,&
                  mynumrows,info)

    CALL pdgemr2d(N,N,EV,1,1,descz,EV_loc,1,1,desc_ev_loc,ictxt)

!   *** START Diagonalising matrix ***

    ! Eigentliche Blockdiagonalisierung

    q_loc = 0

    DO q=start_sec_block,N
      q_loc = q_loc + 1
      DO p=1,(start_sec_block-1)

          IF(ABS(A_loc(p,q_loc))>thresh) THEN

          tau = (eigval(q)-eigval(p))/(2.0_dp*A_loc(p,q_loc))

          tan_theta = SIGN(1.0_dp,tau)/(ABS(tau)+SQRT(1.0_dp+tau*tau))

          ! cos theta
          c = 1.0_dp/SQRT(1.0_dp+tan_theta*tan_theta)
          s = tan_theta*c

          ! Und jetzt noch die Eigenvektoren produzieren:
          ! Q * J
          !  Verstaendliche Version (bevor die BLAS-Aufrufe sie ersetzt haben)
          !  c_ip = c*EV_loc(:,p) - s*EV_loc(:,q)
          !  c_iq = s*EV_loc(:,p) + c*EV_loc(:,q)

          !  EV(:,p)=c_ip
          !  EV(:,q)=c_iq

          CALL dcopy(mynumrows,EV_loc(1,p),1,c_ip(1),1)
          CALL dscal(mynumrows,c,EV_loc(1,p),1)
          CALL daxpy(mynumrows,-s,EV_loc(1,q),1,EV_loc(1,p),1)
          CALL dscal(mynumrows,c,EV_loc(1,q),1)
          CALL daxpy(mynumrows,s,c_ip(1),1,EV_loc(1,q),1)

         END IF

      END DO
    END DO

    ! Nun muessen die Eigenvektoren wieder in die alte Verteilung zurueckverschickt werden.
    ! Verschicke EVs zurueck!'
    CALL pdgemr2d(N,N,EV_loc,1,1,desc_ev_loc,EV,1,1,descz,ictxt)

    ! Speicher freigeben
    DEALLOCATE(A_loc,EV_loc,c_ip)

    CALL cp_blacs_gridexit(ictxt_loc)

#else

    N = matrix%matrix_struct%nrow_global ! Groesse der Matrix A, die bearbeitet werden soll

    ALLOCATE(c_ip(N))   ! Speicher fuer den lokalen Eigenwertvektor

    A => matrix%local_data        ! Contains the Matrix to be worked on
    EV => eigenvectors%local_data ! Contains the eigenvectors up to blocksize: Rest ist Muell

    ! Start diagonalizing matrix

    tan_theta = 0.0_dp
    tau       = 0.0_dp

    DO q=start_sec_block,N
      DO p=1,(start_sec_block-1)

          IF(ABS(A(p,q))>thresh) THEN

          tau = (eigval(q)-eigval(p))/(2.0_dp*A(p,q))

          tan_theta = SIGN(1.0_dp,tau)/(ABS(tau)+SQRT(1.0_dp+tau*tau))

          ! cos theta
          c = 1.0_dp/SQRT(1.0_dp+tan_theta*tan_theta)
          s = tan_theta*c

          ! Und jetzt noch die Eigenvektoren produzieren:
          ! Q * J
          !  Verstaendliche Version (bevor die BLAS-Aufrufe sie ersetzt haben)
          !  c_ip = c*EV(:,p) - s*EV(:,q)
          !  c_iq = s*EV(:,p) + c*EV(:,q)

          !  EV(:,p)=c_ip
          !  EV(:,q)=c_iq

          CALL dcopy(N,EV(1,p),1,c_ip(1),1)
          CALL dscal(N,c,EV(1,p),1)
          CALL daxpy(N,-s,EV(1,q),1,EV(1,p),1)
          CALL dscal(N,c,EV(1,q),1)
          CALL daxpy(N,s,c_ip(1),1,EV(1,q),1)

         END IF

      END DO
    END DO

    ! Release work storage

    DEALLOCATE (c_ip)

#endif

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_block_jacobi

  ! ***************************************************************************

END MODULE cp_fm_diag
