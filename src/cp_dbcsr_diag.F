!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Interface to (sca)lapack for the Cholesky based procedures
!> \author  VW
!> \date    2009-11-09
!> \version 0.8
!>
!> <b>Modification history:</b>
!> - Created 2009-11-09
! *****************************************************************************
MODULE cp_dbcsr_diag

  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_cfm_types,                    ONLY: cp_cfm_release,&
                                             cp_cfm_to_cfm,&
                                             cp_cfm_type
  USE cp_dbcsr_operations,             ONLY: copy_cfm_to_dbcsr,&
                                             copy_dbcsr_to_cfm,&
                                             copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr
  USE cp_fm_types,                     ONLY: cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE dbcsr_methods,                   ONLY: dbcsr_get_info
  USE dbcsr_types,                     ONLY: dbcsr_obj
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_diag'

  ! Public subroutines

  PUBLIC :: dbcsr_syevd, &
            dbcsr_syevr, &
            dbcsr_syevx

CONTAINS

! *****************************************************************************
  SUBROUTINE dbcsr_syevd(matrix,eigenvectors,eigenvalues,error)

    ! Computes all eigenvalues and vectors of a real symmetric matrix
    ! should be quite a bit faster than syevx for that case
    ! especially in parallel with thightly clustered evals
    ! needs more workspace in the worst case, but much better distributed

    TYPE(dbcsr_obj), POINTER                 :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_syevd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info, istat, liwork, &
                                                lwork, mypcol, myprow, n, nmo
    INTEGER, DIMENSION(9)                    :: descm, descv
    INTEGER, DIMENSION(:), POINTER           :: iwork
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eig
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m, v
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_type), POINTER                :: fm_eigenvectors, fm_matrix

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CALL copy_dbcsr_to_fm(matrix, fm_matrix)
    CALL copy_dbcsr_to_fm(eigenvectors, fm_eigenvectors)

    context => fm_matrix%matrix_struct%context
    myprow = context%mepos(1)
    mypcol = context%mepos(2)
    n = fm_matrix%matrix_struct%nrow_global
    nmo = SIZE(eigenvalues,1)
    ALLOCATE(eig(n), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    m => fm_matrix%local_data

#if defined(__SCALAPACK)
    descm(:) = fm_matrix%matrix_struct%descriptor(:)
    v => fm_eigenvectors%local_data
    descv(:) = fm_eigenvectors%matrix_struct%descriptor(:)

    liwork = 7*n + 8*context%num_pe(2) + 2
    ALLOCATE(iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "iwork",liwork*int_size)
    ! work space query

    lwork = -1
    ALLOCATE(work(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work",dp_size)

    CALL pdsyevd('V','U',n,m(1,1),1,1,descm,eig(1),v(1,1),1,1,descv, &
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

    CALL pdsyevd('V','U',n,m(1,1),1,1,descm,eig(1),v(1,1),1,1,descv, &
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
    CALL dsyevd('V','U', n, m(1,1), n, eig(1), work(1), lwork, iwork(1), liwork, info)
    CALL cp_fm_to_fm(fm_matrix,fm_eigenvectors,error=error)
#endif


    IF (info /= 0) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Matrix diagonalization failed")
    END IF

    IF (SIZE(eigenvalues,1) > n) THEN
      eigenvalues(1:n) = eig(1:n)
    ELSE
      eigenvalues(1:nmo) = eig(1:nmo)
    END IF

    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "iwork")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work")

    DEALLOCATE (eig,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "eig")

    CALL copy_fm_to_dbcsr(fm_eigenvectors, eigenvectors)

    CALL cp_fm_release(fm_matrix, error)

    CALL cp_fm_release(fm_eigenvectors, error)

    CALL timestop(handle)

  END SUBROUTINE dbcsr_syevd

! *****************************************************************************
!> \brief   compute eigenvalues and optionally eigenvectors of a real symmetric matrix using scalapack.
!>          If eigenvectors are required this routine will replicate a full matrix on each CPU...
!>          if more than a handful of vectors are needed, use dbcsr_syevd instead
!> \par     matrix is supposed to be in upper triangular form, and overwritten by this routine
!>          neig   is the number of vectors needed (default all)
!>          work_syevx evec calculation only, is the fraction of the working buffer allowed (1.0 use full buffer)
!>                     reducing this saves time, but might cause the routine to fail
! *****************************************************************************
  SUBROUTINE dbcsr_syevx(matrix,eigenvectors,eigenvalues,neig,work_syevx,error)

    ! Diagonalise the symmetric n by n matrix using the LAPACK library.

    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(dbcsr_obj), POINTER, OPTIONAL       :: eigenvectors
    REAL(KIND = dp), OPTIONAL, INTENT(IN)        :: work_syevx
    INTEGER, INTENT(IN), OPTIONAL                :: neig
    REAL(KIND = dp), DIMENSION(:), INTENT(OUT)   :: eigenvalues
    TYPE(cp_error_type), INTENT(inout)  :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "dbcsr_syevx",&
                                   routineP = moduleN//":"//routineN
    REAL(KIND = dp), PARAMETER  :: orfac = -1.0_dp,&
                                   vl = 0.0_dp,&
                                   vu = 0.0_dp

    REAL(KIND = dp) :: abstol, work_syevx_local
    INTEGER  :: handle,info,istat,liwork,lwork,m,mypcol,myprow,n,nb,&
                nn,np0,npcol,npe,nprow,nq0,nz,output_unit,itype, neig_local
    LOGICAL  :: ionode, needs_evecs

    INTEGER, DIMENSION(9) :: desca,descz,descb

    TYPE(cp_blacs_env_type), POINTER           :: context
    TYPE(cp_logger_type), POINTER              :: logger
    REAL(KIND = dp), DIMENSION(:), ALLOCATABLE :: gap,w,work
    INTEGER, DIMENSION(:), ALLOCATABLE         :: iclustr,ifail,iwork
    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: a,z,b
    CHARACTER(LEN=1)                           :: job_type
    TYPE(cp_fm_type), POINTER                  :: fm_matrix
    TYPE(cp_fm_type), POINTER                  :: fm_eigenvectors

    REAL(KIND = dp), EXTERNAL :: dlamch

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL  :: iceil,numroc
#else
    INTEGER, EXTERNAL  :: ilaenv
#endif

    CALL timeset(routineN,handle)

    ! by default all
    CALL dbcsr_get_info(matrix,nblkrows_total=n)
    neig_local=n
    IF (PRESENT(neig)) neig_local=neig
    IF (neig_local == 0) RETURN

    NULLIFY(eigenvectors)
    CALL copy_dbcsr_to_fm(matrix, fm_matrix)
    IF(PRESENT(eigenvectors)) CALL copy_dbcsr_to_fm(eigenvectors, fm_eigenvectors)

    needs_evecs=PRESENT(eigenvectors)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    n = fm_matrix%matrix_struct%nrow_global

    ! by default allocate all needed space
    work_syevx_local=1.0_dp
    IF (PRESENT(work_syevx)) work_syevx_local=work_syevx

    ! set scalapack job type
    IF (needs_evecs) THEN
       job_type="V"
    ELSE
       job_type="N"
    ENDIF

    ! target the most accurate calculation of the eigenvalues
    abstol = 2.0_dp*dlamch("S")


    context =>  fm_matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)


    ALLOCATE (w(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "w",n*dp_size)
    eigenvalues(:) = 0.0_dp
#if defined(__SCALAPACK)

    IF (fm_matrix%matrix_struct%nrow_block /= fm_matrix%matrix_struct%ncol_block) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Invalid blocksize (no square blocks)")
    END IF

    a => fm_matrix%local_data
    desca(:) = fm_matrix%matrix_struct%descriptor(:)

    IF (needs_evecs) THEN
      z => fm_eigenvectors%local_data
      descz(:) = fm_eigenvectors%matrix_struct%descriptor(:)
    ELSE
      ! z will not be referenced
      z => fm_matrix%local_data
      descz=desca
    ENDIF

    ! Get the optimal work storage size

    npe = nprow*npcol
    nb = fm_matrix%matrix_struct%nrow_block
    nn = MAX(n,nb,2)
    np0 = numroc(nn,nb,0,0,nprow)
    nq0 = MAX(numroc(nn,nb,0,0,npcol),nb)

    IF (needs_evecs) THEN
       lwork = 5*n + MAX(5*nn,np0*nq0) + iceil(neig_local,npe)*nn + 2*nb*nb +&
               INT(work_syevx_local*REAL((neig_local - 1)*n,dp)) !!!! allocates a full matrix on every CPU !!!!!
    ELSE
       lwork = 5*n + MAX(5*nn,nb*(np0+1))
    ENDIF
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

    CALL pdsyevx(job_type,"I","U",n,a(1,1),1,1,desca,vl,vu,1,neig_local,abstol,m,nz,w(1),orfac,&
                 z(1,1),1,1,descz,work(1),lwork,iwork(1),liwork,ifail(1),iclustr(1),gap,info)

    ! Error handling

    IF (info /= 0) THEN
      IF (ionode) THEN
        output_unit = cp_logger_get_unit_nr(logger,local=.FALSE.)
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

    a => fm_matrix%local_data
    IF (needs_evecs) THEN
       z => fm_eigenvectors%local_data
    ELSE
       ! z will not be referenced
       z => fm_matrix%local_data
    ENDIF

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

    CALL dsyevx(job_type,"I","U",n,a(1,1),n,vl,vu,1,neig_local,abstol,m,w,z(1,1),n,work(1),lwork,&
                iwork(1),ifail(1),info)

    ! Error handling

    IF (info /= 0) THEN
      output_unit = cp_logger_get_unit_nr(logger,local=.FALSE.)
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
    eigenvalues(1:neig_local) = w(1:neig_local)
    DEALLOCATE (w,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "w")

    IF(PRESENT(eigenvectors)) CALL copy_fm_to_dbcsr(fm_eigenvectors, eigenvectors)

    IF(PRESENT(eigenvectors)) CALL cp_fm_release(fm_eigenvectors, error)

    CALL cp_fm_release(fm_matrix, error)

    CALL timestop(handle)

  END SUBROUTINE dbcsr_syevx

! *****************************************************************************
!> \brief  computes selected eigenvalues and, optionally, eigenvectors of 
!>        a real symmetric matrix A distributed in 2D blockcyclic format by 
!>       calling the recommended sequence of ScaLAPACK routines. 
!>          
!> \par     matrix is supposed to be in upper triangular form, and overwritten by this routine
!>          subsets of eigenvalues/vectors can be selected by
!>          specifying a range of values or a range of indices for the desired eigenvalues.
! *****************************************************************************
  SUBROUTINE dbcsr_syevr(matrix,eigenvectors,eigenvalues,ilow,iup,error)

    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(dbcsr_obj), POINTER, OPTIONAL       :: eigenvectors
    REAL(KIND = dp), DIMENSION(:), INTENT(OUT)   :: eigenvalues
    INTEGER, INTENT(IN), OPTIONAL                :: ilow,iup
    TYPE(cp_error_type), INTENT(inout)  :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "dbcsr_syevr",&
                                   routineP = moduleN//":"//routineN
    REAL(KIND = dp), PARAMETER  :: orfac = -1.0_dp,&
                                   vl = 0.0_dp,&
                                   vu = 0.0_dp

    CHARACTER(LEN=1)                           :: job_type
    INTEGER :: handle, ilow_local, info, istat, iup_local, lwork, liwork, m, mypcol,myprow, n, neig, nb, nz
    INTEGER, DIMENSION(9) :: desca,descz,descb
    INTEGER, DIMENSION(:), ALLOCATABLE  :: ifail, iwork
    LOGICAL :: failure, ionode, needs_evecs
    REAL(dp) :: abstol
    REAL(dp), DIMENSION(:), ALLOCATABLE :: w, work
    REAL(dp), DIMENSION(:,:), POINTER :: a, z

    TYPE(cp_blacs_env_type), POINTER           :: context
    TYPE(cp_logger_type), POINTER              :: logger
    TYPE(cp_fm_type), POINTER                  :: fm_matrix
    TYPE(cp_fm_type), POINTER                  :: fm_eigenvectors

    REAL(KIND = dp), EXTERNAL :: dlamch

#if defined(__SCALAPACK)
#else
    INTEGER, EXTERNAL  :: ilaenv
#endif

    ! -------------------------------------------------------------------------

    ! by default all
    CALL dbcsr_get_info(matrix,nblkrows_total=n)
    neig=n
    iup_local = n
    ilow_local = 1
    IF (PRESENT(ilow) .AND. PRESENT(iup)) THEN
        neig=iup-ilow+1
        iup_local = iup
        ilow_local = ilow
    END IF
    IF (neig <= 0) RETURN

    CALL timeset(routineN,handle)

    failure =.FALSE.
    needs_evecs=PRESENT(eigenvectors)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    n = fm_matrix%matrix_struct%nrow_global

    ! set scalapack job type
    IF (needs_evecs) THEN
       job_type="V"
    ELSE
       job_type="N"
    ENDIF

    context =>  fm_matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)

    ALLOCATE(w(n), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    eigenvalues(:) = 0.0_dp

#if defined(__SCALAPACK)

    IF (fm_matrix%matrix_struct%nrow_block /= fm_matrix%matrix_struct%ncol_block) THEN
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

    a => fm_matrix%local_data
    desca(:) = fm_matrix%matrix_struct%descriptor(:)

    IF (needs_evecs) THEN
      z => fm_eigenvectors%local_data
      descz(:) = fm_eigenvectors%matrix_struct%descriptor(:)
    ELSE
      ! z will not be referenced
      z => fm_matrix%local_data
      descz=desca
    ENDIF

    ! First Call: Determine the needed work_space
    lwork = -1
    ALLOCATE(work(5*n), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(iwork(6*n), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
!    CALL pdsyevr(job_type,'I','U',n,a,1,1,desca,vl,vu,ilow_local,iup_local,m,nz,w(1),z,1,1,descz,work,lwork,iwork,liwork,info)
    lwork = INT(work(1))
    liwork = iwork(1)
    IF(lwork>SIZE(work,1)) THEN
      DEALLOCATE(work,STAT=istat)
      ALLOCATE(work(lwork),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(liwork>SIZE(iwork,1)) THEN
      DEALLOCATE(iwork,STAT=istat)
      ALLOCATE(iwork(liwork),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    !Second call: solve the eigenvalue problem
    info = 0
!    CALL pdsyevr(job_type,'I','U',n,a,1,1,desca,vl,vu,ilow_local,iup_local,m,nz,w(1),z,1,1,descz,work,lwork,iwork,liwork,info)

    IF(info>0) THEN
      WRITE(*,*) 'Processor ', myprow, mypcol, ': Error! INFO code = ', INFO 
    END IF
    CPPrecondition(info==0,cp_failure_level,routineP,error,failure)

    ! Release work storage
    DEALLOCATE (iwork,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (work,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

#else

    a => fm_matrix%local_data
    IF (needs_evecs) THEN
       z => fm_eigenvectors%local_data
    ELSE
       ! z will not be referenced
       z => fm_matrix%local_data
    ENDIF

    ! Get the optimal work storage size

    nb = MAX(ilaenv(1,"DSYTRD","U",n,-1,-1,-1),&
             ilaenv(1,"DORMTR","U",n,-1,-1,-1))

    lwork = MAX((nb + 3)*n,8*n)+n ! sun bug fix
    liwork = 5*n

    ALLOCATE (ifail(n),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ifail = 0

    ALLOCATE (iwork(liwork),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (work(lwork),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    ! target the most accurate calculation of the eigenvalues
    abstol = 2.0_dp*dlamch("S")

    info = 0
    CALL dsyevx(job_type,"I","U",n,a(1,1),n,vl,vu,ilow_local,iup_local,abstol,m,w,z(1,1),n,work(1),lwork,&
                iwork(1),ifail(1),info)

    ! Error handling
    CPPrecondition(info==0,cp_failure_level,routineP,error,failure)

    ! Release work storage
    DEALLOCATE (iwork,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (work,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

#endif

    eigenvalues(ilow_local:iup_local) = w(ilow_local:iup_local)
    DEALLOCATE (w,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    IF(PRESENT(eigenvectors)) CALL copy_fm_to_dbcsr(fm_eigenvectors, eigenvectors)

    IF(PRESENT(eigenvectors)) CALL cp_fm_release(fm_eigenvectors, error)

    CALL cp_fm_release(fm_matrix, error)

    CALL timestop(handle)

  END SUBROUTINE dbcsr_syevr

! *****************************************************************************
  SUBROUTINE dbcsr_heevd(matrix,eigenvectors,eigenvalues,error)

    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(dbcsr_obj), OPTIONAL, POINTER       :: eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_heevd', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: work
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: m, v
    INTEGER                                  :: handle, info, istat, liwork, &
                                                lrwork, lwork, n
    INTEGER, DIMENSION(9)                    :: descm, descv
    INTEGER, DIMENSION(:), POINTER           :: iwork
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rwork
    TYPE(cp_cfm_type), POINTER               :: fm_eigenvectors, fm_matrix

    CALL timeset(routineN,handle)

    failure=.FALSE.

    CALL copy_dbcsr_to_cfm(matrix, fm_matrix)
    CALL copy_dbcsr_to_cfm(eigenvectors, fm_eigenvectors)

    n = fm_matrix%matrix_struct%nrow_global
    m => fm_matrix%local_data
    ALLOCATE(iwork(1),work(1),rwork(1),STAT=istat)
    IF (istat.NE.0) CALL stop_memory(routineN,"iwork")
! work space query
    lwork  = -1
    lrwork = -1
    liwork = -1

#if defined(__SCALAPACK)
    v => fm_eigenvectors%local_data
    descm(:) = fm_matrix%matrix_struct%descriptor(:)
    descv(:) = fm_eigenvectors%matrix_struct%descriptor(:)
    CALL PZHEEVD('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv, &
                  work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
    lwork=work(1)
    lrwork=rwork(1)+1000000 ! needed to correct for a bug in scalapack, unclear how much the right number is
    liwork=iwork(1)
    DEALLOCATE(work,iwork,rwork)
    ALLOCATE(iwork(liwork),work(lwork),rwork(lrwork),STAT=istat)
    IF (istat.NE.0) CALL stop_memory(routineN,"work")
    iwork=0;work=CMPLX(0._dp,0._dp,dp);rwork=0._dp;
    CALL PZHEEVD('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv, &
                  work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
#else
    CALL ZHEEVD('V','U',n,m(1,1),SIZE(m,1),eigenvalues(1), &
                  work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
    lwork=work(1)
    lrwork=rwork(1)
    liwork=iwork(1)
    DEALLOCATE(work,iwork,rwork)
    ALLOCATE(iwork(liwork),work(lwork),rwork(lrwork),STAT=istat)
    IF (istat.NE.0) CALL stop_memory(routineN,"work")
    CALL zheevd('V','U',n,m(1,1),SIZE(m,1),eigenvalues(1), &
                  work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
    CALL cp_cfm_to_cfm(fm_matrix,fm_eigenvectors,error=error)
#endif

    IF (info.NE.0) CALL stop_program(routineN,"unable to diagonalize matrix")
    DEALLOCATE(work,iwork,rwork)


    CALL copy_cfm_to_dbcsr(fm_eigenvectors, eigenvectors)

    CALL cp_cfm_release(fm_matrix, error)

    CALL cp_cfm_release(fm_eigenvectors, error)

    CALL timestop(handle)

  END SUBROUTINE dbcsr_heevd


END MODULE cp_dbcsr_diag
