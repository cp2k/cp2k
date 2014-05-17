!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief used for collecting some of the diagonalization shemes available for
!>      cp_fm_type. cp_fm_power also moved here as it is very related
!> \note
!>      first version : most routines imported
!> \par History
!>      - unused Jacobi routines removed, cosmetics (05.04.06,MK)
!> \author Joost VandeVondele (2003-08)
! *****************************************************************************
MODULE cp_fm_diag

  USE cp_blacs_calls,                  ONLY: cp_blacs_gridexit,&
                                             cp_blacs_gridinit
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release
  USE cp_fm_basic_linalg,              ONLY: cp_fm_syrk,&
                                             cp_fm_upper_to_full
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: default_string_length, dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_comm_free,&
                                             mp_comm_split,&
                                             mp_comm_split_direct,&
                                             mp_sync
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#if defined (__ELPA)
  USE ELPA1
  USE ELPA2
#endif

#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_diag'

  ! these saved variables are DIAGONALIZATION global
   INTEGER, SAVE :: diag_type = 0

  ! Public subroutines

  PUBLIC :: choose_eigv_solver,&
            cp_fm_block_jacobi,&
            cp_fm_power,&
            cp_fm_syevd,&
            cp_fm_syevr,&
            cp_fm_syevx,&
            cp_fm_elpa,&
            diag_init

CONTAINS

! *****************************************************************************
!> \brief Setup the diagonalization library to be used
!         Check of availability not yet fully implemented
!         It should change library to Skalapack if others are not available
!  \par   diag_library flag from GLOBAL section in input
!  \authot  MI 11.2013       
! *****************************************************************************
  SUBROUTINE diag_init(diag_lib,switched,error)
    CHARACTER(LEN=*), INTENT(IN)             :: diag_lib
    LOGICAL, INTENT(INOUT)                   :: switched
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'diag_init', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: try_library

! skalapack is the default and always linked

    IF(diag_lib.EQ."SL") THEN
      try_library = .FALSE.
      diag_type = 1
    ELSE
      try_library =.TRUE.
    END IF

    IF(try_library) THEN
        IF(diag_lib.EQ."ELPA") THEN

#if defined (__ELPA)
            diag_type = 3
#else
          ! ELPA library not linked, switch to SL
             diag_type = 1
             switched = .TRUE.
#endif
         ELSE IF(diag_lib.EQ."SL2") THEN

#if defined (__SCALAPACK2)
            diag_type = 2
#else
          ! SL2 library not linked, switch to SL
             diag_type = 1
             switched = .TRUE.
#endif
         END IF
    END IF

    ! Check that one of the diagonalization type is set
    IF(diag_type<1) THEN
      ! something wrong
      CALL stop_program(routineN,moduleN,__LINE__,&
                         "Unknown DIAG type")
    END IF

  END SUBROUTINE diag_init


! *****************************************************************************
!> \brief   Choose the Eigensolver depending on which library is available
!           ELPA seems to be unstable for small systems
!> \par     info If present returns error code and prevents program stops.
!>               Works currently only for cp_fm_syevd with scalapack.
!>               Other solvers will end the program regardless of PRESENT(info).
! *****************************************************************************
  SUBROUTINE choose_eigv_solver(matrix,eigenvectors,eigenvalues, info, error)

    TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    INTEGER, INTENT(OUT), OPTIONAL           :: info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'choose_eigv_solver', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: myinfo, nmo

    myinfo = 0

    nmo = SIZE(eigenvalues,1)

    IF(diag_type==3) THEN
      CALL cp_fm_elpa(matrix,eigenvectors,eigenvalues,error=error)
    ELSE IF(diag_type==2) THEN
      CALL cp_fm_syevr(matrix,eigenvectors,eigenvalues,1,nmo,error=error)
    ELSE IF(diag_type==1) THEN
      CALL cp_fm_syevd(matrix,eigenvectors,eigenvalues,info=myinfo,error=error)
    END IF

    IF (PRESENT(info)) info = myinfo 
  END SUBROUTINE choose_eigv_solver

! *****************************************************************************
!> \brief   Computes all eigenvalues and vectors of a real symmetric matrix
!>          significantly faster than syevx, scales also much better.
!>          Needs workspace to allocate all the eigenvectors
!> \par     matrix is supposed to be in upper triangular form, and overwritten by this routine
!> \par     info If present returns error code and prevents program stops.
!>               Works currently only for scalapack.
!>               Other solvers will end the program regardless of PRESENT(info).
! *****************************************************************************
  SUBROUTINE cp_fm_syevd(matrix,eigenvectors,eigenvalues,info,error)

    TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    INTEGER, INTENT(OUT), OPTIONAL           :: info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_syevd', &
      routineP = moduleN//':'//routineN

    INTEGER :: fake_descriptor(9), handle, istat, mepos_old, &
      myinfo, n, ngroups, nmo, num_pe_new, num_pe_old, subgroup
    INTEGER, DIMENSION(:), POINTER           :: group_distribution, &
                                                group_partition
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: fake_local_data(1,1)
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eig
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_new
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_new
    TYPE(cp_fm_type), POINTER                :: eigenvectors_new, matrix_new
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_new
#if ! defined(__SCALAPACK)
    INTEGER                                  :: liwork, lwork
    INTEGER, DIMENSION(:), POINTER           :: iwork
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m, v
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
#endif

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.
    myinfo = 0

    n = matrix%matrix_struct%nrow_global
    ALLOCATE(eig(n), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)


#if defined(__SCALAPACK)

    ! first figure out the optimal number of cpus
    ! this is pure heuristics, based on rosa timings
    ! that demonstrate that timings go up sharply if too many tasks are used
    ! we take a multiple of 4, and approximately n/60
    num_pe_old=matrix%matrix_struct%para_env%num_pe
    num_pe_new=((n+60*4-1)/(60*4))*4
    para_env=>matrix%matrix_struct%para_env
    mepos_old=para_env%mepos

    ! if the optimal is smaller than num_pe, we will redistribute the input matrix
    ! it seems sensible to refactor redistributing matrices as this might be useful elsewhere
    IF (num_pe_new<num_pe_old) THEN


       ! split comm, the first num_pe_new tasks will do the work
       ALLOCATE(group_distribution(0:num_pe_old-1))
       ALLOCATE(group_partition(0:1))
       group_partition=(/num_pe_new,num_pe_old-num_pe_new/)
       CALL mp_comm_split(comm=para_env%group,sub_comm=subgroup,&
                          ngroups=ngroups,group_distribution=group_distribution,&
                          n_subgroups=2, group_partition=group_partition)

       IF (group_distribution(mepos_old)==0) THEN

          ! create para_env, might need a proper error, bound to this para_env
          NULLIFY(para_env_new)
          CALL cp_para_env_create(para_env_new,subgroup,error=error)
          ! test a sync
          CALL mp_sync(para_env_new%group)

          ! create blacs, should inherit the preferences for the layout and so on, from the higher level
          NULLIFY(blacs_env_new)
          CALL cp_blacs_env_create(blacs_env=blacs_env_new, para_env=para_env_new, error=error)

          ! create new matrix
          NULLIFY(fm_struct_new)
          CALL cp_fm_struct_create(fmstruct=fm_struct_new, para_env=para_env_new, context=blacs_env_new,&
                                   nrow_global=n, ncol_global=n, error=error)
          CALL cp_fm_create(matrix_new, matrix_struct=fm_struct_new, name="yevd_new_mat", error=error)
          CALL cp_fm_create(eigenvectors_new, matrix_struct=fm_struct_new, name="yevd_new_vec", error=error)

          ! redistribute old
          CALL pdgemr2d(n,n,matrix%local_data(1,1),1,1,matrix%matrix_struct%descriptor, &
                        matrix_new%local_data(1,1),1,1,matrix_new%matrix_struct%descriptor, &
                        matrix%matrix_struct%context%group)

          ! call scalapack
          CALL cp_fm_syevd_base(matrix_new,eigenvectors_new,eig,myinfo,error)

          ! redistribute results
          CALL pdgemr2d(n,n,eigenvectors_new%local_data(1,1),1,1,eigenvectors_new%matrix_struct%descriptor, &
                        eigenvectors%local_data(1,1),1,1,eigenvectors%matrix_struct%descriptor, &
                        eigenvectors%matrix_struct%context%group)

          ! free stuff
          CALL cp_fm_struct_release(fm_struct_new,error=error)
          CALL cp_fm_release(matrix_new,error=error)
          CALL cp_fm_release(eigenvectors_new,error=error)
          CALL cp_blacs_env_release(blacs_env_new,error=error)
          CALL cp_para_env_release(para_env_new,error=error)

       ELSE
          ! these tasks must help redistribute (they own part of the data),
          ! but need fake 'new' data, and their descriptor must indicate this with -1
          ! see also scalapack comments on pdgemr2d
          fake_descriptor=-1
          CALL pdgemr2d(n,n,matrix%local_data(1,1),1,1,matrix%matrix_struct%descriptor, &
                        fake_local_data(1,1),1,1,fake_descriptor, &
                        matrix%matrix_struct%context%group)

          CALL pdgemr2d(n,n,fake_local_data(1,1),1,1,fake_descriptor, &
                        eigenvectors%local_data(1,1),1,1,eigenvectors%matrix_struct%descriptor, &
                        eigenvectors%matrix_struct%context%group)

          ! free stuff
          CALL mp_comm_free(subgroup)
       ENDIF

       ! finally, also the eigenvalues need to end up on the non-group member tasks
       CALL mp_bcast(eig,0,para_env%group)

       ! free more stuff
       DEALLOCATE(group_distribution,group_partition)

    ELSE

        CALL cp_fm_syevd_base(matrix,eigenvectors,eig,myinfo,error)

    ENDIF

#else

    !MK Retrieve the optimal work array sizes first
    myinfo = 0
    lwork = -1
    liwork = -1
    m => matrix%local_data
    eig(:) = 0.0_dp

    ALLOCATE (work(1),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    work(:) = 0.0_dp
    ALLOCATE (iwork(1),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    iwork(:) = 0

    CALL dsyevd('V','U',n,m(1,1),n,eig(1),work(1),lwork,iwork(1),liwork,myinfo)

    IF (myinfo /= 0) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "ERROR in DSYEVD: Could not retrieve work array sizes")
    END IF

    ! Reallocate work arrays and perform diagonalisation
    lwork = INT(work(1))
    DEALLOCATE (work,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(work(lwork),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    liwork = iwork(1)
    DEALLOCATE (iwork,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(iwork(liwork),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    iwork(:) = 0

    CALL dsyevd('V','U',n,m(1,1),n,eig(1),work(1),lwork,iwork(1),liwork,myinfo)

    IF (myinfo /= 0) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Matrix diagonalization failed")
    END IF

    CALL cp_fm_to_fm(matrix,eigenvectors,error=error)

    DEALLOCATE (iwork,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (work,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
#endif

    IF(PRESENT(info)) myinfo = 0

    nmo = SIZE(eigenvalues,1)
    IF (nmo > n) THEN
      eigenvalues(1:n) = eig(1:n)
    ELSE
      eigenvalues(1:nmo) = eig(1:nmo)
    END IF

    DEALLOCATE (eig,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)

  END SUBROUTINE cp_fm_syevd

  SUBROUTINE cp_fm_syevd_base(matrix,eigenvectors,eig,info,error)

    TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:)              :: eig
    INTEGER, INTENT(OUT)                     :: info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_syevd_base', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat, liwork, lwork, &
                                                mypcol, myprow, n
    INTEGER, DIMENSION(9)                    :: descm, descv
    INTEGER, DIMENSION(:), POINTER           :: iwork
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m, v
    TYPE(cp_blacs_env_type), POINTER         :: context

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

#if defined(__SCALAPACK)

    failure = .FALSE.

    n = matrix%matrix_struct%nrow_global
    m => matrix%local_data
    context => matrix%matrix_struct%context
    myprow = context%mepos(1)
    mypcol = context%mepos(2)
    descm(:) = matrix%matrix_struct%descriptor(:)

    v => eigenvectors%local_data
    descv(:) = eigenvectors%matrix_struct%descriptor(:)

    liwork = 7*n + 8*context%num_pe(2) + 2
    ALLOCATE(iwork(liwork),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ! work space query

    lwork = -1
    ALLOCATE(work(1),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL pdsyevd('V','U',n,m(1,1),1,1,descm,eig(1),v(1,1),1,1,descv, &
                 work(1),lwork,iwork(1),liwork,info)

    ! look here for a PDORMTR warning :-)
    ! this routine seems to need more workspace than reported
    ! (? lapack bug ...)
    ! arbitrary additional memory  ... we give 100000 more words
    ! (it seems to depend on the block size used)

    lwork = NINT(work(1)+100000)
!    lwork = NINT(work(1))
    DEALLOCATE (work,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(work(lwork),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL pdsyevd('V','U',n,m(1,1),1,1,descm,eig(1),v(1,1),1,1,descv, &
                 work(1),lwork,iwork(1),liwork,info)

    IF(info /= 0) CALL stop_program &
       (routineN,moduleN,__LINE__, "Matrix diagonalization failed")

    DEALLOCATE (work,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE (iwork,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

#endif


    CALL timestop(handle)

  END SUBROUTINE cp_fm_syevd_base

! *****************************************************************************
!> \brief   compute eigenvalues and optionally eigenvectors of a real symmetric matrix using scalapack.
!>          If eigenvectors are required this routine will replicate a full matrix on each CPU...
!>          if more than a handful of vectors are needed, use cp_fm_syevd instead
!> \par     matrix is supposed to be in upper triangular form, and overwritten by this routine
!>          neig   is the number of vectors needed (default all)
!>          work_syevx evec calculation only, is the fraction of the working buffer allowed (1.0 use full buffer)
!>                     reducing this saves time, but might cause the routine to fail
! *****************************************************************************
  SUBROUTINE cp_fm_syevx(matrix,eigenvectors,eigenvalues,neig,work_syevx,error)

    ! Diagonalise the symmetric n by n matrix using the LAPACK library.

    TYPE(cp_fm_type), POINTER                    :: matrix
    TYPE(cp_fm_type), POINTER, OPTIONAL          :: eigenvectors
    REAL(KIND = dp), OPTIONAL, INTENT(IN)        :: work_syevx
    INTEGER, INTENT(IN), OPTIONAL                :: neig
    REAL(KIND = dp), DIMENSION(:), INTENT(OUT)   :: eigenvalues
    TYPE(cp_error_type), INTENT(inout)  :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "cp_fm_syevx",&
                                   routineP = moduleN//":"//routineN
    REAL(KIND = dp), PARAMETER  :: orfac = -1.0_dp,&
                                   vl = 0.0_dp,&
                                   vu = 0.0_dp

    REAL(KIND = dp) :: abstol, work_syevx_local
    INTEGER  :: handle,info,istat,liwork,lwork,m,mypcol,myprow,n,nb,&
                nn,np0,npcol,npe,nprow,nq0,nz,output_unit,neig_local
    LOGICAL  :: failure, ionode, needs_evecs

    INTEGER, DIMENSION(9) :: desca,descz

    TYPE(cp_blacs_env_type), POINTER           :: context
    TYPE(cp_logger_type), POINTER              :: logger
    REAL(KIND = dp), DIMENSION(:), ALLOCATABLE :: gap,w,work
    INTEGER, DIMENSION(:), ALLOCATABLE         :: iclustr,ifail,iwork
    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: a,z
    CHARACTER(LEN=1)                           :: job_type

    REAL(KIND = dp), EXTERNAL :: dlamch

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL  :: iceil,numroc
#else
    INTEGER, EXTERNAL  :: ilaenv
#endif
    ! -------------------------------------------------------------------------

    failure = .FALSE.
    ! by default all
    n = matrix%matrix_struct%nrow_global
    neig_local=n
    IF (PRESENT(neig)) neig_local=neig
    IF (neig_local == 0) RETURN

    CALL timeset(routineN,handle)

    needs_evecs=PRESENT(eigenvectors)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    n = matrix%matrix_struct%nrow_global

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


    context =>  matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)


    ALLOCATE (w(n),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    eigenvalues(:) = 0.0_dp
#if defined(__SCALAPACK)

    IF (matrix%matrix_struct%nrow_block /= matrix%matrix_struct%ncol_block) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Invalid blocksize (no square blocks)")
    END IF

    a => matrix%local_data
    desca(:) = matrix%matrix_struct%descriptor(:)

    IF (needs_evecs) THEN
      z => eigenvectors%local_data
      descz(:) = eigenvectors%matrix_struct%descriptor(:)
    ELSE
      ! z will not be referenced
      z => matrix%local_data
      descz=desca
    ENDIF

    ! Get the optimal work storage size

    npe = nprow*npcol
    nb = matrix%matrix_struct%nrow_block
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
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    gap = 0.0_dp
    ALLOCATE (iclustr(2*npe),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    iclustr = 0
    ALLOCATE (ifail(n),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ifail = 0
    ALLOCATE (iwork(liwork),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (work(lwork),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

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
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (iclustr,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (ifail,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (iwork,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (work,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

#else

    a => matrix%local_data
    IF (needs_evecs) THEN
       z => eigenvectors%local_data
    ELSE
       ! z will not be referenced
       z => matrix%local_data
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
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (iwork,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (work,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

#endif
    eigenvalues(1:neig_local) = w(1:neig_local)
    DEALLOCATE (w,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_syevx

! *****************************************************************************
!> \brief  computes selected eigenvalues and, optionally, eigenvectors of
!>        a real symmetric matrix A distributed in 2D blockcyclic format by
!>       calling the recommended sequence of ScaLAPACK routines.
!>
!> \par     matrix is supposed to be in upper triangular form, and overwritten by this routine
!>          subsets of eigenvalues/vectors can be selected by
!>          specifying a range of values or a range of indices for the desired eigenvalues.
! *****************************************************************************
  SUBROUTINE cp_fm_syevr(matrix,eigenvectors,eigenvalues,ilow,iup,error)

    TYPE(cp_fm_type), POINTER                    :: matrix
    TYPE(cp_fm_type), POINTER, OPTIONAL          :: eigenvectors
    REAL(KIND = dp), DIMENSION(:), INTENT(OUT)   :: eigenvalues
    INTEGER, INTENT(IN), OPTIONAL                :: ilow,iup
    TYPE(cp_error_type), INTENT(inout)  :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "cp_fm_syevr",&
                                   routineP = moduleN//":"//routineN
    REAL(KIND = dp), PARAMETER  :: orfac = -1.0_dp,&
                                   vl = 0.0_dp,&
                                   vu = 0.0_dp

    CHARACTER(LEN=1)                           :: job_type
    INTEGER :: handle, ilow_local, info, istat, iup_local, lwork, liwork, mypcol,myprow, n, neig
    INTEGER, DIMENSION(9) :: desca, descz
    INTEGER, DIMENSION(:), ALLOCATABLE  :: iwork
    LOGICAL :: failure, ionode, needs_evecs
    REAL(dp), DIMENSION(:), ALLOCATABLE :: w, work
    REAL(dp), DIMENSION(:,:), POINTER :: a, z

    TYPE(cp_blacs_env_type), POINTER           :: context
    TYPE(cp_logger_type), POINTER              :: logger

    REAL(KIND = dp), EXTERNAL :: dlamch

#if defined(__SCALAPACK)
#else
    REAL(dp)                            :: abstol
    INTEGER                             :: m, nb
    INTEGER, DIMENSION(:), ALLOCATABLE  :: ifail
    INTEGER, EXTERNAL                   :: ilaenv
#if defined(__SCALAPACK2)
    INTEGER                             :: nz
#endif
#endif

    ! -------------------------------------------------------------------------


    ! by default all
    n = matrix%matrix_struct%nrow_global
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
    n = matrix%matrix_struct%nrow_global

    ! set scalapack job type
    IF (needs_evecs) THEN
       job_type="V"
    ELSE
       job_type="N"
    ENDIF

    context =>  matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)

    ALLOCATE(w(n), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    eigenvalues(:) = 0.0_dp

#if defined(__SCALAPACK)

    IF (matrix%matrix_struct%nrow_block /= matrix%matrix_struct%ncol_block) THEN
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

    a => matrix%local_data
    desca(:) = matrix%matrix_struct%descriptor(:)

    IF (needs_evecs) THEN
      z => eigenvectors%local_data
      descz(:) = eigenvectors%matrix_struct%descriptor(:)
    ELSE
      ! z will not be referenced
      z => matrix%local_data
      descz=desca
    ENDIF

    ! First Call: Determine the needed work_space
    lwork = -1
    ALLOCATE(work(5*n), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(iwork(6*n), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
#if defined (__SCALAPACK2)
    CALL pdsyevr(job_type,'I','U',n,a,1,1,desca,vl,vu,ilow_local,iup_local,m,nz,w(1),z,1,1,descz,work,lwork,iwork,liwork,info)
#endif
    lwork = INT(work(1))
    lwork = NINT(work(1)+300000)
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
#if defined (__SCALAPACK2)
    CALL pdsyevr(job_type,'I','U',n,a,1,1,desca,vl,vu,ilow_local,iup_local,m,nz,w(1),z,1,1,descz,work,lwork,iwork,liwork,info)
#endif

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

    a => matrix%local_data
    IF (needs_evecs) THEN
       z => eigenvectors%local_data
    ELSE
       ! z will not be referenced
       z => matrix%local_data
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

    CALL timestop(handle)

  END SUBROUTINE cp_fm_syevr

! *****************************************************************************
 SUBROUTINE cp_fm_elpa(matrix,eigenvectors,eigenvalues,error)

    TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:)              :: eigenvalues
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_elpa', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: comm_col, comm_row, group, &
                                                handle, istat, mypcol, &
                                                myprow, n, n_rows, nblk, neig
    INTEGER, DIMENSION(9)                    :: descm, descv
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eval
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m, m1, v
    TYPE(cp_blacs_env_type), POINTER         :: context

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

#if defined(__SCALAPACK)

    failure = .FALSE.

    n = matrix%matrix_struct%nrow_global
    m => matrix%local_data
    context => matrix%matrix_struct%context
    myprow = context%mepos(1)
    mypcol = context%mepos(2)
    group =  matrix%matrix_struct%para_env%group
    descm(:) = matrix%matrix_struct%descriptor(:)


    ! For ELPA, the MPI communicators along rows/cols are sufficient
    ! mpi communicators are created, needed for communicating within
    ! rows or columns of processes
    CALL mp_comm_split_direct( group, comm_row, mypcol, myprow)
    CALL mp_comm_split_direct(group, comm_col, myprow, mypcol)

    v => eigenvectors%local_data
    descv(:) = eigenvectors%matrix_struct%descriptor(:)

    ALLOCATE (m1(SIZE(m,1),SIZE(m,2)),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    m1 = m

    ! elpa needs the full matrix
    CALL cp_fm_upper_to_full(matrix,eigenvectors,error=error)

    n_rows = matrix%matrix_struct%local_leading_dimension
    nblk = matrix%matrix_struct%nrow_block
    neig = SIZE(eigenvalues,1)
    ! the full eigenvalues vector is needed
    ALLOCATE(eval(n),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
#if defined (__ELPA)
    ! Calculate eigenvalues/eigenvectors
    CALL solve_evp_real_2stage(n,neig,m,n_rows, eval,v,n_rows,nblk,comm_row,comm_col,group)
#else
    CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
#endif
    eigenvalues(1:neig) = eval(1:neig)

    DEALLOCATE(eval)
    DEALLOCATE (m1,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    ! mpi communicators are freed
    CALL mp_comm_free(comm_row)
    CALL mp_comm_free(comm_col)

#endif

    CALL timestop(handle)

  END SUBROUTINE cp_fm_elpa


! *****************************************************************************
  SUBROUTINE cp_fm_power(matrix,work,exponent,threshold,n_dependent,verbose,error)

   ! Raise the real symmetric n by n matrix to the power given by
   ! the exponent. All eigenvectors with a corresponding eigenvalue lower
   ! than threshold are quenched. result in matrix

   ! - Creation (29.03.1999, Matthias Krack)
   ! - Parallelised using BLACS and ScaLAPACK (06.06.2001,MK)

   TYPE(cp_fm_type), POINTER                 :: matrix,work
   REAL(KIND = dp), INTENT(IN)               :: exponent,threshold
   INTEGER, INTENT(OUT)                      :: n_dependent
   LOGICAL, INTENT(IN), OPTIONAL             :: verbose
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_power', &
      routineP = moduleN//':'//routineN

   INTEGER         :: handle,icol_global,icol_local,ipcol,iprow,irow_global,&
                      irow_local,istat,mypcol,myprow,ncol_block,ncol_global,npcol,&
                      nprow,nrow_block,nrow_global
   LOGICAL :: failure, my_verbose
   REAL(KIND = dp) :: condition_number,f,p
   REAL(KIND = dp), DIMENSION(:), ALLOCATABLE :: eigenvalues
   REAL(KIND = dp), DIMENSION(:,:), POINTER   :: eigenvectors
   TYPE(cp_blacs_env_type), POINTER           :: context

#if defined(__SCALAPACK)
   INTEGER, EXTERNAL :: indxg2l,indxg2p
#endif
    ! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.
    my_verbose = .FALSE.
    IF(PRESENT(verbose)) my_verbose = verbose

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
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    eigenvalues(:) = 0.0_dp

    ! Compute the eigenvectors and eigenvalues

#if defined (__ELPA)
    CALL cp_fm_elpa(matrix,work,eigenvalues,error=error)
#else
    CALL cp_fm_syevd(matrix,work,eigenvalues,error=error)
#endif

#if defined(__SCALAPACK)
    nrow_block = work%matrix_struct%nrow_block
    ncol_block = work%matrix_struct%ncol_block

    eigenvectors => work%local_data

    ! Build matrix**exponent with eigenvector quenching

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
    CALL cp_fm_syrk("U","N",ncol_global,1.0_dp,work,1,1,0.0_dp,matrix,error=error)
    CALL cp_fm_upper_to_full(matrix,work,error=error)

    ! Print some warnings/notes

    IF (matrix%matrix_struct%para_env%mepos == 0 .AND. my_verbose) THEN
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
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_power

! *****************************************************************************
  SUBROUTINE cp_fm_block_jacobi(matrix,eigenvectors,eigval,thresh,&
                                start_sec_block)

    ! Calculates block diagonalizazion from full symmetric matrix
    ! It has its origin in cp_fm_syevx. This routine rotates only elements
    ! which are larger than a threshold thresh.
    ! start_sec_block is the start of the second block.
    ! IT DOES ONLY ONE SWEEP!

    ! - Creation (07.10.2002, Martin Fengler)
    ! - Cosmetics (05.04.06,MK)

    TYPE(cp_fm_type), POINTER                 :: eigenvectors,matrix
    REAL(KIND = dp), DIMENSION(:), INTENT(IN) :: eigval
    INTEGER, INTENT(IN)                       :: start_sec_block
    REAL(KIND = dp), INTENT(IN)               :: thresh

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_block_jacobi', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle
    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: a,ev

    REAL(KIND = dp) :: tan_theta,tau,c,s
    INTEGER  :: q,p,N,q_loc
    REAL(KIND = dp), DIMENSION (:),ALLOCATABLE :: c_ip

#if defined(__SCALAPACK)
    TYPE(cp_blacs_env_type), POINTER :: context

    INTEGER :: myprow,mypcol,nprow,npcol,ictxt_loc,block_dim_row,block_dim_col,&
               info,ev_row_block_size,iam,source,allgrp,mynumrows,ictxt,mype,npe
    REAL(KIND = dp), DIMENSION(:,:), ALLOCATABLE :: a_loc, ev_loc
    INTEGER, DIMENSION(9)                        :: desca,descz,desc_a_block,&
                                                    desc_ev_loc
    INTEGER, EXTERNAL :: numroc,indxl2g,indxg2l
#endif

    ! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

#if defined(__SCALAPACK)

    ! Umgebung fuer die Kohn-Sham-Matrix und die Eigenvektoren uebernehmen

    context =>  matrix%matrix_struct%context
    source =    matrix%matrix_struct%para_env%source
    allgrp =    matrix%matrix_struct%para_env%group

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

    CALL timestop(handle)

  END SUBROUTINE cp_fm_block_jacobi

END MODULE cp_fm_diag
