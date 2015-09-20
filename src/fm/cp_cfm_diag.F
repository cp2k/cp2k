!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief used for collecting some of the diagonalization shemes available for cp_cfm_type
!>      cp_fm_power also moved here as it is very related
!> \note
!>      first version : only one routine right now
!> \author Joost VandeVondele (2003-09)
! *****************************************************************************
MODULE cp_cfm_diag
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_cholesky_decompose,&
                                             cp_cfm_triangular_invert,&
                                             cp_cfm_triangular_multiply
  USE cp_cfm_types,                    ONLY: cp_cfm_get_info,&
                                             cp_cfm_to_cfm,&
                                             cp_cfm_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#if defined (__HAS_IEEE_EXCEPTIONS)
  USE ieee_exceptions,                 ONLY: ieee_get_halting_mode,&
                                             ieee_set_halting_mode,&
                                             ieee_all
#endif

#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_cfm_diag'

PUBLIC :: cp_cfm_heevd, cp_cfm_geeig

CONTAINS

! *****************************************************************************
!> \brief Perform a diagonalisation of a complex matrix
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \par History
!>      - (De)Allocation checks updated (15.02.2011,MK)
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_cfm_heevd(matrix,eigenvectors,eigenvalues)

    TYPE(cp_cfm_type), POINTER               :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_heevd', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: work
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: m
    INTEGER                                  :: handle, info, istat, liwork, &
                                                lrwork, lwork, n
    INTEGER, DIMENSION(:), POINTER           :: iwork
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rwork
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: descm, descv
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: v
#if defined (__HAS_IEEE_EXCEPTIONS)
    LOGICAL, DIMENSION(5)                    :: halt
#endif
#endif

    CALL timeset(routineN,handle)
    failure = .FALSE.

    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,failure)
    CPPrecondition(ASSOCIATED(eigenvectors),cp_failure_level,routineP,failure)

    n = matrix%matrix_struct%nrow_global
    m => matrix%local_data
    ALLOCATE (iwork(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"iwork",int_size)
    ALLOCATE (rwork(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rwork",dp_size)
    ALLOCATE (work(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"work",2*dp_size)
    ! work space query
    lwork  = -1
    lrwork = -1
    liwork = -1

#if defined(__SCALAPACK)
    v => eigenvectors%local_data
    descm(:) = matrix%matrix_struct%descriptor(:)
    descv(:) = eigenvectors%matrix_struct%descriptor(:)
    CALL PZHEEVD('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv,&
                  work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
    lwork  = INT(REAL(work(1)))
    lrwork = INT(REAL(rwork(1))) + 1000000 ! needed to correct for a bug in scalapack, unclear how much the right number is
    liwork = iwork(1)
#else
    CALL ZHEEVD('V','U',n,m(1,1),SIZE(m,1),eigenvalues(1),&
                work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
    lwork  = INT(REAL(work(1)))
    lrwork = INT(REAL(rwork(1)))
    liwork = iwork(1)
#endif

    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"iwork")
    DEALLOCATE (rwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rwork")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"work")
    ALLOCATE (iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"iwork",int_size*liwork)
    iwork(:) = 0
    ALLOCATE (rwork(lrwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rwork",dp_size*lrwork)
    rwork(:) = 0.0_dp
    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"work",2*dp_size*lwork)
    work(:) = CMPLX(0.0_dp,0.0_dp,KIND=dp)

#if defined(__SCALAPACK)
! Scalapack takes advantage of IEEE754 exceptions for speedup.
! Therefore, we disable floating point traps temporarily.
#if defined (__HAS_IEEE_EXCEPTIONS)
    CALL ieee_get_halting_mode(IEEE_ALL, halt)
    CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif

    CALL PZHEEVD('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv,&
                  work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)

#if defined (__HAS_IEEE_EXCEPTIONS)
    CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
#else
    CALL ZHEEVD('V','U',n,m(1,1),SIZE(m,1),eigenvalues(1), &
                work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
    eigenvectors%local_data = matrix%local_data
#endif

    IF (info /= 0) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"Diagonalisation complex matrix failed")
    END IF
    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"iwork")
    DEALLOCATE (rwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rwork")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"work")

    CALL timestop(handle)

  END SUBROUTINE cp_cfm_heevd

! *****************************************************************************
!> \brief General Eigenvalue Problem  AX = BXE
!>        Single option version: Cholesky decomposition of B
!> \param amatrix ...
!> \param bmatrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param work ...
! *****************************************************************************
 SUBROUTINE cp_cfm_geeig(amatrix,bmatrix,eigenvectors,eigenvalues,work)

    TYPE(cp_cfm_type), POINTER               :: amatrix, bmatrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:)              :: eigenvalues
    TYPE(cp_cfm_type), POINTER               :: work

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_geeig', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat, nao, nmo
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: evals

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CALL cp_cfm_get_info(amatrix,nrow_global=nao)
    ALLOCATE (evals(nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)
    ! Cholesky decompose S=U(T)U
    CALL cp_cfm_cholesky_decompose(bmatrix)
    ! Invert to get U^(-1)
    CALL cp_cfm_triangular_invert(bmatrix)
    ! Reduce to get U^(-T) * H * U^(-1)
    CALL cp_cfm_triangular_multiply(bmatrix,amatrix,side="R")
    CALL cp_cfm_triangular_multiply(bmatrix,amatrix,transa_tr="C")
    ! Diagonalize
    CALL cp_cfm_heevd(matrix=amatrix,eigenvectors=work,eigenvalues=evals)
    ! Restore vectors C = U^(-1) * C*
    CALL cp_cfm_triangular_multiply(bmatrix,work)
    nmo = SIZE(eigenvalues)
    CALL cp_cfm_to_cfm(work,eigenvectors,nmo)
    eigenvalues(1:nmo) = evals(1:nmo)

    DEALLOCATE (evals,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)

    CALL timestop(handle)

  END SUBROUTINE cp_cfm_geeig

END MODULE cp_cfm_diag
