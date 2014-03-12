!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief used for collecting some of the diagonalization shemes available for cp_cfm_type
!>      cp_fm_power also moved here as it is very related
!> \note
!>      first version : only one routine right now
!> \author Joost VandeVondele (2003-09)
! *****************************************************************************
MODULE cp_cfm_diag
  USE cp_cfm_types,                    ONLY: cp_cfm_type
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_cfm_diag'

PUBLIC :: cp_cfm_heevd

CONTAINS

! *****************************************************************************
!> \brief Perform a diagonalisation of a complex matrix
!> \par History
!>      - (De)Allocation checks updated (15.02.2011,MK)
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_cfm_heevd(matrix,eigenvectors,eigenvalues,error)

    TYPE(cp_cfm_type), POINTER               :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_heevd', &
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

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(eigenvectors),cp_failure_level,routineP,error,failure)

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
    lwork  = work(1)
    lrwork = rwork(1) + 1000000 ! needed to correct for a bug in scalapack, unclear how much the right number is
    liwork = iwork(1)
#else
    CALL ZHEEVD('V','U',n,m(1,1),SIZE(m,1),eigenvalues(1),&
                work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
    lwork  = work(1)
    lrwork = rwork(1)
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
    CALL PZHEEVD('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv,&
                  work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
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

END MODULE cp_cfm_diag
