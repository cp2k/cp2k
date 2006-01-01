!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2003 CP2K developers group
!-----------------------------------------------------------------------------!


!!****h* cp2k/cp_cfm_diag [1.0] *
!!
!!   NAME
!!     cp_cfm_diag
!!
!!   FUNCTION
!!     used for collecting some of the diagonalization shemes available for cp_cfm_type
!!     cp_fm_power also moved here as it is very related
!!
!!   NOTES
!!     first version : only one routine right now
!!
!!   AUTHOR
!!     Joost VandeVondele (2003-09)
!!
! ************************************************************************************

MODULE cp_cfm_diag

  USE cp_cfm_types,                    ONLY: cp_cfm_type
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_cfm_diag'

!****************************************************************************

PUBLIC :: cp_cfm_heevd

CONTAINS

  SUBROUTINE cp_cfm_heevd(matrix,eigenvectors,eigenvalues,error)

    TYPE(cp_cfm_type), POINTER               :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    CPPrecondition(associated(matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(eigenvectors),cp_failure_level,routineP,error,failure)

    n = matrix%matrix_struct%nrow_global
    m => matrix%local_data
    ALLOCATE(iwork(1),work(1),rwork(1),STAT=istat)
    IF (istat.NE.0) CALL stop_memory(routineN,"iwork")
! work space query
    lwork  = -1
    lrwork = -1
    liwork = -1

#if defined(__SCALAPACK)
    v => eigenvectors%local_data
    descm(:) = matrix%matrix_struct%descriptor(:)
    descv(:) = eigenvectors%matrix_struct%descriptor(:)
    CALL PZHEEVD('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv, &
                  work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
    lwork=work(1)
    lrwork=rwork(1)+1000000 ! needed to correct for a bug in scalapack, unclear how much the right number is
    liwork=iwork(1)
    DEALLOCATE(work,iwork,rwork)
    ALLOCATE(iwork(liwork),work(lwork),rwork(lrwork),STAT=istat)
    iwork=0;work=CMPLX(0._dp,0._dp,dp);rwork=0._dp;
    IF (istat.NE.0) CALL stop_memory(routineN,"work")
    CALL PZHEEVD('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv, &
                  work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
#else
    CALL ZHEEVD('V','U',n,m(1,1),size(m,1),eigenvalues(1), &
                  work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
    lwork=work(1)
    lrwork=rwork(1)
    liwork=iwork(1)
    DEALLOCATE(work,iwork,rwork)
    ALLOCATE(iwork(liwork),work(lwork),rwork(lrwork),STAT=istat)
    IF (istat.NE.0) CALL stop_memory(routineN,"work")
    CALL ZHEEVD('V','U',n,m(1,1),size(m,1),eigenvalues(1), &
                  work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
    eigenvectors%local_data=matrix%local_data
#endif

    IF (info.NE.0) CALL stop_program("cp_cfm_heevd","unable to diagonalize matrix")
    DEALLOCATE(work,iwork,rwork)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_cfm_heevd

END MODULE cp_cfm_diag
