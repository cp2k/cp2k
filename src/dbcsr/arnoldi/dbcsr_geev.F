!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief provides a unified interface to lapack geev routines 
!> \par History
!>       2014.09 created [Florian Schiffmann]
!> \author Florian Schiffmann
! *****************************************************************************

MODULE dbcsr_geev
  USE kinds,                           ONLY: real_4,&
                                             real_8

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_geev'

  PUBLIC :: dbcsr_general_local_diag

  INTERFACE dbcsr_general_local_diag
     MODULE PROCEDURE dbcsr_sgeev, dbcsr_dgeev, dbcsr_zgeev, dbcsr_cgeev
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param jobvl ...
!> \param jobvr ...
!> \param matrix ...
!> \param ndim ...
!> \param evals ...
!> \param revec ...
!> \param levec ...
! *****************************************************************************
  SUBROUTINE dbcsr_sgeev(jobvl,jobvr,matrix,ndim,evals,revec,levec)
    CHARACTER(1)                             :: jobvl, jobvr
    REAL(real_4), DIMENSION(:, :)            :: matrix
    INTEGER                                  :: ndim
    COMPLEX(real_4), DIMENSION(:)            :: evals
    COMPLEX(real_4), DIMENSION(:, :)         :: revec, levec

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_sgeev', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, info, lwork
    REAL(real_4)                             :: tmp_array(ndim,ndim), &
                                                work(20*ndim)
    REAL(real_4), DIMENSION(ndim)            :: eval1, eval2
    REAL(real_4), DIMENSION(ndim, ndim)      :: evec_l, evec_r

    eval1=REAL(0.0,real_4); eval2=REAL(0.0,real_4)
    ! ask lapack how much space it would like in the work vector, don't ask me why
    lwork=-1
    CALL SGEEV(jobvl, jobvr, ndim, tmp_array, ndim, eval1, eval2, evec_l, ndim, evec_r, ndim, work, lwork, info)
    lwork=MIN(20*ndim, INT(work(1)))

    tmp_array(:,:)=matrix(:,:)
    CALL SGEEV(jobvl, jobvr, ndim, tmp_array, ndim, eval1, eval2, evec_l, ndim, evec_r, ndim, work, lwork, info)

    ! compose the eigenvectors, lapacks way of storing them is a pain
    ! if eval is complex, then the complex conj pair of evec can be constructed from the i and i+1st evec
    i=1
    DO WHILE (i.le.ndim)
      IF(ABS(eval2(i)).LT.EPSILON(REAL(0.0,real_4)))THEN
        revec(:,i)=CMPLX(evec_r(:,i),REAL(0.0,real_4),real_4)
        levec(:,i)=CMPLX(evec_l(:,i),REAL(0.0,real_4),real_4)
        i=i+1
      ELSE IF(eval2(i).GT.EPSILON(REAL(0.0,real_4)))THEN
        revec(:,i)=CMPLX(evec_r(:,i),evec_r(:,i+1),real_4)
        revec(:,i+1)=CMPLX(evec_r(:,i),-evec_r(:,i+1),real_4)
        levec(:,i)=CMPLX(evec_l(:,i),evec_l(:,i+1),real_4)
        levec(:,i+1)=CMPLX(evec_l(:,i),-evec_l(:,i+1),real_4)
        i=i+2
      ELSE
        STOP 'something went wrong while sorting the EV in dbcsr_geev'
      END IF
    END DO
    ! this is to keep the interface consistent with complex geev
    DO i=1,ndim
       evals(i)=CMPLX(eval1(i),eval2(i),real_4)
    END DO

  END SUBROUTINE  dbcsr_sgeev

! *****************************************************************************
!> \brief ...
!> \param jobvl ...
!> \param jobvr ...
!> \param matrix ...
!> \param ndim ...
!> \param evals ...
!> \param revec ...
!> \param levec ...
! *****************************************************************************
  SUBROUTINE dbcsr_dgeev(jobvl,jobvr,matrix,ndim,evals,revec,levec)
    CHARACTER(1)                             :: jobvl, jobvr
    REAL(real_8), DIMENSION(:, :)            :: matrix
    INTEGER                                  :: ndim
    COMPLEX(real_8), DIMENSION(:)            :: evals
    COMPLEX(real_8), DIMENSION(:, :)         :: revec, levec

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_dgeev', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, info, lwork
    REAL(real_8)                             :: tmp_array(ndim,ndim), &
                                                work(20*ndim)
    REAL(real_8), DIMENSION(ndim)            :: eval1, eval2
    REAL(real_8), DIMENSION(ndim, ndim)      :: evec_l, evec_r

    eval1=REAL(0.0,real_8); eval2=REAL(0.0,real_8)
    ! ask lapack how much space it would like in the work vector, don't ask me why
    lwork=-1
    CALL DGEEV(jobvl, jobvr, ndim, tmp_array, ndim, eval1, eval2, evec_l, ndim, evec_r, ndim, work, lwork, info)
    lwork=MIN(20*ndim, INT(work(1)))

    tmp_array(:,:)=matrix(:,:)
    CALL DGEEV(jobvl, jobvr, ndim, tmp_array, ndim, eval1, eval2, evec_l, ndim, evec_r, ndim, work, lwork, info)

    ! compose the eigenvectors, lapacks way of storing them is a pain
    ! if eval is complex, then the complex conj pair of evec can be constructed from the i and i+1st evec
    i=1
    DO WHILE (i.le.ndim)
      IF(ABS(eval2(i)).LT.EPSILON(REAL(0.0,real_8)))THEN
        revec(:,i)=CMPLX(evec_r(:,i),REAL(0.0,real_8),real_8)
        levec(:,i)=CMPLX(evec_l(:,i),REAL(0.0,real_8),real_8)
        i=i+1
      ELSE IF(eval2(i).GT.EPSILON(REAL(0.0,real_8)))THEN
        revec(:,i)=CMPLX(evec_r(:,i),evec_r(:,i+1),real_8)
        revec(:,i+1)=CMPLX(evec_r(:,i),-evec_r(:,i+1),real_8)
        levec(:,i)=CMPLX(evec_l(:,i),evec_l(:,i+1),real_8)
        levec(:,i+1)=CMPLX(evec_l(:,i),-evec_l(:,i+1),real_8)
        i=i+2
      ELSE
        STOP 'something went wrong while sorting the EV in dbcsr_geev'
      END IF
    END DO

    ! this is to keep the interface consistent with complex geev
    DO i=1,ndim
       evals(i)=CMPLX(eval1(i),eval2(i),real_8)
    END DO

  END SUBROUTINE  dbcsr_dgeev 

! *****************************************************************************
!> \brief ...
!> \param jobvl ...
!> \param jobvr ...
!> \param matrix ...
!> \param ndim ...
!> \param evals ...
!> \param revec ...
!> \param levec ...
! *****************************************************************************
  SUBROUTINE dbcsr_zgeev(jobvl,jobvr,matrix,ndim,evals,revec,levec)
    CHARACTER(1)                             :: jobvl, jobvr
    COMPLEX(real_8), DIMENSION(:, :)         :: matrix
    INTEGER                                  :: ndim
    COMPLEX(real_8), DIMENSION(:)            :: evals
    COMPLEX(real_8), DIMENSION(:, :)         :: revec, levec

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_zgeev', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: info, lwork
    COMPLEX(real_8)                          :: work(20*ndim), &
                                                tmp_array(ndim,ndim)
    REAL(real_8)                             :: work2(2*ndim)

    evals=CMPLX(0.0, 0.0, real_8) 
    ! ask lapack how much space it would like in the work vector, don't ask me why
    lwork=-1
    CALL ZGEEV(jobvl, jobvr, ndim, tmp_array, ndim, evals, levec, ndim, revec, ndim, work, lwork, work2, info)
    lwork=MIN(20*ndim, INT(work(1)))

    tmp_array(:,:)=matrix(:,:)
    CALL ZGEEV(jobvl, jobvr, ndim, tmp_array, ndim, evals, levec, ndim, revec, ndim, work, lwork, work2, info)

  END SUBROUTINE  dbcsr_zgeev 

! *****************************************************************************
!> \brief ...
!> \param jobvl ...
!> \param jobvr ...
!> \param matrix ...
!> \param ndim ...
!> \param evals ...
!> \param revec ...
!> \param levec ...
! *****************************************************************************
  SUBROUTINE dbcsr_cgeev(jobvl,jobvr,matrix,ndim,evals,revec,levec)
    CHARACTER(1)                             :: jobvl, jobvr
    COMPLEX(real_4), DIMENSION(:, :)         :: matrix
    INTEGER                                  :: ndim
    COMPLEX(real_4), DIMENSION(:)            :: evals
    COMPLEX(real_4), DIMENSION(:, :)         :: revec, levec

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_cgeev', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: info, lwork
    COMPLEX(real_4)                          :: work(20*ndim), &
                                                tmp_array(ndim,ndim)
    REAL(real_4)                             :: work2(2*ndim)

    evals=CMPLX(0.0, 0.0, real_4)
    ! ask lapack how much space it would like in the work vector, don't ask me why
    lwork=-1
    CALL CGEEV(jobvl, jobvr, ndim, tmp_array, ndim, evals, levec, ndim, revec, ndim, work, lwork, work2, info)
    lwork=MIN(20*ndim, INT(work(1)))

    tmp_array(:,:)=matrix(:,:)
    CALL CGEEV(jobvl, jobvr, ndim, tmp_array, ndim, evals, levec, ndim, revec, ndim, work, lwork, work2, info)

  END SUBROUTINE  dbcsr_cgeev

END MODULE dbcsr_geev 
