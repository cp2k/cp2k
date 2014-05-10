!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Aux types for the ALMO methods
!> \par History
!>       2011.12 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin 
! *****************************************************************************
MODULE almo_scf_aux_types
  USE almo_scf_aux2_methods,           ONLY: copy_object01_gen,&
                                             init_object01_gen,&
                                             op2_object01_gen,&
                                             release_object01_gen,&
                                             set_object01_gen
  USE almo_scf_aux2_types,             ONLY: object01_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_copy,&
                                             cp_dbcsr_create,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_release,&
                                             cp_dbcsr_set,&
                                             cp_dbcsr_trace,&
                                             cp_dbcsr_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  
  INTEGER, PARAMETER :: object00_orthogonal = 1
  
  INTEGER, PARAMETER :: object00_env_dbcsr = 1
  INTEGER, PARAMETER :: object00_env_domain = 2

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf_aux_types'

  PUBLIC :: object00_type,&
            init_object00_gen, release_object00, push_object00,&
            extrapolate_object00

  INTERFACE init_object00_gen
     MODULE PROCEDURE init_object00_dbcsr 
     MODULE PROCEDURE init_object00_domain 
  END INTERFACE
  
  TYPE object00_type
   
    INTEGER :: object00_env_type

    INTEGER :: buffer_length
    INTEGER :: max_buffer_length

   
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: m_var
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: m_err

    TYPE(object01_type), DIMENSION(:,:), ALLOCATABLE :: d_var
    TYPE(object01_type), DIMENSION(:,:), ALLOCATABLE :: d_err
    
    TYPE(object01_type), DIMENSION(:), ALLOCATABLE     :: m_b
    
    INTEGER :: in_point

    INTEGER :: error_type

  END TYPE object00_type

CONTAINS

! *****************************************************************************
!> \par History
!>       2011.12 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE init_object00_dbcsr(object00_env,sample_err,sample_var,error_type,&
    max_length,error)

    TYPE(object00_type), INTENT(INOUT)       :: object00_env
    TYPE(cp_dbcsr_type), INTENT(IN)          :: sample_err, sample_var
    INTEGER, INTENT(IN)                      :: error_type, max_length
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_object00_dbcsr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, im, ndomains
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    
    IF( max_length.le.0 ) THEN
       CPErrorMessage(cp_failure_level,routineP,"max_length is less than zero",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

    object00_env%object00_env_type=object00_env_dbcsr

    object00_env%max_buffer_length=max_length
    object00_env%buffer_length=0
    object00_env%error_type=error_type
    object00_env%in_point=1
    
    ALLOCATE(object00_env%m_err(object00_env%max_buffer_length))
    ALLOCATE(object00_env%m_var(object00_env%max_buffer_length))

    DO im=1,object00_env%max_buffer_length
       CALL cp_dbcsr_init(object00_env%m_err(im),error=error)
       CALL cp_dbcsr_create(object00_env%m_err(im),&
               template=sample_err,error=error)
       CALL cp_dbcsr_init(object00_env%m_var(im),error=error)
       CALL cp_dbcsr_create(object00_env%m_var(im),&
               template=sample_var,error=error)
    ENDDO
    
    ndomains=1
    ALLOCATE(object00_env%m_b(ndomains))
    CALL init_object01_gen(object00_env%m_b,error=error)
    object00_env%m_b(:)%domain=100 
    DO idomain=1, ndomains
       IF (object00_env%m_b(idomain)%domain.gt.0) THEN
          ALLOCATE(object00_env%m_b(idomain)%mdata(1,1))
          object00_env%m_b(idomain)%mdata(:,:)=0.0_dp
       ENDIF
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE init_object00_dbcsr

! *****************************************************************************
!> \par History
!>       2011.12 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE init_object00_domain(object00_env,sample_err,sample_var,error_type,&
    max_length,error)

    TYPE(object00_type), INTENT(INOUT)       :: object00_env
    TYPE(object01_type), DIMENSION(:), &
      INTENT(IN)                             :: sample_err, sample_var
    INTEGER, INTENT(IN)                      :: error_type, max_length
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_object00_domain', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, ndomains
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    
    IF( max_length.le.0 ) THEN
       CPErrorMessage(cp_failure_level,routineP,"max_length is less than zero",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

    object00_env%object00_env_type=object00_env_domain

    object00_env%max_buffer_length=max_length
    object00_env%buffer_length=0
    object00_env%error_type=error_type
    object00_env%in_point=1

    ndomains=SIZE(sample_err)

    ALLOCATE(object00_env%d_err(object00_env%max_buffer_length,ndomains))
    ALLOCATE(object00_env%d_var(object00_env%max_buffer_length,ndomains))
   
    CALL init_object01_gen(object00_env%d_var,error=error)
    CALL init_object01_gen(object00_env%d_err,error=error)
    
    ALLOCATE(object00_env%m_b(ndomains))
    CALL init_object01_gen(object00_env%m_b,error=error)
    object00_env%m_b(:)%domain=sample_err(:)%domain
    DO idomain=1, ndomains
       IF (object00_env%m_b(idomain)%domain.gt.0) THEN
          ALLOCATE(object00_env%m_b(idomain)%mdata(1,1))
          object00_env%m_b(idomain)%mdata(:,:)=0.0_dp
       ENDIF
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE init_object00_domain

! *****************************************************************************
!> \par History
!>       2011.12 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE push_object00(object00_env,var,err,d_var,d_err,error)
    TYPE(object00_type), INTENT(INOUT)       :: object00_env
    TYPE(cp_dbcsr_type), INTENT(IN), &
      OPTIONAL                               :: var, err
    TYPE(object01_type), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: d_var, d_err
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'push_object00', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, in_point, &
                                                irow, ndomains, &
                                                old_buffer_length
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: trace0
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: m_b_tmp

    CALL timeset(routineN,handle)

    IF (object00_env%object00_env_type.eq.object00_env_dbcsr) THEN
       IF ( .NOT.(PRESENT(var).AND.PRESENT(err)) ) THEN
          CPErrorMessage(cp_failure_level,routineP,"provide DBCSR matrices",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
    ELSE IF (object00_env%object00_env_type.eq.object00_env_domain) THEN
       IF ( .NOT.(PRESENT(d_var).AND.PRESENT(d_err)) ) THEN
          CPErrorMessage(cp_failure_level,routineP,"provide matrices",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
    ELSE
       CPErrorMessage(cp_failure_level,routineP,"illegal type",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    in_point=object00_env%in_point

    IF (object00_env%object00_env_type.eq.object00_env_dbcsr) THEN
       CALL cp_dbcsr_copy(object00_env%m_var(in_point),var,error=error) 
       CALL cp_dbcsr_copy(object00_env%m_err(in_point),err,error=error) 
    ELSE IF (object00_env%object00_env_type.eq.object00_env_domain) THEN
       CALL copy_object01_gen(d_var,object00_env%d_var(in_point,:),copy_data=.TRUE.,error=error) 
       CALL copy_object01_gen(d_err,object00_env%d_err(in_point,:),copy_data=.TRUE.,error=error) 
    ENDIF

    old_buffer_length=object00_env%buffer_length
    object00_env%buffer_length=object00_env%buffer_length+1
    IF (object00_env%buffer_length.gt.object00_env%max_buffer_length) &
       object00_env%buffer_length=object00_env%max_buffer_length
    
    ndomains=SIZE(object00_env%m_b)
    IF (old_buffer_length.lt.object00_env%buffer_length) THEN
       ALLOCATE(m_b_tmp(object00_env%buffer_length+1,object00_env%buffer_length+1))
       DO idomain=1,ndomains
          IF (object00_env%m_b(idomain)%domain.gt.0) THEN
             m_b_tmp(:,:)=0.0_dp
             m_b_tmp(1:object00_env%buffer_length,1:object00_env%buffer_length)=&
                object00_env%m_b(idomain)%mdata(:,:)
             DEALLOCATE(object00_env%m_b(idomain)%mdata)
             ALLOCATE(object00_env%m_b(idomain)%mdata(object00_env%buffer_length+1,&
                object00_env%buffer_length+1))
             object00_env%m_b(idomain)%mdata(:,:)=m_b_tmp(:,:)
          ENDIF
       ENDDO
       DEALLOCATE(m_b_tmp)
    ENDIF
    DO idomain=1,ndomains
       IF (object00_env%m_b(idomain)%domain.gt.0) THEN
          object00_env%m_b(idomain)%mdata(1,in_point+1)=-1.0_dp
          object00_env%m_b(idomain)%mdata(in_point+1,1)=-1.0_dp
          DO irow=1,object00_env%buffer_length
             IF (object00_env%object00_env_type.eq.object00_env_dbcsr) THEN
                trace0=object00_overlap(object00_env,&
                   A=object00_env%m_err(irow),B=object00_env%m_err(in_point),&
                   error=error)
             ELSE IF (object00_env%object00_env_type.eq.object00_env_domain) THEN
                trace0=object00_overlap(object00_env,&
                   d_A=object00_env%d_err(irow,idomain),&
                   d_B=object00_env%d_err(in_point,idomain),&
                   error=error)
             ENDIF
             object00_env%m_b(idomain)%mdata(irow+1,in_point+1)=trace0
             object00_env%m_b(idomain)%mdata(in_point+1,irow+1)=trace0
          ENDDO 
       ENDIF
    ENDDO 

    object00_env%in_point=object00_env%in_point+1
    IF (object00_env%in_point.gt.object00_env%max_buffer_length) object00_env%in_point=1
    
    CALL timestop(handle)

  END SUBROUTINE push_object00

! *****************************************************************************
!> \par History
!>       2011.12 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE extrapolate_object00(object00_env,extr_var,d_extr_var,error)
    TYPE(object00_type), INTENT(INOUT)       :: object00_env
    TYPE(cp_dbcsr_type), INTENT(INOUT), &
      OPTIONAL                               :: extr_var
    TYPE(object01_type), DIMENSION(:), &
      INTENT(INOUT), OPTIONAL                :: d_extr_var
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'extrapolate_object00', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, im, INFO, &
                                                LWORK, ndomains, unit_nr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: checksum
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues, WORK
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: coeff, m_b_copy, tmp1
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    IF (object00_env%object00_env_type.eq.object00_env_dbcsr) THEN
       IF ( .NOT.PRESENT(extr_var) ) THEN
          CPErrorMessage(cp_failure_level,routineP,"provide DBCSR matrix",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
    ELSE IF (object00_env%object00_env_type.eq.object00_env_domain) THEN
       IF ( .NOT.PRESENT(d_extr_var) ) THEN
          CPErrorMessage(cp_failure_level,routineP,"provide matrices",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
    ELSE
       CPErrorMessage(cp_failure_level,routineP,"illegal type",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    ALLOCATE(eigenvalues(object00_env%buffer_length+1))
    ALLOCATE(m_b_copy(object00_env%buffer_length+1,object00_env%buffer_length+1))

    ndomains=SIZE(object00_env%m_b)

    DO idomain=1, ndomains

       IF (object00_env%m_b(idomain)%domain.gt.0) THEN

          m_b_copy(:,:)=object00_env%m_b(idomain)%mdata(:,:)

          LWORK = -1
          ALLOCATE(WORK(MAX(1,LWORK)))
          CALL DSYEV('V','L',object00_env%buffer_length+1,m_b_copy,&
                  object00_env%buffer_length+1,eigenvalues,WORK,LWORK,INFO)
          LWORK = INT(WORK( 1 ))
          DEALLOCATE(WORK)
   
          ALLOCATE(WORK(MAX(1,LWORK)))
          CALL DSYEV('V','L',object00_env%buffer_length+1,m_b_copy,&
                  object00_env%buffer_length+1,eigenvalues,WORK,LWORK,INFO)
          IF( INFO.NE.0 ) THEN
             CPErrorMessage(cp_failure_level,routineP,"DSYEV failed",error)
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
          DEALLOCATE(WORK)

          ALLOCATE(tmp1(object00_env%buffer_length+1,1))
          ALLOCATE(coeff(object00_env%buffer_length+1,1))
          tmp1(:,1)=-1.0_dp*m_b_copy(1,:)/eigenvalues(:)
          coeff=MATMUL(m_b_copy,tmp1)
          DEALLOCATE(tmp1)

          checksum=0.0_dp
          IF (object00_env%object00_env_type.eq.object00_env_dbcsr) THEN
             CALL cp_dbcsr_set(extr_var,0.0_dp,error=error)
             DO im=1,object00_env%buffer_length
                CALL cp_dbcsr_add(extr_var,object00_env%m_var(im),&
                        1.0_dp,coeff(im+1,1),error=error)
                checksum=checksum+coeff(im+1,1)
             ENDDO
          ELSE IF (object00_env%object00_env_type.eq.object00_env_domain) THEN
             CALL copy_object01_gen(object00_env%d_var(1,idomain),&
                     d_extr_var(idomain),&
                     copy_data=.FALSE.,error=error)
             CALL set_object01_gen(d_extr_var(idomain),0.0_dp,error=error)
             DO im=1,object00_env%buffer_length
                CALL op2_object01_gen(1.0_dp,d_extr_var(idomain),&
                        coeff(im+1,1),object00_env%d_var(im,idomain),&
                        'N',error=error)
                checksum=checksum+coeff(im+1,1)
             ENDDO
          ENDIF

          DEALLOCATE(coeff)

       ENDIF 

    ENDDO 

    DEALLOCATE(eigenvalues)
    DEALLOCATE(m_b_copy)
    
    CALL timestop(handle)

  END SUBROUTINE extrapolate_object00

! *****************************************************************************
!> \par History
!>       2013.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  FUNCTION object00_overlap(object00_env,A,B,d_A,d_B,error)
    
    TYPE(object00_type), INTENT(INOUT)       :: object00_env
    TYPE(cp_dbcsr_type), INTENT(INOUT), &
      OPTIONAL                               :: A, B
    TYPE(object01_type), INTENT(INOUT), &
      OPTIONAL                               :: d_A, d_B
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    REAL(KIND=dp)                            :: object00_overlap

    CHARACTER(len=*), PARAMETER :: routineN = 'object00_overlap', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: trace

    CALL timeset(routineN,handle)

    IF (object00_env%object00_env_type.eq.object00_env_dbcsr) THEN
       IF ( .NOT.(PRESENT(A).AND.PRESENT(B)) ) THEN
          CPErrorMessage(cp_failure_level,routineP,"provide DBCSR matrices",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
    ELSE IF (object00_env%object00_env_type.eq.object00_env_domain) THEN
       IF ( .NOT.(PRESENT(d_A).AND.PRESENT(d_B)) ) THEN
          CPErrorMessage(cp_failure_level,routineP,"provide matrices",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
    ELSE
       CPErrorMessage(cp_failure_level,routineP,"illegal type",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    SELECT CASE(object00_env%error_type)
    CASE (object00_orthogonal)
       IF (object00_env%object00_env_type.eq.object00_env_dbcsr) THEN
          CALL cp_dbcsr_trace(A, B, trace, 'T', 'N', error=error)
       ELSE IF (object00_env%object00_env_type.eq.object00_env_domain) THEN
          CPPrecondition(SIZE(d_A%mdata,1).eq.SIZE(d_B%mdata,1),cp_failure_level,routineP,error,failure)
          CPPrecondition(SIZE(d_A%mdata,2).eq.SIZE(d_B%mdata,2),cp_failure_level,routineP,error,failure)
          CPPrecondition(d_A%domain.eq.d_B%domain,cp_failure_level,routineP,error,failure)
          CPPrecondition(d_A%domain.gt.0,cp_failure_level,routineP,error,failure)
          CPPrecondition(d_B%domain.gt.0,cp_failure_level,routineP,error,failure)
          trace=SUM(d_A%mdata(:,:)*d_B%mdata(:,:))
       ENDIF
    CASE DEFAULT
       CPErrorMessage(cp_failure_level,routineP,"Vector type is unknown",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    object00_overlap=trace

    CALL timestop(handle)

  END FUNCTION object00_overlap

! *****************************************************************************
!> \par History
!>       2011.12 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE release_object00(object00_env,error)
    TYPE(object00_type), INTENT(INOUT)       :: object00_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_object00', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, im

    CALL timeset(routineN,handle)

    DO im=1,object00_env%max_buffer_length
       IF (object00_env%object00_env_type.eq.object00_env_dbcsr) THEN
          CALL cp_dbcsr_release(object00_env%m_err(im),error=error)
          CALL cp_dbcsr_release(object00_env%m_var(im),error=error)
       ELSE IF (object00_env%object00_env_type.eq.object00_env_domain) THEN
          CALL release_object01_gen(object00_env%d_var(im,:),error=error)
          CALL release_object01_gen(object00_env%d_err(im,:),error=error)
       ENDIF
    ENDDO

    IF (object00_env%object00_env_type.eq.object00_env_domain) THEN
       CALL release_object01_gen(object00_env%m_b(:),error=error)
    ENDIF

    IF (ALLOCATED(object00_env%m_b)) DEALLOCATE(object00_env%m_b)
    IF (ALLOCATED(object00_env%m_err)) DEALLOCATE(object00_env%m_err)
    IF (ALLOCATED(object00_env%m_var)) DEALLOCATE(object00_env%m_var)
    IF (ALLOCATED(object00_env%d_err)) DEALLOCATE(object00_env%d_err)
    IF (ALLOCATED(object00_env%d_var)) DEALLOCATE(object00_env%d_var)

    CALL timestop(handle)

  END SUBROUTINE release_object00

END MODULE almo_scf_aux_types

