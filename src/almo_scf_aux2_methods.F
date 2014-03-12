!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>       2013.01 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
MODULE almo_scf_aux2_methods 

  USE almo_scf_aux2_types,             ONLY: object01_type,&
                                             object02_type,&
                                             select_row_col
  USE array_types,                     ONLY: array_data
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_col_block_sizes, cp_dbcsr_distribution, cp_dbcsr_filter, &
       cp_dbcsr_finalize, cp_dbcsr_get_block_p, cp_dbcsr_get_matrix_type, &
       cp_dbcsr_get_stored_coordinates, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_nblkcols_total, &
       cp_dbcsr_nblkrows_total, cp_dbcsr_reserve_block2d, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_work_create
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_iterator,&
                                             cp_dbcsr_type
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_mp,&
                                             dbcsr_mp_group,&
                                             dbcsr_mp_mynode,&
                                             dbcsr_mp_numnodes
  USE dbcsr_types,                     ONLY: dbcsr_type_antisymmetric,&
                                             dbcsr_type_no_symmetry,&
                                             dbcsr_type_symmetric
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_alltoall,&
                                             mp_sync
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf_aux2_methods'

  PUBLIC :: copy_object01_gen, copy_object01_data,&
            release_object01_gen, op1_object01_gen, op2_object01_gen,&
            construct_object01, init_object01_gen,&
            construct_original_form_object01,&
            op4_object01,&
            set_object01_gen,&
            op5_object01, ij_exists, op3_object01

  INTERFACE init_object01_gen
     MODULE PROCEDURE init_object01_0d 
     MODULE PROCEDURE init_object01_1d
     MODULE PROCEDURE init_object01_2d
  END INTERFACE
  
  INTERFACE set_object01_gen
     MODULE PROCEDURE set_object01_array
     MODULE PROCEDURE set_object01
  END INTERFACE
  
  INTERFACE copy_object01_gen
     MODULE PROCEDURE copy_object01_array
     MODULE PROCEDURE copy_object01
  END INTERFACE
  
  INTERFACE release_object01_gen
     MODULE PROCEDURE release_object01_array
     MODULE PROCEDURE release_object01
  END INTERFACE

  INTERFACE op1_object01_gen
     MODULE PROCEDURE op1_object01_once
     MODULE PROCEDURE op1_object01_array
  END INTERFACE

  INTERFACE op2_object01_gen
     MODULE PROCEDURE op2_object01_once
     MODULE PROCEDURE op2_object01_array
  END INTERFACE
  
CONTAINS 

  SUBROUTINE init_object01_0d(subm,error)

    TYPE(object01_type), INTENT(INOUT)       :: subm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    subm%domain=-1
    subm%nbrows=-1
    subm%nbcols=-1
    subm%nrows=-1
    subm%ncols=-1
    subm%nnodes=-1
    subm%groupid=-1

  END SUBROUTINE init_object01_0d
  
  SUBROUTINE init_object01_1d(subm,error)

    TYPE(object01_type), DIMENSION(:), &
      INTENT(INOUT)                          :: subm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    subm(:)%domain=-1
    subm(:)%nbrows=-1
    subm(:)%nbcols=-1
    subm(:)%nrows=-1
    subm(:)%ncols=-1
    subm(:)%nnodes=-1
    subm(:)%groupid=-1

  END SUBROUTINE init_object01_1d
  
  SUBROUTINE init_object01_2d(subm,error)

    TYPE(object01_type), DIMENSION(:, :), &
      INTENT(INOUT)                          :: subm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    subm(:,:)%domain=-1
    subm(:,:)%nbrows=-1
    subm(:,:)%nbcols=-1
    subm(:,:)%nrows=-1
    subm(:,:)%ncols=-1
    subm(:,:)%nnodes=-1
    subm(:,:)%groupid=-1

  END SUBROUTINE init_object01_2d
  
  SUBROUTINE copy_object01_array(original,copy,copy_data,error)

    TYPE(object01_type), DIMENSION(:), &
      INTENT(IN)                             :: original
    TYPE(object01_type), DIMENSION(:), &
      INTENT(INOUT)                          :: copy
    LOGICAL, INTENT(IN)                      :: copy_data
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'copy_object01_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, ndomains, &
                                                ndomainsB
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    ndomains=SIZE(original)
    ndomainsB=SIZE(copy)
    CPPrecondition(ndomains.eq.ndomainsB,cp_failure_level,routineP,error,failure)
    copy(:)%nnodes=original(:)%nnodes
    copy(:)%groupid=original(:)%groupid
    DO idomain = 1, ndomains
       IF (original(idomain)%domain.gt.0) THEN
          CALL copy_object01(original(idomain),copy(idomain),copy_data,error)
       ENDIF
    ENDDO 

    CALL timestop(handle)

  END SUBROUTINE copy_object01_array
  
  SUBROUTINE copy_object01(original,copy,copy_data,error)

    TYPE(object01_type), INTENT(IN)          :: original
    TYPE(object01_type), INTENT(INOUT)       :: copy
    LOGICAL, INTENT(IN)                      :: copy_data
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'copy_object01', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, icol, irow

    CALL timeset(routineN,handle)

    copy%domain=original%domain
    copy%nnodes=original%nnodes
    copy%groupid=original%groupid

    IF (original%domain.gt.0) THEN

       copy%nbrows=original%nbrows
       copy%nbcols=original%nbcols
       copy%nrows=original%nrows
       copy%ncols=original%ncols
   
       IF (.NOT.ALLOCATED(copy%dbcsr_row)) THEN
          ALLOCATE(copy%dbcsr_row(original%nbrows))
       ELSE
          IF (SIZE(copy%dbcsr_row).ne.SIZE(original%dbcsr_row)) THEN
             DEALLOCATE(copy%dbcsr_row)
             ALLOCATE(copy%dbcsr_row(original%nbrows))
          ENDIF
       ENDIF
       IF (.NOT.ALLOCATED(copy%dbcsr_col)) THEN
          ALLOCATE(copy%dbcsr_col(original%nbcols))
       ELSE
          IF (SIZE(copy%dbcsr_col).ne.SIZE(original%dbcsr_col)) THEN
             DEALLOCATE(copy%dbcsr_col)
             ALLOCATE(copy%dbcsr_col(original%nbcols))
          ENDIF
       ENDIF
       IF (.NOT.ALLOCATED(copy%size_brow)) THEN
          ALLOCATE(copy%size_brow(original%nbrows))
       ELSE
          IF (SIZE(copy%size_brow).ne.SIZE(original%size_brow)) THEN
             DEALLOCATE(copy%size_brow)
             ALLOCATE(copy%size_brow(original%nbrows))
          ENDIF
       ENDIF
       IF (.NOT.ALLOCATED(copy%size_bcol)) THEN
          ALLOCATE(copy%size_bcol(original%nbcols))
       ELSE
          IF (SIZE(copy%size_bcol).ne.SIZE(original%size_bcol)) THEN
             DEALLOCATE(copy%size_bcol)
             ALLOCATE(copy%size_bcol(original%nbcols))
          ENDIF
       ENDIF
   
       DO irow=1, original%nbrows
          copy%dbcsr_row(irow)=original%dbcsr_row(irow)
          copy%size_brow(irow)=original%size_brow(irow)
       ENDDO
   
       DO icol=1, original%nbcols
          copy%dbcsr_col(icol)=original%dbcsr_col(icol)
          copy%size_bcol(icol)=original%size_bcol(icol)
       ENDDO
   
       IF (copy_data) THEN
          CALL copy_object01_data(original%mdata,copy,error)
       ENDIF

    ENDIF 

    CALL timestop(handle)

  END SUBROUTINE copy_object01

  SUBROUTINE copy_object01_data(array,copy,error)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: array
    TYPE(object01_type), INTENT(INOUT)       :: copy
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'copy_object01_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ds1, ds2, handle, ms1, ms2
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    CPPrecondition(copy%domain.gt.0,cp_failure_level,routineP,error,failure)

    ds1=SIZE(array,1)
    ds2=SIZE(array,2)

    IF (.NOT.ALLOCATED(copy%mdata)) THEN
       ALLOCATE(copy%mdata(ds1,ds2))
    ELSE
       ms1=SIZE(copy%mdata,1)
       ms2=SIZE(copy%mdata,2)
       IF ((ds1.ne.ms1).OR.(ds2.ne.ms2)) THEN
          DEALLOCATE(copy%mdata)
          ALLOCATE(copy%mdata(ds1,ds2))
       ENDIF
    ENDIF

    copy%mdata(:,:)=array(:,:)

    CALL timestop(handle)

  END SUBROUTINE copy_object01_data
       
  SUBROUTINE set_object01_array(objects01,scalar,error)

    TYPE(object01_type), DIMENSION(:), &
      INTENT(INOUT)                          :: objects01
    REAL(KIND=dp), INTENT(IN)                :: scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_object01_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, ndomains

    CALL timeset(routineN,handle)

    ndomains=SIZE(objects01)
    DO idomain = 1, ndomains
       IF (objects01(idomain)%domain.gt.0) THEN
          CALL set_object01(objects01(idomain),scalar,error)
       ENDIF
    ENDDO 

    CALL timestop(handle)

  END SUBROUTINE set_object01_array
  
  SUBROUTINE set_object01(object01,scalar,error)

    TYPE(object01_type), INTENT(INOUT)       :: object01
    REAL(KIND=dp), INTENT(IN)                :: scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_object01', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ds1, ds2, handle, ms1, ms2
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    CPPrecondition(object01%domain.gt.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(object01%nrows.gt.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(object01%ncols.gt.0,cp_failure_level,routineP,error,failure)

    ds1=object01%nrows
    ds2=object01%ncols

    IF (.NOT.ALLOCATED(object01%mdata)) THEN
       ALLOCATE(object01%mdata(ds1,ds2))
    ELSE
       ms1=SIZE(object01%mdata,1)
       ms2=SIZE(object01%mdata,2)
       IF ((ds1.ne.ms1).OR.(ds2.ne.ms2)) THEN
          DEALLOCATE(object01%mdata)
          ALLOCATE(object01%mdata(ds1,ds2))
       ENDIF
    ENDIF

    object01%mdata(:,:)=scalar

    CALL timestop(handle)

  END SUBROUTINE set_object01

  SUBROUTINE release_object01_array(subm,error)

    TYPE(object01_type), DIMENSION(:), &
      INTENT(INOUT)                          :: subm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_object01_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, ndomains

    CALL timeset(routineN,handle)

    ndomains=SIZE(subm)
    DO idomain = 1, ndomains
       CALL release_object01(subm(idomain),error)
    ENDDO 

    CALL timestop(handle)

  END SUBROUTINE release_object01_array

  SUBROUTINE release_object01(subm,error)

    TYPE(object01_type), INTENT(INOUT)       :: subm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_object01', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    subm%domain=-1
    subm%nbrows=-1
    subm%nbcols=-1
    subm%nrows=-1
    subm%ncols=-1
    subm%nnodes=-1
    subm%groupid=-1
    
    IF (ALLOCATED(subm%dbcsr_row)) THEN
       DEALLOCATE(subm%dbcsr_row)
    ENDIF
    IF (ALLOCATED(subm%dbcsr_col)) THEN
       DEALLOCATE(subm%dbcsr_col)
    ENDIF
    IF (ALLOCATED(subm%size_brow)) THEN
       DEALLOCATE(subm%size_brow)
    ENDIF
    IF (ALLOCATED(subm%size_bcol)) THEN
       DEALLOCATE(subm%size_bcol)
    ENDIF
    IF (ALLOCATED(subm%mdata)) THEN
       DEALLOCATE(subm%mdata)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE release_object01

  SUBROUTINE op1_object01_once(transA,transB,alpha,A,B,beta,C,error)

    CHARACTER, INTENT(IN)                    :: transA, transB
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(object01_type), INTENT(IN)          :: A, B
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(object01_type), INTENT(INOUT)       :: C
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'op1_object01_once', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cs1, cs2, handle, icol, irow, &
                                                K, K1, LDA, LDB, LDC, M, &
                                                Mblocks, N, Nblocks
    LOGICAL                                  :: failure, NOTA, NOTB

    CALL timeset(routineN,handle)
    
    CPPrecondition(A%domain.gt.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(B%domain.gt.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(C%domain.gt.0,cp_failure_level,routineP,error,failure)

    LDA = SIZE(A%mdata,1)
    LDB = SIZE(B%mdata,1)

    NOTA = (transA.eq.'N').OR.(transA.eq.'n')
    NOTB = (transB.eq.'N').OR.(transB.eq.'n')

    IF (NOTA) THEN
       M = A%nrows
       K = A%ncols
       Mblocks = A%nbrows
    ELSE
       M = A%ncols
       K = A%nrows
       Mblocks = A%nbcols
    ENDIF

    IF (NOTB) THEN
       K1 = B%nrows
       N = B%ncols
       Nblocks = B%nbcols
    ELSE
       K1 = B%ncols
       N = B%nrows
       Nblocks = B%nbrows
    ENDIF

    CPPrecondition(K.eq.K1,cp_failure_level,routineP,error,failure)
   
    C%nrows=M
    C%ncols=N
    C%nbrows=Mblocks
    C%nbcols=Nblocks
    IF (ALLOCATED(C%dbcsr_row)) THEN
       DEALLOCATE(C%dbcsr_row)
    ENDIF
    ALLOCATE(C%dbcsr_row(C%nbrows))
    IF (ALLOCATED(C%dbcsr_col)) THEN
       DEALLOCATE(C%dbcsr_col)
    ENDIF
    ALLOCATE(C%dbcsr_col(C%nbcols))
    IF (ALLOCATED(C%size_brow)) THEN
       DEALLOCATE(C%size_brow)
    ENDIF
    ALLOCATE(C%size_brow(C%nbrows))
    IF (ALLOCATED(C%size_bcol)) THEN
       DEALLOCATE(C%size_bcol)
    ENDIF
    ALLOCATE(C%size_bcol(C%nbcols))

    DO irow=1, C%nbrows
       IF (NOTA) THEN
          C%dbcsr_row(irow)=A%dbcsr_row(irow)
          C%size_brow(irow)=A%size_brow(irow)
       ELSE
          C%dbcsr_row(irow)=A%dbcsr_col(irow)
          C%size_brow(irow)=A%size_bcol(irow)
       ENDIF
    ENDDO

    DO icol=1, C%nbcols
       IF (NOTB) THEN
          C%dbcsr_col(icol)=B%dbcsr_col(icol)
          C%size_bcol(icol)=B%size_bcol(icol)
       ELSE
          C%dbcsr_col(icol)=B%dbcsr_row(icol)
          C%size_bcol(icol)=B%size_brow(icol)
       ENDIF
    ENDDO

    IF (.NOT.ALLOCATED(C%mdata)) THEN
       CPPrecondition(beta.eq.0.0_dp,cp_failure_level,routineP,error,failure)
       ALLOCATE(C%mdata(C%nrows,C%ncols))
    ELSE
       cs1=SIZE(C%mdata,1)
       cs2=SIZE(C%mdata,2)
       IF ((C%nrows.ne.cs1).OR.(C%ncols.ne.cs2)) THEN
          CPPrecondition(beta.eq.0.0_dp,cp_failure_level,routineP,error,failure)
          DEALLOCATE(C%mdata)
          ALLOCATE(C%mdata(C%nrows,C%ncols))
       ENDIF
    ENDIF

    LDC = C%nrows

    CALL DGEMM(transA,transB,M,N,K,alpha,A%mdata,LDA,B%mdata,LDB,beta,C%mdata,LDC)

    C%nnodes=A%nnodes
    C%groupid=A%groupid

    CALL timestop(handle)

  END SUBROUTINE op1_object01_once

  SUBROUTINE op1_object01_array(transA,transB,alpha,A,B,beta,C,error)

    CHARACTER, INTENT(IN)                    :: transA, transB
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(object01_type), DIMENSION(:), &
      INTENT(IN)                             :: A, B
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(object01_type), DIMENSION(:), &
      INTENT(INOUT)                          :: C
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'op1_object01_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, idomainA, &
                                                idomainB, ndomains, &
                                                ndomainsB, ndomainsC
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    ndomains=SIZE(A)
    ndomainsB=SIZE(B)
    ndomainsC=SIZE(C)

    CPPrecondition(ndomains.eq.ndomainsB,cp_failure_level,routineP,error,failure)
    CPPrecondition(ndomainsB.eq.ndomainsC,cp_failure_level,routineP,error,failure)
    
    DO idomain = 1, ndomains

       idomainA = A(idomain)%domain
       idomainB = B(idomain)%domain
       
       CPPrecondition(idomainA.eq.idomainB,cp_failure_level,routineP,error,failure)
       
       C(idomain)%domain = idomainA

       IF (idomainA.gt.0) THEN
          CALL op1_object01_once(transA,transB,alpha,A(idomain),B(idomain),beta,C(idomain),error)
       ENDIF 

    ENDDO 

    CALL timestop(handle)

  END SUBROUTINE op1_object01_array

  SUBROUTINE op2_object01_once(alpha,A,beta,B,transB,error)

    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(object01_type), INTENT(INOUT)       :: A
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(object01_type), INTENT(IN)          :: B
    CHARACTER, INTENT(IN)                    :: transB
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'op2_object01_once', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: C1, C2, handle, icol, R1, R2
    LOGICAL                                  :: failure, NOTB

    CALL timeset(routineN,handle)

    CPPrecondition(A%domain.gt.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(B%domain.gt.0,cp_failure_level,routineP,error,failure)
    
    R1 = A%nrows
    C1 = A%ncols

    NOTB = (transB.eq.'N').OR.(transB.eq.'n')

    IF (NOTB) THEN
       R2 = B%nrows
       C2 = B%ncols
    ELSE
       R2 = B%ncols
       C2 = B%nrows
    ENDIF

    CPPrecondition(C1.eq.C2,cp_failure_level,routineP,error,failure)
    CPPrecondition(R1.eq.R2,cp_failure_level,routineP,error,failure)
   
    IF (NOTB) THEN
       DO icol = 1, C1
          A%mdata(:,icol) = alpha*A%mdata(:,icol) + beta*B%mdata(:,icol)
       ENDDO
    ELSE
       DO icol = 1, C1
          A%mdata(:,icol) = alpha*A%mdata(:,icol) + beta*B%mdata(icol,:)
       ENDDO
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE op2_object01_once

  SUBROUTINE op2_object01_array(alpha,A,beta,B,transB,error)

    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(object01_type), DIMENSION(:), &
      INTENT(INOUT)                          :: A
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(object01_type), DIMENSION(:), &
      INTENT(IN)                             :: B
    CHARACTER, INTENT(IN)                    :: transB
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'op2_object01_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, idomainA, &
                                                idomainB, ndomains, ndomainsB
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    ndomains=SIZE(A)
    ndomainsB=SIZE(B)

    CPPrecondition(ndomains.eq.ndomainsB,cp_failure_level,routineP,error,failure)
    
    DO idomain = 1, ndomains

       idomainA = A(idomain)%domain
       idomainB = B(idomain)%domain
       
       CPPrecondition(idomainA.eq.idomainB,cp_failure_level,routineP,error,failure)
       
       IF (idomainA.gt.0) THEN
          CALL op2_object01_once(alpha,A(idomain),beta,B(idomain),transB,error)
       ENDIF 

    ENDDO 

    CALL timestop(handle)

  END SUBROUTINE op2_object01_array

! *****************************************************************************
!> \par History
!>       2013.03 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE op3_object01(objects01,norm,error)
    
    TYPE(object01_type), DIMENSION(:), &
      INTENT(IN)                             :: objects01
    REAL(KIND=dp), INTENT(OUT)               :: norm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'op3_object01', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, ndomains
    REAL(KIND=dp)                            :: curr_norm, send_norm
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_norm

    CALL timeset(routineN,handle)

    send_norm=0.0_dp
    
    ndomains=SIZE(objects01)

    DO idomain = 1, ndomains

       IF (objects01(idomain)%domain.gt.0) THEN
          curr_norm=MAXVAL(ABS(objects01(idomain)%mdata))
          IF (curr_norm.gt.send_norm) send_norm=curr_norm
       ENDIF

    ENDDO 
    
    ALLOCATE(recv_norm(objects01(1)%nnodes))
    CALL mp_allgather(send_norm,recv_norm,objects01(1)%groupid)
   
    norm=MAXVAL(recv_norm)

    DEALLOCATE(recv_norm)

    CALL timestop(handle)
  
  END SUBROUTINE op3_object01 

! *****************************************************************************
!> \par History
!>       2013.03 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE op4_object01(A,B,trace,error)
    
    TYPE(object01_type), DIMENSION(:), &
      INTENT(IN)                             :: A, B
    REAL(KIND=dp), INTENT(OUT)               :: trace
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'op4_object01', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: domainA, domainB, handle, &
                                                idomain, ndomainsA, ndomainsB
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: curr_trace, send_trace
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_trace

    CALL timeset(routineN,handle)

    send_trace=0.0_dp
    
    ndomainsA=SIZE(A)
    ndomainsB=SIZE(B)

    CPPrecondition(ndomainsA.eq.ndomainsB,cp_failure_level,routineP,error,failure)
    
    DO idomain = 1, ndomainsA

       domainA = A(idomain)%domain
       domainB = B(idomain)%domain
       
       CPPrecondition(domainA.eq.domainB,cp_failure_level,routineP,error,failure)

       IF (domainA.gt.0) THEN

          CPPrecondition(A(idomain)%nrows.eq.B(idomain)%nrows,cp_failure_level,routineP,error,failure)
          CPPrecondition(A(idomain)%ncols.eq.B(idomain)%ncols,cp_failure_level,routineP,error,failure)
          
          curr_trace=SUM(A(idomain)%mdata(:,:)*B(idomain)%mdata(:,:))
          send_trace=send_trace+curr_trace

       ENDIF

    ENDDO 
    
    ALLOCATE(recv_trace(A(1)%nnodes))
    CALL mp_allgather(send_trace,recv_trace,A(1)%groupID)
   
    trace=SUM(recv_trace)

    DEALLOCATE(recv_trace)

    CALL timestop(handle)
  
  END SUBROUTINE op4_object01 

! *****************************************************************************
!> \par History
!>       2013.01 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE construct_object01(matrix,object01,distr_pattern,domain_map,&
     node_of_domain,job_type,error)
    
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(object01_type), DIMENSION(:), &
      INTENT(INOUT)                          :: object01
    TYPE(cp_dbcsr_type), INTENT(IN)          :: distr_pattern
    TYPE(object02_type), INTENT(IN)          :: domain_map
    INTEGER, DIMENSION(:), INTENT(IN)        :: node_of_domain
    INTEGER, INTENT(IN)                      :: job_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'construct_object01', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: matrix_type
    INTEGER :: block_node, block_offset, col, col_offset, col_size, &
      dest_node, GroupID, handle, iBlock, icol, idomain, index_col, index_ec, &
      index_er, index_row, index_sc, index_sr, iNode, ldesc, myNode, &
      nblkcols_tot, nblkrows_tot, ndomains, ndomains2, nNodes, &
      recv_size2_total, recv_size_total, row, row_size, send_size2_total, &
      send_size_total, smcol, smrow, start_data
    INTEGER, ALLOCATABLE, DIMENSION(:) :: first_col, first_row, &
      offset2_block, offset_block, recv_data2, recv_offset2_cpu, &
      recv_offset_cpu, recv_size2_cpu, recv_size_cpu, send_data2, &
      send_offset2_cpu, send_offset_cpu, send_size2_cpu, send_size_cpu
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_descriptor, &
                                                send_descriptor
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: failure, found, transp
    REAL(KIND=dp)                            :: antifactor
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_data, send_data
    REAL(KIND=dp), DIMENSION(:), POINTER     :: block_p

    CALL timeset(routineN,handle)
    
    nblkrows_tot = cp_dbcsr_nblkrows_total(matrix)
    nblkcols_tot = cp_dbcsr_nblkcols_total(matrix)

    ndomains = nblkcols_tot 

    nNodes = dbcsr_mp_numnodes(dbcsr_distribution_mp(&
       cp_dbcsr_distribution(distr_pattern)))
    GroupID = dbcsr_mp_group(dbcsr_distribution_mp(&
       cp_dbcsr_distribution(distr_pattern)))
    myNode = dbcsr_mp_mynode(dbcsr_distribution_mp(&
       cp_dbcsr_distribution(distr_pattern)))

    matrix_type=cp_dbcsr_get_matrix_type(matrix)

    ldesc=2
    ALLOCATE(send_descriptor(ldesc,nNodes))
    ALLOCATE(recv_descriptor(ldesc,nNodes))
    send_descriptor(:,:)=0
    
    DO idomain = 1, ndomains
             
       dest_node=node_of_domain(idomain)

       index_sr=1 
       IF (idomain.gt.1) index_sr=domain_map%index1(idomain-1)
       index_er=domain_map%index1(idomain)-1 

       DO index_row=index_sr, index_er

          row = domain_map%pairs(index_row,1)

          IF (job_type==select_row_col) THEN
             index_sc=1 
             IF (idomain.gt.1) index_sc=domain_map%index1(idomain-1)
             index_ec=domain_map%index1(idomain)-1 
          ELSE
             index_sc=1 
             index_ec=1 
          ENDIF

          DO index_col=index_sc, index_ec

             IF (job_type==select_row_col) THEN
                col = domain_map%pairs(index_col,1)
             ELSE
                col = idomain 
             ENDIF

             transp=.FALSE.
             CALL cp_dbcsr_get_stored_coordinates(matrix,&
                     row, col, transp, block_node)
             IF (block_node.eq.myNode) THEN
                CALL cp_dbcsr_get_block_p(matrix,row,col,block_p,found,row_size,col_size)
                IF (found) THEN
                   send_descriptor(1,dest_node+1)=send_descriptor(1,dest_node+1)+1
                   send_descriptor(2,dest_node+1)=send_descriptor(2,dest_node+1)+&
                      row_size*col_size
                ENDIF
             ENDIF

          ENDDO 

       ENDDO 

    ENDDO

    CALL mp_alltoall(send_descriptor,recv_descriptor,ldesc,GroupID)
    
    ALLOCATE(send_size_cpu(nNodes),send_offset_cpu(nNodes))
    send_offset_cpu(1)=0
    send_size_cpu(1)=send_descriptor(2,1)
    DO iNode=2,nNodes
       send_size_cpu(iNode)=send_descriptor(2,iNode)
       send_offset_cpu(iNode)=send_offset_cpu(iNode-1) + &
          send_size_cpu(iNode-1)
    ENDDO
    send_size_total=send_offset_cpu(nNodes)+send_size_cpu(nNodes)
    
    ALLOCATE(recv_size_cpu(nNodes),recv_offset_cpu(nNodes))
    recv_offset_cpu(1)=0
    recv_size_cpu(1)=recv_descriptor(2,1)
    DO iNode=2,nNodes
       recv_size_cpu(iNode)=recv_descriptor(2,iNode)
       recv_offset_cpu(iNode)=recv_offset_cpu(iNode-1) + &
          recv_size_cpu(iNode-1)
    ENDDO
    recv_size_total=recv_offset_cpu(nNodes)+recv_size_cpu(nNodes)
    
    ALLOCATE(send_size2_cpu(nNodes),send_offset2_cpu(nNodes))
    send_offset2_cpu(1)=0
    send_size2_cpu(1)=2*send_descriptor(1,1)
    DO iNode=2,nNodes
       send_size2_cpu(iNode)=2*send_descriptor(1,iNode)
       send_offset2_cpu(iNode)=send_offset2_cpu(iNode-1) + &
          send_size2_cpu(iNode-1)
    ENDDO
    send_size2_total=send_offset2_cpu(nNodes)+send_size2_cpu(nNodes)
    
    ALLOCATE(recv_size2_cpu(nNodes),recv_offset2_cpu(nNodes))
    recv_offset2_cpu(1)=0
    recv_size2_cpu(1)=2*recv_descriptor(1,1)
    DO iNode=2,nNodes
       recv_size2_cpu(iNode)=2*recv_descriptor(1,iNode)
       recv_offset2_cpu(iNode)=recv_offset2_cpu(iNode-1) + &
          recv_size2_cpu(iNode-1)
    ENDDO
    recv_size2_total=recv_offset2_cpu(nNodes)+recv_size2_cpu(nNodes)
    
    DEALLOCATE(send_descriptor)
    DEALLOCATE(recv_descriptor)
    
    ALLOCATE(send_data(send_size_total))
    ALLOCATE(recv_data(recv_size_total))
    ALLOCATE(send_data2(send_size2_total))
    ALLOCATE(recv_data2(recv_size2_total))
    ALLOCATE(offset_block(nNodes))
    ALLOCATE(offset2_block(nNodes))
    offset_block(:)=0
    offset2_block(:)=0
    DO idomain = 1, ndomains
             
       dest_node=node_of_domain(idomain)

       index_sr=1 
       IF (idomain.gt.1) index_sr=domain_map%index1(idomain-1)
       index_er=domain_map%index1(idomain)-1 

       DO index_row=index_sr, index_er

          row = domain_map%pairs(index_row,1)

          IF (job_type==select_row_col) THEN
             index_sc=1 
             IF (idomain.gt.1) index_sc=domain_map%index1(idomain-1)
             index_ec=domain_map%index1(idomain)-1 
          ELSE
             index_sc=1 
             index_ec=1 
          ENDIF

          DO index_col=index_sc, index_ec

             IF (job_type==select_row_col) THEN
                col = domain_map%pairs(index_col,1)
             ELSE
                col = idomain 
             ENDIF

             transp=.FALSE.
             CALL cp_dbcsr_get_stored_coordinates(matrix,&
                     row, col, transp, block_node)
             IF (block_node.eq.myNode) THEN
                CALL cp_dbcsr_get_block_p(matrix,row,col,block_p,found,row_size,col_size)
                IF (found) THEN
                   col_offset=row_size*col_size
                   start_data=send_offset_cpu(dest_node+1)+&
                                 offset_block(dest_node+1)
                   send_data(start_data+1:start_data+col_offset)=&
                         block_p(1:col_offset)
                   offset_block(dest_node+1)=offset_block(dest_node+1)+col_offset
                   send_data2(send_offset2_cpu(dest_node+1)+&
                         offset2_block(dest_node+1)+1)=row
                   send_data2(send_offset2_cpu(dest_node+1)+&
                         offset2_block(dest_node+1)+2)=col
                   offset2_block(dest_node+1)=offset2_block(dest_node+1)+2
                ENDIF
             ENDIF

          ENDDO 

       ENDDO 

    ENDDO
    
    CALL mp_alltoall(send_data,send_size_cpu,send_offset_cpu,&
            recv_data,recv_size_cpu,recv_offset_cpu,GroupID)
    CALL mp_alltoall(send_data2,send_size2_cpu,send_offset2_cpu,&
            recv_data2,recv_size2_cpu,recv_offset2_cpu,GroupID)
    
    DEALLOCATE(send_size_cpu,send_offset_cpu)
    DEALLOCATE(send_size2_cpu,send_offset2_cpu)
    DEALLOCATE(send_data)
    DEALLOCATE(send_data2)
    DEALLOCATE(offset_block)
    DEALLOCATE(offset2_block)
    
    col_blk_size => array_data(cp_dbcsr_col_block_sizes(matrix))
    row_blk_size => array_data(cp_dbcsr_row_block_sizes(matrix))
    ndomains2=SIZE(object01)
    IF (ndomains2.ne.ndomains) THEN
       CPErrorMessage(cp_failure_level,routineP,"wrong matrix size",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF
    CALL release_object01_gen(object01,error)
    object01(:)%nnodes=nNodes
    object01(:)%groupID=GroupID
    object01(:)%nrows=0
    object01(:)%ncols=0
    
    ALLOCATE(first_row(nblkrows_tot),first_col(nblkcols_tot))
    object01(:)%domain=-1
    DO idomain=1, ndomains 
       dest_node=node_of_domain(idomain)
       IF (dest_node.eq.mynode) THEN
          object01(idomain)%domain=idomain
          object01(idomain)%nbrows=0
          object01(idomain)%nbcols=0
          
          first_row(:)=-1
          index_sr=1 
          IF (idomain.gt.1) index_sr=domain_map%index1(idomain-1)
          index_er=domain_map%index1(idomain)-1 
          DO index_row=index_sr, index_er
             row = domain_map%pairs(index_row,1)
             first_row(row)=object01(idomain)%nrows+1
             object01(idomain)%nrows=object01(idomain)%nrows+row_blk_size(row)
             object01(idomain)%nbrows=object01(idomain)%nbrows+1
          ENDDO
          ALLOCATE(object01(idomain)%dbcsr_row(object01(idomain)%nbrows))
          ALLOCATE(object01(idomain)%size_brow(object01(idomain)%nbrows))
          smrow=1
          index_sr=1 
          IF (idomain.gt.1) index_sr=domain_map%index1(idomain-1)
          index_er=domain_map%index1(idomain)-1 
          DO index_row=index_sr, index_er
             row = domain_map%pairs(index_row,1)
                object01(idomain)%dbcsr_row(smrow)=row
                object01(idomain)%size_brow(smrow)=row_blk_size(row)
                smrow=smrow+1
          ENDDO

          first_col(:)=-1
          IF (job_type==select_row_col) THEN
             index_sc=1 
             IF (idomain.gt.1) index_sc=domain_map%index1(idomain-1)
             index_ec=domain_map%index1(idomain)-1 
          ELSE
             index_sc=1 
             index_ec=1
          ENDIF
          DO index_col=index_sc, index_ec
             IF (job_type==select_row_col) THEN
                col = domain_map%pairs(index_col,1)
             ELSE
                col = idomain 
             ENDIF
                first_col(col)=object01(idomain)%ncols+1
                object01(idomain)%ncols=object01(idomain)%ncols+col_blk_size(col)
                object01(idomain)%nbcols=object01(idomain)%nbcols+1
          ENDDO
          
          ALLOCATE(object01(idomain)%dbcsr_col(object01(idomain)%nbcols))
          ALLOCATE(object01(idomain)%size_bcol(object01(idomain)%nbcols))
          
          smcol=1
          IF (job_type==select_row_col) THEN
             index_sc=1 
             IF (idomain.gt.1) index_sc=domain_map%index1(idomain-1)
             index_ec=domain_map%index1(idomain)-1 
          ELSE
             index_sc=1 
             index_ec=1 
          ENDIF
          DO index_col=index_sc, index_ec
             IF (job_type==select_row_col) THEN
                col = domain_map%pairs(index_col,1)
             ELSE
                col = idomain 
             ENDIF
                object01(idomain)%dbcsr_col(smcol)=col
                object01(idomain)%size_bcol(smcol)=col_blk_size(col)
                smcol=smcol+1
          ENDDO

          ALLOCATE(object01(idomain)%mdata(&
             object01(idomain)%nrows,&
             object01(idomain)%ncols))
          object01(idomain)%mdata(:,:)=0.0_dp
          DO iNode=1, nNodes
             block_offset=0
             DO iBlock=1, recv_size2_cpu(iNode)/2
                row = recv_data2(recv_offset2_cpu(iNode)+(iBlock-1)*2+1)
                col = recv_data2(recv_offset2_cpu(iNode)+(iBlock-1)*2+2)
                IF ((first_col(col).NE.-1).AND.(first_row(row).NE.-1)) THEN
                   start_data=recv_offset_cpu(iNode)+block_offset+1
                   DO icol=0, col_blk_size(col)-1
                      object01(idomain)%mdata(first_row(row):&
                                first_row(row)+row_blk_size(row)-1,&
                                first_col(col)+icol)=&
                                recv_data(start_data:start_data+row_blk_size(row)-1)
                      start_data=start_data+row_blk_size(row)
                   ENDDO
                   IF (job_type==select_row_col) THEN
                      IF (matrix_type==dbcsr_type_symmetric.OR.&
                          matrix_type==dbcsr_type_antisymmetric) THEN
                         antifactor=1.0_dp
                         IF (matrix_type==dbcsr_type_antisymmetric) THEN
                            antifactor=-1.0_dp
                         ENDIF
                         start_data=recv_offset_cpu(iNode)+block_offset+1
                         DO icol=0, col_blk_size(col)-1
                            object01(idomain)%mdata(first_row(col)+icol,&
                                      first_col(row):&
                                      first_col(row)+row_blk_size(row)-1)=&
                                      antifactor*recv_data(start_data:&
                                      start_data+row_blk_size(row)-1)
                            start_data=start_data+row_blk_size(row)
                         ENDDO
                      ELSE IF (matrix_type==dbcsr_type_no_symmetry) THEN
                      ELSE
                         CPErrorMessage(cp_failure_level,routineP,"matrix type is NYI",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                      ENDIF
                   ENDIF
                ENDIF
                block_offset=block_offset+col_blk_size(col)*row_blk_size(row)
             ENDDO
          ENDDO
       ENDIF 
    ENDDO 

    DEALLOCATE(recv_size_cpu,recv_offset_cpu)
    DEALLOCATE(recv_size2_cpu,recv_offset2_cpu)
    DEALLOCATE(recv_data)
    DEALLOCATE(recv_data2)
    DEALLOCATE(first_row,first_col)
       
    CALL timestop(handle)
  
  END SUBROUTINE construct_object01 

! *****************************************************************************
!> \par History
!>       2013.01 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE construct_original_form_object01(matrix,object01,distr_pattern,&
     error)
    
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(object01_type), DIMENSION(:), &
      INTENT(IN)                             :: object01
    TYPE(cp_dbcsr_type), INTENT(IN)          :: distr_pattern
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'construct_original_form_object01', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: matrix_type
    INTEGER :: block_offset, col, col_offset, col_s, colsize, dest_node, &
      GroupID, handle, iBlock, icol, idomain, iNode, irow_subm, ldesc, &
      myNode, nblkcols_tot, nblkrows_tot, ndomains, ndomains2, nNodes, &
      recv_size2_total, recv_size_total, row, row_s, rowsize, &
      send_size2_total, send_size_total, smroff, start_data, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:) :: offset2_block, offset_block, &
      recv_data2, recv_offset2_cpu, recv_offset_cpu, recv_size2_cpu, &
      recv_size_cpu, send_data2, send_offset2_cpu, send_offset_cpu, &
      send_size2_cpu, send_size_cpu
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_descriptor, &
                                                send_descriptor
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: failure, transp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_data, send_data
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_p
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    nblkrows_tot = cp_dbcsr_nblkrows_total(matrix)
    nblkcols_tot = cp_dbcsr_nblkcols_total(matrix)

    ndomains = nblkcols_tot 
    ndomains2 = SIZE(object01)

    IF (ndomains.ne.ndomains2) THEN
       CPErrorMessage(cp_failure_level,routineP,"domain mismatch",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    nNodes = dbcsr_mp_numnodes(dbcsr_distribution_mp(&
       cp_dbcsr_distribution(distr_pattern)))
    GroupID = dbcsr_mp_group(dbcsr_distribution_mp(&
       cp_dbcsr_distribution(distr_pattern)))
    myNode = dbcsr_mp_mynode(dbcsr_distribution_mp(&
       cp_dbcsr_distribution(distr_pattern)))

    matrix_type=cp_dbcsr_get_matrix_type(matrix)
    IF (matrix_type.ne.dbcsr_type_no_symmetry) THEN
       CPErrorMessage(cp_failure_level,routineP,"only non-symmetric matrices so far",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF
   
    CALL cp_dbcsr_iterator_start(iter,matrix)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter,row,col,block_p)
       block_p(:,:)=0.0_dp
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)
    CALL cp_dbcsr_filter(matrix,0.1_dp,error=error)

    CALL cp_dbcsr_work_create(matrix,work_mutable=.TRUE.,&
            error=error)

    ldesc=2
    ALLOCATE(send_descriptor(ldesc,nNodes))
    ALLOCATE(recv_descriptor(ldesc,nNodes))
    send_descriptor(:,:)=0
    
    DO idomain = 1, ndomains

       IF (object01(idomain)%domain.gt.0) THEN

          DO irow_subm = 1, object01(idomain)%nbrows
          
             IF (object01(idomain)%nbcols.ne.1) THEN
                CPErrorMessage(cp_failure_level,routineP,"corrupt object01 structure",error)
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             ENDIF

             row = object01(idomain)%dbcsr_row(irow_subm)
             col = object01(idomain)%dbcsr_col(1)

             IF (col.ne.idomain) THEN
                CPErrorMessage(cp_failure_level,routineP,"corrupt object01 structure",error)
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             ENDIF

             transp=.FALSE.
             row_s = row
             col_s = idomain
             CALL cp_dbcsr_get_stored_coordinates(distr_pattern,&
                     row_s, col_s, transp, dest_node)
             CPPrecondition(row_s==row,cp_failure_level,routineP,error,failure)
             CPPrecondition(col_s==idomain,cp_failure_level,routineP,error,failure)
             
             send_descriptor(1,dest_node+1)=send_descriptor(1,dest_node+1)+1
             send_descriptor(2,dest_node+1)=send_descriptor(2,dest_node+1)+&
                object01(idomain)%size_brow(irow_subm)*&
                object01(idomain)%size_bcol(1)

          ENDDO 

       ENDIF       

    ENDDO 

    CALL mp_alltoall(send_descriptor,recv_descriptor,ldesc,GroupID)
    
    ALLOCATE(send_size_cpu(nNodes),send_offset_cpu(nNodes))
    send_offset_cpu(1)=0
    send_size_cpu(1)=send_descriptor(2,1)
    DO iNode=2,nNodes
       send_size_cpu(iNode)=send_descriptor(2,iNode)
       send_offset_cpu(iNode)=send_offset_cpu(iNode-1) + &
          send_size_cpu(iNode-1)
    ENDDO
    send_size_total=send_offset_cpu(nNodes)+send_size_cpu(nNodes)
    
    ALLOCATE(recv_size_cpu(nNodes),recv_offset_cpu(nNodes))
    recv_offset_cpu(1)=0
    recv_size_cpu(1)=recv_descriptor(2,1)
    DO iNode=2,nNodes
       recv_size_cpu(iNode)=recv_descriptor(2,iNode)
       recv_offset_cpu(iNode)=recv_offset_cpu(iNode-1) + &
          recv_size_cpu(iNode-1)
    ENDDO
    recv_size_total=recv_offset_cpu(nNodes)+recv_size_cpu(nNodes)
    
    ALLOCATE(send_size2_cpu(nNodes),send_offset2_cpu(nNodes))
    send_offset2_cpu(1)=0
    send_size2_cpu(1)=2*send_descriptor(1,1)
    DO iNode=2,nNodes
       send_size2_cpu(iNode)=2*send_descriptor(1,iNode)
       send_offset2_cpu(iNode)=send_offset2_cpu(iNode-1) + &
          send_size2_cpu(iNode-1)
    ENDDO
    send_size2_total=send_offset2_cpu(nNodes)+send_size2_cpu(nNodes)
    
    ALLOCATE(recv_size2_cpu(nNodes),recv_offset2_cpu(nNodes))
    recv_offset2_cpu(1)=0
    recv_size2_cpu(1)=2*recv_descriptor(1,1)
    DO iNode=2,nNodes
       recv_size2_cpu(iNode)=2*recv_descriptor(1,iNode)
       recv_offset2_cpu(iNode)=recv_offset2_cpu(iNode-1) + &
          recv_size2_cpu(iNode-1)
    ENDDO
    recv_size2_total=recv_offset2_cpu(nNodes)+recv_size2_cpu(nNodes)
    
    DEALLOCATE(send_descriptor)
    DEALLOCATE(recv_descriptor)
    
    ALLOCATE(send_data(send_size_total))
    ALLOCATE(recv_data(recv_size_total))
    ALLOCATE(send_data2(send_size2_total))
    ALLOCATE(recv_data2(recv_size2_total))
    ALLOCATE(offset_block(nNodes))
    ALLOCATE(offset2_block(nNodes))
    offset_block(:)=0
    offset2_block(:)=0
    DO idomain = 1, ndomains

       IF (object01(idomain)%domain.gt.0) THEN

          smroff=0
          DO irow_subm = 1, object01(idomain)%nbrows
          
             row = object01(idomain)%dbcsr_row(irow_subm)
             col = object01(idomain)%dbcsr_col(1)

             rowsize = object01(idomain)%size_brow(irow_subm)
             colsize = object01(idomain)%size_bcol(1)

             transp=.FALSE.
             row_s = row
             col_s = idomain
             CALL cp_dbcsr_get_stored_coordinates(distr_pattern,&
                     row_s, col_s, transp, dest_node)
             CPPrecondition(row_s==row,cp_failure_level,routineP,error,failure)
             CPPrecondition(col_s==idomain,cp_failure_level,routineP,error,failure)
            
             col_offset=0
             DO icol=1,colsize
                start_data=send_offset_cpu(dest_node+1)+&
                           offset_block(dest_node+1)+&
                           col_offset
                send_data(start_data+1:start_data+rowsize)=&
                   object01(idomain)%mdata(smroff+1:smroff+rowsize,icol)
                col_offset=col_offset+rowsize
             ENDDO
             offset_block(dest_node+1)=offset_block(dest_node+1)+&
                colsize*rowsize
             send_data2(send_offset2_cpu(dest_node+1)+&
                   offset2_block(dest_node+1)+1)=row
             send_data2(send_offset2_cpu(dest_node+1)+&
                   offset2_block(dest_node+1)+2)=col
             offset2_block(dest_node+1)=offset2_block(dest_node+1)+2

             smroff=smroff+rowsize

          ENDDO

       ENDIF       

    ENDDO 
    
    CALL mp_alltoall(send_data,send_size_cpu,send_offset_cpu,&
            recv_data,recv_size_cpu,recv_offset_cpu,GroupID)
    CALL mp_alltoall(send_data2,send_size2_cpu,send_offset2_cpu,&
            recv_data2,recv_size2_cpu,recv_offset2_cpu,GroupID)
    
    DEALLOCATE(send_size_cpu,send_offset_cpu)
    DEALLOCATE(send_size2_cpu,send_offset2_cpu)
    DEALLOCATE(send_data)
    DEALLOCATE(send_data2)
    DEALLOCATE(offset_block)
    DEALLOCATE(offset2_block)
    
    col_blk_size => array_data(cp_dbcsr_col_block_sizes(matrix))
    row_blk_size => array_data(cp_dbcsr_row_block_sizes(matrix))
    DO iNode=1, nNodes
       block_offset=0
       DO iBlock=1, recv_size2_cpu(iNode)/2
          row = recv_data2(recv_offset2_cpu(iNode)+(iBlock-1)*2+1)
          col = recv_data2(recv_offset2_cpu(iNode)+(iBlock-1)*2+2)
          NULLIFY (block_p)
          CALL cp_dbcsr_reserve_block2d(matrix,row,col,block_p)
          CPPostcondition(ASSOCIATED(block_p),cp_failure_level,routineP,error,failure)
          start_data=recv_offset_cpu(iNode)+block_offset+1
          DO icol=1, col_blk_size(col)
             block_p(:,icol)=&
                recv_data(start_data:start_data+row_blk_size(row)-1)
             start_data=start_data+row_blk_size(row)
          ENDDO
          block_offset=block_offset+col_blk_size(col)*row_blk_size(row)
       ENDDO
    ENDDO

    DEALLOCATE(recv_size_cpu,recv_offset_cpu)
    DEALLOCATE(recv_size2_cpu,recv_offset2_cpu)
    DEALLOCATE(recv_data)
    DEALLOCATE(recv_data2)
       
    CALL cp_dbcsr_finalize(matrix,error=error)
    
    CALL timestop(handle)
  
  END SUBROUTINE construct_original_form_object01 

  SUBROUTINE op5_object01(objects01,mpgroup,error)
    
    TYPE(object01_type), DIMENSION(:), &
      INTENT(IN)                             :: objects01
    INTEGER, INTENT(IN)                      :: mpgroup
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'op5_object01', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=30)                        :: colstr, formatstr
    INTEGER                                  :: handle, i, irow, n, ncols, &
                                                nrows

    CALL timeset(routineN,handle)
   
    n=SIZE(objects01)

    DO i=1, n
       nrows=SIZE(objects01(i)%mdata,1)
       ncols=SIZE(objects01(i)%mdata,2)
       WRITE (colstr,*) ncols
       formatstr='('//TRIM(ADJUSTL(colstr))//'F16.9)'
       IF (objects01(i)%domain.gt.0) THEN
          WRITE (*,*) "object01: ", i, nrows, 'x', ncols
          nrows=SIZE(objects01(i)%mdata,1)
          DO irow=1, nrows
             WRITE (*,formatstr) objects01(i)%mdata(irow,:)
          ENDDO
       ENDIF
       CALL mp_sync(mpgroup)
    ENDDO

    CALL timestop(handle)
  
  END SUBROUTINE op5_object01 

! *****************************************************************************
!> \par History
!>       2013.01 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  FUNCTION ij_exists(map,row,col,error)
  
    TYPE(object02_type), INTENT(IN)          :: map
    INTEGER, INTENT(IN)                      :: row, col
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL                                  :: ij_exists

    CHARACTER(len=*), PARAMETER :: routineN = 'ij_exists', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first, last, mid, ndomains

    ndomains=SIZE(map%index1)

    ij_exists=.FALSE.
    IF (col.lt.1.OR.col.gt.ndomains) RETURN
    first=1
    IF (col.gt.1) first=map%index1(col-1)
    last=map%index1(col)-1

    DO WHILE (last.ge.first)
        mid = first+(last-first)/2
        IF (map%pairs(mid,1).gt.row) THEN
           last = mid-1
        ELSE IF (map%pairs(mid,1).lt.row) THEN
           first = mid+1
        ELSE
            ij_exists=.TRUE.   
            EXIT
        END IF
    END DO
  
    RETURN

  END FUNCTION ij_exists

END MODULE almo_scf_aux2_methods 


