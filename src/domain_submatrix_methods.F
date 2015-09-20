!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Subroutines to handle submatrices 
!> \par History
!>       2013.01 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
MODULE domain_submatrix_methods

  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_col_block_sizes, cp_dbcsr_distribution, cp_dbcsr_filter, &
       cp_dbcsr_finalize, cp_dbcsr_get_block_p, cp_dbcsr_get_matrix_type, &
       cp_dbcsr_get_stored_coordinates, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, &
       cp_dbcsr_nblkcols_total, cp_dbcsr_nblkrows_total, &
       cp_dbcsr_reserve_block2d, cp_dbcsr_row_block_sizes, cp_dbcsr_type, &
       cp_dbcsr_work_create, dbcsr_distribution_mp, dbcsr_mp_group, &
       dbcsr_mp_mynode, dbcsr_mp_numnodes, dbcsr_type_antisymmetric, &
       dbcsr_type_no_symmetry, dbcsr_type_symmetric
  USE domain_submatrix_types,          ONLY: domain_map_type,&
                                             domain_submatrix_type,&
                                             select_row,&
                                             select_row_col
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_alltoall,&
                                             mp_sync
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'domain_submatrix_methods'

  PUBLIC :: copy_submatrices, copy_submatrix_data,&
            release_submatrices, multiply_submatrices, add_submatrices,&
            construct_submatrices, init_submatrices,&
            construct_dbcsr_from_submatrices,&
            set_submatrices,&
            print_submatrices, maxnorm_submatrices

  INTERFACE init_submatrices
     MODULE PROCEDURE init_submatrices_0d 
     MODULE PROCEDURE init_submatrices_1d
     MODULE PROCEDURE init_submatrices_2d
  END INTERFACE
  
  INTERFACE set_submatrices
     MODULE PROCEDURE set_submatrix_array
     MODULE PROCEDURE set_submatrix
  END INTERFACE
  
  INTERFACE copy_submatrices
     MODULE PROCEDURE copy_submatrix_array
     MODULE PROCEDURE copy_submatrix
  END INTERFACE
  
  INTERFACE release_submatrices
     MODULE PROCEDURE release_submatrix_array
     MODULE PROCEDURE release_submatrix
  END INTERFACE

  INTERFACE multiply_submatrices
     MODULE PROCEDURE multiply_submatrices_once
     MODULE PROCEDURE multiply_submatrices_array
  END INTERFACE

  INTERFACE add_submatrices
     MODULE PROCEDURE add_submatrices_once
     MODULE PROCEDURE add_submatrices_array
  END INTERFACE
  
CONTAINS 

! *****************************************************************************
!> \brief ...
!> \param subm ...
! *****************************************************************************
  SUBROUTINE init_submatrices_0d(subm)

    TYPE(domain_submatrix_type), &
      INTENT(INOUT)                          :: subm

    subm%domain=-1
    subm%nbrows=-1
    subm%nbcols=-1
    subm%nrows=-1
    subm%ncols=-1
    subm%nnodes=-1
    subm%groupid=-1

  END SUBROUTINE init_submatrices_0d
  
! *****************************************************************************
!> \brief ...
!> \param subm ...
! *****************************************************************************
  SUBROUTINE init_submatrices_1d(subm)

    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(INOUT)            :: subm

    subm(:)%domain=-1
    subm(:)%nbrows=-1
    subm(:)%nbcols=-1
    subm(:)%nrows=-1
    subm(:)%ncols=-1
    subm(:)%nnodes=-1
    subm(:)%groupid=-1

  END SUBROUTINE init_submatrices_1d
  
! *****************************************************************************
!> \brief ...
!> \param subm ...
! *****************************************************************************
  SUBROUTINE init_submatrices_2d(subm)

    TYPE(domain_submatrix_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: subm

    subm(:,:)%domain=-1
    subm(:,:)%nbrows=-1
    subm(:,:)%nbcols=-1
    subm(:,:)%nrows=-1
    subm(:,:)%ncols=-1
    subm(:,:)%nnodes=-1
    subm(:,:)%groupid=-1

  END SUBROUTINE init_submatrices_2d
  
! *****************************************************************************
!> \brief ...
!> \param original ...
!> \param copy ...
!> \param copy_data ...
! *****************************************************************************
  SUBROUTINE copy_submatrix_array(original,copy,copy_data)

    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(IN)               :: original
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(INOUT)            :: copy
    LOGICAL, INTENT(IN)                      :: copy_data

    CHARACTER(len=*), PARAMETER :: routineN = 'copy_submatrix_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, ndomains, &
                                                ndomainsB

    CALL timeset(routineN,handle)

    ndomains=SIZE(original)
    ndomainsB=SIZE(copy)
    CPASSERT(ndomains.eq.ndomainsB)
    copy(:)%nnodes=original(:)%nnodes
    copy(:)%groupid=original(:)%groupid
    DO idomain = 1, ndomains
       IF (original(idomain)%domain.gt.0) THEN
          CALL copy_submatrix(original(idomain),copy(idomain),copy_data)
       ENDIF
    ENDDO ! loop over domains

    CALL timestop(handle)

  END SUBROUTINE copy_submatrix_array
  
! *****************************************************************************
!> \brief ...
!> \param original ...
!> \param copy ...
!> \param copy_data ...
! *****************************************************************************
  SUBROUTINE copy_submatrix(original,copy,copy_data)

    TYPE(domain_submatrix_type), INTENT(IN)  :: original
    TYPE(domain_submatrix_type), &
      INTENT(INOUT)                          :: copy
    LOGICAL, INTENT(IN)                      :: copy_data

    CHARACTER(len=*), PARAMETER :: routineN = 'copy_submatrix', &
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
          CALL copy_submatrix_data(original%mdata,copy)
       ENDIF

    ENDIF ! do not copy empty submatrix

    CALL timestop(handle)

  END SUBROUTINE copy_submatrix

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param copy ...
! *****************************************************************************
  SUBROUTINE copy_submatrix_data(array,copy)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: array
    TYPE(domain_submatrix_type), &
      INTENT(INOUT)                          :: copy

    CHARACTER(len=*), PARAMETER :: routineN = 'copy_submatrix_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ds1, ds2, handle, ms1, ms2

    CALL timeset(routineN,handle)

    CPASSERT(copy%domain.gt.0)

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

  END SUBROUTINE copy_submatrix_data
       
! *****************************************************************************
!> \brief ...
!> \param submatrices ...
!> \param scalar ...
! *****************************************************************************
  SUBROUTINE set_submatrix_array(submatrices,scalar)

    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(INOUT)            :: submatrices
    REAL(KIND=dp), INTENT(IN)                :: scalar

    CHARACTER(len=*), PARAMETER :: routineN = 'set_submatrix_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, ndomains

    CALL timeset(routineN,handle)

    ndomains=SIZE(submatrices)
    DO idomain = 1, ndomains
       IF (submatrices(idomain)%domain.gt.0) THEN
          CALL set_submatrix(submatrices(idomain),scalar)
       ENDIF
    ENDDO ! loop over domains

    CALL timestop(handle)

  END SUBROUTINE set_submatrix_array
  
! *****************************************************************************
!> \brief ...
!> \param submatrix ...
!> \param scalar ...
! *****************************************************************************
  SUBROUTINE set_submatrix(submatrix,scalar)

    TYPE(domain_submatrix_type), &
      INTENT(INOUT)                          :: submatrix
    REAL(KIND=dp), INTENT(IN)                :: scalar

    CHARACTER(len=*), PARAMETER :: routineN = 'set_submatrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ds1, ds2, handle, ms1, ms2

    CALL timeset(routineN,handle)

    CPASSERT(submatrix%domain.gt.0)
    CPASSERT(submatrix%nrows.gt.0)
    CPASSERT(submatrix%ncols.gt.0)

    ds1=submatrix%nrows
    ds2=submatrix%ncols

    IF (.NOT.ALLOCATED(submatrix%mdata)) THEN
       ALLOCATE(submatrix%mdata(ds1,ds2))
    ELSE
       ms1=SIZE(submatrix%mdata,1)
       ms2=SIZE(submatrix%mdata,2)
       IF ((ds1.ne.ms1).OR.(ds2.ne.ms2)) THEN
          DEALLOCATE(submatrix%mdata)
          ALLOCATE(submatrix%mdata(ds1,ds2))
       ENDIF
    ENDIF

    submatrix%mdata(:,:)=scalar

    CALL timestop(handle)

  END SUBROUTINE set_submatrix

! *****************************************************************************
!> \brief ...
!> \param subm ...
! *****************************************************************************
  SUBROUTINE release_submatrix_array(subm)

    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(INOUT)            :: subm

    CHARACTER(len=*), PARAMETER :: routineN = 'release_submatrix_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, ndomains

    CALL timeset(routineN,handle)

    ndomains=SIZE(subm)
    DO idomain = 1, ndomains
       CALL release_submatrix(subm(idomain))
    ENDDO ! loop over domains

    CALL timestop(handle)

  END SUBROUTINE release_submatrix_array

! *****************************************************************************
!> \brief ...
!> \param subm ...
! *****************************************************************************
  SUBROUTINE release_submatrix(subm)

    TYPE(domain_submatrix_type), &
      INTENT(INOUT)                          :: subm

    CHARACTER(len=*), PARAMETER :: routineN = 'release_submatrix', &
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

  END SUBROUTINE release_submatrix

  ! more complex routine might be necessary if submatrices are distributed
! *****************************************************************************
!> \brief ...
!> \param transA ...
!> \param transB ...
!> \param alpha ...
!> \param A ...
!> \param B ...
!> \param beta ...
!> \param C ...
! *****************************************************************************
  SUBROUTINE multiply_submatrices_once(transA,transB,alpha,A,B,beta,C)

    CHARACTER, INTENT(IN)                    :: transA, transB
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(domain_submatrix_type), INTENT(IN)  :: A, B
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(domain_submatrix_type), &
      INTENT(INOUT)                          :: C

    CHARACTER(len=*), PARAMETER :: routineN = 'multiply_submatrices_once', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cs1, cs2, handle, icol, irow, &
                                                K, K1, LDA, LDB, LDC, M, &
                                                Mblocks, N, Nblocks
    LOGICAL                                  :: NOTA, NOTB

    CALL timeset(routineN,handle)
    
    CPASSERT(A%domain.gt.0)
    CPASSERT(B%domain.gt.0)
    CPASSERT(C%domain.gt.0)

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

    ! these checks are for debugging only
    CPASSERT(K.eq.K1)
   
    ! conform C matrix 
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
       !!! cannot use non-zero beta if C is not allocated
       CPASSERT(beta.eq.0.0_dp)
       ALLOCATE(C%mdata(C%nrows,C%ncols))
    ELSE
       cs1=SIZE(C%mdata,1)
       cs2=SIZE(C%mdata,2)
       IF ((C%nrows.ne.cs1).OR.(C%ncols.ne.cs2)) THEN
          !!! cannot deallocate data if beta is non-zero
          CPASSERT(beta.eq.0.0_dp)
          DEALLOCATE(C%mdata)
          ALLOCATE(C%mdata(C%nrows,C%ncols))
       ENDIF
    ENDIF

    LDC = C%nrows

    CALL DGEMM(transA,transB,M,N,K,alpha,A%mdata,LDA,B%mdata,LDB,beta,C%mdata,LDC)

    C%nnodes=A%nnodes
    C%groupid=A%groupid

    CALL timestop(handle)

  END SUBROUTINE multiply_submatrices_once

! *****************************************************************************
!> \brief ...
!> \param transA ...
!> \param transB ...
!> \param alpha ...
!> \param A ...
!> \param B ...
!> \param beta ...
!> \param C ...
! *****************************************************************************
  SUBROUTINE multiply_submatrices_array(transA,transB,alpha,A,B,beta,C)

    CHARACTER, INTENT(IN)                    :: transA, transB
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(IN)               :: A, B
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(INOUT)            :: C

    CHARACTER(len=*), PARAMETER :: routineN = 'multiply_submatrices_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, idomainA, &
                                                idomainB, ndomains, &
                                                ndomainsB, ndomainsC

    CALL timeset(routineN,handle)

    ndomains=SIZE(A)
    ndomainsB=SIZE(B)
    ndomainsC=SIZE(C)

    CPASSERT(ndomains.eq.ndomainsB)
    CPASSERT(ndomainsB.eq.ndomainsC)
    
    DO idomain = 1, ndomains

       idomainA = A(idomain)%domain
       idomainB = B(idomain)%domain
       
       CPASSERT(idomainA.eq.idomainB)
       
       C(idomain)%domain = idomainA

       ! check if the submatrix exists
       IF (idomainA.gt.0) THEN
          CALL multiply_submatrices_once(transA,transB,alpha,A(idomain),B(idomain),beta,C(idomain))
       ENDIF ! submatrix for the domain exists

    ENDDO ! loop over domains

    CALL timestop(handle)

  END SUBROUTINE multiply_submatrices_array

  ! more complex routine might be necessary if submatrices are distributed
! *****************************************************************************
!> \brief ...
!> \param alpha ...
!> \param A ...
!> \param beta ...
!> \param B ...
!> \param transB ...
! *****************************************************************************
  SUBROUTINE add_submatrices_once(alpha,A,beta,B,transB)

    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(domain_submatrix_type), &
      INTENT(INOUT)                          :: A
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(domain_submatrix_type), INTENT(IN)  :: B
    CHARACTER, INTENT(IN)                    :: transB

    CHARACTER(len=*), PARAMETER :: routineN = 'add_submatrices_once', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: C1, C2, handle, icol, R1, R2
    LOGICAL                                  :: NOTB

    CALL timeset(routineN,handle)

    CPASSERT(A%domain.gt.0)
    CPASSERT(B%domain.gt.0)
    
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

    ! these checks are for debugging only
    CPASSERT(C1.eq.C2)
    CPASSERT(R1.eq.R2)
   
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

  END SUBROUTINE add_submatrices_once

! *****************************************************************************
!> \brief ...
!> \param alpha ...
!> \param A ...
!> \param beta ...
!> \param B ...
!> \param transB ...
! *****************************************************************************
  SUBROUTINE add_submatrices_array(alpha,A,beta,B,transB)

    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(INOUT)            :: A
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(IN)               :: B
    CHARACTER, INTENT(IN)                    :: transB

    CHARACTER(len=*), PARAMETER :: routineN = 'add_submatrices_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, idomainA, &
                                                idomainB, ndomains, ndomainsB

    CALL timeset(routineN,handle)

    ndomains=SIZE(A)
    ndomainsB=SIZE(B)

    CPASSERT(ndomains.eq.ndomainsB)
    
    DO idomain = 1, ndomains

       idomainA = A(idomain)%domain
       idomainB = B(idomain)%domain
       
       CPASSERT(idomainA.eq.idomainB)
       
       ! check if the submatrix exists
       IF (idomainA.gt.0) THEN
          CALL add_submatrices_once(alpha,A(idomain),beta,B(idomain),transB)
       ENDIF ! submatrix for the domain exists

    ENDDO ! loop over domains

    CALL timestop(handle)

  END SUBROUTINE add_submatrices_array

! *****************************************************************************
!> \brief Computes the max norm of the collection of submatrices
!> \param submatrices ...
!> \param norm ...
!> \par History
!>       2013.03 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE maxnorm_submatrices(submatrices,norm)
    
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(IN)               :: submatrices
    REAL(KIND=dp), INTENT(OUT)               :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'maxnorm_submatrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, ndomains
    REAL(KIND=dp)                            :: curr_norm, send_norm
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_norm

!TYPE(cp_dbcsr_type), INTENT(IN)                        :: distr_pattern

    CALL timeset(routineN,handle)

    send_norm=0.0_dp
    
    ndomains=SIZE(submatrices)

    DO idomain = 1, ndomains

       ! check if the submatrix is local
       IF (submatrices(idomain)%domain.gt.0) THEN
          curr_norm=MAXVAL(ABS(submatrices(idomain)%mdata))
          IF (curr_norm.gt.send_norm) send_norm=curr_norm
       ENDIF

    ENDDO ! loop over domains
    
    ! communicate local norm to the other nodes
    ALLOCATE(recv_norm(submatrices(1)%nnodes))
    CALL mp_allgather(send_norm,recv_norm,submatrices(1)%groupid)
   
    norm=MAXVAL(recv_norm)

    DEALLOCATE(recv_norm)

    CALL timestop(handle)
  
  END SUBROUTINE maxnorm_submatrices 

! *****************************************************************************
!> \brief Computes the sum of traces of the submatrix A.tr(B)
!> \param A ...
!> \param B ...
!> \param trace ...
!> \par History
!>       2013.03 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE trace_submatrices(A,B,trace)
    
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(IN)               :: A, B
    REAL(KIND=dp), INTENT(OUT)               :: trace

    CHARACTER(len=*), PARAMETER :: routineN = 'trace_submatrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: domainA, domainB, handle, &
                                                idomain, ndomainsA, ndomainsB
    REAL(KIND=dp)                            :: curr_trace, send_trace
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_trace

!TYPE(cp_dbcsr_type), INTENT(IN)                        :: distr_pattern

    CALL timeset(routineN,handle)

    send_trace=0.0_dp
    
    ndomainsA=SIZE(A)
    ndomainsB=SIZE(B)

    CPASSERT(ndomainsA.eq.ndomainsB)
    
    DO idomain = 1, ndomainsA

       domainA = A(idomain)%domain
       domainB = B(idomain)%domain
       
       CPASSERT(domainA.eq.domainB)

       ! check if the submatrix is local
       IF (domainA.gt.0) THEN

          CPASSERT(A(idomain)%nrows.eq.B(idomain)%nrows)
          CPASSERT(A(idomain)%ncols.eq.B(idomain)%ncols)
          
          curr_trace=SUM(A(idomain)%mdata(:,:)*B(idomain)%mdata(:,:))
          send_trace=send_trace+curr_trace

       ENDIF

    ENDDO ! loop over domains
    
    ! communicate local norm to the other nodes
    ALLOCATE(recv_trace(A(1)%nnodes))
    CALL mp_allgather(send_trace,recv_trace,A(1)%groupID)
   
    trace=SUM(recv_trace)

    DEALLOCATE(recv_trace)

    CALL timestop(handle)
  
  END SUBROUTINE trace_submatrices 

! *****************************************************************************
!> \brief Constructs submatrices for each ALMO domain by collecting distributed
!>        DBCSR blocks to local arrays
!> \param matrix ...
!> \param submatrix ...
!> \param distr_pattern ...
!> \param domain_map ...
!> \param node_of_domain ...
!> \param job_type ...
!> \par History
!>       2013.01 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE construct_submatrices(matrix,submatrix,distr_pattern,domain_map,&
     node_of_domain,job_type)
    
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(INOUT)            :: submatrix
    TYPE(cp_dbcsr_type), INTENT(IN)          :: distr_pattern
    TYPE(domain_map_type), INTENT(IN)        :: domain_map
    INTEGER, DIMENSION(:), INTENT(IN)        :: node_of_domain
    INTEGER, INTENT(IN)                      :: job_type

    CHARACTER(len=*), PARAMETER :: routineN = 'construct_submatrices', &
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
    LOGICAL                                  :: found, transp
    REAL(KIND=dp)                            :: antifactor
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_data, send_data
    REAL(KIND=dp), DIMENSION(:), POINTER     :: block_p

!INTEGER, PARAMETER                         :: select_row_col = 1
!INTEGER, PARAMETER                         :: select_row = 2
!                                                  subm_row_size,&
!                                                  subm_col_size,&

    CALL timeset(routineN,handle)
    
    nblkrows_tot = cp_dbcsr_nblkrows_total(matrix)
    nblkcols_tot = cp_dbcsr_nblkcols_total(matrix)

    ndomains = nblkcols_tot ! RZK-warning not true for atomic distributions

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
    
    ! find: the number of blocks and their sizes that must be sent to each cpu
    ! loop over all domains
    DO idomain = 1, ndomains
             
       dest_node=node_of_domain(idomain)

       ! loop over those rows that have non-zero quencher
       index_sr=1 ! index start row
       IF (idomain.gt.1) index_sr=domain_map%index1(idomain-1)
       index_er=domain_map%index1(idomain)-1 ! index end row

       DO index_row=index_sr, index_er

          row = domain_map%pairs(index_row,1)

          IF (job_type==select_row_col) THEN
             ! loop over those columns that have non-zero quencher
             index_sc=1 ! index start col
             IF (idomain.gt.1) index_sc=domain_map%index1(idomain-1)
             index_ec=domain_map%index1(idomain)-1 ! index end col
          ELSE
             ! fake loop
             index_sc=1 ! index start col
             index_ec=1 ! index end col
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

          ENDDO ! loop over columns

       ENDDO ! loop over rows

    ENDDO

    ! simple but quadratically scaling procedure
    ! loop over local blocks
    !CALL cp_dbcsr_iterator_start(iter,matrix)
    !DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
    !   CALL cp_dbcsr_iterator_next_block(iter,row,col,data_p,&
    !           row_size=row_size,col_size=col_size)
    !   DO idomain = 1, ndomains
    !      IF (job_type==select_row_col) THEN
    !         domain_needs_block=(qblk_exists(domain_map,col,idomain)&
    !            .AND.qblk_exists(domain_map,row,idomain))
    !      ELSE
    !         domain_needs_block=(idomain==col&
    !            .AND.qblk_exists(domain_map,row,idomain))
    !      ENDIF
    !      IF (domain_needs_block) THEN
    !         transp=.FALSE.
    !         dest_node=node_of_domain(idomain)
    !         !CALL cp_dbcsr_get_stored_coordinates(distr_pattern,&
    !         !        idomain, idomain, transp, dest_node)
    !         send_descriptor(1,dest_node+1)=send_descriptor(1,dest_node+1)+1
    !         send_descriptor(2,dest_node+1)=send_descriptor(2,dest_node+1)+&
    !            row_size*col_size
    !      ENDIF
    !   ENDDO
    !ENDDO
    !CALL cp_dbcsr_iterator_stop(iter)
       
    ! communicate number of blocks and their sizes to the other nodes
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
    
    ! collect data from the matrix into send_data
    ALLOCATE(send_data(send_size_total))
    ALLOCATE(recv_data(recv_size_total))
    ALLOCATE(send_data2(send_size2_total))
    ALLOCATE(recv_data2(recv_size2_total))
    ALLOCATE(offset_block(nNodes))
    ALLOCATE(offset2_block(nNodes))
    offset_block(:)=0
    offset2_block(:)=0
    ! loop over all domains
    DO idomain = 1, ndomains
             
       dest_node=node_of_domain(idomain)

       ! loop over those rows that have non-zero quencher
       index_sr=1 ! index start row
       IF (idomain.gt.1) index_sr=domain_map%index1(idomain-1)
       index_er=domain_map%index1(idomain)-1 ! index end row

       DO index_row=index_sr, index_er

          row = domain_map%pairs(index_row,1)

          IF (job_type==select_row_col) THEN
             ! loop over those columns that have non-zero quencher
             index_sc=1 ! index start col
             IF (idomain.gt.1) index_sc=domain_map%index1(idomain-1)
             index_ec=domain_map%index1(idomain)-1 ! index end col
          ELSE
             ! fake loop
             index_sc=1 ! index start col
             index_ec=1 ! index end col
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
                   !col_offset=0
                   !DO icol=1,col_size
                   !   start_data=send_offset_cpu(dest_node+1)+&
                   !              offset_block(dest_node+1)+&
                   !              col_offset
                   !   send_data(start_data+1:start_data+row_size)=&
                   !      data_p(1:row_size,icol)
                   !   col_offset=col_offset+row_size
                   !ENDDO
                   col_offset=row_size*col_size
                   start_data=send_offset_cpu(dest_node+1)+&
                                 offset_block(dest_node+1)
                   send_data(start_data+1:start_data+col_offset)=&
                         block_p(1:col_offset)
                   offset_block(dest_node+1)=offset_block(dest_node+1)+col_offset
                   ! fill out row,col information
                   send_data2(send_offset2_cpu(dest_node+1)+&
                         offset2_block(dest_node+1)+1)=row
                   send_data2(send_offset2_cpu(dest_node+1)+&
                         offset2_block(dest_node+1)+2)=col
                   offset2_block(dest_node+1)=offset2_block(dest_node+1)+2
                ENDIF
             ENDIF

          ENDDO ! loop over columns

       ENDDO ! loop over rows

    ENDDO
    ! more simple but quadratically scaling version
    !CALL cp_dbcsr_iterator_start(iter,matrix)
    !DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
    !   CALL cp_dbcsr_iterator_next_block(iter,row,col,data_p,&
    !           row_size=row_size,col_size=col_size)
    !   DO idomain = 1, ndomains
    !      IF (job_type==select_row_col) THEN
    !         domain_needs_block=(qblk_exists(domain_map,col,idomain)&
    !            .AND.qblk_exists(domain_map,row,idomain))
    !      ELSE
    !         domain_needs_block=(idomain==col&
    !            .AND.qblk_exists(domain_map,row,idomain))
    !      ENDIF
    !      IF (domain_needs_block) THEN
    !         transp=.FALSE.
    !         dest_node=node_of_domain(idomain)
    !         !CALL cp_dbcsr_get_stored_coordinates(distr_pattern,&
    !         !        idomain, idomain, transp, dest_node)
    !         ! place the data appropriately
    !         col_offset=0
    !         DO icol=1,col_size
    !            start_data=send_offset_cpu(dest_node+1)+&
    !                       offset_block(dest_node+1)+&
    !                       col_offset
    !            send_data(start_data+1:start_data+row_size)=&
    !               data_p(1:row_size,icol)
    !            col_offset=col_offset+row_size
    !         ENDDO
    !         offset_block(dest_node+1)=offset_block(dest_node+1)+col_size*row_size
    !         ! fill out row,col information
    !         send_data2(send_offset2_cpu(dest_node+1)+&
    !               offset2_block(dest_node+1)+1)=row
    !         send_data2(send_offset2_cpu(dest_node+1)+&
    !               offset2_block(dest_node+1)+2)=col
    !         offset2_block(dest_node+1)=offset2_block(dest_node+1)+2
    !      ENDIF
    !   ENDDO
    !ENDDO
    !CALL cp_dbcsr_iterator_stop(iter)
       
    ! send-receive all blocks
    CALL mp_alltoall(send_data,send_size_cpu,send_offset_cpu,&
            recv_data,recv_size_cpu,recv_offset_cpu,GroupID)
    ! send-receive rows and cols of the blocks
    CALL mp_alltoall(send_data2,send_size2_cpu,send_offset2_cpu,&
            recv_data2,recv_size2_cpu,recv_offset2_cpu,GroupID)
    
    DEALLOCATE(send_size_cpu,send_offset_cpu)
    DEALLOCATE(send_size2_cpu,send_offset2_cpu)
    DEALLOCATE(send_data)
    DEALLOCATE(send_data2)
    DEALLOCATE(offset_block)
    DEALLOCATE(offset2_block)
    
    ! copy blocks into submatrices
    col_blk_size => cp_dbcsr_col_block_sizes(matrix)
    row_blk_size => cp_dbcsr_row_block_sizes(matrix)
!    ALLOCATE(subm_row_size(ndomains),subm_col_size(ndomains))
!    subm_row_size(:)=0
!    subm_col_size(:)=0
    ndomains2=SIZE(submatrix)
    IF (ndomains2.ne.ndomains) THEN
       CPABORT("wrong submatrix size")
       CPASSERT(.FALSE.)
    ENDIF
    CALL release_submatrices(submatrix)
    submatrix(:)%nnodes=nNodes
    submatrix(:)%groupID=GroupID
    submatrix(:)%nrows=0
    submatrix(:)%ncols=0
    
    ALLOCATE(first_row(nblkrows_tot),first_col(nblkcols_tot))
    submatrix(:)%domain=-1
    DO idomain=1, ndomains 
       dest_node=node_of_domain(idomain)
       !transp=.FALSE.
       !CALL cp_dbcsr_get_stored_coordinates(distr_pattern,&
       !        idomain, idomain, transp, dest_node)
       IF (dest_node.eq.mynode) THEN
          submatrix(idomain)%domain=idomain
          submatrix(idomain)%nbrows=0
          submatrix(idomain)%nbcols=0
          
          ! loop over those rows that have non-zero quencher
          first_row(:)=-1
          index_sr=1 ! index start row
          IF (idomain.gt.1) index_sr=domain_map%index1(idomain-1)
          index_er=domain_map%index1(idomain)-1 ! index end row
          DO index_row=index_sr, index_er
             row = domain_map%pairs(index_row,1)
          !DO row = 1, nblkrows_tot
          !   IF (qblk_exists(domain_map,row,idomain)) THEN
             first_row(row)=submatrix(idomain)%nrows+1
             submatrix(idomain)%nrows=submatrix(idomain)%nrows+row_blk_size(row)
             submatrix(idomain)%nbrows=submatrix(idomain)%nbrows+1
          !   ENDIF
          ENDDO
          ALLOCATE(submatrix(idomain)%dbcsr_row(submatrix(idomain)%nbrows))
          ALLOCATE(submatrix(idomain)%size_brow(submatrix(idomain)%nbrows))
          smrow=1
          ! again loop over those rows that have non-zero quencher
          index_sr=1 ! index start row
          IF (idomain.gt.1) index_sr=domain_map%index1(idomain-1)
          index_er=domain_map%index1(idomain)-1 ! index end row
          DO index_row=index_sr, index_er
             row = domain_map%pairs(index_row,1)
          !DO row = 1, nblkrows_tot
          !   IF (first_row(row).ne.-1) THEN
                submatrix(idomain)%dbcsr_row(smrow)=row
                submatrix(idomain)%size_brow(smrow)=row_blk_size(row)
                smrow=smrow+1
          !   ENDIF
          ENDDO

          ! loop over the necessary columns
          first_col(:)=-1
          IF (job_type==select_row_col) THEN
             ! loop over those columns that have non-zero quencher
             index_sc=1 ! index start col
             IF (idomain.gt.1) index_sc=domain_map%index1(idomain-1)
             index_ec=domain_map%index1(idomain)-1 ! index end col
          ELSE
             ! fake loop
             index_sc=1 ! index start col
             index_ec=1 ! index end col
          ENDIF
          DO index_col=index_sc, index_ec
             IF (job_type==select_row_col) THEN
                col = domain_map%pairs(index_col,1)
             ELSE
                col = idomain 
             ENDIF
          !DO col = 1, nblkcols_tot
          !   IF (job_type==select_row_col) THEN
          !      domain_needs_block=(qblk_exists(domain_map,col,idomain))
          !   ELSE
          !      domain_needs_block=(col==idomain) ! RZK-warning col belongs to the domain
          !   ENDIF
          !   IF (domain_needs_block) THEN
                first_col(col)=submatrix(idomain)%ncols+1
                submatrix(idomain)%ncols=submatrix(idomain)%ncols+col_blk_size(col)
                submatrix(idomain)%nbcols=submatrix(idomain)%nbcols+1
          !   ENDIF
          ENDDO
          
          ALLOCATE(submatrix(idomain)%dbcsr_col(submatrix(idomain)%nbcols))
          ALLOCATE(submatrix(idomain)%size_bcol(submatrix(idomain)%nbcols))
          
          ! loop over the necessary columns again
          smcol=1
          IF (job_type==select_row_col) THEN
             ! loop over those columns that have non-zero quencher
             index_sc=1 ! index start col
             IF (idomain.gt.1) index_sc=domain_map%index1(idomain-1)
             index_ec=domain_map%index1(idomain)-1 ! index end col
          ELSE
             ! fake loop
             index_sc=1 ! index start col
             index_ec=1 ! index end col
          ENDIF
          DO index_col=index_sc, index_ec
             IF (job_type==select_row_col) THEN
                col = domain_map%pairs(index_col,1)
             ELSE
                col = idomain 
             ENDIF
          !DO col = 1, nblkcols_tot
          !   IF (first_col(col).ne.-1) THEN
                submatrix(idomain)%dbcsr_col(smcol)=col
                submatrix(idomain)%size_bcol(smcol)=col_blk_size(col)
                smcol=smcol+1
          !   ENDIF
          ENDDO

          ALLOCATE(submatrix(idomain)%mdata(&
             submatrix(idomain)%nrows,&
             submatrix(idomain)%ncols))
          submatrix(idomain)%mdata(:,:)=0.0_dp
          DO iNode=1, nNodes
             block_offset=0
             DO iBlock=1, recv_size2_cpu(iNode)/2
                ! read the (row,col) of the block
                row = recv_data2(recv_offset2_cpu(iNode)+(iBlock-1)*2+1)
                col = recv_data2(recv_offset2_cpu(iNode)+(iBlock-1)*2+2)
                ! check if this block should be in the submatrix of this domain
                IF ((first_col(col).NE.-1).AND.(first_row(row).NE.-1)) THEN
                   ! copy data from the received array into submatrix
                   start_data=recv_offset_cpu(iNode)+block_offset+1
                   DO icol=0, col_blk_size(col)-1
                      submatrix(idomain)%mdata(first_row(row):&
                                first_row(row)+row_blk_size(row)-1,&
                                first_col(col)+icol)=&
                                recv_data(start_data:start_data+row_blk_size(row)-1)
                      start_data=start_data+row_blk_size(row)
                   ENDDO
                   IF (job_type==select_row_col) THEN
                      IF (matrix_type==dbcsr_type_symmetric.OR.&
                          matrix_type==dbcsr_type_antisymmetric) THEN
                         ! copy data into the transposed block as well
                         antifactor=1.0_dp
                         IF (matrix_type==dbcsr_type_antisymmetric) THEN
                            antifactor=-1.0_dp
                         ENDIF
                         start_data=recv_offset_cpu(iNode)+block_offset+1
                         DO icol=0, col_blk_size(col)-1
                            submatrix(idomain)%mdata(first_row(col)+icol,&
                                      first_col(row):&
                                      first_col(row)+row_blk_size(row)-1)=&
                                      antifactor*recv_data(start_data:&
                                      start_data+row_blk_size(row)-1)
                            start_data=start_data+row_blk_size(row)
                         ENDDO
                      ELSE IF (matrix_type==dbcsr_type_no_symmetry) THEN
                      ELSE
                         CPABORT("matrix type is NYI")
                         CPASSERT(.FALSE.)
                      ENDIF
                   ENDIF
                ENDIF
                block_offset=block_offset+col_blk_size(col)*row_blk_size(row)
             ENDDO
          ENDDO
       ENDIF ! mynode.eq.dest_node
    ENDDO ! loop over domains

    DEALLOCATE(recv_size_cpu,recv_offset_cpu)
    DEALLOCATE(recv_size2_cpu,recv_offset2_cpu)
    DEALLOCATE(recv_data)
    DEALLOCATE(recv_data2)
!    DEALLOCATE(subm_row_size,subm_col_size)
    DEALLOCATE(first_row,first_col)
       
    CALL timestop(handle)
  
  END SUBROUTINE construct_submatrices 

! *****************************************************************************
!> \brief Constructs a DBCSR matrix from submatrices
!> \param matrix ...
!> \param submatrix ...
!> \param distr_pattern ...
!> \par History
!>       2013.01 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE construct_dbcsr_from_submatrices(matrix,submatrix,distr_pattern)
    
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(IN)               :: submatrix
    TYPE(cp_dbcsr_type), INTENT(IN)          :: distr_pattern

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'construct_dbcsr_from_submatrices', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: matrix_type
    INTEGER :: block_offset, col, col_offset, colsize, dest_node, GroupID, &
      handle, iBlock, icol, idomain, iNode, irow_subm, ldesc, myNode, &
      nblkcols_tot, nblkrows_tot, ndomains, ndomains2, nNodes, &
      recv_size2_total, recv_size_total, row, rowsize, send_size2_total, &
      send_size_total, smroff, start_data, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:) :: offset2_block, offset_block, &
      recv_data2, recv_offset2_cpu, recv_offset_cpu, recv_size2_cpu, &
      recv_size_cpu, send_data2, send_offset2_cpu, send_offset_cpu, &
      send_size2_cpu, send_size_cpu
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_descriptor, &
                                                send_descriptor
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: transp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_data, send_data
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_p
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    
    ! get a useful output_unit
    logger => cp_get_default_logger()
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    nblkrows_tot = cp_dbcsr_nblkrows_total(matrix)
    nblkcols_tot = cp_dbcsr_nblkcols_total(matrix)

    ndomains = nblkcols_tot ! RZK-warning not true for atomic distributions
    ndomains2 = SIZE(submatrix)

    IF (ndomains.ne.ndomains2) THEN
       CPABORT("domain mismatch")
       CPASSERT(.FALSE.)
    ENDIF

    nNodes = dbcsr_mp_numnodes(dbcsr_distribution_mp(&
       cp_dbcsr_distribution(distr_pattern)))
    GroupID = dbcsr_mp_group(dbcsr_distribution_mp(&
       cp_dbcsr_distribution(distr_pattern)))
    myNode = dbcsr_mp_mynode(dbcsr_distribution_mp(&
       cp_dbcsr_distribution(distr_pattern)))

    matrix_type=cp_dbcsr_get_matrix_type(matrix)
    IF (matrix_type.ne.dbcsr_type_no_symmetry) THEN
       CPABORT("only non-symmetric matrices so far")
       CPASSERT(.FALSE.)
    ENDIF
   
    ! remove all blocks from the dbcsr matrix   
    CALL cp_dbcsr_iterator_start(iter,matrix)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter,row,col,block_p)
       block_p(:,:)=0.0_dp
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)
    CALL cp_dbcsr_filter(matrix,0.1_dp)

    CALL cp_dbcsr_work_create(matrix,work_mutable=.TRUE.)

    ldesc=2
    ALLOCATE(send_descriptor(ldesc,nNodes))
    ALLOCATE(recv_descriptor(ldesc,nNodes))
    send_descriptor(:,:)=0
    
    ! loop over domains - find how much data to send 
    DO idomain = 1, ndomains

       IF (submatrix(idomain)%domain.gt.0) THEN

          DO irow_subm = 1, submatrix(idomain)%nbrows
          
             IF (submatrix(idomain)%nbcols.ne.1) THEN
                CPABORT("corrupt submatrix structure")
                CPASSERT(.FALSE.)
             ENDIF

             row = submatrix(idomain)%dbcsr_row(irow_subm)
             col = submatrix(idomain)%dbcsr_col(1)

             IF (col.ne.idomain) THEN
                CPABORT("corrupt submatrix structure")
                CPASSERT(.FALSE.)
             ENDIF

             transp=.FALSE.
             CALL cp_dbcsr_get_stored_coordinates(distr_pattern,&
                     row, idomain, transp, dest_node)
             
             send_descriptor(1,dest_node+1)=send_descriptor(1,dest_node+1)+1
             send_descriptor(2,dest_node+1)=send_descriptor(2,dest_node+1)+&
                submatrix(idomain)%size_brow(irow_subm)*&
                submatrix(idomain)%size_bcol(1)

          ENDDO ! loop over submatrix blocks

       ENDIF       

    ENDDO ! loop over domains

    ! communicate number of blocks and their sizes to the other nodes
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
    
    ! collect data from the matrix into send_data
    ALLOCATE(send_data(send_size_total))
    ALLOCATE(recv_data(recv_size_total))
    ALLOCATE(send_data2(send_size2_total))
    ALLOCATE(recv_data2(recv_size2_total))
    ALLOCATE(offset_block(nNodes))
    ALLOCATE(offset2_block(nNodes))
    offset_block(:)=0
    offset2_block(:)=0
    ! loop over domains - collect data to send 
    DO idomain = 1, ndomains

       IF (submatrix(idomain)%domain.gt.0) THEN

          smroff=0
          DO irow_subm = 1, submatrix(idomain)%nbrows
          
             row = submatrix(idomain)%dbcsr_row(irow_subm)
             col = submatrix(idomain)%dbcsr_col(1)

             rowsize = submatrix(idomain)%size_brow(irow_subm)
             colsize = submatrix(idomain)%size_bcol(1)

             transp=.FALSE.
             CALL cp_dbcsr_get_stored_coordinates(distr_pattern,&
                     row, idomain, transp, dest_node)
            
             ! place the data appropriately
             col_offset=0
             DO icol=1,colsize
                start_data=send_offset_cpu(dest_node+1)+&
                           offset_block(dest_node+1)+&
                           col_offset
                send_data(start_data+1:start_data+rowsize)=&
                   submatrix(idomain)%mdata(smroff+1:smroff+rowsize,icol)
                col_offset=col_offset+rowsize
             ENDDO
             offset_block(dest_node+1)=offset_block(dest_node+1)+&
                colsize*rowsize
             ! fill out row,col information
             send_data2(send_offset2_cpu(dest_node+1)+&
                   offset2_block(dest_node+1)+1)=row
             send_data2(send_offset2_cpu(dest_node+1)+&
                   offset2_block(dest_node+1)+2)=col
             offset2_block(dest_node+1)=offset2_block(dest_node+1)+2

             smroff=smroff+rowsize

          ENDDO

       ENDIF       

    ENDDO ! loop over domains
    
    ! send-receive all blocks
    CALL mp_alltoall(send_data,send_size_cpu,send_offset_cpu,&
            recv_data,recv_size_cpu,recv_offset_cpu,GroupID)
    ! send-receive rows and cols of the blocks
    CALL mp_alltoall(send_data2,send_size2_cpu,send_offset2_cpu,&
            recv_data2,recv_size2_cpu,recv_offset2_cpu,GroupID)
    
    DEALLOCATE(send_size_cpu,send_offset_cpu)
    DEALLOCATE(send_size2_cpu,send_offset2_cpu)
    DEALLOCATE(send_data)
    DEALLOCATE(send_data2)
    DEALLOCATE(offset_block)
    DEALLOCATE(offset2_block)
    
    ! copy recieved data into dbcsr matrix
    col_blk_size => cp_dbcsr_col_block_sizes(matrix)
    row_blk_size => cp_dbcsr_row_block_sizes(matrix)
    DO iNode=1, nNodes
       block_offset=0
       DO iBlock=1, recv_size2_cpu(iNode)/2
          ! read the (row,col) of the block
          row = recv_data2(recv_offset2_cpu(iNode)+(iBlock-1)*2+1)
          col = recv_data2(recv_offset2_cpu(iNode)+(iBlock-1)*2+2)
          !CALL cp_dbcsr_get_block_p(matrix,row,col,block_p,found)
          !IF (.NOT.found) THEN
          NULLIFY (block_p)
          CALL cp_dbcsr_reserve_block2d(matrix,row,col,block_p)
          CPASSERT(ASSOCIATED(block_p))
          !ENDIF
          ! copy data from the received array into the matrix block
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
       
    CALL cp_dbcsr_finalize(matrix)
    
    CALL timestop(handle)
  
  END SUBROUTINE construct_dbcsr_from_submatrices 

! *****************************************************************************
!> \brief ...
!> \param submatrices ...
!> \param mpgroup ...
! *****************************************************************************
  SUBROUTINE print_submatrices(submatrices,mpgroup)
    
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(IN)               :: submatrices
    INTEGER, INTENT(IN)                      :: mpgroup

    CHARACTER(len=*), PARAMETER :: routineN = 'print_submatrices', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=30)                        :: colstr, formatstr
    INTEGER                                  :: handle, i, irow, n, ncols, &
                                                nrows

    CALL timeset(routineN,handle)
   
    n=SIZE(submatrices)

    DO i=1, n
       nrows=SIZE(submatrices(i)%mdata,1)
       ncols=SIZE(submatrices(i)%mdata,2)
       WRITE (colstr,*) ncols
       formatstr='('//TRIM(ADJUSTL(colstr))//'F16.9)'
       IF (submatrices(i)%domain.gt.0) THEN
          WRITE (*,*) "SUBMATRIX: ", i, nrows, 'x', ncols
          nrows=SIZE(submatrices(i)%mdata,1)
          DO irow=1, nrows
             WRITE (*,formatstr) submatrices(i)%mdata(irow,:)
          ENDDO
       ENDIF
       CALL mp_sync(mpgroup)
    ENDDO

    CALL timestop(handle)
  
  END SUBROUTINE print_submatrices 

! *****************************************************************************
!> \brief Reports whether the DBCSR block (row,col) exists in the quencher
!> \param map ...
!> \param row ...
!> \param col ...
!> \retval qblk_exists ...
!> \par History
!>       2013.01 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  FUNCTION qblk_exists(map,row,col)
  
    TYPE(domain_map_type), INTENT(IN)        :: map
    INTEGER, INTENT(IN)                      :: row, col
    LOGICAL                                  :: qblk_exists

    CHARACTER(len=*), PARAMETER :: routineN = 'qblk_exists', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first, last, mid, ndomains

!CALL timeset(routineN,handle)

    ndomains=SIZE(map%index1)

    qblk_exists=.FALSE.
    IF (col.lt.1.OR.col.gt.ndomains) RETURN
    first=1
    IF (col.gt.1) first=map%index1(col-1)
    last=map%index1(col)-1

    ! perform binary search within first-last
    DO WHILE (last.ge.first)
        mid = first+(last-first)/2
        IF (map%pairs(mid,1).gt.row) THEN
           last = mid-1
        ELSE IF (map%pairs(mid,1).lt.row) THEN
           first = mid+1
        ELSE
            qblk_exists=.TRUE.   ! SUCCESS!!
            EXIT
        END IF
    END DO

    !CALL timestop(handle)
  
    RETURN

  END FUNCTION qblk_exists

END MODULE domain_submatrix_methods


