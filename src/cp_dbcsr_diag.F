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
  USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_cfm_type
  USE cp_dbcsr_operations,             ONLY: copy_cfm_to_dbcsr,&
                                             copy_dbcsr_to_cfm,&
                                             copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr
  USE cp_fm_diag,                      ONLY: cp_fm_syevd,&
                                             cp_fm_syevr,&
                                             cp_fm_syevx
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE dbcsr_methods,                   ONLY: dbcsr_get_info
  USE dbcsr_types,                     ONLY: dbcsr_obj
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_diag'

  ! Public subroutines

  PUBLIC :: dbcsr_syevd, &
            dbcsr_syevr, &
            dbcsr_syevx, &
            dbcsr_heevd

CONTAINS

! *****************************************************************************
  SUBROUTINE dbcsr_syevd(matrix,eigenvectors,eigenvalues,para_env,blacs_env,error)

    ! Computes all eigenvalues and vectors of a real symmetric matrix
    ! should be quite a bit faster than syevx for that case
    ! especially in parallel with thightly clustered evals
    ! needs more workspace in the worst case, but much better distributed

    TYPE(dbcsr_obj)                          :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_syevd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nfullrows_total
    LOGICAL                                  :: failure
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_eigenvectors, fm_matrix

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(fm_matrix, fm_eigenvectors, fm_struct)
    CALL dbcsr_get_info(matrix,nfullrows_total=nfullrows_total)

    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullrows_total,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrix,fm_struct,name="fm_matrix",error=error)
    CALL cp_fm_create(fm_eigenvectors,fm_struct,name="fm_eigenvectors",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL copy_dbcsr_to_fm(matrix, fm_matrix, error=error)

    CALL cp_fm_syevd(fm_matrix,fm_eigenvectors,eigenvalues,error)

    CALL copy_fm_to_dbcsr(fm_eigenvectors, eigenvectors, error=error)

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
  SUBROUTINE dbcsr_syevx(matrix,eigenvectors,eigenvalues,neig,work_syevx,&
       para_env,blacs_env,error)

    ! Diagonalise the symmetric n by n matrix using the LAPACK library.

    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(dbcsr_obj), OPTIONAL, POINTER       :: eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    INTEGER, INTENT(IN), OPTIONAL            :: neig
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: work_syevx
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_syevx', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n, neig_local
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_eigenvectors, fm_matrix

    CALL timeset(routineN,handle)

    ! by default all
    CALL dbcsr_get_info(matrix,nfullrows_total=n)
    neig_local=n
    IF (PRESENT(neig)) neig_local=neig
    IF (neig_local == 0) RETURN

    NULLIFY(fm_matrix, fm_eigenvectors, fm_struct)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=n,&
         ncol_global=n,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrix,fm_struct,name="fm_matrix",error=error)
    IF(PRESENT(eigenvectors)) &
         CALL cp_fm_create(fm_eigenvectors,fm_struct,name="fm_eigenvectors",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL copy_dbcsr_to_fm(matrix, fm_matrix, error=error)

    CALL cp_fm_syevx(fm_matrix,fm_eigenvectors,eigenvalues,neig,work_syevx,error)

    IF(PRESENT(eigenvectors)) CALL copy_fm_to_dbcsr(fm_eigenvectors, eigenvectors, &
         error=error)
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
  SUBROUTINE dbcsr_syevr(matrix,eigenvectors,eigenvalues,ilow,iup,para_env,blacs_env,error)

    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(dbcsr_obj), OPTIONAL, POINTER       :: eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    INTEGER, INTENT(IN), OPTIONAL            :: ilow, iup
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_syevr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_eigenvectors, fm_matrix

! -------------------------------------------------------------------------

    NULLIFY(fm_matrix, fm_eigenvectors, fm_struct)
    CALL dbcsr_get_info(matrix,nfullrows_total=n)

    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=n,&
         ncol_global=n,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrix,fm_struct,name="fm_matrix",error=error)
    CALL cp_fm_create(fm_eigenvectors,fm_struct,name="fm_eigenvectors",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL copy_dbcsr_to_fm(matrix, fm_matrix, error=error)

    CALL cp_fm_syevr(fm_matrix,fm_eigenvectors,eigenvalues,ilow,iup,error)

    IF(PRESENT(eigenvectors)) CALL copy_fm_to_dbcsr(fm_eigenvectors, eigenvectors, &
         error=error)
    IF(PRESENT(eigenvectors)) CALL cp_fm_release(fm_eigenvectors, error)
    CALL cp_fm_release(fm_matrix, error)

    CALL timestop(handle)

  END SUBROUTINE dbcsr_syevr

! *****************************************************************************
  SUBROUTINE dbcsr_heevd(matrix,eigenvectors,eigenvalues,para_env,blacs_env,error)

    TYPE(dbcsr_obj)                          :: matrix
    TYPE(dbcsr_obj), OPTIONAL, POINTER       :: eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_heevd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nfullrows_total
    LOGICAL                                  :: failure
    TYPE(cp_cfm_type), POINTER               :: fm_eigenvectors, fm_matrix
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct

    CALL timeset(routineN,handle)

    failure=.FALSE.

    NULLIFY(fm_matrix, fm_eigenvectors, fm_struct)
    CALL dbcsr_get_info(matrix,nfullrows_total=nfullrows_total)

    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullrows_total,para_env=para_env,error=error)
    CALL cp_cfm_create(fm_matrix,fm_struct,name="fm_matrix",error=error)
    CALL cp_cfm_create(fm_eigenvectors,fm_struct,name="fm_eigenvectors",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL copy_dbcsr_to_cfm(matrix, fm_matrix, error=error)

    CALL cp_cfm_heevd(fm_matrix,fm_eigenvectors,eigenvalues,error)

    CALL copy_cfm_to_dbcsr(fm_eigenvectors, eigenvectors, error=error)

    CALL cp_cfm_release(fm_matrix, error)
    CALL cp_cfm_release(fm_eigenvectors, error)

    CALL timestop(handle)

  END SUBROUTINE dbcsr_heevd

END MODULE cp_dbcsr_diag
