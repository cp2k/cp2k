!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/kg_gpw_fm_mol_methods [1.0] *
!!
!!   NAME
!!     kg_gpw_fm_mol_methods
!!
!!   FUNCTION
!!     Some utilities to handle the molecular blocks of the full matrix
!!     of the molecular orbitals
!!     
!!   AUTHOR
!!     MI (29.11.2004)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE kg_gpw_fm_mol_utils

  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kg_gpw_fm_mol_types,             ONLY: get_fm_mol_block,&
                                             get_mol_mo_set,&
                                             fm_mol_blocks_type,&
                                             kg_fm_mol_set_type
  USE kinds,                           ONLY: dp,dp_size,int_size
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "kg_gpw_fm_mol_utils"

! *** Public subroutines ***

  PUBLIC :: fm_mol_power, mol_mo_random

! *****************************************************************************

CONTAINS

! *****************************************************************************

!***************************************************************************

  SUBROUTINE fm_mol_power(matrix,eigenvectors,exponent,threshold,n_dependent)

    REAL(dp), DIMENSION(:,:), POINTER            :: matrix, eigenvectors
    REAL(dp), INTENT(IN)                         :: exponent,threshold
    INTEGER , INTENT(OUT)                        :: n_dependent

    CHARACTER(LEN=*), PARAMETER :: routineN ="fm_mol_power",&
         routineP = module_name//':'//routineN

    INTEGER                                      :: handle, icol, istat,&
                                                    ncol,nrow
    REAL(dp)                                     :: f,p
    REAL(dp), DIMENSION(:), ALLOCATABLE          :: eigenvalues

    CALL timeset("fm_mol_power","I","",handle)

    nrow = SIZE(matrix,1)
    ncol = SIZE(matrix,2)

    n_dependent = 0

    ALLOCATE(eigenvalues(ncol),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"eigenvalues",ncol*dp_size)
    eigenvalues(:) = 0.0_dp

    CALL fm_mol_syevx(matrix,eigenvectors,eigenvalues,ncol)

    p = 0.5_dp*exponent

    DO icol = 1,ncol
       IF(eigenvalues(icol) < threshold) THEN 
         n_dependent = n_dependent + 1
         eigenvectors(1:nrow,icol) = 0.0_dp
       ELSE

         f = eigenvalues(icol)**p
         eigenvectors(1:nrow,icol) =&
             f*eigenvectors(1:nrow,icol)
       END IF
    END DO

    CALL dsyrk("U","N",nrow,ncol,1.0_dp,eigenvectors(1,1),nrow,&
                                 0.0_dp,matrix(1,1),nrow)

    DEALLOCATE(eigenvalues, STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineP,"eigenvalues")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE fm_mol_power

!***************************************************************************

  SUBROUTINE fm_mol_syevx(matrix,eigenvetors,eigenvalues,ncol,error)

    REAL(dp), DIMENSION(:,:), POINTER            :: matrix, eigenvetors
    REAL(dp), DIMENSION(:), INTENT(OUT)          :: eigenvalues
    INTEGER, INTENT(IN)                          :: ncol
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN ="fm_mol_syevx",&
         routineP = module_name//':'//routineN

    REAL(KIND = dp), PARAMETER                   :: vl = 0.0_dp,&
                                                    vu = 0.0_dp 

!   *** Local variables ***

    INTEGER  :: info, istat, liwork, lwork, m, nb, nrow
    INTEGER,  DIMENSION(:), ALLOCATABLE          :: ifail,iwork
    REAL(dp)                                     :: abstol
    REAL(dp), DIMENSION(:), ALLOCATABLE          :: w, work
    REAL(dp), DIMENSION(:,:), POINTER            :: a,z

    REAL(KIND = dp), EXTERNAL                    :: dlamch
    INTEGER, EXTERNAL                            :: ilaenv


    abstol = 2.0_dp*dlamch("S")
    nrow = SIZE(matrix,1)

    ALLOCATE(w(nrow), STAT =istat)
    IF (istat /= 0) CALL stop_memory(routineP,"w",nrow*dp_size)
    w(:) = 0.0_dp

    eigenvalues(:) = 0.0_dp

    nb = MAX(ilaenv(1,"DSYTRD","U",nrow,-1,-1,-1),&
             ilaenv(1,"DORMTR","U",nrow,-1,-1,-1))

    lwork = MAX((nb + 3)*nrow,8*nrow)+nrow ! sun bug fix
    liwork = 5*nrow

    ALLOCATE (ifail(nrow),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"ifail",nrow*int_size)
    ALLOCATE (iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iwork",liwork*int_size)
    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work",lwork*dp_size)

    CALL dsyevx("V","I","U",nrow,matrix(1,1),nrow,vl,vu,1,ncol,abstol,&
                m,w,eigenvetors(1,1),nrow,work(1),lwork,&
                iwork(1),ifail(1),info)

!   *** Error handling ***

    IF (info /= 0) CALL stop_program(routineP,"Error in dsyevx")

    eigenvalues(1:ncol) = w(1:ncol)

!   *** Release work storage ***

    DEALLOCATE(w,ifail,iwork,work,STAT =istat)
    IF (istat /= 0) CALL stop_memory(routineP,"w,ifail,iwork,work")

  END SUBROUTINE fm_mol_syevx 

!***************************************************************************

  SUBROUTINE mol_mo_random(mo,mepos,num_pe,ncol,start_col,error)

    REAL(dp), DIMENSION(:,:), POINTER        :: mo
    INTEGER, INTENT(IN)                      :: mepos, num_pe
    INTEGER, INTENT(IN), OPTIONAL            :: start_col, ncol
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN ="mol_mo_random"

    INTEGER                                  :: icol, my_start, my_ncol, nrow
    INTEGER, DIMENSION(4), SAVE              :: ISEED
    LOGICAL, SAVE                            :: FIRST = .TRUE.
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

    failure = .FALSE.

    IF (FIRST) THEN
      ISEED(1)= mepos
      ISEED(2)= mepos + num_pe
      ISEED(3)= mepos + 2*num_pe
      ISEED(4)= mepos + 3*num_pe ! last one has to be odd
      IF(MOD(ISEED(4),2) .EQ. 0) ISEED(4)=ISEED(4)+1
      FIRST=.FALSE.
    ENDIF

    CPPostcondition(ASSOCIATED(mo),cp_failure_level,routineN,error,failure)
    my_start = 1
    IF(PRESENT(start_col)) my_start = start_col
    my_ncol = SIZE(mo,2)
    IF(PRESENT(ncol)) THEN
      CPPostcondition(my_ncol>=ncol,cp_failure_level,routineN,error,failure)
      my_ncol = ncol
    END IF

    nrow = SIZE(mo,1)
    DO icol = my_start ,my_start + my_ncol - 1
      CALL dlarnv(1,iseed,nrow,mo(1,icol))
    END DO 

  END SUBROUTINE mol_mo_random

!***************************************************************************
END MODULE  kg_gpw_fm_mol_utils




