!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Some output routines for sparse matrices
!> \par History
!>      JGH: outsourced to this routine from qs_overlap
!> \author Matthias Krack (03.09.2001,25.06.2003)
! *****************************************************************************
MODULE sparse_matrix_output

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE orbital_pointers,                ONLY: nso
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE sparse_matrix_types,             ONLY: copy_local_sm_to_replicated_fm,&
                                             get_matrix_info,&
                                             real_matrix_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sparse_matrix_output'

  ! *** Public subroutines ***
  PUBLIC ::  write_sparse_matrix

CONTAINS

! *****************************************************************************
!> \brief Write a matrix or a sub-matrix to the output unit (symmetric)
!> \author Creation (01.07.2003,MK)
! *****************************************************************************
  SUBROUTINE write_matrix_sym(matrix,matrix_name,before,after,qs_env,para_env,&
                          first_row,last_row,first_col,last_col,output_unit,error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    CHARACTER(LEN=*), INTENT(IN)             :: matrix_name
    INTEGER, INTENT(IN)                      :: before, after
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: first_row, last_row, &
                                                first_col, last_col, &
                                                output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_matrix_sym', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=25)                        :: fmtstr1
    CHARACTER(LEN=35)                        :: fmtstr2
    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: sgf_symbol
    INTEGER :: from, group, iatom, icol, irow, iset, isgf, ishell, iso, jcol, &
      l, left, natom, ncol, ndigits, nset, nsgf, right, stat, to, width
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: nshell
    INTEGER, DIMENSION(:, :), POINTER        :: lshell
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.
    group = para_env%group

    IF (output_unit>0) THEN
      CALL m_flush(output_unit)

      CALL get_qs_env(qs_env=qs_env,&
                      atomic_kind_set=atomic_kind_set,&
                      particle_set=particle_set,error=error)

      natom = SIZE(particle_set)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

      ALLOCATE (first_sgf(natom),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (last_sgf(natom),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL get_particle_set(particle_set=particle_set,&
                            first_sgf=first_sgf,&
                            last_sgf=last_sgf,error=error)

      ! *** Definition of the variable formats ***
      fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
      fmtstr2 = "(T2,2I5,2X,A2,1X,A8,   (1X,F  .  ))"

      ! *** Write headline ***
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(matrix_name)

      ! *** Write the variable format strings ***
      ndigits = after

      width = before + ndigits + 3
      ncol = INT(56/width)

      right = MAX((ndigits-2),1)
      left =  width - right - 5

      WRITE (UNIT=fmtstr1(11:12),FMT="(I2)") ncol
      WRITE (UNIT=fmtstr1(14:15),FMT="(I2)") left
      WRITE (UNIT=fmtstr1(21:22),FMT="(I2)") right

      WRITE (UNIT=fmtstr2(22:23),FMT="(I2)") ncol
      WRITE (UNIT=fmtstr2(29:30),FMT="(I2)") width - 1
      WRITE (UNIT=fmtstr2(32:33),FMT="(I2)") ndigits

      ! *** Write the matrix in the selected format ***
      DO icol=first_col,last_col,ncol
         from = icol
         to = MIN((from+ncol-1),last_col)
         WRITE (UNIT=output_unit,FMT=fmtstr1) (jcol,jcol=from,to)
         irow = 1
         DO iatom=1,natom
            NULLIFY(orb_basis_set)
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                 element_symbol=element_symbol, orb_basis_set=orb_basis_set)
            IF ( ASSOCIATED(orb_basis_set) ) THEN
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                    nset=nset,nshell=nshell,l=lshell,sgf_symbol=sgf_symbol)
               isgf = 1
               DO iset=1,nset
                  DO ishell=1,nshell(iset)
                     l = lshell(ishell,iset)
                     DO iso=1,nso(l)
                        IF ((irow >= first_row).AND.(irow <= last_row)) THEN
                           WRITE (UNIT=output_unit,FMT=fmtstr2)&
                                irow,iatom,element_symbol,sgf_symbol(isgf),&
                                (matrix(irow,jcol),jcol=from,to)
                        END IF
                        isgf = isgf + 1
                        irow = irow + 1
                     END DO
                  END DO
               END DO
               IF ((irow >= first_row).AND.(irow <= last_row)) THEN
                  WRITE (UNIT=output_unit,FMT="(A)")
               END IF
            ELSE
               DO iso=first_sgf(iatom),last_sgf(iatom)
                  IF ((irow >= first_row).AND.(irow <= last_row)) THEN
                     WRITE (UNIT=output_unit,FMT=fmtstr2)&
                          irow,iatom,element_symbol," ",&
                          (matrix(irow,jcol),jcol=from,to)
                  END IF
                  irow = irow + 1
               END DO
               IF ((irow >= first_row).AND.(irow <= last_row)) THEN
                  WRITE (UNIT=output_unit,FMT="(A)")
               END IF
            END IF
         END DO
      END DO
      
      WRITE (UNIT=output_unit,FMT="(/)")
      DEALLOCATE (first_sgf,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE (last_sgf,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL mp_sync(group)
    IF(output_unit>0) CALL m_flush(output_unit)

  END SUBROUTINE write_matrix_sym

! *****************************************************************************
!> \brief Write a matrix not necessarily symmetric (no index with atomic labels)
!> \author Teodoro Laino [tlaino] - 10.2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE write_matrix_gen(matrix,matrix_name,before,after,para_env,&
       first_row,last_row,first_col,last_col,output_unit,error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    CHARACTER(LEN=*), INTENT(IN)             :: matrix_name
    INTEGER, INTENT(IN)                      :: before, after
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: first_row, last_row, &
                                                first_col, last_col, &
                                                output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_matrix_gen', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=25)                        :: fmtstr1
    CHARACTER(LEN=35)                        :: fmtstr2
    INTEGER                                  :: from, group, icol, irow, &
                                                jcol, left, ncol, ndigits, &
                                                right, to, width
    LOGICAL                                  :: failure

    failure = .FALSE.
    group = para_env%group

    IF (output_unit>0) THEN
      CALL m_flush(output_unit)

      ! *** Definition of the variable formats ***
      fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
      fmtstr2 = "(T2, I5,        18X,   (1X,F  .  ))"

      ! *** Write headline ***
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(matrix_name)

      ! *** Write the variable format strings ***
      ndigits = after

      width = before + ndigits + 3
      ncol = INT(56/width)

      right = MAX((ndigits-2),1)
      left =  width - right - 5

      WRITE (UNIT=fmtstr1(11:12),FMT="(I2)") ncol
      WRITE (UNIT=fmtstr1(14:15),FMT="(I2)") left
      WRITE (UNIT=fmtstr1(21:22),FMT="(I2)") right

      WRITE (UNIT=fmtstr2(22:23),FMT="(I2)") ncol
      WRITE (UNIT=fmtstr2(29:30),FMT="(I2)") width - 1
      WRITE (UNIT=fmtstr2(32:33),FMT="(I2)") ndigits

      ! *** Write the matrix in the selected format ***
      DO icol=first_col,last_col,ncol
         from = icol
         to = MIN((from+ncol-1),last_col)
         WRITE (UNIT=output_unit,FMT=fmtstr1) (jcol,jcol=from,to)
         irow = 1
         DO irow = first_row, last_row
            WRITE (UNIT=output_unit,FMT=fmtstr2)&
                 irow,(matrix(irow,jcol),jcol=from,to)
         END DO
      END DO
      
      WRITE (UNIT=output_unit,FMT="(/)")
    END IF

    CALL mp_sync(group)
    IF(output_unit>0) CALL m_flush(output_unit)

  END SUBROUTINE write_matrix_gen

! *****************************************************************************
!> \brief Print a spherical matrix of sparse_matrix_type.
!> \author Creation (07.06.2000,MK)
!>      Allow for printing of a sub-matrix (01.07.2003,MK)
! *****************************************************************************
  SUBROUTINE write_sparse_matrix(sparse_matrix,before,after,qs_env,para_env,&
                                 first_row,last_row,first_col,last_col,scale,&
                                 output_unit,error)

    TYPE(real_matrix_type), POINTER          :: sparse_matrix
    INTEGER, INTENT(IN)                      :: before, after
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_col, last_col
    REAL(dp), INTENT(IN), OPTIONAL           :: scale
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_sparse_matrix', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=40)                        :: matrix_symmetry
    CHARACTER(LEN=80)                        :: matrix_name
    INTEGER                                  :: col1, col2, dim_col, dim_row, &
                                                group, row1, row2, stat
    LOGICAL                                  :: failure, print_sym
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set

    failure = .FALSE.
    group = para_env%group
    CALL get_matrix_info(matrix=sparse_matrix,&
                         matrix_name=matrix_name,&
                         matrix_symmetry=matrix_symmetry)

    NULLIFY (matrix)

    CALL copy_local_sm_to_replicated_fm(sparse_matrix,matrix)

    IF (matrix_symmetry == "symmetric") THEN
       CALL symmetrize_matrix(matrix,"upper_to_lower")
       print_sym = .TRUE.
    ELSE IF (matrix_symmetry == "antisymmetric") THEN
       CALL symmetrize_matrix(matrix,"anti_upper_to_lower")
       print_sym = .TRUE.
    ELSE IF (matrix_symmetry /= "none") THEN
       print_sym = .FALSE.
    END IF

    CALL mp_sum(matrix,group)

    ! *** Get the matrix dimension and check the optional arguments ***
    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,error=error)
    dim_row = SIZE(matrix,1)
    dim_col = SIZE(matrix,2)

    IF (PRESENT(first_row)) THEN
      row1 = MAX(1,first_row)
    ELSE
      row1 = 1
    END IF

    IF (PRESENT(last_row)) THEN
      row2 = MIN(dim_row,last_row)
    ELSE
      row2 = dim_row
    END IF

    IF (PRESENT(first_col)) THEN
      col1 = MAX(1,first_col)
    ELSE
      col1 = 1
    END IF

    IF (PRESENT(last_col)) THEN
      col2 = MIN(dim_col,last_col)
    ELSE
      col2 = dim_col
    END IF

    IF (PRESENT(scale)) THEN
      matrix=matrix*scale
    END IF

    IF (print_sym) THEN
       CALL write_matrix_sym(matrix,matrix_name,before,after,qs_env,para_env,&
            row1,row2,col1,col2, output_unit,error=error)
    ELSE
       CALL write_matrix_gen(matrix,matrix_name,before,after,para_env,&
            row1,row2,col1,col2, output_unit,error=error)
    END IF

    IF (ASSOCIATED(matrix)) THEN
       DEALLOCATE (matrix,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE write_sparse_matrix

END MODULE sparse_matrix_output
