!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/sparse_matrix_output [1.0] *
!!
!!   NAME
!!     sparse_matrix_output
!!
!!   FUNCTION
!!     Some output routines for sparse matrices
!!
!!   AUTHOR
!!     Matthias Krack (03.09.2001,25.06.2003)
!!
!!   MODIFICATION HISTORY
!!     JGH: outsourced to this routine from qs_overlap
!!
!!***
!******************************************************************************

MODULE sparse_matrix_output

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_should_output
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE orbital_pointers,                ONLY: nso
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE sparse_matrix_types,             ONLY: copy_local_sm_to_replicated_fm,&
                                             get_matrix_info,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sparse_matrix_output'

! *** Public subroutines ***

  PUBLIC :: write_fm_with_basis_info,&
            write_matrix_distribution,&
            write_sparse_matrix

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE write_fm_with_basis_info(blacs_matrix,before,after,qs_env,globenv,&
                                      first_row,last_row,first_col,last_col,output_unit)

!   Purpose: Print a spherical matrix of blacs type.

!   History: - Creation (12.06.2001,MK)
!            - Allow for printing of a sub-matrix (01.07.2003,MK)

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER                :: blacs_matrix
    INTEGER, INTENT(IN)                      :: before, after
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_col, last_col
    INTEGER, INTENT(IN)                      :: output_unit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_fm_with_basis_info'

    CHARACTER(LEN=60)                        :: matrix_name
    INTEGER                                  :: col1, col2, group, istat, &
                                                ncol_global, nrow_global, &
                                                nsgf, row1, row2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set

!   ---------------------------------------------------------------------------

    group = globenv%para_env%group
    IF (.NOT.ASSOCIATED(blacs_matrix)) RETURN
    CALL cp_fm_get_info(blacs_matrix,name=matrix_name,&
            nrow_global=nrow_global,ncol_global=ncol_global)

    ALLOCATE(matrix(nrow_global,ncol_global))
    CALL cp_fm_get_submatrix(blacs_matrix,matrix)


!   *** Get the matrix dimension and check the optional arguments ***

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

    IF (PRESENT(first_row)) THEN
      row1 = MAX(1,first_row)
    ELSE
      row1 = 1
    END IF

    IF (PRESENT(last_row)) THEN
      row2 = MIN(nsgf,last_row)
    ELSE
      row2 = nsgf
    END IF

    IF (PRESENT(first_col)) THEN
      col1 = MAX(1,first_col)
    ELSE
      col1 = 1
    END IF

    IF (PRESENT(last_col)) THEN
      col2 = MIN(nsgf,last_col)
    ELSE
      col2 = nsgf
    END IF

    CALL write_matrix(matrix,matrix_name,before,after,qs_env,globenv,&
                      row1,row2,col1,col2,output_unit)

!   *** Release work storage ***

    IF (ASSOCIATED(matrix)) THEN
      DEALLOCATE (matrix,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "matrix")
    END IF

  END SUBROUTINE write_fm_with_basis_info

! *****************************************************************************

  SUBROUTINE write_matrix(matrix,matrix_name,before,after,qs_env,globenv,&
                          first_row,last_row,first_col,last_col,output_unit)

!   Purpose: Write a matrix or a sub-matrix to the output unit.

!   History: - Creation (01.07.2003,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    CHARACTER(LEN=*), INTENT(IN)             :: matrix_name
    INTEGER, INTENT(IN)                      :: before, after
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER, INTENT(IN)                      :: first_row, last_row, &
                                                first_col, last_col, &
                                                output_unit

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'write_matrix'

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=25)                        :: fmtstr1
    CHARACTER(LEN=35)                        :: fmtstr2
    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: sgf_symbol
    INTEGER :: from, group, iatom, icol, irow, iset, isgf, ishell, iso, jcol, &
      l, left, natom, ncol, ndigits, nset, nsgf, right, to, width
    INTEGER, DIMENSION(:), POINTER           :: nshell
    INTEGER, DIMENSION(:, :), POINTER        :: lshell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!   ---------------------------------------------------------------------------

    group = globenv%para_env%group

    IF (output_unit>0) THEN
      CALL m_flush(output_unit)

      CALL get_qs_env(qs_env=qs_env,&
                      atomic_kind_set=atomic_kind_set,&
                      particle_set=particle_set)

      natom = SIZE(particle_set)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

!     *** Definition of the variable formats ***

      fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
      fmtstr2 = "(T2,2I5,2X,A2,1X,A8,   (1X,F  .  ))"

!     *** Write headline ***

      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(matrix_name)

!     *** Write the variable format strings ***

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

!     *** Write the matrix in the selected format ***

      DO icol=first_col,last_col,ncol
        from = icol
        to = MIN((from+ncol-1),last_col)
        WRITE (UNIT=output_unit,FMT=fmtstr1) (jcol,jcol=from,to)
        irow = 1
        DO iatom=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               element_symbol=element_symbol,&
                               orb_basis_set=orb_basis_set)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                 nset=nset,&
                                 nshell=nshell,&
                                 l=lshell,&
                                 sgf_symbol=sgf_symbol)
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
        END DO
      END DO

      WRITE (UNIT=output_unit,FMT="(/)")

    END IF

    CALL mp_sync(group)
    IF(output_unit>0) CALL m_flush(output_unit)

  END SUBROUTINE write_matrix

! *****************************************************************************

  SUBROUTINE write_matrix_distribution(matrix,force_env_section,globenv,error)

!   Purpose: Print the distribution of a sparse matrix.

!   History: - Creation (25.06.2003,MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_matrix_distribution'

    CHARACTER(LEN=80)                        :: matrix_name, matrix_symmetry
    INTEGER                                  :: group, ipe, istat, mype, &
                                                natom, nblock_sum, &
                                                nblock_tot, nelement_sum, &
                                                npe, nrow, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nblock, nelement
    LOGICAL                                  :: ionode
    REAL(KIND=dp)                            :: occupation
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
!   ---------------------------------------------------------------------------
!   *** Quick return, if no printing is requested ***

    IF (.NOT.BTEST(cp_print_key_should_output(logger%iter_info,&
         force_env_section,"PRINT%DISTRIBUTION",error=error),cp_p_file)) RETURN

    group = globenv%para_env%group
    ionode = globenv%para_env%ionode
    mype = globenv%para_env%mepos + 1
    npe = globenv%para_env%num_pe
    output_unit = globenv%scr

!   *** Allocate work storage ***

    ALLOCATE (nblock(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "nblock",npe*int_size)
    nblock(:) = 0

    ALLOCATE (nelement(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "nelement",npe*int_size)
    nelement(:) = 0

    matrix_symmetry = ""

    CALL get_matrix_info(matrix=matrix,&
                         matrix_name=matrix_name,&
                         matrix_symmetry=matrix_symmetry,&
                         nblock_allocated=nblock(mype),&
                         nelement_allocated=nelement(mype),&
                         nblock_row=natom,&
                         nrow=nrow)

    CALL mp_sum(nblock,group)
    CALL mp_sum(nelement,group)

    nblock_sum = SUM(nblock)
    nelement_sum = SUM(nelement)

    IF (TRIM(matrix_symmetry) == "symmetric") THEN
      nblock_tot = natom*(natom + 1)/2
    ELSE
      nblock_tot = natom*natom
    END IF

    occupation = 100.0_dp*REAL(nblock_sum,dp)/REAL(nblock_tot,dp)

    IF (ionode) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
        "DISTRIBUTION OF THE "//TRIM(matrix_name)
      WRITE (UNIT=output_unit,FMT="(/,T3,A,/,/,(I9,T27,I10,T55,I10))")&
        "Process    Number of matrix blocks   Number of matrix elements",&
        (ipe-1,nblock(ipe),nelement(ipe),ipe=1,npe)
      WRITE (UNIT=output_unit,FMT="(/,T7,A3,T27,I10,T55,I10)")&
        "Sum",nblock_sum,nelement_sum
      WRITE (UNIT=output_unit,FMT="(/,T7,A3,T27,I10,A,F5.1,A,T55,I10,A,F5.1,A)")&
        " of",nblock_tot," (",occupation," % occupation)"
    END IF

!   *** Release work storage ***

    DEALLOCATE (nblock,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "nblock")

    DEALLOCATE (nelement,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "nelement")

  END SUBROUTINE write_matrix_distribution

! *****************************************************************************

  SUBROUTINE write_sparse_matrix(sparse_matrix,before,after,qs_env,globenv,&
                                 first_row,last_row,first_col,last_col,scale,&
                                 output_unit)

!   Purpose: Print a spherical matrix of sparse_matrix_type.

!   History: - Creation (07.06.2000,MK)
!            - Allow for printing of a sub-matrix (01.07.2003,MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: sparse_matrix
    INTEGER, INTENT(IN)                      :: before, after
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_col, last_col
    REAL(dp), INTENT(IN), OPTIONAL           :: scale
    INTEGER, INTENT(IN)                      :: output_unit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_sparse_matrix'

    CHARACTER(LEN=40)                        :: matrix_symmetry
    CHARACTER(LEN=80)                        :: matrix_name
    INTEGER                                  :: col1, col2, group, istat, &
                                                nsgf, row1, row2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set

!   ---------------------------------------------------------------------------

    group = globenv%para_env%group

    CALL get_matrix_info(matrix=sparse_matrix,&
                         matrix_name=matrix_name,&
                         matrix_symmetry=matrix_symmetry)

    NULLIFY (matrix)

    CALL copy_local_sm_to_replicated_fm(sparse_matrix,matrix)

    IF (matrix_symmetry == "symmetric") THEN
      CALL symmetrize_matrix(matrix,"upper_to_lower")
    ELSE IF (matrix_symmetry == "antisymmetric") THEN
      CALL symmetrize_matrix(matrix,"anti_upper_to_lower")
    ELSE IF (matrix_symmetry /= "none") THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Invalid matrix structure",globenv)
    END IF

    CALL mp_sum(matrix,group)

!   *** Get the matrix dimension and check the optional arguments ***

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

    IF (PRESENT(first_row)) THEN
      row1 = MAX(1,first_row)
    ELSE
      row1 = 1
    END IF

    IF (PRESENT(last_row)) THEN
      row2 = MIN(nsgf,last_row)
    ELSE
      row2 = nsgf
    END IF

    IF (PRESENT(first_col)) THEN
      col1 = MAX(1,first_col)
    ELSE
      col1 = 1
    END IF

    IF (PRESENT(last_col)) THEN
      col2 = MIN(nsgf,last_col)
    ELSE
      col2 = nsgf
    END IF

    IF (PRESENT(scale)) THEN
      matrix=matrix*scale
    END IF

    CALL write_matrix(matrix,matrix_name,before,after,qs_env,globenv,&
                      row1,row2,col1,col2, output_unit)

    IF (ASSOCIATED(matrix)) THEN
      DEALLOCATE (matrix,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "matrix")
    END IF

  END SUBROUTINE write_sparse_matrix

! *****************************************************************************

END MODULE sparse_matrix_output
