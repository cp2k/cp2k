!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/om_utilities [1.0] *
!!
!!   NAME
!!     om_utilities
!!
!!   FUNCTION
!!      Utility routines for the printing of the operator matrix
!!
!!   AUTHOR
!!     Matthias Krack (12.07.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE om_utilities

! *****************************************************************************

! Index:

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE message_passing, ONLY: mp_sum

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: write_cartesian_matrix,&
            write_spherical_matrix,&
            write_g_matrix

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE write_cartesian_matrix(block_matrix,before,after,globenv)

!   Purpose: Write a Cartesian matrix.

!   History: - Creation (07.06.2000, Matthias Krack)

!   ***************************************************************************

!   List of variables:

!   after   : Number of digits after point.
!   before  : Number of digits before point.
!   col1    : First column of matrix to be printed.
!   col2    : Last column of matrix to be printed.
!   matrix  : Integral matrix to be printed.
!   row1    : First row of matrix to be printed.

!   ---------------------------------------------------------------------------

    USE atomic_kinds,     ONLY: kind_info
    USE atoms,            ONLY: atom_info,natom,ncgf
    USE global_types,     ONLY: global_environment_type
    USE mathlib,          ONLY: symmetrize_matrix
    USE matrix_types,     ONLY: add_block_node,&
                                allocate_matrix,&
                                copy_matrix,&
                                deallocate_matrix,&
                                get_block_node,&
                                get_matrix_info,&
                                real_matrix_type
    USE orbital_pointers, ONLY: nco

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(real_matrix_type), POINTER           :: block_matrix
    INTEGER, INTENT(IN)                       :: after,before

!   *** Local variables ***

    TYPE(real_matrix_type), POINTER :: full_matrix
    CHARACTER(LEN=60)               :: name
    CHARACTER(LEN=40)               :: symmetry
    CHARACTER(LEN=25)               :: fmtstr1
    CHARACTER(LEN=35)               :: fmtstr2
    INTEGER                         :: from,iatom,icgf,ico,icol,ikind,irow,&
                                       iset,ishell,jcol,l,left,ncol,&
                                       output_unit,right,to,width

    REAL(wp), DIMENSION(:,:), POINTER :: matrix

!   ---------------------------------------------------------------------------

    output_unit = globenv%scr

    CALL get_matrix_info(matrix=block_matrix,&
                         matrix_name=name,&
                         matrix_symmetry=symmetry)

    NULLIFY (full_matrix)

    CALL allocate_matrix(matrix=full_matrix,&
                         nrow=ncgf,&
                         ncol=ncgf,&
                         matrix_name="WORK",&
                         matrix_symmetry=symmetry)

    CALL copy_matrix(block_matrix,full_matrix)

    CALL get_block_node(matrix=full_matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=matrix)

    IF (symmetry == "symmetric") THEN
      CALL symmetrize_matrix(matrix,"upper_to_lower")
    ELSE IF (symmetry == "antisymmetric") THEN
      CALL symmetrize_matrix(matrix,"anti_upper_to_lower")
    END IF

#if defined(__parallel)
    CALL mp_sum(matrix,globenv%group)
#endif

    IF (.NOT.globenv%ionode) RETURN

!   *** Definition of the variable formats ***

    fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
    fmtstr2 = "(T2,2I5,2X,A2,1X,A8,   (1X,F  .  ))"

!   *** Write headline ***

    WRITE (output_unit,"(/,/,T2,A)") TRIM(name)

!   *** Write the variable format strings ***

    width = before + after + 3
    ncol = INT(56/width)

    right = MAX((after-2),1)
    left =  width - right - 5

    WRITE (fmtstr1(11:12),"(I2)") ncol
    WRITE (fmtstr1(14:15),"(I2)") left
    WRITE (fmtstr1(21:22),"(I2)") right

    WRITE (fmtstr2(22:23),"(I2)") ncol
    WRITE (fmtstr2(29:30),"(I2)") width - 1
    WRITE (fmtstr2(32:33),"(I2)") after

!   *** Write the matrix in the selected format ***

    DO icol=1,ncgf,ncol
      from = icol
      to = MIN((from+ncol-1),ncgf)
      WRITE (output_unit,fmtstr1) (jcol,jcol=from,to)
      irow = 1
      DO iatom=1,natom
        IF (iatom /= 1) WRITE (output_unit,"(A)") ""
        ikind = atom_info(iatom)%kind
        icgf = 1
        DO iset=1,kind_info(ikind)%orb_basis_set%nset
          DO ishell=1,kind_info(ikind)%orb_basis_set%nshell(iset)
            l = kind_info(ikind)%orb_basis_set%l(ishell,iset)
            DO ico=1,nco(l)
              WRITE (output_unit,fmtstr2)&
                irow,iatom,kind_info(ikind)%element_symbol,&
                kind_info(ikind)%orb_basis_set%cgf_symbol(icgf),&
                (matrix(irow,jcol),jcol=from,to)
              icgf = icgf + 1
              irow = irow + 1
            END DO
          END DO
        END DO
      END DO
    END DO

    WRITE (output_unit,"(/)")

!   *** Release work storage ***

    CALL deallocate_matrix(full_matrix)

  END SUBROUTINE write_cartesian_matrix

! *****************************************************************************

  SUBROUTINE write_spherical_matrix(block_matrix,before,after,globenv,&
                                    matrix_name)

!   Purpose: Write a spherical matrix.

!   History: - Creation (07.06.2000, Matthias Krack)

!   ***************************************************************************

!   List of variables:

!   after   : Number of digits after point.
!   before  : Number of digits before point.
!   col1    : First column of matrix to be printed.
!   col2    : Last column of matrix to be printed.
!   matrix  : Integral matrix to be printed.
!   row1    : First row of matrix to be printed.

!   ---------------------------------------------------------------------------

    USE atomic_kinds,     ONLY: kind_info
    USE atoms,            ONLY: atom_info,natom,nsgf
    USE global_types,     ONLY: global_environment_type
    USE mathlib,          ONLY: symmetrize_matrix
    USE matrix_types,     ONLY: add_block_node,&
                                allocate_matrix,&
                                deallocate_matrix,&
                                copy_matrix,&
                                get_block_node,&
                                get_matrix_info,&
                                real_matrix_type
    USE orbital_pointers, ONLY: nso

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(real_matrix_type), POINTER           :: block_matrix
    INTEGER, INTENT(IN)                       :: after,before
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: matrix_name

!   *** Local variables ***

    TYPE(real_matrix_type), POINTER :: full_matrix
    CHARACTER(LEN=60)               :: name
    CHARACTER(LEN=40)               :: symmetry
    CHARACTER(LEN=25)               :: fmtstr1
    CHARACTER(LEN=35)               :: fmtstr2
    INTEGER                         :: from,iatom,icol,ikind,irow,iset,isgf,&
                                       iso,ishell,jcol,l,left,ncol,&
                                       output_unit,right,to,width

    REAL(wp), DIMENSION(:,:), POINTER :: matrix

!   ---------------------------------------------------------------------------

    output_unit = globenv%scr

    CALL get_matrix_info(matrix=block_matrix,&
                         matrix_name=name,&
                         matrix_symmetry=symmetry)

    IF (PRESENT(matrix_name)) name = matrix_name

    NULLIFY (full_matrix)

    CALL allocate_matrix(matrix=full_matrix,&
                         nrow=nsgf,&
                         ncol=nsgf,&
                         matrix_name="WORK",&
                         matrix_symmetry=symmetry)

    CALL copy_matrix(block_matrix,full_matrix)

    CALL get_block_node(matrix=full_matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=matrix)

    IF (symmetry == "symmetric") THEN
      CALL symmetrize_matrix(matrix,"upper_to_lower")
    ELSE IF (symmetry == "antisymmetric") THEN
      CALL symmetrize_matrix(matrix,"anti_upper_to_lower")
    END IF

#if defined(__parallel)
    CALL mp_sum(matrix,globenv%group)
#endif

    IF (.NOT.globenv%ionode) RETURN

!   *** Definition of the variable formats ***

    fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
    fmtstr2 = "(T2,2I5,2X,A2,1X,A8,   (1X,F  .  ))"

!   *** Write headline ***

    WRITE (output_unit,"(/,/,T2,A)") TRIM(name)

!   *** Write the variable format strings ***

    width = before + after + 3
    ncol = INT(56/width)

    right = MAX((after-2),1)
    left =  width - right - 5

    WRITE (fmtstr1(11:12),"(I2)") ncol
    WRITE (fmtstr1(14:15),"(I2)") left
    WRITE (fmtstr1(21:22),"(I2)") right

    WRITE (fmtstr2(22:23),"(I2)") ncol
    WRITE (fmtstr2(29:30),"(I2)") width - 1
    WRITE (fmtstr2(32:33),"(I2)") after

!   *** Write the matrix in the selected format ***

    DO icol=1,nsgf,ncol
      from = icol
      to = MIN((from+ncol-1),nsgf)
      WRITE (output_unit,fmtstr1) (jcol,jcol=from,to)
      irow = 1
      DO iatom=1,natom
        IF (iatom /= 1) WRITE (output_unit,"(A)") ""
        ikind = atom_info(iatom)%kind
        isgf = 1
        DO iset=1,kind_info(ikind)%orb_basis_set%nset
          DO ishell=1,kind_info(ikind)%orb_basis_set%nshell(iset)
            l = kind_info(ikind)%orb_basis_set%l(ishell,iset)
            DO iso=1,nso(l)
              WRITE (output_unit,fmtstr2)&
                irow,iatom,kind_info(ikind)%element_symbol,&
                kind_info(ikind)%orb_basis_set%sgf_symbol(isgf),&
                (matrix(irow,jcol),jcol=from,to)
              isgf = isgf + 1
              irow = irow + 1
            END DO
          END DO
        END DO
      END DO
    END DO

    WRITE (output_unit,"(/)")

!   *** Release work storage ***

    CALL deallocate_matrix(full_matrix)

  END SUBROUTINE write_spherical_matrix

! *****************************************************************************

  SUBROUTINE write_g_matrix(block_matrix,before,after,globenv)

!   Purpose: Write the Cartesian G matrix.

!   History: - Creation (07.06.2000, Matthias Krack)

!   ***************************************************************************

!   List of variables:

!   after   : Number of digits after point.
!   before  : Number of digits before point.
!   col1    : First column of matrix to be printed.
!   col2    : Last column of matrix to be printed.
!   matrix  : Integral matrix to be printed.
!   row1    : First row of matrix to be printed.

!   ---------------------------------------------------------------------------

    USE atomic_kinds,     ONLY: kind_info
    USE atoms,            ONLY: atom_info,natom,ncgf_aux
    USE global_types,     ONLY: global_environment_type
    USE mathlib,          ONLY: symmetrize_matrix
    USE matrix_types,     ONLY: add_block_node,&
                                allocate_matrix,&
                                copy_matrix,&
                                deallocate_matrix,&
                                get_block_node,&
                                get_matrix_info,&
                                real_matrix_type
    USE orbital_pointers, ONLY: nco

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(real_matrix_type), POINTER           :: block_matrix
    INTEGER, INTENT(IN)                       :: after,before

!   *** Local variables ***

    TYPE(real_matrix_type), POINTER :: full_matrix
    CHARACTER(LEN=60)               :: name
    CHARACTER(LEN=40)               :: symmetry
    CHARACTER(LEN=25)               :: fmtstr1
    CHARACTER(LEN=35)               :: fmtstr2
    INTEGER                         :: from,iatom,icgf,ico,icol,ikind,irow,&
                                       iset,ishell,jcol,l,left,ncol,&
                                       output_unit,right,to,width

    REAL(wp), DIMENSION(:,:), POINTER :: matrix

!   ---------------------------------------------------------------------------

    output_unit = globenv%scr

    CALL get_matrix_info(matrix=block_matrix,&
                         matrix_name=name,&
                         matrix_symmetry=symmetry)

    NULLIFY (full_matrix)

    CALL allocate_matrix(matrix=full_matrix,&
                         nrow=ncgf_aux,&
                         ncol=ncgf_aux,&
                         matrix_name="WORK",&
                         matrix_symmetry=symmetry)

    CALL copy_matrix(block_matrix,full_matrix)

    CALL get_block_node(matrix=full_matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=matrix)

    IF (symmetry == "symmetric") THEN
      CALL symmetrize_matrix(matrix,"upper_to_lower")
    ELSE IF (symmetry == "antisymmetric") THEN
      CALL symmetrize_matrix(matrix,"anti_upper_to_lower")
    END IF

#if defined(__parallel)
    CALL mp_sum(matrix,globenv%group)
#endif

    IF (.NOT.globenv%ionode) RETURN

!   *** Definition of the variable formats ***

    fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
    fmtstr2 = "(T2,2I5,2X,A2,1X,A8,   (1X,F  .  ))"

!   *** Write headline ***

    WRITE (output_unit,"(/,/,T2,A)") TRIM(name)

!   *** Write the variable format strings ***

    width = before + after + 3
    ncol = INT(56/width)

    right = MAX((after-2),1)
    left =  width - right - 5

    WRITE (fmtstr1(11:12),"(I2)") ncol
    WRITE (fmtstr1(14:15),"(I2)") left
    WRITE (fmtstr1(21:22),"(I2)") right

    WRITE (fmtstr2(22:23),"(I2)") ncol
    WRITE (fmtstr2(29:30),"(I2)") width - 1
    WRITE (fmtstr2(32:33),"(I2)") after

!   *** Write the matrix in the selected format ***

    DO icol=1,ncgf_aux,ncol
      from = icol
      to = MIN((from+ncol-1),ncgf_aux)
      WRITE (output_unit,fmtstr1) (jcol,jcol=from,to)
      irow = 1
      DO iatom=1,natom
        IF (iatom /= 1) WRITE (output_unit,"(A)") ""
        ikind = atom_info(iatom)%kind
        icgf = 1
        DO iset=1,kind_info(ikind)%aux_basis_set%nset
          DO ishell=1,kind_info(ikind)%aux_basis_set%nshell(iset)
            l = kind_info(ikind)%aux_basis_set%l(ishell,iset)
            DO ico=1,nco(l)
              WRITE (output_unit,fmtstr2)&
                irow,iatom,kind_info(ikind)%element_symbol,&
                kind_info(ikind)%aux_basis_set%cgf_symbol(icgf),&
                (matrix(irow,jcol),jcol=from,to)
              icgf = icgf + 1
              irow = irow + 1
            END DO
          END DO
        END DO
      END DO
    END DO

    WRITE (output_unit,"(/)")

!   *** Release work storage ***

    CALL deallocate_matrix(full_matrix)

  END SUBROUTINE write_g_matrix

! *****************************************************************************

END MODULE om_utilities
