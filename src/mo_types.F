!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/mo_types [1.0] *
!!
!!   NAME
!!     mo_types
!!
!!   FUNCTION
!!     Definition and initialisation of the mo data type.
!!
!!   AUTHOR
!!     Matthias Krack (09.05.2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE mo_types

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE global_types,    ONLY: global_environment_type
  USE input_utilities, ONLY: close_file,&
                             open_file
  USE matrix_types,    ONLY: allocate_matrix,&
                             deallocate_matrix,&
                             get_block_node,&
                             real_matrix_set_type

  IMPLICIT NONE

  PRIVATE

  TYPE mo_set_type
    INTEGER                         :: homo,lfomo
    REAL(wp), DIMENSION(:), POINTER :: eigenvalues,occupation_numbers
    TYPE(real_matrix_set_type)      :: eigenvectors
  END TYPE mo_set_type

! *** Public data types ***

  PUBLIC :: mo_set_type

! *** Public subroutines ***

  PUBLIC :: allocate_mo_set,&
            deallocate_mo_set,&
            read_mo_set,&
            write_mo_set

! *****************************************************************************

  INTERFACE write_mo_set
    MODULE PROCEDURE write_mo_set_to_output_unit,write_mo_set_to_restart_unit
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_mo_set(mo_set,nao,nmo)

    USE memory_utilities, ONLY: reallocate

    TYPE(mo_set_type), INTENT(OUT) :: mo_set
    INTEGER, INTENT(IN)            :: nao,nmo

!   ---------------------------------------------------------------------------

    mo_set%homo = 0
    mo_set%lfomo = 0

    NULLIFY (mo_set%eigenvalues)
    mo_set%eigenvalues => reallocate(mo_set%eigenvalues,1,nmo)

    NULLIFY (mo_set%occupation_numbers)
    mo_set%occupation_numbers => reallocate(mo_set%occupation_numbers,1,nmo)

    NULLIFY (mo_set%eigenvectors%matrix)
    CALL allocate_matrix(matrix=mo_set%eigenvectors%matrix,&
                         nrow=nao,&
                         ncol=nmo,&
                         matrix_name="MO EIGENVECTORS",&
                         matrix_symmetry="none")

  END SUBROUTINE allocate_mo_set

! *****************************************************************************

  SUBROUTINE deallocate_mo_set(mo_set)

    TYPE(mo_set_type), INTENT(INOUT) :: mo_set

!   ---------------------------------------------------------------------------

    DEALLOCATE (mo_set%eigenvalues)

    DEALLOCATE (mo_set%occupation_numbers)

    CALL deallocate_matrix(matrix=mo_set%eigenvectors%matrix)

  END SUBROUTINE deallocate_mo_set

! *****************************************************************************

  SUBROUTINE read_mo_set(mo_set,globenv)

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), INTENT(IN)             :: mo_set

!   *** Local variables ***

    INTEGER :: iao,imo,restart_unit

    REAL(wp), DIMENSION(:,:), POINTER :: mo_eigenvectors

!   ---------------------------------------------------------------------------

    CALL get_block_node(matrix=mo_set%eigenvectors%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=mo_eigenvectors)

    CALL open_file(file_name=globenv%restart_file_name,&
                   file_action="READ",&
                   file_form="UNFORMATTED",&
                   file_status="OLD",&
                   unit_number=restart_unit)

    READ (UNIT=restart_unit) ((mo_eigenvectors(iao,imo),&
                               iao=1,SIZE(mo_eigenvectors,1)),&
                               imo=1,SIZE(mo_eigenvectors,2))

    CALL close_file(unit_number=restart_unit)

  END SUBROUTINE read_mo_set

! *****************************************************************************

  SUBROUTINE write_mo_set_to_restart_unit(mo_set,globenv)

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), INTENT(IN)             :: mo_set

!   *** Local variables ***

    INTEGER :: iao,imo,restart_unit

    REAL(wp), DIMENSION(:,:), POINTER :: mo_eigenvectors

!   ---------------------------------------------------------------------------

    CALL get_block_node(matrix=mo_set%eigenvectors%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=mo_eigenvectors)

    CALL open_file(file_name=globenv%restart_file_name,&
                   file_action="WRITE",&
                   file_form="UNFORMATTED",&
                   file_status="REPLACE",&
                   unit_number=restart_unit)

    WRITE (UNIT=restart_unit) ((mo_eigenvectors(iao,imo),&
                                iao=1,SIZE(mo_eigenvectors,1)),&
                                imo=1,SIZE(mo_eigenvectors,2))

    CALL close_file(unit_number=restart_unit)

  END SUBROUTINE write_mo_set_to_restart_unit

! *****************************************************************************

  SUBROUTINE write_mo_set_to_output_unit(mo_set,before,after,globenv)

!   Purpose: Write the MO eigenvalues, MO occupation numbers and
!            MO eigenvectors.

!   History: - Creation (15.05.2001, Matthias Krack)

!   ***************************************************************************

!   List of variables:

!   after : Number of digits after point.
!   before: Number of digits before point.

!   ---------------------------------------------------------------------------

    USE atomic_kinds,     ONLY: kind_info
    USE atoms,            ONLY: atom_info,natom,nsgf
    USE orbital_pointers, ONLY: nso

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), INTENT(IN)             :: mo_set
    INTEGER, INTENT(IN)                       :: after,before

!   *** Local variables ***

    CHARACTER(LEN=60) :: name
    CHARACTER(LEN=25) :: fmtstr1
    CHARACTER(LEN=22) :: fmtstr2
    CHARACTER(LEN=35) :: fmtstr3
    INTEGER           :: from,iatom,icol,ikind,irow,iset,isgf,iso,ishell,jcol,&
                         l,left,ncol,output_unit,right,to,width

    REAL(wp), DIMENSION(:,:), POINTER :: matrix

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%print%mo_eigenvalues) RETURN

    output_unit = globenv%scr

    CALL get_block_node(matrix=mo_set%eigenvectors%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=matrix)

    IF (.NOT.globenv%ionode) RETURN

!   *** Definition of the variable formats ***

    fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
    fmtstr2 = "(T2,23X,  (1X,F  .  ))"
    fmtstr3 = "(T2,2I5,2X,A2,1X,A8,   (1X,F  .  ))"

!   *** Write headline ***

    IF (globenv%print%mo_eigenvectors) THEN
      name = "MO EIGENVALUES, MO OCCUPATION NUMBERS AND MO EIGENVECTORS"
    ELSE IF (globenv%print%mo_occupation_numbers) THEN
      name = "MO EIGENVALUES AND MO OCCUPATION NUMBERS"
    ELSE
      name = "MO EIGENVALUES"
    END IF

    WRITE (output_unit,"(/,/,T2,A)") TRIM(name)

!   *** Write the variable format strings ***

    width = before + after + 3
    ncol = INT(56/width)

    right = MAX((after-2),1)
    left =  width - right - 5

    WRITE (fmtstr1(11:12),"(I2)") ncol
    WRITE (fmtstr1(14:15),"(I2)") left
    WRITE (fmtstr1(21:22),"(I2)") right

    WRITE (fmtstr2(9:10),"(I2)") ncol
    WRITE (fmtstr2(16:17),"(I2)") width - 1
    WRITE (fmtstr2(19:20),"(I2)") after

    WRITE (fmtstr3(22:23),"(I2)") ncol
    WRITE (fmtstr3(29:30),"(I2)") width - 1
    WRITE (fmtstr3(32:33),"(I2)") after

!   *** Write the matrix in the selected format ***

    DO icol=1,nsgf,ncol

      from = icol
      to = MIN((from+ncol-1),nsgf)

      WRITE (output_unit,fmtstr1)&
        (jcol,jcol=from,to)
      WRITE (output_unit,fmtstr2)&
        (mo_set%eigenvalues(jcol),jcol=from,to)
      WRITE (output_unit,"(A)") ""

      IF (globenv%print%mo_occupation_numbers) THEN
        WRITE (output_unit,fmtstr2)&
          (mo_set%occupation_numbers(jcol),jcol=from,to)
        WRITE (output_unit,"(A)") ""
      END IF

      IF (globenv%print%mo_eigenvectors) THEN

        irow = 1

        DO iatom=1,natom

          IF (iatom /= 1) WRITE (output_unit,"(A)") ""

          ikind = atom_info(iatom)%kind
          isgf = 1

          DO iset=1,kind_info(ikind)%orb_basis_set%nset
            DO ishell=1,kind_info(ikind)%orb_basis_set%nshell(iset)
              l = kind_info(ikind)%orb_basis_set%l(ishell,iset)
              DO iso=1,nso(l)
                WRITE (output_unit,fmtstr3)&
                  irow,iatom,kind_info(ikind)%element_symbol,&
                  kind_info(ikind)%orb_basis_set%sgf_symbol(isgf),&
                  (matrix(irow,jcol),jcol=from,to)
                isgf = isgf + 1
                irow = irow + 1
              END DO
            END DO
          END DO

        END DO

      END IF

    END DO

    WRITE (output_unit,"(/)")

  END SUBROUTINE write_mo_set_to_output_unit

! *****************************************************************************

END MODULE mo_types
