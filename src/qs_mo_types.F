!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_mo_types [1.0] *
!!
!!   NAME
!!     qs_mo_types
!!
!!   FUNCTION
!!     Definition and initialisation of the mo data type.
!!
!!   AUTHOR
!!     Matthias Krack (09.05.2001)
!!
!!   MODIFICATION HISTORY
!!     - Adapted to the new QS environment data structure (02.04.2002,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_mo_types

! *****************************************************************************

  USE kinds, ONLY: wp => dp,&
                   wp_size => dp_size

  USE atomic_kind_types,    ONLY: atomic_kind_type,&
                                  get_atomic_kind,&
                                  get_atomic_kind_set
  USE basis_set_types,      ONLY: get_gto_basis_set,&
                                  gto_basis_set_type
  USE qs_blacs,             ONLY: allocate_blacs_matrix,&
                                  blacs_matrix_type,&
                                  copy_blacs_to_full_matrix,&
                                  deallocate_blacs_matrix,&
                                  get_blacs_info,&
                                  read_blacs_matrix,&
                                  write_blacs_matrix
  USE global_types,         ONLY: global_environment_type
  USE particle_types,       ONLY: particle_type
  USE qs_parser,            ONLY: close_file,&
                                  open_file
  USE sparse_matrix_types,  ONLY: allocate_matrix,&
                                  get_block_node,&
                                  real_matrix_type
  USE message_passing,      ONLY: mp_sum
  USE orbital_pointers,     ONLY: nso
  USE qs_environment_types, ONLY: get_qs_env,&
                                  qs_environment_type
  USE termination,          ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

  TYPE mo_set_type
    INTEGER                         :: homo,lfomo
    REAL(wp), DIMENSION(:), POINTER :: eigenvalues,occupation_numbers
    TYPE(blacs_matrix_type)         :: eigenvectors
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

  SUBROUTINE allocate_mo_set(mo_set,nao,nmo,nrow_block,ncol_block,context,&
                             globenv)

!   Purpose: Allocate a wavefunction data structure.

!   History: - Creation (15.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), INTENT(OUT)            :: mo_set
    INTEGER, INTENT(IN)                       :: context,nao,ncol_block,nmo,&
                                                 nrow_block

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_mo_set (MODULE qs_mo_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    mo_set%homo = 0
    mo_set%lfomo = 0

    ALLOCATE (mo_set%eigenvalues(nmo),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"mo_set%eigenvalues",nmo*wp_size)
    END IF
    mo_set%eigenvalues(:) = 0.0_wp

    ALLOCATE (mo_set%occupation_numbers(nmo),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"mo_set%occupation_numbers",nmo*wp_size)
    END IF
    mo_set%occupation_numbers(:) = 0.0_wp

    CALL allocate_blacs_matrix(new_matrix=mo_set%eigenvectors,&
                               nrow_global=nao,&
                               ncol_global=nmo,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="MO EIGENVECTORS",&
                               context=context,&
                               globenv=globenv)

  END SUBROUTINE allocate_mo_set

! *****************************************************************************

  SUBROUTINE deallocate_mo_set(mo_set)

!   Purpose: Deallocate a wavefunction data structure.

!   History: - Creation (15.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(mo_set_type), INTENT(INOUT) :: mo_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_mo_set (MODULE qs_mo_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    DEALLOCATE (mo_set%eigenvalues,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mo_set%eigenvalues")

    DEALLOCATE (mo_set%occupation_numbers,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mo_set%occupation_numbers")

    CALL deallocate_blacs_matrix(mo_set%eigenvectors)

  END SUBROUTINE deallocate_mo_set

! *****************************************************************************

  SUBROUTINE read_mo_set(mo_set,context,globenv)

!   Purpose: Read the MO eigenvectors from the restart file.

!   History: - Creation (15.05.2001, Matthias Krack)
!            - Parallel input (19.05.2001, MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), INTENT(OUT)            :: mo_set
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    CHARACTER(LEN=6)   :: extension
    CHARACTER(LEN=200) :: file_name
    INTEGER            :: mype,npe,restart_unit

!   ---------------------------------------------------------------------------

    CALL get_blacs_info(context=context,&
                        globenv=globenv,&
                        my_process_number=mype,&
                        number_of_processes=npe)

    IF (npe > 1) THEN
      WRITE (UNIT=extension,FMT="(I6)") mype
      file_name = TRIM(globenv%restart_file_name)//"."//ADJUSTL(extension)
    ELSE
      file_name = globenv%restart_file_name
    END IF

    CALL open_file(file_name=file_name,&
                   file_action="READ",&
                   file_form="UNFORMATTED",&
                   file_status="OLD",&
                   unit_number=restart_unit)

    CALL read_blacs_matrix(matrix=mo_set%eigenvectors,&
                           lunit=restart_unit,&
                           context=context,&
                           globenv=globenv)

    CALL close_file(unit_number=restart_unit)

  END SUBROUTINE read_mo_set

! *****************************************************************************

  SUBROUTINE write_mo_set_to_restart_unit(mo_set,context,globenv)

!   Purpose: Write the MO eigenvectors to the restart file.

!   History: - Creation (15.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), INTENT(IN)             :: mo_set
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    CHARACTER(LEN=6)   :: extension
    CHARACTER(LEN=200) :: file_name
    INTEGER            :: mype,npe,restart_unit

!   ---------------------------------------------------------------------------

    CALL get_blacs_info(context=context,&
                        globenv=globenv,&
                        my_process_number=mype,&
                        number_of_processes=npe)

    IF (npe > 1) THEN
      WRITE (UNIT=extension,FMT="(I6)") mype
      file_name = TRIM(globenv%restart_file_name)//"."//ADJUSTL(extension)
    ELSE
      file_name = globenv%restart_file_name
    END IF

    CALL open_file(file_name=file_name,&
                   file_action="WRITE",&
                   file_form="UNFORMATTED",&
                   file_status="REPLACE",&
                   unit_number=restart_unit)

    CALL write_blacs_matrix(matrix=mo_set%eigenvectors,&
                            lunit=restart_unit,&
                            context=context,&
                            globenv=globenv)

    CALL close_file(unit_number=restart_unit)

  END SUBROUTINE write_mo_set_to_restart_unit

! *****************************************************************************

  SUBROUTINE write_mo_set_to_output_unit(mo_set,before,after,context,qs_env,&
                                         globenv)

!   Purpose: Write the MO eigenvalues, MO occupation numbers and
!            MO eigenvectors.

!   History: - Creation (15.05.2001, Matthias Krack)

!   ***************************************************************************

!   List of variables:

!   after : Number of digits after point.
!   before: Number of digits before point.

!   ---------------------------------------------------------------------------

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), INTENT(INOUT)          :: mo_set
    TYPE(qs_environment_type), INTENT(IN)     :: qs_env
    INTEGER, INTENT(IN)                       :: after,before,context

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE write_mo_set_to_output_unit (MODULE qs_mo_types)"

!   *** Local variables ***

    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    CHARACTER(LEN=60) :: name
    CHARACTER(LEN=25) :: fmtstr1
    CHARACTER(LEN=22) :: fmtstr2
    CHARACTER(LEN=35) :: fmtstr3
    CHARACTER(LEN=2)  :: element_symbol
    INTEGER           :: from,iatom,icol,irow,iset,isgf,iso,ishell,istat,jcol,&
                         left,lshell,natom,ncol,nset,nsgf,output_unit,right,&
                         to,width

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set

    CHARACTER(LEN=6), DIMENSION(:), POINTER :: sgf_symbol
    INTEGER, DIMENSION(:), POINTER          :: nshell
    REAL(wp), DIMENSION(:,:), POINTER       :: matrix
    INTEGER, DIMENSION(:,:), POINTER        :: l

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%print%mo_eigenvalues) RETURN

    output_unit = globenv%scr

    IF (globenv%print%mo_eigenvectors) THEN
      NULLIFY (matrix)
      CALL copy_blacs_to_full_matrix(mo_set%eigenvectors,matrix,context,&
                                     globenv)
    END IF

    IF (.NOT.globenv%ionode) RETURN

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom=natom,&
                             nsgf=nsgf)

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

    WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(name)

!   *** Write the variable format strings ***

    width = before + after + 3
    ncol = INT(56/width)

    right = MAX((after-2),1)
    left =  width - right - 5

    WRITE (UNIT=fmtstr1(11:12),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr1(14:15),FMT="(I2)") left
    WRITE (UNIT=fmtstr1(21:22),FMT="(I2)") right

    WRITE (UNIT=fmtstr2(9:10),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr2(16:17),FMT="(I2)") width - 1
    WRITE (UNIT=fmtstr2(19:20),FMT="(I2)") after

    WRITE (UNIT=fmtstr3(22:23),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr3(29:30),FMT="(I2)") width - 1
    WRITE (UNIT=fmtstr3(32:33),FMT="(I2)") after

!   *** Write the matrix in the selected format ***

    DO icol=1,nsgf,ncol

      from = icol
      to = MIN((from+ncol-1),nsgf)

      WRITE (UNIT=output_unit,FMT=fmtstr1)&
        (jcol,jcol=from,to)
      WRITE (UNIT=output_unit,FMT=fmtstr2)&
        (mo_set%eigenvalues(jcol),jcol=from,to)
      WRITE (UNIT=output_unit,FMT="(A)") ""

      IF (globenv%print%mo_occupation_numbers) THEN
        WRITE (UNIT=output_unit,FMT=fmtstr2)&
          (mo_set%occupation_numbers(jcol),jcol=from,to)
        WRITE (UNIT=output_unit,FMT="(A)") ""
      END IF

      IF (globenv%print%mo_eigenvectors) THEN

        irow = 1

        DO iatom=1,natom

          IF (iatom /= 1) WRITE (UNIT=output_unit,FMT="(A)") ""

          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               element_symbol=element_symbol,&
                               orb_basis_set=orb_basis_set)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                 nset=nset,&
                                 nshell=nshell,&
                                 l=l,&
                                 sgf_symbol=sgf_symbol)

          isgf = 1

          DO iset=1,nset
            DO ishell=1,nshell(iset)
              lshell = l(ishell,iset)
              DO iso=1,nso(lshell)
                WRITE (UNIT=output_unit,FMT=fmtstr3)&
                  irow,iatom,element_symbol,sgf_symbol(isgf),&
                  (matrix(irow,jcol),jcol=from,to)
                isgf = isgf + 1
                irow = irow + 1
              END DO
            END DO
          END DO

        END DO

      END IF

    END DO

    WRITE (UNIT=output_unit,FMT="(/)")

!   *** Release work storage ***

    IF (globenv%print%mo_eigenvectors) THEN
      DEALLOCATE (matrix,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"matrix")
    END IF

  END SUBROUTINE write_mo_set_to_output_unit

! *****************************************************************************

END MODULE qs_mo_types
