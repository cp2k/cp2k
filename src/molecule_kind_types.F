!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/molecule_kind_types [1.0] *
!!
!!   NAME
!!     molecule_kind_types
!!
!!   FUNCTION
!!     Define the molecule kind structure types and the corresponding
!!     functionality
!!
!!   AUTHOR
!!     MK (22.08.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE molecule_kind_types

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE force_field_types,               ONLY: bend_kind_type,&
                                             bond_kind_type,&
                                             impr_kind_type,&
                                             onfo_kind_type,&
                                             torsion_kind_type,&
                                             ub_kind_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_size
  USE memory_utilities,                ONLY: reallocate
  USE qs_parser,                       ONLY: finish_parser,&
                                             read_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE string_utilities,                ONLY: ascii_to_string,&
                                             compress,&
                                             string_to_ascii
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "molecule_kind_types"

! *** Define the derived structure types ***

  TYPE atom_type
    TYPE(atomic_kind_type), POINTER :: atomic_kind
    CHARACTER(LEN=default_string_length)  :: name
  END TYPE atom_type

  TYPE bond_type
    TYPE(bond_kind_type), POINTER :: bond_kind
    INTEGER :: a,b
  END TYPE bond_type

  TYPE bend_type
    TYPE(bend_kind_type), POINTER :: bend_kind
    INTEGER :: a,b,c
  END TYPE bend_type

  TYPE ub_type
    TYPE(ub_kind_type), POINTER :: ub_kind
    INTEGER :: a,b,c
  END TYPE ub_type

  TYPE torsion_type
    TYPE(torsion_kind_type), POINTER :: torsion_kind
    INTEGER :: a,b,c,d
  END TYPE torsion_type

  TYPE impr_type
    TYPE(impr_kind_type), POINTER :: impr_kind
    INTEGER :: a,b,c,d
  END TYPE impr_type

  TYPE onfo_type
    TYPE(onfo_kind_type), POINTER :: onfo_kind
    INTEGER :: a,b,c,d
  END TYPE onfo_type

  TYPE dist_constraint_type
    REAL(KIND = dp)               :: dab
    INTEGER                       :: a,b
  END TYPE dist_constraint_type

  TYPE g3x3_constraint_type
    INTEGER                :: a,b,c
    REAL(KIND = dp)        :: dab,dac,dbc
  END TYPE g3x3_constraint_type

  TYPE g4x6_constraint_type
    INTEGER                :: a,b,c,d
    REAL(KIND = dp)        :: dab,dac,dbc,dad,dbd,dcd
  END TYPE g4x6_constraint_type

  TYPE molecule_kind_type
    PRIVATE
    TYPE(atom_type), DIMENSION(:), POINTER            :: atom_list
    TYPE(bend_kind_type), DIMENSION(:), POINTER       :: bend_kind_set
    TYPE(bend_type), DIMENSION(:), POINTER            :: bend_list
    TYPE(bond_kind_type), DIMENSION(:), POINTER       :: bond_kind_set
    TYPE(bond_type), DIMENSION(:), POINTER            :: bond_list
    TYPE(ub_kind_type), DIMENSION(:), POINTER         :: ub_kind_set
    TYPE(ub_type), DIMENSION(:), POINTER              :: ub_list
    TYPE(torsion_kind_type), DIMENSION(:), POINTER    :: torsion_kind_set
    TYPE(torsion_type), DIMENSION(:), POINTER         :: torsion_list
    TYPE(impr_kind_type), DIMENSION(:), POINTER       :: impr_kind_set
    TYPE(impr_type), DIMENSION(:), POINTER            :: impr_list
    TYPE(onfo_kind_type), DIMENSION(:), POINTER       :: onfo_kind_set
    TYPE(onfo_type), DIMENSION(:), POINTER            :: onfo_list
    TYPE(dist_constraint_type), DIMENSION(:), POINTER :: dist_list
    TYPE(g3x3_constraint_type), DIMENSION(:), POINTER :: g3x3_list
    TYPE(g4x6_constraint_type), DIMENSION(:), POINTER :: g4x6_list
    CHARACTER(LEN=default_string_length)                    :: name,&
                                                         nhcopt
    REAL(KIND = dp), DIMENSION(:), POINTER            :: charge,&
                                                         mass
    INTEGER                                           :: kind_number,&
                                                         natom,&
                                                         nbond,&
                                                         nbend,&
                                                         nimpr,&
                                                         nonfo,&
                                                         ntorsion,&
                                                         nub,&
                                                         ndist,&
                                                         ng3x3,&
                                                         ng4x6,&
                                                         nmolecule
    INTEGER, DIMENSION(:), POINTER                    :: molecule_list
  END TYPE molecule_kind_type

! *** Public subroutines ***

  PUBLIC :: allocate_molecule_kind_set,&
            deallocate_molecule_kind_set,&
            get_molecule_kind,&
            get_molecule_kind_set,&
            read_molecule_kind,&
            read_molecule_kind_set,&
            set_molecule_kind,&
            write_molecule_kind,&
            write_molecule_kind_set

! *** Public data types ***

  PUBLIC :: atom_type,&
            bend_type,&
            bond_type,&
            ub_type,&
            torsion_type,&
            impr_type,&
            onfo_type,&
            dist_constraint_type,&
            g3x3_constraint_type,&
            g4x6_constraint_type,&
            molecule_kind_type

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_molecule_kind_set(molecule_kind_set,nmolecule_kind)

!   Purpose: Allocate and initialize a molecule kind set.

!   History: - Creation (22.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    INTEGER, INTENT(IN)                      :: nmolecule_kind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_molecule_kind_set"

    INTEGER                                  :: imolecule_kind, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind_set)) THEN
      CALL deallocate_molecule_kind_set(molecule_kind_set)
    END IF

    ALLOCATE (molecule_kind_set(nmolecule_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "molecule_kind_set",&
                                     nmolecule_kind*int_size)

    DO imolecule_kind=1,nmolecule_kind
      NULLIFY (molecule_kind_set(imolecule_kind)%atom_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%bond_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%bend_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%dist_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%g3x3_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%g4x6_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%torsion_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%bond_kind_set)
      NULLIFY (molecule_kind_set(imolecule_kind)%bend_kind_set)
      NULLIFY (molecule_kind_set(imolecule_kind)%torsion_kind_set)
      NULLIFY (molecule_kind_set(imolecule_kind)%charge)
      NULLIFY (molecule_kind_set(imolecule_kind)%mass)
      molecule_kind_set(imolecule_kind)%name = ""
      molecule_kind_set(imolecule_kind)%nhcopt = ""
      molecule_kind_set(imolecule_kind)%kind_number = imolecule_kind
      molecule_kind_set(imolecule_kind)%natom = 0
      molecule_kind_set(imolecule_kind)%nbend = 0
      molecule_kind_set(imolecule_kind)%nbond = 0
      molecule_kind_set(imolecule_kind)%ndist = 0
      molecule_kind_set(imolecule_kind)%ng3x3 = 0
      molecule_kind_set(imolecule_kind)%ng4x6 = 0
      molecule_kind_set(imolecule_kind)%nmolecule = 0
      molecule_kind_set(imolecule_kind)%ntorsion = 0
      NULLIFY (molecule_kind_set(imolecule_kind)%molecule_list)
    END DO

  END SUBROUTINE allocate_molecule_kind_set

! *****************************************************************************

  SUBROUTINE deallocate_molecule_kind_set(molecule_kind_set)

!   Purpose: Deallocate a molecule kind set.

!   History: - Creation (22.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set

    CHARACTER(LEN=*), PARAMETER :: &
      routine_name = "deallocate_molecule_kind_set"

    INTEGER                                  :: imolecule_kind, istat, &
                                                nmolecule_kind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind_set)) THEN

      nmolecule_kind = SIZE(molecule_kind_set)

      DO imolecule_kind=1,nmolecule_kind

        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%atom_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%atom_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routine_name,module_name,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%atom_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%bend_kind_set)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%bend_kind_set,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routine_name,module_name,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%bend_kind_set")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%bend_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%bend_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routine_name,module_name,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%bend_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%bond_kind_set)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%bond_kind_set,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routine_name,module_name,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%bond_kind_set")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%bond_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%bond_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routine_name,module_name,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%bond_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%dist_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%dist_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routine_name,module_name,__LINE__,&
                           "molecule_kind_set(imolecule_kind)%dist_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%g3x3_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%g3x3_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routine_name,module_name,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%g3x3_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%g4x6_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%g4x6_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routine_name,module_name,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%g4x6_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%torsion_kind_set)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%torsion_kind_set,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routine_name,module_name,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%torsion_kind_set")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%torsion_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%torsion_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routine_name,module_name,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%torsion_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%molecule_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%molecule_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routine_name,module_name,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%molecule_list")
           END IF
        ENDIF
      END DO

      DEALLOCATE (molecule_kind_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "molecule_kind_set")

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule_kind_set is not associated "//&
                        "and cannot be deallocated")

    END IF

  END SUBROUTINE deallocate_molecule_kind_set

! *****************************************************************************

  SUBROUTINE get_molecule_kind(molecule_kind,atom_list,bond_list,bend_list,&
                               ub_list,impr_list,onfo_list,dist_list,&
                               g3x3_list,g4x6_list,torsion_list,name,nhcopt,&
                               mass,charge,kind_number,natom,nbend,nbond,nub,&
                               nimpr,nonfo,nconstraint,&
                               ndist,ng3x3,ng4x6,nmolecule,nsgf,ntorsion,&
                               molecule_list,&
                               bond_kind_set,bend_kind_set,ub_kind_set,&
                               impr_kind_set,onfo_kind_set,torsion_kind_set)

!   Purpose: Get informations about a molecule kind.

!   History: - Creation (27.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atom_list
    TYPE(bond_kind_type), DIMENSION(:), OPTIONAL,POINTER    :: bond_kind_set
    TYPE(bend_kind_type), DIMENSION(:), OPTIONAL,POINTER    :: bend_kind_set
    TYPE(ub_kind_type), DIMENSION(:), OPTIONAL,POINTER      :: ub_kind_set
    TYPE(impr_kind_type), DIMENSION(:), OPTIONAL,POINTER    :: impr_kind_set
    TYPE(onfo_kind_type), DIMENSION(:), OPTIONAL,POINTER    :: onfo_kind_set
    TYPE(torsion_kind_type), DIMENSION(:), OPTIONAL,POINTER :: torsion_kind_set
    TYPE(bond_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bond_list
    TYPE(bend_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bend_list
    TYPE(ub_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ub_list
    TYPE(impr_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: impr_list
    TYPE(onfo_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: onfo_list
    TYPE(dist_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: dist_list
    TYPE(g3x3_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: g3x3_list
    TYPE(g4x6_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: g4x6_list
    TYPE(torsion_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: torsion_list
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name, nhcopt
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: mass, charge
    INTEGER, INTENT(OUT), OPTIONAL           :: kind_number, natom, nbend, &
                                                nbond, nconstraint, ndist, &
                                                ng3x3, ng4x6, nmolecule, nsgf, &
                                                ntorsion, nub, nimpr, nonfo
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: molecule_list


    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_molecule_kind"

    INTEGER                                  :: iatom, n

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind)) THEN

      IF (PRESENT(atom_list)) atom_list => molecule_kind%atom_list
      IF (PRESENT(bend_list)) bend_list => molecule_kind%bend_list
      IF (PRESENT(bond_list)) bond_list => molecule_kind%bond_list
      IF (PRESENT(impr_list)) impr_list => molecule_kind%impr_list
      IF (PRESENT(onfo_list)) onfo_list => molecule_kind%onfo_list
      IF (PRESENT(ub_list)) ub_list => molecule_kind%ub_list
      IF (PRESENT(bond_kind_set)) bond_kind_set => molecule_kind%bond_kind_set
      IF (PRESENT(bend_kind_set)) bend_kind_set => molecule_kind%bend_kind_set
      IF (PRESENT(ub_kind_set)) ub_kind_set => molecule_kind%ub_kind_set
      IF (PRESENT(impr_kind_set)) impr_kind_set => molecule_kind%impr_kind_set
      IF (PRESENT(onfo_kind_set)) onfo_kind_set => molecule_kind%onfo_kind_set
      IF (PRESENT(torsion_kind_set)) torsion_kind_set => molecule_kind%torsion_kind_set
      IF (PRESENT(dist_list)) dist_list => molecule_kind%dist_list
      IF (PRESENT(g3x3_list)) g3x3_list => molecule_kind%g3x3_list
      IF (PRESENT(g4x6_list)) g4x6_list => molecule_kind%g4x6_list
      IF (PRESENT(torsion_list)) torsion_list => molecule_kind%torsion_list
      IF (PRESENT(name)) name = molecule_kind%name
      IF (PRESENT(nhcopt)) nhcopt = molecule_kind%nhcopt
      IF (PRESENT(mass)) mass => molecule_kind%mass
      IF (PRESENT(charge)) charge => molecule_kind%charge
      IF (PRESENT(kind_number)) kind_number = molecule_kind%kind_number
      IF (PRESENT(natom)) natom = molecule_kind%natom
      IF (PRESENT(nbend)) nbend = molecule_kind%nbend
      IF (PRESENT(nbond)) nbond = molecule_kind%nbond
      IF (PRESENT(nub)) nub = molecule_kind%nub
      IF (PRESENT(nimpr)) nimpr = molecule_kind%nimpr
      IF (PRESENT(nonfo)) nonfo = molecule_kind%nonfo
      IF (PRESENT(nconstraint)) nconstraint = molecule_kind%ndist +&
                                              3*molecule_kind%ng3x3 +&
                                              6*molecule_kind%ng4x6
      IF (PRESENT(ndist)) ndist = molecule_kind%ndist
      IF (PRESENT(ng3x3)) ng3x3= molecule_kind%ng3x3
      IF (PRESENT(ng4x6)) ng4x6= molecule_kind%ng4x6
      IF (PRESENT(nmolecule)) nmolecule = molecule_kind%nmolecule
      IF (PRESENT(ntorsion)) ntorsion= molecule_kind%ntorsion
      IF (PRESENT(nsgf)) THEN
        nsgf = 0
        DO iatom=1,SIZE(molecule_kind%atom_list)
          CALL get_atomic_kind(atomic_kind=&
                                 molecule_kind%atom_list(iatom)%atomic_kind,&
                               nsgf=n)
          nsgf = nsgf + n
        END DO
      END IF
      IF (PRESENT(molecule_list)) molecule_list => molecule_kind%molecule_list

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule_kind is not associated")

    END IF

  END SUBROUTINE get_molecule_kind

! *****************************************************************************

  SUBROUTINE get_molecule_kind_set(molecule_kind_set,maxatom,natom,&
                                   nconstraint,nmolecule)

!   Purpose: Get informations about a molecule kind set.

!   History: - Creation (27.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    INTEGER, INTENT(OUT), OPTIONAL           :: maxatom, natom, nconstraint, &
                                                nmolecule

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_molecule_kind_set"

    INTEGER                                  :: imolecule_kind, na, nc, nm, &
                                                nmolecule_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind_set)) THEN

      IF (PRESENT(maxatom)) maxatom = 0
      IF (PRESENT(natom)) natom = 0
      IF (PRESENT(nconstraint)) nconstraint = 0
      IF (PRESENT(nmolecule)) nmolecule = 0

      nmolecule_kind = SIZE(molecule_kind_set)

      DO imolecule_kind=1,nmolecule_kind

        molecule_kind => molecule_kind_set(imolecule_kind)

        CALL get_molecule_kind(molecule_kind=molecule_kind,&
                               natom=na,&
                               nconstraint=nc,&
                               nmolecule=nm)

        IF (PRESENT(maxatom)) maxatom = MAX(maxatom,na)
        IF (PRESENT(natom)) natom = natom + na*nm
        IF (PRESENT(nconstraint)) nconstraint = nconstraint + nc*nm
        IF (PRESENT(nmolecule)) nmolecule = nmolecule + nm

      END DO

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule_kind_set is not associated")

    END IF

  END SUBROUTINE get_molecule_kind_set

! *****************************************************************************

  SUBROUTINE read_molecule_kind(molecule_kind,atomic_kind_set,globenv)

!   Purpose: Read a molecule kind data set from the input file.

!   History: - Creation (02.09.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_molecule_kind", &
      start_section_label = "MOLECULE"
    INTEGER, PARAMETER                       :: nblock = 1

    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=default_string_length)     :: atomic_kind_name, &
                                                current_name, keyword
    INTEGER                                  :: iatom, ikind, istat, natom
    INTEGER, DIMENSION(:), POINTER           :: kind_list
    INTEGER, DIMENSION(:, :), POINTER        :: name_list
    LOGICAL                                  :: atoms_read, name_found, &
                                                section_found
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!   ---------------------------------------------------------------------------

    NULLIFY (kind_list)
    NULLIFY (name_list)

    CALL reallocate(kind_list,1,nblock)
    CALL reallocate(name_list,1,default_string_length,1,nblock)

    natom = 0
    atoms_read = .FALSE.

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label=start_section_label,&
                      start_section_name=molecule_kind%name,&
                      section_found=section_found)

    IF (section_found) THEN
      DO
        test_result = test_object(newline=.TRUE.)
        IF (test_result == "EOS") THEN
          EXIT
        ELSE IF (test_result == "EOF") THEN
          CALL stop_parser(routine_name,module_name,__LINE__,"EOF")
        ELSE
          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE (TRIM(keyword))
          CASE ("ATOMS")
            IF (atoms_read) THEN
              CALL stop_parser(routine_name,module_name,__LINE__,&
                               "A subsection ATOMS was already read for "//&
                               "the current molecule kind")
            END IF
            DO
              test_result = test_object(newline=.TRUE.)
              IF (test_result == "EOS") THEN
                CALL stop_parser(routine_name,module_name,__LINE__,&
                                 "Unexpected end of input section <"//&
                                 start_section_label//" "//&
                                 TRIM(molecule_kind%name)//"> found. "//&
                                 "The subsection ATOMS wasn't terminated by END")
              ELSE IF (test_result == "EOF") THEN
                CALL stop_parser(routine_name,module_name,__LINE__,"EOF")
              ELSE
                CALL read_object(atomic_kind_name,lower_to_upper=.TRUE.,&
                                 newline=.TRUE.)
                IF (TRIM(atomic_kind_name) == "END") EXIT
                natom = natom + 1
                IF (natom > SIZE(kind_list)) THEN
                  CALL reallocate(kind_list,1,natom+nblock)
                  CALL reallocate(name_list,1,default_string_length,1,natom+nblock)
                END IF
                name_found = .FALSE.
                DO ikind=1,SIZE(atomic_kind_set)
                  atomic_kind => atomic_kind_set(ikind)
                  CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                       name=current_name)
                  IF (TRIM(current_name) == TRIM(atomic_kind_name)) THEN
                    name_found = .TRUE.
                    kind_list(natom) = ikind
                    EXIT
                  END IF
                END DO
                IF (.NOT.name_found) THEN
                  CALL stop_parser(routine_name,module_name,__LINE__,&
                                   "An atomic kind with the name <"//&
                                   TRIM(atomic_kind_name)//&
                                   "> could not be found.")
                END IF
                test_result = test_object(newline=.FALSE.)
                IF (test_result /= "EOL") THEN
                  CALL read_object(current_name,lower_to_upper=.FALSE.,&
                                   newline=.FALSE.)
                ELSE
                  WRITE (UNIT=current_name,FMT="(I6)") natom
                END IF
                CALL compress(current_name)
                CALL string_to_ascii(current_name,name_list(:,natom))
              END IF
            END DO
            atoms_read = .TRUE.
          CASE ("BENDS")
            CALL stop_parser(routine_name,module_name,__LINE__,&
                             "Reading BENDS: to be done")
          CASE ("BONDS")
            CALL stop_parser(routine_name,module_name,__LINE__,&
                             "Reading BONDS: to be done")
          CASE ("TORSIONS")
            CALL stop_parser(routine_name,module_name,__LINE__,&
                             "Reading TORSIONS: to be done")
          CASE DEFAULT
            CALL stop_parser(routine_name,module_name,__LINE__,"INVALID_KEYWORD")
          END SELECT
        END IF
      END DO
    ELSE
      CALL stop_parser(routine_name,module_name,__LINE__,"SECTION_NOT_FOUND")
    END IF

    IF (.NOT.atoms_read) THEN
      CALL stop_parser(routine_name,module_name,__LINE__,&
                       "No subsection ATOMS found for the molecule kind <"//&
                       TRIM(molecule_kind%name)//">")
    END IF

    CALL finish_parser()

!   *** Initialize the molecule kind data structure ***

    molecule_kind%natom = natom

    ALLOCATE (molecule_kind%atom_list(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "molecule_kind%atom_list",natom*int_size)

    DO iatom=1,natom
      ikind = kind_list(iatom)
      molecule_kind%atom_list(iatom)%atomic_kind => atomic_kind_set(ikind)
      CALL ascii_to_string(name_list(:,iatom),molecule_kind%atom_list(iatom)%name)
    END DO

!   *** Release work storage ***

    IF (ASSOCIATED(kind_list)) THEN
      DEALLOCATE (kind_list,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "kind_list")
    END IF

    IF (ASSOCIATED(name_list)) THEN
      DEALLOCATE (name_list,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "name_list")
    END IF

  END SUBROUTINE read_molecule_kind

! *****************************************************************************

  SUBROUTINE read_molecule_kind_set(molecule_kind_set,atomic_kind_set,globenv)

!   Purpose: Read all needed molecule kind data sets from the input file.

!   History: - Creation (12.09.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_molecule_kind_set"

    INTEGER                                  :: imolecule_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind_set)) THEN

      DO imolecule_kind=1,SIZE(molecule_kind_set)

        molecule_kind => molecule_kind_set(imolecule_kind)

        CALL read_molecule_kind(molecule_kind,atomic_kind_set,globenv)

      END DO

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule_kind_set is not associated")

    END IF

  END SUBROUTINE read_molecule_kind_set

! *****************************************************************************

  SUBROUTINE set_molecule_kind(molecule_kind,name,mass,charge,kind_number,&
                               molecule_list,atom_list,nhcopt,nbond,bond_list,&
                               nbend,bend_list,nub,ub_list,nimpr,impr_list,&
                               nonfo,onfo_list,ntorsion,torsion_list,&
                               ndist,dist_list,ng3x3,g3x3_list,ng4x6,g4x6_list,&
                               bond_kind_set,bend_kind_set,ub_kind_set,&
                               torsion_kind_set,impr_kind_set,onfo_kind_set)

!   Purpose: Set the components of a molecule kind.

!   History: - Creation (27.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: mass, charge
    INTEGER, INTENT(IN), OPTIONAL            :: kind_number
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: molecule_list
    TYPE(atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atom_list
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: nhcopt
    INTEGER, INTENT(IN), OPTIONAL            :: nbond
    TYPE(bond_kind_type), DIMENSION(:), OPTIONAL,POINTER :: bond_kind_set
    TYPE(bend_kind_type), DIMENSION(:), OPTIONAL,POINTER :: bend_kind_set
    TYPE(ub_kind_type), DIMENSION(:), OPTIONAL,POINTER   :: ub_kind_set
    TYPE(torsion_kind_type), DIMENSION(:), OPTIONAL,POINTER :: torsion_kind_set
    TYPE(impr_kind_type), DIMENSION(:), OPTIONAL,POINTER :: impr_kind_set
    TYPE(onfo_kind_type), DIMENSION(:), OPTIONAL,POINTER :: onfo_kind_set
    TYPE(bond_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bond_list
    INTEGER, INTENT(IN), OPTIONAL            :: nbend
    TYPE(bend_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bend_list
    INTEGER, INTENT(IN), OPTIONAL            :: nub
    TYPE(ub_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ub_list
    INTEGER, INTENT(IN), OPTIONAL            :: ntorsion
    TYPE(torsion_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: torsion_list
    INTEGER, INTENT(IN), OPTIONAL            :: nimpr
    TYPE(impr_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: impr_list
    INTEGER, INTENT(IN), OPTIONAL            :: nonfo
    TYPE(onfo_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: onfo_list
    INTEGER, INTENT(IN), OPTIONAL            :: ndist
    TYPE(dist_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: dist_list
    INTEGER, INTENT(IN), OPTIONAL            :: ng3x3
    TYPE(g3x3_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: g3x3_list
    INTEGER, INTENT(IN), OPTIONAL            :: ng4x6
    TYPE(g4x6_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: g4x6_list


    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_molecule_kind"

    INTEGER                                  :: n

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind)) THEN

      IF (PRESENT(atom_list)) THEN
!        n = SIZE(atom_list)
!        IF (n > 0) THEN
!          IF (ASSOCIATED(molecule_kind%atom_list)) THEN
!            DO i=1,SIZE(molecule_kind%atom_list)
!              NULLIFY (molecule_kind%atom_list(i)%atomic_kind)
!            END DO
!            DEALLOCATE (molecule_kind%atom_list,STAT=istat)
!            IF (istat /= 0) THEN
!              CALL stop_memory(routine_name,module_name,__LINE__,&
!                               "molecule_kind%atom_list")
!!            END IF
!          END IF
!          molecule_kind%natom = n
!          ALLOCATE (molecule_kind%atom_list(n),STAT=istat)
!          IF (istat /= 0) THEN
!            CALL stop_memory(routine_name,module_name,__LINE__,&
!                             "molecule_kind%atom_list",n*int_size)
!          END IF
!          DO i=1,n
!            molecule_kind%atom_list(i)%atomic_kind => atom_list(i)%atomic_kind
!            molecule_kind%atom_list(i)%name = atom_list(i)%name
!          END DO
!        ELSE
!          CALL stop_program(routine_name,module_name,__LINE__,&
!                            "An invalid atom_list was supplied")
!        END IF
        n = SIZE(atom_list)
        molecule_kind%natom = n
        molecule_kind%atom_list => atom_list
      END IF
      IF (PRESENT(name)) molecule_kind%name = name
      IF (PRESENT(nhcopt)) molecule_kind%nhcopt = nhcopt
      IF (PRESENT(mass)) molecule_kind%mass => mass
      IF (PRESENT(charge)) molecule_kind%mass => charge
      IF (PRESENT(kind_number)) molecule_kind%kind_number = kind_number
      IF (PRESENT(nbond)) molecule_kind%nbond = nbond
      IF (PRESENT(bond_list)) molecule_kind%bond_list => bond_list
      IF (PRESENT(nbend)) molecule_kind%nbend = nbend
      IF (PRESENT(bend_list)) molecule_kind%bend_list => bend_list
      IF (PRESENT(nub)) molecule_kind%nub = nub
      IF (PRESENT(ub_list)) molecule_kind%ub_list => ub_list
      IF (PRESENT(ntorsion)) molecule_kind%ntorsion = ntorsion
      IF (PRESENT(torsion_list)) molecule_kind%torsion_list => torsion_list
      IF (PRESENT(nimpr)) molecule_kind%nimpr = nimpr
      IF (PRESENT(impr_list)) molecule_kind%impr_list => impr_list
      IF (PRESENT(nonfo)) molecule_kind%nonfo = nonfo
      IF (PRESENT(onfo_list)) molecule_kind%onfo_list => onfo_list
      IF (PRESENT(ndist)) molecule_kind%ndist = ndist
      IF (PRESENT(dist_list)) molecule_kind%dist_list => dist_list
      IF (PRESENT(ng3x3)) molecule_kind%ng3x3 = ng3x3
      IF (PRESENT(g3x3_list)) molecule_kind%g3x3_list => g3x3_list
      IF (PRESENT(ng4x6)) molecule_kind%ng4x6 = ng4x6
      IF (PRESENT(g4x6_list)) molecule_kind%g4x6_list => g4x6_list
      IF (PRESENT(bond_kind_set)) molecule_kind%bond_kind_set => bond_kind_set
      IF (PRESENT(bend_kind_set)) molecule_kind%bend_kind_set => bend_kind_set
      IF (PRESENT(ub_kind_set)) molecule_kind%ub_kind_set => ub_kind_set
      IF (PRESENT(torsion_kind_set)) molecule_kind%torsion_kind_set => torsion_kind_set
      IF (PRESENT(impr_kind_set)) molecule_kind%impr_kind_set => impr_kind_set
      IF (PRESENT(onfo_kind_set)) molecule_kind%onfo_kind_set => onfo_kind_set
      IF (PRESENT(molecule_list)) THEN
!        n = SIZE(molecule_list)
!        IF (n > 0) THEN
!          IF (ASSOCIATED(molecule_kind%molecule_list)) THEN
!            DEALLOCATE (molecule_kind%molecule_list,STAT=istat)
!            IF (istat /= 0) THEN
!              CALL stop_memory(routine_name,module_name,__LINE__,&
!                               "molecule_kind%molecule_list")
!            END IF
!          END IF
!          molecule_kind%nmolecule = n
!          ALLOCATE (molecule_kind%molecule_list(n),STAT=istat)
!          IF (istat /= 0) THEN
!            CALL stop_memory(routine_name,module_name,__LINE__,&
!                             "molecule_kind%molecule_list",n*int_size)
!          END IF
!          molecule_kind%molecule_list(:) = molecule_list(:)
!        ELSE
!          CALL stop_program(routine_name,module_name,__LINE__,&
!                            "An invalid molecule_list was supplied")
!        END IF
        n = SIZE(molecule_list)
        molecule_kind%nmolecule = n
        molecule_kind%molecule_list => molecule_list
      END IF

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule_kind is not associated")

    END IF

  END SUBROUTINE set_molecule_kind

! *****************************************************************************

  SUBROUTINE write_molecule_kind(molecule_kind,globenv)

!   Purpose: Write a molecule kind data set to the output unit.

!   History: - Creation (24.09.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_molecule_kind"

    CHARACTER(LEN=default_string_length)     :: name
    INTEGER                                  :: iatom, imolecule, natom, &
                                                nmolecule, output_unit
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%kind_information) THEN

      IF (ASSOCIATED(molecule_kind)) THEN

        output_unit = globenv%scr

        natom = SIZE(molecule_kind%atom_list)
        nmolecule = SIZE(molecule_kind%molecule_list)

        WRITE (UNIT=output_unit,FMT="(/,T2,I5,A,T50,A,T75,I6,/,T22,A)")&
          molecule_kind%kind_number,&
          ". Molecule kind: "//TRIM(molecule_kind%name),&
          "Number of atoms:    ",natom,&
          "Atom         Atomic kind name             Atom name (label)"
        DO iatom=1,natom
          atomic_kind => molecule_kind%atom_list(iatom)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,name=name)
          WRITE (UNIT=output_unit,FMT="(T20,I6,(7X,A18,11X,A19))")&
            iatom,TRIM(name),TRIM(molecule_kind%atom_list(iatom)%name)
        END DO
        WRITE (UNIT=output_unit,FMT="(/,T9,A,I6,/,T9,A,(T30,10I5))")&
          "Number of molecules:",nmolecule,"Molecule list:",&
          (molecule_kind%molecule_list(imolecule),imolecule=1,nmolecule)

!MK to be done
!         "Number of bonds:    ",SIZE(molecule_kind%bond_list)
!         "Number of bends:    ",SIZE(molecule_kind%bend_list)
!         "Number of torsions: ",SIZE(molecule_kind%torsion_list)

      ELSE

        CALL stop_program(routine_name,module_name,__LINE__,&
                          "The pointer molecule_kind is not associated",&
                          globenv)

      END IF

    END IF

  END SUBROUTINE write_molecule_kind

! *****************************************************************************

  SUBROUTINE write_molecule_kind_set(molecule_kind_set,globenv)

!   Purpose: Write a moleculeatomic kind set data set to the output unit.

!   History: - Creation (24.09.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_molecule_kind_set"

    INTEGER                                  :: imolecule_kind, &
                                                nmolecule_kind, output_unit
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%kind_information) THEN

      IF (ASSOCIATED(molecule_kind_set)) THEN

        output_unit = globenv%scr

        WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") "MOLECULE KIND INFORMATION"

        nmolecule_kind = SIZE(molecule_kind_set)

        DO imolecule_kind=1,nmolecule_kind
          molecule_kind => molecule_kind_set(imolecule_kind)
          CALL write_molecule_kind(molecule_kind,globenv)
        END DO

      ELSE

        CALL stop_program(routine_name,module_name,__LINE__,&
                          "The pointer molecule_kind_set is not associated",&
                        globenv)

      END IF

    END IF

  END SUBROUTINE write_molecule_kind_set

! *****************************************************************************


END MODULE molecule_kind_types
