!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/atoms [1.0] *
!!
!!   NAME
!!     atoms
!!
!!   FUNCTION
!!     Define the data structure for the atomic information.
!!
!!   AUTHOR
!!     Matthias Krack (22.08.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE atoms

! *****************************************************************************

! Index:

! SUBROUTINE init_atom_info
! SUBROUTINE read_atomic_coordinates(start_section,end_section,globenv)
! SUBROUTINE write_atomic_coordinates(lunit)
! SUBROUTINE write_interatomic_distances(lunit)

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE atomic_kinds, ONLY: allocate_kind_info,&
                          init_kind_info,&
                          kind_info,&
                          max_kind_label_length

  IMPLICIT NONE

  PRIVATE

  TYPE atom_info_type
    INTEGER :: atom,first_cgf,first_cgf_aux,first_sgf,first_sgf_aux,&
               kind,last_cgf,last_cgf_aux,last_sgf,last_sgf_aux
    REAL(wp), DIMENSION(3) :: f,r,r_pbc,v
  END TYPE atom_info_type

  TYPE(atom_info_type), DIMENSION(:), POINTER :: atom_info

  INTEGER :: natom = 0,&
             ncgf = 0,&
             ncgf_aux = 0,&
             nsgf = 0,&
             nsgf_aux = 0

! *** Public subroutines ***

  PUBLIC :: init_atom_info,&
            read_atomic_coordinates,&
            write_atomic_coordinates,&
            write_interatomic_distances

! *** Public data types ***

  PUBLIC :: atom_info_type

! *** Public variables ***

  PUBLIC :: atom_info,natom,ncgf,ncgf_aux,nsgf,nsgf_aux

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE init_atom_info

!   Purpose: Complete the initialization of the atomic data sets.

!   History: - Creation (09.10.2000, Matthias Krack)

!   ***************************************************************************

    USE method_specifications, ONLY: allchem

    INTEGER :: iatom,ikind

!   ---------------------------------------------------------------------------

    ncgf = 0
    nsgf = 0

    DO iatom=1,natom
      ikind = atom_info(iatom)%kind
      atom_info(iatom)%first_cgf = ncgf + 1
      ncgf = ncgf + kind_info(ikind)%orb_basis_set%ncgf
      atom_info(iatom)%last_cgf = ncgf
      atom_info(iatom)%first_sgf = nsgf + 1
      nsgf = nsgf + kind_info(ikind)%orb_basis_set%nsgf
      atom_info(iatom)%last_sgf = nsgf
    END DO

    IF (allchem) THEN

      ncgf_aux = 0
      nsgf_aux = 0

      DO iatom=1,natom
        ikind = atom_info(iatom)%kind
        atom_info(iatom)%first_cgf_aux = ncgf_aux + 1
        ncgf_aux = ncgf_aux + kind_info(ikind)%aux_basis_set%ncgf
        atom_info(iatom)%last_cgf_aux = ncgf_aux
        atom_info(iatom)%first_sgf_aux = nsgf_aux + 1
        nsgf_aux = nsgf_aux + kind_info(ikind)%aux_basis_set%nsgf
        atom_info(iatom)%last_sgf_aux = nsgf_aux
      END DO

    END IF

  END SUBROUTINE init_atom_info

! *****************************************************************************

  SUBROUTINE read_atomic_coordinates(start_section,end_section,globenv)

!   Purpose: Read the input section with the atomic coordinates.

!   History: - Creation (19.05.2000, Matthias Krack)

!   ***************************************************************************

    USE cell_parameters,  ONLY: abc,pbc,scaled_coordinates,unit_of_length
    USE global_types,     ONLY: global_environment_type
    USE input_utilities,  ONLY: finish_input_session,&
                                read_object,&
                                search,&
                                start_input_session,&
                                test_object
    USE memory_utilities, ONLY: reallocate
    USE periodic_table,   ONLY: ptable
    USE string_utilities, ONLY: integer_to_string,&
                                string_to_integer,&
                                uppercase
    USE termination,      ONLY: stop_memory,&
                                stop_program

    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=*), INTENT(IN)              :: start_section,end_section

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_atomic_coordinates (MODULE atoms)"

    INTEGER, PARAMETER :: nblock_atom = 100,&
                          nblock_kind = 10

!   *** Local variables ***

    CHARACTER(LEN=max_kind_label_length) :: current_label,label,string
    CHARACTER(LEN=200)                   :: message
    CHARACTER(LEN=11)                    :: test_result
    CHARACTER(LEN=8)                     :: unit_name
    INTEGER                              :: iatom,ikind,istat,nkind,z
    LOGICAL                              :: found,new_kind

    INTEGER, DIMENSION(:), POINTER :: kind_of,natom_of_kind

    REAL(wp), DIMENSION(:,:), POINTER :: r
    INTEGER, DIMENSION(:,:), POINTER  :: atom_list,kind_label

!   ---------------------------------------------------------------------------

!   *** Deallocate the old data structure for the atomic information ***

    IF (natom > 0) THEN
      DEALLOCATE (atom_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"atom_info")
    END IF

    NULLIFY (atom_list,kind_label,kind_of,natom_of_kind,r)

    CALL start_input_session(globenv%input_file_name,globenv)

!   *** Search for the input section with the atomic coordinates ***

    CALL search(start_section,.TRUE.,found)

    IF (found) THEN

      nkind = 0
      natom = 0

      atom_list => reallocate(atom_list,1,nblock_kind,1,nblock_atom)
      kind_label => reallocate(kind_label,1,max_kind_label_length,&
                                          1,nblock_atom)
      kind_of => reallocate(kind_of,1,nblock_atom)
      natom_of_kind => reallocate(natom_of_kind,1,nblock_kind)
      r => reallocate(r,1,3,1,nblock_atom)

      DO

        CALL test_object(test_result,skip_lines=1)

        SELECT CASE (TRIM(test_result))
        CASE ("string")
          CALL read_object(current_label,skip_lines=1)
          string = current_label
          CALL uppercase(string)
          IF (string == end_section) EXIT
        CASE ("integer")
          CALL read_object(z,skip_lines=1)
          current_label = ptable(z)%symbol
        CASE DEFAULT
          CALL stop_program(routine,&
                            "Floating point number found for the atomic "//&
                            "label or number in section <"//&
                            TRIM(start_section)//">")
        END SELECT

        natom = natom + 1

        IF (natom > SIZE(atom_list,2)) THEN
          atom_list => reallocate(atom_list,1,nkind,1,natom+nblock_atom)
          kind_of => reallocate(kind_of,1,natom+nblock_atom)
          r => reallocate(r,1,3,1,natom+nblock_atom)
        END IF

!       *** Read the atomic coordinates ***

        CALL read_object(r(1,natom))
        CALL read_object(r(2,natom))
        CALL read_object(r(3,natom))

!       *** Check, if we found a new atomic kind ***

        label = current_label
        CALL uppercase(label)

        new_kind = .TRUE.

        DO ikind=1,nkind
          CALL integer_to_string(kind_label(:,ikind),string)
          CALL uppercase(string)
          IF (string == label) THEN
            new_kind = .FALSE.
            EXIT
          END IF
        END DO

        IF (new_kind) THEN
          nkind = nkind + 1
          IF (nkind > SIZE(atom_list,1)) THEN
            atom_list => reallocate(atom_list,1,nkind+nblock_kind,1,natom)
            kind_label => reallocate(kind_label,1,max_kind_label_length,&
                                                1,nkind+nblock_kind)
            natom_of_kind => reallocate(natom_of_kind,1,nkind+nblock_kind)
          END IF
          ikind = nkind
          CALL string_to_integer(current_label,kind_label(:,ikind))
        END IF

        kind_of(natom) = ikind
        natom_of_kind(ikind) = natom_of_kind(ikind) + 1
        atom_list(ikind,natom_of_kind(ikind)) = natom

      END DO

    ELSE

!     *** Stop program, if the end of file is reached ***

      CALL stop_program(routine,&
                        "No input section <"//TRIM(start_section)//&
                        "> found in the input file <"//&
                        TRIM(globenv%input_file_name)//">")

    END IF

    CALL finish_input_session

    atom_list => reallocate(atom_list,1,nkind,1,natom)
    kind_label => reallocate(kind_label,1,max_kind_label_length,1,nkind)
    kind_of => reallocate(kind_of,1,natom)
    natom_of_kind => reallocate(natom_of_kind,1,nkind)
    r => reallocate(r,1,3,1,natom)

!   *** Allocate the data structure for the atomic kind information ***

    CALL allocate_kind_info(nkind)

!   *** Initialize the data structure for atomic kind information ***

    CALL init_kind_info(kind_label,natom_of_kind,atom_list)

!   *** Allocate the data structure for the atomic information ***

    IF (natom > 0) THEN
      ALLOCATE (atom_info(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"atom_info(natom)",0)
    ELSE
      CALL stop_program(routine,&
                        "No atomic coordinates found in the input section <"//&
                        TRIM(start_section))

    END IF

!   *** Initialize the data structure for the atomic information ***

    DO iatom=1,natom
      ikind = kind_of(iatom)
      atom_info(iatom)%kind = ikind
      atom_info(iatom)%atom = iatom
      atom_info(iatom)%f(:) = 0.0_wp
      atom_info(iatom)%first_cgf = 0
      atom_info(iatom)%first_sgf = 0
      atom_info(iatom)%last_cgf = 0
      atom_info(iatom)%last_sgf = 0
      IF (scaled_coordinates) THEN
        atom_info(iatom)%r(:) = r(:,iatom)*abc(:)*unit_of_length
      ELSE
        atom_info(iatom)%r(:) = r(:,iatom)*unit_of_length
      END IF
      atom_info(iatom)%r_pbc(:) = pbc(atom_info(iatom)%r(:))
      atom_info(iatom)%v(:) = 0.0_wp
    END DO

!   *** Release work storage ***

    DEALLOCATE (atom_list,kind_label,kind_of,natom_of_kind,r)

  END SUBROUTINE read_atomic_coordinates

! *****************************************************************************

  SUBROUTINE write_atomic_coordinates(lunit)

!   Purpose: Write the atomic coordinates to the logical unit number lunit.

!   History: - Creation (05.06.2000, Matthias Krack)

!   ***************************************************************************

    USE cell_parameters, ONLY: unit_of_length,unit_of_length_name

    INTEGER, INTENT(IN) :: lunit

!   *** Local variables ***

    INTEGER :: iatom,ikind

!   ---------------------------------------------------------------------------

!   *** Write headlines ***

    WRITE (lunit,"(/,/,T2,A)")&
      "ATOMIC COORDINATES IN "//TRIM(unit_of_length_name)
    WRITE (lunit,"(/,T3,A,7X,2(A1,11X),A1,7X,A6,4X,A4,/)")&
      "Atom  Kind  Element","X","Y","Z","Z(eff)","Mass"

    DO iatom=1,natom
      ikind = atom_info(iatom)%kind
      WRITE (lunit,"(T2,I5,1X,I4,3X,A2,2X,I3,3F12.6,4X,I3,3X,F8.4)")&
        iatom,ikind,kind_info(ikind)%element_symbol,kind_info(ikind)%z,&
        atom_info(iatom)%r(:)/unit_of_length,kind_info(ikind)%zeff,&
        kind_info(ikind)%mass
    END DO

  END SUBROUTINE write_atomic_coordinates

! *****************************************************************************

  SUBROUTINE write_interatomic_distances(lunit)

!   Purpose: Write the matrix of the interatomic distances to the logical
!            output unit lunit.

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    USE cell_parameters,       ONLY: abc,unit_of_length,unit_of_length_name
    USE method_specifications, ONLY: periodic

    INTEGER, INTENT(IN) :: lunit

!   *** Local variables ***

    REAL(wp) :: dab,dab_image
    INTEGER  :: from,iatom,icell,icol,ikind,jatom,jcell,kcell,to

    REAL(wp), DIMENSION(3) :: ra,rab,rab_image,rb,rb_image

    REAL(wp), DIMENSION(natom,natom) :: distance_matrix

!   ---------------------------------------------------------------------------

    distance_matrix(:,:) = 0.0_wp

    DO iatom=1,natom

      ra(:) = atom_info(iatom)%r_pbc(:)

      DO jatom=iatom+1,natom

        rb(:) = atom_info(jatom)%r_pbc(:)

        rab(:) = rb(:) - ra(:)

        dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))

!       *** Search for the image of atom B which is closest to atom A ***

        IF (periodic) THEN
          DO icell=-1,1
            DO jcell=-1,1
              DO kcell=-1,1
                rb_image(:) = rb(:) + REAL((/icell,jcell,kcell/),wp)*abc(:)
                rab_image(:) = rb_image(:) - ra(:)
                dab_image = SQRT(rab_image(1)*rab_image(1) +&
                                 rab_image(2)*rab_image(2) +&
                                 rab_image(3)*rab_image(3))
                IF (dab_image < dab) dab = dab_image
              END DO
            END DO
          END DO
        END IF

        distance_matrix(iatom,jatom) = dab/unit_of_length
        distance_matrix(jatom,iatom) = distance_matrix(iatom,jatom)

      END DO

    END DO

!   *** Print the distance matrix ***

    WRITE (lunit,"(/,/,T2,A)")&
      "INTERATOMIC DISTANCES IN "//TRIM(unit_of_length_name)

    DO jatom=1,natom,5
      from = jatom
      to = MIN(from+4,natom)
      WRITE (lunit,"(/,T2,11X,5(4X,I5,4X))") (icol,icol=from,to)
      DO iatom=1,natom
        ikind = atom_info(iatom)%kind
        WRITE (lunit,"(T2,I5,2X,A2,2X,5(1X,F12.6))")&
          iatom,kind_info(ikind)%element_symbol,&
          (distance_matrix(iatom,icol),icol=from,to)
      END DO
    END DO

  END SUBROUTINE write_interatomic_distances

! *****************************************************************************

END MODULE atoms
