!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/readin [1.0] *
!!
!!   NAME
!!     readin
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     Matthias Krack (18.05.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE readin

! *****************************************************************************

! Index:

! SUBROUTINE read_files(globenv)
! SUBROUTINE write_total_numbers(lunit)

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE global_types, ONLY: global_environment_type

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_files

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE read_files(globenv)

!   Purpose: Read the input and the database files.

!   History: - Creation (22.05.2000, Matthias Krack)

!   ***************************************************************************

    USE atomic_kinds,                    ONLY: all_potential_present,&
                                               kind_info,&
                                               nkind,&
                                               read_kind_info,&
                                               write_kind_info
    USE atoms,                           ONLY: init_atom_info,&
                                               read_atomic_coordinates,&
                                               write_atomic_coordinates,&
                                               write_interatomic_distances
    USE basis_set_types,                 ONLY: init_aux_basis_set,&
                                               init_orb_basis_set,&
                                               maxlcgf,&
                                               write_basis_set
    USE cell_parameters,                 ONLY: read_cell_parameters,&
                                               write_cell_parameters
    USE external_potential_types,        ONLY: init_potential,&
                                               maxlppl,&
                                               maxlppnl
    USE gamma,                           ONLY: init_md_ftable
    USE input_utilities,                 ONLY: get_special_input_character
    USE interactions,                    ONLY: init_interaction_radii,&
                                               read_interaction_thresholds,&
                                               write_interaction_thresholds
    USE method_specifications,           ONLY: maxder,&
                                               read_method_specifications,&
                                               write_method_specifications
    USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
    USE qs_scf,                          ONLY: read_scf_parameters,&
                                               write_scf_parameters
    USE timings,                         ONLY: timeset,timestop

    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local variables ***

    CHARACTER(LEN=4) :: end_section
    CHARACTER(LEN=1) :: section_character
    INTEGER          :: handle,ikind,maxl,output_unit
    LOGICAL          :: ionode

!   ---------------------------------------------------------------------------

    CALL timeset("read_files","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL get_special_input_character(section=section_character)

    end_section = section_character//"END"

!   *** Read the input section with the method specifications ***

    CALL read_method_specifications(section_character//"METHOD",end_section,&
                                    globenv)

!   *** Print the method specifications ***

    IF (ionode.AND.globenv%print%method_specifications) THEN
      CALL write_method_specifications(output_unit)
    END IF

!   *** Read the input section with the periodic cell parameters ***

    CALL read_cell_parameters(section_character//"CELL",end_section,globenv)

!   *** Print the cell parameters ***

    IF (ionode.AND.globenv%print%cell_parameters) THEN
      CALL write_cell_parameters(output_unit)
    END IF

!   *** Read the input section with the atomic coordinates ***

    CALL read_atomic_coordinates(section_character//"CART",end_section,globenv)

!   *** Read the kind information ***

    CALL read_kind_info(section_character//"KIND",end_section,globenv)

!   *** Print the unnormalized basis set information (input data) ***

    IF (ionode.AND.globenv%print%basis_set) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
        "BASIS SET INFORMATION (Unnormalized Gaussian-type functions)"
      DO ikind=1,nkind
        WRITE (UNIT=output_unit,FMT="(/,T2,I2,A,/,/,T6,A,/,/,T6,A)")&
          ikind,". Atomic kind: "//TRIM(kind_info(ikind)%label),&
          "Orbital basis set name: "//&
          TRIM(kind_info(ikind)%orb_basis_set_name),&
          "Set   Shell     n   l            Exponent    Coefficient"
        CALL write_basis_set(kind_info(ikind)%orb_basis_set,output_unit)
        IF (ASSOCIATED(kind_info(ikind)%aux_basis_set)) THEN
          WRITE (UNIT=output_unit,FMT="(/,T6,A)")&
            "Auxiliary basis set name: "//&
            TRIM(kind_info(ikind)%aux_basis_set_name),&
            "Set   Shell     n   l            Exponent    Coefficient"
          CALL write_basis_set(kind_info(ikind)%aux_basis_set,output_unit)
        END IF
      END DO
    END IF

!   *** Initialize the spherical harmonics and ***
!   *** the orbital transformation matrices    ***

    maxl = MAX(maxlcgf+maxder,maxlppl,maxlppnl)

    CALL init_spherical_harmonics(maxl,globenv)

!   *** Initialize the pretabulation for the calculation of the   ***
!   *** incomplete Gamma function F_n(t) after McMurchie-Davidson ***

    maxl = 3*maxlcgf + 1

    CALL init_md_ftable(maxl)

!   *** Initialize the basis set and atomic potential data sets ***

    DO ikind=1,nkind
      CALL init_orb_basis_set(kind_info(ikind)%orb_basis_set)
      IF (ASSOCIATED(kind_info(ikind)%aux_basis_set)) THEN
        CALL init_aux_basis_set(kind_info(ikind)%aux_basis_set)
      END IF
      IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN
        CALL init_potential(kind_info(ikind)%gth_potential)
      END IF
    END DO

!   *** Print the kind information ***

    IF (ionode.AND.globenv%print%basis_set) THEN
      CALL write_kind_info(globenv)
    END IF

!   *** Initialize the atomic data sets completely ***

    CALL init_atom_info

!   *** Print the total number of kinds, atoms, basis functions etc. ***

    IF (ionode.AND.globenv%print%total_numbers) THEN
      CALL write_total_numbers(output_unit)
    END IF

!   *** Print the atomic coordinates ***

    IF (ionode.AND.globenv%print%atomic_coordinates) THEN
      CALL write_atomic_coordinates(output_unit)
    END IF

!   *** Print the interatomic distances ***

    IF (ionode.AND.globenv%print%interatomic_distances) THEN
      CALL write_interatomic_distances(output_unit)
    END IF

!   *** Read the interaction threshold values ***

    CALL read_interaction_thresholds(section_character//"INTER",end_section,&
                                     globenv)

!   *** Print the interaction threshold values ***

    IF (ionode.AND.globenv%print%interaction_thresholds) THEN
      CALL write_interaction_thresholds(output_unit)
    END IF

!   *** Initialize the atomic interaction radii ***

    CALL init_interaction_radii(globenv)

!   *** Read SCF parameters ***

    CALL read_scf_parameters(section_character//"SCF",end_section,globenv)

    IF (ionode.AND.globenv%print%scf) THEN
      CALL write_scf_parameters(output_unit)
    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE read_files

! *****************************************************************************

  SUBROUTINE write_total_numbers(lunit)

!   Purpose: Write the total number of kinds, atoms, etc. to the logical unit
!            number lunit.

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    USE atomic_kinds,             ONLY: gth_potential_present,kind_info,nkind
    USE atoms,                    ONLY: atom_info,natom,ncgf,nsgf
    USE basis_set_types,          ONLY: maxlcgf
    USE external_potential_types, ONLY: maxlppl,maxlppnl

    INTEGER, INTENT(IN) :: lunit

!   *** Local variables ***

    INTEGER :: iatom,ikind,npgf,nset,nshell

!   ---------------------------------------------------------------------------

    npgf = 0
    nset = 0
    nshell = 0

    DO iatom=1,natom
      ikind = atom_info(iatom)%kind
      npgf = npgf + SUM(kind_info(ikind)%orb_basis_set%npgf(:))
      nset = nset + kind_info(ikind)%orb_basis_set%nset
      nshell = nshell + SUM(kind_info(ikind)%orb_basis_set%nshell(:))
    END DO

    WRITE (UNIT=lunit,FMT="(/,/,T2,A)") "TOTAL NUMBERS AND MAXIMUM NUMBERS"

    WRITE (UNIT=lunit,FMT="(/,T3,A,/,/,(T3,A,I6))")&
      "Total number of",&
      "- Atomic kinds:                  ",nkind,&
      "- Atoms:                         ",natom,&
      "- Shell sets:                    ",nset,&
      "- Shells:                        ",nshell,&
      "- Primitive Cartesian functions: ",npgf,&
      "- Cartesian basis functions:     ",ncgf,&
      "- Spherical basis functions:     ",nsgf

    IF (gth_potential_present) THEN
      WRITE (UNIT=lunit,FMT="(/,T3,A,/,/,(T3,A,I6))")&
        "Maximum angular momentum quantum number of the",&
        "- Orbital basis functions:                   ",maxlcgf,&
        "- Local part of the GTH pseudopotential:     ",maxlppl,&
        "- Non-local part of the GTH pseudopotential: ",maxlppnl
    ELSE
      WRITE (UNIT=lunit,FMT="(/,T3,A,/,/,(T3,A,I6))")&
        "Maximum angular momentum quantum number of the",&
        "- Orbital basis functions: ",maxlcgf
    END IF

  END SUBROUTINE write_total_numbers

! *****************************************************************************

END MODULE readin
