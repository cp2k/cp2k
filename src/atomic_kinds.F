!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/atomic_kinds [1.0] *
!!
!!   NAME
!!     atomic_kinds
!!
!!   FUNCTION
!!     Define the data structure for the atomic kind information.
!!
!!   AUTHOR
!!     Matthias Krack (28.06.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE atomic_kinds

! *****************************************************************************

! Index:

! SUBROUTINE allocate_kind_info(number_of_kinds)
! SUBROUTINE deallocate_kind_info(number_of_kinds)
! SUBROUTINE init_kind_info(kind_label,natom_of_kind,atom_list)
! SUBROUTINE read_kind_info(start_section,end_section)
! SUBROUTINE write_kind_info(lunit)

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE basis_set_types,          ONLY: allocate_basis_set,&
                                      deallocate_basis_set,&
                                      gto_basis_set_type,&
                                      maxlcgf,&
                                      read_basis_set
  USE global_types,             ONLY: global_environment_type
  USE method_specifications,    ONLY: allchem
  USE external_potential_types, ONLY: all_potential_type,&
                                      allocate_potential,&
                                      deallocate_potential,&
                                      gth_potential_type,&
                                      maxlppl,&
                                      maxlppnl,&
                                      read_potential

  IMPLICIT NONE

  PRIVATE

! *** Parameters ***

  INTEGER, PARAMETER :: max_kind_label_length = 8,&
                        max_basis_set_name_length = 60

! *** Define the atomic kind info type ***

  TYPE kind_info_type
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    CHARACTER(LEN=max_kind_label_length)     :: label
    CHARACTER(LEN=max_basis_set_name_length) :: aux_basis_set_name,&
                                                orb_basis_set_name,&
                                                potential_name
    CHARACTER(LEN=2)                         :: element_symbol
    REAL(wp)                                 :: alpha_core_charge,&
                                                ccore_charge,&
                                                core_charge_radius,mass,&
                                                prj_paw_radius
    INTEGER                                  :: natom,z,zeff
    INTEGER                                  :: number_of_grid_points
    INTEGER, DIMENSION(:), POINTER           :: atom_list,elec_conf
  END TYPE kind_info_type

  TYPE(kind_info_type), DIMENSION(:), POINTER :: kind_info

  INTEGER :: nkind = 0
  LOGICAL :: all_potential_present = .FALSE.,&
             gth_potential_present = .FALSE.

! *** Public subroutines ***

  PUBLIC :: allocate_kind_info,&
            init_kind_info,&
            read_kind_info,&
            write_kind_info

! *** Public data types ***

  PUBLIC :: kind_info_type

! *** Public parameters ***

  PUBLIC :: max_basis_set_name_length,&
            max_kind_label_length

! *** Public variables ***

  PUBLIC :: all_potential_present,&
            gth_potential_present,&
            kind_info,&
            nkind

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_kind_info(number_of_kinds)

!   Purpose: Allocate the data structure for the atomic kind information.

!   History: - Creation (26.10.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_memory,stop_program

    INTEGER, INTENT(IN) :: number_of_kinds

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_kind_info (MODULE atomic_kinds)"

!   *** Local variables ***

    INTEGER :: ikind,istat

!   ---------------------------------------------------------------------------

!   *** Deallocate the old data structure for the atomic kind information ***

    IF (nkind > 0) THEN
      DO ikind=1,nkind
        IF (ASSOCIATED(kind_info(ikind)%all_potential)) THEN
          CALL deallocate_potential(kind_info(ikind)%all_potential)
        END IF
        IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN
          CALL deallocate_potential(kind_info(ikind)%gth_potential)
        END IF
        IF (ASSOCIATED(kind_info(ikind)%aux_basis_set)) THEN
          CALL deallocate_basis_set(kind_info(ikind)%aux_basis_set)
        END IF
        IF (ASSOCIATED(kind_info(ikind)%orb_basis_set)) THEN
          CALL deallocate_basis_set(kind_info(ikind)%orb_basis_set)
        END IF
        DEALLOCATE (kind_info(ikind)%atom_list,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"kind_info(ikind)%atom_list")
        DEALLOCATE (kind_info(ikind)%elec_conf,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"kind_info(ikind)%elec_conf")
      END DO
      DEALLOCATE (kind_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"kind_info")
    END IF

    nkind = number_of_kinds

    IF (nkind > 0) THEN
      ALLOCATE (kind_info(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"kind_info(nkind)",0)
      DO ikind=1,nkind
        NULLIFY (kind_info(ikind)%all_potential)
        NULLIFY (kind_info(ikind)%gth_potential)
        NULLIFY (kind_info(ikind)%aux_basis_set)
        NULLIFY (kind_info(ikind)%orb_basis_set)
        NULLIFY (kind_info(ikind)%atom_list)
        NULLIFY (kind_info(ikind)%elec_conf)
      END DO
    ELSE
      CALL stop_program(routine,&
                        "An invalid number of atomic kinds was specified")
    END IF

  END SUBROUTINE allocate_kind_info

! *****************************************************************************

  SUBROUTINE init_kind_info(kind_label,natom_of_kind,atom_list)

!   Purpose: Initialize the atomic kind information.

!   History: - Creation (26.10.2000, Matthias Krack)

!   ***************************************************************************

    USE memory_utilities, ONLY: reallocate
    USE string_utilities, ONLY: integer_to_string
    USE termination,      ONLY: stop_program

    INTEGER, DIMENSION(:), INTENT(IN)   :: natom_of_kind
    INTEGER, DIMENSION(:,:), INTENT(IN) :: atom_list,kind_label

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE init_kind_info (MODULE atomic_kinds)"

!   *** Local variables ***

    INTEGER :: ikind

!   ---------------------------------------------------------------------------

!   *** Check the argument dimensions ***

    IF (SIZE(kind_label,1) /= max_kind_label_length) THEN
      CALL stop_program(routine,"The first dimension of the argument "//&
                                "kind_label has an invalid size")
    END IF

    IF (SIZE(kind_label,2) /= nkind) THEN
      CALL stop_program(routine,"The second dimension of the argument "//&
                                "kind_label has an invalid size")
    END IF

    IF (SIZE(natom_of_kind) /= nkind) THEN
      CALL stop_program(routine,"The first dimension of the argument "//&
                                "natom_of_kind has an invalid size")
    END IF

    IF (SIZE(atom_list,1) /= nkind) THEN
      CALL stop_program(routine,"The first dimension of the argument "//&
                                "atom_list has an invalid size")
    END IF

    DO ikind=1,nkind
      CALL integer_to_string(kind_label(:,ikind),kind_info(ikind)%label)
      kind_info(ikind)%natom = natom_of_kind(ikind)
      kind_info(ikind)%atom_list =>&
        reallocate(kind_info(ikind)%atom_list,1,natom_of_kind(ikind))
      kind_info(ikind)%atom_list(1:natom_of_kind(ikind)) =&
        atom_list(ikind,1:natom_of_kind(ikind))
    END DO

  END SUBROUTINE init_kind_info

! *****************************************************************************

  SUBROUTINE read_kind_info(start_section,end_section,globenv)

!   Purpose: Read the kind data from the input file.

!   History: - Creation (05.06.2000, Matthias Krack)

!   ***************************************************************************

    USE input_utilities,  ONLY: finish_input_session,&
                                read_object,&
                                search,&
                                start_input_session
    USE mathconstants,    ONLY: pi
    USE periodic_table,   ONLY: ptable
    USE string_utilities, ONLY: uppercase
    USE termination,      ONLY: stop_program

    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=*), INTENT(IN)              :: start_section,end_section

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_kind_info (MODULE atomic_kinds)"

!   *** Local variables ***

    CHARACTER(LEN=max_basis_set_name_length) :: apname
    CHARACTER(LEN=max_kind_label_length)     :: current_label,requested_label
    CHARACTER(LEN=40)                        :: keyword
    CHARACTER(LEN=2)                         :: current_symbol,element_symbol,&
                                                symbol
    INTEGER                                  :: ikind,ipos,istat,z
    LOGICAL                                  :: found

!   ---------------------------------------------------------------------------

    all_potential_present = .FALSE.
    gth_potential_present = .FALSE.

    maxlcgf = -1
    maxlppl = -1
    maxlppnl = -1

    kind_loop: DO ikind=1,nkind

      CALL start_input_session(globenv%input_file_name,globenv)

      requested_label = kind_info(ikind)%label
      CALL uppercase(requested_label)

!     *** Search for the input sections with ***
!     *** the corresponding kind definition  ***

      search_loop: DO

      CALL search(start_section,.TRUE.,found)

        IF (found) THEN

          CALL read_object(current_label)
          CALL uppercase(current_label)

          IF (current_label == requested_label) THEN

            kind_info(ikind)%aux_basis_set_name = ""
            kind_info(ikind)%orb_basis_set_name = ""
            kind_info(ikind)%potential_name = ""
            kind_info(ikind)%z = -1

            DO
              CALL read_object(keyword,skip_lines=1)
              CALL uppercase(keyword)
              SELECT CASE (TRIM(keyword))
              CASE ("AUXILIARY_BASIS_SET","AUX_BASIS_SET","AUX_BASIS")
                CALL read_object(kind_info(ikind)%aux_basis_set_name)
              CASE ("ORBITAL_BASIS_SET","ORB_BASIS_SET","ORB_BASIS",&
                    "BASIS_SET","BASIS")
                CALL read_object(kind_info(ikind)%orb_basis_set_name)
              CASE ("ELEMENT_SYMBOL","ELEMENT")
                CALL read_object(element_symbol)
                symbol = element_symbol
                CALL uppercase(symbol)
                DO z=LBOUND(ptable,1),UBOUND(ptable,1)
                  current_symbol = ptable(z)%symbol
                  CALL uppercase(current_symbol)
                  IF (current_symbol == symbol) THEN
                    kind_info(ikind)%element_symbol = ptable(z)%symbol
                    kind_info(ikind)%mass = ptable(z)%mass
                    kind_info(ikind)%z = z
                    EXIT
                  END IF
                END DO
              CASE ("MASS","WEIGHT")
                CALL read_object(kind_info(ikind)%mass)
              CASE ("POTENTIAL","POT")
                CALL read_object(kind_info(ikind)%potential_name)
              CASE DEFAULT
                IF (keyword == end_section) THEN
                  EXIT
                ELSE
                  CALL stop_program(routine,&
                                    "Invalid keyword <"//TRIM(keyword)//&
                                    "> found in the input section <"//&
                                    TRIM(start_section)//">")
                END IF
              END SELECT
            END DO

!           *** The corresponding element of the current atomic    ***
!           *** kind is still undefined. Try the kind label itself ***

            IF (kind_info(ikind)%z == -1) THEN
              DO z=LBOUND(ptable,1),UBOUND(ptable,1)
                current_symbol = ptable(z)%symbol
                CALL uppercase(current_symbol)
                IF (current_symbol == current_label) THEN
                  kind_info(ikind)%element_symbol = ptable(z)%symbol
                  kind_info(ikind)%mass = ptable(z)%mass
                  kind_info(ikind)%z = z
                  EXIT
                END IF
              END DO
            END IF

            IF (kind_info(ikind)%orb_basis_set_name == "") THEN
              CALL stop_program(routine,&
                                "No orbital basis set was defined for "//&
                                "the atomic kind <"//TRIM(requested_label)//&
                                ">. Specify a valid orbital basis set "//&
                                "name in the input section <"//&
                                TRIM(start_section)//" "//&
                                TRIM(requested_label)//">")
            END IF

            IF (allchem.AND.(kind_info(ikind)%aux_basis_set_name == "")) THEN
              CALL stop_program(routine,&
                                "No auxiliary basis set was defined for "//&
                                "the atomic kind <"//TRIM(requested_label)//&
                                ">. Specify a valid auxiliary basis set "//&
                                "name in the input section <"//&
                                TRIM(start_section)//" "//&
                                TRIM(requested_label)//">")
            END IF

            IF (kind_info(ikind)%potential_name == "") THEN
              CALL stop_program(routine,&
                                "No atomic potential type was defined "//&
                                "for the atomic kind <"//&
                                TRIM(requested_label)//&
                                ">. Specify a valid atomic potential "//&
                                "type in the input section <"//&
                                TRIM(start_section)//" "//&
                                TRIM(requested_label)//">")
            END IF

            IF (kind_info(ikind)%z == -1) THEN
              CALL stop_program(routine,&
                                "No element could be assign to the "//&
                                "atomic kind <"//TRIM(requested_label)//&
                                ">. Specify a valid element symbol "//&
                                "in the input section <"//&
                                TRIM(start_section)//" "//&
                                TRIM(requested_label)//"> or select an "//&
                                "element symbol as the kind label")
            END IF

            EXIT search_loop

          END IF

        ELSE

          CALL stop_program(routine,&
                            "The input section <"//TRIM(start_section)//&
                            "> for the atomic kind <"//&
                            TRIM(requested_label)//&
                            "> was not found in the input file <"//&
                            TRIM(globenv%input_file_name)//">")
        END IF

      END DO search_loop

      CALL finish_input_session

!     *** Allocate the orbital basis set data structure ***

      CALL allocate_basis_set(kind_info(ikind)%orb_basis_set)

!     *** Read the basis set database file ***

      CALL read_basis_set(kind_info(ikind)%element_symbol,&
                          kind_info(ikind)%orb_basis_set_name,&
                          kind_info(ikind)%orb_basis_set,&
                          globenv)

!     *** Allocate and initialize the auxiliary basis set data structure ***

      IF (allchem) THEN
        CALL allocate_basis_set(kind_info(ikind)%aux_basis_set)
        CALL read_basis_set(kind_info(ikind)%element_symbol,&
                            kind_info(ikind)%aux_basis_set_name,&
                            kind_info(ikind)%aux_basis_set,&
                            globenv)
      END IF

!     *** Read the potential database file ***

      ipos = INDEX(kind_info(ikind)%potential_name,"-")

      IF (ipos > 0) THEN
        apname(1:ipos-1) = kind_info(ikind)%potential_name(1:ipos-1)
      ELSE
        apname = kind_info(ikind)%potential_name
      END IF

      CALL uppercase(apname)

      IF (apname == "ALL") THEN
        CALL allocate_potential(kind_info(ikind)%all_potential)
        CALL read_potential(kind_info(ikind)%element_symbol,apname,&
                            kind_info(ikind)%alpha_core_charge,&
                            kind_info(ikind)%elec_conf,&
                            kind_info(ikind)%zeff,&
                            kind_info(ikind)%all_potential,&
                            globenv)
        all_potential_present = .TRUE.
      ELSE IF (apname == "GTH") THEN
        CALL allocate_potential(kind_info(ikind)%gth_potential)
        CALL read_potential(kind_info(ikind)%element_symbol,apname,&
                            kind_info(ikind)%alpha_core_charge,&
                            kind_info(ikind)%elec_conf,&
                            kind_info(ikind)%zeff,&
                            kind_info(ikind)%gth_potential,&
                            globenv)
        gth_potential_present = .TRUE.
      ELSE
        CALL stop_program(routine,&
                          "An invalid potential type <"//&
                          TRIM(apname)//"> was specified "//&
                          "for the atomic kind <"//&
                          TRIM(kind_info(ikind)%label)//&
                          "> in the input file <"//&
                          TRIM(globenv%input_file_name)//">")
      END IF

      kind_info(ikind)%ccore_charge =&
        REAL(kind_info(ikind)%zeff,wp)*&
        SQRT((kind_info(ikind)%alpha_core_charge/pi)**3)

!     *** Check the number of electrons ***

      z = SUM(kind_info(ikind)%elec_conf)

      IF (z /= kind_info(ikind)%z) THEN
        CALL stop_program(routine,&
                          "The total number of electrons defined by "//&
                          "the electronic configuration of the potential "//&
                          "data set for the atomic kind <"//&
                          TRIM(kind_info(ikind)%label)//&
                          "> in the input file <"//&
                          TRIM(globenv%input_file_name)//&
                          "> does not fit to its atomic number")
      END IF

    END DO kind_loop

  END SUBROUTINE read_kind_info

! *****************************************************************************

  SUBROUTINE write_kind_info(globenv)

!   Purpose: Write the kind information to the logical unit number lunit.

!   History: - Creation (06.06.2000, Matthias Krack)

!   ***************************************************************************

    USE orbital_pointers, ONLY: nco

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    REAL(wp) :: r
    INTEGER  :: i,ico,icgf,ikind,ipgf,iset,ishell,j,l,output_unit

!   ---------------------------------------------------------------------------

    output_unit = globenv%scr

    IF (globenv%print%basis_set.OR.globenv%print%potential) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") "ATOMIC KIND INFORMATION"
    END IF

    DO ikind=1,nkind

      IF (globenv%print%basis_set.OR.globenv%print%potential) THEN
        WRITE (UNIT=output_unit,FMT="(/,T2,I2,A,/,/,T6,A,I6)")&
          ikind,". Atomic kind: "//TRIM(kind_info(ikind)%label),&
          "Number of atoms of this kind in the cell: ",kind_info(ikind)%natom
      END IF

      IF (globenv%print%basis_set) THEN

        WRITE (UNIT=output_unit,FMT="(/,T3,A)")&
          "a) GTO basis set information"

        WRITE (UNIT=output_unit,FMT="(/,T6,A)")&
          "Orbital basis set name: "//TRIM(kind_info(ikind)%orb_basis_set_name)

        WRITE (UNIT=output_unit,FMT="(/,(T6,A,I6))")&
          "Number of orbital shell sets:            ",&
          kind_info(ikind)%orb_basis_set%nset,&
          "Number of orbital shells:                ",&
          SUM(kind_info(ikind)%orb_basis_set%nshell(:)),&
          "Number of primitive Cartesian functions: ",&
          SUM(kind_info(ikind)%orb_basis_set%npgf(:)),&
          "Number of Cartesian basis functions:     ",&
          kind_info(ikind)%orb_basis_set%ncgf,&
          "Number of spherical basis functions:     ",&
          kind_info(ikind)%orb_basis_set%nsgf

        WRITE (UNIT=output_unit,FMT="(/,T6,A,/,/,T6,A)")&
          "Normalized Cartesian orbitals",&
          "Set   Shell   Orbital            Exponent    Coefficient"

        icgf = 0

        DO iset=1,kind_info(ikind)%orb_basis_set%nset
          DO ishell=1,kind_info(ikind)%orb_basis_set%nshell(iset)
            WRITE (UNIT=output_unit,FMT="(A)") ""
            DO ico=1,nco(kind_info(ikind)%orb_basis_set%l(ishell,iset))
              icgf = icgf + 1
              WRITE (UNIT=output_unit,FMT="(T6,I3,4X,I4,3X,A12,(T32,2F15.6))")&
                iset,ishell,&
                kind_info(ikind)%orb_basis_set%cgf_symbol(icgf),&
                (kind_info(ikind)%orb_basis_set%zet(ipgf,iset),&
                 kind_info(ikind)%orb_basis_set%norm_cgf(icgf)*&
                 kind_info(ikind)%orb_basis_set%gcc(ipgf,ishell,iset),&
                 ipgf=1,kind_info(ikind)%orb_basis_set%npgf(iset))
            END DO
          END DO
        END DO

        IF (ASSOCIATED(kind_info(ikind)%aux_basis_set)) THEN

          WRITE (UNIT=output_unit,FMT="(/,T6,A)")&
            "Auxiliary basis set name: "//&
            TRIM(kind_info(ikind)%aux_basis_set_name)

          WRITE (UNIT=output_unit,FMT="(/,(T6,A,I6))")&
            "Number of orbital shell sets:            ",&
            kind_info(ikind)%aux_basis_set%nset,&
            "Number of orbital shells:                ",&
            SUM(kind_info(ikind)%aux_basis_set%nshell(:)),&
            "Number of primitive Cartesian functions: ",&
            SUM(kind_info(ikind)%aux_basis_set%npgf(:)),&
            "Number of Cartesian basis functions:     ",&
            kind_info(ikind)%aux_basis_set%ncgf,&
            "Number of spherical basis functions:     ",&
            kind_info(ikind)%aux_basis_set%nsgf

          WRITE (UNIT=output_unit,FMT="(/,T6,A,/,/,T6,A)")&
            "Normalized Cartesian auxiliary functions",&
            "Set   Shell   Function           Exponent    Coefficient"

          icgf = 0

          DO iset=1,kind_info(ikind)%aux_basis_set%nset
            DO ishell=1,kind_info(ikind)%aux_basis_set%nshell(iset)
              WRITE (UNIT=output_unit,FMT="(A)") ""
              DO ico=1,nco(kind_info(ikind)%aux_basis_set%l(ishell,iset))
                icgf = icgf + 1
                WRITE (UNIT=output_unit,&
                       FMT="(T6,I3,4X,I4,3X,A12,(T32,2F15.6))")&
                  iset,ishell,&
                  kind_info(ikind)%aux_basis_set%cgf_symbol(icgf),&
                  (kind_info(ikind)%aux_basis_set%zet(ipgf,iset),&
                   kind_info(ikind)%aux_basis_set%norm_cgf(icgf)*&
                   kind_info(ikind)%aux_basis_set%gcc(ipgf,ishell,iset),&
                   ipgf=1,kind_info(ikind)%aux_basis_set%npgf(iset))
              END DO
            END DO
          END DO

        END IF

      END IF

      IF (globenv%print%potential) THEN

        WRITE (UNIT=output_unit,FMT="(/,/,T3,A,/,/,T6,A)")&
          "b) Potential information",&
          "Potential type: "//TRIM(kind_info(ikind)%potential_name)

        IF (ASSOCIATED(kind_info(ikind)%all_potential)) THEN
          WRITE (UNIT=output_unit,FMT="(/,(T6,A))")&
            "Description: "//&
            TRIM(kind_info(ikind)%all_potential%description(1)),&
            "             "//&
            TRIM(kind_info(ikind)%all_potential%description(2))
          WRITE (UNIT=output_unit,FMT="(/,T6,A,F12.6)")&
            "Gaussian exponent of the core charge distribution: ",&
            kind_info(ikind)%alpha_core_charge
          WRITE (UNIT=output_unit,FMT="(/,T6,A,10I3)")&
            "Electronic configuration (s p d ...):",&
            kind_info(ikind)%elec_conf
        ELSE IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN
          WRITE (UNIT=output_unit,FMT="(/,(T6,A))")&
            "Description: "//&
            TRIM(kind_info(ikind)%gth_potential%description(1)),&
            "             "//&
            TRIM(kind_info(ikind)%gth_potential%description(2)),&
            "             "//&
            TRIM(kind_info(ikind)%gth_potential%description(3))
          WRITE (UNIT=output_unit,FMT="(/,T6,A,F12.6)")&
            "Gaussian exponent of the core charge distribution: ",&
            kind_info(ikind)%alpha_core_charge
          WRITE (UNIT=output_unit,FMT="(/,T6,A,10I3)")&
            "Electronic configuration (s p d ...):",&
            kind_info(ikind)%gth_potential%elec_conf
          r = 1.0_wp/SQRT(2.0_wp*kind_info(ikind)%gth_potential%alpha_ppl)
          WRITE (UNIT=output_unit,FMT="(/,T6,A,/,/,T8,A,/,T8,4F12.6)")&
            "Parameters of the local part of the GTH pseudopotential:",&
            "      rloc        C1          C2          C3          C4",&
            r,(kind_info(ikind)%gth_potential%cexp_ppl(i)*r**(2*(i-1)),&
               i=1,kind_info(ikind)%gth_potential%nexp_ppl)
          IF (kind_info(ikind)%gth_potential%lppnl > -1) THEN
            WRITE (UNIT=output_unit,FMT="(/,T6,A,/,/,T6,A,/)")&
              "Parameters of the non-local part of the GTH pseudopotential:",&
              " l      r(l)      h(i,j,l)"
            DO l=0,kind_info(ikind)%gth_potential%lppnl
              r = SQRT(0.5_wp/kind_info(ikind)%gth_potential%alpha_ppnl(l))
              WRITE (UNIT=output_unit,FMT="(T6,I2,F12.6,4F12.6)")&
                l,r,(kind_info(ikind)%gth_potential%hprj_ppnl(1,j,l),&
                     j=1,kind_info(ikind)%gth_potential%nprj_ppnl(l))
              DO i=2,kind_info(ikind)%gth_potential%nprj_ppnl(l)
                WRITE (UNIT=output_unit,FMT="(T20,4F12.6)")&
                  (kind_info(ikind)%gth_potential%hprj_ppnl(i,j,l),&
                   j=1,kind_info(ikind)%gth_potential%nprj_ppnl(l))
              END DO
            END DO
          END IF
        END IF

      END IF

    END DO

  END SUBROUTINE write_kind_info

! *****************************************************************************

END MODULE atomic_kinds
