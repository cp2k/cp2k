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

  USE force_field_types,               ONLY: bend_kind_type,&
                                             bond_kind_type,&
                                             torsion_kind_type
  USE kinds,                           ONLY: int_size,&
                                             wp => dp
  USE molecule_types,                  ONLY: molecule_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "molecule_kind_types"
  INTEGER, PARAMETER          :: max_name_length = 60

! *** Define the derived structure types ***

  TYPE bond_type
    TYPE(bond_kind_type), POINTER :: kind
    INTEGER :: a,b
  END TYPE bond_type

  TYPE bend_type
    TYPE(bend_kind_type), POINTER :: kind
    INTEGER :: a,b,c
  END TYPE bend_type

  TYPE torsion_type
    TYPE(torsion_kind_type), POINTER :: kind
    INTEGER :: a,b,c,d
  END TYPE torsion_type

  TYPE molecule_kind_type
    PRIVATE
    TYPE(bend_type), DIMENSION(:), POINTER            :: bend_list
    TYPE(bond_type), DIMENSION(:), POINTER            :: bond_list
!MK TYPE(dist_constraint_type), DIMENSION(:), POINTER :: dist_list
!MK TYPE(g3x3_constraint_type), DIMENSION(:), POINTER :: g3x3_list
    TYPE(torsion_type), DIMENSION(:), POINTER         :: torsion_list
    CHARACTER(LEN=max_name_length)                    :: name,nhcopt
    REAL(wp)                                          :: charge,mass
    INTEGER                                           :: index,natom
    INTEGER, DIMENSION(:), POINTER                    :: molecule_list
  END TYPE molecule_kind_type

! *** Public subroutines ***

  PUBLIC :: allocate_molecule_kind_set,&
            deallocate_molecule_kind_set,&
            get_molecule_kind,&
            get_molecule_kind_set,&
            init_molecule_kind_set,&
            set_molecule_kind

! *** Public data types ***

  PUBLIC :: bend_type,&
            bond_type,&
            molecule_kind_type,&
            torsion_type

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_molecule_kind_set(molecule_kind_set,nkind)

!   Purpose: Allocate and initialize a molecule kind set.

!   History: - Creation (22.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), DIMENSION(:), POINTER :: molecule_kind_set
    INTEGER, INTENT(IN)                             :: nkind

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_molecule_kind_set"

!   *** Local variables ***

    INTEGER :: ikind,istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind_set)) THEN
      CALL deallocate_molecule_kind_set(molecule_kind_set)
    END IF

    ALLOCATE (molecule_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "molecule_kind_set",nkind*int_size)

    DO ikind=1,nkind
      NULLIFY (molecule_kind_set(ikind)%bond_list)
      NULLIFY (molecule_kind_set(ikind)%bend_list)
      NULLIFY (molecule_kind_set(ikind)%torsion_list)
      molecule_kind_set(ikind)%name = ""
      molecule_kind_set(ikind)%nhcopt = ""
      molecule_kind_set(ikind)%charge = 0.0_wp
      molecule_kind_set(ikind)%mass = 0.0_wp
      molecule_kind_set(ikind)%index = ikind
      molecule_kind_set(ikind)%natom = 0
      NULLIFY (molecule_kind_set(ikind)%molecule_list)
    END DO

  END SUBROUTINE allocate_molecule_kind_set

! *****************************************************************************

  SUBROUTINE deallocate_molecule_kind_set(molecule_kind_set)

!   Purpose: Deallocate a molecule kind set.

!   History: - Creation (22.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), DIMENSION(:), POINTER :: molecule_kind_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_molecule_kind_set"

!   *** Local variables ***

    INTEGER :: ikind,istat,nkind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind_set)) THEN

      nkind = SIZE(molecule_kind_set)

      DO ikind=1,nkind
        DEALLOCATE (molecule_kind_set(ikind)%molecule_list,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine_name,module_name,__LINE__,&
                           "molecule_kind_set(ikind)%molecule_list")
        END IF
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

  SUBROUTINE get_molecule_kind(molecule_kind,bond_list,bend_list,torsion_list,&
                               name,nhcopt,mass,index,natom,molecule_list,&
                               nmolecule)

!   Purpose: Get informations about a molecule kind.

!   History: - Creation (27.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), POINTER                     :: molecule_kind
    TYPE(bond_type), DIMENSION(:), OPTIONAL, POINTER      :: bond_list
    TYPE(bend_type), DIMENSION(:), OPTIONAL, POINTER      :: bend_list
    TYPE(torsion_type), DIMENSION(:), OPTIONAL, POINTER   :: torsion_list
    CHARACTER(LEN=max_name_length), OPTIONAL, INTENT(OUT) :: name,nhcopt
    REAL(wp), OPTIONAL, INTENT(OUT)                       :: mass
    INTEGER, OPTIONAL, INTENT(OUT)                        :: index,&
                                                             natom,&
                                                             nmolecule
    INTEGER, DIMENSION(:), OPTIONAL, POINTER              :: molecule_list

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_molecule_kind"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind)) THEN

      IF (PRESENT(bond_list)) bond_list => molecule_kind%bond_list
      IF (PRESENT(bend_list)) bend_list => molecule_kind%bend_list
      IF (PRESENT(torsion_list)) torsion_list => molecule_kind%torsion_list
      IF (PRESENT(name)) name = molecule_kind%name
      IF (PRESENT(nhcopt)) nhcopt = molecule_kind%nhcopt
      IF (PRESENT(mass)) mass = molecule_kind%mass
      IF (PRESENT(index)) index = molecule_kind%index
      IF (PRESENT(natom)) natom = molecule_kind%natom
      IF (PRESENT(molecule_list)) molecule_list => molecule_kind%molecule_list
      IF (PRESENT(nmolecule)) nmolecule = SIZE(molecule_kind%molecule_list)

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule_kind is not associated")

    END IF

  END SUBROUTINE get_molecule_kind

! *****************************************************************************

  SUBROUTINE get_molecule_kind_set(molecule_kind_set,maxatom,natom,nmolecule)

!   Purpose: Get informations about a molecule kind set.

!   History: - Creation (27.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), DIMENSION(:), POINTER :: molecule_kind_set
    INTEGER, OPTIONAL, INTENT(OUT)                  :: maxatom,&
                                                       natom,&
                                                       nmolecule

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_molecule_kind_set"

!   *** Local variables ***

    TYPE(molecule_kind_type), POINTER :: molecule_kind

    INTEGER :: ikind,na,nkind,nm

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind_set)) THEN

      IF (PRESENT(maxatom)) maxatom = 0
      IF (PRESENT(natom)) natom = 0
      IF (PRESENT(nmolecule)) nmolecule = 0

      nkind = SIZE(molecule_kind_set)

      DO ikind=1,nkind

        molecule_kind => molecule_kind_set(ikind)

        CALL get_molecule_kind(molecule_kind=molecule_kind,&
                               natom=na,&
                               nmolecule=nm)

        IF (PRESENT(maxatom)) THEN
          maxatom = MAX(maxatom,na)
        END IF

        IF (PRESENT(natom)) THEN
          natom = natom + na
        END IF

        IF (PRESENT(nmolecule)) THEN
          nmolecule = nmolecule + nm
        END IF

      END DO

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule_kind_set is not associated")

    END IF

  END SUBROUTINE get_molecule_kind_set

! *****************************************************************************

  SUBROUTINE init_molecule_kind_set(molecule_kind_set,mol_setup)

!   Purpose: Initialize a molecule kind set. Currently just an interface to
!            the old mol_setup structure type.

!   History: - Creation (27.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), DIMENSION(:), POINTER :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), INTENT(IN)   :: mol_setup

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_molecule_kind_set"

!   *** Local variables ***

    TYPE(molecule_kind_type), POINTER :: molecule_kind

    INTEGER :: i,imolecule,imolecule_kind,istat,n,nmolecule_kind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind_set)) THEN

      nmolecule_kind = SIZE(molecule_kind_set)

      imolecule = 0

      DO imolecule_kind=1,nmolecule_kind

        molecule_kind => molecule_kind_set(imolecule_kind)

        molecule_kind%index = imolecule_kind
        molecule_kind%charge = SUM(mol_setup(imolecule_kind)%molpar%acharge)
        molecule_kind%mass = SUM(mol_setup(imolecule_kind)%molpar%aweight)
        molecule_kind%name = TRIM(mol_setup(imolecule_kind)%molpar%name)
        molecule_kind%nhcopt = TRIM(mol_setup(imolecule_kind)%nhcopt)
        molecule_kind%natom = mol_setup(imolecule_kind)%molpar%natom

        n = mol_setup(imolecule_kind)%num_mol

        ALLOCATE (molecule_kind%molecule_list(n),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "molecule_kind%molecule_list",&
                                         n*int_size)
        DO i=1,n
          molecule_kind%molecule_list(i) = imolecule + i
        END DO

        imolecule = imolecule + n

        n = mol_setup(imolecule_kind)%molpar%nbonds
        n = mol_setup(imolecule_kind)%molpar%nbends
        n = mol_setup(imolecule_kind)%molpar%ntorsions

      END DO

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule_kind_set is not associated")

    END IF

  END SUBROUTINE init_molecule_kind_set

! *****************************************************************************

  SUBROUTINE set_molecule_kind(molecule_kind,name,mass,index,natom,&
                               molecule_list)

!   Purpose: Set the components of a molecule kind.

!   History: - Creation (27.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), POINTER                    :: molecule_kind
    CHARACTER(LEN=max_name_length), OPTIONAL, INTENT(IN) :: name
    REAL(wp), OPTIONAL, INTENT(IN)                       :: mass
    INTEGER, OPTIONAL, INTENT(IN)                        :: index,&
                                                            natom
    INTEGER, DIMENSION(:), OPTIONAL, POINTER             :: molecule_list

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_molecule_kind"

!   *** Local variables ***

    INTEGER :: istat,n

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind)) THEN

      IF (PRESENT(name)) molecule_kind%name = name
      IF (PRESENT(mass)) molecule_kind%mass = mass
      IF (PRESENT(index)) molecule_kind%index = index
      IF (PRESENT(natom)) molecule_kind%natom = natom
      IF (PRESENT(molecule_list)) THEN
        n = SIZE(molecule_list)
        IF (n > 0) THEN
          molecule_kind%natom = n
          ALLOCATE (molecule_kind%molecule_list(n),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "molecule_kind%molecule_list",&
                                           n*int_size)
          molecule_kind%molecule_list(:) = molecule_list(:)
        ELSE
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "An invalid molecule_list was supplied")
        END IF
      END IF

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule_kind is not associated")

    END IF

  END SUBROUTINE set_molecule_kind

! *****************************************************************************

END MODULE molecule_kind_types
