!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****s* cp2k/molecule_types_new [1.0] *
!!
!!   NAME
!!     molecule_types_new
!!
!!   FUNCTION
!!     Define the data structure for the molecule information.
!!
!!   AUTHOR
!!     MK (29.08.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE molecule_types_new

  USE kinds,                           ONLY: int_size
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "molecule_types_new"

! *** Data types ***

  TYPE local_molecule_type
    INTEGER, DIMENSION(:), POINTER :: states  ! indices of Kohn-Sham states for molecule
    INTEGER                        :: nstates ! Kohn-Sham states for molecule
  END TYPE local_molecule_type

  TYPE molecule_type
    TYPE(molecule_kind_type), POINTER  :: molecule_kind ! pointer to molecule kind information
    TYPE(local_molecule_type), POINTER :: lmi           ! local molecule information
    INTEGER                            :: first_atom,&  ! global index of the first atom in the molecule
                                          nvt_num       ! number of thermostat
  END TYPE molecule_type

! *** Public data types ***

  PUBLIC :: local_molecule_type,&
            molecule_type

! *** Public subroutines ***

  PUBLIC :: allocate_molecule_set,&
            deallocate_molecule_set,&
            get_molecule,&
            set_molecule,&
            set_molecule_set

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_molecule_set(molecule_set,nmolecule)

!   Purpose: Allocate a molecule set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_type), DIMENSION(:), POINTER :: molecule_set
    INTEGER, INTENT(IN)                        :: nmolecule

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_molecule_set"

!   *** Local variables ***

    INTEGER :: imolecule,istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_set)) CALL deallocate_molecule_set(molecule_set)

    ALLOCATE (molecule_set(nmolecule),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "molecule_set",nmolecule*int_size)

    DO imolecule=1,nmolecule
      NULLIFY (molecule_set(imolecule)%molecule_kind)
      NULLIFY (molecule_set(imolecule)%lmi)
      molecule_set(imolecule)%first_atom = 0
      molecule_set(imolecule)%nvt_num = 0
    END DO

  END SUBROUTINE allocate_molecule_set

! *****************************************************************************

  SUBROUTINE deallocate_molecule_set(molecule_set)

!   Purpose: Deallocate a molecule set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_type), DIMENSION(:), POINTER :: molecule_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_molecule_set"

!   *** Local variables ***

    INTEGER :: istat,imolecule

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_set)) THEN

      DO imolecule=1,SIZE(molecule_set)
         IF (ASSOCIATED(molecule_set(imolecule)%lmi)) THEN
             DEALLOCATE (molecule_set(imolecule)%lmi,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                              "molecule_set(i)%lmi")
         ENDIF
      ENDDO
      DEALLOCATE (molecule_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "molecule_set")

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule_set is not associated and "//&
                        "cannot be deallocated")

    END IF

  END SUBROUTINE deallocate_molecule_set

! *****************************************************************************

  SUBROUTINE get_molecule(molecule,molecule_kind,lmi)

!   Purpose: Get components from a molecule data set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_type), POINTER                :: molecule
    TYPE(molecule_kind_type), OPTIONAL, POINTER :: molecule_kind
    TYPE(local_molecule_type), OPTIONAL, POINTER :: lmi

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_molecule"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule)) THEN

      IF (PRESENT(molecule_kind)) molecule_kind => molecule%molecule_kind
      IF (PRESENT(lmi))           lmi           => molecule%lmi

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule is not associated")

    END IF

  END SUBROUTINE get_molecule

! *****************************************************************************

  SUBROUTINE set_molecule(molecule,molecule_kind,lmi)

!   Purpose: Set a molecule data set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_type), POINTER                :: molecule
    TYPE(molecule_kind_type), OPTIONAL, POINTER :: molecule_kind
    TYPE(local_molecule_type), OPTIONAL, POINTER :: lmi

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_molecule"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule)) THEN

      IF (PRESENT(molecule_kind)) molecule%molecule_kind => molecule_kind
      IF (PRESENT(lmi))           molecule%lmi           => lmi

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule is not associated")

    END IF

  END SUBROUTINE set_molecule

! *****************************************************************************

  SUBROUTINE set_molecule_set(molecule_set,first_atom)

!   Purpose: Set a molecule data set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_type), DIMENSION(:), POINTER  :: molecule_set
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: first_atom

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_molecule_set"

!   *** Local variables ***

    INTEGER :: imolecule,istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_set)) THEN

      IF (PRESENT(first_atom)) THEN

        IF (SIZE(first_atom) /= SIZE(molecule_set)) THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "The sizes of first_atom and molecule_set "//&
                            "are different")
        END IF

        DO imolecule=1,SIZE(molecule_set)
          molecule_set(imolecule)%first_atom = first_atom(imolecule)
        END DO
      END IF

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule_set is not associated")

    END IF

  END SUBROUTINE set_molecule_set

! *****************************************************************************

END MODULE molecule_types_new
