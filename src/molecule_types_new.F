!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
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
!!     JGH (22.05.2004) add last_atom information
!!
!!   SOURCE
!******************************************************************************

MODULE molecule_types_new

  USE kinds,                           ONLY: dp,int_size
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "molecule_types_new"

! *** Data types ***
  TYPE local_dist_constraint_type
    REAL(KIND = dp)               :: del_lam,lambda,scale,scale_old,sigma
    REAL(KIND = dp), DIMENSION(3) :: f,ra_old,rb_old
  END TYPE local_dist_constraint_type

  TYPE local_g3x3_constraint_type
    REAL(KIND = dp)               :: scale,scale_old
    REAL(KIND = dp), DIMENSION(3) :: fa,fb,fc,&
                              ra_old,rb_old,rc_old,&
                              va,vb,vc,del_lam,lambda,&
                              va_old,vb_old,vc_old
  END TYPE local_g3x3_constraint_type

  TYPE local_g4x6_constraint_type
    REAL(KIND = dp)               :: scale,scale_old
    REAL(KIND = dp), DIMENSION(3) :: fa,fb,fc,fd,fe,ff,&
                              ra_old,rb_old,rc_old,rd_old,re_old,rf_old,&
                              va,vb,vc,vd,ve,vf,&
                              va_old,vb_old,vc_old,vd_old,ve_old,vf_old
    REAL(KIND = dp), DIMENSION(6) :: del_lam,lambda
  END TYPE local_g4x6_constraint_type

  TYPE local_molecule_type
    INTEGER, DIMENSION(:), POINTER :: states  ! indices of Kohn-Sham states for molecule
    INTEGER                        :: nstates ! Kohn-Sham states for molecule
  END TYPE local_molecule_type

  TYPE local_constraint_type
    TYPE ( local_dist_constraint_type ), DIMENSION ( : ) , POINTER :: ldist
    TYPE ( local_g3x3_constraint_type ), DIMENSION ( : ) , POINTER :: lg3x3
    TYPE ( local_g4x6_constraint_type ), DIMENSION ( : ) , POINTER :: lg4x6
  END TYPE local_constraint_type

  TYPE molecule_type
    TYPE(molecule_kind_type), POINTER  :: molecule_kind ! pointer to molecule kind information
    TYPE(local_molecule_type), POINTER   :: lmi         ! local molecule information
    TYPE(local_constraint_type), POINTER :: lci         ! local molecule constraint info
    INTEGER                              :: first_atom  ! global index of first atom in molecule
    INTEGER                              :: last_atom   ! global index of last atom in molecule
  END TYPE molecule_type

! *** Public data types ***

  PUBLIC :: local_dist_constraint_type,&
            local_g3x3_constraint_type,&
            local_g4x6_constraint_type,&
            local_constraint_type,&
            local_molecule_type,&
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

    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    INTEGER, INTENT(IN)                      :: nmolecule

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_molecule_set"

    INTEGER                                  :: imolecule, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_set)) CALL deallocate_molecule_set(molecule_set)

    ALLOCATE (molecule_set(nmolecule),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "molecule_set",nmolecule*int_size)

    DO imolecule=1,nmolecule
      NULLIFY (molecule_set(imolecule)%molecule_kind)
      NULLIFY (molecule_set(imolecule)%lmi)
      molecule_set(imolecule)%first_atom = 0
      molecule_set(imolecule)%last_atom = 0
    END DO

  END SUBROUTINE allocate_molecule_set

! *****************************************************************************

  SUBROUTINE deallocate_molecule_set(molecule_set)

!   Purpose: Deallocate a molecule set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_molecule_set"

    INTEGER                                  :: imolecule, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_set)) THEN

      DO imolecule=1,SIZE(molecule_set)
         IF (ASSOCIATED(molecule_set(imolecule)%lmi)) THEN
             IF (ASSOCIATED(molecule_set(imolecule)%lmi%states)) THEN
                 DEALLOCATE (molecule_set(imolecule)%lmi%states,STAT=istat)
                 IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                                  "molecule_set(i)%lmi%states")
             ENDIF
             DEALLOCATE (molecule_set(imolecule)%lmi,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                              "molecule_set(i)%lmi")
         ENDIF
         IF (ASSOCIATED(molecule_set(imolecule)%lci)) THEN
             IF (ASSOCIATED(molecule_set(imolecule)%lci%ldist)) THEN
                 DEALLOCATE (molecule_set(imolecule)%lci%ldist,STAT=istat)
                 IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                                  "molecule_set(i)%lci%ldist")
             ENDIF
             IF (ASSOCIATED(molecule_set(imolecule)%lci%lg3x3)) THEN
                 DEALLOCATE (molecule_set(imolecule)%lci%lg3x3,STAT=istat)
                 IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                                  "molecule_set(i)%lci%lg3x3")
             ENDIF
             IF (ASSOCIATED(molecule_set(imolecule)%lci%lg4x6)) THEN
                 DEALLOCATE (molecule_set(imolecule)%lci%lg4x6,STAT=istat)
                 IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                                  "molecule_set(i)%lci%lg4x6")
             ENDIF
             DEALLOCATE (molecule_set(imolecule)%lci,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                              "molecule_set(i)%lci")
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

  SUBROUTINE get_molecule(molecule,molecule_kind,lmi,lg3x3,lg4x6,ldist,&
                          first_atom,last_atom)

!   Purpose: Get components from a molecule data set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(molecule_kind_type), OPTIONAL, &
      POINTER                                :: molecule_kind
    TYPE(local_molecule_type), OPTIONAL, &
      POINTER                                :: lmi
    TYPE(local_dist_constraint_type), OPTIONAL, &
      POINTER                                :: ldist ( : )
    TYPE(local_g3x3_constraint_type), OPTIONAL, &
      POINTER                                :: lg3x3 ( : )
    TYPE(local_g4x6_constraint_type), OPTIONAL, &
      POINTER                                :: lg4x6 ( : )
    INTEGER, OPTIONAL                        :: first_atom
    INTEGER, OPTIONAL                        :: last_atom

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "get_molecule"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule)) THEN

      IF (PRESENT(first_atom)) first_atom = molecule%first_atom
      IF (PRESENT(last_atom)) last_atom = molecule%last_atom
      IF (PRESENT(molecule_kind)) molecule_kind => molecule%molecule_kind
      IF (PRESENT(lmi))           lmi           => molecule%lmi
      IF ( PRESENT ( ldist ) ) THEN
         IF ( ASSOCIATED ( molecule%lci ) ) THEN 
           ldist => molecule%lci%ldist
         ELSE
           CALL stop_program(routine_name,module_name,__LINE__,&
                           "The pointer lci is not associated")
         ENDIF
      ENDIF
      IF ( PRESENT ( lg3x3 ) ) THEN
         IF ( ASSOCIATED ( molecule%lci ) ) THEN 
           lg3x3 => molecule%lci%lg3x3
         ELSE 
           CALL stop_program(routine_name,module_name,__LINE__,&
                           "The pointer lci is not associated")
         ENDIF
      ENDIF
      IF ( PRESENT ( lg4x6 ) ) THEN
         IF ( ASSOCIATED ( molecule%lci ) ) THEN 
           lg4x6 => molecule%lci%lg4x6
         ELSE 
           CALL stop_program(routine_name,module_name,__LINE__,&
                           "The pointer lci is not associated")
         ENDIF
      ENDIF

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer lci is not associated")

    END IF

  END SUBROUTINE get_molecule

! *****************************************************************************

  SUBROUTINE set_molecule(molecule,molecule_kind,lmi,lci,ldist,lg3x3,lg4x6)

!   Purpose: Set a molecule data set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(molecule_kind_type), OPTIONAL, &
      POINTER                                :: molecule_kind
    TYPE(local_molecule_type), OPTIONAL, &
      POINTER                                :: lmi
    TYPE(local_constraint_type), OPTIONAL, &
      POINTER                                :: lci
    TYPE(local_dist_constraint_type), OPTIONAL, &
      POINTER                                :: ldist ( : )
    TYPE(local_g3x3_constraint_type), OPTIONAL, &
      POINTER                                :: lg3x3 ( : )
    TYPE(local_g4x6_constraint_type), OPTIONAL, &
      POINTER                                :: lg4x6 ( : )

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "set_molecule"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule)) THEN

      IF (PRESENT(molecule_kind)) molecule%molecule_kind => molecule_kind
      IF (PRESENT(lmi))           molecule%lmi           => lmi
      IF (PRESENT(lci))           molecule%lci           => lci
      IF (PRESENT(ldist) ) THEN 
         IF ( ASSOCIATED ( molecule%lci ) ) THEN
           molecule%lci%ldist => ldist
         ELSE
           CALL stop_program(routine_name,module_name,__LINE__,&
                           "The pointer lci is not associated")
         ENDIF
      ENDIF
      IF (PRESENT(lg3x3)) THEN 
         IF ( ASSOCIATED ( molecule%lci ) ) THEN
           molecule%lci%lg3x3 => lg3x3
         ELSE
           CALL stop_program(routine_name,module_name,__LINE__,&
                           "The pointer lci is not associated")
         ENDIF
      ENDIF
      IF (PRESENT(lg4x6)) THEN 
         IF ( ASSOCIATED ( molecule%lci ) ) THEN
           molecule%lci%lg4x6 => lg4x6
         ELSE
           CALL stop_program(routine_name,module_name,__LINE__,&
                           "The pointer lci is not associated")
         ENDIF
      ENDIF
    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule is not associated")

    END IF

  END SUBROUTINE set_molecule

! *****************************************************************************

  SUBROUTINE set_molecule_set(molecule_set,first_atom,last_atom)

!   Purpose: Set a molecule data set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: first_atom
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: last_atom

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_molecule_set"

    INTEGER                                  :: imolecule

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

      IF (PRESENT(last_atom)) THEN

        IF (SIZE(last_atom) /= SIZE(molecule_set)) THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "The sizes of last_atom and molecule_set "//&
                            "are different")
        END IF

        DO imolecule=1,SIZE(molecule_set)
          molecule_set(imolecule)%last_atom = last_atom(imolecule)
        END DO

      END IF

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer molecule_set is not associated")

    END IF

  END SUBROUTINE set_molecule_set

! *****************************************************************************

END MODULE molecule_types_new
