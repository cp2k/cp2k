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
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE force_field_types,               ONLY: &
       bend_kind_type, bond_kind_type, impr_kind_dealloc_ref, impr_kind_type, &
       onfo_kind_dealloc_ref, onfo_kind_type, torsion_kind_dealloc_ref, &
       torsion_kind_type, ub_kind_dealloc_ref, ub_kind_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE colvar_types, ONLY: colvar_type

#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'molecule_kind_types'

! *** Define the derived structure types ***

  TYPE atom_type
    TYPE(atomic_kind_type), POINTER :: atomic_kind
    CHARACTER(LEN=default_string_length)  :: name
  END TYPE atom_type

  TYPE bond_type
    TYPE(bond_kind_type), POINTER :: bond_kind
    INTEGER :: a,b
    CHARACTER(LEN=default_string_length)  :: name
  END TYPE bond_type

  TYPE bend_type
    TYPE(bend_kind_type), POINTER :: bend_kind
    INTEGER :: a,b,c
    CHARACTER(LEN=default_string_length)  :: name
  END TYPE bend_type

  TYPE ub_type
    TYPE(ub_kind_type), POINTER :: ub_kind
    INTEGER :: a,b,c
    CHARACTER(LEN=default_string_length)  :: name
  END TYPE ub_type

  TYPE torsion_type
    TYPE(torsion_kind_type), POINTER :: torsion_kind
    INTEGER :: a,b,c,d
    CHARACTER(LEN=default_string_length)  :: name
  END TYPE torsion_type

  TYPE impr_type
    TYPE(impr_kind_type), POINTER :: impr_kind
    INTEGER :: a,b,c,d
    CHARACTER(LEN=default_string_length)  :: name
  END TYPE impr_type

  TYPE onfo_type
    TYPE(onfo_kind_type), POINTER :: onfo_kind
    INTEGER :: a,b
    CHARACTER(LEN=default_string_length)  :: name
  END TYPE onfo_type

  TYPE colvar_constraint_type
     INTEGER                        :: type_id
     REAL(KIND = dp)                :: expected_value
     INTEGER, POINTER, DIMENSION(:) :: i_atoms
  END TYPE colvar_constraint_type

  TYPE g3x3_constraint_type
    INTEGER                :: a,b,c
    REAL(KIND = dp)        :: dab,dac,dbc
  END TYPE g3x3_constraint_type

  TYPE g4x6_constraint_type
    INTEGER                :: a,b,c,d
    REAL(KIND = dp)        :: dab,dac,dbc,dad,dbd,dcd
  END TYPE g4x6_constraint_type

  TYPE fixd_constraint_type
     INTEGER :: fixd
  END TYPE fixd_constraint_type

  TYPE colvar_counters
     INTEGER :: ndist
  END TYPE colvar_counters

  TYPE molecule_kind_type
    PRIVATE
    TYPE(atom_type), DIMENSION(:), POINTER            :: atom_list
    TYPE(bond_kind_type), DIMENSION(:), POINTER       :: bond_kind_set
    TYPE(bond_type), DIMENSION(:), POINTER            :: bond_list
    TYPE(bend_kind_type), DIMENSION(:), POINTER       :: bend_kind_set
    TYPE(bend_type), DIMENSION(:), POINTER            :: bend_list
    TYPE(ub_kind_type), DIMENSION(:), POINTER         :: ub_kind_set
    TYPE(ub_type), DIMENSION(:), POINTER              :: ub_list
    TYPE(torsion_kind_type), DIMENSION(:), POINTER    :: torsion_kind_set
    TYPE(torsion_type), DIMENSION(:), POINTER         :: torsion_list
    TYPE(impr_kind_type), DIMENSION(:), POINTER       :: impr_kind_set
    TYPE(impr_type), DIMENSION(:), POINTER            :: impr_list
    TYPE(onfo_kind_type), DIMENSION(:), POINTER       :: onfo_kind_set
    TYPE(onfo_type), DIMENSION(:), POINTER            :: onfo_list
    TYPE(colvar_constraint_type), DIMENSION(:), POINTER :: colv_list
    TYPE(g3x3_constraint_type), DIMENSION(:), POINTER :: g3x3_list
    TYPE(g4x6_constraint_type), DIMENSION(:), POINTER :: g4x6_list
    TYPE(fixd_constraint_type), DIMENSION(:), POINTER :: fixd_list
    CHARACTER(LEN=default_string_length)              :: name
    INTEGER                                           :: nhcopt
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
                                                         ng3x3,&
                                                         ng4x6,&
                                                         nfixd,&
                                                         nmolecule
    TYPE(colvar_counters)                             :: ncolv
    INTEGER                                           :: nsgf, nelectron
    INTEGER, DIMENSION(:), POINTER                    :: molecule_list
  END TYPE molecule_kind_type

! *** Public subroutines ***

  PUBLIC :: allocate_molecule_kind_set,&
            deallocate_molecule_kind_set,&
            get_molecule_kind,&
            get_molecule_kind_set,&
            num_ao_el_per_molecule,&
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
            colvar_counters,&
            colvar_constraint_type,&
            g3x3_constraint_type,&
            g4x6_constraint_type,&
            fixd_constraint_type,&
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

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_molecule_kind_set'

    INTEGER                                  :: imolecule_kind, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind_set)) THEN
      CALL deallocate_molecule_kind_set(molecule_kind_set)
    END IF

    ALLOCATE (molecule_kind_set(nmolecule_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "molecule_kind_set",&
                                     nmolecule_kind*int_size)

    DO imolecule_kind=1,nmolecule_kind
      NULLIFY (molecule_kind_set(imolecule_kind)%atom_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%bond_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%bend_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%colv_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%ub_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%ub_kind_set)
      NULLIFY (molecule_kind_set(imolecule_kind)%impr_kind_set)      
      NULLIFY (molecule_kind_set(imolecule_kind)%impr_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%onfo_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%onfo_kind_set)
      NULLIFY (molecule_kind_set(imolecule_kind)%g3x3_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%g4x6_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%fixd_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%torsion_list)
      NULLIFY (molecule_kind_set(imolecule_kind)%bond_kind_set)
      NULLIFY (molecule_kind_set(imolecule_kind)%bend_kind_set)
      NULLIFY (molecule_kind_set(imolecule_kind)%torsion_kind_set)
      NULLIFY (molecule_kind_set(imolecule_kind)%charge)
      NULLIFY (molecule_kind_set(imolecule_kind)%mass)
      molecule_kind_set(imolecule_kind)%name = ""
      molecule_kind_set(imolecule_kind)%nhcopt = -1
      molecule_kind_set(imolecule_kind)%kind_number = imolecule_kind
      molecule_kind_set(imolecule_kind)%natom = 0
      molecule_kind_set(imolecule_kind)%nbend = 0
      molecule_kind_set(imolecule_kind)%nbond = 0
      molecule_kind_set(imolecule_kind)%nimpr = 0
      molecule_kind_set(imolecule_kind)%nonfo = 0
      molecule_kind_set(imolecule_kind)%nub = 0
      molecule_kind_set(imolecule_kind)%ncolv%ndist = 0
      molecule_kind_set(imolecule_kind)%ng3x3 = 0
      molecule_kind_set(imolecule_kind)%ng4x6 = 0
      molecule_kind_set(imolecule_kind)%nfixd = 0
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

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_molecule_kind_set'

    INTEGER                                  :: i, imolecule_kind, istat, &
                                                nmolecule_kind, j

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind_set)) THEN

      nmolecule_kind = SIZE(molecule_kind_set)

      DO imolecule_kind=1,nmolecule_kind

        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%atom_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%atom_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%atom_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%bend_kind_set)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%bend_kind_set,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%bend_kind_set")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%bend_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%bend_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%bend_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%ub_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%ub_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%ub_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%ub_kind_set)) THEN
           CALL ub_kind_dealloc_ref(molecule_kind_set(imolecule_kind)%ub_kind_set)
        END IF        
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%impr_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%impr_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%impr_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%impr_kind_set)) THEN
           DO i = 1, SIZE(molecule_kind_set(imolecule_kind)%impr_kind_set)
              CALL impr_kind_dealloc_ref(molecule_kind_set(imolecule_kind)%impr_kind_set(i))
           END DO
           DEALLOCATE (molecule_kind_set(imolecule_kind)%impr_kind_set,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%onfo_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%onfo_list)) THEN
           DEALLOCATE(molecule_kind_set(imolecule_kind)%onfo_list, STAT=istat)
           IF (istat /= 0) THEN
              CALL stop_memory(routineN,moduleN,__LINE__,&
                   "molecule_kind_set(imolecule_kind)%onfo_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%onfo_kind_set)) THEN
          CALL onfo_kind_dealloc_ref(molecule_kind_set(imolecule_kind)%onfo_kind_set)
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%bond_kind_set)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%bond_kind_set,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%bond_kind_set")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%bond_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%bond_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%bond_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%colv_list)) THEN
           DO j = 1, SIZE(molecule_kind_set(imolecule_kind)%colv_list)
              DEALLOCATE (molecule_kind_set(imolecule_kind)%colv_list(j)%i_atoms,STAT=istat)
              IF (istat /= 0) THEN
                 CALL stop_memory(routineN,moduleN,__LINE__,&
                      "molecule_kind_set(imolecule_kind)%colv_list%i_atoms")        
              END IF
           END DO
           DEALLOCATE (molecule_kind_set(imolecule_kind)%colv_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                           "molecule_kind_set(imolecule_kind)%colv_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%g3x3_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%g3x3_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%g3x3_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%g4x6_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%g4x6_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%g4x6_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%fixd_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%fixd_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%fixd_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%torsion_kind_set)) THEN
           DO i = 1, SIZE(molecule_kind_set(imolecule_kind)%torsion_kind_set)
              CALL torsion_kind_dealloc_ref(molecule_kind_set(imolecule_kind)%torsion_kind_set(i))
           END DO
           DEALLOCATE (molecule_kind_set(imolecule_kind)%torsion_kind_set,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%torsion_kind_set")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%torsion_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%torsion_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%torsion_list")
           END IF
        END IF
        IF (ASSOCIATED(molecule_kind_set(imolecule_kind)%molecule_list)) THEN
           DEALLOCATE (molecule_kind_set(imolecule_kind)%molecule_list,STAT=istat)
           IF (istat /= 0) THEN
             CALL stop_memory(routineN,moduleN,__LINE__,&
                              "molecule_kind_set(imolecule_kind)%molecule_list")
           END IF
        ENDIF
      END DO

      DEALLOCATE (molecule_kind_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "molecule_kind_set")

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer molecule_kind_set is not associated "//&
                        "and cannot be deallocated")

    END IF

  END SUBROUTINE deallocate_molecule_kind_set

! *****************************************************************************

  SUBROUTINE get_molecule_kind(molecule_kind,atom_list,bond_list,bend_list,&
                               ub_list,impr_list,onfo_list,colv_list,fixd_list,&
                               g3x3_list,g4x6_list,torsion_list,name,nhcopt,&
                               mass,charge,kind_number,natom,nbend,nbond,nub,&
                               nimpr,nonfo,nconstraint,nfixd,&
                               ncolv,ng3x3,ng4x6,nmolecule,nsgf,ntorsion,&
                               molecule_list,nelectron,&
                               bond_kind_set,bend_kind_set,ub_kind_set,&
                               impr_kind_set,onfo_kind_set,torsion_kind_set)

!   Purpose: Get informations about a molecule kind.

!   History: - Creation (27.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atom_list
    TYPE(bond_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bond_list
    TYPE(bend_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bend_list
    TYPE(ub_type), DIMENSION(:), OPTIONAL, &
      POINTER                                :: ub_list
    TYPE(impr_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: impr_list
    TYPE(onfo_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: onfo_list
    TYPE(colvar_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: colv_list
    TYPE(fixd_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: fixd_list
    TYPE(g3x3_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: g3x3_list
    TYPE(g4x6_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: g4x6_list
    TYPE(torsion_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: torsion_list
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    INTEGER, INTENT(OUT), OPTIONAL           :: nhcopt
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: mass, charge
    INTEGER, INTENT(OUT), OPTIONAL :: kind_number, natom, nbend, nbond, nub, &
      nimpr, nonfo, nconstraint, nfixd, ng3x3, ng4x6, nmolecule, nsgf, &
      ntorsion
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: molecule_list
    INTEGER, INTENT(OUT), OPTIONAL           :: nelectron
    TYPE(colvar_counters), OPTIONAL, INTENT(out) :: ncolv
    TYPE(bond_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bond_kind_set
    TYPE(bend_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bend_kind_set
    TYPE(ub_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ub_kind_set
    TYPE(impr_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: impr_kind_set
    TYPE(onfo_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: onfo_kind_set
    TYPE(torsion_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: torsion_kind_set

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'get_molecule_kind'

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
      IF (PRESENT(colv_list)) colv_list => molecule_kind%colv_list
      IF (PRESENT(g3x3_list)) g3x3_list => molecule_kind%g3x3_list
      IF (PRESENT(g4x6_list)) g4x6_list => molecule_kind%g4x6_list
      IF (PRESENT(fixd_list)) fixd_list => molecule_kind%fixd_list
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
      IF (PRESENT(nconstraint)) nconstraint = molecule_kind%ncolv%ndist +&
                                              3*molecule_kind%ng3x3 +&
                                              6*molecule_kind%ng4x6 
      IF (PRESENT(ncolv)) ncolv = molecule_kind%ncolv
      IF (PRESENT(ng3x3)) ng3x3= molecule_kind%ng3x3
      IF (PRESENT(ng4x6)) ng4x6= molecule_kind%ng4x6
      IF (PRESENT(nfixd)) nfixd= molecule_kind%nfixd
      IF (PRESENT(nmolecule)) nmolecule = molecule_kind%nmolecule
      IF (PRESENT(ntorsion)) ntorsion= molecule_kind%ntorsion
      IF (PRESENT(nsgf)) nsgf = molecule_kind%nsgf
      IF (PRESENT(nelectron)) nelectron = molecule_kind%nelectron
      IF (PRESENT(molecule_list)) molecule_list => molecule_kind%molecule_list

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer molecule_kind is not associated")

    END IF

  END SUBROUTINE get_molecule_kind

! *****************************************************************************

  SUBROUTINE get_molecule_kind_set(molecule_kind_set,maxatom,natom,&
                                   nbond,nbend,nub,ntorsion,nimpr,nonfo,&
                                   nconstraint,nmolecule)

!   Purpose: Get informations about a molecule kind set.

!   History: - Creation (27.08.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    INTEGER, INTENT(OUT), OPTIONAL           :: maxatom, natom, nbond, nbend, &
                                                nub, ntorsion, nimpr, nonfo, &
                                                nconstraint, nmolecule

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_molecule_kind_set'

    INTEGER                                  :: ibend, ibond, iimpr, &
                                                imolecule_kind, ionfo, &
                                                itorsion, iub, na, nc, nfixd, &
                                                nm, nmolecule_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind_set)) THEN

      IF (PRESENT(maxatom)) maxatom = 0
      IF (PRESENT(natom)) natom = 0
      IF (PRESENT(nbond)) nbond = 0
      IF (PRESENT(nbend)) nbend = 0
      IF (PRESENT(nub)) nub = 0
      IF (PRESENT(ntorsion)) ntorsion = 0
      IF (PRESENT(nimpr)) nimpr = 0
      IF (PRESENT(nimpr)) nimpr = 0
      IF (PRESENT(nonfo)) nonfo = 0
      IF (PRESENT(nconstraint)) nconstraint = 0
      IF (PRESENT(nmolecule)) nmolecule = 0

      nmolecule_kind = SIZE(molecule_kind_set)

      DO imolecule_kind=1,nmolecule_kind

        molecule_kind => molecule_kind_set(imolecule_kind)

        CALL get_molecule_kind(molecule_kind=molecule_kind,&
                               natom=na,&
                               nbond=ibond,&
                               nbend=ibend,&
                               nub=iub,&
                               ntorsion=itorsion,&
                               nimpr=iimpr,&
                               nonfo=ionfo,&
                               nconstraint=nc,&
                               nfixd=nfixd,&
                               nmolecule=nm)

        IF (PRESENT(maxatom)) maxatom = MAX(maxatom,na)
        IF (PRESENT(natom)) natom = natom + na*nm
        IF (PRESENT(nbond)) nbond = nbond + ibond*nm
        IF (PRESENT(nbend)) nbend = nbend + ibend*nm
        IF (PRESENT(nub)) nub = nub + iub*nm
        IF (PRESENT(ntorsion)) ntorsion = ntorsion + itorsion*nm
        IF (PRESENT(nimpr)) nimpr = nimpr + iimpr*nm
        IF (PRESENT(nonfo)) nonfo = nonfo + ionfo*nm
        IF (PRESENT(nconstraint)) nconstraint = nconstraint + nc*nm + 3*nfixd
        IF (PRESENT(nmolecule)) nmolecule = nmolecule + nm

      END DO

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer molecule_kind_set is not associated")

    END IF

  END SUBROUTINE get_molecule_kind_set

! *****************************************************************************

  SUBROUTINE num_ao_el_per_molecule(molecule_kind_set)

!   Purpose: Read a molecule kind data set from the input file.

!   History: - Creation (22.11.2004,MI)

!   ***************************************************************************

    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set

    CHARACTER(LEN=*), PARAMETER :: routineN = 'num_ao_el_per_molecule'

    INTEGER                                  :: iatom, imol, n_ao, natom, &
                                                nelectron, nmol_kind, nsgf
    REAL(KIND=dp)                            :: zeff, zeff_correction
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind_set)) THEN

    nmol_kind = SIZE(molecule_kind_set,1)
    natom = 0
!   *** Initialize the molecule kind data structure ***

    DO imol = 1,nmol_kind
      molecule_kind => molecule_kind_set(imol)
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             natom=natom)
      nelectron = 0
      n_ao = 0
      DO iatom=1,natom
        atomic_kind => molecule_kind%atom_list(iatom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             all_potential=all_potential,&
                             gth_potential=gth_potential,&
                             orb_basis_set=orb_basis_set)
        IF (ASSOCIATED(all_potential)) THEN
           CALL get_potential(potential=all_potential,zeff=zeff,&
                zeff_correction=zeff_correction)
        ELSE IF (ASSOCIATED(gth_potential)) THEN
           CALL get_potential(potential=gth_potential,zeff=zeff,&
                zeff_correction=zeff_correction)
        ELSE
           zeff = 0.0_dp
           zeff_correction = 0.0_dp
        END IF
        nelectron = nelectron + INT(zeff-zeff_correction)
  
        IF (ASSOCIATED(orb_basis_set)) THEN
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,nsgf=nsgf)
        ELSE 
          nsgf = 0
        ENDIF
        n_ao = n_ao + nsgf 
      END DO  ! iatom
      CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           nelectron=nelectron, nsgf=n_ao)
    END DO  ! imol
    END IF

  END SUBROUTINE num_ao_el_per_molecule

! *****************************************************************************

  SUBROUTINE set_molecule_kind(molecule_kind,name,mass,charge,kind_number,&
                               molecule_list,atom_list,nhcopt,nbond,bond_list,&
                               nbend,bend_list,nub,ub_list,nimpr,impr_list,&
                               nonfo,onfo_list,ntorsion,torsion_list,fixd_list,&
                               ncolv,colv_list,ng3x3,g3x3_list,ng4x6,nfixd,g4x6_list,&
                               bond_kind_set,bend_kind_set,ub_kind_set,&
                               torsion_kind_set,impr_kind_set,onfo_kind_set,&
                               nelectron,nsgf)

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
    INTEGER, INTENT(IN), OPTIONAL            :: nhcopt, nbond
    TYPE(bond_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bond_list
    INTEGER, INTENT(IN), OPTIONAL            :: nbend
    TYPE(bend_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bend_list
    INTEGER, INTENT(IN), OPTIONAL            :: nub
    TYPE(ub_type), DIMENSION(:), OPTIONAL, &
      POINTER                                :: ub_list
    INTEGER, INTENT(IN), OPTIONAL            :: nimpr
    TYPE(impr_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: impr_list
    INTEGER, INTENT(IN), OPTIONAL            :: nonfo
    TYPE(onfo_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: onfo_list
    INTEGER, INTENT(IN), OPTIONAL            :: ntorsion
    TYPE(torsion_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: torsion_list
    TYPE(fixd_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: fixd_list
    TYPE(colvar_counters), INTENT(IN), OPTIONAL   :: ncolv
    TYPE(colvar_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: colv_list
    INTEGER, INTENT(IN), OPTIONAL            :: ng3x3
    TYPE(g3x3_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: g3x3_list
    INTEGER, INTENT(IN), OPTIONAL            :: ng4x6, nfixd
    TYPE(g4x6_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: g4x6_list
    TYPE(bond_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bond_kind_set
    TYPE(bend_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: bend_kind_set
    TYPE(ub_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ub_kind_set
    TYPE(torsion_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: torsion_kind_set
    TYPE(impr_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: impr_kind_set
    TYPE(onfo_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: onfo_kind_set
    INTEGER, INTENT(IN), OPTIONAL            :: nelectron, nsgf

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'set_molecule_kind'

    INTEGER                                  :: n

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(molecule_kind)) THEN

      IF (PRESENT(atom_list)) THEN
        n = SIZE(atom_list)
        molecule_kind%natom = n
        molecule_kind%atom_list => atom_list
      END IF
      IF (PRESENT(name)) molecule_kind%name = name
      IF (PRESENT(nhcopt)) molecule_kind%nhcopt = nhcopt
      IF (PRESENT(mass)) molecule_kind%mass => mass
      IF (PRESENT(charge)) molecule_kind%charge => charge
      IF (PRESENT(kind_number)) molecule_kind%kind_number = kind_number
      IF (PRESENT(nbond)) molecule_kind%nbond = nbond
      IF (PRESENT(bond_list)) molecule_kind%bond_list => bond_list
      IF (PRESENT(nbend)) molecule_kind%nbend = nbend
      IF (PRESENT(nelectron)) molecule_kind%nelectron = nelectron
      IF (PRESENT(nsgf)) molecule_kind%nsgf = nsgf
      IF (PRESENT(bend_list)) molecule_kind%bend_list => bend_list
      IF (PRESENT(nub)) molecule_kind%nub = nub
      IF (PRESENT(ub_list)) molecule_kind%ub_list => ub_list
      IF (PRESENT(ntorsion)) molecule_kind%ntorsion = ntorsion
      IF (PRESENT(torsion_list)) molecule_kind%torsion_list => torsion_list
      IF (PRESENT(nimpr)) molecule_kind%nimpr = nimpr
      IF (PRESENT(impr_list)) molecule_kind%impr_list => impr_list
      IF (PRESENT(nonfo)) molecule_kind%nonfo = nonfo
      IF (PRESENT(onfo_list)) molecule_kind%onfo_list => onfo_list
      IF (PRESENT(ncolv)) molecule_kind%ncolv = ncolv
      IF (PRESENT(colv_list)) molecule_kind%colv_list => colv_list
      IF (PRESENT(ng3x3)) molecule_kind%ng3x3 = ng3x3
      IF (PRESENT(g3x3_list)) molecule_kind%g3x3_list => g3x3_list
      IF (PRESENT(ng4x6)) molecule_kind%ng4x6 = ng4x6
      IF (PRESENT(nfixd)) molecule_kind%nfixd = nfixd
      IF (PRESENT(g4x6_list)) molecule_kind%g4x6_list => g4x6_list
      IF (PRESENT(fixd_list)) molecule_kind%fixd_list => fixd_list
      IF (PRESENT(bond_kind_set)) molecule_kind%bond_kind_set => bond_kind_set
      IF (PRESENT(bend_kind_set)) molecule_kind%bend_kind_set => bend_kind_set
      IF (PRESENT(ub_kind_set)) molecule_kind%ub_kind_set => ub_kind_set
      IF (PRESENT(torsion_kind_set)) molecule_kind%torsion_kind_set => torsion_kind_set
      IF (PRESENT(impr_kind_set)) molecule_kind%impr_kind_set => impr_kind_set
      IF (PRESENT(onfo_kind_set)) molecule_kind%onfo_kind_set => onfo_kind_set
      IF (PRESENT(molecule_list)) THEN
        n = SIZE(molecule_list)
        molecule_kind%nmolecule = n
        molecule_kind%molecule_list => molecule_list
      END IF

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer molecule_kind is not associated")

    END IF

  END SUBROUTINE set_molecule_kind

! *****************************************************************************

  SUBROUTINE write_molecule_kind(molecule_kind,output_unit,error)

!   Purpose: Write a molecule kind data set to the output unit.

!   History: - Creation (24.09.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    INTEGER, INTENT(in)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_molecule_kind', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name
    INTEGER                                  :: iatom, imolecule, natom, &
                                                nmolecule
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    failure = .FALSE.
    IF (output_unit>0) THEN
      IF (ASSOCIATED(molecule_kind)) THEN
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
        WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
          "Number of bonds:       ",molecule_kind%nbond
        WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
          "Number of bends:       ",molecule_kind%nbend
        WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
          "Number of Urey-Bradley:",molecule_kind%nub
        WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
          "Number of torsions:    ",molecule_kind%ntorsion
        WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
          "Number of improper:    ",molecule_kind%nimpr
        WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
          "Number of one-four:    ",molecule_kind%nonfo

      ELSE         
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
      END IF

    END IF
  END SUBROUTINE write_molecule_kind

! *****************************************************************************

  SUBROUTINE write_molecule_kind_set(molecule_kind_set,force_env_section,error)

!   Purpose: Write a moleculeatomic kind set data set to the output unit.

!   History: - Creation (24.09.2003,MK)

!   ***************************************************************************

    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_molecule_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imolecule_kind, nbend, nbond, &
                                                nimpr, nmolecule_kind, nonfo, &
                                                ntors, nub, output_unit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)       
    output_unit = cp_print_key_unit_nr(logger,force_env_section,&
         "SUBSYS%PRINT%MOLECULES",extension=".Log",error=error)
    IF (output_unit>0) THEN
      IF (ASSOCIATED(molecule_kind_set)) THEN
        WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") "MOLECULE KIND INFORMATION"

        nmolecule_kind = SIZE(molecule_kind_set)

        DO imolecule_kind=1,nmolecule_kind
          molecule_kind => molecule_kind_set(imolecule_kind)
          CALL write_molecule_kind(molecule_kind,output_unit,error)
        END DO


        WRITE(UNIT=output_unit,FMT="(/,/,T2,A)") "MOLECULE KIND SET INFORMATION"
        CALL get_molecule_kind_set(molecule_kind_set=molecule_kind_set,&
                                   nbond=nbond,&
                                   nbend=nbend,&
                                   nub=nub,&
                                   ntorsion=ntors,&
                                   nimpr=nimpr,&
                                   nonfo=nonfo)
        WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
          "Total Number of bonds:       ",nbond
        WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
          "Total Number of bends:       ",nbend
        WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
          "Total Number of Urey-Bradley:",nub
        WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
          "Total Number of torsions:    ",ntors
        WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
          "Total Number of improper:    ",nimpr
        WRITE (UNIT=output_unit,FMT="(1X,A30,I6)")&
          "Total Number of one-four:    ",nonfo
      ELSE
         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
      END IF

    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
         "SUBSYS%PRINT%MOLECULES",error=error)
  END SUBROUTINE write_molecule_kind_set

! *****************************************************************************


END MODULE molecule_kind_types
