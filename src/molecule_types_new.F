!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Define the data structure for the molecule information.
!> \par History
!>      JGH (22.05.2004) add last_atom information
!>      Teodoro Laino [tlaino] 12.2008 - Preparing for VIRTUAL SITE constraints
!>                                       (patch by Marcel Baer)
!> \author MK (29.08.2003)
! *****************************************************************************
MODULE molecule_types_new

  USE colvar_types,                    ONLY: colvar_counters,&
                                             colvar_release,&
                                             colvar_type
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE molecule_kind_types,             ONLY: colvar_constraint_type,&
                                             fixd_constraint_type,&
                                             g3x3_constraint_type,&
                                             g4x6_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type,&
                                             vsite_constraint_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'molecule_types_new'

! *** Data types ***
! *****************************************************************************
  TYPE local_colvar_constraint_type
     LOGICAL                       :: init
     TYPE(colvar_type),POINTER     :: colvar
     TYPE(colvar_type),POINTER     :: colvar_old
     REAL(KIND = dp)               :: lambda, sigma
  END TYPE local_colvar_constraint_type

! *****************************************************************************
  TYPE local_g3x3_constraint_type
    LOGICAL                       :: init
    REAL(KIND = dp)               :: scale,scale_old, imass1, imass2, imass3
    REAL(KIND = dp), DIMENSION(3) :: fa,fb,fc,f_roll1,f_roll2,f_roll3,&
                                     ra_old,rb_old,rc_old,&
                                     va,vb,vc,lambda, del_lambda, lambda_old,&
                                     r0_12, r0_13, r0_23
    REAL(KIND = dp), DIMENSION(3,3) :: amat
  END TYPE local_g3x3_constraint_type

! *****************************************************************************
  TYPE local_g4x6_constraint_type
    LOGICAL                       :: init
    REAL(KIND = dp)               :: scale,scale_old, imass1, imass2, imass3, imass4
    REAL(KIND = dp), DIMENSION(3) :: fa,fb,fc,fd,fe,ff,&
                                     f_roll1, f_roll2, f_roll3, f_roll4, f_roll5, f_roll6,&
                                     ra_old,rb_old,rc_old,rd_old,re_old,rf_old,&
                                     va,vb,vc,vd,ve,vf,&
                                     r0_12, r0_13, r0_14, r0_23, r0_24, r0_34
    REAL(KIND = dp), DIMENSION(6)   :: lambda, del_lambda, lambda_old
    REAL(KIND = dp), DIMENSION(6,6) :: amat
  END TYPE local_g4x6_constraint_type

! *****************************************************************************
  TYPE local_molecule_type
    INTEGER, DIMENSION(:), POINTER :: states  ! indices of Kohn-Sham states for molecule
    INTEGER                        :: nstates ! Kohn-Sham states for molecule
  END TYPE local_molecule_type

! *****************************************************************************
  TYPE local_constraint_type
     TYPE(local_colvar_constraint_type), DIMENSION ( : ) , POINTER :: lcolv
     TYPE(local_g3x3_constraint_type), DIMENSION ( : ) ,   POINTER :: lg3x3
     TYPE(local_g4x6_constraint_type), DIMENSION ( : ) ,   POINTER :: lg4x6
  END TYPE local_constraint_type

! *****************************************************************************
  TYPE global_constraint_type
     TYPE(colvar_counters)                    :: ncolv
     INTEGER                                  :: ntot,  nrestraint
     INTEGER                                  :: ng3x3, ng3x3_restraint
     INTEGER                                  :: ng4x6, ng4x6_restraint
     INTEGER                                  :: nvsite, nvsite_restraint
     TYPE(fixd_constraint_type), DIMENSION(:), POINTER   :: fixd_list
     TYPE(colvar_constraint_type), DIMENSION(:), POINTER :: colv_list
     TYPE(g3x3_constraint_type), DIMENSION(:), POINTER   :: g3x3_list
     TYPE(g4x6_constraint_type), DIMENSION(:), POINTER   :: g4x6_list
     TYPE(vsite_constraint_type), DIMENSION(:), POINTER  :: vsite_list
     TYPE(local_colvar_constraint_type), DIMENSION ( : ) , POINTER :: lcolv
     TYPE(local_g3x3_constraint_type),   DIMENSION ( : ) , POINTER :: lg3x3
     TYPE(local_g4x6_constraint_type),   DIMENSION ( : ) , POINTER :: lg4x6
  END TYPE global_constraint_type

! *****************************************************************************
  TYPE molecule_type
    TYPE(molecule_kind_type), POINTER    :: molecule_kind ! pointer to molecule kind information
    TYPE(local_molecule_type), POINTER   :: lmi         ! local molecule information
    TYPE(local_constraint_type), POINTER :: lci         ! local molecule constraint info
    INTEGER                              :: first_atom  ! global index of first atom in molecule
    INTEGER                              :: last_atom   ! global index of last atom in molecule
    INTEGER                              :: first_shell ! global index of first shell atom in molecule
    INTEGER                              :: last_shell  ! global index of last shell atom in molecule
  END TYPE molecule_type

! *** Public data types ***

  PUBLIC :: local_colvar_constraint_type,&
            local_g3x3_constraint_type,&
            local_g4x6_constraint_type,&
            local_constraint_type,&
            local_molecule_type,&
            global_constraint_type,&
            molecule_type

! *** Public subroutines ***

  PUBLIC :: allocate_molecule_set,&
            deallocate_molecule_set,&
            get_molecule,&
            set_molecule,&
            set_molecule_set,&
            molecule_of_atom,&
            get_molecule_set_info

CONTAINS

! *****************************************************************************
!> \brief   Allocate a molecule set.
!> \author  MK
!> \date    29.08.2003
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_molecule_set(molecule_set,nmolecule,error)
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    INTEGER, INTENT(IN)                      :: nmolecule
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_molecule_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imolecule, istat

    IF (ASSOCIATED(molecule_set)) CALL deallocate_molecule_set(molecule_set,error)

    ALLOCATE (molecule_set(nmolecule),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "molecule_set",nmolecule*int_size)

    DO imolecule=1,nmolecule
      NULLIFY (molecule_set(imolecule)%molecule_kind)
      NULLIFY (molecule_set(imolecule)%lmi)
      NULLIFY (molecule_set(imolecule)%lci)

      molecule_set(imolecule)%first_atom  = 0
      molecule_set(imolecule)%last_atom   = 0
      molecule_set(imolecule)%first_shell = 0
      molecule_set(imolecule)%last_shell  = 0
    END DO

  END SUBROUTINE allocate_molecule_set

! *****************************************************************************
!> \brief   Deallocate a molecule set.
!> \author  MK
!> \date    29.08.2003
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_molecule_set(molecule_set,error)
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_molecule_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imolecule, istat, j

    IF (ASSOCIATED(molecule_set)) THEN

      DO imolecule=1,SIZE(molecule_set)
         IF (ASSOCIATED(molecule_set(imolecule)%lmi)) THEN
             IF (ASSOCIATED(molecule_set(imolecule)%lmi%states)) THEN
                 DEALLOCATE (molecule_set(imolecule)%lmi%states,STAT=istat)
                 IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                                  "molecule_set(i)%lmi%states")
             ENDIF
             DEALLOCATE (molecule_set(imolecule)%lmi,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                              "molecule_set(i)%lmi")
         ENDIF
         IF (ASSOCIATED(molecule_set(imolecule)%lci)) THEN
            IF (ASSOCIATED(molecule_set(imolecule)%lci%lcolv)) THEN
               DO j = 1, SIZE(molecule_set(imolecule)%lci%lcolv)
                  CALL colvar_release(molecule_set(imolecule)%lci%lcolv(j)%colvar,error=error)
                  CALL colvar_release(molecule_set(imolecule)%lci%lcolv(j)%colvar_old,error=error)
                  NULLIFY(molecule_set(imolecule)%lci%lcolv(j)%colvar)
                  NULLIFY(molecule_set(imolecule)%lci%lcolv(j)%colvar_old)
               END DO
               DEALLOCATE (molecule_set(imolecule)%lci%lcolv,STAT=istat)
               IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                                  "molecule_set(i)%lci%lcolv")
             ENDIF
             IF (ASSOCIATED(molecule_set(imolecule)%lci%lg3x3)) THEN
                 DEALLOCATE (molecule_set(imolecule)%lci%lg3x3,STAT=istat)
                 IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                                  "molecule_set(i)%lci%lg3x3")
             ENDIF
             IF (ASSOCIATED(molecule_set(imolecule)%lci%lg4x6)) THEN
                 DEALLOCATE (molecule_set(imolecule)%lci%lg4x6,STAT=istat)
                 IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                                  "molecule_set(i)%lci%lg4x6")
             ENDIF
             DEALLOCATE (molecule_set(imolecule)%lci,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                              "molecule_set(i)%lci")
         ENDIF
      ENDDO
      DEALLOCATE (molecule_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "molecule_set")

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer molecule_set is not associated and "//&
                        "cannot be deallocated")

    END IF

  END SUBROUTINE deallocate_molecule_set

! *****************************************************************************
!> \brief   Get components from a molecule data set.
!> \author  MK
!> \date    29.08.2003
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_molecule(molecule,molecule_kind,lmi,lg3x3,lg4x6,lcolv,&
       first_atom,last_atom,first_shell,last_shell)

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(molecule_kind_type), OPTIONAL, &
      POINTER                                :: molecule_kind
    TYPE(local_molecule_type), OPTIONAL, &
      POINTER                                :: lmi
    TYPE(local_g3x3_constraint_type), &
      OPTIONAL, POINTER                      :: lg3x3( : )
    TYPE(local_g4x6_constraint_type), &
      OPTIONAL, POINTER                      :: lg4x6( : )
    TYPE(local_colvar_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: lcolv
    INTEGER, OPTIONAL                        :: first_atom, last_atom, &
                                                first_shell, last_shell

    CHARACTER(len=*), PARAMETER :: routineN = 'get_molecule', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(molecule)) THEN

      IF (PRESENT(first_atom)) first_atom = molecule%first_atom
      IF (PRESENT(last_atom)) last_atom = molecule%last_atom
      IF (PRESENT(first_shell)) first_shell = molecule%first_shell
      IF (PRESENT(last_shell)) last_shell = molecule%last_shell
      IF (PRESENT(molecule_kind)) molecule_kind => molecule%molecule_kind
      IF (PRESENT(lmi))           lmi           => molecule%lmi
      IF ( PRESENT ( lcolv ) ) THEN
         IF ( ASSOCIATED ( molecule%lci ) ) THEN
           lcolv => molecule%lci%lcolv
         ELSE
           CALL stop_program(routineN,moduleN,__LINE__,&
                           "The pointer lci is not associated")
         ENDIF
      ENDIF
      IF ( PRESENT ( lg3x3 ) ) THEN
         IF ( ASSOCIATED ( molecule%lci ) ) THEN
           lg3x3 => molecule%lci%lg3x3
         ELSE
           CALL stop_program(routineN,moduleN,__LINE__,&
                           "The pointer lci is not associated")
         ENDIF
      ENDIF
      IF ( PRESENT ( lg4x6 ) ) THEN
         IF ( ASSOCIATED ( molecule%lci ) ) THEN
           lg4x6 => molecule%lci%lg4x6
         ELSE
           CALL stop_program(routineN,moduleN,__LINE__,&
                           "The pointer lci is not associated")
         ENDIF
      ENDIF

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer lci is not associated")

    END IF

  END SUBROUTINE get_molecule

! *****************************************************************************
!> \brief   Set a molecule data set.
!> \author  MK
!> \date    29.08.2003
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_molecule(molecule,molecule_kind,lmi,lci,lcolv,lg3x3,lg4x6)
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(molecule_kind_type), OPTIONAL, &
      POINTER                                :: molecule_kind
    TYPE(local_molecule_type), OPTIONAL, &
      POINTER                                :: lmi
    TYPE(local_constraint_type), OPTIONAL, &
      POINTER                                :: lci
    TYPE(local_colvar_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: lcolv
    TYPE(local_g3x3_constraint_type), &
      OPTIONAL, POINTER                      :: lg3x3( : )
    TYPE(local_g4x6_constraint_type), &
      OPTIONAL, POINTER                      :: lg4x6( : )

    CHARACTER(len=*), PARAMETER :: routineN = 'set_molecule', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(molecule)) THEN

      IF (PRESENT(molecule_kind)) molecule%molecule_kind => molecule_kind
      IF (PRESENT(lmi))           molecule%lmi           => lmi
      IF (PRESENT(lci))           molecule%lci           => lci
      IF (PRESENT(lcolv) ) THEN
         IF ( ASSOCIATED ( molecule%lci ) ) THEN
           molecule%lci%lcolv => lcolv
         ELSE
           CALL stop_program(routineN,moduleN,__LINE__,&
                           "The pointer lci is not associated")
         ENDIF
      ENDIF
      IF (PRESENT(lg3x3)) THEN
         IF ( ASSOCIATED ( molecule%lci ) ) THEN
           molecule%lci%lg3x3 => lg3x3
         ELSE
           CALL stop_program(routineN,moduleN,__LINE__,&
                           "The pointer lci is not associated")
         ENDIF
      ENDIF
      IF (PRESENT(lg4x6)) THEN
         IF ( ASSOCIATED ( molecule%lci ) ) THEN
           molecule%lci%lg4x6 => lg4x6
         ELSE
           CALL stop_program(routineN,moduleN,__LINE__,&
                           "The pointer lci is not associated")
         ENDIF
      ENDIF
    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer molecule is not associated")

    END IF

  END SUBROUTINE set_molecule

! *****************************************************************************
!> \brief   Set a molecule data set.
!> \author  MK
!> \date    29.08.2003
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_molecule_set(molecule_set,first_atom,last_atom)
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: first_atom, last_atom

    CHARACTER(len=*), PARAMETER :: routineN = 'set_molecule_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imolecule

    IF (ASSOCIATED(molecule_set)) THEN

      IF (PRESENT(first_atom)) THEN

        IF (SIZE(first_atom) /= SIZE(molecule_set)) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "The sizes of first_atom and molecule_set "//&
                            "are different")
        END IF

        DO imolecule=1,SIZE(molecule_set)
          molecule_set(imolecule)%first_atom = first_atom(imolecule)
        END DO

      END IF

      IF (PRESENT(last_atom)) THEN

        IF (SIZE(last_atom) /= SIZE(molecule_set)) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "The sizes of last_atom and molecule_set "//&
                            "are different")
        END IF

        DO imolecule=1,SIZE(molecule_set)
          molecule_set(imolecule)%last_atom = last_atom(imolecule)
        END DO

      END IF

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer molecule_set is not associated")

    END IF

  END SUBROUTINE set_molecule_set

! *****************************************************************************
!> \brief   finds for each atom the molecule it belongs to
! *****************************************************************************
  SUBROUTINE molecule_of_atom(molecule_set,atom_to_mol,error)
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    INTEGER, DIMENSION(:), INTENT(OUT)       :: atom_to_mol
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'molecule_of_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, iatom, imol, &
                                                last_atom
    TYPE(molecule_type), POINTER             :: molecule

    DO imol = 1,SIZE(molecule_set)
       molecule => molecule_set(imol)
       CALL get_molecule ( molecule=molecule, first_atom = first_atom, last_atom = last_atom )
       DO iatom = first_atom, last_atom
          atom_to_mol(iatom)=imol
       ENDDO  ! iatom
    END DO  ! imol

  END SUBROUTINE molecule_of_atom

! *****************************************************************************
!> \brief returns information about molecules in the set.
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE get_molecule_set_info(molecule_set,atom_to_mol,mol_to_first_atom,&
    mol_to_last_atom,mol_to_nelectrons,mol_to_nbasis,error)
    
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: atom_to_mol, &
                                                mol_to_first_atom, &
                                                mol_to_last_atom
    INTEGER, DIMENSION(:, :), INTENT(OUT), &
      OPTIONAL                               :: mol_to_nelectrons
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: mol_to_nbasis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_molecule_set_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, iatom, imol, &
                                                last_atom, nbasis, nelec
    LOGICAL                                  :: failure
    TYPE(molecule_kind_type), POINTER        :: my_molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    failure=.FALSE.

    DO imol = 1,SIZE(molecule_set)

       molecule => molecule_set(imol)
       CALL get_molecule( molecule=molecule,molecule_kind=my_molecule_kind,&
                           first_atom = first_atom,last_atom = last_atom )

       IF (PRESENT(mol_to_nelectrons)) THEN
         CALL get_molecule_kind(my_molecule_kind,nelectron=nelec)
         IF (MOD(nelec,2).ne.0) THEN
           CALL stop_program(routineN,moduleN,__LINE__,&
                        "ODD NUMBER OF ELECTRONS ON A MOLECULE")
         ELSE
           !ZZZlater number of beta? 
           mol_to_nelectrons(imol,1)=nelec/2
         ENDIF
       ENDIF

       IF (PRESENT(mol_to_nbasis)) THEN
         CALL get_molecule_kind(my_molecule_kind,nsgf=nbasis)
         mol_to_nbasis(imol)=nbasis
       ENDIF

       IF (PRESENT(mol_to_first_atom)) THEN
         mol_to_first_atom(imol)=first_atom
       ENDIF

       IF (PRESENT(mol_to_last_atom)) THEN
         mol_to_last_atom(imol)=last_atom
       ENDIF

       IF (PRESENT(atom_to_mol)) THEN
         DO iatom = first_atom, last_atom
            atom_to_mol(iatom)=imol
         ENDDO  ! iatom
       ENDIF

    END DO  ! imol

  END SUBROUTINE get_molecule_set_info 

END MODULE molecule_types_new
