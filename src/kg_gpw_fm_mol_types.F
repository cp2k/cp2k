!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"


!!****** cp2k/kg_gpw_fm_mol_types [1.0] *
!!
!!   NAME
!!     kg_gpw_fm_mol_types
!!
!!   FUNCTION
!!     Define the full matrix blocks for the KG_GPW method 
!!     Using this method all the operations on the full matrixes can be
!!     splitted inoperations on the single molecular blocks, because
!!     each molecule is treated as independent, i.e. it does not interact
!!     directly with the other molecules 
!!
!!   AUTHOR
!!     MI (20.11.2004)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE kg_gpw_fm_mol_types

  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
   
  USE kinds,                           ONLY: dp,dp_size
  USE qs_mo_types,                     ONLY: set_mo_occupation
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "kg_gpw_fm_mol_types"

! *** Define the derived structure types ***

 TYPE mol_mo_set_type
    ! the actual MO coefficients as a matrix
    REAL(dp), DIMENSION(:,:),  POINTER        :: mo
    ! number of molecular orbitals (# cols in mo_coeff)
    INTEGER                                   :: nmo
    ! number of atomic orbitals (# rows in mo_coeff)
    INTEGER                                   :: nao
    ! occupation - eigenvalues  of the nmo states (if eigenstates)
    REAL(KIND = dp), DIMENSION(:), POINTER    :: eigenvalues,occupation_numbers
    ! maximum allowed occupation number of an MO (1-2)
    REAL(KIND = dp)                           :: maxocc
    ! number of electrons (taking occupation into account)
    INTEGER                                   :: nelectron
    ! highest non-zero occupied orbital
    INTEGER                                   :: homo
    ! lowest non maxocc occupied orbital (e.g. fractional or zero)
    INTEGER                                   :: lfomo
  END TYPE mol_mo_set_type

  TYPE mol_mo_set_p_type
    TYPE(mol_mo_set_type), POINTER :: mo_set
  END TYPE mol_mo_set_p_type

  TYPE fm_mol_blocks_type
    TYPE(mol_mo_set_p_type), DIMENSION(:), &
       POINTER                                :: mos
    REAL(dp), DIMENSION(:,:),  POINTER        :: ortho
    INTEGER, DIMENSION(:), POINTER            :: index_atom
    INTEGER, DIMENSION(:), POINTER            :: index_kind
    INTEGER, DIMENSION(:), POINTER            :: ifirst_ao
    INTEGER, DIMENSION(:), POINTER            :: ilast_ao
  END TYPE fm_mol_blocks_type

  TYPE kg_fm_mol_set_type
    INTEGER  :: multiplicity, &
                n_ao, n_mo(2), natom, &
                nelectron_spin(2),&
                nmolecule_global,&
                nmolecule_local
    REAL(dp) :: maxocc
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
       POINTER                                :: fm_mol_blocks
  END TYPE kg_fm_mol_set_type

  TYPE kg_fm_p_type
    TYPE(kg_fm_mol_set_type), DIMENSION(:),& 
       POINTER                                :: kg_fm_mol_set
    INTEGER  :: ref_count
  END TYPE kg_fm_p_type

! *** Public subroutines ***

  PUBLIC ::  allocate_fm_mol_blocks, allocate_mol_mo_set, &
             init_mol_mo_set, kg_fm_mol_set_create, &
             kg_fm_mol_set_release, &
             kg_fm_mol_set_retain, set_kg_fm_mol_set

! *** Public data types ***

  PUBLIC :: fm_mol_blocks_type, kg_fm_mol_set_type, &
            kg_fm_p_type, mol_mo_set_type, mol_mo_set_p_type 

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE  allocate_fm_mol_blocks(kg_fm_mol_set,nmol,error)

    TYPE(kg_fm_mol_set_type), POINTER         :: kg_fm_mol_set
    INTEGER                                   :: nmol
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_fm_mol_blocks',&
      routineP = module_name//':'//routineN

    TYPE(fm_mol_blocks_type), POINTER         :: fm_mol_block
    INTEGER                                   :: imol, istat
    LOGICAL                                   :: failure

    failure = .FALSE.

    IF(nmol > 0) THEN
      ALLOCATE(kg_fm_mol_set%fm_mol_blocks(nmol), STAT=istat) 
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      DO imol = 1,nmol
        fm_mol_block => kg_fm_mol_set%fm_mol_blocks(imol)
        NULLIFY(fm_mol_block%mos)
        NULLIFY(fm_mol_block%ortho)
        NULLIFY(fm_mol_block%index_atom)
        NULLIFY(fm_mol_block%index_kind)
        NULLIFY(fm_mol_block%ifirst_ao)
        NULLIFY(fm_mol_block%ilast_ao)
      END DO 
    ELSE
      NULLIFY(kg_fm_mol_set%fm_mol_blocks)
    END IF

  END SUBROUTINE allocate_fm_mol_blocks

! *****************************************************************************

  SUBROUTINE allocate_mol_mo_set(mo_set,n_ao,nmo,maxocc,nelectron,error)

    TYPE(mol_mo_set_type), POINTER            :: mo_set
    INTEGER                                   :: n_ao,nmo,nelectron
    REAL(dp)                                  :: maxocc
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_mol_mo_set',&
      routineP = module_name//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE (mo_set,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    mo_set%nmo = nmo
    mo_set%nao = n_ao
    mo_set%maxocc = maxocc
    mo_set%nelectron = nelectron
    mo_set%homo = 0
    mo_set%lfomo = 0

    NULLIFY (mo_set%eigenvalues)
    NULLIFY (mo_set%occupation_numbers)
    NULLIFY (mo_set%mo)

  END SUBROUTINE allocate_mol_mo_set

! *****************************************************************************

  SUBROUTINE init_mol_mo_set(mo_set,n_ao,nmo,error)

    TYPE(mol_mo_set_type), POINTER            :: mo_set
    INTEGER                                   :: n_ao,nmo
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_mol_mo_set',&
      routineP = module_name//':'//routineN

    INTEGER                                  :: homo,istat,lfomo,nelectron
    REAL(dp)                                 :: maxocc
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(mo_set),cp_failure_level,routineP,error,failure)

    ALLOCATE(mo_set%eigenvalues(nmo), STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,"mo_set%eigenvalues",nmo*dp_size)
    END IF
    mo_set%eigenvalues(:) = 0.0_dp

    ALLOCATE (mo_set%occupation_numbers(nmo),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,"mo_set%occupation_numbers",nmo*dp_size)
    END IF

    maxocc = mo_set%maxocc
    nelectron = mo_set%nelectron
    CALL set_mo_occupation(mo_set%occupation_numbers,mo_set%eigenvalues,&
                           homo,lfomo,maxocc,nelectron,smear=0.0_dp)
    mo_set%homo = homo
    mo_set%lfomo = lfomo

    ALLOCATE (mo_set%mo(n_ao,nmo),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,"mo_set%mo",n_ao*nmo*dp_size)
    END IF

  END SUBROUTINE init_mol_mo_set

! *****************************************************************************

  SUBROUTINE  kg_fm_mol_set_create(kg_fm,nmolecule_kind,error)

!   Purpose: 

!   History: - Creation (20.11.2004,MI)

!   ***************************************************************************

    TYPE(kg_fm_p_type), POINTER               :: kg_fm
    INTEGER                                   :: nmolecule_kind
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_fm_mol_set_create',&
      routineP = module_name//':'//routineN

    INTEGER                                   :: i, istat
    LOGICAL                                   :: failure

    failure = .FALSE.

    ALLOCATE(kg_fm,STAT=istat)
    NULLIFY(kg_fm%kg_fm_mol_set)
    ALLOCATE(kg_fm%kg_fm_mol_set(nmolecule_kind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN
      kg_fm%ref_count = 1
      DO i = 1,nmolecule_kind
        kg_fm%kg_fm_mol_set(i)%maxocc = 0.0_dp
        kg_fm%kg_fm_mol_set(i)%n_ao = 0
        kg_fm%kg_fm_mol_set(i)%n_mo = 0
        kg_fm%kg_fm_mol_set(i)%nelectron_spin = 0
        kg_fm%kg_fm_mol_set(i)%nmolecule_global = 0
        kg_fm%kg_fm_mol_set(i)%nmolecule_local = 0
        NULLIFY(kg_fm%kg_fm_mol_set(i)%fm_mol_blocks)
      END DO
    END IF

  END SUBROUTINE kg_fm_mol_set_create

!***************************************************************************

  SUBROUTINE kg_fm_mol_set_distroy(kg_fm_mol_set,error)

    TYPE(kg_fm_mol_set_type) , DIMENSION(:), &
       POINTER                                :: kg_fm_mol_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_fm_mol_set_distroy', &
      routineP = module_name//':'//routineN

    TYPE(fm_mol_blocks_type), POINTER         :: fm_mol_block
    INTEGER                                   :: i, is, istat, j , n, nl, nspins
    LOGICAL                                   :: failure

    failure=.FALSE.

    n = SIZE(kg_fm_mol_set)

    DO i = 1,n
      nl = kg_fm_mol_set(i)%nmolecule_local

      DO j = 1,nl
        fm_mol_block => kg_fm_mol_set(i)%fm_mol_blocks(j)
        nspins = SIZE(fm_mol_block%mos)
        DO is = 1,nspins
          DEALLOCATE(fm_mol_block%mos(is)%mo_set%mo,STAT=istat)
          DEALLOCATE(fm_mol_block%mos(is)%mo_set%eigenvalues,STAT=istat)
          DEALLOCATE(fm_mol_block%mos(is)%mo_set%occupation_numbers,STAT=istat)
          DEALLOCATE(fm_mol_block%mos(is)%mo_set,STAT=istat)
        END DO
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        DEALLOCATE(fm_mol_block%ortho,STAT=istat)
        DEALLOCATE(fm_mol_block%index_atom,STAT=istat)
        DEALLOCATE(fm_mol_block%index_kind,STAT=istat)
        DEALLOCATE(fm_mol_block%ifirst_ao,STAT=istat)
        DEALLOCATE(fm_mol_block%ilast_ao,STAT=istat)
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
      END DO

      DEALLOCATE(kg_fm_mol_set(i)%fm_mol_blocks, STAT=istat)
    END DO 
    DEALLOCATE(kg_fm_mol_set, STAT=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
  END SUBROUTINE

!***************************************************************************

  SUBROUTINE kg_fm_mol_set_release(kg_fm,error)

    TYPE(kg_fm_p_type),  POINTER              :: kg_fm
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_fm_mol_set_release', &
      routineP = module_name//':'//routineN

    INTEGER                                   :: istat
    LOGICAL                                   :: failure

    failure=.FALSE.

    IF (ASSOCIATED(kg_fm)) THEN

      CPPrecondition(kg_fm%ref_count>0,cp_failure_level,routineP,error,failure)
      kg_fm%ref_count=kg_fm%ref_count-1
      IF (kg_fm%ref_count<1) THEN
        CALL kg_fm_mol_set_distroy(kg_fm%kg_fm_mol_set,error=error)
        DEALLOCATE(kg_fm,STAT=istat)
        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

    END IF 
    NULLIFY(kg_fm)

  END SUBROUTINE kg_fm_mol_set_release

!***************************************************************************

  SUBROUTINE kg_fm_mol_set_retain(kg_fm,error)

    TYPE(kg_fm_p_type), POINTER               :: kg_fm
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_fm_mol_set_retain', &
      routineP = module_name//':'//routineN

    LOGICAL                                   :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(kg_fm),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(kg_fm%ref_count>0,cp_failure_level,routineP,error,failure)
     kg_fm%ref_count=kg_fm%ref_count+1
  END IF

  END SUBROUTINE kg_fm_mol_set_retain

!***************************************************************************

  SUBROUTINE set_kg_fm_mol_set(kg_fm_mol_set, nmolecule_global,&
                               nmolecule_local, natom, n_ao, n_mo,&
                               nelectron_spin , maxocc, multiplicity)

    TYPE(kg_fm_mol_set_type), POINTER        :: kg_fm_mol_set

    INTEGER, INTENT(IN), OPTIONAL            :: multiplicity,&
                                                nmolecule_global,&
                                                nmolecule_local,&
                                                natom, n_ao, n_mo(2),&
                                                nelectron_spin(2)
    REAL(dp), INTENT(IN), OPTIONAL           :: maxocc

    CHARACTER(len=*), PARAMETER :: routineN = "set_kg_fm_mol_set", &
      routineP = module_name//':'//routineN

    IF(PRESENT(maxocc)) kg_fm_mol_set%maxocc = maxocc
    IF(PRESENT(multiplicity)) kg_fm_mol_set%multiplicity = multiplicity
    IF(PRESENT(nmolecule_global)) kg_fm_mol_set%nmolecule_global = nmolecule_global
    IF(PRESENT(nmolecule_local)) kg_fm_mol_set%nmolecule_local = nmolecule_local
    IF(PRESENT(natom)) kg_fm_mol_set%natom = natom
    IF(PRESENT(n_ao)) kg_fm_mol_set%n_ao = n_ao
    IF(PRESENT(n_mo)) kg_fm_mol_set%n_mo(1:2) = n_mo(1:2)
    IF(PRESENT(nelectron_spin)) &
               kg_fm_mol_set%nelectron_spin(1:2) = nelectron_spin(1:2)

  END SUBROUTINE set_kg_fm_mol_set

!***************************************************************************

END MODULE  kg_gpw_fm_mol_types
