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

!!****t* kg_gpw_fm_mol_types/mol_mo_set_type [1.0] *
!!
!!   NAME
!!     mol_mo_set_type
!!
!!   FUNCTION
!!     contains MO for one molecule, and relative informations
!!     used in KG_GPW method to separate the block diagonal full matrix of MO
!!     in molecular blocks. Here also eigenvalues and occupation numbers 
!!     relative to one molecule are stored
!!
!!   NOTES
!!      -
!!   
!!   ATTRIBUTES
!!     - mo :  the actual MO coefficients as a matrix for one molecule
!!     - nmo :  number of molecular orbitals (# cols in mo_coeff)
!!     - nao : number of atomic orbitals (# rows in mo_coeff)
!!     - eigenvalues : eigenvalues  of the nmo states 
!!     - occupation_numbers : 
!!     - maxocc : maximum allowed occupation number of an MO (1-2)
!!     - nelectron : number of electrons (taking occupation into account)
!!     - homo:  highest non-zero occupied orbital
!!     - lfomo: lowest non maxocc occupied orbital (e.g. fractional or zero)
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     11.2004 created [MI]
!!
!!*** *********************************************************************

 TYPE mol_mo_set_type
    REAL(dp), DIMENSION(:,:),  POINTER        :: mo
    INTEGER                                   :: nmo
    INTEGER                                   :: nao
    REAL(KIND = dp), DIMENSION(:), POINTER    :: eigenvalues,occupation_numbers
    REAL(KIND = dp)                           :: maxocc
    INTEGER                                   :: nelectron
    INTEGER                                   :: homo
    INTEGER                                   :: lfomo
  END TYPE mol_mo_set_type

!!****t* kg_gpw_fm_mol_types/mol_mo_set_p_type [1.0] *
!!
!!   NAME
!!     mol_mo_set_p_type
!!
!!   FUNCTION
!!     contains array of mol_mo_set_type, needed if nspins =/ 1
!!
!!   NOTES
!!      -
!!   
!!   ATTRIBUTES
!!      mo_set
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     11.2004 created [MI]
!!
!!*** *********************************************************************

  TYPE mol_mo_set_p_type
    TYPE(mol_mo_set_type), POINTER :: mo_set
  END TYPE mol_mo_set_p_type

!!****t* kg_gpw_fm_mol_types/fm_mol_blocks_type [1.0] *
!!
!!   NAME
!!     fm_mol_blocks_type
!!
!!   FUNCTION
!!     contains the info relative to one molecule, which are  
!!     needed to create the corresponding MO,  
!!     in a block-like form. The index and kind of the  atoms are also stored
!!
!!   NOTES
!!      -
!!   
!!   ATTRIBUTES
!!     - mos:  MO orbitals
!!     - ortho : scratch space for the orthogonalization
!!     - work  : scratch space for operation with the MO matrix
!!     - imol_global: index of the molecule
!!     - index_atom: index of atoms belonging to the molecule
!!     - index_kind: index of ithe kinds of the atoms belonging to the molecule
!!     - ifirst_ao: index of the first ao for each atom
!!     - ilast_ao: index of the last ao for each atom
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     11.2004 created [MI]
!!
!!*** *********************************************************************

  TYPE fm_mol_blocks_type
    TYPE(mol_mo_set_p_type), DIMENSION(:), &
       POINTER                                :: mos
    REAL(dp), DIMENSION(:,:),  POINTER        :: ortho
    REAL(dp), DIMENSION(:,:),  POINTER        :: work
    INTEGER                                   :: imol_global
    INTEGER, DIMENSION(:), POINTER            :: index_atom
    INTEGER, DIMENSION(:), POINTER            :: index_kind
    INTEGER, DIMENSION(:), POINTER            :: ifirst_ao
    INTEGER, DIMENSION(:), POINTER            :: ilast_ao
  END TYPE fm_mol_blocks_type

!!****t* kg_gpw_fm_mol_types/kg_fm_mol_set_type [1.0] *
!!
!!   NAME
!!     kg_fm_mol_set_type
!!
!!   FUNCTION
!!     contains the info relative to one molecule kind, which are  
!!     needed to create the corresponding MO,  
!!     in a block-like form. 
!!     The set of molecules that are handled locally by one processor
!!     and the inro related to them are also given here
!!
!!   NOTES
!!      -
!!   
!!   ATTRIBUTES
!!     - multiplicity
!!     - n_ao
!!     - n_mo
!!     - natom
!!     - nelectron_spin
!!     - nmol_kind_global
!!     - nmolecule_local
!!     - maxocc
!!     - fm_mol_blocks
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     11.2004 created [MI]
!!
!!*** *********************************************************************

  TYPE kg_fm_mol_set_type
    INTEGER  :: multiplicity, &
                n_ao, n_mo(2), natom, &
                nelectron_spin(2),&
                nmol_kind_global,&
                nmolecule_local
    REAL(dp) :: maxocc
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
       POINTER                                :: fm_mol_blocks
  END TYPE kg_fm_mol_set_type

!!****t* kg_gpw_fm_mol_types/kg_fm_p_type [1.0] *
!!
!!   NAME
!!     kg_fm_p_type
!!
!!   FUNCTION
!!      global type that drives to  fm_mol_set_type for each molecule kind
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     11.2004 created [MI]
!!
!!*** *********************************************************************

  TYPE kg_fm_p_type
    TYPE(kg_fm_mol_set_type), DIMENSION(:),& 
       POINTER                                :: kg_fm_mol_set
    INTEGER  :: nao_global, nelectron_global, nmolecule_global, ref_count
    INTEGER  :: nao_max, nmo_max
    REAL(dp) :: maxocc_global(2)
  END TYPE kg_fm_p_type

! *** Public subroutines ***

  PUBLIC ::  allocate_fm_mol_blocks, allocate_mol_mo_set, &
             get_fm_mol_block, get_kg_fm_mol_set,&
             get_mol_mo_set,init_mol_mo_set, &
             kg_fm_mol_set_create, kg_fm_mol_set_release, &
             kg_fm_mol_set_retain, kg_gpw_fm_mol_to_fm_mol, & 
             set_kg_fm_mol_set, set_mol_mo_set

! *** Public data types ***

  PUBLIC :: fm_mol_blocks_type, kg_fm_mol_set_distroy, kg_fm_mol_set_type, &
            kg_fm_p_type, mol_mo_set_type, mol_mo_set_p_type 

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE  allocate_fm_mol_blocks(kg_fm_mol_set,nmol,error)

    TYPE(kg_fm_mol_set_type)                  :: kg_fm_mol_set
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
        NULLIFY(fm_mol_block%work)
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
    INTEGER, INTENT(IN)                       :: n_ao,nmo,nelectron
    REAL(dp), INTENT(IN)                      :: maxocc
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

  SUBROUTINE get_mol_mo_set(mo_set,homo,lfomo,maxocc,nmo,nelectron,&
                            mo,eigenvalues,occupation_numbers)

    TYPE(mol_mo_set_type), POINTER            :: mo_set
    INTEGER, INTENT(OUT), OPTIONAL            :: homo,lfomo,nmo,nelectron
    REAL(dp), INTENT(OUT), OPTIONAL           :: maxocc
    REAL(dp), DIMENSION(:,:), POINTER,&
      OPTIONAL                                :: mo
    REAL(dp), DIMENSION(:), POINTER, OPTIONAL :: eigenvalues,&
                                                 occupation_numbers
    IF(PRESENT(mo)) mo  => mo_set%mo
    IF(PRESENT(occupation_numbers)) &
        occupation_numbers  => mo_set%occupation_numbers
    IF(PRESENT(eigenvalues)) eigenvalues => mo_set%eigenvalues
    IF(PRESENT(maxocc)) maxocc  = mo_set%maxocc
    IF(PRESENT(homo)) homo = mo_set%homo
    IF(PRESENT(lfomo)) lfomo = mo_set%lfomo
    IF(PRESENT(nmo)) nmo = mo_set%nmo
    IF(PRESENT(nelectron)) nelectron = mo_set%nelectron

  END SUBROUTINE get_mol_mo_set 

! *****************************************************************************
 
  SUBROUTINE set_mol_mo_set(mo_set,homo,lfomo,maxocc,nmo,nelectron,&
                            mo,eigenvalues,occupation_numbers)

    TYPE(mol_mo_set_type), POINTER            :: mo_set
    INTEGER, INTENT(IN), OPTIONAL             :: homo,lfomo,nmo,nelectron
    REAL(dp), INTENT(In), OPTIONAL            :: maxocc
    REAL(dp), DIMENSION(:,:), POINTER,&
      OPTIONAL                                :: mo
    REAL(dp), DIMENSION(:), POINTER, OPTIONAL :: eigenvalues,&
                                                 occupation_numbers
    CHARACTER(len=*), PARAMETER :: routineN = 'set_mol_mo_set'

    INTEGER                     :: istat

    IF(PRESENT(maxocc)) mo_set%maxocc = maxocc 
    IF(PRESENT(homo)) mo_set%homo = homo
    IF(PRESENT(lfomo)) mo_set%lfomo = lfomo
    IF(PRESENT(nmo)) mo_set%nmo = nmo
    IF(PRESENT(nelectron)) mo_set%nelectron = nelectron
    IF(PRESENT(eigenvalues)) THEN
      IF(ASSOCIATED(mo_set%eigenvalues)) THEN
        DEALLOCATE(mo_set%eigenvalues, STAT =istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,"dealloc. mo_set%eigenvalues")
        END IF
      END IF
      mo_set%eigenvalues => eigenvalues
    END IF
    IF(PRESENT(occupation_numbers)) THEN
      IF(ASSOCIATED(mo_set%occupation_numbers)) THEN
        DEALLOCATE(mo_set%occupation_numbers, STAT =istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,"dealloc. mo_set%occupation_numbers")
        END IF
      END IF
      mo_set%occupation_numbers => occupation_numbers
    END IF

  END SUBROUTINE set_mol_mo_set

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
    LOGICAL                                  :: failure, mo_uocc

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(mo_set),cp_failure_level,routineP,error,failure)

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
    mo_uocc = .TRUE.
    CALL set_mo_occupation(mo_set%occupation_numbers,mo_set%eigenvalues,&
                           homo,lfomo,maxocc,nelectron,mo_uocc)

    mo_set%homo = homo
    mo_set%lfomo = lfomo

    ALLOCATE (mo_set%mo(n_ao,nmo),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,"mo_set%mo",n_ao*nmo*dp_size)
    END IF

  END SUBROUTINE init_mol_mo_set

! *****************************************************************************
!***************************************************************************
!!****f* kg_gpw_fm_mol_types/kg_gpw_fm_mol_to_fm_mol
!!
!!   NAME
!!    kg_gpw_fm_mol_to_fm_mol
!!
!!   FUNCTION
!!     Copy the full set of molecular mos into a new set
!!     This is needed to store the mos when some sort of extrapolation
!!     is used as an initial guess for the next set of MD or GEO_OPT
!!     
!!   NOTE
!!
!!   INPUTS
!!     - fm_mol_set_a : contains the actual  mos 
!!     - fm_mol_set_b : where the mos are stored
!!     - error
!!
!!   AURHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************

  SUBROUTINE kg_gpw_fm_mol_to_fm_mol(fm_mol_set_a,fm_mol_set_b,error)

    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER                                 :: fm_mol_set_a
    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER                                 :: fm_mol_set_b
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_gpw_fm_mol_to_fm_mol',&
      routineP = module_name//':'//routineN

    INTEGER :: i, imol, ispin, istat, n_ao, n_mo(2), natom, nmol, nmolecule_kind, nspins
    LOGICAL :: failure

    TYPE(fm_mol_blocks_type), POINTER         :: mol_block_a,&
                                                 mol_block_b
    TYPE(kg_fm_mol_set_type), POINTER         :: fm_mol

!
!   ---------------------------------------------------------------------------
!

    nmolecule_kind=SIZE(fm_mol_set_a,1)
    CPPostcondition(SIZE(fm_mol_set_b,1)==nmolecule_kind,cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN
      DO i = 1,nmolecule_kind
        fm_mol_set_b(i)%maxocc = fm_mol_set_a(i)%maxocc 
        fm_mol_set_b(i)%multiplicity = fm_mol_set_a(i)%multiplicity
        fm_mol_set_b(i)%natom = fm_mol_set_a(i)%natom
        fm_mol_set_b(i)%n_ao = fm_mol_set_a(i)%n_ao
        fm_mol_set_b(i)%n_mo = fm_mol_set_a(i)%n_mo 
        fm_mol_set_b(i)%nelectron_spin = fm_mol_set_a(i)%nelectron_spin
        fm_mol_set_b(i)%nmol_kind_global = fm_mol_set_a(i)%nmol_kind_global 
        fm_mol_set_b(i)%nmolecule_local = fm_mol_set_a(i)%nmolecule_local
        n_ao = fm_mol_set_a(i)%n_ao
        n_mo = fm_mol_set_a(i)%n_mo
        nmol = SIZE(fm_mol_set_a(i)%fm_mol_blocks,1)
        natom = fm_mol_set_a(i)%natom

        IF(.NOT. ASSOCIATED(fm_mol_set_b(i)%fm_mol_blocks)) THEN
          CALL allocate_fm_mol_blocks(fm_mol_set_b(i),nmol,error)
        END IF

        DO imol = 1,nmol
          mol_block_a => fm_mol_set_a(i)%fm_mol_blocks(imol)
          nspins = SIZE(mol_block_a%mos,1)
          mol_block_b => fm_mol_set_b(i)%fm_mol_blocks(imol)

          IF(.NOT. ASSOCIATED(mol_block_b%mos)) THEN
            ALLOCATE (fm_mol_set_b(i)%fm_mol_blocks(imol)%mos(nspins), STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            DO ispin = 1,nspins
             NULLIFY(fm_mol_set_b(i)%fm_mol_blocks(imol)%mos(ispin)%mo_set)
            END DO 
          END IF
          DO ispin = 1,nspins
            IF(.NOT. ASSOCIATED(mol_block_b%mos(ispin)%mo_set)) THEN
              ALLOCATE (fm_mol_set_b(i)%fm_mol_blocks(imol)%mos(ispin)%mo_set,STAT=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
              NULLIFY(fm_mol_set_b(i)%fm_mol_blocks(imol)%mos(ispin)%mo_set%mo)
              NULLIFY(fm_mol_set_b(i)%fm_mol_blocks(imol)%mos(ispin)%mo_set%eigenvalues)
              NULLIFY(fm_mol_set_b(i)%fm_mol_blocks(imol)%mos(ispin)%mo_set%occupation_numbers)
            END IF

            mol_block_b%mos(ispin)%mo_set%nmo =  n_mo(ispin)
            mol_block_b%mos(ispin)%mo_set%nao = n_ao
            mol_block_b%mos(ispin)%mo_set%maxocc =&
                        mol_block_a%mos(ispin)%mo_set%maxocc
            mol_block_b%mos(ispin)%mo_set%nelectron =&
                        mol_block_a%mos(ispin)%mo_set%nelectron
            mol_block_b%mos(ispin)%mo_set%homo = &
                        mol_block_a%mos(ispin)%mo_set%homo 
            mol_block_b%mos(ispin)%mo_set%lfomo = &
                        mol_block_a%mos(ispin)%mo_set%lfomo 

            IF(.NOT. ASSOCIATED(mol_block_b%mos(ispin)%mo_set%mo)) THEN
              ALLOCATE (fm_mol_set_b(i)%fm_mol_blocks(imol)%mos(ispin)%&
                        mo_set%mo(n_ao,n_mo(ispin)),STAT=istat)
            END IF
            CALL DCOPY(n_ao*n_mo(ispin),mol_block_a%mos(ispin)%mo_set%mo(1,1),1,&
                                        mol_block_b%mos(ispin)%mo_set%mo(1,1),1)
          END DO  ! ispin
          
          IF(.NOT. ASSOCIATED(mol_block_b%ifirst_ao)) THEN
            ALLOCATE(fm_mol_set_b(i)%fm_mol_blocks(imol)%ifirst_ao(natom),STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            ALLOCATE(fm_mol_set_b(i)%fm_mol_blocks(imol)%ilast_ao(natom),STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            ALLOCATE(fm_mol_set_b(i)%fm_mol_blocks(imol)%index_atom(natom),STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF

          fm_mol_set_b(i)%fm_mol_blocks(imol)%ifirst_ao(1:natom) = mol_block_a%ifirst_ao(1:natom)
          fm_mol_set_b(i)%fm_mol_blocks(imol)%ilast_ao(1:natom) = mol_block_a%ilast_ao(1:natom)
          fm_mol_set_b(i)%fm_mol_blocks(imol)%index_atom(1:natom) =&
                                              mol_block_a%index_atom(1:natom)
        END DO  ! imol

      END DO  ! i 
    END IF
  END SUBROUTINE kg_gpw_fm_mol_to_fm_mol

! *****************************************************************************
  SUBROUTINE  kg_fm_mol_set_create(kg_fm,nmolecule_kind,error)

!   Purpose: 

!   History: - Creation (20.11.2004,MI)

!   ***************************************************************************

    TYPE(kg_fm_p_type), POINTER                :: kg_fm
    INTEGER                                    :: nmolecule_kind
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_fm_mol_set_create',&
      routineP = module_name//':'//routineN

    INTEGER                                   :: i, istat
    LOGICAL                                   :: failure

    failure = .FALSE.

    ALLOCATE(kg_fm,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    kg_fm%nelectron_global = 0
    kg_fm%nao_global = 0
    kg_fm%nmolecule_global = 0
    kg_fm%nao_max = 0
    kg_fm%nmo_max = 0
    NULLIFY(kg_fm%kg_fm_mol_set)
    ALLOCATE(kg_fm%kg_fm_mol_set(nmolecule_kind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN
      kg_fm%ref_count = 1
      DO i = 1,nmolecule_kind
        kg_fm%kg_fm_mol_set(i)%maxocc = 0.0_dp
        kg_fm%kg_fm_mol_set(i)%multiplicity = 0
        kg_fm%kg_fm_mol_set(i)%natom = 0
        kg_fm%kg_fm_mol_set(i)%n_ao = 0
        kg_fm%kg_fm_mol_set(i)%n_mo = 0
        kg_fm%kg_fm_mol_set(i)%nelectron_spin = 0
        kg_fm%kg_fm_mol_set(i)%nmol_kind_global = 0
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
          CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
          IF(ASSOCIATED(fm_mol_block%mos(is)%mo_set%eigenvalues)) THEN
            DEALLOCATE(fm_mol_block%mos(is)%mo_set%eigenvalues,STAT=istat)
            CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
          END IF
          IF(ASSOCIATED(fm_mol_block%mos(is)%mo_set%occupation_numbers)) THEN
            DEALLOCATE(fm_mol_block%mos(is)%mo_set%occupation_numbers,STAT=istat)
            CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
          END IF
          DEALLOCATE(fm_mol_block%mos(is)%mo_set,STAT=istat)
          CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        END DO
        DEALLOCATE(fm_mol_block%mos,STAT=istat)
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        IF(ASSOCIATED(fm_mol_block%ortho)) THEN
          DEALLOCATE(fm_mol_block%ortho,STAT=istat)
          CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        END IF
        IF(ASSOCIATED(fm_mol_block%work))  THEN
          DEALLOCATE(fm_mol_block%work,STAT=istat)
          CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        END IF
        IF(ASSOCIATED(fm_mol_block%index_atom)) THEN
          DEALLOCATE(fm_mol_block%index_atom,STAT=istat)
          CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        END IF
        IF(ASSOCIATED(fm_mol_block%index_kind)) THEN
          DEALLOCATE(fm_mol_block%index_kind,STAT=istat)
          CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        END IF
        IF(ASSOCIATED(fm_mol_block%ifirst_ao)) THEN
          DEALLOCATE(fm_mol_block%ifirst_ao,STAT=istat)
          CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        END IF
        IF(ASSOCIATED(fm_mol_block%ilast_ao)) THEN
          DEALLOCATE(fm_mol_block%ilast_ao,STAT=istat)
          CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        END IF
      END DO
      IF(ASSOCIATED(kg_fm_mol_set(i)%fm_mol_blocks)) THEN
        DEALLOCATE(kg_fm_mol_set(i)%fm_mol_blocks, STAT=istat)
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
      END IF
    END DO 
    DEALLOCATE(kg_fm_mol_set, STAT=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

  END SUBROUTINE kg_fm_mol_set_distroy

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

  SUBROUTINE get_kg_fm_mol_set(kg_fm_mol_set, nmol_kind_global,&
                               nmolecule_local, natom, n_ao, n_mo,&
                               nelectron_spin , maxocc, multiplicity,&
                               fm_mol_blocks)

    TYPE(kg_fm_mol_set_type), POINTER         :: kg_fm_mol_set
    INTEGER, INTENT(OUT), OPTIONAL            :: multiplicity,&
                                                nmol_kind_global,&
                                                nmolecule_local,&
                                                natom, n_ao, n_mo(2),&
                                                nelectron_spin(2)
    REAL(dp), INTENT(OUT), OPTIONAL           :: maxocc
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
       POINTER, OPTIONAL                      :: fm_mol_blocks

    IF(PRESENT(maxocc)) maxocc = kg_fm_mol_set%maxocc
    IF(PRESENT(multiplicity)) multiplicity =kg_fm_mol_set%multiplicity
    IF(PRESENT(nmol_kind_global)) nmol_kind_global=kg_fm_mol_set%nmol_kind_global
    IF(PRESENT(nmolecule_local)) nmolecule_local=kg_fm_mol_set%nmolecule_local
    IF(PRESENT(natom)) natom = kg_fm_mol_set%natom
    IF(PRESENT(n_ao)) n_ao = kg_fm_mol_set%n_ao
    IF(PRESENT(n_mo)) n_mo(1:2) = kg_fm_mol_set%n_mo(1:2)
    IF(PRESENT(nelectron_spin)) &
               nelectron_spin(1:2) = kg_fm_mol_set%nelectron_spin(1:2)
    IF(PRESENT(fm_mol_blocks)) fm_mol_blocks => kg_fm_mol_set%fm_mol_blocks

  END SUBROUTINE get_kg_fm_mol_set

!***************************************************************************

  SUBROUTINE get_fm_mol_block(fm_mol_block,mos,ortho,work,&
                              index_atom,index_kind,&
                              ifirst_ao,ilast_ao)

    TYPE(fm_mol_blocks_type), POINTER         :: fm_mol_block
    TYPE(mol_mo_set_p_type), DIMENSION(:), &
      POINTER, OPTIONAL                       :: mos
    REAL(dp), DIMENSION(:,:),  POINTER, &
      OPTIONAL                                :: ortho
    REAL(dp), DIMENSION(:,:),  POINTER, &
      OPTIONAL                                :: work
    INTEGER, DIMENSION(:), POINTER, OPTIONAL  :: index_atom
    INTEGER, DIMENSION(:), POINTER, OPTIONAL  :: index_kind
    INTEGER, DIMENSION(:), POINTER, OPTIONAL  :: ifirst_ao
    INTEGER, DIMENSION(:), POINTER, OPTIONAL  :: ilast_ao

    IF(PRESENT(mos)) mos  => fm_mol_block%mos
    IF(PRESENT(ortho)) ortho  => fm_mol_block%ortho
    IF(PRESENT(work)) work  => fm_mol_block%work
    IF(PRESENT(index_atom)) index_atom  => fm_mol_block%index_atom
    IF(PRESENT(index_kind)) index_kind  => fm_mol_block%index_kind
    IF(PRESENT(ifirst_ao)) ifirst_ao  => fm_mol_block%ifirst_ao
    IF(PRESENT(ilast_ao)) ilast_ao  => fm_mol_block%ilast_ao

  END SUBROUTINE get_fm_mol_block

!***************************************************************************

  SUBROUTINE set_kg_fm_mol_set(kg_fm_mol_set, nmol_kind_global,&
                               nmolecule_local, natom, n_ao, n_mo,&
                               nelectron_spin , maxocc, multiplicity,& 
                               fm_mol_blocks)

    TYPE(kg_fm_mol_set_type), POINTER        :: kg_fm_mol_set

    INTEGER, INTENT(IN), OPTIONAL            :: multiplicity,&
                                                nmol_kind_global,&
                                                nmolecule_local,&
                                                natom, n_ao, n_mo(2),&
                                                nelectron_spin(2)
    REAL(dp), INTENT(IN), OPTIONAL           :: maxocc
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
       POINTER, OPTIONAL                      :: fm_mol_blocks

    CHARACTER(len=*), PARAMETER :: routineN = "set_kg_fm_mol_set", &
      routineP = module_name//':'//routineN

    IF(PRESENT(maxocc)) kg_fm_mol_set%maxocc = maxocc
    IF(PRESENT(multiplicity)) kg_fm_mol_set%multiplicity = multiplicity
    IF(PRESENT(nmol_kind_global)) kg_fm_mol_set%nmol_kind_global = nmol_kind_global
    IF(PRESENT(nmolecule_local)) kg_fm_mol_set%nmolecule_local = nmolecule_local
    IF(PRESENT(natom)) kg_fm_mol_set%natom = natom
    IF(PRESENT(n_ao)) kg_fm_mol_set%n_ao = n_ao
    IF(PRESENT(n_mo)) kg_fm_mol_set%n_mo(1:2) = n_mo(1:2)
    IF(PRESENT(nelectron_spin)) &
               kg_fm_mol_set%nelectron_spin(1:2) = nelectron_spin(1:2)

  END SUBROUTINE set_kg_fm_mol_set

!***************************************************************************

END MODULE  kg_gpw_fm_mol_types
