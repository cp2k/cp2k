!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/kg_gpw_fm_mol_methods [1.0] *
!!
!!   NAME
!!     kg_gpw_fm_mol_methods
!!
!!   FUNCTION
!!     Contruct, initialize and use the full matrix blocks for the KG_GPW method 
!!     Using this method all the operations on the full matrixes can be
!!     splitted inoperations on the single molecular blocks, because
!!     each molecule is treated as independent, i.e. it does not interact
!!     directly with the other molecules 
!!
!!     It should contain also the routines for the :
!!                           initialization of mos (guess or restart)
!!                           orthogonalization
!!                           diagonalization
!!                           transfer of data to and from the sparse matrix
!!     
!!   AUTHOR
!!     MI (20.11.2004)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE kg_gpw_fm_mol_methods

  USE atomic_kind_types,               ONLY: get_atomic_kind
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
  USE distribution_1d_types,           ONLY: distribution_1d_type 
  USE kg_gpw_fm_mol_types,             ONLY: allocate_fm_mol_blocks,&
                                             allocate_mol_mo_set,&
                                             fm_mol_blocks_type,&
                                             init_mol_mo_set,&
                                             kg_fm_mol_set_type, set_kg_fm_mol_set
  USE kinds,                           ONLY: dp,dp_size,int_size
  USE molecule_kind_types,             ONLY: get_molecule_kind,molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "kg_gpw_fm_mol_methods"

! *** Define the derived structure types ***

! *** Public subroutines ***

  PUBLIC ::  build_local_fm_mol

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE  build_local_fm_mol(fm_mol_set,molecule_kind,imolkind,local_molecules,&
                                 molecule_set,particle_set,nspins,&
                                 nmo_eq_nao,added_mos,error)

!   Purpose: 

!   History: - Creation (20.11.2004,MI)

!   ***************************************************************************

    TYPE(kg_fm_mol_set_type), POINTER        :: fm_mol_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind 
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: imolkind, nspins
    LOGICAL                                  :: nmo_eq_nao
    INTEGER, DIMENSION(2)                    :: added_mos
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_local_fm_mol',&
      routineP = module_name//':'//routineN

    TYPE(molecule_type), POINTER              :: molecule
    TYPE(fm_mol_blocks_type), POINTER         :: fm_mol_block

    INTEGER  :: iao, iat, iatom, iglobal, ikind, imol, ispin, istat, &
               multiplicity, n_ao, n_mo(2), natom, nelectron,&
               nelectron_spin(2), nmol_global, nmol_local, nsgf
    REAL(dp) :: maxocc
    INTEGER, DIMENSION(:), POINTER            :: index_mol_local
    INTEGER, DIMENSION(:), POINTER            :: molecule_list
    LOGICAL                                   :: failure

    failure = .FALSE.

    NULLIFY(molecule_list, index_mol_local, fm_mol_block, molecule)

    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                            molecule_list=molecule_list,&
                            natom=natom,&
                            nelectron=nelectron,nsgf=n_ao)
 
!   If this molecule kind is charged this has to appear here.
!   In this case the number of elctrons is modified and we can calculate
!   the right number of molecular orbitals and te occupation numbers

    IF ((MODULO(nelectron,2) /= 0).AND.(nspins == 1)) THEN
      CALL stop_program(routineN,module_name,__LINE__,&
                   "Use the LSD option for an odd number of electrons.")
    END IF

    nmol_local = local_molecules%n_el(imolkind)
    index_mol_local => local_molecules%list(imolkind)%array

    nmol_global = SIZE(molecule_list,1)

    IF (MODULO(nelectron,2) == 0) THEN
      multiplicity = 1
    ELSE
      multiplicity = 2
    END IF


    IF (nspins == 1) THEN

      maxocc = 2.0_dp
      nelectron_spin(1) = nelectron
      nelectron_spin(2) = 0
      n_mo(1) = nelectron/2
      n_mo(2) = 0

    ELSE

      maxocc=1.0_dp

!     *** the simplist spin distribution is written here. Special cases will
!     *** need additional user input

      IF (MODULO(nelectron + multiplicity - 1,2) /= 0) THEN
        CALL stop_program(routineN,module_name,__LINE__,&
                          "LSD: try to use a different multiplicity.")
      END IF

      nelectron_spin(1) = (nelectron + multiplicity - 1)/2
      nelectron_spin(2) = (nelectron - multiplicity + 1)/2

      IF (nelectron_spin(2) < 0) THEN
        CALL stop_program(routineN,module_name,__LINE__,&
                         "LSD: too few electrons for this multiplicity.")
      END IF

      n_mo(1) = nelectron_spin(1)
      n_mo(2) = nelectron_spin(2)

    END IF

!    *** if additional mos need to be taken along, we do it here ***
     n_mo(1) = n_mo(1) + added_mos(1)
     IF (n_mo(2) /= 0) n_mo(2) = n_mo(2) + added_mos(2)

     IF(nmo_eq_nao) THEN
        n_mo(1) = n_ao
        IF (n_mo(2) /= 0) n_mo(2) = n_ao
     END IF
 
     CALL set_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set,&
                           nmolecule_global=nmol_global,&
                           nmolecule_local=nmol_local,&
                           natom = natom, n_ao=n_ao,&
                           n_mo =n_mo , nelectron_spin=nelectron_spin,& 
                           maxocc=maxocc, multiplicity=multiplicity)

     CALL allocate_fm_mol_blocks(fm_mol_set,nmol_local,error=error)

     DO imol = 1,nmol_local
        fm_mol_block => fm_mol_set%fm_mol_blocks(imol)
        iglobal = index_mol_local(imol)
        molecule => molecule_set(iglobal)
        ALLOCATE(fm_mol_block%mos(nspins),STAT=istat)
        DO ispin = 1,nspins
          NULLIFY(fm_mol_block%mos(ispin)%mo_set)
          CALL allocate_mol_mo_set(fm_mol_block%mos(ispin)%mo_set,n_ao,n_mo(ispin),&
                                   maxocc,nelectron_spin(ispin),error=error)
          CALL init_mol_mo_set(fm_mol_block%mos(ispin)%mo_set,n_ao,&
                               n_mo(ispin),error=error)
        END DO
        ALLOCATE(fm_mol_block%ortho(n_ao,n_ao),STAT=istat)
        IF(istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                        "ortho",n_ao*n_ao*dp_size)

        ALLOCATE(fm_mol_block%index_atom(natom),STAT=istat)
        ALLOCATE(fm_mol_block%index_kind(natom),STAT=istat)
        ALLOCATE(fm_mol_block%ifirst_ao(natom),STAT=istat)
        ALLOCATE(fm_mol_block%ilast_ao(natom),STAT=istat)
        IF(istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
             "index_atom,index_kind,ifirst_ao,ilast_ao",natom*4*int_size)

        iao = 0
        DO iat = 1,natom
          iatom = molecule%first_atom +iat -1
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               kind_number=ikind,nsgf=nsgf)
          fm_mol_block%index_atom(iat) = iatom
          fm_mol_block%index_kind(iat) = ikind
          fm_mol_block%ifirst_ao(iat)  = iao + 1
          iao = iao + nsgf
          fm_mol_block%ilast_ao(iat)  = iao
        END DO

     END DO

  END SUBROUTINE build_local_fm_mol

!***************************************************************************

END MODULE  kg_gpw_fm_mol_methods
