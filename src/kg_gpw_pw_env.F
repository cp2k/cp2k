!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"


!!****** cp2k/kg_gpw_pw_env_types [1.0] *
!!
!!   NAME
!!     kg_gpw_pw_env_types
!!
!!   FUNCTION
!!     Define the  PW environment for each molecule kind
!!     A cubic box is given for each molecule kind, the box should
!!     be smaller than the box containing the whole system. The size
!!     of the box is determined by the typical size of the molecule, 
!!     assuming that intra-molecular bonds are not broken during this
!!     kind of simulations. The local grids are defined in order to be consistent 
!!     with the global grid for the entire system.  
!!     For each molecule a rho structure is defined, which is dimensioned on
!!     the size of the local grid in the molecular-box
!!
!!   AUTHOR
!!     MI (20.12.2004)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE kg_gpw_pw_env_types

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
  USE pw_env_types,                    ONLY: pw_env_release,&
                                             pw_env_retain,&
                                             pw_env_type
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_retain,&
                                             qs_rho_type
  USE simulation_cell,                 ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "kg_gpw_pw_env_types"

! *** Define the derived structure types ***

  TYPE rho_mol_blocks_type
    TYPE(qs_rho_type), POINTER                  :: rho_mol
    REAL(dp), DIMENSION(3)                      :: r0_molecule
    REAL(dp), DIMENSION(:,:), POINTER           :: r_in_molbox
    REAL(dp), DIMENSION(:,:), POINTER           :: r_in_totbox
  END TYPE rho_mol_blocks_type

  TYPE kg_molbox_env_type
    INTEGER :: natom, nmolecule_local
    REAL(dp) :: rab_max, rad_max
    TYPE(pw_env_type) , POINTER                 :: pw_env_mol 
    TYPE(cell_type), POINTER                    :: cell_mol
    TYPE(rho_mol_blocks_type), DIMENSION(:),&
      POINTER                                   :: rho_mol_blocks
  END TYPE kg_molbox_env_type

  TYPE kg_sub_pw_env_type
    TYPE(kg_molbox_env_type), DIMENSION(:), POINTER :: molbox_env_set
    INTEGER                                     :: ref_count
  END TYPE kg_sub_pw_env_type


! *** Public data types ***
  PUBLIC :: kg_molbox_env_type, kg_sub_pw_env_type, rho_mol_blocks_type

! *** Public subroutines ***
  PUBLIC :: get_molbox_env, kg_sub_pw_env_create, kg_sub_pw_env_release, &
            kg_sub_pw_env_retain, set_molbox_env

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_molbox_env(molbox_env,cell_mol,pw_env_mol,rho_mol_blocks)

    TYPE(kg_molbox_env_type), POINTER          :: molbox_env
    TYPE(cell_type) , POINTER, OPTIONAL        :: cell_mol
    TYPE(pw_env_type) , POINTER, OPTIONAL      :: pw_env_mol
    TYPE(rho_mol_blocks_type), DIMENSION(:),&
      POINTER, OPTIONAL                        :: rho_mol_blocks

    IF(PRESENT(cell_mol)) cell_mol => molbox_env%cell_mol
    IF(PRESENT(pw_env_mol)) pw_env_mol => molbox_env%pw_env_mol
    IF(PRESENT(rho_mol_blocks)) rho_mol_blocks => molbox_env%rho_mol_blocks

  END SUBROUTINE get_molbox_env

! *****************************************************************************

  SUBROUTINE kg_sub_pw_env_create(kg_sub_pw_env,nmolecule_kind,error)

    TYPE(kg_sub_pw_env_type), POINTER          :: kg_sub_pw_env
    INTEGER                                    :: nmolecule_kind
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_sub_pw_env_create',&
      routineP = module_name//':'//routineN

    INTEGER                                    :: i, istat
    LOGICAL                                    :: failure

    failure = .FALSE.

    ALLOCATE(kg_sub_pw_env,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(kg_sub_pw_env%molbox_env_set)
    ALLOCATE(kg_sub_pw_env%molbox_env_set(nmolecule_kind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN
      kg_sub_pw_env%ref_count = 1
      DO i = 1,nmolecule_kind
        kg_sub_pw_env%molbox_env_set(i)%natom = 0
        kg_sub_pw_env%molbox_env_set(i)%nmolecule_local = 0
        NULLIFY(kg_sub_pw_env%molbox_env_set(i)%cell_mol)
        NULLIFY(kg_sub_pw_env%molbox_env_set(i)%pw_env_mol)
        NULLIFY(kg_sub_pw_env%molbox_env_set(i)%rho_mol_blocks)
      END DO

    END IF
  
  END SUBROUTINE kg_sub_pw_env_create

! *****************************************************************************

  SUBROUTINE kg_sub_pw_env_destroy(molbox_env_set,error)

    TYPE(kg_molbox_env_type), DIMENSION(:), &
      POINTER                                  :: molbox_env_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_sub_pw_env_destroy',&
      routineP = module_name//':'//routineN

    INTEGER                                    :: ikind, imol, istat, n, nmol
    LOGICAL                                    :: failure
    TYPE(rho_mol_blocks_type), POINTER         :: rho_mol_block

    failure=.FALSE.

    n = SIZE(molbox_env_set,1)

    DO ikind = 1,n
      nmol = molbox_env_set(ikind)%nmolecule_local
 
      CALL cell_release(molbox_env_set(ikind)%cell_mol)
      CALL pw_env_release(molbox_env_set(ikind)%pw_env_mol,error=error)

      DO imol = 1, nmol
        rho_mol_block => molbox_env_set(ikind)%rho_mol_blocks(imol)
        CALL qs_rho_release(rho_mol_block%rho_mol,error=error)
        DEALLOCATE(rho_mol_block%r_in_molbox,&
                   rho_mol_block%r_in_totbox, STAT = istat)
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
      END DO  ! imol

      DEALLOCATE(molbox_env_set(ikind)%rho_mol_blocks, STAT = istat)
      CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    END DO  ! ikind
    DEALLOCATE(molbox_env_set, STAT = istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

  END SUBROUTINE kg_sub_pw_env_destroy

! *****************************************************************************

  SUBROUTINE kg_sub_pw_env_release(kg_sub_pw_env,error)

    TYPE(kg_sub_pw_env_type), POINTER          :: kg_sub_pw_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_sub_pw_env_release',&
      routineP = module_name//':'//routineN

    INTEGER                                    :: istat
    LOGICAL                                    :: failure

    failure=.FALSE.

    IF (ASSOCIATED(kg_sub_pw_env)) THEN

      CPPrecondition(kg_sub_pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
      kg_sub_pw_env%ref_count = kg_sub_pw_env%ref_count - 1
      IF(kg_sub_pw_env%ref_count<1) THEN
        CALL kg_sub_pw_env_destroy(kg_sub_pw_env%molbox_env_set, error=error)
        DEALLOCATE(kg_sub_pw_env, STAT = istat)
        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

    END IF
    NULLIFY(kg_sub_pw_env)

  END SUBROUTINE kg_sub_pw_env_release

! *****************************************************************************

  SUBROUTINE kg_sub_pw_env_retain(kg_sub_pw_env,error)

    TYPE(kg_sub_pw_env_type), POINTER          :: kg_sub_pw_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_sub_pw_env_retain',&
      routineP = module_name//':'//routineN

    LOGICAL                                   :: failure

    failure=.FALSE.
  
    CPPrecondition(ASSOCIATED(kg_sub_pw_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CPPrecondition(kg_sub_pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
      kg_sub_pw_env%ref_count = kg_sub_pw_env%ref_count + 1
    END IF

  END SUBROUTINE kg_sub_pw_env_retain

! *****************************************************************************

  SUBROUTINE set_molbox_env(molbox_env,cell_mol,pw_env_mol,rho_mol_blocks,error)

    TYPE(kg_molbox_env_type), POINTER          :: molbox_env
    TYPE(cell_type) , POINTER, OPTIONAL        :: cell_mol
    TYPE(pw_env_type) , POINTER, OPTIONAL      :: pw_env_mol
    TYPE(rho_mol_blocks_type), DIMENSION(:),&
      POINTER, OPTIONAL                        :: rho_mol_blocks
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routineN = "set_molbox_env", &
      routineP = module_name//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(molbox_env),cp_failure_level,routineP,error,failure)

    IF(PRESENT(cell_mol)) THEN
      CALL cell_retain(cell_mol, error=error)
      CALL cell_release(molbox_env%cell_mol, error=error)
      molbox_env%cell_mol => cell_mol
    END IF
    IF(PRESENT(pw_env_mol)) THEN
      CALL pw_env_retain(pw_env_mol,error=error)
      CALL pw_env_release(molbox_env%pw_env_mol,error=error)
      molbox_env%pw_env_mol => pw_env_mol
    END IF
    IF(PRESENT(rho_mol_blocks)) THEN
      molbox_env%rho_mol_blocks => rho_mol_blocks
    END IF

  END SUBROUTINE set_molbox_env

! *****************************************************************************

END MODULE kg_gpw_pw_env_types
