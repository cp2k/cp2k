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
  USE cube_utils,                      ONLY: destroy_cube_info
  USE kg_gpw_fm_mol_types,             ONLY: fm_mol_blocks_type
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
  USE green_types,                     ONLY: greens_function_type,&
                                             pw_green_release,&
                                             pw_green_retain,&
                                             pw_green_create
  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "kg_gpw_pw_env_types"

! *** Define the derived structure types ***

  TYPE rho_mol_blocks_type
    INTEGER, DIMENSION(:), POINTER              :: index_atom
    REAL(dp), DIMENSION(3)                      :: r0_molecule
    REAL(dp), DIMENSION(:,:), POINTER           :: r_in_molbox
    REAL(dp), DIMENSION(:,:), POINTER           :: r_in_totbox
    TYPE(qs_rho_type), POINTER                  :: rho_mol
  END TYPE rho_mol_blocks_type

  TYPE kg_molbox_env_type
    INTEGER :: natom, nmolecule_local
    REAL(dp) :: rab_max, rad_max
    TYPE(pw_env_type) , POINTER                 :: pw_env_mol 
    TYPE(greens_function_type), POINTER         :: green
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
  PUBLIC :: get_molbox_env, get_rho_mol_block, &
            kg_sub_pw_env_create, kg_sub_pw_env_release, &
            kg_sub_pw_env_retain, rho_mol_blocks_create, &
            set_molbox_env, set_rho_mol_block

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE deallocate_rho_mol_blocks(rho_mol_blocks,error)

    TYPE(rho_mol_blocks_type), DIMENSION(:),POINTER :: rho_mol_blocks
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_rho_mol_blocks',&
      routineP = module_name//':'//routineN

    INTEGER :: imol, istat, nmol
    LOGICAL :: failure
    TYPE(rho_mol_blocks_type), POINTER         :: rho_mol_iblock

    failure = .FALSE.

    IF(ASSOCIATED(rho_mol_blocks)) THEN
      nmol = SIZE(rho_mol_blocks,1)
      DO imol = 1, nmol
        rho_mol_iblock => rho_mol_blocks(imol)
        CALL qs_rho_release(rho_mol_iblock%rho_mol,error=error)
        DEALLOCATE(rho_mol_iblock%r_in_molbox,&
                   rho_mol_iblock%r_in_totbox, STAT = istat)
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        DEALLOCATE(rho_mol_iblock%index_atom,STAT=istat)
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
      END DO  ! imol

      DEALLOCATE(rho_mol_blocks, STAT = istat)
      CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    ELSE
      NULLIFY(rho_mol_blocks)
    END IF

  END SUBROUTINE deallocate_rho_mol_blocks

! *****************************************************************************

  SUBROUTINE get_molbox_env(molbox_env, natom, nmolecule_local,&
                            cell_mol,green,pw_env_mol,rho_mol_blocks)

    TYPE(kg_molbox_env_type), POINTER          :: molbox_env
    INTEGER, INTENT(OUT), OPTIONAL             :: natom, nmolecule_local
    TYPE(cell_type) , POINTER, OPTIONAL        :: cell_mol
    TYPE(greens_function_type), POINTER ,OPTIONAL     :: green
    TYPE(pw_env_type) , POINTER, OPTIONAL      :: pw_env_mol
    TYPE(rho_mol_blocks_type), DIMENSION(:),&
      POINTER, OPTIONAL                        :: rho_mol_blocks

    IF(PRESENT(natom))  natom = molbox_env%natom
    IF(PRESENT(nmolecule_local)) nmolecule_local =  molbox_env%nmolecule_local
    IF(PRESENT(cell_mol)) cell_mol => molbox_env%cell_mol
    IF(PRESENT(green)) green => molbox_env%green
    IF(PRESENT(pw_env_mol)) pw_env_mol => molbox_env%pw_env_mol
    IF(PRESENT(rho_mol_blocks)) rho_mol_blocks => molbox_env%rho_mol_blocks

  END SUBROUTINE get_molbox_env

! *****************************************************************************

  SUBROUTINE get_rho_mol_block(rho_block,rho_mol,index_atom,r_in_molbox,r_in_totbox)

    TYPE(rho_mol_blocks_type), POINTER          :: rho_block
    TYPE(qs_rho_type), POINTER, OPTIONAL        :: rho_mol
    INTEGER, DIMENSION(:), POINTER, OPTIONAL    :: index_atom
    REAL(dp), DIMENSION(:,:), POINTER, OPTIONAL :: r_in_molbox, r_in_totbox

    IF(PRESENT(rho_mol)) rho_mol => rho_block%rho_mol
    IF(PRESENT(index_atom)) index_atom => rho_block%index_atom
    IF(PRESENT(r_in_molbox)) r_in_molbox => rho_block%r_in_molbox
    IF(PRESENT(r_in_totbox)) r_in_totbox => rho_block%r_in_totbox

  END SUBROUTINE get_rho_mol_block

! *****************************************************************************

  SUBROUTINE set_rho_mol_block(rho_block,r0_molecule,rho_mol,error)

    TYPE(rho_mol_blocks_type), POINTER          :: rho_block
    REAL(dp), INTENT(IN), OPTIONAL              :: r0_molecule(3)
    TYPE(qs_rho_type), POINTER, OPTIONAL        :: rho_mol
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                 :: error

    IF(PRESENT(r0_molecule)) rho_block%r0_molecule(1:3) = r0_molecule(1:3)
    IF (PRESENT(rho_mol)) THEN ! accepts also null pointers !
      IF (ASSOCIATED(rho_mol)) CALL qs_rho_retain(rho_mol,error=error)
      CALL qs_rho_release(rho_block%rho_mol,error=error)
      rho_block%rho_mol => rho_mol
    END IF

  END SUBROUTINE set_rho_mol_block

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
        NULLIFY(kg_sub_pw_env%molbox_env_set(i)%green)
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

    INTEGER                                    :: igrid_level, ikind, imol, istat, n, nmol
    LOGICAL                                    :: failure
    TYPE(rho_mol_blocks_type), POINTER         :: rho_mol_block

    failure=.FALSE.

    n = SIZE(molbox_env_set,1)

    DO ikind = 1,n
      nmol = molbox_env_set(ikind)%nmolecule_local
      
      CALL pw_green_release(molbox_env_set(ikind)%green)
      CALL cell_release(molbox_env_set(ikind)%cell_mol)

      IF (ASSOCIATED(molbox_env_set(ikind)%pw_env_mol%cube_info)) THEN

         NULLIFY(molbox_env_set(ikind)%pw_env_mol%cube_info(1)%lb_cube)
         NULLIFY(molbox_env_set(ikind)%pw_env_mol%cube_info(1)%ub_cube)
         NULLIFY(molbox_env_set(ikind)%pw_env_mol%cube_info(1)%sphere_bounds)
         NULLIFY(molbox_env_set(ikind)%pw_env_mol%cube_info(1)%sphere_bounds_count)

         DO igrid_level=2,SIZE(molbox_env_set(ikind)%pw_env_mol%cube_info)
            CALL destroy_cube_info(molbox_env_set(ikind)%pw_env_mol%cube_info(igrid_level))
         END DO
         DEALLOCATE(molbox_env_set(ikind)%pw_env_mol%cube_info,stat=istat)
         CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
      END IF

      CALL pw_env_release(molbox_env_set(ikind)%pw_env_mol,error=error)

      CALL deallocate_rho_mol_blocks(molbox_env_set(ikind)%rho_mol_blocks,error=error)

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

  SUBROUTINE rho_mol_blocks_create(rho_mol_blocks,nmol,natom,fm_mol_blocks,error)

    TYPE(rho_mol_blocks_type), DIMENSION(:), POINTER :: rho_mol_blocks
    INTEGER, INTENT(IN)                        :: natom, nmol
    TYPE(fm_mol_blocks_type), DIMENSION(:),&
      POINTER                                  :: fm_mol_blocks
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routineN = "rho_mol_blocks_create", &
      routineP = module_name//':'//routineN

    INTEGER :: imol, istat
    LOGICAL :: failure

    failure = .FALSE.

    IF(ASSOCIATED(rho_mol_blocks)) &
        CALL deallocate_rho_mol_blocks(rho_mol_blocks,error=error)

    ALLOCATE(rho_mol_blocks(nmol), STAT = istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO imol = 1, nmol
      NULLIFY(rho_mol_blocks(imol)%rho_mol)
      ALLOCATE(rho_mol_blocks(imol)%index_atom(natom), STAT = istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      rho_mol_blocks(imol)%index_atom(1:natom) = fm_mol_blocks(imol)%index_atom(1:natom)
      rho_mol_blocks(imol)%r0_molecule(1:3) = 0.0_dp
      ALLOCATE(rho_mol_blocks(imol)%r_in_molbox(3,natom), STAT =istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      rho_mol_blocks(imol)%r_in_molbox = 0.0_dp
      ALLOCATE(rho_mol_blocks(imol)%r_in_totbox(3,natom), STAT =istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      rho_mol_blocks(imol)%r_in_totbox = 0.0_dp
    END DO

  END SUBROUTINE rho_mol_blocks_create

! *****************************************************************************

  SUBROUTINE set_molbox_env(molbox_env,cell_mol,green,pw_env_mol,rho_mol_blocks,error)

    TYPE(kg_molbox_env_type), POINTER          :: molbox_env
    TYPE(cell_type) , POINTER, OPTIONAL        :: cell_mol
    TYPE(greens_function_type), POINTER, OPTIONAL      :: green
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

    IF(PRESENT(green)) THEN
      CALL pw_green_retain(green, error=error)
      CALL pw_green_release(molbox_env%green, error=error)
      molbox_env%green => green
    END IF
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
      CALL deallocate_rho_mol_blocks( molbox_env%rho_mol_blocks,error=error)
      molbox_env%rho_mol_blocks => rho_mol_blocks
    END IF

  END SUBROUTINE set_molbox_env

! *****************************************************************************

END MODULE kg_gpw_pw_env_types
