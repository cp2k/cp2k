!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002-2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/kg_diis_types [1.0] *
!!
!!   NAME
!!     kg_diis_types
!!
!!   FUNCTION
!!     buffer for the diis of the scf for the KG_GPW method
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     12.2004 
!!
!!   SOURCE
!****************************************************************************
module kg_diis_types
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: module_name='kg_diis_types'

! Publilc Types
  public :: kg_diis_buffer_type, mol_diis_block_type, mol_diis_buffer_set_type

! Publilc Procedures
  public :: get_mol_diis_set, init_mol_diis_block, kg_diis_b_retain, kg_diis_b_release, set_mol_diis_set

!****************************************************************************

  TYPE mol_diis_block_type
     TYPE(cp_2d_r_p_type), DIMENSION(:,:), POINTER :: ks_buffer    
     TYPE(cp_2d_r_p_type), DIMENSION(:,:), POINTER :: er_buffer    
     REAL(dp) :: mol_error_max
     REAL(dp), DIMENSION(:,:), POINTER             :: b_buffer    
  END TYPE mol_diis_block_type

  TYPE mol_diis_buffer_set_type
    INTEGER :: nmol_global, nmol_local
    TYPE(mol_diis_block_type) , DIMENSION(:), POINTER :: mol_diis_block
  END TYPE mol_diis_buffer_set_type

  TYPE kg_diis_buffer_type
    INTEGER                                          :: id_nr,nbuffer,&
                                                        ncall,nspins,ref_count
    TYPE(mol_diis_buffer_set_type), DIMENSION(:),&
       POINTER                                       :: mol_diis_buffer_set
  END TYPE kg_diis_buffer_type


  TYPE kg_diis_buffer_p_type
    TYPE(kg_diis_buffer_type), DIMENSION(:), POINTER :: kg_diis_buffer
  END TYPE kg_diis_buffer_p_type

!****************************************************************************

CONTAINS

!!****f* kg_diis_types/kg_diis_b_retain [1.0] *
!!
!!   NAME
!!     kg_diis_b_retain
!!
!!   FUNCTION
!!     retains a diis buffer (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - diis_buffer: the buffer to retain
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     12.2003 created 
!!
!!*** **********************************************************************
  SUBROUTINE kg_diis_b_retain(diis_buffer,error)

    TYPE(kg_diis_buffer_type), POINTER       :: diis_buffer
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), parameter :: routineN = 'kg_diis_b_retain', &
      routineP = module_name//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(associated(diis_buffer),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CPPreconditionNoFail(diis_buffer%ref_count>0,cp_failure_level,routineP,error)
      diis_buffer%ref_count=diis_buffer%ref_count+1
    END IF
 
  END SUBROUTINE kg_diis_b_retain

!***************************************************************************

!!****f* kg_diis_types/kg_diis_b_release [1.0] *
!!
!!   NAME
!!     kg_diis_b_release
!!
!!   FUNCTION
!!     releases a diis buffer (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - diis_buffer: the buffer to release
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     12.2003 created 
!!
!!*** **********************************************************************
  SUBROUTINE kg_diis_b_release(diis_buffer,error)

    TYPE(kg_diis_buffer_type), POINTER       :: diis_buffer
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_diis_b_release', &
      routineP = module_name//':'//routineN

    INTEGER                                  :: ibuffer, ikind, imol , ispin, istat, &
                                                nbuffer, nmol_kind, nmol_local,&
                                                nspin
    LOGICAL                                  :: failure
    TYPE(mol_diis_block_type), DIMENSION(:),&
       POINTER                               :: mol_diis_block

    failure=.FALSE.

    IF(ASSOCIATED(diis_buffer)) THEN

      CPPreconditionNoFail(diis_buffer%ref_count>0,cp_failure_level,routineP,error)
      diis_buffer%ref_count=diis_buffer%ref_count-1
      IF (diis_buffer%ref_count<1) THEN

        nbuffer = diis_buffer%nbuffer
        nspin = diis_buffer%nspins
        IF(ASSOCIATED(diis_buffer%mol_diis_buffer_set)) THEN
          nmol_kind = SIZE(diis_buffer%mol_diis_buffer_set,1)

          ! For each molecule kind
          DO ikind = 1,nmol_kind
            IF(ASSOCIATED(diis_buffer%mol_diis_buffer_set(ikind)%mol_diis_block)) THEN
              nmol_local = &
                  SIZE(diis_buffer%mol_diis_buffer_set(ikind)%mol_diis_block,1)

              ! For each molecule
              mol_diis_block=>diis_buffer%mol_diis_buffer_set(ikind)%mol_diis_block
              CALL distroy_mol_diis_block(mol_diis_block)

            END IF 
          END DO  ! ikind

          DEALLOCATE(diis_buffer%mol_diis_buffer_set, STAT = istat)
          CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        END IF
      END IF 
      DEALLOCATE(diis_buffer, STAT = istat)
    END IF

  END SUBROUTINE kg_diis_b_release

!***************************************************************************
  SUBROUTINE get_mol_diis_set(mol_diis_buffer_set,mol_diis_block,imol,&
                              mol_diis_iblock,error)

    TYPE(mol_diis_buffer_set_type), POINTER    :: mol_diis_buffer_set
    TYPE(mol_diis_block_type), DIMENSION(:),&
      POINTER, OPTIONAL                        :: mol_diis_block
    INTEGER, INTENT(IN), OPTIONAL              :: imol
    TYPE(mol_diis_block_type), POINTER, &
      OPTIONAL                                 :: mol_diis_iblock
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                 :: error


    CHARACTER(len=*), PARAMETER :: routine_name ="get_mol_diis_set"

    INTEGER :: nmol_local


    IF(PRESENT(mol_diis_block)) mol_diis_block => mol_diis_buffer_set%mol_diis_block
    IF(PRESENT(imol)) THEN
      nmol_local =  mol_diis_buffer_set%nmol_local
      CPPostconditionNoFail(imol<=nmol_local,cp_failure_level,routine_name,error)
      CPPostconditionNoFail(ASSOCIATED(mol_diis_buffer_set%mol_diis_block),cp_failure_level,routine_name,error)
      IF(PRESENT(mol_diis_iblock)) mol_diis_iblock => &
                            mol_diis_buffer_set%mol_diis_block(imol)
    END IF

  END SUBROUTINE get_mol_diis_set

!***************************************************************************

  SUBROUTINE set_mol_diis_set(mol_diis_buffer_set,nmol_global,nmol_local,&
                              mol_diis_block)

    TYPE(mol_diis_buffer_set_type), POINTER    :: mol_diis_buffer_set
    INTEGER, INTENT(IN), OPTIONAL              :: nmol_global,nmol_local
    TYPE(mol_diis_block_type), DIMENSION(:),&
      POINTER, OPTIONAL                        :: mol_diis_block

    CHARACTER(len=*), PARAMETER :: routine_name ="set_mol_diis_set"

    IF(PRESENT(nmol_global)) mol_diis_buffer_set%nmol_global = nmol_global
    IF(PRESENT(nmol_local)) mol_diis_buffer_set%nmol_local = nmol_local
    IF(PRESENT(mol_diis_block)) THEN
      IF(ASSOCIATED(mol_diis_buffer_set%mol_diis_block )) &
         CALL distroy_mol_diis_block(mol_diis_buffer_set%mol_diis_block)
      mol_diis_buffer_set%mol_diis_block => mol_diis_block
    END IF

  END SUBROUTINE set_mol_diis_set

!***************************************************************************

  SUBROUTINE distroy_mol_diis_block(mol_diis_block)

    TYPE(mol_diis_block_type), DIMENSION(:), POINTER  :: mol_diis_block 

    CHARACTER(len=*), PARAMETER :: routine_name ="distroy_mol_diis_block"

    INTEGER :: ibuffer, imol, ispin, istat, nbuffer, nmol_local, nspins

    nmol_local = SIZE(mol_diis_block,1)

    DO imol = 1,nmol_local
      IF(ASSOCIATED(mol_diis_block(imol)%ks_buffer)) THEN
        nbuffer = SIZE(mol_diis_block(imol)%ks_buffer,1)
        nspins = SIZE(mol_diis_block(imol)%ks_buffer,2)
        DO ispin = 1,nspins
          DO ibuffer = 1,nbuffer
            DEALLOCATE(mol_diis_block(imol)%ks_buffer(ibuffer,ispin)%array, &
                       STAT = istat)
            IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                               "deall. mol_diis_block%ks_buffer%array")
          END DO
        END DO
        DEALLOCATE(mol_diis_block(imol)%ks_buffer, STAT = istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                               "deall. mol_diis_block%ks_buffer%")
      END IF
      IF(ASSOCIATED(mol_diis_block(imol)%er_buffer)) THEN
        nbuffer = SIZE(mol_diis_block(imol)%er_buffer,1)
        nspins = SIZE(mol_diis_block(imol)%er_buffer,2)
        DO ispin = 1,nspins
          DO ibuffer = 1,nbuffer
            DEALLOCATE(mol_diis_block(imol)%er_buffer(ibuffer,ispin)%array, &
                       STAT = istat)
            IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                               "deall. mol_diis_block%er_buffer%array")
          END DO
        END DO
        DEALLOCATE(mol_diis_block(imol)%er_buffer, STAT = istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                               "deall. mol_diis_block%er_buffer%")
      END IF
      IF(ASSOCIATED(mol_diis_block(imol)%b_buffer)) THEN
        DEALLOCATE(mol_diis_block(imol)%b_buffer, STAT = istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                               "deall. mol_diis_block%b_buffer%")
      END IF
    END DO
    DEALLOCATE(mol_diis_block, STAT = istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                               "mol_diis_block")

  END SUBROUTINE distroy_mol_diis_block

!***************************************************************************

  SUBROUTINE init_mol_diis_block(mol_diis_buffer_set,nmol_local,error)

    TYPE(mol_diis_buffer_set_type), POINTER    :: mol_diis_buffer_set
    INTEGER, INTENT(IN)                        :: nmol_local
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routine_name ="set_mol_diis_set"

    INTEGER                                    :: imol, nsize
    LOGICAL                                    :: failure
    TYPE(mol_diis_block_type), DIMENSION(:),&
      POINTER                                  :: mol_diis_block

    failure = .FALSE.

    mol_diis_block => mol_diis_buffer_set%mol_diis_block
    nsize = SIZE(mol_diis_block,1)
    CPPostcondition(nsize/=nmol_local,cp_failure_level,routine_name,error,failure)

    DO imol = 1, nmol_local
      NULLIFY(mol_diis_block(imol)%ks_buffer)
      NULLIFY(mol_diis_block(imol)%er_buffer)
      NULLIFY(mol_diis_block(imol)%b_buffer)
    END DO
    mol_diis_block%mol_error_max = 0.0_dp

  END SUBROUTINE init_mol_diis_block

!***************************************************************************

!***************************************************************************
END MODULE kg_diis_types
!***************************************************************************
