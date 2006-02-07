!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002-2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


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
MODULE kg_diis_types
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_diis_types'

! Publilc Types
  PUBLIC :: kg_diis_buffer_type, mol_diis_block_type, mol_diis_buffer_set_type

! Publilc Procedures
  PUBLIC :: allocate_diis_block, get_mol_diis_set, init_mol_diis_set,&
            kg_diis_b_retain, kg_diis_b_release, set_mol_diis_set

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

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_diis_b_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(diis_buffer),cp_failure_level,routineP,error,failure)
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
!!*************************************************************************
  SUBROUTINE kg_diis_b_release(diis_buffer,error)

    TYPE(kg_diis_buffer_type), POINTER       :: diis_buffer
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_diis_b_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, nbuffer, &
                                                nmol_kind, nmol_local, nspin
    LOGICAL                                  :: failure
    TYPE(mol_diis_block_type), &
      DIMENSION(:), POINTER                  :: mol_diis_block

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

    TYPE(mol_diis_buffer_set_type), POINTER  :: mol_diis_buffer_set
    TYPE(mol_diis_block_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: mol_diis_block
    INTEGER, INTENT(IN), OPTIONAL            :: imol
    TYPE(mol_diis_block_type), OPTIONAL, &
      POINTER                                :: mol_diis_iblock
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER              :: routineN = 'get_mol_diis_set'

    INTEGER                                  :: nmol_local

    IF(PRESENT(mol_diis_block)) mol_diis_block => mol_diis_buffer_set%mol_diis_block
    IF(PRESENT(imol)) THEN
      nmol_local =  mol_diis_buffer_set%nmol_local
      CPPostconditionNoFail(imol<=nmol_local,cp_failure_level,routineN,error)
      CPPostconditionNoFail(ASSOCIATED(mol_diis_buffer_set%mol_diis_block),cp_failure_level,routineN,error)
      IF(PRESENT(mol_diis_iblock)) mol_diis_iblock => &
                            mol_diis_buffer_set%mol_diis_block(imol)
    END IF

  END SUBROUTINE get_mol_diis_set

!***************************************************************************
  SUBROUTINE allocate_diis_block(mol_diis_buffer_set,nmol_global,nmol_local,error)

    TYPE(mol_diis_buffer_set_type), POINTER  :: mol_diis_buffer_set
    INTEGER, INTENT(IN)                      :: nmol_global, nmol_local
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_diis_block'

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(mol_diis_buffer_set),cp_failure_level,routineN,error,failure)

    IF(ASSOCIATED(mol_diis_buffer_set%mol_diis_block)) &
       CALL distroy_mol_diis_block(mol_diis_buffer_set%mol_diis_block)

    ALLOCATE(mol_diis_buffer_set%mol_diis_block(nmol_local), STAT = istat)
    CPPostcondition(istat==0,cp_failure_level,routineN,error,failure)

    mol_diis_buffer_set%nmol_global = nmol_global
    mol_diis_buffer_set%nmol_local = nmol_local
    
    CALL init_mol_diis_block(mol_diis_buffer_set,nmol_local,error)

  END SUBROUTINE allocate_diis_block

!***************************************************************************

  SUBROUTINE init_mol_diis_set(mol_diis_buffer_set,nmol_kind,error)

    TYPE(mol_diis_buffer_set_type), &
      DIMENSION(:), POINTER                  :: mol_diis_buffer_set
    INTEGER, INTENT(IN)                      :: nmol_kind
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER              :: routineN = 'init_mol_diis_set'

    INTEGER                                  :: ikind
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(mol_diis_buffer_set),cp_failure_level,routineN,error,failure)
    DO ikind = 1,nmol_kind

      NULLIFY(mol_diis_buffer_set(ikind)%mol_diis_block)

    END DO

  END SUBROUTINE init_mol_diis_set

!***************************************************************************

  SUBROUTINE set_mol_diis_set(mol_diis_buffer_set,nmol_global,nmol_local,&
                              mol_diis_block)

    TYPE(mol_diis_buffer_set_type), POINTER  :: mol_diis_buffer_set
    INTEGER, INTENT(IN), OPTIONAL            :: nmol_global, nmol_local
    TYPE(mol_diis_block_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: mol_diis_block

    CHARACTER(len=*), PARAMETER              :: routineN = 'set_mol_diis_set'

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

    TYPE(mol_diis_block_type), &
      DIMENSION(:), POINTER                  :: mol_diis_block

    CHARACTER(len=*), PARAMETER :: routineN = 'distroy_mol_diis_block'

    INTEGER                                  :: ibuffer, imol, ispin, istat, &
                                                nbuffer, nmol_local, nspins

    nmol_local = SIZE(mol_diis_block,1)

    DO imol = 1,nmol_local
      IF(ASSOCIATED(mol_diis_block(imol)%ks_buffer)) THEN
        nbuffer = SIZE(mol_diis_block(imol)%ks_buffer,1)
        nspins = SIZE(mol_diis_block(imol)%ks_buffer,2)
        DO ispin = 1,nspins
          DO ibuffer = 1,nbuffer
            DEALLOCATE(mol_diis_block(imol)%ks_buffer(ibuffer,ispin)%array, &
                       STAT = istat)
            IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                               "deall. mol_diis_block%ks_buffer%array")
          END DO
        END DO
        DEALLOCATE(mol_diis_block(imol)%ks_buffer, STAT = istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                               "deall. mol_diis_block%ks_buffer%")
      END IF
      IF(ASSOCIATED(mol_diis_block(imol)%er_buffer)) THEN
        nbuffer = SIZE(mol_diis_block(imol)%er_buffer,1)
        nspins = SIZE(mol_diis_block(imol)%er_buffer,2)
        DO ispin = 1,nspins
          DO ibuffer = 1,nbuffer
            DEALLOCATE(mol_diis_block(imol)%er_buffer(ibuffer,ispin)%array, &
                       STAT = istat)
            IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                               "deall. mol_diis_block%er_buffer%array")
          END DO
        END DO
        DEALLOCATE(mol_diis_block(imol)%er_buffer, STAT = istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                               "deall. mol_diis_block%er_buffer%")
      END IF
      IF(ASSOCIATED(mol_diis_block(imol)%b_buffer)) THEN
        DEALLOCATE(mol_diis_block(imol)%b_buffer, STAT = istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                               "deall. mol_diis_block%b_buffer%")
      END IF
    END DO
    DEALLOCATE(mol_diis_block, STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                               "mol_diis_block")

  END SUBROUTINE distroy_mol_diis_block

!***************************************************************************

  SUBROUTINE init_mol_diis_block(mol_diis_buffer_set,nmol_local,error)

    TYPE(mol_diis_buffer_set_type), POINTER  :: mol_diis_buffer_set
    INTEGER, INTENT(IN)                      :: nmol_local
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_mol_diis_block'

    INTEGER                                  :: imol, nsize
    LOGICAL                                  :: failure
    TYPE(mol_diis_block_type), &
      DIMENSION(:), POINTER                  :: mol_diis_block

    failure = .FALSE.

    mol_diis_block => mol_diis_buffer_set%mol_diis_block
    nsize = SIZE(mol_diis_block,1)
    CPPostcondition(nsize==nmol_local,cp_failure_level,routineN,error,failure)

    DO imol = 1, nmol_local
      NULLIFY(mol_diis_block(imol)%ks_buffer)
      NULLIFY(mol_diis_block(imol)%er_buffer)
      NULLIFY(mol_diis_block(imol)%b_buffer)
      mol_diis_block(imol)%mol_error_max = 0.0_dp
    END DO

  END SUBROUTINE init_mol_diis_block

!***************************************************************************

!***************************************************************************
END MODULE kg_diis_types
!***************************************************************************
