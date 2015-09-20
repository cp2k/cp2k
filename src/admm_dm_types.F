!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types and set/get functions for auxiliary denisty matrix methods
!> \par History
!>      11.2014 created [Ole Schuett]
!> \author Ole Schuett
! *****************************************************************************
MODULE admm_dm_types
  USE cp_control_types,                ONLY: admm_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_release,&
                                             cp_dbcsr_type
  USE input_constants,                 ONLY: do_admm_basis_projection,&
                                             do_admm_purify_mcweeny
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: admm_dm_type, mcweeny_history_type
  PUBLIC :: admm_dm_create, admm_dm_release

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'admm_dm_types'

  TYPE mcweeny_history_type
    TYPE(cp_dbcsr_type)                        :: m
    INTEGER                                    :: count = -1
    TYPE(mcweeny_history_type), POINTER        :: next => Null()
  END TYPE mcweeny_history_type

  TYPE mcweeny_history_p_type
    TYPE(mcweeny_history_type), POINTER        :: p => Null()
  END TYPE mcweeny_history_p_type

  TYPE admm_dm_type
    LOGICAL                                  :: purify = .FALSE.
    INTEGER                                  :: method = -1
    TYPE(cp_dbcsr_type), POINTER             :: matrix_a => Null()
    REAL(KIND=dp)                            :: eps_filter = 1e-20_dp
    INTEGER                                  :: mcweeny_max_steps = 100
    INTEGER, DIMENSION(:,:), POINTER         :: block_map => Null()
    TYPE(mcweeny_history_p_type),&
      DIMENSION(:), POINTER                  :: mcweeny_history => Null()
  END TYPE


  CONTAINS

! *****************************************************************************
!> \brief Create a new admm_dm type
!> \param admm_dm ...
!> \param admm_control ...
!> \param nspins ...
!> \param natoms ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE admm_dm_create(admm_dm, admm_control, nspins, natoms)
    TYPE(admm_dm_type), POINTER              :: admm_dm
    TYPE(admm_control_type), POINTER         :: admm_control
    INTEGER, INTENT(IN)                      :: nspins, natoms

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_dm_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, iblock, j, jatom

    ALLOCATE(admm_dm)
    ! copy settings from admm_control
    admm_dm%purify = (admm_control%purification_method == do_admm_purify_mcweeny)
    admm_dm%method = admm_control%method
    admm_dm%eps_filter = admm_control%eps_filter

    ALLOCATE(admm_dm%mcweeny_history(nspins))

    IF(admm_dm%method /= do_admm_basis_projection) THEN
       ! create block map
       ALLOCATE(admm_dm%block_map(natoms,natoms))
       admm_dm%block_map(:,:) = 0
       DO iblock = 1, SIZE(admm_control%blocks)
         DO i = 1,SIZE(admm_control%blocks(iblock)%list)
           iatom = admm_control%blocks(iblock)%list(i)
           DO j = 1,SIZE(admm_control%blocks(iblock)%list)
             jatom = admm_control%blocks(iblock)%list(j)
             admm_dm%block_map(iatom,jatom) = 1
           END DO
         END DO
       END DO
    ENDIF
  END SUBROUTINE admm_dm_create


! *****************************************************************************
!> \brief Release a admm_dm type
!> \param admm_dm ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE admm_dm_release(admm_dm)
    TYPE(admm_dm_type), POINTER              :: admm_dm

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_dm_release', &
      routineP = moduleN//':'//routineN

    IF(.NOT. ASSOCIATED(admm_dm)) RETURN

    IF(ASSOCIATED(admm_dm%matrix_a)) THEN
      CALL cp_dbcsr_release(admm_dm%matrix_a)
      DEALLOCATE(admm_dm%matrix_a)
    ENDIF

    IF(ASSOCIATED(admm_dm%block_map))&
      DEALLOCATE(admm_dm%block_map)

    DEALLOCATE(admm_dm%mcweeny_history)
    DEALLOCATE(admm_dm)

  END SUBROUTINE admm_dm_release

END MODULE admm_dm_types

