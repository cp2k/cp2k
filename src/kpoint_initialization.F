!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
!!****** cp2k/kpoint_initialization [1.0] *
!!
!!   NAME
!!     kpoint_initialization
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
!
    MODULE kpoint_initialization
!
!------------------------------------------------------------------------------!

      USE kinds, ONLY : dbl

      USE brillouin, ONLY : kpoint_type
      USE global_types, ONLY : global_environment_type
      USE k290, ONLY : set_k290, set_k290_atoms, kp_sym_gen
      USE particle_types, ONLY : particle_type
      USE termination, ONLY : stop_memory
!
      IMPLICIT NONE
!
      PRIVATE
!
      PUBLIC :: initialize_kpoints

!!*****
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!!****** kpoint_initialization/initialize_kpoints [1.0] *
!!
!!   NAME
!!     initialize_kpoints
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

      SUBROUTINE initialize_kpoints(inpar,kp,symm,hmat,part)

        IMPLICIT NONE

        TYPE (global_environment_type) :: inpar
        TYPE (kpoint_type) :: kp
        LOGICAL, INTENT (IN) :: symm
        REAL (dbl), INTENT (IN) :: hmat(:,:)
        TYPE (particle_type), INTENT (IN) :: part(:)

        REAL (dbl), ALLOCATABLE :: coor(:,:)
        INTEGER, ALLOCATABLE :: types(:)
        INTEGER :: i, nat, isos

        IF (symm .OR. kp%scheme=='MONKHORST-PACK' .OR. kp%scheme=='MACDONALD') &
            THEN
          CALL set_k290(hmat,kp%nk,shift=kp%shift,symm=kp%symmetry, &
            unit=inpar%scr)
          nat = size(part)
          ALLOCATE (coor(3,nat),STAT=isos)
          IF (isos/=0) CALL stop_memory('initialize_kpoints','coor',3*nat)
          ALLOCATE (types(nat),STAT=isos)
          IF (isos/=0) CALL stop_memory('initialize_kpoints','types',nat)
          DO i = 1, nat
            coor(:,i) = part(i) %r
            types(i) = part(i) %prop%ptype
          END DO
          CALL set_k290_atoms(coor,types)
          DEALLOCATE (coor,STAT=isos)
          IF (isos/=0) CALL stop_memory('initialize_kpoints','coor')
          DEALLOCATE (types,STAT=isos)
          IF (isos/=0) CALL stop_memory('initialize_kpoints','types')
          CALL kp_sym_gen
        END IF

      END SUBROUTINE initialize_kpoints

!!*****
!------------------------------------------------------------------------------!

    END MODULE kpoint_initialization

!------------------------------------------------------------------------------!
