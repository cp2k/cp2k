!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Add CP2K error reporting, new add_force routine [07.2014,JGH]
!> \author MK (03.06.2002)
! *****************************************************************************
MODULE qs_force_types

  !USE cp_control_types,                ONLY: qs_control_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_force_types'
  PRIVATE

  TYPE qs_force_type
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: all_potential,&
                                         core_overlap,&
                                         gth_ppl,&
                                         gth_nlcc,&
                                         gth_ppnl,&
                                         kinetic,&
                                         overlap,&
                                         overlap_admm,&
                                         rho_core,&
                                         rho_elec,&
                                         rho_lri_elec,&
                                         vhxc_atom,&
                                         g0s_Vh_elec,&
                                         repulsive,&
                                         dispersion,&
                                         other,&
                                         ch_pulay,&
                                         fock_4c,&
                                         ehrenfest,&
                                         efield,&
                                         eev,&
                                         mp2_sep,&
                                         mp2_non_sep,&
                                         total
  END TYPE qs_force_type

  PUBLIC :: qs_force_type

  PUBLIC :: allocate_qs_force,&
            add_qs_force,&
            deallocate_qs_force,&
            zero_qs_force

CONTAINS

! *****************************************************************************
!> \brief   Allocate a Quickstep force data structure.
!> \param qs_force ...
!> \param natom_of_kind ...
!> \param error ...
!> \date    05.06.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_qs_force(qs_force,natom_of_kind,error)

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force
    INTEGER, DIMENSION(:), INTENT(IN)        :: natom_of_kind
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_qs_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, n, nkind
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(qs_force)) CALL deallocate_qs_force(qs_force,error)

    nkind = SIZE(natom_of_kind)

    ALLOCATE (qs_force(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO ikind=1,nkind
      n = natom_of_kind(ikind)
      ALLOCATE (qs_force(ikind)%all_potential(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%core_overlap(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%gth_ppl(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%gth_nlcc(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%gth_ppnl(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%kinetic(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%overlap(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%overlap_admm(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%rho_core(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%rho_elec(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%rho_lri_elec(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%vhxc_atom(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%g0s_Vh_elec(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%repulsive(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%dispersion(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%other(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%ch_pulay(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%ehrenfest(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%efield(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%eev(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ! Always initialize ch_pulay to zero..
      qs_force(ikind)%ch_pulay  = 0.0_dp
      ALLOCATE (qs_force(ikind)%fock_4c(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%mp2_sep(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%mp2_non_sep(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (qs_force(ikind)%total(3,n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END DO

  END SUBROUTINE allocate_qs_force

! *****************************************************************************
!> \brief   Deallocate a Quickstep force data structure.
!> \param qs_force ...
!> \param error ...
!> \date    05.06.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_qs_force(qs_force,error)

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_qs_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, nkind
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(qs_force),cp_failure_level,routineP,error,failure)

    nkind = SIZE(qs_force)

    DO ikind=1,nkind

      IF (ASSOCIATED(qs_force(ikind)%all_potential)) THEN
        DEALLOCATE (qs_force(ikind)%all_potential,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%core_overlap)) THEN
        DEALLOCATE (qs_force(ikind)%core_overlap,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%gth_ppl)) THEN
        DEALLOCATE (qs_force(ikind)%gth_ppl,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%gth_nlcc)) THEN
        DEALLOCATE (qs_force(ikind)%gth_nlcc,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%gth_ppnl)) THEN
        DEALLOCATE (qs_force(ikind)%gth_ppnl,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%kinetic)) THEN
        DEALLOCATE (qs_force(ikind)%kinetic,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%overlap)) THEN
        DEALLOCATE (qs_force(ikind)%overlap,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%overlap_admm)) THEN
        DEALLOCATE (qs_force(ikind)%overlap_admm,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%rho_core)) THEN
        DEALLOCATE (qs_force(ikind)%rho_core,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%rho_elec)) THEN
        DEALLOCATE (qs_force(ikind)%rho_elec,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      IF (ASSOCIATED(qs_force(ikind)%rho_lri_elec)) THEN
        DEALLOCATE (qs_force(ikind)%rho_lri_elec,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%vhxc_atom)) THEN
        DEALLOCATE (qs_force(ikind)%vhxc_atom,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%g0s_Vh_elec)) THEN
        DEALLOCATE (qs_force(ikind)%g0s_Vh_elec,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%repulsive)) THEN
        DEALLOCATE (qs_force(ikind)%repulsive,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%dispersion)) THEN
        DEALLOCATE (qs_force(ikind)%dispersion,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%other)) THEN
        DEALLOCATE (qs_force(ikind)%other,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%total)) THEN
        DEALLOCATE (qs_force(ikind)%total,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%ch_pulay)) THEN
        DEALLOCATE (qs_force(ikind)%ch_pulay,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%fock_4c)) THEN
        DEALLOCATE (qs_force(ikind)%fock_4c,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%mp2_sep)) THEN
        DEALLOCATE (qs_force(ikind)%mp2_sep,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%mp2_non_sep)) THEN
        DEALLOCATE (qs_force(ikind)%mp2_non_sep,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%ehrenfest)) THEN
        DEALLOCATE (qs_force(ikind)%ehrenfest,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%efield)) THEN
        DEALLOCATE (qs_force(ikind)%efield,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%eev)) THEN
        DEALLOCATE (qs_force(ikind)%eev,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
    END DO

    DEALLOCATE (qs_force,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE deallocate_qs_force

! *****************************************************************************
!> \brief    Initialize a Quickstep force data structure.
!> \param qs_force ...
!> \param error ...
!> \date    15.07.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE zero_qs_force(qs_force,error)

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'zero_qs_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(qs_force),cp_failure_level,routineP,error,failure)

    DO ikind=1,SIZE(qs_force)
      qs_force(ikind)%all_potential(:,:) = 0.0_dp
      qs_force(ikind)%core_overlap(:,:) = 0.0_dp
      qs_force(ikind)%gth_ppl(:,:) = 0.0_dp
      qs_force(ikind)%gth_nlcc(:,:) = 0.0_dp
      qs_force(ikind)%gth_ppnl(:,:) = 0.0_dp
      qs_force(ikind)%kinetic(:,:) = 0.0_dp
      qs_force(ikind)%overlap(:,:) = 0.0_dp
      qs_force(ikind)%overlap_admm(:,:) = 0.0_dp
      qs_force(ikind)%rho_core(:,:) = 0.0_dp
      qs_force(ikind)%rho_elec(:,:) = 0.0_dp
      qs_force(ikind)%rho_lri_elec(:,:) = 0.0_dp
      qs_force(ikind)%vhxc_atom(:,:) = 0.0_dp
      qs_force(ikind)%g0s_Vh_elec(:,:) = 0.0_dp
      qs_force(ikind)%repulsive(:,:) = 0.0_dp
      qs_force(ikind)%dispersion(:,:) = 0.0_dp
      qs_force(ikind)%other(:,:) = 0.0_dp
      qs_force(ikind)%fock_4c(:,:) = 0.0_dp
      qs_force(ikind)%ehrenfest(:,:) = 0.0_dp
      qs_force(ikind)%efield(:,:) = 0.0_dp
      qs_force(ikind)%eev(:,:) = 0.0_dp
      qs_force(ikind)%mp2_non_sep(:,:) = 0.0_dp
      qs_force(ikind)%mp2_sep(:,:) = 0.0_dp
      qs_force(ikind)%total(:,:) = 0.0_dp
    END DO

  END SUBROUTINE zero_qs_force

! *****************************************************************************
!> \brief Duplicates a force type variable.
!> \param qs_force_input The force type variable to be duplicated
!> \param qs_force_output The force type duplicate
!> \param natom_of_kind ...
!> \param error ...
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!> \note
!>      The qs_force_variable is expected of the same size as
!>      qs_force_input if present.
! *****************************************************************************
  SUBROUTINE duplicate_qs_force(qs_force_input, qs_force_output, &
                                natom_of_kind,error)

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force_input, &
                                                qs_force_output
    INTEGER, DIMENSION(:), INTENT(IN)        :: natom_of_kind
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'duplicate_qs_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(qs_force_input),cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
      IF (.NOT. ASSOCIATED(qs_force_output)) THEN
        CALL allocate_qs_force(qs_force_output, natom_of_kind, error)
      ELSE
        CALL zero_qs_force(qs_force_output, error)
      END IF

      DO ikind=1,SIZE(natom_of_kind)
        qs_force_output(ikind)%all_potential(:,:) = &
           qs_force_input(ikind)%all_potential(:,:)
        qs_force_output(ikind)%core_overlap(:,:) = &
           qs_force_input(ikind)%core_overlap(:,:)
        qs_force_output(ikind)%gth_ppl(:,:) = qs_force_input(ikind)%gth_ppl(:,:)
        qs_force_output(ikind)%gth_nlcc(:,:) = qs_force_input(ikind)%gth_nlcc(:,:)
        qs_force_output(ikind)%gth_ppnl(:,:) = qs_force_input(ikind)%gth_ppnl(:,:)
        qs_force_output(ikind)%kinetic(:,:) = qs_force_input(ikind)%kinetic(:,:)
        qs_force_output(ikind)%overlap(:,:) = qs_force_input(ikind)%overlap(:,:)
        qs_force_output(ikind)%overlap_admm(:,:) = qs_force_input(ikind)%overlap_admm(:,:)
        qs_force_output(ikind)%rho_core(:,:) = qs_force_input(ikind)%rho_core(:,:)
        qs_force_output(ikind)%rho_elec(:,:) = qs_force_input(ikind)%rho_elec(:,:)
        qs_force_output(ikind)%rho_lri_elec(:,:) = qs_force_input(ikind)%rho_lri_elec(:,:)
        qs_force_output(ikind)%vhxc_atom(:,:) = qs_force_input(ikind)%vhxc_atom(:,:)
        qs_force_output(ikind)%g0s_Vh_elec(:,:) = &
           qs_force_input(ikind)%g0s_Vh_elec(:,:)
        qs_force_output(ikind)%repulsive(:,:) = &
           qs_force_input(ikind)%repulsive(:,:)
        qs_force_output(ikind)%dispersion(:,:) = &
           qs_force_input(ikind)%dispersion(:,:)
        qs_force_output(ikind)%other(:,:) = qs_force_input(ikind)%other(:,:)
        qs_force_output(ikind)%fock_4c(:,:) = qs_force_input(ikind)%fock_4c(:,:)
        qs_force_output(ikind)%ehrenfest(:,:) = qs_force_input(ikind)%ehrenfest(:,:)
        qs_force_output(ikind)%efield(:,:) = qs_force_input(ikind)%efield(:,:)
        qs_force_output(ikind)%eev(:,:) = qs_force_input(ikind)%eev(:,:)
        qs_force_output(ikind)%mp2_non_sep(:,:) = qs_force_input(ikind)%mp2_non_sep(:,:)
        qs_force_output(ikind)%mp2_sep(:,:) = qs_force_input(ikind)%mp2_sep(:,:)
        qs_force_output(ikind)%total(:,:) = qs_force_input(ikind)%total(:,:)
      END DO
    END IF

  END SUBROUTINE duplicate_qs_force

! *****************************************************************************
!> \brief Add force to a force_type  variable.
!> \param force Input force, dimension (3,natom)
!> \param qs_force The force type variable to be used
!> \param forcetype ...
!> \param atomic_kind_set ...
!> \param error ...
!> \par History
!>      07.2014 JGH
!> \author JGH
! *****************************************************************************
  SUBROUTINE add_qs_force(force, qs_force, forcetype, atomic_kind_set, error)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: force
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force
    CHARACTER(LEN=*), INTENT(IN)             :: forcetype
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_qs_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, iatom, ikind, natom_kind
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!   ------------------------------------------------------------------------

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(qs_force),cp_failure_level,routineP,error,failure)

    SELECT CASE (forcetype)
       CASE ("overlap_admm")
         DO ikind=1,SIZE(atomic_kind_set,1)
            atomic_kind => atomic_kind_set(ikind)
            CALL get_atomic_kind(atomic_kind=atomic_kind,natom=natom_kind)
            DO ia=1,natom_kind
              iatom = atomic_kind%atom_list(ia)
              qs_force(ikind)%overlap_admm(:,ia) = qs_force(ikind)%overlap_admm(:,ia) + force(:,iatom)
            END DO
         END DO
       CASE DEFAULT
         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE add_qs_force

END MODULE qs_force_types
