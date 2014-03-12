!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      none
!> \author MK (03.06.2002)
! *****************************************************************************
MODULE qs_force_types

  !USE cp_control_types,                ONLY: qs_control_type
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
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
                                         vhxc_atom,&
                                         g0s_Vh_elec,&
                                         repulsive,&
                                         dispersion,&
                                         other,&
                                         ch_pulay,&
                                         fock_4c,&
                                         hfx_ri,&
                                         ehrenfest,&
                                         efield,&
                                         eev,&
                                         mp2_sep,&
                                         mp2_non_sep,&
                                         total
  END TYPE qs_force_type

  PUBLIC :: qs_force_type

  PUBLIC :: allocate_qs_force,&
            deallocate_qs_force,&
            zero_qs_force, &
            duplicate_qs_force

CONTAINS

! *****************************************************************************
!> \brief   Allocate a Quickstep force data structure.
!> \author  MK
!> \date    05.06.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_qs_force(qs_force,natom_of_kind)

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force
    INTEGER, DIMENSION(:), INTENT(IN)        :: natom_of_kind

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_qs_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, n, nkind

    IF (ASSOCIATED(qs_force)) THEN
      CALL deallocate_qs_force(qs_force)
    END IF

    nkind = SIZE(natom_of_kind)

    ALLOCATE (qs_force(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "qs_force",nkind*int_size)

    DO ikind=1,nkind
      n = natom_of_kind(ikind)
      ALLOCATE (qs_force(ikind)%all_potential(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%all_potential",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%core_overlap(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%core_overlap",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%gth_ppl(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%gth_ppl",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%gth_nlcc(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%gth_nlcc",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%gth_ppnl(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%gth_ppnl",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%kinetic(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%kinetic",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%overlap(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%overlap",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%overlap_admm(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%overlap_admm",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%rho_core(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%rho_core",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%rho_elec(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%rho_elec",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%vhxc_atom(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%vhxc_atom",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%g0s_Vh_elec(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%g0s_Vh_elec",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%repulsive(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%repulsive",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%dispersion(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%dispersion",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%other(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%other",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%ch_pulay(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%other",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%ehrenfest(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%ehrenfest",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%efield(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%efield",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%eev(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%eev",3*n*dp_size)
      END IF
      ! Always initialize ch_pulay to zero..
      qs_force(ikind)%ch_pulay  = 0.0_dp
      ALLOCATE (qs_force(ikind)%fock_4c(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%fock_4c",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%hfx_ri(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%hfx_ri",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%mp2_sep(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%mp2_sep",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%mp2_non_sep(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%mp2_non_sep",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%total(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "qs_force(ikind)%total",3*n*dp_size)
      END IF
    END DO

  END SUBROUTINE allocate_qs_force

! *****************************************************************************
!> \brief   Deallocate a Quickstep force data structure.
!> \author  MK
!> \date    05.06.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_qs_force(qs_force)

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_qs_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, nkind

    IF (.NOT.ASSOCIATED(qs_force)) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The qs_force pointer is not associated "//&
                        "and cannot be deallocated")
    END IF

    nkind = SIZE(qs_force)

    DO ikind=1,nkind

      IF (ASSOCIATED(qs_force(ikind)%all_potential)) THEN
        DEALLOCATE (qs_force(ikind)%all_potential,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%all_potential")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%core_overlap)) THEN
        DEALLOCATE (qs_force(ikind)%core_overlap,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%core_overlap")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%gth_ppl)) THEN
        DEALLOCATE (qs_force(ikind)%gth_ppl,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%gth_ppl")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%gth_nlcc)) THEN
        DEALLOCATE (qs_force(ikind)%gth_nlcc,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%gth_nlcc")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%gth_ppnl)) THEN
        DEALLOCATE (qs_force(ikind)%gth_ppnl,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%gth_ppnl")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%kinetic)) THEN
        DEALLOCATE (qs_force(ikind)%kinetic,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%kinetic")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%overlap)) THEN
        DEALLOCATE (qs_force(ikind)%overlap,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%overlap")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%overlap_admm)) THEN
        DEALLOCATE (qs_force(ikind)%overlap_admm,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%overlap_admm")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%rho_core)) THEN
        DEALLOCATE (qs_force(ikind)%rho_core,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%rho_core")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%rho_elec)) THEN
        DEALLOCATE (qs_force(ikind)%rho_elec,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%rho_elec")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%vhxc_atom)) THEN
        DEALLOCATE (qs_force(ikind)%vhxc_atom,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%vhxc_atom")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%g0s_Vh_elec)) THEN
        DEALLOCATE (qs_force(ikind)%g0s_Vh_elec,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%g0s_Vh_elec")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%repulsive)) THEN
        DEALLOCATE (qs_force(ikind)%repulsive,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%repulsive")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%dispersion)) THEN
        DEALLOCATE (qs_force(ikind)%dispersion,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%dispersion")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%other)) THEN
        DEALLOCATE (qs_force(ikind)%other,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%other")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%total)) THEN
        DEALLOCATE (qs_force(ikind)%total,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%total")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%ch_pulay)) THEN
         DEALLOCATE (qs_force(ikind)%ch_pulay,STAT=istat)
         IF (istat /= 0) THEN
            CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%total")
         END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%fock_4c)) THEN
        DEALLOCATE (qs_force(ikind)%fock_4c,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%fock_4c")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%mp2_sep)) THEN
        DEALLOCATE (qs_force(ikind)%mp2_sep,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%mp2_sep")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%mp2_non_sep)) THEN
        DEALLOCATE (qs_force(ikind)%mp2_non_sep,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%mp2_non_sep")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%hfx_ri)) THEN
        DEALLOCATE (qs_force(ikind)%hfx_ri,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%hfx_ri")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%ehrenfest)) THEN
        DEALLOCATE (qs_force(ikind)%ehrenfest,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%ehrenfest")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%efield)) THEN
        DEALLOCATE (qs_force(ikind)%efield,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%efield")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%eev)) THEN
        DEALLOCATE (qs_force(ikind)%eev,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"qs_force%eev")
        END IF
      END IF
    END DO

    DEALLOCATE (qs_force,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"qs_force")

  END SUBROUTINE deallocate_qs_force

! *****************************************************************************
!> \brief    Initialize a Quickstep force data structure.
!> \author  MK
!> \date    15.07.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE zero_qs_force(qs_force)

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force

    CHARACTER(len=*), PARAMETER :: routineN = 'zero_qs_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind

    IF (.NOT.ASSOCIATED(qs_force)) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The qs_force pointer is not associated "//&
                        "and cannot be initialized")
    END IF

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
      qs_force(ikind)%vhxc_atom(:,:) = 0.0_dp
      qs_force(ikind)%g0s_Vh_elec(:,:) = 0.0_dp
      qs_force(ikind)%repulsive(:,:) = 0.0_dp
      qs_force(ikind)%dispersion(:,:) = 0.0_dp
      qs_force(ikind)%other(:,:) = 0.0_dp
      qs_force(ikind)%fock_4c(:,:) = 0.0_dp
      qs_force(ikind)%hfx_ri(:,:) = 0.0_dp
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
!> \param natom_of_kind
!> \note
!>      The qs_force_variable is expected of the same size as
!>      qs_force_input if present.
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE duplicate_qs_force(qs_force_input, qs_force_output, &
                                natom_of_kind)

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force_input, &
                                                qs_force_output
    INTEGER, DIMENSION(:), INTENT(IN)        :: natom_of_kind

    CHARACTER(len=*), PARAMETER :: routineN = 'duplicate_qs_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    IF (.NOT.ASSOCIATED(qs_force_input)) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The qs_force_input pointer is not associated "//&
                        "and cannot be initialized")
    END IF

    IF (.NOT. failure) THEN
      IF (.NOT. ASSOCIATED(qs_force_output)) THEN
        CALL allocate_qs_force(qs_force_output, natom_of_kind)
      ELSE
        CALL zero_qs_force(qs_force_output)
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
        qs_force_output(ikind)%vhxc_atom(:,:) = qs_force_input(ikind)%vhxc_atom(:,:)
        qs_force_output(ikind)%g0s_Vh_elec(:,:) = &
           qs_force_input(ikind)%g0s_Vh_elec(:,:)
        qs_force_output(ikind)%repulsive(:,:) = &
           qs_force_input(ikind)%repulsive(:,:)
        qs_force_output(ikind)%dispersion(:,:) = &
           qs_force_input(ikind)%dispersion(:,:)
        qs_force_output(ikind)%other(:,:) = qs_force_input(ikind)%other(:,:)
        qs_force_output(ikind)%fock_4c(:,:) = qs_force_input(ikind)%fock_4c(:,:)
        qs_force_output(ikind)%hfx_ri(:,:) = qs_force_input(ikind)%hfx_ri(:,:)
        qs_force_output(ikind)%ehrenfest(:,:) = qs_force_input(ikind)%ehrenfest(:,:)
        qs_force_output(ikind)%efield(:,:) = qs_force_input(ikind)%efield(:,:)
        qs_force_output(ikind)%eev(:,:) = qs_force_input(ikind)%eev(:,:)
        qs_force_output(ikind)%mp2_non_sep(:,:) = qs_force_input(ikind)%mp2_non_sep(:,:)
        qs_force_output(ikind)%mp2_sep(:,:) = qs_force_input(ikind)%mp2_sep(:,:)
        qs_force_output(ikind)%total(:,:) = qs_force_input(ikind)%total(:,:)
      END DO
    END IF

  END SUBROUTINE duplicate_qs_force

END MODULE qs_force_types
