!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_force_types [1.0] *
!!
!!   NAME
!!     qs_force_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MK (03.06.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_force_types

  !USE cp_control_types,                ONLY: qs_control_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  TYPE qs_force_type
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: all_potential,&
                                         core_overlap,&
                                         gth_ppl,&
                                         gth_ppnl,&
                                         kinetic,&
                                         overlap,&
                                         rho_core,&
                                         rho_elec,&
                                         vhxc_atom,&
                                         g0s_Vh_elec,&
                                         hartree_2c,&
                                         hartree_3c,&
                                         kg_gpw_ekin_mol,&
                                         other,&
                                         ch_pulay,&
                                         total
  END TYPE qs_force_type

! *** Public data types ***

  PUBLIC :: qs_force_type

! *** Public subroutines ***

  PUBLIC :: allocate_qs_force,&
            deallocate_qs_force,&
            zero_qs_force, &
            duplicate_qs_force

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_qs_force(qs_force,natom_of_kind)

!   Purpose: Allocate a Quickstep force data structure.

!   History: - Creation (05.06.2002,MK)

!   ***************************************************************************

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force
    INTEGER, DIMENSION(:), INTENT(IN)        :: natom_of_kind

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE allocate_qs_force (MODULE qs_force_types)"

    INTEGER                                  :: ikind, istat, n, nkind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(qs_force)) THEN
      CALL deallocate_qs_force(qs_force)
    END IF

    nkind = SIZE(natom_of_kind)

    ALLOCATE (qs_force(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"qs_force",nkind*int_size)

    DO ikind=1,nkind
      n = natom_of_kind(ikind)
      ALLOCATE (qs_force(ikind)%all_potential(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%all_potential",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%core_overlap(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%core_overlap",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%gth_ppl(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%gth_ppl",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%gth_ppnl(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%gth_ppnl",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%kinetic(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%kinetic",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%overlap(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%overlap",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%rho_core(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%rho_core",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%rho_elec(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%rho_elec",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%vhxc_atom(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%vhxc_atom",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%g0s_Vh_elec(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%g0s_Vh_elec",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%hartree_2c(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%hartree_2c",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%hartree_3c(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%hartree_3c",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%kg_gpw_ekin_mol(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%kg_gpw_ekin_mol",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%other(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%other",3*n*dp_size)
      END IF
      ALLOCATE (qs_force(ikind)%ch_pulay(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%other",3*n*dp_size)
      END IF
      ! Always initialize ch_pulay to zero..
      qs_force(ikind)%ch_pulay  = 0.0_dp
      ALLOCATE (qs_force(ikind)%total(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"qs_force(ikind)%total",3*n*dp_size)
      END IF
    END DO

  END SUBROUTINE allocate_qs_force

! *****************************************************************************

  SUBROUTINE deallocate_qs_force(qs_force)

!   Purpose: Deallocate a Quickstep force data structure.

!   History: - Creation (05.06.2002,MK)

!   ***************************************************************************

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE deallocate_qs_force (MODULE qs_force_types)"

    INTEGER                                  :: ikind, istat, nkind

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(qs_force)) THEN
      CALL stop_program(routine,"The qs_force pointer is not associated "//&
                                "and cannot be deallocated")
    END IF

    nkind = SIZE(qs_force)

    DO ikind=1,nkind

      IF (ASSOCIATED(qs_force(ikind)%all_potential)) THEN
        DEALLOCATE (qs_force(ikind)%all_potential,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%all_potential")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%core_overlap)) THEN
        DEALLOCATE (qs_force(ikind)%core_overlap,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%core_overlap")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%gth_ppl)) THEN
        DEALLOCATE (qs_force(ikind)%gth_ppl,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%gth_ppl")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%gth_ppnl)) THEN
        DEALLOCATE (qs_force(ikind)%gth_ppnl,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%gth_ppnl")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%kinetic)) THEN
        DEALLOCATE (qs_force(ikind)%kinetic,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%kinetic")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%overlap)) THEN
        DEALLOCATE (qs_force(ikind)%overlap,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%overlap")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%rho_core)) THEN
        DEALLOCATE (qs_force(ikind)%rho_core,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%rho_core")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%rho_elec)) THEN
        DEALLOCATE (qs_force(ikind)%rho_elec,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%rho_elec")
        END IF
      END IF
 
      IF (ASSOCIATED(qs_force(ikind)%vhxc_atom)) THEN
        DEALLOCATE (qs_force(ikind)%vhxc_atom,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%vhxc_atom")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%g0s_Vh_elec)) THEN
        DEALLOCATE (qs_force(ikind)%g0s_Vh_elec,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%g0s_Vh_elec")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%hartree_2c)) THEN
        DEALLOCATE (qs_force(ikind)%hartree_2c,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%hartree_2c")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%hartree_3c)) THEN
        DEALLOCATE (qs_force(ikind)%hartree_3c,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%hartree_3c")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%kg_gpw_ekin_mol)) THEN
        DEALLOCATE (qs_force(ikind)%kg_gpw_ekin_mol,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%kg_gpw_ekin_mol")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%other)) THEN
        DEALLOCATE (qs_force(ikind)%other,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%other")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%total)) THEN
        DEALLOCATE (qs_force(ikind)%total,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%total")
        END IF
      END IF

      IF (ASSOCIATED(qs_force(ikind)%ch_pulay)) THEN
         DEALLOCATE (qs_force(ikind)%ch_pulay,STAT=istat)
         IF (istat /= 0) THEN
            CALL stop_memory(routine,"qs_force%total")
         END IF
      END IF

    END DO

    DEALLOCATE (qs_force,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"qs_force")

  END SUBROUTINE deallocate_qs_force

! *****************************************************************************

  SUBROUTINE zero_qs_force(qs_force)

!   Purpose: Initialize a Quickstep force data structure.

!   History: - Creation (15.07.2002,MK)

!   ***************************************************************************

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE init_qs_force (MODULE qs_force_types)"

    INTEGER                                  :: ikind

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(qs_force)) THEN
      CALL stop_program(routine,"The qs_force pointer is not associated "//&
                                "and cannot be initialized")
    END IF

    DO ikind=1,SIZE(qs_force)
      qs_force(ikind)%all_potential(:,:) = 0.0_dp
      qs_force(ikind)%core_overlap(:,:) = 0.0_dp
      qs_force(ikind)%gth_ppl(:,:) = 0.0_dp
      qs_force(ikind)%gth_ppnl(:,:) = 0.0_dp
      qs_force(ikind)%kinetic(:,:) = 0.0_dp
      qs_force(ikind)%overlap(:,:) = 0.0_dp
      qs_force(ikind)%rho_core(:,:) = 0.0_dp
      qs_force(ikind)%rho_elec(:,:) = 0.0_dp
      qs_force(ikind)%vhxc_atom(:,:) = 0.0_dp
      qs_force(ikind)%g0s_Vh_elec(:,:) = 0.0_dp
      qs_force(ikind)%hartree_2c(:,:) = 0.0_dp
      qs_force(ikind)%hartree_3c(:,:) = 0.0_dp
      qs_force(ikind)%kg_gpw_ekin_mol(:,:) = 0.0_dp
      qs_force(ikind)%other(:,:) = 0.0_dp
      qs_force(ikind)%total(:,:) = 0.0_dp
    END DO

  END SUBROUTINE zero_qs_force

! *****************************************************************************

!!****f* qs_force_types/duplicate_qs_force [1.0] *
!!
!!   NAME
!!     duplicate_qs_force
!!
!!   SYNOPSIS
!!     SUBROUTINE duplicate_qs_force(qs_force_input, qs_force_output, natom_of_kind)
!!       TYPE(qs_force_type), DIMENSION(:), POINTER :: qs_force_input, qs_force_output
!!       INTEGER, DIMENSION(:), INTENT(IN) :: natom_of_kind
!!     END SUBROUTINE duplicate_qs_force
!!
!!   FUNCTION
!!     Duplicates a force type variable.
!!
!!   NOTES
!!     The qs_force_variable is expected of the same size as 
!!     qs_force_input if present.
!!
!!   ARGUMENTS
!!     - qs_force_input: The force type variable to be duplicated
!!     - qs_force_output: The force type duplicate
!!     - natom_of_kind: 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE duplicate_qs_force(qs_force_input, qs_force_output, &
                                natom_of_kind)

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force_input, &
                                                qs_force_output
    INTEGER, DIMENSION(:), INTENT(IN)        :: natom_of_kind

    CHARACTER(len=*), PARAMETER :: &
      routine = "SUBROUTINE duplicate_qs_force (MODULE qs_force_types)"

    INTEGER                                  :: ikind
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    IF (.NOT.ASSOCIATED(qs_force_input)) THEN
      CALL stop_program(routine,"The qs_force_input pointer is not associated "//&
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
        qs_force_output(ikind)%gth_ppnl(:,:) = qs_force_input(ikind)%gth_ppnl(:,:)
        qs_force_output(ikind)%kinetic(:,:) = qs_force_input(ikind)%kinetic(:,:)
        qs_force_output(ikind)%overlap(:,:) = qs_force_input(ikind)%overlap(:,:)
        qs_force_output(ikind)%rho_core(:,:) = qs_force_input(ikind)%rho_core(:,:)
        qs_force_output(ikind)%rho_elec(:,:) = qs_force_input(ikind)%rho_elec(:,:)
        qs_force_output(ikind)%vhxc_atom(:,:) = qs_force_input(ikind)%vhxc_atom(:,:)
        qs_force_output(ikind)%g0s_Vh_elec(:,:) = &
           qs_force_input(ikind)%g0s_Vh_elec(:,:)
        qs_force_output(ikind)%hartree_2c(:,:) = qs_force_input(ikind)%hartree_2c(:,:)
        qs_force_output(ikind)%hartree_3c(:,:) = qs_force_input(ikind)%hartree_3c(:,:)
        qs_force_output(ikind)%kg_gpw_ekin_mol(:,:) = &
           qs_force_input(ikind)%kg_gpw_ekin_mol(:,:)
        qs_force_output(ikind)%other(:,:) = qs_force_input(ikind)%other(:,:)
        qs_force_output(ikind)%total(:,:) = qs_force_input(ikind)%total(:,:)
      END DO
    END IF

  END SUBROUTINE duplicate_qs_force
!***************************************************************************

END MODULE qs_force_types
