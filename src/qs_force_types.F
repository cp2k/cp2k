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

  USE kinds, ONLY: int_size,&
                   wp => dp,&
                   wp_size => dp_size

  USE global_types, ONLY: global_environment_type
  USE termination,  ONLY: stop_memory,&
                          stop_program

  IMPLICIT NONE

  PRIVATE

  TYPE qs_force_type
    REAL(wp), DIMENSION(:,:), POINTER :: core_overlap,&
                                         gth_ppl,&
                                         gth_ppnl,&
                                         kinetic,&
                                         overlap,&
                                         rho_core,&
                                         rho_elec,&
                                         total
  END TYPE qs_force_type

! *** Public data types ***

  PUBLIC :: qs_force_type

! *** Public subroutines ***

  PUBLIC :: allocate_qs_force,&
            deallocate_qs_force,&
            write_qs_force

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_qs_force(qs_force,nkind)

!   Purpose: Allocate a Quickstep force data structure.

!   History: - Creation (05.06.2002,MK)

!   ***************************************************************************

    TYPE(qs_force_type), DIMENSION(:), POINTER :: qs_force
    INTEGER, INTENT(IN)                        :: nkind

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_qs_force (MODULE qs_force_types)"

!   *** Local variables ***

    INTEGER :: ikind,istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(qs_force)) THEN
      CALL deallocate_qs_force(qs_force)
    END IF

    ALLOCATE (qs_force(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"qs_force",nkind*int_size)

    DO ikind=1,nkind
      NULLIFY (qs_force(ikind)%core_overlap)
      NULLIFY (qs_force(ikind)%gth_ppl)
      NULLIFY (qs_force(ikind)%gth_ppnl)
      NULLIFY (qs_force(ikind)%kinetic)
      NULLIFY (qs_force(ikind)%overlap)
      NULLIFY (qs_force(ikind)%rho_core)
      NULLIFY (qs_force(ikind)%rho_elec)
      NULLIFY (qs_force(ikind)%total)
    END DO

  END SUBROUTINE allocate_qs_force

! *****************************************************************************

  SUBROUTINE deallocate_qs_force(qs_force)

!   Purpose: Deallocate a Quickstep force data structure.

!   History: - Creation (05.06.2002,MK)

!   ***************************************************************************

    TYPE(qs_force_type), DIMENSION(:), POINTER :: qs_force

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_qs_force (MODULE qs_force_types)"

!   *** Local variables ***

    INTEGER :: ikind,istat,nkind

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(qs_force)) THEN
      CALL stop_program(routine,"The qs_force pointer is not associated "//&
                                "and cannot be deallocated")
    END IF

    nkind = SIZE(qs_force)

    DO ikind=1,nkind

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

      IF (ASSOCIATED(qs_force(ikind)%total)) THEN
        DEALLOCATE (qs_force(ikind)%total,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"qs_force%total")
        END IF
      END IF

    END DO

    DEALLOCATE (qs_force,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"qs_force")

  END SUBROUTINE deallocate_qs_force

! *****************************************************************************

  SUBROUTINE write_qs_force(qs_force,globenv)

!   Purpose: Write a Quickstep force data structure to 

!   History: - Creation (05.06.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN)  :: globenv
    TYPE(qs_force_type), DIMENSION(:), POINTER :: qs_force

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE write_qs_force (MODULE qs_force_types)"

!   *** Local variables ***

    INTEGER :: iatom,ikind,istat,output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (.NOT.globenv%print%forces) RETURN

    IF (.NOT.ASSOCIATED(qs_force)) THEN
      CALL stop_program(routine,"The qs_force pointer is not associated "//&
                                "and cannot be printed")
    END IF

    output_unit = globenv%scr

    WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") "FORCES"

    DO ikind=1,SIZE(qs_force)
      WRITE (UNIT=output_unit,FMT="(/,T2,A,I5,/)") "ATOMIC KIND",ikind
      IF (ASSOCIATED(qs_force(ikind)%core_overlap)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%core_overlap,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "core_overlap:",iatom,qs_force(ikind)%core_overlap(1:3,iatom)
        END DO
      END IF
      IF (ASSOCIATED(qs_force(ikind)%gth_ppl)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%gth_ppl,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "gth_ppl:",iatom,qs_force(ikind)%gth_ppl(1:3,iatom)
        END DO
      END IF
      IF (ASSOCIATED(qs_force(ikind)%gth_ppnl)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%gth_ppnl,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "gth_ppnl:",iatom,qs_force(ikind)%gth_ppnl(1:3,iatom)
        END DO
      END IF
      IF (ASSOCIATED(qs_force(ikind)%kinetic)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%kinetic,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "kinetic:",iatom,qs_force(ikind)%kinetic(1:3,iatom)
        END DO
      END IF
      IF (ASSOCIATED(qs_force(ikind)%overlap)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%overlap,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "overlap:",iatom,qs_force(ikind)%overlap(1:3,iatom)
        END DO
      END IF
      IF (ASSOCIATED(qs_force(ikind)%rho_core)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%rho_core,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "rho_core:",iatom,qs_force(ikind)%rho_core(1:3,iatom)
        END DO
      END IF
      IF (ASSOCIATED(qs_force(ikind)%rho_elec)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%rho_elec,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "rho_elec:",iatom,qs_force(ikind)%rho_elec(1:3,iatom)
        END DO
      END IF
!      IF (ASSOCIATED(qs_force(ikind)%total)) THEN
!        DO iatom=1,SIZE(qs_force(ikind)%total,2)
!          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
!            "total:",iatom,qs_force(ikind)%rho_total(1:3,iatom)
!        END DO
!      END IF
    END DO

  END SUBROUTINE write_qs_force

! *****************************************************************************

END MODULE qs_force_types
