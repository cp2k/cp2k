!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Define all structures types related to force_fields
!> \par History
!>      10.2008 Teodoro Laino [tlaino] - University of Zurich
!>      added native support for AMBER forcefield
!>      10.2014 Moved kind types into force_field_kind_types.F [Ole Schuett]
!> \author MK (29.08.2003)
! *****************************************************************************
MODULE force_field_types

  USE damping_dipole_types,            ONLY: damping_info_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE pair_potential_types,            ONLY: pair_potential_p_release,&
                                             pair_potential_p_type
  USE shell_potential_types,           ONLY: shell_p_release,&
                                             shell_p_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_field_types'

! *** Define the derived structure types ***

! *****************************************************************************
  TYPE input_info_type
    CHARACTER (LEN=default_string_length), POINTER :: charge_atm(:)
    REAL(KIND=dp), POINTER                         :: charge(:)
    CHARACTER (LEN=default_string_length), POINTER :: apol_atm(:)
    REAL(KIND=dp), POINTER                         :: apol(:)
    CHARACTER (LEN=default_string_length), POINTER :: cpol_atm(:)
    REAL(KIND=dp), POINTER                         :: cpol(:)
    INTEGER, POINTER                               :: bond_kind(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_b(:)
    REAL(KIND=dp), POINTER                         :: bond_k(:,:)
    REAL(KIND=dp), POINTER                         :: bond_r0(:)
    REAL(KIND=dp), POINTER                         :: bond_cs(:)
    INTEGER, POINTER                               :: bend_kind(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_c(:)
    REAL(KIND=dp), POINTER                         :: bend_k(:)
    REAL(KIND=dp), POINTER                         :: bend_theta0(:)
    REAL(KIND=dp), POINTER                         :: bend_cb(:)
    REAL(KIND=dp), POINTER, DIMENSION(:)           :: bend_r012,&
                                                      bend_r032,&
                                                      bend_kbs12,&
                                                      bend_kbs32,&
                                                      bend_kss
    INTEGER, POINTER                               :: ub_kind(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_c(:)
    REAL(KIND=dp), POINTER                         :: ub_k(:,:)
    REAL(KIND=dp), POINTER                         :: ub_r0(:)
    INTEGER, POINTER                               :: torsion_kind(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_d(:)
    REAL(KIND=dp), POINTER                         :: torsion_k(:)
    INTEGER, POINTER                               :: torsion_m(:)
    REAL(KIND=dp), POINTER                         :: torsion_phi0(:)
    INTEGER, POINTER                               :: impr_kind(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_d(:)
    REAL(KIND=dp), POINTER                         :: impr_k(:)
    REAL(KIND=dp), POINTER                         :: impr_phi0(:)
    INTEGER, POINTER                               :: opbend_kind(:)
    CHARACTER (LEN=default_string_length), POINTER :: opbend_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: opbend_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: opbend_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: opbend_d(:)
    REAL(KIND=dp), POINTER                         :: opbend_k(:)
    REAL(KIND=dp), POINTER                         :: opbend_phi0(:)
    TYPE (pair_potential_p_type), POINTER          :: nonbonded
    TYPE (pair_potential_p_type), POINTER          :: nonbonded14
    TYPE (shell_p_type), DIMENSION(:), POINTER     :: shell_list
    TYPE (damping_info_type),DIMENSION(:), POINTER :: damping_list
  END TYPE input_info_type

! *****************************************************************************
  TYPE charmm_info_type
    CHARACTER (LEN=default_string_length), POINTER :: bond_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_b(:)
    REAL(KIND=dp), POINTER                         :: bond_k(:)
    REAL(KIND=dp), POINTER                         :: bond_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_c(:)
    REAL(KIND=dp), POINTER                         :: bend_k(:)
    REAL(KIND=dp), POINTER                         :: bend_theta0(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_c(:)
    REAL(KIND=dp), POINTER                         :: ub_k(:)
    REAL(KIND=dp), POINTER                         :: ub_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_d(:)
    REAL(KIND=dp), POINTER                         :: torsion_k(:)
    INTEGER, POINTER                               :: torsion_m(:)
    REAL(KIND=dp), POINTER                         :: torsion_phi0(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_d(:)
    REAL(KIND=dp), POINTER                         :: impr_k(:)
    REAL(KIND=dp), POINTER                         :: impr_phi0(:)
    CHARACTER (LEN=default_string_length), POINTER :: nonbond_a(:)
    REAL (KIND=dp), POINTER                        :: nonbond_eps(:)
    REAL (KIND=dp), POINTER                        :: nonbond_rmin2(:)
    CHARACTER (LEN=default_string_length), POINTER :: nonbond_a_14(:)
    REAL (KIND=dp), POINTER                        :: nonbond_eps_14(:)
    REAL (KIND=dp), POINTER                        :: nonbond_rmin2_14(:)
  END TYPE charmm_info_type

! *****************************************************************************
  TYPE amber_info_type
    CHARACTER (LEN=default_string_length), POINTER :: bond_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_b(:)
    REAL(KIND=dp), POINTER                         :: bond_k(:)
    REAL(KIND=dp), POINTER                         :: bond_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_c(:)
    REAL(KIND=dp), POINTER                         :: bend_k(:)
    REAL(KIND=dp), POINTER                         :: bend_theta0(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_d(:)
    REAL(KIND=dp), POINTER                         :: torsion_k(:)
    INTEGER, POINTER                               :: torsion_m(:)
    REAL(KIND=dp), POINTER                         :: torsion_phi0(:)
    CHARACTER (LEN=default_string_length), POINTER :: nonbond_a(:)
    REAL (KIND=dp), POINTER                        :: nonbond_eps(:)
    REAL (KIND=dp), POINTER                        :: nonbond_rmin2(:)
  END TYPE amber_info_type

! *****************************************************************************
  TYPE gromos_info_type
    INTEGER                                        :: ff_gromos_type
    REAL(KIND=dp), POINTER                         :: solvent_k(:)
    REAL(KIND=dp), POINTER                         :: solvent_r0(:)
    REAL(KIND=dp), POINTER                         :: bond_k(:)
    REAL(KIND=dp), POINTER                         :: bond_r0(:)
    REAL(KIND=dp), POINTER                         :: bend_k(:)
    REAL(KIND=dp), POINTER                         :: bend_theta0(:)
    REAL(KIND=dp), POINTER                         :: torsion_k(:)
    INTEGER, POINTER                               :: torsion_m(:)
    REAL(KIND=dp), POINTER                         :: torsion_phi0(:)
    REAL(KIND=dp), POINTER                         :: impr_k(:)
    REAL(KIND=dp), POINTER                         :: impr_phi0(:)
    CHARACTER (LEN=default_string_length), POINTER :: nonbond_a(:)
    REAL (KIND=dp), POINTER                        :: nonbond_c6(:,:)
    REAL (KIND=dp), POINTER                        :: nonbond_c12(:,:)
    CHARACTER (LEN=default_string_length), POINTER :: nonbond_a_14(:)
    REAL (KIND=dp), POINTER                        :: nonbond_c6_14(:,:)
    REAL (KIND=dp), POINTER                        :: nonbond_c12_14(:,:)
  END TYPE gromos_info_type

! *****************************************************************************
  TYPE force_field_type
    LOGICAL                                   :: shift_cutoff,&
                                                 do_nonbonded,&
                                                 electrostatics,&
                                                 multiple_potential,&
                                                 ignore_missing_critical
    INTEGER                                   :: ff_type
    REAL(KIND=dp)                             :: rcut_nb
    REAL(KIND=dp)                             :: ei_scale14
    REAL(KIND=dp)                             :: vdw_scale14
    REAL(KIND=dp)                             :: eps_spline,&
                                                 max_energy,&
                                                 emax_spline,&
                                                 rlow_nb
    INTEGER                                   :: npoints
    CHARACTER (LEN=default_path_length)       :: ff_file_name
    TYPE (input_info_type),POINTER            :: inp_info
    TYPE (charmm_info_type),POINTER           :: chm_info
    TYPE (gromos_info_type),POINTER           :: gro_info
    TYPE (amber_info_type), POINTER           :: amb_info
  END TYPE force_field_type

! *** Public subroutines ***

  PUBLIC :: init_ff_type,&
            deallocate_ff_type,&
            init_inp_info,&
            deallocate_inp_info

! *** Public data types ***

  PUBLIC :: charmm_info_type,&
            amber_info_type,&
            gromos_info_type,&
            input_info_type,&
            force_field_type

CONTAINS

! *****************************************************************************
!> \brief 1. Just NULLIFY and zero all the stuff
!> \param ff_type ...
!> \author ikuo
! *****************************************************************************
  SUBROUTINE init_ff_type (ff_type)
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_ff_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    !-----------------------------------------------------------------------------
    ! 1. Initialize
    !-----------------------------------------------------------------------------
    NULLIFY(ff_type%inp_info,ff_type%chm_info,ff_type%amb_info,ff_type%gro_info)
    ALLOCATE(ff_type%inp_info,STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(ff_type%chm_info,STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(ff_type%gro_info,STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(ff_type%amb_info,STAT=stat)
    CPASSERT(stat==0)

    !-----------------------------------------------------------------------------
    ! 2. Initialize and Nullify things in ff_type%inp_info
    !-----------------------------------------------------------------------------
    CALL init_inp_info(ff_type%inp_info)

    !-----------------------------------------------------------------------------
    ! 3. Initialize and Nullify things in ff_type%chm_info
    !-----------------------------------------------------------------------------
    CALL init_chm_info(ff_type%chm_info)

    !-----------------------------------------------------------------------------
    ! 4. Initialize and Nullify things in ff_type%gro_info
    !-----------------------------------------------------------------------------
    CALL init_gromos_info(ff_type%gro_info)

    !-----------------------------------------------------------------------------
    ! 5. Initialize and Nullify things in ff_type%amb_info
    !-----------------------------------------------------------------------------
    CALL init_amber_info(ff_type%amb_info)

  END SUBROUTINE init_ff_type

! *****************************************************************************
!> \brief Specific function to initialize the gro_info
!> \param gro_info ...
!> \author ikuo
! *****************************************************************************
  SUBROUTINE init_gromos_info(gro_info)
    TYPE(gromos_info_type), POINTER          :: gro_info

    NULLIFY(gro_info%solvent_k)
    NULLIFY(gro_info%solvent_r0)
    NULLIFY(gro_info%bond_k)
    NULLIFY(gro_info%bond_r0)
    NULLIFY(gro_info%bend_k)
    NULLIFY(gro_info%bend_theta0)
    NULLIFY(gro_info%torsion_k)
    NULLIFY(gro_info%torsion_m)
    NULLIFY(gro_info%torsion_phi0)
    NULLIFY(gro_info%impr_k)
    NULLIFY(gro_info%impr_phi0)
    NULLIFY(gro_info%nonbond_a)
    NULLIFY(gro_info%nonbond_c6)
    NULLIFY(gro_info%nonbond_c12)
    NULLIFY(gro_info%nonbond_a_14)
    NULLIFY(gro_info%nonbond_c6_14)
    NULLIFY(gro_info%nonbond_c12_14)

  END SUBROUTINE init_gromos_info

! *****************************************************************************
!> \brief Specific function to initialize the chm_info
!> \param chm_info ...
!> \par History
!>      none
!> \author teo
! *****************************************************************************
  SUBROUTINE init_chm_info(chm_info)
    TYPE(charmm_info_type), POINTER          :: chm_info

    NULLIFY(chm_info%bond_a)
    NULLIFY(chm_info%bond_b)
    NULLIFY(chm_info%bond_k)
    NULLIFY(chm_info%bond_r0)
    NULLIFY(chm_info%bend_a)
    NULLIFY(chm_info%bend_b)
    NULLIFY(chm_info%bend_c)
    NULLIFY(chm_info%bend_k)
    NULLIFY(chm_info%bend_theta0)
    NULLIFY(chm_info%ub_a)
    NULLIFY(chm_info%ub_b)
    NULLIFY(chm_info%ub_c)
    NULLIFY(chm_info%ub_k)
    NULLIFY(chm_info%ub_r0)
    NULLIFY(chm_info%torsion_a)
    NULLIFY(chm_info%torsion_b)
    NULLIFY(chm_info%torsion_c)
    NULLIFY(chm_info%torsion_d)
    NULLIFY(chm_info%torsion_k)
    NULLIFY(chm_info%torsion_m)
    NULLIFY(chm_info%torsion_phi0)
    NULLIFY(chm_info%impr_a)
    NULLIFY(chm_info%impr_b)
    NULLIFY(chm_info%impr_c)
    NULLIFY(chm_info%impr_d)
    NULLIFY(chm_info%impr_k)
    NULLIFY(chm_info%impr_phi0)
    NULLIFY(chm_info%nonbond_a)
    NULLIFY(chm_info%nonbond_eps)
    NULLIFY(chm_info%nonbond_rmin2)
    NULLIFY(chm_info%nonbond_a_14)
    NULLIFY(chm_info%nonbond_eps_14)
    NULLIFY(chm_info%nonbond_rmin2_14)

  END SUBROUTINE init_chm_info

! *****************************************************************************
!> \brief Specific function to initialize the chm_info
!> \param amb_info ...
!> \par History
!>      none
!> \author teo
! *****************************************************************************
  SUBROUTINE init_amber_info(amb_info)
    TYPE(amber_info_type), POINTER           :: amb_info

    NULLIFY(amb_info%bond_a)
    NULLIFY(amb_info%bond_b)
    NULLIFY(amb_info%bond_k)
    NULLIFY(amb_info%bond_r0)
    NULLIFY(amb_info%bend_a)
    NULLIFY(amb_info%bend_b)
    NULLIFY(amb_info%bend_c)
    NULLIFY(amb_info%bend_k)
    NULLIFY(amb_info%bend_theta0)
    NULLIFY(amb_info%torsion_a)
    NULLIFY(amb_info%torsion_b)
    NULLIFY(amb_info%torsion_c)
    NULLIFY(amb_info%torsion_d)
    NULLIFY(amb_info%torsion_k)
    NULLIFY(amb_info%torsion_m)
    NULLIFY(amb_info%torsion_phi0)
    NULLIFY(amb_info%nonbond_a)
    NULLIFY(amb_info%nonbond_eps)
    NULLIFY(amb_info%nonbond_rmin2)

  END SUBROUTINE init_amber_info

! *****************************************************************************
!> \brief Specific function to initialize the inp_info
!> \param inp_info ...
!> \par History
!>      none
!> \author teo
! *****************************************************************************
  SUBROUTINE init_inp_info(inp_info)
    TYPE(input_info_type), POINTER           :: inp_info

    NULLIFY(inp_info%charge_atm)
    NULLIFY(inp_info%charge)
    NULLIFY(inp_info%apol_atm)
    NULLIFY(inp_info%apol)
    NULLIFY(inp_info%cpol_atm)
    NULLIFY(inp_info%cpol)
    NULLIFY(inp_info%bond_kind)
    NULLIFY(inp_info%bond_a)
    NULLIFY(inp_info%bond_b)
    NULLIFY(inp_info%bond_k)
    NULLIFY(inp_info%bond_r0)
    NULLIFY(inp_info%bond_cs)
    NULLIFY(inp_info%bend_kind)
    NULLIFY(inp_info%bend_a)
    NULLIFY(inp_info%bend_b)
    NULLIFY(inp_info%bend_c)
    NULLIFY(inp_info%bend_k)
    NULLIFY(inp_info%bend_theta0)
    NULLIFY(inp_info%bend_cb)
    NULLIFY(inp_info%bend_r012)
    NULLIFY(inp_info%bend_r032)
    NULLIFY(inp_info%bend_kbs12)
    NULLIFY(inp_info%bend_kbs32)
    NULLIFY(inp_info%bend_kss)
    NULLIFY(inp_info%ub_kind)
    NULLIFY(inp_info%ub_a)
    NULLIFY(inp_info%ub_b)
    NULLIFY(inp_info%ub_c)
    NULLIFY(inp_info%ub_k)
    NULLIFY(inp_info%ub_r0)
    NULLIFY(inp_info%torsion_kind)
    NULLIFY(inp_info%torsion_a)
    NULLIFY(inp_info%torsion_b)
    NULLIFY(inp_info%torsion_c)
    NULLIFY(inp_info%torsion_d)
    NULLIFY(inp_info%torsion_k)
    NULLIFY(inp_info%torsion_m)
    NULLIFY(inp_info%torsion_phi0)
    NULLIFY(inp_info%impr_kind)
    NULLIFY(inp_info%impr_a)
    NULLIFY(inp_info%impr_b)
    NULLIFY(inp_info%impr_c)
    NULLIFY(inp_info%impr_d)
    NULLIFY(inp_info%impr_k)
    NULLIFY(inp_info%impr_phi0)
    NULLIFY(inp_info%opbend_kind)
    NULLIFY(inp_info%opbend_a)
    NULLIFY(inp_info%opbend_b)
    NULLIFY(inp_info%opbend_c)
    NULLIFY(inp_info%opbend_d)
    NULLIFY(inp_info%opbend_k)
    NULLIFY(inp_info%opbend_phi0)
    NULLIFY(inp_info%nonbonded)
    NULLIFY(inp_info%nonbonded14)
    NULLIFY(inp_info%shell_list)
    NULLIFY(inp_info%damping_list)

  END SUBROUTINE init_inp_info

! *****************************************************************************
!> \brief 1. Just DEALLOCATE all the stuff
!> \param ff_type ...
! *****************************************************************************
  SUBROUTINE deallocate_ff_type (ff_type)
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_ff_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    !-----------------------------------------------------------------------------
    ! 1. DEALLOCATE things in ff_type%inp_info
    !-----------------------------------------------------------------------------
    CALL deallocate_inp_info(ff_type%inp_info)

    !-----------------------------------------------------------------------------
    ! 2. DEALLOCATE things in ff_type%chm_info
    !-----------------------------------------------------------------------------
    CALL deallocate_chm_info(ff_type%chm_info)

    !-----------------------------------------------------------------------------
    ! 3. DEALLOCATE things in ff_type%gro_info
    !-----------------------------------------------------------------------------
    CALL deallocate_gromos_info(ff_type%gro_info)

    !-----------------------------------------------------------------------------
    ! 4. DEALLOCATE things in ff_type%amb_info
    !-----------------------------------------------------------------------------
    CALL deallocate_amb_info(ff_type%amb_info)

    !-----------------------------------------------------------------------------
    ! 5. DEALLOCATE things in ff_type
    !-----------------------------------------------------------------------------
    IF(ASSOCIATED(ff_type%inp_info)) THEN
       DEALLOCATE(ff_type%inp_info,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(ff_type%chm_info)) THEN
       DEALLOCATE(ff_type%chm_info,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(ff_type%gro_info)) THEN
       DEALLOCATE(ff_type%gro_info,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(ff_type%amb_info)) THEN
       DEALLOCATE(ff_type%amb_info,STAT=stat)
       CPASSERT(stat==0)
    END IF

  END SUBROUTINE deallocate_ff_type

! *****************************************************************************
!> \brief Specific function to deallocate the gro_info
!> \param gro_info ...
!> \author ikuo
! *****************************************************************************
  SUBROUTINE deallocate_gromos_info(gro_info)
    TYPE(gromos_info_type), POINTER          :: gro_info

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_gromos_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF(ASSOCIATED(gro_info%solvent_k)) THEN
       DEALLOCATE(gro_info%solvent_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%solvent_r0)) THEN
       DEALLOCATE(gro_info%solvent_r0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%bond_k)) THEN
       DEALLOCATE(gro_info%bond_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%bond_r0)) THEN
       DEALLOCATE(gro_info%bond_r0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%bend_k)) THEN
       DEALLOCATE(gro_info%bend_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%bend_theta0)) THEN
       DEALLOCATE(gro_info%bend_theta0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%torsion_k)) THEN
       DEALLOCATE(gro_info%torsion_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%torsion_m)) THEN
       DEALLOCATE(gro_info%torsion_m,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%torsion_phi0)) THEN
       DEALLOCATE(gro_info%torsion_phi0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%impr_k)) THEN
       DEALLOCATE(gro_info%impr_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%impr_phi0)) THEN
       DEALLOCATE(gro_info%impr_phi0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%nonbond_a)) THEN
       DEALLOCATE(gro_info%nonbond_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%nonbond_c6)) THEN
       DEALLOCATE(gro_info%nonbond_c6,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%nonbond_c12)) THEN
       DEALLOCATE(gro_info%nonbond_c12,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%nonbond_a_14)) THEN
       DEALLOCATE(gro_info%nonbond_a_14,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%nonbond_c6_14)) THEN
       DEALLOCATE(gro_info%nonbond_c6_14,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(gro_info%nonbond_c12_14)) THEN
       DEALLOCATE(gro_info%nonbond_c12_14,STAT=stat)
       CPASSERT(stat==0)
    END IF

  END SUBROUTINE deallocate_gromos_info

! *****************************************************************************
!> \brief Specific function to deallocate the chm_info
!> \param chm_info ...
!> \par History
!>      none
!> \author teo
! *****************************************************************************
  SUBROUTINE deallocate_chm_info(chm_info)
    TYPE(charmm_info_type), POINTER          :: chm_info

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_chm_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF(ASSOCIATED(chm_info%bond_a)) THEN
       DEALLOCATE(chm_info%bond_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%bond_b)) THEN
       DEALLOCATE(chm_info%bond_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%bond_k)) THEN
       DEALLOCATE(chm_info%bond_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%bond_r0)) THEN
       DEALLOCATE(chm_info%bond_r0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%bend_a)) THEN
       DEALLOCATE(chm_info%bend_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%bend_b)) THEN
       DEALLOCATE(chm_info%bend_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%bend_c)) THEN
       DEALLOCATE(chm_info%bend_c,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%bend_k)) THEN
       DEALLOCATE(chm_info%bend_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%bend_theta0)) THEN
       DEALLOCATE(chm_info%bend_theta0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%ub_a)) THEN
       DEALLOCATE(chm_info%ub_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%ub_b)) THEN
       DEALLOCATE(chm_info%ub_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%ub_c)) THEN
       DEALLOCATE(chm_info%ub_c,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%ub_k)) THEN
       DEALLOCATE(chm_info%ub_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%ub_r0)) THEN
       DEALLOCATE(chm_info%ub_r0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%torsion_a)) THEN
       DEALLOCATE(chm_info%torsion_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%torsion_b)) THEN
       DEALLOCATE(chm_info%torsion_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%torsion_c)) THEN
       DEALLOCATE(chm_info%torsion_c,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%torsion_d)) THEN
       DEALLOCATE(chm_info%torsion_d,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%torsion_k)) THEN
       DEALLOCATE(chm_info%torsion_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%torsion_m)) THEN
       DEALLOCATE(chm_info%torsion_m,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%torsion_phi0)) THEN
       DEALLOCATE(chm_info%torsion_phi0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%impr_a)) THEN
       DEALLOCATE(chm_info%impr_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%impr_b)) THEN
       DEALLOCATE(chm_info%impr_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%impr_c)) THEN
       DEALLOCATE(chm_info%impr_c,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%impr_d)) THEN
       DEALLOCATE(chm_info%impr_d,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%impr_k)) THEN
       DEALLOCATE(chm_info%impr_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%impr_phi0)) THEN
       DEALLOCATE(chm_info%impr_phi0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%nonbond_a)) THEN
       DEALLOCATE(chm_info%nonbond_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%nonbond_eps)) THEN
       DEALLOCATE(chm_info%nonbond_eps,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%nonbond_rmin2)) THEN
       DEALLOCATE(chm_info%nonbond_rmin2,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%nonbond_a_14)) THEN
       DEALLOCATE(chm_info%nonbond_a_14,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%nonbond_eps_14)) THEN
       DEALLOCATE(chm_info%nonbond_eps_14,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(chm_info%nonbond_rmin2_14)) THEN
       DEALLOCATE(chm_info%nonbond_rmin2_14,STAT=stat)
       CPASSERT(stat==0)
    END IF

  END SUBROUTINE deallocate_chm_info

! *****************************************************************************
!> \brief Specific function to deallocate the chm_info
!> \param amb_info ...
!> \par History
!>      none
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  SUBROUTINE deallocate_amb_info(amb_info)
    TYPE(amber_info_type), POINTER           :: amb_info

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_amb_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF(ASSOCIATED(amb_info%bond_a)) THEN
       DEALLOCATE(amb_info%bond_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%bond_b)) THEN
       DEALLOCATE(amb_info%bond_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%bond_k)) THEN
       DEALLOCATE(amb_info%bond_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%bond_r0)) THEN
       DEALLOCATE(amb_info%bond_r0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%bend_a)) THEN
       DEALLOCATE(amb_info%bend_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%bend_b)) THEN
       DEALLOCATE(amb_info%bend_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%bend_c)) THEN
       DEALLOCATE(amb_info%bend_c,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%bend_k)) THEN
       DEALLOCATE(amb_info%bend_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%bend_theta0)) THEN
       DEALLOCATE(amb_info%bend_theta0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%torsion_a)) THEN
       DEALLOCATE(amb_info%torsion_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%torsion_b)) THEN
       DEALLOCATE(amb_info%torsion_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%torsion_c)) THEN
       DEALLOCATE(amb_info%torsion_c,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%torsion_d)) THEN
       DEALLOCATE(amb_info%torsion_d,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%torsion_k)) THEN
       DEALLOCATE(amb_info%torsion_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%torsion_m)) THEN
       DEALLOCATE(amb_info%torsion_m,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%torsion_phi0)) THEN
       DEALLOCATE(amb_info%torsion_phi0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%nonbond_a)) THEN
       DEALLOCATE(amb_info%nonbond_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%nonbond_eps)) THEN
       DEALLOCATE(amb_info%nonbond_eps,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(amb_info%nonbond_rmin2)) THEN
       DEALLOCATE(amb_info%nonbond_rmin2,STAT=stat)
       CPASSERT(stat==0)
    END IF

  END SUBROUTINE deallocate_amb_info

! *****************************************************************************
!> \brief Specific function to deallocate the inp_info
!> \param inp_info ...
!> \par History
!>      none
!> \author teo
! *****************************************************************************
  SUBROUTINE deallocate_inp_info(inp_info)
    TYPE(input_info_type), POINTER           :: inp_info

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_inp_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF(ASSOCIATED(inp_info%charge_atm)) THEN
       DEALLOCATE(inp_info%charge_atm,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%charge)) THEN
       DEALLOCATE(inp_info%charge,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%apol_atm)) THEN
       DEALLOCATE(inp_info%apol_atm,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%apol)) THEN
       DEALLOCATE(inp_info%apol,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%cpol_atm)) THEN
       DEALLOCATE(inp_info%cpol_atm,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%cpol)) THEN
       DEALLOCATE(inp_info%cpol,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bond_kind)) THEN
       DEALLOCATE(inp_info%bond_kind,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bond_a)) THEN
       DEALLOCATE(inp_info%bond_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bond_b)) THEN
       DEALLOCATE(inp_info%bond_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bond_k)) THEN
       DEALLOCATE(inp_info%bond_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bond_r0)) THEN
       DEALLOCATE(inp_info%bond_r0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bond_cs)) THEN
       DEALLOCATE(inp_info%bond_cs,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bend_kind)) THEN
       DEALLOCATE(inp_info%bend_kind,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bend_a)) THEN
       DEALLOCATE(inp_info%bend_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bend_b)) THEN
       DEALLOCATE(inp_info%bend_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bend_c)) THEN
       DEALLOCATE(inp_info%bend_c,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bend_k)) THEN
       DEALLOCATE(inp_info%bend_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bend_theta0)) THEN
       DEALLOCATE(inp_info%bend_theta0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bend_cb)) THEN
       DEALLOCATE(inp_info%bend_cb,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bend_r012)) THEN
       DEALLOCATE(inp_info%bend_r012,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bend_r032)) THEN
       DEALLOCATE(inp_info%bend_r032,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bend_kbs12)) THEN
       DEALLOCATE(inp_info%bend_kbs12,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bend_kbs32)) THEN
       DEALLOCATE(inp_info%bend_kbs32,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%bend_kss)) THEN
       DEALLOCATE(inp_info%bend_kss,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%ub_kind)) THEN
       DEALLOCATE(inp_info%ub_kind,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%ub_a)) THEN
       DEALLOCATE(inp_info%ub_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%ub_b)) THEN
       DEALLOCATE(inp_info%ub_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%ub_c)) THEN
       DEALLOCATE(inp_info%ub_c,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%ub_k)) THEN
       DEALLOCATE(inp_info%ub_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%ub_r0)) THEN
       DEALLOCATE(inp_info%ub_r0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%torsion_kind)) THEN
       DEALLOCATE(inp_info%torsion_kind,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%torsion_a)) THEN
       DEALLOCATE(inp_info%torsion_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%torsion_b)) THEN
       DEALLOCATE(inp_info%torsion_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%torsion_c)) THEN
       DEALLOCATE(inp_info%torsion_c,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%torsion_d)) THEN
       DEALLOCATE(inp_info%torsion_d,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%torsion_k)) THEN
       DEALLOCATE(inp_info%torsion_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%torsion_m)) THEN
       DEALLOCATE(inp_info%torsion_m,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%torsion_phi0)) THEN
       DEALLOCATE(inp_info%torsion_phi0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%impr_kind)) THEN
       DEALLOCATE(inp_info%impr_kind,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%impr_a)) THEN
       DEALLOCATE(inp_info%impr_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%impr_b)) THEN
       DEALLOCATE(inp_info%impr_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%impr_c)) THEN
       DEALLOCATE(inp_info%impr_c,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%impr_d)) THEN
       DEALLOCATE(inp_info%impr_d,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%impr_k)) THEN
       DEALLOCATE(inp_info%impr_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%impr_phi0)) THEN
       DEALLOCATE(inp_info%impr_phi0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%opbend_kind)) THEN
       DEALLOCATE(inp_info%opbend_kind,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%opbend_a)) THEN
       DEALLOCATE(inp_info%opbend_a,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%opbend_b)) THEN
       DEALLOCATE(inp_info%opbend_b,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%opbend_c)) THEN
       DEALLOCATE(inp_info%opbend_c,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%opbend_d)) THEN
       DEALLOCATE(inp_info%opbend_d,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%opbend_k)) THEN
       DEALLOCATE(inp_info%opbend_k,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%opbend_phi0)) THEN
       DEALLOCATE(inp_info%opbend_phi0,STAT=stat)
       CPASSERT(stat==0)
    END IF
    IF(ASSOCIATED(inp_info%nonbonded)) THEN
       CALL pair_potential_p_release(inp_info%nonbonded)
    END IF
    IF(ASSOCIATED(inp_info%nonbonded14)) THEN
       CALL pair_potential_p_release(inp_info%nonbonded14)
    END IF
    IF(ASSOCIATED(inp_info%shell_list)) THEN
       CALL shell_p_release(inp_info%shell_list)
    END IF
    IF(ASSOCIATED(inp_info%damping_list)) THEN
       DEALLOCATE(inp_info%damping_list,STAT=stat)
       CPASSERT(stat==0)
    END IF

  END SUBROUTINE deallocate_inp_info

END MODULE force_field_types
