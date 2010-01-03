!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition of the atomic potential types.
!> \par History
!>      GT, 22.09.2002: added elp_potential_types
!> \author Matthias Krack (04.07.2000)
! *****************************************************************************
MODULE external_potential_types
  USE bibliography,                    ONLY: Goedecker1996,&
                                             Hartwigsen1998,&
                                             Krack2000,&
                                             Krack2005,&
                                             cite_reference
  USE cp_linked_list_val,              ONLY: cp_sll_val_next,&
                                             cp_sll_val_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_object,&
                                             parser_search_string,&
                                             parser_test_next_token
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_list_get,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE input_val_types,                 ONLY: val_get,&
                                             val_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: dfac,&
                                             fac,&
                                             pi,&
                                             rootpi
  USE mathlib,                         ONLY: symmetrize_matrix
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: co,&
                                             coset,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset,&
                                             nso
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE periodic_table,                  ONLY: ptable
  USE string_utilities,                ONLY: remove_word,&
                                             uppercase
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'external_potential_types'

  ! Define the all-electron potential type

  ! Literature: M. Krack and M. Parrinello,
  !             Phys. Chem. Chem. Phys. 2, 2105 (2000)

! *****************************************************************************
  TYPE all_potential_type
     !MK PRIVATE
     CHARACTER(LEN=default_string_length)   :: name
     CHARACTER(LEN=default_string_length),&
          DIMENSION(2)                      :: description
     REAL(KIND = dp)                        :: alpha_core_charge,&
                                               ccore_charge,&
                                               core_charge_radius,&
                                               zeff, zeff_correction
     INTEGER                                :: z
     INTEGER, DIMENSION(:), POINTER         :: elec_conf
     REAL(KIND = dp),DIMENSION(:),POINTER   :: orb_radius_3c
     REAL(KIND = dp),DIMENSION(:,:),POINTER :: set_radius_3c
  END TYPE all_potential_type

  ! Define the GTH potential type

  ! Literature: - S. Goedecker, M. Teter and J. Hutter,
  !               Phys. Rev. B 54, 1703 (1996)
  !             - C. Hartwigsen, S. Goedecker and J. Hutter,
  !               Phys. Rev. B 58, 3641 (1998)
  !             - M. Krack,
  !               Theor. Chem. Acc. 114, 145 (2005)

! *****************************************************************************
  TYPE gth_potential_type
     PRIVATE
     CHARACTER(LEN=default_string_length)       :: name
     CHARACTER(LEN=default_string_length),&
          DIMENSION(4)                          :: description
     REAL(KIND = dp)                            :: alpha_core_charge,&
          alpha_ppl,ccore_charge,&
          cerf_ppl,zeff,&
          core_charge_radius,&
          ppl_radius,ppnl_radius,&
          zeff_correction
     INTEGER                                    :: lppnl,lprj_ppnl_max,&
                                                   nexp_ppl,nppnl,&
                                                   nprj_ppnl_max,z
     REAL(KIND = dp), DIMENSION(:), POINTER     :: alpha_ppnl,cexp_ppl
     INTEGER, DIMENSION(:), POINTER             :: elec_conf,nprj_ppnl
     REAL(KIND = dp), DIMENSION(:,:), POINTER   :: cprj,cprj_ppnl,vprj_ppnl
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: hprj_ppnl
  END TYPE gth_potential_type

  ! Define the effective point charge potential type
! *****************************************************************************
  TYPE elp_potential_type
     PRIVATE
     CHARACTER(LEN=default_string_length)   :: name
     CHARACTER(LEN=default_string_length),&
          DIMENSION(2)                      :: description
     REAL(KIND = dp)                        :: alpha_core_charge,&
                                               ccore_charge,&
                                               cerf_ppl,&
                                               core_charge_radius,&
                                               ccore_charge_radius,&
                                               ppl_radius,qeff,a_pol
     INTEGER                                :: nalpha_ppl
     INTEGER,DIMENSION (:), POINTER         :: nexp_ppl
     REAL(KIND = dp),DIMENSION(:), POINTER  :: alpha_ppl
     REAL(KIND = dp),DIMENSION(:,:),POINTER :: cexp_ppl
  END TYPE elp_potential_type

  ! Public subroutines
  PUBLIC :: allocate_potential,&
            deallocate_potential,&
            get_potential,&
            init_potential,&
            read_potential,&
            set_potential,&
            set_default_all_potential,&
            write_potential

  ! Public data types

  PUBLIC :: all_potential_type,&
            elp_potential_type,&
            gth_potential_type

  INTERFACE allocate_potential
     MODULE PROCEDURE allocate_all_potential,&
                      allocate_elp_potential,&
                      allocate_gth_potential
  END INTERFACE

  INTERFACE deallocate_potential
     MODULE PROCEDURE deallocate_all_potential,&
                      deallocate_elp_potential,&
                      deallocate_gth_potential
  END INTERFACE

  INTERFACE get_potential
     MODULE PROCEDURE get_all_potential,&
                      get_elp_potential,&
                      get_gth_potential
  END INTERFACE

  INTERFACE init_potential
     MODULE PROCEDURE init_all_potential,&
                      init_gth_potential
  END INTERFACE

  INTERFACE read_potential
     MODULE PROCEDURE read_all_potential,&
                      read_elp_potential,&
                      read_gth_potential
  END INTERFACE

  INTERFACE set_potential
     MODULE PROCEDURE set_all_potential,&
                      set_elp_potential,&
                      set_gth_potential
  END INTERFACE

  INTERFACE write_potential
     MODULE PROCEDURE write_all_potential,&
                      write_elp_potential,&
                      write_gth_potential
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief   Allocate an atomic all-electron potential data set.
!> \author  MK
!> \date    25.07.2000,
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_all_potential(potential,error)
    TYPE(all_potential_type), POINTER        :: potential
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_all_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(potential)) CALL deallocate_potential(potential,error)

    ALLOCATE (potential,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY (potential%elec_conf)
    NULLIFY (potential%orb_radius_3c)
    NULLIFY (potential%set_radius_3c)

    potential%description(1) = "All-electron potential"
    potential%description(2) = "Krack, Parrinello, PCCP 2, 2105 (2000)"

  END SUBROUTINE allocate_all_potential

! *****************************************************************************
!> \brief   Allocate an effective charge potential data set.
!> \author  gt
!> \date    22.09.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_elp_potential(potential,error)
    TYPE(elp_potential_type), POINTER        :: potential
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_elp_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(potential)) CALL deallocate_potential(potential,error)

    ALLOCATE (potential,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    potential%a_pol = 0.0_dp
    NULLIFY (potential%alpha_ppl)
    NULLIFY (potential%cexp_ppl)
    NULLIFY (potential%nexp_ppl)

    potential%description(1) = "Local pseudopotential for KG-MD"

  END SUBROUTINE allocate_elp_potential

! *****************************************************************************
!> \brief   Allocate an atomic GTH potential data set.
!> \author  MK
!> \date    25.07.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_gth_potential(potential,error)
    TYPE(gth_potential_type), POINTER        :: potential
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_gth_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(potential)) CALL deallocate_potential(potential,error)

    ALLOCATE (potential,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY (potential%alpha_ppnl)
    NULLIFY (potential%cexp_ppl)
    NULLIFY (potential%elec_conf)
    NULLIFY (potential%nprj_ppnl)
    NULLIFY (potential%cprj)
    NULLIFY (potential%cprj_ppnl)
    NULLIFY (potential%vprj_ppnl)
    NULLIFY (potential%hprj_ppnl)

    potential%description(1) = "Goedecker-Teter-Hutter pseudopotential"
    potential%description(2) = "Goedecker et al., PRB 54, 1703 (1996)"
    potential%description(3) = "Hartwigsen et al., PRB 58, 3641 (1998)"
    potential%description(4) = "Krack, TCA 114, 145 (2005)"

  END SUBROUTINE allocate_gth_potential

! *****************************************************************************
!> \brief   Deallocate an atomic all-electron potential data set.
!> \author  MK
!> \date    03.11.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_all_potential(potential,error)
    TYPE(all_potential_type), POINTER        :: potential
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_all_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(potential)) THEN
       DEALLOCATE (potential%elec_conf,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (potential%orb_radius_3c,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (potential%set_radius_3c,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (potential,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
       CALL stop_program(routineN,moduleN,&
            __LINE__,&
            "The pointer potential is not associated")
    END IF

  END SUBROUTINE deallocate_all_potential

! *****************************************************************************
!> \brief   deallocate an effective charge potential data set.
!> \author  gt
!> \date    22.09.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_elp_potential(potential,error)
    TYPE(elp_potential_type), POINTER        :: potential
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_elp_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(potential)) THEN
       IF (ASSOCIATED(potential%nexp_ppl)) THEN
          DEALLOCATE (potential%nexp_ppl,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(potential%cexp_ppl)) THEN
          DEALLOCATE (potential%cexp_ppl,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(potential%alpha_ppl)) THEN
          DEALLOCATE (potential%alpha_ppl,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       DEALLOCATE (potential,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
       CALL stop_program(routineN,moduleN,&
            __LINE__,&
            "The pointer potential is not associated")
    END IF

  END SUBROUTINE deallocate_elp_potential

! *****************************************************************************
!> \brief   Deallocate an atomic GTH potential data set.
!> \author  MK
!> \date    03.11.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_gth_potential(potential,error)
    TYPE(gth_potential_type), POINTER        :: potential
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_gth_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(potential)) THEN

       DEALLOCATE (potential%elec_conf,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       !     *** Deallocate the parameters of the local part ***

       IF (ASSOCIATED(potential%cexp_ppl)) THEN
          DEALLOCATE (potential%cexp_ppl,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       !     *** Deallocate the parameters of the non-local part ***
       IF (ASSOCIATED(potential%alpha_ppnl)) THEN
          DEALLOCATE (potential%alpha_ppnl,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (potential%cprj,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (potential%cprj_ppnl,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (potential%hprj_ppnl,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (potential%nprj_ppnl,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (potential%vprj_ppnl,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       DEALLOCATE (potential,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE

       CALL stop_program(routineN,moduleN,&
            __LINE__,&
            "The pointer potential is not associated")
    END IF

  END SUBROUTINE deallocate_gth_potential

! *****************************************************************************
!> \brief   Get informations about an all-electron potential data set.
!> \author  MK
!> \date    11.01.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_all_potential(potential,name,alpha_core_charge,&
       ccore_charge,core_charge_radius,z,zeff,&
       zeff_correction, elec_conf)
    TYPE(all_potential_type), POINTER        :: potential
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: alpha_core_charge, &
                                                ccore_charge, &
                                                core_charge_radius
    INTEGER, INTENT(OUT), OPTIONAL           :: z
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: zeff, zeff_correction
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: elec_conf

    CHARACTER(len=*), PARAMETER :: routineN = 'get_all_potential', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(potential)) THEN

       IF (PRESENT(name)) name = potential%name
       IF (PRESENT(alpha_core_charge))&
            alpha_core_charge = potential%alpha_core_charge
       IF (PRESENT(ccore_charge)) ccore_charge = potential%ccore_charge
       IF (PRESENT(core_charge_radius))&
            core_charge_radius = potential%core_charge_radius
       IF (PRESENT(z)) z = potential%z
       IF (PRESENT(zeff)) zeff = potential%zeff
       IF (PRESENT(zeff_correction)) zeff_correction = potential%zeff_correction
       IF (PRESENT(elec_conf)) elec_conf => potential%elec_conf

    ELSE

       CALL stop_program(routineN,moduleN,&
            __LINE__,&
            "The pointer potential is not associated")

    END IF

  END SUBROUTINE get_all_potential

! *****************************************************************************
!> \brief   Get informations about effective pint charge potential.
!> \author  GT
!> \date    22-09-2002
!> \par History
!>          28-11-2003 gt: supports a local gth-like analytical form
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_elp_potential(potential,name,alpha_core_charge,&
       alpha_ppl,cerf_ppl,cexp_ppl,nalpha_ppl,&
       ccore_charge,core_charge_radius,&
       ccore_charge_radius,&
       nexp_ppl,ppl_radius,ppl_present,a_pol,qeff)
    TYPE(elp_potential_type), POINTER        :: potential
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: alpha_core_charge
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: alpha_ppl
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: cerf_ppl
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: cexp_ppl
    INTEGER, OPTIONAL                        :: nalpha_ppl
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: ccore_charge, &
                                                core_charge_radius, &
                                                ccore_charge_radius
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: nexp_ppl
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: ppl_radius
    LOGICAL, INTENT(OUT), OPTIONAL           :: ppl_present
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: a_pol, qeff

    CHARACTER(len=*), PARAMETER :: routineN = 'get_elp_potential', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(potential)) THEN

       IF (PRESENT(name)) name = potential%name
       IF (PRESENT(qeff)) qeff = potential%qeff
       IF (PRESENT(a_pol)) a_pol = potential%a_pol
       IF (PRESENT(nalpha_ppl)) nalpha_ppl=potential%nalpha_ppl
       IF (PRESENT(alpha_ppl)) alpha_ppl => potential%alpha_ppl
       IF (PRESENT(alpha_core_charge))&
            alpha_core_charge = potential%alpha_core_charge
       IF (PRESENT(cerf_ppl)) cerf_ppl=potential%cerf_ppl
       IF (PRESENT(cexp_ppl)) cexp_ppl => potential%cexp_ppl
       IF (PRESENT(ccore_charge)) ccore_charge = potential%ccore_charge
       IF (PRESENT(core_charge_radius))&
            core_charge_radius = potential%core_charge_radius
       IF (PRESENT(ccore_charge_radius))&
            ccore_charge_radius = potential%ccore_charge_radius
       IF (PRESENT(ppl_radius)) ppl_radius=potential%ppl_radius
       IF (PRESENT(nexp_ppl)) nexp_ppl => potential%nexp_ppl
       IF (PRESENT(ppl_present)) ppl_present = (ASSOCIATED(potential%nexp_ppl) )

    ELSE

       CALL stop_program(routineN,moduleN,&
            __LINE__,&
            "The pointer potential is not associated")

    END IF

  END SUBROUTINE get_elp_potential

! *****************************************************************************
!> \brief   Get informations about a GTH potential data set.
!> \author  MK
!> \date    11.01.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_gth_potential(potential,name,alpha_core_charge,&
       alpha_ppl,ccore_charge,cerf_ppl,&
       core_charge_radius,ppl_radius,ppnl_radius,&
       lppnl,lprj_ppnl_max,nexp_ppl,nppnl,&
       nprj_ppnl_max,z,zeff,zeff_correction,&
       ppl_present,ppnl_present,&
       alpha_ppnl,cexp_ppl,elec_conf,nprj_ppnl,cprj,&
       cprj_ppnl,vprj_ppnl,hprj_ppnl)

    TYPE(gth_potential_type), POINTER        :: potential
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: alpha_core_charge, alpha_ppl, &
      ccore_charge, cerf_ppl, core_charge_radius, ppl_radius, ppnl_radius
    INTEGER, INTENT(OUT), OPTIONAL           :: lppnl, lprj_ppnl_max, &
                                                nexp_ppl, nppnl, &
                                                nprj_ppnl_max, z
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: zeff, zeff_correction
    LOGICAL, INTENT(OUT), OPTIONAL           :: ppl_present, ppnl_present
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: alpha_ppnl, cexp_ppl
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: elec_conf, nprj_ppnl
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: cprj, cprj_ppnl, vprj_ppnl
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: hprj_ppnl

    CHARACTER(len=*), PARAMETER :: routineN = 'get_gth_potential', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(potential)) THEN

       IF (PRESENT(name)) name = potential%name
       IF (PRESENT(alpha_core_charge))&
            alpha_core_charge = potential%alpha_core_charge
       IF (PRESENT(alpha_ppl)) alpha_ppl = potential%alpha_ppl
       IF (PRESENT(ccore_charge)) ccore_charge = potential%ccore_charge
       IF (PRESENT(cerf_ppl)) cerf_ppl = potential%cerf_ppl
       IF (PRESENT(core_charge_radius))&
            core_charge_radius = potential%core_charge_radius
       IF (PRESENT(ppl_radius)) ppl_radius = potential%ppl_radius
       IF (PRESENT(ppnl_radius)) ppnl_radius = potential%ppnl_radius
       IF (PRESENT(lppnl)) lppnl = potential%lppnl
       IF (PRESENT(lprj_ppnl_max)) lprj_ppnl_max = potential%lprj_ppnl_max
       IF (PRESENT(nexp_ppl)) nexp_ppl = potential%nexp_ppl
       IF (PRESENT(nppnl)) nppnl = potential%nppnl
       IF (PRESENT(nprj_ppnl_max)) nprj_ppnl_max = potential%nprj_ppnl_max
       IF (PRESENT(z)) z = potential%z
       IF (PRESENT(zeff)) zeff = potential%zeff
       IF (PRESENT(zeff_correction)) zeff_correction = potential%zeff_correction
       IF (PRESENT(ppl_present)) ppl_present = (potential%nexp_ppl > 0)
       IF (PRESENT(ppnl_present)) ppnl_present = (potential%nppnl > 0)
       IF (PRESENT(alpha_ppnl)) alpha_ppnl => potential%alpha_ppnl
       IF (PRESENT(cexp_ppl)) cexp_ppl => potential%cexp_ppl
       IF (PRESENT(elec_conf)) elec_conf => potential%elec_conf
       IF (PRESENT(nprj_ppnl)) nprj_ppnl => potential%nprj_ppnl
       IF (PRESENT(cprj)) cprj => potential%cprj
       IF (PRESENT(cprj_ppnl)) cprj_ppnl => potential%cprj_ppnl
       IF (PRESENT(vprj_ppnl)) vprj_ppnl => potential%vprj_ppnl
       IF (PRESENT(hprj_ppnl)) hprj_ppnl => potential%hprj_ppnl

    ELSE

       CALL stop_program(routineN,moduleN,&
            __LINE__,&
            "The pointer potential is not associated")

    END IF

  END SUBROUTINE get_gth_potential

! *****************************************************************************
!> \brief   Initialise the coefficients of the projectors of the non-local
!>          part of the GTH pseudopotential and the transformation matrices
!>          for Cartesian overlap integrals between the orbital basis
!>          functions and the projector functions.
!> \author  MK
!> \date    16.10.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_cprj_ppnl(potential,error)
    TYPE(gth_potential_type), POINTER        :: potential
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_cprj_ppnl', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cpx, cpy, cpz, cx, cy, cz, &
                                                ico, iprj, iprj_ppnl, l, lp, &
                                                lprj_ppnl, nprj, px, py, pz
    REAL(KIND=dp)                            :: alpha_ppnl, cp

    nprj = 0

    DO l=0,potential%lppnl
       alpha_ppnl = potential%alpha_ppnl(l)
       DO iprj_ppnl=1,potential%nprj_ppnl(l)
          lp = iprj_ppnl - 1
          lprj_ppnl = l + 2*lp
          cp = SQRT(2.0_dp**(2.0_dp*REAL(lprj_ppnl,dp) + 3.5_dp)*&
               alpha_ppnl**(REAL(lprj_ppnl,dp) + 1.5_dp)/&
               (rootpi*dfac(2*lprj_ppnl + 1)))
          potential%cprj_ppnl(iprj_ppnl,l) = cp
          DO cx=0,l
             DO cy=0,l-cx
                cz = l - cx - cy
                iprj = nprj + co(cx,cy,cz)
                DO px=0,lp
                   DO py=0,lp-px
                      pz = lp - px - py
                      cpx = cx + 2*px
                      cpy = cy + 2*py
                      cpz = cz + 2*pz
                      ico = coset(cpx,cpy,cpz)
                      potential%cprj(ico,iprj) = cp*fac(lp)/(fac(px)*fac(py)*fac(pz))
                   END DO
                END DO
             END DO
          END DO
          nprj = nprj + nco(l)
       END DO
    END DO

  END SUBROUTINE init_cprj_ppnl

! *****************************************************************************
!> \brief   Initialise a GTH potential data set structure. 
!> \author  MK
!> \date    27.10.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_gth_potential(potential,error)
    TYPE(gth_potential_type), POINTER        :: potential
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_gth_potential', &
      routineP = moduleN//':'//routineN

    IF (.NOT.ASSOCIATED(potential)) RETURN

    IF (potential%nppnl > 0) THEN

       !     *** Initialise the projector coefficients of the    ***
       !     *** non-local part of the GTH pseudopotential and   ***
       !     *** the transformation matrices "pgf" -> "prj_ppnl" ***

       CALL init_cprj_ppnl(potential,error)

       !     *** Initialise the h(i,j) projector coefficients of ***
       !     *** the non-local part of the GTH pseudopotential   ***

       CALL init_vprj_ppnl(potential,error)

    END IF

  END SUBROUTINE init_gth_potential

! *****************************************************************************
!> \brief   Initialise the h(i,j) projector coefficients of the non-local part
!>          of the GTH pseudopotential.
!> \author  MK 
!> \date    24.10.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_vprj_ppnl(potential,error)
    TYPE(gth_potential_type), POINTER        :: potential
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_vprj_ppnl', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ico, iprj, iprj_ppnl, iso, &
                                                j, jco, jprj, jprj_ppnl, l, &
                                                nprj

    nprj = 0

    DO l=0,potential%lppnl
       DO iprj_ppnl=1,potential%nprj_ppnl(l)
          iprj = nprj + (iprj_ppnl - 1)*nco(l)
          DO jprj_ppnl=1,potential%nprj_ppnl(l)
             jprj = nprj + (jprj_ppnl - 1)*nco(l)
             DO ico=1,nco(l)
                i = iprj + ico
                DO jco=1,nco(l)
                   j = jprj + jco
                   DO iso=1,nso(l)
                      potential%vprj_ppnl(i,j) = potential%vprj_ppnl(i,j) +&
                           orbtramat(l)%slm(iso,ico)*&
                           potential%hprj_ppnl(iprj_ppnl,&
                           jprj_ppnl,l)*&
                           orbtramat(l)%slm(iso,jco)
                   END DO
                END DO
             END DO
          END DO
       END DO
       nprj = nprj + potential%nprj_ppnl(l)*nco(l)
    END DO

  END SUBROUTINE init_vprj_ppnl

! *****************************************************************************
  SUBROUTINE init_all_potential(potential,itype,zeff,zeff_correction,error)

    TYPE(all_potential_type), POINTER        :: potential
    CHARACTER(LEN=*), OPTIONAL               :: itype
    REAL(KIND=dp), OPTIONAL                  :: zeff, zeff_correction
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_all_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dz

    IF (.NOT.ASSOCIATED(potential)) RETURN

    IF ( PRESENT (zeff) ) potential%zeff = zeff
    IF ( PRESENT (zeff_correction) ) potential%zeff_correction = zeff_correction
    dz = potential%z - INT(potential%zeff - potential%zeff_correction)
    SELECT CASE (dz)
    CASE DEFAULT
    CASE (2)
       potential%elec_conf(0) = potential%elec_conf(0) - 2
    CASE (10)
       potential%elec_conf(0) = potential%elec_conf(0) - 4
       potential%elec_conf(1) = potential%elec_conf(1) - 6
    CASE (18)
       potential%elec_conf(0) = potential%elec_conf(0) - 6
       potential%elec_conf(1) = potential%elec_conf(1) - 12
    CASE (28)
       potential%elec_conf(0) = potential%elec_conf(0) - 6
       potential%elec_conf(1) = potential%elec_conf(1) - 12
       potential%elec_conf(2) = potential%elec_conf(2) - 10
    CASE (30)
       potential%elec_conf(0) = potential%elec_conf(0) - 8
       potential%elec_conf(1) = potential%elec_conf(1) - 12
       potential%elec_conf(2) = potential%elec_conf(2) - 10
    CASE (36)
       potential%elec_conf(0) = potential%elec_conf(0) - 8
       potential%elec_conf(1) = potential%elec_conf(1) - 18
       potential%elec_conf(2) = potential%elec_conf(2) - 10
    CASE (46)
       potential%elec_conf(0) = potential%elec_conf(0) - 8
       potential%elec_conf(1) = potential%elec_conf(1) - 18
       potential%elec_conf(2) = potential%elec_conf(2) - 20
    CASE (48)
       potential%elec_conf(0) = potential%elec_conf(0) - 10
       potential%elec_conf(1) = potential%elec_conf(1) - 18
       potential%elec_conf(2) = potential%elec_conf(2) - 20
    CASE (54)
       potential%elec_conf(0) = potential%elec_conf(0) - 10
       potential%elec_conf(1) = potential%elec_conf(1) - 24
       potential%elec_conf(2) = potential%elec_conf(2) - 20
    CASE (68)
       potential%elec_conf(0) = potential%elec_conf(0) - 10
       potential%elec_conf(1) = potential%elec_conf(1) - 24
       potential%elec_conf(2) = potential%elec_conf(2) - 20
       potential%elec_conf(3) = potential%elec_conf(3) - 14
    CASE (78)
       potential%elec_conf(0) = potential%elec_conf(0) - 10
       potential%elec_conf(1) = potential%elec_conf(1) - 24
       potential%elec_conf(2) = potential%elec_conf(2) - 30
       potential%elec_conf(3) = potential%elec_conf(3) - 14
    CASE (80)
       potential%elec_conf(0) = potential%elec_conf(0) - 12
       potential%elec_conf(1) = potential%elec_conf(1) - 24
       potential%elec_conf(2) = potential%elec_conf(2) - 30
       potential%elec_conf(3) = potential%elec_conf(3) - 14
    CASE (86)
       potential%elec_conf(0) = potential%elec_conf(0) - 12
       potential%elec_conf(1) = potential%elec_conf(1) - 30
       potential%elec_conf(2) = potential%elec_conf(2) - 30
       potential%elec_conf(3) = potential%elec_conf(3) - 14
    CASE (100)
       potential%elec_conf(0) = potential%elec_conf(0) - 12
       potential%elec_conf(1) = potential%elec_conf(1) - 30
       potential%elec_conf(2) = potential%elec_conf(2) - 30
       potential%elec_conf(3) = potential%elec_conf(3) - 28
    END SELECT

    IF ( PRESENT (itype) ) THEN
       IF ( itype == "BARE" ) THEN
          potential%description(1) = "Bare Coulomb Potential"
          IF ( dz > 0 ) THEN
             potential%description(2) = "Valence charge only"
          ELSE
             potential%description(2) = "Full atomic charge"
          END IF
       ENDIF
    END IF

  END SUBROUTINE init_all_potential

! *****************************************************************************
!> \brief   Read an atomic all-electron potential data set.
!> \author  MK
!> \date    14.05.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE read_all_potential(element_symbol,potential_name,potential,zeff_correction,&
       para_env, dft_section, potential_section, error)

    CHARACTER(LEN=*), INTENT(IN)             :: element_symbol, potential_name
    TYPE(all_potential_type), POINTER        :: potential
    REAL(KIND=dp)                            :: zeff_correction
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: dft_section, potential_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_all_potential', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=240)                       :: line
    CHARACTER(LEN=242)                       :: line2
    CHARACTER(len=5*default_string_length)   :: line_att
    CHARACTER(len=default_path_length)       :: potential_file_name
    CHARACTER(LEN=LEN(element_symbol))       :: symbol
    CHARACTER(LEN=LEN(element_symbol)+2)     :: symbol2
    CHARACTER(LEN=LEN(potential_name))       :: apname
    CHARACTER(LEN=LEN(potential_name)+2)     :: apname2
    INTEGER                                  :: irep, l, stat, strlen1, &
                                                strlen2
    INTEGER, DIMENSION(:), POINTER           :: elec_conf
    LOGICAL                                  :: failure, found, is_ok, match, &
                                                read_from_input
    REAL(KIND=dp)                            :: alpha, r
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_sll_val_type), POINTER           :: list
    TYPE(val_type), POINTER                  :: val

    line2 = ""
    symbol2 = ""
    apname2 = ""
    NULLIFY(parser)
    failure = .FALSE.
    CALL cite_reference(Krack2000)

    potential%name = potential_name
    read_from_input = .FALSE.
    CALL section_vals_get(potential_section,explicit=read_from_input, error=error)
    IF (.NOT.read_from_input) THEN
       CALL section_vals_val_get(dft_section,"POTENTIAL_FILE_NAME",c_val=potential_file_name,&
            error=error)
       CALL parser_create(parser,potential_file_name,para_env=para_env,error=error)
    END IF

    !   *** Search for the requested potential in the potential file   ***
    !   *** until the potential is found or the end of file is reached ***

    apname = potential_name
    symbol = element_symbol
    irep   = 0
    search_loop: DO
       IF (read_from_input) THEN
          NULLIFY(list, val)
          found = .TRUE.
          CALL section_vals_list_get(potential_section,"_DEFAULT_KEYWORD_",list=list,error=error)
       ELSE
          CALL parser_search_string(parser,TRIM(apname),.TRUE.,found,line,error=error)
       END IF
       IF (found) THEN
          CALL uppercase(symbol)
          CALL uppercase(apname)

          IF (read_from_input) THEN
             match = .TRUE.
          ELSE
             !       *** Check both the element symbol and the atomic potential name ***
             match = .FALSE.
             CALL uppercase(line)
             line2 = " "//line//" "
             symbol2 = " "//TRIM(symbol)//" "
             apname2 = " "//TRIM(apname)//" "
             strlen1 = LEN_TRIM(symbol2) + 1
             strlen2 = LEN_TRIM(apname2) + 1

             IF ( (INDEX(line2,symbol2(:strlen1)) > 0).AND.&
                  (INDEX(line2,apname2(:strlen2)) > 0)) match = .TRUE.
          END IF
          IF (match) THEN
             !         *** Read the electronic configuration ***
             NULLIFY (elec_conf)
             l = 0
             CALL reallocate(elec_conf,0,l)
             IF (read_from_input) THEN
                is_ok=cp_sll_val_next(list,val,error=error)
                IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
                     __LINE__,&
                     "Error reading the Potential from input file!!")
                CALL val_get(val,c_val=line_att,error=error)
                READ(line_att,*)elec_conf(l)
                CALL remove_word(line_att)
                DO WHILE (LEN_TRIM(line_att) /= 0)
                   l = l + 1
                   CALL reallocate(elec_conf,0,l)
                   READ(line_att,*)elec_conf(l)
                   CALL remove_word(line_att)
                END DO
             ELSE
                CALL parser_get_object(parser,elec_conf(l),newline=.TRUE.,error=error)
                DO WHILE (parser_test_next_token(parser,error=error) == "INT")
                   l = l + 1
                   CALL reallocate(elec_conf,0,l)
                   CALL parser_get_object(parser,elec_conf(l),error=error)
                END DO
                irep = irep + 1
                WRITE(line_att,'(100(1X,I0))')elec_conf
                CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
                     c_val=TRIM(line_att), error=error)              
             END IF

             CALL reallocate(potential%elec_conf,0,l)
             potential%elec_conf(:) = elec_conf(:)

             potential%zeff_correction = zeff_correction
             potential%zeff = REAL(SUM(elec_conf),dp)+zeff_correction

             DEALLOCATE (elec_conf,STAT=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

             !         *** Read r(loc) to define the exponent of the core charge    ***
             !         *** distribution and calculate the corresponding coefficient ***

             IF (read_from_input) THEN
                is_ok=cp_sll_val_next(list,val,error=error)
                IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
                     __LINE__,&
                     "Error reading the Potential from input file!!")
                CALL val_get(val,c_val=line_att,error=error)
                READ(line_att,*)r                          
             ELSE
                CALL parser_get_object(parser,r,newline=.TRUE.,error=error)
                irep = irep + 1
                WRITE(line_att,'(E24.16)')r
                CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
                     c_val=TRIM(line_att), error=error)
             END IF
             alpha = 1.0_dp/(2.0_dp*r**2)

             potential%alpha_core_charge = alpha
             potential%ccore_charge = potential%zeff*SQRT((alpha/pi)**3)

             EXIT search_loop
          END IF
       ELSE
          !       *** Stop program, if the end of file is reached ***
          CALL stop_program(routineN,moduleN,&
               __LINE__,&
               "The requested atomic potential <"//TRIM(potential_name)//&
               "> was not found in the potential file <"//TRIM(potential_file_name)//">")
       END IF
    END DO search_loop

    IF (.NOT.read_from_input) THEN
       ! Dump the potential info the in potential section
       IF (match) THEN
          irep = irep + 1
          WRITE(line_att,'(A)')"         # Potential name: "//apname2(:strlen2)//" for symbol: "//symbol2(:strlen1)
          CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
               c_val=TRIM(line_att), error=error)
          irep = irep + 1
          WRITE(line_att,'(A)')"         # Potential read from the potential filename: "//TRIM(potential_file_name)
          CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
               c_val=TRIM(line_att), error=error)
       END IF
       CALL parser_release(parser,error=error)
    END IF
  END SUBROUTINE read_all_potential

! *****************************************************************************
!> \brief   Read an effective charge potential data set.
!> \author  GT
!> \date    22-09-2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE read_elp_potential(element_symbol,potential_name,potential,para_env,&
       dft_section,potential_section,error)
    CHARACTER(LEN=*), INTENT(IN)             :: element_symbol, potential_name
    TYPE(elp_potential_type), POINTER        :: potential
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: dft_section, potential_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_elp_potential', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=240)                       :: line
    CHARACTER(LEN=242)                       :: line2
    CHARACTER(LEN=5*default_string_length)   :: line_att
    CHARACTER(LEN=default_path_length)       :: potential_file_name
    CHARACTER(LEN=LEN(element_symbol))       :: symbol
    CHARACTER(LEN=LEN(element_symbol)+2)     :: symbol2
    CHARACTER(LEN=LEN(potential_name))       :: apname
    CHARACTER(LEN=LEN(potential_name)+2)     :: apname2
    INTEGER                                  :: i, irep, istr, j, maxlppl, n, &
                                                nalpha_ppl, nmax, strlen1, &
                                                strlen2
    LOGICAL                                  :: found, is_ok, match, &
                                                read_from_input
    REAL(KIND=dp)                            :: alpha, qeff
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_sll_val_type), POINTER           :: list
    TYPE(val_type), POINTER                  :: val

    line2 = ""
    symbol2 = ""
    apname2 = ""
    nalpha_ppl = 0
    nmax = 0
    NULLIFY(parser)

    potential%name = potential_name
    read_from_input = .FALSE.
    CALL section_vals_get(potential_section,explicit=read_from_input, error=error)
    IF (.NOT.read_from_input) THEN
       CALL section_vals_val_get(dft_section,"POTENTIAL_FILE_NAME",c_val=potential_file_name,&
            error=error)
       CALL parser_create(parser,potential_file_name,para_env=para_env,error=error)
    END IF

    !   *** Search for the requested potential in the potential file   ***
    !   *** until the potential is found or the end of file is reached ***

    apname = potential_name
    symbol = element_symbol
    irep   = 0
    search_loop: DO
       IF (read_from_input) THEN
          NULLIFY(list, val)
          found = .TRUE.
          CALL section_vals_list_get(potential_section,"_DEFAULT_KEYWORD_",list=list,error=error)
       ELSE
          CALL parser_search_string(parser,TRIM(apname),.TRUE.,found,line,error=error)
       END IF
       IF (found) THEN
          CALL uppercase(symbol)
          CALL uppercase(apname)
          IF (read_from_input) THEN
             match = .TRUE.
          ELSE
             !       *** Check both the element symbol and the atomic potential name ***
             match = .FALSE.
             CALL uppercase(line)
             line2 = " "//line//" "
             symbol2 = " "//TRIM(symbol)//" "
             apname2 = " "//TRIM(apname)//" "
             strlen1 = LEN_TRIM(symbol2) + 1
             strlen2 = LEN_TRIM(apname2) + 1

             IF ( (INDEX(line2,symbol2(:strlen1)) > 0).AND.&
                  (INDEX(line2,apname2(:strlen2)) > 0) ) match = .TRUE.
          END IF
          IF (match) THEN
             !         *** Read the effective charge qeff ***
             IF (read_from_input) THEN
                is_ok=cp_sll_val_next(list,val,error=error)
                IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
                     __LINE__,&
                     "Error reading the Potential from input file!!")
                CALL val_get(val,c_val=line_att,error=error)
                READ(line_att,*)qeff                
             ELSE
                CALL parser_get_object(parser,qeff,newline=.TRUE.,error=error)
                irep = irep + 1
                WRITE(line_att,'(E24.16)')qeff
                CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
                     c_val=TRIM(line_att), error=error)                
             END IF
             potential%qeff = qeff

             !         *** read alpha and calculate the corresponding coefficient ***
             IF (read_from_input) THEN
                is_ok=cp_sll_val_next(list,val,error=error)
                IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
                     __LINE__,&
                     "Error reading the Potential from input file!!")
                CALL val_get(val,c_val=line_att,error=error)
                READ(line_att,*)alpha
             ELSE
                CALL parser_get_object(parser,alpha,newline=.TRUE.,error=error)
                irep = irep + 1
                WRITE(line_att,'(E24.16)')alpha
                CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
                     c_val=TRIM(line_att), error=error) 
             END IF
             potential%alpha_core_charge = alpha
             potential%ccore_charge = potential%qeff*SQRT((alpha/pi)**3)
             potential%cerf_ppl = potential%qeff*SQRT((alpha/pi)**3)

             IF (read_from_input) THEN
                is_ok=cp_sll_val_next(list,val,error=error)
                IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
                     __LINE__,&
                     "Error reading the Potential from input file!!")
                CALL val_get(val,c_val=line_att,error=error)
                READ(line_att,*)nalpha_ppl
             ELSE
                CALL parser_get_object(parser,nalpha_ppl,newline=.TRUE.,error=error)
                irep = irep + 1
                WRITE(line_att,'(1X,I0)')nalpha_ppl
                CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
                     c_val=TRIM(line_att), error=error) 
             END IF
             potential%nalpha_ppl = nalpha_ppl

             IF (nalpha_ppl > 0) THEN
                CALL reallocate(potential%alpha_ppl,1,nalpha_ppl)
                CALL reallocate(potential%nexp_ppl,1,nalpha_ppl)
                DO i=1,nalpha_ppl
                   IF (read_from_input) THEN
                      is_ok=cp_sll_val_next(list,val,error=error)
                      IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
                           __LINE__,&
                           "Error reading the Potential from input file!!")
                      CALL val_get(val,c_val=line_att,error=error)
                      READ(line_att,*)potential%alpha_ppl(i)
                      CALL remove_word(line_att)
                      READ(line_att,*)n
                      CALL remove_word(line_att)
                   ELSE
                      line_att = ""
                      CALL parser_get_object(parser,potential%alpha_ppl(i),newline=.TRUE.,error=error)
                      CALL parser_get_object(parser,n,error=error)
                      WRITE(line_att,'(E24.16,1X,I0)')potential%alpha_ppl(i),n
                   END IF
                   potential%nexp_ppl(i)=n
                   IF (n > nmax) THEN
                      nmax = n
                      CALL reallocate(potential%cexp_ppl,1,nalpha_ppl,1,nmax)
                   END IF
                   DO j=1,n
                      IF (read_from_input) THEN
                         READ(line_att,*)potential%cexp_ppl(i,j)
                         CALL remove_word(line_att)
                      ELSE
                         CALL parser_get_object(parser,potential%cexp_ppl(i,j),error=error)
                         istr=LEN_TRIM(line_att)+1
                         WRITE(line_att(istr:),'(100E24.16)')potential%cexp_ppl(i,j)
                      END IF
                      potential%cexp_ppl(i,j) = (2.0_dp*potential%alpha_ppl(i))**(j - 1)*&
                           potential%cexp_ppl(i,j)
                   END DO
                   IF (.NOT.read_from_input) THEN
                      irep = irep + 1
                      CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
                           c_val=TRIM(line_att), error=error)
                   ELSE
                      IF (LEN_TRIM(line_att)/=0) CALL stop_program(routineN,moduleN,&
                           __LINE__,&
                           "Error reading the Potential from input file!!")
                   END IF
                END DO

                !           *** Check if alpha_core_charge = alpha_ppl(1) ***
                IF (potential%alpha_core_charge /= potential%alpha_ppl(1)) THEN
                   CALL stop_program(routineN,moduleN,&
                        __LINE__,&
                        "exponent of the core charge and "//&
                        "alpha_ppl(1) must be equal")
                END IF

                maxlppl = 2*(nmax - 1)
                IF (maxlppl > -1)  CALL init_orbital_pointers(maxlppl)
             END IF
             EXIT search_loop
          END IF
       ELSE
          !       *** Stop program, if the end of file is reached ***
          CALL stop_program(routineN,moduleN,&
               __LINE__,&
               "The requested atomic potential <"//TRIM(potential_name)//&
               "> was not found in the potential file <"//TRIM(potential_file_name)//">",para_env)
       END IF
    END DO search_loop
    IF (.NOT.read_from_input) THEN
       ! Dump the potential info the in potential section
       IF (match) THEN
          irep = irep + 1
          WRITE(line_att,'(A)')"         # Potential name: "//apname2(:strlen2)//" for symbol: "//symbol2(:strlen1)
          CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
               c_val=TRIM(line_att), error=error)
          irep = irep + 1
          WRITE(line_att,'(A)')"         # Potential read from the potential filename: "//TRIM(potential_file_name)
          CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
               c_val=TRIM(line_att), error=error)
       END IF
       CALL parser_release(parser,error=error)
    END IF
  END SUBROUTINE read_elp_potential

! *****************************************************************************
!> \brief   Read an atomic GTH potential data set.
!> \author  MK
!> \date    14.05.2000
!> \par  Literature
!>         - S. Goedecker, M. Teter and J. Hutter,
!>                Phys. Rev. B 54, 1703 (1996)
!>         - C. Hartwigsen, S. Goedecker and J. Hutter,
!>                Phys. Rev. B 58, 3641 (1998) 
!> \version 1.0
! *****************************************************************************
  SUBROUTINE read_gth_potential(element_symbol,potential_name,potential,zeff_correction,&
       para_env, dft_section, potential_section, error)

    CHARACTER(LEN=*), INTENT(IN)             :: element_symbol, potential_name
    TYPE(gth_potential_type), POINTER        :: potential
    REAL(KIND=dp), INTENT(IN)                :: zeff_correction
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: dft_section, potential_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_gth_potential', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=240)                       :: line
    CHARACTER(LEN=242)                       :: line2
    CHARACTER(len=5*default_string_length)   :: line_att
    CHARACTER(len=default_path_length)       :: potential_file_name
    CHARACTER(LEN=LEN(element_symbol))       :: symbol
    CHARACTER(LEN=LEN(element_symbol)+2)     :: symbol2
    CHARACTER(LEN=LEN(potential_name))       :: apname
    CHARACTER(LEN=LEN(potential_name)+2)     :: apname2
    INTEGER :: i, irep, istr, j, l, lppnl, lprj_ppnl_max, maxlppl, n, nppnl, &
      nprj_ppnl, nprj_ppnl_max, stat, strlen1, strlen2
    INTEGER, DIMENSION(:), POINTER           :: elec_conf
    LOGICAL                                  :: failure, found, is_ok, match, &
                                                read_from_input
    REAL(KIND=dp)                            :: alpha, ci, r, rc2
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: hprj_ppnl
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_sll_val_type), POINTER           :: list
    TYPE(val_type), POINTER                  :: val

    line2 = ""
    symbol2 = ""
    apname2 = ""
    NULLIFY(parser)
    failure = .FALSE.
    CALL cite_reference(Goedecker1996)
    CALL cite_reference(Hartwigsen1998)
    CALL cite_reference(Krack2005)

    potential%name = potential_name
    read_from_input = .FALSE.
    CALL section_vals_get(potential_section,explicit=read_from_input, error=error)
    IF (.NOT.read_from_input) THEN
       CALL section_vals_val_get(dft_section,"POTENTIAL_FILE_NAME",c_val=potential_file_name,&
            error=error)
       CALL parser_create(parser,potential_file_name,para_env=para_env,error=error)
    END IF

    !   *** Search for the requested potential in the potential file   ***
    !   *** until the potential is found or the end of file is reached ***

    apname = potential_name
    symbol = element_symbol
    irep   = 0
    search_loop: DO
       IF (read_from_input) THEN
          NULLIFY(list, val)
          found = .TRUE.
          CALL section_vals_list_get(potential_section,"_DEFAULT_KEYWORD_",list=list,error=error)
       ELSE
          CALL parser_search_string(parser,TRIM(apname),.TRUE.,found,line,error=error)
       END IF
       IF (found) THEN
          CALL uppercase(symbol)
          CALL uppercase(apname)
          IF (read_from_input) THEN
             match = .TRUE.
          ELSE
             !       *** Check both the element symbol and the atomic potential name ***
             match = .FALSE.
             CALL uppercase(line)
             line2 = " "//line//" "
             symbol2 = " "//TRIM(symbol)//" "
             apname2 = " "//TRIM(apname)//" "
             strlen1 = LEN_TRIM(symbol2) + 1
             strlen2 = LEN_TRIM(apname2) + 1

             IF ( (INDEX(line2,symbol2(:strlen1)) > 0).AND.&
                  (INDEX(line2,apname2(:strlen2)) > 0) ) match = .TRUE.
          END IF
          IF (match) THEN
             !         *** Read the electronic configuration ***
             NULLIFY (elec_conf)
             l = 0
             CALL reallocate(elec_conf,0,l)
             IF (read_from_input) THEN
                is_ok=cp_sll_val_next(list,val,error=error)
                IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
                     __LINE__,&
                     "Error reading the Potential from input file!!")
                CALL val_get(val,c_val=line_att,error=error)
                READ(line_att,*)elec_conf(l)
                CALL remove_word(line_att)
                DO WHILE (LEN_TRIM(line_att) /= 0)
                   l = l + 1
                   CALL reallocate(elec_conf,0,l)
                   READ(line_att,*)elec_conf(l)
                   CALL remove_word(line_att)
                END DO
             ELSE            
                CALL parser_get_object(parser,elec_conf(l),newline=.TRUE.,error=error)
                DO WHILE (parser_test_next_token(parser,error=error) == "INT")
                   l = l + 1
                   CALL reallocate(elec_conf,0,l)
                   CALL parser_get_object(parser,elec_conf(l),error=error)
                END DO
                irep = irep + 1
                WRITE(line_att,'(100(1X,I0))')elec_conf
                CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
                     c_val=TRIM(line_att), error=error)              
             END IF

             CALL reallocate(potential%elec_conf,0,l)
             potential%elec_conf(:) = elec_conf(:)

             potential%zeff_correction = zeff_correction
             potential%zeff = REAL(SUM(elec_conf),dp)+zeff_correction

             DEALLOCATE (elec_conf,STAT=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

             !         *** Read r(loc) to define the exponent of the core charge    ***
             !         *** distribution and calculate the corresponding coefficient ***

             IF (read_from_input) THEN
                is_ok=cp_sll_val_next(list,val,error=error)
                IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
                     __LINE__,&
                     "Error reading the Potential from input file!!")
                CALL val_get(val,c_val=line_att,error=error)
                READ(line_att,*)r
                CALL remove_word(line_att)
             ELSE
                line_att = ""
                CALL parser_get_object(parser,r,newline=.TRUE.,error=error) 
                istr=LEN_TRIM(line_att)+1
                WRITE(line_att(istr:),'(E24.16)')r
             END IF
             alpha = 1.0_dp/(2.0_dp*r**2)

             potential%alpha_core_charge = alpha
             potential%ccore_charge = potential%zeff*SQRT((alpha/pi)**3)

             potential%alpha_ppl = alpha
             potential%cerf_ppl = potential%zeff*SQRT((alpha/pi)**3)

             !         *** Read the parameters for the local part ***
             !         *** of the GTH pseudopotential (ppl)       ***

             IF (read_from_input) THEN
                READ(line_att,*)n
                CALL remove_word(line_att)
             ELSE
                CALL parser_get_object(parser,n,error=error)
                istr=LEN_TRIM(line_att)+1
                WRITE(line_att(istr:),'(1X,I0)')n
             END IF
             potential%nexp_ppl = n
             CALL reallocate(potential%cexp_ppl,1,n)

             DO i=1,n
                IF (read_from_input) THEN
                   READ(line_att,*) ci
                   CALL remove_word(line_att)
                ELSE
                   CALL parser_get_object(parser,ci,error=error)
                   istr=LEN_TRIM(line_att)+1
                   WRITE(line_att(istr:),'(E24.16)') ci
                END IF
                rc2 = (2.0_dp*potential%alpha_ppl)
                potential%cexp_ppl(i) = rc2**(i - 1)*ci
             END DO

             IF (.NOT.read_from_input) THEN
                irep = irep + 1
                CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
                     c_val=TRIM(line_att), error=error)
                line_att = ""
             ELSE
                IF (LEN_TRIM(line_att)/=0) CALL stop_program(routineN,moduleN,&
                     __LINE__,&
                     "Error reading the Potential from input file!!")
             END IF
             maxlppl = 2*(n - 1)

             IF (maxlppl > -1) CALL init_orbital_pointers(maxlppl)

             !         *** Read the parameters for the non-local  ***
             !         *** part of the GTH pseudopotential (ppnl) ***
             IF (read_from_input) THEN
                is_ok=cp_sll_val_next(list,val,error=error)
                IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
                     __LINE__,&
                     "Error reading the Potential from input file!!")
                CALL val_get(val,c_val=line_att,error=error)
                READ(line_att,*)n
                CALL remove_word(line_att)               
             ELSE
                CALL parser_get_object(parser,n,newline=.TRUE.,error=error)
                irep = irep + 1
                WRITE(line_att,'(1X,I0)')n
                CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
                     c_val=TRIM(line_att), error=error)
             END IF
             potential%lppnl = n - 1
             potential%nppnl = 0

             potential%lprj_ppnl_max = n - 1
             potential%nprj_ppnl_max = 0

             IF (n > 0) THEN

                lppnl = potential%lppnl
                nppnl = potential%nppnl

                CALL init_orbital_pointers(lppnl)

                NULLIFY (hprj_ppnl)

                !           *** Load the parameter for n non-local projectors ***

                CALL reallocate(potential%alpha_ppnl,0,lppnl)
                CALL reallocate(potential%nprj_ppnl,0,lppnl)

                lprj_ppnl_max = -1
                nprj_ppnl_max = 0

                DO l=0,lppnl
                   IF (read_from_input) THEN
                      is_ok=cp_sll_val_next(list,val,error=error)
                      IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
                           __LINE__,&
                           "Error reading the Potential from input file!!")
                      CALL val_get(val,c_val=line_att,error=error)
                      READ(line_att,*)r
                      CALL remove_word(line_att)
                      READ(line_att,*)nprj_ppnl
                      CALL remove_word(line_att)
                   ELSE
                      line_att =""
                      CALL parser_get_object(parser,r,newline=.TRUE.,error=error)
                      CALL parser_get_object(parser,nprj_ppnl,error=error)
                      istr=LEN_TRIM(line_att)+1
                      WRITE(line_att(istr:),'(E24.16,1X,I0)')r, nprj_ppnl
                   END IF
                   IF (r==0.0_dp.AND.nprj_ppnl/=0) THEN 
                      CALL stop_program(routineN,moduleN,&
                           __LINE__,&
                           "An error was detected in the atomic potential <"//&
                           TRIM(potential_name)//&
                           "> potential file <"//&
                           TRIM(potential_file_name)//">",para_env)
                   END IF
                   potential%alpha_ppnl(l) = 0.0_dp
                   IF (r/=0.0_dp.AND.n/=0)  potential%alpha_ppnl(l) = 1.0_dp/(2.0_dp*r**2)
                   potential%nprj_ppnl(l) = nprj_ppnl
                   nppnl = nppnl + nprj_ppnl*nco(l)
                   IF (nprj_ppnl > nprj_ppnl_max) THEN
                      nprj_ppnl_max = nprj_ppnl
                      CALL reallocate(hprj_ppnl,1,nprj_ppnl_max,&
                           1,nprj_ppnl_max,&
                           0,lppnl)
                   END IF
                   DO i=1,nprj_ppnl
                      IF (i == 1) THEN
                         IF (read_from_input) THEN
                            READ(line_att,*)hprj_ppnl(i,i,l)
                            CALL remove_word(line_att) 
                         ELSE
                            CALL parser_get_object(parser,hprj_ppnl(i,i,l),error=error)
                            istr=LEN_TRIM(line_att)+1
                            WRITE(line_att(istr:),'(E24.16)')hprj_ppnl(i,i,l)
                         END IF
                      ELSE
                         IF (read_from_input) THEN
                            IF (LEN_TRIM(line_att)/=0) CALL stop_program(routineN,moduleN,&
                                 __LINE__,&
                                 "Error reading the Potential from input file!!")
                            is_ok=cp_sll_val_next(list,val,error=error)
                            IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
                                 __LINE__,&
                                 "Error reading the Potential from input file!!")
                            CALL val_get(val,c_val=line_att,error=error)
                            READ(line_att,*)hprj_ppnl(i,i,l)
                            CALL remove_word(line_att)
                         ELSE
                            irep = irep + 1
                            CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
                                 c_val=TRIM(line_att), error=error)                 
                            line_att =""
                            CALL parser_get_object(parser,hprj_ppnl(i,i,l),newline=.TRUE.,error=error)  
                            istr=LEN_TRIM(line_att)+1
                            WRITE(line_att(istr:),'(E24.16)')hprj_ppnl(i,i,l)
                         END IF
                      END IF
                      DO j=i+1,nprj_ppnl
                         IF (read_from_input) THEN
                            READ(line_att,*)hprj_ppnl(i,j,l)
                            CALL remove_word(line_att)
                         ELSE
                            CALL parser_get_object(parser,hprj_ppnl(i,j,l),error=error)
                            istr=LEN_TRIM(line_att)+1
                            WRITE(line_att(istr:),'(E24.16)')hprj_ppnl(i,j,l)
                         END IF
                      END DO
                   END DO
                   IF (.NOT.read_from_input) THEN
                      irep = irep + 1
                      CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
                           c_val=TRIM(line_att), error=error)                 
                      line_att =""
                   ELSE
                      IF (LEN_TRIM(line_att)/=0) CALL stop_program(routineN,moduleN,&
                           __LINE__,&
                           "Error reading the Potential from input file!!")
                   END IF
                   IF (nprj_ppnl > 1) THEN
                      CALL symmetrize_matrix(hprj_ppnl(:,:,l),"upper_to_lower")
                   END IF
                   lprj_ppnl_max = MAX(lprj_ppnl_max,l + 2*(nprj_ppnl - 1))
                END DO

                potential%nppnl = nppnl
                CALL init_orbital_pointers(lprj_ppnl_max)

                potential%lprj_ppnl_max = lprj_ppnl_max
                potential%nprj_ppnl_max = nprj_ppnl_max
                CALL reallocate(potential%hprj_ppnl,1,nprj_ppnl_max,&
                     1,nprj_ppnl_max,&
                     0,lppnl)
                potential%hprj_ppnl(:,:,:) = hprj_ppnl(:,:,:)

                CALL reallocate(potential%cprj,1,ncoset(lprj_ppnl_max),1,nppnl)
                CALL reallocate(potential%cprj_ppnl,1,nprj_ppnl_max,0,lppnl)
                CALL reallocate(potential%vprj_ppnl,1,nppnl,1,nppnl)

                DEALLOCATE (hprj_ppnl,STAT=stat)
                CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             EXIT search_loop
          END IF
       ELSE
          !       *** Stop program, if the end of file is reached ***
          CALL stop_program(routineN,moduleN,&
               __LINE__,&
               "The requested atomic potential <"//&
               TRIM(potential_name)//&
               "> was not found in the potential file <"//&
               TRIM(potential_file_name)//">",para_env)
       END IF
    END DO search_loop
    IF (.NOT.read_from_input) THEN
       ! Dump the potential info the in potential section
       IF (match) THEN
          irep = irep + 1
          WRITE(line_att,'(A)')"         # Potential name: "//apname2(:strlen2)//" for symbol: "//symbol2(:strlen1)             
          CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
               c_val=TRIM(line_att), error=error)
          irep = irep + 1
          WRITE(line_att,'(A)')"         # Potential read from the potential filename: "//TRIM(potential_file_name)
          CALL section_vals_val_set(potential_section,"_DEFAULT_KEYWORD_",i_rep_val=irep,&
               c_val=TRIM(line_att), error=error)
       END IF
       CALL parser_release(parser,error=error)
    END IF
  END SUBROUTINE read_gth_potential

! *****************************************************************************
  SUBROUTINE set_default_all_potential(potential,z,zeff_correction,error)

    TYPE(all_potential_type), POINTER        :: potential
    INTEGER, INTENT(IN)                      :: z
    REAL(KIND=dp)                            :: zeff_correction
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_default_all_potential', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name
    INTEGER, DIMENSION(:), POINTER           :: elec_conf
    REAL(KIND=dp)                            :: alpha, alpha_core_charge, &
                                                ccore_charge, &
                                                core_charge_radius, r, zeff

    ALLOCATE ( elec_conf(0:3) )
    elec_conf(0:3) = ptable(z)%e_conv(0:3)
    zeff = REAL(SUM ( elec_conf ),dp)+zeff_correction
    name = ptable(z)%name

    r = ptable(z)%covalent_radius * 0.5_dp
    r = MAX ( r, 0.2_dp )
    r = MIN ( r, 1.0_dp )
    alpha = 1.0_dp/(2.0_dp*r**2)

    core_charge_radius = r
    alpha_core_charge = alpha
    ccore_charge = zeff*SQRT((alpha/pi)**3)

    CALL set_all_potential(potential,&
         name=name,&
         alpha_core_charge=alpha_core_charge,&
         ccore_charge=ccore_charge,&
         core_charge_radius=core_charge_radius,&
         z=z,&
         zeff=zeff,&
         zeff_correction=zeff_correction,&
         elec_conf=elec_conf)

    DEALLOCATE ( elec_conf )

  END SUBROUTINE set_default_all_potential

! *****************************************************************************
!> \brief   Set the components of an all-electron potential data set.
!> \author  MK
!> \date    11.01.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_all_potential(potential,name,alpha_core_charge,&
       ccore_charge,core_charge_radius,z,zeff,&
       zeff_correction,elec_conf)

    TYPE(all_potential_type), POINTER        :: potential
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: name
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: alpha_core_charge, &
                                                ccore_charge, &
                                                core_charge_radius
    INTEGER, INTENT(IN), OPTIONAL            :: z
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: zeff, zeff_correction
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: elec_conf

    CHARACTER(len=*), PARAMETER :: routineN = 'set_all_potential', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(potential)) THEN

       IF (PRESENT(name)) potential%name = name
       IF (PRESENT(alpha_core_charge))&
            potential%alpha_core_charge = alpha_core_charge
       IF (PRESENT(ccore_charge)) potential%ccore_charge = ccore_charge
       IF (PRESENT(core_charge_radius))&
            potential%core_charge_radius = core_charge_radius
       IF (PRESENT(z)) potential%z = z
       IF (PRESENT(zeff)) potential%zeff = zeff
       IF (PRESENT(zeff_correction)) potential%zeff_correction = zeff_correction
       IF (PRESENT(elec_conf)) THEN
          IF ( .NOT. ASSOCIATED( potential%elec_conf ) ) THEN
             CALL reallocate(potential%elec_conf,0,SIZE(elec_conf)-1)
          END IF
          potential%elec_conf(:) = elec_conf(:)
       END IF

    ELSE

       CALL stop_program(routineN,moduleN,&
            __LINE__,&
            "The pointer potential is not associated")

    END IF

  END SUBROUTINE set_all_potential

! *****************************************************************************
!> \brief   Set the effective charge in the elp potential data set.
!> \author  GT
!> \date    22-09-2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_elp_potential(potential,name,alpha_core_charge,ccore_charge,&
       core_charge_radius,ccore_charge_radius,cexp_ppl,cerf_ppl,qeff,a_pol,&
       ppl_radius,alpha_ppl,nexp_ppl,nalpha_ppl,description)
    TYPE(elp_potential_type), POINTER        :: potential
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: name
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: alpha_core_charge, &
                                                ccore_charge, &
                                                core_charge_radius, &
                                                ccore_charge_radius
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: cexp_ppl
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: cerf_ppl, qeff, a_pol, &
                                                ppl_radius
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: alpha_ppl
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: nexp_ppl
    INTEGER, OPTIONAL                        :: nalpha_ppl
    CHARACTER(LEN=default_string_length), &
      DIMENSION(2), OPTIONAL                 :: description

    CHARACTER(len=*), PARAMETER :: routineN = 'set_elp_potential', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(potential)) THEN

       IF (PRESENT(name)) potential%name = name
       IF (PRESENT(description)) potential%description = description
       IF (PRESENT(nalpha_ppl)) potential%nalpha_ppl = nalpha_ppl
       IF (PRESENT(alpha_core_charge))&
            potential%alpha_core_charge = alpha_core_charge
       IF (PRESENT(ccore_charge)) potential%ccore_charge = ccore_charge
       IF (PRESENT(core_charge_radius))&
            potential%core_charge_radius = core_charge_radius
       IF (PRESENT(ccore_charge_radius))&
            potential%ccore_charge_radius = ccore_charge_radius
       IF (PRESENT(qeff)) potential%qeff = qeff
       IF (PRESENT(a_pol)) potential%a_pol = a_pol
       IF (PRESENT(alpha_ppl)) potential%alpha_ppl => alpha_ppl
       IF (PRESENT(cexp_ppl)) potential%cexp_ppl => cexp_ppl
       IF (PRESENT(nexp_ppl)) potential%nexp_ppl => nexp_ppl
       IF (PRESENT(cerf_ppl)) potential%cerf_ppl = cerf_ppl
       IF (PRESENT(ppl_radius)) potential%ppl_radius = ppl_radius

    ELSE

       CALL stop_program(routineN,moduleN,&
            __LINE__,&
            "The pointer potential is not associated")

    END IF

  END SUBROUTINE set_elp_potential

! *****************************************************************************
!> \brief   Get the components of a GTH potential data set.
!> \author  MK
!> \date    11.01.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_gth_potential(potential,name,alpha_core_charge,alpha_ppl,&
       ccore_charge,cerf_ppl,core_charge_radius,&
       ppl_radius,ppnl_radius,lppnl,lprj_ppnl_max,&
       nexp_ppl,nppnl,nprj_ppnl_max,z,zeff,zeff_correction,&
       alpha_ppnl,cexp_ppl,elec_conf,nprj_ppnl,cprj,cprj_ppnl,&
       vprj_ppnl,hprj_ppnl)
    TYPE(gth_potential_type), POINTER        :: potential
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: name
    REAL(KIND=dp), INTENT(IN), OPTIONAL :: alpha_core_charge, alpha_ppl, &
      ccore_charge, cerf_ppl, core_charge_radius, ppl_radius, ppnl_radius
    INTEGER, INTENT(IN), OPTIONAL            :: lppnl, lprj_ppnl_max, &
                                                nexp_ppl, nppnl, &
                                                nprj_ppnl_max, z
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: zeff, zeff_correction
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: alpha_ppnl, cexp_ppl
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: elec_conf, nprj_ppnl
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: cprj, cprj_ppnl, vprj_ppnl
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: hprj_ppnl

    CHARACTER(len=*), PARAMETER :: routineN = 'set_gth_potential', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(potential)) THEN

       IF (PRESENT(name)) potential%name = name
       IF (PRESENT(alpha_core_charge))&
            potential%alpha_core_charge = alpha_core_charge
       IF (PRESENT(alpha_ppl)) potential%alpha_ppl = alpha_ppl
       IF (PRESENT(ccore_charge)) potential%ccore_charge = ccore_charge
       IF (PRESENT(cerf_ppl)) potential%cerf_ppl = cerf_ppl
       IF (PRESENT(core_charge_radius))&
            potential%core_charge_radius = core_charge_radius
       IF (PRESENT(ppl_radius)) potential%ppl_radius = ppl_radius
       IF (PRESENT(ppnl_radius)) potential%ppnl_radius = ppnl_radius
       IF (PRESENT(lppnl)) potential%lppnl = lppnl
       IF (PRESENT(lprj_ppnl_max)) potential%lprj_ppnl_max = lprj_ppnl_max
       IF (PRESENT(nexp_ppl)) potential%nexp_ppl = nexp_ppl
       IF (PRESENT(nppnl)) potential%nppnl = nppnl
       IF (PRESENT(nprj_ppnl_max)) potential%nprj_ppnl_max = nprj_ppnl_max
       IF (PRESENT(z)) potential%z = z
       IF (PRESENT(zeff)) potential%zeff = zeff
       IF (PRESENT(zeff_correction)) potential%zeff_correction = zeff_correction
       IF (PRESENT(alpha_ppnl)) potential%alpha_ppnl => alpha_ppnl
       IF (PRESENT(cexp_ppl)) potential%cexp_ppl => cexp_ppl
       IF (PRESENT(elec_conf)) potential%elec_conf => elec_conf
       IF (PRESENT(nprj_ppnl)) potential%nprj_ppnl => nprj_ppnl
       IF (PRESENT(cprj)) potential%cprj => cprj
       IF (PRESENT(cprj_ppnl)) potential%cprj_ppnl => cprj_ppnl
       IF (PRESENT(vprj_ppnl)) potential%vprj_ppnl => vprj_ppnl
       IF (PRESENT(hprj_ppnl)) potential%hprj_ppnl => hprj_ppnl

    ELSE

       CALL stop_program(routineN,moduleN,&
            __LINE__,&
            "The pointer potential is not associated")

    END IF

  END SUBROUTINE set_gth_potential

! *****************************************************************************
  SUBROUTINE write_all_potential(potential,output_unit,error)

    ! Write an atomic all-electron potential data set to the output unit.

    ! - Creation (09.02.2002,MK)

    TYPE(all_potential_type), POINTER        :: potential
    INTEGER, INTENT(in)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=20)                        :: string

    IF (output_unit > 0.AND.ASSOCIATED(potential)) THEN
       WRITE (UNIT=output_unit,FMT="(/,T6,A,T41,A40,/)")&
            "Potential information for",ADJUSTR(TRIM(potential%name))
       WRITE (UNIT=output_unit,FMT="(T8,A,T41,A40)")&
            "Description: ",TRIM(potential%description(1)),&
            "             ",TRIM(potential%description(2))
       WRITE (UNIT=output_unit,FMT="(/,T8,A,T69,F12.6)")&
            "Gaussian exponent of the core charge distribution: ",&
            potential%alpha_core_charge
       WRITE (UNIT=string,FMT="(5I4)") potential%elec_conf
       WRITE (UNIT=output_unit,FMT="(T8,A,T61,A20)")&
            "Electronic configuration (s p d ...):",&
            ADJUSTR(TRIM(string))
    END IF

  END SUBROUTINE write_all_potential

! *****************************************************************************
  SUBROUTINE write_elp_potential(potential,output_unit,error)

    ! Write an effective-point-charge potential data set to the output unit.

    ! - Creation (24.09.2002,GT)

    TYPE(elp_potential_type), POINTER        :: potential
    INTEGER, INTENT(in)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    IF (output_unit>0.AND.ASSOCIATED(potential)) THEN
       WRITE (UNIT=output_unit,FMT="(/,T6,A,T41,A40,/)")&
            "Potential information for",ADJUSTR(TRIM(potential%name))
       WRITE (UNIT=output_unit,FMT="(T8,A,T41,A40)")&
            "Description: ",TRIM(potential%description(1))
       WRITE (UNIT=output_unit,FMT="(/,T8,A,T69,F12.6)")&
            "Atomic Effective Charge:",potential%qeff,&
            "Gaussian exponent of the core charge distribution: ",&
            potential%alpha_core_charge
    END IF

  END SUBROUTINE write_elp_potential

! *****************************************************************************
  SUBROUTINE write_gth_potential(potential,output_unit,error)

    ! Write an atomic GTH potential data set to the output unit.
    ! - Creation (09.02.2002,MK)

    TYPE(gth_potential_type), POINTER        :: potential
    INTEGER, INTENT(in)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=20)                        :: string
    INTEGER                                  :: i, j, l
    REAL(KIND=dp)                            :: r

    IF (output_unit>0.AND.ASSOCIATED(potential)) THEN
       WRITE (UNIT=output_unit,FMT="(/,T6,A,T41,A40,/)")&
            "Potential information for",ADJUSTR(TRIM(potential%name))
       WRITE (UNIT=output_unit,FMT="(T8,A,T41,A40)")&
            "Description: ",ADJUSTR(TRIM(potential%description(1))),&
            "             ",ADJUSTR(TRIM(potential%description(2))),&
            "             ",ADJUSTR(TRIM(potential%description(3))),&
            "             ",ADJUSTR(TRIM(potential%description(4)))
       WRITE (UNIT=output_unit,FMT="(/,T8,A,T69,F12.6)")&
            "Gaussian exponent of the core charge distribution: ",&
            potential%alpha_core_charge
       WRITE (UNIT=string,FMT="(5I4)") potential%elec_conf
       WRITE (UNIT=output_unit,FMT="(T8,A,T61,A20)")&
            "Electronic configuration (s p d ...):",&
            ADJUSTR(TRIM(string))

       r = 1.0_dp/SQRT(2.0_dp*potential%alpha_ppl)

       WRITE (UNIT=output_unit,FMT="(/,T8,A,/,/,T27,A,/,T21,5F12.6)")&
            "Parameters of the local part of the GTH pseudopotential:",&
            "rloc        C1          C2          C3          C4",&
            r,(potential%cexp_ppl(i)*r**(2*(i-1)),i=1,potential%nexp_ppl)

       IF (potential%lppnl > -1) THEN
          WRITE (UNIT=output_unit,FMT="(/,T8,A,/,/,T20,A,/)")&
               "Parameters of the non-local part of the GTH pseudopotential:",&
               "l      r(l)      h(i,j,l)"
          DO l=0,potential%lppnl
             r = SQRT(0.5_dp/potential%alpha_ppnl(l))
             WRITE (UNIT=output_unit,FMT="(T19,I2,5F12.6)")&
                  l,r,(potential%hprj_ppnl(1,j,l),j=1,potential%nprj_ppnl(l))
             DO i=2,potential%nprj_ppnl(l)
                WRITE (UNIT=output_unit,FMT="(T33,4F12.6)")&
                     (potential%hprj_ppnl(i,j,l),j=1,potential%nprj_ppnl(l))
             END DO
          END DO
       END IF
    END IF

  END SUBROUTINE write_gth_potential

END MODULE external_potential_types
