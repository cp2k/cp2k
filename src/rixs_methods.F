!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods for Resonant Inelastic XRAY Scattering (RIXS) calculations
!> \author BSG (02.2025)
! **************************************************************************************************
MODULE rixs_methods
    USE bibliography,                   ONLY: cite_reference, &
                                              VazdaCruz2021
    USE input_section_types,            ONLY: section_release,&
                                              section_type,&
                                              section_vals_create,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
    USE qs_environment_types,           ONLY: get_qs_env,&
                                              qs_environment_type
    USE cp_control_types,               ONLY: dft_control_type, rixs_control_type, &
                                              rixs_control_create, &
                                              rixs_control_release
    USE cp_control_utils,               ONLY: read_rixs_control
    USE cp_log_handling,                ONLY: cp_logger_get_default_io_unit
    ! USE rixs_types,                     ONLY: read_rixs_control !,&
    USE qs_tddfpt2_methods,             ONLY: tddfpt
    USE xas_tdp_methods,                ONLY: xas_tdp
#include "./base/base_uses.f90"
    IMPLICIT NONE
    PRIVATE

    CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rixs_methods'

    PUBLIC :: rixs, rixs_core

CONTAINS

! **************************************************************************************************
!> \brief Driver for RIXS calculations.
!> \param qs_env the inherited qs_environment
!> \author BSG
!> \note Empty for now...
! **************************************************************************************************

    SUBROUTINE rixs(qs_env)

        TYPE(qs_environment_type), POINTER                 :: qs_env
        TYPE(dft_control_type), POINTER                    :: dft_control

        CHARACTER(len=*), PARAMETER                        :: routineN = 'rixs'

        TYPE(section_vals_type), POINTER                   :: rixs_section, tddfp2_section, xas_tdp_section
        INTEGER                                            :: handle, output_unit

        CALL timeset(routineN, handle)

        NULLIFY(rixs_section)
        rixs_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%RIXS")
        output_unit = cp_logger_get_default_io_unit()

        qs_env%do_rixs = .TRUE.

        IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A,/,T3,A,/,T3,A,/,T3,A,/)") &
            "!===========================================================================!", &
            "!                                RIXS                                       !", &
            "!             Resonant Inelastic X-RAY Scattering calculations              !", &
            "!===========================================================================!"
        END IF

        CALL cite_reference(VazdaCruz2021)

        CALL get_qs_env(qs_env, dft_control=dft_control)

        ! TODO: is it possible to restart?

        ! check if tddfp2 and xas_tdp are "explicit"ily defined in this section,
        xas_tdp_section => section_vals_get_subs_vals(rixs_section, "XAS_TDP")
        tddfp2_section => section_vals_get_subs_vals(rixs_section, "TDDFPT")

        IF (.NOT. ASSOCIATED(xas_tdp_section)) THEN
            CPABORT("XAS_TDP calculation missing")
        END IF
        IF (.NOT. ASSOCIATED(tddfp2_section)) THEN
            CPABORT("TDDFPT calculation missing")
        END IF

        CALL rixs_core(rixs_section, qs_env)

        ! stuff

        ! IF (output_unit > 0) THEN
        !     WRITE (UNIT=output_unit, FMT="(/,T3,A,/,T3,A,/,T3,A,/,T3,A,/)") &
        !     "!===========================================================================!", &
        !     "!           End of Resonant Inelastic X-RAY Scattering calculations         !", &
        !     "!===========================================================================!"
        ! END IF

        CALL timestop(handle)

    END SUBROUTINE rixs

! **************************************************************************************************
!> \brief Perform RIXS calculation.
! **************************************************************************************************
    SUBROUTINE rixs_core(rixs_section, qs_env)

        TYPE(section_vals_type), POINTER                   :: rixs_section
        TYPE(qs_environment_type), POINTER                 :: qs_env

        CHARACTER(len=*), PARAMETER                        :: routineN = 'rixs_core'

        INTEGER                                            :: handle
        TYPE(dft_control_type), POINTER                    :: dft_control
        TYPE(rixs_control_type), POINTER                   :: rixs_control

        CALL timeset(routineN, handle)

        NULLIFY(rixs_control, dft_control)

        CALL get_qs_env(qs_env, dft_control=dft_control)
        CALL rixs_control_create(rixs_control)
        ! dft_control needed for tddft_control
        CALL read_rixs_control(rixs_control, rixs_section, dft_control%qs_control) 

        ! create rixs_env

        ! first, xas_tdp calculation
        CALL xas_tdp(qs_env)

        ! then tddpt2 calculation
        CALL tddfpt(qs_env, calc_forces=.FALSE., rixs_control=rixs_control) ! hardcode no forces for now

        ! the necessary variables are stored in rixs_env%xas_tdp_env%bla and rixs_env%valence_state%bla
        ! calculate rixs here

        ! should not nullify rixs_env because it will contain structure from xas_tdp and tddfpt2
        ! rixs_env%core should be allocated in xas_tdp
        ! rixs_env%valence should be allocated in tddpt2
        ! maybe for this one has to call rixs_env_create in xas, or maybe one can copy the structure to rixs_env
        ! both should be deallocated at the end of this section

        ! print stuff

        ! nullify rixs_control, rixs_env
        CALL rixs_control_release(rixs_control)

        CALL timestop(handle)            

    END SUBROUTINE rixs_core

END MODULE rixs_methods