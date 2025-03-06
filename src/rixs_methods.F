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
   USE input_section_types,             ONLY: section_release,&
                                              section_type,&
                                              section_vals_create,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
    USE qs_environment_types,           ONLY: get_qs_env,&
                                              qs_environment_type
    USE cp_control_types,               ONLY: dft_control_type
    USE cp_log_handling,                ONLY: cp_logger_get_default_io_unit
    USE rixs_types,                     ONLY: rixs_control_type,&
                                              rixs_control_create,&
                                              read_rixs_control
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

        TYPE(section_vals_type), POINTER                   :: rixs_section
        INTEGER                                            :: handle, output_unit

        CALL timeset(routineN, handle)

        NULLIFY(rixs_section)
        rixs_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%RIXS")
        output_unit = cp_logger_get_default_io_unit()

        IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A,/,T3,A,/,T3,A,/,T3,A,/)") &
            "!===========================================================================!", &
            "!                                RIXS                                       !", &
            "!             Resonant Inelastic X-RAY Scattering calculations              !", &
            "!===========================================================================!"
        END IF

        CALL get_qs_env(qs_env, dft_control=dft_control)

        ! TODO: is it possible to restart?

        IF (.NOT. dft_control%do_xas_tdp_calculation) THEN
            CPABORT("XAS_TDP calculation missing")
        END IF
        IF (.NOT. dft_control%tddfpt2_control%enabled) THEN
            CPABORT("TDDFPT calculation missing")
        END IF

        CALL rixs_core(rixs_section, qs_env)

        ! stuff?

        ! print stuff

        IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A,/,T3,A,/,T3,A,/,T3,A,/)") &
            "!===========================================================================!", &
            "!           End of Resonant Inelastic X-RAY Scattering calculations         !", &
            "!===========================================================================!"
        END IF

        CALL timestop(handle)

    END SUBROUTINE rixs

! **************************************************************************************************
!> \brief Perform RIXS calculation.
! **************************************************************************************************
    SUBROUTINE rixs_core(rixs_section, qs_env)

        TYPE(section_vals_type), POINTER                   :: rixs_section
        TYPE(qs_environment_type), POINTER                 :: qs_env

        CHARACTER(len=*), PARAMETER                        :: routineN = 'rixs_init'

        INTEGER                                            :: handle
        TYPE(dft_control_type), POINTER                    :: dft_control
        TYPE(rixs_control_type), POINTER                   :: rixs_control

        CALL timeset(routineN, handle)

        NULLIFY(rixs_section, rixs_control, dft_control)
        rixs_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%RIXS")

        CALL get_qs_env(qs_env, dft_control=dft_control)
        CALL rixs_control_create(rixs_control)
        CALL read_rixs_control(rixs_control, rixs_section)

        ! should not nullify rixs_env because it will contain structure from xas_tdp and tddpt2
        ! rixs_env%core should be allocated in xas_tdp
        ! rixs_env%valence should be allocated in tddpt2
        ! maybe for this one has to call rixs_env_create in xas, or maybe one can copy the structure to rixs_env
        ! both should be deallocated at the end of this section


        CALL timestop(handle)            

    END SUBROUTINE rixs_core

END MODULE rixs_methods