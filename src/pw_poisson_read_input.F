!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Reading of input parameters for the pw_poisson-modules.
!> \par History
!>      01.2014 Code moved into separate module to make pw_poisson-modules
!>              independet from input_section_types and input_constants.
!> \author Ole Schuett
! *****************************************************************************
MODULE pw_poisson_read_input
  USE f77_blas
  USE input_constants,                 ONLY: use_perd_none,&
                                             use_perd_x,&
                                             use_perd_xy,&
                                             use_perd_xyz,&
                                             use_perd_xz,&
                                             use_perd_y,&
                                             use_perd_yz,&
                                             use_perd_z
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE ps_wavelet_types,                ONLY: WAVELET0D,&
                                             WAVELET2D,&
                                             WAVELET3D
  USE pw_poisson_types,                ONLY: do_ewald_none,&
                                             pw_poisson_analytic,&
                                             pw_poisson_mt,&
                                             pw_poisson_multipole,&
                                             pw_poisson_none,&
                                             pw_poisson_parameter_type,&
                                             pw_poisson_periodic,&
                                             pw_poisson_wavelet
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_poisson_read_input'

  PUBLIC :: pw_poisson_read_parameters

CONTAINS

! *****************************************************************************
!> \brief Reads the POISSON input-section and into pw_poisson_parameter_type.
!> \param poisson_section ...
!> \param params ...
!> \param error ...
!> \par History
!>      01.2014 Code moved into separate module from pw_poisson_types,
!>              pw_poisson_methods and ps_wavelet_types.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE pw_poisson_read_parameters(poisson_section, params, error)
    TYPE(section_vals_type), POINTER         :: poisson_section
    TYPE(pw_poisson_parameter_type), &
      INTENT(INOUT)                          :: params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_poisson_read_parameters', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: periodic
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: ewald_section, mt_section, &
                                                wavelet_section

    NULLIFY(ewald_section, mt_section, wavelet_section)
    failure=.FALSE.

    CALL section_vals_val_get(poisson_section,"POISSON_SOLVER",i_val=params%solver,error=error)

    ! Decoding PERIODIC depending on chosen solver,
    ! because not all solvers support every possible periodicity
    CALL section_vals_val_get(poisson_section,"PERIODIC",i_val=periodic,error=error)
    SELECT CASE (params%solver)
      CASE(pw_poisson_periodic,pw_poisson_analytic,pw_poisson_mt,pw_poisson_multipole)
         CALL decode_periodic_green(periodic, params, error)
      CASE(pw_poisson_wavelet)
         CALL decode_periodic_wavelet(periodic, params, error)
      CASE(pw_poisson_none)
      CASE default
         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    ! parsing EWALD subsection
    params%ewald_type = do_ewald_none
    ewald_section => section_vals_get_subs_vals(poisson_section,"EWALD", can_return_null=.TRUE., error=error)
    IF ( ASSOCIATED(ewald_section) ) THEN
         CALL section_vals_val_get(ewald_section,"EWALD_TYPE",i_val=params%ewald_type,error=error)
         CALL section_vals_val_get(ewald_section,"o_spline", i_val=params%ewald_o_spline,error=error)
         CALL section_vals_val_get(ewald_section,"alpha",r_val=params%ewald_alpha,error=error)
    ENDIF

    ! parsing MT subsection
    mt_section => section_vals_get_subs_vals(poisson_section,"MT",error=error)
    CALL section_vals_val_get(mt_section,"REL_CUTOFF",r_val=params%mt_rel_cutoff,error=error)
    CALL section_vals_val_get(mt_section,"ALPHA",r_val=params%mt_alpha,error=error)

    ! parsing WAVELET subsection
    wavelet_section => section_vals_get_subs_vals(poisson_section,"WAVELET",error=error)
    CALL section_vals_val_get(wavelet_section,"SCF_TYPE",i_val=params%wavelet_scf_type,error=error)
  END SUBROUTINE pw_poisson_read_parameters


! *****************************************************************************
!> \brief Helper routien for pw_poisson_read_parameters
!> \param periodic ...
!> \param params ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE decode_periodic_green(periodic, params, error)
    INTEGER, INTENT(IN)                      :: periodic
    TYPE(pw_poisson_parameter_type), &
      INTENT(INOUT)                          :: params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'decode_periodic_green', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    SELECT CASE(periodic)
      CASE(use_perd_x)
         params%periodic = (/1,0,0/)
      CASE(use_perd_y)
         params%periodic = (/0,1,0/)
      CASE(use_perd_z)
         params%periodic = (/0,0,1/)
      CASE(use_perd_xy)
         params%periodic = (/1,1,0/)
      CASE(use_perd_xz)
         params%periodic = (/1,0,1/)
      CASE(use_perd_yz)
         params%periodic = (/0,1,1/)
      CASE(use_perd_xyz)
         params%periodic = (/1,1,1/)
      CASE(use_perd_none)
         params%periodic = (/0,0,0/)
      CASE DEFAULT
         CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
    ! check for consistent use of periodicity (cell <-> Poisson solver)
    !CPPostcondition(ALL(perd == cell%perd),cp_fatal_level,routineP,error,failure)

  END SUBROUTINE decode_periodic_green


! *****************************************************************************
!> \brief Helper routien for pw_poisson_read_parameters
!> \param periodic ...
!> \param params ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE decode_periodic_wavelet(periodic, params, error)
    INTEGER, INTENT(IN)                      :: periodic
    TYPE(pw_poisson_parameter_type), &
      INTENT(INOUT)                          :: params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'decode_periodic_wavelet', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    params%wavelet_special_dimension=0

    SELECT CASE (periodic)
    CASE (use_perd_none)
       params%wavelet_geocode="F"
       params%wavelet_method=WAVELET0D
    CASE (use_perd_xz)
       params%wavelet_geocode="S"
       params%wavelet_method=WAVELET2D
       params%wavelet_special_dimension=2
    CASE (use_perd_xyz)
       params%wavelet_geocode="P"
       params%wavelet_method=WAVELET3D
    CASE(use_perd_x,use_perd_y,use_perd_z,use_perd_xy,use_perd_yz)
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "Poisson solver for this periodicity not yet implemented",&
            error=error,failure=failure)
    CASE DEFAULT
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE decode_periodic_wavelet

END MODULE pw_poisson_read_input
