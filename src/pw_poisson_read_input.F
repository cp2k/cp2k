!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Reading of input parameters for the pw_poisson-modules.
!> \par History
!>      01.2014 Code moved into separate module to make pw_poisson-modules
!>              independet from input_section_types and input_constants.
!> \author Ole Schuett
! *****************************************************************************
MODULE pw_poisson_read_input
  USE cell_types,                      ONLY: use_perd_none,&
                                             use_perd_x,&
                                             use_perd_xy,&
                                             use_perd_xyz,&
                                             use_perd_xz,&
                                             use_perd_y,&
                                             use_perd_yz,&
                                             use_perd_z
  USE dirichlet_bc_types,              ONLY: xy_aligned_rectangle,&
                                             xz_aligned_rectangle,&
                                             yz_aligned_rectangle
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE ps_wavelet_types,                ONLY: WAVELET0D,&
                                             WAVELET2D,&
                                             WAVELET3D
  USE pw_poisson_types,                ONLY: &
       do_ewald_none, pw_poisson_analytic, pw_poisson_implicit, &
       pw_poisson_mt, pw_poisson_multipole, pw_poisson_none, &
       pw_poisson_parameter_type, pw_poisson_periodic, pw_poisson_wavelet
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
      CASE(pw_poisson_implicit)
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

    ! parsing IMPLICIT subsection
    CALL ps_implicit_read_parameters(poisson_section, params, error)

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

! *****************************************************************************
!> \brief Reads the subsection IMPLICIT and initializes corresponding parameters in
!>        pw_poisson_parameter_type
!> \param poisson_section poisson section to be read from input
!> \param params poisson_env parameters 
!> \param error  cp2k error 
!> \par History
!>      08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE ps_implicit_read_parameters(poisson_section, params, error)
    TYPE(section_vals_type), POINTER         :: poisson_section
    TYPE(pw_poisson_parameter_type), &
      INTENT(INOUT)                          :: params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ps_implicit_read_parameters', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, n_aac_rep, n_aap_rep, n_p_rep, n_xaac_rep, parallel_plane, &
      xaa_cylindrical_apxtyp, xaa_cylindrical_nsides
    INTEGER, DIMENSION(:), POINTER           :: aa_cuboidal_nprtn, &
                                                aa_planar_nprtn, &
                                                planar_nprtn, &
                                                xaa_cylindrical_nprtn
    LOGICAL :: aa_cuboidal_explicit, aa_planar_explicit, has_dielectric, &
      planar_explicit, smooth, xaa_cylindrical_explicit
    REAL(dp)                                 :: intercept, v_D, &
                                                xaa_cylindrical_brad, zeta
    REAL(dp), DIMENSION(:), POINTER :: aa_cuboidal_xxtnt, aa_cuboidal_yxtnt, &
      aa_cuboidal_zxtnt, aa_planar_xxtnt, aa_planar_yxtnt, aa_planar_zxtnt, &
      planar_Avtx, planar_Bvtx, planar_Cvtx, xaa_cylindrical_bctr, &
      xaa_cylindrical_xxtnt
    TYPE(section_vals_type), POINTER :: aa_cuboidal_section, &
      aa_planar_section, dbc_section, dielectric_section, planar_section, &
      ps_implicit_section, xaa_cylindrical_section

    NULLIFY(ps_implicit_section, dielectric_section, dbc_section)

    ! parsing IMPLICIT subsection
    ps_implicit_section => section_vals_get_subs_vals(poisson_section, "IMPLICIT", error=error)
    CALL section_vals_val_get(ps_implicit_section, "BOUNDARY_CONDITIONS", &
                              i_val=params%ps_implicit_params%boundary_condition, error=error)
    CALL section_vals_val_get(ps_implicit_section, "ZERO_INITIAL_GUESS", &
                              l_val=params%ps_implicit_params%zero_initial_guess, error=error)
    CALL section_vals_val_get(ps_implicit_section, "max_iter", i_val=params%ps_implicit_params%max_iter, error=error)
    CALL section_vals_val_get(ps_implicit_section, "tol", r_val=params%ps_implicit_params%tol, error=error)
    CALL section_vals_val_get(ps_implicit_section, "omega", r_val=params%ps_implicit_params%omega, error=error)

    ! parsing DIELECTRIC subsection
    dielectric_section => section_vals_get_subs_vals(ps_implicit_section, "DIELECTRIC", error=error)
    CALL section_vals_get(dielectric_section, explicit=has_dielectric, error=error)
    CALL section_vals_val_get(dielectric_section, "epsilon", r_val=params%dielectric_params%eps0, error=error)
    CALL section_vals_val_get(dielectric_section, "rho_min", r_val=params%dielectric_params%rho_min, error=error)
    CALL section_vals_val_get(dielectric_section, "rho_max", r_val=params%dielectric_params%rho_max, error=error)
    CALL section_vals_val_get(dielectric_section, "DERIVATIVE_METHOD", &
                              i_val=params%dielectric_params%derivative_method, error=error)
    params%has_dielectric = has_dielectric

    ! parsing DIRICHLET_BC subsection
    dbc_section => section_vals_get_subs_vals(ps_implicit_section, "DIRICHLET_BC", error=error)
    CALL section_vals_val_get(dbc_section, "DEBUG", l_val=params%dbc_params%debug, error=error)
    aa_planar_section => section_vals_get_subs_vals(dbc_section, "AA_PLANAR", error=error)
    planar_section => section_vals_get_subs_vals(dbc_section, "PLANAR", error=error)
    xaa_cylindrical_section => section_vals_get_subs_vals(dbc_section, "XAA_CYLINDRICAL", error=error)
    aa_cuboidal_section => section_vals_get_subs_vals(dbc_section, "AA_CUBOIDAL", error=error)
    CALL section_vals_get(aa_planar_section, explicit=aa_planar_explicit, n_repetition=n_aap_rep, error=error)
    CALL section_vals_get(planar_section, explicit=planar_explicit, n_repetition=n_p_rep, error=error)
    CALL section_vals_get(xaa_cylindrical_section, explicit=xaa_cylindrical_explicit, n_repetition=n_xaac_rep, error=error)
    CALL section_vals_get(aa_cuboidal_section, explicit=aa_cuboidal_explicit, n_repetition=n_aac_rep, error=error)
    IF (params%solver .EQ. pw_poisson_implicit) THEN
       IF (aa_planar_explicit) THEN
          params%dbc_params%n_aa_planar = n_aap_rep
          ALLOCATE(params%dbc_params%aa_planar_xxtnt(2, n_aap_rep), &
                   params%dbc_params%aa_planar_yxtnt(2, n_aap_rep), &
                   params%dbc_params%aa_planar_zxtnt(2, n_aap_rep), &
                   params%dbc_params%aa_planar_nprtn(3, n_aap_rep), &
                   params%dbc_params%aa_planar_vD(n_aap_rep),&
                   params%dbc_params%aa_planar_zeta(n_aap_rep),&
                   params%dbc_params%aa_planar_smooth(n_aap_rep))
          NULLIFY(aa_planar_xxtnt, aa_planar_yxtnt, aa_planar_zxtnt, aa_planar_nprtn)
          DO i = 1, n_aap_rep
             CALL section_vals_val_get(aa_planar_section,"v_D"   ,&
                                       i_rep_section=i,r_val=v_D,error=error)
             CALL section_vals_val_get(aa_planar_section,"zeta"  ,&
                                       i_rep_section=i,r_val=zeta,error=error)
             CALL section_vals_val_get(aa_planar_section,"smooth",&
                                       i_rep_section=i,l_val=smooth,error=error)

             params%dbc_params%aa_planar_vD(i) = v_D
             params%dbc_params%aa_planar_zeta(i) = zeta
             params%dbc_params%aa_planar_smooth(i) = smooth

             CALL section_vals_val_get(aa_planar_section,"PARALLEL_PLANE",&
                                       i_rep_section=i,i_val=parallel_plane,error=error)
             CALL section_vals_val_get(aa_planar_section,"INTERCEPT",&
                                       i_rep_section=i,r_val=intercept,error=error)
             SELECT CASE (parallel_plane)
             CASE (xy_aligned_rectangle)
                CALL section_vals_val_get(aa_planar_section,"X_xtnt",&
                                       i_rep_section=i,r_vals=aa_planar_xxtnt,error=error)
                CALL section_vals_val_get(aa_planar_section,"Y_xtnt",&
                                       i_rep_section=i,r_vals=aa_planar_yxtnt,error=error)

                params%dbc_params%aa_planar_xxtnt(:,i) = aa_planar_xxtnt
                params%dbc_params%aa_planar_yxtnt(:,i) = aa_planar_yxtnt
                params%dbc_params%aa_planar_zxtnt(:,i) = intercept

                CALL section_vals_val_get(aa_planar_section,"n_prtn",&
                                       i_rep_section=i,i_vals=aa_planar_nprtn,error=error)

                params%dbc_params%aa_planar_nprtn(1,i) = aa_planar_nprtn(1)
                params%dbc_params%aa_planar_nprtn(2,i) = aa_planar_nprtn(2)
                params%dbc_params%aa_planar_nprtn(3,i) = 1
             CASE (yz_aligned_rectangle)
                CALL section_vals_val_get(aa_planar_section,"Y_xtnt",&
                                       i_rep_section=i,r_vals=aa_planar_yxtnt,error=error)
                CALL section_vals_val_get(aa_planar_section,"Z_xtnt",&
                                       i_rep_section=i,r_vals=aa_planar_zxtnt,error=error)

                params%dbc_params%aa_planar_xxtnt(:,i) = intercept
                params%dbc_params%aa_planar_yxtnt(:,i) = aa_planar_yxtnt
                params%dbc_params%aa_planar_zxtnt(:,i) = aa_planar_zxtnt

                CALL section_vals_val_get(aa_planar_section,"n_prtn",&
                                       i_rep_section=i,i_vals=aa_planar_nprtn,error=error)

                params%dbc_params%aa_planar_nprtn(1,i) = 1
                params%dbc_params%aa_planar_nprtn(2,i) = aa_planar_nprtn(1)
                params%dbc_params%aa_planar_nprtn(3,i) = aa_planar_nprtn(2)
             CASE (xz_aligned_rectangle)
                CALL section_vals_val_get(aa_planar_section,"X_xtnt",&
                                       i_rep_section=i,r_vals=aa_planar_xxtnt,error=error)
                CALL section_vals_val_get(aa_planar_section,"Z_xtnt",&
                                       i_rep_section=i,r_vals=aa_planar_zxtnt,error=error)

                params%dbc_params%aa_planar_xxtnt(:,i) = aa_planar_xxtnt
                params%dbc_params%aa_planar_yxtnt(:,i) = intercept
                params%dbc_params%aa_planar_zxtnt(:,i) = aa_planar_zxtnt

                CALL section_vals_val_get(aa_planar_section,"n_prtn",&
                                       i_rep_section=i,i_vals=aa_planar_nprtn,error=error)

                params%dbc_params%aa_planar_nprtn(1,i) = aa_planar_nprtn(1)
                params%dbc_params%aa_planar_nprtn(2,i) = 1
                params%dbc_params%aa_planar_nprtn(3,i) = aa_planar_nprtn(2)
             END SELECT

          END DO
       ELSE
         params%dbc_params%n_aa_planar = 0
       END IF
       IF (planar_explicit) THEN
          params%dbc_params%n_planar = n_p_rep
          ALLOCATE(params%dbc_params%planar_Avtx(3, n_p_rep), &
                   params%dbc_params%planar_Bvtx(3, n_p_rep), &
                   params%dbc_params%planar_Cvtx(3, n_p_rep), &
                   params%dbc_params%planar_nprtn(3, n_p_rep), &
                   params%dbc_params%planar_vD(n_p_rep), &
                   params%dbc_params%planar_zeta(n_p_rep), &
                   params%dbc_params%planar_smooth(n_p_rep))
          NULLIFY(planar_Avtx, planar_Bvtx, planar_Cvtx, planar_nprtn)
          DO i = 1, n_p_rep
             CALL section_vals_val_get(planar_section,"v_D"   ,&
                                       i_rep_section=i,r_val=v_D,error=error)
             CALL section_vals_val_get(planar_section,"zeta"  ,&
                                       i_rep_section=i,r_val=zeta,error=error)
             CALL section_vals_val_get(planar_section,"smooth",&
                                       i_rep_section=i,l_val=smooth,error=error)
             CALL section_vals_val_get(planar_section,"A"     ,&
                                       i_rep_section=i,r_vals=planar_Avtx,error=error)
             CALL section_vals_val_get(planar_section,"B"     ,&
                                       i_rep_section=i,r_vals=planar_Bvtx,error=error)
             CALL section_vals_val_get(planar_section,"C"     ,&
                                       i_rep_section=i,r_vals=planar_Cvtx,error=error)
             CALL section_vals_val_get(planar_section,"n_prtn",&
                                       i_rep_section=i,i_vals=planar_nprtn,error=error)
             params%dbc_params%planar_vD(i) = v_D
             params%dbc_params%planar_zeta(i) = zeta
             params%dbc_params%planar_smooth(i) = smooth
             params%dbc_params%planar_Avtx(:,i) = planar_Avtx
             params%dbc_params%planar_Bvtx(:,i) = planar_Bvtx
             params%dbc_params%planar_Cvtx(:,i) = planar_Cvtx
             params%dbc_params%planar_nprtn(:,i) = planar_nprtn
          END DO
       ELSE
          params%dbc_params%n_planar = 0
       END IF
       IF (xaa_cylindrical_explicit) THEN
          params%dbc_params%n_xaa_cylindrical = n_xaac_rep
          ALLOCATE(params%dbc_params%xaa_cylindrical_xxtnt(2, n_xaac_rep), &
                   params%dbc_params%xaa_cylindrical_nprtn(3, n_xaac_rep), &
                   params%dbc_params%xaa_cylindrical_bctr(2, n_xaac_rep), &
                   params%dbc_params%xaa_cylindrical_brad(n_xaac_rep), &
                   params%dbc_params%xaa_cylindrical_nsides(n_xaac_rep), &
                   params%dbc_params%xaa_cylindrical_apxtyp(n_xaac_rep), &
                   params%dbc_params%xaa_cylindrical_vD(n_xaac_rep), &
                   params%dbc_params%xaa_cylindrical_zeta(n_xaac_rep), &
                   params%dbc_params%xaa_cylindrical_smooth(n_xaac_rep))
          NULLIFY(xaa_cylindrical_xxtnt, xaa_cylindrical_bctr, xaa_cylindrical_nprtn)
          DO i = 1, n_xaac_rep
             CALL section_vals_val_get(xaa_cylindrical_section,"v_D"        ,&
                                       i_rep_section=i,r_val=v_D,error=error)
             CALL section_vals_val_get(xaa_cylindrical_section,"zeta"       ,&
                                       i_rep_section=i,r_val=zeta,error=error)
             CALL section_vals_val_get(xaa_cylindrical_section,"smooth"     ,&
                                       i_rep_section=i,l_val=smooth,error=error)
             CALL section_vals_val_get(xaa_cylindrical_section,"X_xtnt"  ,&
                                       i_rep_section=i,r_vals=xaa_cylindrical_xxtnt,error=error)
             CALL section_vals_val_get(xaa_cylindrical_section,"base_center",&
                                       i_rep_section=i,r_vals=xaa_cylindrical_bctr,error=error)
             CALL section_vals_val_get(xaa_cylindrical_section,"base_radius",&
                                       i_rep_section=i,r_val=xaa_cylindrical_brad,error=error)
             CALL section_vals_val_get(xaa_cylindrical_section,"n_sides"    ,&
                                       i_rep_section=i,i_val=xaa_cylindrical_nsides,error=error)
             CALL section_vals_val_get(xaa_cylindrical_section,"apx_type"   ,&
                                       i_rep_section=i,i_val=xaa_cylindrical_apxtyp,error=error)
             CALL section_vals_val_get(xaa_cylindrical_section,"n_prtn"     ,&
                                       i_rep_section=i,i_vals=xaa_cylindrical_nprtn,error=error)
             params%dbc_params%xaa_cylindrical_vD(i) = v_D
             params%dbc_params%xaa_cylindrical_zeta(i) = zeta
             params%dbc_params%xaa_cylindrical_smooth(i) = smooth
             params%dbc_params%xaa_cylindrical_xxtnt(:,i) = xaa_cylindrical_xxtnt
             params%dbc_params%xaa_cylindrical_bctr(:,i)  = xaa_cylindrical_bctr
             params%dbc_params%xaa_cylindrical_brad(i)    = xaa_cylindrical_brad
             params%dbc_params%xaa_cylindrical_nsides(i)  = xaa_cylindrical_nsides
             params%dbc_params%xaa_cylindrical_apxtyp(i)  = xaa_cylindrical_apxtyp
             params%dbc_params%xaa_cylindrical_nprtn(:,i) = xaa_cylindrical_nprtn
          END DO
       ELSE
          params%dbc_params%n_xaa_cylindrical = 0
          ALLOCATE(params%dbc_params%xaa_cylindrical_nsides(n_xaac_rep))
       END IF
       IF (aa_cuboidal_explicit) THEN
          params%dbc_params%n_aa_cuboidal = n_aac_rep
          ALLOCATE(params%dbc_params%aa_cuboidal_xxtnt(2, n_aac_rep), &
                   params%dbc_params%aa_cuboidal_yxtnt(2, n_aac_rep), &
                   params%dbc_params%aa_cuboidal_zxtnt(2, n_aac_rep), &
                   params%dbc_params%aa_cuboidal_nprtn(3, n_aac_rep), &
                   params%dbc_params%aa_cuboidal_vD(n_aac_rep),&
                   params%dbc_params%aa_cuboidal_zeta(n_aac_rep),&
                   params%dbc_params%aa_cuboidal_smooth(n_aac_rep))
          NULLIFY(aa_cuboidal_xxtnt, aa_cuboidal_yxtnt, aa_cuboidal_zxtnt, aa_cuboidal_nprtn)
          DO i = 1, n_aac_rep
             CALL section_vals_val_get(aa_cuboidal_section,"v_D"   ,&
                                       i_rep_section=i,r_val=v_D,error=error)
             CALL section_vals_val_get(aa_cuboidal_section,"zeta"  ,&
                                       i_rep_section=i,r_val=zeta,error=error)
             CALL section_vals_val_get(aa_cuboidal_section,"smooth",&
                                       i_rep_section=i,l_val=smooth,error=error)
             CALL section_vals_val_get(aa_cuboidal_section,"X_xtnt",&
                                       i_rep_section=i,r_vals=aa_cuboidal_xxtnt,error=error)
             CALL section_vals_val_get(aa_cuboidal_section,"Y_xtnt",&
                                       i_rep_section=i,r_vals=aa_cuboidal_yxtnt,error=error)
             CALL section_vals_val_get(aa_cuboidal_section,"Z_xtnt",&
                                       i_rep_section=i,r_vals=aa_cuboidal_zxtnt,error=error)
             CALL section_vals_val_get(aa_cuboidal_section,"n_prtn",&
                                       i_rep_section=i,i_vals=aa_cuboidal_nprtn,error=error)
             params%dbc_params%aa_cuboidal_vD(i) = v_D
             params%dbc_params%aa_cuboidal_zeta(i) = zeta
             params%dbc_params%aa_cuboidal_smooth(i) = smooth
             params%dbc_params%aa_cuboidal_xxtnt(:,i) = aa_cuboidal_xxtnt
             params%dbc_params%aa_cuboidal_yxtnt(:,i) = aa_cuboidal_yxtnt
             params%dbc_params%aa_cuboidal_zxtnt(:,i) = aa_cuboidal_zxtnt
             params%dbc_params%aa_cuboidal_nprtn(:,i) = aa_cuboidal_nprtn
          END DO
       ELSE
         params%dbc_params%n_aa_cuboidal = 0
       END IF
    END IF

  END SUBROUTINE ps_implicit_read_parameters

END MODULE pw_poisson_read_input
