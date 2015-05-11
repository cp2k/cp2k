!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief function that build the poisson section of the input
!> \par History
!>      03.2006 fusing of poisson_dft and poisson_mm
!> \author fawzi
! *****************************************************************************
MODULE input_cp2k_poisson
  USE bibliography,                    ONLY: &
       Aguado2003, Blochl1995, Darden1993, Essmann1995, Ewald1921, &
       Genovese2006, Genovese2007, Laino2008, Martyna1999, Toukmaji1996
  USE cell_types,                      ONLY: use_perd_none,&
                                             use_perd_x,&
                                             use_perd_xy,&
                                             use_perd_xyz,&
                                             use_perd_xz,&
                                             use_perd_y,&
                                             use_perd_yz,&
                                             use_perd_z
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             low_print_level,&
                                             medium_print_level
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE dielectric_types,                ONLY: derivative_cd3,&
                                             derivative_cd5,&
                                             derivative_cd7,&
                                             derivative_fft,&
                                             derivative_fft_use_deps,&
                                             derivative_fft_use_drho
  USE dirichlet_bc_types,              ONLY: CIRCUMSCRIBED,&
                                             INSCRIBED,&
                                             xy_aligned_rectangle,&
                                             xz_aligned_rectangle,&
                                             yz_aligned_rectangle
  USE input_constants,                 ONLY: do_fist_pol_cg,&
                                             do_fist_pol_none,&
                                             do_fist_pol_sc
  USE input_cp2k_rsgrid,               ONLY: create_rsgrid_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: enum_t,&
                                             integer_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE multipole_types,                 ONLY: do_multipole_charge,&
                                             do_multipole_dipole,&
                                             do_multipole_none,&
                                             do_multipole_quadrupole
  USE ps_implicit_types,               ONLY: MIXED_BC,&
                                             MIXED_PERIODIC_BC,&
                                             NEUMANN_BC,&
                                             PERIODIC_BC
  USE pw_poisson_types,                ONLY: &
       do_ewald_ewald, do_ewald_none, do_ewald_pme, do_ewald_spme, &
       pw_poisson_analytic, pw_poisson_implicit, pw_poisson_mt, &
       pw_poisson_multipole, pw_poisson_periodic, pw_poisson_wavelet
  USE pw_spline_utils,                 ONLY: no_precond,&
                                             precond_spl3_1,&
                                             precond_spl3_2,&
                                             precond_spl3_3,&
                                             precond_spl3_aint,&
                                             precond_spl3_aint2
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_poisson'

  PUBLIC :: create_poisson_section,&
            create_gspace_interp_section,&
            create_ewald_section
!***
CONTAINS

! *****************************************************************************
!> \brief Creates the Poisson section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_poisson_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_poisson_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="poisson",&
            description="Sets up the poisson resolutor.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)
       CALL keyword_create(keyword, name="POISSON_SOLVER",&
            variants=(/"POISSON", "PSOLVER"/),&
            description="Specify which kind of solver to use to solve the Poisson equation.",&
            usage="POISSON_SOLVER char",&
            enum_c_vals=s2a( "PERIODIC", "ANALYTIC", "MT", "MULTIPOLE", "WAVELET", "IMPLICIT" ),&
            enum_i_vals=(/ pw_poisson_periodic, pw_poisson_analytic, pw_poisson_mt, pw_poisson_multipole, &
                           pw_poisson_wavelet, pw_poisson_implicit/),&
            citations=(/Blochl1995,Martyna1999,Genovese2006,Genovese2007/),&
            default_i_val=pw_poisson_periodic, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PERIODIC",&
            description="Specify the directions on wich apply PBC. Important notice, "//&
                        " this only applies to the electrostatics."//&
                        " See the CELL section to specify the periodicity used for e.g. the pair lists."//&
                        " Typically the settings should be the same.",&
            usage="PERIODIC (x|y|z|xy|xz|yz|xyz|none)",&
            enum_c_vals=s2a( "x","y","z","xy","xz","yz","xyz","none"),&
            enum_i_vals=(/ use_perd_x,  use_perd_y,   use_perd_z,&
                           use_perd_xy, use_perd_xz, use_perd_yz,&
                           use_perd_xyz, use_perd_none /),&
            default_i_val=use_perd_xyz, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_mt_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_wavelet_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_multipole_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_ewald_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_implicit_ps_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_poisson_section

! *****************************************************************************
!> \brief Section to set-up parameters for decoupling using the Bloechl scheme
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_multipole_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_multipole_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="MULTIPOLE",&
            description="This section is used to set up the decoupling of QM periodic images with "//&
            "the use of density derived atomic point charges.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)
       CALL keyword_create(keyword, name="RCUT",&
            description="Real space cutoff for the Ewald sum.",&
            usage="RCUT {real}", n_var=1, type_of_var=real_t,&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EWALD_PRECISION",&
            description="Precision achieved in the Ewald sum.",&
            usage="EWALD_PRECISION {real}", n_var=1, type_of_var=real_t,&
            unit_str="hartree",default_r_val=1.0E-6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ANALYTICAL_GTERM",&
            description="Evaluates the Gterm in the Ewald Scheme analytically instead of using Splines.",&
            usage="ANALYTICAL_GTERM <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NGRIDS",&
            description="Specifies the number of grid points used for the Interpolation of the G-space term",&
            usage="NGRIDS <integer> <iteger> <integer> ",n_var=3,default_i_vals=(/50,50,50/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_gspace_interp_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL cp_print_key_section_create(subsection,"check_spline",&
            description="Controls the checking of the G-space term Spline Interpolation.",&
            print_level=medium_print_level,filename="GSpace-SplInterp",&
            error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL cp_print_key_section_create(subsection,"program_run_info",&
            description="Controls the printing of basic information during the run", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_multipole_section

! *****************************************************************************
!> \brief Creates the Martyna-Tuckerman section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_mt_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mt_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="mt",&
            description="Sets up parameters of  Martyna-Tuckerman poisson solver. "//&
                        "Note that exact results are only guaranteed if the unit cell is "//&
                        "twice as large as charge density (and serious artefacts can result "//&
                        "if the cell is much smaller).",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            citations=(/Martyna1999/),&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ALPHA",&
            description="Convergence parameter ALPHA*RMIN. Default value 7.0",&
            usage="ALPHA real",&
            n_var=1,default_r_val=7.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="REL_CUTOFF",&
            description="Specify the multiplicative factor for the CUTOFF keyword in MULTI_GRID "//&
            " section. The result gives the cutoff at which the 1/r non-periodic FFT3D is evaluated."//&
            "Default is 2.0",&
            usage="REL_CUTOFF real",&
            n_var=1,default_r_val=2.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_mt_section

! *****************************************************************************
!> \brief ...
!> \param section will contain the ewald section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE create_ewald_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ewald_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

  failure=.FALSE.

  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
    CALL section_create(section,name="ewald",&
         description="Ewald parameters controlling electrostatic only for CLASSICAL MM.",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         citations=(/Ewald1921,Darden1993,Essmann1995,Toukmaji1996,Laino2008/),&
         error=error)

    NULLIFY(keyword,print_key,subsection)
    CALL keyword_create(keyword, name="EWALD_TYPE",&
         description="The type of ewald you want to perform.",&
         citations=(/Ewald1921,Essmann1995,Darden1993/),&
         usage="EWALD_TYPE (NONE|EWALD|PME|SPME)",&
         default_i_val=do_ewald_ewald,&
         enum_c_vals=(/"none     ",&
                       "ewald    ",&
                       "pme      ",&
                       "spme     " /),&
         enum_i_vals=(/do_ewald_none,&
                       do_ewald_ewald,&
                       do_ewald_pme,&
                       do_ewald_spme/),&
         enum_desc=s2a("NONE standard real-space coulomb potential is computed together with the non-bonded contributions",&
                       "EWALD is the standard non-fft based ewald",&
                       "PME is the particle mesh using fft interpolation",&
                       "SPME is the smooth particle mesh using beta-Euler splines (recommended)"),&
                       error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EWALD_ACCURACY",&
         description="Expected accuracy in the Ewald sum. This number affects only the calculation of "//&
         "the cutoff for the real-space term of the ewald summation (EWALD|PME|SPME) as well as the "//&
         "construction of the neighbor lists (if the cutoff for non-bonded terms is smaller than the "//&
         "value employed to compute the EWALD real-space term). This keyword has no "//&
         "effect on the reciprocal space term (which can be tuned independently).",&
         usage="EWALD_ACCURACY {real}", n_var=1, type_of_var=real_t,&
         unit_str="hartree",default_r_val=1.0E-6_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RCUT",&
         description="Explicitly provide the real-space cutoff of the ewald summation (EWALD|PME|SPME). "//&
                     "If present, overwrites the estimate of EWALD_ACCURACY and may affect the "//&
                     "construction of the neighbor lists for non-bonded terms (in FIST), if the value "//&
                     "specified is larger than the cutoff for non-bonded interactions.",&
         usage="RCUT 5.0", n_var=1, type_of_var=real_t, unit_str="angstrom",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="alpha",&
         description="alpha parameter associated with Ewald (EWALD|PME|SPME). "//&
                     "Recommended for small systems is is alpha = 3.5 / r_cut. "//&
                     "Tuning alpha, r_cut and gmax is needed to obtain O(N**1.5) scaling for ewald.",&
         usage="alpha .30",&
         default_r_val=cp_unit_to_cp2k(value=0.35_dp,unit_str="angstrom^-1", error=error), &
         unit_str='angstrom^-1',error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="gmax",&
         description="number of grid points (SPME and EWALD). If a single number is specified,"// &
                     "the same number of points is used for all three directions on the grid."// &
                     "If three numbers are given, each direction can have a different number of points."// &
                     "The number of points needs to be FFTable (which depends on the library used) and odd for EWALD."// &
                     "The optimal number depends e.g. on alpha and the size of the cell. 1 point per Angstrom is common.",&
         usage="gmax 25 25 25", n_var=-1, type_of_var=integer_t, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ns_max",&
         description="number of grid points on small mesh (PME only), should be odd.",&
         usage="ns_max 11", default_i_val=11,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="o_spline",&
         description="order of the beta-Euler spline (SPME only)",&
         usage="o_spline 6", default_i_val=6,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="epsilon",&
         description="tolerance of gaussians for fft interpolation (PME only)",&
         usage="epsilon 1e-6", default_r_val=1.e-6_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    NULLIFY(subsection)
    CALL create_rsgrid_section(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    NULLIFY(subsection)
    CALL section_create(subsection,name="MULTIPOLES",&
         description="Enables the use of multipoles in the treatment of the electrostatics.",&
         n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
         citations=(/Aguado2003,Laino2008/), error=error)

    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="Controls the activation of the Multipoles",&
         usage="&MULTIPOLES T",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_MULTIPOLE_EXPANSION",&
         description="Specify the maximum level of multipoles expansion used "//&
         " for the electrostatics.",&
         usage="MAX_MULTIPOLE_EXPANSION DIPOLE",&
         enum_c_vals=s2a("NONE","CHARGE","DIPOLE","QUADRUPOLE"),&
         enum_desc=s2a("No multipolar terms! Check the codes providing a zero contribution.",&
                       "Use up to the Charge term", &
                       "Use up to the Dipole term",&
                       "Use up to the Quadrupole term"),&
         enum_i_vals=(/ do_multipole_none, do_multipole_charge, do_multipole_dipole,&
         do_multipole_quadrupole/), type_of_var=enum_t, error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="POL_SCF", &
         description="Specify the method to obtain self consistent induced "// &
                     "multipole moments.", &
         usage="POL_SCF CONJUGATE_GRADIENT", &
         enum_c_vals=s2a("NONE","SELF_CONSISTENT","CONJUGATE_GRADIENT"), &
         enum_desc=s2a("No inducible multipoles.", &
                       "Conventional self-consistent iteration.", &
                       "Linear conjugate-gradient optimization of the sum "//&
                       "of the electrostatic and induction energy. This "//&
                       "method does not support non-linear polarization "//&
                       "but is sometimes faster."), &
         enum_i_vals=(/ do_fist_pol_none, do_fist_pol_sc, do_fist_pol_cg/), &
         type_of_var=enum_t, default_i_val=do_fist_pol_none, error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_IPOL_ITER",&
         description="Specify the maximum number of iterations for induced "//&
         "dipoles",&
         usage="MAX_IPOL_ITER {int}",required=.FALSE., type_of_var=integer_t,&
            n_var=1, default_i_val=0,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_POL",&
         description="Specify the rmsd threshold for the derivatives "//&
         "of the energy towards the Cartesian dipoles components",&
         usage="EPS_POL {real}",required=.FALSE., type_of_var=real_t,&
            n_var=1, default_r_val=0.5e-07_dp,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    NULLIFY(subsection)
    CALL section_create(subsection,name="PRINT",&
         description="Controls printing of Ewald properties",&
         n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
         error=error)
    NULLIFY(print_key)
    CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
         description="controls the printing of ewald setup",&
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
         error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

  END IF
END SUBROUTINE create_ewald_section

! *****************************************************************************
!> \brief creates the interpolation section for the periodic QM/MM
!> \param section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_gspace_interp_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_gspace_interp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="interpolator",&
            description="controls the interpolation for the G-space term",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, print_key)

       CALL keyword_create(keyword, name="aint_precond",&
            description="the approximate inverse to use to get the starting point"//&
            " for the linear solver of the spline3 methods",&
            usage="kind spline3",&
            default_i_val=precond_spl3_aint,&
            enum_c_vals=s2a( "copy","spl3_nopbc_aint1","spl3_nopbc_precond1",&
            "spl3_nopbc_aint2","spl3_nopbc_precond2","spl3_nopbc_precond3"),&
            enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_1,&
            precond_spl3_aint2, precond_spl3_2, precond_spl3_3/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="precond",&
            description="The preconditioner used"//&
            " for the linear solver of the spline3 methods",&
            usage="kind spline3",&
            default_i_val=precond_spl3_3,&
            enum_c_vals=s2a("copy","spl3_nopbc_aint1","spl3_nopbc_precond1",&
            "spl3_nopbc_aint2","spl3_nopbc_precond2","spl3_nopbc_precond3"),&
            enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_1,&
            precond_spl3_aint2, precond_spl3_2, precond_spl3_3/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_x",&
            description="accuracy on the solution for spline3 the interpolators",&
            usage="eps_x 1.e-15", default_r_val=1.e-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_r",&
            description="accuracy on the residual for spline3 the interpolators",&
            usage="eps_r 1.e-15", default_r_val=1.e-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="max_iter",&
            variants=(/'maxiter'/),&
            description="the maximum number of iterations",&
            usage="max_iter 200", default_i_val=100, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"conv_info",&
            description="if convergence information about the linear solver"//&
            " of the spline methods should be printed", &
            print_level=medium_print_level,each_iter_names=s2a("SPLINE_FIND_COEFFS"),&
            each_iter_values=(/10/),filename="__STD_OUT__",&
            add_last=add_last_numeric,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_gspace_interp_section

! *****************************************************************************
!> \brief Creates the wavelet section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fschiff
!> \note
!>      this approach is based on the development of T. Deutsch and S. Goedecker
! *****************************************************************************
  SUBROUTINE create_wavelet_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_wavelet_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="wavelet",&
            description="Sets up parameters of  wavelet based poisson solver."//&
                  "This solver allows for non-periodic (PERIODIC NONE) boundary conditions and slab-boundary conditions "//&
                  "(but only PERIODIC XZ)."//&
                  "It does not require very large unit cells, only that the density goes to zero on the faces of the cell."//&
                  "The use of PREFERRED_FFT_LIBRARY FFTSG is required",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            citations=(/Genovese2006,Genovese2007/),&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="SCF_TYPE",&
            description="Type of scaling function used in the wavelet approach, the total energy depends on this choice,"//&
                        "and the convergence with respect to cutoff depends on the selected scaling functions."//&
                        "Possible values are 8,14,16,20,24,30,40,50,60,100 ",&
            usage="SCF_TYPE integer",&
            n_var=1,default_i_val=40,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_wavelet_section

! *****************************************************************************
!> \brief Creates the section for the generalized (iterative) poisson solver
!> \param section the section to be created
!> \param error cp2k error
!>
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE create_implicit_ps_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_implicit_ps_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="IMPLICIT",&
            description="Parameters for the implicit (generalized) Poisson solver.",&
            n_keywords=4, n_subsections=2, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(subsection, keyword)

       CALL create_dielectric_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection, error=error)

       CALL create_dbc_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection, error=error)

       CALL keyword_create(keyword, name="BOUNDARY_CONDITIONS", &
            enum_c_vals=s2a( 'PERIODIC', 'MIXED', 'MIXED_PERIODIC', 'NEUMANN' ), &
            enum_desc=s2a( 'periodic boundary conditions', 'Dirichlet + homogeneous Neumann boundary conditions', & 
                           'Dirichlet + periodic boundary conditions', 'homogeneous Neumann BC (zero-average solution)' ), &
            enum_i_vals=(/ PERIODIC_BC, MIXED_BC, MIXED_PERIODIC_BC, NEUMANN_BC /), &
            description="Specifies the type of boundary conditions. Dirichlet=fixed value, Neumann=zero normal deriv. "//&
                        "Mixed boundaries essentially requires FFTW3 so that all grid sizes are FFT-able.", &
            usage="BOUNDARY_CONDITIONS <bc_type>", default_i_val=PERIODIC_BC, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="ZERO_INITIAL_GUESS", &
            description="Whether or not to use zero potential as initial guess.", &
            usage="ZERO_INITIAL_GUESS <logical>", default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="max_iter", &
            description="Maximum number of iterations.", &
            usage="max_iter <integer>", default_i_val=30, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="tol", &
            description="Stopping tolerance.", &
            usage="tol <real>", default_r_val=1.0E-8_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="OR_PARAMETER", variants=s2a('omega'), &
            description="Over-relaxation parameter (large epsilon requires smaller omega ~0.1).", &
            usage="OR_PARAMETER <real>", default_r_val=1.0_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

    END IF
  END SUBROUTINE create_implicit_ps_section

! *****************************************************************************
!> \brief Creates the dielectric constant section.
!>  The dielectric constant is defined as a function of electronic density.
!>  [see O. Andreussi, I. Dabo, and N. Marzari, J. Chem. Phys., 136, 064102(2012)]
!> \param section the section to be created
!> \param error cp2k error
!>
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE create_dielectric_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_dielectric_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="DIELECTRIC", &
            description="Parameters for the dielectric constant function as proposed in the SCCS model.", &
            n_keywords=4, n_subsections=0, repeats=.FALSE., required=.FALSE., &
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="dielectric_constant", variants=s2a('epsilon'), &
            description="Dielectric constant in the bulk of the solvent.", &
            usage="dielectric_constant <real>", default_r_val=80.0_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="rho_min", &
            description="Lower density threshold.", &
            usage="rho_min <real>", default_r_val=1.0E-4_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="rho_max", &
            description="Upper density threshold.", &
            usage="rho_max <real>", default_r_val=1.0E-3_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="DERIVATIVE_METHOD", &
            enum_c_vals=s2a('fft', 'fft_use_deps', 'fft_use_drho', 'cd3', 'cd5', 'cd7'), &
            enum_i_vals=(/ derivative_fft, derivative_fft_use_deps, derivative_fft_use_drho, &
                           derivative_cd3, derivative_cd5, derivative_cd7 /), &
           enum_desc=s2a("FFT based deriv of epsilon, without correction (high cutoff needed).",&
                         "FFT based deriv of epsilon, with correction using gradient of epsilon (high cutoff needed).",&
                         "FFT based deriv of epsilon, with correction using gradient of rho (high cutoff needed).",&
                         "3-point central difference derivative.",&
                         "5-point central difference derivative.",&
                         "7-point central difference derivative (recommended)."),&
            description="Preferred method for evaluating the gradient of ln(eps).", &
            usage="DERIVATIVE_METHOD  <method>", default_i_val=derivative_cd7, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)
    END IF

  END SUBROUTINE create_dielectric_section

! *****************************************************************************
!> \brief Creates the section for Dirichlet boundary conditions
!> \param section the section to be created
!> \param error cp2k error
!>
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE create_dbc_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_dbc_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="DIRICHLET_BC", &
            description="Parameters for creating Dirichlet type boundary conditions.",&
            n_keywords=1, n_subsections=4, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="DEBUG", &
            description="Whether or not to print out the values of the electrostatic potential "// &
            "at Dirichlet regions, as well as the coordinates of the vertices of the region and "// &
            "its tessellation.", &
            usage="DEBUG <logical>", default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       NULLIFY(subsection)

       CALL create_aa_planar_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection, error=error)

       CALL create_planar_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection, error=error)

       CALL create_xaa_cylindrical_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection, error=error)

       CALL create_aa_cuboidal_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection, error=error)

    END IF
  END SUBROUTINE create_dbc_section

! *****************************************************************************
!> \brief Creates the section for creating axis-aligned planar Dirichlet BC.
!> \param section the section to be created
!> \param error cp2k error
!>
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE create_aa_planar_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_aa_planar_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="AA_PLANAR", &
            description="Parameters for creating axis-aligned planar (rectangular) Dirichlet boundary region.",&
            n_keywords=9, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="v_D", &
            description="value of the potential at the Dirichlet boundary (e.g. gate voltage).", &
            usage="v_D <real>", unit_str="volt", type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="PARALLEL_PLANE", &
            enum_c_vals=s2a( 'XY', 'YZ' , 'XZ'), &
            enum_i_vals=(/ xy_aligned_rectangle, yz_aligned_rectangle, xz_aligned_rectangle /), &
            description="The coordinate plane that the region is parallel to.", &
            usage="PARALLEL_PLANE <plane>", required=.TRUE., &
            type_of_var=enum_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="INTERCEPT", &
            description="The intercept of the rectangle's plane.", &
            usage="INTERCEPT <real>", unit_str="angstrom", required=.TRUE., &
            type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="X_xtnt", &
            description="The X extents of the plane.", &
            usage="X_xtnt <xmin(real)> <xmax(real)>", unit_str="angstrom", &
            n_var=2, type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="Y_xtnt", &
            description="The Y extents of the plane.", &
            usage="Y_xtnt <ymin(real)> <ymax(real)>", unit_str="angstrom", &
            n_var=2, type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="Z_xtnt", &
            description="The Z extents of the plane.", &
            usage="Z_xtnt <zmin(real)> <zmax(real)>", unit_str="angstrom", &
            n_var=2, type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="N_PRTN", &
            description="The number of partitions in xy, yz or xz directions for tiling the plane.", &
            usage="N_PRTN <integer> <integer>", &
            n_var=2, default_i_vals=(/1, 1/), error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="SMOOTH", &
            description="Whether or not to smoothen (mollify) the dirichlet boundary region to avoid "//&
            "possible Gibbs phenomenon.", &
            usage="SMOOTH <logical>", default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="smoothing_width", variants = s2a('zeta'),  &
            description="The width of the standard mollifier.", &
            usage="smoothing_width <real>", unit_str="angstrom", default_r_val=0.2_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

    END IF
  END SUBROUTINE create_aa_planar_section

! *****************************************************************************
!> \brief Creates the section for creating axis-aligned planar Dirichlet BC.
!> \param section the section to be created
!> \param error cp2k error
!>
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE create_planar_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_planar_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="PLANAR", &
            description="Parameters for creating an arbitrary planar Dirichlet boundary region.",&
            n_keywords=7, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="v_D", &
            description="value of the potential at the Dirichlet boundary (e.g. gate voltage).", &
            usage="v_D <real>", unit_str="volt", type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="A", &
            description="Coordinates of the vertex A.", &
            usage="A <x(real)> <y(real)> <z(real)>", unit_str="angstrom", &
            n_var=3, type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="B", &
            description="Coordinates of the vertex B.", &
            usage="B <x(real)> <y(real)> <z(real)>", unit_str="angstrom", &
            n_var=3, type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="C", &
            description="Coordinates of the vertex C.", &
            usage="C <x(real)> <y(real)> <z(real)>", unit_str="angstrom", &
            n_var=3, type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="N_PRTN", &
            description="The number of partitions in the x, y and z directions for tiling the plane.", &
            usage="N_PRTN <integer> <integer> <integer>", &
            n_var=3, default_i_vals=(/1, 1, 1/), error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="SMOOTH", &
            description="Whether or not to smoothen (mollify) the dirichlet boundary region to avoid "//&
            "possible Gibbs phenomenon.", &
            usage="SMOOTH <logical>", default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="smoothing_width", variants = s2a('zeta'),  &
            description="The width of the standard mollifier.", &
            usage="smoothing_width <real>", unit_str="angstrom", default_r_val=0.2_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

    END IF
  END SUBROUTINE create_planar_section

! *****************************************************************************
!> \brief Creates the section for creating x-axis-aligned cylindrical Dirichlet BC.
!> \param section the section to be created
!> \param error cp2k error
!>
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE create_xaa_cylindrical_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_xaa_cylindrical_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="XAA_CYLINDRICAL", &
            description="Parameters for creating x-axis-aligned cylindrical Dirichlet boundary region.",&
            n_keywords=9, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="v_D", &
            description="value of the potential at the Dirichlet boundary (e.g. gate voltage).", &
            usage="v_D <real>", unit_str="volt", type_of_var=real_t,  error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="X_xtnt", &
            description="The X extents of the cylinder.", &
            usage="X_xtnt <xmin(real)> <xmax(real)>", unit_str="angstrom", &
            n_var=2, type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="base_center", &
            description="The y and z coordinates of the cylinder's base center.", &
            usage="base_center <y(real)> <z(real)>", unit_str="angstrom", &
            n_var=2, type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="base_radius", &
            description="The base radius of the cylinder.", &
            usage="base_radius <real>", default_r_val=1.0_dp, unit_str="angstrom", error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="n_sides", &
            description="The number of sides os the n-gonal prism approximating the cylinder.", &
            usage="n_sides <integer>", default_i_val=5, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="apx_type", &
            enum_c_vals=s2a( 'CIRCUMSCRIBED', 'INSCRIBED' ), &
            enum_i_vals=(/ CIRCUMSCRIBED, INSCRIBED /), &
            description="Specifies the type of the n-gonal prism approximating the cylinder.", &
            usage="apx_type <apx_type>", default_i_val=CIRCUMSCRIBED, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="N_PRTN", &
            description="The number of partitions in the x, y and z directions for tiling the plane.", &
            usage="N_PRTN <integer> <integer> <integer>", &
            n_var=3, default_i_vals=(/1, 1, 1/), error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="SMOOTH", &
            description="Whether or not to smoothen (mollify) the dirichlet boundary region to avoid "//&
            "possible Gibbs phenomenon.", &
            usage="SMOOTH <logical>", default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="smoothing_width", variants = s2a('zeta'),  &
            description="The width of the standard mollifier.", &
            usage="smoothing_width <real>", unit_str="angstrom", default_r_val=0.2_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

    END IF
  END SUBROUTINE create_xaa_cylindrical_section

! *****************************************************************************
!> \brief Creates the section for creating axis-aligned cuboidal Dirichlet region.
!> \param section the section to be created
!> \param error cp2k error
!>
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE create_aa_cuboidal_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_aa_cuboidal_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="AA_CUBOIDAL", &
            description="Parameters for creating axis-aligned cuboidal region where a fixed voltage is applied.",&
            n_keywords=7, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="v_D", &
            description="value of the potential at the Dirichlet boundary (e.g. gate voltage).", &
            usage="v_D <real>", unit_str="volt", type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="X_xtnt", &
            description="The X extents of the cuboid.", &
            usage="X_xtnt <xmin(real)> <xmax(real)>", unit_str="angstrom", &
            n_var=2, type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="Y_xtnt", &
            description="The Y extents of the cuboid.", &
            usage="Y_xtnt <ymin(real)> <ymax(real)>", unit_str="angstrom", &
            n_var=2, type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="Z_xtnt", &
            description="The Z extents of the cuboid.", &
            usage="Z_xtnt <zmin(real)> <zmax(real)>", unit_str="angstrom", &
            n_var=2, type_of_var=real_t, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="N_PRTN", &
            description="The number of partitions in the x, y and z directions for partitioning the cuboid.", &
            usage="N_PRTN <integer> <integer> <integer>", &
            n_var=3, default_i_vals=(/1, 1, 1/), error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="SMOOTH", &
            description="Whether or not to smoothen (mollify) the dirichlet region to avoid "//&
            "possible Gibbs phenomenon.", &
            usage="SMOOTH <logical>", default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="smoothing_width", variants = s2a('zeta'),  &
            description="The width of the standard mollifier.", &
            usage="smoothing_width <real>", unit_str="angstrom", default_r_val=0.2_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

    END IF
  END SUBROUTINE create_aa_cuboidal_section

END MODULE input_cp2k_poisson
