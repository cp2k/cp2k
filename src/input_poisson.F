!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_poisson *
!!
!!   NAME
!!     input_poisson
!!
!!   FUNCTION
!!     function that build the poisson section of the input
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     03.2006 fusing of poisson_dft and poisson_mm
!!
!!   SOURCE
!****************************************************************************
MODULE input_poisson
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_constants
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: integer_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_poisson'

  PUBLIC :: create_poisson_section,&
            create_gspace_interp_section,&
            create_multipole_section
!***
!****************************************************************************
CONTAINS

!!****f* input_poisson/create_poisson_section [1.0] *
!!
!!   NAME
!!     create_poisson_section
!!
!!   FUNCTION
!!     Creates the Poisson section
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_poisson_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
            enum_c_vals=s2a( "PERIODIC", "ANALYTIC", "MT", "MULTIPOLE"),&
            enum_i_vals=(/ use_periodic, use_analytic, use_mt, use_multipole/),&
            default_i_val=use_periodic, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL keyword_create(keyword, name="PERIODIC",&
            description="Specify the directions on wich apply PBC ",&
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

       CALL create_multipole_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_screening_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_ewald_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_poisson_section
!***************************************************************************

!***************************************************************************
!!****f* input_poisson/create_multipole_section [1.0] *
!!
!!   NAME
!!     create_multipole_section
!!
!!   FUNCTION
!!     Section to set-up parameters for decoupling using the Bloechl scheme
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** ********************************************************************** 
  SUBROUTINE create_multipole_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

       CALL keyword_create(keyword, name="NUM_GAUSS",&
            description="Specifies the numbers of gaussian used to fit the QM density for each atomic site.",&
            usage="NUM_GAUSS {integer}", required=.FALSE.,&
            n_var=1, type_of_var=integer_t, default_i_val=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PFACTOR",&
            description="Specifies the progression factor for the gaussian exponent for each atomic site.",&
            usage="PFACTOR {real}", required=.FALSE.,&
            n_var=1, type_of_var=real_t, default_r_val=1.5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MIN_RADIUS",&
            description="Specifies the smallest radius of the gaussian used in the fit. All other radius are"//&
            " obtained with the progression factor.",&
            usage="MIN_RADIUS {real}", required=.FALSE.,&
            unit_str="angstrom",n_var=1, type_of_var=real_t, default_r_val=0.5_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="RADII",&
            description="Specifies all the radius of the gaussian used in the fit for each atomic site. The use"//&
            " of this keyword disables all other keywords of this section.",&
            usage="RADII {real} {real} .. {real}", required=.FALSE.,&
            unit_str="angstrom",n_var=-1, type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GCUT",&
            description="Cutoff for charge fit in G-space.",&
            usage="GCUT {real}", required=.FALSE.,&
            n_var=1, type_of_var=real_t, default_r_val=SQRT(6.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

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
            print_level=medium_print_level,each=(/1/),filename="GSpace-SplInterp",&
            error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)       

       CALL cp_print_key_section_create(subsection,"program_run_info",&
            description="Controls the printing of basic information during the run", &
            print_level=low_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_multipole_section

!***************************************************************************
!!****f* input_poisson/create_screening_section [1.0] *
!!
!!   NAME
!!     create_screening_section
!!
!!   FUNCTION
!!     Section to set-up parameters for screening of short range Green
!!     function
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     jgh
!!
!!*** **********************************************************************
  SUBROUTINE create_screening_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_screening_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="SCREENING",&
            description="This section is used to set up a screening function for the short "//&
            "range part of the interaction potential.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="SR_SCREENING",&
            description="Use a short range screening function.",&
            usage="SR_SCREENING", required=.FALSE.,&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RC_INTERACTION",&
            description="The range of the excluded interaction.",&
            usage="RC_INTERACTION {real}", required=.FALSE.,&
            n_var=1, type_of_var=real_t, default_r_val=0._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RC_RANGE",&
            description="The range of the switching function.",&
            usage="RC_INTERACTION {real}", required=.FALSE.,&
            n_var=1, type_of_var=real_t, default_r_val=8._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_screening_section

!***************************************************************************

!!****f* input_cp2k_dft/create_mt_section [1.0] *
!!
!!   NAME
!!     create_mt_section
!!
!!   FUNCTION
!!     Creates the Martyna-Tukerman section
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_mt_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mt_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="mt",&
            description="Sets up parameters of  Martyna-Tuckerman poisson solver",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
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
!***************************************************************************

!!****f* input_cp2k_mm/create_ewald_section *
!!
!!   NAME
!!     create_ewald_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the ewald section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************

SUBROUTINE create_ewald_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ewald_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

  failure=.FALSE.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
    CALL section_create(section,name="ewald",&
         description="ewald parameters",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="EWALD_TYPE",&
         description="The type of ewald you want to perform."//&
         " EWALD is the standard non-fft based ewald."//&
         " SPME is the smooth particle mesh using beta-Euler splines."//&
         " PME is the particle mesh using fft interpolation.",&
         usage="EWALD_TYPE {EWALD|SPME|PME}",&
         default_i_val=do_ewald_ewald,&
         enum_c_vals=ewald_name_types,&
         enum_i_vals=ewald_name_index,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="alpha",&
         description="alpha parameter associated with Ewald (EWALD|PME|SPME)",&
         usage="alpha .44", default_r_val=.44_dp, unit_str='angstrom^-1',error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="gmax",&
         description="maximum number of grid points (SPME and EWALD)",&
         usage="gmax [integer]", default_i_val=25,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ns_max",&
         description="number of grid points on small mesh (PME only)",&
         usage="ns_max 12", default_i_val=12,error=error)
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

  END IF
END SUBROUTINE create_ewald_section
!**************************************************************************
!!****f* input_cp2k_qmmm/create_gspace_interp_section *
!!
!!   NAME
!!     create_gspace_interp_section
!!
!!   FUNCTION
!!     creates the interpolation section for the periodic QM/MM 
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     tlaino
!!
!!*** **********************************************************************
  SUBROUTINE create_gspace_interp_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
            print_level=medium_print_level,each=(/10/),filename="__STD_OUT__",&
            add_last=(/1/),error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_gspace_interp_section
!**************************************************************************

END MODULE input_poisson
