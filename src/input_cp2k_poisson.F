!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_poisson *
!!
!!   NAME
!!     input_cp2k_poisson
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
MODULE input_cp2k_poisson
  USE bibliography,                    ONLY: Darden1993,&
                                             Essmann1995,&
                                             Ewald1921,&
                                             Genovese2006,&
                                             Genovese2007,&
                                             Martyna1999,&
                                             Toukmaji1996
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
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_poisson'

  PUBLIC :: create_poisson_section,&
            create_gspace_interp_section,&
            create_multipole_section,&
            create_ewald_section
!***
!****************************************************************************
CONTAINS

!!****f* input_cp2k_poisson/create_poisson_section [1.0] *
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
            enum_c_vals=s2a( "PERIODIC", "ANALYTIC", "MT", "MULTIPOLE","WAVELET"),&
            enum_i_vals=(/ use_periodic, use_analytic, use_mt, use_multipole,use_wavelet/),&
            citations=(/Martyna1999/),&
            default_i_val=use_periodic, error=error)
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

    END IF
  END SUBROUTINE create_poisson_section
!***************************************************************************

!***************************************************************************
!!****f* input_cp2k_poisson/create_multipole_section [1.0] *
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

!!****f* input_cp2k_dft/create_mt_section [1.0] *
!!
!!   NAME
!!     create_mt_section
!!
!!   FUNCTION
!!     Creates the Martyna-Tuckerman section
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
         citations=(/Ewald1921,Darden1993,Essmann1995,Toukmaji1996/),&
         error=error)

    NULLIFY(keyword,print_key,subsection)
    CALL keyword_create(keyword, name="EWALD_TYPE",&
         description="The type of ewald you want to perform."//&
         " NONE standard real-space coulomb potential is computed together"//&
         " with the non-bonded contributions."//&
         " EWALD is the standard non-fft based ewald."//&
         " SPME is the smooth particle mesh using beta-Euler splines (recommended)."//&
         " PME is the particle mesh using fft interpolation.",&
         citations=(/Ewald1921,Essmann1995,Darden1993/),&
         usage="EWALD_TYPE {EWALD|SPME|PME|NONE}",&
         default_i_val=do_ewald_ewald,&
         enum_c_vals=(/"none     ",&
                       "ewald    ",&
                       "pme      ",&
                       "spme     " /),&
         enum_i_vals=(/do_ewald_none,&
                       do_ewald_ewald,&
                       do_ewald_pme,&
                       do_ewald_spme/),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="alpha",&
         description="alpha parameter associated with Ewald (EWALD|PME|SPME). Recommended is alpha = 3.5 / r_cut.",&
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

    CALL keyword_create(keyword, name="RS_GRID",&
         description="Parallellisation strategy of the realspace grids",&
         usage="RS_GRID DISTRIBUTED",&
         enum_c_vals=s2a("DISTRIBUTED","REPLICATED"),&
         enum_desc=s2a("Better memory usage and parallel scalability.", &
                       "Has been tested more extensively (old default)."),&
         enum_i_vals=(/ do_rsgrid_distributed, do_rsgrid_replicated/),&
         default_i_val=do_rsgrid_distributed, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    NULLIFY(subsection)
    CALL section_create(subsection,name="PRINT",&
         description="Controls printing of Ewald properties",&
         n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
         error=error)
    NULLIFY(print_key)
    CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
         description="controls the printing of ewald setup",&
         print_level=low_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
         error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

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
            print_level=medium_print_level,each=(/10/),filename="__STD_OUT__",&
            add_last=(/1/),error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_gspace_interp_section

!***************************************************************************

!!****f* input_cp2k_poisson/create_wavelet_section [1.0] *
!!
!!   NAME
!!     create_wavelet_section
!!
!!   FUNCTION
!!     Creates the wavelet section
!!
!!   NOTES
!!     this approach is based on the development of T. Deutsch and S. Goedecker
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fschiff
!!
!!*** **********************************************************************
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
                  "It does not require very large unit cells, only that the density goes to zero on the faces of the cell.",&
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
!**************************************************************************

END MODULE input_cp2k_poisson
