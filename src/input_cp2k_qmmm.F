!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_qmmm *
!!
!!   NAME
!!     input_cp2k_qmmm
!!
!!   FUNCTION
!!     creates the qmmm section of the input
!!
!!   NOTES
!!     moved out of input_cp2k
!!
!!   AUTHOR
!!     teo & fawzi
!!
!!   MODIFICATION HISTORY
!!     10.2005 split out of input_cp2k
!!
!!   SOURCE
!****************************************************************************
MODULE input_cp2k_qmmm
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE input_constants
  USE input_poisson,                   ONLY: create_gspace_interp_section,&
                                             create_multipole_section
  USE input_cp2k_mm,                   ONLY: create_Goodwin_section,&
                                             create_LJ_section,&
                                             create_NONBONDED14_section,&
                                             create_Williams_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: char_t,&
                                             integer_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_qmmm'
  
  PUBLIC :: create_qmmm_section

!***
!****************************************************************************
CONTAINS
!!****f* input_cp2k_qmmm/create_qmmm_section [1.0] *
!!
!!   NAME
!!     create_qmmm_section
!!
!!   FUNCTION
!!     Creates the QM/MM section
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
  SUBROUTINE create_qmmm_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qmmm_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="qmmm",&
            description="Input for QM/MM calculations.",&
            n_keywords=6, n_subsections=3, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, subsection )

       CALL keyword_create(keyword, name="MM_POTENTIAL_FILE_NAME",&
            description="Name of the file containing the potential expansion in gaussians.",&
            usage="MM_POTENTIAL_FILE_NAME {filename}",&
            default_lc_val="MM_POTENTIAL",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="A_QM",&
            description="Set the A vector of the QM Cell to the Reals value. The other values will be"//&
            "referred to the MM cell. It Changes the shape of the QM cell in one dimension.",&
            usage="A_QM real real real",&
            n_var=3,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B_QM",&
            description="Set the B vector of the QM Cell to the Reals value. The other values will be"//&
            "referred to the MM cell. It Changes the shape of the QM cell in one dimension.",&
            usage="B_QM real real real",&
            n_var=3,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="C_QM",&
            description="Set the C vector of the QM Cell to the Reals value. The other values will be"//&
            "referred to the MM cell. It Changes the shape of the QM cell in one dimension.",&
            usage="C_QM real real real",&
            n_var=3,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ABC_QM",&
            description="Set the A,B,C vectors of  QM Cell to the Reals value."//&
            "It changes the overall shape of the QM cell.",&
            usage="ABC_QM real real real",&
            n_var=3,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_mm_rspace",&
            description="Set the threshold g(r) - threshold = 0. Where g(r) is a gaussian function."//&
            "Standard value set to 1.0E-10_dp",&
            usage="eps_mm_rspace real",&
            default_r_val=1.0E-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)     

       CALL keyword_create(keyword, name="nocenter",&
            description="This keyword disables the automatic centering of the qm system every MD step."//&
            "It centers the system only for the first step. Default .FALSE. (i.e. the QM system is centered every step)",&
            usage="nocenter LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="nocenter0",&
            description="This keyword disables the automatic centering of the qm system every MD step "//&
            " even for the first step. Default .FALSE. (i.e. the QM system is centered every step)",&
            usage="nocenter LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="nocompatibility",&
            description="This keyword disables the compatibility of QM/MM "//&
            "potential between CPMD and CP2K implementations. The compatibility"//&
            " is achieved using an MM potential of the form: Erf[x/rc]/x + (1/rc -2/(pi^1/2*rc))*Exp[-(x/rc)^2] .",&
            usage="nocompatibility LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="use_geep_lib",&
            description=" This keyword enables the use of the internal GEEP library to generate the "//&
            "gaussian expansion of the MM potential. Using this keyword there's no need to provide "//&
            "the MM_POTENTIAL_FILENAME. It expects a number from 2 to 15 (the number of gaussian funtions"//&
            " to be used in the expansion.",&
            usage="use_geep_lib INTEGER",&
            default_i_val=0, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="typ_center",&
            description="This keyword specifies whether the centering is performes continuosly (NONE) "//&
            "or by unit grid steps (GRID).",&
            usage="typ_center (NONE|GRID)",&
            enum_c_vals=s2a( "NONE","GRID"),&
            enum_i_vals=(/ do_center_none, do_center_grid /),&
            default_i_val=do_center_none, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)     

       CALL keyword_create(keyword, name="E_COUPL",&
            variants=(/"QMMM_COUPLING",&
            "ECOUPL       "/),& 
            description="Set the Coupling's type between QM and MM system. Possible Options:"//&
            "(1) NONE: Just Mechanical Coupling. (2) GRID: Computing the electrostatic"//&
            "potential analytically on the Grid. (3) GAUSS: Mapping the MM electrostatic"//&
            "potential on the Grid with Gaussian.",&
            usage="E_COUPL char",&
            enum_c_vals=s2a( "NONE","GRID","GAUSS"),&
            enum_i_vals=(/ do_qmmm_none, do_qmmm_grid, do_qmmm_gauss/),&
            default_i_val=do_qmmm_none, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL section_create(subsection,name="QM_KIND",&
            description="Information about the qm kind in the qm/mm scheme",&
            n_keywords=3, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="The qm kind",&
            usage="O",&
            n_var=1,type_of_var=char_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MM_INDEX",&
            description="The indexes of the mm atoms that have this kind. This keyword can be repeated several times "//&
            "(useful if you have to specify many indexes).",&
            usage="MM_INDEX 1 2",&
            n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SUBSYS",&
            description="Specifies which subsystem the qm atom belongs to.",&
            usage="SUBSYS integer",n_var=1,type_of_var=integer_t,&
            default_i_val=1,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)    

       CALL section_create(subsection,name="MM_KIND",&
            description="Information about the mm kind in the qm/mm scheme",&
            n_keywords=2, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="The mm kind",&
            usage="O",&
            n_var=1,type_of_var=char_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RADIUS",&
            description="Specifies the radius of the atomic kinds [angstrom]",&
            usage="RADIUS real",n_var=1,type_of_var=real_t,&
            default_r_val=RADIUS_QMMM_DEFAULT,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORR_RADIUS",&
            description="Specifies the correction radius of the atomic kinds [angstrom]",&
            usage="RADIUS real",n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)    

       CALL create_qmmm_periodic_section(subsection, error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)       

       CALL create_qmmm_link_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_qmmm_interp_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_qmmm_forcefield_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)       

       CALL create_print_qmmm_section(subsection, error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_qmmm_section

!!****f* input_cp2k_qmmm/create_qmmm_periodic_section *
!!
!!   NAME
!!     create_qmmm_periodic_section
!!
!!   FUNCTION
!!     Input section to set QM/MM periodic boundary conditions
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
  SUBROUTINE create_qmmm_periodic_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qmmm_periodic_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    NULLIFY(keyword, subsection)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="PERIODIC",&
            description="Specify parameters for QM/MM periodic boundary conditions calculations",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)    

       CALL keyword_create(keyword, name="GMAX",&
            description="Specifies the maximum value of G in the reciprocal space over which perform the Ewald sum.",&
            usage="GMAX <real>",n_var=1,default_r_val=1.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="REPLICA",&
            description="Specifies the number of replica to take into consideration for the real part of the "//&
            "calculation. Default is letting the qmmm module decide how many replica you really need.",&
            usage="REPLICA <integer>",n_var=1,default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       
       CALL keyword_create(keyword, name="NGRIDS",&
            description="Specifies the number of grid points used for the Interpolation of the G-space term",&
            usage="NGRIDS <integer> <iteger> <integer> ",n_var=3,default_i_vals=(/50,50,50/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL create_multipole_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)       

       CALL create_gspace_interp_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)       

       CALL cp_print_key_section_create(subsection,"check_spline",&
            description="Controls the checking of the G-space term Spline Interpolation.",&
            print_level=medium_print_level,each=(/1/),filename="GSpace-SplInterp",&
            error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)       

    END IF
  END SUBROUTINE create_qmmm_periodic_section

!!****f* input_cp2k_qmmm/create_qmmm_forcefield_section *
!!
!!   NAME
!!     create_qmmm_forcefield_section
!!
!!   FUNCTION
!!     creates the qm/mm forcefield section to override to the FF specification
!!     given in the FIST input
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
  SUBROUTINE create_qmmm_forcefield_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_qmmm_forcefield_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    NULLIFY(subsection )
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="FORCEFIELD",&
            description="Specify information on the QM/MM forcefield",&
            n_keywords=0, n_subsections=2, repeats=.TRUE., required=.TRUE.,&
            error=error)    

       CALL create_qmmm_ff_nb_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_NONBONDED14_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_qmmm_forcefield_section

!!****f* input_cp2k_qmmm/create_qmmm_ff_nb_section *
!!
!!   NAME
!!     create_qmmm_ff_nb_section
!!
!!   FUNCTION
!!     creates the qm/mm forcefield section to override to the FF specification
!!     given in the FIST input - NONBONDED PART
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
  SUBROUTINE create_qmmm_ff_nb_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qmmm_ff_nb_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    NULLIFY(subsection )
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="NONBONDED",&
            description="Specify information on the QM/MM non-bonded forcefield",&
            n_keywords=0, n_subsections=2, repeats=.TRUE., required=.TRUE.,&
            error=error)    

       CALL create_LJ_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
       
       CALL create_Williams_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Goodwin_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_qmmm_ff_nb_section

!!****f* input_cp2k_qmmm/create_qmmm_link_section *
!!
!!   NAME
!!     create_qmmm_link_section
!!
!!   FUNCTION
!!     creates the qm/mm link section
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
  SUBROUTINE create_qmmm_link_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qmmm_link_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    NULLIFY(keyword, subsection )
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="LINK",&
            description="Specify information on the QM/MM link treatment",&
            n_keywords=7, n_subsections=2, repeats=.TRUE., required=.TRUE.,&
            error=error)    
       
       CALL keyword_create(keyword, name="QM_INDEX",&
            variants=(/ "QM" /),&
            description="Specifies the index of the QM atom involved in the QM/MM link",&
            usage="QM_INDEX integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="QM_KIND",&
            description="Specifies the element of the QM capping atom involved in the QM/MM link",&
            usage="QM_KIND char",n_var=1,type_of_var=char_t,&
            default_c_val="H",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MM_INDEX",&
            variants=(/ "MM" /),&
            description="Specifies the index of the MM atom involved in the QM/MM link, Default hydrogen.",&
            usage="MM_INDEX integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RADIUS",&
            description="Overwrite the specification of the radius only for the MM atom involved in the link."//&
            "Default is to use the same radius as for the specified type.",&
            usage="RADIUS real",n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORR_RADIUS",&
            description="Overwrite the specification of the correction radius only for the MM atom involved in the link."//&
            "Default is to use the same correction radius as for the specified type.",&
            usage="RADIUS real",n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

        CALL keyword_create(keyword, name="LINK_TYPE",&
            variants=(/ "LINK ","TYPE ","LTYPE"/),&
            description="Specifies the method to use to treat the defined QM/MM link",&
            usage="LINK_TYPE char",&
            enum_c_vals=s2a( "IMOMM","GHO","PSEUDO"),&
            enum_i_vals=(/ do_qmmm_link_imomm, do_qmmm_link_gho, do_qmmm_link_pseudo/),&
            default_i_val=do_qmmm_link_imomm, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)      

       CALL keyword_create(keyword, name="ALPHA_IMOMM",&
            variants=s2a("ALPHA"),&
            description="Specifies the scaling factor to be used in the IMOMM QM/MM link scheme",&
            usage="ALPHA_IMOMM real",n_var=1,type_of_var=real_t,&
            default_r_val=ALPHA_IMOMM_DEFAULT,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)      

       CALL keyword_create(keyword, name="QMMM_SCALE_FACTOR",&
            variants=(/ "QMMM_CHARGE_SCALE ",&
                        "QMMM_CHARGE_FACTOR",&
                        "QMMM_SCALE_CHARGE " /),&
            description="Specifies the scaling factor for the MM charge involved in the link QM/MM."//&
            "This keyword affects only the QM/MM potential, it doesn't affect the electrostatic in "//&
            "the classical part of the code"//&
            "Default 1.0 i.e. delete the charge of the MM atom of the QM/MM link bond.",&
            usage="SCALE_FACTOR real",n_var=1,type_of_var=real_t,&
            default_r_val=CHARGE_SCALE_FACTOR,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FIST_SCALE_FACTOR",&
            variants=(/ "FIST_CHARGE_SCALE ",&
                        "FIST_CHARGE_FACTOR",&
                        "FIST_SCALE_CHARGE " /),&
            description="Specifies the scaling factor for the MM charge involved in the link QM/MM."//&
            "This keyword modifies the MM charge in FIST. The modified charge will be used then also"//&
            "for the generation of the QM/MM potential."//&
            "Default 1.0 i.e. delete the charge of the MM atom of the QM/MM link bond.",&
            usage="SCALE_FACTOR real",n_var=1,type_of_var=real_t,&
            default_r_val=CHARGE_SCALE_FACTOR,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="SUBSYS",&
            description="Specifies which subsystem the qm/mm link atoms belong to.",&
            usage="SUBSYS integer",n_var=1,type_of_var=integer_t,&
            default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL section_create(subsection,name="MOVE_MM_CHARGE",&
            description="Specify information to move a classical charge before the"//&
            " QM/MM energies and forces evaluation",&
            n_keywords=4, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)      

       CALL keyword_create(keyword, name="ATOM_INDEX_1",&
            variants=(/ "MM1" /),&
            description="Specifies the index of the MM atom involved in the QM/MM link to be moved",&
            usage="ATOM_INDEX_1 integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOM_INDEX_2",&
            variants=(/ "MM2" /),&
            description="Specifies the index of the second atom defining the direction along which "//&
            " the atom will be moved",&
            usage="ATOM_INDEX_2 integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALPHA",&
            description="Specifies the scaling factor that defines the movement along the defined direction",&
            usage="ALPHA real",n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)      

       CALL keyword_create(keyword, name="RADIUS",&
            description="Specifies the radius used for the QM/MM electrostatic coupling after movement",&
            usage="RADIUS real",n_var=1,type_of_var=real_t, default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORR_RADIUS",&
            description="Specifies the correction radius used for the QM/MM electrostatic coupling after movement",&
            usage="RADIUS real",n_var=1,type_of_var=real_t, default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)


       CALL section_create(subsection,name="ADD_MM_CHARGE",&
            description="Specify information to add a classical charge before the"//&
            " QM/MM energies and forces evaluation",&
            n_keywords=5, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)      

       CALL keyword_create(keyword, name="ATOM_INDEX_1",&
            variants=(/ "MM1" /),&
            description="Specifies the index of the first atom defining the direction along which"//&
            " the atom will be added",&
            usage="ATOM_INDEX_1 integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOM_INDEX_2",&
            variants=(/ "MM2" /),&
            description="Specifies the index of the second atom defining the direction along which "//&
            " the atom will be added",&
            usage="ATOM_INDEX_2 integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALPHA",&
            description="Specifies the scaling factor that defines the movement along the defined direction",&
            usage="ALPHA real",n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)      

       CALL keyword_create(keyword, name="RADIUS",&
            description="Specifies the radius used for the QM/MM electrostatic coupling for the added source",&
            usage="RADIUS real",default_r_val=RADIUS_QMMM_DEFAULT,n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORR_RADIUS",&
            description="Specifies the correction radius used for the QM/MM electrostatic coupling for the added source",&
            usage="RADIUS real",default_r_val=RADIUS_QMMM_DEFAULT,n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE",&
            description="Specifies the charge for the added source of QM/MM potential",&
            usage="CHARGE real",default_r_val=0.0_dp,n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_qmmm_link_section

!!****f* input_cp2k_qmmm/create_qmmm_interp_section *
!!
!!   NAME
!!     create_qmmm_interp_section
!!
!!   FUNCTION
!!     creates the interpolation section
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
  SUBROUTINE create_qmmm_interp_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qmmm_interp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="interpolator",&
            description="kind of interpolation used between the multigrids",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, print_key)

       CALL keyword_create(keyword, name="kind",&
            description="the interpolator to use",&
            usage="kind spline3",&
            default_i_val=spline3_nopbc_interp,&
            enum_c_vals=s2a("linear","spline3_nopbc"),&
            enum_i_vals=(/linear_interp,&
            spline3_nopbc_interp/),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="safe_computation",&
            description="if a non unrolled calculation is to be performed in parallel",&
            usage="safe_computation OFF",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

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

       CALL cp_print_key_section_create(print_key,"spl_coeffs",&
            description="outputs a cube with the coefficents calculated for "//&
            "the spline interpolation", &
            print_level=debug_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
    END IF
  END SUBROUTINE create_qmmm_interp_section

!***************************************************************************
!!****f* input_cp2k_qmmm/create_print_qmmm_section [1.0] *
!!
!!   NAME
!!     create_print_qmmm_section
!!
!!   FUNCTION
!!     Create the print qmmm section
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
  SUBROUTINE create_print_qmmm_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_qmmm_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options specific of the QMMM code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key)

       CALL cp_print_key_section_create(print_key,"PGF",&
            description="Controls the printing of the gaussian expansion basis set of the"//&
            " electrostatic potential", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"POTENTIAL",&
            description="Controls the printing of the QMMM  potential",&
            print_level=high_print_level,each=(/1/),filename="MM_ELPOT_QMMM",&
            common_iter_levels=1,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MM_POTENTIAL",&
            description="Controls the printing of the MM unidimensional potential on file",&
            print_level=high_print_level,each=(/1/),filename="MM_ELPOT.data",&
            common_iter_levels=1,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_BANNER",&
            description="Controls the printing of the banner of the MM program",&
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of information regarding the run.",&
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PERIODIC_INFO",&
            description="Controls the printing of information regarding the periodic boundary condition.",&
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"derivatives",&
            description="Print all derivatives after QM/MM calculation", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"qmmm_charges",&
            description="Print all charges generating the QM/MM potential", &
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"qmmm_link_info",&
            description="Print all information on QM/MM links", &
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
 
       CALL cp_print_key_section_create(print_key,"qs_derivatives",&
            description="Print QM derivatives after QS calculation", &
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_qmmm_section
!***************************************************************************

END MODULE input_cp2k_qmmm
