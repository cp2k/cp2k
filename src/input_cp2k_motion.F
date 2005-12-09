!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_motion *
!!
!!   NAME
!!     input_cp2k_motion
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     teo & fawzi
!!
!!   MODIFICATION HISTORY
!!     10.2005 split input_cp2k into smaller modules [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE input_cp2k_motion
  use kinds, only: dp
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_vals_type,&
       section_vals_get_subs_vals, section_vals_get
  USE input_constants
  USE input_keyword_types, ONLY: keyword_type, keyword_create,&
       keyword_release
  USE input_val_types, ONLY: val_type, no_t, logical_t, integer_t, real_t, char_t,&
       enum_t
  USE cp_output_handling, only: cp_print_key_section_create
  USE cp_units, only: cp_unit_to_cp2k
  USE string_utilities, ONLY: s2a
#include "cp_common_uses.h"
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='input_cp2k_motion'

PUBLIC :: create_motion_section

!***
!****************************************************************************
contains

!!****f* input_cp2k/create_metadyn_section *
!!
!!   NAME
!!     create_metadyn_section
!!
!!   FUNCTION
!!     creates the metadynamics section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_motion_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_motion_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="motion",&
            description="This section defines a set of tool connected with the motion of the nuclei.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)

       CALL create_geoopt_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_md_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_metadyn_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
       
       CALL create_mc_section(subsection, error)
       call section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
       
    END IF

  END SUBROUTINE create_motion_section

!!****f* input_cp2k/create_mc_section *
!!
!!   NAME
!!     create_mc_section
!!
!!   FUNCTION
!!     creates the metadynamics section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_mc_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_mc_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="mc",&
            description="This section sets parameters to set up a MonteCarlo calculation.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)
        
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="NSTEP",&
            description="Specifies the number of MC cycles.",&
            usage="NSTEP {integer}",&
            default_i_val=100,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NMOVES",&
            description="Specifies the number of classical moves between energy evaluations. ",&
            usage="NMOVES {integer}",&
            default_i_val=4,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="NVOLMOVES",&
            description="How many classical volume moves to bias with.",&
            usage="NVOLMOVES {integer}",&
            default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="NSWAPMOVES",&
            description="How many classical swap moves to try per QS swap move.",&
            usage="NSWAPMOVES {integer}",&
            default_i_val=16,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)      

       CALL keyword_create(keyword, name="LBIAS",&
            description="Dictates if we presample moves with a different potential.",&
            usage="LBIAS {logical}",&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LSTOP",&
            description="Makes nstep in terms of steps, instead of cycles.",&
            usage="LSTOP {logical}",&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LCLASSICAL",&
            description="Bias the moves with an inexpensive classical potential.",&
            usage="LCLASSICAL {logical}",&
            default_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="RESTART",&
            description="Read initial configuration from restart file.",&
            usage="RESTART {logical}",&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="IUPVOLUME",&
            description="Every iupvolume steps update maximum volume displacement.",&
            usage="IUPVOLUME {integer}",&
            default_i_val=10000,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="IUPTRANS",&
            description="Every iuptrans steps update maximum translation.",&
            usage="IUPTRANS {integer}",&
            default_i_val=10000,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="PROGRAM",&
            description="Specify the type of simulation",&
            usage="PROGRAM (TRADITIONAL|GEMC)",&
            enum_c_vals=s2a( "TRADITIONAL","GEMC"),&
            enum_i_vals=(/do_mc_traditional,do_mc_gemc/),&
            default_i_val=do_mc_traditional,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       
       CALL keyword_create(keyword, name="RESTART_FILE_NAME",&
            description="Name of the restart file for MC information.",&
            usage="RESTART_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COORDINATE_FILE_NAME",&
            description="The file to print the current coordinates to.",&
            usage="COORDINATE_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ENERGY_FILE_NAME",&
            description="The file to print current energies to.",&
            usage="ENERGY_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DATA_FILE_NAME",&
            description="The file to print current configurational info to.",&
            usage="DATA_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CELL_FILE_NAME",&
            description="The file to print current cell length info to.",&
            usage="CELL_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_DISP_FILE_NAME",&
            description="The file to print current maximum displacement info to.",&
            usage="MAX_DISP_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BOX2_FILE_NAME",&
            description="For GEMC, the name of the input file for the other box.",&
            usage="BOX2_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PRESSURE",&
            description="The pressure for NpT simulations, in bar.",&
            usage="PRESSURE {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TEMP",&
            description="The temperature of the simulation, in Kelvin.",&
            usage="TEMP {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMSWAP",&
            description="The probability of attempting a swap move.",&
            usage="PMSWAP {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMVOLUME",&
            description="The probability of attempting a volume move.",&
            usage="PMVOLUME {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMTRAION",&
            description="The probability of attempting a conformational change.",&
            usage="PMTRAION {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMTRANS",&
            description="The probability of attempting a molecule translation.",&
            usage="PMTRANS {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMVOLUME",&
            description="Maximum volume displacement, in angstrom**3.",&
            usage="RMVOLUME {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMBOND",&
            description="Maximum bond length displacement, in angstrom.",&
            usage="RMBOND {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMANGLE",&
            description="Maximum bond angle distplacement, in degrees.",&
            usage="RMANGLE {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMROT",&
            description="Maximum rotational displacement, in degrees.",&
            usage="RMROT {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMTRANS",&
            description="Maximum translational displacement, in angstroms.",&
            usage="RMTRANS {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_mc_section

!!****f* input_cp2k/create_geoopt_section *
!!
!!   NAME
!!     create_geoopt_section
!!
!!   FUNCTION
!!     creates the geometry optimization section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_geoopt_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_geoopt_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="geoopt",&
            description="This section sets the environment of the geometry optimizer.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="OPTIMIZER",&
            variants=(/"MINIMIZER"/),&
            description="Specify which method to use to perform a geometry optimization.",&
            usage="OPTIMIZER {BFGS|LBFGS|CG}",&
            enum_c_vals=s2a( "BFGS","LBFGS","CG"),&
            enum_i_vals=(/default_bfgs_method_id,default_lbfgs_method_id,default_cg_method_id/),&
            default_i_val=default_bfgs_method_id,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_ITER",&
            description="Specifies the maximum number of geometry optimization steps.",&
            usage="MAX_ITER {integer}",&
            default_i_val=200,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="MAX_H_RANK",&
            description="Maximum rank (and consequently size) of the "//&
            "approximate hessian matrix used by the LBFGS optimizer",&
            usage="MAX_H_RANK {integer}",&
            default_i_val=5,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_F_PER_ITER",&
            description="Maximum number of force evaluations per iteration"// &
            "(used for the line search)",&
            usage="MAX_F_PER_ITER {integer}",&
            default_i_val=3,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WANTED_PROJ_GRADIENT",&
            description="Requested norm thresold of the gradient multiplied "// &
            "by the approximate hessian",&
            usage="WANTED_PROJ_GRADIENT {real}",&
            default_r_val=1.0E-3_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WANTED_REL_F_ERROR",&
            description="Requested relative error on the objective function"//&
            "of the optimizer (the energy)",&
            usage="WANTED_REL_F_ERROR {real}",&
            default_r_val=1.0E-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_DR",&
            description="Tolerance on the maximum value of the displacement.",&
            usage="MAX_DR {real}",&
            default_r_val=0.0030_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_FORCE",&
            description="Tolerance on the maximum value of Forces.",&
            usage="MAX_FORCE {real}",&
            default_r_val=0.00045_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMS_DR",&
            description="Tolerance on RMS displacements.",&
            usage="RMS_DR {real}",&
            default_r_val=0.0015_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMS_FORCE",&
            description="Tolerance on RMS Forces.",&
            usage="RMS_FORCE {real}",&
            default_r_val=0.00030_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BRENT_TOL",&
            description="Tolerance requested during Brent line search in Conjugate Gradients Optimization.",&
            usage="BRENT_TOL {real}",&
            default_r_val=0.0001_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="BRACK_LIMIT",&
            description="Limit in 1D bracketing during line search in Conjugate Gradients Optimization.",&
            usage="BRACK_LIMIT {real}",&
            default_r_val=100.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="BRENT_MAX_ITER",&
            description="Maximum number of iterations in brent algorithm "// &
            "(used for the line search in Conjugated Gradients Optimization)",&
            usage="BRENT_MAX_ITER {integer}",&
            default_i_val=100,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="MAX_STEEP_STEPS",&
            description="Maximum number of steepest descent steps before starting the"//&
            " conjugate gradients optimization.",&
            usage="MAX_STEEP_STEPS {integer}",&
            default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="CG_RESTART_LIMIT",&
            description="Cosine of the angle between two consecutive searching directions."//&
            " If the angle during a CG optimization is less than the one corresponding to "//&
            " to the CG_RESTART_LIMIT the CG is reset and one step of steepest descent is "//&
            " performed.",&
            usage="CG_RESTART_LIMIT {real}",&
            default_r_val=0.1_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="USE_FIT_LINE_SEARCH",&
            description="Use instead of bracketing the minima and than applying the Brent"//&
            " algorithm, a least mean square to interpolate the minima in the line search.",&
            usage="USE_FIT_LINE_SEARCH",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART",&
            description="The guess starting HESSIAN is read from the RESTART.",&
            usage="RESTART",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FLETCHER_REEVES",&
            description="Uses FLETCHER-REEVES instead of POLAK-RIBIERE when using Conjugate Gradients",&
            usage="FLETCHER-REEVES",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_geoopt_section

!!****f* input_cp2k/create_metadyn_section *
!!
!!   NAME
!!     create_metadyn_section
!!
!!   FUNCTION
!!     creates the metadynamics section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_metadyn_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_metadyn_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="metadyn",&
            description="This section sets parameters to set up a calculation of metadynamics.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)
        
       NULLIFY(subsection,keyword)

       CALL keyword_create(keyword, name="NCOLVAR",&
            variants=(/"N_COLVAR"/),&
            description="Specifies the number of collective variables to use in the metadynamics."//&
            "Defaults 1 variable.",&
            usage="NCOLVAR {integer}",&
            default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NT_HILLS",&
            description="Specify the number the interval of time step to spawn the hills.",&
            usage="NT_HILLS {integer}",&
            default_i_val=30,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TEMPERATURE",&
            variants=(/"TEMPCONTROL"/),&
            description="If a Lagrangian scheme is used the temperature control for the collective "//&
            "variables is specified. Default 100.0 K.",&
            usage="TEMPERATURE {real}",&
            default_r_val=100.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WW",&
            description="Specifies the width of the gaussian to spawn. Default 0.1 .",&
            usage="WW {real}",&
            default_r_val=0.1_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DO_HILLS",&
            description="This keyword enables the spawning of the hills. Default .FALSE.",&
            usage="DO_HILLS",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LAGRANGE",&
            description="Specifies whether an extended-lagrangian should be used. Default .FALSE.",&
            usage="LAGRANGE",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="RESTART_HILLS",&
            variants=(/"RESTART"/),&
            description="Read hills from the RESTART file.",&
            usage="RESTART_HILLS",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_colvar_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_metadyn_section


!!****f* input_cp2k/create_colvar_section *
!!
!!   NAME
!!     create_colvar_section
!!
!!   FUNCTION
!!     creates the colvar section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="colvar",&
            description="This section specify the nature of the collective variables.",&
            n_keywords=1, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
            error=error)
        
       NULLIFY(subsection,keyword)

       CALL keyword_create(keyword, name="LAMBDA",&
            description="Specifies the lambda parameter of the collective variable in the"//&
            " extended lagrangian scheme.",&
            usage="LAMBDA {real}",&
            default_r_val=-HUGE(0.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MASS",&
            description="Specifies the mass parameter of the collective variable in the"//&
            " extended lagrangian scheme. Default 0.0 .",&
            usage="MASS {real}",&
            default_r_val=-HUGE(0.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SCALE",&
            description="Specifies the scale factor for the following collective variable.",&
            usage="SCALE {real}",&
            default_r_val=-HUGE(0.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WALL_PLUS",&
            variants=(/"WALLP"/),&
            description="Activates the reflective wall on the upper limit of the collective variable.",&
            usage="WALL_PLUS {real}",&
            default_r_val=-HUGE(0.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WALL_MINUS",&
            variants=(/"WALLM"/),&
            description="Activates the reflective wall on the lower limit of the collective variable.",&
            usage="WALL_MINUS {real}",&
            default_r_val=-HUGE(0.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SUBSYS",&
            description="Specifies the subsystem on which to apply the collective variables. Default 1.",&
            usage="SUBSYS {integer}",&
            default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
    
       CALL create_colvar_var_section(subsection, section, error)
       
    END IF

  END SUBROUTINE create_colvar_section


!!****f* input_cp2k/create_colvar_var_section *
!!
!!   NAME
!!     create_colvar_var_section
!!
!!   FUNCTION
!!     creates the collective variables for the colvar section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_var_section(subsection, section, error)
    TYPE(section_type), POINTER :: section, subsection
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_var_section',&
         routineP=moduleN//':'//routineN

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(subsection),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)

       CALL create_colvar_dist_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)      

       CALL create_colvar_angle_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)      

       CALL create_colvar_torsion_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)      

       CALL create_colvar_coord_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)      

    END IF
  END SUBROUTINE create_colvar_var_section

!!****f* input_cp2k/create_colvar_coord_section *
!!
!!   NAME
!!     create_colvar_coord_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables dist
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_coord_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_coord_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="coordination",&
            description="Section to define the coordination as a collective variables.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
        
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS_FROM",&
            description="Specify indexes of atoms building the coordination variable. ",&
            usage="ATOMS_FROM {integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_TO",&
            description="Specify indexes of atoms building the coordination variable. ",&
            usage="ATOMS_TO {integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_FROM",&
            description="Specify kinds of atoms building the coordination variable. ",&
            usage="KINDS_FROM {integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_TO",&
            description="Specify kinds of atoms building the coordination variable. ",&
            usage="KINDS_TO {integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="R0",&
            variants=(/"R_0"/),&
            description="Specify the R0 parameter in the coordination function.",&
            usage="R0 {real}",default_r_val=3.0_dp,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NN",&
            variants=(/"EXPON_NUMERATOR"/),&
            description="Sets the value of the numerator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="NN {integer}",default_i_val=6,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ND",&
            variants=(/"EXPON_DENOMINATOR"/),&
            description="Sets the value of the denominator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="ND {integer}",default_i_val=12,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_colvar_coord_section


!!****f* input_cp2k/create_colvar_torsion_section *
!!
!!   NAME
!!     create_colvar_torsion_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables dist
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_torsion_section(section,error,description)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: description
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_torsion_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       IF (PRESENT(description)) THEN
          CALL section_create(section,name="torsion",&
               description=description,&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       ELSE
          CALL section_create(section,name="torsion",&
               description="Section to define the torsion as a collective variables.",&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       END IF
        
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the atom indexes defining the torsion.",&
            usage="ATOMS {integer} {integer} {integer} {integer}",required=.TRUE.,&
            n_var=4, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_colvar_torsion_section

!!****f* input_cp2k/create_colvar_angle_section *
!!
!!   NAME
!!     create_colvar_angle_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables dist
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_angle_section(section,error,description)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: description
    
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_angle_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       IF (PRESENT(description)) THEN
          CALL section_create(section,name="angle",&
               description=description,&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       ELSE
          CALL section_create(section,name="angle",&
               description="Section to define the angle as a collective variables.",&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       END IF

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the atom indexes defining the angle.",&
            usage="ATOMS {integer} {integer} {integer}",required=.TRUE.,&
            n_var=3, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_colvar_angle_section
!****************************************************************************


!!****f* input_cp2k/create_colvar_dist_section *
!!
!!   NAME
!!     create_colvar_dist_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables dist
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_dist_section(section,error,description)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: description

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_dist_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       IF (PRESENT(description)) THEN
          CALL section_create(section,name="distance",&
               description=description,&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       ELSE
          CALL section_create(section,name="distance",&
               description="Section to define the distance as a collective variables.",&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       END IF

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the atom indexes defining the distance.",&
            usage="ATOMS {integer} {integer}",required=.TRUE.,&
            n_var=2, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_colvar_dist_section
!***************************************************************************

!!****f* input_cp2k/create_md_section *
!!
!!   NAME
!!     create_md_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the md section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE create_md_section(section,error)
    TYPE(section_type), pointer :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_md_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword
    TYPE(section_type), POINTER :: coeff_section, nose_section,uni_section, baro_section,&
                                   hmc_section, subsection, print_key

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="md",&
            description="parameters needed perform an md run",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, coeff_section, nose_section, uni_section, baro_section, subsection)
       CALL keyword_create(keyword, name="ensemble",&
            description="The ensemble that you want to sample with the md."//&
            " npt_i uses an isotropic resize of the cell, npt_f a flexible cell.",&
            usage="ensemble nve",&
            default_i_val=nve_ensemble,&
            enum_c_vals=ensemble_name_types,&
            enum_i_vals=ensemble_name_index,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="steps",&
            description="The number of steps to be performed",&
            usage="steps 100", default_i_val=3,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="timestep",&
            description="The length of an integration step.",&
            usage="timestep 1.0", default_r_val=0.5_dp,unit_str='fs',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="restart",&
            description="The way the trajectory should be done (init: from the "//&
            "initial (input file) positions, with velocities from a Boltzmann "//&
            "distribution; pos: using the latest positions stored in the restart"//&
            " file, with with velocities from a Boltzmann distribution"//&
            "; posvel: using both the latest positions and velocites stored in "//&
            "the restart file; all: restarts all.", usage="RESTART posvel", &
            default_i_val=md_init_start,&
            enum_c_vals=restart_name_types,&
            enum_i_vals=restart_name_index,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="temperature",&
            description="The temperature in K used to initialize the velocities"//&
            "with init and pos restart, and in the npt/nvt simulations", &
            usage="temperature 325.0",&
            default_r_val=300.0_dp, unit_str='K',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="temp_tol",&
            variants=(/"temp_to              ", "temperature_tolerance"/),&
            description="The maximum accepted oscillations of the temperature"//&
            "from the wanted temperature before rescaling the velocites."//&
            "If it is 0 no rescaling is performed", &
            usage="temp_tol 0.0", default_r_val=0.0_dp, unit_str='K', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="annealing",&
            description="Specifies the rescaling factor for velocities after the integration"//&
            "scheme used. Automatically enables the annealing procedure.", &
            usage="annealing 0.0", default_r_val=1.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="gamma",&
            description="Gamma parameter for the Langevin dynamics",&
            usage="gamma 0.001",&
            default_r_val=0.0_dp, unit_str='fs^-1',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="pv_ava",&
            description="Chris has to provide docs",&
            usage="pv_ava LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="pv_num",&
            description="Chris has to provide docs",&
            usage="pv_num LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_uniaxial_section(uni_section, error=error)
       CALL section_add_subsection(section, uni_section, error=error)
       CALL section_release(uni_section,error=error)       

       CALL create_barostat_section(baro_section, error=error)
       CALL section_add_subsection(section, baro_section, error=error)
       CALL section_release(baro_section,error=error)

       CALL create_coeff_section(coeff_section, error=error)
       CALL section_add_subsection(section, coeff_section, error=error)
       CALL section_release(coeff_section,error=error)

       CALL create_nose_section(nose_section, error=error)
       CALL section_add_subsection(section, nose_section, error=error)
       CALL section_release(nose_section,error=error)

       NULLIFY(hmc_section)
       CALL create_hmc_section(hmc_section, error=error)
       CALL section_add_subsection(section, hmc_section, error=error)
       CALL section_release(hmc_section, error=error)


       CALL create_qs_md_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)       


       CALL section_create(subsection,name="print",&
            description="Controls the printing properties during an MD run",&
            n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
            error=error)

       
       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"FILES",&
            description="Controls the printing of MD information on all files.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SCREEN",&
            description="Controls the printing of MD information on screen.", &
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"COORDINATE",&
            description="Controls the printing of coordinates during an MD run.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"VELOCITIES",&
            description="Controls the printing of velocities during an MD run.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"COEFFICIENTS",&
            description="Controls the printing of coefficients during an MD run.", &
            print_level=medium_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PRESSURE",&
            description="Controls the printing of the Pressure info during an MD run.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"ENERGIES",&
            description="Controls the printing of energies during an MD run.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"TEMPERATURE",&
            description="Controls the printing of temperature info during an MD run.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESTART",&
            description="Controls the dumping of the restart file during an MD run.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
       
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection, error=error)

    END IF
END SUBROUTINE create_md_section

!***************************************************************************
!!****f* input_cp2k/create_qs_md_section *
!!
!!   NAME
!!     create_qs_md_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the coeff section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************

SUBROUTINE create_qs_md_section(section, error)
  IMPLICIT NONE
  TYPE(section_type), POINTER :: section
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='create_qs_md_section',&
       routineP=moduleN//':'//routineN
  TYPE(keyword_type), POINTER :: keyword

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN  

     CALL section_create(section,name="QS_MD",&
          description="QS specific MD parameters.",&
          n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
       
     NULLIFY(keyword)

     CALL keyword_create(keyword, name="neighbor_lists_from_scratch",&
          variants=(/"NEIGH"/),&
          description="This keyword enables the building of the neighbouring list from scratch."//&
          " Default .TRUE.",&
          usage="NEIGH logical",&
          default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     
  END IF
END SUBROUTINE create_qs_md_section

!***************************************************************************
!!****f* input_cp2k/create_uniaxial_section *
!!
!!   NAME
!!     create_uniaxial_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the coeff section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************

SUBROUTINE create_uniaxial_section(section, error)
  IMPLICIT NONE
  TYPE(section_type), POINTER :: section
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='create_uniaxial_section',&
       routineP=moduleN//':'//routineN
  TYPE(keyword_type), POINTER :: keyword
  TYPE(section_type), POINTER :: subsection

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN  

       CALL section_create(section,name="uniaxial",&
            description="Parameters for uniaxial shock.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword )
       CALL keyword_create(keyword, name="PRESSURE",&
            description="Initial pressure",&
            usage="PRESSURE real",&
            default_r_val=0._dp, n_var=1, unit_str='bar', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)            

       CALL keyword_create(keyword, name="ENERGY",&
            description="Initial energy",&
            usage="ENERGY real",&
            default_r_val=0._dp, n_var=1, unit_str='hartree',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)         

       CALL keyword_create(keyword, name="VOLUME",&
            description="Initial compression",&
            usage="VOLUME real",&
            default_r_val=0._dp, n_var=1, unit_str='angstrom^3',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)             

       CALL keyword_create(keyword, name="CMASS",&
            description="Barostat mass",&
            usage="CMASS real",&
            default_r_val=0._dp, n_var=1,unit_str='au_m',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)             

       CALL keyword_create(keyword, name="VSHOCK",variants=(/"V_SHOCK"/),&
            description="Velocity of shock front",&
            usage="VSHOCK real",&
            default_r_val=0._dp, n_var=1,unit_str='m/s',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)             

       CALL keyword_create(keyword, name="GAMMA",&
            description="Damping coefficient for barostat and particles",&
            usage="GAMMA real",&
            unit_str='fs^-1', &
            default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)             

    END IF
  END SUBROUTINE create_uniaxial_section

!***************************************************************************
!!****f* input_cp2k/create_barostat_section *
!!
!!   NAME
!!     create_barostat_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the coeff section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************

SUBROUTINE create_barostat_section(section, error)
  IMPLICIT NONE
  TYPE(section_type), POINTER :: section
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='create_barostat_section',&
       routineP=moduleN//':'//routineN
  TYPE(keyword_type), POINTER :: keyword

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN  

       CALL section_create(section,name="barostat",&
            description="Parameters of barostat.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="PRESSURE",&
            description="Initial pressure",&
            usage="PRESSURE real",&
            default_r_val=0._dp,unit_str='bar',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)      

       CALL keyword_create(keyword, name="TIMECON",&
            description="Barostat time constant",&
            usage="TIMECON real",&
            default_r_val=cp_unit_to_cp2k(1000._dp,"fs",error=error),&
            unit_str='fs',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)        

    END IF
  END SUBROUTINE create_barostat_section

!***************************************************************************
!!****f* input_cp2k/create_coeff_section *
!!
!!   NAME
!!     create_coeff_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the coeff section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************

SUBROUTINE create_coeff_section(coeff_section, error)
  IMPLICIT NONE
  TYPE(section_type), POINTER :: coeff_section
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='create_coeff_section',&
       routineP=moduleN//':'//routineN
  TYPE(keyword_type), POINTER :: keyword
  TYPE(section_type), POINTER :: nose_section

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(coeff_section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN  

       CALL section_create(coeff_section,name="kg_coeff",&
            description="coefficent parameters needed to perform a kg-pol md run",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       
       NULLIFY(nose_section, keyword)
       CALL keyword_create(keyword, name="ensemble",&
            description="md ensemble for KG-POL coefficients", usage="ENSEMBLE NVT", &
            default_i_val=nvt_ensemble,&
            enum_c_vals=(/"nve  ","nvt  "/),&
            enum_i_vals=(/nve_ensemble,nvt_ensemble/),&
            error=error)
       CALL section_add_keyword(coeff_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       call keyword_create(keyword, name="nhcopt",&
            description="global or massive thermostats for KG-POL coefficients", &
            usage="NHCOPT GLOBAL", &
            default_i_val=nhc_global,&
            enum_c_vals=s2a("global","massive"),&
            enum_desc=s2a("just one thermostat for all the coefficents",&
            "each coefficient has its own thermostat"),&
            enum_i_vals=(/nhc_global,nhc_massive/),&
            error=error)
       CALL section_add_keyword(coeff_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       call keyword_create(keyword, name="storage",&
            description="method to store the effective masses for"//&
            "KG-POL coefficients", usage="storage global", &
            default_i_val=dyn_coeff_mass_global,&
            enum_c_vals=s2a("global","kind","atom","coeff"),&
            enum_desc=s2a("just one mass for all the coefficents",&
            "each atomic kind has its own effective mass",&
            "coefficients of each atom have their own effective mass",&
            "each coefficient has its own effective mass"),&
            enum_i_vals=(/dyn_coeff_mass_global,dyn_coeff_mass_kind,&
                          dyn_coeff_mass_atom,dyn_coeff_mass_coeff/),&
            error=error)
       CALL section_add_keyword(coeff_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="temperature",&
            description="target temperature for the KG-POL coeffs",&
            default_r_val=1.0_dp, usage="temperature 1.0",&
            unit_str='K',error=error)
       CALL section_add_keyword(coeff_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="mass",&
            description="mass of the KG-POL coeffs if storage is global",&
            default_r_val=1.0_dp, usage="mass 1.0",unit_str='au_m',&
            error=error)
       CALL section_add_keyword(coeff_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="restart",&
            description="The way the coefficents should initialized.", usage="RESTART posvel", &
            default_i_val=kg_coeff_init_start,&
            enum_c_vals=s2a("init","pos","posvel"),&
            enum_desc=s2a("from scratch",&
            "position from restart, velocities from Boltzmann distribution",&
            "positions and velocities from restart"),&
            enum_i_vals=(/kg_coeff_init_start,kg_coeff_pos_restart, &
            kg_coeff_pos_vel_restart/),&
            error=error)
       CALL section_add_keyword(coeff_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_nose_section(nose_section, error=error)
       CALL section_add_subsection(coeff_section, nose_section, error=error)
       CALL section_release(nose_section,error=error)
       
    END IF
  END SUBROUTINE create_coeff_section
!***************************************************************************
!!****f* input_cp2k/create_nose_section *
!!
!!   NAME
!!     create_nose_section
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
!!     gloria
!!
!!*** **********************************************************************

SUBROUTINE create_nose_section(section, error)
  IMPLICIT NONE
  TYPE(section_type), POINTER :: section
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='create_nose_section',&
       routineP=moduleN//':'//routineN
  TYPE(keyword_type), POINTER :: keyword

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN  
       CALL section_create(section,name="nose",&
            description="paramameters of the Nose Hoover thermostat chain for KG coeffs",&
            required=.FALSE.,error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="length",&
            description="length of the NoseHoover chain", usage="length integer", &
            default_i_val=3,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="Yoshida",&
            description="order of the yoshida integretor usef for the thermostat",&
            usage="Yoshida integer", &
            default_i_val=3,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="timecon",&
            description="timeconstant of the thermostat chain",&
            usage="timeconstant real", &
            default_r_val=cp_unit_to_cp2k(1000._dp,'fs',error=error),&
            unit_str='fs',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="mts", &
            variants=s2a("multiple_time_steps","mult_t_steps"),&
            description="number of multiple timesteps to be used for the NoseHoover chain",&
            usage="mts integer", &
            default_i_val=2,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

  END IF

END SUBROUTINE create_nose_section

!***************************************************************************

!!****f* input_cp2k/create_hmc_section [1.0] *
!!
!!   NAME
!!     create_hmc_section
!!
!!   FUNCTION
!!     This section specifies the parameters for the Hybrid Monte Carlo simulation
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
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!*** **********************************************************************
  SUBROUTINE create_hmc_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_hmc_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="HMC",&
            description="This section specifies the parameters for the "//&
                        "Hybrid Metropolis Monte Carlo simulation",&
            n_keywords=14, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       !CALL keyword_create(keyword, name="gamma",&
       !     description="Gamma parameter for the Langevin dynamics",&
       !     usage="gamma 0.001",&
       !     default_r_val=0.0_dp, unit_str='fs^-1',error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="ld_steps", &
            description="Number of LD steps per MC cycle.", &
            usage="nsteps 20", required=.TRUE.,&
            default_i_val=10, lone_keyword_i_val=20,&
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="temp_alpha", &
            description="Additive temperature correction", &
            usage="temp_alpha -10.0", required=.FALSE., &
            default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="temp_beta", &
            description="Mutiplicative temperature correction", &
            usage="temp_beta 0.75", required=.FALSE., &
            default_r_val=1.0_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="tot_energy", &
            description="Indicates if the total energy differences "//&
            "or the potential energy differences are used", &
            usage="tot_energy TRUE", required=.FALSE., &
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="ext_temp", &
            description="Indicates if the external or the computed "//&
            "current temperature is used", &
            usage="ext_temp TRUE", required=.FALSE., &
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="qs_accept", &
            description="Denotes if the QS or the in the force "//&
            "environment specified energies are used", &
            usage="qs_accept TRUE", required=.FALSE., &
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="rnd_velocities", &
            description="Activates random velocity initialization", &
            usage="rnd_velocities TRUE", required=.FALSE., &
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="vel_scaling", &
            description="The velocities are scaled, such that the current"//&
            "temperature is equivalent with the external temperature", &
            usage="vel scaling TRUE", required=.FALSE., &
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="qs_start", &
            description="Denotes if a fully converged QS energy "//&
            "evaluation is performed at the start-up or not", &
            usage="qs_start TRUE", required=.FALSE., &
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="restore_history", &
            description="If activated the WF history is fully restored "//&
            "in case a MD approach is rejected", &
            usage="restore_history TRUE", required=.FALSE., &
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="recover_momenta", &
            description="Causes that the momentas are recovered in "//&
            "CASE a MD approach is rejected", &
            usage="recover_momenta TRUE", required=.FALSE., &
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="semi_hybrid", &
            description="Results that the momenta are choosen at "//&
            "random, only in case a MD approach is rejected.", &
            usage="semi_hybrid TRUE", required=.FALSE., &
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)
       
       CALL keyword_create(keyword, name="MD_method", &
            description="Indicates which method is used for the MD approach", &
            enum_c_vals=s2a("FIST", "KG", "QS", "HARRIS") , &
            enum_i_vals=(/FIST_MD_method, KG_MD_method, QS_MD_method, Harris_MD_method/),&
            usage="MD_method KG", required=.FALSE., &
            default_i_val=KG_MD_method, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="MC_method", &
            description="Indicates according to which method the MC criterion is evaluated", &
            enum_c_vals=s2a("KG", "QS", "HARRIS") , &
            enum_i_vals=(/KG_MC_method, QS_MC_method, Harris_MC_method/),&
            usage="MC_method Harris", required=.FALSE., &
            default_i_val=Harris_MC_method, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

    END IF
  END SUBROUTINE create_hmc_section

!***************************************************************************

END MODULE input_cp2k_motion
