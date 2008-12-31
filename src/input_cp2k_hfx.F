!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief function that builds the hartree fock exchange section of the input
!> \par History
!>      09.2007 created
!> \author Manuel Guidon
! *****************************************************************************
MODULE input_cp2k_hfx
  USE bibliography,                    ONLY: Guidon2008
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE f77_blas
  USE input_constants
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_hfx'

  PUBLIC :: create_hfx_section

CONTAINS

! *****************************************************************************
!> \brief creates the input section for the hf part
!> \param section the section to create
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE create_hfx_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_hfx_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"HF",&
            description="Sets up the Hartree-Fock parameters if requested ",&
            n_keywords=5, n_subsections=2, repeats=.TRUE., required=.TRUE.,&
            citations=(/Guidon2008/),&
            error=error)

       NULLIFY(keyword, print_key, subsection)

       CALL keyword_create(keyword, name="FRACTION",&
            description="The fraction of Hartree-Fock to add to the total energy. "//&
                        "1.0 implies standard Hartree-Fock if used with XC_FUNCTIONAL NONE. "//&
                        "NOTE: In a mixed potential calculation this should be set to 1.0, otherwise "//&
                        "all parts are multiplied with this factor. ",&
            usage="FRACTION 1.0", default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PW_HFX",&
            description="Compute the Hartree-Fock energy also in the plane wave basis."//&
                        "The value is ignored, and intended for debugging only.",&
            usage="PW_HFX FALSE", default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PW_HFX_BLOCKSIZE",&
            description="Improve the performance of pw_hfx at the cost of some additional memory "//&
                        "by storing the realspace representation of PW_HFX_BLOCKSIZE states.",&
            usage="PW_HFX_BLOCKSIZE 20", default_i_val=20,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"HF_INFO",&
            description="Controls the printing basic info about hf method", &
            print_level=medium_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
       
       CALL create_hf_pbc_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_hf_screening_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_hf_potential_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_hf_load_balance_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
      
       CALL create_hf_memory_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
 
  END SUBROUTINE create_hfx_section

! *****************************************************************************
!> \brief !****f* input_cp2k_dft/create_hf_load_balance_section [1.0] *
!> 
!>      creates the input section for the hf potential part
!> \param section the section to create
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE create_hf_load_balance_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_hf_load_balance_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"LOAD_BALANCE",&
            description="Sets up interaction potential if requested ",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            citations=(/guidon2008/),&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(&
            keyword=keyword,&
            name="NBINS",&
            description="Number of bins to be used on each process to distribute "//&
                        "in a first step.",&
            usage="NBINS 32",&
            default_i_val=16,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOVES",&
            description="Total number of Monte Carlo moves to use in the optimization",&
            usage="MOVES 1000000 ",&
            n_var=1,default_i_val=100000,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="CYCLES",&
            description="How many annealing cycles to use",&
            usage="CYCLES 10", default_i_val=10, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="REDUCTION_STEPS",&
            description="How many temperature reduction steps have to be performed per cycle",&
            usage="REDUCTION_STEPS 5", default_i_val=5, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="REDUCTION_FACTOR",&
            description="How to scale the temperature at each reduction step",&
            usage="REDUCTION_FACTOR 0.8", default_r_val=0.8_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="TEMPERATURE",&
            description="Highest temperature to be used (in units of the maximum cost of a non-zero 2D block)",&
            usage="TEMPERATURE 0.5", default_r_val=0.5_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="TOLERANCE",&
            description="If the found distribution is within a fraction of the theoretical optimal one, the optimization stops",&
            usage="TOLERANCE 0.01", default_r_val=0.02_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="SWAP_PROBABILITY",&
            description="Probability of a col/row assignment swap",&
            usage="SWAP_PROBABILITY 0.8", default_r_val=0.9_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of basic information during the run", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section, print_key, error=error)
       CALL section_release(print_key,error=error)

    END IF
 
  END SUBROUTINE create_hf_load_balance_section

! *****************************************************************************
!> \brief !****f* input_cp2k_dft/create_hf_potential_section [1.0] *
!> 
!>      creates the input section for the hf potential part
!> \param section the section to create
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE create_hf_potential_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_hf_potential_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"INTERACTION_POTENTIAL",&
            description="Sets up interaction potential if requested ",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            citations=(/guidon2008/),&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(&
            keyword=keyword,&
            name="POTENTIAL_TYPE",&
            description="Which interaction potential should be used "//&
                         "(Coulomb, longrange or shortrange).",&
            usage="POTENTIAL_TYPE SHORTRANGE",&
            enum_c_vals=s2a("COULOMB","SHORTRANGE","LONGRANGE","MIX_CL","GAUSSIAN","MIX_LG","IDENTITY","TRUNCATED"),&
            enum_i_vals=(/do_hfx_potential_coulomb, do_hfx_potential_short, do_hfx_potential_long, &
                          do_hfx_potential_mix_cl, do_hfx_potential_gaussian, do_hfx_potential_mix_lg,&
                          do_hfx_potential_id, do_hfx_potential_truncated/),&
            enum_desc=s2a("Coulomb potential: 1/r", &
                          "Shortrange potential: erfc(omega*r)/r",&
                          "Longrange potential: erf(omega*r)/r",&
                          "Mix coulomb and longrange potential: 1/r + erf(omega*r)/r",&
                          "Damped Gaussian potential: exp(-omega^2*r^2)",&
                          "Mix Gaussian and longrange potential: erf(omega*r)/r + exp(-omega^2*r^2)",&
                          "Overlap",&
                          "Truncated coulomb potential: if(r<R_c) 1/r else 0"),&
            default_i_val=do_hfx_potential_coulomb, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(keyword)
       CALL keyword_create(&
            keyword=keyword,&
            name="OMEGA",&
            description="Parameter for short/longrange interaction",&
            usage="OMEGA 0.5",&
            default_r_val=0.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SCALE_COULOMB",&
            description="Scales Hartree-Fock contribution arising from a coulomb potential. "//&
                        "Only valid when doing a mixed potential calculation",&
            usage="SCALE_COULOMB 1.0", default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SCALE_LONGRANGE",&
            description="Scales Hartree-Fock contribution arising from a longrange potential. "//&
                        "Only valid when doing a mixed potential calculation",&
            usage="SCALE_LONGRANGE 1.0", default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SCALE_GAUSSIAN",&
            description="Scales Hartree-Fock contribution arising from a gaussian potential. "//&
                        "Only valid when doing a mixed potential calculation",&
            usage="SCALE_GAUSSIAN 1.0", default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CUTOFF_RADIUS",&
            description="Determines cutoff radius for the truncated 1/r potential. "//&
                        "Only valid when doing truncated calculation",&
            usage="CUTOFF_RADIUS 10.0",required=.TRUE.,type_of_var=real_t,&! default_r_val=10.0_dp,&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="T_C_G_DATA",&
            description="Location of the file t_c_g.dat that contains the data for the "//&
                        "evaluation of the truncated gamma function ",&
            usage="T_C_G_DATA /data/t_c_g.dat",&
            default_c_val="../../t_c_g.dat",&  ! this is tests/QS i.e. regtests are possible
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
 
  END SUBROUTINE create_hf_potential_section

!****f* input_cp2k_dft/create_hf_screening_section [1.0] *

! *****************************************************************************
!> \brief creates the input section for the hf screening part
!> \param section the section to create
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE create_hf_screening_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_hf_screening_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"SCREENING",&
            description="Sets up screening parameters if requested ",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/guidon2008/),&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(&
            keyword=keyword,&
            name="EPS_SCHWARZ",&
            description="Screens the near field part of the electronic repulsion "//&
                        "integrals using the Schwarz inequality for the given "//&
                        "threshold.",&  
            usage="EPS_SCHWARZ 1.0E-6",&
            default_r_val=1.0E-10_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       NULLIFY(keyword)
       CALL keyword_create(&
            keyword=keyword,&
            name="EPS_BOX",&
            description="Screens the far field part of the electronic repulsion "//&
                        "integrals using a box screening for the given "//&
                        "threshold.",&
            usage="EPS_BOX 1.0E-6",&
            default_r_val=1.0E-10_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(keyword)
       CALL keyword_create(&
            keyword=keyword,&
            name="EPS_SCHWARZ_FORCES",&
            description="Screens the near field part of the electronic repulsion "//&
                        "integrals using the Schwarz inequality for the given "//&
                        "threshold. This will be approximately the accuracy of the forces, "//&
                        " and should normally be similar to EPS_SCF",&  
            usage="EPS_SCHWARZ_FORCES 1.0E-5",&
            default_r_val=1.0E-6_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       NULLIFY(keyword)
       CALL keyword_create(&
            keyword=keyword,&
            name="EPS_BOX_FORCES",&
            description="Screens the far field part of the electronic repulsion "//&
                        "integrals using a box screening for the given "//&
                        "threshold. This will be approximately the accuracy of the forces, "//&
                        "and should normally be similar to EPS_SCF",&
            usage="EPS_BOX_FORCES 1.0E-5",&
            default_r_val=1.0E-6_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(keyword)
       CALL keyword_create(&
            keyword=keyword,&
            name="SCREEN_P_FORCES",&
            description="Screens the electronic repulsion integrals for the forces "//&
                        "using the density matrix. "//&
                        "This results in a significant speedup for large systems, "//&
                        "but might require a somewhat tigher EPS_SCHWARZ_FORCES.", &
            usage="SCREEN_P_FORCES TRUE",&
            default_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="SCREEN_ON_INITIAL_P",&
            description="Screen on an initial density matrix. For the first MD step"//&
                        " this matrix must be provided by a Restart File.",&
            usage="SCREEN_ON_INITIAL_P  TRUE", default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="P_SCREEN_CORRECTION_FACTOR",&
            description="Recalculates integrals on the fly if the actual density matrix is"//&
                        " larger by a given factor than the initial one. If the factor is set"//&
                        " to 0.0_dp, this feature is disbled.",&
            usage="P_SCREEN_CORRECTION_FACTOR  0.0_dp", default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
 
  END SUBROUTINE create_hf_screening_section

! *****************************************************************************
!> \brief creates the input section for the hf-pbc part
!> \param section the section to create
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE create_hf_pbc_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_hf_pbc_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"PERIODIC",&
            description="Sets up periodic boundary condition parameters if requested ",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            citations=(/guidon2008/),&
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(&
            keyword=keyword,&
            name="NUMBER_OF_SHELLS",&
            description="Number of shells taken into account for periodicity. "//&
            "By default, cp2k tries to automatically evaluate this number. "//&
            "This algorithm might be to conservative, resulting in some overhead. "//&
            "You can try to adjust this number in order to make a calculation cheaper. ", &
            usage="NUMBER_OF_SHELLS 2",&
            default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
 
  END SUBROUTINE create_hf_pbc_section

! *****************************************************************************
!> \brief creates the input section for the hf-memory part
!> \param section the section to create
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE create_hf_memory_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_hf_memory_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"MEMORY",&
            description="Sets up memory parameters for the storage of the ERI's if requested ",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            citations=(/guidon2008/),&
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(&
            keyword=keyword,&
            name="EPS_STORAGE_SCALING",&
            variants=(/"EPS_STORAGE"/),&
            description="Scaling factor to scale eps_schwarz. Storage threshold for compression "//&
                        "will be EPS_SCHWARZ*EPS_STORAGE_SCALING.",&
            usage="EPS_STORAGE 1.0E-2",&
            default_r_val=1.0E0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="MAX_MEMORY",&
            description="Defines the maximum amount of memory [MB] used to store precomputed "//&
                         "compressed four-center integrals",&
            usage="MAX_MEMORY 256",&
            default_i_val=512,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(&
            keyword=keyword,&
            name="STORAGE_LOCATION",&
            description="Loaction where ERI's are stored if MAX_DISK_SPACE /=0 "//&
                        "Expects a path to a directory. ",&
            usage="STORAGE_LOCATION /data/scratch",&
            default_c_val=".",&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="MAX_DISK_SPACE",&
            description="Defines the maximum amount of disk space [MB] used to store precomputed "//&
                         "compressed four-center integrals. If 0, nothing is stored to disk",&
            usage="MAX_DISK_SPACE 256",&
            default_i_val=0,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

 
    END IF
 
  END SUBROUTINE create_hf_memory_section

END MODULE input_cp2k_hfx
