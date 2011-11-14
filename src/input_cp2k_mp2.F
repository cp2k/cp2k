!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief input section for MP2
!> \par History
!>      05.2011 created
!> \author MDB
! *****************************************************************************
MODULE input_cp2k_mp2
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

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_mp2'

  PUBLIC :: create_mp2_section

CONTAINS

! *****************************************************************************
!> \brief creates the input section for the mp2 part
!> \param section the section to create
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author MDB
! *****************************************************************************
  SUBROUTINE create_mp2_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mp2_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"MP2",&
            description="Sets up the MP2 parameters if requested ",&
            n_keywords=6, n_subsections=3, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, print_key, subsection)

       CALL keyword_create(&
            keyword=keyword,&
            name="METHOD",&
            description="Which method should be used to compute the MP2 energy",&
            usage="METHOD LAPLACE",&
            enum_c_vals=s2a("NONE","DIRECT_CANONICAL","MP2_GPW"),&
            enum_i_vals=(/mp2_method_none,mp2_method_direct,mp2_method_gpw/),&
            enum_desc=s2a("Skip MP2 calculation",&
                          "Use the direct mp2 canonical approach",&
                          "Use a GPW approach"),&
            default_i_val=mp2_method_direct, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="MEMORY",&
            description="Memory usage in direct mp2 method",&
            usage="MEMORY 1024.0",&
            default_r_val=1.024E+3_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="GROUP_SIZE",&
            variants=(/"NUMBER_PROC"/),&
            description="Group size used in the computation of the integrals. Default is to use all processors (GROUP_SIZE=-1)."//&
             "A smaller group size (for example the node size), might a better choice if the actual MP2 time is large "//&
             "compared to integral computation time. This is usually the case if the total number of processors is not too large.",&
            usage="GROUP_SIZE 32",&
            default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"MP2_INFO",&
            description="Controls the printing basic info about MP2 method", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
       
       CALL create_mp2_direct(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_mp2_gpw(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_mp2_potential(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
 
  END SUBROUTINE create_mp2_section

  SUBROUTINE create_mp2_direct(section,error)                                 
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mp2_direct', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"DIRECT_CANONICAL",&
            description="Parameters influencing the direct canonical method",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(&
            keyword=keyword,&
            name="BIG_SEND",&
            description="Send big messages between processes (useful for >48 processors).",&
            usage="BIG_SEND",&
            default_l_val=.TRUE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_mp2_direct

  SUBROUTINE create_mp2_GPW(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mp2_GPW', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"MP2_GPW",&
            description="Parameters for the GPW approach",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="EPS_GRID",&
            description="Determines a threshold for the GPW based integration",&
            usage="EPS_GRID 1.0E-5 ",required=.FALSE.,type_of_var=real_t,&
            default_r_val=1.0E-6_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_FILTER",&
            description="Determines a threshold for the DBCSR based multiply (usually 10 times smaller than EPS_GRID).",&
            usage="EPS_FILTER 1.0E-6 ",required=.FALSE.,type_of_var=real_t,&
            default_r_val=1.0E-7_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_mp2_GPW


  SUBROUTINE create_mp2_potential(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mp2_potential', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"INTERACTION_POTENTIAL",&
            description="Parameters the interaction potential in computing the biel integrals",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(&
            keyword=keyword,&
            name="POTENTIAL_TYPE",&
            description="Which interaction potential should be used "//&
                         "(Coulomb, TShPSC operator).",&
            usage="POTENTIAL_TYPE TSHPSC",&
            enum_c_vals=s2a("COULOMB","TShPSC"),&
            enum_i_vals=(/do_hfx_potential_coulomb,& 
                          do_mp2_potential_TShPSC/),&
            enum_desc=s2a("Coulomb potential: 1/r", &
                          "TShPSC: 1/x - s/Rc  ->  for x<=Rc "//&
                          "(1-s)/Rc - "//&
                          "(x-Rc)/Rc^2 + "//&
                          "(x-Rc)^2/Rc^3 - "//&
                          "(2*n^2-7*n+9-4*s)*(x-Rc)^3/(Rc^4*(n^2-2*n+1)*(n-1)) + "//&
                          "(6-3*s-4*n+n^2)*(x-Rc)^4/(Rc^5*(n^4-4*n^3+6*n^2-4*n+1))   -> for Rc<x<=n*Rc (4th order polynomial) "//&
                          " 0  -> for x>n*Rc"),&
            default_i_val=do_hfx_potential_coulomb, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
      
       CALL keyword_create(keyword, name="TRUNCATION_RADIUS",&
            description="Determines truncation radius for the truncated TShPSC potential. "//&
                        "Only valid when doing truncated calculation",&
            usage="TRUNCATION_RADIUS 10.0",required=.FALSE.,type_of_var=real_t,&  
            default_r_val=10.0_dp,&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="TShPSC_DATA",&
            description="Location of the file TShPSC.dat that contains the data for the "//&
                        "evaluation of the TShPSC G0 ",&
            usage="TShPSC_DATA /data/t_sh_p_s_c.dat",&
            default_c_val="../../t_sh_p_s_c.dat",&  ! this is tests/QS i.e. regtests are possible
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_mp2_potential
   
END MODULE input_cp2k_mp2

