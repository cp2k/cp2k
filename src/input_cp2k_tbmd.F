!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_tbmd *
!!
!!   NAME
!!     input_cp2k_tbmd
!!
!!   FUNCTION
!!     Creates the TBMD section of the input
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 created
!!
!!   SOURCE
!****************************************************************************
MODULE input_cp2k_tbmd
  USE input_constants
  USE input_keyword_types, ONLY: keyword_type, keyword_create,&
       keyword_release
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release
  USE input_val_types,     ONLY: char_t, real_t, enum_t, logical_t
  USE kinds,               ONLY: dp
  USE string_utilities,    ONLY: s2a
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='input_cp2k_tbmd'

  PUBLIC :: create_tbmd_section, create_tb_model_section
!!***
!****************************************************************************
CONTAINS

!!****f* input_cp2k_tbmd/create_tbmd_section [1.0] *
!!
!!   NAME
!!     create_tbmd_section
!!
!!   FUNCTION
!!     Create the input section for TBMD
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
!!   MODIFICATION HISTORY
!!     01.2006 created
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE create_tbmd_section(section, error)
    TYPE(section_type), POINTER                  :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_tbmd_section', &
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    TYPE(keyword_type), POINTER :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="TBMD", &
            description="This section contains all information to run a TBMD calculation.", &
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE., &
            error=error)

       NULLIFY(subsection)

       CALL create_tb_model_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection, error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="Noisy-TB", &
            description="Turns on Noisy-TightBinding", &
            usage="Noisy-TB TRUE", required=.TRUE., type_of_var=logical_t, &
            n_var=1, repeats=.FALSE., variants=(/"Noisy_TB"/), &
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

    END IF

  END SUBROUTINE create_tbmd_section
!***************************************************************************

!!****f* input_cp2k_tbmd/create_tb_model_section [1.0] *
!!
!!   NAME
!!     create_tb_model_section
!!
!!   FUNCTION
!!     This section specifies the Tight-Binding parametrization model
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
!!   MODIFICATION HISTORY
!!     01.2006 created
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE create_tb_model_section(section, error)
    TYPE(section_type), POINTER                  :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_tb_model_section', &
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="TB_Model", &
            description="This section specifies the Tight-Binding parametrization model.", &
            n_keywords=18, n_subsections=0, repeats=.FALSE., required=.TRUE., &
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="TB-Model", &
            description="Defines the TB parametrization model.", &
            usage="TB-Model GSP", required=.TRUE., type_of_var=enum_t, &
            n_var=1, repeats=.FALSE., variants=(/"TB_Model"/), &
            enum_c_vals=s2a("GSP", "Sawada", "Kwon", "Lenosky", "UDP"), &
            enum_i_vals=(/use_gsp, use_sawada, use_kwon, use_lenosky, use_udp/), &
            default_i_val=use_gsp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="r_cutoff", &
            description="Cutoff length: More distant interactions are neglected.", &
            usage="r_cutoff 3.75", required=.TRUE., type_of_var=real_t, &
            n_var=1, repeats=.FALSE., default_r_val=3.75_dp, &
            unit_str="angstrom", error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       ! E_tot = Sum_n[2*<Psi_n|H_TB|Psi_n> + E_rep + E_0*N_at]
       ! My general form for hopping parameters and pair potentials: 
       ! h_alpha(r) = h_alpha(r_0)(r_0/r)**n * exp{n*[-(r/r_calpha)**n_calpha + (r_0/r_calpha)**n_calpha]}
       ! E_rep = Sum_i[f*Sum_j[Phi(r_ij)]], with f(x)=C_1*x + C_2*x**2 + C_3*x**3 + C_4*x**4 and
       ! Phi(r) = Phi(1)*(r_0/r)**m * exp{m*[-(r/d_c)**m_c + (r_0/d_c)**m_c]}

       CALL keyword_create(keyword, name="h_ssS", &
            description="The value of the hopping coefficient h(1)_ssS", &
            usage="h_ssS -1.82", required=.FALSE., type_of_var=real_t, &
            n_var=1, repeats=.FALSE., variants=(/"h_ssSigma"/), &
            default_r_val=-1.82_dp, unit_str="eV", error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="h_spS", &
            description="The value of the hopping coefficient h(1)_spS", &
            usage="h_spS 1.96", required=.FALSE., type_of_var=real_t, &
            n_var=1, repeats=.FALSE., variants=(/"h_spSigma"/), &
            default_r_val=+1.96_dp, unit_str="eV", error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="h_ppS", &
            description="The value of the hopping coefficient h(1)_ppS", &
            usage="h_ppS 3.06", required=.FALSE., type_of_var=real_t, &
            n_var=1, repeats=.FALSE., variants=(/"h_ppSigma"/), &
            default_r_val=+3.06_dp, unit_str="eV", error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="h_ppP", &
            description="The value of the hopping coefficient h(1)_ppP", &
            usage="h_ppP -0.87", required=.FALSE., type_of_var=real_t, &
            n_var=1, repeats=.FALSE., variants=(/"h_ppPi"/), &
            default_r_val=-0.87_dp, unit_str="eV", error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="n", &
            description="Hopping exponent", &
            usage="n 2", required=.FALSE., type_of_var=real_t, &
            n_var=1, repeats=.FALSE., default_r_val=2.0_dp, &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="r_0", &
            description="Nearest neighbor distance", &
            usage="r_0 2.360352", required=.FALSE., type_of_var=real_t, &
            n_var=1, repeats=.FALSE., default_r_val=2.35126_dp, &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="r_calpha", &
            description="Individual cutoffs for each hopping parameter", &
            usage="r_calpha <Real r_c1> <Real r_c2> <Real r_c3> <Real r_c4>", &
            required=.FALSE., type_of_var=real_t, n_var=4, repeats=.FALSE., &
            default_r_vals=(/3.67_dp, 3.67_dp, 3.67_dp, 3.67_dp/), &
            unit_str="angstrom", error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="n_calpha", &
            description="Individual exponents for each hopping parameter", &
            usage="n_calpha <Real n_c1> <Real n_c2> <Real n_c3> <Real n_c4>", &
            required=.FALSE., type_of_var=real_t, n_var=4, repeats=.FALSE., &
            default_r_vals=(/6.48_dp, 6.48_dp, 6.48_dp, 6.48_dp/), error=error)            
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)


       CALL keyword_create(keyword, name="Phi(1)", &
            description="Phi(1) pair-potential parameter", &
            usage="Phi(1) 3.4581", required=.FALSE., type_of_var=real_t, &
            n_var=1, repeats=.FALSE., default_r_val=3.4581_dp, &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="m", &
            description="The value of the pair-potential exponent", &
            usage="m 4.54", required=.FALSE., type_of_var=real_t, &
            n_var=1, repeats=.FALSE., default_r_val=4.54_dp, &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="d_c", &
            description="Pair-potential cutoff", &
            usage="d_c 3.67", required=.FALSE., type_of_var=real_t, &
            n_var=1, repeats=.FALSE., default_r_val=3.67_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="m_c", &
            description="Pair-potential exponent m_c", &
            usage="m_c", required=.FALSE., type_of_var=real_t, n_var=1, &
            repeats=.FALSE., default_r_val=6.48_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="C14", &
            description="Parameters C1-C4 of Kwon's 4th order polynomial for the pair-repulsion", &
            usage="C14 <Real C1> <Real C2> <Real C3> <Real C4>", &
            required=.FALSE., type_of_var=real_t, n_var=4, repeats=.FALSE., &
            default_r_vals=(/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/), error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)


       CALL keyword_create(keyword, name="dEpsilon", &
            description="s-p energy splitting dEpsilon = Epsilon_p-Epsilon_s", &
            usage="dEpsilon 8.295", required=.FALSE., type_of_var=real_t, &
            n_var=1, repeats=.FALSE., default_r_val=8.295_dp, &
            unit_str="eV", error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="E_s", &
            description="Energy of the s-Orbital", &
            usage="E_s 0.0", required=.FALSE., type_of_var=real_t, &
            n_var=1, repeats=.FALSE., default_r_val=0.0_dp, &
            unit_str="eV", error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="E_0", &
            description="Energy shift constant", &
            usage="E_0 0.0", required=.FALSE., type_of_var=real_t, &
            n_var=1, repeats=.FALSE., default_r_val=0.0_dp, &
            unit_str="eV", error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

    END IF
  END SUBROUTINE create_tb_model_section
!***************************************************************************

END MODULE input_cp2k_tbmd
