!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief function that build the kpoints section of the input
!> \par History
!>      init [07.2014]
!> \author JGH
! *****************************************************************************
MODULE input_cp2k_kpoints
  
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_create,&
                                             section_type
  USE input_val_types,                 ONLY: char_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_kpoints'

  INTEGER, PARAMETER                      :: use_real_wfn=101,&
                                             use_complex_wfn=100

  PUBLIC :: create_kpoints_section
  PUBLIC :: use_real_wfn, use_complex_wfn

CONTAINS

! *****************************************************************************
!> \brief Creates the Kpoints section
!>  SECTION: &kpoint... &end
!>      SCHEME [None, Gamma, Monkhorst-Pack, MacDonald, General]
!>        { nx ny nz }
!>        { nx ny nz  sx sy sz }
!>      KPOINT           x1 y1 z1 w1
!>      SYMMETRY         [on, off]
!>      WAVEFUNCTION     [real, complex]
!>      FULL_GRID        [on, off]        
!>      VERBOSE          [on, off]        
!>      EPS_GEO          value
!>      PARALLEL_GROUP_SIZE   [-1,0,n]
!>
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author JGH
! *****************************************************************************
  SUBROUTINE create_kpoints_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_kpoints_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="KPOINTS",&
            description="Sets up the kpoints.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY (keyword)
       CALL keyword_create(keyword, name="SCHEME",&
            description="Kpoint scheme to be used. ",&
            usage="SCHEME {KPMETHOD}{integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KPOINT",&
            description="Specify kpoint coordinates and weight. ",&
            usage="KPOINT x  y  z  w",repeats=.TRUE.,&
            n_var=4, type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SYMMETRY",&
            description="Use symmetry to reduce the number of kpoints.",&
            usage="SYMMETRY <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FULL_GRID",&
            description="Use full non-reduced kpoint grid.",&
            usage="FULL_GRID <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VERBOSE",&
            description="Verbose output information.",&
            usage="VERBOSE <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_GEO",&
            description="Accuracy in symmetry determination.",&
            usage="EPS_GEO <real>",&
            default_r_val=1.0e-6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARALLEL_GROUP_SIZE",&
            description="Number of processors to be used for a single kpoint."//&
                        " Value=-1 (minimum number of processes)."//&
                        " Value=0 (maximum number of processes)."//&
                        " Value=n (exactly n processes).",&
            usage="PARALLEL_GROUP_SIZE <integer>",&
            default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WAVEFUNCTIONS",&
            description="Use real/complex wavefunctions if possible.",&
            usage="WAVEFUNCTIONS REAL",&
            default_i_val=use_complex_wfn,&
            enum_c_vals=s2a("REAL","COMPLEX"),&
            enum_desc=s2a("Use real wavefunctions (if possible by kpoints specified).",&
                          "Use complex wavefunctions (default)."),&
            enum_i_vals=(/ use_real_wfn, use_complex_wfn/), error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_kpoints_section

END MODULE input_cp2k_kpoints
