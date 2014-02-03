!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH FEB-13-2007 : Distributed/replicated realspace grids
!>      Teodoro Laino [tlaino] - University of Zurich - 12.2007
!> \author CJM NOV-30-2003
! *****************************************************************************
MODULE ewald_environment_types
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE f77_blas
  USE input_cp2k_poisson,              ONLY: create_ewald_section
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_get,&
                                             keyword_type
  USE input_section_types,             ONLY: section_get_keyword,&
                                             section_release,&
                                             section_type,&
                                             section_vals_get_subs_vals,&
                                             section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                             do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

! *****************************************************************************
!> \brief to build arrays of pointers
!> \param ewald_env the pointer to the ewald_env
!> \par History
!>      11/03
!> \author CJM
! *****************************************************************************
  TYPE ewald_environment_type
    PRIVATE
    INTEGER   :: id_nr, ref_count
    LOGICAL   :: do_multipoles            ! Flag for using the multipole code
    INTEGER   :: do_ipol                  ! Solver for induced dipoles
    INTEGER   :: max_multipole            ! max expansion in the multipoles
    INTEGER   :: max_ipol_iter            ! max number of interation for induced dipoles
    INTEGER   :: ewald_type               ! type of ewald
    INTEGER   :: gmax(3)                  ! max Miller index
    INTEGER   :: ns_max                   ! # grid points for small grid (PME)
    INTEGER   :: o_spline                 ! order of spline (SPME)
    REAL (KIND=dp) :: precs               ! precision achieved when evaluating the real-space part
    REAL (KIND=dp) :: alpha, rcut         ! ewald alpha and real-space cutoff
    REAL (KIND=dp) :: epsilon             ! tolerance for small grid (PME)
    REAL (KIND=dp) :: eps_pol             ! tolerance for convergence of induced dipoles
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: poisson_section
    ! interaction cutoff is required to make the electrostatic interaction
    ! continuous at a pair distance equal to rcut. this is ignored by the
    ! multipole code and is otherwise only active when SHIFT_CUTOFF is used.
    REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: interaction_cutoffs
  END TYPE ewald_environment_type

! *****************************************************************************
!> \brief to build arrays of pointers
!> \param ewald_env the pointer to the ewald_env
!> \par History
!>      11/03
!> \author CJM
! *****************************************************************************
  TYPE ewald_environment_p_type
     TYPE(ewald_environment_type), POINTER :: ewald_env
  END TYPE ewald_environment_p_type

! *** Public data types ***
  PUBLIC :: ewald_environment_type,        &
            ewald_environment_p_type

! *** Public subroutines ***
  PUBLIC :: ewald_env_get,&
            ewald_env_set,&
            ewald_env_create, &
            ewald_env_retain, &
            ewald_env_release, &
            read_ewald_section

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewald_environment_types'
  INTEGER, PRIVATE, SAVE :: last_ewald_env_id_nr=0

CONTAINS

! *****************************************************************************
!> \brief Purpose: Get the EWALD environment.
!> \param ewald_env the pointer to the ewald_env
!> \par History
!>      11/03
!> \author CJM
! *****************************************************************************
  SUBROUTINE ewald_env_get(ewald_env, ewald_type, alpha, eps_pol, epsilon, &
       gmax, ns_max, o_spline, group, para_env, id_nr, poisson_section, precs, &
       rcut, do_multipoles, max_multipole, do_ipol, max_ipol_iter, &
       interaction_cutoffs, error)
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    INTEGER, OPTIONAL                        :: ewald_type
    REAL(KIND=dp), OPTIONAL                  :: alpha, eps_pol, epsilon
    INTEGER, OPTIONAL                        :: gmax(3), ns_max, o_spline, &
                                                group
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    INTEGER, INTENT(OUT), OPTIONAL           :: id_nr
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: poisson_section
    REAL(KIND=dp), OPTIONAL                  :: precs, rcut
    LOGICAL, INTENT(OUT), OPTIONAL           :: do_multipoles
    INTEGER, INTENT(OUT), OPTIONAL           :: max_multipole, do_ipol, &
                                                max_ipol_iter
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: interaction_cutoffs
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(ewald_env),cp_failure_level,routineP,error,failure)

    IF (PRESENT(id_nr)) id_nr = ewald_env%id_nr
    IF (PRESENT(ewald_type)) ewald_type = ewald_env%ewald_type
    IF (PRESENT(do_multipoles)) do_multipoles = ewald_env%do_multipoles
    IF (PRESENT(do_ipol)) do_ipol = ewald_env%do_ipol
    IF (PRESENT(max_multipole)) max_multipole = ewald_env%max_multipole
    IF (PRESENT(max_ipol_iter)) max_ipol_iter = ewald_env%max_ipol_iter
    IF (PRESENT(alpha)) alpha = ewald_env%alpha
    IF (PRESENT(precs)) precs = ewald_env%precs
    IF (PRESENT(rcut)) rcut = ewald_env%rcut
    IF (PRESENT(epsilon)) epsilon = ewald_env%epsilon
    IF (PRESENT(eps_pol)) eps_pol = ewald_env%eps_pol
    IF (PRESENT(gmax)) gmax = ewald_env%gmax
    IF (PRESENT(ns_max)) ns_max = ewald_env%ns_max
    IF (PRESENT(o_spline)) o_spline = ewald_env%o_spline
    IF (PRESENT(group)) group = ewald_env%para_env%group
    IF (PRESENT(para_env)) para_env => ewald_env%para_env
    IF (PRESENT(poisson_section)) poisson_section => ewald_env%poisson_section
    IF (PRESENT(interaction_cutoffs)) interaction_cutoffs => &
      ewald_env%interaction_cutoffs

  END SUBROUTINE ewald_env_get

! *****************************************************************************
!> \brief Purpose: Set the EWALD environment.
!> \param ewald_env the pointer to the ewald_env
!> \par History
!>      11/03
!> \author CJM
! *****************************************************************************
  SUBROUTINE ewald_env_set(ewald_env, ewald_type, alpha, epsilon, eps_pol, &
       gmax, ns_max, precs, o_spline, para_env, id_nr, poisson_section, &
       interaction_cutoffs, error)

    TYPE(ewald_environment_type), POINTER    :: ewald_env
    INTEGER, OPTIONAL                        :: ewald_type
    REAL(KIND=dp), OPTIONAL                  :: alpha, epsilon, eps_pol
    INTEGER, OPTIONAL                        :: gmax(3), ns_max
    REAL(KIND=dp), OPTIONAL                  :: precs
    INTEGER, OPTIONAL                        :: o_spline
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    INTEGER, INTENT(IN), OPTIONAL            :: id_nr
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: poisson_section
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: interaction_cutoffs
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_env_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(ewald_env),cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       IF (PRESENT(id_nr)) ewald_env%id_nr = id_nr
       IF (PRESENT(ewald_type)) ewald_env%ewald_type = ewald_type
       IF (PRESENT(alpha)) ewald_env%alpha = alpha
       IF (PRESENT(precs)) ewald_env%precs = precs
       IF (PRESENT(epsilon)) ewald_env%epsilon = epsilon
       IF (PRESENT(eps_pol)) ewald_env%eps_pol = eps_pol
       IF (PRESENT(gmax)) ewald_env%gmax = gmax
       IF (PRESENT(ns_max)) ewald_env%ns_max = ns_max
       IF (PRESENT(o_spline)) ewald_env%o_spline = o_spline
       IF (PRESENT(para_env))ewald_env% para_env => para_env
       IF (PRESENT(poisson_section)) THEN
          CALL section_vals_retain(poisson_section,error=error)
          CALL section_vals_release(ewald_env%poisson_section,error=error)
          ewald_env%poisson_section => poisson_section
       END IF
       IF (PRESENT(interaction_cutoffs)) ewald_env%interaction_cutoffs => &
         interaction_cutoffs
    ENDIF
  END SUBROUTINE ewald_env_set

! *****************************************************************************
!> \brief allocates and intitializes a ewald_env
!> \param ewald_env the object to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE ewald_env_create ( ewald_env, para_env, error )
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE( ewald_env, stat=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       ewald_env%ref_count=1
       last_ewald_env_id_nr=last_ewald_env_id_nr+1
       ewald_env%id_nr=last_ewald_env_id_nr
       NULLIFY(ewald_env%poisson_section)
       CALL cp_para_env_retain(para_env, error=error)
       ewald_env%para_env => para_env
       NULLIFY(ewald_env%interaction_cutoffs) ! allocated and initialized later
    END IF
  END SUBROUTINE ewald_env_create

! *****************************************************************************
!> \brief retains the given ewald_env (see doc/ReferenceCounting.html)
!> \param ewald_env the object to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE ewald_env_retain(ewald_env,error)
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(ewald_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(ewald_env%ref_count>0,cp_failure_level,routineP,error,failure)
       ewald_env%ref_count=ewald_env%ref_count+1
    END IF
  END SUBROUTINE ewald_env_retain

! *****************************************************************************
!> \brief releases the given ewald_env (see doc/ReferenceCounting.html)
!> \param ewald_env the object to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE ewald_env_release(ewald_env,error)
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(ewald_env)) THEN
       CPPrecondition(ewald_env%ref_count>0,cp_failure_level,routineP,error,failure)
       ewald_env%ref_count=ewald_env%ref_count-1
       IF (ewald_env%ref_count<1) THEN
          CALL cp_para_env_release ( ewald_env%para_env, error )
          CALL section_vals_release(ewald_env%poisson_section,error=error)
          IF (ASSOCIATED(ewald_env%interaction_cutoffs)) THEN
             DEALLOCATE(ewald_env%interaction_cutoffs, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          DEALLOCATE (  ewald_env, stat = stat )
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF
    END IF
    NULLIFY(ewald_env)
  END SUBROUTINE ewald_env_release

! *****************************************************************************
!> \brief Purpose: read the EWALD section
!> \param ewald_env the pointer to the ewald_env
!> \author Teodoro Laino [tlaino] -University of Zurich - 2005
! *****************************************************************************
  SUBROUTINE read_ewald_section ( ewald_env, ewald_section, error )
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(section_vals_type), POINTER         :: ewald_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_ewald_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iw
    INTEGER, DIMENSION(:), POINTER           :: gmax_read
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: dummy
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: multipole_section

    NULLIFY(enum, keyword, section, multipole_section)
    logger=>cp_error_get_logger(error)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       CALL section_vals_val_get(ewald_section,"EWALD_TYPE",i_val=ewald_env%ewald_type,error=error)
       CALL section_vals_val_get(ewald_section,"ALPHA",r_val=ewald_env%alpha,error=error)
       CALL section_vals_val_get(ewald_section,"EWALD_ACCURACY",r_val=ewald_env%precs,error=error)

       IF (ewald_env%ewald_type==do_ewald_none) THEN
          ewald_env%rcut = 0.0_dp
       ELSE
          CALL section_vals_val_get(ewald_section,"RCUT",explicit=explicit,error=error)
          IF (explicit) THEN
             CALL section_vals_val_get(ewald_section,"RCUT",r_val=ewald_env%rcut,error=error)
          ELSE
             ewald_env%rcut = find_ewald_optimal_value(ewald_env%precs,error)/ewald_env%alpha
          ENDIF
       END IF
       ! we have no defaults for gmax, gmax is only needed for ewald and spme
       SELECT CASE ( ewald_env%ewald_type )
       CASE (do_ewald_ewald, do_ewald_spme)
          CALL section_vals_val_get(ewald_section,"GMAX",i_vals=gmax_read,error=error)
          SELECT CASE (SIZE(gmax_read,1))
          CASE (1)
             ewald_env%gmax = gmax_read(1)
          CASE (3)
             ewald_env%gmax = gmax_read
          CASE DEFAULT
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
          IF (ewald_env%ewald_type==do_ewald_spme) THEN
             CALL section_vals_val_get(ewald_section,"O_SPLINE",i_val=ewald_env%o_spline,error=error)
          END IF
       CASE (do_ewald_pme)
          CALL section_vals_val_get(ewald_section,"NS_MAX",i_val=ewald_env%ns_max,error=error)
          CALL section_vals_val_get(ewald_section,"EPSILON",r_val=ewald_env%epsilon,error=error)
       CASE DEFAULT
          ! this should not be used for do_ewald_none
          ewald_env%gmax   = HUGE(0)
          ewald_env%ns_max = HUGE(0)
       END SELECT

       ! Multipoles
       multipole_section => section_vals_get_subs_vals(ewald_section,"MULTIPOLES",error=error)
       CALL section_vals_val_get(multipole_section,"_SECTION_PARAMETERS_",l_val=ewald_env%do_multipoles,error=error)
       CALL section_vals_val_get(multipole_section,"POL_SCF",i_val=ewald_env%do_ipol,error=error)
       CALL section_vals_val_get(multipole_section,"EPS_POL",r_val=ewald_env%eps_pol,error=error)
       IF (ewald_env%do_multipoles) THEN
          SELECT CASE(ewald_env%ewald_type)
          CASE(do_ewald_ewald)
             CALL section_vals_val_get(multipole_section,"MAX_MULTIPOLE_EXPANSION",i_val=ewald_env%max_multipole,error=error)
             CALL section_vals_val_get(multipole_section,"MAX_IPOL_ITER",i_val=ewald_env%max_ipol_iter,error=error)
          CASE DEFAULT
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Multipole code works at the moment only with standard EWALD sums.", &
                  error=error, error_level=cp_failure_level)
          END SELECT
       END IF

       iw=cp_print_key_unit_nr(logger,ewald_section,"PRINT%PROGRAM_RUN_INFO",&
            extension=".log",error=error)
       IF ( iw > 0 ) THEN
          NULLIFY(keyword, enum)
          CALL create_ewald_section(section,error=error)
          IF ( ewald_env%ewald_type /= do_ewald_none ) THEN
             keyword => section_get_keyword(section,"EWALD_TYPE",error=error)
             CALL keyword_get(keyword,enum=enum,error=error)
             WRITE ( iw, '(/,T2,"EWALD| ",A,T67,A14 )')'Summation is done by:', &
                  ADJUSTR(TRIM(enum_i2c(enum,ewald_env%ewald_type,error=error)))
             IF (ewald_env%do_multipoles) THEN
                NULLIFY(keyword, enum)
                keyword => section_get_keyword(section,"MULTIPOLES%MAX_MULTIPOLE_EXPANSION",error=error)
                CALL keyword_get(keyword,enum=enum,error=error)
                WRITE ( iw, '( T2,"EWALD| ",A )' ) 'Enabled Multipole Method'
                WRITE ( iw, '( T2,"EWALD| ",A,T67,A14 )' ) 'Max Term in Multipole Expansion :',&
                     ADJUSTR(TRIM(enum_i2c(enum,ewald_env%max_multipole,error=error)))
                WRITE ( iw, '( T2,"EWALD| ",A,T67,3I10 )' ) 'Max number Iterations for IPOL :',&
                     ewald_env%max_ipol_iter
             END IF
             dummy = cp_unit_from_cp2k(ewald_env%alpha,"angstrom^-1",error=error)
             WRITE ( iw, '( T2,"EWALD| ",A,A18,A,T71,F10.4 )' ) &
                  'Alpha parameter [','ANGSTROM^-1',']',dummy
             dummy = cp_unit_from_cp2k(ewald_env%rcut,"angstrom",error=error)
             WRITE ( iw, '( T2,"EWALD| ",A,A18,A,T71,F10.4 )' ) &
                  'Real Space Cutoff [','ANGSTROM',']',dummy

             SELECT CASE ( ewald_env%ewald_type )
             CASE ( do_ewald_ewald)
                WRITE ( iw, '( T2,"EWALD| ",A,T51,3I10 )' ) &
                     'G-space max. Miller index', ewald_env%gmax
             CASE ( do_ewald_pme)
                WRITE ( iw, '( T2,"EWALD| ",A,T71,I10 )' ) &
                     'Max small-grid points (input) ', ewald_env%ns_max
                WRITE ( iw, '( T2,"EWALD| ",A,T71,E10.4 )' ) &
                     'Gaussian tolerance (input) ', ewald_env%epsilon
             CASE ( do_ewald_spme )
                WRITE ( iw, '( T2,"EWALD| ",A,T51,3I10 )' ) &
                     'G-space max. Miller index', ewald_env%gmax
                WRITE ( iw, '( T2,"EWALD| ",A,T71,I10 )' ) &
                     'Spline interpolation order ', ewald_env%o_spline
             CASE DEFAULT
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
          ELSE
             WRITE ( iw, '( T2,"EWALD| ",T73, A )' ) 'not used'
          END IF
          CALL section_release(section,error=error)
       END IF
       CALL cp_print_key_finished_output(iw,logger,ewald_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)
    END IF

  END SUBROUTINE read_ewald_section

! *****************************************************************************
!> \brief triggers (by bisection) the optimal value for EWALD parameter x
!>      EXP(-x^2)/x^2 = EWALD_ACCURACY
!> \param ewald_env the pointer to the ewald_env
!> \author Teodoro Laino [tlaino] - University of Zurich - 12.2007
! *****************************************************************************
  FUNCTION find_ewald_optimal_value(precs,error) RESULT(value)
    REAL(KIND=dp)                            :: precs
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: value

    CHARACTER(len=*), PARAMETER :: routineN = 'find_ewald_optimal_value', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: func, func1, func2, s, s1, s2

    failure = .FALSE.
    s    =  0.1_dp
    func =  EXP(-s**2)/s**2 - precs
    CPPrecondition(func>0.0_dp,cp_failure_level,routineP,error,failure)
    DO WHILE (func>0.0_dp)
       s = s + 0.1_dp
       func =  EXP(-s**2)/s**2 - precs
    END DO
    s2 = s
    s1 = s - 0.1_dp
    ! Start bisection
    DO WHILE (.TRUE.)
       func2 = EXP(-s2**2)/s2**2 - precs
       func1 = EXP(-s1**2)/s1**2 - precs
       CPPrecondition(func1>=0,cp_failure_level,routineP,error,failure)
       CPPrecondition(func2<=0,cp_failure_level,routineP,error,failure)
       s    = 0.5_dp * (s1 + s2)
       func = EXP(-s**2)/s**2 - precs
       IF      (func > 0.0_dp) THEN
          s1 = s
       ELSE IF (func < 0.0_dp) THEN
          s2 = s
       END IF
       IF (ABS(func)<100.0_dp*EPSILON(0.0_dp)) EXIT
    END DO
    value = s
  END FUNCTION find_ewald_optimal_value

END MODULE ewald_environment_types

