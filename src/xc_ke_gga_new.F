!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/xc_ke_gga [1.0] *
!!
!!   NAME
!!     xc_ke_gga
!!
!!   FUNCTION
!!     Calculate the several different kinetic energy functionals
!!     with a GGA form 
!!
!!   AUTHOR
!!     JGH (20.02.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!     fawzi (04.2004)  : adapted to the new xc interface
!!
!******************************************************************************

MODULE xc_ke_gga_new
! *****************************************************************************
  USE cp_array_r_utils,                ONLY: cp_3d_r_p_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_functionals_utilities,        ONLY: calc_wave_vector,&
                                             set_util
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_type,&
                                             xc_rho_set_get,&
                                             xc_rho_set_type

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f43 = 4._dbl*f13, &
                          f53 = 5._dbl*f13

  PUBLIC :: ke_gga_info, ke_gga_lda_eval, ke_gga_lsd_eval

  INTEGER, PARAMETER, PUBLIC :: ke_ol1=1, ke_ol2=2,ke_llp=3,ke_pw86=4,&
       ke_pw91=5, ke_lc=6,ke_t92=7

  REAL(dbl) :: cf, b, flda, flsd, sfac, t13
  REAL(dbl) :: fact, tact
  REAL(dbl) :: eps_rho
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="xc_ke_gga_new"
!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE ke_gga_init ( cutoff )

    REAL(dbl), INTENT(IN)                    :: cutoff

    eps_rho = cutoff
    CALL set_util ( cutoff )

    cf = 0.3_dbl*(3._dbl*pi*pi)**f23
    flda = cf
    flsd = flda * 2._dbl**f23
    b = 2._dbl**f43 * (3._dbl*pi*pi)**(-f13)
    sfac = 1._dbl / (2._dbl*(3._dbl*pi*pi)**f13)
    t13 = 2._dbl**f13

  END SUBROUTINE ke_gga_init

! *****************************************************************************

  SUBROUTINE ke_gga_info ( functional,lsd,reference, shortform, needs, max_deriv, error )
    INTEGER, INTENT(in)                      :: functional
    LOGICAL, INTENT(in)                      :: lsd
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "ke_gga_info", &
      routineP = moduleN//":"//routineN

    IF ( PRESENT ( reference ) ) THEN
       SELECT CASE(functional)
       CASE(ke_ol1)
          reference = "H. Ou-Yang and M. Levy, "//&
               "Intl. J. Quant. Chem. 40, 379 (1991); Functional 1"
       CASE(ke_ol2)
          reference = "H. Ou-Yang and M. Levy, "//&
               "Intl. J. Quant. Chem. 40, 379 (1991); Functional 2"
       CASE(ke_llp)
          reference = "H. Lee, C. Lee, R.G. Parr, Phys. Rev. A, 44, 768 (1991)"
       CASE(ke_pw86)
          reference = "J.P. Perdew and Y. Wang, Phys. Rev. B, 33, 8800 (1986)"
       CASE(ke_pw91)
          reference = "J.P. Perdew and Y. Wang, Electronic Structure of Solids 91"
       CASE(ke_lc)
          reference = "A. Lembarki and H. Chermette, Phys. Rev. A, 50, 5328 (1994)"
       CASE(ke_t92)
          reference = "A.J. Thakkar, Phys. Rev. A, 46, 6920 (1992)"
       END SELECT
       IF (.not.lsd) THEN
          IF (LEN_TRIM(reference)+6<LEN(reference)) THEN
             reference(LEN_TRIM(reference):LEN_TRIM(reference)+6)=' {LDA}'
          END IF
       END IF
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       SELECT CASE(functional)
       CASE(ke_ol1)
          shortform = "Ou-Yang-Levy Functional 1"
       CASE(ke_ol2)
          shortform = "Ou-Yang-Levy Functional 2"
       CASE(ke_llp)
          shortform = "Lee-Lee-Parr Functional"
       CASE(ke_pw86)
          shortform = "Perdew-Wang 1986 Functional (kinetic energy)"
       CASE(ke_pw91)
          shortform = "Perdew-Wang 1991 Functional (kinetic energy)"
       CASE(ke_lc)
          shortform = "Lembarki-Chermette kinetic energy functional"
       CASE(ke_t92)
          shortform = "Thakkar 1992 Functional"
       END SELECT
       IF (.not.lsd) THEN
          IF (LEN_TRIM(shortform)+6<LEN(shortform)) THEN
             shortform(LEN_TRIM(shortform):LEN_TRIM(shortform)+6)=' {LDA}'
          END IF
       END IF
    END IF
    IF (PRESENT(needs)) THEN
       IF (lsd) THEN
          needs%rho_spin=.TRUE.
          needs%rho_spin_1_3=.TRUE.
          needs%norm_drho=.TRUE.
       ELSE
          needs%rho=.TRUE.
          needs%rho_1_3=.TRUE.
          needs%norm_drho=.TRUE.
       END IF
    END IF
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE ke_gga_info

! *****************************************************************************

  SUBROUTINE ke_gga_lda_eval ( functional, rho_set, deriv_set, order, error )

    !   ---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                      :: functional
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: order
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "ke_gga_lda_eval", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, m, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(dbl)                                :: drho_cutoff, rho_cutoff
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: s
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: fs
    REAL(dbl), DIMENSION(:, :, :), POINTER :: e_0, e_ndrho, e_ndrho_ndrho, &
      e_ndrho_ndrho_ndrho, e_rho, e_rho_ndrho, e_rho_ndrho_ndrho, e_rho_rho, &
      e_rho_rho_ndrho, e_rho_rho_rho, grho, rho, rho13
    TYPE(xc_derivative_type), POINTER        :: deriv

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(bo,rho, rho13,e_0,e_rho,e_ndrho,&
         e_rho_rho, e_rho_ndrho, e_ndrho_ndrho,&
         e_rho_rho_rho, e_rho_rho_ndrho, e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho)
    m = ABS(order)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       CALL xc_rho_set_get(rho_set,rho_1_3=rho13,rho=rho,&
            norm_drho=grho,local_bounds=bo,rho_cutoff=rho_cutoff,&
            drho_cutoff=drho_cutoff,error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)
       CALL ke_gga_init(rho_cutoff)

       ALLOCATE ( s(npoints), STAT=stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( fs(npoints,m+1), STAT=stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       CALL calc_wave_vector ( "u", rho, grho, s )

       SELECT CASE (functional)
       CASE (ke_ol1)
          CALL efactor_ol1 ( s, fs, m )
       CASE (ke_ol2)
          CALL efactor_ol2 ( s, fs, m )
       CASE (ke_llp)
          CALL efactor_llp ( s, fs, m )
       CASE (ke_pw86)
          CALL efactor_pw86 ( s, fs, m )
       CASE (ke_pw91)
          CALL efactor_pw91 ( s, fs, m, 1 )
       CASE (ke_lc)
          CALL efactor_pw91 ( s, fs, m, 2 )
       CASE (ke_t92)
          CALL efactor_t92 ( s, fs, m )
       CASE DEFAULT
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT

       fact = flda
       tact = t13

       IF ( order>=0 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE., error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)

          CALL kex_p_0 ( rho, rho13, fs, e_0, npoints, error=error )
       END IF

       IF ( order>=1 .OR. order==-1 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)

          CALL kex_p_1 ( rho, grho, rho13, s, fs, e_rho=e_rho, e_ndrho=e_ndrho,&
               npoints=npoints, error=error)
       END IF
       IF ( order>=2 .OR. order==-2 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(rho)(norm_drho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_ndrho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drho)(norm_drho)", allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho,error=error)

          CALL kex_p_2 ( rho, grho, rho13, s, fs, e_rho_rho=e_rho_rho,&
               e_rho_ndrho=e_rho_ndrho, e_ndrho_ndrho=e_ndrho_ndrho, npoints=npoints,&
               error=error )
       END IF
       IF ( order>=3 .OR. order==-3 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(rho)(rho)(norm_drho)",allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_ndrho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(rho)(norm_drho)(norm_drho)",allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_ndrho_ndrho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drho)(norm_drho)(norm_drho)", allocate_deriv=.TRUE.,&
               error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho,error=error)

          CALL kex_p_3 ( rho, grho, rho13, s, fs, e_rho_rho_rho=e_rho_rho_rho,&
               e_rho_rho_ndrho=e_rho_rho_ndrho, e_rho_ndrho_ndrho=e_rho_ndrho_ndrho,&
               e_ndrho_ndrho_ndrho=e_ndrho_ndrho_ndrho,npoints=npoints,error=error)
       END IF
       IF ( order>3.OR.order<-3) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 3 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       DEALLOCATE ( s, STAT=stat )
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       DEALLOCATE ( fs, STAT=stat )
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE ke_gga_lda_eval

! *****************************************************************************

  SUBROUTINE ke_gga_lsd_eval ( functional, rho_set, deriv_set, order, &
       error )

    !   ---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                      :: functional
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN), OPTIONAL            :: order
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "ke_gga_lsd_eval", &
      routineP = moduleN//":"//routineN

    CHARACTER(len=12), DIMENSION(2) :: &
      norm_drho_spin_name = (/"(norm_drhoa)","(norm_drhob)"/)
    CHARACTER(len=6), DIMENSION(2) :: rho_spin_name = (/"(rhoa)","(rhob)"/)
    INTEGER                                  :: handle, ispin, m, npoints, &
                                                stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(dbl)                                :: rho_cutoff
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: s
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: fs
    REAL(dbl), DIMENSION(:, :, :), POINTER :: e_0, e_ndrho, e_ndrho_ndrho, &
      e_ndrho_ndrho_ndrho, e_rho, e_rho_ndrho, e_rho_ndrho_ndrho, e_rho_rho, &
      e_rho_rho_ndrho, e_rho_rho_rho
    TYPE(cp_3d_r_p_type), DIMENSION(2)       :: norm_drho, rho, rho_1_3
    TYPE(xc_derivative_type), POINTER        :: deriv

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(e_0,e_ndrho, e_ndrho_ndrho, e_ndrho_ndrho_ndrho, e_rho_ndrho_ndrho, &
         e_rho_ndrho, e_rho_rho_ndrho, e_rho, e_rho_rho, e_rho_rho_rho)
    DO ispin=1,2
       NULLIFY(norm_drho(ispin)%array, rho(ispin)%array, rho_1_3(ispin)%array)
    END DO

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       CALL xc_rho_set_get(rho_set,rhoa_1_3=rho_1_3(1)%array,&
            rhob_1_3=rho_1_3(2)%array,rhoa=rho(1)%array,&
            rhob=rho(2)%array,norm_drhoa=norm_drho(1)%array, &
            norm_drhob=norm_drho(2)%array,rho_cutoff=rho_cutoff,&
            local_bounds=bo, error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)
       m=ABS(order)
       CALL ke_gga_init(rho_cutoff)

       ALLOCATE ( s(npoints), STAT=stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( fs(npoints,m+1), STAT=stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       fact = flsd
       tact = 1._dbl

       DO ispin=1,2
          CALL calc_wave_vector ( "p", rho(ispin)%array, norm_drho(ispin)%array, s )

          SELECT CASE (functional)
          CASE (ke_ol1)
             CALL efactor_ol1 ( s, fs, m )
          CASE (ke_ol2)
             CALL efactor_ol2 ( s, fs, m )
          CASE (ke_llp)
             CALL efactor_llp ( s, fs, m )
          CASE (ke_pw86)
             CALL efactor_pw86 ( s, fs, m )
          CASE (ke_pw91)
             CALL efactor_pw91 ( s, fs, m, 1 )
          CASE (ke_lc)
             CALL efactor_pw91 ( s, fs, m, 2 )
          CASE (ke_t92)
             CALL efactor_t92 ( s, fs, m )
          CASE DEFAULT
             CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT

          IF ( order>=0 ) THEN
             deriv => xc_dset_get_derivative(deriv_set,"",&
                  allocate_deriv=.TRUE., error=error)
             CALL xc_derivative_get(deriv, deriv_data=e_0,error=error)

             CALL kex_p_0 ( rho(ispin)%array, rho_1_3(ispin)%array, fs, &
                  e_0=e_0, npoints=npoints, error=error )
          END IF
          IF ( order>=1.or.order==-1 ) THEN
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin),&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,norm_drho_spin_name(ispin),&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)

             CALL kex_p_1 ( rho=rho(ispin)%array, grho=norm_drho(ispin)%array, &
                  r13=rho_1_3(ispin)%array, s=s, fs=fs, e_rho=e_rho, &
                  e_ndrho=e_ndrho, npoints=npoints, error=error )
          END IF
          IF ( order>=2.OR.order==-2 ) THEN
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin)//&
                  rho_spin_name(ispin),allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin)//&
                  norm_drho_spin_name(ispin),allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_ndrho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,norm_drho_spin_name(ispin)//&
                  norm_drho_spin_name(ispin), allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho,error=error)

             CALL kex_p_2 ( rho(ispin)%array, norm_drho(ispin)%array, rho_1_3(ispin)%array,&
                  s, fs, e_rho_rho,e_rho_ndrho,e_ndrho_ndrho,npoints,error=error )
          END IF
          IF ( order>=3 .OR. order==-3 ) THEN
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin)//&
                  rho_spin_name(ispin)//rho_spin_name(ispin),&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin)//&
                  rho_spin_name(ispin)//norm_drho_spin_name(ispin),&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_ndrho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin)//&
                  norm_drho_spin_name(ispin)//norm_drho_spin_name(ispin), &
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_ndrho_ndrho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,norm_drho_spin_name(ispin)//&
                  norm_drho_spin_name(ispin)//norm_drho_spin_name(ispin),&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho,error=error)

             CALL kex_p_3 ( rho(ispin)%array, norm_drho(ispin)%array, &
                  rho_1_3(ispin)%array, s, fs, e_rho_rho_rho, e_rho_rho_ndrho,&
                  e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho, npoints, error)
          END IF
          IF ( order>3.OR.order<-3) THEN
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="derivatives bigger than 3 not implemented", &
                  error=error, error_level=cp_failure_level)
          END IF

       END DO

       DEALLOCATE ( s, STAT=stat )
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       DEALLOCATE ( fs, STAT=stat )
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    CALL timestop(0.0_dbl,handle)
  END SUBROUTINE ke_gga_lsd_eval

! *****************************************************************************

  SUBROUTINE kex_p_0 ( rho, r13, fs, e_0, npoints, error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, r13
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(*), INTENT(INOUT)   :: e_0
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ip

!   ---------------------------------------------------------------------------

    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         e_0(ip) = e_0(ip) + fact*r13(ip)*r13(ip)*rho(ip) * fs(ip,1)
      END IF

    END DO

  END SUBROUTINE kex_p_0

! *****************************************************************************

  SUBROUTINE kex_p_1 ( rho, grho, r13, s, fs, e_rho, e_ndrho, npoints, error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(*), INTENT(INOUT)   :: e_rho, e_ndrho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: a0, a1, sx, sy

!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,a0,a1,sx,sy)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*r13(ip)*rho(ip)
         a1 = f53*fact*r13(ip)*r13(ip)
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         e_rho(ip) = e_rho(ip) + a1*fs(ip,1) + a0*fs(ip,2)*sx
         e_ndrho(ip) = e_ndrho(ip) + a0*fs(ip,2)*sy

      END IF

    END DO

  END SUBROUTINE kex_p_1

! *****************************************************************************

  SUBROUTINE kex_p_2 ( rho, grho, r13, s, fs, e_rho_rho, e_rho_ndrho,e_ndrho_ndrho,&
      npoints,error)

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(*), INTENT(INOUT)   :: e_rho_rho, e_rho_ndrho, &
                                                e_ndrho_ndrho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: a0, a1, a2, sx, sxx, sxy, sy

!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,a0,a1,a2,sx,sy,sxx,sxy)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*r13(ip)*rho(ip)
         a1 = f53*fact*r13(ip)*r13(ip)
         a2 = f23*f53*fact/r13(ip)
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         sxx= 28._dbl/9._dbl*s(ip)/(rho(ip)*rho(ip))
         sxy= -f43*sfac*tact/(r13(ip)*rho(ip)*rho(ip))
         e_rho_rho(ip) = e_rho_rho(ip) + a2*fs(ip,1) +2._dbl*a1*fs(ip,2)*sx + &
                     a0*fs(ip,3)*sx*sx + a0*fs(ip,2)*sxx
         e_rho_ndrho(ip) = e_rho_ndrho(ip) + a1*fs(ip,2)*sy + a0*fs(ip,3)*sx*sy + &
                     a0*fs(ip,2)*sxy
         e_ndrho_ndrho(ip) = e_ndrho_ndrho(ip) + a0*fs(ip,3)*sy*sy

      END IF

    END DO

  END SUBROUTINE kex_p_2

! *****************************************************************************

  SUBROUTINE kex_p_3 ( rho, grho, r13, s, fs, e_rho_rho_rho,e_rho_rho_ndrho,&
      e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho, npoints, error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(*), INTENT(inout)   :: e_rho_rho_rho, &
                                                e_rho_rho_ndrho, &
                                                e_rho_ndrho_ndrho, &
                                                e_ndrho_ndrho_ndrho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: a0, a1, a2, a3, sx, sxx, &
                                                sxxx, sxxy, sxy, sy

!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,a0,a1,a2,a3,sx,sy,sxx,sxy,sxxx,sxxy)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*r13(ip)*rho(ip)
         a1 = f53*fact*r13(ip)*r13(ip)
         a2 = f23*f53*fact/r13(ip)
         a3 = -f13*f23*f53*fact/(r13(ip)*rho(ip))
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         sxx= 28._dbl/9._dbl*s(ip)/(rho(ip)*rho(ip))
         sxy= -f43*sfac*tact/(r13(ip)*rho(ip)*rho(ip))
         sxxx= -280._dbl/27._dbl*s(ip)/(rho(ip)*rho(ip)*rho(ip))
         sxxy= 28._dbl/9._dbl*sfac*tact/(r13(ip)*rho(ip)*rho(ip)*rho(ip))
         e_rho_rho_rho(ip) = e_rho_rho_rho(ip) + a3*fs(ip,1) + 3._dbl*a2*fs(ip,2)*sx + &
                     3._dbl*a1*fs(ip,3)*sx*sx + 3._dbl*a1*fs(ip,2)*sxx + &
                     a0*fs(ip,4)*sx*sx*sx + 3._dbl*a0*fs(ip,3)*sx*sxx + &
                     a0*fs(ip,2)*sxxx
         e_rho_rho_ndrho(ip) = e_rho_rho_ndrho(ip) + a2*fs(ip,2)*sy + 2._dbl*a1*fs(ip,3)*sx*sy + &
                     2._dbl*a1*fs(ip,2)*sxy + a0*fs(ip,4)*sx*sx*sy + &
                     2._dbl*a0*fs(ip,3)*sx*sxy + a0*fs(ip,3)*sxx*sy + &
                     a0*fs(ip,2)*sxxy
         e_rho_ndrho_ndrho(ip) = e_rho_ndrho_ndrho(ip) + a1*fs(ip,3)*sy*sy + a0*fs(ip,4)*sx*sy*sy + &
                     2._dbl*a0*fs(ip,3)*sxy*sy
         e_ndrho_ndrho_ndrho(ip) = e_ndrho_ndrho_ndrho(ip) + a0*fs(ip,4)*sy*sy*sy

      END IF

    END DO

  END SUBROUTINE kex_p_3

! *****************************************************************************
! Enhancement Factors
! *****************************************************************************
  SUBROUTINE efactor_ol1 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl)                                :: t1, t2

    t1 = b*b/(72._dbl*cf)
    t2 = 0.001878_dbl*b
!$omp parallel do private(ip)
    DO ip=1,SIZE(s)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)
      CASE (1)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)
         fs(ip,2) = 2._dbl*t1*s(ip) + t2
      CASE (2:3)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)
         fs(ip,2) = 2._dbl*t1*s(ip) + t2
         fs(ip,3) = 2._dbl*t1
      CASE DEFAULT
         CALL stop_program("efactor_ol1","Illegal order")
      END SELECT
    END DO
    IF(m==3) fs(:,4)=0._dbl

  END SUBROUTINE efactor_ol1
! *****************************************************************************
  SUBROUTINE efactor_ol2 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl)                                :: t1, t2, t3, y

    t1 = b*b/(72._dbl*cf)
    t2 = 0.0245_dbl*b
    t3 = 2._dbl**f53*b
!$omp parallel do private(ip,y)
    DO ip=1,SIZE(s)
      y = 1._dbl/(1._dbl+t3*s(ip))
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)*y
      CASE (1)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)*y
         fs(ip,2) = 2._dbl*t1*s(ip) + t2*y*y
      CASE (2)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)*y
         fs(ip,2) = 2._dbl*t1*s(ip) + t2*y*y
         fs(ip,3) = 2._dbl*(t1 - t2*t3*y*y*y)
      CASE (3)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)*y
         fs(ip,2) = 2._dbl*t1*s(ip) + t2*y*y
         fs(ip,3) = 2._dbl*(t1 - t2*t3*y*y*y)
         fs(ip,4) = 6._dbl*t2*t3*t3*y*y*y*y
      CASE DEFAULT
         CALL stop_program("efactor_ol2","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_ol2
! *****************************************************************************
  SUBROUTINE efactor_llp ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl) :: as, asp, asp2, asp3, bs, p, q, sas, sbs, sbs3, sbs5, t1, &
      t11, t12, t13, t16, t19, t22, t24, t25, t3, t33, t34, t36, t39, t4, &
      t43, t44, t60, t62, w1, x, ys

    p = 0.0045_dbl*b*b
    q = 0.0253_dbl*b
!$omp parallel do private(ip,x,bs,sbs,as,sas,ys,asp,sbs3,asp2,asp3,sbs5,w1) &
!$omp             private(t1,t3,t4,t11,t12,t13,t16,t19,t22,t24,t25,t33,t34,t36,t39) &
!$omp             private(t43,t44,t60,t62)
    DO ip=1,SIZE(s)
      x = s(ip)
      bs = b*x
      sbs = SQRT(bs*bs+1._dbl)
      as = LOG(bs+sbs)
      sas = x*as
      ys =  1._dbl/(1._dbl+q*sas)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1._dbl + p*x*x*ys
      CASE (1)
         asp = as + bs / sbs
         fs(ip,1) = 1._dbl + p*x*x*ys
         fs(ip,2) = 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys
      CASE (2)
         asp = as + bs / sbs
         sbs3 = sbs*sbs*sbs
         asp2 = 2._dbl*b / sbs - b*bs*bs/ sbs3
         fs(ip,1) = 1._dbl + p*x*x*ys
         fs(ip,2) = 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys
         fs(ip,3) = 2._dbl*p*ys - p*q*x*(4._dbl*asp+x*asp2)*ys*ys + &
                    2._dbl*p*q*q*x*x*asp*asp*ys*ys*ys
      CASE (3)
         asp = as + bs / sbs
         sbs3 = sbs*sbs*sbs
         sbs5 = sbs3*sbs*sbs
         asp2 = 2._dbl*b / sbs - b*bs*bs/ sbs3
         asp3 = -4._dbl*b*b*bs / sbs3 + 3._dbl*b*b*bs*bs*bs/ sbs5
         fs(ip,1) = 1._dbl + p*x*x*ys
         fs(ip,2) = 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys
         w1 = (4._dbl*asp+x*asp2)
         fs(ip,3) = 2._dbl*p*ys - p*q*x*w1*ys*ys + &
                    2._dbl*p*q*q*x*x*asp*asp*ys*ys*ys
         t1 = q*x; t3 = b**2; t4 = x**2; t11 = 1+t1*as; t12 = t11**2; t13 = 1/t12
         t16 = 1/sbs; t19 = q*as+t1*b*t16; t22 = p*x; t24 = 1/t12/t11; t25 = t19**2
         t33 = t3*b; t34 = sbs**2; t36 = 1/t34/sbs; t39 = 2*q*b*t16-q*t4*t33*t36
         t43 = p*t4; t44 = t12**2; t60 = t3**2; t62 = t34**2
         fs(ip,4) = &
            -6*p*t13*t19+12*t22*t24*t25-6*t22*t13*t39-6*t43/t44*t25*t19+ &
            6*t43*t24*t19*t39-t43*t13*(-4*q*t33*t36*x+3*q*t4*x*t60*b/t62/sbs)

      CASE DEFAULT
         CALL stop_program("efactor_llp","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_llp
! *****************************************************************************
  SUBROUTINE efactor_pw86 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl)                                :: f15, p0, p1, p15, s2, s4, s6, &
                                                t1, t10, t12, t13, t14, t19, &
                                                t2, t25, t3, t8, t9

    t1 = 1.296_dbl
    t2 = 14._dbl
    t3 = 0.2_dbl
    f15 = 1._dbl/15._dbl
!$omp parallel do private(ip,s2,s4,s6,p0,p1,p15,t8,t9,t10,t12,t13,t14,t19,t25)
    DO ip=1,SIZE(s)
      s2 = s(ip)*s(ip)
      s4 = s2*s2
      s6 = s2*s4
      SELECT CASE (m)
      CASE (0)
         p0 = 1._dbl + t1*s2 + t2*s4 + t3*s6
         fs(ip,1) = p0**f15
      CASE (1)
         p0 = 1._dbl + t1*s2 + t2*s4 + t3*s6
         p1 = s(ip)*(2._dbl*t1 + 4._dbl*t2*s2 + 6._dbl*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
      CASE (2)
         p0 = 1._dbl + t1*s2 + t2*s4 + t3*s6
         p1 = s(ip)*(2._dbl*t1 + 4._dbl*t2*s2 + 6._dbl*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
         t9 = p15**2; t10 = t9**2; t12 = t10**2; t13 = t12*t10*t9
         t25 = p1*p1
         fs(ip,3) = -14._dbl/225._dbl/t13/p0*t25 + &
                     1._dbl/t13*(2._dbl*t1+12*t2*s2+30._dbl*t3*s4)/15._dbl
      CASE (3)
         p0 = 1._dbl + t1*s2 + t2*s4 + t3*s6
         p1 = s(ip)*(2._dbl*t1 + 4._dbl*t2*s2 + 6._dbl*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
         t9 = p15**2; t10 = t9**2; t12 = t10**2; t13 = t12*t10*t9
         t25 = p1*p1
         fs(ip,3) = -14._dbl/225._dbl/t13/p0*t25 + &
                     1._dbl/t13*(2._dbl*t1+12*t2*s2+30._dbl*t3*s4)/15._dbl
         t8 = p0**2; t9 = p0**f15; t14 = p0/t9; t19 = s2*s(ip)
         fs(ip,4) = 406._dbl/3375._dbl/t14/t8*p1*p1*p1-14._dbl/&
                    75._dbl/t14/p0*p1*(2*t1+12*t2*s2+30*t3*s4)+&
                    1/t14*(24*t2*s(ip)+120*t3*t19)*f15
      CASE DEFAULT
         CALL stop_program("efactor_ol1","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_pw86
! *****************************************************************************
  SUBROUTINE efactor_t92 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl)                                :: a1, a2, as, asp, asp2, asp3, &
                                                bs, p, q, s1, s2, sas, sbs, &
                                                sbs3, sbs5, t0, w1, x, ys

    p = 0.0055_dbl*b*b
    q = 0.0253_dbl*b
    a1 = 0.072_dbl * b
    a2 = 2._dbl**f53 * b
!$omp parallel do private(ip,x,as,sas,asp,asp2,asp3,bs,ys,sbs,sbs3,sbs5,w1,t0,s1,s2)
    DO ip=1,SIZE(s)
      x = s(ip)
      bs = b*x
      sbs = SQRT(bs*bs+1._dbl)
      as = LOG(bs+sbs)
      sas = x*as
      ys =  1._dbl/(1._dbl+q*sas)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1._dbl + p*x*x*ys - a1*x/(1+a2*x)
      CASE (1)
         asp = as + bs / sbs
         fs(ip,1) = 1._dbl + p*x*x*ys - a1*x/(1+a2*x)
         fs(ip,2) = 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys - a1/(1+a2*x)**2
      CASE (2)
         asp = as + bs / sbs
         sbs3 = sbs*sbs*sbs
         asp2 = 2._dbl*b / sbs - b*bs*bs/ sbs3
         fs(ip,1) = 1._dbl + p*x*x*ys - a1*x/(1+a2*x)
         fs(ip,2) = 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys - a1/(1+a2*x)**2
         fs(ip,3) = 2._dbl*p*ys - p*q*x*(4._dbl*asp+x*asp2)*ys*ys + &
                    2._dbl*p*q*q*x*x*asp*asp*ys*ys*ys + 2._dbl*a1*a2/(1+a2*x)**3
      CASE (3)
         asp = as + bs / sbs
         sbs3 = sbs*sbs*sbs
         sbs5 = sbs3*sbs*sbs
         asp2 = 2._dbl*b / sbs - b*bs*bs/ sbs3
         asp3 = -4._dbl*b*b*bs / sbs3 + 3._dbl*b*b*bs*bs*bs/ sbs5
         w1 = (4._dbl*asp+x*asp2)
         fs(ip,1) = 1._dbl + p*x*x*ys - a1*x/(1+a2*x)
         fs(ip,2) = 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys - a1/(1+a2*x)**2
         fs(ip,3) = 2._dbl*p*ys - p*q*x*w1*ys*ys + &
                    2._dbl*p*q*q*x*x*asp*asp*ys*ys*ys + 2._dbl*a1*a2/(1+a2*x)**3

      s2 = -6*p/(1+q*x*LOG(b*x+SQRT(1+b**2*x**2)))**2*(q*LOG(b*x+SQRT(1+b**2*x**2))+&
           q*x*(b+1/SQRT(1+b**2*x**2)*b**2*x)/(b*x+SQRT(1+b**2*x**2)))+12*p*x/&
           (1+q*x*LOG(b*x+SQRT(1+b**2*x**2)))**3*(q*LOG(b*x+SQRT(1+b**2*x**2))+&
           q*x*(b+1/SQRT(1+b**2*x**2)*b**2*x)/(b*x+SQRT(1+b**2*x**2)))**2
      s1 = s2-6*p*x/(1+q*x*LOG(b*x+SQRT(1+b**2*x**2)))**2*(2*q*(b+1/SQRT(1+b**2*x**2)*b**2*x)/&
           (b*x+SQRT(1+b**2*x**2))+q*x*(-1/SQRT(1+b**2*x**2)**3*b**4*x**2+1/SQRT(1+b**2*x**2)*b**2)/&
           (b*x+SQRT(1+b**2*x**2))-q*x*(b+1/SQRT(1+b**2*x**2)*b**2*x)**2/&
           (b*x+SQRT(1+b**2*x**2))**2)-6*p*x**2/(1+q*x*LOG(b*x+SQRT(1+b**2*x**2)))**4&
           *(q*LOG(b*x+SQRT(1+b**2*x**2))+q*x*(b+1/SQRT(1+b**2*x**2)*b**2*x)/(b*x+SQRT(1+b**2*x**2)))**3
      s2 = s1+6*p*x**2/(1+q*x*LOG(b*x+SQRT(1+b**2*x**2)))**3*(q*LOG(b*x+SQRT(1+b**2*x**2))+&
           q*x*(b+1/SQRT(1+b**2*x**2)*b**2*x)/(b*x+SQRT(1+b**2*x**2)))*(2*q*(b+1/SQRT(1+b**2*x**2)*b**2*x)&
           /(b*x+SQRT(1+b**2*x**2))+q*x*(-1/SQRT(1+b**2*x**2)**3*b**4*x**2+1/SQRT(1+b**2*x**2)*&
           b**2)/(b*x+SQRT(1+b**2*x**2))-q*x*(b+1/SQRT(1+b**2*x**2)*b**2*x)**2/(b*x+SQRT(1+b**2*x**2))**2)
      t0 = s2-p*x**2/(1+q*x*LOG(b*x+SQRT(1+b**2*x**2)))**2*(3*q*(-1/SQRT(1+b**2*x**2)**3*b**4*x**2+&
           1/SQRT(1+b**2*x**2)*b**2)/(b*x+SQRT(1+b**2*x**2))-3*q*(b+1/SQRT(1+b**2*x**2)*b**2*x)**2/&
           (b*x+SQRT(1+b**2*x**2))**2+q*x*(3/SQRT(1+b**2*x**2)**5*b**6*x**3-3/SQRT(1+b**2*x**2)**3*b**4*x)/&
           (b*x+SQRT(1+b**2*x**2))-3*q*x*(-1/SQRT(1+b**2*x**2)**3*b**4*x**2+1/SQRT(1+b**2*x**2)*b**2)/&
           (b*x+SQRT(1+b**2*x**2))**2*(b+1/SQRT(1+b**2*x**2)*b**2*x)+2*q*x*(b+1/SQRT(1+b**2*x**2)*&
           b**2*x)**3/(b*x+SQRT(1+b**2*x**2))**3)-6*a1/(1+a2*x)**3*a2**2+6*a1*x/(1+a2*x)**4*a2**3

         fs(ip,4) = t0
!d       fs(ip,4) = -2._dbl*p*q*asp*ys*ys+p*q*x*(5._dbl*asp2-x*asp3)*ys*ys - &
!d                  p*q*w1*ys*ys + 2._dbl*p*q*q*x*asp*w1*ys*ys*ys - &
!d                  6._dbl*p*q*q*q*x*x*asp*asp*asp*ys*ys*ys*ys + &
!d                  4._dbl*p*q*q*ys*ys*ys*(x*asp*asp+x*x*asp*asp2) - &
!d                  6._dbl * a1*a2*a2/(1+a2*x)**4
      CASE DEFAULT
         CALL stop_program("efactor_t92","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_t92
! *****************************************************************************
  SUBROUTINE efactor_pw91 ( s, fs, m, pset )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m, pset

    INTEGER                                  :: ip
    REAL(dbl)                                :: a, a1, a2, a3, a4, b1, o, &
                                                pa(6,2), s1, s2, s3, x

! parameter set 1: Perdew-Wang
! parameter set 2: Lembarki-Chermette

    pa(1:6,1)=(/0.19645_dbl,    0.2743_dbl,  &
                0.1508_dbl,   100._dbl,      &
                0.8145161_dbl,  0.004_dbl/)
    pa(1:6,2)=(/0.093907_dbl,   0.26608_dbl, &
                0.0809615_dbl,100._dbl,      &
               76.320_dbl,      0.57767e-4_dbl/)
    o=1._dbl
    a1 = pa(1,pset)
    a2 = pa(2,pset)
    a3 = pa(3,pset)
    a4 = pa(4,pset)
    a  = pa(5,pset)
    b1 = pa(6,pset)
    IF ( m >= 0 ) THEN
!$omp parallel do private(ip,x)
      DO ip=1,SIZE(s)
        x = s(ip)
        fs(ip,1) = (o+a1*x*LOG(a*x+SQRT(o+a**2*x**2))+a2-a3*EXP(-a4*x**2)*x**2)&
                   /(o+a1*x*LOG(a*x+SQRT(o+a**2*x**2))+b1*x**4)
      END DO
    END IF
    IF ( m >= 1 ) THEN
!$omp parallel do private(ip,x)
      DO ip=1,SIZE(s)
        x = s(ip)
        fs(ip,2) =&
           (a1*LOG(a*x+SQRT(o+a**2*x**2))+a1*x*a/SQRT(o+a**2*x**2)+2*a3*       &
           a4*x**3*EXP(-a4*x**2)-2*a3*EXP(-a4*x**2)*x)/(o+a1*x*LOG(a*x+SQRT(o  &
           +a**2*x**2))+b1*x**4)-(o+a1*x*LOG(a*x+SQRT(o+a**2*x**2))+a2-a3*exp  &
           (-a4*x**2)*x**2)/(o+a1*x*LOG(a*x+SQRT(o+a**2*x**2))+b1*x**4)**2*(a1 &
           *LOG(a*x+SQRT(o+a**2*x**2))+a1*x*a/SQRT(o+a**2*x**2)+4*b1*x**3)
      END DO
    END IF
    IF ( m >= 2 ) THEN
!$omp parallel do private(ip,x,s1)
      DO ip=1,SIZE(s)
        x = s(ip)
        s1 = (2*a1*a/SQRT(1+a**2*x**2)-a1*x**2*a**3/SQRT(1+a**2*x**2)**3+  &
             10*a3*a4*x**2*EXP(-a4*x**2)-4*a3*a4**2*x**4*EXP(-a4*x**2)     &
             -2*a3*EXP(-a4*x**2))/(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+      &
             b1*x**4)-2*(a1*LOG(a*x+SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2) &
             +2*a3*a4*x**3*EXP(-a4*x**2)-2*a3*EXP(-a4*x**2)*x)/(1+a1*x*    &
             LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**2*(a1*LOG(a*x+SQRT(1+a**2*x**2))&
             +a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)
        fs(ip,3) =&
             s1+2*(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-a3*EXP(-a4*x**2)*x**2)&
             /(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**3*(a1*LOG(a*x+     &
             SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)**2-       &
             (1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-a3*EXP(-a4*x**2)*x**2)     &
             /(1+a1*x*LOG(a*x +SQRT(1+a**2*x**2))+b1*x**4)**2*(2*a1*a/        &
             SQRT(1+a**2*x**2)-a1*x**2*a**3/SQRT(1+a**2*x**2)**3+12*b1*x**2)
      END DO
    END IF
    IF ( m >= 2 ) THEN
!$omp parallel do private(ip,x,s1,s2,s3)
      DO ip=1,SIZE(s)
        x = s(ip)
        s2 = (-4*a1*a**3/SQRT(1+a**2*x**2)**3*x+3*a1*x**3*a**5/         &
             SQRT(1+a**2*x**2)**5+24*a3*a4*x*EXP(-a4*x**2)-36*a3*a4**2* &
             x**3*EXP(-a4*x**2)+8*a3*a4**3*x**5*EXP(-a4*x**2))/(1+a1*x* &
             LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)
        s3 = -3*(2*a1*a/SQRT(1+a**2*x**2)-a1*x**2*a**3/SQRT(1+a**2*     &
             x**2)**3+10*a3*a4*x**2*EXP(-a4*x**2)-4*a3*a4**2*x**4*      &
             EXP(-a4*x**2)-2*a3*EXP(-a4*x**2))/(1+a1*x*LOG(a*x+         &
             SQRT(1+a**2*x**2))+b1*x**4)**2*(a1*LOG(a*x+SQRT(1+a**2*    &
             x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)+6*(a1*LOG(a*x+  &
             SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2)+2*a3*a4*x**3*  &
             EXP(-a4*x**2)-2*a3*EXP(-a4*x**2)*x)/(1+a1*x*LOG(a*x+       &
             SQRT(1+a**2*x**2))+b1*x**4)**3*(a1*LOG(a*x+SQRT(1+a**2*    &
             x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)**2
        s1 = s2+s3
        s2 = s1-3*(a1*LOG(a*x+SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2) &
             +2*a3*a4*x**3*EXP(-a4*x**2)-2*a3*EXP(-a4*x**2)*x)/(1+a1*x*   &
             LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**2*(2*a1*a/SQRT(1+a**2*  &
             x**2)-a1*x**2*a**3/SQRT(1+a**2*x**2)**3+12*b1*x**2)
        s3 = s2-6*(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-a3*EXP(-a4*x**2)* &
             x**2)/(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**4*(a1*    &
             LOG(a*x+SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)**3
        fs(ip,4) =&
             s3+6*(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-a3*EXP(-a4*x**2)  &
             *x**2)/(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**3*(a1*   &
             LOG(a*x+SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*    &
             x**3)*(2*a1*a/SQRT(1+a**2*x**2)-a1*x**2*a**3/SQRT(1+a**2*    &
             x**2)**3+12*b1*x**2)-(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-  &
             a3*EXP(-a4*x**2)*x**2)/(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+   &
             b1*x**4)**2*(-4*a1*a**3/SQRT(1+a**2*x**2)**3*x+3*a1*x**3*    &
             a**5/SQRT(1+a**2*x**2)**5+24*b1*x)
      END DO
    END IF

  END SUBROUTINE efactor_pw91

! *****************************************************************************

END MODULE xc_ke_gga_new

! *****************************************************************************
