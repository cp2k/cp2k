!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate the several different kinetic energy functionals
!>      with a GGA form
!> \par History
!>      JGH (26.02.2003) : OpenMP enabled
!>      fawzi (04.2004)  : adapted to the new xc interface
!> \author JGH (20.02.2002)
! *****************************************************************************
MODULE xc_ke_gga

  USE cp_array_utils_r,                ONLY: cp_3d_r_p_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_functionals_utilities,        ONLY: calc_wave_vector,&
                                             set_util
  USE xc_input_constants,              ONLY: &
       ke_lc, ke_llp, ke_ol1, ke_ol2, ke_pbe, ke_pw86, ke_pw91, ke_revpbe, &
       ke_t92
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  REAL(KIND=dp), PARAMETER :: f13 = 1.0_dp/3.0_dp, &
                              f23 = 2.0_dp*f13, &
                              f43 = 4.0_dp*f13, &
                              f53 = 5.0_dp*f13

  PUBLIC :: ke_gga_info, ke_gga_lda_eval, ke_gga_lsd_eval

  REAL(KIND=dp) :: cf, b, b_lda, b_lsd, flda, flsd, sfac, t13
  REAL(KIND=dp) :: fact, tact
  REAL(KIND=dp) :: eps_rho

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_ke_gga'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param cutoff ...
! *****************************************************************************
  SUBROUTINE ke_gga_init ( cutoff )

    REAL(KIND=dp), INTENT(IN)                :: cutoff

    eps_rho = cutoff
    CALL set_util ( cutoff )

    cf = 0.3_dp*(3.0_dp*pi*pi)**f23
    flda = cf
    flsd = flda * 2.0_dp**f23
!   the_factor 2^(1/3) for LDA is here
    b_lda = 2.0_dp**f43 * (3.0_dp*pi*pi)**(f13)
    b_lsd = 2.0_dp * (3.0_dp*pi*pi)**(f13)
    sfac = 1.0_dp / (2.0_dp*(3.0_dp*pi*pi)**f13)
    t13 = 2.0_dp**f13

  END SUBROUTINE ke_gga_init

! *****************************************************************************
!> \brief ...
!> \param functional ...
!> \param lsd ...
!> \param reference ...
!> \param shortform ...
!> \param needs ...
!> \param max_deriv ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE ke_gga_info ( functional,lsd,reference, shortform, needs, max_deriv, error )
    INTEGER, INTENT(in)                      :: functional
    LOGICAL, INTENT(in)                      :: lsd
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ke_gga_info', &
      routineP = moduleN//':'//routineN

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
       CASE(ke_pbe)
          reference = " "
       CASE(ke_revpbe)
          reference = " "
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
       CASE(ke_pbe)
          shortform = " "
       CASE(ke_revpbe)
          shortform = " "
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
          needs%norm_drho_spin=.TRUE.
       ELSE
          needs%rho=.TRUE.
          needs%rho_1_3=.TRUE.
          needs%norm_drho=.TRUE.
       END IF
    END IF
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE ke_gga_info

! *****************************************************************************
!> \brief ...
!> \param functional ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param order ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE ke_gga_lda_eval ( functional, rho_set, deriv_set, order, error )


    INTEGER, INTENT(IN)                      :: functional
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: order
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ke_gga_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, m, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: drho_cutoff, rho_cutoff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: s
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: fs
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: e_0, e_ndrho, &
      e_ndrho_ndrho, e_ndrho_ndrho_ndrho, e_rho, e_rho_ndrho, &
      e_rho_ndrho_ndrho, e_rho_rho, e_rho_rho_ndrho, e_rho_rho_rho, grho, &
      rho, rho13
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
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

!      s = norm_drho/(rho^(4/3)*2*(pi*pi*3)^(1/3))
       CALL calc_wave_vector ( "p", rho, grho, s )

       fact = flda
!      Definition of s has changed
       b = b_lda
!       tact = t13
       tact =  1.0_dp

       SELECT CASE (functional)
       CASE (ke_ol1)
          CALL efactor_ol1 ( s, fs, m )
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"OL1 functional currently not working properly",&
               error,failure)
       CASE (ke_ol2)
          CALL efactor_ol2 ( s, fs, m )
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"OL2 functional currently not working properly",&
               error,failure)
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
       CASE (ke_pbe)
          CALL efactor_pbex (s, fs, m, 1)
       CASE (ke_revpbe)
          CALL efactor_pbex (s, fs, m, 2)
       CASE DEFAULT
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT

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

          CALL kex_p_1 ( rho, rho13, s, fs, e_rho=e_rho, e_ndrho=e_ndrho,&
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

          CALL kex_p_2 ( rho, rho13, s, fs, e_rho_rho=e_rho_rho,&
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

          CALL kex_p_3 ( rho, rho13, s, fs, e_rho_rho_rho=e_rho_rho_rho,&
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

    CALL timestop(handle)

  END SUBROUTINE ke_gga_lda_eval

! *****************************************************************************
!> \brief ...
!> \param functional ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param order ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE ke_gga_lsd_eval ( functional, rho_set, deriv_set, order, &
       error )


    INTEGER, INTENT(IN)                      :: functional
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN), OPTIONAL            :: order
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ke_gga_lsd_eval', &
      routineP = moduleN//':'//routineN
    CHARACTER(len=12), DIMENSION(2), PARAMETER :: &
      norm_drho_spin_name = (/"(norm_drhoa)","(norm_drhob)"/)
    CHARACTER(len=6), DIMENSION(2), PARAMETER :: &
      rho_spin_name = (/"(rhoa)","(rhob)"/)

    INTEGER                                  :: handle, ispin, m, npoints, &
                                                stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rho_cutoff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: s
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: fs
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: e_0, e_ndrho, &
      e_ndrho_ndrho, e_ndrho_ndrho_ndrho, e_rho, e_rho_ndrho, &
      e_rho_ndrho_ndrho, e_rho_rho, e_rho_rho_ndrho, e_rho_rho_rho
    TYPE(cp_3d_r_p_type), DIMENSION(2)       :: norm_drho, rho, rho_1_3
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
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
       b = b_lsd
       tact = 1.0_dp

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
             tact = (1.0_dp/2.0_dp)**(1.0_dp/3.0_dp)
             CALL efactor_pw86 ( s, fs, m, f2_lsd=tact )
             tact = 1.0_dp
          CASE (ke_pbe)
             tact = (1.0_dp/2.0_dp)**(1.0_dp/3.0_dp)
             CALL efactor_pbex ( s, fs, m, 1, f2_lsd=tact )
             tact = 1.0_dp
          CASE (ke_revpbe)
             tact = (1.0_dp/2.0_dp)**(1.0_dp/3.0_dp)
             CALL efactor_pbex ( s, fs, m, 2, f2_lsd=tact )
             tact = 1.0_dp
          CASE (ke_pw91)
             tact = (1.0_dp/2.0_dp)**(1.0_dp/3.0_dp)
             CALL efactor_pw91 ( s, fs, m, 1, f2_lsd=tact )
             tact = 1.0_dp
          CASE (ke_lc)
             tact = (1.0_dp/2.0_dp)**(1.0_dp/3.0_dp)
             CALL efactor_pw91 ( s, fs, m, 2 , f2_lsd=tact )
             tact = 1.0_dp
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

             CALL kex_p_1 ( rho=rho(ispin)%array, &
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

             CALL kex_p_2 ( rho(ispin)%array, rho_1_3(ispin)%array,&
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

             CALL kex_p_3 ( rho(ispin)%array, &
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
    CALL timestop(handle)
  END SUBROUTINE ke_gga_lsd_eval

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param r13 ...
!> \param fs ...
!> \param e_0 ...
!> \param npoints ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE kex_p_0 ( rho, r13, fs, e_0, npoints, error )

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, r13
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: fs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_0
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip

    !$omp parallel do default(none) &
    !$omp             shared(npoints, rho, e_0, fact, r13, fs, eps_rho) &
    !$omp             private(ip)

    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         e_0(ip) = e_0(ip) + fact*r13(ip)*r13(ip)*rho(ip) * fs(ip,1)
      END IF

    END DO

    !$omp end parallel do

  END SUBROUTINE kex_p_0

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param r13 ...
!> \param s ...
!> \param fs ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param npoints ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE kex_p_1 ( rho, r13, s, fs, e_rho, e_ndrho, npoints, error )

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, r13, s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: fs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_rho, e_ndrho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: a0, a1, sx, sy

    !$omp parallel do default(none) &
    !$omp             shared(npoints, rho, eps_rho, fact, r13, sfac, tact) &
    !$omp             shared(fs, e_rho, e_ndrho, s) &
    !$omp             private(ip,a0,a1,sx,sy)

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

    !$omp end parallel do

  END SUBROUTINE kex_p_1

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param r13 ...
!> \param s ...
!> \param fs ...
!> \param e_rho_rho ...
!> \param e_rho_ndrho ...
!> \param e_ndrho_ndrho ...
!> \param npoints ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE kex_p_2 ( rho, r13, s, fs, e_rho_rho, e_rho_ndrho,e_ndrho_ndrho,&
      npoints,error)

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, r13, s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: fs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_rho_rho, e_rho_ndrho, &
                                                e_ndrho_ndrho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: a0, a1, a2, sx, sxx, sxy, sy

    !$omp parallel do default(none) &
    !$omp             shared (npoints, rho, eps_rho, fact, r13) &
    !$omp             shared (e_rho_rho, e_rho_ndrho, e_ndrho_ndrho, fs) &
    !$omp             shared(s, sfac, tact) &
    !$omp             private(ip,a0,a1,a2,sx,sy,sxx,sxy)

    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*r13(ip)*rho(ip)
         a1 = f53*fact*r13(ip)*r13(ip)
         a2 = f23*f53*fact/r13(ip)
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         sxx= 28.0_dp/9.0_dp*s(ip)/(rho(ip)*rho(ip))
         sxy= -f43*sfac*tact/(r13(ip)*rho(ip)*rho(ip))
         e_rho_rho(ip) = e_rho_rho(ip) + a2*fs(ip,1) +2.0_dp*a1*fs(ip,2)*sx + &
                     a0*fs(ip,3)*sx*sx + a0*fs(ip,2)*sxx
         e_rho_ndrho(ip) = e_rho_ndrho(ip) + a1*fs(ip,2)*sy + a0*fs(ip,3)*sx*sy + &
                     a0*fs(ip,2)*sxy
         e_ndrho_ndrho(ip) = e_ndrho_ndrho(ip) + a0*fs(ip,3)*sy*sy

      END IF

    END DO

    !$omp end parallel do

  END SUBROUTINE kex_p_2

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param r13 ...
!> \param s ...
!> \param fs ...
!> \param e_rho_rho_rho ...
!> \param e_rho_rho_ndrho ...
!> \param e_rho_ndrho_ndrho ...
!> \param e_ndrho_ndrho_ndrho ...
!> \param npoints ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE kex_p_3 ( rho, r13, s, fs, e_rho_rho_rho,e_rho_rho_ndrho,&
      e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho, npoints, error )

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, r13, s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: fs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(inout)                          :: e_rho_rho_rho, &
                                                e_rho_rho_ndrho, &
                                                e_rho_ndrho_ndrho, &
                                                e_ndrho_ndrho_ndrho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: a0, a1, a2, a3, sx, sxx, &
                                                sxxx, sxxy, sxy, sy

    !$omp parallel do default(none) &
    !$omp             shared(npoints, rho, eps_rho, fact, r13) &
    !$omp             shared(s, sfac, tact, fs, e_rho_rho_rho) &
    !$omp             shared(e_rho_rho_ndrho, e_rho_ndrho_ndrho) &
    !$omp             shared(e_ndrho_ndrho_ndrho) &
    !$omp             private(ip,a0,a1,a2,a3,sx,sy,sxx,sxy,sxxx,sxxy)

    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*r13(ip)*rho(ip)
         a1 = f53*fact*r13(ip)*r13(ip)
         a2 = f23*f53*fact/r13(ip)
         a3 = -f13*f23*f53*fact/(r13(ip)*rho(ip))
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         sxx= 28.0_dp/9.0_dp*s(ip)/(rho(ip)*rho(ip))
         sxy= -f43*sfac*tact/(r13(ip)*rho(ip)*rho(ip))
         sxxx= -280.0_dp/27.0_dp*s(ip)/(rho(ip)*rho(ip)*rho(ip))
         sxxy= 28.0_dp/9.0_dp*sfac*tact/(r13(ip)*rho(ip)*rho(ip)*rho(ip))
         e_rho_rho_rho(ip) = e_rho_rho_rho(ip) + a3*fs(ip,1) + 3.0_dp*a2*fs(ip,2)*sx + &
                     3.0_dp*a1*fs(ip,3)*sx*sx + 3.0_dp*a1*fs(ip,2)*sxx + &
                     a0*fs(ip,4)*sx*sx*sx + 3.0_dp*a0*fs(ip,3)*sx*sxx + &
                     a0*fs(ip,2)*sxxx
         e_rho_rho_ndrho(ip) = e_rho_rho_ndrho(ip) + a2*fs(ip,2)*sy + 2.0_dp*a1*fs(ip,3)*sx*sy + &
                     2.0_dp*a1*fs(ip,2)*sxy + a0*fs(ip,4)*sx*sx*sy + &
                     2.0_dp*a0*fs(ip,3)*sx*sxy + a0*fs(ip,3)*sxx*sy + &
                     a0*fs(ip,2)*sxxy
         e_rho_ndrho_ndrho(ip) = e_rho_ndrho_ndrho(ip) + a1*fs(ip,3)*sy*sy + a0*fs(ip,4)*sx*sy*sy + &
                     2.0_dp*a0*fs(ip,3)*sxy*sy
         e_ndrho_ndrho_ndrho(ip) = e_ndrho_ndrho_ndrho(ip) + a0*fs(ip,4)*sy*sy*sy

      END IF

    END DO

    !$omp end parallel do

  END SUBROUTINE kex_p_3

! Enhancement Factors
! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param fs ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE efactor_ol1 ( s, fs, m )
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: fs
    INTEGER, INTENT(IN)                      :: m

    CHARACTER(len=*), PARAMETER :: routineN = 'efactor_ol1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: t1, t2

    t1 = b*b/(72.0_dp*cf)
    t2 = 0.001878_dp*b

    !$omp parallel do default(none) &
    !$omp             shared(s, m, fs, t1, t2) &
    !$omp             private(ip)

    DO ip=1,SIZE(s)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1.0_dp + t1*s(ip)*s(ip) + t2*s(ip)
      CASE (1)
         fs(ip,1) = 1.0_dp + t1*s(ip)*s(ip) + t2*s(ip)
         fs(ip,2) = 2.0_dp*t1*s(ip) + t2
      CASE (2:3)
         fs(ip,1) = 1.0_dp + t1*s(ip)*s(ip) + t2*s(ip)
         fs(ip,2) = 2.0_dp*t1*s(ip) + t2
         fs(ip,3) = 2.0_dp*t1
      CASE DEFAULT
         CALL stop_program(routineN,moduleN,__LINE__,"Illegal order.")
      END SELECT
    END DO

    !$omp end parallel do

    IF(m==3) fs(:,4)=0.0_dp

  END SUBROUTINE efactor_ol1
! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param fs ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE efactor_ol2 ( s, fs, m )
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: fs
    INTEGER, INTENT(IN)                      :: m

    CHARACTER(len=*), PARAMETER :: routineN = 'efactor_ol2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: t1, t2, t3, y

    t1 = b*b/(72.0_dp*cf)
    t2 = 0.0245_dp*b
    t3 = 2.0_dp**f53*b

    !$omp parallel do default(none) &
    !$omp             shared(s, m, t1, t2, t3, fs) &
    !$omp             private(ip,y)

    DO ip=1,SIZE(s)
      y = 1.0_dp/(1.0_dp+t3*s(ip))
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1.0_dp + t1*s(ip)*s(ip) + t2*s(ip)*y
      CASE (1)
         fs(ip,1) = 1.0_dp + t1*s(ip)*s(ip) + t2*s(ip)*y
         fs(ip,2) = 2.0_dp*t1*s(ip) + t2*y*y
      CASE (2)
         fs(ip,1) = 1.0_dp + t1*s(ip)*s(ip) + t2*s(ip)*y
         fs(ip,2) = 2.0_dp*t1*s(ip) + t2*y*y
         fs(ip,3) = 2.0_dp*(t1 - t2*t3*y*y*y)
      CASE (3)
         fs(ip,1) = 1.0_dp + t1*s(ip)*s(ip) + t2*s(ip)*y
         fs(ip,2) = 2.0_dp*t1*s(ip) + t2*y*y
         fs(ip,3) = 2.0_dp*(t1 - t2*t3*y*y*y)
         fs(ip,4) = 6.0_dp*t2*t3*t3*y*y*y*y
      CASE DEFAULT
         CALL stop_program(routineN,moduleN,__LINE__,"Illegal order.")
      END SELECT
    END DO

    !$omp end parallel do

  END SUBROUTINE efactor_ol2
! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param fs ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE efactor_llp ( s, fs, m )
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: fs
    INTEGER, INTENT(IN)                      :: m

    CHARACTER(len=*), PARAMETER :: routineN = 'efactor_llp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(KIND=dp) :: as, bs, p, q, sas, sbs, t1, t10, t11, t12, t133, t16, &
      t17, t19, t2, t20, t22, t23, t26, t28, t29, t3, t30, t33, t36, t39, t4, &
      t40, t42, t43, t45, t46, t47, t49, t5, t50, t54, t55, t7, t71, t8, t9, &
      x, ys

    p = 0.0044188_dp*b*b
    q = 0.0253_dp*b

    !$omp parallel do default(none) &
    !$omp             shared(s, m, fs, p, q, b) &
    !$omp             private(ip,x,bs,sbs,as,sas,ys,t2,t4,t5,t8,t9,t10,t12) &
    !$omp             private(t1,t3,t7,t11,t16,t17,t20,t22,t23,t26,t30,t33) &
    !$omp             private(t40,t42,t43,t45,t46,t47,t49,t50,t71,t133) &
    !$omp             private(t19,t28,t29,t36,t39,t54,t55)

    DO ip=1,SIZE(s)
      x = s(ip)
      bs = b*x
      sbs = SQRT(bs*bs+1.0_dp)
      as = LOG(bs+sbs)
      sas = x*as
      ys =  1.0_dp/(1.0_dp+q*sas)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1.0_dp + p*x*x*ys
      CASE (1)
         fs(ip,1) = 1.0_dp + p*x*x*ys
         t2 = q*x
         t4 = b**2
         t5 = x**2
         t8 = SQRT(1.0_dp+t4*t5)
         t9 = b*x+t8
         t10 = LOG(t9)
         t12 = 1.0_dp+t2*t10
         t17 = t12**2
         fs(ip,2) = 2.0_dp*p*x/t12-p*t5/t17*(q*t10+t2*(b+1.0_dp/t8*t4*x)/t9)

      CASE (2)
         fs(ip,1) = 1.0_dp + p*x*x*ys
         ! first der
         t2 = q*x
         t4 = b**2
         t5 = x**2
         t8 = SQRT(1.0_dp+t4*t5)
         t9 = b*x+t8
         t10 = LOG(t9)
         t12 = 1.0_dp+t2*t10
         t17 = t12**2
         fs(ip,2) = 2.0_dp*p*x/t12-p*t5/t17*(q*t10+t2*(b+1.0_dp/t8*t4*x)/t9)

         ! second der
         t1 = q*x
         t3 = b**2
         t4 = x**2
         t7 = SQRT(1.0_dp+t3*t4)
         t8 = b*x+t7
         t9 = LOG(t8)
         t11 = 1.0_dp+t1*t9
         t16 = t11**2
         t17 = 1.0_dp/t16
         t20 = 1.0_dp/t7*t3
         t22 = b+t20*x
         t23 = 1/t8
         t26 = q*t9+t1*t22*t23
         t30 = p*t4
         t33 = t26**2
         t40 = t7**2
         t43 = t3**2
         t49 = t22**2
         t50 = t8**2
         fs(ip,3) = 2.0_dp*p/t11-4.0_dp*p*x*t17*t26+2.0_dp*t30/t16/&
               t11*t33-t30*t17*(2.0_dp*q*t22*t23+t1*&
              (-1.0_dp/t40/t7*t43*t4+t20)*t23-t1*t49/t50)

      CASE (3)

         fs(ip,1) = 1.0_dp + p*x*x*ys
         ! first der
         t2 = q*x
         t4 = b**2
         t5 = x**2
         t8 = SQRT(1.0_dp+t4*t5)
         t9 = b*x+t8
         t10 = LOG(t9)
         t12 = 1.0_dp+t2*t10
         t17 = t12**2
         fs(ip,2) = 2.0_dp*p*x/t12-p*t5/t17*(q*t10+t2*(b+1.0_dp/t8*t4*x)/t9)

         ! second der
         t1 = q*x
         t3 = b**2
         t4 = x**2
         t7 = SQRT(1.0_dp+t3*t4)
         t8 = b*x+t7
         t9 = LOG(t8)
         t11 = 1.0_dp+t1*t9
         t16 = t11**2
         t17 = 1.0_dp/t16
         t20 = 1.0_dp/t7*t3
         t22 = b+t20*x
         t23 = 1/t8
         t26 = q*t9+t1*t22*t23
         t30 = p*t4
         t33 = t26**2
         t40 = t7**2
         t43 = t3**2
         t49 = t22**2
         t50 = t8**2
         fs(ip,3) = 2.0_dp*p/t11-4.0_dp*p*x*t17*t26+2.0_dp*t30/t16/&
               t11*t33-t30*t17*(2.0_dp*q*t22*t23+t1*&
              (-1.0_dp/t40/t7*t43*t4+t20)*t23-t1*t49/t50)

         t1 = q*x
         t3 = b**2
         t4 = x**2
         t7 = SQRT(1+t3*t4)
         t8 = b*x+t7
         t9 = LOG(t8)
         t11 = 1.0_dp+t1*t9
         t12 = t11**2
         t133 = 1.0_dp/t12
         t17 = 1.0_dp/t7*t3
         t19 = b+t17*x
         t20 = 1.0_dp/t8
         t23 = q*t9+t1*t19*t20
         t26 = p*x
         t28 = 1.0_dp/t12/t11
         t29 = t23**2
         t36 = t7**2
         t39 = t3**2
         t40 = 1.0_dp/t36/t7*t39
         t42 = -t40*t4+t17
         t45 = t19**2
         t46 = t8**2
         t47 = 1.0_dp/t46
         t50 = 2.0_dp*q*t19*t20+t1*t42*t20-t1*t45*t47
         t54 = p*t4
         t55 = t12**2
         t71 = t36**2
         fs(ip,4) = &
             -6.0_dp*p*t133*t23+12.0_dp*t26*t28*t29-&
              6.0_dp*t26*t133*t50-6.0_dp*t54/t55*t29*t23+&
              6.0_dp*t54*t28*t23*t50-t54*t133*&
              (3.0_dp*q*t42*t20-3.0_dp*q*t45*t47+3.0_dp*t1*&
              (1.0_dp/t71/t7*t39*t3*t4*x-t40*x)*t20-&
              3.0_dp*t1*t42*t47*t19+2.0_dp*t1*t45*t19/t46/t8)

      CASE DEFAULT
         CALL stop_program(routineN,moduleN,__LINE__,"Illegal order.")
      END SELECT
    END DO

    !$omp end parallel do

  END SUBROUTINE efactor_llp
! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param fs ...
!> \param m ...
!> \param f2_lsd ...
! *****************************************************************************
  SUBROUTINE efactor_pw86 ( s, fs, m , f2_lsd)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: fs
    INTEGER, INTENT(IN)                      :: m
    REAL(dp), OPTIONAL                       :: f2_lsd

    CHARACTER(len=*), PARAMETER :: routineN = 'efactor_pw86', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: f15, ff, p0, p1, p15, p2, p3, &
                                                s1, s2, s4, s6, t1, t2, t3

    t1 = 1.296_dp
    t2 = 14.0_dp
    t3 = 0.2_dp
    f15 = 1.0_dp/15.0_dp
    ff = 1.0_dp
    IF(PRESENT(f2_lsd)) ff = f2_lsd

    !$omp parallel do default(none) &
    !$omp             shared(s, fs, m, t1, t2, t3, f15, ff) &
    !$omp             private(ip, s1, s2, s4, s6, p0, p1, p2, p3, p15)

    DO ip=1,SIZE(s)
      s1 = s(ip)*ff
      s2 = s1 * s1
      s4 = s2*s2
      s6 = s2*s4
      SELECT CASE (m)
      CASE (0)
         p0 = 1.0_dp + t1*s2 + t2*s4 + t3*s6
         fs(ip,1) = p0**f15
      CASE (1)
         p0 = 1.0_dp + t1*s2 + t2*s4 + t3*s6
         p1 = s1*ff*(2.0_dp*t1 + 4.0_dp*t2*s2 + 6.0_dp*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
      CASE (2)
         p0 = 1.0_dp + t1*s2 + t2*s4 + t3*s6
         p1 = s1*ff*(2.0_dp*t1 + 4.0_dp*t2*s2 + 6.0_dp*t3*s4)
         p2 = ff*ff*(2.0_dp*t1 + 12.0_dp*t2*s2 + 30.0_dp*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
          fs(ip,3) = f15*p15/p0*(p2-14.0_dp/15.0_dp*p1*p1/p0)
      CASE (3)
         p0 = 1.0_dp + t1*s2 + t2*s4 + t3*s6
         p1 =s1*ff*(2.0_dp*t1 + 4.0_dp*t2*s2 + 6.0_dp*t3*s4)
         p2 = ff*ff*(2.0_dp*t1 + 12.0_dp*t2*s2 + 30.0_dp*t3*s4)
         p3 = s1*ff*ff*ff*(24.0_dp*t2+120.0_dp*t3*s2)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
          fs(ip,3) = f15*p15/p0*(p2-14.0_dp/15.0_dp*p1*p1/p0)
          fs(ip,4) = f15*p15/p0*(-14.0_dp*f15*p1*p1/p0+14.0_dp*14.0_dp*f15*p1*p1*p1/p0/p0+&
                     p3-14.0_dp*p2*p1/p0+14.0_dp*p1*p1/p0/p0)
      CASE DEFAULT
         CALL stop_program(routineN,moduleN,__LINE__,"Illegal order.")
      END SELECT
    END DO

    !$omp end parallel do

  END SUBROUTINE efactor_pw86
! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param fs ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE efactor_t92 ( s, fs, m )
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: fs
    INTEGER, INTENT(IN)                      :: m

    CHARACTER(len=*), PARAMETER :: routineN = 'efactor_t92', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: a1, a2, as, asp, asp2, asp3, &
                                                bs, p, q, s1, s2, sas, sbs, &
                                                sbs3, sbs5, t0, w1, x, ys

    p = 0.0055_dp*b*b
    q = 0.0253_dp*b
    a1 = 0.072_dp * b
    a2 = 2.0_dp**f53 * b

    !$omp parallel do default(none) &
    !$omp             shared(s, fs, m, p, q, a1, a2, b) &
    !$omp             private(ip, x, bs, sbs, sas, ys, asp, sbs3, asp2, sbs5) &
    !$omp             private(asp3, as, s2, s1, t0, w1)

    DO ip=1,SIZE(s)
      x = s(ip)
      bs = b*x
      sbs = SQRT(bs*bs+1.0_dp)
      as = LOG(bs+sbs)
      sas = x*as
      ys =  1.0_dp/(1.0_dp+q*sas)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1.0_dp + p*x*x*ys - a1*x/(1+a2*x)
      CASE (1)
         asp = as + bs / sbs
         fs(ip,1) = 1.0_dp + p*x*x*ys - a1*x/(1+a2*x)
         fs(ip,2) = 2.0_dp*p*x*ys - p*q*x*x*asp*ys*ys - a1/(1+a2*x)**2
      CASE (2)
         asp = as + bs / sbs
         sbs3 = sbs*sbs*sbs
         asp2 = 2.0_dp*b / sbs - b*bs*bs/ sbs3
         fs(ip,1) = 1.0_dp + p*x*x*ys - a1*x/(1+a2*x)
         fs(ip,2) = 2.0_dp*p*x*ys - p*q*x*x*asp*ys*ys - a1/(1+a2*x)**2
         fs(ip,3) = 2.0_dp*p*ys - p*q*x*(4.0_dp*asp+x*asp2)*ys*ys + &
                    2.0_dp*p*q*q*x*x*asp*asp*ys*ys*ys + 2.0_dp*a1*a2/(1+a2*x)**3
      CASE (3)
         asp = as + bs / sbs
         sbs3 = sbs*sbs*sbs
         sbs5 = sbs3*sbs*sbs
         asp2 = 2.0_dp*b / sbs - b*bs*bs/ sbs3
         asp3 = -4.0_dp*b*b*bs / sbs3 + 3.0_dp*b*b*bs*bs*bs/ sbs5
         w1 = (4.0_dp*asp+x*asp2)
         fs(ip,1) = 1.0_dp + p*x*x*ys - a1*x/(1+a2*x)
         fs(ip,2) = 2.0_dp*p*x*ys - p*q*x*x*asp*ys*ys - a1/(1+a2*x)**2
         fs(ip,3) = 2.0_dp*p*ys - p*q*x*w1*ys*ys + &
                    2.0_dp*p*q*q*x*x*asp*asp*ys*ys*ys + 2.0_dp*a1*a2/(1+a2*x)**3

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

      CASE DEFAULT
         CALL stop_program(routineN,moduleN,__LINE__,"Illegal order")
      END SELECT
    END DO

    !$omp end parallel do

  END SUBROUTINE efactor_t92
! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param fs ...
!> \param m ...
!> \param pset ...
!> \param f2_lsd ...
! *****************************************************************************
  SUBROUTINE efactor_pbex ( s, fs, m, pset,f2_lsd )

    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: fs
    INTEGER, INTENT(IN)                      :: m, pset
    REAL(dp), OPTIONAL                       :: f2_lsd

    CHARACTER(len=*), PARAMETER :: routineN = 'efactor_pbex', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: kappa1 = 0.804_dp, &
                                                kappa2 = 1.245_dp, &
                                                mu = 0.2195149727645171_dp

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: f0, mk, x, x2, y

    IF( pset==1 ) mk=mu/kappa1
    IF( pset==2 ) mk=mu/kappa2

    f0 = 1.0_dp/tact
    IF(PRESENT(f2_lsd)) f0 = f2_lsd

    !$omp parallel do default(none) &
    !$omp             shared(s, m, fs, f0, mk) &
    !$omp             private(ip,x,x2,y)

    DO ip=1,SIZE(s)
      x = s(ip)*f0
      x2 = x*x
      y = 1.0_dp/(1.0_dp+mk*x2)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1.0_dp + mu*x2*y
      CASE (1)
         fs(ip,1) = 1.0_dp + mu*x2*y
         fs(ip,2) = 2.0_dp*mu*x*y*y*f0
      CASE (2)
         fs(ip,1) = 1.0_dp + mu*x2*y
         fs(ip,2) = 2.0_dp*mu*x*y*y*f0
         fs(ip,3) = -2.0_dp*mu*(3.0_dp*mk*x2-1.0_dp)*y*y*y*f0*f0
      CASE (3)
         fs(ip,1) = 1.0_dp + mu*x2*y
         fs(ip,2) = 2.0_dp*mu*x*y*y*f0
         fs(ip,3) = -2.0_dp*mu*(3.0_dp*mk*x2-1.0_dp)*y*y*y*f0*f0
         fs(ip,4) = 24.0_dp*mu*mk*x*(mk*x2-1.0_dp)*y*y*y*y*f0*f0*f0
      CASE DEFAULT
         CALL stop_program(routineN,moduleN,__LINE__,"Illegal order.")
      END SELECT
    END DO

    !$omp end parallel do

  END SUBROUTINE efactor_pbex

! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param fs ...
!> \param m ...
!> \param pset ...
!> \param f2_lsd ...
! *****************************************************************************
  SUBROUTINE efactor_pw91 ( s, fs, m, pset , f2_lsd)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: fs
    INTEGER, INTENT(IN)                      :: m, pset
    REAL(dp), OPTIONAL                       :: f2_lsd

    INTEGER                                  :: ip
    REAL(dp) :: ff, t1, t10, t101, t106, t109, t111, t113, t119, t12, t123, &
      t124, t13, t14, t15, t16, t17, t18, t19, t2, t20, t21, t22, t23, t25, &
      t26, t27, t28, t29, t3, t30, t31, t33, t35, t37, t38, t39, t4, t40, &
      t44, t47, t48, t5, t50, t51, t53, t55, t56, t57, t58, t59, t6, t60, &
      t64, t65, t69, t7, t70, t71, t73, t77, t78, t8, t80, t82, t9, t90, t93, &
      t94, t96, t98
    REAL(KIND=dp)                            :: a1, a2, a3, a4, a5, bb, o, &
                                                pa(6,2), x

! parameter set 1: Perdew-Wang
! parameter set 2: Lembarki-Chermette

    pa(1:6,1)=(/0.19645_dp,    0.2743_dp,  &
                0.1508_dp,   100.0_dp,      &
                7.7956_dp,  0.004_dp/)
    pa(1:6,2)=(/0.093907_dp,   0.26608_dp, &
                0.0809615_dp,100.0_dp,      &
               76.320_dp,      0.57767e-4_dp/)
    o=1.0_dp
    ff = 1.0_dp
    IF(PRESENT(f2_lsd)) ff = f2_lsd
    a1 = pa(1,pset)*FF
    a2 = pa(2,pset)*ff*ff
    a3 = pa(3,pset)*ff*ff
    a4 = pa(4,pset)*ff*ff
    bb = pa(5,pset)*ff
!   it should be valid also for lsd
    a5 = pa(6,pset)*ff*ff*ff*ff

    !$omp parallel default(none) &
    !$omp          shared(s, fs, m, a1, a2, a3, a4, a5, bb, pa, o, ff) &
    !$omp          private(x,t1,t10,t101,t106,t109,t111,t113,t119,t12,t123) &
    !$omp          private(t124,t13,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22) &
    !$omp          private(t23,t25,t26,t27,t28,t29,t3,t30,t31,t33,t35,t37) &
    !$omp          private(t38,t39,t4,t40,t44,t47,t48,t5,t50,t51,t53,t55) &
    !$omp          private(t56,t57,t58,t59,t6,t60,t64,t65,t69,t7,t70,t71) &
    !$omp          private(t73,t77,t78,t8,t80,t82,t9,t90,t93,t94,t96,t98,ip)

    IF ( m >= 0 ) THEN
      !$omp do
      DO ip=1,SIZE(s)
        x = s(ip)
        t3 = bb**2
        t4 = x**2
        t7 = SQRT(o+t3*t4)
        t9 = LOG(bb*x+t7)
        t10 = a1*x*t9
        t12 = EXP(-a4*t4)
        t17 = t4**2
        fs(ip,1) = (o+t10+(a2-a3*t12)*t4)/(o+t10+a5*t17)
      END DO
      !$omp end do
    END IF
    IF ( m >= 1 ) THEN
      !$omp do
      DO ip=1,SIZE(s)
        x = s(ip)
        t2 = bb**2
        t3 = x**2
        t6 = SQRT(o+t2*t3)
        t7 = bb*x+t6
        t8 = LOG(t7)
        t9 = a1*t8
        t10 = a1*x
        t17 = t10*(bb+1/t6*t2*x)/t7
        t19 = t3*x
        t21 = EXP(-a4*t3)
        t26 = a2-a3*t21
        t30 = t10*t8
        t31 = t3**2
        t33 = o+t30+a5*t31
        t38 = t33**2
        fs(ip,2) =&
           (t9+t17+2._dp*a3*a4*t19*t21+2._dp*t26*x)/&
           t33-(o+t30+t26*t3)/t38*(t9+t17+4._dp*a5*t19)
      END DO
      !$omp end do
    END IF
    IF ( m >= 2 ) THEN
      !$omp do
      DO ip=1,SIZE(s)
        x = s(ip)
        t1 = bb**2
        t2 = x**2
        t5 = SQRT(o+t1*t2)
        t7 = o/t5*t1
        t9 = bb+t7*x
        t12 = bb*x+t5
        t13 = o/t12
        t15 = 2._dp*a1*t9*t13
        t16 = a1*x
        t17 = t5**2
        t20 = t1**2
        t25 = t16*(-o/t17/t5*t20*t2+t7)*t13
        t26 = t9**2
        t27 = t12**2
        t30 = t16*t26/t27
        t31 = a3*a4
        t33 = EXP(-a4*t2)
        t37 = a4**2
        t39 = t2**2
        t44 = a3*t33
        t47 = LOG(t12)
        t48 = t16*t47
        t50 = o+t48+a5*t39
        t53 = a1*t47
        t55 = t16*t9*t13
        t56 = t2*x
        t60 = a2-t44
        t64 = t50**2
        t65 = o/t64
        t69 = t53+t55+4._dp*a5*t56
        t73 = o+t48+t60*t2
        t77 = t69**2
        fs(ip,3) = &
         (t15+t25-t30+10._dp*t31*t2*t33-4._dp*a3*t37*t39*t33+&
          2._dp*a2-2._dp*t44)/t50-2._dp*&
         (t53+t55+2._dp*t31*t56*t33+2._dp*t60*x)*&
         t65*t69+2._dp*t73/t64/t50*t77-t73*t65*(t15+t25-t30+12._dp*a5*t2)
      END DO
      !$omp end do
    END IF
    IF ( m >= 3 ) THEN
      !$omp do
      DO ip=1,SIZE(s)
        x = s(ip)
          t1 = bb**2
          t2 = x**2
          t5 = SQRT(0.1e1_dp+t1*t2)
          t6 = t5**2
          t9 = t1**2
          t10 = 1/t6/t5*t9
          t13 = 1/t5*t1
          t14 = -t10*t2+t13
          t17 = bb*x+t5
          t18 = 1/t17
          t20 = 3*a1*t14*t18
          t22 = bb+t13*x
          t23 = t22**2
          t25 = t17**2
          t26 = 1/t25
          t28 = 3*a1*t23*t26
          t29 = a1*x
          t30 = t6**2
          t35 = t2*x
          t40 = 3*t29*(1/t30/t5*t1*t9*t35-t10*x)*t18
          t44 = 3*t29*t14*t26*t22
          t50 = 2*t29*t23*t22/t25/t17
          t51 = a3*a4
          t53 = EXP(-a4*t2)
          t57 = a4**2
          t58 = a3*t57
          t59 = t35*t53
          t64 = t2**2
          t70 = LOG(t17)
          t71 = t29*t70
          t73 = 0.1e1_dp+t71+a5*t64
          t78 = 2*a1*t22*t18
          t80 = t29*t14*t18
          t82 = t29*t23*t26
          t90 = a3*t53
          t93 = t73**2
          t94 = 1/t93
          t96 = a1*t70
          t98 = t29*t18*t22
          t101 = t96+t98+4*a5*t35
          t106 = a2-t90
          t109 = t96+t98+2*t51*t59+2*t106*x
          t111 = 1/t93/t73
          t113 = t101**2
          t119 = t78+t80-t82+12*a5*t2
          t123 = 0.1e1_dp+t71+t106*t2
          t124 = t93**2
          fs(ip,4) = &
            (t20-t28+t40-t44+t50+24*t51*x*t53-36._dp*t58*t59+8._dp*a3*t57*a4*t64*&
            x*t53)/t73-3._dp*(t78+t80-t82+10._dp*t51*t2*t53-&
            4._dp*t58*t64*t53+2._dp*a2-2._dp*t90)*t94*t101+&
            6._dp*t109*t111*t113-3._dp*t109*t94*t119-6*t123/t124*t113*t101+&
            6._dp*t123*t111*t101*t119-t123*t94*(t20-t28+t40-t44+t50+24._dp*a5*x)
      END DO
      !$omp end do
    END IF

    !$omp end parallel

  END SUBROUTINE efactor_pw91

END MODULE xc_ke_gga

