!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/pws [1.0] *
!!
!!   NAME
!!     pws
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (4-Jan-2001) : fft_wrap moved to coefficient_types and pw_types
!!     JGH (13-Mar-2001) : new pw_poisson_solver, delete pw_greens_fn
!!
!!   SOURCE
!******************************************************************************

MODULE pws
  USE coefficient_types,               ONLY: SQUARE,&
                                             SQUAREROOT,&
                                             coeff_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE greens_fn,                       ONLY: ANALYTIC0D,&
                                             ANALYTIC1D,&
                                             ANALYTIC2D,&
                                             HOCKNEY0D,&
                                             HOCKNEY1D,&
                                             HOCKNEY2D,&
                                             MT0D,&
                                             MT1D,&
                                             MT2D,&
                                             PERIODIC3D,&
                                             greens_function_type,&
                                             pw_green_fn
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: fourpi
  USE message_passing,                 ONLY: mp_sum
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             PW_MODE_DISTRIBUTED,&
                                             pw_grid_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA1D,&
                                             REALDATA3D,&
                                             RECIPROCALSPACE,&
                                             pw_allocate,&
                                             pw_copy,&
                                             pw_deallocate,&
                                             pw_derive,&
                                             pw_fft_wrap,&
                                             pw_transfer,&
                                             pw_type
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pw_poisson_solver, pw_kinetic_evaluate, init_pw_poisson_solver,&
       cleanup_pw_poisson_solver
  PUBLIC :: pw_integral_aa, pw_integral_ab, pw_integral_a2b, pw_h_psi, &
            pw_collect

! Global variable
  TYPE ( greens_function_type ), POINTER, PRIVATE :: gf
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN='pws'

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pws/init_pw_poisson_solver [1.0] *
!!
!!   NAME
!!     init_pw_poisson_solver
!!
!!   FUNCTION
!!     Initialize the module variable gf, used in the poisson solver
!!     You should call this just before calling the work routine
!!     pw_poisson_solver (don't assume that no other piece of the program
!!     uses the poisson solver in another context).
!!     Call cleanup_pw_poisson_solver when you have finished
!!
!!   AUTHOR
!!     JGH (12-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************
SUBROUTINE init_pw_poisson_solver ( gftype, error )


    TYPE(greens_function_type), INTENT(IN), &
      TARGET                                 :: gftype
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_pw_poisson_solver', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    LOGICAL, SAVE                            :: first_time = .TRUE.

!------------------------------------------------------------------------------

  failure=.FALSE.
  IF (first_time) THEN
     first_time=.FALSE.
  ELSE
     CPPrecondition(.not.ASSOCIATED(gf),cp_failure_level,routineP,error,failure)
  END IF
  gf => gftype

END SUBROUTINE init_pw_poisson_solver
!******************************************************************************

!!****** pws/cleanup_pw_poisson_solver [1.0] *
!!
!!   NAME
!!     cleanup_pw_poisson_solver
!!
!!   FUNCTION
!!     Frees again the module variable gf, used in the poisson solver.
!!     call this as soon as you don't need the poisson solver anymore
!!     (so that other pieces of the program can use it).
!!
!!   AUTHOR
!!     fawzi (11.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************
SUBROUTINE cleanup_pw_poisson_solver ( gftype, error)


    TYPE(greens_function_type), INTENT(IN), &
      OPTIONAL, TARGET                       :: gftype
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cleanup_pw_poisson_solver', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!------------------------------------------------------------------------------

  failure=.FALSE.
  IF (PRESENT(gftype)) THEN
     CPPrecondition(ASSOCIATED(gf),cp_failure_level,routineP,error,failure)
     IF (.not.failure) THEN
        CPPrecondition(gf%id_nr==gftype%id_nr,cp_failure_level,routineP,error,failure)
     END IF
  END IF
  NULLIFY(gf)

END SUBROUTINE cleanup_pw_poisson_solver
!******************************************************************************

!!****** pws/pw_poisson_solver [1.0] *
!!
!!   NAME
!!     pw_poisson_solver
!!
!!   FUNCTION
!!     Solve Poisson equation in a plane wave basis set
!!     Obtains electrostatic potential and its derivatives with respect to r
!!     from the density
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (13-Mar-2001) : completely revised
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_poisson_solver ( density, ehartree, vhartree, dvhartree, h_stress )


    TYPE(coeff_type), INTENT(IN)             :: density
    REAL(dbl), INTENT(OUT), OPTIONAL         :: ehartree
    TYPE(coeff_type), INTENT(INOUT), &
      OPTIONAL, TARGET                       :: vhartree
    TYPE(coeff_type), DIMENSION(3), &
      INTENT(INOUT), OPTIONAL                :: dvhartree
    REAL(dbl), DIMENSION(3, 3), &
      INTENT(OUT), OPTIONAL                  :: h_stress

    INTEGER                                  :: alpha, beta, i, ng
    INTEGER, DIMENSION(3)                    :: n
    REAL(dbl)                                :: ffa
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_type)                            :: dvg( 3 ), rhog

!------------------------------------------------------------------------------

  IF ( .NOT. ASSOCIATED ( gf ) ) THEN
    CALL stop_program ( "pw_poisson_solver", "init_pw_poisson_solver "// &
      "has to be used before pw_poisson_solver" )
  END IF

  NULLIFY (pw_grid)

! point pw
  IF ( PRESENT ( vhartree ) .OR. PRESENT ( dvhartree ) ) &
  pw_grid => vhartree % pw % pw_grid
  IF ( PRESENT ( dvhartree ) .AND. .NOT. PRESENT ( vhartree ) ) &
  CALL stop_program ( 'pw_poisson_solver', 'Need vhartree for dvhartree calculation' )

  IF ( PRESENT ( vhartree ) ) THEN
! density in G space
    CALL pw_allocate ( rhog, pw_grid, COMPLEXDATA1D )
    rhog % in_space = RECIPROCALSPACE

! update the greens function
    CALL pw_green_fn ( gf )
! apply the greens function
    ng = SIZE ( pw_grid % gsq )
    SELECT CASE ( gf % method )
    CASE ( PERIODIC3D, ANALYTIC2D, ANALYTIC1D, ANALYTIC0D, MT2D, MT1D, MT0D )
      CALL pw_transfer ( density % pw, rhog )
      rhog % cc ( 1 : ng )  = rhog % cc ( 1 : ng ) *  &
                              gf % influence_function % cr  ( 1 : ng )
      CALL pw_transfer ( rhog, vhartree % pw )
    CASE ( HOCKNEY2D, HOCKNEY1D, HOCKNEY0D )
       CALL stop_program ( "pw_poisson_solver", " Method not implemented ")
    CASE DEFAULT
       CALL stop_program ( "pw_poisson_solver", " Method not implemented ")
    END SELECT
  ENDIF

  ! do we need to calculate the energy
  IF ( PRESENT ( ehartree ) .AND. PRESENT ( vhartree ) ) THEN
    ehartree = 0.5_dbl * pw_integral_ab ( density % pw, vhartree % pw )
  ELSE  IF ( PRESENT ( ehartree ) .AND. .NOT. PRESENT ( vhartree ) ) THEN
    CALL stop_program ( 'pw_poisson_solver', 'Need vhartree for ehartree calculation' )
  ENDIF

  ! do we need to calculate the derivative of the potential?
  IF ( PRESENT ( h_stress ) .OR. PRESENT ( dvhartree ) ) THEN
    IF ( gf % method /= PERIODIC3D ) THEN
      CALL stop_program ( "pw_poisson_solver", &
        "Pressure calculation works ONLY for fully priodic systems" )
    END IF
    DO i = 1, 3
      CALL pw_allocate ( dvg ( i ), pw_grid, COMPLEXDATA1D )
      dvg ( i ) % in_space = RECIPROCALSPACE
      n = 0
      n ( i ) = 1
      CALL pw_copy ( rhog, dvg ( i ) )
      CALL pw_derive ( dvg ( i ), n )
    END DO
    ! save the derivatives
    IF ( PRESENT ( dvhartree ) ) THEN
      DO i = 1, 3
        CALL pw_transfer ( dvg ( i ), dvhartree ( i ) % pw )
      END DO
    END IF
    ! calculate the contribution to the stress tensor
    ! this is only the contribution from the Greens function
    ! and the volume factor of the plane waves
    IF ( PRESENT ( h_stress ) ) THEN

      ffa = -1.0_dbl / fourpi
      h_stress = 0._dbl
      DO alpha = 1, 3
        h_stress ( alpha, alpha ) = ehartree
        DO beta = alpha, 3
          h_stress ( alpha, beta ) = h_stress ( alpha, beta ) &
             + ffa * pw_integral_ab ( dvg ( alpha ), dvg ( beta ) )
          h_stress ( beta, alpha ) = h_stress ( alpha, beta )
        END DO
      END DO
   
    END IF
    DO i = 1, 3
      CALL pw_deallocate ( dvg ( i ) )
    END DO
  END IF

  IF ( PRESENT ( vhartree ) ) &
  CALL pw_deallocate ( rhog )

END SUBROUTINE pw_poisson_solver

!!*****
!******************************************************************************
!!****** pws/pw_integral [1.0] *
!!
!!   NAME
!!     pw_integral
!!
!!   FUNCTION
!!     Calculate integral over unit cell for functions in plane wave basis
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (14-Mar-2001) : Parallel sum and some tests, HALFSPACE case
!!
!!   SOURCE
!******************************************************************************

FUNCTION pw_integral_ab ( pw1, pw2, flag ) RESULT ( integral_value )

    TYPE(pw_type), INTENT(IN)                :: pw1, pw2
    INTEGER, INTENT(IN), OPTIONAL            :: flag
    REAL(dbl)                                :: integral_value

!------------------------------------------------------------------------------

  IF ( pw1 % pw_grid % identifier /= pw2 % pw_grid % identifier ) THEN
     CALL stop_program ( "integral_ab", "grids incompatible" )
  END IF

  IF ( pw1 % in_use == REALDATA3D .AND. pw2 % in_use == REALDATA3D ) THEN
     integral_value = SUM ( pw1 % cr3d ( :, :, : ) &
          * pw2 % cr3d ( :, :, : ) )
  ELSE IF ( pw1 % in_use == REALDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     integral_value = SUM ( pw1 % cr3d ( :, :, : ) &
          * pw2 % cc3d ( :, :, : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == REALDATA3D ) THEN
     integral_value = SUM ( pw1 % cc3d ( :, :, : ) &
          * pw2 % cr3d ( :, :, : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     integral_value = SUM ( CONJG ( pw1 % cc3d ( :, :, : ) ) &
          * pw2 % cc3d ( :, :, : ) )

  ELSE IF ( pw1 % in_use == REALDATA1D &
       .AND. pw2 % in_use == REALDATA1D ) THEN
     integral_value = SUM ( pw1 % cr ( : ) * pw2 % cr ( : ) )
  ELSE IF ( pw1 % in_use == REALDATA1D &
       .AND. pw2 % in_use == COMPLEXDATA1D ) THEN
     integral_value = SUM ( pw1 % cr ( : ) * pw2 % cc ( : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D &
       .AND. pw2 % in_use == REALDATA1D ) THEN
     integral_value = SUM ( pw1 % cc ( : ) * pw2 % cr ( : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D &
       .AND. pw2 % in_use == COMPLEXDATA1D ) THEN
     integral_value = SUM ( CONJG ( pw1 % cc ( : ) ) * pw2 % cc ( : ) )
  ELSE
     CALL stop_program ( "integral_ab", "no possible DATA" )
  END IF

  IF ( pw1 % in_use == REALDATA3D .OR. pw1 % in_use == COMPLEXDATA3D ) THEN
     integral_value = integral_value * pw1 % pw_grid % dvol
  ELSE
     integral_value = integral_value * pw1 % pw_grid % vol
  ENDIF
  IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
     IF ( pw1 % pw_grid % grid_span == HALFSPACE ) THEN
        integral_value = 2._dbl * integral_value
        IF ( pw1 % pw_grid % have_g0 ) integral_value = integral_value - &
             CONJG ( pw1 % cc ( 1 ) ) * pw2 % cc ( 1 )
     END IF
  END IF

  IF ( pw1 % pw_grid % para % mode == PW_MODE_DISTRIBUTED ) &
     CALL mp_sum ( integral_value, pw1 % pw_grid % para % group )

END FUNCTION pw_integral_ab

!******************************************************************************

FUNCTION pw_integral_aa ( pw1, flag ) RESULT ( integral_value )

    TYPE(pw_type), INTENT(IN)                :: pw1
    INTEGER, INTENT(IN), OPTIONAL            :: flag
    REAL(dbl)                                :: integral_value

!------------------------------------------------------------------------------

  IF ( PRESENT ( flag ) ) THEN
     IF ( flag == SQUARE ) THEN
        IF ( pw1 % in_use == REALDATA3D ) THEN
           integral_value = SUM ( pw1 % cr3d ( :, :, : ) ** 2 )
        ELSE IF ( pw1 % in_use == COMPLEXDATA3D ) THEN
           integral_value = SUM ( CONJG ( pw1 % cc3d ( :, :, : ) ) &
                * pw1 % cc3d ( :, :, : ) )
        ELSE IF ( pw1 % in_use == REALDATA1D ) THEN
           integral_value = SUM ( pw1 % cr ( : ) ** 2 )
        ELSE IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
           integral_value = SUM ( CONJG ( pw1 % cc ( : ) ) &
                * pw1 % cc ( : ) )
           IF ( pw1 % pw_grid % grid_span == HALFSPACE ) THEN
              integral_value = 2._dbl * integral_value
              IF ( pw1 % pw_grid % have_g0 ) integral_value = integral_value - &
                   CONJG ( pw1 % cc ( 1 ) ) * pw1 % cc ( 1 )
           END IF
        ELSE
           CALL stop_program ( "integral_aa", "no possible SQUARE DATA" )
        END IF

     ELSE IF ( flag == SQUAREROOT ) THEN
        CALL stop_program ( "integral_aa", "no SQUAREROOT defined" )
     END IF

  ELSE
     IF ( pw1 % in_use == REALDATA3D ) THEN
        integral_value = SUM ( pw1 % cr3d ( :, :, : ) )
     ELSE IF ( pw1 % in_use == COMPLEXDATA3D ) THEN
        integral_value = SUM ( pw1 % cc3d ( :, :, : ) )
     ELSE IF ( pw1 % in_use == REALDATA1D ) THEN
        integral_value = SUM ( pw1 % cr ( : ) )
     ELSE IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
        integral_value = SUM ( pw1 % cc ( : ) )
     ELSE
        CALL stop_program ( "integral_aa", "no possible DATA" )
     END IF
  END IF

  IF ( pw1 % in_use == REALDATA3D .OR. pw1 % in_use == COMPLEXDATA3D ) THEN
     integral_value = integral_value * pw1 % pw_grid % dvol
  ELSE
     integral_value = integral_value * pw1 % pw_grid % vol
  END IF

  IF ( pw1 % pw_grid % para % mode == PW_MODE_DISTRIBUTED ) &
     CALL mp_sum ( integral_value, pw1 % pw_grid % para % group )

END FUNCTION pw_integral_aa

!******************************************************************************

FUNCTION pw_integral_a2b ( pw1, pw2 ) RESULT ( integral_value )

    TYPE(pw_type), INTENT(IN)                :: pw1, pw2
    REAL(dbl)                                :: integral_value

!------------------------------------------------------------------------------

  IF ( pw1 % pw_grid % identifier /= pw2 % pw_grid % identifier ) THEN
     CALL stop_program ( "integral_a2b", "grids incompatible" )
  END IF
  IF ( pw1 % in_use == REALDATA1D .AND. &
       pw2 % in_use == REALDATA1D ) THEN
     integral_value = SUM ( pw1 % cr ( : ) * pw2 % cr ( : ) &
                            * pw1 % pw_grid % gsq ( : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D .AND. &
            pw2 % in_use == COMPLEXDATA1D ) THEN
     integral_value = SUM ( REAL ( CONJG ( pw1 % cc ( : ) ) &
          * pw2 % cc ( : ), dbl ) * pw1 % pw_grid % gsq ( : ) )
     IF ( pw1 % pw_grid % grid_span == HALFSPACE ) THEN
        integral_value = 2._dbl * integral_value
     END IF
  ELSE
     CALL stop_program ( "integral_a2b", "no possible DATA" )
  END IF

  IF ( pw1 % in_use == REALDATA3D .OR. pw1 % in_use == COMPLEXDATA3D ) THEN
     integral_value = integral_value * pw1 % pw_grid % dvol
  ELSE
     integral_value = integral_value * pw1 % pw_grid % vol
  END IF

  IF ( pw1 % pw_grid % para % mode == PW_MODE_DISTRIBUTED ) &
     CALL mp_sum ( integral_value, pw1 % pw_grid % para % group )

END FUNCTION pw_integral_a2b

!!*****
!******************************************************************************
!!****** pws/pw_kinetic_evaluate [1.0] *
!!
!!   NAME
!!     pw_kinetic_evaluate
!!
!!   FUNCTION
!!     Calculate kinetic energy, potential, derivative of potential and
!!     stress tensor
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_kinetic_evaluate ( psi, xkpt, e_kinetic, kinetic_oper_psi, &
     stress_kinetic )
  

    TYPE(coeff_type), INTENT(IN)             :: psi
    REAL(dbl), INTENT(IN)                    :: xkpt( 3 )
    REAL(dbl), INTENT(OUT), OPTIONAL         :: e_kinetic
    TYPE(coeff_type), INTENT(INOUT), &
      OPTIONAL                               :: kinetic_oper_psi
    REAL(dbl), INTENT(INOUT), OPTIONAL       :: stress_kinetic( 3, 3 )

    INTEGER                                  :: gpt

!------------------------------------------------------------------------------

  IF ( .NOT. PRESENT ( kinetic_oper_psi ) ) THEN
     CALL stop_program ( "pw_kinetic_evaluate", &
          "currently kinetic_oper_psi must be present" )
  END IF
  
  IF ( psi % pw % in_use /= COMPLEXDATA1D .OR. &
       kinetic_oper_psi % pw % in_use /= COMPLEXDATA1D ) THEN
     CALL stop_program ( "pw_kinetic_evaluate", "wrong DATA in USE" )
  END IF

  !OpenMP: parallel do
  DO gpt = 1, psi % pw % pw_grid % ngpts_cut
     kinetic_oper_psi % pw % cc ( gpt ) = kinetic_oper_psi % pw % cc ( gpt ) &
          + 0.5_dbl &
          * SUM ( ( psi % pw % pw_grid % g ( :, gpt ) + xkpt ( : ) ) ** 2 ) &
          * psi % pw % cc ( gpt )
  END DO
  
  IF ( PRESENT ( e_kinetic ) ) &
       e_kinetic = pw_integral_ab ( psi % pw, kinetic_oper_psi % pw )
  
  IF ( PRESENT ( stress_kinetic ) ) THEN
     CALL stop_program ( "pw_kinetic_evaluate", "stress not yet implemented" )
  END IF
  
END SUBROUTINE pw_kinetic_evaluate

!!*****
!******************************************************************************
!!****** pws/pw_h_psi [1.0] *
!!
!!   NAME
!!     pw_h_psi
!!
!!   FUNCTION
!!     Apply local potential to a wavefunction in plane wave basis
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_h_psi ( psi_iks, xkpt, vloc, h_psi, psi_iks_r )


    TYPE(coeff_type), INTENT(IN)             :: psi_iks
    REAL(dbl), INTENT(IN)                    :: xkpt( 3 )
    TYPE(coeff_type), INTENT(IN)             :: vloc
    TYPE(coeff_type), INTENT(INOUT)          :: h_psi
    TYPE(coeff_type), INTENT(IN), OPTIONAL   :: psi_iks_r

    INTEGER                                  :: gpt, gpt2, i1, i2, i3, nhits
    REAL(dbl)                                :: ekinetic
    TYPE(pw_type)                            :: cfft

!TMPTMPTMP
!TMPTMPTMP
!  COMPLEX ( dbl ) :: foo !TMPTMPTMP
!------------------------------------------------------------------------------

  CALL pw_allocate ( cfft, psi_iks % pw % pw_grid, &
       use_data = COMPLEXDATA3D )

#if defined ( DEBUG )
  nhits = 0
  OPEN ( UNIT = 99, FILE = "psik", FORM = "FORMATTED" )
  DO gpt = 1, psi_iks % grid % ngpts
     READ ( UNIT = 99, FMT = * ) i1, i2, i3, foo
     DO gpt2 = 1, psi_iks % grid % ngpts
        IF ( psi_iks % grid % g_hat ( 1, gpt2 ) == i1 .AND. &
             psi_iks % grid % g_hat ( 2, gpt2 ) == i2 .AND. &
             psi_iks % grid % g_hat ( 3, gpt2 ) == i3 ) THEN
           psi_iks % cc ( gpt2 ) = foo
           nhits = nhits + 1
           EXIT
        END IF
     END DO
  END DO
  CLOSE ( UNIT = 99 )
  WRITE ( 6, * ) "DDD: #hits = ", nhits
#endif

  IF ( PRESENT ( psi_iks_r ) ) THEN
     cfft % cr3d = vloc % pw % cr3d * psi_iks_r % pw % cc3d
  ELSE
     CALL pw_fft_wrap ( psi_iks%pw, cfft )

     cfft % cc3d = cfft % cc3d * vloc % pw % cr3d
  END IF

  CALL pw_fft_wrap ( cfft, h_psi%pw )

!! This must be done after the local potential
  CALL pw_kinetic_evaluate ( psi_iks, xkpt, ekinetic, h_psi )

  CALL pw_deallocate ( cfft )

END SUBROUTINE pw_h_psi

!******************************************************************************

SUBROUTINE pw_stateset_density ( psi_iks, n_s )
  
  
    TYPE(coeff_type), DIMENSION(:), &
      INTENT(IN)                             :: psi_iks
    TYPE(coeff_type), INTENT(INOUT)          :: n_s

    INTEGER                                  :: nstates, state
    REAL(dbl)                                :: wkptf
    TYPE(pw_type)                            :: cfft

!------------------------------------------------------------------------------

  CALL pw_allocate ( cfft, psi_iks ( 1 ) % pw % pw_grid, &
       use_data = COMPLEXDATA3D )
  
  nstates = SIZE ( psi_iks )
  
  !OpenMP parallel do ?
  DO state = 1, nstates
     CALL pw_fft_wrap ( psi_iks ( state )%pw , cfft )
     
     wkptf = 1.0_dbl / psi_iks ( state ) % pw % pw_grid % vol ! TMPTMPTMP
     
     n_s % pw % cr3d = n_s % pw % cr3d + wkptf &
          * ( REAL ( cfft % cc3d ) ** 2 + AIMAG ( cfft % cc3d ) ** 2 )
  END DO
  
  CALL pw_deallocate ( cfft )
  
END SUBROUTINE pw_stateset_density

!!*****
!******************************************************************************

SUBROUTINE pw_collect ( pw, label, op )


    TYPE(pw_type), INTENT(IN)                :: pw
    CHARACTER(LEN=*)                         :: label
    CHARACTER(LEN=*), OPTIONAL               :: op

    REAL(dbl)                                :: xsum

!------------------------------------------------------------------------------

  IF ( .NOT. PRESENT ( op ) ) THEN
    IF ( pw % in_use == COMPLEXDATA1D ) THEN
      xsum = ABS ( SUM ( pw % cc ) )
    ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
      xsum = ABS ( SUM ( pw % cc3d ) )
    ELSEIF ( pw % in_use == REALDATA1D ) THEN
      xsum = SUM ( pw % cr )
    ELSEIF ( pw % in_use == REALDATA3D ) THEN
      xsum = SUM ( pw % cr3d )
    END IF
  ELSE
    IF ( pw % in_use == COMPLEXDATA1D ) THEN
      xsum = SUM ( ABS ( pw % cc ) ) 
    ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
      xsum = SUM ( ABS ( pw % cc3d ) )
    ELSEIF ( pw % in_use == REALDATA1D ) THEN
      xsum = SUM ( ABS ( pw % cr ) )
    ELSEIF ( pw % in_use == REALDATA3D ) THEN
      xsum = SUM ( ABS ( pw % cr3d ) )
    END IF
  END IF

  IF ( pw % pw_grid % para % mode == PW_MODE_DISTRIBUTED ) THEN
    CALL mp_sum ( xsum, pw % pw_grid % para % group )
    IF ( pw % pw_grid % para % group_head ) &
      WRITE ( *, '(A,A10,G20.8)' ) " COLLECT ",label,xsum
  ELSE
    WRITE ( *, '(A,A10,G20.8)' ) " COLLECT ",label,xsum
  END IF

END SUBROUTINE pw_collect

!******************************************************************************

END MODULE pws

!******************************************************************************
