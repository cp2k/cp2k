!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/wfs [1.0] *
!!
!!   NAME
!!     wfs
!!
!!   FUNCTION
!!     operations with the wave functions
!!
!!   AUTHOR
!!     JGH (12-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE wfs
  
  USE coefficient_types, ONLY : coeff_type, &
       coeff_allocate, coeff_deallocate, PLANEWAVES, SQUARE, PW_COMPLEXDATA1D
  USE coefficients, ONLY : integral
  USE cp_control_types, ONLY : dft_control_type
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE kinetic, ONLY : kinetic_evaluate
  USE message_passing, ONLY : mp_sum
  USE occupations, ONLY : occupations_set, occupations_fsm_set
  USE potential_types, ONLY : potential_type
  USE pws, ONLY : pw_h_psi
  USE termination, ONLY : stop_program
  USE wf_types, ONLY : wf_type, wf_spinset_type, wf_kset_type, &
       wf_allocate_eigenvalues
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: wf_update, wf_randomise, wf_read_restart, wf_write_restart
  
  INTERFACE wf_innerproduct
     MODULE PROCEDURE wf_innerproduct_s1s2
  END INTERFACE
  
CONTAINS

!******************************************************************************

SUBROUTINE wf_update ( vloc, psi, dft_control, wavepar, tprint )
  IMPLICIT NONE
  
! Arguments
  TYPE ( potential_type ), DIMENSION ( : ), INTENT ( IN ) :: vloc
  TYPE ( dft_control_type ), INTENT ( IN ) :: dft_control
  TYPE ( global_environment_type ), INTENT ( IN ) :: wavepar
  LOGICAL, INTENT ( IN ) :: tprint
  
  TYPE ( wf_type ), DIMENSION ( : ), INTENT ( INOUT ) :: psi
  
!------------------------------------------------------------------------------
  
  IF ( wavepar % print_level > 5 ) &
       WRITE ( wavepar % scr, '( A )' ) "Entering wf_update..."
  
  IF ( wavepar % program_name == "WAVE" .AND. &
       dft_control % wave_control % wf_optimisation_method == "SD" ) THEN
     IF ( .NOT. ASSOCIATED ( &
          psi ( 1 ) % psi_s ( 1 ) % psi_ks ( 1 ) % eigenvalue ) ) THEN
        CALL wf_allocate_eigenvalues ( psi )
     END IF
     
     CALL wf_steepest_descent ( vloc, psi, wavepar )
  ELSE
     CALL stop_program ( "wf_update", "no suitable algorithm implemented" )
  END IF
  
  IF ( wavepar % print_level > 5 ) &
       WRITE ( wavepar % scr, '( A )' ) "Entering occupation_update..."
  
  CALL wf_occupation_update ( psi, dft_control, wavepar )
  
END SUBROUTINE wf_update

!******************************************************************************

SUBROUTINE wf_gramschmidt ( psi_ks )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( wf_kset_type ) :: psi_ks
  
! Locals
  INTEGER :: state, state2, nstates
  REAL ( dbl ) :: anorm
  COMPLEX ( dbl ) :: cscv
  
!------------------------------------------------------------------------------
  
  nstates = SIZE ( psi_ks % psi_iks )
  
  DO state = 1, nstates
     
     DO state2 = 1, state - 1
        cscv = integral ( psi_ks % psi_iks ( state2 ), &
             psi_ks % psi_iks ( state ) )
        
! Subtract from the current state: |c_i^> = |c_i> - <c_j|c_i> |c_j>
        IF ( psi_ks % psi_iks ( state ) % in_use == PLANEWAVES ) THEN
           psi_ks % psi_iks ( state ) % pw % cc &
                = psi_ks % psi_iks ( state ) % pw % cc &
                - cscv * psi_ks % psi_iks ( state2 ) % pw % cc
        END IF
     END DO
     
! Normalise
     anorm = integral ( psi_ks % psi_iks ( state ), flag = SQUARE )
     
     anorm = 1._dbl / SQRT ( anorm )
     
     IF ( psi_ks % psi_iks ( state ) % in_use == PLANEWAVES ) THEN
        psi_ks % psi_iks ( state ) % pw % cc &
             = psi_ks % psi_iks ( state ) % pw % cc * anorm
     END IF
     
  END DO
  
END SUBROUTINE wf_gramschmidt

!******************************************************************************

SUBROUTINE wf_randomise ( psi, wavepar )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( IN ) :: wavepar
  
  TYPE ( wf_type ), INTENT ( IN ) :: psi
  
! Locals
  INTEGER :: spin, nspins, kpt, nkpts
  TYPE ( wf_kset_type ), POINTER :: psi_ks
  
!------------------------------------------------------------------------------
  
  nspins = SIZE ( psi % psi_s )
  
  DO spin = 1, nspins
     nkpts = SIZE ( psi % psi_s ( spin ) % psi_ks )
     
     DO kpt = 1, nkpts
        psi_ks => psi % psi_s ( spin ) % psi_ks ( kpt )
        
        CALL wf_kset_randomise ( psi_ks, wavepar )
     END DO
  END DO
  
END SUBROUTINE wf_randomise

!******************************************************************************

SUBROUTINE wf_kset_randomise ( psi_ks, wavepar )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( wf_kset_type ) :: psi_ks
  TYPE ( global_environment_type ), INTENT ( IN ) :: wavepar
  
! Locals
  INTEGER :: state, nstates, gpt, ngpts, state2
  REAL ( dbl ) :: rand1, rand2
  TYPE ( coeff_type ), POINTER :: psi_iks
!  COMPLEX ( dbl ), DIMENSION ( : ) :: get_random_numbers
  COMPLEX ( dbl ) :: innerproduct
  
!------------------------------------------------------------------------------
  
  nstates = SIZE ( psi_ks % psi_iks )
  
  DO state = 1, nstates
     psi_iks => psi_ks % psi_iks ( state )
     IF ( psi_iks % in_use == PLANEWAVES ) THEN
        ngpts = SIZE ( psi_ks % psi_iks ( state ) % pw % cc )
        
        psi_iks % pw % cc = get_random_numbers ( ngpts )
     ELSE
        stop "wf_kset_randomise, no valid datatype"
     END IF
     
  END DO
  
  CALL wf_gramschmidt ( psi_ks )
  
  WRITE ( wavepar % scr, '( A )' ) "Orthonormality check"
  DO state = 1, nstates
     DO state2 = 1, state
        innerproduct = wf_innerproduct ( psi_ks % psi_iks ( state ), &
             psi_ks % psi_iks ( state2 ) )
        IF ( state2 == state ) THEN
           IF ( ABS ( REAL  ( innerproduct ) - 1.0_dbl ) > 1.0E-12_dbl .OR. &
                AIMAG ( innerproduct ) > 1.0E-12_dbl ) THEN
              WRITE ( wavepar % scr, '( "WARNING!!!" )' )
              WRITE ( wavepar % scr, '( A, 2I4, 2F10.6 )' ) &
                   "OVL: ", state, state2, innerproduct
           END IF
        ELSE
           IF ( REAL  ( innerproduct ) > 1.0E-12_dbl .OR. &
                AIMAG ( innerproduct ) > 1.0E-12_dbl ) THEN
              WRITE ( wavepar % scr, '( "WARNING!!!" )' )
              WRITE ( wavepar % scr, '( A, 2I4, 2F10.6 )' ) &
                   "OVL: ", state, state2, innerproduct
           END IF
        END IF
     END DO
  END DO
  WRITE ( wavepar % scr, '( )' )
  
!------------------------------------------------------------------------------
  
CONTAINS
  
  FUNCTION get_random_numbers ( npts ) RESULT ( random_number_array )
    IMPLICIT NONE
    INTEGER, INTENT ( IN ) :: npts
    COMPLEX ( dbl ), DIMENSION ( npts ) :: random_number_array
    REAL ( dbl ), DIMENSION ( npts ) :: rand1, rand2
    
    CALL RANDOM_NUMBER ( rand1 ( : ) )
    CALL RANDOM_NUMBER ( rand2 ( : ) )
    random_number_array = CMPLX ( rand1, rand2, dbl )
  END FUNCTION get_random_numbers
  
END SUBROUTINE wf_kset_randomise

!******************************************************************************

SUBROUTINE wf_h_psi ( psi_iks, xkpt, vloc, h_psi, psi_iks_r )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( IN ) :: psi_iks
  TYPE ( coeff_type ), INTENT ( IN ), OPTIONAL :: psi_iks_r
  TYPE ( coeff_type ), INTENT ( IN ) :: vloc
  REAL ( dbl ), INTENT ( IN ) :: xkpt ( 3 )
  
  TYPE ( coeff_type ), INTENT ( INOUT ) :: h_psi
  
! Locals
  
!------------------------------------------------------------------------------
  
  IF ( psi_iks % in_use == PLANEWAVES ) THEN
! This is so heavily dependent on the details that it is better to
!    separate it
     CALL pw_h_psi ( psi_iks, xkpt, vloc, h_psi, psi_iks_r )
  END IF
  
END SUBROUTINE wf_h_psi

!******************************************************************************

SUBROUTINE wf_steepest_descent ( vloc, psi, wavepar )
  
  USE kpoints, ONLY : kset
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( potential_type ), DIMENSION ( : ), INTENT ( IN ) :: vloc
!  TYPE ( coeff_type ), INTENT ( IN ), OPTIONAL :: psi_iks_r
  TYPE ( global_environment_type ), INTENT ( IN ) :: wavepar
  
  TYPE ( wf_type ), DIMENSION ( : ), INTENT ( INOUT ) :: psi
  
! Locals
  INTEGER :: kpt, nkpt, spin, nspins, wfset, nwfsets, part
  TYPE ( wf_spinset_type ), POINTER :: psi_s
  REAL ( dbl ), DIMENSION ( : ), POINTER :: xkpt
  
!------------------------------------------------------------------------------
  
  part = 1 !TMPTMPTMP
  
  IF ( wavepar % print_level > 5 ) &
       WRITE ( wavepar % scr, '( A )' ) "Entering SD..."
  
  nwfsets = SIZE ( psi )
  
  DO wfset = 1, nwfsets
     nspins = SIZE ( psi ( wfset ) % psi_s )
     
     DO spin = 1, nspins
        psi_s => psi ( wfset ) % psi_s ( spin )
        nkpt = SIZE ( psi_s % psi_ks )
        
        DO kpt = 1, nkpt
           xkpt => kset % xkpt ( :, kpt )
           CALL wf_stateset_steepest_descent ( &
                vloc ( wfset ) % pot_part ( spin, part ), &
                xkpt, psi_s % psi_ks ( kpt ), wavepar )
        END DO
     END DO
  END DO
  
END SUBROUTINE wf_steepest_descent

!******************************************************************************

SUBROUTINE wf_stateset_steepest_descent ( vloc, xkpt, psi_ks, wavepar )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( IN ) :: vloc
  TYPE ( global_environment_type ), INTENT ( IN ) :: wavepar
  REAL ( dbl ), INTENT ( IN ) :: xkpt ( 3 )
  
  TYPE ( wf_kset_type ), INTENT ( INOUT ) :: psi_ks
  
! Locals
  INTEGER :: state, nstates, gpt
  REAL ( dbl ) :: dt_by_e = 0.5_dbl !TMPTMPTMP
  TYPE ( coeff_type ), POINTER :: psi_iks
  TYPE ( coeff_type ) :: h_psi
  
!------------------------------------------------------------------------------
  
  IF ( wavepar % print_level > 5 ) &
       WRITE ( wavepar % scr, '( A )' ) "Entering stateset SD..."
  
  nstates = SIZE ( psi_ks % psi_iks )
  
  CALL coeff_allocate ( h_psi, psi_ks % psi_iks ( 1 ), &
       use_basis = PLANEWAVES )
  
!OpenMP parallel do ?
  DO state = 1, nstates
     psi_iks => psi_ks % psi_iks ( state )
     
     IF ( wavepar % print_level > 7 ) &
          WRITE ( wavepar % scr, '( A )' ) "Going for h|psi>..."
     
     CALL wf_h_psi ( psi_iks, xkpt, vloc, h_psi )
     
     psi_ks % eigenvalue ( state ) = integral ( psi_iks, h_psi )
     
     IF ( wavepar % print_level > 6 ) THEN
        WRITE ( wavepar % scr, '( "eig", I4, F12.6 )' ) state, &
             psi_ks % eigenvalue ( state ) * 27.2116168391_dbl
     END IF
     
     IF ( psi_iks % in_use == PLANEWAVES ) THEN
        CALL wf_precondition ( h_psi )
        DO gpt = 1, SIZE ( psi_iks % pw % cc )
           psi_iks % pw % cc ( gpt ) = psi_iks % pw % cc ( gpt ) &
                - dt_by_e * h_psi % pw % cc ( gpt )
        END DO
     ELSE
        CALL stop_program ( "wf_stateset_steepest_descent", "no suitable field" )
     END IF
  END DO
  
  CALL coeff_deallocate ( h_psi )
  
  CALL wf_gramschmidt ( psi_ks )
  
  CALL wf_subspace_rotate ( psi_ks, vloc, xkpt, wavepar )
  
END SUBROUTINE wf_stateset_steepest_descent

!******************************************************************************

SUBROUTINE wf_precondition ( array )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( INOUT ) :: array
  
! Locals
  INTEGER :: gpt, ngpts
  REAL ( dbl ) :: glen, glen0 = 5.0_dbl !TMPTMPTMP
  
!------------------------------------------------------------------------------
  
  IF ( array % in_use == PLANEWAVES ) THEN
     IF ( array % pw % in_use /= PW_COMPLEXDATA1D ) THEN
        CALL stop_program ( "wf_precondition", "wrong pw part in use" )
     END IF
     
     ngpts = array % pw % pw_grid % ngpts_cut
     DO gpt = 1, ngpts
        glen = SUM ( array % pw % pw_grid % g ( :, gpt ) ** 2 )
        IF ( glen > glen0 ) THEN
           array % pw % cc ( gpt ) = array % pw % cc ( gpt ) * glen0 / glen
        ELSE
           array % pw % cc ( gpt ) = array % pw % cc ( gpt )
        END IF
     END DO
  ELSE
     CALL stop_program ( "wf_precondition", "no basis in use" )
  END IF
  
END SUBROUTINE wf_precondition

!******************************************************************************

FUNCTION wf_innerproduct_s1s2 ( psi_iks1, psi_iks2 ) RESULT ( innerproduct )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ) :: psi_iks1, psi_iks2
  COMPLEX ( dbl ) :: innerproduct
  
! Locals
  
!------------------------------------------------------------------------------
  
  innerproduct = integral ( psi_iks1, psi_iks2 )
  
END FUNCTION wf_innerproduct_s1s2

!******************************************************************************

SUBROUTINE wf_subspace_rotate ( psi_ks, vloc, xkpt, wavepar )
  
  USE eigenvalueproblems, ONLY : diagonalise
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( IN ) :: vloc
  REAL ( dbl ), INTENT ( IN ) :: xkpt ( 3 )
  TYPE ( global_environment_type ), INTENT ( IN ) :: wavepar
  
  TYPE ( wf_kset_type ), INTENT ( INOUT ) :: psi_ks
  
! Locals
  INTEGER :: state, state2, nstates, gpt, ngpts, index
  REAL ( dbl ) :: eigenvalues_k ( SIZE ( psi_ks % psi_iks ) )
  COMPLEX ( dbl ) :: eigenvectors_k ( SIZE ( psi_ks % psi_iks ), &
       SIZE ( psi_ks % psi_iks ) )
  COMPLEX ( dbl ) :: cphp ( SIZE ( psi_ks % psi_iks ) &
       * ( SIZE ( psi_ks % psi_iks ) + 1 ) / 2 )
  TYPE ( coeff_type ), POINTER :: psi_iks
  TYPE ( coeff_type ) :: h_psi
  
  INTEGER :: allocstat
  COMPLEX ( dbl ), DIMENSION ( : ), ALLOCATABLE :: cwork
  
!------------------------------------------------------------------------------
  
  IF ( wavepar % print_level > 5 ) &
       WRITE ( wavepar % scr, '( A )' ) "SUBROT"
  
  nstates = SIZE ( psi_ks % psi_iks )
  
  CALL coeff_allocate ( h_psi, psi_ks % psi_iks ( 1 ), &
       use_basis = PLANEWAVES )
  
  index = 0
  DO state = 1, nstates
     psi_iks => psi_ks % psi_iks ( state )
     
! Calculate H |psi>
     CALL wf_h_psi ( psi_iks, xkpt, vloc, h_psi )
     
! Calculate <psi|H|psi> and store in 'cphp'
     DO state2 = state, nstates
        index = index + 1
        cphp ( index ) = integral ( psi_ks % psi_iks ( state2 ), h_psi )
     END DO
     
  END DO ! i; states
  
  CALL diagonalise ( cphp, nstates, "Lower", eigenvalues_k, eigenvectors_k )
  
  IF ( wavepar % print_level > 3 ) &
       WRITE ( wavepar % scr, '( "EIG: ", 8F9.4 )' ) &
       eigenvalues_k * 27.2116169301_dbl
  
  IF ( psi_ks % psi_iks ( 1 ) % in_use /= PLANEWAVES ) THEN
     CALL stop_program ( "wf_subspace_rotate", &
          "sub-space rotation still plain for plane waves..." )
  END IF
  
!TMPTMPTMP >>>
  ngpts = psi_ks % psi_iks ( 1 ) % pw % pw_grid % ngpts_cut
  
  ALLOCATE ( cwork ( nstates ), STAT = allocstat )
  IF ( allocstat /= 0 ) THEN
     CALL stop_program ( "wf_subspace_rotate", "error allocating cwork" )
  END IF
  
  DO gpt = 1, ngpts
     DO state = 1, nstates
        cwork ( state ) = 0.0_dbl
        DO state2 = 1, nstates
           cwork ( state ) = cwork ( state ) &
                + eigenvectors_k ( state2, state ) &
                * psi_ks % psi_iks ( state2 ) % pw % cc ( gpt )
        END DO
        
     END DO
     
     DO state = 1, nstates
        psi_ks % psi_iks ( state ) % pw % cc ( gpt ) = cwork ( state )
     END DO
  END DO
!TMPTMPTMP <<<
  
  DEALLOCATE ( cwork, STAT = allocstat )
  IF ( allocstat /= 0 ) THEN
     CALL stop_program ( "wf_subspace_rotate", "error deallocating cwork" )
  END IF
  
  CALL coeff_deallocate ( h_psi )
  
END SUBROUTINE wf_subspace_rotate

!******************************************************************************

SUBROUTINE wf_read_restart ( psi, io_unit )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: io_unit
  TYPE ( wf_type ), DIMENSION ( : ), INTENT ( OUT ) :: psi
  
! Locals
  INTEGER :: state, nstates, kpt, nkpts, spin, nspins, wfset, nwfsets
  
!------------------------------------------------------------------------------
  
  nwfsets = SIZE ( psi )
  DO wfset = 1, nwfsets
     nspins = SIZE ( psi ( wfset ) % psi_s )
     DO spin = 1, nspins
        nkpts = SIZE ( psi ( wfset ) % psi_s ( spin ) % psi_ks )
        DO kpt = 1, nkpts
           nstates = SIZE ( psi ( wfset ) % psi_s ( spin ) &
                % psi_ks ( kpt ) % psi_iks )
           DO state = 1, nstates
              READ ( UNIT = io_unit ) &
                   psi ( wfset ) % psi_s ( spin ) % psi_ks ( kpt ) &
                   % psi_iks ( state ) % pw % cc
           END DO
        END DO
     END DO
  END DO
  
END SUBROUTINE wf_read_restart

!******************************************************************************

SUBROUTINE wf_write_restart ( psi, io_unit )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: io_unit
  TYPE ( wf_type ), DIMENSION ( : ), INTENT ( IN ) :: psi
  
! Locals
  INTEGER :: state, nstates, kpt, nkpts, spin, nspins, wfset, nwfsets
  
!------------------------------------------------------------------------------
  
  nwfsets = SIZE ( psi )
  DO wfset = 1, nwfsets
     nspins = SIZE ( psi ( wfset ) % psi_s )
     DO spin = 1, nspins
        nkpts = SIZE ( psi ( wfset ) % psi_s ( spin ) % psi_ks )
        DO kpt = 1, nkpts
           nstates = SIZE ( psi ( wfset ) % psi_s ( spin ) &
                % psi_ks ( kpt ) % psi_iks )
           DO state = 1, nstates
              WRITE ( UNIT = io_unit ) &
                   psi ( wfset ) % psi_s ( spin ) % psi_ks ( kpt ) &
                   % psi_iks ( state ) % pw % cc
           END DO
        END DO
     END DO
  END DO
  
END SUBROUTINE wf_write_restart

!******************************************************************************

SUBROUTINE wf_occupation_update ( psi, dft_control, wavepar )
  
  USE kpoints, ONLY : kset
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( dft_control_type ), INTENT ( IN ) :: dft_control
  TYPE ( global_environment_type ), INTENT ( IN ) :: wavepar
  
  TYPE ( wf_type ), DIMENSION ( : ), INTENT ( INOUT ) :: psi
  
! Locals
  REAL ( dbl ) :: efermi, efermis ( 2 )
  TYPE ( wf_kset_type ), POINTER :: psi_ks
  TYPE ( wf_spinset_type ), POINTER :: psi_s
  
!------------------------------------------------------------------------------
  
  IF ( wavepar % print_level > 5 ) &
       WRITE ( wavepar % scr, '( A )' ) "Entering wf_occupation_update..."
  
  efermi = 0.0_dbl
  efermis ( : ) = 0.0_dbl
  
!TMPTMPTMP
  IF ( .NOT. ASSOCIATED ( psi ( 1 ) % psi_s ( 1 ) % psi_ks ( 1 ) % eigenvalue ) ) STOP
  
  IF ( dft_control % fixedspinmoment_amount /= 0.0_dbl ) THEN
     CALL occupations_fsm_set ( psi, efermis, dft_control )
  ELSE
     CALL occupations_set ( psi, efermi, dft_control )
  END IF
  
#if defined ( FOO )
  IF ( dft_control % fixedspinmoment_amount /= 0.0_dbl ) THEN
!! fixed spin moment
     
     nwfsets = SIZE ( psi )
     wfset_loop: DO wfset = 1, nwfsets
        nspins = SIZE ( psi ( wfset ) % psi_s )
        
        spin_loop: DO spin = 1, nspins
           psi_s => psi ( wfset ) % psi_s ( spin )
           
           IF ( ASSOCIATED ( psi_s % psi_ks ( 1 ) % eigenvalue ) ) THEN
              ...
           ELSE
              nkpts = SIZE ( psi_s % psi_ks )
              nelcs_occupied = 0
              DO kpt = 1, nkpts
                 psi_ks => psi_s % psi_ks ( kpt )
                 
                 psi_ks % occupation ( : ) = 2.0_dbl / REAL ( nspins, dbl )
                 nelcs_occupied = nelcs_occupied &
                      + kset % wkpt ( kpt ) * SUM ( psi_ks % occupation ( : ) )
              END DO
              WRITE ( wavepar % scr, '( A, I4 )' ) "DDD: ", nelcs_occupied
              
! Check that all the states are fully occupied
              IF ( nspins == 1 ) THEN
                 IF ( ABS ( 0.5_dbl * ( nelcs_occupied &
                      - dft_control % nelectrons ) ) > 1.0E-8_dbl ) THEN
                    CALL stop_program ( "wf_occupation_update", &
                         "can only use full occupations without eigenvalues" )
                 END IF
                 
              ELSE IF ( nspins == 2 ) THEN
              IF ( ( spin == 1 .AND. ABS ( 0.5_dbl * &
                   ( nelcs_occupied + dft_control % fixedspinmoment_amount ) &
                   - dft_control % nelectrons ) > 1.0E-8_dbl ) .OR. &
                   ( spin == 2 .AND. ABS ( 0.5_dbl * &
                   ( nelcs_occupied - dft_control % fixedspinmoment_amount ) &
                   - dft_control % nelectrons ) > 1.0E-8_dbl ) ) THEN
                 CALL stop_program ( "wf_occupation_update", &
                      "can only use full occupations without eigenvalues" )
              END IF
           END IF
        END IF
        
     END DO spin_loop
     
  END DO wfset_loop
  
#endif

END SUBROUTINE wf_occupation_update

!******************************************************************************

END MODULE wfs
