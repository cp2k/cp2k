!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!

!!*apsi* Deleted units

MODULE force_fit

  USE amoeba, ONLY : amoeba_evaluate
  USE global_types, ONLY : globenv
! USE io_parameters, ONLY : scr
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type, thermodynamic_type
  USE pair_potential, ONLY : potentialparm_type, spline_nonbond_control
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       cfield, p_error, get_real, get_int, search_label
  USE structure_types, ONLY : structure_type
  USE termination, ONLY : stop_program

!  USE fist_input, ONLY : atom_names
!  USE force_env_types, ONLY : force
!  USE md, ONLY : simulation_parameters
!  USE pair_potential, ONLY : potentialparm, spline_nonbond_control

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fit_control

  TYPE linklist_real_type
     REAL ( dbl ), POINTER :: real_value
     TYPE ( linklist_real_type ), POINTER :: next
  END TYPE linklist_real_type

  TYPE linklist_parameters_type
     TYPE ( linklist_real_type ), POINTER :: lleps
     TYPE ( linklist_real_type ), POINTER :: llsig
     TYPE ( linklist_real_type ), POINTER :: lla
     TYPE ( linklist_real_type ), POINTER :: llb
     TYPE ( linklist_real_type ), POINTER :: llc
  END TYPE linklist_parameters_type

  TYPE fit_parameter_type
     CHARACTER ( len = 10 ) :: type, datatype, restart
     INTEGER :: max_iter, nsteps, iskip, iwrite, idump
     REAL ( dbl ) :: max_tol
  END TYPE fit_parameter_type

  TYPE fit_info_type
     REAL ( dbl ), DIMENSION ( 3 ) :: f_ref, f_fit, r
  END TYPE fit_info_type

! dbg
!  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: v_ref, v_fit
! dbg
!  TYPE ( linklist_parameters ), ALLOCATABLE, DIMENSION ( :, : ) :: param
!  TYPE ( fit_info ), ALLOCATABLE, DIMENSION ( :, : ) :: finfo
!  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: ave_sig, ave_eps
!  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: ave_a, ave_b, ave_c
!  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: y
!  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE, TARGET :: p
!  REAL ( dbl ) :: ftol
!  INTEGER :: iter, iwrite

CONTAINS

!******************************************************************************

!!>---------------------------------------------------------------------------!
!! SECTION: &fit ... &end                                                     !
!!                                                                            !
!!  type:    [LJ, WILL]                  type of fit                          !
!!  restart: [INITIAL, ALL]              restart option                       !
!!  data type: [POS_FORCE, POS_VEL]      type of data to read                 !
!!  out unit: [ATOMIC, KELVIN]           output unit                          !
!!  tol:      max_tol                    maximum tolerance (fractional ) !
!!  iter:     max_iter                   maximum number of iterations         !
!!  nsteps:   nsteps                     number of trajectory steps           !
!!  iskip:    iskip                      analyse every iskip step             !
!!  idump:    idump                      write to dump file every idump steps !
!!                                                                            !
!!<---------------------------------------------------------------------------!

SUBROUTINE read_fit_section ( fit_param )

  IMPLICIT NONE

! Arguments
  TYPE ( fit_parameter_type ), INTENT ( OUT ) :: fit_param

! Locals
  INTEGER :: ierror, ilen, scr
  CHARACTER ( len = 20 ) :: string
  CHARACTER ( len = 5 ) :: label

!------------------------------------------------------------------------------
! this is just a temporary fix to get rid of io_parameters
  scr = 6

!..defaults
  fit_param % max_tol = 1.0E-6_dbl
  fit_param % max_iter = 10000
  fit_param % nsteps = 1
  fit_param % iskip = 1
  fit_param % type = 'LJ'
  fit_param % datatype = 'POS_VEL'
  fit_param % idump = 10
  fit_param % restart = 'INITIAL'

!..parse the input section


  label = '&FIT'

  CALL parser_init ( globenv % input_file_name )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )

  IF ( ierror /= 0 ) THEN

     IF ( globenv % ionode ) THEN
        CALL stop_program ( 'read_fit_section', 'no input section &FIT found' )
     END IF

  ELSE

     CALL read_line()

     DO WHILE ( test_next ( ) /= 'X' )

        ilen = 6
        CALL cfield ( string, ilen )

        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error ( )
           CALL stop_program ( 'read_fit_section', 'unknown option' )

        CASE ( 'TYPE' )
           ilen = 6
           CALL cfield ( string, ilen )
           fit_param % type = string ( 1:ilen )

        CASE ( 'DATA' )
           ilen = 10
           CALL cfield ( string, ilen )
           fit_param % datatype = string ( 1:ilen )

        CASE ( 'RESTAR' )
           ilen = 7
           CALL cfield ( string, ilen )
           fit_param % restart = string ( 1:ilen )

        CASE ( 'TOL' )
           fit_param % max_tol = get_real ( )

        CASE ( 'ITER' )
           fit_param % max_iter = get_int ( )

        CASE ( 'NSTEPS' )
           fit_param % nsteps = get_int ( )

        CASE ( 'IDUMP' )
           fit_param % idump = get_int ( )

        CASE ( 'ISKIP' )
           fit_param % iskip = get_int ( )

        END SELECT

! Get the next line
        CALL read_line

     END DO

  END IF
  CALL parser_end

!..end of parsing the input section
!..write some information to output
  IF ( globenv % print_level >= 0 ) THEN
     WRITE ( scr, '( A )' ) ' FIT| Fit Protocol '
     WRITE ( scr, '( A, T71, A )' ) &
          ' FIT| Fit type ', ADJUSTR ( fit_param % type )
     WRITE ( scr, '( A, T71, A )' ) &
          ' FIT| Data type ', ADJUSTR ( fit_param % datatype )
     WRITE ( scr, '( A, T71, A )' ) &
          ' FIT| Restart type', ADJUSTR ( fit_param % restart )
     WRITE ( scr, '( A, T71, I10 )' ) &
          ' FIT| Number of steps ', fit_param % nsteps
     WRITE ( scr, '( A, T71, I10 )' ) &
          ' FIT| Iskip ', fit_param % iskip
     WRITE ( scr, '( A, T71, I10 )' ) &
          ' FIT| Idump ', fit_param % idump
     WRITE ( scr, '( A, T71, G10.4 )' ) &
          ' FIT| Max tolerance ', fit_param % max_tol
     WRITE ( scr, '( A, T71, I10 )' ) &
          ' FIT| Max iterations ', fit_param % max_iter
  END IF

END SUBROUTINE read_fit_section

!******************************************************************************

SUBROUTINE fit_control ( struc, potparm, thermo, simpar )

  IMPLICIT NONE

! Arguments
  TYPE ( structure_type ), INTENT ( IN ) :: struc
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), INTENT ( IN ) :: potparm
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( thermodynamic_type ), INTENT ( IN ) :: thermo

! Locals
  LOGICAL :: conv_flag
  INTEGER :: isos
  TYPE ( fit_parameter_type ) :: fit_param

!------------------------------------------------------------------------------

  ALLOCATE ( param ( SIZE ( potparm, 1 ), SIZE ( potparm, 2 ) ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_program ( 'fit_control',&
                                       'failed to allocate param' )

  CALL read_fit_section ( fit_param )

! allocating the arrays for the analysis
  SELECT CASE ( fit_param % type )
  CASE ( 'LJ' )
     ALLOCATE ( ave_sig ( SIZE ( potparm, 1 ), SIZE ( potparm, 2 )  ), &
          STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'fit_control', &
          'failed to allocate LJ ave_sig' )
     ALLOCATE ( ave_eps ( SIZE ( potparm, 1 ), SIZE ( potparm, 2 )  ), &
          STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'fit_control', &
          'failed to allocate LJ ave_eps' )

  CASE ( 'WILL' )
     ALLOCATE ( ave_a ( SIZE ( potparm, 1 ), SIZE ( potparm, 2 )  ), &
          STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'fit_control', &
          'failed to allocate WILL ave_a' )
     ALLOCATE ( ave_b ( SIZE ( potparm, 1 ), SIZE ( potparm, 2 )  ), &
          STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'fit_control', &
          'failed to allocate WILL ave_b' )
     ALLOCATE ( ave_c ( SIZE ( potparm, 1 ), SIZE ( potparm, 2 )  ), &
          STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'fit_control', &
          'failed to allocate WILL ave_c' )

  CASE DEFAULT
     CALL stop_program ( 'fit_control', 'no valid potential type' )
  END SELECT

  conv_flag = .FALSE.

  CALL fit ( conv_flag, fit_param, struc, simpar )

  IF ( conv_flag ) THEN
     WRITE ( scr, '( )' )
     WRITE ( scr, '( T23, A )' ) '******************************'
     WRITE ( scr, '( T23, A )' ) '*                            *'
     WRITE ( scr, '( T23, A )' ) '*    CONVERGENCE ACHIEVED    *'
     WRITE ( scr, '( T23, A )' ) '*                            *'
     WRITE ( scr, '( T23, A )' ) '******************************'
     WRITE ( scr, '( )' )
  ELSE
     WRITE ( scr, '( )' )
     WRITE ( scr, '( T23, A )' ) '******************************'
     WRITE ( scr, '( T23, A )' ) '*                            *'
     WRITE ( scr, '( T23, A )' ) '*  WARNING: NO CONVERGENCE   *'
     WRITE ( scr, '( T23, A )' ) '*                            *'
     WRITE ( scr, '( T23, A )' ) '******************************'
     WRITE ( scr, '( )' )
  END IF

END SUBROUTINE fit_control

!******************************************************************************

!! Given the reference force, and the derived force, optimize
!!    paramters with a least squares fit.

SUBROUTINE fit ( conv_flag, fit_param, struc, simpar )

  IMPLICIT NONE

! Arguments
  LOGICAL, INTENT ( INOUT ) :: conv_flag
  TYPE ( fit_parameter_type ), INTENT ( INOUT ) :: fit_param
  TYPE ( structure_type ), INTENT ( IN ) :: struc
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar

! Locals
  INTEGER :: i, j
  REAL ( dbl ) :: tol
  REAL ( dbl ), DIMENSION ( : ), POINTER :: y
  REAL ( dbl ), DIMENSION ( :, : ), POINTER :: p

  INTERFACE
     FUNCTION func ( x )
       USE kinds, ONLY : dbl
       IMPLICIT NONE
       REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: x
       REAL ( dbl ) :: func
     END FUNCTION func
  END INTERFACE

!------------------------------------------------------------------------------

! Allocate and intialize arrays for amoeba
  CALL copy_params ( fit_param, potparm )

! Get reference force
  CALL force_ref ( y, p, fit_param, struc, simpar )

! initialize y
  SELECT CASE ( fit_param % restart )
  CASE ( 'ALL' )

     CALL fit_read_restart ( y, p )

  CASE ( 'INITIAL' )

     DO i = 1, SIZE ( y )
        y ( i ) = func ( p ( i, : )  )
     END DO

  CASE DEFAULT
     CALL stop_program ( "fit", "no valid restart option" )

  END SELECT

! optimize
  DO iter = 1, fit_param % max_iter, fit_param % iwrite

     CALL amoeba ( p, y, fit_param % max_tol, tol, func, fit_param % iwrite )

     SELECT CASE ( fit_param % type )
     CASE ( 'LJ' )
        CALL get_lj_params ( potparm )
     CASE ( 'WILL' )
        CALL get_will_params ( potparm )
     CASE DEFAULT
        CALL stop_program ( "fit", "no valid potential option" )
     END SELECT

     CALL fit_print_info ( fit_param, tol, potparm )

! write restart file
     IF ( MOD ( iter, fit_param % idump ) == 0 ) &
          CALL fit_write_restart ( y, p )

! convergence achieved
     IF ( tol <= fit_param % max_tol ) THEN
        conv_flag = .TRUE.

        CALL fit_write_restart ( y, p )

        EXIT
     END IF

  END DO

END SUBROUTINE fit

!******************************************************************************

SUBROUTINE copy_params ( fit_param, potparm, x )

  IMPLICIT NONE

! Arguments
  TYPE ( fit_parameter_type ), INTENT ( IN ) :: fit_param
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), INTENT ( INOUT ) :: potparm
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ), OPTIONAL :: x

! Locals

!------------------------------------------------------------------------------

  SELECT CASE ( fit_param % type )
  CASE ( 'LJ' )
     CALL lj_fit ( potparm, x )
  CASE ( 'WILL' )
     CALL will_fit ( potparm, x )
  CASE ( 'TB' )
     CALL stop_program ( "copy_params", "case TB not implemented" )
  CASE DEFAULT
     CALL stop_program ( "copy_params", "no valid potential type" )
  END SELECT

END SUBROUTINE copy_params

!******************************************************************************

SUBROUTINE lj_fit ( potparm, x, y, p )

  IMPLICIT NONE

! Arguments
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), INTENT ( INOUT ) :: potparm
  REAL ( dbl ), DIMENSION ( : ), POINTER, OPTIONAL :: y
  REAL ( dbl ), DIMENSION ( :, : ), POINTER, OPTIONAL :: p
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ), OPTIONAL :: x
  TYPE ( linklist_real_type ), POINTER :: llsig
  TYPE ( linklist_real_type ), POINTER :: lleps

! Locals
  INTEGER :: i, icount, iatom, jatom, numpar, isos
  REAL ( dbl ) :: sig, eps

!------------------------------------------------------------------------------

  IF ( PRESENT ( x )  ) THEN
     icount = 0
     DO iatom = 1, SIZE ( potparm, 1 )
        DO jatom = iatom, SIZE ( potparm, 2 )
           icount = icount + 1
           potparm ( iatom, jatom ) % lj % epsilon = x ( icount ) ** 2

! symmetrizing
           potparm ( jatom, iatom ) % lj % epsilon &
                = potparm ( iatom, jatom ) % lj % epsilon
           icount = icount + 1
           potparm ( iatom, jatom ) % lj % sigma6 = x ( icount )
           potparm ( iatom, jatom ) % lj % sigma12 = x ( icount ) ** 4

! symmetrizing
           potparm ( jatom, iatom ) % lj % sigma6 &
                = potparm ( iatom, jatom ) % lj % sigma6
           potparm ( jatom, iatom ) % lj % sigma12 &
                = potparm ( iatom, jatom ) % lj % sigma12
        END DO
     END DO

! initialize splines
     CALL spline_nonbond_control ( potparm, pstat, 2500, cff % ewald_type )

  ELSE

! number of potential parameters to be optimized is 2*number
! of interactions ( remember: we only need the upper triangle )
! e.g. for LJ each interaction has 2 parameters, sigma and epsilon
     numpar = 2*SIZE ( potparm, 1 ) * ( SIZE ( potparm, 1 ) +1 ) / 2
     ALLOCATE ( y ( numpar+1 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'lj_fit', 'failed to allocate y' )
     ALLOCATE ( p ( numpar+1, numpar ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'lj_fit', 'failed to allocate p' )

! i = 0, numpar is the numpar+1 vector in p ( i, j ) .  It is also
! used to generate the initial guesses for the verticies
! of the amoeba
     ALLOCATE ( lleps, STAT = isos )
     ALLOCATE ( llsig, STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'lj_fit', &
          'failed to allocate lleps and llsig' )

     icount = 0
     DO iatom = 1, SIZE ( potparm, 1 )
        DO jatom = iatom, SIZE ( potparm, 2 )
           param ( iatom, jatom ) % lleps => lleps
           param ( iatom, jatom ) % llsig => llsig

           icount = icount + 1

! loop for eps
           DO i = 1, numpar + 1
              eps = potparm ( iatom, jatom ) % lj % epsilon
              p ( i, icount ) = SQRT ( eps )
              lleps % eps => p ( i, icount )
              ALLOCATE ( lleps % next, STAT = isos )
              IF ( isos /= 0 ) CALL stop_program ( 'lj_fit', &
                   'failed to allocate lleps % next' )
              lleps => lleps % next
           END DO

           icount = icount + 1

! loop for sig
           DO i = 1, numpar + 1
              sig = potparm ( iatom, jatom ) % lj % sigma6
              p ( i, icount ) = SQRT ( sig )
              llsig % sig => p ( i, icount )
              ALLOCATE ( llsig % next, STAT = isos )
              IF ( isos /= 0 ) CALL stop_program ( 'lj_fit', &
                   'failed to allocate llsig % next' )
              llsig => llsig % next
           END DO
        END DO
     END DO

! scaling the appropriate interaction parameter by 1.1
! to construct the amoeba
     DO i = 0, numpar
        IF ( i>0 ) p ( i + 1, i ) = p ( i + 1, i ) * 1.1_dbl
     END DO

  END IF

END SUBROUTINE lj_fit

!******************************************************************************

SUBROUTINE will_fit ( potparm, x, y, p )

  IMPLICIT NONE

! Arguments
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), INTENT ( INOUT ) :: potparm
  REAL ( dbl ), DIMENSION ( : ), POINTER, OPTIONAL :: y
  REAL ( dbl ), DIMENSION ( :, : ), POINTER, OPTIONAL :: p
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ), OPTIONAL :: x
  TYPE ( linklist_a ), POINTER :: lla
  TYPE ( linklist_b ), POINTER :: llb
  TYPE ( linklist_c ), POINTER :: llc

! Locals
  INTEGER :: i, icount, iatom, jatom, numpar, isos
  REAL ( dbl ) a, b, c

!------------------------------------------------------------------------------

  IF ( PRESENT ( x )  ) THEN
     icount = 0
     DO iatom = 1, SIZE ( potparm, 1 )
        DO jatom = iatom, SIZE ( potparm, 2 )
           icount = icount + 1
           potparm ( iatom, jatom ) % willis % a &
                = x ( icount ) ** 2
! symmetrizing
           potparm ( jatom, iatom ) % willis % a &
                = potparm ( iatom, jatom ) % willis % a
           icount = icount + 1
           potparm ( iatom, jatom ) % willis % b &
                = x ( icount ) ** 2
! symmetrizing
           potparm ( jatom, iatom ) % willis % b &
                = potparm ( iatom, jatom ) % willis % b
           icount = icount + 1
           potparm ( iatom, jatom ) % willis % c &
                = x ( icount ) ** 2
! symmetrizing
           potparm ( jatom, iatom ) % willis % c &
                = potparm ( iatom, jatom ) % willis % c
        END DO
     END DO

! initialize splines
     CALL spline_nonbond_control ( potparm, pstat, 2500, cff % ewald_type )

  ELSE

! number of potential parameters to be optimized is 3*number
! of interactions ( remember: we only need the upper triangle )
! e.g. for WILL each interaction has 3 parameters, a, b, and c.
     numpar = 3.0_dbl &
          * SIZE ( potparm, 1 ) * ( SIZE ( potparm, 1 ) + 1 ) / 2

     ALLOCATE ( y ( numpar+1 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'will_fit', 'failed to allocate y' )
     ALLOCATE ( p ( numpar+1, numpar ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'will_fit', 'failed to allocate p' )

! i = 0, numpar is the numpar+1 vector in p ( i, j ) .  It is also
! used to generate the initial guesses for the verticies
! of the amoeba
     ALLOCATE ( lla, STAT = isos )
     ALLOCATE ( llb, STAT = isos )
     ALLOCATE ( llc, STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'will_fit', 'failed to allocate lla' )

     icount = 0
     DO iatom = 1, SIZE ( potparm, 1 )
        DO jatom = iatom, SIZE ( potparm, 2 )
           param ( iatom, jatom ) % lla => lla
           param ( iatom, jatom ) % llb => llb
           param ( iatom, jatom ) % llc => llc

           icount = icount + 1

! loop for a
           DO i = 1, numpar + 1
              a = potparm ( iatom, jatom ) % willis % a
              p ( i, icount ) = SQRT ( a )
              lla % a => p ( i, icount )

              ALLOCATE ( lla % next, STAT = isos )
              IF ( isos /= 0 ) CALL stop_program ( 'will_fit', &
                   'failed to allocate lla % next' )
              lla => lla % next
           END DO

           icount = icount + 1

! loop for b
           DO i = 1, numpar + 1
              b = potparm ( iatom, jatom ) % willis % b
              p ( i, icount ) = SQRT ( b )
              llb % b => p ( i, icount )
              ALLOCATE ( llb % next, STAT = isos )
              IF ( isos /= 0 ) CALL stop_program ( 'will_fit', &
                   'failed to allocate llb % next' )
              llb => llb % next
           END DO

           icount = icount + 1

! loop for c
           DO i = 1, numpar + 1
              c = potparm ( iatom, jatom ) % willis % c
              p ( i, icount ) = SQRT ( c )
              llc % c => p ( i, icount )
              ALLOCATE ( llc % next, STAT = isos )
              IF ( isos /= 0 ) CALL stop_program ( 'will_fit', &
                   'failed to allocate llc % next' )
              llc => llc % next
           END DO
        END DO
     END DO

! scaling the appropriate interaction parameter by 1.1
!    to construct the amoeba
     DO i = 1, numpar
        p ( i + 1, i ) = p ( i + 1, i ) * 1.1_dbl
     END DO

  END IF

END SUBROUTINE will_fit

!******************************************************************************

!! From the first principles trajectory file obtain the total force
!!    on the particle by inverting the velcocity Verlet equations

SUBROUTINE force_ref ( y, p, fit_param, struc, simpar )

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), DIMENSION ( : ), INTENT ( OUT ) :: y
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ) :: p
  TYPE ( fit_parameter_type ), INTENT ( INOUT ) :: fit_param
  TYPE ( structure_type ), INTENT ( IN ) :: struc
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar

! Locals
  INTEGER :: npart, ipart, itime, ii, it, isos
  REAL ( dbl ), ALLOCATABLE, DIMENSION ( :, :, : ) :: r, v
  REAL ( dbl ), DIMENSION ( 3 ) :: fref
  LOGICAL :: copy
! dbg
  REAL ( dbl ) :: dummy
! dbg

!------------------------------------------------------------------------------

  npart = SIZE ( struc % pnode )
  itime = fit_param % nsteps / fit_param % iskip

  IF ( fit_param % datatype == 'POS_FORCE' ) THEN
     ALLOCATE ( finfo ( itime, npart ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'force_ref',&
          'failed to allocate finfo' )
! dbg
     ALLOCATE ( v_ref ( itime ), STAT = isos )
     ALLOCATE ( v_fit ( itime ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'force_ref',&
          'failed to allocate v_ref' )
! dbg

  ELSE IF ( fit_param % datatype == 'POS_VEL' ) THEN

! Can only compute force up to nsteps-1 step: thus allocated so
     IF ( MOD ( fit_param % nsteps, itime ) == 0 ) THEN
        ALLOCATE ( finfo ( itime-1, npart ), STAT = isos )
        IF ( isos /= 0 ) CALL stop_program ( 'force_ref', &
             'failed to allocate finfo' )
     ELSE
        ALLOCATE ( finfo ( itime, npart ), STAT = isos )
        IF ( isos /= 0 ) CALL stop_program ( 'force_ref', &
             'failed to allocate finfo' )
     END IF
  END IF

  SELECT CASE ( fit_param % datatype )
! FOR DATA TYPE FORCES
  CASE ( 'POS_FORCE' )
     ALLOCATE ( r ( 1, 1, 3 ), STAT = isos )
     ALLOCATE ( v ( 1, 1, 3 ), STAT = isos ) !  reads force data
     IF ( isos /= 0 ) CALL stop_program ( 'force_ref', &
          'failed to allocate r and v' )

     copy = .FALSE.
     it = 0
     DO itime = 1, fit_param % nsteps
        IF ( MOD ( itime, fit_param % iskip ) == 0 ) THEN
           copy = .TRUE.
           it = it + 1
        END IF

        DO ipart = 1, npart
           READ ( 250, FMT = * ) ii, r ( 1, 1, : ), v ( 1, 1, : )
           IF ( copy ) THEN
              finfo ( it, ipart ) % r ( : ) = r ( 1, 1, : )
              finfo ( it, ipart ) % f_ref ( : ) = v ( 1, 1, : )
           END IF
        END DO
        READ ( 250, FMT = * ) ii, dummy
! dbg
!           IF ( copy ) v_ref ( it ) = dummy
! dbg
        copy = .FALSE.
     END DO

! Reassign fit_param % nsteps = it
     fit_param % nsteps = it

     DEALLOCATE ( r, STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'force_ref',&
          'failed to deallocate r' )
     DEALLOCATE ( v, STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'force_ref',&
          'failed to deallocate v' )

! FOR DATA TYPE TRAJECTORY
  CASE ( 'POS_VEL' )
     ALLOCATE ( r ( 2, npart, 3 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'force_ref',&
          'failed to allocate r' )
     ALLOCATE ( v ( 2, npart, 3 ), STAT = isos ) ! reads velocity data
     IF ( isos /= 0 ) CALL stop_program ( 'force_ref', 'failed to allocate v' )

! Read in particles from t = 1
     DO ipart = 1, npart
        READ ( 200, FMT = * ) ii, r ( 1, ipart, : ), v ( 1, ipart, : )
     END DO

! initialize counter and logical
     copy = .FALSE.

     it = 1
     DO itime = 2, fit_param % nsteps

! set flag
        IF ( MOD ( itime-1, fit_param % iskip ) == 0 ) copy = .TRUE.

        DO ipart = 1, npart
           READ ( 200, FMT = * ) ii, r ( 2, ipart, : ), v ( 2, ipart, : )

! Invert velocity verlet equations
           fref ( : ) = 2.0_dbl &
                * struc % pnode ( ipart ) % p % prop % mass &
                * ( r ( 2, ipart, : ) -r ( 1, ipart, : ) &
                - simpar % dt * v ( 1, ipart, : ) ) &
                / ( simpar % dt ** 2 )

! reassigning r ( 1 ) to r ( 2 ) and v ( 1 ) to v ( 2 )
!    for next force calculation
           IF ( copy ) THEN
              finfo ( it, ipart ) % r ( : ) = r ( 1, ipart, : )
              finfo ( it, ipart ) % f_ref ( : ) = fref ( : )
           END IF
        END DO

! increment the counter
        IF (  ( itime /= fit_param % nsteps ) .AND. copy ) it = it + 1

! reset logical
        copy = .FALSE.
        r ( 1, :, : ) = r ( 2, :, : )
        v ( 1, :, : ) = v ( 2, :, : )
     END DO

! Reassign fit_param % nsteps = it-1
     fit_param % nsteps = it - 1

     DEALLOCATE ( r, STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'force_ref',&
          'failed to deallocate r' )

     DEALLOCATE ( v, STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'force_ref',&
          'failed to deallocate v' )

  CASE DEFAULT
     CALL stop_program ( 'force_ref', 'no valid type' )
  END SELECT

END SUBROUTINE force_ref

!******************************************************************************

!! From the first principles trajectory file obtain the fit force
!!    using a CP2K force field.

SUBROUTINE forcefit ( fit_param, struc )

  IMPLICIT NONE

! Arguments
  TYPE ( fit_parameter_type ), INTENT ( INOUT ) :: fit_param
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc

! Locals
  INTEGER :: itime, ipart, npart

!------------------------------------------------------------------------------

  npart = SIZE ( struc % pnode )

  DO itime = 1, fit_param % nsteps
     DO ipart = 1, npart
        struc % pnode ( ipart ) % p % r ( : ) &
             = finfo ( itime, ipart ) % r ( : )
     END DO

     CALL force ( struc, thermo, simpar )

     DO ipart = 1, npart
        finfo ( itime, ipart ) % f_fit ( : ) &
             = struc % pnode ( ipart ) % p % f ( : )
     END DO

! dbg
!     v_fit ( itime ) = thermo % pot
! dbg

  END DO

END SUBROUTINE forcefit

!******************************************************************************

FUNCTION local ( x, fit_param, struc, potparm )

  IMPLICIT NONE

! Return type
  REAL ( dbl ) :: local

! Arguments
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: x
  TYPE ( fit_parameter_type ), INTENT ( INOUT ) :: fit_param
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), INTENT ( INOUT ) :: potparm

! Locals
  REAL ( dbl ), DIMENSION ( 3 ) :: vec
! dbg
!        REAL ( dbl ) :: vdiff
! dbg
  INTEGER :: itime, ipart, npart

!------------------------------------------------------------------------------

  npart = SIZE ( struc % pnode )
  CALL copy_params ( fit_param, potparm, x )
  CALL forcefit ( fit_param, struc )

  local = 0.0_dbl
! dbg
!  vdiff = 0.0_dbl
! dbg
  DO itime = 1, fit_param % nsteps
     DO ipart = 1, npart
        vec ( : ) = finfo ( itime, ipart ) % f_fit ( : ) &
             - finfo ( itime, ipart ) % f_ref ( : )
        local = local + DOT_PRODUCT ( vec, vec )
     END DO
! dbg
!     vdiff = vdiff + ( v_ref ( itime ) - v_fit ( itime )  ) ** 2
! dbg
  END DO

! The addition of the 1 is a safety mechanism so
! that one can define the FRACTIONAL tolerance for
! systems whos func value may be converged at zero.
! It does not effect other converged func values.
! dbg
!        local = local / ( 3._dbl*npart*fit_param % nsteps ) + 1._dbl
  local = local / ( 3.0_dbl * npart * fit_param % nsteps )
! dbg

! dbg
!     local = local+vdiff / fit_param % nsteps
! dbg

END FUNCTION local

!******************************************************************************

SUBROUTINE fit_print_info ( fit_param, tol, potparm )

  IMPLICIT NONE

! Arguments
  TYPE ( fit_parameter_type ), INTENT ( IN ) :: fit_param
  REAL ( dbl ), INTENT ( IN ) :: tol
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), INTENT ( INOUT ) :: potparm

! Locals
  INTEGER :: iatom, jatom

!------------------------------------------------------------------------------

  WRITE ( scr, '( )' )
  WRITE ( scr, '( A, A )' ) ' **************************************', &
       '*****************************************'
  WRITE ( scr, '( A, T71, A )' ) ' FIT TYPE = ', &
       ADJUSTR ( fit_param % type )
  WRITE ( scr, '( A, T71, I10 )' ) ' ITERATION NUMBER = ', iter
! dbg
  WRITE ( scr, '( A, T66, E15.7 )' ) ' | ( F_REF-F_FIT ) |^2 = ', &
       sum ( y ) / SIZE ( y )
!        WRITE ( scr,  '( A, T66, E15.7 )' ) ' | ( F_REF-F_FIT ) |^2 = ', &
!         sum ( y ) / SIZE ( y ) -1._dbl
! dbg
  WRITE ( scr, '( A, T66, E15.7 )' ) ' TOLERANCE = ', tol
  WRITE ( scr, '( )' )

  IF ( fit_param % type == 'LJ' ) THEN
     WRITE ( scr, '( T5, A, T15, A, T50, A, T73, A )' ) &
          'ATOMI', 'ATOMJ', 'EPSILON', 'SIGMA'
     WRITE ( scr, '( )' )
     DO iatom = 1, SIZE ( potparm, 1 )
        DO jatom = iatom, SIZE ( potparm, 2 )
           WRITE ( scr, '( T5, A, T15, A, T45, E15.8, T66, E15.8 )' ) &
                atom_names ( iatom ), atom_names ( jatom ), &
                ave_eps ( iatom, jatom ), ave_sig ( iatom, jatom )
        END DO
     END DO

  ELSE IF ( fit_param % type == 'WILL' ) THEN

     WRITE ( scr, '( T5, A, T15, A, T44, A, T60, A, T75, A )' ) &
          'ATOMI', 'ATOMJ', 'A', 'B', 'C'
     WRITE ( scr, '( )' )
     DO iatom = 1, SIZE ( potparm, 1 )
        DO jatom = iatom, SIZE ( potparm, 2 )
           WRITE ( scr, &
                '( T5, A, T15, A, T39, E12.6, T54, E12.6, T69, E12.6 )' ) &
                atom_names ( iatom ), atom_names ( jatom ), &
                ave_a ( iatom, jatom ), ave_b ( iatom, jatom ), &
                ave_c ( iatom, jatom )
        END DO
     END DO

  ELSE
     CALL stop_program ( "fit_print_info", "no valid potential type" )
  END IF

  WRITE ( scr, '( A, A )' ) ' **************************************', &
       '*****************************************'

END SUBROUTINE fit_print_info

!******************************************************************************

SUBROUTINE get_lj_params ( potparm )

  IMPLICIT NONE

! Arguments
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), INTENT ( INOUT ) :: potparm

! Locals
  INTEGER :: iatom, jatom, i
  TYPE ( linklist_real_type ), POINTER :: lleps
  TYPE ( linklist_real_type ), POINTER :: llsig

!------------------------------------------------------------------------------

  ave_sig = 0.0_dbl
  ave_eps = 0.0_dbl
  DO iatom = 1, SIZE ( potparm, 1 )
     DO jatom = iatom, SIZE ( potparm, 2 )
        lleps => param ( iatom, jatom ) % lleps
        llsig => param ( iatom, jatom ) % llsig
        DO i = 1, SIZE ( y )
           ave_sig ( iatom, jatom ) = ave_sig ( iatom, jatom ) &
                + llsig % sig ** 2
           ave_eps ( iatom, jatom ) = ave_eps ( iatom, jatom ) &
                + lleps % eps ** 2
           llsig => llsig % next
           lleps => lleps % next
        END DO
        ave_sig ( iatom, jatom ) = ave_sig ( iatom, jatom ) / SIZE ( y )
        ave_eps ( iatom, jatom ) = ave_eps ( iatom, jatom ) / SIZE ( y )
     END DO
  END DO

  ave_sig = ave_sig ** ( 1.0_dbl / 6.0_dbl )

END SUBROUTINE get_lj_params

!******************************************************************************

SUBROUTINE get_will_params ( potparm )

  IMPLICIT NONE

! Arguments
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), INTENT ( INOUT ) :: potparm

! Locals
  INTEGER :: iatom, jatom, i
  TYPE ( linklist_real_type ), POINTER :: lla
  TYPE ( linklist_real_type ), POINTER :: llb
  TYPE ( linklist_real_type ), POINTER :: llc

!------------------------------------------------------------------------------

  ave_a = 0.0_dbl
  ave_b = 0.0_dbl
  ave_c = 0.0_dbl
  DO iatom = 1, SIZE ( potparm, 1 )
     DO jatom = iatom, SIZE ( potparm, 2 )
        lla => param ( iatom, jatom ) % lla
        llb => param ( iatom, jatom ) % llb
        llc => param ( iatom, jatom ) % llc
        DO i = 1, SIZE ( y )
           ave_a ( iatom, jatom ) = ave_a ( iatom, jatom ) + lla % a * lla % a
           ave_b ( iatom, jatom ) = ave_b ( iatom, jatom ) + llb % b * llb % b
           ave_c ( iatom, jatom ) = ave_c ( iatom, jatom ) + llc % c * llc % c
           lla => lla % next
           llb => llb % next
           llc => llc % next
        END DO
        ave_a ( iatom, jatom ) = ave_a ( iatom, jatom ) / SIZE ( y )
        ave_b ( iatom, jatom ) = ave_b ( iatom, jatom ) / SIZE ( y )
        ave_c ( iatom, jatom ) = ave_c ( iatom, jatom ) / SIZE ( y )
     END DO
  END DO

END SUBROUTINE get_will_params

!******************************************************************************

SUBROUTINE fit_write_restart ( y, p )

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: y
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: p

! Locals
  INTEGER :: i, j, iochannel = 99

!------------------------------------------------------------------------------

! write to the restart file
  OPEN ( UNIT = iochannel, FILE = 'amoeba.rst', STATUS = 'UNKNOWN' )

  DO i = 1, SIZE ( y )
     WRITE ( iochannel, FMT = * ) y ( i )
  END DO

  DO i = 1, SIZE ( y )
     WRITE ( iochannel, FMT = * ) ( p ( i, j ), j = 1, SIZE ( p, 2 )  )
  END DO

  CLOSE ( UNIT = iochannel )

END SUBROUTINE fit_write_restart

!******************************************************************************

SUBROUTINE fit_read_restart ( y, p )

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), DIMENSION ( : ), INTENT ( OUT ) :: y
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ) :: p

! Locals
  INTEGER :: i, j, iochannel = 99

!------------------------------------------------------------------------------

! read from the restart file
  OPEN ( UNIT = iochannel, FILE = 'amoeba.rst', STATUS = 'UNKNOWN' )

  DO i = 1, SIZE ( y )
     READ ( iochannel, FMT = * ) y ( i )
  END DO

  DO i = 1, SIZE ( y )
     READ ( iochannel, FMT = * ) ( p ( i, j ), j = 1, SIZE ( p, 2 )  )
  END DO

  CLOSE ( UNIT = iochannel )

END SUBROUTINE fit_read_restart

!******************************************************************************

END MODULE force_fit
