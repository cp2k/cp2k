!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/initialize_extended_types [1.0] *
!!
!!   NAME
!!     initialize_extended_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     cjm, FEB 20 2001:  added subroutine initialize_extended_parameters
!!     cjm, MAY 03 2001:  reorganized and added separtate routines for
!!                        nhc_part, nhc_baro, nhc_ao, npt
!!   SOURCE
!******************************************************************************

MODULE initialize_extended_types

  USE coefficient_types, ONLY : coeff_type
  USE distribution_1d_types, ONLY : distribution_1d_type
  USE extended_system_initialize, ONLY : initialize_nhc, initialize_npt
  USE extended_system_types, ONLY :  lnhc_parameters_type, npt_info_type
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type
  USE molecule_types_new, ONLY : molecule_type
  USE molecule_kind_types, ONLY: get_molecule_kind_set,&
                                             molecule_kind_type
  USE message_passing, ONLY : mp_sum
  USE nhc, ONLY : nhc_particle_setup, nhc_barostat_setup!, nhc_coef_setup
  USE simulation_cell, ONLY : cell_type
  USE termination, ONLY : stop_program

  IMPLICIT NONE

  PUBLIC :: initialize_nhc_part, initialize_nhc_baro, initialize_npt_type, &
            initialize_nhc_forces
  PRIVATE

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** initialize_extended_types/initialize_npt_type [1.0] *
!!
!!   NAME
!!     initialize_npt_type
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!    
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE initialize_npt_type ( simpar, globenv, npt_info, box )

  IMPLICIT NONE

! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( npt_info_type ), DIMENSION ( :, : ), POINTER :: npt_info
  TYPE ( cell_type ), INTENT (inout) :: box
! locals
  INTEGER :: isos

!------------------------------------------------------------------------------

! first allocating the npt_info_type
  SELECT CASE ( simpar % ensemble )
  CASE ( 'NPT_I')
     ALLOCATE ( npt_info ( 1, 1 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'initial_npt', &
          'failed to allocate npt_info')
     npt_info  ( :, : ) % eps = LOG ( box % deth ) / 3.0_dbl

  CASE ( 'NPT_F')
     ALLOCATE ( npt_info ( 3, 3 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'initial_npt', &
          'failed to allocate npt_info')

  CASE  DEFAULT
     ALLOCATE ( npt_info ( 0, 0 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'initial_npt', &
          'failed to allocate npt_info')
  END SELECT

  IF ( simpar % read_type /= 'ALL' ) &
  CALL initialize_npt ( npt_info, simpar % tau_cell, simpar % temp_ext,  &
                       box % deth, simpar % nfree, simpar % ensemble, globenv )

END SUBROUTINE initialize_npt_type

!******************************************************************************
!!****** initialize_extended_types/initialize_nhc_baro [1.0] *
!!
!!   NAME
!!     initialize_nhc_baro
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!    
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE initialize_nhc_baro ( simpar, globenv, nhc )

  IMPLICIT NONE

! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( lnhc_parameters_type ),  INTENT ( INOUT ) :: nhc
! locals
  INTEGER :: isos

!------------------------------------------------------------------------------

! fire up the thermostats, if NPT
  CALL nhc_barostat_setup ( simpar, nhc, globenv )

  IF ( simpar % read_type /= 'ALL' ) THEN

! Initializing thermostat forces and velocities for the Nose-Hoover
! Chain variables
      IF ( nhc % nhc_len /= 0 ) &
       CALL initialize_nhc ( nhc, simpar % tau_nhc,  &
                       simpar % temp_ext, globenv )
  END IF

END SUBROUTINE initialize_nhc_baro

!******************************************************************************
!!****** initialize_extended_types/initialize_nhc_part [1.0] *
!!
!!   NAME
!!     initialize_nhc_part
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE initialize_nhc_part ( box, simpar, local_molecules, molecule, &
                                 molecule_kind_set, globenv, nhc )

  IMPLICIT NONE

! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( distribution_1d_type ), POINTER  :: local_molecules
  TYPE ( molecule_type ), POINTER :: molecule ( : )
  TYPE ( molecule_kind_type ), POINTER :: molecule_kind_set ( : )
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( lnhc_parameters_type ),  INTENT ( INOUT ) :: nhc
! locals
  INTEGER :: isos, nfree, nconstraint, iw, natom 

!------------------------------------------------------------------------------

! Initialize simpar%nfree to dim (because at this point no
! assumptions have been made about thermostatting options
! and COM is ALWAYS zeroed in this code.  If there are
! multiple thermostats etc...simpar % nfree will be set
! to zero in nhc_particle_setup.f)

  simpar % nfree = 6 - SUM ( box % perd ( 1:3 ) )

! fire up the thermostats, if not NVE
  CALL nhc_particle_setup ( simpar, local_molecules, molecule, &
                            molecule_kind_set, nhc, globenv )

  IF ( simpar % read_type /= 'ALL' ) THEN

! Initializing thermostat forces and velocities for the Nose-Hoover
! Chain variables
     IF ( nhc % nhc_len /= 0 ) &
       CALL initialize_nhc ( nhc, simpar % tau_nhc,  &
                             simpar % temp_ext, globenv )

  ENDIF

! compute the total number of degrees of freedom for temperature
  CALL get_molecule_kind_set(molecule_kind_set=molecule_kind_set,&
                             natom=natom,&
                             nconstraint=nconstraint)

  nfree = 3*natom - nconstraint - simpar%nfree

  IF (globenv%ionode.AND.globenv%print_level > 2) THEN
    iw = globenv % scr
    WRITE ( iw, '( /, A )' ) ' Calculation of degrees of freedom'
    WRITE ( iw, '( T48, A, T71, I10 )' ) '       Number of atoms:', natom
    WRITE ( iw, '( T48, A, T71, I10 )' ) ' Number of constraints:', nconstraint
    WRITE ( iw, '( T34, A, T71, I10 )' ) &
         ' Invariants(translation + rotations):', simpar % nfree
    WRITE ( iw, '( T48, A, T71, I10 )' ) '    Degrees of freedom:', nfree
    WRITE ( iw, '()' )
  END IF

  simpar%nfree = nfree


END SUBROUTINE initialize_nhc_part

!******************************************************************************
!!****** initialize_extended_types/initialize_nhc_coef [1.0] *
!!
!!   NAME
!!     initialize_nhc_coef
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************

!SUBROUTINE initialize_nhc_coef ( simpar, coeff, globenv, nhc )

!  IMPLICIT NONE

! Arguments
!  TYPE ( coeff_type ), INTENT ( IN ) :: coeff
!  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
!  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
!  TYPE ( lnhc_parameters_type ),  INTENT ( INOUT ) :: nhc
! locals
!  INTEGER :: isos

!------------------------------------------------------------------------------

! fire up the thermostats
!  CALL nhc_coef_setup ( simpar, coeff, nhc, globenv )

!  IF ( simpar % read_type_coef /= 'ALL' ) THEN

! Initializing thermostat forces and velocities for the Nose-Hoover
! Chain variables
!      IF ( nhc % nhc_len /= 0 ) CALL initialize_nhc ( nhc, simpar % tau_nhc_coef,  &
!                                      simpar % temp_coef, globenv )
!  END IF

!END SUBROUTINE initialize_nhc_coef

!!*****
!******************************************************************************
!!****** initialize_extended_types/initialize_extended_parameters [1.0] *
!!
!!   NAME
!!     initialize_extended_parameters
!!
!!   FUNCTION
!!    Assigns extended parameters from the restart file.
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE initialize_nhc_forces ( nhc_part, nhc_baro )
  IMPLICIT NONE
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_baro
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_part
! locals
  INTEGER :: i, j
  
! assign the forces
  DO i = 1, SIZE ( nhc_part % nvt, 2 )
     DO j = 2, SIZE ( nhc_part % nvt, 1 )
        nhc_part % nvt ( j, i ) % f = nhc_part % nvt ( j-1, i ) % mass * &
                                      nhc_part % nvt ( j-1, i ) % v *  &
                                      nhc_part % nvt ( j-1, i ) % v -  &
                                      nhc_part % nvt ( j, i ) % nkt
        nhc_part % nvt ( j, i ) % f = nhc_part % nvt ( j, i ) % f &
                                     /nhc_part % nvt ( j, i ) % mass
     END DO
  END DO

! assign the forces
  DO i = 1, SIZE ( nhc_baro % nvt, 2 )
     DO j = 2, SIZE ( nhc_baro % nvt, 1 )
        nhc_baro % nvt ( j, i ) % f = nhc_baro % nvt ( j-1, i ) % mass * &
                                      nhc_baro % nvt ( j-1, i ) % v *  &
                                      nhc_baro % nvt ( j-1, i ) % v -  &
                                      nhc_baro % nvt ( j, i ) % nkt
        nhc_baro % nvt ( j, i ) % f = nhc_baro % nvt ( j, i ) % f &
                                     /nhc_baro % nvt ( j, i ) % mass
     END DO
  END DO



END SUBROUTINE initialize_nhc_forces

!!*****
!******************************************************************************

END MODULE initialize_extended_types

!******************************************************************************
