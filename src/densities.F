!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/densities [1.0] *
!!
!!   NAME
!!     densities
!!
!!   FUNCTION
!!     definition of the electronic density type
!!
!!   AUTHOR
!!     JGH (12-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE densities
  
  USE coefficient_types, ONLY : coeff_type,  &
       coeff_zero, PLANEWAVES
  USE coefficients, ONLY : integral
  USE dft_types, ONLY : dft_control_type
  USE density_types, ONLY : density_type
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
!  USE pws, ONLY : pw_stateset_density
  USE termination, ONLY : stop_program
  USE wf_types, ONLY : wf_type, wf_spinset_type, wf_kset_type
  
  PRIVATE
  PUBLIC :: density_update, zero
  
  INTERFACE zero
     MODULE PROCEDURE density_zero_spinset, density_zero
  END INTERFACE

!!***

CONTAINS

!******************************************************************************

SUBROUTINE density_update ( psi, n, dft_control, tprint, wavepar )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( wf_type ), DIMENSION ( : ), INTENT ( IN ) :: psi
  TYPE ( dft_control_type ), INTENT ( IN ) :: dft_control
  LOGICAL, INTENT ( IN ) :: tprint
  TYPE ( density_type ), DIMENSION ( : ), INTENT ( OUT ) :: n
  TYPE ( global_environment_type ), INTENT ( IN ) :: wavepar
  
! Local
  INTEGER :: kpt, nkpts, spin, nspins, wfset, nwfsets, part
  REAL ( dbl ) :: rho_total, rho_spin
  TYPE ( coeff_type ), POINTER :: n_s
  TYPE ( wf_spinset_type ), POINTER :: psi_s
  
!------------------------------------------------------------------------------
  
  nwfsets = SIZE ( psi )
  part = 1 !TMPTMPTMP
  
  IF ( tprint .AND. wavepar % ionode &
       .AND. wavepar % print_level > 2 ) WRITE ( wavepar % scr, '()' )
  
  rho_total = 0.0_dbl
  
  DO wfset = 1, nwfsets
     nspins = SIZE ( psi ( wfset ) % psi_s )
     
     DO spin = 1, nspins
        psi_s => psi ( wfset ) % psi_s ( spin )
        n_s => n ( wfset ) % dpart ( part, spin )
        CALL zero ( n_s )
        
        nkpts = SIZE ( psi_s % psi_ks )
        
        DO kpt = 1, nkpts
           CALL stateset_density ( psi_s % psi_ks ( kpt ), n_s )
        END DO
        
        IF ( tprint .AND. wavepar % print_level > 2 ) THEN
           rho_spin = integral ( n_s )
           rho_total = rho_total + rho_spin
           WRITE ( wavepar % scr, '( "DENSITY = ", F12.6 )' ) rho_spin
        END IF
     END DO
  END DO
  
  IF ( tprint .AND. wavepar % print_level > 2 ) THEN
     WRITE ( wavepar % scr, '( "TOTAL DENSITY = ", F12.6 )' ) rho_total
     WRITE ( wavepar % scr, '()' )
  END IF
  
  ! Augment / Vb
  
  ! Symmetrise
  
END SUBROUTINE density_update

!******************************************************************************

SUBROUTINE stateset_density ( psi_ks, n_s )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( wf_kset_type ), INTENT ( IN ) :: psi_ks
  TYPE ( coeff_type ), INTENT ( OUT ) :: n_s
  
!------------------------------------------------------------------------------
  
  IF ( psi_ks % in_use == PLANEWAVES ) THEN
!     CALL pw_stateset_density ( psi_ks, n_s )
     CALL stop_program ( "stateset_density", &
          "pw_stateset_density not availabe" )
  ELSE
     CALL stop_program ( "stateset_density", &
          "no way available to obtain density" )
  END IF
  
END SUBROUTINE stateset_density

!******************************************************************************

SUBROUTINE density_zero ( n )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( density_type ), INTENT ( INOUT ) :: n
  
! Locals
  INTEGER :: spin, nspins, part
  
!------------------------------------------------------------------------------
  
  nspins = n % nspin
  part = 1 !TMPTMPTMP
  
  DO spin = 1, nspins
     CALL coeff_zero ( n % dpart ( part, spin ) )
  END DO
  
END SUBROUTINE density_zero

!******************************************************************************

SUBROUTINE density_zero_spinset ( n_s )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( INOUT ) :: n_s
  
!------------------------------------------------------------------------------
  
  CALL coeff_zero ( n_s )
  
END SUBROUTINE density_zero_spinset

!******************************************************************************

END MODULE densities
