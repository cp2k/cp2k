!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/potentials [1.0] *
!!
!!   NAME
!!     potentials
!!
!!   FUNCTION
!!     operations on the local potential type
!!
!!   AUTHOR
!!     JGH (13-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE potentials
  
  USE coefficient_types,               ONLY: coeff_add,&
                                             coeff_allocate,&
                                             coeff_copy,&
                                             coeff_deallocate,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_control_types,                ONLY: dft_control_type
  USE density_types,                   ONLY: density_type
  USE hartree,                         ONLY: calculate_hartree
  USE kinds,                           ONLY: dbl
  USE potential_types,                 ONLY: potential_type
!JH  USE xc, ONLY : calculate_xc


  PRIVATE
  PUBLIC :: potential_update, zero
  
  INTERFACE zero
     MODULE PROCEDURE potential_zero_spinset, potential_zero, &
          potential_set_zero
  END INTERFACE
  
!!***
CONTAINS

!******************************************************************************

SUBROUTINE potential_update ( n, vloc, dft_control, tprint )
  
  
    TYPE(density_type), DIMENSION(:), &
      INTENT(IN)                             :: n
    TYPE(potential_type), DIMENSION(:), &
      INTENT(INOUT), TARGET                  :: vloc
    TYPE(dft_control_type), INTENT(IN)       :: dft_control
    LOGICAL, INTENT(IN)                      :: tprint

    INTEGER                                  :: densityset, ndensitysets, &
                                                nspins, spin
    REAL(dbl)                                :: e_hartree
    TYPE(coeff_type)                         :: n_total
    TYPE(potential_type), POINTER            :: v_ext, v_hartree, v_xc

! TMPTMPTMP
!------------------------------------------------------------------------------

  ndensitysets = SIZE ( n )
  
  DO densityset = 1, ndensitysets
     nspins = n ( densityset ) % nspin
     
     DO spin = 1, nspins
        CALL zero ( vloc )
     END DO
     
     v_xc => vloc ( densityset )
     v_hartree => vloc ( densityset )
     v_ext => vloc ( densityset )
     
! Sum up the spin densities into the total electron density
     CALL coeff_allocate ( n_total, n ( densityset ) % dpart ( 1, 1 ) )
     IF ( nspins == 2 ) THEN
        CALL coeff_add ( n ( densityset ) % dpart ( 1, 1 ), &
             n ( densityset ) % dpart ( nspins, 1 ), n_total )
     ELSE
        CALL coeff_copy ( n ( densityset ) % dpart ( 1, 1 ), n_total )
     END IF
     
! Calculate the Hartree energy and potential
!   (depends only on the total density) into the spin up part
     CALL calculate_hartree ( n_total, e_hartree, &
          v_hartree % pot_part ( 1, 1 ) )
     
     CALL coeff_deallocate ( n_total )
     
! Because the potential is the same for spin up and dn components
     IF ( nspins == 2 ) THEN
        CALL coeff_copy ( v_hartree % pot_part ( 1, 1 ), &
             v_hartree % pot_part ( 1, nspins ) )
     END IF
     
! Calculate the exchange-correlation energy and potential
!JH     CALL calculate_xc ( n ( densityset ), e_xc, xc_pot = v_xc )
     
! Calculate the contribution due to the external potential
!*apsi*TMPTMPTMP     CALL potential_external_evaluate ( n ( densityset ), e_ext, v_ext )
     
  END DO
  
END SUBROUTINE potential_update

!******************************************************************************

SUBROUTINE potential_set_zero ( vloc )
  
  
    TYPE(potential_type), DIMENSION(:), &
      INTENT(INOUT)                          :: vloc

    INTEGER                                  :: npotsets, nspins, potset, spin

!------------------------------------------------------------------------------

  npotsets = SIZE ( vloc )
  
  DO potset = 1, npotsets
     nspins = vloc ( potset ) % nspin
     
     DO spin = 1, nspins
        CALL coeff_zero ( vloc ( potset ) % pot_part ( spin, 1 ) )
     END DO
  END DO
  
END SUBROUTINE potential_set_zero

!******************************************************************************

SUBROUTINE potential_zero ( vloc )
  
  
    TYPE(potential_type), INTENT(INOUT)      :: vloc

    INTEGER                                  :: nspins, spin

!------------------------------------------------------------------------------

  nspins = vloc % nspin
  
  DO spin = 1, nspins
     CALL coeff_zero ( vloc % pot_part ( spin, 1 ) )
  END DO
  
END SUBROUTINE potential_zero

!******************************************************************************

SUBROUTINE potential_zero_spinset ( vloc_s )
  
  
    TYPE(coeff_type), INTENT(INOUT)          :: vloc_s

!------------------------------------------------------------------------------

  CALL coeff_zero ( vloc_s )
  
END SUBROUTINE potential_zero_spinset

!******************************************************************************

END MODULE potentials
