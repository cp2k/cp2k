!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!


MODULE energy_output

  USE energy_input,                    ONLY: min_info_type
  USE global_types,                    ONLY: global_environment_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: print_energy_info, print_convergence                 

CONTAINS

!******************************************************************************
!
  SUBROUTINE print_energy_info ( iter, min_info, min_type, globenv )

    INTEGER, INTENT(IN)                      :: iter
    TYPE(min_info_type), INTENT(IN)          :: min_info
    CHARACTER(len=*), INTENT(IN)             :: min_type
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    integer                                  :: iw

    iw = globenv % scr
    
    WRITE (iw,*)
    WRITE (iw,'(A,A)') ' **************************************', &
      '*****************************************'
    WRITE (iw,'(A,T71,I10)') ' ITERATION NUMBER        = ', iter
    WRITE (iw,'(A,T66,E15.7)') ' MAX FORCE ON COEFFICENTS = ', min_info%tol
    IF ( min_type == 'DIIS' ) &
    WRITE (iw,'(A,T66,E15.7)') ' ODIIS ERROR              = ', min_info%err
    WRITE (iw,'(A,T66,E15.7)') ' CLASSIC  ENERGY    = ', min_info%e_class
    WRITE (iw,'(A,T66,E15.7)') ' TOTAL OVERLAP ENERGY   = ', min_info%e_over
    WRITE (iw,'(A,T66,E15.7)') ' <i|i> OVERLAP ENERGY   = ', min_info%e_over_ii
    WRITE (iw,'(A,T66,E15.7)') ' <i|j> OVERLAP ENERGY   = ', min_info%e_over_ij
    WRITE (iw,'(A,T66,E15.7)') ' TOTAL ELECTROSTATIC ENERGY  = ', min_info%e_electro
    WRITE (iw,'(A,T66,E15.7)') ' RHO0-RHO0 ELECTR ENERGY  = ', min_info%e_electro_0
    WRITE (iw,'(A,T66,E15.7)') ' RHO0-DRHO ELECTR ENERGY  = ', min_info%e_electro_1
    WRITE (iw,'(A,T66,E15.7)') ' DRHO-DRHO ELECTR ENERGY  = ', min_info%e_electro_2
    WRITE (iw,'(A,T66,E15.7)') ' HK CONTRIBUTION     = ', min_info%e_thk
    WRITE (iw,'(A,T66,E15.7)') ' TOTAL POLARIZATION ENERGY = ', min_info%e_pol
    WRITE (iw,'(A,T66,E15.7)') ' TOTAL ENERGY       = ', min_info%e
    WRITE (iw,'(A,T66,E15.7)') ' DELTA E           = ', min_info%de
    WRITE (iw,'(A,A)') ' **************************************', &
      '*****************************************'
!
  END SUBROUTINE print_energy_info
!
!******************************************************************************

  SUBROUTINE print_convergence ( conv_flag, globenv )

    LOGICAL, INTENT(IN)                      :: conv_flag
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    integer                                  :: iw

    iw = globenv % scr

     IF (conv_flag) THEN
       WRITE (iw,*)
       WRITE (iw,'(T23,A)') '******************************'
       WRITE (iw,'(T23,A)') '*    CONVERGENCE ACHIEVED    *'
       WRITE (iw,'(T23,A)') '*    IN COEF MINIMIZATION    *'
       WRITE (iw,'(T23,A)') '******************************'
       WRITE (iw,*)
     ELSE
       WRITE (iw,*)
       WRITE (iw,'(T23,A)') '******************************'
       WRITE (iw,'(T23,A)') '*  WARNING: NO CONVERGENCE!  *'
       WRITE (iw,'(T23,A)') '*  RESTART WITH LAST COEF    *'
       WRITE (iw,'(T23,A)') '******************************'
       WRITE (iw,*)
     END IF

  END SUBROUTINE print_convergence 

!******************************************************************************

END MODULE energy_output

!******************************************************************************
