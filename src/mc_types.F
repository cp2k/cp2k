!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** mc_types.F *********************************************************
!!
!!   NAME
!!     mc_types
!!
!!   FUNCTION
!!     holds all the structure types needed for Monte Carlo, except
!!     the mc_environment_type
!!
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE mc_types
  USE kinds,                           ONLY: dp

  IMPLICIT NONE

  PRIVATE 
! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_types"

  PUBLIC :: mc_simulation_parameters_type,&
            mc_simulation_parameters_p_type,&
            mc_averages_type,mc_averages_p_type,&
            mc_moves_type, mc_moves_p_type,accattempt


!****************************************************************************

  TYPE mc_simulation_parameters_type
     INTEGER :: nstep
     INTEGER :: iprint
     INTEGER :: iupvolume
     INTEGER :: iuptrans
     INTEGER :: iblock
     INTEGER :: nbox
     INTEGER :: nchain
     INTEGER :: nunits
     INTEGER :: nunits_tot
     INTEGER :: nmolty
     INTEGER :: nmoves
     INTEGER :: nvolmoves
     INTEGER :: crd
     INTEGER :: ene
     INTEGER :: rm
     INTEGER :: cl
     INTEGER :: diff
     INTEGER :: data
     INTEGER :: nstart
     INTEGER :: source
     INTEGER :: group
     LOGICAL :: lclassical
     LOGICAL :: ionode
     LOGICAL :: restart
     CHARACTER ( LEN = 3 ) :: lstop
     CHARACTER ( LEN = 20 ) :: ensemble
     CHARACTER ( LEN = 40 ) :: restart_file
     CHARACTER ( LEN = 40 ) :: molecules_file
     CHARACTER ( LEN = 40 ) :: moves_file
     CHARACTER ( LEN = 40 ) :: coords_file
     CHARACTER ( LEN = 40 ) :: energy_file
     CHARACTER ( LEN = 40 ) :: displacement_file
     CHARACTER ( LEN = 40 ) :: cell_file
     CHARACTER ( LEN = 40 ) :: dat_file
     CHARACTER ( LEN = 40 ) :: data_file
     CHARACTER ( LEN = 200 ) :: box2_file
     CHARACTER ( LEN = 50 ) :: program
     REAL ( dp ) :: rmvolume
     REAL ( dp ) :: rmbond
     REAL ( dp ) :: rmangle
     REAL ( dp ) :: rmrot
     REAL ( dp ) :: rmtrans
     REAL ( dp ) :: temp
     REAL ( dp ) :: pressure
     REAL ( dp ) :: pmswap
     REAL ( dp ) :: pmvolume
     REAL ( dp ) :: pmtraion
     REAL ( dp ) :: pmtrans
     REAL ( dp ) :: BETA
     REAL ( dp ) :: epsilon
     REAL ( dp ) :: sigma
     REAL ( dp ) :: charge
  END TYPE mc_simulation_parameters_type

  TYPE mc_simulation_parameters_p_type
      TYPE (mc_simulation_parameters_type),POINTER :: mc_par
  END TYPE mc_simulation_parameters_p_type

  TYPE mc_averages_type
      REAL(KIND = dp) :: ave_energy
      REAL(KIND = dp) :: ave_energy_squared
      REAL(KIND = dp) :: ave_cell_length
      REAL(KIND = dp) :: bond1
      REAL(KIND = dp) :: bond2
      REAL(KIND = dp) :: angle
      REAL(KIND = dp) :: molecules
  END TYPE mc_averages_type

  TYPE mc_averages_p_type
      TYPE (mc_averages_type),POINTER :: averages
  END TYPE mc_averages_p_type

! *** Types ***
  TYPE mc_moves_type
      TYPE ( accattempt ), POINTER :: classical_bond
      TYPE ( accattempt ), POINTER :: classical_angle
      TYPE ( accattempt ), POINTER :: classical_trans
      TYPE ( accattempt ), POINTER :: classical_rot
      TYPE ( accattempt ), POINTER :: classical_volume
      TYPE ( accattempt ), POINTER :: bond
      TYPE ( accattempt ), POINTER :: angle
      TYPE ( accattempt ), POINTER :: trans
      TYPE ( accattempt ), POINTER :: rot
      TYPE ( accattempt ), POINTER :: swap
      TYPE ( accattempt ), POINTER :: volume
      TYPE ( accattempt ), POINTER :: Quickstep
      REAL(KIND = dp) :: classical_energy,trans_dis,qtrans_dis
      integer :: empty,grown
  END TYPE mc_moves_type

  TYPE accattempt
      INTEGER :: successes
      INTEGER :: qsuccesses
      INTEGER :: attempts
  END TYPE accattempt

  TYPE mc_moves_p_type
      TYPE(mc_moves_type), POINTER :: moves
  END TYPE mc_moves_p_type

END MODULE mc_types

!******************************************************************************
