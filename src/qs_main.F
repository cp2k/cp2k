!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_main [1.0] *
!!
!!   NAME
!!     qs_main
!!
!!   FUNCTION
!!     Perform a LCGTO-SCF-DF calculation using the Gaussian plane wave (GPW)
!!     method
!!
!!   AUTHOR
!!     Matthias Krack (03.04.2001)
!!
!!   LITERATURE
!!     G. Lippert, J. Hutter, M. Parrinello, Mol. Phys. 92, 477 (1997)
!!     G. Lippert, J. Hutter, M. Parrinello, Theor. Chem. Acc. 103, 124 (1999)
!!     M. Krack, M. Parrinello, Phys. Chem. Chem. Phys. 2, 2105 (2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_main

  USE kinds, ONLY: wp => dp

  USE global_types,         ONLY: global_environment_type
  USE qs_core_hamiltonian,  ONLY: build_core_hamiltonian_matrix
  USE qs_environment_types, ONLY: qs_environment_type
  USE qs_init,              ONLY: read_qs_env
  USE qs_neighbor_lists,    ONLY: build_qs_neighbor_lists
  USE qs_overlap,           ONLY: distribute_overlap_matrix
  USE qs_scf,               ONLY: scf
  USE timings,              ONLY: timeset,&
                                  timestop

  IMPLICIT NONE

  PRIVATE

  TYPE(qs_environment_type) :: qs_env

  PUBLIC :: quickstep

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE quickstep(globenv)

!   Purpose: Main driver routine for a Quickstep calculation.

!   History: - Creation (23.01.2002, MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: handle

!   ---------------------------------------------------------------------------

    CALL timeset("Quickstep","I","",handle)

!   *** Read the input and the database files and perform further  ***
!   *** initializations for the setup of the QUICKSTEP environment ***

    CALL read_qs_env(qs_env,globenv)

!   *** Build neighbor list sets ***

    CALL build_qs_neighbor_lists(qs_env,globenv)

!   *** Allocate the overlap matrix and distribute the atomic blocks ***

    CALL distribute_overlap_matrix(qs_env,globenv)

!   *** Calculate the overlap and the core Hamiltonian integral matrix ***

    CALL build_core_hamiltonian_matrix(qs_env,globenv,.FALSE.)

!   *** Perform a SCF run ***

    CALL scf(qs_env,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE quickstep

! *****************************************************************************

END MODULE qs_main
