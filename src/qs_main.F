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
!!     Perform a LCGTO-SCF-DF calculation using the Gaussian and plane wave
!!     (GPW) method
!!
!!   AUTHOR
!!     Matthias Krack (03.04.2001)
!!
!!   LITERATURE
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_main

! *****************************************************************************

! Index:

! SUBROUTINE qs_main

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE global_types, ONLY: global_environment_type

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: quickstep

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE quickstep(globenv)

!   Purpose: Main driver routine for a Quickstep calculation.

!   History: - Creation (04.00.2001, Matthias Krack)

!   ***************************************************************************

    USE core_hamiltonian,      ONLY: build_core_hamiltonian_matrix
    USE nl,                    ONLY: build_neighbor_list_sets
    USE qs_scf,                ONLY: scf
    USE readin,                ONLY: read_files
    USE start_program_run,     ONLY: start_qs_program_run
    USE timings,               ONLY: timeset,timestop

    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE quickstep (MODULE qs_main)"

!   *** Local variables ***

    INTEGER  :: handle

!   ---------------------------------------------------------------------------

    CALL timeset("Quickstep","I","",handle)

!   *** Initialize the Quickstep program run ***

    CALL start_qs_program_run(globenv)

!   *** Read the input and the database files ***
!   *** and perform further initializations   ***

    CALL read_files(globenv)

!   *** Build neighbor list sets ***

    CALL build_neighbor_list_sets(globenv)

!   *** Calculate the core Hamiltonian integral matrix ***

    CALL build_core_hamiltonian_matrix(globenv)

!   *** Perform SCF run ***

    CALL scf(globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE quickstep

! *****************************************************************************

END MODULE qs_main
