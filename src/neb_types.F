!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/neb_types [1.0] *
!!
!!   NAME
!!     neb_types
!!
!!   FUNCTION
!!     Typo for Nudged Elastic Band Calculation
!!
!!   NOTES
!!     Numerical accuracy for parallel runs:
!!      Each replica starts the SCF run from the one optimized
!!      in a previous run. It may happen then energies and derivatives
!!      of a serial run and a parallel run could be slightly different
!!      'cause of a different starting density matrix.
!!      Exact results are obtained using:
!!         EXTRAPOLATION USE_GUESS in QS section (Teo 09.2006)
!!
!!   AUTHOR
!!     Teodoro Laino 10.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
MODULE neb_types
  USE force_env_types,                 ONLY: force_env_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'neb_types'

  PUBLIC :: neb_type

  TYPE neb_type
     ! NEB parameters
     LOGICAL       :: rotate_frames, optimize_end_points
     INTEGER       :: id_type
     INTEGER       :: opt_type
     INTEGER       :: number_of_replica
     INTEGER       :: nsteps_it, istep
     INTEGER       :: nr_HE_image
     REAL(KIND=dp) :: K, spring_energy, avg_distance
     CHARACTER(LEN=default_string_length) :: opt_type_label
     ! Section used for restart
     TYPE(force_env_type), POINTER        :: force_env
     TYPE(section_vals_type), POINTER     :: root_section
     TYPE(section_vals_type), POINTER     :: motion_print_section
     TYPE(section_vals_type), POINTER     :: force_env_section
     TYPE(section_vals_type), POINTER     :: neb_section
  END TYPE neb_type

END MODULE neb_types
