!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qmmm_gaussian_input [1.0] *
!!
!!   NAME
!!     qmmm_gaussian_input
!!
!!   FUNCTION
!!     Set the QMMM Gaussian Input Environment
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     6.2004 created [tlaino]
!!
!!   SOURCE
!****************************************************************************
MODULE qmmm_gaussian_input
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dbl
  USE qs_parser,                       ONLY: finish_parser,&
                                             open_file,&
                                             read_object,&
                                             search_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE memory_utilities,                ONLY: reallocate
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE global_types,                    ONLY: global_environment_type
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_program
  USE string_utilities,                ONLY: uppercase
  USE physcon,                         ONLY: bohr
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_gaussian_input'

  PUBLIC  :: read_mm_potential
!***
!****************************************************************************
CONTAINS

  !!****f* qmmm_gaussian_input/read_mm_potential [1.0] *
  !!
  !!   NAME
  !!     read_mm_potential
  !!
  !!   FUNCTION
  !!     read MM_POTENTIAL file
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     06.2004 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE read_mm_potential(globenv, qmmm_gaussian_fns, qm_cell_small,&
                               error)
    IMPLICIT NONE
    ! Arguments
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(qmmm_gaussian_p_type), &
         DIMENSION(:), POINTER               :: qmmm_gaussian_fns
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(cell_type), POINTER                 :: qm_cell_small
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'read_mm_potential',&
         routineP = moduleN//':'//routineN
    INTEGER :: j, IRad, Nog, stat, Ival , Nval
    LOGICAL ::  Found, Found_Radius, failure
    CHARACTER(LEN=240)                       :: line
    CHARACTER(len=10)                        :: Ftarget, Units
    REAL(KIND=dbl)                           :: Radius, fconv


    Nval = size(qmmm_gaussian_fns)
    Ival = 0
    CALL start_parser(globenv%mm_potential_file_name,globenv)

    search_loop: DO
       Ftarget = "&MM_FIT_POT"
       IF (Ival.EQ.Nval) EXIT search_loop
       CALL search_object(Ftarget,.TRUE.,found,line)
       IF (Found) THEN
          Ival = Ival + 1
!
! Structure example of the MM fit potential file:
!
!           &MM_FIT_POT 
!           RADIUS  0.4 Angstrom
!           7   
!           0.223396   0.811453  Bohr
!           0.306814   1.01988   Bohr
!           0.254879   1.37404   Bohr
!           0.188293   1.87929   Bohr
!           0.136391   2.56745   Bohr
!           0.100305   3.50033   Bohr
!           0.0790169  4.82046   Bohr
!           &END
!

          CALL read_object  (Ftarget,newline=.TRUE.)
          CPPostcondition(TRIM(Ftarget)=="RADIUS",cp_failure_level,routineP,error,Failure)
          CALL read_object  (radius)
          CALL read_object  (units)
          CALL uppercase (units)
          fconv = 1.0_dbl
          IF (TRIM(units).EQ."ANGSTROM") fconv = bohr 
          Found_Radius = .FALSE.
          radius = radius * fconv
          Radius_Loop: DO J=1,SIZE(qmmm_gaussian_fns)
             IF (ABS(radius-qmmm_gaussian_fns(J)%pgf%Elp_Radius).LT.EPSILON(0.0_dbl)) THEN
                Found_Radius = .TRUE.
                EXIT Radius_Loop
             END IF
          END DO Radius_Loop
          IF (.NOT.Found_Radius) THEN
             WRITE(*,'(A,F12.6,A)')"Discrepancies in Radius: ",Radius," definition." 
             CALL stop_parser(routine=routineN,message="ERROR IN READ_MM_POTENTIAL: "&
                  //"RADIUS FIELD IN QMMM_GAUSSIAN_TYPE NOT FOUND!")   
          END IF
          IRad = J          
          ! Read  Rmin, Rmax
          CALL read_object(qmmm_gaussian_fns(J)%pgf%Number_of_Gaussians,newline=.TRUE.)
          ! Allocate Vectors
          NOG = qmmm_gaussian_fns(IRad)%pgf%Number_of_Gaussians
          ALLOCATE(qmmm_gaussian_fns(IRad)%pgf%Ak(NOG),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
          ALLOCATE(qmmm_gaussian_fns(IRad)%pgf%Gk(NOG),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
          DO J = 1, qmmm_gaussian_fns(IRad)%pgf%Number_of_Gaussians
             CALL read_object(qmmm_gaussian_fns(IRad)%pgf%Ak(J),newline=.TRUE.)
             CALL read_object(qmmm_gaussian_fns(IRad)%pgf%Gk(J))
             CALL read_object(units)
             CALL uppercase (units)
             fconv = 1.0_dbl
             IF (TRIM(units).EQ."ANGSTROM") fconv = bohr
             qmmm_gaussian_fns(IRad)%pgf%Ak(J) = qmmm_gaussian_fns(IRad)%pgf%Ak(J) * fconv
             qmmm_gaussian_fns(IRad)%pgf%Gk(J) = qmmm_gaussian_fns(IRad)%pgf%Gk(J) * fconv
          END DO
       ELSE
!       *** Stop program, if the end of file is reached ***
          CALL stop_program(routineN,moduleN,__LINE__,"Error in MM_POTENTIAL file",globenv)
       END IF
       
    END DO search_loop

    CALL finish_parser()

  END SUBROUTINE read_mm_potential

END MODULE qmmm_gaussian_input

