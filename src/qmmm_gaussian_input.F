!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


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
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: rootpi
  USE cp_parser,                       ONLY: parser_release,&
                                             parser_get_object,&
                                             parser_search_string,&
                                             parser_create,&
                                             cp_parser_type
  USE memory_utilities,                ONLY: reallocate
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE global_types,                    ONLY: global_environment_type
  USE cell_types,                 ONLY: cell_type
  USE termination,                     ONLY: stop_program
  USE string_utilities,                ONLY: uppercase
  USE physcon,                         ONLY: bohr
  USE qmmm_gaussian_data
  USE input_section_types, ONLY: section_vals_val_get,section_vals_type
#include "cp_common_uses.h"
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_gaussian_input'

  PUBLIC  :: read_mm_potential, set_mm_potential
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
                               compatibility, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(qmmm_gaussian_p_type), &
         DIMENSION(:), POINTER               :: qmmm_gaussian_fns
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(cell_type), POINTER                 :: qm_cell_small
    LOGICAL, INTENT(IN)                      :: compatibility
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'read_mm_potential',&
         routineP = moduleN//':'//routineN
    INTEGER :: j, IRad, Nog, stat, Ival , Nval, istart
    LOGICAL ::  Found, Found_Radius, failure
    CHARACTER(LEN=240)                       :: line
    CHARACTER(len=10)                        :: Ftarget, Units
    REAL(KIND=dp)                           :: Radius, fconv, my_radius
    TYPE(cp_parser_type), pointer :: parser
    CHARACTER(default_string_length) :: mm_potential_file_name


    Nval = size(qmmm_gaussian_fns)
    Ival = 0
    CALL section_vals_val_get(globenv%input_file,"FORCE_EVAL%QMMM%MM_POTENTIAL_FILE_NAME",&
         c_val=mm_potential_file_name,error=error)

    nullify(parser)
    CALL parser_create(parser,mm_potential_file_name,para_env=globenv%para_env)

    search_loop: DO
       Ftarget = "&MM_FIT_POT"
       IF (Ival.EQ.Nval) EXIT search_loop
       CALL parser_search_string(parser,Ftarget,.TRUE.,found,line)
       IF (Found) THEN
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
          CALL parser_get_object  (parser,Ftarget,newline=.TRUE.)
          CPPostcondition(TRIM(Ftarget)=="RADIUS",cp_failure_level,routineP,error,Failure)
          CALL parser_get_object  (parser,radius)
          CALL parser_get_object  (parser,units)
          CALL uppercase (units)
          fconv = 1.0_dp
          IF (TRIM(units).EQ."ANGSTROM") fconv = bohr 
          Found_Radius = .FALSE.
          radius = radius * fconv
          Radius_Loop: DO J=1,SIZE(qmmm_gaussian_fns)
             IF (ABS(radius-qmmm_gaussian_fns(J)%pgf%Elp_Radius).LT.EPSILON(0.0_dp)) THEN
                Found_Radius = .TRUE.
                EXIT Radius_Loop
             END IF
          END DO Radius_Loop
          IF (.NOT.Found_Radius) THEN
             CYCLE search_loop
          END IF
          Ival = Ival + 1
          IRad = J          
          ! Read  Rmin, Rmax
          CALL parser_get_object(parser,qmmm_gaussian_fns(J)%pgf%Number_of_Gaussians,newline=.TRUE.)
          ! Allocate Vectors
          istart = 1
          IF (compatibility) THEN 
             qmmm_gaussian_fns(J)%pgf%Number_of_Gaussians = qmmm_gaussian_fns(J)%pgf%Number_of_Gaussians + 1
             istart = 2             
          END IF
          NOG = qmmm_gaussian_fns(IRad)%pgf%Number_of_Gaussians
          ALLOCATE(qmmm_gaussian_fns(IRad)%pgf%Ak(NOG),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
          ALLOCATE(qmmm_gaussian_fns(IRad)%pgf%Gk(NOG),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
          IF (compatibility) THEN
             my_radius = qmmm_gaussian_fns(J)%pgf%Elp_Radius_corr
             qmmm_gaussian_fns(IRad)%pgf%Ak(1) = 1.0_dp/radius - 2.0_dp/(rootpi*radius)
             qmmm_gaussian_fns(IRad)%pgf%Gk(1) = my_radius
          END IF
          DO J = istart, qmmm_gaussian_fns(IRad)%pgf%Number_of_Gaussians
             CALL parser_get_object(parser,qmmm_gaussian_fns(IRad)%pgf%Ak(J),newline=.TRUE.)
             CALL parser_get_object(parser,qmmm_gaussian_fns(IRad)%pgf%Gk(J))
             CALL parser_get_object(parser,units)
             CALL uppercase (units)
             fconv = 1.0_dp
             IF (TRIM(units).EQ."ANGSTROM") fconv = bohr
             qmmm_gaussian_fns(IRad)%pgf%Ak(J) = qmmm_gaussian_fns(IRad)%pgf%Ak(J) * fconv
             qmmm_gaussian_fns(IRad)%pgf%Gk(J) = qmmm_gaussian_fns(IRad)%pgf%Gk(J) * fconv
          END DO
       ELSE
!       *** Stop program, if the end of file is reached ***
          WRITE(*,'(A,F12.6,A)')"Discrepancies in Radius: ",Radius," definition." 
          CALL stop_program(routineN,moduleN,__LINE__,"Radius Value not found in MM_POTENTIAL file",globenv)
       END IF
       
    END DO search_loop

    CALL parser_release(parser,error=error)

  END SUBROUTINE read_mm_potential

  !!****f* qmmm_gaussian_input/set_mm_potential [1.0] *
  !!
  !!   NAME
  !!     set_mm_potential
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
  !!     07.2005 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE set_mm_potential(globenv, qmmm_gaussian_fns, qm_cell_small,&
                               compatibility, num_geep_gauss, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(qmmm_gaussian_p_type), &
         DIMENSION(:), POINTER               :: qmmm_gaussian_fns
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(cell_type), POINTER                 :: qm_cell_small
    LOGICAL, INTENT(IN)                      :: compatibility
    INTEGER, INTENT(IN)                      :: num_geep_gauss
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'set_mm_potential',&
         routineP = moduleN//':'//routineN
    INTEGER :: j, IRad, Nog, stat, Ival , Nval, istart
    LOGICAL ::  Found, Found_Radius, failure
    REAL(KIND=dp)                            :: rc, radius, my_radius

    failure = .FALSE.
    Nval = size(qmmm_gaussian_fns)
    DO IRad = 1, Nval
       qmmm_gaussian_fns(IRad)%pgf%Number_of_Gaussians = num_geep_gauss
       radius = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius
       istart = 0
       ! Allocate Vectors
       IF (compatibility) THEN 
          qmmm_gaussian_fns(IRad)%pgf%Number_of_Gaussians = qmmm_gaussian_fns(IRad)%pgf%Number_of_Gaussians + 1
          istart = 1
       END IF
       SELECT CASE(num_geep_gauss)
       CASE(2)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g2_rc*bohr)
       CASE(3)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g3_rc*bohr)
       CASE(4)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g4_rc*bohr)
       CASE(5)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g5_rc*bohr)
       CASE(6)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g6_rc*bohr)
       CASE(7)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g7_rc*bohr)
       CASE(8)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g8_rc*bohr)
       CASE(9)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g9_rc*bohr)
       CASE(10)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g10_rc*bohr)
       CASE(11)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g11_rc*bohr)
       CASE(12)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g12_rc*bohr)
       CASE(13)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g13_rc*bohr)
       CASE(14)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g14_rc*bohr)
       CASE(15)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g15_rc*bohr)
       CASE(16)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g16_rc*bohr)
       CASE(17)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g16_rc*bohr)
       CASE(18)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (g16_rc*bohr)
       END SELECT             
       NOG = qmmm_gaussian_fns(IRad)%pgf%Number_of_Gaussians
       ALLOCATE(qmmm_gaussian_fns(IRad)%pgf%Ak(NOG),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       ALLOCATE(qmmm_gaussian_fns(IRad)%pgf%Gk(NOG),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       IF (compatibility) THEN
          my_radius = qmmm_gaussian_fns(IRad)%pgf%Elp_Radius_corr
          qmmm_gaussian_fns(IRad)%pgf%Ak(1) = 1.0_dp/radius - 2.0_dp/(rootpi*radius)
          qmmm_gaussian_fns(IRad)%pgf%Gk(1) = my_radius
       END IF
       SELECT CASE(num_geep_gauss)
       CASE(2)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g2_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g2_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g2_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g2_b2
       CASE(3)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g3_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g3_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g3_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g3_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g3_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g3_b3
       CASE(4)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g4_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g4_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g4_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g4_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g4_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g4_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g4_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g4_b4
       CASE(5)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g5_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g5_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g5_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g5_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g5_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g5_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g5_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g5_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g5_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g5_b5
       CASE(6)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g6_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g6_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g6_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g6_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g6_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g6_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g6_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g6_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g6_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g6_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = g6_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = g6_b6
       CASE(7)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g7_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g7_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g7_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g7_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g7_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g7_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g7_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g7_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g7_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g7_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = g7_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = g7_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = g7_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = g7_b7
       CASE(8)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g8_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g8_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g8_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g8_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g8_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g8_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g8_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g8_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g8_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g8_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = g8_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = g8_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = g8_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = g8_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = g8_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = g8_b8
       CASE(9)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g9_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g9_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g9_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g9_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g9_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g9_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g9_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g9_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g9_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g9_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = g9_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = g9_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = g9_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = g9_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = g9_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = g9_b8
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = g9_a9
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = g9_b9
       CASE(10)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g10_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g10_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g10_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g10_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g10_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g10_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g10_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g10_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g10_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g10_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = g10_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = g10_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = g10_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = g10_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = g10_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = g10_b8
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = g10_a9
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = g10_b9
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= g10_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= g10_b10
       CASE(11)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g11_a1 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g11_b1 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g11_a2 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g11_b2 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g11_a3 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g11_b3 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g11_a4 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g11_b4 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g11_a5 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g11_b5 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = g11_a6 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = g11_b6 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = g11_a7 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = g11_b7 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = g11_a8 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = g11_b8 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = g11_a9 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = g11_b9 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= g11_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= g11_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= g11_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= g11_b11
       CASE(12)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g12_a1 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g12_b1 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g12_a2 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g12_b2 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g12_a3 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g12_b3 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g12_a4 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g12_b4 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g12_a5 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g12_b5 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = g12_a6 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = g12_b6 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = g12_a7 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = g12_b7 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = g12_a8 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = g12_b8 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = g12_a9 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = g12_b9 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= g12_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= g12_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= g12_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= g12_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= g12_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= g12_b12
       CASE(13)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g13_a1 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g13_b1 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g13_a2 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g13_b2 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g13_a3 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g13_b3 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g13_a4 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g13_b4 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g13_a5 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g13_b5 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = g13_a6 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = g13_b6 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = g13_a7 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = g13_b7 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = g13_a8 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = g13_b8 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = g13_a9 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = g13_b9 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= g13_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= g13_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= g13_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= g13_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= g13_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= g13_b12
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+13)= g13_a13
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+13)= g13_b13
       CASE(14)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g14_a1 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g14_b1 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g14_a2 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g14_b2 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g14_a3 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g14_b3 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g14_a4 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g14_b4 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g14_a5 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g14_b5 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = g14_a6 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = g14_b6 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = g14_a7 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = g14_b7 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = g14_a8 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = g14_b8 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = g14_a9 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = g14_b9 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= g14_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= g14_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= g14_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= g14_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= g14_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= g14_b12
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+13)= g14_a13
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+13)= g14_b13
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+14)= g14_a14
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+14)= g14_b14
       CASE(15)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g15_a1 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g15_b1 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g15_a2 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g15_b2 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g15_a3 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g15_b3 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g15_a4 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g15_b4 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g15_a5 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g15_b5 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = g15_a6 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = g15_b6 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = g15_a7 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = g15_b7 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = g15_a8 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = g15_b8 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = g15_a9 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = g15_b9 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= g15_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= g15_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= g15_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= g15_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= g15_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= g15_b12
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+13)= g15_a13
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+13)= g15_b13
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+14)= g15_a14
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+14)= g15_b14
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+15)= g15_a15
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+15)= g15_b15
       CASE(16)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g16_a1 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g16_b1 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g16_a2 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g16_b2 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g16_a3 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g16_b3 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g16_a4 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g16_b4 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g16_a5 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g16_b5 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = g16_a6 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = g16_b6 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = g16_a7 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = g16_b7 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = g16_a8 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = g16_b8 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = g16_a9 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = g16_b9 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= g16_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= g16_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= g16_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= g16_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= g16_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= g16_b12
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+13)= g16_a13
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+13)= g16_b13
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+14)= g16_a14
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+14)= g16_b14
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+15)= g16_a15
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+15)= g16_b15
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+16)= g16_a16
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+16)= g16_b16
       CASE(17)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g17_a1 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g17_b1 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g17_a2 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g17_b2 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g17_a3 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g17_b3 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g17_a4 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g17_b4 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g17_a5 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g17_b5 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = g17_a6 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = g17_b6 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = g17_a7 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = g17_b7 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = g17_a8 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = g17_b8 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = g17_a9 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = g17_b9 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= g17_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= g17_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= g17_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= g17_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= g17_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= g17_b12
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+13)= g17_a13
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+13)= g17_b13
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+14)= g17_a14
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+14)= g17_b14
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+15)= g17_a15
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+15)= g17_b15
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+16)= g17_a16
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+16)= g17_b16
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+17)= g17_a17
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+17)= g17_b17
       CASE(18)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = g18_a1 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = g18_b1 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = g18_a2 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = g18_b2 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = g18_a3 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = g18_b3 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = g18_a4 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = g18_b4 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = g18_a5 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = g18_b5 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = g18_a6 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = g18_b6 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = g18_a7 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = g18_b7 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = g18_a8 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = g18_b8 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = g18_a9 
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = g18_b9 
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= g18_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= g18_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= g18_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= g18_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= g18_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= g18_b12
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+13)= g18_a13
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+13)= g18_b13
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+14)= g18_a14
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+14)= g18_b14
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+15)= g18_a15
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+15)= g18_b15
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+16)= g18_a16
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+16)= g18_b16
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+17)= g18_a17
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+17)= g18_b17
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+18)= g18_a18
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+18)= g18_b18
       END SELECT
       qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1:) = qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1:) / rc
       qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1:) = qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1:) * rc
    END DO
  END SUBROUTINE set_mm_potential


END MODULE qmmm_gaussian_input

