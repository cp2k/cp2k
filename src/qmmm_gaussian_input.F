!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Set the QMMM Gaussian Input Environment
!> \par History
!>      6.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE qmmm_gaussian_input
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_object,&
                                             parser_search_string
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: rootpi
  USE physcon,                         ONLY: bohr
  USE qmmm_gaussian_data,              ONLY: &
       g10_a1, g10_a10, g10_a2, g10_a3, g10_a4, g10_a5, g10_a6, g10_a7, &
       g10_a8, g10_a9, g10_b1, g10_b10, g10_b2, g10_b3, g10_b4, g10_b5, &
       g10_b6, g10_b7, g10_b8, g10_b9, g10_rc, g11_a1, g11_a10, g11_a11, &
       g11_a2, g11_a3, g11_a4, g11_a5, g11_a6, g11_a7, g11_a8, g11_a9, &
       g11_b1, g11_b10, g11_b11, g11_b2, g11_b3, g11_b4, g11_b5, g11_b6, &
       g11_b7, g11_b8, g11_b9, g11_rc, g12_a1, g12_a10, g12_a11, g12_a12, &
       g12_a2, g12_a3, g12_a4, g12_a5, g12_a6, g12_a7, g12_a8, g12_a9, &
       g12_b1, g12_b10, g12_b11, g12_b12, g12_b2, g12_b3, g12_b4, g12_b5, &
       g12_b6, g12_b7, g12_b8, g12_b9, g12_rc, g13_a1, g13_a10, g13_a11, &
       g13_a12, g13_a13, g13_a2, g13_a3, g13_a4, g13_a5, g13_a6, g13_a7, &
       g13_a8, g13_a9, g13_b1, g13_b10, g13_b11, g13_b12, g13_b13, g13_b2, &
       g13_b3, g13_b4, g13_b5, g13_b6, g13_b7, g13_b8, g13_b9, g13_rc, &
       g14_a1, g14_a10, g14_a11, g14_a12, g14_a13, g14_a14, g14_a2, g14_a3, &
       g14_a4, g14_a5, g14_a6, g14_a7, g14_a8, g14_a9, g14_b1, g14_b10, &
       g14_b11, g14_b12, g14_b13, g14_b14, g14_b2, g14_b3, g14_b4, g14_b5, &
       g14_b6, g14_b7, g14_b8, g14_b9, g14_rc, g15_a1, g15_a10, g15_a11, &
       g15_a12, g15_a13, g15_a14, g15_a15, g15_a2, g15_a3, g15_a4, g15_a5, &
       g15_a6, g15_a7, g15_a8, g15_a9, g15_b1, g15_b10, g15_b11, g15_b12, &
       g15_b13, g15_b14, g15_b15, g15_b2, g15_b3, g15_b4, g15_b5, g15_b6, &
       g15_b7, g15_b8, g15_b9, g15_rc, g16_a1, g16_a10, g16_a11, g16_a12, &
       g16_a13, g16_a14, g16_a15, g16_a16, g16_a2, g16_a3, g16_a4, g16_a5, &
       g16_a6, g16_a7, g16_a8, g16_a9, g16_b1, g16_b10, g16_b11, g16_b12, &
       g16_b13, g16_b14, g16_b15, g16_b16, g16_b2, g16_b3, g16_b4, g16_b5, &
       g16_b6, g16_b7, g16_b8, g16_b9, g16_rc, g17_a1, g17_a10, g17_a11, &
       g17_a12, g17_a13, g17_a14, g17_a15, g17_a16, g17_a17, g17_a2, g17_a3, &
       g17_a4, g17_a5, g17_a6, g17_a7, g17_a8, g17_a9, g17_b1, g17_b10, &
       g17_b11, g17_b12, g17_b13, g17_b14, g17_b15, g17_b16, g17_b17, g17_b2, &
       g17_b3, g17_b4, g17_b5, g17_b6, g17_b7, g17_b8, g17_b9, g18_a1, &
       g18_a10, g18_a11, g18_a12, g18_a13, g18_a14, g18_a15, g18_a16, &
       g18_a17, g18_a18, g18_a2, g18_a3, g18_a4, g18_a5, g18_a6, g18_a7, &
       g18_a8, g18_a9, g18_b1, g18_b10, g18_b11, g18_b12, g18_b13, g18_b14, &
       g18_b15, g18_b16, g18_b17, g18_b18, g18_b2, g18_b3, g18_b4, g18_b5, &
       g18_b6, g18_b7, g18_b8, g18_b9, g2_a1, g2_a2, g2_b1, g2_b2, g2_rc, &
       g3_a1, g3_a2, g3_a3, g3_b1, g3_b2, g3_b3, g3_rc, g4_a1, g4_a2, g4_a3, &
       g4_a4, g4_b1, g4_b2, g4_b3, g4_b4, g4_rc, g5_a1, g5_a2, g5_a3, g5_a4, &
       g5_a5, g5_b1, g5_b2, g5_b3, g5_b4, g5_b5, g5_rc, g6_a1, g6_a2, g6_a3, &
       g6_a4, g6_a5, g6_a6, g6_b1, g6_b2, g6_b3, g6_b4, g6_b5, g6_b6, g6_rc, &
       g7_a1, g7_a2, g7_a3, g7_a4, g7_a5, g7_a6, g7_a7, g7_b1, g7_b2, g7_b3, &
       g7_b4, g7_b5, g7_b6, g7_b7, g7_rc, g8_a1, g8_a2, g8_a3, g8_a4, g8_a5, &
       g8_a6, g8_a7, g8_a8, g8_b1, g8_b2, g8_b3, g8_b4, g8_b5, g8_b6, g8_b7, &
       g8_b8, g8_rc, g9_a1, g9_a2, g9_a3, g9_a4, g9_a5, g9_a6, g9_a7, g9_a8, &
       g9_a9, g9_b1, g9_b2, g9_b3, g9_b4, g9_b5, g9_b6, g9_b7, g9_b8, g9_b9, &
       g9_rc, s10_a1, s10_a10, s10_a2, s10_a3, s10_a4, s10_a5, s10_a6, &
       s10_a7, s10_a8, s10_a9, s10_b1, s10_b10, s10_b2, s10_b3, s10_b4, &
       s10_b5, s10_b6, s10_b7, s10_b8, s10_b9, s10_rc, s11_a1, s11_a10, &
       s11_a11, s11_a2, s11_a3, s11_a4, s11_a5, s11_a6, s11_a7, s11_a8, &
       s11_a9, s11_b1, s11_b10, s11_b11, s11_b2, s11_b3, s11_b4, s11_b5, &
       s11_b6, s11_b7, s11_b8, s11_b9, s11_rc, s12_a1, s12_a10, s12_a11, &
       s12_a12, s12_a2, s12_a3, s12_a4, s12_a5, s12_a6, s12_a7, s12_a8, &
       s12_a9, s12_b1, s12_b10, s12_b11, s12_b12, s12_b2, s12_b3, s12_b4, &
       s12_b5, s12_b6, s12_b7, s12_b8, s12_b9, s12_rc, s13_a1, s13_a10, &
       s13_a11, s13_a12, s13_a13, s13_a2, s13_a3, s13_a4, s13_a5, s13_a6, &
       s13_a7, s13_a8, s13_a9, s13_b1, s13_b10, s13_b11, s13_b12, s13_b13, &
       s13_b2, s13_b3, s13_b4, s13_b5, s13_b6, s13_b7, s13_b8, s13_b9, &
       s13_rc, s14_a1, s14_a10, s14_a11, s14_a12, s14_a13, s14_a14, s14_a2, &
       s14_a3, s14_a4, s14_a5, s14_a6, s14_a7, s14_a8, s14_a9, s14_b1, &
       s14_b10, s14_b11, s14_b12, s14_b13, s14_b14, s14_b2, s14_b3, s14_b4, &
       s14_b5, s14_b6, s14_b7, s14_b8, s14_b9, s14_rc, s15_a1, s15_a10, &
       s15_a11, s15_a12, s15_a13, s15_a14, s15_a15, s15_a2, s15_a3, s15_a4, &
       s15_a5, s15_a6, s15_a7, s15_a8, s15_a9, s15_b1, s15_b10, s15_b11, &
       s15_b12, s15_b13, s15_b14, s15_b15, s15_b2, s15_b3, s15_b4, s15_b5, &
       s15_b6, s15_b7, s15_b8, s15_b9, s15_rc, s16_a1, s16_a10, s16_a11, &
       s16_a12, s16_a13, s16_a14, s16_a15, s16_a16, s16_a2, s16_a3, s16_a4, &
       s16_a5, s16_a6, s16_a7, s16_a8, s16_a9, s16_b1, s16_b10, s16_b11, &
       s16_b12, s16_b13, s16_b14, s16_b15, s16_b16, s16_b2, s16_b3, s16_b4, &
       s16_b5, s16_b6, s16_b7, s16_b8, s16_b9, s16_rc, s17_a1, s17_a10, &
       s17_a11, s17_a12, s17_a13, s17_a14, s17_a15, s17_a16, s17_a17, s17_a2, &
       s17_a3, s17_a4, s17_a5, s17_a6, s17_a7, s17_a8, s17_a9, s17_b1, &
       s17_b10, s17_b11, s17_b12, s17_b13, s17_b14, s17_b15, s17_b16, &
       s17_b17, s17_b2, s17_b3, s17_b4, s17_b5, s17_b6, s17_b7, s17_b8, &
       s17_b9, s17_rc, s18_a1, s18_a10, s18_a11, s18_a12, s18_a13, s18_a14, &
       s18_a15, s18_a16, s18_a17, s18_a18, s18_a2, s18_a3, s18_a4, s18_a5, &
       s18_a6, s18_a7, s18_a8, s18_a9, s18_b1, s18_b10, s18_b11, s18_b12, &
       s18_b13, s18_b14, s18_b15, s18_b16, s18_b17, s18_b18, s18_b2, s18_b3, &
       s18_b4, s18_b5, s18_b6, s18_b7, s18_b8, s18_b9, s18_rc, s2_a1, s2_a2, &
       s2_b1, s2_b2, s2_rc, s3_a1, s3_a2, s3_a3, s3_b1, s3_b2, s3_b3, s3_rc, &
       s4_a1, s4_a2, s4_a3, s4_a4, s4_b1, s4_b2, s4_b3, s4_b4, s4_rc, s5_a1, &
       s5_a2, s5_a3, s5_a4, s5_a5, s5_b1, s5_b2, s5_b3, s5_b4, s5_b5, s5_rc, &
       s6_a1, s6_a2, s6_a3, s6_a4, s6_a5, s6_a6, s6_b1, s6_b2, s6_b3, s6_b4, &
       s6_b5, s6_b6, s6_rc, s7_a1, s7_a2, s7_a3, s7_a4, s7_a5, s7_a6, s7_a7, &
       s7_b1, s7_b2, s7_b3, s7_b4, s7_b5, s7_b6, s7_b7, s7_rc, s8_a1, s8_a2, &
       s8_a3, s8_a4, s8_a5, s8_a6, s8_a7, s8_a8, s8_b1, s8_b2, s8_b3, s8_b4, &
       s8_b5, s8_b6, s8_b7, s8_b8, s8_rc, s9_a1, s9_a2, s9_a3, s9_a4, s9_a5, &
       s9_a6, s9_a7, s9_a8, s9_a9, s9_b1, s9_b2, s9_b3, s9_b4, s9_b5, s9_b6, &
       s9_b7, s9_b8, s9_b9, s9_rc
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_program
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_gaussian_input'

  PUBLIC  :: read_mm_potential,&
             set_mm_potential_swave,&
             set_mm_potential_erf
!***
CONTAINS

! *****************************************************************************
!> \brief read MM_POTENTIAL file
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE read_mm_potential(para_env, qmmm_gaussian_fns, qm_cell_small,&
                               compatibility, qmmm_section, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: qmmm_gaussian_fns
    TYPE(cell_type), POINTER                 :: qm_cell_small
    LOGICAL, INTENT(IN)                      :: compatibility
    TYPE(section_vals_type), POINTER         :: qmmm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_mm_potential', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=240)                       :: line
    CHARACTER(LEN=default_string_length)     :: Ftarget, &
                                                mm_potential_file_name, Units
    INTEGER                                  :: IRad, istart, Ival, j, Nog, &
                                                Nval, output_unit, stat
    LOGICAL                                  :: failure, Found, Found_Radius
    REAL(KIND=dp)                            :: fconv, my_radius, Radius
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    Nval = SIZE(qmmm_gaussian_fns)
    Ival = 0
    CALL section_vals_val_get(qmmm_section,"MM_POTENTIAL_FILE_NAME",&
         c_val=mm_potential_file_name,error=error)

    NULLIFY(parser)
    CALL parser_create(parser,mm_potential_file_name,para_env=para_env,error=error)

    search_loop: DO
       Ftarget = "&MM_FIT_POT"
       IF (Ival.EQ.Nval) EXIT search_loop
       CALL parser_search_string(parser,Ftarget,.TRUE.,found,line,error=error)
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
          CALL parser_get_object  (parser,Ftarget,newline=.TRUE.,error=error)
          CPPostcondition(TRIM(Ftarget)=="RADIUS",cp_failure_level,routineP,error,Failure)
          CALL parser_get_object  (parser,radius,error=error)
          CALL parser_get_object  (parser,units,error=error)
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
          CALL parser_get_object(parser,qmmm_gaussian_fns(J)%pgf%Number_of_Gaussians,newline=.TRUE.,error=error)
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
             CALL parser_get_object(parser,qmmm_gaussian_fns(IRad)%pgf%Ak(J),newline=.TRUE.,error=error)
             CALL parser_get_object(parser,qmmm_gaussian_fns(IRad)%pgf%Gk(J),error=error)
             CALL parser_get_object(parser,units,error=error)
             CALL uppercase (units)
             fconv = 1.0_dp
             IF (TRIM(units).EQ."ANGSTROM") fconv = bohr
             qmmm_gaussian_fns(IRad)%pgf%Ak(J) = qmmm_gaussian_fns(IRad)%pgf%Ak(J) * fconv
             qmmm_gaussian_fns(IRad)%pgf%Gk(J) = qmmm_gaussian_fns(IRad)%pgf%Gk(J) * fconv
          END DO
       ELSE
!       *** Stop program, if the end of file is reached ***
          IF (output_unit>0) WRITE(output_unit,'(A,F12.6,A)')"Discrepancies in Radius: ",Radius," definition."
          CALL stop_program(routineN,moduleN,__LINE__,"Radius Value not found in MM_POTENTIAL file",para_env)
       END IF

    END DO search_loop

    CALL parser_release(parser,error=error)

  END SUBROUTINE read_mm_potential

! *****************************************************************************
!> \brief set the GEEP information for Erf(r/rc)/r
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE set_mm_potential_erf(qmmm_gaussian_fns, qm_cell_small,&
                               compatibility, num_geep_gauss, error)
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: qmmm_gaussian_fns
    TYPE(cell_type), POINTER                 :: qm_cell_small
    LOGICAL, INTENT(IN)                      :: compatibility
    INTEGER, INTENT(IN)                      :: num_geep_gauss
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_mm_potential_erf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: IRad, istart, Nog, Nval, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_radius, radius, rc

    failure = .FALSE.
    Nval = SIZE(qmmm_gaussian_fns)
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
  END SUBROUTINE set_mm_potential_erf

! *****************************************************************************
!> \brief set the GEEP information for the S-WAVE expansion
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2007 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE set_mm_potential_swave(qmmm_gaussian_fns, qm_cell_small,&
                                    num_geep_gauss, error)
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: qmmm_gaussian_fns
    TYPE(cell_type), POINTER                 :: qm_cell_small
    INTEGER, INTENT(IN)                      :: num_geep_gauss
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_mm_potential_swave', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: IRad, istart, Nog, Nval, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: radius, rc

    failure = .FALSE.
    Nval = SIZE(qmmm_gaussian_fns)
    DO IRad = 1, Nval
       qmmm_gaussian_fns(IRad)%pgf%Number_of_Gaussians = num_geep_gauss
       radius = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius
       istart = 0
       ! Allocate Vectors
       SELECT CASE(num_geep_gauss)
       CASE(2)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s2_rc*bohr)
       CASE(3)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s3_rc*bohr)
       CASE(4)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s4_rc*bohr)
       CASE(5)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s5_rc*bohr)
       CASE(6)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s6_rc*bohr)
       CASE(7)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s7_rc*bohr)
       CASE(8)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s8_rc*bohr)
       CASE(9)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s9_rc*bohr)
       CASE(10)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s10_rc*bohr)
       CASE(11)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s11_rc*bohr)
       CASE(12)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s12_rc*bohr)
       CASE(13)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s13_rc*bohr)
       CASE(14)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s14_rc*bohr)
       CASE(15)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s15_rc*bohr)
       CASE(16)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s16_rc*bohr)
       CASE(17)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s17_rc*bohr)
       CASE(18)
          rc  = qmmm_gaussian_fns(Irad)%pgf%Elp_Radius / (s18_rc*bohr)
       END SELECT
       NOG = qmmm_gaussian_fns(IRad)%pgf%Number_of_Gaussians
       ALLOCATE(qmmm_gaussian_fns(IRad)%pgf%Ak(NOG),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       ALLOCATE(qmmm_gaussian_fns(IRad)%pgf%Gk(NOG),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       SELECT CASE(num_geep_gauss)
       CASE(2)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s2_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s2_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s2_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s2_b2
       CASE(3)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s3_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s3_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s3_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s3_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s3_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s3_b3
       CASE(4)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s4_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s4_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s4_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s4_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s4_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s4_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s4_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s4_b4
       CASE(5)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s5_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s5_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s5_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s5_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s5_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s5_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s5_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s5_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s5_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s5_b5
       CASE(6)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s6_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s6_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s6_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s6_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s6_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s6_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s6_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s6_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s6_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s6_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = s6_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = s6_b6
       CASE(7)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s7_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s7_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s7_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s7_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s7_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s7_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s7_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s7_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s7_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s7_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = s7_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = s7_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = s7_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = s7_b7
       CASE(8)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s8_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s8_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s8_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s8_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s8_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s8_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s8_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s8_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s8_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s8_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = s8_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = s8_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = s8_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = s8_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = s8_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = s8_b8
       CASE(9)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s9_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s9_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s9_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s9_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s9_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s9_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s9_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s9_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s9_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s9_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = s9_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = s9_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = s9_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = s9_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = s9_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = s9_b8
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = s9_a9
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = s9_b9
       CASE(10)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s10_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s10_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s10_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s10_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s10_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s10_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s10_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s10_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s10_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s10_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = s10_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = s10_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = s10_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = s10_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = s10_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = s10_b8
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = s10_a9
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = s10_b9
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= s10_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= s10_b10
       CASE(11)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s11_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s11_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s11_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s11_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s11_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s11_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s11_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s11_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s11_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s11_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = s11_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = s11_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = s11_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = s11_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = s11_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = s11_b8
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = s11_a9
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = s11_b9
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= s11_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= s11_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= s11_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= s11_b11
       CASE(12)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s12_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s12_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s12_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s12_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s12_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s12_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s12_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s12_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s12_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s12_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = s12_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = s12_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = s12_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = s12_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = s12_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = s12_b8
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = s12_a9
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = s12_b9
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= s12_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= s12_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= s12_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= s12_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= s12_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= s12_b12
       CASE(13)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s13_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s13_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s13_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s13_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s13_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s13_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s13_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s13_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s13_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s13_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = s13_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = s13_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = s13_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = s13_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = s13_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = s13_b8
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = s13_a9
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = s13_b9
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= s13_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= s13_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= s13_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= s13_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= s13_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= s13_b12
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+13)= s13_a13
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+13)= s13_b13
       CASE(14)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s14_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s14_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s14_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s14_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s14_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s14_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s14_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s14_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s14_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s14_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = s14_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = s14_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = s14_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = s14_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = s14_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = s14_b8
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = s14_a9
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = s14_b9
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= s14_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= s14_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= s14_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= s14_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= s14_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= s14_b12
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+13)= s14_a13
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+13)= s14_b13
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+14)= s14_a14
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+14)= s14_b14
       CASE(15)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s15_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s15_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s15_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s15_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s15_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s15_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s15_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s15_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s15_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s15_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = s15_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = s15_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = s15_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = s15_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = s15_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = s15_b8
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = s15_a9
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = s15_b9
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= s15_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= s15_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= s15_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= s15_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= s15_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= s15_b12
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+13)= s15_a13
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+13)= s15_b13
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+14)= s15_a14
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+14)= s15_b14
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+15)= s15_a15
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+15)= s15_b15
       CASE(16)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s16_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s16_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s16_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s16_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s16_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s16_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s16_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s16_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s16_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s16_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = s16_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = s16_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = s16_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = s16_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = s16_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = s16_b8
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = s16_a9
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = s16_b9
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= s16_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= s16_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= s16_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= s16_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= s16_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= s16_b12
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+13)= s16_a13
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+13)= s16_b13
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+14)= s16_a14
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+14)= s16_b14
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+15)= s16_a15
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+15)= s16_b15
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+16)= s16_a16
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+16)= s16_b16
       CASE(17)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s17_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s17_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s17_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s17_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s17_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s17_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s17_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s17_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s17_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s17_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = s17_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = s17_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = s17_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = s17_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = s17_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = s17_b8
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = s17_a9
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = s17_b9
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= s17_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= s17_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= s17_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= s17_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= s17_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= s17_b12
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+13)= s17_a13
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+13)= s17_b13
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+14)= s17_a14
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+14)= s17_b14
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+15)= s17_a15
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+15)= s17_b15
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+16)= s17_a16
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+16)= s17_b16
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+17)= s17_a17
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+17)= s17_b17
       CASE(18)
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1) = s18_a1
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1) = s18_b1
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+2) = s18_a2
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+2) = s18_b2
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+3) = s18_a3
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+3) = s18_b3
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+4) = s18_a4
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+4) = s18_b4
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+5) = s18_a5
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+5) = s18_b5
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+6) = s18_a6
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+6) = s18_b6
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+7) = s18_a7
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+7) = s18_b7
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+8) = s18_a8
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+8) = s18_b8
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+9) = s18_a9
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+9) = s18_b9
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+10)= s18_a10
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+10)= s18_b10
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+11)= s18_a11
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+11)= s18_b11
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+12)= s18_a12
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+12)= s18_b12
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+13)= s18_a13
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+13)= s18_b13
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+14)= s18_a14
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+14)= s18_b14
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+15)= s18_a15
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+15)= s18_b15
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+16)= s18_a16
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+16)= s18_b16
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+17)= s18_a17
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+17)= s18_b17
          qmmm_gaussian_fns(IRad)%pgf%Ak(istart+18)= s18_a18
          qmmm_gaussian_fns(IRad)%pgf%Gk(istart+18)= s18_b18
       END SELECT
       qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1:) = qmmm_gaussian_fns(IRad)%pgf%Ak(istart+1:) / rc
       qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1:) = qmmm_gaussian_fns(IRad)%pgf%Gk(istart+1:) * rc
    END DO
  END SUBROUTINE set_mm_potential_swave

END MODULE qmmm_gaussian_input

