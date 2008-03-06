!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Output Utilities for MOTION_SECTION
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   02.2008 
! *****************************************************************************
MODULE motion_utils

  USE cell_types,                      ONLY: cell_type
  USE cp2k_units,                      ONLY: unit_convert_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_constants,                 ONLY: dump_atomic,&
                                             dump_dcd,&
                                             dump_xmol
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             sp
  USE machine,                         ONLY: m_flush
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_coordinates
  USE physcon,                         ONLY: angstrom
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'motion_utils'

  PUBLIC :: write_trajectory, write_stress_tensor, write_simulation_cell,&
            get_output_format

CONTAINS

! *****************************************************************************
!> \brief   Prints the information controlled by the TRAJECTORY section
!> \author  Teodoro Laino [tlaino] - University of Zurich
!> \date    02.2008
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_trajectory (force_env, root_section, it, time, etot, pk_name,&
       pos, act, middle_name, particles, extended_xmol_title, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    INTEGER, INTENT(IN)                      :: it
    REAL(KIND=dp), INTENT(IN)                :: time, etot
    CHARACTER(LEN=*), OPTIONAL               :: pk_name
    CHARACTER(LEN=default_string_length), &
      OPTIONAL                               :: pos, act
    CHARACTER(LEN=*), OPTIONAL               :: middle_name
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: particles
    LOGICAL, INTENT(IN), OPTIONAL            :: extended_xmol_title
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_trajectory', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=4)                         :: id_dcd
    CHARACTER(LEN=default_string_length) :: id_label, id_wpc, my_act, my_ext, &
      my_form, my_middle, my_pk_name, my_pos, remark1, remark2, section_ref, &
      title
    INTEGER                                  :: iskip, nat, outformat, &
                                                traj_unit
    LOGICAL                                  :: failure, &
                                                my_extended_xmol_title, &
                                                new_file
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: my_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.
    NULLIFY(logger, cell, subsys, my_particles, particle_set)
    logger => cp_error_get_logger(error)
    id_label = logger%iter_info%level_name(logger%iter_info%n_rlevel)
    my_pos     = "APPEND"
    my_act     = "WRITE"
    my_middle  = "pos"
    my_pk_name = "TRAJECTORY"
    IF (PRESENT(middle_name))  my_middle  = middle_name
    IF (PRESENT(pos))          my_pos     = pos
    IF (PRESENT(act))          my_act     = act
    IF (PRESENT(pk_name))      my_pk_name = pk_name
    
    SELECT CASE(TRIM(my_pk_name))
    CASE("TRAJECTORY","SHELL_TRAJECTORY","CORE_TRAJECTORY")
       id_dcd = "CORD"
       id_wpc = "POS"
    CASE("VELOCITIES","SHELL_VELOCITIES","CORE_VELOCITIES")
       id_dcd = "VEL "
       id_wpc = "VEL"
    CASE("FORCES","SHELL_FORCES","CORE_FORCES")
       id_dcd = "FRC "
       id_wpc = "FORCE"
    CASE DEFAULT
       CPPostcondition(.FALSE.,cp_fatal_level,routineP,error,failure)
    END SELECT

    CALL force_env_get(force_env, cell=cell, subsys=subsys, error=error)
    IF (PRESENT(particles)) THEN
       CPPostcondition(ASSOCIATED(particles),cp_fatal_level,routineP,error,failure)
       my_particles => particles
    ELSE
       CALL cp_subsys_get(subsys=subsys, particles=my_particles,error=error)
    END IF
    particle_set => my_particles%els
    nat=my_particles%n_els

    CALL get_output_format(root_section,"MOTION%PRINT%"//TRIM(my_pk_name),my_form,my_ext,error)
    traj_unit = cp_print_key_unit_nr(logger,root_section,"MOTION%PRINT%"//TRIM(my_pk_name),&
            extension=my_ext, file_position=my_pos, file_action=my_act,&
            file_form=my_form, middle_name=TRIM(my_middle), is_new_file=new_file,&
            error=error)
    IF (traj_unit > 0) THEN
       CALL section_vals_val_get(root_section,"MOTION%PRINT%"//TRIM(my_pk_name)//"%FORMAT",&
            i_val=outformat,error=error)
       title = ""
       IF  (outformat==dump_dcd) THEN
          IF (new_file) THEN
             !Lets write the header for the coordinate dcd
             section_ref = "MOTION%PRINT%"//TRIM(my_pk_name)//"%EACH%"//TRIM(id_label)
             iskip = section_get_ival(root_section,TRIM(section_ref),error=error)
             WRITE(traj_unit) id_dcd,0,it,iskip,0,0,0,0,0,0,REAL(time,KIND=sp),1,0,0,0,0,0,0,0,0,24
             remark1= "REMARK FILETYPE "//id_dcd//" DCD GENERATED BY CP2K"
             remark2= "REMARK Support new DCD format with cell information"
             WRITE(traj_unit) 2,remark1,remark2
             WRITE(traj_unit) nat
             CALL m_flush(traj_unit)
          ENDIF
       ELSE IF(outformat==dump_xmol) THEN
          my_extended_xmol_title = .FALSE.
          IF (PRESENT(extended_xmol_title)) my_extended_xmol_title = extended_xmol_title
          ! This information can be digested by Molden
          IF (my_extended_xmol_title) THEN
            WRITE (UNIT=title,FMT="(A,I8,A,F12.3,A,F20.10)")&
              " i = ",it,", time = ",time*REAL(it,KIND=dp),", E = ",etot
          ELSE
            WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")" i = ",it,", E = ",etot
          END IF
       ELSE IF(outformat==dump_atomic) THEN
          ! Do Nothing..
       ELSE
          CPPostcondition(.FALSE.,cp_fatal_level,routineP,error,failure)
       END IF
       CALL write_particle_coordinates(particle_set,traj_unit,outformat,TRIM(id_wpc),TRIM(title),cell)
    END IF
    CALL cp_print_key_finished_output(traj_unit,logger,root_section,"MOTION%PRINT%"//TRIM(my_pk_name),error=error)

  END SUBROUTINE write_trajectory

! *****************************************************************************
!> \brief Info on the unit to be opened to dump MD informations
!> \author Teodoro Laino - University of Zurich - 07.2007
! *****************************************************************************
  SUBROUTINE get_output_format(section,path,my_form,my_ext,error)
    TYPE(section_vals_type), POINTER         :: section
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: path
    CHARACTER(LEN=*), INTENT(OUT)            :: my_form, my_ext
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_output_format', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: output_format

    my_form = "FORMATTED"
    my_ext  = ".xyz"
    IF (PRESENT(path)) THEN
       CALL section_vals_val_get(section,TRIM(path)//"%FORMAT",i_val=output_format,error=error)
    ELSE
       CALL section_vals_val_get(section,"FORMAT",i_val=output_format,error=error)
    END IF
    IF(output_format==dump_dcd) THEN
       my_form = "UNFORMATTED"
       my_ext  = ".dcd"
    END IF
  END SUBROUTINE get_output_format

! *****************************************************************************
!> \brief   Prints the Stress Tensor
!> \author  Teodoro Laino [tlaino] - University of Zurich
!> \date    02.2008
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_stress_tensor(virial, cell, motion_section, unit, itimes, time, pos, act, error)
    TYPE(virial_type), POINTER               :: virial
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(unit_convert_type), POINTER         :: unit
    INTEGER, INTENT(IN)                      :: itimes
    REAL(KIND=dp), INTENT(IN)                :: time
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: pos, act
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_stress_tensor', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_act, my_pos
    INTEGER                                  :: str
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    
    IF( virial % pv_availability ) THEN
      my_pos = "APPEND"
      my_act = "WRITE"
      IF (PRESENT(pos)) my_pos = pos
      IF (PRESENT(act)) my_act = act
      str = cp_print_key_unit_nr(logger,motion_section,"PRINT%STRESS",&
           extension=".stress", file_position=my_pos, file_action=my_act,&
           file_form="FORMATTED", error=error)
    ELSE
      str = 0
    END IF
    IF ( str > 0 ) THEN
       WRITE(UNIT=str,FMT="(I10,F12.3,9F20.10)") itimes,time,&
                     virial%pv_total(1,1)/cell%deth*unit%bar,&
                     virial%pv_total(1,2)/cell%deth*unit%bar,&
                     virial%pv_total(1,3)/cell%deth*unit%bar,&
                     virial%pv_total(2,1)/cell%deth*unit%bar,&
                     virial%pv_total(2,2)/cell%deth*unit%bar,&
                     virial%pv_total(2,3)/cell%deth*unit%bar,& 
                     virial%pv_total(3,1)/cell%deth*unit%bar,&
                     virial%pv_total(3,2)/cell%deth*unit%bar,&
                     virial%pv_total(3,3)/cell%deth*unit%bar
       CALL m_flush(str)
    END IF
    IF( virial % pv_availability ) THEN
      CALL cp_print_key_finished_output(str,logger,motion_section,&
           "PRINT%STRESS", error=error)
    END IF
  END SUBROUTINE write_stress_tensor

! *****************************************************************************
!> \brief   Prints the Simulation Cell
!> \author  Teodoro Laino [tlaino] - University of Zurich
!> \date    02.2008
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_simulation_cell(cell, motion_section, itimes, time, pos, act, error)
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: motion_section
    INTEGER, INTENT(IN)                      :: itimes
    REAL(KIND=dp), INTENT(IN)                :: time
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: pos, act
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_simulation_cell', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_act, my_pos
    INTEGER                                  :: box
    LOGICAL                                  :: new_file
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    
    my_pos = "APPEND"
    my_act = "WRITE"
    IF (PRESENT(pos)) my_pos = pos
    IF (PRESENT(act)) my_act = act
    box = cp_print_key_unit_nr(logger,motion_section,"PRINT%CELL",&
         extension=".cell", file_position=my_pos, file_action=my_act,&
         file_form="FORMATTED", is_new_file=new_file, error=error)
    IF ( box > 0 ) THEN
       IF (new_file) THEN
          WRITE (box,'("#",3X,A,2X,A,2X,9(6X,A,6X),4X,A)')"Step Nr.","Time[fs]","XX[Ang.]","YX[Ang.]",&
               "ZX[Ang.]","YX[Ang.]","YY[Ang.]","YZ{Ang.]","ZX[Ang.]","ZY[Ang.]","ZZ[Ang.]","Vol[Ang.^3]"
       END IF
       WRITE(UNIT=box,FMT="(I8, F12.3,10F20.10)") itimes,time,&
            cell%hmat(1,1)*angstrom,cell%hmat(2,1)*angstrom,cell%hmat(3,1)*angstrom,&
            cell%hmat(1,2)*angstrom,cell%hmat(2,2)*angstrom,cell%hmat(3,2)*angstrom,&
            cell%hmat(1,3)*angstrom,cell%hmat(2,3)*angstrom,cell%hmat(3,3)*angstrom,&
            cell%deth*angstrom*angstrom*angstrom
       CALL m_flush(box)
    END IF
    CALL cp_print_key_finished_output(box,logger,motion_section,"PRINT%CELL", error=error)

  END SUBROUTINE write_simulation_cell

END MODULE motion_utils
