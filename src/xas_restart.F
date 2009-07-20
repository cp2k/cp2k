!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Initialize the XAS orbitals for specific core excitations
!>       Either the GS orbitals are used as initial guess, or the
!>       xas mos are read from a previous calculation.
!>       In the latter case, the core-hole potetial should be the same.
!> \note
!>       The restart with the same core-hole potential should be checked
!>       and a wrong restart should stop the program
!> \par History
!>      created 09.2006
!> \author MI (09.2006)
! *****************************************************************************
MODULE xas_restart

  USE cp_control_types,                ONLY: dft_control_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_types,                     ONLY: cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_bcast
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             make_basis_sm
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             read_mos_restart_low,&
                                             wfn_restart_file_name,&
                                             write_mo_set_low
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE string_utilities,                ONLY: xstring
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_control,                     ONLY: xas_control_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_restart'

! *** Public subroutines ***

  PUBLIC ::  xas_read_restart, xas_write_restart, xas_initialize_rho

CONTAINS

! *****************************************************************************
!> \brief Set up for reading the restart
!>      corresponing to the excitation of iatom
!>      If the corresponding restart file does not exist
!>      the GS orbitals are used as initial guess
!> \param xas_control control variables
!> \param xas_section input section for XAS calculations
!>      qs_env:
!> \param iatom index of the absorbing atom
!> \param estate index of the core-hole orbital
!>      error:
!> \par History
!>      09.2006 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE xas_read_restart(xas_control,xas_section,qs_env,iatom,estate,error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(section_vals_type), POINTER         :: xas_section
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom
    INTEGER, INTENT(INOUT)                   :: estate
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_read_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: group, handle, homo, ia, ie, &
                                                ispin, istate, natom, nspin, &
                                                output_unit, rst_unit, source
    LOGICAL                                  :: failure, file_exists
    REAL(dp), DIMENSION(:), POINTER          :: occupation_numbers
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s

    CALL timeset(routineN,handle)

    failure = .FALSE.
    file_exists = .FALSE.

    NULLIFY(matrix_s,mos,occupation_numbers,para_env)
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    output_unit = cp_print_key_unit_nr(logger,xas_section,&
         "PRINT%PROGRAM_RUN_INFO", extension=".Log",error=error)

    CALL get_qs_env( qs_env=qs_env, para_env=para_env ,error=error)
    group = para_env%group
    source = para_env%source

    IF (para_env%ionode) THEN
      CALL wfn_restart_file_name(filename,file_exists,xas_section,logger,&
           xas=.TRUE.,error=error)

      CALL xstring (filename, ia, ie )
      filename = filename(ia:ie)//'-at'//&
                 TRIM(ADJUSTL(cp_to_string(iatom)))//'.rst'

      INQUIRE (FILE=filename,EXIST=file_exists)
      ! open file
      IF (file_exists) THEN

         CALL open_file(file_name=TRIM(filename),&
                     file_action="READ",&
                     file_form="UNFORMATTED",&
                     file_position="REWIND",&
                     file_status="OLD",&
                     unit_number=rst_unit)

         IF (output_unit > 0) WRITE (UNIT=output_unit,FMT="(/,T20,A,I5,/)")&
             "Read restart file for atom ", iatom

      ELSE IF (.NOT.file_exists) THEN
          IF(output_unit > 0) WRITE (UNIT=output_unit,FMT="(/,T10,A,I5,A,/)")&
               "Restart file for atom ", iatom,&
               " not available. Initialization done with GS orbitals"
      END IF
    END IF
    CALL mp_bcast(file_exists,source,group)

    IF(file_exists) THEN
      CALL get_qs_env(qs_env=qs_env, mos=mos, matrix_s=matrix_s,&
           para_env=para_env, particle_set=particle_set,error=error)

      natom = SIZE(particle_set,1)

      CALL read_mos_restart_low (mos, para_env, particle_set, natom, &
           rst_unit, xas=.TRUE., error=error)

      CALL get_mo_set(mo_set=mos(1)%mo_set,occupation_numbers=occupation_numbers,&
           homo=homo)
      estate = 0
      DO istate = 1,xas_control%nexc_search
        IF(occupation_numbers(istate)==xas_control%occ_estate) THEN
          estate = istate
          EXIT
        END IF
      END DO
      IF(estate==0) THEN
        CALL stop_program(routineP,"Core-hole orbital not found,"//&
        " probably the core-hole potential of the restart file is different")
      END IF
      occupation_numbers(homo) =  xas_control%occ_homo

      nspin = SIZE(mos,1)
      DO ispin = 1,nspin
      ! ortho so that one can restart for different positions (basis sets?)
         NULLIFY(mo_coeff)
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff,homo=homo)
         CALL make_basis_sm(mo_coeff,homo,matrix_s(1)%matrix,error=error)
      END DO
    END IF

    IF (para_env%ionode) THEN
      IF (file_exists) CALL close_file(unit_number=rst_unit)
    END IF

    CALL timestop(handle)

  END SUBROUTINE xas_read_restart

! *****************************************************************************
  SUBROUTINE xas_write_restart(xas_control,xas_section,qs_env,iatom,error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(section_vals_type), POINTER         :: xas_section
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_write_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    CHARACTER(LEN=default_string_length)     :: my_middle
    INTEGER                                  :: handle, output_unit, rst_unit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(mos,para_env,particle_set,logger,print_key)
    logger => cp_error_get_logger(error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
              xas_section,"PRINT%RESTART",used_print_key=print_key,error=error),&
              cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger,xas_section,&
            "PRINT%PROGRAM_RUN_INFO",extension=".Log",error=error)
       
       CALL get_qs_env( qs_env=qs_env, mos=mos, &
            para_env=para_env, particle_set=particle_set ,error=error)
       
       ! Open file
       rst_unit = -1
       my_middle= 'at'//TRIM(ADJUSTL(cp_to_string(iatom)))
       rst_unit = cp_print_key_unit_nr(logger,xas_section,"PRINT%RESTART",&
                  extension=".rst", file_status="REPLACE", file_action="WRITE",&
                  file_form="UNFORMATTED",middle_name=TRIM(my_middle), error=error)

       filename = cp_print_key_generate_filename(logger,print_key,&
                  middle_name=TRIM(my_middle),extension=".rst",&
                  my_local=.FALSE.,error=error)

       IF(output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,T10,A,I5,A,A,/)")&
               "Xas orbitals  for the absorbing atom ", iatom,&
               " are written in ", filename
          
       END IF

       ! Write mos
       CALL write_mo_set_low(mos, particle_set, rst_unit, error=error)
       ! Close file
       CALL cp_print_key_finished_output(rst_unit,logger,xas_section,&
                  "PRINT%RESTART", error=error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE xas_write_restart

!****f* xas_restart/xas_initialize_rho [1.0] *

! *****************************************************************************
!> \brief Once the mos and the occupation numbers are initialized
!>      the electronic density of the excited state can be calclated
!> \par History
!>      09-2006 MI created
!> \author MI
! *****************************************************************************
  SUBROUTINE xas_initialize_rho(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_initialize_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(dft_control,mos)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         mos=mos,rho=rho,error=error)

    DO ispin=1,SIZE(mos)
       CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set,&
            density_matrix=rho%rho_ao(ispin)%matrix,error=error)
    END DO

    CALL qs_rho_update_rho(rho,qs_env=qs_env, error=error)

    CALL timestop(handle)

  END SUBROUTINE xas_initialize_rho

END MODULE xas_restart
