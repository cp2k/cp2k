!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  I/O subroutines for pint_env
!> \author Lukasz Walewski
!> \date   2009-06-04
! *****************************************************************************
MODULE pint_io

  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE f77_interface,                   ONLY: f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get
  USE input_constants,                 ONLY: dump_atomic,&
                                             dump_dcd,&
                                             dump_xmol
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE pint_public,                     ONLY: pint_com_pos
  USE pint_transformations,            ONLY: pint_u2x
  USE pint_types,                      ONLY: e_conserved_id,&
                                             e_kin_thermo_id,&
                                             e_kin_virial_id,&
                                             e_potential_id,&
                                             pint_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_io'

  PUBLIC :: pint_write_line
  PUBLIC :: pint_write_centroids
  PUBLIC :: pint_write_trajectory
  PUBLIC :: pint_write_com
  PUBLIC :: pint_write_ener
  PUBLIC :: pint_write_step_info
  PUBLIC :: pint_write_rgyr

  CONTAINS

  ! ***************************************************************************
  !> \brief  Writes out a line of text to the default output unit.
  !> \author Lukasz Walewski
  !> \date   2009-07-10
  ! ***************************************************************************
  SUBROUTINE pint_write_line(line,error)

    CHARACTER(len=*), INTENT(IN)             :: line
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_write_line', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: my_label
    INTEGER                                  :: unit_nr
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    my_label = "PINT|"

    IF (logger%para_env%ionode) THEN
      unit_nr = cp_logger_get_default_unit_nr(logger)
      WRITE(unit_nr,'(T2,A)') TRIM(my_label)//" "//TRIM(line)
    END IF

    RETURN
  END SUBROUTINE pint_write_line


  ! ***************************************************************************
  !> \brief Write out the trajectory of the centroid (positions and velocities)
  !> \author fawzi
  !> \par History
  !>      various bug fixes - hforbert
  !>      2010-11-25 rewritten, added support for velocity printing,
  !>                 calc of the stddev of the beads turned off [lwalewski]
  ! ***************************************************************************
  SUBROUTINE pint_write_centroids(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_write_centroids', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: n_ids = 2, pos_id = 1, &
                                                vel_id = 2

    CHARACTER(len=default_string_length)     :: ext, form, my_middle_name, &
                                                unit_str
    CHARACTER(len=default_string_length), &
      DIMENSION(2)                           :: content_id, middle_name, &
                                                sect_path, title
    INTEGER                                  :: handle, handle1, iat, ib, id, &
                                                idim, idir, ierr, outformat, &
                                                should_output, unit_nr
    LOGICAL                                  :: failure, new_file
    REAL(kind=dp)                            :: nb, ss, unit_conv, vv
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle1)

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (failure) THEN
      CALL timestop(handle1)
      RETURN
    END IF

    sect_path(pos_id) = "MOTION%PINT%PRINT%CENTROID_POS"
    sect_path(vel_id) = "MOTION%PINT%PRINT%CENTROID_VEL"
    middle_name(pos_id) = "centroid-pos"
    middle_name(vel_id) = "centroid-vel"
    content_id(pos_id) = "POS"
    content_id(vel_id) = "VEL"
    WRITE(UNIT=title(pos_id),FMT="(A,I8,A,F20.10)")&
          " i =",pint_env%iter,&
          ", E =",SUM(pint_env%e_pot_bead)/REAL(pint_env%p,dp)
    WRITE(UNIT=title(vel_id),FMT="(A,I8,A,F20.10,A,F20.10)")&
          " i =",pint_env%iter,&
          ", E_trm =",pint_env%energy(e_kin_thermo_id),&
          ", E_vir =",pint_env%energy(e_kin_virial_id)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    CALL pint_u2x(pint_env,ux=pint_env%uv,x=pint_env%v,error=error)

    ! iterate over the properties that we know how to print
    ! (currently positions and velocities)
    DO id = 1, n_ids

      print_key => section_vals_get_subs_vals(pint_env%input,&
                   TRIM(sect_path(id)), error=error)

      should_output = cp_print_key_should_output(&
                      iteration_info=logger%iter_info,&
                      basis_section=print_key, error=error)
      IF ( .NOT. BTEST(should_output,cp_p_file) ) CONTINUE

      ! get units of measure for output (if available)
      CALL section_vals_val_get(print_key,"UNIT",&
           c_val=unit_str,error=error)
      unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)

      ! get the format for output
      CALL section_vals_val_get(print_key,"FORMAT",i_val=outformat,&
           error=error)

      SELECT CASE(outformat)
        CASE (dump_dcd)
          form="UNFORMATTED"
          ext=".dcd"
        CASE (dump_atomic)
          form="FORMATTED"
          ext=""
        CASE (dump_xmol)
          form="FORMATTED"
          ext=".xyz"
        CASE default
          CPPostcondition(.FALSE.,cp_fatal_level,routineP,error,failure)
      END SELECT

      NULLIFY(f_env,cell,subsys)
      CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id,&
           f_env=f_env,new_error=new_error, failure=failure, handle=handle)
      CALL force_env_get(force_env=f_env%force_env,&
           cell=cell, subsys=subsys, error=new_error)
      CALL cp_subsys_get(subsys,particles=particles,&
           error=new_error)

      ! calculate and copy the requested property
      ! to the particles structure
      nb=REAL(pint_env%p,dp)
      idim=0
      DO iat=1,pint_env%ndim/3
        DO idir=1,3
          idim=idim+1
          ss=0.0_dp
          vv=0.0_dp
!          ss2=0.0_dp
          DO ib=1,pint_env%p
            ss=ss+pint_env%x(ib,idim)
            vv=vv+pint_env%v(ib,idim)
!            ss2=ss2+pint_env%x(ib,idim)**2
          END DO
          particles%els(iat)%r(idir)=ss/nb
          particles%els(iat)%v(idir)=vv/nb
!          particles%els(iat)%v(idir)=SQRT(ss2/nb-(ss/nb)**2)
        END DO
      END DO

      ! set up the output unit number and file name
      ! for the current property
      my_middle_name = TRIM(middle_name(id))
      unit_nr = cp_print_key_unit_nr(logger=logger, &
                basis_section=print_key, print_key_path="", &
                extension=TRIM(ext), middle_name=TRIM(my_middle_name), &
                local=.FALSE., file_form=form, is_new_file=new_file, &
                error=new_error)

      ! don't write the 0-th frame if the file already exists
      IF ( .NOT. new_file .AND. (pint_env%iter .LE. pint_env%first_step) ) THEN
        CALL cp_print_key_finished_output(unit_nr,logger,&
             print_key, error=new_error)
        CONTINUE
      END IF

      ! actually perform the i/o - on the ionode only
      IF (unit_nr>0) THEN

        CALL write_particle_coordinates( &
             particles%els, &
             iunit=unit_nr, &
             output_format=outformat, &
             content=content_id(id), &
             title=title(id), &
             cell=cell, &
             unit_conv=unit_conv, &
             error=error)

        CALL cp_print_key_finished_output(unit_nr,logger,&
             print_key, "", local=.FALSE., error=new_error)

      END IF

      CALL f_env_rm_defaults(f_env,new_error,ierr,handle)
      CPAssert(ierr==0,cp_failure_level,routineP,error,failure)

    END DO

    CALL timestop(handle1)
    RETURN
  END SUBROUTINE pint_write_centroids

  ! ***************************************************************************
  !> \brief  Write out the trajectory of the beads (positions and velocities)
  !> \author hforbert
  !> \par    History
  !>         2010-11-25 added support for velocity printing [lwalewski]
  ! ***************************************************************************
  SUBROUTINE pint_write_trajectory(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_write_trajectory', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: n_ids = 2, pos_id = 1, &
                                                vel_id = 2

    CHARACTER(len=default_string_length)     :: ext, form, ib_str, &
                                                my_middle_name, title, &
                                                unit_str
    CHARACTER(len=default_string_length), &
      DIMENSION(2)                           :: content_id, middle_name, &
                                                sect_path
    INTEGER                                  :: handle, handle1, iat, ib, id, &
                                                idim, idir, ierr, outformat, &
                                                should_output, unit_nr
    LOGICAL                                  :: failure, new_file
    REAL(kind=dp)                            :: unit_conv
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle1)

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (failure) THEN
      CALL timestop(handle1)
      RETURN
    END IF

    sect_path(pos_id) = "MOTION%PRINT%TRAJECTORY"
    sect_path(vel_id) = "MOTION%PRINT%VELOCITIES"
    middle_name(pos_id) = "pos-"
    middle_name(vel_id) = "vel-"
    content_id(pos_id) = "POS"
    content_id(vel_id) = "VEL"

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    CALL pint_u2x(pint_env,ux=pint_env%uv,x=pint_env%v,error=error)

    ! iterate over the properties that we know how to print
    ! (currently positions and velocities)
    DO id = 1, n_ids

      print_key => section_vals_get_subs_vals(pint_env%input,&
                   TRIM(sect_path(id)), error=error)

      should_output = cp_print_key_should_output(&
                      iteration_info=logger%iter_info,&
                      basis_section=print_key, error=error)
      IF ( .NOT. BTEST(should_output,cp_p_file) ) CONTINUE

      ! get units of measure for output (if available)
      CALL section_vals_val_get(print_key,"UNIT",&
           c_val=unit_str,error=error)
      unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)

      ! get the format for output
      CALL section_vals_val_get(print_key,"FORMAT",i_val=outformat,&
           error=error)

      SELECT CASE(outformat)
        CASE (dump_dcd)
          form="UNFORMATTED"
          ext=".dcd"
        CASE (dump_atomic)
          form="FORMATTED"
          ext=""
        CASE (dump_xmol)
          form="FORMATTED"
          ext=".xyz"
        CASE default
          CPPostcondition(.FALSE.,cp_fatal_level,routineP,error,failure)
      END SELECT

      NULLIFY(f_env,cell,subsys)
      CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id,&
           f_env=f_env,new_error=new_error, failure=failure, handle=handle)
      CALL force_env_get(force_env=f_env%force_env,&
           cell=cell, subsys=subsys, error=new_error)
      CALL cp_subsys_get(subsys,particles=particles,&
           error=new_error)

      ! iterate over beads
      DO ib = 1, pint_env%p

        ! copy the requested property of the current bead
        ! to the particles structure
        idim=0
        DO iat=1,pint_env%ndim/3
          DO idir=1,3
            idim=idim+1
            particles%els(iat)%r(idir)=pint_env%x(ib,idim)
            particles%els(iat)%v(idir)=pint_env%v(ib,idim)
          END DO
        END DO

        ! set up the output unit number and file name
        ! for the current property and bead
        ib_str = ""
        WRITE(ib_str,*) ib
        my_middle_name = TRIM(middle_name(id)) // TRIM(ADJUSTL(ib_str))
        unit_nr = cp_print_key_unit_nr(logger=logger, &
                  basis_section=print_key, print_key_path="", &
                  extension=TRIM(ext), middle_name=TRIM(my_middle_name), &
                  local=.FALSE., file_form=form, is_new_file=new_file, &
                  error=new_error)

        ! don't write the 0-th frame if the file already exists
        IF ( .NOT. new_file .AND. (pint_env%iter .LE. pint_env%first_step) ) THEN
          CALL cp_print_key_finished_output(unit_nr,logger,&
               print_key, error=new_error)
          CONTINUE
        END IF

        ! actually perform the i/o - on the ionode only
        IF (unit_nr>0) THEN

          IF (outformat == dump_xmol) THEN
            WRITE(UNIT=title,FMT="(A,I8,A,F20.10)")&
                  " i =",pint_env%iter,&
                  ", E =",pint_env%e_pot_bead(ib)
          END IF

          CALL write_particle_coordinates( &
               particles%els, &
               iunit=unit_nr, &
               output_format=outformat, &
               content=content_id(id), &
               title=title, &
               cell=cell, &
               unit_conv=unit_conv, &
               error=error)

          CALL cp_print_key_finished_output(unit_nr,logger,&
               print_key, "", local=.FALSE., error=new_error)

        END IF

      END DO

      CALL f_env_rm_defaults(f_env,new_error,ierr,handle)
      CPAssert(ierr==0,cp_failure_level,routineP,error,failure)

    END DO

    CALL timestop(handle1)
    RETURN
  END SUBROUTINE pint_write_trajectory

  ! ***************************************************************************
  !> \brief  Write center of mass (COM) position according to PINT%PRINT%COM
  !> \author Lukasz Walewski
  !> \date   2010-02-17
  ! ***************************************************************************
  SUBROUTINE pint_write_com( pint_env, error )

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_write_com', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: stmp1, stmp2
    INTEGER                                  :: ic, unit_nr
    LOGICAL                                  :: failure, new_file, &
                                                should_output
    REAL(kind=dp), DIMENSION(3)              :: com_r
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    ! decide whether to write anything or not
    NULLIFY(print_key)
    print_key => section_vals_get_subs_vals(pint_env%input, &
      "MOTION%PINT%PRINT%COM", error=error)
    should_output = BTEST(cp_print_key_should_output( &
      iteration_info=logger%iter_info,&
      basis_section=print_key,error=error),cp_p_file)
    IF ( .NOT. should_output ) THEN
      RETURN
    END IF

    com_r = pint_com_pos( pint_env, error )
    DO ic = 1, 3
      com_r(ic) = cp_unit_from_cp2k(com_r(ic), "angstrom", error=error)
    END DO

    unit_nr=cp_print_key_unit_nr(logger, print_key, is_new_file=new_file,&
            middle_name="com-pos",extension=".xyz",error=error)

    ! don't write the 0-th frame if the file already exists
    IF ( .NOT. new_file .AND. (pint_env%iter .LE. pint_env%first_step) ) THEN
      CALL cp_print_key_finished_output(unit_nr,logger,&
           print_key, error=error)
      RETURN
    END IF

    ! actually perform the i/o - on the ionode only
    IF (unit_nr>0) THEN

      WRITE(unit_nr,'(I2)') 1
      WRITE(stmp1,*) pint_env%iter
      WRITE(stmp2,'(F20.10)') pint_env%energy(e_conserved_id)
      WRITE(unit_nr,'(4A)') " Iteration = ", TRIM(ADJUSTL(stmp1)), &
                            ", E_conserved = ", TRIM(ADJUSTL(stmp2))
      WRITE(unit_nr,'(A2,3(1X,F20.10))') "X ", (com_r(ic),ic=1,3)

      CALL m_flush(unit_nr)

    END IF

    CALL cp_print_key_finished_output(unit_nr,logger,print_key,error=error)

    RETURN
  END SUBROUTINE pint_write_com

  ! ***************************************************************************
  !> \brief  Writes out the energies according to PINT%PRINT%ENERGY
  !> \author fawzi
  !> \param  pint_env path integral environment
  !> \param  error variable to control error logging, stopping,...
  !>         see module cp_error_handling
  !> \par    History
  !>           various bug fixes [hforbert]
  !>           2009-11-16 energy components calc moved out of here [lwalewski]
  ! ***************************************************************************
  SUBROUTINE pint_write_ener(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_write_ener', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ndof, unit_nr
    LOGICAL                                  :: failure, file_is_new
    REAL(kind=dp)                            :: t, temp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF ( failure ) RETURN

    NULLIFY(print_key,logger)
    print_key => section_vals_get_subs_vals(pint_env%input, &
      "MOTION%PINT%PRINT%ENERGY", error=error)
    logger => cp_error_get_logger(error)
    IF ( BTEST(cp_print_key_should_output(iteration_info=logger%iter_info,&
      basis_section=print_key,error=error),cp_p_file) ) THEN

      unit_nr=cp_print_key_unit_nr(logger,print_key,middle_name="energy",&
        extension=".dat",is_new_file=file_is_new,error=error )

      ! don't write the 0-th frame if the file already exists
      IF ( .NOT. file_is_new .AND. (pint_env%iter .LE. pint_env%first_step) ) THEN
        CALL cp_print_key_finished_output(unit_nr,logger,&
             print_key, error=error)
        RETURN
      END IF

      ! cp_print_key_unit_nr returns -1 on nodes other than logger%para_env%ionode
      IF (unit_nr>0) THEN

        ! please keep the format explaination up to date
        ! keep the constant of motion the true constant of motion !
        IF ( file_is_new ) THEN
          WRITE(unit_nr,"(A8,1X,A12,1X,5(A20,1X),A12)")&
                            "# StepNr",&
                        "   Time [fs]",&
                "      Kinetic [a.u.]",&
                "    VirialKin [a.u.]",&
                "     Temperature [K]",&
                "    Potential [a.u.]",&
                "      ConsQty [a.u.]",&
                        "     CPU [s]"
        END IF

        t=cp_unit_from_cp2k(pint_env%t,"fs",error=error)

        ndof = pint_env%p
        IF ( pint_env%first_propagated_mode .EQ. 2 ) THEN
          ndof = ndof - 1
        END IF
        temp = cp_unit_from_cp2k(2.0_dp*pint_env%e_kin_beads/&
               REAL(ndof,dp)/REAL(pint_env%ndim,dp),&
               "K",error=error)

        WRITE (unit_nr,"(I8,1X,F12.3,1X,5(F20.9,1X),F12.1)")&
               pint_env%iter,&
               t,&
               pint_env%energy(e_kin_thermo_id),&
               pint_env%energy(e_kin_virial_id),&
               temp,&
               pint_env%energy(e_potential_id),&
               pint_env%energy(e_conserved_id),&
               pint_env%time_per_step
        CALL m_flush(unit_nr)

      END IF

      CALL cp_print_key_finished_output(unit_nr,logger,print_key,error=error)
    END IF

    RETURN
  END SUBROUTINE pint_write_ener


  ! ***************************************************************************
  !> \brief  Write step info to the output file.
  !> \date   2009-11-16
  !> \author Lukasz Walewski
  !> \par History
  !>      2010-01-27 getting default unit nr now only on ionode [lwalewski]
  ! ***************************************************************************
  SUBROUTINE pint_write_step_info( pint_env, error )
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_write_step_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: msgstr, stmp, time_unit
    INTEGER                                  :: unit_nr
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: time_used
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    time_used = pint_env%time_per_step
    time_unit = "sec"
    IF ( time_used .GE. 60.0_dp ) THEN
      time_used = time_used / 60.0_dp
      time_unit = "min"
    END IF
    IF ( time_used .GE. 60.0_dp ) THEN
      time_used = time_used / 60.0_dp
      time_unit = "hours"
    END IF
    msgstr = "MD step"
    stmp = ""
    WRITE(stmp, *) pint_env%iter
    msgstr = TRIM(ADJUSTL(msgstr)) // " " // TRIM(ADJUSTL(stmp)) // " of"
    stmp = ""
    WRITE(stmp, *) pint_env%last_step
    msgstr = TRIM(ADJUSTL(msgstr)) // " " // TRIM(ADJUSTL(stmp)) // " in"
    stmp = ""
    WRITE(stmp, '(F20.1)') time_used
    msgstr = TRIM(ADJUSTL(msgstr)) // " " // TRIM(ADJUSTL(stmp))
    msgstr = TRIM(ADJUSTL(msgstr)) // " " // TRIM(ADJUSTL(time_unit)) // "."

    IF (logger%para_env%ionode) THEN
      unit_nr = cp_logger_get_default_unit_nr(logger)
      WRITE(unit_nr,'(T2,A)') "PINT| " // TRIM(ADJUSTL(msgstr))
    END IF

    ! print out the total energy - for regtest evaluation
    stmp = ""
    WRITE(stmp, *) pint_env%energy(e_conserved_id)
    msgstr = "Total energy = " // TRIM(ADJUSTL(stmp))
    IF (logger%para_env%ionode) THEN
      WRITE(unit_nr,'(T2,A)') "PINT| " // TRIM(ADJUSTL(msgstr))
    END IF

    RETURN
  END SUBROUTINE pint_write_step_info

  ! ***************************************************************************
  !> \brief  Write radii of gyration according to PINT%PRINT%CENTROID_GYR
  !> \author Lukasz Walewski
  !> \date   2011-01-07
  ! ***************************************************************************
  SUBROUTINE pint_write_rgyr( pint_env, error )

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_write_rgyr', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: unit_str
    INTEGER                                  :: ia, ib, ic, idim, unit_nr
    LOGICAL                                  :: failure, new_file, &
                                                should_output
    REAL(kind=dp)                            :: nb, ss, unit_conv
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    ! decide whether to write anything or not
    NULLIFY(print_key)
    print_key => section_vals_get_subs_vals(pint_env%input, &
      "MOTION%PINT%PRINT%CENTROID_GYR", error=error)
    should_output = BTEST(cp_print_key_should_output( &
      iteration_info=logger%iter_info,&
      basis_section=print_key,error=error),cp_p_file)
    IF ( .NOT. should_output ) THEN
      RETURN
    END IF

    ! get the units conversion factor
    CALL section_vals_val_get(print_key,"UNIT",c_val=unit_str,error=error)
    unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)

    ! calculate the centroid positions
    nb=REAL(pint_env%p,dp)
    idim=0
    DO ia = 1, pint_env%ndim/3
      DO ic = 1, 3
        idim=idim+1
        ss=0.0_dp
        DO ib = 1, pint_env%p
          ss=ss+pint_env%x(ib,idim)
        END DO
        pint_env%rtmp_ndim(idim) = ss/nb
      END DO
    END DO

    ! calculate the radii of gyration
    idim=0
    DO ia = 1, pint_env%ndim/3
      ss=0.0_dp
      DO ic = 1, 3
        idim=idim+1
        DO ib = 1, pint_env%p
          ss=ss+(pint_env%x(ib,idim)-pint_env%rtmp_ndim(idim))**2
        END DO
      END DO
      pint_env%rtmp_natom(ia) = SQRT(ss/nb) * unit_conv
    END DO

    unit_nr=cp_print_key_unit_nr(logger, print_key, is_new_file=new_file, &
            middle_name="centroid-gyr",extension=".dat",error=error)

    ! don't write the 0-th frame if the file already exists
    IF ( .NOT. new_file .AND. (pint_env%iter .LE. pint_env%first_step) ) THEN
      CALL cp_print_key_finished_output(unit_nr,logger,&
           print_key, error=error)
      RETURN
    END IF

    ! actually perform the i/o - on the ionode only
    IF (unit_nr>0) THEN

      DO ia = 1, pint_env%ndim/3
        WRITE(unit_nr,'(F20.10,1X)',ADVANCE='NO') pint_env%rtmp_natom(ia)
      END DO
      WRITE(unit_nr,'(A)') ""

      CALL m_flush(unit_nr)

    END IF

    CALL cp_print_key_finished_output(unit_nr,logger,print_key,error=error)

    RETURN
  END SUBROUTINE pint_write_rgyr

END MODULE pint_io
