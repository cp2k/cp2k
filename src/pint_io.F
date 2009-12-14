!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
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
  USE input_section_types,             ONLY: section_vals_add_values,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set,&
                                             section_vals_write
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE pint_transformations,            ONLY: pint_u2x
  USE pint_types,                      ONLY: e_conserved_id,&
                                             e_kin_thermo_id,&
                                             e_kin_virial_id,&
                                             e_potential_id,&
                                             pint_env_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_io'

  PUBLIC :: pint_write_centroids
  PUBLIC :: pint_write_trajectory
  PUBLIC :: pint_write_ener
  PUBLIC :: pint_write_restart
  PUBLIC :: pint_write_step_info

  CONTAINS

  ! ***************************************************************************
  !> \brief Writes out the trajectory of the centroids (and their spread)
  !> \author fawzi
  !> \param pint_env path integral environment
  !> \param error variable to control error logging, stopping,...
  !>        see module cp_error_handling
  !> \par History
  !>      various bug fixes - hforbert
  ! ***************************************************************************
  SUBROUTINE pint_write_centroids(pint_env,error)

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_write_centroids', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: ext, form, my_middle_name, &
                                                title, unit_str
    INTEGER                                  :: handle, iat, ib, idim, idir, &
                                                ierr, outformat, unit_nr
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: nb, ss, ss2, unit_conv
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: my_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: print_key

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF ( failure ) RETURN

    NULLIFY(print_key,logger,f_env,cell,subsys)
    print_key => section_vals_get_subs_vals(pint_env%input, &
      "MOTION%PINT%PRINT%CENTROIDS", error=error)
    logger => cp_error_get_logger(error)
    IF ( BTEST(cp_print_key_should_output(iteration_info=logger%iter_info,&
      basis_section=print_key,error=error),cp_p_file) ) THEN

      ! Gather units of measure for output (if available)
      CALL section_vals_val_get(print_key,"UNIT",&
           c_val=unit_str,error=error)
      unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)

      ! Get the format for output
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
      my_middle_name="centroid"

      unit_nr=cp_print_key_unit_nr(logger, print_key, &
              middle_name=TRIM(my_middle_name), extension=TRIM(ext), &
              local=.FALSE., file_form=form, error=error)

      ! cp_print_key_unit_nr returns -1 on nodes other than logger%para_env%ionode
      IF (unit_nr>0) THEN

        ! FIXME!!!!
        ! This call should be splitted in one for position and in one for velocities
        ! at the moment it will print only positions. Please, fix it if really used..

        CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id,&
             f_env=f_env,new_error=my_error, failure=failure, handle=handle)
        CALL force_env_get(force_env=f_env%force_env,cell=cell,subsys=subsys,error=my_error)
        CALL cp_subsys_get(subsys,particles=particles,error=my_error)

        nb=REAL(pint_env%p,dp)
        idim=0
        DO iat=1,pint_env%ndim/3
          DO idir=1,3
            idim=idim+1
            ss=0._dp
            ss2=0._dp
            DO ib=1,pint_env%p
              ss=ss+pint_env%x(ib,idim)
              ss2=ss2+pint_env%x(ib,idim)**2
            END DO
            particles%els(iat)%r(idir)=ss/nb
            particles%els(iat)%v(idir)=SQRT(ss2/nb-(ss/nb)**2)
          END DO
        END DO

        IF (outformat == dump_xmol) THEN
          WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")&
                " i =",pint_env%iter,&
                ", E =",SUM(pint_env%e_pot_bead)/REAL(pint_env%p,dp)
        END IF

        CALL write_particle_coordinates(particles%els,&
             iunit=unit_nr, output_format=outformat, content="POS",&
             title=title, cell=cell, unit_conv=unit_conv, error=error)
        CALL m_flush(unit_nr)

        CALL f_env_rm_defaults(f_env,my_error,ierr,handle)
        CPAssert(ierr==0,cp_failure_level,routineP,error,failure)

      END IF

      CALL cp_print_key_finished_output( unit_nr, logger, print_key, &
           error=error )
    END IF

    RETURN
  END SUBROUTINE pint_write_centroids

  ! ***************************************************************************
  !> \brief writes out the trajectory of the beads
  !> \param pint_env path integral environment
  !> \param error variable to control error logging, stopping,...
  !>        see module cp_error_handling
  !> \author hforbert
  ! ***************************************************************************
  SUBROUTINE pint_write_trajectory(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_write_trajectory', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: ext, form, ib_str, &
                                                my_middle_name, title, &
                                                unit_str
    INTEGER                                  :: handle, iat, ib, idim, idir, &
                                                ierr, outformat, unit_nr
    LOGICAL                                  :: failure, new_file
    REAL(kind=dp)                            :: unit_conv
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: print_key

    failure=.FALSE.

    NULLIFY(f_env,cell,subsys)
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
     CALL pint_u2x(pint_env,ux=pint_env%uv,x=pint_env%v,error=error)
     IF (pint_env%replicas%para_env%ionode) THEN

        print_key => section_vals_get_subs_vals(pint_env%input,&
           "MOTION%PRINT%TRAJECTORY", error=error)

        ! Gather units of measure for output (if available)
        CALL section_vals_val_get(print_key,"UNIT",&
           c_val=unit_str,error=error)
        unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)

        ! Get the format for output
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

        CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id,&
             f_env=f_env,new_error=new_error, failure=failure, handle=handle)
        logger => cp_error_get_logger(error)
        IF (BTEST(cp_print_key_should_output(iteration_info=logger%iter_info,&
             basis_section=print_key, error=new_error),cp_p_file)) THEN

           ! iterate over beads
           DO ib=1,pint_env%p

              ib_str = ""
              WRITE(ib_str,*) ib
              my_middle_name = "pos-" // TRIM(ADJUSTL(ib_str))

              unit_nr=cp_print_key_unit_nr(logger=logger, &
                     basis_section=print_key, print_key_path="", &
                     extension=TRIM(ext), middle_name=TRIM(my_middle_name), &
                     local=.FALSE., file_form=form, is_new_file=new_file, &
                     error=new_error)

              IF (unit_nr>0) THEN
                 CALL force_env_get(force_env=f_env%force_env,&
                      cell=cell, subsys=subsys, error=new_error)
                 CALL cp_subsys_get(subsys,particles=particles,&
                      error=new_error)

                 idim=0
                 DO iat=1,pint_env%ndim/3
                    DO idir=1,3
                       idim=idim+1
                       particles%els(iat)%r(idir)=pint_env%x(ib,idim)
                       particles%els(iat)%v(idir)=pint_env%v(ib,idim)
                    END DO
                 END DO

                IF (outformat == dump_xmol) THEN
                   WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")&
                        " i =",pint_env%iter,&
                       ", E =",pint_env%e_pot_bead(ib)
                END IF

!TODO/FIXME separate printkey for velocities? since both is not
! possible in one file atm it seems
                 CALL write_particle_coordinates(particles%els, &
                      iunit=unit_nr, output_format=outformat, content="POS", &
                      title=title, cell=cell, unit_conv=unit_conv, error=error)
                 CALL cp_print_key_finished_output(unit_nr,logger,&
                      print_key,"", local=.FALSE., error=new_error)
              END IF
           END DO
        END IF
        CALL f_env_rm_defaults(f_env,new_error,ierr,handle)
        CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
     END IF
    END IF
  END SUBROUTINE pint_write_trajectory

  ! ***************************************************************************
  !> \brief  Writes out the energies according to PINT%PRINT%ENER_INFO
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

    INTEGER                                  :: iter, unit_nr
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
      "MOTION%PINT%PRINT%ENER_INFO", error=error)
    logger => cp_error_get_logger(error)
    IF ( BTEST(cp_print_key_should_output(iteration_info=logger%iter_info,&
      basis_section=print_key,error=error),cp_p_file) ) THEN

      unit_nr=cp_print_key_unit_nr(logger,print_key,middle_name="energy",&
        extension=".dat",is_new_file=file_is_new,error=error )

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

        temp = cp_unit_from_cp2k(2.0_dp*pint_env%e_kin_beads/&
               REAL(pint_env%p,dp)/REAL(pint_env%ndim,dp),&
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
  !> \brief Writes the restart file for PINT
  !> \author Fawzi Mohamed
  !> \par History
  !>      11.2003 created [fawzi]
  !>      2009-11-25 major cleanup, logic changed so that one doesn't need to
  !>                 include empty &BEADS and &NOSE sections in the input file
  !>                 any longer, see changes in input_cp2k_motion:create_pint_
  !>                 section as well [lwalewski]
  ! ***************************************************************************
  SUBROUTINE pint_write_restart(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_write_restart', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ibead, idim, inos, status, &
                                                unit_nr
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_vals
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key, tmpsec

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    NULLIFY(print_key)
    print_key => section_vals_get_subs_vals(pint_env%input, &
                 "MOTION%PINT%PRINT%RESTART",error=error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
      basis_section=print_key,error=error),cp_p_file) ) THEN

      ! update iteration count
      CALL section_vals_val_set(pint_env%input, &
           "MOTION%PINT%ITERATION", &
           i_val=pint_env%iter,error=error)

      ! allocate memory for COORDs and VELOCITYs if the BEADS section was not
      ! explicitly given in the input (this is actually done only once since
      ! after section_vals_add_values section becomes explict)
      NULLIFY(tmpsec)
      tmpsec => section_vals_get_subs_vals(pint_env%input,"MOTION%PINT%BEADS",&
                error=error)
      CALL section_vals_get(tmpsec,explicit=explicit,error=error)
      IF ( .NOT. explicit ) THEN
        CALL section_vals_add_values(tmpsec,error)
      END IF

      ! update bead coordinates in the global input structure
      NULLIFY(r_vals)
      ALLOCATE(r_vals(pint_env%p*pint_env%ndim),STAT=status)
      CPAssert(status==0, cp_failure_level, routineP, error, failure)
      i=1
      CALL pint_u2x(pint_env,error=error)
      DO idim=1, pint_env%ndim
        DO ibead=1, pint_env%p
          r_vals(i)=pint_env%x(ibead,idim)
          i = i+1
        END DO
      END DO
      CALL section_vals_val_set(pint_env%input, &
           "MOTION%PINT%BEADS%COORD%_DEFAULT_KEYWORD_", &
           r_vals_ptr=r_vals,error=error)

      ! update bead velocities in the global input structure
      NULLIFY(r_vals)
      ALLOCATE(r_vals(pint_env%p*pint_env%ndim),STAT=status)
      CPAssert(status==0, cp_failure_level, routineP, error, failure)
      i=1
      CALL pint_u2x(pint_env,ux=pint_env%uv,x=pint_env%v,error=error)
      DO idim=1, pint_env%ndim
        DO ibead=1, pint_env%p
          r_vals(i)=pint_env%v(ibead,idim)
          i = i+1
        END DO
      END DO
      CALL section_vals_val_set(pint_env%input, &
           "MOTION%PINT%BEADS%VELOCITY%_DEFAULT_KEYWORD_", &
           r_vals_ptr=r_vals,error=error)

      IF (pint_env%nnos > 0) THEN

        ! allocate memory for COORDs and VELOCITYs if the NOSE section was not
        ! explicitly given in the input (this is actually done only once since
        ! after section_vals_add_values section becomes explict)
        NULLIFY(tmpsec)
        tmpsec => section_vals_get_subs_vals(pint_env%input,"MOTION%PINT%NOSE",&
                  error=error)
        CALL section_vals_get(tmpsec,explicit=explicit,error=error)
        IF ( .NOT. explicit ) THEN
          CALL section_vals_add_values(tmpsec,error)
        END IF

        ! update thermostat coordinates in the global input structure
        NULLIFY(r_vals)
        ALLOCATE(r_vals(pint_env%p*pint_env%ndim*pint_env%nnos),STAT=status)
        CPAssert(status==0, cp_failure_level, routineP, error, failure)
        i=1
        DO idim=1, pint_env%ndim
          DO ibead=1, pint_env%p
            DO inos=1, pint_env%nnos
              r_vals(i)=pint_env%tx(inos,ibead,idim)
              i = i+1
            END DO
          END DO
        END DO
        CALL section_vals_val_set(pint_env%input, &
             "MOTION%PINT%NOSE%COORD%_DEFAULT_KEYWORD_", &
             r_vals_ptr=r_vals,error=error)

        ! update thermostat velocities in the global input structure
        NULLIFY(r_vals)
        ALLOCATE(r_vals(pint_env%p*pint_env%ndim*pint_env%nnos),STAT=status)
        CPAssert(status==0, cp_failure_level, routineP, error, failure)
        i=1
        DO idim=1, pint_env%ndim
          DO ibead=1, pint_env%p
            DO inos=1, pint_env%nnos
              r_vals(i)=pint_env%tv(inos,ibead,idim)
              i = i+1
            END DO
          END DO
        END DO
        CALL section_vals_val_set(pint_env%input, &
             "MOTION%PINT%NOSE%VELOCITY%_DEFAULT_KEYWORD_", &
             r_vals_ptr=r_vals,error=error)
        NULLIFY(r_vals)

      END IF

      ! actually write the restart file
      unit_nr = cp_print_key_unit_nr(logger,print_key,&
                extension=".restart", do_backup=.TRUE., error=error)
      IF (unit_nr>0) THEN
        CALL section_vals_write(pint_env%input,unit_nr,hide_root=.TRUE.,&
             error=error)
        CALL m_flush(unit_nr)
      END IF
      CALL cp_print_key_finished_output(unit_nr,logger,print_key,error=error)

    END IF

    RETURN
  END SUBROUTINE pint_write_restart

  ! ***************************************************************************
  !> \brief  Write step info to the output file.
  !> \date   2009-11-16
  !> \author Lukasz Walewski
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
    unit_nr = cp_logger_get_default_unit_nr(logger)

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

END MODULE pint_io
