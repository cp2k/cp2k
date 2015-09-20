!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  I/O subroutines for helium
!> \author Lukasz Walewski
!> \date   2009-06-08
! *****************************************************************************
MODULE helium_io

  USE cell_types,                      ONLY: get_cell
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE helium_common,                   ONLY: helium_cycle_number,&
                                             helium_is_winding,&
                                             helium_path_length,&
                                             helium_pbc
  USE helium_types,                    ONLY: e_id_interact,&
                                             e_id_kinetic,&
                                             e_id_potential,&
                                             e_id_thermo,&
                                             e_id_total,&
                                             e_id_virial,&
                                             helium_solvent_type
  USE input_constants,                 ONLY: helium_cell_shape_cube,&
                                             helium_cell_shape_octahedron
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_gather,&
                                             mp_sum
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_io'

  PUBLIC :: helium_write_rdf
  PUBLIC :: helium_write_rho
  PUBLIC :: helium_write_line
  PUBLIC :: helium_write_setup
  PUBLIC :: helium_write_energy
  PUBLIC :: helium_write_sdensity
  PUBLIC :: helium_write_wnumber
  PUBLIC :: helium_write_plength
  PUBLIC :: helium_write_coordinates
  PUBLIC :: helium_write_force
  PUBLIC :: helium_write_accepts
  PUBLIC :: helium_write_perm

  CONTAINS


! ***************************************************************************
!> \brief  Write helium parameters to the output unit.
!> \param helium ...
!> \date   2009-06-03
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_write_setup( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_setup', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: my_label, stmp, stmp1, stmp2, &
                                                unit_str
    INTEGER                                  :: i, itmp, j, unit_nr
    INTEGER(kind=int_8)                      :: i8tmp
    LOGICAL                                  :: failure, first
    REAL(kind=dp)                            :: rtmp, v1, v2, v3
    REAL(kind=dp), DIMENSION(3)              :: my_abc
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    CPASSERT(ASSOCIATED(helium))

    NULLIFY(logger)
    logger => cp_get_default_logger()
    my_label = "HELIUM| "

    IF (logger%para_env%ionode) THEN
      unit_nr = cp_logger_get_default_unit_nr(logger)

      WRITE(unit_nr,*)
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Number of helium environments:     ", helium%num_env

      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Number of solvent atoms:           ", helium%atoms
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Number of solvent beads:           ", helium%beads
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Total number of solvent particles: ", helium%atoms*helium%beads

      unit_str = "angstrom^-3"
      rtmp = cp_unit_from_cp2k(helium%density, &
        unit_str)
      WRITE(unit_nr,'(T2,A,F12.6)') TRIM(my_label)//" Density   ["// &
        TRIM(unit_str)//"]:", rtmp

      unit_str = "angstrom"
      rtmp = cp_unit_from_cp2k(helium%cell_size, &
        unit_str)
      WRITE(unit_nr,'(T2,A,F12.6)') TRIM(my_label)//" Cell size ["// &
        TRIM(unit_str)//"]:   ", rtmp

      IF ( helium%periodic ) THEN
        IF ( helium%cell_shape .EQ. helium_cell_shape_cube ) THEN
          CALL helium_write_line("PBC cell shape: CUBE.")
        ELSE IF ( helium%cell_shape .EQ. helium_cell_shape_octahedron ) THEN
          CALL helium_write_line("PBC cell shape: TRUNCATED OCTAHEDRON.")
        ELSE
          CALL helium_write_line("*** Warning: unknown cell shape.")
        END IF
      ELSE
        CALL helium_write_line("PBC turned off.")
      END IF

      ! first step gets incremented during first iteration
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " First MC step                      :", helium%first_step + 1
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Last MC step                       :", helium%last_step
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Total number of MC steps           :", helium%num_steps
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Number of outer MC trials per step :", helium%iter_rot
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Number of inner MC trials per step :", helium%iter_norot
      i8tmp = helium%iter_rot
      i8tmp = i8tmp * helium%iter_norot
      stmp = ""
      WRITE(stmp, *) i8tmp
      WRITE(unit_nr,'(T2,A)') TRIM(my_label)//&
        " Total number of MC trials per step : " // TRIM(ADJUSTL(stmp))
      i8tmp = helium%num_steps
      i8tmp = i8tmp * helium%iter_rot
      i8tmp = i8tmp * helium%iter_norot
      stmp = ""
      WRITE(stmp, *) i8tmp
      WRITE(unit_nr,'(T2,A)') TRIM(my_label) //&
        " Total number of MC trials          : " // TRIM(ADJUSTL(stmp))

      ! permutation cycle length sampling
      stmp = ""
      CALL helium_write_line(stmp)
      WRITE(stmp,*) helium%maxcycle
      stmp2 = ""
      WRITE(stmp2,*) "Using maximum permutation cycle length: " //&
        TRIM(ADJUSTL(stmp))
      CALL helium_write_line(stmp2)
      stmp = ""
      stmp1 = ""
      WRITE(stmp1,*) helium%m_ratio
      stmp2 = ""
      WRITE(stmp2,*) helium%m_value
      WRITE(stmp,*) "Using ratio " // TRIM(ADJUSTL(stmp1)) // " for M = " // TRIM(ADJUSTL(stmp2))
      CALL helium_write_line(stmp)
      stmp = ""
      CALL helium_write_line(stmp)

      IF (helium%solute_present) THEN
        WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
          " Number of solute atoms:            ", helium%solute_atoms
        WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
          " Number of solute beads:            ", helium%solute_beads
        WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
          " Total number of solute particles:  ", helium%solute_atoms*&
          helium%solute_beads
        DO i = 1, 3
          WRITE(unit_nr,'(T2,A,1X,I0,A,I0,A)', ADVANCE='NO') TRIM(my_label)//&
            " Solute atom type: ", i,&
            "(",helium%solute_number(i),")"
          first = .TRUE.
          DO j = 1, helium%solute_number(i)
            IF (first) THEN
              WRITE(unit_nr,'(A)', ADVANCE='NO')", indices: "
              first = .FALSE.
            ELSE
              WRITE(unit_nr,'(A)', ADVANCE='NO') ", "
            END IF
            WRITE(unit_nr,'(I0)', ADVANCE='NO') helium%solute_index(i,j)
          END DO
          WRITE(unit_nr,'(1X)')
        END DO
        CALL get_cell(helium%solute_cell, abc=my_abc)
        unit_str = "angstrom"
        v1 = cp_unit_from_cp2k(my_abc(1), unit_str)
        v2 = cp_unit_from_cp2k(my_abc(2), unit_str)
        v3 = cp_unit_from_cp2k(my_abc(3), unit_str)
        WRITE(unit_nr,'(T2,A,F12.6,1X,F12.6,1X,F12.6)') &
          TRIM(my_label)//" Solute cell size ["// &
          TRIM(unit_str)//"]:   ", v1, v2, v3
      ELSE
        WRITE(unit_nr,'(T2,A)') TRIM(my_label)//" Solute is NOT present"
      END IF
    END IF

    ! radial distribution function related settings
    rtmp = cp_unit_from_cp2k(helium%rdf_delr, "angstrom")
    WRITE(stmp, '(1X,F12.6)') rtmp
    CALL helium_write_line("RDF| delr [angstrom]: "//TRIM(stmp))
    rtmp = cp_unit_from_cp2k(helium%rdf_maxr, "angstrom")
    WRITE(stmp, '(1X,F12.6)') rtmp
    CALL helium_write_line("RDF| maxr [angstrom]: "//TRIM(stmp))
    itmp = helium%rdf_nbin
    WRITE(stmp, '(I6)') itmp
    CALL helium_write_line("RDF| nbin           : "//TRIM(stmp))

    CALL helium_write_line("")

    ! density related
    IF (helium%rho_present) THEN
      CALL helium_write_line("RHO| Calculating densities on the grid")
      itmp = helium%rho_nbin
      stmp = ""
      WRITE(stmp, '(I6)') itmp
      CALL helium_write_line("RHO| "//TRIM(ADJUSTL(stmp))//"x"//&
            TRIM(ADJUSTL(stmp))//"x"//TRIM(ADJUSTL(stmp)))
      CALL helium_write_line("")
    END IF

    RETURN
  END SUBROUTINE helium_write_setup

! ***************************************************************************
!> \brief  Writes out a line of text to the default output unit.
!> \param line ...
!> \date   2009-07-10
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_write_line(line)

    CHARACTER(len=*), INTENT(IN)             :: line

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_line', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: my_label
    INTEGER                                  :: unit_nr
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_get_default_logger()
    my_label = "HELIUM|"

    IF (logger%para_env%ionode) THEN
      unit_nr = cp_logger_get_default_unit_nr(logger)
      WRITE(unit_nr,'(T2,A)') TRIM(my_label)//" "//TRIM(line)
    END IF

    RETURN
  END SUBROUTINE helium_write_line

! ***************************************************************************
!> \brief  Writes out helium energies according to HELIUM%PRINT%ENERGY
!> \param helium ...
!> \date   2009-06-08
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_write_energy( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iteration, m, unit_nr
    LOGICAL                                  :: failure, file_is_new
    REAL(kind=dp)                            :: naccptd
    REAL(kind=dp), DIMENSION(:), POINTER     :: my_energy
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPASSERT(ASSOCIATED(helium))
    CPASSERT(ASSOCIATED(helium%input))

    my_energy => helium%energy_avrg

    NULLIFY(print_key,logger)
    print_key => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM%PRINT%ENERGY")
    logger => cp_get_default_logger()
    iteration = logger%iter_info%iteration(2)
    IF ( BTEST(cp_print_key_should_output(iteration_info=logger%iter_info,&
      basis_section=print_key),cp_p_file) ) THEN
      unit_nr=cp_print_key_unit_nr(logger,print_key, &
        middle_name="helium-energy",extension=".dat",&
        is_new_file=file_is_new)
      ! cp_print_key_unit_nr returns -1 on nodes other than logger%para_env%ionode

      naccptd = 0.0_dp
      DO m = 1, helium%maxcycle
        naccptd = naccptd + helium%num_accepted(helium%bisctlog2+2,m)
      END DO
      CALL mp_sum(naccptd,logger%para_env%group)
      naccptd = naccptd / REAL(logger%para_env%num_pe,dp)

      IF (unit_nr>0) THEN

        IF ( file_is_new ) THEN
          WRITE(unit_nr,'(A9,1X,A12,6(1X,A20))')&
                           "#    Step",&
                        "     Naccptd",&
                "               E_pot",&
                "               E_kin",&
                "            E_thermo",&
                "            E_virial",&
                "             E_inter",&
                "               E_tot"
        END IF

        WRITE (unit_nr,"(I9,1X,F12.1,6(1X,F20.9))") &
          iteration, &
          naccptd,&
          my_energy(e_id_potential), &
          my_energy(e_id_kinetic), &
          my_energy(e_id_thermo), &
          my_energy(e_id_virial), &
          my_energy(e_id_interact), &
          my_energy(e_id_total)
        CALL m_flush(unit_nr)
      END IF
      CALL cp_print_key_finished_output(unit_nr,logger,print_key)
    END IF

    CALL timestop(handle)
    RETURN
  END SUBROUTINE helium_write_energy

! ***************************************************************************
!> \brief  Writes out helium energies according to HELIUM%PRINT%SDENSITY
!> \param helium ...
!> \date   2010-06-15
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_write_sdensity( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_sdensity', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iteration, unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPASSERT(ASSOCIATED(helium))
    CPASSERT(ASSOCIATED(helium%input))

    NULLIFY(print_key,logger)
    print_key => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM%PRINT%SDENSITY")
    logger => cp_get_default_logger()
    iteration = logger%iter_info%iteration(2)
    IF ( BTEST(cp_print_key_should_output(iteration_info=logger%iter_info,&
      basis_section=print_key),cp_p_file) ) THEN
      unit_nr=cp_print_key_unit_nr(logger,print_key, &
        middle_name="helium-sdensity",extension=".dat")
      IF (unit_nr>0) THEN
        WRITE (unit_nr,"(F20.9)") helium%sdensity_avrg
        CALL m_flush(unit_nr)
      END IF
      CALL cp_print_key_finished_output(unit_nr,logger,print_key)
    END IF

    CALL timestop(handle)
    RETURN
  END SUBROUTINE helium_write_sdensity

! ***************************************************************************
!> \brief  Writes out helium winding number according to HELIUM%PRINT%WNUMBER
!> \param helium ...
!> \date   2009-10-19
!> \par    History
!>         2010-06-15 output W for each He environment/processor [lwalewski]
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_write_wnumber( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_wnumber', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPASSERT(ASSOCIATED(helium))
    CPASSERT(ASSOCIATED(helium%input))

    NULLIFY(print_key,logger)
    print_key => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM%PRINT%WNUMBER")
    logger => cp_get_default_logger()
    IF ( BTEST(cp_print_key_should_output(iteration_info=logger%iter_info,&
      basis_section=print_key),cp_p_file) ) THEN
      unit_nr=cp_print_key_unit_nr(logger,print_key, &
        middle_name="helium-wnumber",extension=".dat")

      ! gather winding number from all processors to logger%para_env%source
      helium%rtmp_3_np_1d(:) = 0
      CALL mp_gather( helium%wnumber_avrg, &
                      helium%rtmp_3_np_1d, &
                      logger%para_env%source, logger%para_env%group )

      IF (unit_nr>0) THEN

        DO i = 1, 3 * helium%num_env
          WRITE(unit_nr,'(F20.9)',ADVANCE='NO') helium%rtmp_3_np_1d(i)
          IF ( i .LT. 3 * helium%num_env ) THEN
            WRITE(unit_nr,'(1X)',ADVANCE='NO')
          END IF
        END DO
        WRITE(unit_nr, '(A)') ""

        CALL m_flush(unit_nr)
      END IF
      CALL cp_print_key_finished_output(unit_nr,logger,print_key)
    END IF

    CALL timestop(handle)
    RETURN
  END SUBROUTINE helium_write_wnumber

! ***************************************************************************
!> \brief  Writes out acceptance counts according to HELIUM%PRINT%ACCEPTS
!> \param helium ...
!> \date   2010-05-27
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_write_accepts( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_accepts', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iteration, j, &
                                                unit_nr
    LOGICAL                                  :: failure, file_is_new
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPASSERT(ASSOCIATED(helium))
    CPASSERT(ASSOCIATED(helium%input))

    NULLIFY(print_key,logger)
    print_key => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM%PRINT%ACCEPTS")
    logger => cp_get_default_logger()
    iteration = logger%iter_info%iteration(2)
    IF ( BTEST(cp_print_key_should_output(iteration_info=logger%iter_info,&
      basis_section=print_key),cp_p_file) ) THEN
      unit_nr=cp_print_key_unit_nr(logger,print_key, &
        middle_name="helium-accepts",extension=".dat",&
        is_new_file=file_is_new)
      IF (unit_nr>0) THEN

        IF ( file_is_new ) THEN
          WRITE(unit_nr,'(A8,1X,A15,1X,A20)',ADVANCE='NO')&
                            "# Length",&
                     "         Trials",&
                "            Selected"
          DO j = 1, helium%bisctlog2
            WRITE(unit_nr,'(A17,1X,I3)',ADVANCE='NO') "            Level", j
          END DO
          WRITE(unit_nr, '(A)') ""
        END IF

        DO i = 1, helium%maxcycle
          WRITE(unit_nr, '(I3)',ADVANCE='NO') i
          DO j = 1, helium%bisctlog2 + 2
            WRITE(unit_nr,'(1X,F20.2)',ADVANCE='NO') helium%num_accepted(j,i)
          END DO
          WRITE(unit_nr, '(A)') ""
        END DO
        WRITE(unit_nr, '(A)') "&"

        CALL m_flush(unit_nr)
      END IF
      CALL cp_print_key_finished_output(unit_nr,logger,print_key)
    END IF

    CALL timestop(handle)
    RETURN
  END SUBROUTINE helium_write_accepts

! ***************************************************************************
!> \brief  Writes out permutation state according to HELIUM%PRINT%PERM
!> \param helium ...
!> \date   2010-06-07
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_write_perm( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_perm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iteration, j, &
                                                offset, unit_nr
    LOGICAL                                  :: failure, file_is_new
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPASSERT(ASSOCIATED(helium))
    CPASSERT(ASSOCIATED(helium%input))

    NULLIFY(print_key,logger)
    print_key => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM%PRINT%PERM")
    logger => cp_get_default_logger()
    iteration = logger%iter_info%iteration(2)
    IF ( BTEST(cp_print_key_should_output(iteration_info=logger%iter_info,&
      basis_section=print_key),cp_p_file) ) THEN
      unit_nr=cp_print_key_unit_nr(logger,print_key, &
        middle_name="helium-perm",extension=".dat",&
        is_new_file=file_is_new)

      ! gather permutation state from all processors to logger%para_env%source
      helium%itmp_atoms_np_1d(:) = 0
      CALL mp_gather( helium%permutation, &
                      helium%itmp_atoms_np_1d, &
                      logger%para_env%source, logger%para_env%group )

      IF (unit_nr>0) THEN

        DO i = 1, helium%atoms
          DO j = 1, helium%num_env
            offset = (j-1) * helium%atoms
            WRITE(unit_nr,'(I6)',ADVANCE='NO') helium%itmp_atoms_np_1d(offset+i)
            IF ( j .LT. helium%num_env ) THEN
              WRITE(unit_nr,'(1X)',ADVANCE='NO')
            END IF
          END DO
          WRITE(unit_nr, '(A)') ""
        END DO
        WRITE(unit_nr, '(A)') "&"

        CALL m_flush(unit_nr)
      END IF
      CALL cp_print_key_finished_output(unit_nr,logger,print_key)
    END IF

    CALL timestop(handle)
    RETURN
  END SUBROUTINE helium_write_perm

! ***************************************************************************
!> \brief  Writes helium configuration according to HELIUM%PRINT%COORDINATES
!> \param helium ...
!> \date   2009-07-16
!> \par    History
!>         2010-02-15 output from all processors added [lwalewski]
!> \author Lukasz Walewski
!> \note   particle_types->write_particle_coordinates is of no use here,
!>         since it does not support atom connectivity information needed
!>         for helium paths
! *****************************************************************************
  SUBROUTINE helium_write_coordinates( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_coordinates', &
      routineP = moduleN//':'//routineN

    CHARACTER(3)                             :: resName
    CHARACTER(len=default_string_length)     :: fmt_string, my_middle_name, &
                                                stmp
    INTEGER                                  :: handle, ia, ib, ib1, ib2, ic, &
                                                icycle, irank, msglen, &
                                                offset, tmp1, tmp2, unit_nr
    INTEGER, DIMENSION(:), POINTER           :: my_perm
    LOGICAL                                  :: are_connected, failure, &
                                                is_winding, ltmp, &
                                                should_output
    REAL(kind=dp)                            :: xtmp, ytmp, ztmp
    REAL(kind=dp), DIMENSION(3)              :: r0, r1, r2
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPASSERT(ASSOCIATED(helium))

    NULLIFY(logger,print_key)
    logger => cp_get_default_logger()

    ! decide whether to write anything or not
    print_key => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM%PRINT%COORDINATES")
    should_output = BTEST(cp_print_key_should_output( &
      iteration_info=logger%iter_info,&
      basis_section=print_key),cp_p_file)
    IF ( .NOT. should_output ) THEN
      CALL timestop(handle)
      RETURN
    END IF

    ! prepare the coordinates for output (print one image of the periodic
    ! space that is centered around r0)
    r0(:) = helium%origin(:)
    DO ia = 1, helium%atoms
      DO ib = 1, helium%beads
        r1(:) = helium%pos(:,ia,ib) - r0(:)
        r2(:) = helium%pos(:,ia,ib) - r0(:)
        CALL helium_pbc( helium, r2 )
        ltmp = .FALSE.
        DO ic = 1, 3
          IF ( ABS(r1(ic)-r2(ic)) .GT. 100.0_dp*EPSILON(0.0_dp) ) THEN
            ltmp = .TRUE.
            CYCLE
          END IF
        END DO
        IF ( ltmp ) THEN
          helium%work(:,ia,ib) = r0(:) + r2(:)
        ELSE
          helium%work(:,ia,ib) = helium%pos(:,ia,ib)
        END IF
      END DO
    END DO

    ! gather positions from all processors to logger%para_env%source
    helium%rtmp_3_atoms_beads_1d(:) = PACK( helium%work, .TRUE. )
    CALL mp_gather( helium%rtmp_3_atoms_beads_1d, &
                    helium%rtmp_3_atoms_beads_np_1d, &
                    logger%para_env%source, logger%para_env%group )

    ! gather permutation state from all processors to logger%para_env%source
    CALL mp_gather( helium%permutation, &
                    helium%itmp_atoms_np_1d, &
                    logger%para_env%source, logger%para_env%group )

    ! set logical mask for unpacking coordinates gathered from other ranks
    helium%ltmp_3_atoms_beads_3d(:,:,:) = .TRUE.

    ! I/O only on the ionode
    IF (logger%para_env%ionode) THEN

      ! iterate over processors/helium environments
      DO irank = 1, helium%num_env

        ! generate one file per processor
        stmp = ""
        WRITE(stmp,*) irank
        my_middle_name = "helium-pos-" // TRIM(ADJUSTL(stmp))
        unit_nr=cp_print_key_unit_nr( logger, print_key, &
          middle_name=TRIM(my_middle_name), extension=".pdb")

        ! write out the unit cell parameters
        fmt_string = "(A6,3F9.3,3F7.2,1X,A11,1X,I3)"
        xtmp = helium%cell_size
        xtmp = cp_unit_from_cp2k(xtmp, "angstrom")
        IF ( helium%cell_shape .EQ. helium_cell_shape_octahedron ) THEN
          stmp = "O          "
        ELSE
          stmp = "C          "
        END IF
        WRITE(unit_nr,fmt_string) "CRYST1", &
          xtmp, xtmp, xtmp, &
          90.0_dp, 90.0_dp, 90.0_dp, &
          stmp, helium%beads

        ! unpack coordinates
        msglen = SIZE(helium%rtmp_3_atoms_beads_1d)
        offset = (irank-1) * msglen
        helium%work(:,:,:) = &
          UNPACK(helium%rtmp_3_atoms_beads_np_1d(offset+1:offset+msglen), &
          MASK=helium%ltmp_3_atoms_beads_3d, FIELD=0.0_dp )

        ! unpack permutation state (actually point to the right section only)
        msglen = SIZE(helium%permutation)
        offset = (irank-1) * msglen
        my_perm => helium%itmp_atoms_np_1d(offset+1:offset+msglen)

        ! write out coordinates
        fmt_string = &
          "(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)"
        DO ia = 1, helium%atoms
          icycle = helium_cycle_number(helium, ia, my_perm)
          is_winding = helium_is_winding(helium,ia,helium%work,my_perm)
          IF ( is_winding ) THEN
            resName = "SPR"
          ELSE
            resName = "NRM"
          END IF
          DO ib = 1, helium%beads
            xtmp = helium%work(1,ia,ib)
            xtmp = cp_unit_from_cp2k(xtmp, "angstrom")
            ytmp = helium%work(2,ia,ib)
            ytmp = cp_unit_from_cp2k(ytmp, "angstrom")
            ztmp = helium%work(3,ia,ib)
            ztmp = cp_unit_from_cp2k(ztmp, "angstrom")
            WRITE(unit_nr,fmt_string) "ATOM  ", &
              (ia-1)*helium%beads+ib, &
              " He ", " ", resName, "X", &
              icycle, &
              " ", &
              xtmp, ytmp, ztmp, &
              1.0_dp, 0.0_dp, "HE", "  "
          END DO
        END DO

        ! write out the bead connectivity information
        DO ia = 1, helium%atoms

          ! write connectivity records for this atom only if the path
          ! it belongs to is longer than 1.
          IF ( helium_path_length(helium, ia, my_perm) .LE. 1 ) THEN
            CYCLE
          END IF

          DO ib = 1, helium%beads-1
            ! check wheather the consecutive beads belong to the same box
            r1(:) = helium%work(:,ia,ib) - helium%work(:,ia,ib+1)
            r2(:) = r1(:)
            CALL helium_pbc( helium, r2 )
            are_connected = .TRUE.
            DO ic = 1, 3
              IF ( ABS(r1(ic)-r2(ic)) .GT. 100.0_dp*EPSILON(0.0_dp) ) THEN
                ! if the distance betw ib and ib+1 changes upon applying
                ! PBC do not connect them
                are_connected = .FALSE.
                CYCLE
              END IF
            END DO
            IF ( are_connected ) THEN
              tmp1 = (ia-1)*helium%beads+ib
              tmp2 = (ia-1)*helium%beads+ib+1
              ! smaller value has to go first
              IF (tmp1 .LT. tmp2) THEN
                ib1 = tmp1
                ib2 = tmp2
              ELSE
                ib1 = tmp2
                ib2 = tmp1
              END IF
              WRITE(unit_nr,'(A6,2I5)') "CONECT", ib1, ib2
            END IF
          END DO

          ! last bead of atom <ia> connects to the first bead
          ! of the next atom in the permutation cycle
          r1(:) = helium%work(:,ia,helium%beads) - helium%work(:,my_perm(ia),1)
          r2(:) = r1(:)
          CALL helium_pbc( helium, r2 )
          are_connected = .TRUE.
          DO ic = 1, 3
            IF ( ABS(r1(ic)-r2(ic)) .GT. 100.0_dp*EPSILON(0.0_dp) ) THEN
              ! if the distance betw ib and ib+1 changes upon applying
              ! PBC do not connect them
              are_connected = .FALSE.
              CYCLE
            END IF
          END DO
          IF ( are_connected ) THEN
            tmp1 = ia*helium%beads
            tmp2 = (my_perm(ia)-1)*helium%beads+1
            IF (tmp1 .LT. tmp2) THEN
              ib1 = tmp1
              ib2 = tmp2
            ELSE
              ib1 = tmp2
              ib2 = tmp1
            END IF
            WRITE(unit_nr,'(A6,2I5)') "CONECT", ib1, ib2
          END IF
        END DO
        WRITE(unit_nr,'(A)') "END"

        CALL m_flush(unit_nr)
        CALL cp_print_key_finished_output(unit_nr,logger,print_key)

      END DO

    END IF

    CALL timestop(handle)
    RETURN
  END SUBROUTINE helium_write_coordinates

! ***************************************************************************
!> \brief  Write helium RDF according to HELIUM%PRINT%RDF
!> \param helium ...
!> \date   2009-07-23
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_write_rdf( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_rdf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, unit_nr
    LOGICAL                                  :: failure, is_new, should_output
    REAL(kind=dp)                            :: rtmp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPASSERT(ASSOCIATED(helium))
    CPASSERT(ASSOCIATED(helium%input))

    NULLIFY(logger,print_key)
    logger => cp_get_default_logger()
    print_key => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM%PRINT%RDF")
    should_output = BTEST(cp_print_key_should_output( &
      iteration_info=logger%iter_info,&
      basis_section=print_key),cp_p_file)
    IF (should_output) THEN
      IF (logger%para_env%ionode) THEN
        unit_nr=cp_print_key_unit_nr(logger,print_key, &
          middle_name="helium-hst",extension=".dat",&
          is_new_file=is_new)
        IF (.NOT. is_new) THEN
          WRITE(unit_nr,'(A1)') "&"
        END IF
        DO i = 1, helium%rdf_nbin
          rtmp = ( REAL(i) - 0.5_dp ) * helium%rdf_delr
          rtmp = cp_unit_from_cp2k(rtmp, "angstrom")
          WRITE(unit_nr,'(2F20.10)') rtmp,&
            helium%rdf_avrg(i)
        END DO
        CALL m_flush(unit_nr)
        CALL cp_print_key_finished_output(unit_nr,logger,print_key)
      END IF
    END IF

    CALL timestop(handle)
    RETURN
  END SUBROUTINE helium_write_rdf


! ***************************************************************************
!> \brief  Write helium densities according to HELIUM%PRINT%RHO
!> \param helium ...
!> \date   2011-06-21
!> \par    History
!>         2011-11-11 output of any number of density estimators to different
!>                    files according to helium%rho_num
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_write_rho( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_rho', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: stmp
    INTEGER                                  :: handle, i, itmp, iweight, ix, &
                                                iy, iz, N, nsteps, unit_nr
    LOGICAL                                  :: failure, should_output
    REAL(kind=dp)                            :: inv_norm, rtmp, rx, ry, rz
    REAL(kind=dp), DIMENSION(3)              :: r0
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: message
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPASSERT(ASSOCIATED(helium))
    CPASSERT(ASSOCIATED(helium%input))

    NULLIFY(logger,print_key)
    logger => cp_get_default_logger()
    print_key => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM%PRINT%RHO")
    should_output = BTEST(cp_print_key_should_output( &
      iteration_info=logger%iter_info,&
      basis_section=print_key),cp_p_file)
    IF (should_output) THEN

      ! work on the temporary array so that accumulated data remains intact
      helium%rho_inst(:,:,:,:) = helium%rho_avrg(:,:,:,:)

      ! average over helium environments
      DO i = 1, helium%rho_num
        NULLIFY(message)
        message => helium%rho_inst(i,1::1,1::1,1::1)
        CALL mp_sum(message,logger%para_env%group)
      END DO
      itmp = logger%para_env%num_pe
      CPASSERT(itmp>0)
      inv_norm = 1.0_dp / REAL(itmp,dp)
      helium%rho_inst(:,:,:,:) = helium%rho_inst(:,:,:,:) * inv_norm

      ! average over steps performed so far in this run
      nsteps = helium%current_step-helium%first_step
      CPASSERT(nsteps>0)
      inv_norm = 1.0_dp / REAL(nsteps,dp)
      helium%rho_inst(:,:,:,:) = helium%rho_inst(:,:,:,:) * inv_norm

      iweight = helium%rho_iweight
      IF ( helium%rho_restart ) THEN
        ! average over the old and the current density (observe the weights!)
        helium%rho_inst(:,:,:,:) = nsteps * helium%rho_inst(:,:,:,:) + &
                                   iweight * helium%rho_rstr(:,:,:,:)
        helium%rho_inst(:,:,:,:) = helium%rho_inst(:,:,:,:) / ( nsteps + iweight )
      END IF

      IF (logger%para_env%ionode) THEN

        N = helium%rho_nbin
        rtmp = cp_unit_from_cp2k(helium%rho_delr, "bohr")
        r0(:) = helium%origin(:)
        rx = r0(1) - helium%rho_maxr / 2.0_dp
        ry = r0(2) - helium%rho_maxr / 2.0_dp
        rz = r0(3) - helium%rho_maxr / 2.0_dp

        DO i = 1, helium%rho_num
        stmp = ""
        WRITE(stmp,*) i
        unit_nr=cp_print_key_unit_nr(logger,print_key, &
          middle_name="helium-rho-"//TRIM(ADJUSTL(stmp)), &
          extension=".cube", &
          file_position="REWIND", do_backup=.FALSE.)

        WRITE(unit_nr,'(A)') "Helium density " // TRIM(ADJUSTL(stmp))
        WRITE(unit_nr,'(A)') "CP2K"
        WRITE(unit_nr,'(I5,3F12.8)') 1, rx, ry, rz
        WRITE(unit_nr,'(I5,3F12.8)') N, rtmp, 0.0_dp, 0.0_dp
        WRITE(unit_nr,'(I5,3F12.8)') N, 0.0_dp, rtmp, 0.0_dp
        WRITE(unit_nr,'(I5,3F12.8)') N, 0.0_dp, 0.0_dp, rtmp
        WRITE(unit_nr,'(I5,4F12.8)') 1, 0.0_dp, r0(1), r0(2), r0(3)
        DO ix = 1, N
          DO iy = 1, N
            DO iz = 1, N
              WRITE(unit_nr,'(E13.5,1X)', ADVANCE='NO') helium%rho_inst(i,ix,iy,iz)
              IF ( MOD(iz,6) .EQ. 5 ) THEN
                WRITE(unit_nr,*)
              END IF
            END DO
            WRITE(unit_nr,*)
          END DO
        END DO

        CALL m_flush(unit_nr)
        CALL cp_print_key_finished_output(unit_nr,logger,print_key)

        END DO

      END IF
    END IF

    CALL timestop(handle)
  END SUBROUTINE helium_write_rho


! ***************************************************************************
!> \brief  Write helium permutation length according to HELIUM%PRINT%PLENGTH
!> \param helium ...
!> \date   2010-06-07
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_write_plength( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_plength', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, unit_nr
    LOGICAL                                  :: failure, is_new, should_output
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPASSERT(ASSOCIATED(helium))
    CPASSERT(ASSOCIATED(helium%input))

    NULLIFY(logger,print_key)
    logger => cp_get_default_logger()
    print_key => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM%PRINT%PLENGTH")
    should_output = BTEST(cp_print_key_should_output( &
      iteration_info=logger%iter_info,&
      basis_section=print_key),cp_p_file)
    IF (should_output) THEN
      IF (logger%para_env%ionode) THEN
        unit_nr=cp_print_key_unit_nr(logger,print_key, &
          middle_name="helium-plength",extension=".dat",&
          is_new_file=is_new)

        DO i = 1, helium%atoms
          WRITE(unit_nr,'(F20.10)',ADVANCE='NO') helium%plength_avrg(i)
            IF ( i .LT. helium%atoms ) THEN
              WRITE(unit_nr,'(1X)',ADVANCE='NO')
            END IF
        END DO
        WRITE(unit_nr, '(A)') ""

        CALL m_flush(unit_nr)
        CALL cp_print_key_finished_output(unit_nr,logger,print_key)
      END IF
    END IF

    CALL timestop(handle)
    RETURN
  END SUBROUTINE helium_write_plength

! ***************************************************************************
!> \brief  Write helium force according to HELIUM%PRINT%FORCE
!> \param helium ...
!> \date   2010-01-27
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_write_force( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_force', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: msgstr
    INTEGER                                  :: handle, ia, ib, ic, idim, &
                                                unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPASSERT(ASSOCIATED(helium))

    NULLIFY(logger)
    logger => cp_get_default_logger()

    ! decide whether to write anything or not
    NULLIFY(print_key)
    print_key => section_vals_get_subs_vals( helium%input, &
      "MOTION%PINT%HELIUM%PRINT%FORCES")
    IF ( .NOT. BTEST(cp_print_key_should_output(logger%iter_info, &
        basis_section=print_key),cp_p_file) ) THEN
      CALL timestop(handle)
      RETURN
    END IF

    ! check if there is anything to be printed out
    IF ( .NOT. helium%solute_present ) THEN
      msgstr = "Warning: force printout requested but there is no solute!"
      CALL helium_write_line( msgstr)
      CALL timestop(handle)
      RETURN
    END IF

    ! I/O only on the ionode
    IF (logger%para_env%ionode) THEN

      unit_nr=cp_print_key_unit_nr(logger, print_key, &
              middle_name="helium-force",extension=".dat")

      ! print all force components in one line
      DO ib = 1, helium%solute_beads
        idim = 0
        DO ia = 1, helium%solute_atoms
          DO ic = 1, 3
            idim = idim + 1
            WRITE(unit_nr,'(F20.10)',ADVANCE='NO') helium%force_avrg(ib,idim)
          END DO
        END DO
      END DO
      WRITE(unit_nr,*)

      ! finalize the printout
      CALL m_flush(unit_nr)
      CALL cp_print_key_finished_output(unit_nr,logger,print_key)

    END IF

    CALL timestop(handle)
    RETURN
  END SUBROUTINE helium_write_force

! ***************************************************************************
!> \brief  Write instantaneous helium forces
!> \param helium ...
!> \date   2010-01-29
!> \author Lukasz Walewski
!> \note Collects instantaneous helium forces from all processors on
!>         logger%para_env%source and writes them to files - one file per processor.
!>         This subroutine does message passing, frequent calls can slow down your
!>         code significantly.
! *****************************************************************************
  SUBROUTINE helium_write_force_inst( helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_force_inst', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: my_middle_name, stmp
    INTEGER                                  :: handle, ia, ib, ic, idim, &
                                                irank, offset, unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPASSERT(ASSOCIATED(helium))

    NULLIFY(logger)
    logger => cp_get_default_logger()

    ! decide whether to write anything or not
    NULLIFY(print_key)
    print_key => section_vals_get_subs_vals( helium%input, &
      "MOTION%PINT%HELIUM%PRINT%FORCES_INST")
    IF ( .NOT. BTEST(cp_print_key_should_output(logger%iter_info, &
        basis_section=print_key),cp_p_file) ) THEN
      CALL timestop(handle)
      RETURN
    END IF

    ! check if there is anything to be printed out
    IF ( .NOT. helium%solute_present ) THEN
      stmp = "Warning: force printout requested but there is no solute!"
      CALL helium_write_line( stmp)
      CALL timestop(handle)
      RETURN
    END IF

    ! fill the tmp buffer with instantaneous helium forces at each proc
    helium%rtmp_p_ndim_1d(:) = PACK( helium%force_inst, .TRUE. )

    ! pass the message from all processors to logger%para_env%source
    helium%rtmp_p_ndim_np_1d(:) = 0.0_dp
    CALL mp_gather( helium%rtmp_p_ndim_1d, helium%rtmp_p_ndim_np_1d, &
         logger%para_env%source, logger%para_env%group )

    ! I/O only on the ionode
    IF (logger%para_env%ionode) THEN

      ! iterate over processors/helium environments
      DO irank = 1, helium%num_env

        ! generate one file per processor
        stmp = ""
        WRITE(stmp,*) irank
        my_middle_name = "helium-force-inst-" // TRIM(ADJUSTL(stmp))
        unit_nr=cp_print_key_unit_nr( logger, print_key, &
          middle_name=TRIM(my_middle_name), extension=".dat")

        ! unpack and actually print the forces - all components in one line
        offset = (irank-1) * SIZE(helium%rtmp_p_ndim_1d)
        idim = 0
        DO ib = 1, helium%solute_beads
          DO ia = 1, helium%solute_atoms
            DO ic = 1, 3
              idim = idim + 1
              WRITE(unit_nr,'(F20.10)',ADVANCE='NO') helium%rtmp_p_ndim_np_1d(offset+idim)
            END DO
          END DO
        END DO
        WRITE(unit_nr,*)

        ! finalize the printout
        CALL m_flush(unit_nr)
        CALL cp_print_key_finished_output(unit_nr,logger,print_key)

      END DO

    END IF

    CALL timestop(handle)
    RETURN
  END SUBROUTINE helium_write_force_inst

END MODULE helium_io
