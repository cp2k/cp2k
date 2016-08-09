!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  I/O subroutines for helium
!> \author Lukasz Walewski
!> \date   2009-06-08
! **************************************************************************************************
MODULE helium_io

   USE cell_types,                      ONLY: get_cell
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                              parser_get_object
   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release
   USE cp_units,                        ONLY: cp_unit_from_cp2k,&
                                              cp_unit_to_cp2k
   USE helium_common,                   ONLY: helium_cycle_number,&
                                              helium_cycle_of,&
                                              helium_is_winding,&
                                              helium_path_length,&
                                              helium_pbc
   USE helium_interactions,             ONLY: helium_total_inter_action,&
                                              helium_total_link_action,&
                                              helium_total_pair_action
   USE helium_types,                    ONLY: &
        e_id_interact, e_id_kinetic, e_id_potential, e_id_thermo, e_id_total, e_id_virial, &
        helium_solvent_p_type, helium_solvent_type, rho_num
   USE input_constants,                 ONLY: helium_cell_shape_cube,&
                                              helium_cell_shape_octahedron,&
                                              helium_sampling_ceperley,&
                                              helium_solute_intpot_mwater,&
                                              helium_solute_intpot_none,&
                                              perm_cycle,&
                                              perm_plain
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              int_8
   USE machine,                         ONLY: m_flush
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_gather,&
                                              mp_sum
   USE physcon,                         ONLY: angstrom
   USE pint_types,                      ONLY: pint_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_io'

   CHARACTER(len=*), DIMENSION(6), PARAMETER, PRIVATE :: m_dist_name = (/ &
                                                         "SINGLEV    ", &
                                                         "UNIFORM    ", &
                                                         "LINEAR     ", &
                                                         "QUADRATIC  ", &
                                                         "EXPONENTIAL", &
                                                         "GAUSSIAN   "/)

   PUBLIC :: helium_read_xyz
   PUBLIC :: helium_print_rdf
   PUBLIC :: helium_print_rho
   PUBLIC :: helium_write_line
   PUBLIC :: helium_write_setup
   PUBLIC :: helium_print_energy
   PUBLIC :: helium_print_vector
   PUBLIC :: helium_print_plength
   PUBLIC :: helium_print_coordinates
   PUBLIC :: helium_print_force
   PUBLIC :: helium_print_accepts
   PUBLIC :: helium_print_perm
   PUBLIC :: helium_print_action
   PUBLIC :: helium_write_cubefile

CONTAINS

! ***************************************************************************
!> \brief  Read XYZ coordinates from file
!> \param coords ...
!> \param file_name ...
!> \param para_env ...
!> \date   2009-06-03
!> \author Lukasz Walewski
!> \note   This is a parallel routine, all ranks get coords defined
! **************************************************************************************************
   SUBROUTINE helium_read_xyz(coords, file_name, para_env)

      REAL(KIND=dp), DIMENSION(:), POINTER               :: coords
      CHARACTER(LEN=default_path_length)                 :: file_name
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_read_xyz', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: strtmp
      INTEGER                                            :: frame, handle, istat, j, natom
      LOGICAL                                            :: found, my_end
      TYPE(cp_parser_type), POINTER                      :: parser

      CALL timeset(routineN, handle)

      ! check if the file can be accessed
      INQUIRE (FILE=file_name, EXIST=found, IOSTAT=istat)
      IF (istat /= 0) THEN
         WRITE (UNIT=strtmp, FMT="(A,I0,A)") &
            "An error occurred inquiring the file <"// &
            TRIM(file_name)//">"
         CPABORT(TRIM(strtmp))
      END IF
      IF (.NOT. found) THEN
         CPABORT("Coordinate file <"//TRIM(file_name)//"> not found.")
      END IF

      NULLIFY (parser)
      CALL parser_create( &
         parser, &
         file_name, &
         para_env=para_env, &
         parse_white_lines=.TRUE.)

      natom = 0
      frame = 0
      CALL parser_get_next_line(parser, 1)
      Frames: DO
         ! Atom number
         CALL parser_get_object(parser, natom)
         frame = frame+1
         IF (frame == 1) THEN
            ALLOCATE (coords(3*natom))
         ELSE
            strtmp = "Warning: more than one frame on file <"//TRIM(file_name)//">"
            CALL helium_write_line(strtmp)
            strtmp = "Warning: using the first frame only!"
            CALL helium_write_line(strtmp)
            EXIT Frames
         END IF
         ! Dummy line
         CALL parser_get_next_line(parser, 2)
         DO j = 1, natom
            ! Atom coordinates
            READ (parser%input_line, *) strtmp, &
               coords(3*(j-1)+1), &
               coords(3*(j-1)+2), &
               coords(3*(j-1)+3)
            coords(3*(j-1)+1) = cp_unit_to_cp2k(coords(3*(j-1)+1), "angstrom")
            coords(3*(j-1)+2) = cp_unit_to_cp2k(coords(3*(j-1)+2), "angstrom")
            coords(3*(j-1)+3) = cp_unit_to_cp2k(coords(3*(j-1)+3), "angstrom")
            ! If there's a white line or end of file exit.. otherwise go on
            CALL parser_get_next_line(parser, 1, at_end=my_end)
            my_end = my_end .OR. (LEN_TRIM(parser%input_line) == 0)
            IF (my_end) THEN
               IF (j /= natom) THEN
                  CPABORT("Error in XYZ format.")
               END IF
               EXIT Frames
            END IF
         END DO
      END DO Frames
      CALL parser_release(parser)

      CALL timestop(handle)

      RETURN
   END SUBROUTINE helium_read_xyz

! ***************************************************************************
!> \brief  Write helium parameters to the output unit
!> \param helium ...
!> \date   2009-06-03
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_write_setup(helium)

      TYPE(helium_solvent_type), POINTER                 :: helium

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_setup', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: msg_str, my_label, stmp, stmp1, stmp2, &
                                                            unit_str
      INTEGER                                            :: i, itmp, j, unit_nr
      INTEGER(KIND=int_8)                                :: i8tmp
      LOGICAL                                            :: first
      REAL(KIND=dp)                                      :: rtmp, v1, v2, v3
      REAL(KIND=dp), DIMENSION(3)                        :: my_abc
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()
      my_label = "HELIUM| "

      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger)

         WRITE (unit_nr, *)
         WRITE (unit_nr, '(T2,A,1X,I0)') TRIM(my_label)// &
            " Number of helium environments:     ", helium%num_env

         WRITE (unit_nr, '(T2,A,1X,I0)') TRIM(my_label)// &
            " Number of solvent atoms:           ", helium%atoms
         WRITE (unit_nr, '(T2,A,1X,I0)') TRIM(my_label)// &
            " Number of solvent beads:           ", helium%beads
         WRITE (unit_nr, '(T2,A,1X,I0)') TRIM(my_label)// &
            " Total number of solvent particles: ", helium%atoms*helium%beads

         unit_str = "angstrom^-3"
         rtmp = cp_unit_from_cp2k(helium%density, &
                                  unit_str)
         WRITE (unit_nr, '(T2,A,F12.6)') TRIM(my_label)//" Density   ["// &
            TRIM(unit_str)//"]:", rtmp

         unit_str = "angstrom"
         rtmp = cp_unit_from_cp2k(helium%cell_size, &
                                  unit_str)
         WRITE (unit_nr, '(T2,A,F12.6)') TRIM(my_label)//" Cell size ["// &
            TRIM(unit_str)//"]:   ", rtmp

         IF (helium%periodic) THEN
            SELECT CASE (helium%cell_shape)
            CASE (helium_cell_shape_cube)
               CALL helium_write_line("PBC cell shape: CUBE.")
            CASE (helium_cell_shape_octahedron)
               CALL helium_write_line("PBC cell shape: TRUNCATED OCTAHEDRON.")
            CASE DEFAULT
               CALL helium_write_line("*** Warning: unknown cell shape.")
            END SELECT
         ELSE
            CALL helium_write_line("PBC turned off.")
         END IF

         IF (helium%droplet_radius .LT. HUGE(1.0_dp)) THEN
            WRITE (stmp, *) helium%droplet_radius*angstrom
            CALL helium_write_line("Droplet radius: "//TRIM(ADJUSTL(stmp))//" [angstrom]")
         END IF

         ! first step gets incremented during first iteration
         WRITE (unit_nr, '(T2,A,1X,I0)') TRIM(my_label)// &
            " First MC step                      :", helium%first_step+1
         WRITE (unit_nr, '(T2,A,1X,I0)') TRIM(my_label)// &
            " Last MC step                       :", helium%last_step
         WRITE (unit_nr, '(T2,A,1X,I0)') TRIM(my_label)// &
            " Total number of MC steps           :", helium%num_steps
         WRITE (unit_nr, '(T2,A,1X,I0)') TRIM(my_label)// &
            " Number of outer MC trials per step :", helium%iter_rot
         WRITE (unit_nr, '(T2,A,1X,I0)') TRIM(my_label)// &
            " Number of inner MC trials per step :", helium%iter_norot
         i8tmp = INT(helium%iter_rot, int_8)
         i8tmp = i8tmp*INT(helium%iter_norot, int_8)
         stmp = ""
         WRITE (stmp, *) i8tmp
         WRITE (unit_nr, '(T2,A)') TRIM(my_label)// &
            " Total number of MC trials per step : "//TRIM(ADJUSTL(stmp))
         i8tmp = INT(helium%num_steps, int_8)
         i8tmp = i8tmp*INT(helium%iter_rot, int_8)
         i8tmp = i8tmp*INT(helium%iter_norot, int_8)
         stmp = ""
         WRITE (stmp, *) i8tmp
         WRITE (unit_nr, '(T2,A)') TRIM(my_label)// &
            " Total number of MC trials          : "//TRIM(ADJUSTL(stmp))

         SELECT CASE (helium%sampling_method)

         CASE (helium_sampling_ceperley)

            ! permutation cycle length sampling
            stmp = ""
            CALL helium_write_line(stmp)
            WRITE (stmp, *) helium%maxcycle
            stmp2 = ""
            WRITE (stmp2, *) "Using maximum permutation cycle length: "// &
               TRIM(ADJUSTL(stmp))
            CALL helium_write_line(stmp2)
            stmp = ""
            WRITE (stmp, *) "Permutation cycle length distribution: "// &
               TRIM(ADJUSTL(m_dist_name(helium%m_dist_type)))
            CALL helium_write_line(stmp)
            stmp = ""
            stmp1 = ""
            WRITE (stmp1, *) helium%m_ratio
            stmp2 = ""
            WRITE (stmp2, *) helium%m_value
            WRITE (stmp, *) "Using ratio "//TRIM(ADJUSTL(stmp1))// &
               " for M = "//TRIM(ADJUSTL(stmp2))
            CALL helium_write_line(stmp)
            stmp = ""
            CALL helium_write_line(stmp)

         CASE DEFAULT
            WRITE (msg_str, *) helium%sampling_method
            msg_str = "Sampling method ("//TRIM(ADJUSTL(msg_str))//") not supported."
            CPABORT(msg_str)

         END SELECT

         ! solute related data
         IF (helium%solute_present) THEN
            WRITE (unit_nr, '(T2,A,1X,I0)') TRIM(my_label)// &
               " Number of solute atoms:            ", helium%solute_atoms
            WRITE (unit_nr, '(T2,A,1X,I0)') TRIM(my_label)// &
               " Number of solute beads:            ", helium%solute_beads
            WRITE (unit_nr, '(T2,A,1X,I0)') TRIM(my_label)// &
               " Total number of solute particles:  ", helium%solute_atoms* &
               helium%solute_beads
            DO i = 1, helium%enum
               WRITE (unit_nr, '(T2,A,1X,A2,A,I0,A,I0,A)', ADVANCE='NO') &
                  TRIM(my_label)//" Solute atom type: ", &
                  helium%ename(i), &
                  "(", &
                  helium%eid(i), &
                  "), count: ", &
                  SIZE(helium%solute_i(helium%eid(i))%iap), &
                  ", "
               first = .TRUE.
               DO j = 1, SIZE(helium%solute_i(helium%eid(i))%iap)
                  IF (first) THEN
                     WRITE (unit_nr, '(A)', ADVANCE='NO') "indices: "
                     first = .FALSE.
                  ELSE
                     WRITE (unit_nr, '(A)', ADVANCE='NO') ", "
                  END IF
                  WRITE (unit_nr, '(I0)', ADVANCE='NO') helium%solute_i(helium%eid(i))%iap(j)
               END DO
               WRITE (unit_nr, '(1X)')
            END DO

            stmp1 = ""
            SELECT CASE (helium%solute_interaction)
            CASE (helium_solute_intpot_none)
               WRITE (stmp1, *) "NONE"
            CASE (helium_solute_intpot_mwater)
               WRITE (stmp1, *) "MWATER"
            CASE DEFAULT
               WRITE (stmp1, *) "***UNKNOWN***"
            END SELECT
            WRITE (unit_nr, '(T2,A,A,A)') &
               TRIM(my_label), &
               " Solute interaction type: ", &
               TRIM(ADJUSTL(stmp1))

            CALL get_cell(helium%solute_cell, abc=my_abc)
            unit_str = "angstrom"
            v1 = cp_unit_from_cp2k(my_abc(1), unit_str)
            v2 = cp_unit_from_cp2k(my_abc(2), unit_str)
            v3 = cp_unit_from_cp2k(my_abc(3), unit_str)
            WRITE (unit_nr, '(T2,A,F12.6,1X,F12.6,1X,F12.6)') &
               TRIM(my_label)//" Solute cell size ["// &
               TRIM(unit_str)//"]:   ", v1, v2, v3
         ELSE
            WRITE (unit_nr, '(T2,A)') TRIM(my_label)//" Solute is NOT present"
         END IF
      END IF
      CALL helium_write_line("")

      ! RDF related settings
      IF (helium%rdf_present) THEN
         WRITE (stmp, '(I6)') helium%rdf_num_ctr
         CALL helium_write_line("RDF| number of centers: "//TRIM(stmp))
         rtmp = cp_unit_from_cp2k(helium%rdf_delr, "angstrom")
         WRITE (stmp, '(1X,F12.6)') rtmp
         CALL helium_write_line("RDF| delr [angstrom]  : "//TRIM(stmp))
         rtmp = cp_unit_from_cp2k(helium%rdf_maxr, "angstrom")
         WRITE (stmp, '(1X,F12.6)') rtmp
         CALL helium_write_line("RDF| maxr [angstrom]  : "//TRIM(stmp))
         itmp = helium%rdf_nbin
         WRITE (stmp, '(I6)') itmp
         CALL helium_write_line("RDF| nbin             : "//TRIM(stmp))
         CALL helium_write_line("")
      ELSE
         CALL helium_write_line("RDF| radial distributions will NOT be calculated.")
         CALL helium_write_line("")
      END IF

      ! RHO related settings
      IF (helium%rho_present) THEN
         CALL helium_write_line('RHO| The following densities will be calculated:')
         DO i = 1, rho_num
            IF (helium%rho_property(i)%is_calculated) THEN
               WRITE (stmp, '(A)') 'RHO|    '//TRIM(helium%rho_property(i)%name)
               CALL helium_write_line(TRIM(stmp))
            END IF
         END DO
         CALL helium_write_line('RHO|')
         rtmp = cp_unit_from_cp2k(helium%rho_delr, "angstrom")
         WRITE (stmp, '(1X,F12.6)') rtmp
         CALL helium_write_line("RHO| delr [angstrom]  : "//TRIM(stmp))
         rtmp = cp_unit_from_cp2k(helium%rho_maxr, "angstrom")
         WRITE (stmp, '(1X,F12.6)') rtmp
         CALL helium_write_line("RHO| maxr [angstrom]  : "//TRIM(stmp))
         itmp = helium%rho_nbin
         WRITE (stmp, '(I6)') itmp
         CALL helium_write_line("RHO| nbin             : "//TRIM(stmp))
         CALL helium_write_line("")
      ELSE
         CALL helium_write_line("RHO| Density distributions will NOT be calculated.")
         CALL helium_write_line("")
      END IF

      RETURN
   END SUBROUTINE helium_write_setup

! ***************************************************************************
!> \brief  Writes out a line of text to the default output unit
!> \param line ...
!> \date   2009-07-10
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_write_line(line)

      CHARACTER(len=*), INTENT(IN)                       :: line

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_line', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: my_label
      INTEGER                                            :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()
      my_label = "HELIUM|"

      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger)
         WRITE (unit_nr, '(T2,A)') TRIM(my_label)//" "//TRIM(line)
      END IF

      RETURN
   END SUBROUTINE helium_write_line

! ***************************************************************************
!> \brief  Print energies according to HELIUM%PRINT%ENERGY
!> \param helium_env ...
!> \date   2009-06-08
!> \par    History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_print_energy(helium_env)

      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_print_energy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iteration, k, m, unit_nr
      LOGICAL                                            :: file_is_new, should_output
      REAL(KIND=dp)                                      :: naccptd
      REAL(KIND=dp), DIMENSION(:), POINTER               :: my_energy
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      IF (logger%para_env%ionode) THEN
         print_key => section_vals_get_subs_vals(helium_env(1)%helium%input, &
                                                 "MOTION%PINT%HELIUM%PRINT%ENERGY")
         should_output = BTEST(cp_print_key_should_output( &
                               iteration_info=logger%iter_info, &
                               basis_section=print_key), cp_p_file)
      END IF
      CALL mp_bcast(should_output, logger%para_env%source, helium_env(1)%comm)

      IF (should_output) THEN

         naccptd = 0.0_dp
         DO k = 1, SIZE(helium_env)
            DO m = 1, helium_env(k)%helium%maxcycle
               naccptd = naccptd+helium_env(k)%helium%num_accepted(helium_env(k)%helium%bisctlog2+2, m)
            END DO
         END DO
         CALL mp_sum(naccptd, helium_env(1)%comm)
         naccptd = naccptd/REAL(helium_env(1)%helium%num_env, dp)

         IF (logger%para_env%ionode) THEN
            my_energy => helium_env(1)%helium%energy_avrg
            iteration = logger%iter_info%iteration(2)

            unit_nr = cp_print_key_unit_nr( &
                      logger, &
                      print_key, &
                      middle_name="helium-energy", &
                      extension=".dat", &
                      is_new_file=file_is_new)

            IF (file_is_new) THEN
               WRITE (unit_nr, '(A9,1X,A12,6(1X,A20))') &
                  "#    Step", &
                  "     Naccptd", &
                  "               E_pot", &
                  "               E_kin", &
                  "            E_thermo", &
                  "            E_virial", &
                  "             E_inter", &
                  "               E_tot"
            END IF

            WRITE (unit_nr, "(I9,1X,F12.1,6(1X,F20.9))") &
               iteration, &
               naccptd, &
               my_energy(e_id_potential), &
               my_energy(e_id_kinetic), &
               my_energy(e_id_thermo), &
               my_energy(e_id_virial), &
               my_energy(e_id_interact), &
               my_energy(e_id_total)

            CALL m_flush(unit_nr)
            CALL cp_print_key_finished_output(unit_nr, logger, print_key)

         END IF
      END IF

      CALL timestop(handle)

      RETURN
   END SUBROUTINE helium_print_energy

! ***************************************************************************
!> \brief  Print a 3D real vector accoring to printkey <pkey>
!> \param helium_env ...
!> \param pkey ...
!> \param DATA ...
!> \param uconv ...
!> \param col_label ...
!> \param cmmnt ...
!> \param fname ...
!> \param fpos ...
!> \param avg ...
!> \date   2014-09-09
!> \par    History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_print_vector(helium_env, pkey, DATA, uconv, col_label, cmmnt, fname, fpos, avg)

      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env
      CHARACTER(len=*)                                   :: pkey
      REAL(KIND=dp), DIMENSION(:), POINTER               :: DATA
      REAL(KIND=dp)                                      :: uconv
      CHARACTER(len=*), DIMENSION(3)                     :: col_label
      CHARACTER(len=*)                                   :: cmmnt, fname
      CHARACTER(len=*), OPTIONAL                         :: fpos
      LOGICAL, OPTIONAL                                  :: avg

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_print_vector', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: my_fpos
      INTEGER                                            :: handle, i, irank, msglen, nenv, offset, &
                                                            unit_nr
      LOGICAL                                            :: is_new, my_avg, should_output
      REAL(KIND=dp), DIMENSION(3)                        :: avg_data
      REAL(KIND=dp), DIMENSION(:), POINTER               :: data_p
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: psection

      CALL timeset(routineN, handle)

      IF (PRESENT(avg)) THEN
         my_avg = avg
      ELSE
         my_avg = .FALSE.
      END IF

      IF (PRESENT(fpos)) THEN
         my_fpos = fpos
      ELSE
         my_fpos = "APPEND"
      END IF

      NULLIFY (logger, psection)
      logger => cp_get_default_logger()

      psection => section_vals_get_subs_vals(helium_env(1)%helium%input, pkey)
      should_output = BTEST(cp_print_key_should_output( &
                            iteration_info=logger%iter_info, &
                            basis_section=psection), cp_p_file)

      IF (.NOT. should_output) THEN
         CALL timestop(handle)
         RETURN
      END IF

      IF (my_avg) THEN
         ! average data over all processors and environments
         avg_data(:) = 0.0_dp
         msglen = SIZE(avg_data)
         DO i = 0, SIZE(helium_env)-1
            avg_data(:) = avg_data(:)+DATA(msglen*i+1:msglen*(i+1))
         END DO
         CALL mp_sum(avg_data, helium_env(1)%comm)
         nenv = helium_env(1)%helium%num_env
         avg_data(:) = avg_data(:)/REAL(nenv, dp)
      ELSE
         ! gather data from all processors
         offset = 0
         DO i = 1, logger%para_env%mepos
            offset = offset+helium_env(1)%env_all(i)
         END DO

         helium_env(1)%helium%rtmp_3_np_1d = 0.0_dp
         msglen = SIZE(avg_data)
         DO i = 0, SIZE(helium_env)-1
            helium_env(1)%helium%rtmp_3_np_1d(msglen*(offset+i)+1:msglen*(offset+i+1)) = DATA(msglen*i+1:msglen*(i+1))
         END DO
         CALL mp_sum(helium_env(1)%helium%rtmp_3_np_1d, &
                     helium_env(1)%comm)
      END IF

      unit_nr = cp_print_key_unit_nr( &
                logger, &
                psection, &
                middle_name=fname, &
                extension=".dat", &
                file_position=my_fpos, &
                is_new_file=is_new)

      IF (logger%para_env%ionode) THEN

         IF (is_new) THEN
            ! comment
            IF (cmmnt .NE. "") THEN
               WRITE (unit_nr, '(A)') "# "//cmmnt
            END IF
            ! column labels
            WRITE (unit_nr, '(A2,A18,1X,A20,1X,A20)') &
               "# ", &
               col_label(1), &
               col_label(2), &
               col_label(3)
         END IF

         IF (my_avg) THEN
            DO i = 1, 3
               WRITE (unit_nr, '(E27.20)', ADVANCE='NO') uconv*avg_data(i)
               IF (i .LT. 3) THEN
                  WRITE (unit_nr, '(1X)', ADVANCE='NO')
               END IF
            END DO
            WRITE (unit_nr, '(A)') ""
         ELSE
            ! iterate over processors/helium environments
            DO irank = 1, helium_env(1)%helium%num_env
               ! unpack data (actually point to the right fragment only)
               msglen = SIZE(avg_data)
               offset = (irank-1)*msglen
               data_p => helium_env(1)%helium%rtmp_3_np_1d(offset+1:offset+msglen)
               ! write out the data
               DO i = 1, 3
                  WRITE (unit_nr, '(E27.20)', ADVANCE='NO') uconv*data_p(i)
                  IF (i .LT. 3) THEN
                     WRITE (unit_nr, '(1X)', ADVANCE='NO')
                  END IF
               END DO
               WRITE (unit_nr, '(A)') ""
            END DO
         END IF

         CALL m_flush(unit_nr)
         CALL cp_print_key_finished_output(unit_nr, logger, psection)

      END IF

      CALL timestop(handle)

      RETURN
   END SUBROUTINE helium_print_vector

! ***************************************************************************
!> \brief  Print acceptance counts according to HELIUM%PRINT%ACCEPTS
!> \param helium_env ...
!> \date   2010-05-27
!> \par    History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_print_accepts(helium_env)

      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_print_accepts', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, unit_nr
      LOGICAL                                            :: file_is_new, should_output
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      IF (logger%para_env%ionode) THEN
         print_key => section_vals_get_subs_vals(helium_env(1)%helium%input, &
                                                 "MOTION%PINT%HELIUM%PRINT%ACCEPTS")
         should_output = BTEST(cp_print_key_should_output( &
                               iteration_info=logger%iter_info, &
                               basis_section=print_key), cp_p_file)

         IF (should_output) THEN
            unit_nr = cp_print_key_unit_nr( &
                      logger, &
                      print_key, &
                      middle_name="helium-accepts", &
                      extension=".dat", &
                      is_new_file=file_is_new)

            IF (file_is_new) THEN
               WRITE (unit_nr, '(A8,1X,A15,1X,A20)', ADVANCE='NO') &
                  "# Length", &
                  "         Trials", &
                  "            Selected"
               DO j = 1, helium_env(1)%helium%bisctlog2
                  WRITE (unit_nr, '(A17,1X,I3)', ADVANCE='NO') "            Level", j
               END DO
               WRITE (unit_nr, '(A)') ""
            END IF

            DO i = 1, helium_env(1)%helium%maxcycle
               WRITE (unit_nr, '(I3)', ADVANCE='NO') i
               DO j = 1, helium_env(1)%helium%bisctlog2+2
                  WRITE (unit_nr, '(1X,F20.2)', ADVANCE='NO') helium_env(1)%helium%num_accepted(j, i)
               END DO
               WRITE (unit_nr, '(A)') ""
            END DO
            WRITE (unit_nr, '(A)') "&"

            CALL m_flush(unit_nr)
            CALL cp_print_key_finished_output(unit_nr, logger, print_key)

         END IF
      END IF

      CALL timestop(handle)
      RETURN
   END SUBROUTINE helium_print_accepts

! ***************************************************************************
!> \brief  Print permutation state according to HELIUM%PRINT%PERM
!> \param helium_env ...
!> \date   2010-06-07
!> \par    History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_print_perm(helium_env)

      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_print_perm', &
         routineP = moduleN//':'//routineN

      CHARACTER                                          :: left_delim, right_delim
      CHARACTER(len=default_string_length)               :: msg_str, my_middle_name, stmp
      INTEGER                                            :: curat, handle, i, irank, j, lb, msglen, &
                                                            nused, offset, outformat, ub, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: my_cycle, my_perm, used_indices
      LOGICAL                                            :: first, first_cycle, should_output
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      ! determine offset for arrays
      offset = 0
      DO i = 1, logger%para_env%mepos
         offset = offset+helium_env(1)%env_all(i)
      END DO

      print_key => section_vals_get_subs_vals(helium_env(1)%helium%input, &
                                              "MOTION%PINT%HELIUM%PRINT%PERM")
      should_output = BTEST(cp_print_key_should_output( &
                            iteration_info=logger%iter_info, &
                            basis_section=print_key), cp_p_file)

      IF (.NOT. should_output) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! get the output type
      CALL section_vals_val_get(print_key, "FORMAT", i_val=outformat)
      IF (outformat .EQ. perm_cycle) THEN
         ! gather positions from all processors
         helium_env(1)%helium%rtmp_3_atoms_beads_np_1d = 0.0_dp
         j = SIZE(helium_env(1)%helium%rtmp_3_atoms_beads_1d)
         DO i = 1, SIZE(helium_env)
            helium_env(1)%helium%rtmp_3_atoms_beads_np_1d(j*(offset+i-1):) = &
               PACK(helium_env(i)%helium%pos(:, :, 1:helium_env(1)%helium%beads), .TRUE.)
         END DO
         CALL mp_sum(helium_env(1)%helium%rtmp_3_atoms_beads_np_1d, &
                     helium_env(1)%comm)
         ! set logical mask for unpacking coordinates gathered from other ranks
         helium_env(1)%helium%ltmp_3_atoms_beads_3d(:, :, :) = .TRUE.
      END IF

      ! gather permutation state from all processors to logger%para_env%source
      helium_env(1)%helium%itmp_atoms_np_1d(:) = 0
      DO i = 1, SIZE(helium_env)
         helium_env(1)%helium%itmp_atoms_np_1d(offset+i:) = helium_env(i)%helium%permutation
      END DO

      CALL mp_sum(helium_env(1)%helium%itmp_atoms_np_1d, &
                  helium_env(1)%comm)

      IF (logger%para_env%ionode) THEN

         ! iterate over helium environments
         DO irank = 1, helium_env(1)%helium%num_env

            ! generate one file per environment
            stmp = ""
            WRITE (stmp, *) irank
            my_middle_name = "helium-perm-"//TRIM(ADJUSTL(stmp))
            unit_nr = cp_print_key_unit_nr( &
                      logger, &
                      print_key, &
                      middle_name=TRIM(my_middle_name), &
                      extension=".dat")

            ! unpack permutation state (actually point to the right section only)
            lb = (irank-1)*helium_env(1)%helium%atoms+1
            ub = irank*helium_env(1)%helium%atoms
            my_perm => helium_env(1)%helium%itmp_atoms_np_1d(lb:ub)

            SELECT CASE (outformat)

            CASE (perm_cycle)
               ! write the permutation grouped into cycles

               ! unpack coordinates (necessary only for winding path delimiters)
               msglen = SIZE(helium_env(1)%helium%rtmp_3_atoms_beads_1d)
               offset = (irank-1)*msglen
               helium_env(1)%helium%work(:, :, 1:helium_env(1)%helium%beads) = &
                  UNPACK(helium_env(1)%helium%rtmp_3_atoms_beads_np_1d(offset+1:offset+msglen), &
                         MASK=helium_env(1)%helium%ltmp_3_atoms_beads_3d, FIELD=0.0_dp)

               curat = 1
               nused = 0
               first_cycle = .TRUE.
               DO WHILE (curat .LE. helium_env(1)%helium%atoms)

                  ! get the permutation cycle the current atom belongs to
                  my_cycle => helium_cycle_of(curat, my_perm)

                  ! include the current cycle in the pool of "used" indices
                  nused = nused+SIZE(my_cycle)
                  CALL reallocate(used_indices, 1, nused)
                  used_indices(nused-SIZE(my_cycle)+1:nused) = my_cycle(1:SIZE(my_cycle))

                  ! select delimiters accoring to the cycle's winding state
                  IF (helium_is_winding(helium_env(1)%helium, curat, helium_env(1)%helium%work, my_perm)) THEN
                     left_delim = "["
                     right_delim = "]"
                  ELSE
                     left_delim = "("
                     right_delim = ")"
                  END IF

                  ! cycle delimiter
                  IF (first_cycle) THEN
                     first_cycle = .FALSE.
                  ELSE
                     WRITE (unit_nr, '(1X)', ADVANCE='NO')
                  END IF

                  ! write out the current cycle
                  WRITE (unit_nr, '(A1)', ADVANCE='NO') left_delim
                  first = .TRUE.
                  DO i = 1, SIZE(my_cycle)
                     IF (first) THEN
                        first = .FALSE.
                     ELSE
                        WRITE (unit_nr, '(1X)', ADVANCE='NO')
                     END IF
                     WRITE (unit_nr, '(I0)', ADVANCE='NO') my_cycle(i)
                  END DO
                  WRITE (unit_nr, '(A1)', ADVANCE='NO') right_delim

                  ! clean up
                  DEALLOCATE (my_cycle)
                  NULLIFY (my_cycle)

                  ! try to increment the current atom index
                  DO WHILE (ANY(used_indices .EQ. curat))
                     curat = curat+1
                  END DO

               END DO
               WRITE (unit_nr, *)

               DEALLOCATE (used_indices)
               NULLIFY (used_indices)

            CASE (perm_plain)
               ! write the plain permutation

               first = .TRUE.
               DO i = 1, helium_env(1)%helium%atoms
                  IF (first) THEN
                     first = .FALSE.
                  ELSE
                     WRITE (unit_nr, '(1X)', ADVANCE='NO')
                  END IF
                  WRITE (unit_nr, '(I0)', ADVANCE='NO') my_perm(i)
               END DO
               WRITE (unit_nr, '(A)') ""

            CASE default

               msg_str = "Format for permutation output unknown."
               CPABORT(msg_str)

            END SELECT

            CALL m_flush(unit_nr)
            CALL cp_print_key_finished_output(unit_nr, logger, print_key)

         END DO

      END IF

      CALL timestop(handle)

      RETURN
   END SUBROUTINE helium_print_perm

! **************************************************************************************************
!> \brief Print helium action file to HELIUM%PRINT%ACTION
!> \param pint_env ...
!> \param helium_env ...
!> \date   2016-06-07
!> \par    History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Felix Uhl
! **************************************************************************************************
   SUBROUTINE helium_print_action(pint_env, helium_env)

      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_print_action', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: my_middle_name, stmp
      INTEGER                                            :: handle, i, irank, iteration, k, offset, &
                                                            unit_nr
      LOGICAL                                            :: file_is_new, should_output
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tmp_inter_action, tmp_link_action, &
                                                            tmp_pair_action
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      ! determine offset for arrays
      offset = 0
      DO i = 1, logger%para_env%mepos
         offset = offset+helium_env(1)%env_all(i)
      END DO

      print_key => section_vals_get_subs_vals(helium_env(1)%helium%input, &
                                              "MOTION%PINT%HELIUM%PRINT%ACTION")
      should_output = BTEST(cp_print_key_should_output( &
                            iteration_info=logger%iter_info, &
                            basis_section=print_key), cp_p_file)

      IF (.NOT. should_output) THEN
         CALL timestop(handle)
         RETURN
      END IF

      DO k = 1, SIZE(helium_env)
         helium_env(k)%helium%link_action = helium_total_link_action(helium_env(k)%helium)
         helium_env(k)%helium%pair_action = helium_total_pair_action(helium_env(k)%helium)
         helium_env(k)%helium%inter_action = helium_total_inter_action(pint_env, helium_env(k)%helium)
      END DO

      NULLIFY (tmp_link_action)
      NULLIFY (tmp_pair_action)
      NULLIFY (tmp_inter_action)
      ALLOCATE (tmp_link_action(helium_env(1)%helium%num_env))
      ALLOCATE (tmp_pair_action(helium_env(1)%helium%num_env))
      ALLOCATE (tmp_inter_action(helium_env(1)%helium%num_env))
      tmp_link_action(:) = 0.0_dp
      tmp_pair_action(:) = 0.0_dp
      tmp_inter_action(:) = 0.0_dp
      ! gather Action from all processors to logger%para_env%source
      DO k = 1, SIZE(helium_env)
         tmp_link_action(offset+k) = helium_env(k)%helium%link_action
         tmp_pair_action(offset+k) = helium_env(k)%helium%pair_action
         tmp_inter_action(offset+k) = helium_env(k)%helium%inter_action
      END DO
      CALL mp_sum(tmp_link_action, helium_env(1)%comm)
      CALL mp_sum(tmp_pair_action, helium_env(1)%comm)
      CALL mp_sum(tmp_inter_action, helium_env(1)%comm)

      IF (logger%para_env%ionode) THEN
         iteration = logger%iter_info%iteration(2)
         ! iterate over processors/helium environments
         DO irank = 1, helium_env(1)%helium%num_env

            ! generate one file per helium_env
            stmp = ""
            WRITE (stmp, *) irank
            my_middle_name = "helium-action-"//TRIM(ADJUSTL(stmp))
            unit_nr = cp_print_key_unit_nr( &
                      logger, &
                      print_key, &
                      middle_name=TRIM(my_middle_name), &
                      extension=".dat", &
                      is_new_file=file_is_new)

            IF (file_is_new) THEN
               WRITE (unit_nr, '(A9,3(1X,A25))') &
                  "#    Step", &
                  "     He_Total_Link_Action", &
                  "     He_Total_Pair_Action", &
                  "     He_Total_Interaction"
            END IF

            WRITE (unit_nr, "(I9,3(1X,F25.14))") &
               iteration, &
               tmp_link_action(irank), &
               tmp_pair_action(irank), &
               tmp_inter_action(irank)

            CALL m_flush(unit_nr)
            CALL cp_print_key_finished_output(unit_nr, logger, print_key)

         END DO
      END IF

      DEALLOCATE (tmp_link_action)
      DEALLOCATE (tmp_pair_action)
      DEALLOCATE (tmp_inter_action)
      NULLIFY (tmp_link_action)
      NULLIFY (tmp_pair_action)
      NULLIFY (tmp_inter_action)

      CALL timestop(handle)

      RETURN
   END SUBROUTINE helium_print_action

! ***************************************************************************
!> \brief  Print coordinates according to HELIUM%PRINT%COORDINATES
!> \param helium_env ...
!> \date   2009-07-16
!> \par    History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_print_coordinates(helium_env)

      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_print_coordinates', &
         routineP = moduleN//':'//routineN

      CHARACTER(3)                                       :: resName
      CHARACTER(len=default_string_length)               :: fmt_string, my_middle_name, stmp
      INTEGER                                            :: handle, i, ia, ib, ib1, ib2, ic, icycle, &
                                                            irank, j, k, msglen, offset, tmp1, &
                                                            tmp2, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: my_perm
      LOGICAL                                            :: are_connected, is_winding, ltmp, &
                                                            should_output
      REAL(KIND=dp)                                      :: xtmp, ytmp, ztmp
      REAL(KIND=dp), DIMENSION(3)                        :: r0, r1, r2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      ! determine offset for arrays
      offset = 0
      DO i = 1, logger%para_env%mepos
         offset = offset+helium_env(1)%env_all(i)
      END DO

      print_key => section_vals_get_subs_vals(helium_env(1)%helium%input, &
                                              "MOTION%PINT%HELIUM%PRINT%COORDINATES")
      should_output = BTEST(cp_print_key_should_output( &
                            iteration_info=logger%iter_info, &
                            basis_section=print_key), cp_p_file)

      IF (.NOT. should_output) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! prepare the coordinates for output (use unit cell centered around r0)
      DO k = 1, SIZE(helium_env)
         r0(:) = helium_env(k)%helium%center(:)
         DO ia = 1, helium_env(k)%helium%atoms
            DO ib = 1, helium_env(k)%helium%beads
               r1(:) = helium_env(k)%helium%pos(:, ia, ib)-r0(:)
               r2(:) = helium_env(k)%helium%pos(:, ia, ib)-r0(:)
               CALL helium_pbc(helium_env(k)%helium, r2)
               ltmp = .FALSE.
               DO ic = 1, 3
                  IF (ABS(r1(ic)-r2(ic)) .GT. 100.0_dp*EPSILON(0.0_dp)) THEN
                     ltmp = .TRUE.
                     CYCLE
                  END IF
               END DO
               IF (ltmp) THEN
                  helium_env(k)%helium%work(:, ia, ib) = r0(:)+r2(:)
               ELSE
                  helium_env(k)%helium%work(:, ia, ib) = helium_env(k)%helium%pos(:, ia, ib)
               END IF
            END DO
         END DO
      END DO

      ! gather positions from all processors to logger%para_env%source
      helium_env(1)%helium%rtmp_3_atoms_beads_np_1d = 0.0_dp
      j = SIZE(helium_env(1)%helium%rtmp_3_atoms_beads_1d)
      DO i = 1, SIZE(helium_env)
         helium_env(1)%helium%rtmp_3_atoms_beads_np_1d(j*(offset+i-1)+1:j*(offset+i)) = &
            PACK(helium_env(i)%helium%pos(:, :, 1:helium_env(1)%helium%beads), .TRUE.)
      END DO
      CALL mp_sum(helium_env(1)%helium%rtmp_3_atoms_beads_np_1d, &
                  helium_env(1)%comm)

      ! gather permutation state from all processors to logger%para_env%source
      helium_env(1)%helium%itmp_atoms_np_1d(:) = 0
      DO i = 1, SIZE(helium_env)
         helium_env(1)%helium%itmp_atoms_np_1d((offset+i-1):) = helium_env(i)%helium%permutation
      END DO

      CALL mp_sum(helium_env(1)%helium%itmp_atoms_np_1d, &
                  helium_env(1)%comm)

      ! set logical mask for unpacking coordinates gathered from other ranks
      helium_env(1)%helium%ltmp_3_atoms_beads_3d(:, :, :) = .TRUE.

      IF (logger%para_env%ionode) THEN

         ! iterate over helium environments
         DO irank = 1, helium_env(1)%helium%num_env

            ! generate one file per environment
            stmp = ""
            WRITE (stmp, *) irank
            my_middle_name = "helium-pos-"//TRIM(ADJUSTL(stmp))
            unit_nr = cp_print_key_unit_nr( &
                      logger, &
                      print_key, &
                      middle_name=TRIM(my_middle_name), &
                      extension=".pdb")

            ! write out the unit cell parameters
            fmt_string = "(A6,3F9.3,3F7.2,1X,A11,1X,I3)"
            xtmp = helium_env(1)%helium%cell_size
            xtmp = cp_unit_from_cp2k(xtmp, "angstrom")
            SELECT CASE (helium_env (1)%helium%cell_shape)
            CASE (helium_cell_shape_cube)
               stmp = "C          "
            CASE (helium_cell_shape_octahedron)
               stmp = "O          "
            CASE DEFAULT
               stmp = "UNKNOWN    "
            END SELECT
            WRITE (unit_nr, fmt_string) "CRYST1", &
               xtmp, xtmp, xtmp, &
               90.0_dp, 90.0_dp, 90.0_dp, &
               stmp, helium_env(1)%helium%beads

            ! unpack coordinates
            msglen = SIZE(helium_env(1)%helium%rtmp_3_atoms_beads_1d)
            offset = (irank-1)*msglen
            helium_env(1)%helium%work(:, :, 1:helium_env(1)%helium%beads) = &
               UNPACK(helium_env(1)%helium%rtmp_3_atoms_beads_np_1d(offset+1:offset+msglen), &
                      MASK=helium_env(1)%helium%ltmp_3_atoms_beads_3d, FIELD=0.0_dp)

            ! unpack permutation state (actually point to the right section only)
            msglen = SIZE(helium_env(1)%helium%permutation)
            offset = (irank-1)*msglen
            my_perm => helium_env(1)%helium%itmp_atoms_np_1d(offset+1:offset+msglen)

            ! write out coordinates
            fmt_string = &
               "(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)"
            DO ia = 1, helium_env(1)%helium%atoms
               icycle = helium_cycle_number(helium_env(1)%helium, ia, my_perm)
               is_winding = helium_is_winding(helium_env(1)%helium, ia, helium_env(1)%helium%work, my_perm)
               IF (is_winding) THEN
                  resName = "SPR"
               ELSE
                  resName = "NRM"
               END IF
               DO ib = 1, helium_env(1)%helium%beads
                  xtmp = helium_env(1)%helium%work(1, ia, ib)
                  xtmp = cp_unit_from_cp2k(xtmp, "angstrom")
                  ytmp = helium_env(1)%helium%work(2, ia, ib)
                  ytmp = cp_unit_from_cp2k(ytmp, "angstrom")
                  ztmp = helium_env(1)%helium%work(3, ia, ib)
                  ztmp = cp_unit_from_cp2k(ztmp, "angstrom")
                  WRITE (unit_nr, fmt_string) "ATOM  ", &
                     (ia-1)*helium_env(1)%helium%beads+ib, &
                     " He ", " ", resName, "X", &
                     icycle, &
                     " ", &
                     xtmp, ytmp, ztmp, &
                     1.0_dp, 0.0_dp, "HE", "  "
               END DO
            END DO

            ! write out the bead connectivity information
            DO ia = 1, helium_env(1)%helium%atoms

               ! write connectivity records for this atom only if the path
               ! it belongs to is longer than 1.
               IF (helium_path_length(helium_env(1)%helium, ia, my_perm) .LE. 1) THEN
                  CYCLE
               END IF

               DO ib = 1, helium_env(1)%helium%beads-1
                  ! check wheather the consecutive beads belong to the same box
                  r1(:) = helium_env(1)%helium%work(:, ia, ib)-helium_env(1)%helium%work(:, ia, ib+1)
                  r2(:) = r1(:)
                  CALL helium_pbc(helium_env(1)%helium, r2)
                  are_connected = .TRUE.
                  DO ic = 1, 3
                     IF (ABS(r1(ic)-r2(ic)) .GT. 100.0_dp*EPSILON(0.0_dp)) THEN
                        ! if the distance betw ib and ib+1 changes upon applying
                        ! PBC do not connect them
                        are_connected = .FALSE.
                        CYCLE
                     END IF
                  END DO
                  IF (are_connected) THEN
                     tmp1 = (ia-1)*helium_env(1)%helium%beads+ib
                     tmp2 = (ia-1)*helium_env(1)%helium%beads+ib+1
                     ! smaller value has to go first
                     IF (tmp1 .LT. tmp2) THEN
                        ib1 = tmp1
                        ib2 = tmp2
                     ELSE
                        ib1 = tmp2
                        ib2 = tmp1
                     END IF
                     WRITE (unit_nr, '(A6,2I5)') "CONECT", ib1, ib2
                  END IF
               END DO

               ! last bead of atom <ia> connects to the first bead
               ! of the next atom in the permutation cycle
               r1(:) = helium_env(1)%helium%work(:, ia, helium_env(1)%helium%beads)-helium_env(1)%helium%work(:, my_perm(ia), 1)
               r2(:) = r1(:)
               CALL helium_pbc(helium_env(1)%helium, r2)
               are_connected = .TRUE.
               DO ic = 1, 3
                  IF (ABS(r1(ic)-r2(ic)) .GT. 100.0_dp*EPSILON(0.0_dp)) THEN
                     ! if the distance betw ib and ib+1 changes upon applying
                     ! PBC do not connect them
                     are_connected = .FALSE.
                     CYCLE
                  END IF
               END DO
               IF (are_connected) THEN
                  tmp1 = ia*helium_env(1)%helium%beads
                  tmp2 = (my_perm(ia)-1)*helium_env(1)%helium%beads+1
                  IF (tmp1 .LT. tmp2) THEN
                     ib1 = tmp1
                     ib2 = tmp2
                  ELSE
                     ib1 = tmp2
                     ib2 = tmp1
                  END IF
                  WRITE (unit_nr, '(A6,2I5)') "CONECT", ib1, ib2
               END IF
            END DO
            WRITE (unit_nr, '(A)') "END"

            CALL m_flush(unit_nr)
            CALL cp_print_key_finished_output(unit_nr, logger, print_key)

         END DO

      END IF

      CALL timestop(handle)

      RETURN
   END SUBROUTINE helium_print_coordinates

! ***************************************************************************
!> \brief  Print radial distribution functions according to HELIUM%PRINT%RDF
!> \param helium_env ...
!> \date   2009-07-23
!> \par    History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_print_rdf(helium_env)

      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_print_rdf', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: comment, stmp
      INTEGER                                            :: handle, ia, ic, id, itmp, iweight, k, &
                                                            nsteps, unit_nr
      LOGICAL                                            :: should_output
      REAL(KIND=dp)                                      :: inv_norm, rtmp
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: message
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      IF (logger%para_env%ionode) THEN
         print_key => section_vals_get_subs_vals(helium_env(1)%helium%input, &
                                                 "MOTION%PINT%HELIUM%PRINT%RDF")
         should_output = BTEST(cp_print_key_should_output( &
                               iteration_info=logger%iter_info, &
                               basis_section=print_key), cp_p_file)
      END IF
      CALL mp_bcast(should_output, logger%para_env%source, helium_env(1)%comm)

      IF (should_output) THEN
         ! work on the temporary array so that accumulated data remains intact
         ! save accumulated data of different env on same core in first temp
         helium_env(1)%helium%rdf_inst(:, :, :) = 0.0_dp
         DO k = 1, SIZE(helium_env)
            helium_env(1)%helium%rdf_inst(:, :, :) = helium_env(1)%helium%rdf_inst(:, :, :)+ &
                                                     helium_env(k)%helium%rdf_accu(:, :, :)
         END DO

         ! average over processors
         NULLIFY (message)
         message => helium_env(1)%helium%rdf_inst(:, :, :)
         CALL mp_sum(message, helium_env(1)%comm)
         itmp = helium_env(1)%helium%num_env
         inv_norm = 1.0_dp/REAL(itmp, dp)
         helium_env(1)%helium%rdf_inst(:, :, :) = helium_env(1)%helium%rdf_inst(:, :, :)*inv_norm

         ! average over steps performed so far in this run
         nsteps = helium_env(1)%helium%current_step-helium_env(1)%helium%first_step
         inv_norm = 1.0_dp/REAL(nsteps, dp)
         helium_env(1)%helium%rdf_inst(:, :, :) = helium_env(1)%helium%rdf_inst(:, :, :)*inv_norm

         iweight = helium_env(1)%helium%rdf_iweight
         ! average over the old and the current density (observe the weights!)
         helium_env(1)%helium%rdf_inst(:, :, :) = nsteps*helium_env(1)%helium%rdf_inst(:, :, :)+ &
                                                  iweight*helium_env(1)%helium%rdf_rstr(:, :, :)
         helium_env(1)%helium%rdf_inst(:, :, :) = helium_env(1)%helium%rdf_inst(:, :, :)/REAL(nsteps+iweight, dp)

         IF (logger%para_env%ionode) THEN

            DO id = 1, helium_env(1)%helium%rdf_num ! loop over estimators ---

               ! overwrite RDF file each time it is written
               stmp = ""
               WRITE (stmp, *) id
               unit_nr = cp_print_key_unit_nr( &
                         logger, &
                         print_key, &
                         middle_name="helium-rdf-"//TRIM(ADJUSTL(stmp)), &
                         extension=".dat", &
                         file_position="REWIND", &
                         do_backup=.FALSE.)

               comment = "Radial distribution function"

               DO ic = 1, helium_env(1)%helium%rdf_nbin
                  rtmp = (REAL(ic, dp)-0.5_dp)*helium_env(1)%helium%rdf_delr
                  rtmp = cp_unit_from_cp2k(rtmp, "angstrom")
                  WRITE (unit_nr, '(F20.10)', ADVANCE='NO') rtmp
                  DO ia = 1, helium_env(1)%helium%rdf_num_ctr
                     WRITE (unit_nr, '(F20.10)', ADVANCE='NO') helium_env(1)%helium%rdf_inst(id, ic, ia)
                  END DO
                  WRITE (unit_nr, *)
               END DO

               CALL m_flush(unit_nr)
               CALL cp_print_key_finished_output(unit_nr, logger, print_key)

            END DO ! loop over estimators ---

         END IF
      END IF

      CALL timestop(handle)

      RETURN
   END SUBROUTINE helium_print_rdf

! ***************************************************************************
!> \brief  Print densities according to HELIUM%PRINT%RHO
!> \param helium_env ...
!> \date   2011-06-21
!> \par History
!>      08.2015 cleaned code from unneded arrays
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_print_rho(helium_env)

      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_print_rho', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: comment, fname
      INTEGER                                            :: handle, ic, id, itmp, iweight, k, &
                                                            nsteps, unit_nr
      LOGICAL                                            :: should_output
      REAL(KIND=dp)                                      :: inv_norm, invproc
      REAL(KIND=dp), DIMENSION(3)                        :: center
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: cubdata
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      print_key => section_vals_get_subs_vals(helium_env(1)%helium%input, &
                                              "MOTION%PINT%HELIUM%PRINT%RHO")
      should_output = BTEST(cp_print_key_should_output( &
                            iteration_info=logger%iter_info, &
                            basis_section=print_key), cp_p_file)

      IF (.NOT. should_output) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! work on temporary array so that the average remains intact
      ! save accumulated data of different env on same core in first temp
      helium_env(1)%helium%rho_inst(:, :, :, :) = 0.0_dp
      DO k = 1, SIZE(helium_env)
         helium_env(1)%helium%rho_inst(:, :, :, :) = helium_env(1)%helium%rho_inst(:, :, :, :)+ &
                                                     helium_env(k)%helium%rho_accu(:, :, :, :)
      END DO

      ! average over processors
      CALL mp_sum(helium_env(1)%helium%rho_inst, helium_env(1)%comm)
      itmp = helium_env(1)%helium%num_env
      invproc = 1.0_dp/REAL(itmp, dp)
      helium_env(1)%helium%rho_inst(:, :, :, :) = helium_env(1)%helium%rho_inst(:, :, :, :)*invproc

      ! average over steps performed so far in this run
      nsteps = helium_env(1)%helium%current_step-helium_env(1)%helium%first_step
      inv_norm = 1.0_dp/REAL(nsteps, dp)
      helium_env(1)%helium%rho_inst(:, :, :, :) = helium_env(1)%helium%rho_inst(:, :, :, :)*inv_norm

      iweight = helium_env(1)%helium%rho_iweight
      ! average over the old and the current density (observe the weights!)
      helium_env(1)%helium%rho_inst(:, :, :, :) = nsteps*helium_env(1)%helium%rho_inst(:, :, :, :)+ &
                                                  iweight*helium_env(1)%helium%rho_rstr(:, :, :, :)
      helium_env(1)%helium%rho_inst(:, :, :, :) = helium_env(1)%helium%rho_inst(:, :, :, :)/REAL(nsteps+iweight, dp)

      ! set center of the cubefile
      IF (helium_env(1)%helium%solute_present) THEN
         ! should be set to solute's COM
         center(:) = helium_env(1)%helium%center(:)
      ELSE
         ! regardless of whether we are periodic or not use the origin, since
         ! pure cluster's COM might drift, but we want the cube fixed (note that
         ! the densities are correctly calculated wrt to COM in such case)
         center(:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
      END IF

      DO id = 1, rho_num ! loop over densities ---

         IF (.NOT. helium_env(1)%helium%rho_property(id)%is_calculated) THEN
            ! skip densities that are not requested by the user
            CYCLE
         END IF

         DO ic = 1, helium_env(1)%helium%rho_property(id)%num_components ! loop over components

            WRITE (fname, '(A)') "helium-rho-"// &
               TRIM(ADJUSTL(helium_env(1)%helium%rho_property(id)%filename_suffix(ic)))
            IF (helium_env(1)%helium%rho_property(id)%component_name(ic) .EQ. "") THEN
               WRITE (comment, '(A)') TRIM(helium_env(1)%helium%rho_property(id)%name)
            ELSE
               WRITE (comment, '(A)') TRIM(helium_env(1)%helium%rho_property(id)%name)// &
                  ", "// &
                  TRIM(helium_env(1)%helium%rho_property(id)%component_name(ic))
            END IF
            cubdata => helium_env(1)%helium%rho_inst(helium_env(1)%helium%rho_property(id)%component_index(ic), :, :, :)

            unit_nr = cp_print_key_unit_nr( &
                      logger, &
                      print_key, &
                      middle_name=TRIM(ADJUSTL(fname)), &
                      extension=".cube", &
                      file_position="REWIND", &
                      do_backup=.FALSE.)

            IF (logger%para_env%ionode) THEN
               CALL helium_write_cubefile( &
                  unit_nr, &
                  comment, &
                  center-0.5_dp*(helium_env(1)%helium%rho_maxr-helium_env(1)%helium%rho_delr), &
                  helium_env(1)%helium%rho_delr, &
                  helium_env(1)%helium%rho_nbin, &
                  cubdata)
               CALL m_flush(unit_nr)
               CALL cp_print_key_finished_output(unit_nr, logger, print_key)
            END IF

         END DO ! loop over components

      END DO ! loop over densities ---

      CALL timestop(handle)

   END SUBROUTINE helium_print_rho

! ***************************************************************************
!> \brief Write volumetric data to an orthorhombic cubefile
!> \param   unit - unit number to which output will be sent
!> \param   comment - description of the data stored in the cubefile
!> \param   origin - position of the cubefile origin
!> \param   deltar - voxel side length
!> \param   ndim - number of voxels in each dimension
!> \param   DATA - array (ndim x ndim x ndim) with the data for output
!> \date 2013-11-25
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_write_cubefile(unit, comment, origin, deltar, ndim, DATA)

      INTEGER, INTENT(IN)                                :: unit
      CHARACTER(len=default_string_length), INTENT(IN)   :: comment
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: origin
      REAL(KIND=dp), INTENT(IN)                          :: deltar
      INTEGER, INTENT(IN)                                :: ndim
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         POINTER                                         :: DATA

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_cubefile', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ix, iy, iz, nw
      REAL(KIND=dp)                                      :: delr, inva3
      REAL(KIND=dp), DIMENSION(3)                        :: orig

      CALL timeset(routineN, handle)

      ! convert cubefile data to the proper units of measure
      delr = angstrom*deltar
      orig(:) = angstrom*origin(:)

      ! output cube file header
      WRITE (unit, '(A)') comment
      WRITE (unit, '(A)') "Volumetric data in cubefile format generated by CP2K"
      WRITE (unit, '(I5,3(1X,F12.8))') 0, orig(1), orig(2), orig(3)
      WRITE (unit, '(I5,3(1X,F12.8))') ndim, delr, 0.0_dp, 0.0_dp
      WRITE (unit, '(I5,3(1X,F12.8))') ndim, 0.0_dp, delr, 0.0_dp
      WRITE (unit, '(I5,3(1X,F12.8))') ndim, 0.0_dp, 0.0_dp, delr

      ! output cube file data
      nw = 0
      inva3 = 1.0_dp/(angstrom*angstrom*angstrom)
      DO ix = 1, ndim
         DO iy = 1, ndim
            DO iz = 1, ndim
               WRITE (unit, '(1X,E13.5)', ADVANCE='NO') inva3*DATA(ix, iy, iz)
               nw = nw+1
               IF (MOD(nw, 6) .EQ. 0) THEN
                  nw = 0
                  WRITE (unit, *)
               END IF
            END DO
         END DO
      END DO
      ! some compilers write over the first entry on a line loosing all previous
      ! values written on that line unless line terminator is written at the end
      ! so make sure that the last WRITE statement runs without ADVANCE='NO'
      ! (observed for ifort 14.0.1 and 14.0.2 but not for gfortran 4.8.2)
      IF (nw .NE. 0) THEN
         WRITE (unit, *)
      END IF

      CALL timestop(handle)

   END SUBROUTINE helium_write_cubefile

! ***************************************************************************
!> \brief  Print permutation length according to HELIUM%PRINT%PLENGTH
!> \param helium_env ...
!> \date   2010-06-07
!> \par    History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_print_plength(helium_env)

      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_print_plength', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, unit_nr
      LOGICAL                                            :: is_new, should_output
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      IF (logger%para_env%ionode) THEN
         print_key => section_vals_get_subs_vals(helium_env(1)%helium%input, &
                                                 "MOTION%PINT%HELIUM%PRINT%PLENGTH")
         should_output = BTEST(cp_print_key_should_output( &
                               iteration_info=logger%iter_info, &
                               basis_section=print_key), cp_p_file)

         IF (should_output) THEN

            unit_nr = cp_print_key_unit_nr( &
                      logger, &
                      print_key, &
                      middle_name="helium-plength", &
                      extension=".dat", &
                      is_new_file=is_new)

            DO i = 1, helium_env(1)%helium%atoms
               WRITE (unit_nr, '(F20.10)', ADVANCE='NO') helium_env(1)%helium%plength_avrg(i)
               IF (i .LT. helium_env(1)%helium%atoms) THEN
                  WRITE (unit_nr, '(1X)', ADVANCE='NO')
               END IF
            END DO
            WRITE (unit_nr, '(A)') ""

            CALL m_flush(unit_nr)
            CALL cp_print_key_finished_output(unit_nr, logger, print_key)

         END IF
      END IF

      CALL timestop(handle)

      RETURN
   END SUBROUTINE helium_print_plength

! ***************************************************************************
!> \brief  Print helium force according to HELIUM%PRINT%FORCE
!> \param helium_env ...
!> \date   2010-01-27
!> \par    History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_print_force(helium_env)

      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_print_force', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: msgstr
      INTEGER                                            :: handle, ia, ib, ic, idim, unit_nr
      LOGICAL                                            :: should_output
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      print_key => section_vals_get_subs_vals(helium_env(1)%helium%input, &
                                              "MOTION%PINT%HELIUM%PRINT%FORCES")
      should_output = BTEST(cp_print_key_should_output( &
                            logger%iter_info, &
                            basis_section=print_key), cp_p_file)

      IF (.NOT. should_output) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! check if there is anything to be printed out
      IF (.NOT. helium_env(1)%helium%solute_present) THEN
         msgstr = "Warning: force printout requested but there is no solute!"
         CALL helium_write_line(msgstr)
         CALL timestop(handle)
         RETURN
      END IF

      IF (logger%para_env%ionode) THEN

         unit_nr = cp_print_key_unit_nr( &
                   logger, &
                   print_key, &
                   middle_name="helium-force", &
                   extension=".dat")

         ! print all force components in one line
         DO ib = 1, helium_env(1)%helium%solute_beads
            idim = 0
            DO ia = 1, helium_env(1)%helium%solute_atoms
               DO ic = 1, 3
                  idim = idim+1
                  WRITE (unit_nr, '(F20.10)', ADVANCE='NO') helium_env(1)%helium%force_avrg(ib, idim)
               END DO
            END DO
         END DO
         WRITE (unit_nr, *)

         CALL m_flush(unit_nr)
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)

      END IF

      CALL timestop(handle)

      RETURN
   END SUBROUTINE helium_print_force

#if 0

! ***************************************************************************
!> \brief  Print instantaneous force according to HELIUM%PRINT%FORCES_INST
!> \param helium ...
!> \date   2010-01-29
!> \author Lukasz Walewski
!> \note   Collects instantaneous helium forces from all processors on
!> logger%para_env%source and writes them to files - one file per processor.
!>         !TODO: Generalize for helium_env
! **************************************************************************************************
   SUBROUTINE helium_print_force_inst(helium)

      TYPE(helium_solvent_type), POINTER                 :: helium

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_print_force_inst', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: my_middle_name, stmp
      INTEGER                                            :: handle, ia, ib, ic, idim, irank, offset, &
                                                            unit_nr
      LOGICAL                                            :: should_output
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()
      print_key => section_vals_get_subs_vals(helium%input, &
                                              "MOTION%PINT%HELIUM%PRINT%FORCES_INST")
      should_output = BTEST(cp_print_key_should_output( &
                            logger%iter_info, &
                            basis_section=print_key), cp_p_file)

      IF (should_output) THEN

         ! check if there is anything to be printed out
         IF (.NOT. helium%solute_present) THEN
            stmp = "Warning: force printout requested but there is no solute!"
            CALL helium_write_line(stmp)
            CALL timestop(handle)
            RETURN
         END IF

         ! fill the tmp buffer with instantaneous helium forces at each proc
         helium%rtmp_p_ndim_1d(:) = PACK(helium%force_inst, .TRUE.)

         ! pass the message from all processors to logger%para_env%source
         helium%rtmp_p_ndim_np_1d(:) = 0.0_dp
         CALL mp_gather(helium%rtmp_p_ndim_1d, helium%rtmp_p_ndim_np_1d, &
                        logger%para_env%source, logger%para_env%group)

         IF (logger%para_env%ionode) THEN

            ! iterate over processors/helium environments
            DO irank = 1, helium%num_env

               ! generate one file per processor
               stmp = ""
               WRITE (stmp, *) irank
               my_middle_name = "helium-force-inst-"//TRIM(ADJUSTL(stmp))
               unit_nr = cp_print_key_unit_nr( &
                         logger, &
                         print_key, &
                         middle_name=TRIM(my_middle_name), &
                         extension=".dat")

               ! unpack and actually print the forces - all components in one line
               offset = (irank-1)*SIZE(helium%rtmp_p_ndim_1d)
               idim = 0
               DO ib = 1, helium%solute_beads
                  DO ia = 1, helium%solute_atoms
                     DO ic = 1, 3
                        idim = idim+1
                        WRITE (unit_nr, '(F20.10)', ADVANCE='NO') helium%rtmp_p_ndim_np_1d(offset+idim)
                     END DO
                  END DO
               END DO
               WRITE (unit_nr, *)

               CALL m_flush(unit_nr)
               CALL cp_print_key_finished_output(unit_nr, logger, print_key)

            END DO

         END IF
      END IF

      CALL timestop(handle)

      RETURN
   END SUBROUTINE helium_print_force_inst

#endif

END MODULE helium_io
