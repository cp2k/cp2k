!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for Molecular Dynamics
!> \author Teodoro Laino [tlaino] - University of Zurich - 09.2007
! **************************************************************************************************
MODULE md_util

   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_generate_filename
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE input_cp2k_restarts,             ONLY: write_restart
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE md_energies,                     ONLY: md_write_output
   USE md_environment_types,            ONLY: md_environment_type
   USE message_passing,                 ONLY: mp_bcast
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_util'

   PUBLIC :: md_output, &
             read_vib_eigs_unformatted

CONTAINS

! **************************************************************************************************
!> \brief collects the part of the MD that, basically, does the output
!> \param md_env ...
!> \param md_section ...
!> \param root_section ...
!> \param forced_io ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE md_output(md_env, md_section, root_section, forced_io)
      TYPE(md_environment_type), POINTER                 :: md_env
      TYPE(section_vals_type), POINTER                   :: md_section, root_section
      LOGICAL, INTENT(IN)                                :: forced_io

      CHARACTER(LEN=*), PARAMETER :: routineN = 'md_output', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: do_print
      TYPE(section_vals_type), POINTER                   :: print_section

      CALL timeset(routineN, handle)
      do_print = .TRUE.
      IF (forced_io) THEN
         print_section => section_vals_get_subs_vals(md_section, "PRINT")
         CALL section_vals_val_get(print_section, "FORCE_LAST", l_val=do_print)
      END IF
      IF (do_print) THEN
         ! Dumps all files related to the MD run
         CALL md_write_output(md_env)
         CALL write_restart(md_env=md_env, root_section=root_section)
      END IF
      CALL timestop(handle)

   END SUBROUTINE md_output

! **************************************************************************************************
!> \brief read eigenvalues and eigenvectors of Hessian from vibrational analysis results, for use
!>        of initialising MD simulations. Expects to read an unformatted binary file
!> \param md_section : input section object containing MD subsections and keywords. This should
!>                     provide the filename to read vib analysis eigenvalues and eigenvectors.
!>                     If the filename is not explicitly specified by the user in the input, then
!>                     it will use the default CARTESIAN_EIGS print key filename defined in the
!>                     vibrational analysis input section as the filename.
!> \param vib_section : input section object containing vibrational analysis subsections
!>                      and keywords
!> \param para_env : cp2k mpi environment object, needed for IO in parallel computations
!> \param dof : outputs the total number of eigenvalues (no. degrees of freedom) read from the file
!> \param eigenvalues : outputs the eigenvalues (Cartesian frequencies) read from the file
!> \param eigenvectors : outputs the corresponding eigenvectors read from the file
!> \author Lianheng Tong, lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   SUBROUTINE read_vib_eigs_unformatted(md_section, &
                                        vib_section, &
                                        para_env, &
                                        dof, &
                                        eigenvalues, &
                                        eigenvectors)
      TYPE(section_vals_type), POINTER                   :: md_section, vib_section
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(OUT)                               :: dof
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: eigenvectors

      CHARACTER(len=*), PARAMETER :: routineN = 'read_vib_eigs_unformatted', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: jj, n_rep_val, unit_nr
      LOGICAL                                            :: exist
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      logger => cp_get_default_logger()
      dof = 0
      eigenvalues = 0.0_dp
      eigenvectors = 0.0_dp
      ! obtain file name
      CALL section_vals_val_get(md_section, "INITIAL_VIBRATION%VIB_EIGS_FILE_NAME", &
                                n_rep_val=n_rep_val)
      IF (n_rep_val > 0) THEN
         CALL section_vals_val_get(md_section, "VIB_EIGS_FILE_NAME", c_val=filename)
      ELSE
         print_key => section_vals_get_subs_vals(vib_section, "PRINT%CARTESIAN_EIGS")
         filename = cp_print_key_generate_filename(logger, print_key, extension="eig", &
                                                   my_local=.FALSE.)
      END IF
      ! read file
      IF (para_env%ionode) THEN
         INQUIRE (FILE=filename, exist=exist)
         IF (.NOT. exist) THEN
            CPABORT("File "//filename//" is not found.")
         END IF
         CALL open_file(file_name=filename, &
                        file_action="READ", &
                        file_form="UNFORMATTED", &
                        file_status="OLD", &
                        unit_number=unit_nr)
         ! the first record contains one integer giving degrees of freedom
         READ (unit_nr) dof
         IF (dof .GT. SIZE(eigenvalues)) THEN
            CPABORT("Too many DoFs found in "//filename)
         END IF
         ! the second record contains the eigenvalues
         READ (unit_nr) eigenvalues(1:dof)
         ! the rest of the records contain the eigenvectors
         DO jj = 1, dof
            READ (unit_nr) eigenvectors(1:dof, jj)
         END DO
      END IF
      ! broadcast to all compulational nodes. note that it is assumed
      ! that source is the ionode
      CALL mp_bcast(dof, para_env%source, para_env%group)
      CALL mp_bcast(eigenvalues, para_env%source, para_env%group)
      CALL mp_bcast(eigenvectors, para_env%source, para_env%group)
      ! close file
      IF (para_env%ionode) THEN
         CALL close_file(unit_number=unit_nr)
      END IF
   END SUBROUTINE read_vib_eigs_unformatted

END MODULE md_util
