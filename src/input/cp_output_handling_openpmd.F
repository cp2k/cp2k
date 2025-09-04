!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines to handle the output, The idea is to remove the
!>      decision of wheter to output and what to output from the code
!>      that does the output, and centralize it here.
!> \note
!>      These were originally together with the log handling routines,
!>      but have been spawned off. Some dependencies are still there,
!>      and some of the comments about log handling also applies to output
!>      handling: @see cp_log_handling
!> \par History
!>      12.2001 created [fawzi]
!>      08.2002 updated to new logger [fawzi]
!>      10.2004 big rewrite of the output methods, connected to the new
!>              input, and iteration_info [fawzi]
!>      08.2005 property flags [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE cp_output_handling_openpmd
   USE cp_output_handling, only: cp_print_key_should_output, cp_p_file
   USE cp_files, ONLY: close_file, &
                       open_file
   USE cp_iter_types, ONLY: cp_iteration_info_release, &
                            cp_iteration_info_retain, &
                            cp_iteration_info_type, &
                            each_desc_labels, &
                            each_possible_labels
   USE cp_log_handling, ONLY: cp_logger_generate_filename, &
                              cp_logger_get_default_unit_nr, &
                              cp_logger_get_unit_nr, &
                              cp_logger_type, &
                              cp_to_string
   USE input_keyword_types, ONLY: keyword_create, &
                                  keyword_release, &
                                  keyword_type
   USE input_section_types, ONLY: section_add_keyword, &
                                  section_add_subsection, &
                                  section_create, &
                                  section_release, &
                                  section_type, &
                                  section_vals_get_subs_vals, &
                                  section_vals_type, &
                                  section_vals_val_get
   USE kinds, ONLY: default_path_length, &
                    default_string_length
   USE machine, ONLY: m_mov
   USE memory_utilities, ONLY: reallocate
   USE message_passing, ONLY: mp_file_delete, &
                              mp_file_get_amode, &
                              mp_file_type
   USE openpmd_api, ONLY: &
      openpmd_access_create, openpmd_attributable_set_attribute_vec_int, &
      openpmd_iteration_get_mesh, &
      openpmd_iteration_get_particle_species, openpmd_iteration_type, openpmd_mesh_type, &
      openpmd_particle_species_get_record, openpmd_particle_species_type, &
      openpmd_record_component_make_empty, openpmd_record_type, &
      openpmd_series_create, openpmd_series_present, openpmd_series_type, &
      openpmd_series_write_iteration, openpmd_type_int, openpmd_json_merge
   USE string_utilities, ONLY: compress, &
                               s2a
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_output_handling_openpmd'
   PUBLIC :: cp_openpmd_print_key_unit_nr, cp_openpmd_print_key_finished_output
   PUBLIC :: cp_openpmd_get_value_unit_nr, cp_openpmd_per_call_value_type
   PUBLIC :: cp_openpmd_output_finalize

   TYPE :: cp_openpmd_per_call_value_type
      TYPE(openpmd_series_type) :: series
      TYPE(openpmd_iteration_type) :: iteration
      ! TYPE(openpmd_mesh_type) :: mesh
      ! TYPE(openpmd_particle_species_type) :: particle_species
      CHARACTER(len=default_string_length) :: name_prefix ! e.g. 'WFN_00008_1'
   END TYPE

   TYPE :: cp_openpmd_per_call_type
      INTEGER :: key ! unit_nr
      TYPE(cp_openpmd_per_call_value_type) :: value
   END TYPE cp_openpmd_per_call_type

   TYPE :: cp_current_iteration_counter_type
      INTEGER :: flat_iteration = 0
      INTEGER, ALLOCATABLE :: complex_iteration(:)
      INTEGER :: complex_iteration_depth = 0
   END TYPE cp_current_iteration_counter_type

   TYPE :: cp_openpmd_per_callsite_value_type
      ! openPMD output Series.
      TYPE(openpmd_series_type) :: output_series
      ! Information on the last Iteration that was written to, including
      ! CP2Ks complex Iteration number and its associated contiguous scalar
      ! openPMD Iteration number.
      TYPE(cp_current_iteration_counter_type) :: iteration_counter
   END TYPE cp_openpmd_per_callsite_value_type

   TYPE :: cp_openpmd_per_callsite_type
      CHARACTER(len=default_string_length) :: key ! openpmd_basename
      TYPE(cp_openpmd_per_callsite_value_type) :: value
   END TYPE cp_openpmd_per_callsite_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin data members for openPMD output. !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Map that associates opened unit numbers with their associated openPMD content.
   ! Since CP2K logically opens a new file for every single dataset, multiple
   ! unit numbers may point to the same openPMD Series.
   TYPE(cp_openpmd_per_call_type), ALLOCATABLE :: cp_openpmd_per_call(:)
   INTEGER :: cp_num_openpmd_per_call = 0
   INTEGER :: cp_capacity_openpmd_per_call = 0

   ! Map that associates callsites from which functions of this module may be invoked
   ! to their associated openPMD content.
   ! This stores the actual output Series (which stays open across calls from the
   ! same callsite) and the Iteration counter (which associates complex CP2k
   ! Iterations with flattened scalar Iteration indexes in the openPMD output).
   TYPE(cp_openpmd_per_callsite_type), ALLOCATABLE :: cp_openpmd_per_callsite(:)
   INTEGER :: cp_num_openpmd_per_callsite = 0
   INTEGER :: cp_capacity_openpmd_per_callsite = 0

   INTEGER, PARAMETER :: cp_allocation_size = 100
   CHARACTER(len=*), PARAMETER :: cp_default_backend_config = &
                                  "[hdf5]"//new_line('a')// &
                                  "# will be overridden by particle flushes"//new_line('a')// &
                                  "independent_stores = false"//new_line('a')// &
                                  "dont_warn_unused_keys = ['independent_stores']"//new_line('a')// &
                                  ""//new_line('a')// &
                                  "[adios2]"//new_line('a')// &
                                  "# discard any attributes written on ranks other than 0"//new_line('a')// &
                                  "attribute_writing_ranks = 0"//new_line('a')// &
                                  "[adios2.engine]"//new_line('a')// &
                                  "# CP2K generally has many small IO operations, "//new_line('a')// &
                                  "# so stage IO memory to the buffer first and then "//new_line('a')// &
                                  "# run it all at once, instead of writing to disk directly."//new_line('a')// &
                                  "# Save memory by specifying 'disk' here instead."//new_line('a')// &
                                  "# TODO: In future, maybe implement some input variable"//new_line('a')// &
                                  "#       to specify intervals at which to flush to disk."//new_line('a')// &
                                  "preferred_flush_target = 'buffer'"//new_line('a')// &
#ifndef _WIN32
                                  "# Raise the BufferChunkSize to the maximum (2GB), since large operations"//new_line('a')// &
                                  "# improve IO performance and the allocation overhead only cuts into"//new_line('a')// &
                                  "# virtual memory (except on Windows, hence do not do that there)"//new_line('a')// &
                                  "[adios2.engine.parameters]"//new_line('a')// &
                                  "BufferChunkSize = 2147381248"//new_line('a')
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! End data members for openPMD output. !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

   #:set name_suffixes = ['unit_nr', 'filedata']
   #:set key_types = ['INTEGER', 'CHARACTER(len=default_string_length)']
   #:set value_types = ['cp_openpmd_per_call_value_type', 'cp_openpmd_per_callsite_value_type']
   #:set map_storages = ['cp_openpmd_per_call', 'cp_openpmd_per_callsite']
   #:set map_counters = ['cp_num_openpmd_per_call', 'cp_num_openpmd_per_callsite']
   #:set map_capacities = ['cp_capacity_openpmd_per_call', 'cp_capacity_openpmd_per_callsite']

   #:for name_suffix, key_type, value_type, map_storage, map_counter, map_capacity in zip(name_suffixes, key_types, value_types, map_storages, map_counters, map_capacities)

      ! TODO No reallocation support for now, change cp_allocation_size if larger sizes needed.

      ! **************************************************************************************************
      !> \brief ...
      !> \param key ...
      !> \param value ...
      ! **************************************************************************************************
      FUNCTION cp_openpmd_add_${name_suffix}$ (key, value) RESULT(index)
         ${key_type}$, INTENT(in)                           :: key
         TYPE(${value_type}$), INTENT(in)                   :: value
         INTEGER                                            :: index

         INTEGER                                            :: i

         ! Check if the key already exists
         DO i = 1, ${map_counter}$
            IF (${map_storage}$ (i)%key == key) THEN
               ${map_storage}$ (i)%value = value
               index = i
               RETURN
            END IF
         END DO

         IF (${map_capacity}$ == 0) THEN
            ALLOCATE (${map_storage}$ (cp_allocation_size))
            ${map_capacity}$ = cp_allocation_size
         END IF

         ! No idea how to do reallocations, so for now just assert that they're not needed
         CPASSERT(${map_counter}$ < ${map_capacity}$)

         ! Add a new entry
         ${map_counter}$ = ${map_counter}$+1
         ${map_storage}$ (${map_counter}$)%key = key
         ${map_storage}$ (${map_counter}$)%value = value
         index = ${map_counter}$
      END FUNCTION cp_openpmd_add_${name_suffix}$

      ! **************************************************************************************************
      !> \brief ...
      !> \param key ...
      !> \return ...
      ! **************************************************************************************************
      FUNCTION cp_openpmd_get_index_${name_suffix}$ (key) RESULT(index)
         ${key_type}$, INTENT(in)                           :: key
         INTEGER                                            :: index

         INTEGER                                            :: i

         index = -1

         DO i = 1, ${map_counter}$
            IF (${map_storage}$ (i)%key == key) THEN
               index = i
               RETURN
            END IF
         END DO
      END FUNCTION cp_openpmd_get_index_${name_suffix}$

      FUNCTION cp_openpmd_get_value_${name_suffix}$ (key) RESULT(value)
         ${key_type}$, INTENT(in)                           :: key
         TYPE(${value_type}$)                               :: value

         INTEGER                                            :: i

         i = cp_openpmd_get_index_${name_suffix}$ (key)
         IF (i == -1) RETURN

         value = ${map_storage}$ (i)%value
      END FUNCTION cp_openpmd_get_value_${name_suffix}$

      ! **************************************************************************************************
      !> \brief ...
      !> \param key ...
      !> \return ...
      ! **************************************************************************************************
      FUNCTION cp_openpmd_remove_${name_suffix}$ (key) RESULT(was_found)
         ${key_type}$, INTENT(in)                           :: key
         LOGICAL                                            :: was_found

         INTEGER                                            :: i

         was_found = .FALSE.

         DO i = 1, ${map_counter}$
            IF (${map_storage}$ (i)%key == key) THEN
               was_found = .TRUE.
               IF (i /= ${map_counter}$) THEN
                  ! Swap last element to now freed place
                  ${map_storage}$ (i) = ${map_storage}$ (${map_counter}$)
               END IF

               ${map_counter}$ = ${map_counter}$-1
               IF (${map_counter}$ == 0) THEN
                  DEALLOCATE (${map_storage}$)
                  ${map_capacity}$ = 0
               END IF
               RETURN
            END IF
         END DO
      END FUNCTION cp_openpmd_remove_${name_suffix}$

   #:endfor

! **************************************************************************************************
!> \brief Simplified version of cp_print_key_generate_filename. Since an openPMD Series encompasses
!         multiple datasets that would be separate outputs in e.g. .cube files, this needs not
!         consider dataset names for creation of a filename.
!> \param logger ...
!> \param print_key ...
!> \param openpmd_basename ...
!> \param extension ...
!> \param my_local ...
!> \return ...
! **************************************************************************************************
   FUNCTION cp_print_key_generate_openpmd_filename(logger, print_key, openpmd_basename, extension) RESULT(filename)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key
      CHARACTER(len=*), INTENT(IN)                       :: openpmd_basename, extension
      CHARACTER(len=default_path_length)                 :: filename

      CHARACTER(len=default_path_length)                 :: outPath, root
      CHARACTER(len=default_string_length)               :: outName
      INTEGER                                            :: my_ind1, my_ind2
      LOGICAL                                            :: has_root

      CALL section_vals_val_get(print_key, "FILENAME", c_val=outPath)
      IF (outPath(1:1) == '=') THEN
         CPASSERT(LEN(outPath) - 1 <= LEN(filename))
         filename = outPath(2:)
         RETURN
      END IF
      IF (outPath == "__STD_OUT__") outPath = ""
      outName = outPath
      has_root = .FALSE.
      my_ind1 = INDEX(outPath, "/")
      my_ind2 = LEN_TRIM(outPath)
      IF (my_ind1 /= 0) THEN
         has_root = .TRUE.
         DO WHILE (INDEX(outPath(my_ind1 + 1:my_ind2), "/") /= 0)
            my_ind1 = INDEX(outPath(my_ind1 + 1:my_ind2), "/") + my_ind1
         END DO
         IF (my_ind1 == my_ind2) THEN
            outName = ""
         ELSE
            outName = outPath(my_ind1 + 1:my_ind2)
         END IF
      END IF

      IF (.NOT. has_root) THEN
         root = TRIM(logger%iter_info%project_name)
      ELSE IF (outName == "") THEN
         root = outPath(1:my_ind1)//TRIM(logger%iter_info%project_name)
      ELSE
         root = outPath(1:my_ind1)
      END IF

      filename = ADJUSTL(TRIM(root)//"_"//TRIM(openpmd_basename)//TRIM(extension))

   END FUNCTION cp_print_key_generate_openpmd_filename

! **************************************************************************************************
!> \brief CP2K Iteration numbers are n-dimensional while openPMD Iteration numbers are scalars.
!         This checks if the Iteration number has changed from the previous call (stored in
!         openpmd_file%iteration_counter) and updates it if needed.
!> \param logger ...
!> \return ...
! **************************************************************************************************
   FUNCTION cp_advance_iteration_number(logger, openpmd_file) RESULT(did_advance_iteration)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_openpmd_per_callsite_value_type)           :: openpmd_file
      LOGICAL                                            :: did_advance_iteration

      INTEGER                                            :: len

      #:set ic = 'openpmd_file%iteration_counter'

      did_advance_iteration = .FALSE.
      len = SIZE(logger%iter_info%iteration)
      IF (len /= ${ic}$%complex_iteration_depth) THEN
         did_advance_iteration = .TRUE.
         ${ic}$%complex_iteration_depth = len
         ALLOCATE (${ic}$%complex_iteration(len))
      ELSE
         did_advance_iteration &
            = ANY(${ic}$%complex_iteration(1:len) &
                  /= logger%iter_info%iteration(1:len))
      END IF

      IF (.NOT. did_advance_iteration) RETURN

      ${ic}$%flat_iteration = ${ic}$%flat_iteration + 1
      ${ic}$%complex_iteration(1:len) &
         = logger%iter_info%iteration(1:len)

   END FUNCTION cp_advance_iteration_number

! **************************************************************************************************
!> \brief CP2K deals with output handles in terms of unit numbers.
!         The openPMD output logic does not change this association.
!         For this, we need to emulate unit numbers as (1) they are not native to openPMD and
!         (2) a single openPMD Series might contain multiple datasets treated logically by CP2K
!         as distinct outputs. As a result, a single unit number is resolved by the openPMD logic
!         to the values represented by the cp_print_key_generate_openpmd_filename struct,
!         containing the output Series and the referred datasets therein (Iteration number, mesh
!         name, particle species name, ...).
!> \param series ...
!> \param middle_name ...
!> \param logger ...
!> \return ...
! **************************************************************************************************
   FUNCTION cp_openpmd_create_unit_nr_entry(openpmd_file_index, middle_name, logger) RESULT(res)
      INTEGER                                            :: openpmd_file_index
      CHARACTER(len=*), INTENT(IN)                       :: middle_name
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_openpmd_per_call_value_type)               :: res

      LOGICAL                                            :: opened_new_iteration = .FALSE.

      #:set opmd = 'cp_openpmd_per_callsite(openpmd_file_index)%value'

      res%series = ${opmd}$%output_series

      opened_new_iteration = cp_advance_iteration_number(logger, ${opmd}$)

      res%iteration = openpmd_series_write_iteration( &
                      ${opmd}$%output_series, ${opmd}$%iteration_counter%flat_iteration)

      ! res%mesh = openpmd_iteration_get_mesh(res%iteration, middle_name)
      ! CALL openpmd_record_component_make_empty(res%mesh%openpmd_mesh_as_record_component(), openpmd_type_int, 3)
      ! res%particle_species = openpmd_iteration_get_particle_species(res%iteration, middle_name)
      ! position = openpmd_particle_species_get_record(res%particle_species, "position")
      ! CALL openpmd_record_component_make_empty_int(position%openpmd_record_as_record_component(), 3)
      ! positionOffset = openpmd_particle_species_get_record(res%particle_species, "positionOffset")
      ! CALL openpmd_record_component_make_empty_int(positionOffset%openpmd_record_as_record_component(), 3)
      res%name_prefix = TRIM(middle_name)

      IF (opened_new_iteration) THEN
         CALL openpmd_attributable_set_attribute_vec_int( &
            res%iteration%openpmd_iteration_as_attributable(), &
            "ndim_iteration_index", &
            ${opmd}$%iteration_counter%complex_iteration)
      END IF
   END FUNCTION

! **************************************************************************************************
!> \brief Check if there is already an output Series created for the callsite identified
!         by openpmd_basename. If so, then return it (by index), otherwise open the Series now
!         and return the index then.
   FUNCTION cp_openpmd_get_openpmd_file_entry(openpmd_basename, filename, openpmd_config, logger, use_mpi) RESULT(file_index)
      CHARACTER(len=*), INTENT(IN)             :: openpmd_basename, filename, openpmd_config
      TYPE(cp_logger_type), POINTER            :: logger
      LOGICAL                                  :: use_mpi
      INTEGER                                  :: file_index
      CHARACTER(:), ALLOCATABLE                :: merged_config

      CHARACTER(len=default_string_length)     :: basename_copied = ' '
      TYPE(cp_openpmd_per_callsite_value_type) :: emplace_new

      basename_copied(1:LEN_TRIM(openpmd_basename)) = TRIM(openpmd_basename)

      file_index = cp_openpmd_get_index_filedata(basename_copied)
      IF (file_index /= -1) RETURN

      IF (use_mpi) THEN
         merged_config = openpmd_json_merge(cp_default_backend_config, openpmd_config, logger%para_env)
         emplace_new%output_series = openpmd_series_create( &
                                     filename, openpmd_access_create, logger%para_env, merged_config)
      ELSE
         merged_config = openpmd_json_merge(cp_default_backend_config, openpmd_config)
         emplace_new%output_series = openpmd_series_create( &
                                     filename, openpmd_access_create, config=merged_config)
      END IF
      DEALLOCATE (merged_config)
      file_index = cp_openpmd_add_filedata(basename_copied, emplace_new)
   END FUNCTION

! **************************************************************************************************
!> \brief Close all outputs.
! **************************************************************************************************
   SUBROUTINE cp_openpmd_output_finalize()
      INTEGER :: i
      DO i = 1, cp_num_openpmd_per_callsite
         DEALLOCATE (cp_openpmd_per_callsite(i)%value%iteration_counter%complex_iteration)
         CALL cp_openpmd_per_callsite(i)%value%output_series%close()
      END DO
      IF (ALLOCATED(cp_openpmd_per_callsite)) THEN
         DEALLOCATE (cp_openpmd_per_callsite)
      END IF
      cp_num_openpmd_per_callsite = 0
   END SUBROUTINE cp_openpmd_output_finalize

! **************************************************************************************************
!> \brief ...
!> \param logger ...
!> \param basis_section ...
!> \param print_key_path ...
!> \param extension ...
!> \param middle_name ...
!> \param local ...
!> \param log_filename ...
!> \param ignore_should_output ...
!> \param do_backup ...
!> \param is_new_file true if this rank created a new (or rewound) file, false otherwise
!> \param mpi_io True if the file should be opened in parallel on all processors belonging to
!>               the communicator group. Automatically disabled if the file form or access mode
!>               is unsuitable for MPI IO. Return value indicates whether MPI was actually used
!>               and therefore the flag must also be passed to the file closing directive.
!> \param fout   Name of the actual file where the output will be written. Needed mainly for MPI IO
!>               because inquiring the filename from the MPI filehandle does not work across
!>               all MPI libraries.
!> \param openpmd_basename Used to associate an identifier to each callsite of this module
!> \param use_openpmd ...
!> \return ...
! **************************************************************************************************
   FUNCTION cp_openpmd_print_key_unit_nr(logger, basis_section, print_key_path, &
                                         middle_name, log_filename, ignore_should_output, &
                                         mpi_io, &
                                         fout, openpmd_basename) RESULT(res)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), INTENT(IN)                :: basis_section
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: print_key_path
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: middle_name
      LOGICAL, INTENT(IN), OPTIONAL                      :: log_filename, ignore_should_output
      LOGICAL, INTENT(INOUT), OPTIONAL                   :: mpi_io
      CHARACTER(len=default_path_length), INTENT(OUT), &
         OPTIONAL                                        :: fout
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: openpmd_basename
      INTEGER                                            :: res

      CHARACTER(len=default_path_length)                 :: filename

      CHARACTER(len=default_string_length)               :: openpmd_config, outPath, file_extension
      LOGICAL                                            :: do_log, found, &
                                                            my_mpi_io, &
                                                            my_should_output, &
                                                            replace
      INTEGER                                            :: openpmd_file_index, openpmd_call_index
      TYPE(section_vals_type), POINTER                   :: print_key

      my_mpi_io = .FALSE.
      replace = .FALSE.
      found = .FALSE.
      res = -1
      IF (PRESENT(mpi_io)) THEN
#if defined(__parallel)
         IF (logger%para_env%num_pe > 1 .AND. mpi_io) THEN
            my_mpi_io = .TRUE.
         ELSE
            my_mpi_io = .FALSE.
         END IF
#else
         my_mpi_io = .FALSE.
#endif
         ! Set return value
         mpi_io = my_mpi_io
      END IF
      NULLIFY (print_key)
      CPASSERT(ASSOCIATED(logger))
      CPASSERT(basis_section%ref_count > 0)
      CPASSERT(logger%ref_count > 0)
      my_should_output = BTEST(cp_print_key_should_output(logger%iter_info, &
                                                          basis_section, print_key_path, used_print_key=print_key), cp_p_file)
      IF (PRESENT(ignore_should_output)) my_should_output = my_should_output .OR. ignore_should_output
      IF (.NOT. my_should_output) RETURN
      IF (logger%para_env%is_source() .OR. my_mpi_io) THEN

         CALL section_vals_val_get(print_key, "FILENAME", c_val=outPath)
         CALL section_vals_val_get(print_key, "OPENPMD_EXTENSION", c_val=file_extension)
         CALL section_vals_val_get(print_key, "OPENPMD_CFG_FILE", c_val=openpmd_config)
         IF (LEN_TRIM(openpmd_config) == 0) THEN
            CALL section_vals_val_get(print_key, "OPENPMD_CFG", c_val=openpmd_config)
         ELSE
            openpmd_config = "@"//openpmd_config
         END IF
         filename = cp_print_key_generate_openpmd_filename(logger, print_key, openpmd_basename, file_extension)

         ! Check is we have to log any operation performed on the file..
         do_log = .FALSE.
         IF (PRESENT(log_filename)) THEN
            do_log = log_filename
         ELSE
            CALL section_vals_val_get(print_key, "LOG_PRINT_KEY", l_val=do_log)
         END IF

         IF (PRESENT(fout)) THEN
            fout = filename
         END IF

         openpmd_file_index = cp_openpmd_get_openpmd_file_entry( &
                              openpmd_basename, filename, openpmd_config, logger, my_mpi_io)

         OPEN (newunit=res, status='scratch', action='write')
         openpmd_call_index = cp_openpmd_add_unit_nr( &
                              res, &
                              cp_openpmd_create_unit_nr_entry( &
                              openpmd_file_index, middle_name, logger))

      ELSE
         res = -1
      END IF
   END FUNCTION cp_openpmd_print_key_unit_nr

! **************************************************************************************************
!> \brief should be called after you finish working with a unit obtained with
!>      cp_openpmd_print_key_unit_nr, so that the file that might have been opened
!>      can be closed.
!>
!>      the inputs should be exactly the same of the corresponding
!>      cp_openpmd_print_key_unit_nr
!> \param unit_nr ...
!> \param logger ...
!> \param basis_section ...
!> \param print_key_path ...
!> \param local ...
!> \param ignore_should_output ...
!> \param mpi_io True if file was opened in parallel with MPI
!> \param use_openpmd ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      closes if the corresponding filename of the printkey is
!>      not __STD_OUT__
! **************************************************************************************************
   SUBROUTINE cp_openpmd_print_key_finished_output(unit_nr, logger, basis_section, &
                                                   print_key_path, local, ignore_should_output, &
                                                   mpi_io)
      INTEGER, INTENT(INOUT)                             :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), INTENT(IN)                :: basis_section
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: print_key_path
      LOGICAL, INTENT(IN), OPTIONAL                      :: local, ignore_should_output, &
                                                            mpi_io

      CHARACTER(len=default_string_length)               :: outPath
      LOGICAL                                            :: my_local, my_mpi_io, &
                                                            my_should_output
      TYPE(section_vals_type), POINTER                   :: print_key

      my_local = .FALSE.
      my_mpi_io = .FALSE.
      NULLIFY (print_key)
      IF (PRESENT(local)) my_local = local
      IF (PRESENT(mpi_io)) my_mpi_io = mpi_io
      CPASSERT(ASSOCIATED(logger))
      CPASSERT(basis_section%ref_count > 0)
      CPASSERT(logger%ref_count > 0)
      my_should_output = BTEST(cp_print_key_should_output(logger%iter_info, basis_section, &
                                                          print_key_path, used_print_key=print_key), cp_p_file)
      IF (PRESENT(ignore_should_output)) my_should_output = my_should_output .OR. ignore_should_output
      IF (my_should_output .AND. (my_local .OR. &
                                  logger%para_env%is_source() .OR. &
                                  my_mpi_io)) THEN
         CALL section_vals_val_get(print_key, "FILENAME", c_val=outPath)
         IF (cp_openpmd_remove_unit_nr(unit_nr)) THEN
            CLOSE (unit_nr)
         END IF

         unit_nr = -1
      END IF
      CPASSERT(unit_nr == -1)
      unit_nr = -1
   END SUBROUTINE cp_openpmd_print_key_finished_output

END MODULE cp_output_handling_openpmd
