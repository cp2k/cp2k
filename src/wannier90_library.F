!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief In-process Wannier90 localization with distributed k-point overlaps.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE wannier90_library
   USE bibliography, ONLY: cite_reference, Pizzi2020, Marzari1997, Souza2001
   USE cp_files, ONLY: close_file, open_file
   USE ieee_arithmetic, ONLY: ieee_is_finite
   USE input_section_types, ONLY: section_vals_type, section_vals_val_get
   USE kinds, ONLY: dp, default_string_length
   USE message_passing, ONLY: mp_para_env_type
   USE physcon, ONLY: angstrom, evolt
#if defined(__WANNIER90)
   USE w90_library, ONLY: lib_common_type, w90_get_centres, w90_get_gkpb, w90_get_nn, &
                          w90_get_nnkp, w90_get_spreads, w90_input_setopt, w90_project_overlap, &
                          w90_set_comm, w90_set_eigval, w90_set_m_local, w90_set_option, &
                          w90_set_u_matrix, w90_set_u_opt, w90_wannierise, w90_disentangle, w90_plot
#endif
#include "./base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'wannier90_library'
   CHARACTER(LEN=12), PARAMETER :: window_keywords(4) = &
                                   [CHARACTER(LEN=12) :: 'dis_win_min', 'dis_win_max', 'dis_froz_min', 'dis_froz_max']

   TYPE w90_library_type
      PRIVATE
#if defined(__WANNIER90)
      TYPE(lib_common_type), ALLOCATABLE :: data
#endif
      INTEGER :: output_unit = -1, num_iter = 0, conv_window = 0, num_cg_steps = 5, num_print_cycles = 10
      LOGICAL :: keep_output = .FALSE., inputs_consumed = .FALSE.
      LOGICAL :: write_hr = .FALSE., write_u_matrices = .FALSE., converged = .FALSE.
      REAL(KIND=dp) :: conv_tol = 0.0_dp
      INTEGER :: dis_num_iter = 1000, dis_conv_window = 5
      REAL(KIND=dp) :: dis_conv_tol = 1.0e-10_dp, dis_mix_ratio = 0.5_dp, energy_window(4) = 0.0_dp
      LOGICAL :: window_set(4) = .FALSE.
      INTEGER, ALLOCATABLE :: local_index(:)
      COMPLEX(KIND=dp), POINTER :: overlaps(:, :, :, :) => NULL()
      COMPLEX(KIND=dp), POINTER :: projections(:, :, :) => NULL(), rotations(:, :, :) => NULL()
      REAL(KIND=dp), POINTER :: eigenvalues(:, :) => NULL()
   END TYPE w90_library_type

   PUBLIC :: w90_library_type, w90_library_init, w90_library_overlap, w90_library_eigenvalues, &
             w90_library_run, w90_library_release, w90_library_write_options, w90_library_projections, &
             w90_library_copy_matrices, w90_library_promote_output, w90_write_amn, w90_library_write_input, &
             w90_library_project_trials, w90_library_write_results

CONTAINS

! **************************************************************************************************
!> \brief Initialize the library and obtain its finite-difference mesh.
!> \param context Library state and matrix storage.
!> \param input LIBRARY input section.
!> \param para_env Communicator shared by all export ranks.
!> \param seed_name Output file prefix.
!> \param lattice Real lattice vectors in columns, in angstrom.
!> \param kpoints Fractional full-mesh k-points.
!> \param mp_grid Mesh dimensions.
!> \param num_bands Number of exported bands.
!> \param num_wann Number of Wannier functions.
!> \param nnlist Neighbour k-point indices.
!> \param nncell Reciprocal-cell shifts for the neighbours.
!> \param log_suffix Optional trial suffix for log files, without changing the library seed name.
!> \param num_cg_steps Optional optimizer override for an independent trial.
! **************************************************************************************************
   SUBROUTINE w90_library_init(context, input, para_env, seed_name, lattice, kpoints, mp_grid, &
                               num_bands, num_wann, nnlist, nncell, log_suffix, num_cg_steps)
      TYPE(w90_library_type), INTENT(INOUT) :: context
      TYPE(section_vals_type), POINTER :: input
      TYPE(mp_para_env_type), POINTER :: para_env
      CHARACTER(LEN=*), INTENT(IN) :: seed_name
      REAL(KIND=dp), INTENT(IN) :: lattice(3, 3), kpoints(:, :)
      INTEGER, INTENT(IN) :: mp_grid(3), num_bands, num_wann
      INTEGER, ALLOCATABLE, INTENT(OUT) :: nnlist(:, :), nncell(:, :, :)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: log_suffix
      INTEGER, INTENT(IN), OPTIONAL :: num_cg_steps
#if defined(__WANNIER90)
      INTEGER :: ierr, ik, iwann, nk, nn, nlocal, iwindow
      INTEGER, ALLOCATABLE :: distk(:)
      CHARACTER(LEN=16) :: rank_suffix
      CHARACTER(LEN=default_string_length) :: log_seed

      CPASSERT(.NOT. ALLOCATED(context%data))
      CALL cite_reference(Pizzi2020)
      CALL cite_reference(Marzari1997)
      context%inputs_consumed = .FALSE.
      context%converged = .FALSE.
      ALLOCATE (context%data)
      IF (LEN_TRIM(seed_name) > LEN(context%data%seedname)) THEN
         CPABORT("SEED_NAME exceeds the filename prefix length supported by the Wannier90 library.")
      END IF
      nk = SIZE(kpoints, 2)
      CALL section_vals_val_get(input, "NUM_ITER", i_val=context%num_iter)
      CALL section_vals_val_get(input, "CONV_WINDOW", i_val=context%conv_window)
      CALL section_vals_val_get(input, "CONV_TOL", r_val=context%conv_tol)
      CALL section_vals_val_get(input, "NUM_CG_STEPS", i_val=context%num_cg_steps)
      CALL section_vals_val_get(input, "NUM_PRINT_CYCLES", i_val=context%num_print_cycles)
      CALL section_vals_val_get(input, 'WRITE_HR', l_val=context%write_hr)
      CALL section_vals_val_get(input, 'WRITE_U_MATRICES', l_val=context%write_u_matrices)
      IF (PRESENT(num_cg_steps)) context%num_cg_steps = num_cg_steps
      IF (context%num_cg_steps < 0) CPABORT("Wannier90 NUM_CG_STEPS must be nonnegative.")
      IF (context%num_print_cycles < 1) CPABORT("Wannier90 NUM_PRINT_CYCLES must be positive.")
      IF (context%num_iter < 1 .OR. context%conv_window < 2 .OR. context%conv_tol <= 0.0_dp) THEN
         CPABORT("WANNIER90 LIBRARY requires NUM_ITER > 0, CONV_WINDOW >= 2 and CONV_TOL > 0.")
      END IF
      IF (num_wann < 1 .OR. num_wann > num_bands) THEN
         CPABORT("WANNIER90 LIBRARY requires 1 <= WANNIER_FUNCTIONS <= number of exported bands.")
      END IF
      CALL section_vals_val_get(input, 'DIS_NUM_ITER', i_val=context%dis_num_iter)
      CALL section_vals_val_get(input, 'DIS_CONV_WINDOW', i_val=context%dis_conv_window)
      CALL section_vals_val_get(input, 'DIS_CONV_TOL', r_val=context%dis_conv_tol)
      CALL section_vals_val_get(input, 'DIS_MIX_RATIO', r_val=context%dis_mix_ratio)
      IF (context%dis_num_iter < 1 .OR. context%dis_conv_window < 2 .OR. context%dis_conv_tol <= 0.0_dp) THEN
         CPABORT('Invalid Wannier90 disentanglement convergence settings.')
      END IF
      IF (context%dis_mix_ratio <= 0.0_dp .OR. context%dis_mix_ratio > 1.0_dp) THEN
         CPABORT('Wannier90 DIS_MIX_RATIO must be in (0, 1].')
      END IF
      DO iwindow = 1, SIZE(window_keywords)
         CALL section_vals_val_get(input, TRIM(window_keywords(iwindow)), explicit=context%window_set(iwindow))
         IF (context%window_set(iwindow)) THEN
            CALL section_vals_val_get(input, TRIM(window_keywords(iwindow)), r_val=context%energy_window(iwindow))
            context%energy_window(iwindow) = context%energy_window(iwindow)*evolt
         END IF
      END DO
      IF (ANY(context%window_set) .AND. num_wann == num_bands) THEN
         CPABORT('Wannier90 energy windows require more exported bands than Wannier functions.')
      END IF
      ! Each overlap block has exactly one owner, including when some ranks have no k-points.
      ALLOCATE (distk(nk), context%local_index(nk))
      context%local_index(:) = 0
      nlocal = 0
      DO ik = 1, nk
         distk(ik) = MOD(ik - 1, para_env%num_pe)
         IF (distk(ik) == para_env%mepos) THEN
            nlocal = nlocal + 1
            context%local_index(ik) = nlocal
         END IF
      END DO
      log_seed = seed_name
      IF (PRESENT(log_suffix)) log_seed = TRIM(seed_name)//log_suffix
      IF (para_env%is_source()) THEN
         context%keep_output = .TRUE.
         CALL open_file(TRIM(log_seed)//".library.wout", unit_number=context%output_unit, &
                        file_status="REPLACE", file_action="READWRITE")
      ELSE
         WRITE (rank_suffix, '(I0)') para_env%mepos
         CALL open_file(TRIM(log_seed)//".library-rank-"//TRIM(rank_suffix)//".wout", &
                        unit_number=context%output_unit, file_status="REPLACE", file_action="WRITE")
      END IF
#if defined(__parallel)
      CALL w90_set_comm(context%data, para_env%get_handle())
#endif
      CALL w90_set_option(context%data, "kpoints", kpoints)
      CALL w90_set_option(context%data, "mp_grid", mp_grid)
      CALL w90_set_option(context%data, "unit_cell_cart", lattice)
      CALL w90_set_option(context%data, "num_bands", num_bands)
      CALL w90_set_option(context%data, "num_wann", num_wann)
      CALL w90_set_option(context%data, "distk", distk)
      CALL w90_set_option(context%data, "num_iter", context%num_iter)
      CALL w90_set_option(context%data, "conv_tol", context%conv_tol)
      CALL w90_set_option(context%data, "conv_window", context%conv_window)
      CALL w90_set_option(context%data, "num_cg_steps", context%num_cg_steps)
      CALL w90_set_option(context%data, "num_print_cycles", context%num_print_cycles)
      CALL w90_set_option(context%data, 'write_hr', context%write_hr)
      CALL w90_set_option(context%data, 'write_u_matrices', context%write_u_matrices)
      CALL w90_set_option(context%data, 'dis_num_iter', context%dis_num_iter)
      CALL w90_set_option(context%data, 'dis_conv_window', context%dis_conv_window)
      CALL w90_set_option(context%data, 'dis_conv_tol', context%dis_conv_tol)
      CALL w90_set_option(context%data, 'dis_mix_ratio', context%dis_mix_ratio)
      DO iwindow = 1, SIZE(window_keywords)
         IF (.NOT. context%window_set(iwindow)) CYCLE
         CALL w90_set_option(context%data, TRIM(window_keywords(iwindow)), context%energy_window(iwindow))
      END DO
      CALL w90_input_setopt(context%data, seed_name, context%output_unit, context%output_unit, ierr)
      CALL check_w90_error(ierr, "initialization")
      CALL w90_get_nn(context%data, nn, context%output_unit, context%output_unit, ierr)
      CALL check_w90_error(ierr, "neighbour count")
      ALLOCATE (nnlist(nk, nn), nncell(3, nk, nn))
      CALL w90_get_nnkp(context%data, nnlist, context%output_unit, context%output_unit, ierr)
      CALL check_w90_error(ierr, "neighbour indices")
      CALL w90_get_gkpb(context%data, nncell, context%output_unit, context%output_unit, ierr)
      CALL check_w90_error(ierr, "reciprocal-cell shifts")
      ALLOCATE (context%overlaps(num_bands, num_bands, nn, nlocal))
      ALLOCATE (context%projections(num_bands, num_wann, nk))
      ALLOCATE (context%rotations(num_wann, num_wann, nk), context%eigenvalues(num_bands, nk))
      context%overlaps(:, :, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      context%projections(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      context%rotations(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      context%eigenvalues(:, :) = 0.0_dp
      DO iwann = 1, num_wann
         context%projections(iwann, iwann, :) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      END DO
      CALL w90_set_m_local(context%data, context%overlaps)
      CALL w90_set_u_opt(context%data, context%projections)
      CALL w90_set_u_matrix(context%data, context%rotations)
      CALL w90_set_eigval(context%data, context%eigenvalues)
      DEALLOCATE (distk)
#else
      CPABORT("WANNIER90 LIBRARY requires CP2K compiled with CP2K_USE_WANNIER90=ON.")
#endif
   END SUBROUTINE w90_library_init

! **************************************************************************************************
!> \brief Store a replicated overlap block on its owning rank.
!> \param context Library state.
!> \param ik Full-mesh k-point index.
!> \param inn Neighbour index in the library ordering.
!> \param matrix Complex Bloch-gauge overlap matrix.
! **************************************************************************************************
   SUBROUTINE w90_library_overlap(context, ik, inn, matrix)
      TYPE(w90_library_type), INTENT(INOUT)              :: context
      INTEGER, INTENT(IN)                                :: ik, inn
      COMPLEX(KIND=dp), INTENT(IN)                       :: matrix(:, :)

      INTEGER                                            :: ilocal

      ilocal = context%local_index(ik)
      IF (ilocal > 0) context%overlaps(:, :, inn, ilocal) = matrix(:, :)
   END SUBROUTINE w90_library_overlap

! **************************************************************************************************
!> \brief Store eigenvalues in the energy units used by Wannier90.
!> \param context Library state.
!> \param ik Full-mesh k-point index.
!> \param eigenvalues Eigenvalues in electronvolts.
! **************************************************************************************************
   SUBROUTINE w90_library_eigenvalues(context, ik, eigenvalues)
      TYPE(w90_library_type), INTENT(INOUT)              :: context
      INTEGER, INTENT(IN)                                :: ik
      REAL(KIND=dp), INTENT(IN)                          :: eigenvalues(:)

      context%eigenvalues(:, ik) = eigenvalues(:)
   END SUBROUTINE w90_library_eigenvalues

! **************************************************************************************************
!> \brief Restrict AO trial selection to the subspace admitted by explicit energy windows.
!> \param context Library state with the exported eigenvalues already populated.
!> \param ao Band/AO/k-point projections, empty on non-source ranks.
! **************************************************************************************************
   SUBROUTINE w90_library_project_trials(context, ao)
      TYPE(w90_library_type), INTENT(IN)                 :: context
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: ao(:, :, :)

      INTEGER                                            :: ib, ik, nb, nw
      LOGICAL, ALLOCATABLE                               :: allowed(:), frozen(:)

      nb = SIZE(ao, 1)
      IF (nb == 0 .OR. .NOT. ANY(context%window_set)) RETURN
      nw = SIZE(context%rotations, 1)
      CPASSERT(nb == SIZE(context%eigenvalues, 1))
      CPASSERT(SIZE(ao, 3) == SIZE(context%eigenvalues, 2))
      ALLOCATE (allowed(nb), frozen(nb))
      DO ik = 1, SIZE(ao, 3)
         allowed(:) = .TRUE.
         IF (context%window_set(1)) allowed(:) = allowed .AND. context%eigenvalues(:, ik) >= context%energy_window(1)
         IF (context%window_set(2)) allowed(:) = allowed .AND. context%eigenvalues(:, ik) <= context%energy_window(2)
         IF (COUNT(allowed) < nw) CPABORT('Wannier90 outer energy window contains too few bands.')
         IF (context%window_set(4)) THEN
            frozen(:) = context%eigenvalues(:, ik) <= context%energy_window(4)
            IF (context%window_set(3)) frozen(:) = frozen .AND. context%eigenvalues(:, ik) >= context%energy_window(3)
            IF (COUNT(frozen) > nw) CPABORT('Wannier90 frozen window contains too many bands.')
            IF (ANY(frozen .AND. .NOT. allowed)) THEN
               CPABORT('Wannier90 frozen states lie outside the outer window.')
            END IF
            ! If every target state is frozen there is no remaining subspace freedom:
            ! scoring trials against unoccupied states would select irrelevant AOs.
            IF (COUNT(frozen) == nw) allowed(:) = frozen(:)
         END IF
         DO ib = 1, nb
            IF (.NOT. allowed(ib)) ao(ib, :, ik) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         END DO
      END DO
      DEALLOCATE (allowed, frozen)
   END SUBROUTINE w90_library_project_trials

! **************************************************************************************************
!> \brief Copy immutable input matrices into an independently initialized library instance.
!> \param source Unused library instance holding the original overlaps and eigenvalues.
!> \param destination Fresh library instance with the same mesh, distribution and dimensions.
! **************************************************************************************************
   SUBROUTINE w90_library_copy_matrices(source, destination)
      TYPE(w90_library_type), INTENT(IN)                 :: source
      TYPE(w90_library_type), INTENT(INOUT)              :: destination

      CPASSERT(.NOT. source%inputs_consumed .AND. .NOT. destination%inputs_consumed)
      CPASSERT(ALL(SHAPE(source%overlaps) == SHAPE(destination%overlaps)))
      CPASSERT(ALL(source%local_index == destination%local_index))
      CPASSERT(ALL(SHAPE(source%eigenvalues) == SHAPE(destination%eigenvalues)))
      destination%overlaps = source%overlaps
      destination%eigenvalues = source%eigenvalues
   END SUBROUTINE w90_library_copy_matrices

! **************************************************************************************************
!> \brief Set the replicated initial projections and optionally export them before optimization.
!> \param context Library state.
!> \param amn Band/trial/k-point projections in the same Bloch gauge as the overlap matrices.
!> \param para_env Global parallel environment.
!> \param seed_name Output filename prefix.
!> \param write_files Whether to write a matching external projection file.
! **************************************************************************************************
   SUBROUTINE w90_library_projections(context, amn, para_env, seed_name, write_files)
      TYPE(w90_library_type), INTENT(INOUT)              :: context
      COMPLEX(KIND=dp), INTENT(IN)                       :: amn(:, :, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      CHARACTER(LEN=*), INTENT(IN)                       :: seed_name
      LOGICAL, INTENT(IN)                                :: write_files

      CPASSERT(.NOT. context%inputs_consumed)
      context%projections(:, :, :) = amn(:, :, :)
      IF (write_files) CALL w90_write_amn(amn, para_env, seed_name)
   END SUBROUTINE w90_library_projections

! **************************************************************************************************
!> \brief Write initial projections without overwriting the optimized library state.
!> \param amn Initial band/trial/k-point projections.
!> \param para_env Global parallel environment.
!> \param seed_name Output filename prefix.
! **************************************************************************************************
   SUBROUTINE w90_write_amn(amn, para_env, seed_name)
      COMPLEX(KIND=dp), INTENT(IN)                       :: amn(:, :, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      CHARACTER(LEN=*), INTENT(IN)                       :: seed_name

      INTEGER                                            :: ib, ik, iunit, iwann, nb, nk, nw

      IF (para_env%is_source()) THEN
         nb = SIZE(amn, 1)
         nw = SIZE(amn, 2)
         nk = SIZE(amn, 3)
         CALL open_file(TRIM(seed_name)//'.amn', unit_number=iunit, file_status='REPLACE', file_action='WRITE')
         WRITE (iunit, '(A)') '! Initial projections supplied to the CP2K Wannier90 library'
         WRITE (iunit, '(3I10)') nb, nk, nw
         DO ik = 1, nk
            DO iwann = 1, nw
               DO ib = 1, nb
                  WRITE (iunit, '(3I10,2ES24.16)') ib, iwann, ik, REAL(amn(ib, iwann, ik), KIND=dp), &
                     AIMAG(amn(ib, iwann, ik))
               END DO
            END DO
         END DO
         CALL close_file(iunit)
      END IF
   END SUBROUTINE w90_write_amn

! **************************************************************************************************
!> \brief Publish the winning trial log under the canonical output name.
!> \param context Winning library state; the original canonical log must already be closed.
!> \param seed_name Canonical output filename prefix.
! **************************************************************************************************
   SUBROUTINE w90_library_promote_output(context, seed_name)
      TYPE(w90_library_type), INTENT(INOUT)              :: context
      CHARACTER(LEN=*), INTENT(IN)                       :: seed_name

      CHARACTER(LEN=default_string_length)               :: line
      INTEGER                                            :: ios, output_unit

      IF (.NOT. context%keep_output) RETURN
      CALL open_file(TRIM(seed_name)//'.library.wout', unit_number=output_unit, &
                     file_status='REPLACE', file_action='READWRITE')
      REWIND (context%output_unit)
      DO
         READ (context%output_unit, '(A)', IOSTAT=ios) line
         IF (IS_IOSTAT_END(ios)) EXIT
         CPASSERT(ios == 0)
         WRITE (output_unit, '(A)') TRIM(line)
      END DO
      CALL close_file(context%output_unit)
      context%output_unit = output_unit
   END SUBROUTINE w90_library_promote_output

! **************************************************************************************************
!> \brief Write options matching the in-process calculation for external validation.
!> \param context Library state.
!> \param iunit Open Wannier90 input file.
! **************************************************************************************************
   SUBROUTINE w90_library_write_options(context, iunit)
      TYPE(w90_library_type), INTENT(IN)                 :: context
      INTEGER, INTENT(IN)                                :: iunit

      INTEGER                                            :: iwindow

      WRITE (iunit, '(A,I0)') "num_iter = ", context%num_iter
      WRITE (iunit, '(A,I0)') "num_cg_steps = ", context%num_cg_steps
      WRITE (iunit, '(A,I0)') "num_print_cycles = ", context%num_print_cycles
      WRITE (iunit, '(A,I0)') "conv_window = ", context%conv_window
      WRITE (iunit, '(A,ES24.16)') "conv_tol = ", context%conv_tol
      WRITE (iunit, '(A,L1)') 'write_hr = ', context%write_hr
      WRITE (iunit, '(A,L1)') 'write_u_matrices = ', context%write_u_matrices
      IF (SIZE(context%eigenvalues, 1) > SIZE(context%rotations, 1)) THEN
         WRITE (iunit, '(A,I0)') 'dis_num_iter = ', context%dis_num_iter
         WRITE (iunit, '(A,I0)') 'dis_conv_window = ', context%dis_conv_window
         WRITE (iunit, '(A,ES24.16)') 'dis_conv_tol = ', context%dis_conv_tol
         WRITE (iunit, '(A,ES24.16)') 'dis_mix_ratio = ', context%dis_mix_ratio
         DO iwindow = 1, SIZE(window_keywords)
            IF (.NOT. context%window_set(iwindow)) CYCLE
            WRITE (iunit, '(2A,ES24.16)') TRIM(window_keywords(iwindow)), ' = ', context%energy_window(iwindow)
         END DO
      END IF
   END SUBROUTINE w90_library_write_options

! **************************************************************************************************
!> \brief Write a complete external input matching this library instance and its optimizer.
!> \param context Library state providing dimensions and localization options.
!> \param para_env Global parallel environment.
!> \param seed_name Output filename prefix.
!> \param lattice Real lattice vectors in columns, in angstrom.
!> \param kpoints Full mesh in fractional coordinates.
!> \param mp_grid Mesh dimensions.
!> \param atoms_cart Atomic positions in bohr.
!> \param atom_symbols Element symbols.
! **************************************************************************************************
   SUBROUTINE w90_library_write_input(context, para_env, seed_name, lattice, kpoints, mp_grid, atoms_cart, atom_symbols)
      TYPE(w90_library_type), INTENT(IN)                 :: context
      TYPE(mp_para_env_type), POINTER                    :: para_env
      CHARACTER(LEN=*), INTENT(IN)                       :: seed_name
      REAL(KIND=dp), INTENT(IN)                          :: lattice(3, 3), kpoints(:, :)
      INTEGER, INTENT(IN)                                :: mp_grid(3)
      REAL(KIND=dp), INTENT(IN)                          :: atoms_cart(:, :)
      CHARACTER(LEN=*), INTENT(IN)                       :: atom_symbols(:)

      INTEGER                                            :: i, iunit

      IF (.NOT. para_env%is_source()) RETURN
      CALL open_file(TRIM(seed_name)//'.win', unit_number=iunit, file_status='REPLACE', file_action='WRITE')
      WRITE (iunit, '(A)') '! CP2K Wannier90 library input'
      WRITE (iunit, '(A,I0)') 'num_wann = ', SIZE(context%rotations, 1)
      WRITE (iunit, '(A,I0)') 'num_bands = ', SIZE(context%eigenvalues, 1)
      CALL w90_library_write_options(context, iunit)
      WRITE (iunit, '(A)') 'length_unit = ang', 'begin unit_cell_cart', 'ang'
      DO i = 1, 3
         WRITE (iunit, '(3ES24.16)') lattice(:, i)
      END DO
      WRITE (iunit, '(A)') 'end unit_cell_cart', 'begin atoms_cart', 'bohr'
      DO i = 1, SIZE(atom_symbols)
         WRITE (iunit, '(A,3ES24.16)') atom_symbols(i), atoms_cart(:, i)
      END DO
      WRITE (iunit, '(A)') 'end atoms_cart'
      WRITE (iunit, '(A,3I6)') 'mp_grid = ', mp_grid
      WRITE (iunit, '(A)') 'begin kpoints'
      DO i = 1, SIZE(kpoints, 2)
         WRITE (iunit, '(3ES24.16)') kpoints(:, i)
      END DO
      WRITE (iunit, '(A)') 'end kpoints'
      CALL close_file(iunit)
   END SUBROUTINE w90_library_write_input

! **************************************************************************************************
!> \brief Localize the supplied subspace and return centres and spreads in atomic units.
!> \param context Library state.
!> \param para_env Communicator shared by all export ranks.
!> \param centres Wannier centres in bohr.
!> \param spreads Individual quadratic spreads in bohr squared.
!> \param converged Whether Wannier90 satisfied its localization convergence criterion.
! **************************************************************************************************
   SUBROUTINE w90_library_run(context, para_env, centres, spreads, converged)
      TYPE(w90_library_type), INTENT(INOUT) :: context
      TYPE(mp_para_env_type), POINTER :: para_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: centres(:, :), spreads(:)
      LOGICAL, INTENT(OUT) :: converged
#if defined(__WANNIER90)
      INTEGER :: ierr, nw
      LOGICAL :: dis_converged
      COMPLEX(KIND=dp), POINTER :: reduced_overlaps(:, :, :, :)

      CPASSERT(.NOT. context%inputs_consumed)
      context%inputs_consumed = .TRUE.
      nw = SIZE(context%rotations, 1)
      dis_converged = .TRUE.
      IF (SIZE(context%eigenvalues, 1) > nw) THEN
         CALL cite_reference(Souza2001)
         CALL w90_disentangle(context%data, context%output_unit, context%output_unit, ierr)
         CALL check_w90_error(ierr, 'disentanglement')
         dis_converged = w90_library_report_contains(context, 'Disentanglement convergence criteria satisfied')
         CALL para_env%bcast(dis_converged)
         ! The v4 API leaves the optimized overlaps in the leading nw-by-nw block
         ! of nb-by-nb storage. Wannierisation's BLAS calls require compact nw strides.
         ALLOCATE (reduced_overlaps(nw, nw, SIZE(context%overlaps, 3), SIZE(context%overlaps, 4)))
         reduced_overlaps(:, :, :, :) = context%overlaps(1:nw, 1:nw, :, :)
         CALL w90_set_m_local(context%data, reduced_overlaps)
         DEALLOCATE (context%overlaps)
         context%overlaps => reduced_overlaps
         NULLIFY (reduced_overlaps)
      END IF
      CALL w90_project_overlap(context%data, context%output_unit, context%output_unit, ierr)
      CALL check_w90_error(ierr, "projection")
      CALL w90_wannierise(context%data, context%output_unit, context%output_unit, ierr)
      CALL check_w90_error(ierr, "localization")
      ! Wannier90 4.0.2 does not return nonconvergence in ierr. Fail closed if its
      ! explicit convergence report is absent; do not select an unconverged minimum.
      converged = w90_library_report_contains(context, 'Wannierisation convergence criteria satisfied')
      converged = converged .AND. dis_converged
      CALL para_env%bcast(converged)
      ALLOCATE (centres(3, nw), spreads(nw))
      CALL w90_get_centres(context%data, centres)
      CALL w90_get_spreads(context%data, spreads)
      centres(:, :) = centres(:, :)/angstrom
      spreads(:) = spreads(:)/angstrom**2
      IF (para_env%is_source()) THEN
         converged = converged .AND. ALL(ieee_is_finite(centres)) .AND. ALL(ieee_is_finite(spreads))
      END IF
      CALL para_env%bcast(converged)
      context%converged = converged
#else
      CPABORT("WANNIER90 library is not available.")
#endif
   END SUBROUTINE w90_library_run

! **************************************************************************************************
!> \brief Write native optimized matrices only after selecting a converged library result.
!> \param context Converged state, including the winning trial for AUTO.
!> \param para_env Communicator shared by all export ranks.
! **************************************************************************************************
   SUBROUTINE w90_library_write_results(context, para_env)
      TYPE(w90_library_type), INTENT(INOUT)              :: context
      TYPE(mp_para_env_type), POINTER                    :: para_env
#if defined(__WANNIER90)
      INTEGER                                            :: ierr, iunit, ik, ib, packed_band, mapped_band

      CPASSERT(context%converged)
      IF (.NOT. (context%write_hr .OR. context%write_u_matrices)) RETURN
      CALL w90_plot(context%data, context%output_unit, context%output_unit, ierr)
      CALL check_w90_error(ierr, 'optimized output')
      IF (context%write_u_matrices .AND. para_env%is_source()) THEN
         ! Native U_dis/V rows follow the packed outer window, not the original
         ! band indices. Keep Wannier90's format and publish that mapping explicitly.
         CALL open_file(TRIM(context%data%seedname)//'_band_map.dat', unit_number=iunit, &
                        file_status='REPLACE', file_action='WRITE')
         WRITE (iunit, '(A)') '# k-point, exported band, matrix row (0 outside window), eigenvalue [eV]'
         WRITE (iunit, '(2I10)') SIZE(context%eigenvalues, 2), SIZE(context%eigenvalues, 1)
         DO ik = 1, SIZE(context%eigenvalues, 2)
            packed_band = 0
            DO ib = 1, SIZE(context%eigenvalues, 1)
               mapped_band = ib
               IF (context%data%have_disentangled) THEN
                  mapped_band = 0
                  IF (context%data%dis_manifold%lwindow(ib, ik)) THEN
                     packed_band = packed_band + 1
                     mapped_band = packed_band
                  END IF
               END IF
               WRITE (iunit, '(3I10,ES24.16)') ik, ib, mapped_band, context%eigenvalues(ib, ik)
            END DO
         END DO
         CALL close_file(iunit)
      END IF
#else
      CPABORT('WANNIER90 library is not available.')
#endif
   END SUBROUTINE w90_library_write_results

! **************************************************************************************************
!> \brief Release a library instance and the arrays whose pointers it borrowed.
!> \param context Library state.
! **************************************************************************************************
   SUBROUTINE w90_library_release(context)
      TYPE(w90_library_type), INTENT(INOUT) :: context

#if defined(__WANNIER90)
      DEALLOCATE (context%data)
#endif
      DEALLOCATE (context%overlaps, context%projections, context%rotations, context%eigenvalues)
      NULLIFY (context%overlaps, context%projections, context%rotations, context%eigenvalues)
      DEALLOCATE (context%local_index)
      IF (context%keep_output) THEN
         CALL close_file(context%output_unit)
      ELSE
         CALL close_file(context%output_unit, file_status="DELETE")
      END IF
      context%output_unit = -1
   END SUBROUTINE w90_library_release

#if defined(__WANNIER90)
! **************************************************************************************************
!> \brief Read a library convergence report and restore its position for subsequent output.
!> \param context Library instance, whose root log is readable.
!> \param marker Exact success marker emitted by the library.
!> \return Whether the source-rank log contains the marker; false on other ranks.
! **************************************************************************************************
   FUNCTION w90_library_report_contains(context, marker) RESULT(found)
      TYPE(w90_library_type), INTENT(IN)                 :: context
      CHARACTER(LEN=*), INTENT(IN)                       :: marker
      LOGICAL                                            :: found

      CHARACTER(LEN=default_string_length)               :: line
      INTEGER                                            :: ios

      found = .FALSE.
      IF (.NOT. context%keep_output) RETURN
      FLUSH (context%output_unit)
      REWIND (context%output_unit)
      DO
         READ (context%output_unit, '(A)', IOSTAT=ios) line
         IF (IS_IOSTAT_END(ios)) EXIT
         CPASSERT(ios == 0)
         IF (INDEX(line, marker) > 0) found = .TRUE.
      END DO
      BACKSPACE (context%output_unit)
   END FUNCTION w90_library_report_contains

! **************************************************************************************************
!> \brief Translate library errors to CP2K's collective error handling.
!> \param ierr Library error code.
!> \param operation Operation being checked.
! **************************************************************************************************
   SUBROUTINE check_w90_error(ierr, operation)
      INTEGER, INTENT(IN)                                :: ierr
      CHARACTER(LEN=*), INTENT(IN)                       :: operation

      IF (ierr /= 0) THEN
         CPABORT("Wannier90 library failed during "//operation//"; see the .library.wout file.")
      END IF
   END SUBROUTINE check_w90_error
#endif

END MODULE wannier90_library
