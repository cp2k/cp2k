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
   USE cp_files, ONLY: close_file, open_file
   USE input_section_types, ONLY: section_vals_type, section_vals_val_get
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_para_env_type
   USE physcon, ONLY: angstrom
#if defined(__WANNIER90)
   USE w90_library, ONLY: lib_common_type, w90_get_centres, w90_get_gkpb, w90_get_nn, &
                          w90_get_nnkp, w90_get_spreads, w90_input_setopt, w90_project_overlap, &
                          w90_set_comm, w90_set_eigval, w90_set_m_local, w90_set_option, &
                          w90_set_u_matrix, w90_set_u_opt, w90_wannierise
#endif
#include "./base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'wannier90_library'

   TYPE w90_library_type
      PRIVATE
#if defined(__WANNIER90)
      TYPE(lib_common_type), ALLOCATABLE :: data
#endif
      INTEGER :: output_unit = -1, num_iter = 0, conv_window = 0
      LOGICAL :: keep_output = .FALSE.
      REAL(KIND=dp) :: conv_tol = 0.0_dp
      INTEGER, ALLOCATABLE :: local_index(:)
      COMPLEX(KIND=dp), POINTER :: overlaps(:, :, :, :) => NULL()
      COMPLEX(KIND=dp), POINTER :: projections(:, :, :) => NULL(), rotations(:, :, :) => NULL()
      REAL(KIND=dp), POINTER :: eigenvalues(:, :) => NULL()
   END TYPE w90_library_type

   PUBLIC :: w90_library_type, w90_library_init, w90_library_overlap, w90_library_eigenvalues, &
             w90_library_run, w90_library_release, w90_library_write_options, w90_library_projections

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
! **************************************************************************************************
   SUBROUTINE w90_library_init(context, input, para_env, seed_name, lattice, kpoints, mp_grid, &
                               num_bands, num_wann, nnlist, nncell)
      TYPE(w90_library_type), INTENT(INOUT) :: context
      TYPE(section_vals_type), POINTER :: input
      TYPE(mp_para_env_type), POINTER :: para_env
      CHARACTER(LEN=*), INTENT(IN) :: seed_name
      REAL(KIND=dp), INTENT(IN) :: lattice(3, 3), kpoints(:, :)
      INTEGER, INTENT(IN) :: mp_grid(3), num_bands, num_wann
      INTEGER, ALLOCATABLE, INTENT(OUT) :: nnlist(:, :), nncell(:, :, :)
#if defined(__WANNIER90)
      INTEGER :: ierr, ik, iwann, nk, nn, nlocal
      INTEGER, ALLOCATABLE :: distk(:)
      CHARACTER(LEN=16) :: rank_suffix

      CPASSERT(.NOT. ALLOCATED(context%data))
      ALLOCATE (context%data)
      IF (LEN_TRIM(seed_name) > LEN(context%data%seedname)) THEN
         CPABORT("SEED_NAME exceeds the filename prefix length supported by the Wannier90 library.")
      END IF
      nk = SIZE(kpoints, 2)
      CALL section_vals_val_get(input, "NUM_ITER", i_val=context%num_iter)
      CALL section_vals_val_get(input, "CONV_WINDOW", i_val=context%conv_window)
      CALL section_vals_val_get(input, "CONV_TOL", r_val=context%conv_tol)
      IF (context%num_iter < 1 .OR. context%conv_window < 2 .OR. context%conv_tol <= 0.0_dp) THEN
         CPABORT("WANNIER90 LIBRARY requires NUM_ITER > 0, CONV_WINDOW >= 2 and CONV_TOL > 0.")
      END IF
      IF (num_wann < 1 .OR. num_wann /= num_bands) THEN
         CPABORT("WANNIER90 LIBRARY identity projections require a complete band subspace.")
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
      IF (para_env%is_source()) THEN
         context%keep_output = .TRUE.
         CALL open_file(TRIM(seed_name)//".library.wout", unit_number=context%output_unit, &
                        file_status="REPLACE", file_action="READWRITE")
      ELSE
         WRITE (rank_suffix, '(I0)') para_env%mepos
         CALL open_file(TRIM(seed_name)//".library-rank-"//TRIM(rank_suffix)//".wout", &
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

      INTEGER                                            :: ib, ik, iunit, iwann, nb, nk, nw

      context%projections(:, :, :) = amn(:, :, :)
      IF (write_files .AND. para_env%is_source()) THEN
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
   END SUBROUTINE w90_library_projections

! **************************************************************************************************
!> \brief Write options matching the in-process calculation for external validation.
!> \param context Library state.
!> \param iunit Open Wannier90 input file.
! **************************************************************************************************
   SUBROUTINE w90_library_write_options(context, iunit)
      TYPE(w90_library_type), INTENT(IN)                 :: context
      INTEGER, INTENT(IN)                                :: iunit

      WRITE (iunit, '(A,I0)') "num_iter = ", context%num_iter
      WRITE (iunit, '(A,I0)') "conv_window = ", context%conv_window
      WRITE (iunit, '(A,ES24.16)') "conv_tol = ", context%conv_tol
   END SUBROUTINE w90_library_write_options

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
      CHARACTER(LEN=256) :: line
      INTEGER :: ierr, ios, nw

      CALL w90_project_overlap(context%data, context%output_unit, context%output_unit, ierr)
      CALL check_w90_error(ierr, "projection")
      CALL w90_wannierise(context%data, context%output_unit, context%output_unit, ierr)
      CALL check_w90_error(ierr, "localization")
      ! Wannier90 4.0.2 does not return nonconvergence in ierr. Fail closed if its
      ! explicit convergence report is absent; do not select an unconverged minimum.
      converged = .FALSE.
      IF (para_env%is_source()) THEN
         FLUSH (context%output_unit)
         REWIND (context%output_unit)
         DO
            READ (context%output_unit, '(A)', IOSTAT=ios) line
            IF (IS_IOSTAT_END(ios)) EXIT
            CPASSERT(ios == 0)
            IF (INDEX(line, 'Wannierisation convergence criteria satisfied') > 0) converged = .TRUE.
         END DO
      END IF
      CALL para_env%bcast(converged)
      nw = SIZE(context%rotations, 1)
      ALLOCATE (centres(3, nw), spreads(nw))
      CALL w90_get_centres(context%data, centres)
      CALL w90_get_spreads(context%data, spreads)
      centres(:, :) = centres(:, :)/angstrom
      spreads(:) = spreads(:)/angstrom**2
#else
      CPABORT("WANNIER90 library is not available.")
#endif
   END SUBROUTINE w90_library_run

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
