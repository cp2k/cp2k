!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the calculation of wannier states
!> \author Alin M Elena
! **************************************************************************************************
MODULE wannier_states
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_element,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indco,&
                                              nco,&
                                              nso
   USE orbital_symbols,                 ONLY: cgf_symbol,&
                                              sgf_symbol
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
   USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE wannier_states_types,            ONLY: wannier_centres_type
!!!! this ones are needed to mapping
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'wannier_states'

   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .TRUE.

! *** Public subroutines ***

   PUBLIC :: construct_wannier_states

CONTAINS

! **************************************************************************************************
!> \brief constructs wannier states. mo_localized should not be overwritten!
!> \param mo_localized ...
!> \param Hks ...
!> \param qs_env ...
!> \param loc_print_section ...
!> \param WannierCentres ...
!> \param ns ...
!> \param states ...
!> \par History
!>      - Printout Wannier states in AO basis (11.09.2025, H. Elgabarty)
! **************************************************************************************************
   SUBROUTINE construct_wannier_states(mo_localized, &
                                       Hks, qs_env, loc_print_section, WannierCentres, ns, states)

      TYPE(cp_fm_type), INTENT(in)                       :: mo_localized
      TYPE(dbcsr_type), POINTER                          :: Hks
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: loc_print_section
      TYPE(wannier_centres_type), INTENT(INOUT)          :: WannierCentres
      INTEGER, INTENT(IN)                                :: ns
      INTEGER, INTENT(IN), POINTER                       :: states(:)

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_wannier_states'

      CHARACTER(default_string_length)                   :: unit_str
      CHARACTER(LEN=12)                                  :: symbol
      CHARACTER(LEN=12), DIMENSION(:), POINTER           :: bcgf_symbol
      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=40)                                  :: fmtstr1, fmtstr2, fmtstr3
      CHARACTER(LEN=6), DIMENSION(:), POINTER            :: bsgf_symbol
      INTEGER :: after, before, from, handle, i, iatom, icgf, ico, icol, ikind, iproc, irow, iset, &
         isgf, ishell, iso, jcol, left, lmax, lshell, natom, ncgf, ncol, ncol_global, nrow_global, &
         nset, nsgf, nstates(2), output_unit, right, to, unit_mat
      INTEGER, DIMENSION(:), POINTER                     :: nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: l
      LOGICAL                                            :: print_cartesian
      REAL(KIND=dp)                                      :: unit_conv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cmatrix, smatrix
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type)                                   :: b, c, d
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_parameter
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: print_key

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      CALL timeset(routineN, handle)
      NULLIFY (logger, para_env)

      CALL get_qs_env(qs_env, para_env=para_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set)

      logger => cp_get_default_logger()

      output_unit = cp_logger_get_default_io_unit(logger)

      CALL cp_fm_get_info(mo_localized, &
                          ncol_global=ncol_global, &
                          nrow_global=nrow_global)

      nstates(1) = ns
      nstates(2) = para_env%mepos
      iproc = nstates(2)
      NULLIFY (fm_struct_tmp, print_key)

      print_key => section_vals_get_subs_vals(loc_print_section, "WANNIER_CENTERS")
      CALL section_vals_val_get(print_key, "UNIT", c_val=unit_str)
      unit_conv = cp_unit_from_cp2k(1.0_dp, TRIM(unit_str))

      print_key => section_vals_get_subs_vals(loc_print_section, "WANNIER_STATES")
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nrow_global, &
                               ncol_global=1, &
                               para_env=mo_localized%matrix_struct%para_env, &
                               context=mo_localized%matrix_struct%context)

      CALL cp_fm_create(b, fm_struct_tmp, name="b")
      CALL cp_fm_create(c, fm_struct_tmp, name="c")

      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=1, ncol_global=1, &
                               para_env=mo_localized%matrix_struct%para_env, &
                               context=mo_localized%matrix_struct%context)

      CALL cp_fm_create(d, fm_struct_tmp, name="d")
      CALL cp_fm_struct_release(fm_struct_tmp)

      WannierCentres%WannierHamDiag = 0.0_dp

      unit_mat = cp_print_key_unit_nr(logger, loc_print_section, &
                                      "WANNIER_STATES", extension=".whks", &
                                      ignore_should_output=.FALSE.)
      IF (unit_mat > 0) THEN
         WRITE (unit_mat, '(a16,1(i0,1x))') "Wannier states: ", ns
         WRITE (unit_mat, '(a16)') "#No x y z energy "
      END IF
      DO i = 1, ns
         CALL cp_fm_to_fm(mo_localized, b, 1, states(i), 1)
         CALL cp_dbcsr_sm_fm_multiply(Hks, b, c, 1)
         CALL parallel_gemm('T', 'N', 1, 1, nrow_global, 1.0_dp, &
                            b, c, 0.0_dp, d)
         CALL cp_fm_get_element(d, 1, 1, WannierCentres%WannierHamDiag(i))
         !               if (iproc==para_env%mepos) WRITE(unit_mat,'(f16.8,2x)', advance='no')WannierCentres%WannierHamDiag(i)
         IF (unit_mat > 0) WRITE (unit_mat, '(i0,1x,4(f16.8,2x))') states(i), &
            WannierCentres%centres(1:3, states(i))*unit_conv, WannierCentres%WannierHamDiag(states(i))
      END DO

      IF (unit_mat > 0) WRITE (unit_mat, *)

      IF (output_unit > 0) THEN
         WRITE (output_unit, *) ""
         WRITE (output_unit, *) "NUMBER OF Wannier STATES  ", ns
         WRITE (output_unit, *) "ENERGY      original MO-index"
         DO i = 1, ns
            WRITE (output_unit, '(f16.8,2x,i0)') WannierCentres%WannierHamDiag(i), states(i)
         END DO
      END IF

      CALL cp_fm_release(b)
      CALL cp_fm_release(c)
      CALL cp_fm_release(d)

      ! Print the states in AO basis
      CALL section_vals_val_get(print_key, "CARTESIAN", l_val=print_cartesian)

      ALLOCATE (smatrix(nrow_global, ncol_global))
      CALL cp_fm_get_submatrix(mo_localized, smatrix(1:nrow_global, 1:ncol_global))

      IF (unit_mat > 0) THEN

         NULLIFY (nshell)
         NULLIFY (bsgf_symbol)
         NULLIFY (l)
         NULLIFY (bsgf_symbol)

         CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
         CALL get_qs_kind_set(qs_kind_set, ncgf=ncgf, nsgf=nsgf)

         ! Print header, define column widths and string templates
         after = 6
         before = 4
         ncol = INT(56/(before + after + 3))

         fmtstr1 = "(T2,A,21X,  (  X,I5,  X))"
         fmtstr2 = "(T2,A,9X,  (1X,F  .  ))"
         fmtstr3 = "(T2,A,I5,1X,I5,1X,A,1X,A6,  (1X,F  .  ))"

         right = MAX((after - 2), 1)
         left = (before + after + 3) - right - 5

         IF (print_cartesian) THEN
            WRITE (UNIT=unit_mat, FMT="(T2,A,16X,A)") "WS|", "Wannier states in the cartesian AO basis"
         ELSE
            WRITE (UNIT=unit_mat, FMT="(T2,A,16X,A)") "WS|", "Wannier states in the spherical AO basis"
         END IF
         WRITE (UNIT=fmtstr1(11:12), FMT="(I2)") ncol
         WRITE (UNIT=fmtstr1(14:15), FMT="(I2)") left
         WRITE (UNIT=fmtstr1(21:22), FMT="(I2)") right

         WRITE (UNIT=fmtstr2(10:11), FMT="(I2)") ncol
         WRITE (UNIT=fmtstr2(17:18), FMT="(I2)") before + after + 2
         WRITE (UNIT=fmtstr2(20:21), FMT="(I2)") after

         WRITE (UNIT=fmtstr3(27:28), FMT="(I2)") ncol
         WRITE (UNIT=fmtstr3(34:35), FMT="(I2)") before + after + 2
         WRITE (UNIT=fmtstr3(37:38), FMT="(I2)") after

         ! get MO coefficients in terms of contracted cartesian functions
         IF (print_cartesian) THEN

            ALLOCATE (cmatrix(ncgf, ncgf))
            cmatrix = 0.0_dp

            ! Transform spherical to Cartesian AO basis
            icgf = 1
            isgf = 1
            DO iatom = 1, natom
               NULLIFY (orb_basis_set, dftb_parameter)
               CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), &
                                basis_set=orb_basis_set, &
                                dftb_parameter=dftb_parameter)
               IF (ASSOCIATED(orb_basis_set)) THEN
                  CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                         nset=nset, &
                                         nshell=nshell, &
                                         l=l)
                  DO iset = 1, nset
                     DO ishell = 1, nshell(iset)
                        lshell = l(ishell, iset)
                        CALL dgemm("T", "N", nco(lshell), ncol_global, nso(lshell), 1.0_dp, &
                                   orbtramat(lshell)%c2s, nso(lshell), &
                                   smatrix(isgf, 1), nsgf, 0.0_dp, &
                                   cmatrix(icgf, 1), ncgf)
                        icgf = icgf + nco(lshell)
                        isgf = isgf + nso(lshell)
                     END DO
                  END DO
               ELSE IF (ASSOCIATED(dftb_parameter)) THEN
                  CALL get_dftb_atom_param(dftb_parameter, lmax=lmax)
                  DO ishell = 1, lmax + 1
                     lshell = ishell - 1
                     CALL dgemm("T", "N", nco(lshell), nsgf, nso(lshell), 1.0_dp, &
                                orbtramat(lshell)%c2s, nso(lshell), &
                                smatrix(isgf, 1), nsgf, 0.0_dp, &
                                cmatrix(icgf, 1), ncgf)
                     icgf = icgf + nco(lshell)
                     isgf = isgf + nso(lshell)
                  END DO
               ELSE
                  ! assume atom without basis set
               END IF
            END DO ! iatom

         END IF

         ! Print to file
         DO icol = 1, ncol_global, ncol

            from = icol
            to = MIN((from + ncol - 1), ncol_global)

            WRITE (UNIT=unit_mat, FMT="(T2,A)") "WS|"
            WRITE (UNIT=unit_mat, FMT=fmtstr1) "WS|", (jcol, jcol=from, to)
            WRITE (UNIT=unit_mat, FMT=fmtstr2) "WS|    Energies", &
               (WannierCentres%WannierHamDiag(states(jcol)), jcol=from, to)
            WRITE (UNIT=unit_mat, FMT="(T2,A)") "WS|"

            irow = 1

            DO iatom = 1, natom

               IF (iatom /= 1) WRITE (UNIT=unit_mat, FMT="(T2,A)") "WS|"

               NULLIFY (orb_basis_set, dftb_parameter)
               CALL get_atomic_kind(particle_set(iatom)%atomic_kind, &
                                    element_symbol=element_symbol, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, &
                                dftb_parameter=dftb_parameter)

               IF (print_cartesian) THEN

                  NULLIFY (bcgf_symbol)
                  IF (ASSOCIATED(orb_basis_set)) THEN
                     CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                            nset=nset, &
                                            nshell=nshell, &
                                            l=l, &
                                            cgf_symbol=bcgf_symbol)

                     icgf = 1
                     DO iset = 1, nset
                        DO ishell = 1, nshell(iset)
                           lshell = l(ishell, iset)
                           DO ico = 1, nco(lshell)
                              WRITE (UNIT=unit_mat, FMT=fmtstr3) &
                                 "WS|", irow, iatom, ADJUSTR(element_symbol), bcgf_symbol(icgf), &
                                 (cmatrix(irow, jcol), jcol=from, to)
                              icgf = icgf + 1
                              irow = irow + 1
                           END DO
                        END DO
                     END DO
                  ELSE IF (ASSOCIATED(dftb_parameter)) THEN
                     CALL get_dftb_atom_param(dftb_parameter, lmax=lmax)
                     icgf = 1
                     DO ishell = 1, lmax + 1
                        lshell = ishell - 1
                        DO ico = 1, nco(lshell)
                           symbol = cgf_symbol(1, indco(1:3, icgf))
                           symbol(1:2) = "  "
                           WRITE (UNIT=unit_mat, FMT=fmtstr3) &
                              "WS|", irow, iatom, ADJUSTR(element_symbol), symbol, &
                              (cmatrix(irow, jcol), jcol=from, to)
                           icgf = icgf + 1
                           irow = irow + 1
                        END DO
                     END DO
                  ELSE
                     ! assume atom without basis set
                  END IF

               ELSE !print in spherical AO basis

                  IF (ASSOCIATED(orb_basis_set)) THEN
                     CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                            nset=nset, &
                                            nshell=nshell, &
                                            l=l, &
                                            sgf_symbol=bsgf_symbol)
                     isgf = 1
                     DO iset = 1, nset
                        DO ishell = 1, nshell(iset)
                           lshell = l(ishell, iset)
                           DO iso = 1, nso(lshell)
                              WRITE (UNIT=unit_mat, FMT=fmtstr3) &
                                 "WS|", irow, iatom, ADJUSTR(element_symbol), bsgf_symbol(isgf), &
                                 (smatrix(irow, jcol), jcol=from, to)
                              isgf = isgf + 1
                              irow = irow + 1
                           END DO
                        END DO
                     END DO
                  ELSE IF (ASSOCIATED(dftb_parameter)) THEN
                     CALL get_dftb_atom_param(dftb_parameter, lmax=lmax)
                     isgf = 1
                     DO ishell = 1, lmax + 1
                        lshell = ishell - 1
                        DO iso = 1, nso(lshell)
                           symbol = sgf_symbol(1, lshell, -lshell + iso - 1)
                           symbol(1:2) = "  "
                           WRITE (UNIT=unit_mat, FMT=fmtstr3) &
                              "WS|", irow, iatom, ADJUSTR(element_symbol), symbol, &
                              (smatrix(irow, jcol), jcol=from, to)
                           isgf = isgf + 1
                           irow = irow + 1
                        END DO
                     END DO
                  ELSE
                     ! assume atom without basis set
                  END IF

               END IF ! print cartesian

            END DO ! iatom

         END DO ! icol

         WRITE (UNIT=unit_mat, FMT="(T2,A)") "WS|"

         IF (print_cartesian) THEN
            DEALLOCATE (cmatrix)
         END IF

      END IF ! output Wannier states in AO
      DEALLOCATE (smatrix)

      CALL cp_print_key_finished_output(unit_mat, logger, loc_print_section, &
                                        "WANNIER_STATES")

      CALL timestop(handle)
   END SUBROUTINE construct_wannier_states

END MODULE wannier_states

