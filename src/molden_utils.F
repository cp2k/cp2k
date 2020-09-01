!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Functions handling the MOLDEN format. Split from mode_selective.
!> \author Teodoro Laino, 03.2009
! **************************************************************************************************
MODULE molden_utils
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_get_submatrix
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE input_constants,                 ONLY: gto_cartesian,&
                                              gto_spherical
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: nco,&
                                              nso
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: get_ptable_info
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: mo_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'molden_utils'
   LOGICAL, PARAMETER                   :: debug_this_module = .FALSE.

   INTEGER, PARAMETER                   :: molden_lmax = 4
   INTEGER, PARAMETER                   :: molden_ncomax = (molden_lmax + 1)*(molden_lmax + 2)/2 ! 15

   PUBLIC :: write_vibrations_molden, write_mos_molden

CONTAINS

! **************************************************************************************************
!> \brief Write out the MOs in molden format for visualisation
!> \param mos the set of MOs (both spins, if UKS)
!> \param qs_kind_set for basis set info
!> \param particle_set particles data structure, for positions and kinds
!> \param print_section input section containing relevant print key
!> \author MattW, IainB
! **************************************************************************************************
   SUBROUTINE write_mos_molden(mos, qs_kind_set, particle_set, print_section)
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: print_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'write_mos_molden'
      CHARACTER(LEN=molden_lmax+1), PARAMETER            :: angmom = "spdfg"

      CHARACTER(LEN=15)                                  :: fmtstr1
      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER :: gto_kind, handle, i, iatom, icgf, icol, ikind, ipgf, irow, irow_in, iset, isgf, &
         ishell, ispin, iw, lshell, ncgf, ncol_global, ndigits, nrow_global, nset, nsgf, z
      INTEGER, DIMENSION(:), POINTER                     :: npgf, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: l
      INTEGER, DIMENSION(molden_ncomax, 0:molden_lmax)   :: orbmap
      LOGICAL                                            :: print_warn
      REAL(KIND=dp)                                      :: expzet, prefac
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cmatrix, smatrix
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_section, ""), cp_p_file)) THEN

         iw = cp_print_key_unit_nr(logger, print_section, "", &
                                   extension=".molden", file_status='REPLACE')

         print_warn = .TRUE.

         CALL section_vals_val_get(print_section, "NDIGITS", i_val=ndigits)
         IF (ndigits < 1) ndigits = 3
         IF (ndigits > 92) ndigits = 92
         WRITE (UNIT=fmtstr1, FMT='("(I6,1X,ES",I0,".",I0,")")') ndigits + 7, ndigits

         CALL section_vals_val_get(print_section, "GTO_KIND", i_val=gto_kind)

         IF (mos(1)%mo_set%use_mo_coeff_b) THEN
            ! we are using the dbcsr mo_coeff
            ! we copy it to the fm anyway
            DO ispin = 1, SIZE(mos)
               IF (.NOT. ASSOCIATED(mos(ispin)%mo_set%mo_coeff_b)) THEN
                  CPASSERT(.FALSE.)
               ENDIF
               CALL copy_dbcsr_to_fm(mos(ispin)%mo_set%mo_coeff_b, &
                                     mos(ispin)%mo_set%mo_coeff) !fm->dbcsr
            ENDDO
         ENDIF

         IF (iw > 0) THEN
            WRITE (iw, '(T2,A)') "[Molden Format]"
            WRITE (iw, '(T2,A)') "[Atoms] AU"
            DO i = 1, SIZE(particle_set)
               CALL get_atomic_kind(atomic_kind=particle_set(i)%atomic_kind, &
                                    element_symbol=element_symbol)
               CALL get_ptable_info(element_symbol, number=z)

               WRITE (iw, '(T2,A2,I8,I8,3X,3(F12.6,3X))') &
                  element_symbol, i, z, particle_set(i)%r(:)
            END DO

            WRITE (iw, '(T2,A)') "[GTO]"

            DO i = 1, SIZE(particle_set)
               CALL get_atomic_kind(atomic_kind=particle_set(i)%atomic_kind, kind_number=ikind, &
                                    element_symbol=element_symbol)
               CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
               IF (ASSOCIATED(orb_basis_set)) THEN
                  WRITE (iw, '(T2,I8,I8)') i, 0
                  CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                         nset=nset, &
                                         npgf=npgf, &
                                         nshell=nshell, &
                                         l=l, &
                                         zet=zet, &
                                         gcc=gcc)

                  DO iset = 1, nset
                     DO ishell = 1, nshell(iset)
                        lshell = l(ishell, iset)
                        IF (lshell <= molden_lmax) THEN
                           WRITE (UNIT=iw, &
                                  FMT="(T25,A2,4X,I4,4X,F4.2)") &
                              angmom(lshell + 1:lshell + 1), npgf(iset), 1.0
                           ! MOLDEN expects the contraction coefficient of spherical NOT CARTESIAN NORMALISED
                           ! functions. So we undo the normalisation factors included in the gccs
                           ! Reverse engineered from basis_set_types, normalise_gcc_orb
                           prefac = 2_dp**lshell*(2/pi)**0.75_dp
                           expzet = 0.25_dp*(2*lshell + 3.0_dp)
                           WRITE (UNIT=iw, &
                                  FMT="((T51,2F15.6))") &
                              (zet(ipgf, iset), gcc(ipgf, ishell, iset)/(prefac*zet(ipgf, iset)**expzet), &
                               ipgf=1, npgf(iset))
                        ELSE
                           IF (print_warn) THEN
                              CALL cp_warn(__LOCATION__, &
                                           "MOLDEN format does not support Gaussian orbitals with l > 4.")
                              print_warn = .FALSE.
                           END IF
                        END IF
                     END DO
                  END DO

                  WRITE (iw, '(A4)') "    "

               END IF

            END DO

            IF (gto_kind == gto_spherical) THEN
               WRITE (iw, '(T2,A)') "[5D7F]"
               WRITE (iw, '(T2,A)') "[9G]"
            END IF

            WRITE (iw, '(T2,A)') "[MO]"
         END IF

         !------------------------------------------------------------------------
         ! convert from CP2K to MOLDEN format ordering
         ! http://www.cmbi.ru.nl/molden/molden_format.html
         !"The following order of D, F and G functions is expected:
         !
         !   5D: D 0, D+1, D-1, D+2, D-2
         !   6D: xx, yy, zz, xy, xz, yz
         !
         !   7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
         !  10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
         !
         !   9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
         !  15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
         !       xxyy xxzz yyzz xxyz yyxz zzxy
         !"
         ! CP2K has x in the outer (slower loop), so
         ! xx, xy, xz, yy, yz,zz for l=2, for instance
         !
         ! iorb_cp2k = orbmap(iorb_molden, l), l = 0 .. 4
         ! -----------------------------------------------------------------------
         IF (iw > 0) THEN
            IF (gto_kind == gto_cartesian) THEN
               ! -----------------------------------------------------------------
               ! Use cartesian (6D, 10F, 15G) representation.
               ! This is only format VMD can process.
               ! -----------------------------------------------------------------
               orbmap = RESHAPE((/1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                                  1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                                  1, 4, 6, 2, 3, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                                  1, 7, 10, 4, 2, 3, 6, 9, 8, 5, 0, 0, 0, 0, 0, &
                                  1, 11, 15, 2, 3, 7, 12, 10, 14, 4, 6, 13, 5, 8, 9/), &
                                (/molden_ncomax, molden_lmax + 1/))
            ELSE IF (gto_kind == gto_spherical) THEN
               ! -----------------------------------------------------------------
               ! Use spherical (5D, 7F, 9G) representation.
               ! -----------------------------------------------------------------
               orbmap = RESHAPE((/1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                                  3, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                                  3, 4, 2, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                                  4, 5, 3, 6, 2, 7, 1, 0, 0, 0, 0, 0, 0, 0, 0, &
                                  5, 6, 4, 7, 3, 8, 2, 9, 1, 0, 0, 0, 0, 0, 0/), &
                                (/molden_ncomax, molden_lmax + 1/))
            END IF
         END IF

         DO ispin = 1, SIZE(mos)
            CALL cp_fm_get_info(mos(ispin)%mo_set%mo_coeff, &
                                nrow_global=nrow_global, &
                                ncol_global=ncol_global)
            ALLOCATE (smatrix(nrow_global, ncol_global))
            CALL cp_fm_get_submatrix(mos(ispin)%mo_set%mo_coeff, smatrix)

            IF (iw > 0) THEN
               IF (gto_kind == gto_cartesian) THEN
                  CALL get_qs_kind_set(qs_kind_set, ncgf=ncgf, nsgf=nsgf)

                  ALLOCATE (cmatrix(ncgf, ncgf))

                  cmatrix = 0.0_dp

                  ! Transform spherical MOs to Cartesian MOs

                  icgf = 1
                  isgf = 1
                  DO iatom = 1, SIZE(particle_set)
                     NULLIFY (orb_basis_set)
                     CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
                     CALL get_qs_kind(qs_kind_set(ikind), &
                                      basis_set=orb_basis_set)
                     IF (ASSOCIATED(orb_basis_set)) THEN
                        CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                               nset=nset, &
                                               nshell=nshell, &
                                               l=l)
                        DO iset = 1, nset
                           DO ishell = 1, nshell(iset)
                              lshell = l(ishell, iset)
                              CALL dgemm("T", "N", nco(lshell), mos(ispin)%mo_set%nmo, nso(lshell), 1.0_dp, &
                                         orbtramat(lshell)%s2c, nso(lshell), &
                                         smatrix(isgf, 1), nsgf, 0.0_dp, &
                                         cmatrix(icgf, 1), ncgf)
                              icgf = icgf + nco(lshell)
                              isgf = isgf + nso(lshell)
                           END DO
                        END DO
                     END IF
                  END DO ! iatom
               END IF

               DO icol = 1, mos(ispin)%mo_set%nmo
                  ! index of the first basis function for the given atom, set, and shell
                  irow = 1

                  ! index of the first basis function in MOLDEN file.
                  ! Due to limitation of the MOLDEN format, basis functions with l > molden_lmax
                  ! cannot be exported, so we need to renumber atomic orbitals
                  irow_in = 1

                  WRITE (iw, '(A,ES20.10)') 'Ene=', mos(ispin)%mo_set%eigenvalues(icol)
                  IF (ispin < 2) THEN
                     WRITE (iw, '(A)') 'Spin= Alpha'
                  ELSE
                     WRITE (iw, '(A)') 'Spin= Beta'
                  END IF
                  WRITE (iw, '(A,F12.7)') 'Occup=', mos(ispin)%mo_set%occupation_numbers(icol)

                  DO iatom = 1, SIZE(particle_set)
                     NULLIFY (orb_basis_set)
                     CALL get_atomic_kind(particle_set(iatom)%atomic_kind, &
                                          element_symbol=element_symbol, kind_number=ikind)
                     CALL get_qs_kind(qs_kind_set(ikind), &
                                      basis_set=orb_basis_set)
                     IF (ASSOCIATED(orb_basis_set)) THEN
                        CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                               nset=nset, &
                                               nshell=nshell, &
                                               l=l)

                        IF (gto_kind == gto_cartesian) THEN
                           ! ----------------------------------------------
                           ! Use cartesian (6D, 10F, 15G) representation.
                           ! ----------------------------------------------
                           icgf = 1
                           DO iset = 1, nset
                              DO ishell = 1, nshell(iset)
                                 lshell = l(ishell, iset)

                                 IF (lshell <= molden_lmax) THEN
                                    CALL print_coeffs(iw, fmtstr1, ndigits, irow_in, orbmap(:, lshell), &
                                                      cmatrix(irow:irow + nco(lshell) - 1, icol))
                                    irow_in = irow_in + nco(lshell)
                                 END IF

                                 irow = irow + nco(lshell)
                              END DO ! ishell
                           END DO

                        ELSE IF (gto_kind == gto_spherical) THEN
                           ! ----------------------------------------------
                           ! Use spherical (5D, 7F, 9G) representation.
                           ! ----------------------------------------------
                           DO iset = 1, nset
                              DO ishell = 1, nshell(iset)
                                 lshell = l(ishell, iset)

                                 IF (lshell <= molden_lmax) THEN
                                    CALL print_coeffs(iw, fmtstr1, ndigits, irow_in, orbmap(:, lshell), &
                                                      smatrix(irow:irow + nso(lshell) - 1, icol))
                                    irow_in = irow_in + nso(lshell)
                                 END IF

                                 irow = irow + nso(lshell)
                              END DO
                           END DO
                        END IF

                     END IF
                  END DO ! iatom
               END DO
            END IF

            IF (ALLOCATED(cmatrix)) DEALLOCATE (cmatrix)
            IF (ALLOCATED(smatrix)) DEALLOCATE (smatrix)
         END DO

         CALL cp_print_key_finished_output(iw, logger, print_section, "")

      END IF

      CALL timestop(handle)

   END SUBROUTINE write_mos_molden

! **************************************************************************************************
!> \brief Output MO coefficients formatted correctly for MOLDEN, omitting those <= 1E(-digits)
!> \param iw       output file unit
!> \param fmtstr1  format string
!> \param ndigits  number of significant digits in MO coefficients
!> \param irow_in  index of the first atomic orbital: mo_coeff(orbmap(1))
!> \param orbmap   array to map Gaussian functions from MOLDEN to CP2K ordering
!> \param mo_coeff MO coefficients
! **************************************************************************************************
   SUBROUTINE print_coeffs(iw, fmtstr1, ndigits, irow_in, orbmap, mo_coeff)
      INTEGER, INTENT(in)                                :: iw
      CHARACTER(LEN=*), INTENT(in)                       :: fmtstr1
      INTEGER, INTENT(in)                                :: ndigits, irow_in
      INTEGER, DIMENSION(molden_ncomax), INTENT(in)      :: orbmap
      REAL(KIND=dp), DIMENSION(:), INTENT(in)            :: mo_coeff

      INTEGER                                            :: orbital

      DO orbital = 1, molden_ncomax
         IF (orbmap(orbital) /= 0) THEN
            IF (ABS(mo_coeff(orbmap(orbital))) >= 10.0_dp**(-ndigits)) THEN
               WRITE (iw, fmtstr1) irow_in + orbital - 1, mo_coeff(orbmap(orbital))
            END IF
         END IF
      END DO

   END SUBROUTINE print_coeffs

! **************************************************************************************************
!> \brief writes the output for vibrational analysis in MOLDEN format
!> \param input ...
!> \param particles ...
!> \param freq ...
!> \param eigen_vec ...
!> \param intensities ...
!> \param calc_intens ...
!> \param dump_only_positive ...
!> \param logger ...
!> \param list array of mobile atom indices
!> \author Florian Schiffmann 11.2007
! **************************************************************************************************
   SUBROUTINE write_vibrations_molden(input, particles, freq, eigen_vec, intensities, calc_intens, &
                                      dump_only_positive, logger, list)

      TYPE(section_vals_type), POINTER                   :: input
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      REAL(KIND=dp), DIMENSION(:)                        :: freq
      REAL(KIND=dp), DIMENSION(:, :)                     :: eigen_vec
      REAL(KIND=dp), DIMENSION(:), POINTER               :: intensities
      LOGICAL, INTENT(in)                                :: calc_intens, dump_only_positive
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: list

      CHARACTER(len=*), PARAMETER :: routineN = 'write_vibrations_molden'

      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: handle, i, iw, j, k, l, z
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_list

      CALL timeset(routineN, handle)

      iw = cp_print_key_unit_nr(logger, input, "VIBRATIONAL_ANALYSIS%PRINT%MOLDEN_VIB", &
                                extension=".mol", file_status='REPLACE')

      IF (iw .GT. 0) THEN
         CPASSERT(MOD(SIZE(eigen_vec, 1), 3) == 0)
         CPASSERT(SIZE(freq, 1) == SIZE(eigen_vec, 2))
         ALLOCATE (my_list(SIZE(particles)))
         ! Either we have a list of the subset of mobile atoms,
         ! Or the eigenvectors must span the full space (all atoms)
         IF (PRESENT(list)) THEN
            my_list(:) = 0
            DO i = 1, SIZE(list)
               my_list(list(i)) = i
            END DO
         ELSE
            CPASSERT(SIZE(particles) == SIZE(eigen_vec, 1)/3)
            DO i = 1, SIZE(my_list)
               my_list(i) = i
            END DO
         END IF
         WRITE (iw, '(T2,A)') "[Molden Format]"
         WRITE (iw, '(T2,A)') "[Atoms] AU"
         DO i = 1, SIZE(particles)
            CALL get_atomic_kind(atomic_kind=particles(i)%atomic_kind, &
                                 element_symbol=element_symbol)
            CALL get_ptable_info(element_symbol, number=z)

            WRITE (iw, '(T2,A2,I8,I8,3X,3(F12.6,3X))') &
               element_symbol, i, z, particles(i)%r(:)

         END DO
         WRITE (iw, '(T2,A)') "[FREQ]"
         DO i = 1, SIZE(freq, 1)
            IF ((.NOT. dump_only_positive) .OR. (freq(i) >= 0._dp)) WRITE (iw, '(T5,F12.6)') freq(i)
         END DO
         WRITE (iw, '(T2,A)') "[FR-COORD]"
         DO i = 1, SIZE(particles)
            CALL get_atomic_kind(atomic_kind=particles(i)%atomic_kind, &
                                 element_symbol=element_symbol)
            WRITE (iw, '(T2,A2,3X,3(F12.6,3X))') &
               element_symbol, particles(i)%r(:)
         END DO
         WRITE (iw, '(T2,A)') "[FR-NORM-COORD]"
         l = 0
         DO i = 1, SIZE(eigen_vec, 2)
            IF ((.NOT. dump_only_positive) .OR. (freq(i) >= 0._dp)) THEN
               l = l + 1
               WRITE (iw, '(T2,A,1X,I6)') "vibration", l
               DO j = 1, SIZE(particles)
                  IF (my_list(j) .NE. 0) THEN
                     k = (my_list(j) - 1)*3
                     WRITE (iw, '(T2,3(F12.6,3X))') eigen_vec(k + 1, i), eigen_vec(k + 2, i), eigen_vec(k + 3, i)
                  ELSE
                     WRITE (iw, '(T2,3(F12.6,3X))') 0.0_dp, 0.0_dp, 0.0_dp
                  END IF
               END DO
            END IF
         END DO
         IF (calc_intens) THEN
            WRITE (iw, '(T2,A)') "[INT]"
            DO i = 1, SIZE(intensities)
               IF ((.NOT. dump_only_positive) .OR. (freq(i) >= 0._dp)) WRITE (iw, '(3X,F18.6)') intensities(i)
            END DO
         END IF
         DEALLOCATE (my_list)
      END IF
      CALL cp_print_key_finished_output(iw, logger, input, "VIBRATIONAL_ANALYSIS%PRINT%MOLDEN_VIB")

      CALL timestop(handle)

   END SUBROUTINE write_vibrations_molden

END MODULE molden_utils
