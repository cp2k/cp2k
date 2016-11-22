!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
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
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: nco,&
                                              nso
   USE orbital_symbols,                 ONLY: cgf_symbol
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

   PUBLIC :: write_vibrations_molden, write_mos_molden

CONTAINS

! **************************************************************************************************
!> \brief Write out the MOs in molden format for visualisation
!> \param mos the set of MOs (both spins, if UKS)
!> \param qs_kind_set for basis set info
!> \param particle_set particles data structure, for positions and kinds
!> \param dft_section input section containing relevant print key
!> \author MattW, IainB
! **************************************************************************************************
   SUBROUTINE write_mos_molden(mos, qs_kind_set, particle_set, dft_section)

      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(LEN=1)                                   :: angmom
      CHARACTER(LEN=12), DIMENSION(15)                   :: shell_symbol
      CHARACTER(LEN=12), DIMENSION(:), POINTER           :: bcgf_symbol
      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER :: i, iatom, icgf, ico, icol, ikind, ipgf, irow, irow_in, iset, isgf, ishell, ispin, &
         iw, lshell, ncgf, ncol_global, nrow_global, nset, nsgf, shellgf, z
      INTEGER, DIMENSION(:), POINTER                     :: nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: l
      REAL(KIND=dp)                                      :: expzet, prefac
      REAL(KIND=dp), DIMENSION(15)                       :: mo_coeff
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cmatrix, smatrix
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      NULLIFY (nshell)
      NULLIFY (bcgf_symbol)

      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger, dft_section, "SCF%PRINT%MOS_MOLDEN", &
                                extension=".molden", file_status='REPLACE')

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

      IF (iw .GT. 0) THEN
         WRITE (iw, '(T2,A)') "[Molden Format]"
         WRITE (iw, '(T2,A)') "[Atoms] AU"
         DO i = 1, SIZE(particle_set)
            CALL get_atomic_kind(atomic_kind=particle_set(i)%atomic_kind, &
                                 element_symbol=element_symbol)
            CALL get_ptable_info(element_symbol, number=z)

            WRITE (iw, '(T2,A2,I8,I8,3X,3(F12.6,3X))') &
               element_symbol, i, z, particle_set(i)%r(:)

         END DO

         WRITE (iw, '(T2,A)') "[GTO] "

         DO i = 1, SIZE(particle_set)
            CALL get_atomic_kind(atomic_kind=particle_set(i)%atomic_kind, kind_number=ikind, &
                                 element_symbol=element_symbol)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
            IF (ASSOCIATED(orb_basis_set)) THEN
               WRITE (iw, '(T2,I8,I8)') i, 0
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                      nset=nset, &
                                      nshell=nshell, &
                                      l=l)

               DO iset = 1, orb_basis_set%nset
                  DO ishell = 1, orb_basis_set%nshell(iset)
                     IF (orb_basis_set%l(ishell, iset) == 0) angmom = 's'
                     IF (orb_basis_set%l(ishell, iset) == 1) angmom = 'p'
                     IF (orb_basis_set%l(ishell, iset) == 2) angmom = 'd'
                     IF (orb_basis_set%l(ishell, iset) == 3) angmom = 'f'
                     IF (orb_basis_set%l(ishell, iset) == 4) angmom = 'g'
                     WRITE (UNIT=iw, &
                            FMT="(T25,A2,4X,I4,4X,F4.2)") &
                        angmom, &
                        orb_basis_set%npgf(iset), 1.0
                     ! MOLDEN expects the contaction coefficient of spherical NOT CARTESIAN NORMALISED
                     ! functions. So we undo the normalisation factors included in the gccs
                     ! Reverse engineered from basis_set_types, normalise_gcc_orb
                     prefac = 2_dp**orb_basis_set%l(ishell, iset)*(2/pi)**0.75_dp
                     expzet = 0.25_dp*(2*orb_basis_set%l(ishell, iset)+3.0_dp)
                     WRITE (UNIT=iw, &
                            FMT="((T51,2F15.6))") &
                        (orb_basis_set%zet(ipgf, iset), &
                         orb_basis_set%gcc(ipgf, ishell, iset)/(prefac*orb_basis_set%zet(ipgf, iset)**expzet), &
                         ipgf=1, orb_basis_set%npgf(iset))
                  END DO
               END DO

               WRITE (iw, '(A4)') "    "

            END IF

         END DO

         WRITE (iw, '(T2,A)') "[MO] "

      END IF

      DO ispin = 1, SIZE(mos)
         CALL cp_fm_get_info(mos(ispin)%mo_set%mo_coeff, &
                             nrow_global=nrow_global, &
                             ncol_global=ncol_global)
         ALLOCATE (smatrix(nrow_global, ncol_global))
         CALL cp_fm_get_submatrix(mos(ispin)%mo_set%mo_coeff, smatrix)
         IF (.NOT. iw > 0) THEN
            DEALLOCATE (smatrix)
         END IF

         IF (iw > 0) THEN
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
                        icgf = icgf+nco(lshell)
                        isgf = isgf+nso(lshell)
                     END DO
                  END DO
               END IF
            END DO ! iatom

            DO icol = 1, mos(ispin)%mo_set%nmo

               irow = 1
               IF (iw > 0) WRITE (iw, *) 'Ene=', mos(ispin)%mo_set%eigenvalues(icol)
               IF (ispin < 2) THEN
                  IF (iw > 0) WRITE (iw, *) 'Spin= Alpha'
               ELSE
                  IF (iw > 0) WRITE (iw, *) 'Spin= Beta'
               END IF
               IF (iw > 0) WRITE (iw, *) 'Occup=', mos(ispin)%mo_set%occupation_numbers(icol)

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
                                            l=l, &
                                            cgf_symbol=bcgf_symbol)
                     icgf = 1
                     DO iset = 1, nset
                        DO ishell = 1, nshell(iset)
                           lshell = l(ishell, iset)
                           irow_in = irow
                           shellgf = 1
                           DO ico = 1, nco(lshell)
                              mo_coeff(shellgf) = cmatrix(irow, icol)
                              shell_symbol(shellgf) = bcgf_symbol(icgf)
                              icgf = icgf+1
                              shellgf = shellgf+1
                              irow = irow+1
                           END DO ! ico
                           !------------------------------------------------------------------------
                           ! convert from CP2K MOLDEN format ordering
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
                           ! Always use cartesian (6D, 10F, 15G) as this is only format
                           ! VMD can process anyway
                           !XXX can't be bother making this more elegant with mask etc
                           ! -----------------------------------------------------------------------
                           IF (lshell .EQ. 0) THEN
                              WRITE (iw, *) irow_in+0, mo_coeff(1) !, shell_symbol(1)
                           ELSE IF (lshell .EQ. 1) THEN
                              WRITE (iw, *) irow_in+0, mo_coeff(1) !, shell_symbol(1)
                              WRITE (iw, *) irow_in+1, mo_coeff(2) !, shell_symbol(2)
                              WRITE (iw, *) irow_in+2, mo_coeff(3) !, shell_symbol(3)
                           ELSE IF (lshell .EQ. 2) THEN
                              WRITE (iw, *) irow_in+0, mo_coeff(1) !, shell_symbol(1)
                              WRITE (iw, *) irow_in+1, mo_coeff(4) !, shell_symbol(4)
                              WRITE (iw, *) irow_in+2, mo_coeff(6) !, shell_symbol(6)
                              WRITE (iw, *) irow_in+3, mo_coeff(2) !, shell_symbol(2)
                              WRITE (iw, *) irow_in+4, mo_coeff(3) !, shell_symbol(3)
                              WRITE (iw, *) irow_in+5, mo_coeff(5) !, shell_symbol(5)
                           ELSE IF (lshell .EQ. 3) THEN
                              WRITE (iw, *) irow_in+0, mo_coeff(1) !, shell_symbol(1)
                              WRITE (iw, *) irow_in+1, mo_coeff(7) !, shell_symbol(7)
                              WRITE (iw, *) irow_in+2, mo_coeff(10) !, shell_symbol(10)
                              WRITE (iw, *) irow_in+3, mo_coeff(4) !, shell_symbol(4)
                              WRITE (iw, *) irow_in+4, mo_coeff(2) !, shell_symbol(2)
                              WRITE (iw, *) irow_in+5, mo_coeff(3) !, shell_symbol(3)
                              WRITE (iw, *) irow_in+6, mo_coeff(6) !, shell_symbol(6)
                              WRITE (iw, *) irow_in+7, mo_coeff(9) !, shell_symbol(9)
                              WRITE (iw, *) irow_in+8, mo_coeff(8) !, shell_symbol(8)
                              WRITE (iw, *) irow_in+9, mo_coeff(5) !, shell_symbol(5)
                           ELSE IF (lshell .EQ. 4) THEN
                              WRITE (iw, *) irow_in+0, mo_coeff(1) !, shell_symbol(1)
                              WRITE (iw, *) irow_in+1, mo_coeff(11) !, shell_symbol(11)
                              WRITE (iw, *) irow_in+2, mo_coeff(15) !, shell_symbol(15)
                              WRITE (iw, *) irow_in+3, mo_coeff(2) !, shell_symbol(2)
                              WRITE (iw, *) irow_in+4, mo_coeff(3) !, shell_symbol(3)
                              WRITE (iw, *) irow_in+5, mo_coeff(7) !, shell_symbol(7)
                              WRITE (iw, *) irow_in+6, mo_coeff(12) !, shell_symbol(12)
                              WRITE (iw, *) irow_in+7, mo_coeff(10) !, shell_symbol(10)
                              WRITE (iw, *) irow_in+8, mo_coeff(14) !, shell_symbol(14)
                              WRITE (iw, *) irow_in+9, mo_coeff(4) !, shell_symbol(4)
                              WRITE (iw, *) irow_in+10, mo_coeff(6) !, shell_symbol(6)
                              WRITE (iw, *) irow_in+11, mo_coeff(13) !, shell_symbol(13)
                              WRITE (iw, *) irow_in+12, mo_coeff(5) !, shell_symbol(5)
                              WRITE (iw, *) irow_in+13, mo_coeff(8) !, shell_symbol(8)
                              WRITE (iw, *) irow_in+14, mo_coeff(9) !, shell_symbol(9)
                           ELSE
                              ! ignore larger than g - VMD doesn't render anyway
                              ! and MOLDEN format not defined
                              CONTINUE
                           END IF
                        END DO ! ishell
                     END DO
                  END IF
               END DO ! iatom
            END DO
            DEALLOCATE (smatrix, cmatrix)
         END IF
      END DO

   END SUBROUTINE

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
      LOGICAL                                            :: calc_intens, dump_only_positive
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: list

      CHARACTER(len=*), PARAMETER :: routineN = 'write_vibrations_molden', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: i, iw, j, k, l
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_list

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
               l = l+1
               WRITE (iw, '(T2,A,1X,I6)') "vibration", l
               DO j = 1, SIZE(particles)
                  IF (my_list(j) .NE. 0) THEN
                     k = (my_list(j)-1)*3
                     WRITE (iw, '(T2,3(F12.6,3X))') eigen_vec(k+1, i), eigen_vec(k+2, i), eigen_vec(k+3, i)
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
   END SUBROUTINE write_vibrations_molden

END MODULE molden_utils
