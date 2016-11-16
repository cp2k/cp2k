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
   USE cp_fm_types,                     ONLY: cp_fm_get_element
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: get_ptable_info
   USE qs_kind_types,                   ONLY: get_qs_kind,&
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
      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: i, iao, ikind, imo, ipgf, iset, ishell, &
                                                            ispin, iw, nmo, nset, z
      INTEGER, DIMENSION(:), POINTER                     :: nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: l
      REAL(KIND=dp)                                      :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

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
                     WRITE (UNIT=iw, &
                            FMT="((T51,2F15.6))") &
                        (orb_basis_set%zet(ipgf, iset), &
                         orb_basis_set%gcc(ipgf, ishell, iset), &
                         ipgf=1, orb_basis_set%npgf(iset))
                  END DO
               END DO

               WRITE (iw, '(A4)') "    "

            END IF

         END DO

         WRITE (iw, '(T2,A)') "[MO] "

      END IF

      DO ispin = 1, SIZE(mos)
         nmo = mos(ispin)%mo_set%nmo
         DO imo = 1, nmo
            IF (iw > 0) WRITE (iw, *) 'Ene=', mos(ispin)%mo_set%eigenvalues(imo)
            IF (ispin < 2) THEN
               IF (iw > 0) WRITE (iw, *) 'Spin= Alpha'
            ELSE
               IF (iw > 0) WRITE (iw, *) 'Spin= Beta'
            END IF
            IF (iw > 0) WRITE (iw, *) 'Occup=', mos(ispin)%mo_set%occupation_numbers(imo)
            DO iao = 1, mos(ispin)%mo_set%nao
               CALL cp_fm_get_element(mos(ispin)%mo_set%mo_coeff, iao, imo, mo_coeff)
               !Hardcoded threshold for compressing
               IF (iw > 0 .AND. mo_coeff > 1.0E-6) WRITE (iw, '(I6, 2X, ES14.7)') iao, mo_coeff
            END DO
         END DO
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
               element_symbol, particles((i))%r(:)
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
