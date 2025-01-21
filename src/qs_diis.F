!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Apply the direct inversion in the iterative subspace (DIIS) of Pulay
!>      in the framework of an SCF iteration for convergence acceleration
!> \par Literature
!>      - P. Pulay, Chem. Phys. Lett. 73, 393 (1980)
!>      - P. Pulay, J. Comput. Chem. 3, 556 (1982)
!> \par History
!>      - Changed to BLACS matrix usage (08.06.2001,MK)
!>      - rewritten to include LSD (1st attempt) (01.2003, Joost VandeVondele)
!>      - DIIS for ROKS (05.04.06,MK)
!>      - DIIS for k-points (04.2023, Augustin Bussy)
!> \author Matthias Krack (28.06.2000)
! **************************************************************************************************
MODULE qs_diis
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_trace
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_dot, dbcsr_maxabs, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_release, dbcsr_set, dbcsr_transposed, dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_symm,&
                                              cp_fm_trace
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_maxabsval,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathlib,                         ONLY: diag_complex,&
                                              diamat_all
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_diis_types,                   ONLY: qs_diis_buffer_type,&
                                              qs_diis_buffer_type_kp,&
                                              qs_diis_buffer_type_sparse
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE string_utilities,                ONLY: compress
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_diis'

   ! Public subroutines

   PUBLIC :: qs_diis_b_clear, &
             qs_diis_b_create, &
             qs_diis_b_step
   PUBLIC :: qs_diis_b_clear_sparse, &
             qs_diis_b_create_sparse, &
             qs_diis_b_step_4lscf
   PUBLIC :: qs_diis_b_clear_kp, &
             qs_diis_b_create_kp, &
             qs_diis_b_step_kp, &
             qs_diis_b_calc_err_kp, &
             qs_diis_b_info_kp

CONTAINS

! **************************************************************************************************
!> \brief Allocates an SCF DIIS buffer
!> \param diis_buffer the buffer to create
!> \param nbuffer ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE qs_diis_b_create(diis_buffer, nbuffer)

      TYPE(qs_diis_buffer_type), INTENT(OUT)             :: diis_buffer
      INTEGER, INTENT(in)                                :: nbuffer

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_diis_b_create'

      INTEGER                                            :: handle

! -------------------------------------------------------------------------

      CALL timeset(routineN, handle)

      NULLIFY (diis_buffer%b_matrix)
      NULLIFY (diis_buffer%error)
      NULLIFY (diis_buffer%param)
      diis_buffer%nbuffer = nbuffer
      diis_buffer%ncall = 0

      CALL timestop(handle)

   END SUBROUTINE qs_diis_b_create

! **************************************************************************************************
!> \brief Allocate and initialize a DIIS buffer for nao*nao parameter
!>      variables and with a buffer size of nbuffer.
!> \param diis_buffer the buffer to initialize
!> \param matrix_struct the structure for the matrix of the buffer
!> \param nspin ...
!> \param scf_section ...
!> \par History
!>      - Creation (07.05.2001, Matthias Krack)
!>      - Changed to BLACS matrix usage (08.06.2001,MK)
!>      - DIIS for ROKS (05.04.06,MK)
!> \author Matthias Krack
!> \note
!>      check to allocate matrixes only when needed, using a linked list?
! **************************************************************************************************
   SUBROUTINE qs_diis_b_check_i_alloc(diis_buffer, matrix_struct, nspin, &
                                      scf_section)

      TYPE(qs_diis_buffer_type), INTENT(INOUT)           :: diis_buffer
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      INTEGER, INTENT(IN)                                :: nspin
      TYPE(section_vals_type), POINTER                   :: scf_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_diis_b_check_i_alloc'

      INTEGER                                            :: handle, ibuffer, ispin, nbuffer, &
                                                            output_unit
      TYPE(cp_logger_type), POINTER                      :: logger

! -------------------------------------------------------------------------

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()

      nbuffer = diis_buffer%nbuffer

      IF (.NOT. ASSOCIATED(diis_buffer%error)) THEN
         ALLOCATE (diis_buffer%error(nbuffer, nspin))

         DO ispin = 1, nspin
            DO ibuffer = 1, nbuffer
               CALL cp_fm_create(diis_buffer%error(ibuffer, ispin), &
                                 name="qs_diis_b%error("// &
                                 TRIM(ADJUSTL(cp_to_string(ibuffer)))//","// &
                                 TRIM(ADJUSTL(cp_to_string(ibuffer)))//")", &
                                 matrix_struct=matrix_struct)
            END DO
         END DO
      END IF

      IF (.NOT. ASSOCIATED(diis_buffer%param)) THEN
         ALLOCATE (diis_buffer%param(nbuffer, nspin))

         DO ispin = 1, nspin
            DO ibuffer = 1, nbuffer
               CALL cp_fm_create(diis_buffer%param(ibuffer, ispin), &
                                 name="qs_diis_b%param("// &
                                 TRIM(ADJUSTL(cp_to_string(ibuffer)))//","// &
                                 TRIM(ADJUSTL(cp_to_string(ibuffer)))//")", &
                                 matrix_struct=matrix_struct)
            END DO
         END DO
      END IF

      IF (.NOT. ASSOCIATED(diis_buffer%b_matrix)) THEN
         ALLOCATE (diis_buffer%b_matrix(nbuffer + 1, nbuffer + 1))
         diis_buffer%b_matrix = 0.0_dp
         output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DIIS_INFO", &
                                            extension=".scfLog")
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T9,A)") &
               "DIIS | The SCF DIIS buffer was allocated and initialized"
         END IF
         CALL cp_print_key_finished_output(output_unit, logger, scf_section, &
                                           "PRINT%DIIS_INFO")
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_diis_b_check_i_alloc

! **************************************************************************************************
!> \brief Update the SCF DIIS buffer, and if appropriate does a diis step.
!> \param diis_buffer ...
!> \param mo_array ...
!> \param kc ...
!> \param sc ...
!> \param delta ...
!> \param error_max ...
!> \param diis_step ...
!> \param eps_diis ...
!> \param nmixing ...
!> \param s_matrix ...
!> \param scf_section ...
!> \param roks ...
!> \par History
!>      - Creation (07.05.2001, Matthias Krack)
!>      - Changed to BLACS matrix usage (08.06.2001, MK)
!>      - 03.2003 rewamped [fawzi]
!>      - Adapted for high-spin ROKS (08.04.06,MK)
!> \author Matthias Krack
! **************************************************************************************************
   SUBROUTINE qs_diis_b_step(diis_buffer, mo_array, kc, sc, delta, error_max, &
                             diis_step, eps_diis, nmixing, s_matrix, scf_section, roks)

      TYPE(qs_diis_buffer_type), POINTER                 :: diis_buffer
      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mo_array
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: kc
      TYPE(cp_fm_type), INTENT(IN)                       :: sc
      REAL(KIND=dp), INTENT(IN)                          :: delta
      REAL(KIND=dp), INTENT(OUT)                         :: error_max
      LOGICAL, INTENT(OUT)                               :: diis_step
      REAL(KIND=dp), INTENT(IN)                          :: eps_diis
      INTEGER, INTENT(IN), OPTIONAL                      :: nmixing
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: s_matrix
      TYPE(section_vals_type), POINTER                   :: scf_section
      LOGICAL, INTENT(IN), OPTIONAL                      :: roks

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_diis_b_step'
      REAL(KIND=dp), PARAMETER :: eigenvalue_threshold = 1.0E-12_dp

      CHARACTER(LEN=2*default_string_length)             :: message
      INTEGER                                            :: handle, homo, ib, imo, ispin, jb, &
                                                            my_nmixing, nao, nb, nb1, nmo, nspin, &
                                                            output_unit
      LOGICAL                                            :: eigenvectors_discarded, my_roks
      REAL(KIND=dp)                                      :: maxocc, tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ev, occ
      REAL(KIND=dp), DIMENSION(:), POINTER               :: occa, occb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: a, b
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), POINTER                          :: c, new_errors, old_errors, parameters
      TYPE(cp_logger_type), POINTER                      :: logger

! -------------------------------------------------------------------------

      CALL timeset(routineN, handle)

      nspin = SIZE(mo_array)
      diis_step = .FALSE.

      IF (PRESENT(roks)) THEN
         my_roks = .TRUE.
         nspin = 1
      ELSE
         my_roks = .FALSE.
      END IF

      my_nmixing = 2
      IF (PRESENT(nmixing)) my_nmixing = nmixing

      NULLIFY (c, new_errors, old_errors, parameters, matrix_struct, a, b, occa, occb)
      logger => cp_get_default_logger()

      ! Quick return, if no DIIS is requested

      IF (diis_buffer%nbuffer < 1) THEN
         CALL timestop(handle)
         RETURN
      END IF

      CALL cp_fm_get_info(kc(1), &
                          matrix_struct=matrix_struct)
      CALL qs_diis_b_check_i_alloc(diis_buffer, &
                                   matrix_struct=matrix_struct, &
                                   nspin=nspin, &
                                   scf_section=scf_section)

      error_max = 0.0_dp

      ib = MODULO(diis_buffer%ncall, diis_buffer%nbuffer) + 1
      diis_buffer%ncall = diis_buffer%ncall + 1
      nb = MIN(diis_buffer%ncall, diis_buffer%nbuffer)

      DO ispin = 1, nspin

         CALL get_mo_set(mo_set=mo_array(ispin), &
                         nao=nao, &
                         nmo=nmo, &
                         homo=homo, &
                         mo_coeff=c, &
                         occupation_numbers=occa, &
                         maxocc=maxocc)

         new_errors => diis_buffer%error(ib, ispin)
         parameters => diis_buffer%param(ib, ispin)

         ! Copy the Kohn-Sham matrix K to the DIIS buffer

         CALL cp_fm_to_fm(kc(ispin), parameters)

         IF (my_roks) THEN

            ALLOCATE (occ(nmo))

            CALL get_mo_set(mo_set=mo_array(2), &
                            occupation_numbers=occb)

            DO imo = 1, nmo
               occ(imo) = SQRT(occa(imo) + occb(imo))
            END DO

            CALL cp_fm_to_fm(c, sc)
            CALL cp_fm_column_scale(sc, occ(1:homo))

            ! KC <- K*C
            CALL cp_fm_symm("L", "U", nao, homo, 1.0_dp, parameters, sc, 0.0_dp, kc(ispin))

            IF (PRESENT(s_matrix)) THEN
               CALL copy_dbcsr_to_fm(s_matrix(1)%matrix, new_errors)
               ! SC <- S*C
               CALL cp_fm_symm("L", "U", nao, homo, 1.0_dp, new_errors, c, 0.0_dp, sc)
               CALL cp_fm_column_scale(sc, occ(1:homo))
            END IF

            ! new_errors <- KC*(SC)^T - (SC)*(KC)^T = K*P*S - S*P*K
            ! or for an orthogonal basis
            ! new_errors <- KC*C^T - C*(KC)^T = K*P - P*K with S = I
            CALL parallel_gemm("N", "T", nao, nao, homo, 1.0_dp, sc, kc(ispin), 0.0_dp, new_errors)
            CALL parallel_gemm("N", "T", nao, nao, homo, 1.0_dp, kc(ispin), sc, -1.0_dp, new_errors)

            DEALLOCATE (occ)

         ELSE

            ! KC <- K*C
            CALL cp_fm_symm("L", "U", nao, homo, maxocc, parameters, c, 0.0_dp, kc(ispin))

            IF (PRESENT(s_matrix)) THEN
               ! I guess that this copy can be avoided for LSD
               CALL copy_dbcsr_to_fm(s_matrix(1)%matrix, new_errors)
               ! sc <- S*C
               CALL cp_fm_symm("L", "U", nao, homo, 2.0_dp, new_errors, c, 0.0_dp, sc)
               ! new_errors <- KC*(SC)^T - (SC)*(KC)^T = K*P*S - S*P*K
               CALL parallel_gemm("N", "T", nao, nao, homo, 1.0_dp, sc, kc(ispin), 0.0_dp, new_errors)
               CALL parallel_gemm("N", "T", nao, nao, homo, 1.0_dp, kc(ispin), sc, -1.0_dp, new_errors)
            ELSE
               ! new_errors <- KC*(C)^T - C*(KC)^T = K*P - P*K
               CALL parallel_gemm("N", "T", nao, nao, homo, 1.0_dp, c, kc(ispin), 0.0_dp, new_errors)
               CALL parallel_gemm("N", "T", nao, nao, homo, 1.0_dp, kc(ispin), c, -1.0_dp, new_errors)
            END IF

         END IF

         CALL cp_fm_maxabsval(new_errors, tmp)
         error_max = MAX(error_max, tmp)

      END DO

      ! Check, if a DIIS step is appropriate

      diis_step = ((diis_buffer%ncall >= my_nmixing) .AND. (delta < eps_diis))

      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DIIS_INFO", &
                                         extension=".scfLog")
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T9,A,I4,/,(T9,A,ES12.3))") &
            "DIIS | Current SCF DIIS buffer size:         ", nb, &
            "DIIS | Maximum SCF DIIS error vector element:", error_max, &
            "DIIS | Current SCF convergence:              ", delta, &
            "DIIS | Threshold value for a DIIS step:      ", eps_diis
         IF (error_max < eps_diis) THEN
            WRITE (UNIT=output_unit, FMT="(T9,A)") &
               "DIIS | => The SCF DIIS buffer will be updated"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T9,A)") &
               "DIIS | => No update of the SCF DIIS buffer"
         END IF
         IF (diis_step .AND. (error_max < eps_diis)) THEN
            WRITE (UNIT=output_unit, FMT="(T9,A,/)") &
               "DIIS | => A SCF DIIS step will be performed"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T9,A,/)") &
               "DIIS | => No SCF DIIS step will be performed"
         END IF
      END IF

      ! Update the SCF DIIS buffer

      IF (error_max < eps_diis) THEN

         b => diis_buffer%b_matrix

         DO jb = 1, nb
            b(jb, ib) = 0.0_dp
            DO ispin = 1, nspin
               old_errors => diis_buffer%error(jb, ispin)
               new_errors => diis_buffer%error(ib, ispin)
               CALL cp_fm_trace(old_errors, new_errors, tmp)
               b(jb, ib) = b(jb, ib) + tmp
            END DO
            b(ib, jb) = b(jb, ib)
         END DO

      ELSE

         diis_step = .FALSE.

      END IF

      ! Perform DIIS step

      IF (diis_step) THEN

         nb1 = nb + 1

         ALLOCATE (a(nb1, nb1))
         ALLOCATE (b(nb1, nb1))
         ALLOCATE (ev(nb1))

         ! Set up the linear DIIS equation system

         b(1:nb, 1:nb) = diis_buffer%b_matrix(1:nb, 1:nb)

         b(1:nb, nb1) = -1.0_dp
         b(nb1, 1:nb) = -1.0_dp
         b(nb1, nb1) = 0.0_dp

         ! Solve the linear DIIS equation system

         ev(1:nb1) = 0.0_dp
         CALL diamat_all(b(1:nb1, 1:nb1), ev(1:nb1))

         a(1:nb1, 1:nb1) = b(1:nb1, 1:nb1)

         eigenvectors_discarded = .FALSE.

         DO jb = 1, nb1
            IF (ABS(ev(jb)) < eigenvalue_threshold) THEN
               IF (output_unit > 0) THEN
                  IF (.NOT. eigenvectors_discarded) THEN
                     WRITE (UNIT=output_unit, FMT="(T9,A)") &
                        "DIIS | Checking eigenvalues of the DIIS error matrix"
                  END IF
                  WRITE (UNIT=message, FMT="(T9,A,I6,A,ES10.1,A,ES10.1)") &
                     "DIIS | Eigenvalue ", jb, " = ", ev(jb), " is smaller than "// &
                     "threshold ", eigenvalue_threshold
                  CALL compress(message)
                  WRITE (UNIT=output_unit, FMT="(T9,A)") TRIM(message)
                  eigenvectors_discarded = .TRUE.
               END IF
               a(1:nb1, jb) = 0.0_dp
            ELSE
               a(1:nb1, jb) = a(1:nb1, jb)/ev(jb)
            END IF
         END DO

         IF ((output_unit > 0) .AND. eigenvectors_discarded) THEN
            WRITE (UNIT=output_unit, FMT="(T9,A,/)") &
               "DIIS | The corresponding eigenvectors were discarded"
         END IF

         ev(1:nb) = MATMUL(a(1:nb, 1:nb1), b(nb1, 1:nb1))

         ! Update Kohn-Sham matrix

         DO ispin = 1, nspin
            CALL cp_fm_set_all(kc(ispin), 0.0_dp)
            DO jb = 1, nb
               parameters => diis_buffer%param(jb, ispin)
               CALL cp_fm_scale_and_add(1.0_dp, kc(ispin), -ev(jb), parameters)
            END DO
         END DO

         DEALLOCATE (a)
         DEALLOCATE (b)
         DEALLOCATE (ev)

      ELSE

         DO ispin = 1, nspin
            parameters => diis_buffer%param(ib, ispin)
            CALL cp_fm_to_fm(parameters, kc(ispin))
         END DO

      END IF

      CALL cp_print_key_finished_output(output_unit, logger, scf_section, &
                                        "PRINT%DIIS_INFO")

      CALL timestop(handle)

   END SUBROUTINE qs_diis_b_step

! **************************************************************************************************
!> \brief clears the buffer
!> \param diis_buffer the buffer to clear
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   PURE SUBROUTINE qs_diis_b_clear(diis_buffer)

      TYPE(qs_diis_buffer_type), INTENT(INOUT)           :: diis_buffer

      diis_buffer%ncall = 0

   END SUBROUTINE qs_diis_b_clear

! **************************************************************************************************
!> \brief Update the SCF DIIS buffer in linear scaling SCF (LS-SCF),
!>        and if appropriate does a diis step.
!> \param diis_buffer ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param unit_nr ...
!> \param iscf ...
!> \param diis_step ...
!> \param eps_diis ...
!> \param nmixing ...
!> \param s_matrix ...
!> \param threshold ...
!> \par History
!>      - Adapted for LS-SCF (10-11-14) from qs_diis_b_step
!> \author Fredy W. Aquino
! **************************************************************************************************

   SUBROUTINE qs_diis_b_step_4lscf(diis_buffer, qs_env, ls_scf_env, unit_nr, iscf, &
                                   diis_step, eps_diis, nmixing, s_matrix, threshold)
! Note.- Input: ls_scf_env%matrix_p(ispin) , Density Matrix
!               matrix_ks (from qs_env)    , Kohn-Sham Matrix  (IN/OUT)

      TYPE(qs_diis_buffer_type_sparse), POINTER          :: diis_buffer
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type)                              :: ls_scf_env
      INTEGER, INTENT(IN)                                :: unit_nr, iscf
      LOGICAL, INTENT(OUT)                               :: diis_step
      REAL(KIND=dp), INTENT(IN)                          :: eps_diis
      INTEGER, INTENT(IN), OPTIONAL                      :: nmixing
      TYPE(dbcsr_type), OPTIONAL                         :: s_matrix
      REAL(KIND=dp), INTENT(IN)                          :: threshold

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_diis_b_step_4lscf'
      REAL(KIND=dp), PARAMETER :: eigenvalue_threshold = 1.0E-12_dp

      INTEGER                                            :: handle, ib, ispin, jb, my_nmixing, nb, &
                                                            nb1, nspin
      REAL(KIND=dp)                                      :: error_max, tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ev
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: a, b
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dbcsr_type)                                   :: matrix_KSerr_t, matrix_tmp
      TYPE(dbcsr_type), POINTER                          :: new_errors, old_errors, parameters
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      nspin = ls_scf_env%nspins
      diis_step = .FALSE.
      my_nmixing = 2
      IF (PRESENT(nmixing)) my_nmixing = nmixing
      NULLIFY (new_errors, old_errors, parameters, a, b)
      logger => cp_get_default_logger()
      ! Quick return, if no DIIS is requested
      IF (diis_buffer%nbuffer < 1) THEN
         CALL timestop(handle)
         RETURN
      END IF

! Getting current Kohn-Sham matrix from qs_env
      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      matrix_ks=matrix_ks)
      CALL qs_diis_b_check_i_alloc_sparse( &
         diis_buffer, &
         ls_scf_env, &
         nspin)
      error_max = 0.0_dp

      ib = MODULO(diis_buffer%ncall, diis_buffer%nbuffer) + 1
      diis_buffer%ncall = diis_buffer%ncall + 1
      nb = MIN(diis_buffer%ncall, diis_buffer%nbuffer)
! Create scratch arrays
      CALL dbcsr_create(matrix_tmp, &
                        template=ls_scf_env%matrix_ks(1), &
                        matrix_type='N')
      CALL dbcsr_set(matrix_tmp, 0.0_dp) ! reset matrix
      CALL dbcsr_create(matrix_KSerr_t, &
                        template=ls_scf_env%matrix_ks(1), &
                        matrix_type='N')
      CALL dbcsr_set(matrix_KSerr_t, 0.0_dp) ! reset matrix

      DO ispin = 1, nspin ! ------ Loop-ispin----START

         new_errors => diis_buffer%error(ib, ispin)%matrix
         parameters => diis_buffer%param(ib, ispin)%matrix
         ! Copy the Kohn-Sham matrix K to the DIIS buffer
         CALL dbcsr_copy(parameters, & ! out
                         matrix_ks(ispin)%matrix) ! in

         IF (PRESENT(s_matrix)) THEN ! if-s_matrix ---------- START
! Calculate Kohn-Sham error (non-orthogonal)= K*P*S-(K*P*S)^T
! matrix_tmp = P*S
            CALL dbcsr_multiply("N", "N", &
                                1.0_dp, ls_scf_env%matrix_p(ispin), &
                                s_matrix, &
                                0.0_dp, matrix_tmp, &
                                filter_eps=threshold)
! new_errors= K*P*S
            CALL dbcsr_multiply("N", "N", &
                                1.0_dp, matrix_ks(ispin)%matrix, &
                                matrix_tmp, &
                                0.0_dp, new_errors, &
                                filter_eps=threshold)
! matrix_KSerr_t= transpose(K*P*S)
            CALL dbcsr_transposed(matrix_KSerr_t, &
                                  new_errors)
! new_errors=K*P*S-transpose(K*P*S)
            CALL dbcsr_add(new_errors, &
                           matrix_KSerr_t, &
                           1.0_dp, -1.0_dp)
         ELSE ! if-s_matrix ---------- MID
! Calculate Kohn-Sham error (orthogonal)= K*P - P*K
! new_errors=K*P
            CALL dbcsr_multiply("N", "N", &
                                1.0_dp, matrix_ks(ispin)%matrix, &
                                ls_scf_env%matrix_p(ispin), &
                                0.0_dp, new_errors, &
                                filter_eps=threshold)
! matrix_KSerr_t= transpose(K*P)
            CALL dbcsr_transposed(matrix_KSerr_t, &
                                  new_errors)
! new_errors=K*P-transpose(K*P)
            CALL dbcsr_add(new_errors, &
                           matrix_KSerr_t, &
                           1.0_dp, -1.0_dp)
         END IF ! if-s_matrix ---------- END

         tmp = dbcsr_maxabs(new_errors)
         error_max = MAX(error_max, tmp)

      END DO ! ------ Loop-ispin----END

      ! Check, if a DIIS step is appropriate

      diis_step = (diis_buffer%ncall >= my_nmixing)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(A29,I3,A3,4(I3,A1))') &
            "DIIS: (ncall,nbuffer,ib,nb)=(", iscf, ")=(", &
            diis_buffer%ncall, ",", diis_buffer%nbuffer, ",", ib, ",", nb, ")"
         WRITE (unit_nr, '(A57,I3,A3,L1,A1,F10.8,A1,F4.2,A1,L1,A1)') &
            "DIIS: (diis_step,error_max,eps_diis,error_max<eps_diis)=(", &
            iscf, ")=(", diis_step, ",", error_max, ",", eps_diis, ",", &
            (error_max < eps_diis), ")"
         WRITE (unit_nr, '(A75)') &
            "DIIS: diis_step=T : Perform DIIS  error_max<eps_diis=T : Update DIIS buffer"
      END IF

      ! Update the SCF DIIS buffer
      IF (error_max < eps_diis) THEN
         b => diis_buffer%b_matrix
         DO jb = 1, nb
            b(jb, ib) = 0.0_dp
            DO ispin = 1, nspin
               old_errors => diis_buffer%error(jb, ispin)%matrix
               new_errors => diis_buffer%error(ib, ispin)%matrix
               CALL dbcsr_dot(old_errors, &
                              new_errors, &
                              tmp) ! out : < f_i | f_j >
               b(jb, ib) = b(jb, ib) + tmp
            END DO ! end-loop-ispin
            b(ib, jb) = b(jb, ib)
         END DO ! end-loop-jb
      ELSE
         diis_step = .FALSE.
      END IF

      ! Perform DIIS step
      IF (diis_step) THEN
         nb1 = nb + 1
         ALLOCATE (a(nb1, nb1))
         ALLOCATE (b(nb1, nb1))
         ALLOCATE (ev(nb1))
         ! Set up the linear DIIS equation system
         b(1:nb, 1:nb) = diis_buffer%b_matrix(1:nb, 1:nb)
         b(1:nb, nb1) = -1.0_dp
         b(nb1, 1:nb) = -1.0_dp
         b(nb1, nb1) = 0.0_dp
         ! Solve the linear DIIS equation system
         CALL diamat_all(b(1:nb1, 1:nb1), ev(1:nb1))
         a(1:nb1, 1:nb1) = b(1:nb1, 1:nb1)
         DO jb = 1, nb1
            IF (ABS(ev(jb)) < eigenvalue_threshold) THEN
               a(1:nb1, jb) = 0.0_dp
            ELSE
               a(1:nb1, jb) = a(1:nb1, jb)/ev(jb)
            END IF
         END DO ! end-loop-jb

         ev(1:nb) = MATMUL(a(1:nb, 1:nb1), b(nb1, 1:nb1))

         ! Update Kohn-Sham matrix
         IF (iscf .GE. ls_scf_env%iter_ini_diis) THEN ! if-iscf-to-updateKS------ START

            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(A40,I3)') 'DIIS: Updating Kohn-Sham matrix at iscf=', iscf
            END IF

            DO ispin = 1, nspin
               CALL dbcsr_set(matrix_ks(ispin)%matrix, & ! reset matrix
                              0.0_dp)
               DO jb = 1, nb
                  parameters => diis_buffer%param(jb, ispin)%matrix
                  CALL dbcsr_add(matrix_ks(ispin)%matrix, parameters, &
                                 1.0_dp, -ev(jb))
               END DO ! end-loop-jb
            END DO ! end-loop-ispin
         END IF ! if-iscf-to-updateKS------ END

         DEALLOCATE (a)
         DEALLOCATE (b)
         DEALLOCATE (ev)

      ELSE
         DO ispin = 1, nspin
            parameters => diis_buffer%param(ib, ispin)%matrix
            CALL dbcsr_copy(parameters, & ! out
                            matrix_ks(ispin)%matrix) ! in
         END DO ! end-loop-ispin
      END IF
      CALL dbcsr_release(matrix_tmp)
      CALL dbcsr_release(matrix_KSerr_t)
      CALL timestop(handle)

   END SUBROUTINE qs_diis_b_step_4lscf

! **************************************************************************************************
!> \brief Allocate and initialize a DIIS buffer with a buffer size of nbuffer.
!> \param diis_buffer the buffer to initialize
!> \param ls_scf_env ...
!> \param nspin ...
!> \par History
!>      - Adapted from qs_diis_b_check_i_alloc for sparse matrices and
!>        used in LS-SCF module (ls_scf_main) (10-11-14)
!> \author Fredy W. Aquino
!> \note
!>      check to allocate matrices only when needed
! **************************************************************************************************

   SUBROUTINE qs_diis_b_check_i_alloc_sparse(diis_buffer, ls_scf_env, &
                                             nspin)

      TYPE(qs_diis_buffer_type_sparse), INTENT(INOUT)    :: diis_buffer
      TYPE(ls_scf_env_type)                              :: ls_scf_env
      INTEGER, INTENT(IN)                                :: nspin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_diis_b_check_i_alloc_sparse'

      INTEGER                                            :: handle, ibuffer, ispin, nbuffer
      TYPE(cp_logger_type), POINTER                      :: logger

! -------------------------------------------------------------------------

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()

      nbuffer = diis_buffer%nbuffer

      IF (.NOT. ASSOCIATED(diis_buffer%error)) THEN
         ALLOCATE (diis_buffer%error(nbuffer, nspin))

         DO ispin = 1, nspin
            DO ibuffer = 1, nbuffer
               ALLOCATE (diis_buffer%error(ibuffer, ispin)%matrix)

               CALL dbcsr_create(diis_buffer%error(ibuffer, ispin)%matrix, &
                                 template=ls_scf_env%matrix_ks(1), &
                                 matrix_type='N')
            END DO
         END DO
      END IF

      IF (.NOT. ASSOCIATED(diis_buffer%param)) THEN
         ALLOCATE (diis_buffer%param(nbuffer, nspin))

         DO ispin = 1, nspin
            DO ibuffer = 1, nbuffer
               ALLOCATE (diis_buffer%param(ibuffer, ispin)%matrix)
               CALL dbcsr_create(diis_buffer%param(ibuffer, ispin)%matrix, &
                                 template=ls_scf_env%matrix_ks(1), &
                                 matrix_type='N')
            END DO
         END DO
      END IF

      IF (.NOT. ASSOCIATED(diis_buffer%b_matrix)) THEN
         ALLOCATE (diis_buffer%b_matrix(nbuffer + 1, nbuffer + 1))

         diis_buffer%b_matrix = 0.0_dp
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_diis_b_check_i_alloc_sparse

! **************************************************************************************************
!> \brief clears the DIIS buffer in LS-SCF calculation
!> \param diis_buffer the buffer to clear
!> \par History
!>      10-11-14 created [FA] modified from qs_diis_b_clear
!> \author Fredy W. Aquino
! **************************************************************************************************

   PURE SUBROUTINE qs_diis_b_clear_sparse(diis_buffer)

      TYPE(qs_diis_buffer_type_sparse), INTENT(INOUT)    :: diis_buffer

      diis_buffer%ncall = 0

   END SUBROUTINE qs_diis_b_clear_sparse

! **************************************************************************************************
!> \brief Allocates an SCF DIIS buffer for LS-SCF calculation
!> \param diis_buffer the buffer to create
!> \param nbuffer ...
!> \par History
!>      10-11-14 created [FA] modified from qs_diis_b_create
!> \author Fredy W. Aquino
! **************************************************************************************************
   PURE SUBROUTINE qs_diis_b_create_sparse(diis_buffer, nbuffer)

      TYPE(qs_diis_buffer_type_sparse), INTENT(OUT)      :: diis_buffer
      INTEGER, INTENT(in)                                :: nbuffer

      NULLIFY (diis_buffer%b_matrix)
      NULLIFY (diis_buffer%error)
      NULLIFY (diis_buffer%param)
      diis_buffer%nbuffer = nbuffer
      diis_buffer%ncall = 0

   END SUBROUTINE qs_diis_b_create_sparse

! **************************************************************************************************
!> \brief Allocates an SCF DIIS buffer for k-points
!> \param diis_buffer the buffer to create
!> \param nbuffer ...
! **************************************************************************************************
   SUBROUTINE qs_diis_b_create_kp(diis_buffer, nbuffer)

      TYPE(qs_diis_buffer_type_kp), INTENT(OUT)          :: diis_buffer
      INTEGER, INTENT(in)                                :: nbuffer

      NULLIFY (diis_buffer%b_matrix)
      NULLIFY (diis_buffer%error)
      NULLIFY (diis_buffer%param)
      NULLIFY (diis_buffer%smat)
      diis_buffer%nbuffer = nbuffer
      diis_buffer%ncall = 0

   END SUBROUTINE qs_diis_b_create_kp

! **************************************************************************************************
!> \brief Allocate and initialize a DIIS buffer for nao*nao parameter
!>      variables and with a buffer size of nbuffer, in the k-point case
!> \param diis_buffer the buffer to initialize
!> \param matrix_struct the structure for the matrix of the buffer note: this is in the kp subgroup
!> \param nspin ...
!> \param nkp ...
!> \param scf_section ...
! **************************************************************************************************
   SUBROUTINE qs_diis_b_check_i_alloc_kp(diis_buffer, matrix_struct, nspin, nkp, scf_section)

      TYPE(qs_diis_buffer_type_kp), INTENT(INOUT)        :: diis_buffer
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      INTEGER, INTENT(IN)                                :: nspin, nkp
      TYPE(section_vals_type), POINTER                   :: scf_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_diis_b_check_i_alloc_kp'

      INTEGER                                            :: handle, ibuffer, ikp, ispin, nbuffer, &
                                                            output_unit
      TYPE(cp_logger_type), POINTER                      :: logger

! -------------------------------------------------------------------------

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()

      nbuffer = diis_buffer%nbuffer

      IF (.NOT. ASSOCIATED(diis_buffer%error)) THEN
         ALLOCATE (diis_buffer%error(nbuffer, nspin, nkp))

         DO ikp = 1, nkp
            DO ispin = 1, nspin
               DO ibuffer = 1, nbuffer
                  CALL cp_cfm_create(diis_buffer%error(ibuffer, ispin, ikp), &
                                     name="qs_diis_b%error("// &
                                     TRIM(ADJUSTL(cp_to_string(ibuffer)))//","// &
                                     TRIM(ADJUSTL(cp_to_string(ibuffer)))//")", &
                                     matrix_struct=matrix_struct)
               END DO
            END DO
         END DO
      END IF

      IF (.NOT. ASSOCIATED(diis_buffer%param)) THEN
         ALLOCATE (diis_buffer%param(nbuffer, nspin, nkp))

         DO ikp = 1, nkp
            DO ispin = 1, nspin
               DO ibuffer = 1, nbuffer
                  CALL cp_cfm_create(diis_buffer%param(ibuffer, ispin, ikp), &
                                     name="qs_diis_b%param("// &
                                     TRIM(ADJUSTL(cp_to_string(ibuffer)))//","// &
                                     TRIM(ADJUSTL(cp_to_string(ibuffer)))//")", &
                                     matrix_struct=matrix_struct)
               END DO
            END DO
         END DO
      END IF

      IF (.NOT. ASSOCIATED(diis_buffer%smat)) THEN
         ALLOCATE (diis_buffer%smat(nkp))
         DO ikp = 1, nkp
            CALL cp_cfm_create(diis_buffer%smat(ikp), &
                               name="kp_cfm_smat("// &
                               TRIM(ADJUSTL(cp_to_string(ibuffer)))//","// &
                               TRIM(ADJUSTL(cp_to_string(ibuffer)))//")", &
                               matrix_struct=matrix_struct)
         END DO
      END IF

      IF (.NOT. ASSOCIATED(diis_buffer%b_matrix)) THEN
         ALLOCATE (diis_buffer%b_matrix(nbuffer + 1, nbuffer + 1))
         diis_buffer%b_matrix = 0.0_dp
         output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DIIS_INFO", &
                                            extension=".scfLog")
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T9,A)") &
               "DIIS | The SCF DIIS buffer was allocated and initialized"
         END IF
         CALL cp_print_key_finished_output(output_unit, logger, scf_section, &
                                           "PRINT%DIIS_INFO")
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_diis_b_check_i_alloc_kp

! **************************************************************************************************
!> \brief clears the buffer
!> \param diis_buffer the buffer to clear
! **************************************************************************************************
   PURE SUBROUTINE qs_diis_b_clear_kp(diis_buffer)

      TYPE(qs_diis_buffer_type_kp), INTENT(INOUT)        :: diis_buffer

      diis_buffer%ncall = 0

   END SUBROUTINE qs_diis_b_clear_kp

! **************************************************************************************************
!> \brief Update info about the current buffer step ib and the current number of buffers nb
!> \param diis_buffer ...
!> \param ib ...
!> \param nb ...
! **************************************************************************************************
   SUBROUTINE qs_diis_b_info_kp(diis_buffer, ib, nb)
      TYPE(qs_diis_buffer_type_kp), POINTER              :: diis_buffer
      INTEGER, INTENT(OUT)                               :: ib, nb

      ib = MODULO(diis_buffer%ncall, diis_buffer%nbuffer) + 1
      diis_buffer%ncall = diis_buffer%ncall + 1
      nb = MIN(diis_buffer%ncall, diis_buffer%nbuffer)

   END SUBROUTINE qs_diis_b_info_kp

! **************************************************************************************************
!> \brief Calculate and store the error for a given k-point
!> \param diis_buffer ...
!> \param ib ...
!> \param mos ...
!> \param kc ...
!> \param sc ...
!> \param ispin ...
!> \param ikp ...
!> \param nkp_local ...
!> \param scf_section ...
!> \note We assume that we always have an overlap matrix and complex matrices
!> TODO: do we need to pass the kp weight for the back Fourier transform?
! **************************************************************************************************
   SUBROUTINE qs_diis_b_calc_err_kp(diis_buffer, ib, mos, kc, sc, ispin, ikp, nkp_local, scf_section)
      TYPE(qs_diis_buffer_type_kp), POINTER              :: diis_buffer
      INTEGER, INTENT(IN)                                :: ib
      TYPE(mo_set_type), DIMENSION(:, :), POINTER        :: mos
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: kc, sc
      INTEGER, INTENT(IN)                                :: ispin, ikp, nkp_local
      TYPE(section_vals_type), POINTER                   :: scf_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_diis_b_calc_err_kp'

      INTEGER                                            :: handle, homo, nao, nmo, nspin
      REAL(dp)                                           :: maxocc
      TYPE(cp_cfm_type)                                  :: cmos
      TYPE(cp_cfm_type), POINTER                         :: new_errors, parameters, smat
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), POINTER                          :: imos, rmos

      NULLIFY (matrix_struct, imos, rmos, parameters, new_errors, smat)

      CALL timeset(routineN, handle)

      !Calculate the error for this given k-point, store the KS matrix as well as the ovlp matrix
      !All of this happens within the kp subgroups

      ! Quick return, if no DIIS is requested
      IF (diis_buffer%nbuffer < 1) THEN
         CALL timestop(handle)
         RETURN
      END IF
      nspin = SIZE(mos, 2)

      CALL cp_cfm_get_info(kc, matrix_struct=matrix_struct)
      CALL qs_diis_b_check_i_alloc_kp(diis_buffer, &
                                      matrix_struct=matrix_struct, &
                                      nspin=nspin, nkp=nkp_local, &
                                      scf_section=scf_section)

      !We calculate: e(ikp) = F(ikp)*P(ikp)*S(ikp) - S(ikp)*P(ikp)*F(ikp)
      CALL get_mo_set(mos(1, ispin), nao=nao, nmo=nmo, homo=homo, mo_coeff=rmos, maxocc=maxocc)
      CALL get_mo_set(mos(2, ispin), mo_coeff=imos)
      NULLIFY (matrix_struct)
      CALL cp_fm_get_info(rmos, matrix_struct=matrix_struct)
      CALL cp_cfm_create(cmos, matrix_struct)
      CALL cp_fm_to_cfm(rmos, imos, cmos)

      new_errors => diis_buffer%error(ib, ispin, ikp)
      parameters => diis_buffer%param(ib, ispin, ikp)
      smat => diis_buffer%smat(ikp)

      !copy the KS and overlap matrices to the DIIS buffer
      CALL cp_cfm_to_cfm(kc, parameters)
      CALL cp_cfm_to_cfm(sc, smat)

      ! KC <- K*C
      CALL parallel_gemm("N", "N", nao, homo, nao, CMPLX(maxocc, KIND=dp), parameters, cmos, (0.0_dp, 0.0_dp), kc)
      ! SC <- S*C
      CALL parallel_gemm("N", "N", nao, homo, nao, (2.0_dp, 0.0_dp), smat, cmos, (0.0_dp, 0.0_dp), sc)

      ! new_errors <- KC*(SC)^T - (SC)*(KC)^T = K*P*S - S*P*K
      CALL parallel_gemm("N", "T", nao, nao, homo, (1.0_dp, 0.0_dp), sc, kc, (0.0_dp, 0.0_dp), new_errors)
      CALL parallel_gemm("N", "T", nao, nao, homo, (1.0_dp, 0.0_dp), kc, sc, (-1.0_dp, 0.0_dp), new_errors)

      !clean-up
      CALL cp_cfm_release(cmos)

      CALL timestop(handle)

   END SUBROUTINE qs_diis_b_calc_err_kp

! **************************************************************************************************
!> \brief Update the SCF DIIS buffer, and if appropriate does a diis step, for k-points
!> \param diis_buffer ...
!> \param coeffs ...
!> \param ib ...
!> \param nb ...
!> \param delta ...
!> \param error_max ...
!> \param diis_step ...
!> \param eps_diis ...
!> \param nspin ...
!> \param nkp ...
!> \param nkp_local ...
!> \param nmixing ...
!> \param scf_section ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE qs_diis_b_step_kp(diis_buffer, coeffs, ib, nb, delta, error_max, diis_step, eps_diis, &
                                nspin, nkp, nkp_local, nmixing, scf_section, para_env)

      TYPE(qs_diis_buffer_type_kp), POINTER              :: diis_buffer
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(INOUT)      :: coeffs
      INTEGER, INTENT(IN)                                :: ib, nb
      REAL(KIND=dp), INTENT(IN)                          :: delta
      REAL(KIND=dp), INTENT(OUT)                         :: error_max
      LOGICAL, INTENT(OUT)                               :: diis_step
      REAL(KIND=dp), INTENT(IN)                          :: eps_diis
      INTEGER, INTENT(IN)                                :: nspin, nkp, nkp_local
      INTEGER, INTENT(IN), OPTIONAL                      :: nmixing
      TYPE(section_vals_type), POINTER                   :: scf_section
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_diis_b_step_kp'
      REAL(KIND=dp), PARAMETER :: eigenvalue_threshold = 1.0E-12_dp

      CHARACTER(LEN=2*default_string_length)             :: message
      COMPLEX(KIND=dp)                                   :: tmp
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: a, b
      INTEGER                                            :: handle, ikp, ispin, jb, my_nmixing, nb1, &
                                                            output_unit
      LOGICAL                                            :: eigenvectors_discarded
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ev
      TYPE(cp_cfm_type)                                  :: old_errors
      TYPE(cp_cfm_type), POINTER                         :: new_errors
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type)                                   :: ierr, rerr
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (matrix_struct, new_errors, logger)

      CALL timeset(routineN, handle)

      diis_step = .FALSE.

      my_nmixing = 2
      IF (PRESENT(nmixing)) my_nmixing = nmixing

      logger => cp_get_default_logger()

      ! Quick return, if no DIIS is requested
      IF (diis_buffer%nbuffer < 1) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! Check, if a DIIS step is appropriate
      diis_step = ((diis_buffer%ncall >= my_nmixing) .AND. (delta < eps_diis))

      ! Calculate the DIIS buffer, and update it if max_error < eps_diis
      CALL cp_cfm_get_info(diis_buffer%error(ib, 1, 1), matrix_struct=matrix_struct)
      CALL cp_fm_create(ierr, matrix_struct)
      CALL cp_fm_create(rerr, matrix_struct)
      CALL cp_cfm_create(old_errors, matrix_struct)
      ALLOCATE (b(nb, nb))
      b = 0.0_dp
      DO jb = 1, nb
         DO ikp = 1, nkp_local
            DO ispin = 1, nspin
               new_errors => diis_buffer%error(ib, ispin, ikp)
               CALL cp_cfm_to_fm(diis_buffer%error(jb, ispin, ikp), rerr, ierr)
               CALL cp_fm_scale(-1.0_dp, ierr)
               CALL cp_fm_to_cfm(rerr, ierr, old_errors)
               CALL cp_cfm_trace(old_errors, new_errors, tmp)
               b(jb, ib) = b(jb, ib) + 1.0_dp/REAL(nkp, dp)*tmp
            END DO
         END DO
         b(ib, jb) = CONJG(b(jb, ib))
      END DO
      CALL cp_fm_release(ierr)
      CALL cp_fm_release(rerr)
      CALL cp_cfm_release(old_errors)
      CALL para_env%sum(b)

      error_max = SQRT(REAL(b(ib, ib))**2 + AIMAG(b(ib, ib))**2)

      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DIIS_INFO", &
                                         extension=".scfLog")
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T9,A,I4,/,(T9,A,ES12.3))") &
            "DIIS | Current SCF DIIS buffer size:         ", nb, &
            "DIIS | Maximum SCF DIIS error at last step:  ", error_max, &
            "DIIS | Current SCF convergence:              ", delta, &
            "DIIS | Threshold value for a DIIS step:      ", eps_diis
         IF (error_max < eps_diis) THEN
            WRITE (UNIT=output_unit, FMT="(T9,A)") &
               "DIIS | => The SCF DIIS buffer will be updated"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T9,A)") &
               "DIIS | => No update of the SCF DIIS buffer"
         END IF
         IF (diis_step .AND. (error_max < eps_diis)) THEN
            WRITE (UNIT=output_unit, FMT="(T9,A,/)") &
               "DIIS | => A SCF DIIS step will be performed"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T9,A,/)") &
               "DIIS | => No SCF DIIS step will be performed"
         END IF
      END IF

      ! Update the SCF DIIS buffer
      IF (error_max < eps_diis) THEN
         DO jb = 1, nb
            diis_buffer%b_matrix(ib, jb) = b(ib, jb)
            diis_buffer%b_matrix(jb, ib) = b(jb, ib)
         END DO
      ELSE

         diis_step = .FALSE.
      END IF
      DEALLOCATE (b)

      ! Perform DIIS step
      IF (diis_step) THEN

         nb1 = nb + 1

         ALLOCATE (a(nb1, nb1))
         ALLOCATE (b(nb1, nb1))
         ALLOCATE (ev(nb1))

         ! Set up the linear DIIS equation system
         b(1:nb, 1:nb) = diis_buffer%b_matrix(1:nb, 1:nb)

         b(1:nb, nb1) = -1.0_dp
         b(nb1, 1:nb) = -1.0_dp
         b(nb1, nb1) = 0.0_dp

         ! Solve the linear DIIS equation system
         ev(1:nb1) = 0.0_dp !eigenvalues
         a(1:nb1, 1:nb1) = 0.0_dp !eigenvectors
         CALL diag_complex(b(1:nb1, 1:nb1), a(1:nb1, 1:nb1), ev(1:nb1))
         b(1:nb1, 1:nb1) = a(1:nb1, 1:nb1)

         eigenvectors_discarded = .FALSE.

         DO jb = 1, nb1
            IF (ABS(ev(jb)) < eigenvalue_threshold) THEN
               IF (output_unit > 0) THEN
                  IF (.NOT. eigenvectors_discarded) THEN
                     WRITE (UNIT=output_unit, FMT="(T9,A)") &
                        "DIIS | Checking eigenvalues of the DIIS error matrix"
                  END IF
                  WRITE (UNIT=message, FMT="(T9,A,I6,A,ES10.1,A,ES10.1)") &
                     "DIIS | Eigenvalue ", jb, " = ", ev(jb), " is smaller than "// &
                     "threshold ", eigenvalue_threshold
                  CALL compress(message)
                  WRITE (UNIT=output_unit, FMT="(T9,A)") TRIM(message)
                  eigenvectors_discarded = .TRUE.
               END IF
               a(1:nb1, jb) = 0.0_dp
            ELSE
               a(1:nb1, jb) = a(1:nb1, jb)/ev(jb)
            END IF
         END DO

         IF ((output_unit > 0) .AND. eigenvectors_discarded) THEN
            WRITE (UNIT=output_unit, FMT="(T9,A,/)") &
               "DIIS | The corresponding eigenvectors were discarded"
         END IF

         coeffs(1:nb) = -MATMUL(a(1:nb, 1:nb1), CONJG(b(nb1, 1:nb1)))
      ELSE

         coeffs(:) = 0.0_dp
         coeffs(ib) = 1.0_dp
      END IF

      CALL cp_print_key_finished_output(output_unit, logger, scf_section, &
                                        "PRINT%DIIS_INFO")

      CALL timestop(handle)

   END SUBROUTINE qs_diis_b_step_kp
END MODULE qs_diis
