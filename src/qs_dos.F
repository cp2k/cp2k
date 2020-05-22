!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation and writing of density of  states
!> \par History
!>      -
!> \author JGH
! **************************************************************************************************
MODULE qs_dos
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE message_passing,                 ONLY: mp_max,&
                                              mp_min
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dos'

   PUBLIC :: calculate_dos, calculate_dos_kp

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Compute and write density of states
!> \param mos ...
!> \param qs_env ...
!> \param dft_section ...
!> \date    26.02.2008
!> \par History:
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE calculate_dos(mos, qs_env, dft_section)

      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_dos', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: my_act, my_pos
      INTEGER                                            :: handle, i, iounit, ispin, iterstep, iv, &
                                                            iw, nhist, nmo(2), nspins
      LOGICAL                                            :: append, ionode, should_output
      REAL(KIND=dp)                                      :: de, e1, e2, e_fermi(2), emax, emin, eval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ehist, hist, occval
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_set_type), POINTER                         :: mo_set

      CPASSERT(ASSOCIATED(qs_env))
      MARK_USED(qs_env)

      NULLIFY (logger)
      logger => cp_get_default_logger()
      ionode = logger%para_env%ionode
      should_output = BTEST(cp_print_key_should_output(logger%iter_info, dft_section, &
                                                       "PRINT%DOS"), cp_p_file)
      iounit = cp_logger_get_default_io_unit(logger)
      IF ((.NOT. should_output)) RETURN

      CALL timeset(routineN, handle)
      iterstep = logger%iter_info%iteration(logger%iter_info%n_rlevel)

      IF (iounit > 0) WRITE (UNIT=iounit, FMT='(/,(T3,A,T61,I10))') &
         " Calculate DOS at iteration step ", iterstep

      CALL section_vals_val_get(dft_section, "PRINT%DOS%DELTA_E", r_val=de)
      CALL section_vals_val_get(dft_section, "PRINT%PDOS%APPEND", l_val=append)
      IF (append .AND. iterstep > 1) THEN
         my_pos = "APPEND"
      ELSE
         my_pos = "REWIND"
      END IF

      emin = 1.e10_dp
      emax = -1.e10_dp
      nspins = SIZE(mos)
      DO ispin = 1, nspins
         mo_set => mos(ispin)%mo_set
         CALL get_mo_set(mo_set=mo_set, nmo=nmo(ispin), mu=e_fermi(ispin))
         eigenvalues => mo_set%eigenvalues
         e1 = MINVAL(eigenvalues(1:nmo(ispin)))
         e2 = MAXVAL(eigenvalues(1:nmo(ispin)))
         emin = MIN(emin, e1)
         emax = MAX(emax, e2)
      END DO

      IF (de > 0.0_dp) THEN
         nhist = NINT((emax - emin)/de) + 1
         ALLOCATE (hist(nhist, nspins), occval(nhist, nspins), ehist(nhist, nspins))
         hist = 0.0_dp
         occval = 0.0_dp
         ehist = 0.0_dp
         DO ispin = 1, nspins
            mo_set => mos(ispin)%mo_set
            occupation_numbers => mo_set%occupation_numbers
            eigenvalues => mo_set%eigenvalues
            DO i = 1, nmo(ispin)
               eval = eigenvalues(i) - emin
               iv = NINT(eval/de) + 1
               CPASSERT((iv > 0) .AND. (iv <= nhist))
               hist(iv, ispin) = hist(iv, ispin) + 1.0_dp
               occval(iv, ispin) = occval(iv, ispin) + occupation_numbers(i)
            END DO
            hist(:, ispin) = hist(:, ispin)/REAL(nmo(ispin), KIND=dp)
         END DO
         DO i = 1, nhist
            ehist(i, 1:nspins) = emin + (i - 1)*de
         END DO
      ELSE
         nhist = MAXVAL(nmo)
         ALLOCATE (hist(nhist, nspins), occval(nhist, nspins), ehist(nhist, nspins))
         hist = 0.0_dp
         occval = 0.0_dp
         ehist = 0.0_dp
         DO ispin = 1, nspins
            mo_set => mos(ispin)%mo_set
            occupation_numbers => mo_set%occupation_numbers
            eigenvalues => mo_set%eigenvalues
            DO i = 1, nmo(ispin)
               ehist(i, ispin) = eigenvalues(i)
               hist(i, ispin) = 1.0_dp
               occval(i, ispin) = occupation_numbers(i)
            END DO
            hist(:, ispin) = hist(:, ispin)/REAL(nmo(ispin), KIND=dp)
         END DO
      END IF

      my_act = "WRITE"
      iw = cp_print_key_unit_nr(logger, dft_section, "PRINT%DOS", &
                                extension=".dos", file_position=my_pos, file_action=my_act, &
                                file_form="FORMATTED")
      IF (iw > 0) THEN
         IF (nspins == 2) THEN
            WRITE (UNIT=iw, FMT="(T2,A,I0,A,2F12.6)") &
               "# DOS at iteration step i = ", iterstep, ", E_Fermi[a.u.] = ", e_fermi(1:2)
            WRITE (UNIT=iw, FMT="(T2,A, A)") "#  Energy[a.u.]  Alpha_Density     Occupation", &
               "   Energy[a.u.]  Beta_Density      Occupation"
         ELSE
            WRITE (UNIT=iw, FMT="(T2,A,I0,A,F12.6)") &
               "# DOS at iteration step i = ", iterstep, ", E_Fermi[a.u.] = ", e_fermi(1)
            WRITE (UNIT=iw, FMT="(T2,A)") "#  Energy[a.u.]       Density     Occupation"
         END IF
         DO i = 1, nhist
            IF (nspins == 2) THEN
               e1 = ehist(i, 1)
               e2 = ehist(i, 2)
               WRITE (UNIT=iw, FMT="(2(F15.8,2F15.4))") e1, hist(i, 1), occval(i, 1), &
                  e2, hist(i, 2), occval(i, 2)
            ELSE
               eval = ehist(i, 1)
               WRITE (UNIT=iw, FMT="(F15.8,2F15.4)") eval, hist(i, 1), occval(i, 1)
            END IF
         END DO
      END IF
      CALL cp_print_key_finished_output(iw, logger, dft_section, "PRINT%DOS")
      DEALLOCATE (hist, occval, ehist)

      CALL timestop(handle)

   END SUBROUTINE calculate_dos

! **************************************************************************************************
!> \brief Compute and write density of states (kpoints)
!> \param kpoints ...
!> \param qs_env ...
!> \param dft_section ...
!> \date    26.02.2008
!> \par History:
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE calculate_dos_kp(kpoints, qs_env, dft_section)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_dos_kp', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: my_act, my_pos
      INTEGER                                            :: handle, i, ik, iounit, ispin, iterstep, &
                                                            iv, iw, nhist, nmo(2), nspins
      LOGICAL                                            :: append, ionode, should_output
      REAL(KIND=dp)                                      :: de, e1, e2, emax, emin, eval, wkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ehist, hist, occval
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:, :), POINTER      :: mos
      TYPE(mo_set_type), POINTER                         :: mo_set

      NULLIFY (logger)
      logger => cp_get_default_logger()
      ionode = logger%para_env%ionode
      should_output = BTEST(cp_print_key_should_output(logger%iter_info, dft_section, &
                                                       "PRINT%DOS"), cp_p_file)
      iounit = cp_logger_get_default_io_unit(logger)
      IF ((.NOT. should_output)) RETURN

      CALL timeset(routineN, handle)
      iterstep = logger%iter_info%iteration(logger%iter_info%n_rlevel)

      IF (iounit > 0) WRITE (UNIT=iounit, FMT='(/,(T3,A,T61,I10))') &
         " Calculate DOS at iteration step ", iterstep

      CALL section_vals_val_get(dft_section, "PRINT%DOS%DELTA_E", r_val=de)
      CALL section_vals_val_get(dft_section, "PRINT%DOS%APPEND", l_val=append)
      ! ensure a lower value for the histogram width
      de = MAX(de, 0.00001_dp)
      IF (append .AND. iterstep > 1) THEN
         my_pos = "APPEND"
      ELSE
         my_pos = "REWIND"
      END IF

      CALL get_qs_env(qs_env, dft_control=dft_control, para_env=para_env)
      nspins = dft_control%nspins

      emin = 1.e10_dp
      emax = -1.e10_dp
      IF (kpoints%nkp /= 0) THEN
         DO ik = 1, SIZE(kpoints%kp_env)
            mos => kpoints%kp_env(ik)%kpoint_env%mos
            CPASSERT(ASSOCIATED(mos))
            DO ispin = 1, nspins
               mo_set => mos(1, ispin)%mo_set
               CALL get_mo_set(mo_set=mo_set, nmo=nmo(ispin))
               eigenvalues => mo_set%eigenvalues
               e1 = MINVAL(eigenvalues(1:nmo(ispin)))
               e2 = MAXVAL(eigenvalues(1:nmo(ispin)))
               emin = MIN(emin, e1)
               emax = MAX(emax, e2)
            END DO
         END DO
      END IF
      CALL mp_min(emin, para_env%group)
      CALL mp_max(emax, para_env%group)
      CALL mp_max(nmo, para_env%group)

      nhist = NINT((emax - emin)/de) + 1
      ALLOCATE (hist(nhist, nspins), occval(nhist, nspins), ehist(nhist, nspins))
      hist = 0.0_dp
      occval = 0.0_dp
      ehist = 0.0_dp

      IF (kpoints%nkp /= 0) THEN
         DO ik = 1, SIZE(kpoints%kp_env)
            mos => kpoints%kp_env(ik)%kpoint_env%mos
            wkp = kpoints%kp_env(ik)%kpoint_env%wkp
            DO ispin = 1, nspins
               mo_set => mos(1, ispin)%mo_set
               occupation_numbers => mo_set%occupation_numbers
               eigenvalues => mo_set%eigenvalues
               DO i = 1, nmo(ispin)
                  eval = eigenvalues(i) - emin
                  iv = NINT(eval/de) + 1
                  CPASSERT((iv > 0) .AND. (iv <= nhist))
                  hist(iv, ispin) = hist(iv, ispin) + wkp
                  occval(iv, ispin) = occval(iv, ispin) + wkp*occupation_numbers(i)
               END DO
            END DO
         END DO
      END IF
      DO ispin = 1, nspins
         hist(:, ispin) = hist(:, ispin)/REAL(nmo(ispin), KIND=dp)
      END DO
      DO i = 1, nhist
         ehist(i, 1:nspins) = emin + (i - 1)*de
      END DO

      my_act = "WRITE"
      iw = cp_print_key_unit_nr(logger, dft_section, "PRINT%DOS", &
                                extension=".dos", file_position=my_pos, file_action=my_act, &
                                file_form="FORMATTED")
      IF (iw > 0) THEN
         IF (nspins == 2) THEN
            WRITE (UNIT=iw, FMT="(T2,A,I0)") "# DOS at iteration step i = ", iterstep
            WRITE (UNIT=iw, FMT="(T2,A, A)") "   Energy[a.u.]  Alpha_Density     Occupation", &
               "   Beta_Density      Occupation"
         ELSE
            WRITE (UNIT=iw, FMT="(T2,A,I0)") "# DOS at iteration step i = ", iterstep
            WRITE (UNIT=iw, FMT="(T2,A)") "   Energy[a.u.]       Density     Occupation"
         END IF
         DO i = 1, nhist
            eval = emin + (i - 1)*de
            IF (nspins == 2) THEN
               WRITE (UNIT=iw, FMT="(F15.8,4F15.4)") eval, hist(i, 1), occval(i, 1), &
                  hist(i, 2), occval(i, 2)
            ELSE
               WRITE (UNIT=iw, FMT="(F15.8,2F15.4)") eval, hist(i, 1), occval(i, 1)
            END IF
         END DO
      END IF
      CALL cp_print_key_finished_output(iw, logger, dft_section, "PRINT%DOS")
      DEALLOCATE (hist, occval)

      CALL timestop(handle)

   END SUBROUTINE calculate_dos_kp

END MODULE qs_dos

