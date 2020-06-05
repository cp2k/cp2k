!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines saving to files various data related to POD method
!> \author Z. Futera (03.2019)
! **************************************************************************************************
MODULE et_coupling_pod_save

   USE cp_blacs_calls,                  ONLY: cp_blacs_gridinfo,&
                                              cp_blacs_gridinit,&
                                              cp_blacs_gridexit
   USE cp_files,                        ONLY: open_file
   USE cp_fm_types,                     ONLY: cp_fm_type,&
                                              cp_fm_write_unformatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_should_output,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE et_coupling_pod_mo,              ONLY: pod_mo_save
   USE et_coupling_pod_types,           ONLY: et_cpl
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE message_passing,                 ONLY: mp_recv,&
                                              mp_send
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_mo_types,                     ONLY: mo_set_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_pod_save'

   PUBLIC :: pod_save_coupl_close,&
             pod_save_coupl_elements,&
             pod_save_coupl_header,&
             pod_save_coupl_states,&
             pod_save_state_close,&
             pod_save_state_coeffs,&
             pod_save_state_cube,&
             pod_save_state_header,&
             pod_save_transf_mat

CONTAINS

! **************************************************************************************************
!> \brief save header of binary POD coupling file
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param n_spins number of spin components
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_coupl_header(qs_env, ec, n_spins, unit_nr)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: n_spins
      INTEGER, INTENT(OUT)                               :: unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_coupl_header',&
                                     routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: i
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      ! Initialization
      NULLIFY(logger)
      NULLIFY(print_key)

      ! Output settings
      logger => cp_get_default_logger()
      print_key => section_vals_get_subs_vals(qs_env%input, &
                     'PROPERTIES%ET_COUPLING%PROJECTION%PRINT')

      ! Open the file
      unit_nr = cp_print_key_unit_nr(logger,print_key,'HAMILTONIAN', &
                   middle_name='RESTART',extension='.hab',file_form='UNFORMATTED', &
                   file_position='REWIND',file_action='WRITE',file_status='REPLACE')

      IF (unit_nr > 0) THEN

         ! Number of spin components
         WRITE(UNIT=unit_nr) n_spins
         ! Number of blocks
         WRITE(UNIT=unit_nr) ec%n_blocks
         ! Block dimensions
         WRITE(UNIT=unit_nr) (ec%block(i)%n_ao,i=1,ec%n_blocks)

         FLUSH(unit_nr)

      END IF

   END SUBROUTINE pod_save_coupl_header

! **************************************************************************************************
!> \brief save POD block energies to binary data file
!> \param ec electronic coupling data structure
!> \param ib POD block ID
!> \param is spin component
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_coupl_states(ec,ib,is,unit_nr)

      ! Routine arguments
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: ib
      INTEGER                                            :: is
      INTEGER                                            :: unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_coupl_states',&
                                     routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: i

      IF (unit_nr > 0) THEN

         ! Save block energies
         WRITE(UNIT=unit_nr) (ec%block(ib)%mo(is)%mo_set%eigenvalues(i), &
                              i=1,ec%block(ib)%mo(is)%mo_set%nmo)

         FLUSH(unit_nr)

      END IF


   END SUBROUTINE pod_save_coupl_states

! **************************************************************************************************
!> \brief save POD coupling elements from one off-diagonal block to binary data file
!> \param ec electronic coupling data structure
!> \param i1 ID of the first POD block
!> \param i2 ID of the second POD block
!> \param is spin component
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_coupl_elements(ec,i1,i2,is,unit_nr)

      ! Routine arguments
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: i1
      INTEGER                                            :: i2
      INTEGER                                            :: is
      INTEGER                                            :: unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_coupl_elements',&
                                     routineP = moduleN//':'//routineN
 
      ! Save data to file
      CALL cp_fm_write_unformatted(ec%block(i1)%hab(is,i2)%matrix,unit_nr)

      IF (unit_nr > 0) &
         FLUSH(unit_nr)

   END SUBROUTINE pod_save_coupl_elements

! **************************************************************************************************
!> \brief close output stream for saving the couplings
!> \param qs_env QuickStep environment containing all system data
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_coupl_close(qs_env, unit_nr)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(OUT)                               :: unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_coupl_close',&
                                     routineP = moduleN//':'//routineN

      ! Local variables
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      ! Initialization
      NULLIFY(logger)
      NULLIFY(print_key)

      ! Output settings
      logger => cp_get_default_logger()
      print_key => section_vals_get_subs_vals(qs_env%input, &
                     'PROPERTIES%ET_COUPLING%PROJECTION%PRINT')

      ! Close the file
      CALL cp_print_key_finished_output(unit_nr,logger,print_key,'HAMILTONIAN')

   END SUBROUTINE pod_save_coupl_close

! **************************************************************************************************
!> \brief save specified electronic states to cube files
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param n_spins number of spin states
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_save_state_cube(qs_env, ec, n_spins)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: n_spins

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_state_cube',&
                                     routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: i, j, k, l, n
      INTEGER, DIMENSION(:), POINTER                     :: list
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_set_type), POINTER                         :: mo
      TYPE(section_vals_type), POINTER                   :: block_sec, mo_sec, print_sec

      NULLIFY (logger)
      NULLIFY (block_sec)
      NULLIFY (print_sec)
      NULLIFY (mo_sec)

      ! Output logger
      logger => cp_get_default_logger()
      block_sec => section_vals_get_subs_vals(qs_env%input, &
                                              'PROPERTIES%ET_COUPLING%PROJECTION%BLOCK')

      ! Print states of all blocks
      DO i = 1, ec%n_blocks

         print_sec => section_vals_get_subs_vals(block_sec, 'PRINT', i_rep_section=i)

         ! Check if the print input section is active
         IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                              print_sec, 'MO_CUBES'), cp_p_file)) THEN

            mo_sec => section_vals_get_subs_vals(print_sec, 'MO_CUBES')

            ! Spin states
            DO j = 1, n_spins

               mo => ec%block(i)%mo(j)%mo_set

               CALL section_vals_val_get(mo_sec, keyword_name='MO_LIST', n_rep_val=n)

               ! List of specific MOs
               IF (n > 0) THEN

                  DO k = 1, n
                     NULLIFY (list)
                     CALL section_vals_val_get(mo_sec, keyword_name='MO_LIST', &
                                               i_rep_val=k, i_vals=list)
                     IF (ASSOCIATED(list)) THEN
                        DO l = 1, SIZE(list)
                           CALL pod_mo_save(qs_env, logger, print_sec, mo, i, list(l), j)
                        END DO
                     END IF
                  END DO

                  ! Frontier MOs
               ELSE

                  ! Occupied states
                  CALL section_vals_val_get(mo_sec, keyword_name='NHOMO', i_val=n)

                  IF (n > 0) THEN
                     DO k = MAX(1, mo%homo-n+1), mo%homo
                        CALL pod_mo_save(qs_env, logger, print_sec, mo, i, k, j)
                     END DO
                  END IF

                  ! Unoccupied states
                  CALL section_vals_val_get(mo_sec, keyword_name='NLUMO', i_val=n)

                  IF (n > 0) THEN
                     DO k = mo%lfomo, MIN(mo%lfomo+n-1, mo%nmo)
                        CALL pod_mo_save(qs_env, logger, print_sec, mo, i, k, j)
                     END DO
                  END IF

               END IF

            END DO

         END IF

      END DO

   END SUBROUTINE pod_save_state_cube

! **************************************************************************************************
!> \brief save header of binary POD wavefunction file
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param n_spins number of spin components
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_state_header(qs_env, ec, n_spins, unit_nr)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: n_spins
      INTEGER, INTENT(OUT)                               :: unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_state_header',&
                                     routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: i
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      ! Initialization
      NULLIFY(logger)
      NULLIFY(print_key)

      ! Output settings
      logger => cp_get_default_logger()
      print_key => section_vals_get_subs_vals(qs_env%input, &
                     'PROPERTIES%ET_COUPLING%PROJECTION%PRINT')

      ! Open the file
      unit_nr = cp_print_key_unit_nr(logger,print_key,'WAVEFUNCTION', &
                   middle_name='RESTART',extension='.wfn',file_form='UNFORMATTED', &
                   file_position='REWIND',file_action='WRITE',file_status='REPLACE')

      IF (unit_nr > 0) THEN

         ! Number of spin components
         WRITE(UNIT=unit_nr) n_spins
         ! Number of blocks
         WRITE(UNIT=unit_nr) ec%n_blocks
         ! Block dimensions
         WRITE(UNIT=unit_nr) (ec%block(i)%n_ao,i=1,ec%n_blocks)

         FLUSH(unit_nr)

      END IF

   END SUBROUTINE pod_save_state_header

! **************************************************************************************************
!> \brief save POD orbital coefficients to binary data file
!> \param ec electronic coupling data structure
!> \param ib POD block ID
!> \param is spin component
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_state_coeffs(ec,ib,is,unit_nr)

      ! Routine arguments
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: ib
      INTEGER                                            :: is
      INTEGER                                            :: unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_state_coeffs',&
                                     routineP = moduleN//':'//routineN

      ! Save data to file
      CALL cp_fm_write_unformatted(ec%block(ib)%mo(is)%mo_set%mo_coeff,unit_nr)

      IF (unit_nr > 0) &
         FLUSH(unit_nr)
      
   END SUBROUTINE pod_save_state_coeffs

! **************************************************************************************************
!> \brief close output stream for saving the localized wavefunctions
!> \param qs_env QuickStep environment containing all system data
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_state_close(qs_env, unit_nr)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(OUT)                               :: unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_state_close',&
                                     routineP = moduleN//':'//routineN

      ! Local variables
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      ! Initialization
      NULLIFY(logger)
      NULLIFY(print_key)

      ! Output settings
      logger => cp_get_default_logger()
      print_key => section_vals_get_subs_vals(qs_env%input, &
                     'PROPERTIES%ET_COUPLING%PROJECTION%PRINT')

      ! Close the file
      CALL cp_print_key_finished_output(unit_nr,logger,print_key,'WAVEFUNCTION')

   END SUBROUTINE pod_save_state_close

! **************************************************************************************************
!> \brief save transformation matrix S^{+/-1/2} to binary data file
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param key 'F' for forward, 'B' for backward
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_transf_mat(qs_env,ec,key)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      CHARACTER                                          :: key

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_transf_mat',&
                                     routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      ! Initialization
      NULLIFY(logger)
      NULLIFY(print_key)

      ! Output settings
      logger => cp_get_default_logger()
      print_key => section_vals_get_subs_vals(qs_env%input, &
                     'PROPERTIES%ET_COUPLING%PROJECTION%PRINT')

      ! Forward transformation matrix
      IF (key == 'F') THEN

         unit_nr = cp_print_key_unit_nr(logger,print_key,'TRANSF_MATRIX_F', &
                      middle_name='RESTART',extension='.tmf',file_form='UNFORMATTED', &
                      file_position='REWIND',file_action='WRITE',file_status='REPLACE')

         CALL cp_fm_write_unformatted(ec%tmatf,unit_nr)

         CALL cp_print_key_finished_output(unit_nr,logger,print_key,'TRANSF_MATRIX_F')

      ! Backward transformation matrix
      ELSE IF (key == 'B') THEN

         unit_nr = cp_print_key_unit_nr(logger,print_key,'TRANSF_MATRIX_R', &
                      middle_name='RESTART',extension='.tmr',file_form='UNFORMATTED', &
                      file_position='REWIND',file_action='WRITE',file_status='REPLACE')

         CALL cp_fm_write_unformatted(ec%tmatr,unit_nr)

         CALL cp_print_key_finished_output(unit_nr,logger,print_key,'TRANSF_MATRIX_R')

      END IF

   END SUBROUTINE pod_save_transf_mat

END MODULE et_coupling_pod_save
