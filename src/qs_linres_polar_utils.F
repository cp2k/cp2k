!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Polarizability calculation by dfpt
!>      Initialization of the polar_env,
!>      Perturbation Hamiltonian by application of the Berry phase operator to psi0
!>      Write output
!>      Deallocate everything
!> periodic Raman SL February 2013
!> \note
! **************************************************************************************************
MODULE qs_linres_polar_utils
   USE bibliography,                    ONLY: Luber2014,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_trace
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_result_methods,               ONLY: cp_results_erase,&
                                              put_results
   USE cp_result_types,                 ONLY: cp_result_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE machine,                         ONLY: m_flush
   USE mathconstants,                   ONLY: twopi
   USE physcon,                         ONLY: angstrom
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_linres_methods,               ONLY: linres_read_restart,&
                                              linres_solver,&
                                              linres_write_restart
   USE qs_linres_types,                 ONLY: get_polar_env,&
                                              linres_control_type,&
                                              polar_env_create,&
                                              polar_env_type,&
                                              set_polar_env
   USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_p_env_types,                  ONLY: qs_p_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: polar_env_init, polar_polar, polar_print, polar_response, write_polarisability_tensor

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_polar_utils'

CONTAINS

! **************************************************************************************************
!> \brief Initialize the polar environment
!> \param qs_env ...
!> \par History
!>      06.2018 polar_env integrated into qs_env (MK)
! **************************************************************************************************
   SUBROUTINE polar_env_init(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'polar_env_init'

      INTEGER                                            :: handle, idir, iounit, ispin, m, nao, &
                                                            nmo, nspins
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(polar_env_type), POINTER                      :: polar_env
      TYPE(section_vals_type), POINTER                   :: lr_section, polar_section

      CALL timeset(routineN, handle)

      NULLIFY (dft_control)
      NULLIFY (linres_control)
      NULLIFY (logger)
      NULLIFY (matrix_s)
      NULLIFY (mos)
      NULLIFY (polar_env)
      NULLIFY (lr_section, polar_section)

      logger => cp_get_default_logger()
      lr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES")

      iounit = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", &
                                    extension=".linresLog")

      IF (iounit > 0) THEN
         WRITE (iounit, "(/,(T2,A))") "POLAR| Starting polarizability calculation", &
            "POLAR| Initialization of the polar environment"
      ENDIF

      polar_section => section_vals_get_subs_vals(qs_env%input, &
                                                  "PROPERTIES%LINRES%POLAR")

      CALL get_qs_env(qs_env=qs_env, &
                      polar_env=polar_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      linres_control=linres_control, &
                      mos=mos)

      ! Create polar environment if needed
      IF (.NOT. ASSOCIATED(polar_env)) THEN
         CALL polar_env_create(polar_env)
         CALL set_qs_env(qs_env=qs_env, polar_env=polar_env)
      ENDIF

      nspins = dft_control%nspins

      CALL section_vals_val_get(polar_section, "DO_RAMAN", l_val=polar_env%do_raman)
      CALL section_vals_val_get(polar_section, "PERIODIC_DIPOLE_OPERATOR", l_val=polar_env%do_periodic)

      ! Allocate components of the polar environment if needed
      IF (.NOT. ASSOCIATED(polar_env%polar)) THEN
         ALLOCATE (polar_env%polar(3, 3))
         polar_env%polar(:, :) = 0.0_dp
      ENDIF
      IF (.NOT. ASSOCIATED(polar_env%dBerry_psi0)) THEN
         ALLOCATE (polar_env%dBerry_psi0(3, nspins))
      ELSE
         ! Remove previous matrices
         DO ispin = 1, nspins
            DO idir = 1, 3
               CALL cp_fm_release(polar_env%dBerry_psi0(idir, ispin)%matrix)
            ENDDO
         ENDDO
      ENDIF
      IF (.NOT. ASSOCIATED(polar_env%psi1_dBerry)) THEN
         ALLOCATE (polar_env%psi1_dBerry(3, nspins))
      ELSE
         ! Remove previous matrices
         DO ispin = 1, nspins
            DO idir = 1, 3
               CALL cp_fm_release(polar_env%psi1_dBerry(idir, ispin)%matrix)
            ENDDO
         ENDDO
      ENDIF
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao, nmo=nmo)
         CALL cp_fm_get_info(mo_coeff, ncol_global=m, nrow_global=nao)
         NULLIFY (tmp_fm_struct)
         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nao, &
                                  ncol_global=m, &
                                  context=mo_coeff%matrix_struct%context)
         DO idir = 1, 3
            CALL cp_fm_create(polar_env%dBerry_psi0(idir, ispin)%matrix, tmp_fm_struct)
            CALL cp_fm_create(polar_env%psi1_dBerry(idir, ispin)%matrix, tmp_fm_struct)
         ENDDO
         CALL cp_fm_struct_release(tmp_fm_struct)
      END DO

      CALL cp_print_key_finished_output(iounit, logger, lr_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      CALL timestop(handle)

   END SUBROUTINE polar_env_init

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \par History
!>      06.2018 polar_env integrated into qs_env (MK)
! **************************************************************************************************
   SUBROUTINE polar_polar(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'polar_polar'

      INTEGER                                            :: handle, i, iounit, ispin, nspins, z
      LOGICAL                                            :: do_periodic, do_raman, run_stopped
      REAL(dp)                                           :: ptmp
      REAL(dp), DIMENSION(:, :), POINTER                 :: polar, polar_tmp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: dBerry_psi0, psi1_dBerry
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(polar_env_type), POINTER                      :: polar_env

      CALL timeset(routineN, handle)

      NULLIFY (cell, dft_control, polar, psi1_dBerry, logger)
      NULLIFY (mos, dBerry_psi0)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      mos=mos, &
                      polar_env=polar_env)

      nspins = dft_control%nspins

      CALL get_polar_env(polar_env=polar_env, &
                         do_raman=do_raman, &
                         run_stopped=run_stopped)

      IF (.NOT. run_stopped .AND. do_raman) THEN

         CALL cite_reference(Luber2014)

         CALL get_polar_env(polar_env=polar_env, &
                            do_periodic=do_periodic, &
                            dBerry_psi0=dBerry_psi0, &
                            polar=polar, &
                            psi1_dBerry=psi1_dBerry)

         ! Initialize
         ALLOCATE (polar_tmp(3, 3))
         polar_tmp(:, :) = 0.0_dp

         DO i = 1, 3 ! directions of electric field
            DO z = 1, 3 !dipole directions
               DO ispin = 1, dft_control%nspins
                  !SL compute trace
                  ptmp = 0.0_dp
                  CALL cp_fm_trace(psi1_dBerry(i, ispin)%matrix, dBerry_psi0(z, ispin)%matrix, ptmp)
                  polar_tmp(i, z) = polar_tmp(i, z) - 2.0_dp*ptmp
               END DO
            END DO
         END DO !spin

         IF (do_periodic) THEN
            polar(:, :) = MATMUL(MATMUL(cell%hmat, polar_tmp), TRANSPOSE(cell%hmat))/(twopi*twopi)
         ELSE
            polar(:, :) = polar_tmp(:, :)
         END IF
         !SL evtl maxocc instead?
         IF (dft_control%nspins == 1) THEN
            polar(:, :) = 2.0_dp*polar(:, :)
         END IF

         IF (ASSOCIATED(polar_tmp)) THEN
            DEALLOCATE (polar_tmp)
         END IF

      ENDIF ! do_raman

      CALL timestop(handle)

   END SUBROUTINE polar_polar

! **************************************************************************************************
!> \brief Print information related to the polarisability tensor
!> \param qs_env ...
!> \par History
!>      06.2018 polar_env integrated into qs_env (MK)
! **************************************************************************************************
   SUBROUTINE polar_print(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'polar_print'

      CHARACTER(LEN=default_string_length)               :: description
      INTEGER                                            :: iounit, unit_p
      LOGICAL                                            :: do_raman, run_stopped
      REAL(dp), DIMENSION(:, :), POINTER                 :: polar
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(polar_env_type), POINTER                      :: polar_env
      TYPE(section_vals_type), POINTER                   :: polar_section

      NULLIFY (logger, dft_control, para_env, results)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      polar_env=polar_env, &
                      results=results, &
                      para_env=para_env)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      polar_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES%POLAR")

      CALL get_polar_env(polar_env=polar_env, &
                         polar=polar, &
                         do_raman=do_raman, &
                         run_stopped=run_stopped)

      IF (.NOT. run_stopped .AND. do_raman) THEN

         description = "[POLAR]"
         CALL cp_results_erase(results, description=description)
         CALL put_results(results, description=description, values=polar(:, :))

         IF (BTEST(cp_print_key_should_output(logger%iter_info, polar_section, &
                                              "PRINT%POLAR_MATRIX"), cp_p_file)) THEN

            unit_p = cp_print_key_unit_nr(logger, polar_section, "PRINT%POLAR_MATRIX", &
                                          extension=".data", middle_name="raman", log_filename=.FALSE.)
            IF (unit_p > 0) THEN
               IF (unit_p /= iounit) THEN
                  WRITE (unit_p, *)
                  WRITE (unit_p, '(T10,A)') 'POLARIZABILITY TENSOR (atomic units):'
                  WRITE (unit_p, '(T10,A,3F15.5)') "xx,yy,zz", polar(1, 1), polar(2, 2), polar(3, 3)
                  WRITE (unit_p, '(T10,A,3F15.5)') "xy,xz,yz", polar(1, 2), polar(1, 3), polar(2, 3)
                  WRITE (unit_p, '(T10,A,3F15.5)') "yx,zx,zy", polar(2, 1), polar(3, 1), polar(3, 2)
                  WRITE (unit_p, '(T10,A)') 'POLARIZABILITY TENSOR (Angstrom^3):'
                  WRITE (unit_p, '(T10,A,3F15.5)') "xx,yy,zz", polar(1, 1)*angstrom**3, &
                     polar(2, 2)*angstrom**3, polar(3, 3)*angstrom**3
                  WRITE (unit_p, '(T10,A,3F15.5)') "xy,xz,yz", polar(1, 2)*angstrom**3, &
                     polar(1, 3)*angstrom**3, polar(2, 3)*angstrom**3
                  WRITE (unit_p, '(T10,A,3F15.5)') "yx,zx,zy", polar(2, 1)*angstrom**3, &
                     polar(3, 1)*angstrom**3, polar(3, 2)*angstrom**3
                  CALL cp_print_key_finished_output(unit_p, logger, polar_section, &
                                                    "PRINT%POLAR_MATRIX")
               ENDIF
            ENDIF
         ENDIF
         IF (iounit > 0) THEN
            WRITE (iounit, *)
            WRITE (iounit, '(T2,A)') 'POLAR| POLARIZABILITY TENSOR (atomic units):'
            WRITE (iounit, '(T2,A,3F15.5)') "POLAR| xx,yy,zz", polar(1, 1), polar(2, 2), polar(3, 3)
            WRITE (iounit, '(T2,A,3F15.5)') "POLAR| xy,xz,yz", polar(1, 2), polar(1, 3), polar(2, 3)
            WRITE (iounit, '(T2,A,3F15.5)') "POLAR| yx,zx,zy", polar(2, 1), polar(3, 1), polar(3, 2)
            WRITE (iounit, '(T2,A)') 'POLAR| POLARIZABILITY TENSOR (Angstrom^3):'
            WRITE (iounit, '(T2,A,3F15.5)') "POLAR| xx,yy,zz", polar(1, 1)*angstrom**3, &
               polar(2, 2)*angstrom**3, polar(3, 3)*angstrom**3
            WRITE (iounit, '(T2,A,3F15.5)') "POLAR| xy,xz,yz", polar(1, 2)*angstrom**3, &
               polar(1, 3)*angstrom**3, polar(2, 3)*angstrom**3
            WRITE (iounit, '(T2,A,3F15.5)') "POLAR| yx,zx,zy", polar(2, 1)*angstrom**3, &
               polar(3, 1)*angstrom**3, polar(3, 2)*angstrom**3
         END IF
      END IF

   END SUBROUTINE polar_print

! **************************************************************************************************
!> \brief Calculate the polarisability tensor using response theory
!> \param p_env ...
!> \param qs_env ...
!> \par History
!>      06.2018 polar_env integrated into qs_env (MK)
! **************************************************************************************************
   SUBROUTINE polar_response(p_env, qs_env)

      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'polar_response'

      INTEGER                                            :: handle, idir, iounit, ispin, nao, nmo, &
                                                            nspins
      LOGICAL                                            :: do_periodic, do_raman, should_stop
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: h1_psi0, psi0_order, psi1, psi1_ptr
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: dBerry_psi0, psi1_dBerry
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(polar_env_type), POINTER                      :: polar_env
      TYPE(qs_matrix_pools_type), POINTER                :: mpools
      TYPE(section_vals_type), POINTER                   :: lr_section, polar_section

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, linres_control, lr_section, polar_section)
      NULLIFY (logger, mpools, psi1, h1_psi0, mo_coeff, para_env)
      NULLIFY (tmp_fm_struct, psi1_dBerry, dBerry_psi0)

      logger => cp_get_default_logger()
      lr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES")
      polar_section => section_vals_get_subs_vals(qs_env%input, &
                                                  "PROPERTIES%LINRES%POLAR")

      iounit = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", &
                                    extension=".linresLog")
      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(T2,A,/)") &
            "POLAR| Self consistent optimization of the response wavefunctions"
      ENDIF

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      mpools=mpools, &
                      linres_control=linres_control, &
                      mos=mos, &
                      polar_env=polar_env, &
                      para_env=para_env)

      nspins = dft_control%nspins

      CALL get_polar_env(polar_env=polar_env, do_raman=do_raman, do_periodic=do_periodic)

      ! Allocate the vectors
      ALLOCATE (psi0_order(nspins))
      ALLOCATE (psi1(nspins), h1_psi0(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
         psi0_order(ispin)%matrix => mo_coeff
         CALL cp_fm_get_info(mo_coeff, ncol_global=nmo, nrow_global=nao)
         NULLIFY (tmp_fm_struct, psi1(ispin)%matrix, h1_psi0(ispin)%matrix)
         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nao, &
                                  ncol_global=nmo, &
                                  context=mo_coeff%matrix_struct%context)
         CALL cp_fm_create(psi1(ispin)%matrix, tmp_fm_struct)
         CALL cp_fm_create(h1_psi0(ispin)%matrix, tmp_fm_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
      ENDDO

      IF (do_raman) THEN
         CALL get_polar_env(polar_env=polar_env, &
                            psi1_dBerry=psi1_dBerry, &
                            dBerry_psi0=dBerry_psi0)
         ! Restart
         IF (linres_control%linres_restart) THEN
            DO idir = 1, 3
               psi1_ptr => psi1_dBerry(idir, :)
               CALL linres_read_restart(qs_env, lr_section, psi1_ptr, idir, "psi1_dBerry")
            ENDDO
         ELSE
            DO idir = 1, 3
               DO ispin = 1, nspins
                  CALL cp_fm_set_all(psi1_dBerry(idir, ispin)%matrix, 0.0_dp)
               ENDDO
            ENDDO
         ENDIF
         loop_idir: DO idir = 1, 3
            IF (iounit > 0) THEN
               IF (do_periodic) THEN
                  WRITE (iounit, "(/,T2,A)") &
                     "POLAR| Response to the perturbation operator Berry phase_"//ACHAR(idir + 119)
               ELSE
                  WRITE (iounit, "(/,T2,A)") &
                     "POLAR| Response to the perturbation operator R_"//ACHAR(idir + 119)
               END IF
            ENDIF
            ! Do scf cycle to optimize psi1
            DO ispin = 1, nspins
               CALL cp_fm_to_fm(psi1_dBerry(idir, ispin)%matrix, psi1(ispin)%matrix)
               CALL cp_fm_to_fm(dBerry_psi0(idir, ispin)%matrix, h1_psi0(ispin)%matrix)
            ENDDO
            !
            linres_control%lr_triplet = .FALSE. ! we do singlet response
            linres_control%do_kernel = .TRUE. ! we do coupled response
            linres_control%converged = .FALSE.
            CALL linres_solver(p_env, qs_env, psi1, h1_psi0, psi0_order, iounit, should_stop)

            ! Copy the response
            DO ispin = 1, nspins
               CALL cp_fm_to_fm(psi1(ispin)%matrix, psi1_dBerry(idir, ispin)%matrix)
            ENDDO
            !
            ! Write the new result to the restart file
            IF (linres_control%linres_restart) THEN
               psi1_ptr => psi1_dBerry(idir, :)
               CALL linres_write_restart(qs_env, lr_section, psi1_ptr, idir, "psi1_dBerry")
            ENDIF
         ENDDO loop_idir
      ENDIF ! do_raman

      CALL set_polar_env(polar_env, run_stopped=should_stop)

      ! Clean up
      DO ispin = 1, nspins
         CALL cp_fm_release(psi1(ispin)%matrix)
         CALL cp_fm_release(h1_psi0(ispin)%matrix)
      ENDDO
      DEALLOCATE (psi1, h1_psi0, psi0_order)

      CALL cp_print_key_finished_output(iounit, logger, lr_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      CALL timestop(handle)

   END SUBROUTINE polar_response

! **************************************************************************************************
!> \brief Prints the polarisability tensor to a file during MD runs
!> \param force_env ...
!> \param motion_section ...
!> \param itimes ...
!> \param time ...
!> \param pos ...
!> \param act ...
!> \par History
!>      06.2018 Creation (MK)
!> \author Matthias Krack (MK)
! **************************************************************************************************
   SUBROUTINE write_polarisability_tensor(force_env, motion_section, itimes, time, pos, act)

      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(section_vals_type), POINTER                   :: motion_section
      INTEGER, INTENT(IN)                                :: itimes
      REAL(KIND=dp), INTENT(IN)                          :: time
      CHARACTER(LEN=default_string_length), INTENT(IN), &
         OPTIONAL                                        :: pos, act

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_polarisability_tensor', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: my_act, my_pos
      INTEGER                                            :: iounit
      LOGICAL                                            :: do_raman, new_file, run_stopped
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: polar
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(polar_env_type), POINTER                      :: polar_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      NULLIFY (qs_env)

      CALL force_env_get(force_env, qs_env=qs_env)
      IF (ASSOCIATED(qs_env)) THEN
         NULLIFY (polar_env)
         CALL get_qs_env(qs_env=qs_env, polar_env=polar_env)
         IF (ASSOCIATED(polar_env)) THEN
            CALL get_polar_env(polar_env=polar_env, &
                               polar=polar, &
                               do_raman=do_raman, &
                               run_stopped=run_stopped)
            IF (.NOT. run_stopped .AND. do_raman) THEN
               NULLIFY (logger)
               logger => cp_get_default_logger()
               my_pos = "APPEND"
               my_act = "WRITE"
               IF (PRESENT(pos)) my_pos = pos
               IF (PRESENT(act)) my_act = act
               iounit = cp_print_key_unit_nr(logger, motion_section, "PRINT%POLAR_MATRIX", &
                                             extension=".polar", file_position=my_pos, &
                                             file_action=my_act, file_form="FORMATTED", &
                                             is_new_file=new_file)
            ELSE
               iounit = 0
            END IF
            IF (iounit > 0) THEN
               IF (new_file) THEN
                  WRITE (UNIT=iounit, FMT='(A,9(11X,A2," [a.u.]"),6X,A)') &
                     "#   Step   Time [fs]", "xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"
               END IF
               WRITE (UNIT=iounit, FMT='(I8,F12.3,9(1X,F19.8))') itimes, time, &
                  polar(1, 1), polar(1, 2), polar(1, 3), &
                  polar(2, 1), polar(2, 2), polar(2, 3), &
                  polar(3, 1), polar(3, 2), polar(3, 3)
               CALL m_flush(iounit)
               CALL cp_print_key_finished_output(iounit, logger, motion_section, "PRINT%POLAR_MATRIX")
            END IF
         END IF ! polar_env
      END IF ! qs_env

   END SUBROUTINE write_polarisability_tensor

END MODULE qs_linres_polar_utils
