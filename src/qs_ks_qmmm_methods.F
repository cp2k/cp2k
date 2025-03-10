!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \author T.Laino
! **************************************************************************************************
MODULE qs_ks_qmmm_methods
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cube_utils,                      ONLY: cube_info_type,&
                                              init_cube_info
   USE d3_poly,                         ONLY: init_d3_poly_module
   USE input_constants,                 ONLY: do_qmmm_gauss,&
                                              do_qmmm_swave
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_retain
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integral_ab
   USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qmmm_types_low,                  ONLY: qmmm_env_qm_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_qmmm_methods'

   PUBLIC :: ks_qmmm_env_rebuild, qmmm_calculate_energy, &
             qmmm_modify_hartree_pot

CONTAINS

! **************************************************************************************************
!> \brief Initialize the ks_qmmm_env
!> \param qs_env ...
!> \param qmmm_env ...
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE ks_qmmm_env_rebuild(qs_env, qmmm_env)
      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env
      TYPE(qmmm_env_qm_type), POINTER                    :: qmmm_env

      TYPE(qs_ks_qmmm_env_type), POINTER                 :: ks_qmmm_env

      NULLIFY (ks_qmmm_env)
      CALL get_qs_env(qs_env=qs_env, &
                      ks_qmmm_env=ks_qmmm_env)

      !   *** allocate the ks_qmmm env if not allocated yet!**
      IF (.NOT. ASSOCIATED(ks_qmmm_env)) THEN
         ALLOCATE (ks_qmmm_env)
         CALL qs_ks_qmmm_create(ks_qmmm_env=ks_qmmm_env, qs_env=qs_env, &
                                qmmm_env=qmmm_env)
         CALL set_qs_env(qs_env=qs_env, ks_qmmm_env=ks_qmmm_env)
      END IF
   END SUBROUTINE ks_qmmm_env_rebuild

! **************************************************************************************************
!> \brief allocates and initializes the given ks_qmmm_env.
!> \param ks_qmmm_env the ks_qmmm env to be initialized
!> \param qs_env the qs environment
!> \param qmmm_env ...
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE qs_ks_qmmm_create(ks_qmmm_env, qs_env, qmmm_env)
      TYPE(qs_ks_qmmm_env_type), INTENT(OUT)             :: ks_qmmm_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qmmm_env_qm_type), POINTER                    :: qmmm_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_ks_qmmm_create'

      INTEGER                                            :: handle, igrid
      TYPE(cube_info_type), DIMENSION(:), POINTER        :: cube_info
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      NULLIFY (ks_qmmm_env%pw_env, &
               ks_qmmm_env%cube_info)
      NULLIFY (auxbas_pw_pool)
      CALL get_qs_env(qs_env=qs_env, &
                      pw_env=ks_qmmm_env%pw_env)
      CALL pw_env_get(ks_qmmm_env%pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL pw_env_retain(ks_qmmm_env%pw_env)

      ks_qmmm_env%pc_ener = 0.0_dp
      ks_qmmm_env%n_evals = 0

      CALL auxbas_pw_pool%create_pw(ks_qmmm_env%v_qmmm_rspace)

      IF ((qmmm_env%qmmm_coupl_type == do_qmmm_gauss) .OR. (qmmm_env%qmmm_coupl_type == do_qmmm_swave)) THEN
         CALL init_d3_poly_module() ! a fairly arbitrary but sufficient spot to do this
         CALL pw_env_get(ks_qmmm_env%pw_env, pw_pools=pools)
         ALLOCATE (cube_info(SIZE(pools)))
         DO igrid = 1, SIZE(pools)
            CALL init_cube_info(cube_info(igrid), &
                                pools(igrid)%pool%pw_grid%dr(:), &
                                pools(igrid)%pool%pw_grid%dh(:, :), &
                                pools(igrid)%pool%pw_grid%dh_inv(:, :), &
                                pools(igrid)%pool%pw_grid%orthorhombic, &
                                qmmm_env%maxRadius(igrid))
         END DO
         ks_qmmm_env%cube_info => cube_info
      END IF
      NULLIFY (ks_qmmm_env%matrix_h)
      !

      CALL timestop(handle)

   END SUBROUTINE qs_ks_qmmm_create

! **************************************************************************************************
!> \brief Computes the contribution to the total energy of the QM/MM
!>      electrostatic coupling
!> \param qs_env ...
!> \param rho ...
!> \param v_qmmm ...
!> \param qmmm_energy ...
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE qmmm_calculate_energy(qs_env, rho, v_qmmm, qmmm_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), DIMENSION(:), INTENT(IN)     :: rho
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: v_qmmm
      REAL(KIND=dp), INTENT(INOUT)                       :: qmmm_energy

      CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_calculate_energy'

      INTEGER                                            :: handle, ispin, output_unit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_r3d_rs_type), POINTER                      :: rho0_s_rs
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)
      NULLIFY (dft_control, input, logger)
      logger => cp_get_default_logger()

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      input=input)

      output_unit = cp_print_key_unit_nr(logger, input, "QMMM%PRINT%PROGRAM_RUN_INFO", &
                                         extension=".qmmmLog")
      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(T3,A)") &
         "Adding QM/MM electrostatic potential to the Kohn-Sham potential."

      qmmm_energy = 0.0_dp
      DO ispin = 1, dft_control%nspins
         qmmm_energy = qmmm_energy + pw_integral_ab(rho(ispin), v_qmmm)
      END DO
      IF (dft_control%qs_control%gapw) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         rho0_s_rs=rho0_s_rs)
         CPASSERT(ASSOCIATED(rho0_s_rs))
         qmmm_energy = qmmm_energy + pw_integral_ab(rho0_s_rs, v_qmmm)
      END IF

      CALL cp_print_key_finished_output(output_unit, logger, input, &
                                        "QMMM%PRINT%PROGRAM_RUN_INFO")

      CALL timestop(handle)
   END SUBROUTINE qmmm_calculate_energy

! **************************************************************************************************
!> \brief Modify the hartree potential in order to include the QM/MM correction
!> \param v_hartree ...
!> \param v_qmmm ...
!> \param scale ...
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE qmmm_modify_hartree_pot(v_hartree, v_qmmm, scale)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: v_hartree
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: v_qmmm
      REAL(KIND=dp), INTENT(IN)                          :: scale

      CALL pw_axpy(v_qmmm, v_hartree, v_qmmm%pw_grid%dvol*scale)

   END SUBROUTINE qmmm_modify_hartree_pot

END MODULE qs_ks_qmmm_methods
