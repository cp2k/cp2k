!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief pw_types
!> \author CJM
! **************************************************************************************************
MODULE ewald_pw_types
   USE ao_util,                         ONLY: exp_radius
   USE cell_types,                      ONLY: cell_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_init,          ONLY: init_input_type
   USE dg_types,                        ONLY: dg_create,&
                                              dg_release,&
                                              dg_type
   USE dgs,                             ONLY: dg_pme_grid_setup
   USE ewald_environment_types,         ONLY: ewald_env_get,&
                                              ewald_environment_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_comm_self,&
                                              mp_para_env_type
   USE pw_grid_types,                   ONLY: HALFSPACE,&
                                              pw_grid_type
   USE pw_grids,                        ONLY: pw_grid_create,&
                                              pw_grid_release,&
                                              pw_grid_setup
   USE pw_poisson_methods,              ONLY: pw_poisson_set
   USE pw_poisson_read_input,           ONLY: pw_poisson_read_parameters
   USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                              do_ewald_none,&
                                              do_ewald_pme,&
                                              do_ewald_spme,&
                                              pw_poisson_parameter_type,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create,&
                                              pw_pool_p_type,&
                                              pw_pool_release,&
                                              pw_pool_type
   USE realspace_grid_types,            ONLY: &
        realspace_grid_desc_type, realspace_grid_input_type, realspace_grid_type, rs_grid_create, &
        rs_grid_create_descriptor, rs_grid_print, rs_grid_release, rs_grid_release_descriptor, &
        rs_grid_retain_descriptor
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewald_pw_types'
   PUBLIC :: ewald_pw_type, ewald_pw_release, &
             ewald_pw_create, &
             ewald_pw_get, ewald_pw_set

! **************************************************************************************************
   TYPE ewald_pw_type
      PRIVATE
      TYPE(pw_pool_type), POINTER       :: pw_small_pool
      TYPE(pw_pool_type), POINTER       :: pw_big_pool
      TYPE(realspace_grid_desc_type), POINTER    :: rs_desc
      TYPE(pw_poisson_type), POINTER    :: poisson_env
      TYPE(dg_type), POINTER            :: dg
   END TYPE ewald_pw_type

CONTAINS

! **************************************************************************************************
!> \brief creates the structure ewald_pw_type
!> \param ewald_pw ...
!> \param ewald_env ...
!> \param cell ...
!> \param cell_ref ...
!> \param print_section ...
! **************************************************************************************************
   SUBROUTINE ewald_pw_create(ewald_pw, ewald_env, cell, cell_ref, print_section)
      TYPE(ewald_pw_type), INTENT(OUT)                   :: ewald_pw
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(cell_type), POINTER                           :: cell, cell_ref
      TYPE(section_vals_type), POINTER                   :: print_section

      NULLIFY (ewald_pw%pw_big_pool)
      NULLIFY (ewald_pw%pw_small_pool)
      NULLIFY (ewald_pw%rs_desc)
      NULLIFY (ewald_pw%poisson_env)
      ALLOCATE (ewald_pw%dg)
      CALL dg_create(ewald_pw%dg)
      CALL ewald_pw_init(ewald_pw, ewald_env, cell, cell_ref, print_section)
   END SUBROUTINE ewald_pw_create

! **************************************************************************************************
!> \brief releases the memory used by the ewald_pw
!> \param ewald_pw ...
! **************************************************************************************************
   SUBROUTINE ewald_pw_release(ewald_pw)
      TYPE(ewald_pw_type), INTENT(INOUT)                 :: ewald_pw

      CALL pw_pool_release(ewald_pw%pw_small_pool)
      CALL pw_pool_release(ewald_pw%pw_big_pool)
      CALL rs_grid_release_descriptor(ewald_pw%rs_desc)
      IF (ASSOCIATED(ewald_pw%poisson_env)) THEN
         CALL ewald_pw%poisson_env%release()
         DEALLOCATE (ewald_pw%poisson_env)
      END IF
      CALL dg_release(ewald_pw%dg)
      DEALLOCATE (ewald_pw%dg)

   END SUBROUTINE ewald_pw_release

! **************************************************************************************************
!> \brief ...
!> \param ewald_pw ...
!> \param ewald_env ...
!> \param cell ...
!> \param cell_ref ...
!> \param print_section ...
!> \par History
!>      JGH (12-Jan-2001): Added SPME part
!>      JGH (15-Mar-2001): Work newly distributed between initialize, setup,
!>                         and force routine
!> \author CJM
! **************************************************************************************************
   SUBROUTINE ewald_pw_init(ewald_pw, ewald_env, cell, cell_ref, print_section)
      TYPE(ewald_pw_type), INTENT(INOUT)                 :: ewald_pw
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(cell_type), POINTER                           :: cell, cell_ref
      TYPE(section_vals_type), POINTER                   :: print_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ewald_pw_init'

      INTEGER                                            :: bo(2, 3), ewald_type, gmax(3), handle, &
                                                            npts_s(3), ns_max, o_spline, &
                                                            output_unit
      REAL(KIND=dp)                                      :: alpha, alphasq, cutoff_radius, epsilon, &
                                                            norm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_grid_type), POINTER                        :: pw_big_grid, pw_small_grid
      TYPE(pw_poisson_parameter_type)                    :: poisson_params
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(realspace_grid_desc_type), POINTER            :: rs_desc
      TYPE(realspace_grid_input_type)                    :: input_settings
      TYPE(section_vals_type), POINTER                   :: poisson_section, rs_grid_section

      CALL timeset(routineN, handle)

      NULLIFY (pw_big_grid)
      NULLIFY (pw_small_grid, poisson_section)

      CPASSERT(ASSOCIATED(ewald_env))
      CPASSERT(ASSOCIATED(cell))
      CALL ewald_env_get(ewald_env=ewald_env, &
                         para_env=para_env, &
                         gmax=gmax, alpha=alpha, &
                         ns_max=ns_max, &
                         ewald_type=ewald_type, &
                         o_spline=o_spline, &
                         poisson_section=poisson_section, &
                         epsilon=epsilon)

      rs_grid_section => section_vals_get_subs_vals(poisson_section, "EWALD%RS_GRID")

      SELECT CASE (ewald_type)
      CASE (do_ewald_ewald)
         ! set up Classic EWALD sum
         logger => cp_get_default_logger()
         output_unit = cp_print_key_unit_nr(logger, print_section, "", extension=".Log")
         CALL pw_grid_create(pw_big_grid, mp_comm_self)

         IF (ANY(gmax == 2*(gmax/2))) THEN
            CPABORT("gmax has to be odd.")
         END IF
         bo(1, :) = -gmax/2
         bo(2, :) = +gmax/2
         CALL pw_grid_setup(cell_ref%hmat, pw_big_grid, grid_span=HALFSPACE, bounds=bo, spherical=.TRUE., &
                            fft_usage=.FALSE., iounit=output_unit)
         NULLIFY (pw_pool)
         CALL pw_pool_create(pw_pool, pw_grid=pw_big_grid)
         ewald_pw%pw_big_pool => pw_pool
         CALL pw_grid_release(pw_big_grid)
         CALL cp_print_key_finished_output(output_unit, logger, print_section, "")

      CASE (do_ewald_pme)
         ! set up Particle-Mesh EWALD sum
         logger => cp_get_default_logger()
         output_unit = cp_print_key_unit_nr(logger, print_section, "", extension=".Log")
         IF (.NOT. ASSOCIATED(ewald_pw%poisson_env)) THEN
            ALLOCATE (ewald_pw%poisson_env)
            CALL ewald_pw%poisson_env%create()
         END IF
         CALL pw_grid_create(pw_small_grid, mp_comm_self)
         CALL pw_grid_create(pw_big_grid, para_env)
         IF (ns_max == 2*(ns_max/2)) THEN
            CPABORT("ns_max has to be odd.")
         END IF
         npts_s(:) = ns_max
         ! compute cut-off radius
         alphasq = alpha**2
         norm = (2.0_dp*alphasq/pi)**(1.5_dp)
         cutoff_radius = exp_radius(0, 2.0_dp*alphasq, epsilon, norm)

         CALL dg_pme_grid_setup(cell_ref%hmat, npts_s, cutoff_radius, &
                                pw_small_grid, pw_big_grid, rs_dims=(/para_env%num_pe, 1/), &
                                iounit=output_unit, fft_usage=.TRUE.)
         ! Write some useful info
         IF (output_unit > 0) THEN
            WRITE (output_unit, '( A,T71,E10.4 )') &
               ' EWALD| Gaussian tolerance (effective) ', epsilon
            WRITE (output_unit, '( A,T63,3I6 )') &
               ' EWALD| Small box grid ', pw_small_grid%npts
            WRITE (output_unit, '( A,T63,3I6 )') &
               ' EWALD| Full box grid ', pw_big_grid%npts
         END IF

         ! pw pools initialized
         NULLIFY (pw_pool)
         CALL pw_pool_create(pw_pool, pw_grid=pw_big_grid)
         ewald_pw%pw_big_pool => pw_pool

         NULLIFY (pw_pool)
         CALL pw_pool_create(pw_pool, pw_grid=pw_small_grid)
         ewald_pw%pw_small_pool => pw_pool

         NULLIFY (rs_desc)
         CALL init_input_type(input_settings, nsmax=MAXVAL(pw_small_grid%npts(1:3)), &
                              rs_grid_section=rs_grid_section, ilevel=1, &
                              higher_grid_layout=(/-1, -1, -1/))
         CALL rs_grid_create_descriptor(rs_desc, pw_big_grid, input_settings)

         BLOCK
            TYPE(realspace_grid_type) :: rs
            CALL rs_grid_create(rs, rs_desc)
            CALL rs_grid_print(rs, output_unit)
            CALL rs_grid_release(rs)
         END BLOCK

         CALL cp_print_key_finished_output(output_unit, logger, print_section, "")

         ewald_pw%rs_desc => rs_desc

         CALL rs_grid_retain_descriptor(ewald_pw%rs_desc)
         CALL rs_grid_release_descriptor(rs_desc)

         CALL pw_grid_release(pw_small_grid)
         CALL pw_grid_release(pw_big_grid)

      CASE (do_ewald_spme)
         ! set up the Smooth-Particle-Mesh EWALD sum
         logger => cp_get_default_logger()
         output_unit = cp_print_key_unit_nr(logger, print_section, "", extension=".Log")
         IF (.NOT. ASSOCIATED(ewald_pw%poisson_env)) THEN
            ALLOCATE (ewald_pw%poisson_env)
            CALL ewald_pw%poisson_env%create()
         END IF
         CALL pw_grid_create(pw_big_grid, para_env)
         npts_s = gmax
         CALL pw_grid_setup(cell_ref%hmat, pw_big_grid, grid_span=HALFSPACE, npts=npts_s, spherical=.TRUE., &
                            rs_dims=(/para_env%num_pe, 1/), iounit=output_unit, fft_usage=.TRUE.)

         ! pw pools initialized
         NULLIFY (pw_pool)
         CALL pw_pool_create(pw_pool, pw_grid=pw_big_grid)
         ewald_pw%pw_big_pool => pw_pool

         NULLIFY (rs_desc)
         CALL init_input_type(input_settings, nsmax=o_spline, &
                              rs_grid_section=rs_grid_section, ilevel=1, &
                              higher_grid_layout=(/-1, -1, -1/))
         CALL rs_grid_create_descriptor(rs_desc, pw_big_grid, input_settings)

         BLOCK
            TYPE(realspace_grid_type) :: rs

            CALL rs_grid_create(rs, rs_desc)
            CALL rs_grid_print(rs, output_unit)
            CALL rs_grid_release(rs)
         END BLOCK
         CALL cp_print_key_finished_output(output_unit, logger, print_section, "")

         ewald_pw%rs_desc => rs_desc

         CALL rs_grid_retain_descriptor(ewald_pw%rs_desc)
         CALL rs_grid_release_descriptor(rs_desc)

         CALL pw_grid_release(pw_big_grid)
      CASE (do_ewald_none)
         ! No EWALD sums..
      CASE default
         CPABORT("")
      END SELECT
      ! Poisson Environment
      IF (ASSOCIATED(ewald_pw%poisson_env)) THEN
         ALLOCATE (pw_pools(1))
         pw_pools(1)%pool => ewald_pw%pw_big_pool
         CALL pw_poisson_read_parameters(poisson_section, poisson_params)
         poisson_params%ewald_type = ewald_type
         poisson_params%ewald_o_spline = o_spline
         poisson_params%ewald_alpha = alpha
         CALL pw_poisson_set(ewald_pw%poisson_env, cell_hmat=cell%hmat, parameters=poisson_params, &
                             use_level=1, pw_pools=pw_pools)
         DEALLOCATE (pw_pools)
      END IF
      CALL timestop(handle)
   END SUBROUTINE ewald_pw_init

! **************************************************************************************************
!> \brief get the ewald_pw environment to the correct program.
!> \param ewald_pw ...
!> \param pw_big_pool ...
!> \param pw_small_pool ...
!> \param rs_desc ...
!> \param poisson_env ...
!> \param dg ...
!> \author CJM
! **************************************************************************************************
   SUBROUTINE ewald_pw_get(ewald_pw, pw_big_pool, pw_small_pool, rs_desc, poisson_env, dg)

      TYPE(ewald_pw_type), INTENT(IN)                    :: ewald_pw
      TYPE(pw_pool_type), OPTIONAL, POINTER              :: pw_big_pool, pw_small_pool
      TYPE(realspace_grid_desc_type), OPTIONAL, POINTER  :: rs_desc
      TYPE(pw_poisson_type), OPTIONAL, POINTER           :: poisson_env
      TYPE(dg_type), OPTIONAL, POINTER                   :: dg

      IF (PRESENT(poisson_env)) poisson_env => ewald_pw%poisson_env
      IF (PRESENT(pw_big_pool)) pw_big_pool => ewald_pw%pw_big_pool
      IF (PRESENT(pw_small_pool)) pw_small_pool => ewald_pw%pw_small_pool
      IF (PRESENT(rs_desc)) rs_desc => ewald_pw%rs_desc
      IF (PRESENT(dg)) dg => ewald_pw%dg

   END SUBROUTINE ewald_pw_get

! **************************************************************************************************
!> \brief set the ewald_pw environment to the correct program.
!> \param ewald_pw ...
!> \param pw_big_pool ...
!> \param pw_small_pool ...
!> \param rs_desc ...
!> \param dg ...
!> \param poisson_env ...
!> \author CJM
! **************************************************************************************************
   SUBROUTINE ewald_pw_set(ewald_pw, pw_big_pool, pw_small_pool, rs_desc, dg, &
                           poisson_env)

      TYPE(ewald_pw_type), INTENT(INOUT)                 :: ewald_pw
      TYPE(pw_pool_type), OPTIONAL, POINTER              :: pw_big_pool, pw_small_pool
      TYPE(realspace_grid_desc_type), OPTIONAL, POINTER  :: rs_desc
      TYPE(dg_type), OPTIONAL, POINTER                   :: dg
      TYPE(pw_poisson_type), OPTIONAL, POINTER           :: poisson_env

      IF (PRESENT(pw_big_pool)) THEN
         CALL pw_big_pool%retain()
         CALL pw_pool_release(ewald_pw%pw_big_pool)
         ewald_pw%pw_big_pool => pw_big_pool
      END IF
      IF (PRESENT(pw_small_pool)) THEN
         CALL pw_small_pool%retain()
         CALL pw_pool_release(ewald_pw%pw_small_pool)
         ewald_pw%pw_small_pool => pw_small_pool
      END IF
      IF (PRESENT(rs_desc)) THEN
         CALL rs_grid_retain_descriptor(rs_desc)
         CALL rs_grid_release_descriptor(ewald_pw%rs_desc)
         ewald_pw%rs_desc => rs_desc
      END IF
      IF (PRESENT(dg)) THEN
         CALL dg_release(ewald_pw%dg)
         ewald_pw%dg => dg
      END IF
      IF (PRESENT(poisson_env)) THEN
         IF (ASSOCIATED(ewald_pw%poisson_env)) THEN
         IF (.NOT. ASSOCIATED(ewald_pw%poisson_env, poisson_env)) THEN
            CALL ewald_pw%poisson_env%release()
            DEALLOCATE (ewald_pw%poisson_env)
         END IF
         END IF
         ewald_pw%poisson_env => poisson_env
      END IF

   END SUBROUTINE ewald_pw_set

END MODULE ewald_pw_types
