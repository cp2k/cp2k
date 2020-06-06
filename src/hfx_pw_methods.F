!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Test routines for HFX caclulations using PW
!>
!>
!> \par History
!>     refactoring 03-2011 [MI]
!>     Made GAPW compatible 12.2019 (A. Bussy)
!>     Refactored from hfx_admm_utils (JGH)
!> \author MI
! **************************************************************************************************
MODULE hfx_pw_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE dbcsr_api,                       ONLY: dbcsr_type
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_copy,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_create,&
                                              pw_p_type,&
                                              pw_release,&
                                              pw_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: pw_hfx

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_pw_methods'

CONTAINS

! **************************************************************************************************
!> \brief computes the Hartree-Fock energy brute force in a pw basis
!> \param qs_env ...
!> \param ehfx ...
!> \param hfx_section ...
!> \param poisson_env ...
!> \param auxbas_pw_pool ...
!> \param irep ...
!> \par History
!>      12.2007 created [Joost VandeVondele]
!> \note
!>      only computes the HFX energy, no derivatives as yet
! **************************************************************************************************
   SUBROUTINE pw_hfx(qs_env, ehfx, hfx_section, poisson_env, auxbas_pw_pool, irep)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: ehfx
      TYPE(section_vals_type), POINTER                   :: hfx_section
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      INTEGER                                            :: irep

      CHARACTER(*), PARAMETER :: routineN = 'pw_hfx', routineP = moduleN//':'//routineN

      INTEGER                                            :: blocksize, handle, ig, iloc, iorb, &
                                                            iorb_block, ispin, iw, jloc, jorb, &
                                                            jorb_block, norb, potential_type
      LOGICAL                                            :: do_pw_hfx, explicit
      REAL(KIND=dp)                                      :: exchange_energy, fraction, g2, g3d, gg, &
                                                            omega, pair_energy, rcut, scaling
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_b
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_p_type)                                    :: pot_g, rho_g, rho_r
      TYPE(pw_p_type), ALLOCATABLE, DIMENSION(:)         :: rho_i, rho_j
      TYPE(pw_type), POINTER                             :: greenfn
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: ip_section

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      CALL section_vals_val_get(hfx_section, "PW_HFX", l_val=do_pw_hfx, i_rep_section=irep)

      IF (do_pw_hfx) THEN
         CALL section_vals_val_get(hfx_section, "FRACTION", r_val=fraction)
         CALL section_vals_val_get(hfx_section, "PW_HFX_BLOCKSIZE", i_val=blocksize)

         CALL get_qs_env(qs_env, mos=mo_array, pw_env=pw_env, dft_control=dft_control, &
                         cell=cell, particle_set=particle_set, &
                         atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)

         ! limit the blocksize by the number of orbitals
         CALL get_mo_set(mo_set=mo_array(1)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, ncol_global=norb)
         blocksize = MAX(1, MIN(blocksize, norb))

         CALL pw_pool_create_pw(auxbas_pw_pool, rho_r%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_g%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, pot_g%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)

         ALLOCATE (rho_i(blocksize))
         ALLOCATE (rho_j(blocksize))

         NULLIFY (greenfn)
         CALL pw_pool_create_pw(auxbas_pw_pool, greenfn, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         ip_section => section_vals_get_subs_vals(hfx_section, "INTERACTION_POTENTIAL")
         CALL section_vals_get(ip_section, explicit=explicit)
         potential_type = do_potential_coulomb
         IF (explicit) THEN
            CALL section_vals_val_get(ip_section, "POTENTIAL_TYPE", i_val=potential_type)
         ENDIF
         IF (potential_type == do_potential_coulomb) THEN
            CALL pw_copy(poisson_env%green_fft%influence_fn, greenfn)
         ELSEIF (potential_type == do_potential_truncated) THEN
            CALL section_vals_val_get(ip_section, "CUTOFF_RADIUS", r_val=rcut)
            grid => poisson_env%green_fft%influence_fn%pw_grid
            DO ig = grid%first_gne0, grid%ngpts_cut_local
               g2 = grid%gsq(ig)
               gg = SQRT(g2)
               g3d = fourpi/g2
               greenfn%cc(ig) = g3d*(1.0_dp - COS(rcut*gg))
            END DO
            IF (grid%have_g0) &
               greenfn%cc(1) = 0.5_dp*fourpi*rcut*rcut
         ELSEIF (potential_type == do_potential_short) THEN
            CALL section_vals_val_get(ip_section, "OMEGA", r_val=omega)
            IF (omega > 0.0_dp) omega = 0.25_dp/(omega*omega)
            grid => poisson_env%green_fft%influence_fn%pw_grid
            DO ig = grid%first_gne0, grid%ngpts_cut_local
               g2 = grid%gsq(ig)
               gg = -omega*g2
               g3d = fourpi/g2
               greenfn%cc(ig) = g3d*(1.0_dp - EXP(gg))
            END DO
            IF (grid%have_g0) greenfn%cc(1) = 0.0_dp
         ELSE
            CPWARN("PW_SCF: Potential type not supported, calculation uses Coulomb potential.")
         END IF

         DO iorb_block = 1, blocksize
            NULLIFY (rho_i(iorb_block)%pw)
            CALL pw_create(rho_i(iorb_block)%pw, rho_r%pw%pw_grid, &
                           use_data=REALDATA3D, &
                           in_space=REALSPACE)
            NULLIFY (rho_j(iorb_block)%pw)
            CALL pw_create(rho_j(iorb_block)%pw, rho_r%pw%pw_grid, &
                           use_data=REALDATA3D, &
                           in_space=REALSPACE)
         ENDDO

         exchange_energy = 0.0_dp

         DO ispin = 1, SIZE(mo_array)
            CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, mo_coeff=mo_coeff, mo_coeff_b=mo_coeff_b)

            IF (mo_array(ispin)%mo_set%use_mo_coeff_b) THEN !fm->dbcsr
               CALL copy_dbcsr_to_fm(mo_coeff_b, mo_coeff) !fm->dbcsr
            ENDIF !fm->dbcsr

            CALL cp_fm_get_info(mo_coeff, ncol_global=norb)

            DO iorb_block = 1, norb, blocksize

               DO iorb = iorb_block, MIN(iorb_block + blocksize - 1, norb)

                  iloc = iorb - iorb_block + 1
                  CALL calculate_wavefunction(mo_coeff, iorb, rho_i(iloc), rho_g, &
                                              atomic_kind_set, qs_kind_set, cell, dft_control, particle_set, &
                                              pw_env)

               ENDDO

               DO jorb_block = iorb_block, norb, blocksize

                  DO jorb = jorb_block, MIN(jorb_block + blocksize - 1, norb)

                     jloc = jorb - jorb_block + 1
                     CALL calculate_wavefunction(mo_coeff, jorb, rho_j(jloc), rho_g, &
                                                 atomic_kind_set, qs_kind_set, cell, dft_control, particle_set, &
                                                 pw_env)

                  ENDDO

                  DO iorb = iorb_block, MIN(iorb_block + blocksize - 1, norb)
                     iloc = iorb - iorb_block + 1
                     DO jorb = jorb_block, MIN(jorb_block + blocksize - 1, norb)
                        jloc = jorb - jorb_block + 1
                        IF (jorb < iorb) CYCLE

                        ! compute the pair density
                        rho_r%pw%cr3d = rho_i(iloc)%pw%cr3d*rho_j(jloc)%pw%cr3d

                        ! go the g-space and compute hartree energy
                        CALL pw_transfer(rho_r%pw, rho_g%pw)
                        CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, pot_g%pw, &
                                              greenfn=greenfn)

                        ! sum up to the full energy
                        scaling = fraction
                        IF (SIZE(mo_array) == 1) scaling = scaling*2.0_dp
                        IF (iorb /= jorb) scaling = scaling*2.0_dp

                        exchange_energy = exchange_energy - scaling*pair_energy

                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO

         DO iorb_block = 1, blocksize
            CALL pw_release(rho_i(iorb_block)%pw)
            CALL pw_release(rho_j(iorb_block)%pw)
         ENDDO

         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, pot_g%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, greenfn)

         iw = cp_print_key_unit_nr(logger, hfx_section, "HF_INFO", &
                                   extension=".scfLog")

         IF (iw > 0) THEN
            WRITE (UNIT=iw, FMT="((T3,A,T61,F20.10))") &
               "HF_PW_HFX| PW exchange energy:", exchange_energy
            WRITE (UNIT=iw, FMT="((T3,A,T61,F20.10),/)") &
               "HF_PW_HFX| Gaussian exchange energy:", ehfx
         ENDIF

         CALL cp_print_key_finished_output(iw, logger, hfx_section, "HF_INFO")

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE pw_hfx

END MODULE hfx_pw_methods
