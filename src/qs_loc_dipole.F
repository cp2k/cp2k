!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
! **************************************************************************************************
MODULE qs_loc_dipole
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE cp_output_handling,              ONLY: cp_iter_string,&
                                              cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_result_methods,               ONLY: cp_results_erase,&
                                              get_results,&
                                              put_results
   USE cp_result_types,                 ONLY: cp_result_type
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: twopi
   USE moments_utils,                   ONLY: get_reference_point
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: debye
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_loc_types,                    ONLY: qs_loc_env_new_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Global parameters
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_loc_dipole'
   PUBLIC :: loc_dipole

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Computes and prints the Dipole (using localized charges)
!> \param input ...
!> \param dft_control ...
!> \param qs_loc_env ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
! **************************************************************************************************
   SUBROUTINE loc_dipole(input, dft_control, qs_loc_env, logger, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'loc_dipole'

      CHARACTER(LEN=default_string_length)               :: description, descriptionThisDip, iter
      COMPLEX(KIND=dp)                                   :: zeta
      COMPLEX(KIND=dp), DIMENSION(3)                     :: ggamma, zphase
      INTEGER                                            :: handle, i, ikind, ispins, j, n_rep, &
                                                            reference, unit_nr
      LOGICAL                                            :: do_berry, first_time, floating, ghost
      REAL(KIND=dp)                                      :: charge_tot, theta, zeff, zwfc
      REAL(KIND=dp), DIMENSION(3)                        :: ci, dipole, dipole_old, gvec, rcc, ria
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      print_key => section_vals_get_subs_vals(input, "DFT%LOCALIZE%PRINT%TOTAL_DIPOLE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key, first_time=first_time) &
                , cp_p_file)) THEN
         NULLIFY (cell, particle_set, qs_kind_set, ref_point, results)
         CALL get_qs_env(qs_env=qs_env, &
                         cell=cell, &
                         particle_set=particle_set, &
                         qs_kind_set=qs_kind_set, &
                         results=results)

         reference = section_get_ival(print_key, keyword_name="REFERENCE")
         CALL section_vals_val_get(print_key, "REF_POINT", r_vals=ref_point)
         CALL section_vals_val_get(print_key, "PERIODIC", l_val=do_berry)
         description = '[DIPOLE]'
         descriptionThisDip = '[TOTAL_DIPOLE]'
         CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)

         dipole = 0.0_dp
         IF (do_berry) THEN
            rcc = pbc(rcc, cell)
            charge_tot = REAL(dft_control%charge, KIND=dp)
            ria = twopi*MATMUL(cell%h_inv, rcc)
            zphase = CMPLX(COS(ria), SIN(ria), KIND=dp)**charge_tot
            ggamma = CMPLX(1.0_dp, 0.0_dp, KIND=dp)

            ! Nuclear charges
            DO i = 1, SIZE(particle_set)
               CALL get_atomic_kind(particle_set(i)%atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), ghost=ghost, floating=floating)
               IF (.NOT. ghost .AND. .NOT. floating) THEN
                  CALL get_qs_kind(qs_kind_set(ikind), core_charge=zeff)
                  ria = pbc(particle_set(i)%r, cell)
                  DO j = 1, 3
                     gvec = twopi*cell%h_inv(j, :)
                     theta = SUM(ria(:)*gvec(:))
                     zeta = CMPLX(COS(theta), SIN(theta), KIND=dp)**(zeff)
                     ggamma(j) = ggamma(j)*zeta
                  END DO
               END IF
            END DO

            ! Charges of the wfc involved
            ! Warning, this assumes the same occupation for all states
            zwfc = 3.0_dp - REAL(dft_control%nspins, dp)

            DO ispins = 1, dft_control%nspins
               DO i = 1, SIZE(qs_loc_env%localized_wfn_control%centers_set(ispins)%array, 2)
                  ria = pbc(qs_loc_env%localized_wfn_control%centers_set(ispins)%array(1:3, i), cell)
                  DO j = 1, 3
                     gvec = twopi*cell%h_inv(j, :)
                     theta = SUM(ria(:)*gvec(:))
                     zeta = CMPLX(COS(theta), SIN(theta), KIND=dp)**(-zwfc)
                     ggamma(j) = ggamma(j)*zeta
                  END DO
               ENDDO
            ENDDO
            ggamma = ggamma*zphase
            ci = AIMAG(LOG(ggamma))/twopi
            dipole = MATMUL(cell%hmat, ci)
         ELSE
            ! Charges of the atoms involved
            DO i = 1, SIZE(particle_set)
               CALL get_atomic_kind(particle_set(i)%atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), ghost=ghost, floating=floating)
               IF (.NOT. ghost .AND. .NOT. floating) THEN
                  CALL get_qs_kind(qs_kind_set(ikind), core_charge=zeff)
                  ria = pbc(particle_set(i)%r, cell)
                  dipole = dipole + zeff*(ria - rcc)
               END IF
            END DO

            ! Charges of the wfc involved
            ! Warning, this assumes the same occupation for all states
            zwfc = 3.0_dp - REAL(dft_control%nspins, dp)

            DO ispins = 1, dft_control%nspins
               DO i = 1, SIZE(qs_loc_env%localized_wfn_control%centers_set(ispins)%array, 2)
                  ria = pbc(qs_loc_env%localized_wfn_control%centers_set(ispins)%array(1:3, i), cell)
                  dipole = dipole - zwfc*(ria - rcc)
               ENDDO
            ENDDO
         END IF

         ! Print and possibly store results
         unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".Dipole", &
                                        middle_name="TOTAL_DIPOLE")
         IF (unit_nr > 0) THEN
            IF (first_time) THEN
               WRITE (unit=unit_nr, fmt="(A,T31,A,T88,A,T136,A)") &
                  "# iter_level", "dipole(x,y,z)[atomic units]", &
                  "dipole(x,y,z)[debye]", &
                  "delta_dipole(x,y,z)[atomic units]"
            END IF
            iter = cp_iter_string(logger%iter_info)
            CALL get_results(results, descriptionThisDip, n_rep=n_rep)
            IF (n_rep == 0) THEN
               dipole_old = 0._dp
            ELSE
               CALL get_results(results, descriptionThisDip, dipole_old, nval=n_rep)
            END IF
            IF (do_berry) THEN
               WRITE (unit=unit_nr, fmt="(a,9(es18.8))") &
                  iter(1:15), dipole, dipole*debye, pbc(dipole - dipole_old, cell)
            ELSE
               WRITE (unit=unit_nr, fmt="(a,9(es18.8))") &
                  iter(1:15), dipole, dipole*debye, (dipole - dipole_old)
            END IF
         END IF
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
         CALL cp_results_erase(results, description)
         CALL put_results(results, description, dipole)
         CALL cp_results_erase(results, descriptionThisDip)
         CALL put_results(results, descriptionThisDip, dipole)
      END IF

      CALL timestop(handle)

   END SUBROUTINE loc_dipole

END MODULE qs_loc_dipole
