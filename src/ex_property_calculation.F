!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for property calculations of excited states
!> \par History
!>       02.2020 Adapted from ec_properties
!> \author JGH
! **************************************************************************************************
MODULE ex_property_calculation
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_result_methods,               ONLY: cp_results_erase,&
                                              put_results
   USE cp_result_types,                 ONLY: cp_result_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_dot,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_lval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_sum
   USE moments_utils,                   ONLY: get_reference_point
   USE mulliken,                        ONLY: mulliken_charges
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: debye
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_p_env_types,                  ONLY: qs_p_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ex_property_calculation'

   PUBLIC :: ex_properties

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_pe ...
!> \param p_env ...
! **************************************************************************************************
   SUBROUTINE ex_properties(qs_env, matrix_pe, p_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pe
      TYPE(qs_p_env_type), POINTER                       :: p_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ex_properties', routineP = moduleN//':'//routineN

      CHARACTER(LEN=8), DIMENSION(3)                     :: rlab
      CHARACTER(LEN=default_path_length)                 :: filename
      CHARACTER(LEN=default_string_length)               :: description
      INTEGER                                            :: akind, handle, i, ia, iatom, idir, &
                                                            ikind, iounit, ispin, maxmom, natom, &
                                                            nspins, reference, unit_nr
      LOGICAL                                            :: magnetic, periodic, tb
      REAL(KIND=dp)                                      :: charge, dd, q, tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mcharge
      REAL(KIND=dp), DIMENSION(3)                        :: cdip, pdip, pedip, rcc, rdip, ria, tdip
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p, matrix_s, moments
      TYPE(dbcsr_type), POINTER                          :: matrix_pall
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      rlab(1) = "X"
      rlab(2) = "Y"
      rlab(3) = "Z"

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      tb = (dft_control%qs_control%dftb .OR. dft_control%qs_control%xtb)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      END IF

      print_key => section_vals_get_subs_vals(section_vals=qs_env%input, &
                                              subsection_name="DFT%PRINT%MOMENTS")

      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN

         maxmom = section_get_ival(section_vals=qs_env%input, &
                                   keyword_name="DFT%PRINT%MOMENTS%MAX_MOMENT")
         periodic = section_get_lval(section_vals=qs_env%input, &
                                     keyword_name="DFT%PRINT%MOMENTS%PERIODIC")
         reference = section_get_ival(section_vals=qs_env%input, &
                                      keyword_name="DFT%PRINT%MOMENTS%REFERENCE")
         magnetic = section_get_lval(section_vals=qs_env%input, &
                                     keyword_name="DFT%PRINT%MOMENTS%MAGNETIC")
         NULLIFY (ref_point)
         CALL section_vals_val_get(qs_env%input, "DFT%PRINT%MOMENTS%REF_POINT", r_vals=ref_point)
         unit_nr = cp_print_key_unit_nr(logger=logger, basis_section=qs_env%input, &
                                        print_key_path="DFT%PRINT%MOMENTS", extension=".dat", &
                                        middle_name="moments", log_filename=.FALSE.)

         IF (iounit > 0) THEN
            IF (unit_nr /= iounit .AND. unit_nr > 0) THEN
               INQUIRE (UNIT=unit_nr, NAME=filename)
               WRITE (UNIT=iounit, FMT="(/,T2,A,2(/,T3,A),/)") &
                  "MOMENTS", "The electric/magnetic moments are written to file:", &
                  TRIM(filename)
            ELSE
               WRITE (UNIT=iounit, FMT="(/,T2,A)") "ELECTRIC/MAGNETIC MOMENTS"
            END IF
         END IF

         IF (periodic) THEN
            CPABORT("Periodic moments not implemented with TDDFT")
         ELSE
            CPASSERT(maxmom < 2)
            CPASSERT(.NOT. magnetic)
            IF (maxmom == 1) THEN
               CALL get_qs_env(qs_env=qs_env, cell=cell, para_env=para_env)
               ! reference point
               CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)
               ! nuclear contribution
               cdip = 0.0_dp
               CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, &
                               qs_kind_set=qs_kind_set, local_particles=local_particles)
               DO ikind = 1, SIZE(local_particles%n_el)
                  DO ia = 1, local_particles%n_el(ikind)
                     iatom = local_particles%list(ikind)%array(ia)
                     ! fold atomic positions back into unit cell
                     ria = pbc(particle_set(iatom)%r - rcc, cell) + rcc
                     ria = ria - rcc
                     atomic_kind => particle_set(iatom)%atomic_kind
                     CALL get_atomic_kind(atomic_kind, kind_number=akind)
                     CALL get_qs_kind(qs_kind_set(akind), core_charge=charge)
                     cdip(1:3) = cdip(1:3) - charge*ria(1:3)
                  END DO
               END DO
               CALL mp_sum(cdip, para_env%group)
               !
               ! electronic contribution
               CALL get_qs_env(qs_env=qs_env, rho=rho, matrix_s=matrix_s)
               CALL qs_rho_get(rho, rho_ao=matrix_p)
               nspins = SIZE(matrix_p, 1)
               IF (tb) THEN
                  ALLOCATE (matrix_pall)
                  CALL dbcsr_create(matrix_pall, template=matrix_s(1)%matrix)
                  CALL dbcsr_copy(matrix_pall, matrix_s(1)%matrix, "Moments")
                  CALL dbcsr_set(matrix_pall, 0.0_dp)
                  DO ispin = 1, nspins
                     CALL dbcsr_add(matrix_pall, matrix_p(ispin)%matrix, 1.0_dp, 1.0_dp)
                     CALL dbcsr_add(matrix_pall, matrix_pe(ispin)%matrix, 1.0_dp, 1.0_dp)
                     CALL dbcsr_add(matrix_pall, p_env%p1(ispin)%matrix, 1.0_dp, 1.0_dp)
                  END DO
                  CALL get_qs_env(qs_env=qs_env, natom=natom)
                  ! Mulliken charges
                  ALLOCATE (mcharge(natom))
                  !
                  CALL mulliken_charges(matrix_pall, matrix_s(1)%matrix, para_env, mcharge)
                  !
                  rdip = 0.0_dp
                  pdip = 0.0_dp
                  pedip = 0.0_dp
                  DO i = 1, SIZE(particle_set)
                     ria = pbc(particle_set(i)%r - rcc, cell) + rcc
                     ria = ria - rcc
                     q = mcharge(i)
                     rdip = rdip + q*ria
                  END DO
                  CALL dbcsr_release(matrix_pall)
                  DEALLOCATE (matrix_pall)
                  DEALLOCATE (mcharge)
               ELSE
                  ! KS-DFT
                  NULLIFY (moments)
                  CALL dbcsr_allocate_matrix_set(moments, 4)
                  DO i = 1, 4
                     ALLOCATE (moments(i)%matrix)
                     CALL dbcsr_copy(moments(i)%matrix, matrix_s(1)%matrix, "Moments")
                     CALL dbcsr_set(moments(i)%matrix, 0.0_dp)
                  END DO
                  CALL build_local_moment_matrix(qs_env, moments, 1, ref_point=rcc)
                  !
                  rdip = 0.0_dp
                  pdip = 0.0_dp
                  pedip = 0.0_dp
                  DO ispin = 1, nspins
                     DO idir = 1, 3
                        CALL dbcsr_dot(matrix_pe(ispin)%matrix, moments(idir)%matrix, tmp)
                        pedip(idir) = pedip(idir) + tmp
                        CALL dbcsr_dot(matrix_p(ispin)%matrix, moments(idir)%matrix, tmp)
                        pdip(idir) = pdip(idir) + tmp
                        CALL dbcsr_dot(p_env%p1(ispin)%matrix, moments(idir)%matrix, tmp)
                        rdip(idir) = rdip(idir) + tmp
                     END DO
                  END DO
                  CALL dbcsr_deallocate_matrix_set(moments)
               END IF

               tdip = -(rdip + pedip + pdip + cdip)

               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, "(T3,A)") "Dipoles are based on the traditional operator."
                  dd = SQRT(SUM(tdip(1:3)**2))*debye
                  WRITE (unit_nr, "(T3,A)") "Dipole moment [Debye]"
                  WRITE (unit_nr, "(T5,3(A,A,F14.8,1X),T60,A,T67,F14.8)") &
                     (TRIM(rlab(i)), "=", tdip(i)*debye, i=1, 3), "Total=", dd
                  WRITE (unit_nr, FMT="(T2,A,T61,E20.12)") ' DIPOLE : CheckSum  =', SUM(ABS(tdip))
               END IF
            END IF
         END IF

         CALL get_qs_env(qs_env=qs_env, results=results)
         description = "[DIPOLE]"
         CALL cp_results_erase(results=results, description=description)
         CALL put_results(results=results, description=description, values=tdip(1:3))

         CALL cp_print_key_finished_output(unit_nr=unit_nr, logger=logger, &
                                           basis_section=qs_env%input, print_key_path="DFT%PRINT%MOMENTS")
      END IF

      CALL timestop(handle)

   END SUBROUTINE ex_properties

END MODULE ex_property_calculation
