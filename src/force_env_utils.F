!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Util force_env module
!> \author Teodoro Laino [tlaino] - 02.2011
! **************************************************************************************************
MODULE force_env_utils

   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE cell_types,                      ONLY: cell_type
   USE constraint,                      ONLY: rattle_control,&
                                              shake_control
   USE constraint_util,                 ONLY: getold
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE molecule_kind_list_types,        ONLY: molecule_kind_list_type
   USE molecule_list_types,             ONLY: molecule_list_type
   USE molecule_types,                  ONLY: global_constraint_type
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: update_particle_set
   USE physcon,                         ONLY: angstrom
   USE string_utilities,                ONLY: lowercase,&
                                              uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_env_utils'

   PUBLIC :: force_env_shake, &
             force_env_rattle, &
             rescale_forces, &
             write_atener, &
             write_forces

CONTAINS

! **************************************************************************************************
!> \brief perform shake (enforcing of constraints)
!> \param force_env the force env to shake
!> \param dt the dt for shake (if you are not interested in the velocities
!>        it can be any positive number)
!> \param shake_tol the tolerance for shake
!> \param log_unit if >0 then some information on the shake is printed,
!>        defaults to -1
!> \param lagrange_mult ...
!> \param dump_lm ...
!> \param pos ...
!> \param vel ...
!> \param compold ...
!> \param reset ...
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE force_env_shake(force_env, dt, shake_tol, log_unit, lagrange_mult, dump_lm, &
                              pos, vel, compold, reset)

      TYPE(force_env_type), POINTER                      :: force_env
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: dt
      REAL(kind=dp), INTENT(IN)                          :: shake_tol
      INTEGER, INTENT(in), OPTIONAL                      :: log_unit, lagrange_mult
      LOGICAL, INTENT(IN), OPTIONAL                      :: dump_lm
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL, TARGET                                :: pos, vel
      LOGICAL, INTENT(IN), OPTIONAL                      :: compold, reset

      CHARACTER(len=*), PARAMETER                        :: routineN = 'force_env_shake'

      INTEGER :: handle, i, iparticle, iparticle_kind, iparticle_local, j, my_lagrange_mult, &
         my_log_unit, nparticle_kind, nparticle_local
      LOGICAL                                            :: has_pos, has_vel, my_dump_lm
      REAL(KIND=dp)                                      :: mydt
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: my_pos, my_vel
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_molecules, local_particles
      TYPE(global_constraint_type), POINTER              :: gci
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(particle_list_type), POINTER                  :: particles

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(force_env))
      CPASSERT(force_env%ref_count > 0)
      my_log_unit = -1
      IF (PRESENT(log_unit)) my_log_unit = log_unit
      my_lagrange_mult = -1
      IF (PRESENT(lagrange_mult)) my_lagrange_mult = lagrange_mult
      my_dump_lm = .FALSE.
      IF (PRESENT(dump_lm)) my_dump_lm = dump_lm
      NULLIFY (subsys, cell, molecules, molecule_kinds, local_molecules, particles, &
               my_pos, my_vel, gci)
      IF (PRESENT(pos)) my_pos => pos
      IF (PRESENT(vel)) my_vel => vel
      mydt = 0.1_dp
      IF (PRESENT(dt)) mydt = dt
      CALL force_env_get(force_env, subsys=subsys, cell=cell)
      CALL cp_subsys_get(subsys, &
                         atomic_kinds=atomic_kinds, &
                         local_molecules=local_molecules, &
                         local_particles=local_particles, &
                         molecules=molecules, &
                         molecule_kinds=molecule_kinds, &
                         particles=particles, &
                         gci=gci)
      nparticle_kind = atomic_kinds%n_els
      IF (PRESENT(compold)) THEN
         IF (compold) THEN
            CALL getold(gci, local_molecules, molecules%els, molecule_kinds%els, &
                        particles%els, cell)
         END IF
      END IF
      has_pos = .FALSE.
      IF (.NOT. ASSOCIATED(my_pos)) THEN
         has_pos = .TRUE.
         ALLOCATE (my_pos(3, particles%n_els))
         my_pos = 0.0_dp
         DO iparticle_kind = 1, nparticle_kind
            nparticle_local = local_particles%n_el(iparticle_kind)
            DO iparticle_local = 1, nparticle_local
               iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
               my_pos(:, iparticle) = particles%els(iparticle)%r(:)
            END DO
         END DO
      END IF
      has_vel = .FALSE.
      IF (.NOT. ASSOCIATED(my_vel)) THEN
         has_vel = .TRUE.
         ALLOCATE (my_vel(3, particles%n_els))
         my_vel = 0.0_dp
         DO iparticle_kind = 1, nparticle_kind
            nparticle_local = local_particles%n_el(iparticle_kind)
            DO iparticle_local = 1, nparticle_local
               iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
               my_vel(:, iparticle) = particles%els(iparticle)%v(:)
            END DO
         END DO
      END IF

      CALL shake_control(gci=gci, local_molecules=local_molecules, &
                         molecule_set=molecules%els, molecule_kind_set=molecule_kinds%els, &
                         particle_set=particles%els, pos=my_pos, vel=my_vel, dt=mydt, &
                         shake_tol=shake_tol, log_unit=my_log_unit, lagrange_mult=my_lagrange_mult, &
                         dump_lm=my_dump_lm, cell=cell, group=force_env%para_env, &
                         local_particles=local_particles)

      ! Possibly reset the lagrange multipliers
      IF (PRESENT(reset)) THEN
         IF (reset) THEN
            ! Reset Intramolecular constraints
            DO i = 1, SIZE(molecules%els)
               IF (ASSOCIATED(molecules%els(i)%lci%lcolv)) THEN
                  DO j = 1, SIZE(molecules%els(i)%lci%lcolv)
                     ! Reset langrange multiplier
                     molecules%els(i)%lci%lcolv(j)%lambda = 0.0_dp
                  END DO
               END IF
               IF (ASSOCIATED(molecules%els(i)%lci%lg3x3)) THEN
                  DO j = 1, SIZE(molecules%els(i)%lci%lg3x3)
                     ! Reset langrange multiplier
                     molecules%els(i)%lci%lg3x3(j)%lambda = 0.0_dp
                  END DO
               END IF
               IF (ASSOCIATED(molecules%els(i)%lci%lg4x6)) THEN
                  DO j = 1, SIZE(molecules%els(i)%lci%lg4x6)
                     ! Reset langrange multiplier
                     molecules%els(i)%lci%lg4x6(j)%lambda = 0.0_dp
                  END DO
               END IF
            END DO
            ! Reset Intermolecular constraints
            IF (ASSOCIATED(gci)) THEN
               IF (ASSOCIATED(gci%lcolv)) THEN
                  DO j = 1, SIZE(gci%lcolv)
                     ! Reset langrange multiplier
                     gci%lcolv(j)%lambda = 0.0_dp
                  END DO
               END IF
               IF (ASSOCIATED(gci%lg3x3)) THEN
                  DO j = 1, SIZE(gci%lg3x3)
                     ! Reset langrange multiplier
                     gci%lg3x3(j)%lambda = 0.0_dp
                  END DO
               END IF
               IF (ASSOCIATED(gci%lg4x6)) THEN
                  DO j = 1, SIZE(gci%lg4x6)
                     ! Reset langrange multiplier
                     gci%lg4x6(j)%lambda = 0.0_dp
                  END DO
               END IF
            END IF
         END IF
      END IF

      IF (has_pos) THEN
         CALL update_particle_set(particles%els, force_env%para_env, pos=my_pos)
         DEALLOCATE (my_pos)
      END IF
      IF (has_vel) THEN
         CALL update_particle_set(particles%els, force_env%para_env, vel=my_vel)
         DEALLOCATE (my_vel)
      END IF
      CALL timestop(handle)
   END SUBROUTINE force_env_shake

! **************************************************************************************************
!> \brief perform rattle (enforcing of constraints on velocities)
!>      This routine can be easily adapted to performe rattle on whatever
!>      other vector different from forces..
!> \param force_env the force env to shake
!> \param dt the dt for shake (if you are not interested in the velocities
!>        it can be any positive number)
!> \param shake_tol the tolerance for shake
!> \param log_unit if >0 then some information on the shake is printed,
!>        defaults to -1
!> \param lagrange_mult ...
!> \param dump_lm ...
!> \param vel ...
!> \param reset ...
!> \author tlaino
! **************************************************************************************************
   SUBROUTINE force_env_rattle(force_env, dt, shake_tol, log_unit, lagrange_mult, dump_lm, &
                               vel, reset)

      TYPE(force_env_type), POINTER                      :: force_env
      REAL(kind=dp), INTENT(in), OPTIONAL                :: dt
      REAL(kind=dp), INTENT(in)                          :: shake_tol
      INTEGER, INTENT(in), OPTIONAL                      :: log_unit, lagrange_mult
      LOGICAL, INTENT(IN), OPTIONAL                      :: dump_lm
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL, TARGET                                :: vel
      LOGICAL, INTENT(IN), OPTIONAL                      :: reset

      CHARACTER(len=*), PARAMETER                        :: routineN = 'force_env_rattle'

      INTEGER :: handle, i, iparticle, iparticle_kind, iparticle_local, j, my_lagrange_mult, &
         my_log_unit, nparticle_kind, nparticle_local
      LOGICAL                                            :: has_vel, my_dump_lm
      REAL(KIND=dp)                                      :: mydt
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: my_vel
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_molecules, local_particles
      TYPE(global_constraint_type), POINTER              :: gci
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(particle_list_type), POINTER                  :: particles

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(force_env))
      CPASSERT(force_env%ref_count > 0)
      my_log_unit = -1
      IF (PRESENT(log_unit)) my_log_unit = log_unit
      my_lagrange_mult = -1
      IF (PRESENT(lagrange_mult)) my_lagrange_mult = lagrange_mult
      my_dump_lm = .FALSE.
      IF (PRESENT(dump_lm)) my_dump_lm = dump_lm
      NULLIFY (subsys, cell, molecules, molecule_kinds, local_molecules, particles, &
               my_vel)
      IF (PRESENT(vel)) my_vel => vel
      mydt = 0.1_dp
      IF (PRESENT(dt)) mydt = dt
      CALL force_env_get(force_env, subsys=subsys, cell=cell)
      CALL cp_subsys_get(subsys, &
                         atomic_kinds=atomic_kinds, &
                         local_molecules=local_molecules, &
                         local_particles=local_particles, &
                         molecules=molecules, &
                         molecule_kinds=molecule_kinds, &
                         particles=particles, &
                         gci=gci)
      nparticle_kind = atomic_kinds%n_els
      has_vel = .FALSE.
      IF (.NOT. ASSOCIATED(my_vel)) THEN
         has_vel = .TRUE.
         ALLOCATE (my_vel(3, particles%n_els))
         my_vel = 0.0_dp
         DO iparticle_kind = 1, nparticle_kind
            nparticle_local = local_particles%n_el(iparticle_kind)
            DO iparticle_local = 1, nparticle_local
               iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
               my_vel(:, iparticle) = particles%els(iparticle)%v(:)
            END DO
         END DO
      END IF

      CALL rattle_control(gci=gci, local_molecules=local_molecules, &
                          molecule_set=molecules%els, molecule_kind_set=molecule_kinds%els, &
                          particle_set=particles%els, vel=my_vel, dt=mydt, &
                          rattle_tol=shake_tol, log_unit=my_log_unit, lagrange_mult=my_lagrange_mult, &
                          dump_lm=my_dump_lm, cell=cell, group=force_env%para_env, &
                          local_particles=local_particles)

      ! Possibly reset the lagrange multipliers
      IF (PRESENT(reset)) THEN
         IF (reset) THEN
            ! Reset Intramolecular constraints
            DO i = 1, SIZE(molecules%els)
               IF (ASSOCIATED(molecules%els(i)%lci%lcolv)) THEN
                  DO j = 1, SIZE(molecules%els(i)%lci%lcolv)
                     ! Reset langrange multiplier
                     molecules%els(i)%lci%lcolv(j)%lambda = 0.0_dp
                  END DO
               END IF
               IF (ASSOCIATED(molecules%els(i)%lci%lg3x3)) THEN
                  DO j = 1, SIZE(molecules%els(i)%lci%lg3x3)
                     ! Reset langrange multiplier
                     molecules%els(i)%lci%lg3x3(j)%lambda = 0.0_dp
                  END DO
               END IF
               IF (ASSOCIATED(molecules%els(i)%lci%lg4x6)) THEN
                  DO j = 1, SIZE(molecules%els(i)%lci%lg4x6)
                     ! Reset langrange multiplier
                     molecules%els(i)%lci%lg4x6(j)%lambda = 0.0_dp
                  END DO
               END IF
            END DO
            ! Reset Intermolecular constraints
            IF (ASSOCIATED(gci)) THEN
               IF (ASSOCIATED(gci%lcolv)) THEN
                  DO j = 1, SIZE(gci%lcolv)
                     ! Reset langrange multiplier
                     gci%lcolv(j)%lambda = 0.0_dp
                  END DO
               END IF
               IF (ASSOCIATED(gci%lg3x3)) THEN
                  DO j = 1, SIZE(gci%lg3x3)
                     ! Reset langrange multiplier
                     gci%lg3x3(j)%lambda = 0.0_dp
                  END DO
               END IF
               IF (ASSOCIATED(gci%lg4x6)) THEN
                  DO j = 1, SIZE(gci%lg4x6)
                     ! Reset langrange multiplier
                     gci%lg4x6(j)%lambda = 0.0_dp
                  END DO
               END IF
            END IF
         END IF
      END IF

      IF (has_vel) THEN
         CALL update_particle_set(particles%els, force_env%para_env, vel=my_vel)
      END IF
      DEALLOCATE (my_vel)
      CALL timestop(handle)
   END SUBROUTINE force_env_rattle

! **************************************************************************************************
!> \brief Rescale forces if requested
!> \param force_env the force env to shake
!> \author tlaino
! **************************************************************************************************
   SUBROUTINE rescale_forces(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'rescale_forces'

      INTEGER                                            :: handle, iparticle
      LOGICAL                                            :: explicit
      REAL(KIND=dp)                                      :: force(3), max_value, mod_force
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(section_vals_type), POINTER                   :: rescale_force_section

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(force_env))
      CPASSERT(force_env%ref_count > 0)
      rescale_force_section => section_vals_get_subs_vals(force_env%force_env_section, "RESCALE_FORCES")
      CALL section_vals_get(rescale_force_section, explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(rescale_force_section, "MAX_FORCE", r_val=max_value)
         CALL force_env_get(force_env, subsys=subsys)
         CALL cp_subsys_get(subsys, particles=particles)
         DO iparticle = 1, SIZE(particles%els)
            force = particles%els(iparticle)%f(:)
            mod_force = SQRT(DOT_PRODUCT(force, force))
            IF ((mod_force > max_value) .AND. (mod_force /= 0.0_dp)) THEN
               force = force/mod_force*max_value
               particles%els(iparticle)%f(:) = force
            END IF
         END DO
      END IF
      CALL timestop(handle)
   END SUBROUTINE rescale_forces

! **************************************************************************************************
!> \brief Write forces either to the screen or to a file.
!> \param particles ...
!> \param iw ...
!> \param label ...
!> \param ndigits ...
!> \param unit_string ...
!> \param total_force ...
!> \param grand_total_force ...
!> \param zero_force_core_shell_atom ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE write_forces(particles, iw, label, ndigits, unit_string, total_force, &
                           grand_total_force, zero_force_core_shell_atom)

      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, INTENT(IN)                                :: iw
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      INTEGER, INTENT(IN)                                :: ndigits
      CHARACTER(LEN=default_string_length), INTENT(IN)   :: unit_string
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: total_force
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT), &
         OPTIONAL                                        :: grand_total_force
      LOGICAL, INTENT(IN), OPTIONAL                      :: zero_force_core_shell_atom

      IF (iw == cp_logger_get_default_io_unit()) THEN
         CALL write_forces_to_screen(particles, iw, label, ndigits, unit_string, total_force, &
                                     grand_total_force, zero_force_core_shell_atom)
      ELSE
         CALL write_forces_to_file(particles, iw, label, ndigits, total_force, &
                                   grand_total_force, zero_force_core_shell_atom)
      END IF
   END SUBROUTINE write_forces

! **************************************************************************************************
!> \brief Write forces to the screen.
!>
!> \param particles ...
!> \param iw ...
!> \param label ...
!> \param ndigits ...
!> \param unit_string ...
!> \param total_force ...
!> \param grand_total_force ...
!> \param zero_force_core_shell_atom ...
!> \author MK (06.09.2010)
! **************************************************************************************************
   SUBROUTINE write_forces_to_screen(particles, iw, label, ndigits, unit_string, total_force, &
                                     grand_total_force, zero_force_core_shell_atom)

      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, INTENT(IN)                                :: iw
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      INTEGER, INTENT(IN)                                :: ndigits
      CHARACTER(LEN=default_string_length), INTENT(IN)   :: unit_string
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: total_force
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT), &
         OPTIONAL                                        :: grand_total_force
      LOGICAL, INTENT(IN), OPTIONAL                      :: zero_force_core_shell_atom

      CHARACTER(LEN=18)                                  :: fmtstr4
      CHARACTER(LEN=29)                                  :: fmtstr3
      CHARACTER(LEN=38)                                  :: fmtstr2
      CHARACTER(LEN=49)                                  :: fmtstr1
      CHARACTER(LEN=7)                                   :: tag
      CHARACTER(LEN=LEN_TRIM(label))                     :: lc_label
      INTEGER                                            :: i, iparticle, n
      LOGICAL                                            :: zero_force
      REAL(KIND=dp)                                      :: fconv
      REAL(KIND=dp), DIMENSION(3)                        :: f

      IF (iw > 0) THEN
         CPASSERT(ASSOCIATED(particles))
         tag = "FORCES|"
         lc_label = TRIM(label)
         CALL lowercase(lc_label)
         n = MIN(MAX(1, ndigits), 20)
         fmtstr1 = "(/,T2,A,1X,A,/,T2,A,3X,A,T20,A3,2(  X,A3),  X,A3)"
         WRITE (UNIT=fmtstr1(35:36), FMT="(I2)") n + 5
         WRITE (UNIT=fmtstr1(43:44), FMT="(I2)") n + 6
         fmtstr2 = "(T2,A,I7,T16,3(1X,ES  .  ),2X,ES  .  )"
         WRITE (UNIT=fmtstr2(21:22), FMT="(I2)") n + 7
         WRITE (UNIT=fmtstr2(24:25), FMT="(I2)") n
         WRITE (UNIT=fmtstr2(33:34), FMT="(I2)") n + 7
         WRITE (UNIT=fmtstr2(36:37), FMT="(I2)") n
         fmtstr3 = "(T2,A,T16,3(1X,ES  .  ))"
         WRITE (UNIT=fmtstr3(18:19), FMT="(I2)") n + 7
         WRITE (UNIT=fmtstr3(21:22), FMT="(I2)") n
         fmtstr4 = "(T2,A,T  ,ES  .  )"
         WRITE (UNIT=fmtstr4(8:9), FMT="(I2)") 3*(n + 8) + 18
         WRITE (UNIT=fmtstr4(13:14), FMT="(I2)") n + 7
         WRITE (UNIT=fmtstr4(16:17), FMT="(I2)") n
         IF (PRESENT(zero_force_core_shell_atom)) THEN
            zero_force = zero_force_core_shell_atom
         ELSE
            zero_force = .FALSE.
         END IF
         fconv = cp_unit_from_cp2k(1.0_dp, TRIM(unit_string))
         WRITE (UNIT=iw, FMT=fmtstr1) &
            tag, label//" forces ["//TRIM(ADJUSTL(unit_string))//"]", &
            tag, "Atom", " x ", " y ", " z ", "|f|"
         total_force(1:3) = 0.0_dp
         DO iparticle = 1, particles%n_els
            IF (particles%els(iparticle)%atom_index /= 0) THEN
               i = particles%els(iparticle)%atom_index
            ELSE
               i = iparticle
            END IF
            IF (zero_force .AND. (particles%els(iparticle)%shell_index /= 0)) THEN
               f(1:3) = 0.0_dp
            ELSE
               f(1:3) = particles%els(iparticle)%f(1:3)*fconv
            END IF
            WRITE (UNIT=iw, FMT=fmtstr2) &
               tag, i, f(1:3), SQRT(SUM(f(1:3)**2))
            total_force(1:3) = total_force(1:3) + f(1:3)
         END DO
         WRITE (UNIT=iw, FMT=fmtstr3) &
            tag//" Sum", total_force(1:3)
         WRITE (UNIT=iw, FMT=fmtstr4) &
            tag//" Total "//TRIM(ADJUSTL(lc_label))//" force", &
            SQRT(SUM(total_force(1:3)**2))
      END IF

      IF (PRESENT(grand_total_force)) THEN
         grand_total_force(1:3) = grand_total_force(1:3) + total_force(1:3)
         WRITE (UNIT=iw, FMT="(A)") ""
         WRITE (UNIT=iw, FMT=fmtstr4) &
            tag//" Grand total force ["//TRIM(ADJUSTL(unit_string))//"]", &
            SQRT(SUM(grand_total_force(1:3)**2))
      END IF

   END SUBROUTINE write_forces_to_screen

! **************************************************************************************************
!> \brief Write forces to a file using our stable legacy format. Please don't change the format.
!>
!> \param particles ...
!> \param iw ...
!> \param label ...
!> \param ndigits ...
!> \param total_force ...
!> \param grand_total_force ...
!> \param zero_force_core_shell_atom ...
!> \author MK (06.09.2010)
! **************************************************************************************************
   SUBROUTINE write_forces_to_file(particles, iw, label, ndigits, total_force, &
                                   grand_total_force, zero_force_core_shell_atom)

      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, INTENT(IN)                                :: iw
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      INTEGER, INTENT(IN)                                :: ndigits
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: total_force
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT), &
         OPTIONAL                                        :: grand_total_force
      LOGICAL, INTENT(IN), OPTIONAL                      :: zero_force_core_shell_atom

      CHARACTER(LEN=23)                                  :: fmtstr3
      CHARACTER(LEN=36)                                  :: fmtstr2
      CHARACTER(LEN=46)                                  :: fmtstr1
      CHARACTER(LEN=LEN_TRIM(label))                     :: uc_label
      INTEGER                                            :: i, ikind, iparticle, n
      LOGICAL                                            :: zero_force
      REAL(KIND=dp), DIMENSION(3)                        :: f

      IF (iw > 0) THEN
         CPASSERT(ASSOCIATED(particles))
         n = MIN(MAX(1, ndigits), 20)
         fmtstr1 = "(/,T2,A,/,/,T2,A,T11,A,T18,A,T35,A1,2(  X,A1))"
         WRITE (UNIT=fmtstr1(39:40), FMT="(I2)") n + 6
         fmtstr2 = "(T2,I6,1X,I6,T21,A,T28,3(1X,F  .  ))"
         WRITE (UNIT=fmtstr2(33:34), FMT="(I2)") n
         WRITE (UNIT=fmtstr2(30:31), FMT="(I2)") n + 6
         fmtstr3 = "(T2,A,T28,4(1X,F  .  ))"
         WRITE (UNIT=fmtstr3(20:21), FMT="(I2)") n
         WRITE (UNIT=fmtstr3(17:18), FMT="(I2)") n + 6
         IF (PRESENT(zero_force_core_shell_atom)) THEN
            zero_force = zero_force_core_shell_atom
         ELSE
            zero_force = .FALSE.
         END IF
         uc_label = TRIM(label)
         CALL uppercase(uc_label)
         WRITE (UNIT=iw, FMT=fmtstr1) &
            uc_label//" FORCES in [a.u.]", "# Atom", "Kind", "Element", "X", "Y", "Z"
         total_force(1:3) = 0.0_dp
         DO iparticle = 1, particles%n_els
            ikind = particles%els(iparticle)%atomic_kind%kind_number
            IF (particles%els(iparticle)%atom_index /= 0) THEN
               i = particles%els(iparticle)%atom_index
            ELSE
               i = iparticle
            END IF
            IF (zero_force .AND. (particles%els(iparticle)%shell_index /= 0)) THEN
               f(1:3) = 0.0_dp
            ELSE
               f(1:3) = particles%els(iparticle)%f(1:3)
            END IF
            WRITE (UNIT=iw, FMT=fmtstr2) &
               i, ikind, particles%els(iparticle)%atomic_kind%element_symbol, f(1:3)
            total_force(1:3) = total_force(1:3) + f(1:3)
         END DO
         WRITE (UNIT=iw, FMT=fmtstr3) &
            "SUM OF "//uc_label//" FORCES", total_force(1:3), SQRT(SUM(total_force(:)**2))
      END IF

      IF (PRESENT(grand_total_force)) THEN
         grand_total_force(1:3) = grand_total_force(1:3) + total_force(1:3)
         WRITE (UNIT=iw, FMT="(A)") ""
         WRITE (UNIT=iw, FMT=fmtstr3) &
            "GRAND TOTAL FORCE", grand_total_force(1:3), SQRT(SUM(grand_total_force(:)**2))
      END IF

   END SUBROUTINE write_forces_to_file

! **************************************************************************************************
!> \brief Write the atomic coordinates, types, and energies
!> \param iounit ...
!> \param particles ...
!> \param atener ...
!> \param label ...
!> \date    05.06.2023
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE write_atener(iounit, particles, atener, label)

      INTEGER, INTENT(IN)                                :: iounit
      TYPE(particle_list_type), POINTER                  :: particles
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: atener
      CHARACTER(LEN=*), INTENT(IN)                       :: label

      INTEGER                                            :: i, natom

      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") TRIM(label)
         WRITE (UNIT=iounit, FMT="(T4,A,T30,A,T42,A,T54,A,T69,A)") &
            "Atom  Element", "X", "Y", "Z", "Energy[a.u.]"
         natom = particles%n_els
         DO i = 1, natom
            WRITE (iounit, "(I6,T12,A2,T24,3F12.6,F21.10)") i, &
               TRIM(ADJUSTL(particles%els(i)%atomic_kind%element_symbol)), &
               particles%els(i)%r(1:3)*angstrom, atener(i)
         END DO
         WRITE (iounit, "(A)") ""
      END IF

   END SUBROUTINE write_atener

END MODULE force_env_utils
