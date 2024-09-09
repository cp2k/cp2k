!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Harris method environment setup and handling
!> \par History
!>       2024.07 created
!> \author JGH
! **************************************************************************************************
MODULE qs_harris_utils
   USE atom_kind_orbitals,              ONLY: calculate_atomic_density
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE ec_methods,                      ONLY: create_kernel
   USE input_constants,                 ONLY: hden_atomic,&
                                              hfun_harris,&
                                              horb_default
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integral_ab,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec,&
                                              collocate_function
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_harris_types,                 ONLY: harris_energy_type,&
                                              harris_print_energy,&
                                              harris_rhoin_type,&
                                              harris_type
   USE qs_integrate_potential,          ONLY: integrate_function,&
                                              integrate_v_core_rspace,&
                                              integrate_v_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_harris_utils'

   PUBLIC :: harris_env_create, harris_write_input, harris_density_update, calculate_harris_density, &
             harris_energy_correction, harris_set_potentials

CONTAINS

! **************************************************************************************************
!> \brief Allocates and intitializes harris_env
!> \param qs_env The QS environment
!> \param harris_env The Harris method environment (the object to create)
!> \param harris_section The Harris method input section
!> \par History
!>       2024.07 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE harris_env_create(qs_env, harris_env, harris_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(harris_type), POINTER                         :: harris_env
      TYPE(section_vals_type), OPTIONAL, POINTER         :: harris_section

      CPASSERT(.NOT. ASSOCIATED(harris_env))
      ALLOCATE (harris_env)
      CALL init_harris_env(qs_env, harris_env, harris_section)

   END SUBROUTINE harris_env_create

! **************************************************************************************************
!> \brief Initializes Harris method environment
!> \param qs_env The QS environment
!> \param harris_env The Harris method environment
!> \param harris_section The Harris method input section
!> \par History
!>       2024.07 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE init_harris_env(qs_env, harris_env, harris_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(harris_type), POINTER                         :: harris_env
      TYPE(section_vals_type), OPTIONAL, POINTER         :: harris_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'init_harris_env'

      INTEGER                                            :: handle, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      IF (qs_env%harris_method) THEN

         CPASSERT(PRESENT(harris_section))
         ! get a useful output_unit
         logger => cp_get_default_logger()
         IF (logger%para_env%is_source()) THEN
            unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
         ELSE
            unit_nr = -1
         END IF

         CALL section_vals_val_get(harris_section, "ENERGY_FUNCTIONAL", &
                                   i_val=harris_env%energy_functional)
         CALL section_vals_val_get(harris_section, "DENSITY_SOURCE", &
                                   i_val=harris_env%density_source)
         CALL section_vals_val_get(harris_section, "ORBITAL_BASIS", &
                                   i_val=harris_env%orbital_basis)
         !
         CALL section_vals_val_get(harris_section, "DEBUG_FORCES", &
                                   l_val=harris_env%debug_forces)
         CALL section_vals_val_get(harris_section, "DEBUG_STRESS", &
                                   l_val=harris_env%debug_stress)

      END IF

      CALL timestop(handle)

   END SUBROUTINE init_harris_env

! **************************************************************************************************
!> \brief Print out the Harris method input section
!>
!> \param harris_env ...
!> \par History
!>       2024.07 created [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE harris_write_input(harris_env)
      TYPE(harris_type), POINTER                         :: harris_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'harris_write_input'

      INTEGER                                            :: handle, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, '(/,T2,A)') &
            "!"//REPEAT("-", 29)//"   Harris Model    "//REPEAT("-", 29)//"!"

         ! Type of energy functional
         SELECT CASE (harris_env%energy_functional)
         CASE (hfun_harris)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Energy Functional: ", "Harris"
         END SELECT
         ! density source
         SELECT CASE (harris_env%density_source)
         CASE (hden_atomic)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Harris model density: Type", " Atomic kind density"
         END SELECT
         WRITE (unit_nr, '(T2,A,T71,A10)') "Harris model density: Basis type", &
            ADJUSTR(TRIM(harris_env%rhoin%basis_type))
         WRITE (unit_nr, '(T2,A,T71,I10)') "Harris model density: Number of basis functions", &
            harris_env%rhoin%nbas
         ! orbital basis
         SELECT CASE (harris_env%orbital_basis)
         CASE (horb_default)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Harris model basis: ", "Atomic kind orbitals"
         END SELECT

         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
         WRITE (unit_nr, '()')

      END IF ! unit_nr

      CALL timestop(handle)

   END SUBROUTINE harris_write_input

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param harris_env ...
! **************************************************************************************************
   SUBROUTINE harris_density_update(qs_env, harris_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(harris_type), POINTER                         :: harris_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'harris_density_update'

      INTEGER                                            :: handle, i, ikind, ngto, nkind, nset, nsgf
      INTEGER, DIMENSION(:), POINTER                     :: lmax, npgf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coef
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: density
      REAL(KIND=dp), DIMENSION(:), POINTER               :: norm
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      SELECT CASE (harris_env%density_source)
      CASE (hden_atomic)
         IF (.NOT. harris_env%rhoin%frozen) THEN
            CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                            nkind=nkind)
            DO ikind = 1, nkind
               atomic_kind => atomic_kind_set(ikind)
               qs_kind => qs_kind_set(ikind)
               CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, &
                                basis_type=harris_env%rhoin%basis_type)
               CALL get_gto_basis_set(gto_basis_set=basis_set, nset=nset, lmax=lmax, nsgf=nsgf, &
                                      npgf=npgf, norm_cgf=norm, zet=zet, gcc=gcc)
               IF (nset /= 1 .OR. lmax(1) /= 0 .OR. npgf(1) /= nsgf) THEN
                  CPABORT("RHOIN illegal basis type")
               END IF
               DO i = 1, npgf(1)
                  IF (SUM(ABS(gcc(1:npgf(1), i, 1))) /= MAXVAL(ABS(gcc(1:npgf(1), i, 1)))) THEN
                     CPABORT("RHOIN illegal basis type")
                  END IF
               END DO
               !
               ngto = npgf(1)
               ALLOCATE (density(ngto, 2))
               density(1:ngto, 1) = zet(1:ngto, 1)
               density(1:ngto, 2) = 0.0_dp
               CALL calculate_atomic_density(density, atomic_kind, qs_kind, ngto, &
                                             optbasis=.FALSE., confine=.TRUE.)
               ALLOCATE (coef(ngto))
               DO i = 1, ngto
                  coef(i) = density(i, 2)/gcc(i, i, 1)/norm(i)
               END DO
               IF (harris_env%rhoin%nspin == 2) THEN
                  DO i = 1, SIZE(harris_env%rhoin%rhovec(ikind, 1)%rvecs, 2)
                     harris_env%rhoin%rhovec(ikind, 1)%rvecs(1:ngto, i) = coef(1:ngto)*0.5_dp
                     harris_env%rhoin%rhovec(ikind, 2)%rvecs(1:ngto, i) = coef(1:ngto)*0.5_dp
                  END DO
               ELSE
                  DO i = 1, SIZE(harris_env%rhoin%rhovec(ikind, 1)%rvecs, 2)
                     harris_env%rhoin%rhovec(ikind, 1)%rvecs(1:ngto, i) = coef(1:ngto)
                  END DO
               END IF
               DEALLOCATE (density, coef)
            END DO
         END IF
      CASE DEFAULT
         CPABORT("Illeagal value of harris_env%density_source")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE harris_density_update

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rhoin ...
!> \param rho_struct ...
! **************************************************************************************************
   SUBROUTINE calculate_harris_density(qs_env, rhoin, rho_struct)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(harris_rhoin_type), INTENT(IN)                :: rhoin
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_harris_density'

      INTEGER                                            :: handle, i1, i2, iatom, ikind, ilocal, &
                                                            ispin, n, nkind, nlocal, nspin
      REAL(KIND=dp)                                      :: eps_rho_rspace
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vector
      REAL(KIND=dp), DIMENSION(:), POINTER               :: total_rho
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_gspace
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_rspace
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control, para_env=para_env)
      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
      CALL get_qs_env(qs_env, &
                      atomic_kind_set=atomic_kind_set, particle_set=particle_set, &
                      local_particles=local_particles, &
                      qs_kind_set=qs_kind_set, cell=cell, pw_env=pw_env)

      CALL qs_rho_get(rho_struct, rho_r=rho_rspace, rho_g=rho_gspace, &
                      tot_rho_r=total_rho)

      ALLOCATE (vector(rhoin%nbas))

      nkind = SIZE(rhoin%rhovec, 1)
      nspin = SIZE(rhoin%rhovec, 2)

      DO ispin = 1, nspin
         vector = 0.0_dp
         DO ikind = 1, nkind
            nlocal = local_particles%n_el(ikind)
            DO ilocal = 1, nlocal
               iatom = local_particles%list(ikind)%array(ilocal)
               i1 = rhoin%basptr(iatom, 1)
               i2 = rhoin%basptr(iatom, 2)
               n = i2 - i1 + 1
               vector(i1:i2) = rhoin%rhovec(ikind, ispin)%rvecs(1:n, ilocal)
            END DO
         END DO
         CALL para_env%sum(vector)
         !
         CALL collocate_function(vector, rho_rspace(ispin), rho_gspace(ispin), &
                                 atomic_kind_set, qs_kind_set, cell, particle_set, pw_env, &
                                 eps_rho_rspace, rhoin%basis_type)
         total_rho(ispin) = pw_integrate_function(rho_rspace(ispin), isign=-1)
      END DO

      DEALLOCATE (vector)

      CALL timestop(handle)

   END SUBROUTINE calculate_harris_density

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rhoin ...
!> \param v_rspace ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE calculate_harris_integrals(qs_env, rhoin, v_rspace, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(harris_rhoin_type), INTENT(INOUT)             :: rhoin
      TYPE(pw_r3d_rs_type), DIMENSION(:), INTENT(IN)     :: v_rspace
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_harris_integrals'

      INTEGER                                            :: handle, i1, i2, iatom, ikind, ilocal, &
                                                            ispin, n, nkind, nlocal, nspin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: integral, vector
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, para_env=para_env, local_particles=local_particles)

      ALLOCATE (vector(rhoin%nbas))
      ALLOCATE (integral(rhoin%nbas))

      nkind = SIZE(rhoin%rhovec, 1)
      nspin = SIZE(rhoin%rhovec, 2)

      DO ispin = 1, nspin
         vector = 0.0_dp
         integral = 0.0_dp
         DO ikind = 1, nkind
            nlocal = local_particles%n_el(ikind)
            DO ilocal = 1, nlocal
               iatom = local_particles%list(ikind)%array(ilocal)
               i1 = rhoin%basptr(iatom, 1)
               i2 = rhoin%basptr(iatom, 2)
               n = i2 - i1 + 1
               vector(i1:i2) = rhoin%rhovec(ikind, ispin)%rvecs(1:n, ilocal)
            END DO
         END DO
         CALL para_env%sum(vector)
         !
         CALL integrate_function(qs_env, v_rspace(ispin), vector, integral, &
                                 calculate_forces, rhoin%basis_type)
         DO ikind = 1, nkind
            nlocal = local_particles%n_el(ikind)
            DO ilocal = 1, nlocal
               iatom = local_particles%list(ikind)%array(ilocal)
               i1 = rhoin%basptr(iatom, 1)
               i2 = rhoin%basptr(iatom, 2)
               n = i2 - i1 + 1
               rhoin%intvec(ikind, ispin)%rvecs(1:n, ilocal) = integral(i1:i2)
            END DO
         END DO
      END DO

      DEALLOCATE (vector, integral)

      CALL timestop(handle)

   END SUBROUTINE calculate_harris_integrals

! **************************************************************************************************
!> \brief ...
!> \param harris_env ...
!> \param vh_rspace ...
!> \param vxc_rspace ...
! **************************************************************************************************
   SUBROUTINE harris_set_potentials(harris_env, vh_rspace, vxc_rspace)
      TYPE(harris_type), POINTER                         :: harris_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vh_rspace
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: vxc_rspace

      INTEGER                                            :: iab, ispin, nspins
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      ! release possible old potentials
      IF (ASSOCIATED(harris_env%vh_rspace%pw_grid)) THEN
         CALL harris_env%vh_rspace%release()
      END IF
      IF (ASSOCIATED(harris_env%vxc_rspace)) THEN
         DO iab = 1, SIZE(harris_env%vxc_rspace)
            CALL harris_env%vxc_rspace(iab)%release()
         END DO
         DEALLOCATE (harris_env%vxc_rspace)
      END IF

      ! generate new potential data structures
      nspins = harris_env%rhoin%nspin
      ALLOCATE (harris_env%vxc_rspace(nspins))

      pw_grid => vh_rspace%pw_grid
      CALL harris_env%vh_rspace%create(pw_grid)
      DO ispin = 1, nspins
         CALL harris_env%vxc_rspace(ispin)%create(pw_grid)
      END DO

      ! copy potentials
      CALL pw_transfer(vh_rspace, harris_env%vh_rspace)
      IF (ASSOCIATED(vxc_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_transfer(vxc_rspace(ispin), harris_env%vxc_rspace(ispin))
            CALL pw_scale(harris_env%vxc_rspace(ispin), vxc_rspace(ispin)%pw_grid%dvol)
         END DO
      ELSE
         DO ispin = 1, nspins
            CALL pw_zero(harris_env%vxc_rspace(ispin))
         END DO
      END IF

   END SUBROUTINE harris_set_potentials

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE harris_energy_correction(qs_env, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'harris_energy_correction'

      INTEGER                                            :: handle, iounit, ispin, nspins
      REAL(KIND=dp)                                      :: dvol, ec, eh, exc, vxc
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(harris_energy_type), POINTER                  :: energy
      TYPE(harris_type), POINTER                         :: harris_env
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_core
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: core_rspace
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(qs_energy_type), POINTER                      :: ks_energy
      TYPE(qs_rho_type), POINTER                         :: rho

      MARK_USED(calculate_forces)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, harris_env=harris_env, energy=ks_energy)
      energy => harris_env%energy
      energy%eband = ks_energy%band
      energy%ewald_correction = ks_energy%core_overlap + ks_energy%core_self
      energy%dispersion = ks_energy%dispersion

      nspins = harris_env%rhoin%nspin

      CALL get_qs_env(qs_env, rho=rho, rho_core=rho_core)
      CALL qs_rho_get(rho, rho_r=rho_r)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL auxbas_pw_pool%create_pw(core_rspace)
      CALL pw_transfer(rho_core, core_rspace)

      dvol = harris_env%vh_rspace%pw_grid%dvol
      eh = 0.0_dp
      DO ispin = 1, nspins
         eh = eh + pw_integral_ab(rho_r(ispin), harris_env%vh_rspace)/dvol
      END DO
      ec = pw_integral_ab(core_rspace, harris_env%vh_rspace)/dvol
      eh = 0.5_dp*(eh + ec)
      energy%eh_correction = ec - eh

      exc = ks_energy%exc
      vxc = 0.0_dp
      IF (ASSOCIATED(harris_env%vxc_rspace)) THEN
         DO ispin = 1, nspins
            vxc = vxc + pw_integral_ab(rho_r(ispin), harris_env%vxc_rspace(ispin))/ &
                  harris_env%vxc_rspace(ispin)%pw_grid%dvol
         END DO
      END IF
      energy%exc_correction = exc - vxc

      ! Total Harris model energy
      energy%eharris = energy%eband + energy%eh_correction + energy%exc_correction + &
                       energy%ewald_correction + energy%dispersion

      CALL auxbas_pw_pool%give_back_pw(core_rspace)

      ks_energy%total = ks_energy%total + ks_energy%core
      ks_energy%nonscf_correction = energy%eharris - ks_energy%total
      ks_energy%total = energy%eharris

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CALL harris_print_energy(iounit, energy)

      IF (calculate_forces) THEN
         CALL harris_forces(qs_env, iounit)
      END IF

      CALL timestop(handle)

   END SUBROUTINE harris_energy_correction

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param iounit ...
! **************************************************************************************************
   SUBROUTINE harris_forces(qs_env, iounit)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: iounit

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'harris_forces'
      LOGICAL, PARAMETER                                 :: debug_forces = .TRUE.

      INTEGER                                            :: handle, ispin, nspins
      REAL(KIND=dp)                                      :: ehartree
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      TYPE(dbcsr_p_type)                                 :: scrm
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rhoh_ao, smat
      TYPE(harris_type), POINTER                         :: harris_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_c1d_gs_type)                               :: rhoh_tot_gspace, vhout_gspace
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g, rhoh_g
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_core
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: vhout_rspace, vhxc_rspace
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fhxc_rspace, ftau, fxc, rho_r, rhoh_r, &
                                                            tauh_r
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: xc_section

      CALL timeset(routineN, handle)

      IF (debug_forces) THEN
         IF (iounit > 0) WRITE (iounit, "(/,T3,A)") &
            "DEBUG:: Harris Method Forces (density dependent)"
      END IF

      CALL get_qs_env(qs_env, harris_env=harris_env, force=force, para_env=para_env)
      nspins = harris_env%rhoin%nspin

      CALL get_qs_env(qs_env, rho=rho, rho_core=rho_core, matrix_s=smat)
      ! Warning: rho_ao = output DM; rho_r = rhoin
      CALL qs_rho_get(rho, rho_ao=rhoh_ao, rho_r=rho_r, rho_g=rho_g)
      ALLOCATE (scrm%matrix)
      CALL dbcsr_create(scrm%matrix, template=rhoh_ao(1)%matrix)
      CALL dbcsr_copy(scrm%matrix, smat(1)%matrix)
      CALL dbcsr_set(scrm%matrix, 0.0_dp)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, ks_env=ks_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL auxbas_pw_pool%create_pw(vhxc_rspace)

      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      DO ispin = 1, nspins
         CALL pw_copy(harris_env%vh_rspace, vhxc_rspace)
         CALL pw_axpy(harris_env%vxc_rspace(ispin), vhxc_rspace)
         CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                 hmat=scrm, pmat=rhoh_ao(ispin), &
                                 qs_env=qs_env, calculate_forces=.TRUE.)
      END DO
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*(Vh[in]+Vxc)", fodeb
      END IF

      CALL dbcsr_release(scrm%matrix)
      DEALLOCATE (scrm%matrix)
      CALL auxbas_pw_pool%give_back_pw(vhxc_rspace)

      ALLOCATE (rhoh_r(nspins), rhoh_g(nspins))
      DO ispin = 1, nspins
         CALL auxbas_pw_pool%create_pw(rhoh_r(ispin))
         CALL auxbas_pw_pool%create_pw(rhoh_g(ispin))
      END DO
      CALL auxbas_pw_pool%create_pw(rhoh_tot_gspace)
      CALL pw_copy(rho_core, rhoh_tot_gspace)
      DO ispin = 1, nspins
         CALL calculate_rho_elec(ks_env=ks_env, matrix_p=rhoh_ao(ispin)%matrix, &
                                 rho=rhoh_r(ispin), rho_gspace=rhoh_g(ispin))
         CALL pw_axpy(rhoh_g(ispin), rhoh_tot_gspace)
      END DO
      ! no meta functionals here
      NULLIFY (tauh_r)

      CALL auxbas_pw_pool%create_pw(vhout_rspace)
      CALL auxbas_pw_pool%create_pw(vhout_gspace)
      CALL pw_env_get(pw_env, poisson_env=poisson_env)
      !
      CALL pw_poisson_solve(poisson_env, rhoh_tot_gspace, ehartree, vhout_gspace)
      !
      CALL pw_transfer(vhout_gspace, vhout_rspace)
      CALL pw_scale(vhout_rspace, vhout_rspace%pw_grid%dvol)

      IF (debug_forces) fodeb(1:3) = force(1)%rho_core(1:3, 1)
      CALL integrate_v_core_rspace(vhout_rspace, qs_env)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_core(1:3, 1) - fodeb(1:3)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Vh[out]*dncore ", fodeb
      END IF

      ALLOCATE (fhxc_rspace(nspins))
      DO ispin = 1, nspins
         CALL auxbas_pw_pool%create_pw(fhxc_rspace(ispin))
      END DO
      ! vh = vh[out] - vh[in]
      CALL pw_axpy(harris_env%vh_rspace, vhout_rspace, alpha=-1._dp, beta=1.0_dp)
      ! kernel fxc
      ! drho = rho[out] - rho[in]
      DO ispin = 1, nspins
         CALL pw_axpy(rho_r(ispin), rhoh_r(ispin), alpha=-1._dp, beta=1.0_dp)
         CALL pw_axpy(rho_g(ispin), rhoh_g(ispin), alpha=-1._dp, beta=1.0_dp)
      END DO
      xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      NULLIFY (fxc, ftau)
      CALL create_kernel(qs_env, vxc=fxc, vxc_tau=ftau, &
                         rho=rho, rho1_r=rhoh_r, rho1_g=rhoh_g, tau1_r=tauh_r, &
                         xc_section=xc_section)
      CPASSERT(.NOT. ASSOCIATED(ftau))

      DO ispin = 1, nspins
         CALL pw_copy(vhout_rspace, fhxc_rspace(ispin))
         IF (ASSOCIATED(fxc)) THEN
            CALL pw_scale(fxc(ispin), fxc(ispin)%pw_grid%dvol)
            CALL pw_axpy(fxc(ispin), fhxc_rspace(ispin))
         END IF
      END DO

      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      CALL calculate_harris_integrals(qs_env, harris_env%rhoin, fhxc_rspace, .TRUE.)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: (dVh+fxc)*dn[in] ", fodeb
      END IF

      IF (ASSOCIATED(fxc)) THEN
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%give_back_pw(fxc(ispin))
         END DO
         DEALLOCATE (fxc)
      END IF
      IF (ASSOCIATED(ftau)) THEN
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%give_back_pw(ftau(ispin))
         END DO
         DEALLOCATE (ftau)
      END IF

      CALL auxbas_pw_pool%give_back_pw(rhoh_tot_gspace)
      CALL auxbas_pw_pool%give_back_pw(vhout_rspace)
      CALL auxbas_pw_pool%give_back_pw(vhout_gspace)

      DO ispin = 1, nspins
         CALL auxbas_pw_pool%give_back_pw(rhoh_r(ispin))
         CALL auxbas_pw_pool%give_back_pw(rhoh_g(ispin))
         CALL auxbas_pw_pool%give_back_pw(fhxc_rspace(ispin))
      END DO
      DEALLOCATE (rhoh_r, rhoh_g, fhxc_rspace)

      CALL timestop(handle)

   END SUBROUTINE harris_forces

END MODULE qs_harris_utils
