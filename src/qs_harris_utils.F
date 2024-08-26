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
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_constants,                 ONLY: hden_atomic,&
                                              hfun_harris,&
                                              horb_default
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_type
   USE pw_methods,                      ONLY: pw_integrate_function
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: collocate_function
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_harris_types,                 ONLY: harris_rhoin_type,&
                                              harris_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_harris_utils'

   PUBLIC :: harris_env_create, harris_write_input, harris_density_update, calculate_harris_density

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

END MODULE qs_harris_utils
