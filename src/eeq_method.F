!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of charge equilibration method
!> \author JGH
! **************************************************************************************************
MODULE eeq_method
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc,&
                                              plane_distance
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert,&
                                              cp_fm_matvec,&
                                              cp_fm_solve
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_unit_nr
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE eeq_data,                        ONLY: get_eeq_data
   USE eeq_input,                       ONLY: eeq_solver_type
   USE ewald_environment_types,         ONLY: ewald_env_create,&
                                              ewald_env_get,&
                                              ewald_env_release,&
                                              ewald_env_set,&
                                              ewald_environment_type,&
                                              read_ewald_section_tb
   USE ewald_pw_types,                  ONLY: ewald_pw_create,&
                                              ewald_pw_release,&
                                              ewald_pw_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: oorootpi,&
                                              twopi
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_types,                  ONLY: molecule_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: bohr
   USE pw_poisson_types,                ONLY: do_ewald_spme
   USE qs_dispersion_cnum,              ONLY: cnumber_init,&
                                              cnumber_release,&
                                              dcnum_type
   USE qs_dispersion_types,             ONLY: qs_dispersion_release,&
                                              qs_dispersion_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE spme,                            ONLY: spme_forces,&
                                              spme_potential,&
                                              spme_virial
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eeq_method'

   INTEGER, PARAMETER                                    :: maxElem = 86
   ! covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
   ! values for metals decreased by 10 %
   REAL(KIND=dp), PARAMETER :: rcov(1:maxElem) = [&
      & 0.32_dp, 0.46_dp, 1.20_dp, 0.94_dp, 0.77_dp, 0.75_dp, 0.71_dp, 0.63_dp, &
      & 0.64_dp, 0.67_dp, 1.40_dp, 1.25_dp, 1.13_dp, 1.04_dp, 1.10_dp, 1.02_dp, &
      & 0.99_dp, 0.96_dp, 1.76_dp, 1.54_dp, 1.33_dp, 1.22_dp, 1.21_dp, 1.10_dp, &
      & 1.07_dp, 1.04_dp, 1.00_dp, 0.99_dp, 1.01_dp, 1.09_dp, 1.12_dp, 1.09_dp, &
      & 1.15_dp, 1.10_dp, 1.14_dp, 1.17_dp, 1.89_dp, 1.67_dp, 1.47_dp, 1.39_dp, &
      & 1.32_dp, 1.24_dp, 1.15_dp, 1.13_dp, 1.13_dp, 1.08_dp, 1.15_dp, 1.23_dp, &
      & 1.28_dp, 1.26_dp, 1.26_dp, 1.23_dp, 1.32_dp, 1.31_dp, 2.09_dp, 1.76_dp, &
      & 1.62_dp, 1.47_dp, 1.58_dp, 1.57_dp, 1.56_dp, 1.55_dp, 1.51_dp, 1.52_dp, &
      & 1.51_dp, 1.50_dp, 1.49_dp, 1.49_dp, 1.48_dp, 1.53_dp, 1.46_dp, 1.37_dp, &
      & 1.31_dp, 1.23_dp, 1.18_dp, 1.16_dp, 1.11_dp, 1.12_dp, 1.13_dp, 1.32_dp, &
      & 1.30_dp, 1.30_dp, 1.36_dp, 1.31_dp, 1.38_dp, 1.42_dp]

   PUBLIC :: eeq_solver, eeq_print, eeq_charges, eeq_forces, &
             eeq_efield_energy, eeq_efield_pot, eeq_efield_force_loc, eeq_efield_force_periodic

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param iounit ...
!> \param print_level ...
!> \param ext ...
! **************************************************************************************************
   SUBROUTINE eeq_print(qs_env, iounit, print_level, ext)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: iounit, print_level
      LOGICAL, INTENT(IN)                                :: ext

      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: enshift_type, iatom, ikind, natom
      REAL(KIND=dp), DIMENSION(:), POINTER               :: charges
      TYPE(cell_type), POINTER                           :: cell
      TYPE(eeq_solver_type)                              :: eeq_sparam
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      MARK_USED(print_level)

      CALL get_qs_env(qs_env, natom=natom, particle_set=particle_set, cell=cell)
      IF (ext) THEN
         NULLIFY (charges)
         CALL get_qs_env(qs_env, eeq=charges)
         CPASSERT(ASSOCIATED(charges))
         enshift_type = 0
      ELSE
         ALLOCATE (charges(natom))
         ! enforce en shift method 1 (original/molecular)
         ! method 2 from paper on PBC seems not to work
         enshift_type = 1
         !IF (ALL(cell%perd == 0)) enshift_type = 1
         CALL eeq_charges(qs_env, charges, eeq_sparam, 2, enshift_type)
      END IF

      IF (iounit > 0) THEN

         IF (enshift_type == 0) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A)") "EEQ Charges (External)"
         ELSE IF (enshift_type == 1) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A)") "EEQ Charges (Parametrization 2019 (Molecules))"
         ELSE IF (enshift_type == 2) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A)") "EEQ Charges (Parametrization 2019 (Crystals))"
         ELSE
            CPABORT("Unknown enshift_type")
         END IF
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "#     Atom  Element     Kind            Atomic Charge"

         DO iatom = 1, natom
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                 element_symbol=element_symbol, &
                                 kind_number=ikind)
            WRITE (UNIT=iounit, FMT="(T4,I8,T18,A2,I10,T43,F12.4)") &
               iatom, element_symbol, ikind, charges(iatom)
         END DO

      END IF

      IF (.NOT. ext) DEALLOCATE (charges)

   END SUBROUTINE eeq_print

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
!> \param eeq_sparam ...
!> \param eeq_model ...
!> \param enshift_type ...
! **************************************************************************************************
   SUBROUTINE eeq_charges(qs_env, charges, eeq_sparam, eeq_model, enshift_type)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: charges
      TYPE(eeq_solver_type), INTENT(IN)                  :: eeq_sparam
      INTEGER, INTENT(IN)                                :: eeq_model, enshift_type

      CHARACTER(len=*), PARAMETER                        :: routineN = 'eeq_charges'

      INTEGER                                            :: handle, iatom, ikind, iunit, jkind, &
                                                            natom, nkind, za, zb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: do_ewald
      REAL(KIND=dp)                                      :: ala, alb, eeq_energy, esg, kappa, &
                                                            lambda, scn, sgamma, totalcharge, xi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: chia, cnumbers, efr, gam
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell, cell_ref
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dcnum_type), ALLOCATABLE, DIMENSION(:)        :: dcnum
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: ewald_section, poisson_section, &
                                                            print_section

      CALL timeset(routineN, handle)

      iunit = cp_logger_get_default_unit_nr()

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      para_env=para_env, &
                      blacs_env=blacs_env, &
                      cell=cell, &
                      dft_control=dft_control)
      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)

      totalcharge = dft_control%charge

      CALL get_cnumbers(qs_env, cnumbers, dcnum, .FALSE.)

      ! gamma[a,b]
      ALLOCATE (gab(nkind, nkind), gam(nkind))
      gab = 0.0_dp
      gam = 0.0_dp
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), zatom=za)
         CALL get_eeq_data(za, eeq_model, eta=gam(ikind), rad=ala)
         DO jkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(jkind), zatom=zb)
            CALL get_eeq_data(zb, eeq_model, rad=alb)
            !
            gab(ikind, jkind) = SQRT(1._dp/(ala*ala + alb*alb))
            !
         END DO
      END DO

      ! Chi[a,a]
      sgamma = 8.0_dp ! see D4 for periodic systems paper
      esg = 1.0_dp + EXP(sgamma)
      ALLOCATE (chia(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), zatom=za)
         CALL get_eeq_data(za, eeq_model, chi=xi, kcn=kappa)
         !
         IF (enshift_type == 1) THEN
            scn = cnumbers(iatom)/SQRT(cnumbers(iatom) + 1.0e-14_dp)
         ELSE IF (enshift_type == 2) THEN
            scn = LOG(esg/(esg - cnumbers(iatom)))
         ELSE
            CPABORT("Unknown enshift_type")
         END IF
         chia(iatom) = xi - kappa*scn
         !
      END DO

      ! efield
      IF (dft_control%apply_period_efield .OR. dft_control%apply_efield .OR. &
          dft_control%apply_efield_field) THEN
         ALLOCATE (efr(natom))
         efr(1:natom) = 0.0_dp
         CALL eeq_efield_pot(qs_env, efr)
         chia(1:natom) = chia(1:natom) + efr(1:natom)
         DEALLOCATE (efr)
      END IF

      CALL cnumber_release(cnumbers, dcnum, .FALSE.)

      CALL get_cell(cell, periodic=periodic)
      do_ewald = .NOT. ALL(periodic == 0)
      IF (do_ewald) THEN
         ALLOCATE (ewald_env)
         CALL ewald_env_create(ewald_env, para_env)
         poisson_section => section_vals_get_subs_vals(qs_env%input, "DFT%POISSON")
         CALL ewald_env_set(ewald_env, poisson_section=poisson_section)
         ewald_section => section_vals_get_subs_vals(poisson_section, "EWALD")
         print_section => section_vals_get_subs_vals(qs_env%input, "PRINT%GRID_INFORMATION")
         CALL get_qs_env(qs_env, cell_ref=cell_ref)
         CALL read_ewald_section_tb(ewald_env, ewald_section, cell_ref%hmat, &
                                    silent=.TRUE., pset="EEQ")
         ALLOCATE (ewald_pw)
         CALL ewald_pw_create(ewald_pw, ewald_env, cell, cell_ref, print_section=print_section)
         !
         CALL eeq_solver(charges, lambda, eeq_energy, &
                         particle_set, kind_of, cell, chia, gam, gab, &
                         para_env, blacs_env, dft_control, eeq_sparam, &
                         totalcharge=totalcharge, ewald=do_ewald, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw, iounit=iunit)
         !
         CALL ewald_env_release(ewald_env)
         CALL ewald_pw_release(ewald_pw)
         DEALLOCATE (ewald_env, ewald_pw)
      ELSE
         CALL eeq_solver(charges, lambda, eeq_energy, &
                         particle_set, kind_of, cell, chia, gam, gab, &
                         para_env, blacs_env, dft_control, eeq_sparam, &
                         totalcharge=totalcharge, iounit=iunit)
      END IF

      DEALLOCATE (gab, gam, chia)

      CALL timestop(handle)

   END SUBROUTINE eeq_charges

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
!> \param dcharges ...
!> \param gradient ...
!> \param stress ...
!> \param eeq_sparam ...
!> \param eeq_model ...
!> \param enshift_type ...
!> \param response_only ...
! **************************************************************************************************
   SUBROUTINE eeq_forces(qs_env, charges, dcharges, gradient, stress, &
                         eeq_sparam, eeq_model, enshift_type, response_only)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: charges, dcharges
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: gradient
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: stress
      TYPE(eeq_solver_type), INTENT(IN)                  :: eeq_sparam
      INTEGER, INTENT(IN)                                :: eeq_model, enshift_type
      LOGICAL, INTENT(IN)                                :: response_only

      CHARACTER(len=*), PARAMETER                        :: routineN = 'eeq_forces'

      INTEGER                                            :: handle, i, ia, iatom, ikind, iunit, &
                                                            jatom, jkind, katom, natom, nkind, za, &
                                                            zb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: do_ewald, use_virial
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: default_present
      REAL(KIND=dp)                                      :: ala, alb, alpha, cn, ctot, dr, dr2, drk, &
                                                            elag, esg, fe, gam2, gama, grc, kappa, &
                                                            qlam, qq, qq1, qq2, rcut, scn, sgamma, &
                                                            subcells, totalcharge, xi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: c_radius, cnumbers, gam, qlag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: epforce, gab, pair_radius
      REAL(KIND=dp), DIMENSION(3)                        :: fdik, ri, rij, rik, rj
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pvir
      REAL(KIND=dp), DIMENSION(:), POINTER               :: chrgx, dchia
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell, cell_ref
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dcnum_type), ALLOCATABLE, DIMENSION(:)        :: dcnum
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: distribution_1d, local_particles
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_ew
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: ewald_section, poisson_section, &
                                                            print_section
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      iunit = cp_logger_get_default_unit_nr()

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      para_env=para_env, &
                      blacs_env=blacs_env, &
                      cell=cell, &
                      force=force, &
                      virial=virial, &
                      atprop=atprop, &
                      dft_control=dft_control)
      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

      totalcharge = dft_control%charge

      CALL get_cnumbers(qs_env, cnumbers, dcnum, .TRUE.)

      ! gamma[a,b]
      ALLOCATE (gab(nkind, nkind), gam(nkind))
      gab = 0.0_dp
      gam = 0.0_dp
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), zatom=za)
         CALL get_eeq_data(za, eeq_model, eta=gam(ikind), rad=ala)
         DO jkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(jkind), zatom=zb)
            CALL get_eeq_data(zb, eeq_model, rad=alb)
            !
            gab(ikind, jkind) = SQRT(1._dp/(ala*ala + alb*alb))
            !
         END DO
      END DO

      ALLOCATE (qlag(natom))

      CALL get_cell(cell, periodic=periodic)
      do_ewald = .NOT. ALL(periodic == 0)
      IF (do_ewald) THEN
         ALLOCATE (ewald_env)
         CALL ewald_env_create(ewald_env, para_env)
         poisson_section => section_vals_get_subs_vals(qs_env%input, "DFT%POISSON")
         CALL ewald_env_set(ewald_env, poisson_section=poisson_section)
         ewald_section => section_vals_get_subs_vals(poisson_section, "EWALD")
         print_section => section_vals_get_subs_vals(qs_env%input, "PRINT%GRID_INFORMATION")
         CALL get_qs_env(qs_env, cell_ref=cell_ref)
         CALL read_ewald_section_tb(ewald_env, ewald_section, cell_ref%hmat, &
                                    silent=.TRUE., pset="EEQ")
         ALLOCATE (ewald_pw)
         CALL ewald_pw_create(ewald_pw, ewald_env, cell, cell_ref, print_section=print_section)
         !
         CALL eeq_solver(qlag, qlam, elag, &
                         particle_set, kind_of, cell, -dcharges, gam, gab, &
                         para_env, blacs_env, dft_control, eeq_sparam, &
                         ewald=do_ewald, ewald_env=ewald_env, ewald_pw=ewald_pw, iounit=iunit)
      ELSE
         CALL eeq_solver(qlag, qlam, elag, &
                         particle_set, kind_of, cell, -dcharges, gam, gab, &
                         para_env, blacs_env, dft_control, eeq_sparam, iounit=iunit)
      END IF

      sgamma = 8.0_dp ! see D4 for periodic systems paper
      esg = 1.0_dp + EXP(sgamma)
      ALLOCATE (chrgx(natom), dchia(natom))
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), zatom=za)
         CALL get_eeq_data(za, eeq_model, chi=xi, kcn=kappa)
         !
         IF (response_only) THEN
            ctot = -0.5_dp*qlag(iatom)
         ELSE
            ctot = 0.5_dp*(charges(iatom) - qlag(iatom))
         END IF
         IF (enshift_type == 1) THEN
            scn = SQRT(cnumbers(iatom)) + 1.0e-14_dp
            dchia(iatom) = -ctot*kappa/scn
         ELSE IF (enshift_type == 2) THEN
            cn = cnumbers(iatom)
            scn = 1.0_dp/(esg - cn)
            dchia(iatom) = -ctot*kappa*scn
         ELSE
            CPABORT("Unknown enshift_type")
         END IF
      END DO

      ! Efield
      IF (dft_control%apply_period_efield) THEN
         CALL eeq_efield_force_periodic(qs_env, charges, qlag)
      ELSE IF (dft_control%apply_efield) THEN
         CALL eeq_efield_force_loc(qs_env, charges, qlag)
      ELSE IF (dft_control%apply_efield_field) THEN
         CPABORT("apply field")
      END IF

      ! Forces from q*X
      CALL get_qs_env(qs_env=qs_env, local_particles=local_particles)
      DO ikind = 1, nkind
         DO ia = 1, local_particles%n_el(ikind)
            iatom = local_particles%list(ikind)%array(ia)
            DO i = 1, dcnum(iatom)%neighbors
               katom = dcnum(iatom)%nlist(i)
               rik = dcnum(iatom)%rik(:, i)
               drk = SQRT(SUM(rik(:)**2))
               IF (drk > 1.e-3_dp) THEN
                  fdik(:) = dchia(iatom)*dcnum(iatom)%dvals(i)*rik(:)/drk
                  gradient(:, iatom) = gradient(:, iatom) - fdik(:)
                  gradient(:, katom) = gradient(:, katom) + fdik(:)
                  IF (use_virial) THEN
                     CALL virial_pair_force(stress, 1._dp, fdik, rik)
                  END IF
               END IF
            END DO
         END DO
      END DO

      ! Forces from (0.5*q+l)*dA/dR*q
      IF (do_ewald) THEN

         ! Build the neighbor lists for the CN
         CALL get_qs_env(qs_env, &
                         distribution_2d=distribution_2d, &
                         local_particles=distribution_1d, &
                         molecule_set=molecule_set)
         subcells = 2.0_dp
         CALL ewald_env_get(ewald_env, alpha=alpha, rcut=rcut)
         rcut = 2.0_dp*rcut
         NULLIFY (sab_ew)
         ALLOCATE (c_radius(nkind), default_present(nkind), pair_radius(nkind, nkind))
         c_radius(:) = rcut
         default_present = .TRUE.
         ALLOCATE (atom2d(nkind))
         CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                           molecule_set, .FALSE., particle_set=particle_set)
         CALL pair_radius_setup(default_present, default_present, c_radius, c_radius, pair_radius)
         CALL build_neighbor_lists(sab_ew, particle_set, atom2d, cell, pair_radius, &
                                   subcells=subcells, operator_type="PP", nlname="sab_ew")
         DEALLOCATE (c_radius, pair_radius, default_present)
         CALL atom2d_cleanup(atom2d)
         !
         CALL neighbor_list_iterator_create(nl_iterator, sab_ew)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                   iatom=iatom, jatom=jatom, r=rij)
            !
            dr2 = SUM(rij**2)
            dr = SQRT(dr2)
            IF (dr > rcut .OR. dr < 1.E-6_dp) CYCLE
            fe = 1.0_dp
            IF (iatom == jatom) fe = 0.5_dp
            IF (response_only) THEN
               qq = -qlag(iatom)*charges(jatom)
            ELSE
               qq = (0.5_dp*charges(iatom) - qlag(iatom))*charges(jatom)
            END IF
            gama = gab(ikind, jkind)
            gam2 = gama*gama
            grc = 2._dp*gama*EXP(-gam2*dr2)*oorootpi/dr - erf(gama*dr)/dr2 &
                  - 2._dp*alpha*EXP(-alpha**2*dr2)*oorootpi/dr + erf(alpha*dr)/dr2
            IF (response_only) THEN
               qq1 = -qlag(iatom)*charges(jatom)
               qq2 = -qlag(jatom)*charges(iatom)
            ELSE
               qq1 = (0.5_dp*charges(iatom) - qlag(iatom))*charges(jatom)
               qq2 = (0.5_dp*charges(jatom) - qlag(jatom))*charges(iatom)
            END IF
            fdik(:) = -qq1*grc*rij(:)/dr
            gradient(:, iatom) = gradient(:, iatom) + fdik(:)
            gradient(:, jatom) = gradient(:, jatom) - fdik(:)
            IF (use_virial) THEN
               CALL virial_pair_force(stress, -fe, fdik, rij)
            END IF
            fdik(:) = qq2*grc*rij(:)/dr
            gradient(:, iatom) = gradient(:, iatom) - fdik(:)
            gradient(:, jatom) = gradient(:, jatom) + fdik(:)
            IF (use_virial) THEN
               CALL virial_pair_force(stress, fe, fdik, rij)
            END IF
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)
         !
         CALL release_neighbor_list_sets(sab_ew)
      ELSE
         DO ikind = 1, nkind
            DO ia = 1, local_particles%n_el(ikind)
               iatom = local_particles%list(ikind)%array(ia)
               ri(1:3) = particle_set(iatom)%r(1:3)
               DO jatom = 1, natom
                  IF (iatom == jatom) CYCLE
                  jkind = kind_of(jatom)
                  IF (response_only) THEN
                     qq = -qlag(iatom)*charges(jatom)
                  ELSE
                     qq = (0.5_dp*charges(iatom) - qlag(iatom))*charges(jatom)
                  END IF
                  rj(1:3) = particle_set(jatom)%r(1:3)
                  rij(1:3) = ri(1:3) - rj(1:3)
                  rij = pbc(rij, cell)
                  dr2 = SUM(rij**2)
                  dr = SQRT(dr2)
                  gama = gab(ikind, jkind)
                  gam2 = gama*gama
                  grc = 2._dp*gama*EXP(-gam2*dr2)*oorootpi/dr - erf(gama*dr)/dr2
                  fdik(:) = qq*grc*rij(:)/dr
                  gradient(:, iatom) = gradient(:, iatom) + fdik(:)
                  gradient(:, jatom) = gradient(:, jatom) - fdik(:)
               END DO
            END DO
         END DO
      END IF

      ! Forces from Ewald potential: (q+l)*A*q
      IF (do_ewald) THEN
         ALLOCATE (epforce(3, natom))
         epforce = 0.0_dp
         IF (response_only) THEN
            dchia(1:natom) = qlag(1:natom)
         ELSE
            dchia(1:natom) = -charges(1:natom) + qlag(1:natom)
         END IF
         chrgx(1:natom) = charges(1:natom)
         CALL spme_forces(ewald_env, ewald_pw, cell, particle_set, chrgx, &
                          particle_set, dchia, epforce)
         dchia(1:natom) = charges(1:natom)
         chrgx(1:natom) = qlag(1:natom)
         CALL spme_forces(ewald_env, ewald_pw, cell, particle_set, chrgx, &
                          particle_set, dchia, epforce)
         gradient(1:3, 1:natom) = gradient(1:3, 1:natom) + epforce(1:3, 1:natom)
         DEALLOCATE (epforce)

         ! virial
         IF (use_virial) THEN
            chrgx(1:natom) = charges(1:natom) - qlag(1:natom)
            CALL spme_virial(ewald_env, ewald_pw, particle_set, cell, chrgx, pvir)
            stress = stress - pvir
            chrgx(1:natom) = qlag(1:natom)
            CALL spme_virial(ewald_env, ewald_pw, particle_set, cell, chrgx, pvir)
            stress = stress + pvir
            IF (response_only) THEN
               chrgx(1:natom) = charges(1:natom)
               CALL spme_virial(ewald_env, ewald_pw, particle_set, cell, chrgx, pvir)
               stress = stress + pvir
            END IF
         END IF
         !
         CALL ewald_env_release(ewald_env)
         CALL ewald_pw_release(ewald_pw)
         DEALLOCATE (ewald_env, ewald_pw)
      END IF

      CALL cnumber_release(cnumbers, dcnum, .TRUE.)

      DEALLOCATE (gab, gam, qlag, chrgx, dchia)

      CALL timestop(handle)

   END SUBROUTINE eeq_forces

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cnumbers ...
!> \param dcnum ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE get_cnumbers(qs_env, cnumbers, dcnum, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cnumbers
      TYPE(dcnum_type), ALLOCATABLE, DIMENSION(:)        :: dcnum
      LOGICAL, INTENT(IN)                                :: calculate_forces

      INTEGER                                            :: ikind, natom, nkind, za
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: default_present
      REAL(KIND=dp)                                      :: subcells
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: c_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_cn
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_dispersion_type), POINTER                  :: disp
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      cell=cell, &
                      distribution_2d=distribution_2d, &
                      local_particles=distribution_1d, &
                      molecule_set=molecule_set)
      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)

      ! Check for dispersion_env and sab_cn needed for cnumbers
      ALLOCATE (disp)
      disp%k1 = 16.0_dp
      disp%k2 = 4._dp/3._dp
      disp%eps_cn = 1.E-6_dp
      disp%max_elem = maxElem
      ALLOCATE (disp%rcov(maxElem))
      disp%rcov(1:maxElem) = bohr*disp%k2*rcov(1:maxElem)
      subcells = 2.0_dp
      ! Build the neighbor lists for the CN
      NULLIFY (sab_cn)
      ALLOCATE (c_radius(nkind), default_present(nkind), pair_radius(nkind, nkind))
      c_radius(:) = 0.0_dp
      default_present = .TRUE.
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), z=za)
         c_radius(ikind) = 4._dp*rcov(za)*bohr
      END DO
      ALLOCATE (atom2d(nkind))
      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, .FALSE., particle_set=particle_set)
      CALL pair_radius_setup(default_present, default_present, c_radius, c_radius, pair_radius)
      CALL build_neighbor_lists(sab_cn, particle_set, atom2d, cell, pair_radius, &
                                subcells=subcells, operator_type="PP", nlname="sab_cn")
      disp%sab_cn => sab_cn
      DEALLOCATE (c_radius, pair_radius, default_present)
      CALL atom2d_cleanup(atom2d)

      ! Calculate coordination numbers
      CALL cnumber_init(qs_env, cnumbers, dcnum, 2, calculate_forces, disp_env=disp)

      CALL qs_dispersion_release(disp)

   END SUBROUTINE get_cnumbers

! **************************************************************************************************
!> \brief ...
!> \param charges ...
!> \param lambda ...
!> \param eeq_energy ...
!> \param particle_set ...
!> \param kind_of ...
!> \param cell ...
!> \param chia ...
!> \param gam ...
!> \param gab ...
!> \param para_env ...
!> \param blacs_env ...
!> \param dft_control ...
!> \param eeq_sparam ...
!> \param totalcharge ...
!> \param ewald ...
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param iounit ...
! **************************************************************************************************
   SUBROUTINE eeq_solver(charges, lambda, eeq_energy, particle_set, kind_of, cell, &
                         chia, gam, gab, para_env, blacs_env, dft_control, eeq_sparam, &
                         totalcharge, ewald, ewald_env, ewald_pw, iounit)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: charges
      REAL(KIND=dp), INTENT(INOUT)                       :: lambda, eeq_energy
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: chia, gam
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: gab
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(eeq_solver_type), INTENT(IN)                  :: eeq_sparam
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: totalcharge
      LOGICAL, INTENT(IN), OPTIONAL                      :: ewald
      TYPE(ewald_environment_type), OPTIONAL, POINTER    :: ewald_env
      TYPE(ewald_pw_type), OPTIONAL, POINTER             :: ewald_pw
      INTEGER, INTENT(IN), OPTIONAL                      :: iounit

      CHARACTER(len=*), PARAMETER                        :: routineN = 'eeq_solver'

      INTEGER                                            :: handle, ierror, iunit, natom, nkind, ns
      LOGICAL                                            :: do_direct, do_ewald, do_sparse
      REAL(KIND=dp)                                      :: alpha, deth, ftime, qtot
      TYPE(cp_fm_struct_type), POINTER                   :: mat_struct
      TYPE(cp_fm_type)                                   :: eeq_mat

      CALL timeset(routineN, handle)

      do_ewald = .FALSE.
      IF (PRESENT(ewald)) do_ewald = ewald
      !
      qtot = 0.0_dp
      IF (PRESENT(totalcharge)) qtot = totalcharge
      !
      iunit = -1
      IF (PRESENT(iounit)) iunit = iounit

      ! EEQ solver parameters
      do_direct = eeq_sparam%direct
      do_sparse = eeq_sparam%sparse

      ! sparse option NYA
      IF (do_sparse) THEN
         CALL cp_abort(__LOCATION__, "EEQ: Sparse option not yet available")
      END IF
      !
      natom = SIZE(particle_set)
      nkind = SIZE(gam)
      ns = natom + 1
      CALL cp_fm_struct_create(mat_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ns, ncol_global=ns)
      CALL cp_fm_create(eeq_mat, mat_struct)
      CALL cp_fm_set_all(eeq_mat, 0.0_dp, 0.0_dp)
      !
      IF (do_ewald) THEN
         CPASSERT(PRESENT(ewald_env))
         CPASSERT(PRESENT(ewald_pw))
         IF (do_direct) THEN
            IF (dft_control%apply_period_efield .AND. dft_control%period_efield%displacement_field) THEN
               CPABORT("NYA")
            ELSE
               CALL fpbc_solver(charges, lambda, eeq_energy, eeq_mat, particle_set, &
                                kind_of, cell, chia, gam, gab, qtot, &
                                ewald_env, ewald_pw, iounit)
            END IF
         ELSE
            IF (dft_control%apply_period_efield .AND. dft_control%period_efield%displacement_field) THEN
               CPABORT("NYA")
            ELSE
               ierror = 0
               CALL pbc_solver(charges, lambda, eeq_energy, eeq_mat, particle_set, &
                               kind_of, cell, chia, gam, gab, qtot, &
                               ewald_env, ewald_pw, eeq_sparam, ierror, iounit)
               IF (ierror /= 0) THEN
                  ! backup to non-iterative method
                  CALL fpbc_solver(charges, lambda, eeq_energy, eeq_mat, particle_set, &
                                   kind_of, cell, chia, gam, gab, qtot, &
                                   ewald_env, ewald_pw, iounit)
               END IF
            END IF
         END IF
         IF (qtot /= 0._dp) THEN
            CALL get_cell(cell=cell, deth=deth)
            CALL ewald_env_get(ewald_env, alpha=alpha)
            eeq_energy = eeq_energy - 0.5_dp*qtot**2/alpha**2/deth
         END IF
      ELSE
         CALL mi_solver(charges, lambda, eeq_energy, eeq_mat, particle_set, kind_of, &
                        cell, chia, gam, gab, qtot, ftime)
         IF (iounit > 0) THEN
            WRITE (iunit, '(A,T67,F14.3)') " EEQ| Molecular solver time[s]", ftime
         END IF
      END IF
      CALL cp_fm_struct_release(mat_struct)
      CALL cp_fm_release(eeq_mat)

      CALL timestop(handle)

   END SUBROUTINE eeq_solver

! **************************************************************************************************
!> \brief ...
!> \param charges ...
!> \param lambda ...
!> \param eeq_energy ...
!> \param eeq_mat ...
!> \param particle_set ...
!> \param kind_of ...
!> \param cell ...
!> \param chia ...
!> \param gam ...
!> \param gab ...
!> \param qtot ...
!> \param ftime ...
! **************************************************************************************************
   SUBROUTINE mi_solver(charges, lambda, eeq_energy, eeq_mat, particle_set, kind_of, cell, &
                        chia, gam, gab, qtot, ftime)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: charges
      REAL(KIND=dp), INTENT(INOUT)                       :: lambda, eeq_energy
      TYPE(cp_fm_type)                                   :: eeq_mat
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: chia, gam
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: gab
      REAL(KIND=dp), INTENT(IN)                          :: qtot
      REAL(KIND=dp), INTENT(OUT)                         :: ftime

      CHARACTER(len=*), PARAMETER                        :: routineN = 'mi_solver'

      INTEGER                                            :: handle, ia, iac, iar, ic, ikind, ir, &
                                                            jkind, natom, ncloc, ncvloc, nkind, &
                                                            nrloc, nrvloc, ns
      INTEGER, DIMENSION(:), POINTER                     :: cind, cvind, rind, rvind
      REAL(KIND=dp)                                      :: dr, grc, te, ti, xr
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rj
      TYPE(cp_fm_struct_type), POINTER                   :: mat_struct, vec_struct
      TYPE(cp_fm_type)                                   :: rhs_vec
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      ti = m_walltime()

      natom = SIZE(particle_set)
      nkind = SIZE(gam)
      !
      ns = natom + 1
      CALL cp_fm_get_info(eeq_mat, matrix_struct=mat_struct, para_env=para_env)
      CALL cp_fm_get_info(eeq_mat, nrow_local=nrloc, ncol_local=ncloc, &
                          row_indices=rind, col_indices=cind)
      CALL cp_fm_struct_create(vec_struct, template_fmstruct=mat_struct, &
                               nrow_global=ns, ncol_global=1)
      CALL cp_fm_create(rhs_vec, vec_struct)
      CALL cp_fm_get_info(rhs_vec, nrow_local=nrvloc, ncol_local=ncvloc, &
                          row_indices=rvind, col_indices=cvind)
      !
      ! set up matrix
      CALL cp_fm_set_all(eeq_mat, 1.0_dp, 0.0_dp)
      CALL cp_fm_set_all(rhs_vec, 0.0_dp)
      DO ir = 1, nrloc
         iar = rind(ir)
         IF (iar > natom) CYCLE
         ikind = kind_of(iar)
         ri(1:3) = particle_set(iar)%r(1:3)
         DO ic = 1, ncloc
            iac = cind(ic)
            IF (iac > natom) CYCLE
            jkind = kind_of(iac)
            rj(1:3) = particle_set(iac)%r(1:3)
            IF (iar == iac) THEN
               grc = gam(ikind) + 2.0_dp*gab(ikind, ikind)*oorootpi
            ELSE
               rij(1:3) = ri(1:3) - rj(1:3)
               rij = pbc(rij, cell)
               dr = SQRT(SUM(rij**2))
               grc = erf(gab(ikind, jkind)*dr)/dr
            END IF
            eeq_mat%local_data(ir, ic) = grc
         END DO
      END DO
      ! set up rhs vector
      DO ir = 1, nrvloc
         iar = rvind(ir)
         DO ic = 1, ncvloc
            iac = cvind(ic)
            ia = MAX(iar, iac)
            IF (ia > natom) THEN
               xr = qtot
            ELSE
               xr = -chia(ia)
            END IF
            rhs_vec%local_data(ir, ic) = xr
         END DO
      END DO
      !
      CALL cp_fm_solve(eeq_mat, rhs_vec)
      !
      charges = 0.0_dp
      lambda = 0.0_dp
      DO ir = 1, nrvloc
         iar = rvind(ir)
         DO ic = 1, ncvloc
            iac = cvind(ic)
            ia = MAX(iar, iac)
            IF (ia <= natom) THEN
               xr = rhs_vec%local_data(ir, ic)
               charges(ia) = xr
            ELSE
               lambda = rhs_vec%local_data(ir, ic)
            END IF
         END DO
      END DO
      CALL para_env%sum(lambda)
      CALL para_env%sum(charges)
      !
      ! energy:   0.5*(q^T.X - lambda*totalcharge)
      eeq_energy = 0.5*SUM(charges(1:natom)*chia(1:natom)) - 0.5_dp*lambda*qtot

      CALL cp_fm_struct_release(vec_struct)
      CALL cp_fm_release(rhs_vec)

      te = m_walltime()
      ftime = te - ti
      CALL timestop(handle)

   END SUBROUTINE mi_solver

! **************************************************************************************************
!> \brief ...
!> \param charges ...
!> \param lambda ...
!> \param eeq_energy ...
!> \param eeq_mat ...
!> \param particle_set ...
!> \param kind_of ...
!> \param cell ...
!> \param chia ...
!> \param gam ...
!> \param gab ...
!> \param qtot ...
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param eeq_sparam ...
!> \param ierror ...
!> \param iounit ...
! **************************************************************************************************
   SUBROUTINE pbc_solver(charges, lambda, eeq_energy, eeq_mat, particle_set, &
                         kind_of, cell, chia, gam, gab, qtot, &
                         ewald_env, ewald_pw, eeq_sparam, ierror, iounit)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: charges
      REAL(KIND=dp), INTENT(INOUT)                       :: lambda, eeq_energy
      TYPE(cp_fm_type)                                   :: eeq_mat
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: chia, gam
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: gab
      REAL(KIND=dp), INTENT(IN)                          :: qtot
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(eeq_solver_type), INTENT(IN)                  :: eeq_sparam
      INTEGER, INTENT(OUT)                               :: ierror
      INTEGER, OPTIONAL                                  :: iounit

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pbc_solver'

      INTEGER :: ewald_type, handle, i, iac, iar, ic, ikind, info, ir, iunit, iv, ix, iy, iz, &
         jkind, max_diis, mdiis, natom, ncloc, ndiis, nkind, now, nrloc, ns, sdiis
      INTEGER, DIMENSION(3)                              :: cvec, ncell, periodic
      INTEGER, DIMENSION(:), POINTER                     :: cind, rind
      REAL(KIND=dp)                                      :: ad, alpha, astep, deth, dr, eeqn, &
                                                            eps_diis, ftime, grc1, grc2, rcut, &
                                                            res, resin, rmax, te, ti
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: bvec, dvec
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dmat, fvec, vmat, xvec
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rijl, rj
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rhs, rv0, xv0
      TYPE(cp_fm_struct_type), POINTER                   :: mat_struct
      TYPE(cp_fm_type)                                   :: mmat, pmat
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      ti = m_walltime()

      ierror = 0

      iunit = -1
      IF (PRESENT(iounit)) iunit = iounit

      natom = SIZE(particle_set)
      nkind = SIZE(gam)
      !
      CALL get_cell(cell=cell, deth=deth)
      CALL ewald_env_get(ewald_env, alpha=alpha, rcut=rcut, ewald_type=ewald_type)
      ad = 2.0_dp*alpha*oorootpi
      IF (ewald_type /= do_ewald_spme) THEN
         CALL cp_abort(__LOCATION__, "Only SPME Ewald method available with EEQ.")
      END IF
      !
      rmax = 2.0_dp*rcut
      ! max cells used
      CALL get_cell(cell, h=hmat, periodic=periodic)
      ncell(1) = CEILING(rmax/plane_distance(1, 0, 0, cell))
      ncell(2) = CEILING(rmax/plane_distance(0, 1, 0, cell))
      ncell(3) = CEILING(rmax/plane_distance(0, 0, 1, cell))
      IF (periodic(1) == 0) ncell(1) = 0
      IF (periodic(2) == 0) ncell(2) = 0
      IF (periodic(3) == 0) ncell(3) = 0
      !
      CALL mi_solver(charges, lambda, eeqn, eeq_mat, particle_set, kind_of, cell, &
                     chia, gam, gab, qtot, ftime)
      IF (iunit > 0) THEN
         WRITE (iunit, '(A,T67,F14.3)') " EEQ| Iterative PBC guess time[s]", ftime
      END IF
      CALL cp_fm_get_info(eeq_mat, matrix_struct=mat_struct, para_env=para_env)
      CALL cp_fm_create(pmat, mat_struct)
      CALL cp_fm_create(mmat, mat_struct)
      !
      ! response matrix
      CALL cp_fm_get_info(eeq_mat, nrow_local=nrloc, ncol_local=ncloc, &
                          row_indices=rind, col_indices=cind)
      CALL cp_fm_set_all(eeq_mat, 0.0_dp, 0.0_dp)
      DO ir = 1, nrloc
         iar = rind(ir)
         ri = 0.0_dp
         IF (iar <= natom) THEN
            ikind = kind_of(iar)
            ri(1:3) = particle_set(iar)%r(1:3)
         END IF
         DO ic = 1, ncloc
            iac = cind(ic)
            IF (iac > natom .AND. iar > natom) THEN
               eeq_mat%local_data(ir, ic) = 0.0_dp
               CYCLE
            ELSE IF ((iac > natom) .OR. (iar > natom)) THEN
               eeq_mat%local_data(ir, ic) = 1.0_dp
               CYCLE
            END IF
            jkind = kind_of(iac)
            rj(1:3) = particle_set(iac)%r(1:3)
            rij(1:3) = ri(1:3) - rj(1:3)
            rij = pbc(rij, cell)
            DO ix = -ncell(1), ncell(1)
               DO iy = -ncell(2), ncell(2)
                  DO iz = -ncell(3), ncell(3)
                     cvec = [ix, iy, iz]
                     rijl = rij + MATMUL(hmat, cvec)
                     dr = SQRT(SUM(rijl**2))
                     IF (dr > rmax) CYCLE
                     IF (iar == iac .AND. dr < 0.00001_dp) THEN
                        grc1 = gam(ikind) + 2.0_dp*gab(ikind, ikind)*oorootpi - ad
                     ELSE
                        grc1 = erf(gab(ikind, jkind)*dr)/dr - erf(alpha*dr)/dr
                     END IF
                     eeq_mat%local_data(ir, ic) = eeq_mat%local_data(ir, ic) + grc1
                  END DO
               END DO
            END DO
         END DO
      END DO
      !
      ! preconditioner
      CALL cp_fm_get_info(pmat, nrow_local=nrloc, ncol_local=ncloc, &
                          row_indices=rind, col_indices=cind)
      CALL cp_fm_set_all(pmat, 0.0_dp, 0.0_dp)
      DO ir = 1, nrloc
         iar = rind(ir)
         ri = 0.0_dp
         IF (iar <= natom) THEN
            ikind = kind_of(iar)
            ri(1:3) = particle_set(iar)%r(1:3)
         END IF
         DO ic = 1, ncloc
            iac = cind(ic)
            IF (iac > natom .AND. iar > natom) THEN
               pmat%local_data(ir, ic) = 0.0_dp
               CYCLE
            ELSE IF ((iac > natom) .OR. (iar > natom)) THEN
               pmat%local_data(ir, ic) = 1.0_dp
               CYCLE
            END IF
            jkind = kind_of(iac)
            rj(1:3) = particle_set(iac)%r(1:3)
            rij(1:3) = ri(1:3) - rj(1:3)
            rij = pbc(rij, cell)
            IF (iar == iac) THEN
               grc2 = gam(ikind) + 2.0_dp*gab(ikind, ikind)*oorootpi
            ELSE
               grc2 = erf(gab(ikind, jkind)*dr)/dr
            END IF
            pmat%local_data(ir, ic) = grc2
         END DO
      END DO
      CALL cp_fm_set_all(mmat, 0.0_dp, 0.0_dp)
      ! preconditioner invers
      CALL cp_fm_invert(pmat, mmat)
      !
      ! rhs
      ns = natom + 1
      ALLOCATE (rhs(ns))
      rhs(1:natom) = chia(1:natom)
      rhs(ns) = -qtot
      !
      ALLOCATE (xv0(ns), rv0(ns))
      ! initial guess
      xv0(1:natom) = charges(1:natom)
      xv0(ns) = 0.0_dp
      ! DIIS optimizer
      max_diis = eeq_sparam%max_diis
      mdiis = eeq_sparam%mdiis
      sdiis = eeq_sparam%sdiis
      eps_diis = eeq_sparam%eps_diis
      astep = eeq_sparam%alpha
      ALLOCATE (xvec(ns, mdiis), fvec(ns, mdiis), bvec(ns))
      xvec = 0.0_dp; fvec = 0.0_dp
      ALLOCATE (vmat(mdiis, mdiis), dmat(mdiis + 1, mdiis + 1), dvec(mdiis + 1))
      dmat = 0.0_dp; dvec = 0.0_dp
      ndiis = 1
      now = 1
      CALL get_energy_gradient(eeqn, eeq_mat, mmat, ewald_env, ewald_pw, &
                               cell, particle_set, xv0, rhs, rv0)
      resin = SQRT(SUM(rv0*rv0))
      !
      DO iv = 1, max_diis
         res = SQRT(SUM(rv0*rv0))
         IF (res > 10._dp*resin) EXIT
         IF (res < eps_diis) EXIT
         !
         now = MOD(iv - 1, mdiis) + 1
         ndiis = MIN(iv, mdiis)
         xvec(1:ns, now) = xv0(1:ns)
         fvec(1:ns, now) = rv0(1:ns)
         DO i = 1, ndiis
            vmat(now, i) = SUM(fvec(:, now)*fvec(:, i))
            vmat(i, now) = vmat(now, i)
         END DO
         IF (ndiis < sdiis) THEN
            xv0(1:ns) = xv0(1:ns) + astep*rv0(1:ns)
         ELSE
            dvec = 0.0_dp
            dvec(ndiis + 1) = 1.0_dp
            dmat(1:ndiis, 1:ndiis) = vmat(1:ndiis, 1:ndiis)
            dmat(ndiis + 1, 1:ndiis) = 1.0_dp
            dmat(1:ndiis, ndiis + 1) = 1.0_dp
            dmat(ndiis + 1, ndiis + 1) = 0.0_dp
            CALL invmat(dmat(1:ndiis + 1, 1:ndiis + 1), info)
            dvec(1:ndiis + 1) = MATMUL(dmat(1:ndiis + 1, 1:ndiis + 1), dvec(1:ndiis + 1))
            xv0(1:ns) = MATMUL(xvec(1:ns, 1:ndiis), dvec(1:ndiis))
            xv0(1:ns) = xv0(1:ns) + MATMUL(fvec(1:ns, 1:ndiis), dvec(1:ndiis))
         END IF
         !
         CALL get_energy_gradient(eeqn, eeq_mat, mmat, ewald_env, ewald_pw, &
                                  cell, particle_set, xv0, rhs, rv0)
      END DO
      charges(1:natom) = xv0(1:natom)
      lambda = xv0(ns)
      eeq_energy = eeqn
      IF (res > eps_diis) ierror = 1
      !
      DEALLOCATE (xvec, fvec, bvec)
      DEALLOCATE (vmat, dmat, dvec)
      DEALLOCATE (xv0, rv0)
      DEALLOCATE (rhs)
      CALL cp_fm_release(pmat)
      CALL cp_fm_release(mmat)

      te = m_walltime()
      IF (iunit > 0) THEN
         IF (ierror == 1) THEN
            WRITE (iunit, '(A)') " EEQ| PBC solver failed to converge "
         ELSE
            WRITE (iunit, '(A,T50,I4,T61,E20.5)') " EEQ| PBC solver: iterations/accuracy ", iv, res
         END IF
         WRITE (iunit, '(A,T67,F14.3)') " EEQ| Iterative PBC solver: time[s]", te - ti
      END IF
      CALL timestop(handle)

   END SUBROUTINE pbc_solver

! **************************************************************************************************
!> \brief ...
!> \param charges ...
!> \param lambda ...
!> \param eeq_energy ...
!> \param eeq_mat ...
!> \param particle_set ...
!> \param kind_of ...
!> \param cell ...
!> \param chia ...
!> \param gam ...
!> \param gab ...
!> \param qtot ...
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param iounit ...
! **************************************************************************************************
   SUBROUTINE fpbc_solver(charges, lambda, eeq_energy, eeq_mat, particle_set, &
                          kind_of, cell, chia, gam, gab, qtot, ewald_env, ewald_pw, iounit)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: charges
      REAL(KIND=dp), INTENT(INOUT)                       :: lambda, eeq_energy
      TYPE(cp_fm_type)                                   :: eeq_mat
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: chia, gam
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: gab
      REAL(KIND=dp), INTENT(IN)                          :: qtot
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      INTEGER, INTENT(IN), OPTIONAL                      :: iounit

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fpbc_solver'

      INTEGER                                            :: ewald_type, handle, ia, iac, iar, ic, &
                                                            ikind, ir, iunit, ix, iy, iz, jkind, &
                                                            natom, ncloc, ncvloc, nkind, nrloc, &
                                                            nrvloc, ns
      INTEGER, DIMENSION(3)                              :: cvec, ncell, periodic
      INTEGER, DIMENSION(:), POINTER                     :: cind, cvind, rind, rvind
      REAL(KIND=dp)                                      :: ad, alpha, deth, dr, grc1, rcut, rmax, &
                                                            te, ti, xr
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rijl, rj
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: pval, xval
      TYPE(cp_fm_struct_type), POINTER                   :: mat_struct, vec_struct
      TYPE(cp_fm_type)                                   :: rhs_vec
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      ti = m_walltime()

      iunit = -1
      IF (PRESENT(iounit)) iunit = iounit

      natom = SIZE(particle_set)
      nkind = SIZE(gam)
      ns = natom + 1
      !
      CALL get_cell(cell=cell, deth=deth)
      CALL ewald_env_get(ewald_env, alpha=alpha, rcut=rcut, ewald_type=ewald_type)
      ad = 2.0_dp*alpha*oorootpi
      IF (ewald_type /= do_ewald_spme) THEN
         CALL cp_abort(__LOCATION__, "Only SPME Ewald method available with EEQ.")
      END IF
      !
      rmax = 2.0_dp*rcut
      ! max cells used
      CALL get_cell(cell, h=hmat, periodic=periodic)
      ncell(1) = CEILING(rmax/plane_distance(1, 0, 0, cell))
      ncell(2) = CEILING(rmax/plane_distance(0, 1, 0, cell))
      ncell(3) = CEILING(rmax/plane_distance(0, 0, 1, cell))
      IF (periodic(1) == 0) ncell(1) = 0
      IF (periodic(2) == 0) ncell(2) = 0
      IF (periodic(3) == 0) ncell(3) = 0
      !
      CALL cp_fm_get_info(eeq_mat, matrix_struct=mat_struct, para_env=para_env)
      CALL cp_fm_set_all(eeq_mat, 0.0_dp, 0.0_dp)
      CALL cp_fm_get_info(eeq_mat, nrow_local=nrloc, ncol_local=ncloc, &
                          row_indices=rind, col_indices=cind)
      CALL cp_fm_struct_create(vec_struct, template_fmstruct=mat_struct, &
                               nrow_global=ns, ncol_global=1)
      CALL cp_fm_create(rhs_vec, vec_struct)
      CALL cp_fm_get_info(rhs_vec, nrow_local=nrvloc, ncol_local=ncvloc, &
                          row_indices=rvind, col_indices=cvind)
      ! response matrix
      DO ir = 1, nrloc
         iar = rind(ir)
         ri = 0.0_dp
         IF (iar <= natom) THEN
            ikind = kind_of(iar)
            ri(1:3) = particle_set(iar)%r(1:3)
         END IF
         DO ic = 1, ncloc
            iac = cind(ic)
            IF (iac > natom .AND. iar > natom) THEN
               eeq_mat%local_data(ir, ic) = 0.0_dp
               CYCLE
            ELSE IF ((iac > natom) .OR. (iar > natom)) THEN
               eeq_mat%local_data(ir, ic) = 1.0_dp
               CYCLE
            END IF
            jkind = kind_of(iac)
            rj(1:3) = particle_set(iac)%r(1:3)
            rij(1:3) = ri(1:3) - rj(1:3)
            rij = pbc(rij, cell)
            DO ix = -ncell(1), ncell(1)
               DO iy = -ncell(2), ncell(2)
                  DO iz = -ncell(3), ncell(3)
                     cvec = [ix, iy, iz]
                     rijl = rij + MATMUL(hmat, cvec)
                     dr = SQRT(SUM(rijl**2))
                     IF (dr > rmax) CYCLE
                     IF (iar == iac .AND. dr < 0.0001_dp) THEN
                        grc1 = gam(ikind) + 2.0_dp*gab(ikind, ikind)*oorootpi - ad
                     ELSE
                        grc1 = erf(gab(ikind, jkind)*dr)/dr - erf(alpha*dr)/dr
                     END IF
                     eeq_mat%local_data(ir, ic) = eeq_mat%local_data(ir, ic) + grc1
                  END DO
               END DO
            END DO
         END DO
      END DO
      !
      ALLOCATE (xval(natom), pval(natom))
      DO ia = 1, natom
         xval = 0.0_dp
         xval(ia) = 1.0_dp
         CALL apply_potential(ewald_env, ewald_pw, cell, particle_set, xval, pval)
         !
         DO ir = 1, nrloc
            iar = rind(ir)
            IF (iar /= ia) CYCLE
            DO ic = 1, ncloc
               iac = cind(ic)
               IF (iac > natom) CYCLE
               eeq_mat%local_data(ir, ic) = eeq_mat%local_data(ir, ic) + pval(iac)
            END DO
         END DO
      END DO
      DEALLOCATE (xval, pval)
      !
      ! set up rhs vector
      DO ir = 1, nrvloc
         iar = rvind(ir)
         DO ic = 1, ncvloc
            iac = cvind(ic)
            ia = MAX(iar, iac)
            IF (ia > natom) THEN
               xr = qtot
            ELSE
               xr = -chia(ia)
            END IF
            rhs_vec%local_data(ir, ic) = xr
         END DO
      END DO
      !
      CALL cp_fm_solve(eeq_mat, rhs_vec)
      !
      charges = 0.0_dp
      lambda = 0.0_dp
      DO ir = 1, nrvloc
         iar = rvind(ir)
         DO ic = 1, ncvloc
            iac = cvind(ic)
            ia = MAX(iar, iac)
            IF (ia <= natom) THEN
               xr = rhs_vec%local_data(ir, ic)
               charges(ia) = xr
            ELSE
               lambda = rhs_vec%local_data(ir, ic)
            END IF
         END DO
      END DO
      CALL para_env%sum(lambda)
      CALL para_env%sum(charges)
      !
      ! energy:   0.5*(q^T.X - lambda*totalcharge)
      eeq_energy = 0.5*SUM(charges(1:natom)*chia(1:natom)) - 0.5_dp*lambda*qtot

      CALL cp_fm_struct_release(vec_struct)
      CALL cp_fm_release(rhs_vec)

      te = m_walltime()
      IF (iunit > 0) THEN
         WRITE (iunit, '(A,T67,F14.3)') " EEQ| Direct PBC solver: time[s]", te - ti
      END IF
      CALL timestop(handle)

   END SUBROUTINE fpbc_solver

! **************************************************************************************************
!> \brief ...
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param cell ...
!> \param particle_set ...
!> \param charges ...
!> \param potential ...
! **************************************************************************************************
   SUBROUTINE apply_potential(ewald_env, ewald_pw, cell, particle_set, charges, potential)
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), POINTER   :: charges
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: potential

      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL ewald_env_get(ewald_env, para_env=para_env)
      potential = 0.0_dp
      CALL spme_potential(ewald_env, ewald_pw, cell, particle_set, charges, &
                          particle_set, potential)
      CALL para_env%sum(potential)

   END SUBROUTINE apply_potential

! **************************************************************************************************
!> \brief ...
!> \param eeqn ...
!> \param fm_mat ...
!> \param mmat ...
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param cell ...
!> \param particle_set ...
!> \param charges ...
!> \param rhs ...
!> \param potential ...
! **************************************************************************************************
   SUBROUTINE get_energy_gradient(eeqn, fm_mat, mmat, ewald_env, ewald_pw, &
                                  cell, particle_set, charges, rhs, potential)
      REAL(KIND=dp), INTENT(INOUT)                       :: eeqn
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat, mmat
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), POINTER   :: charges
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rhs
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: potential

      INTEGER                                            :: na, ns
      REAL(KIND=dp)                                      :: lambda
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mvec
      TYPE(mp_para_env_type), POINTER                    :: para_env

      ns = SIZE(charges)
      na = ns - 1
      CALL ewald_env_get(ewald_env, para_env=para_env)
      potential = 0.0_dp
      CALL spme_potential(ewald_env, ewald_pw, cell, particle_set, charges(1:na), &
                          particle_set, potential(1:na))
      CALL para_env%sum(potential(1:na))
      CALL cp_fm_matvec(fm_mat, charges, potential, alpha=1.0_dp, beta=1.0_dp)
      eeqn = 0.5_dp*SUM(charges(1:na)*potential(1:na)) + SUM(charges(1:na)*rhs(1:na))
      potential(1:ns) = potential(1:ns) + rhs(1:ns)
      ALLOCATE (mvec(ns))
      CALL cp_fm_matvec(mmat, potential, mvec, alpha=-1.0_dp, beta=0.0_dp)
      lambda = -SUM(mvec(1:na))/REAL(na, KIND=dp)
      potential(1:na) = mvec(1:na) + lambda
      DEALLOCATE (mvec)

   END SUBROUTINE get_energy_gradient

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
!> \param ef_energy ...
! **************************************************************************************************
   SUBROUTINE eeq_efield_energy(qs_env, charges, ef_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: charges
      REAL(KIND=dp), INTENT(OUT)                         :: ef_energy

      COMPLEX(KIND=dp)                                   :: zdeta
      COMPLEX(KIND=dp), DIMENSION(3)                     :: zi
      INTEGER                                            :: ia, idir, natom
      LOGICAL                                            :: dfield
      REAL(KIND=dp)                                      :: kr, omega, q
      REAL(KIND=dp), DIMENSION(3)                        :: ci, dfilter, fieldpol, fpolvec, kvec, &
                                                            qi, ria
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL get_qs_env(qs_env, natom=natom, dft_control=dft_control)
      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set)

      IF (dft_control%apply_period_efield) THEN
         dfield = dft_control%period_efield%displacement_field

         fieldpol = dft_control%period_efield%polarisation
         fieldpol = fieldpol/SQRT(DOT_PRODUCT(fieldpol, fieldpol))
         fieldpol = -fieldpol*dft_control%period_efield%strength
         hmat = cell%hmat(:, :)/twopi
         DO idir = 1, 3
            fpolvec(idir) = fieldpol(1)*hmat(1, idir) + fieldpol(2)*hmat(2, idir) &
                            + fieldpol(3)*hmat(3, idir)
         END DO

         zi(:) = CMPLX(1._dp, 0._dp, dp)
         DO ia = 1, natom
            q = charges(ia)
            ria = particle_set(ia)%r
            ria = pbc(ria, cell)
            DO idir = 1, 3
               kvec(:) = twopi*cell%h_inv(idir, :)
               kr = SUM(kvec(:)*ria(:))
               zdeta = CMPLX(COS(kr), SIN(kr), KIND=dp)**q
               zi(idir) = zi(idir)*zdeta
            END DO
         END DO
         qi = AIMAG(LOG(zi))
         IF (dfield) THEN
            dfilter(1:3) = dft_control%period_efield%d_filter(1:3)
            omega = cell%deth
            ci = MATMUL(hmat, qi)/omega
            ef_energy = 0.0_dp
            DO idir = 1, 3
               ef_energy = ef_energy + dfilter(idir)*(fieldpol(idir) - 2._dp*twopi*ci(idir))**2
            END DO
            ef_energy = -0.25_dp*omega/twopi*ef_energy
         ELSE
            ef_energy = SUM(fpolvec(:)*qi(:))
         END IF

      ELSE IF (dft_control%apply_efield) THEN

         fieldpol = dft_control%efield_fields(1)%efield%polarisation* &
                    dft_control%efield_fields(1)%efield%strength

         ef_energy = 0.0_dp
         DO ia = 1, natom
            ria = particle_set(ia)%r
            ria = pbc(ria, cell)
            q = charges(ia)
            ef_energy = ef_energy - q*SUM(fieldpol*ria)
         END DO

      ELSE
         CPABORT("apply field")
      END IF

   END SUBROUTINE eeq_efield_energy

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param efr ...
! **************************************************************************************************
   SUBROUTINE eeq_efield_pot(qs_env, efr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: efr

      INTEGER                                            :: ia, idir, natom
      LOGICAL                                            :: dfield
      REAL(KIND=dp)                                      :: kr
      REAL(KIND=dp), DIMENSION(3)                        :: fieldpol, fpolvec, kvec, ria
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL get_qs_env(qs_env, natom=natom, dft_control=dft_control)
      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set)

      IF (dft_control%apply_period_efield) THEN
         dfield = dft_control%period_efield%displacement_field

         fieldpol = dft_control%period_efield%polarisation
         fieldpol = fieldpol/SQRT(DOT_PRODUCT(fieldpol, fieldpol))
         fieldpol = -fieldpol*dft_control%period_efield%strength
         hmat = cell%hmat(:, :)/twopi
         DO idir = 1, 3
            fpolvec(idir) = fieldpol(1)*hmat(1, idir) + fieldpol(2)*hmat(2, idir) &
                            + fieldpol(3)*hmat(3, idir)
         END DO

         IF (dfield) THEN
            ! dE/dq depends on q, postpone calculation
            efr = 0.0_dp
         ELSE
            efr = 0.0_dp
            DO ia = 1, natom
               ria = particle_set(ia)%r
               ria = pbc(ria, cell)
               DO idir = 1, 3
                  kvec(:) = twopi*cell%h_inv(idir, :)
                  kr = SUM(kvec(:)*ria(:))
                  efr(ia) = efr(ia) + kr*fpolvec(idir)
               END DO
            END DO
         END IF

      ELSE IF (dft_control%apply_efield) THEN

         fieldpol = dft_control%efield_fields(1)%efield%polarisation* &
                    dft_control%efield_fields(1)%efield%strength

         DO ia = 1, natom
            ria = particle_set(ia)%r
            ria = pbc(ria, cell)
            efr(ia) = -SUM(fieldpol*ria)
         END DO

      ELSE
         CPABORT("apply field")
      END IF

   END SUBROUTINE eeq_efield_pot

! **************************************************************************************************
!> \brief ...
!> \param charges ...
!> \param dft_control ...
!> \param particle_set ...
!> \param cell ...
!> \param efr ...
! **************************************************************************************************
   SUBROUTINE eeq_dfield_pot(charges, dft_control, particle_set, cell, efr)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: charges
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: efr

      COMPLEX(KIND=dp)                                   :: zdeta
      COMPLEX(KIND=dp), DIMENSION(3)                     :: zi
      INTEGER                                            :: ia, idir, natom
      REAL(KIND=dp)                                      :: kr, omega, q
      REAL(KIND=dp), DIMENSION(3)                        :: ci, dfilter, fieldpol, kvec, qi, ria
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat

      natom = SIZE(particle_set)

      dfilter(1:3) = dft_control%period_efield%d_filter(1:3)
      fieldpol = dft_control%period_efield%polarisation
      fieldpol = fieldpol/SQRT(DOT_PRODUCT(fieldpol, fieldpol))
      fieldpol = -fieldpol*dft_control%period_efield%strength
      hmat = cell%hmat(:, :)/twopi
      omega = cell%deth
      !
      zi(:) = CMPLX(1._dp, 0._dp, dp)
      DO ia = 1, natom
         q = charges(ia)
         ria = particle_set(ia)%r
         ria = pbc(ria, cell)
         DO idir = 1, 3
            kvec(:) = twopi*cell%h_inv(idir, :)
            kr = SUM(kvec(:)*ria(:))
            zdeta = CMPLX(COS(kr), SIN(kr), KIND=dp)**q
            zi(idir) = zi(idir)*zdeta
         END DO
      END DO
      qi = AIMAG(LOG(zi))
      ci = MATMUL(hmat, qi)/omega
      ci = dfilter*(fieldpol - 2._dp*twopi*ci)
      DO ia = 1, natom
         ria = particle_set(ia)%r
         ria = pbc(ria, cell)
         efr(ia) = efr(ia) - SUM(ci*ria)
      END DO

   END SUBROUTINE eeq_dfield_pot

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
!> \param qlag ...
! **************************************************************************************************
   SUBROUTINE eeq_efield_force_loc(qs_env, charges, qlag)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: charges, qlag

      INTEGER                                            :: atom_a, ia, iatom, ikind, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind
      REAL(KIND=dp)                                      :: q
      REAL(KIND=dp), DIMENSION(3)                        :: fieldpol
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      cell=cell, particle_set=particle_set, &
                      nkind=nkind, natom=natom, &
                      para_env=para_env, &
                      local_particles=local_particles)

      fieldpol = dft_control%efield_fields(1)%efield%polarisation* &
                 dft_control%efield_fields(1)%efield%strength

      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind)
      CALL get_qs_env(qs_env=qs_env, force=force)

      DO ikind = 1, nkind
         force(ikind)%efield = 0.0_dp
         DO ia = 1, local_particles%n_el(ikind)
            iatom = local_particles%list(ikind)%array(ia)
            q = charges(iatom) - qlag(iatom)
            atom_a = atom_of_kind(iatom)
            force(ikind)%efield(1:3, atom_a) = -q*fieldpol(1:3)
         END DO
         CALL para_env%sum(force(ikind)%efield)
      END DO

   END SUBROUTINE eeq_efield_force_loc

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
!> \param qlag ...
! **************************************************************************************************
   SUBROUTINE eeq_efield_force_periodic(qs_env, charges, qlag)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: charges, qlag

      INTEGER                                            :: atom_a, ia, iatom, ikind, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind
      LOGICAL                                            :: dfield, use_virial
      REAL(KIND=dp)                                      :: q
      REAL(KIND=dp), DIMENSION(3)                        :: fa, fieldpol, ria
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pve
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(virial_type), POINTER                         :: virial

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      cell=cell, particle_set=particle_set, &
                      virial=virial, &
                      nkind=nkind, natom=natom, &
                      para_env=para_env, &
                      local_particles=local_particles)

      dfield = dft_control%period_efield%displacement_field
      CPASSERT(.NOT. dfield)

      fieldpol = dft_control%period_efield%polarisation
      fieldpol = fieldpol/SQRT(DOT_PRODUCT(fieldpol, fieldpol))
      fieldpol = -fieldpol*dft_control%period_efield%strength

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind)
      CALL get_qs_env(qs_env=qs_env, force=force)

      pve = 0.0_dp
      DO ikind = 1, nkind
         force(ikind)%efield = 0.0_dp
         DO ia = 1, local_particles%n_el(ikind)
            iatom = local_particles%list(ikind)%array(ia)
            q = charges(iatom) - qlag(iatom)
            fa(1:3) = q*fieldpol(1:3)
            atom_a = atom_of_kind(iatom)
            force(ikind)%efield(1:3, atom_a) = fa
            IF (use_virial) THEN
               ria = particle_set(ia)%r
               ria = pbc(ria, cell)
               CALL virial_pair_force(pve, -0.5_dp, fa, ria)
               CALL virial_pair_force(pve, -0.5_dp, ria, fa)
            END IF
         END DO
         CALL para_env%sum(force(ikind)%efield)
      END DO
      virial%pv_virial = virial%pv_virial + pve

   END SUBROUTINE eeq_efield_force_periodic

END MODULE eeq_method
