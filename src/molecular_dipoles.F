!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Set of routines handling the localization for molecular properties
! **************************************************************************************************
MODULE molecular_dipoles
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_sum
   USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                              molecule_kind_type
   USE molecule_types,                  ONLY: molecule_type
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

   ! *** Public ***
   PUBLIC :: calculate_molecular_dipole

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'molecular_dipoles'

CONTAINS

! **************************************************************************************************
!> \brief maps wfc's to molecules and also prints molecular dipoles
!> \param qs_env the qs_env in which the qs_env lives
!> \param qs_loc_env ...
!> \param loc_print_key ...
!> \param molecule_set ...
! **************************************************************************************************
   SUBROUTINE calculate_molecular_dipole(qs_env, qs_loc_env, loc_print_key, molecule_set)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_loc_env_new_type), INTENT(IN)              :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: loc_print_key
      TYPE(molecule_type), POINTER                       :: molecule_set(:)

      COMPLEX(KIND=dp)                                   :: zeta
      COMPLEX(KIND=dp), DIMENSION(3)                     :: ggamma, zphase
      INTEGER                                            :: akind, first_atom, i, iatom, ikind, &
                                                            imol, imol_now, iounit, ispin, istate, &
                                                            j, natom, nkind, nmol, nspins, nstate, &
                                                            reference
      LOGICAL                                            :: do_berry, floating, ghost
      REAL(KIND=dp)                                      :: charge_tot, dipole(3), ria(3), theta, &
                                                            zeff, zwfc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: charge_set
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dipole_set
      REAL(KIND=dp), DIMENSION(3)                        :: ci, gvec, rcc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: center(:, :)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      logger => cp_get_default_logger()

      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins

      ! Setup reference point
      reference = section_get_ival(loc_print_key, keyword_name="MOLECULAR_DIPOLES%REFERENCE")
      CALL section_vals_val_get(loc_print_key, "MOLECULAR_DIPOLES%REF_POINT", r_vals=ref_point)
      CALL section_vals_val_get(loc_print_key, "MOLECULAR_DIPOLES%PERIODIC", l_val=do_berry)

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, cell=cell)
      particle_set => qs_loc_env%particle_set
      para_env => qs_loc_env%para_env
      local_molecules => qs_loc_env%local_molecules
      nkind = SIZE(local_molecules%n_el)
      zwfc = 3.0_dp - REAL(nspins, KIND=dp)

      ALLOCATE (dipole_set(3, SIZE(molecule_set)))
      ALLOCATE (charge_set(SIZE(molecule_set)))
      dipole_set = 0.0_dp
      charge_set = 0.0_dp

      DO ispin = 1, nspins
         center => qs_loc_env%localized_wfn_control%centers_set(ispin)%array
         nstate = SIZE(center, 2)
         DO ikind = 1, nkind ! loop over different molecules
            nmol = SIZE(local_molecules%list(ikind)%array)
            DO imol = 1, nmol ! all the molecules of the kind
               imol_now = local_molecules%list(ikind)%array(imol) ! index in the global array
               IF (.NOT. ASSOCIATED(molecule_set(imol_now)%lmi(ispin)%states)) CYCLE
               molecule_kind => molecule_set(imol_now)%molecule_kind
               first_atom = molecule_set(imol_now)%first_atom
               CALL get_molecule_kind(molecule_kind=molecule_kind, natom=natom)

               ! Get reference point for this molecule
               CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, &
                                        ref_point=ref_point, ifirst=first_atom, &
                                        ilast=first_atom + natom - 1)

               dipole = 0.0_dp
               IF (do_berry) THEN
                  rcc = pbc(rcc, cell)
                  ! Find out the total charge of the molecule
                  DO iatom = 1, natom
                     i = first_atom + iatom - 1
                     atomic_kind => particle_set(i)%atomic_kind
                     CALL get_atomic_kind(atomic_kind, kind_number=akind)
                     CALL get_qs_kind(qs_kind_set(akind), ghost=ghost, floating=floating)
                     IF (.NOT. ghost .AND. .NOT. floating) THEN
                        CALL get_qs_kind(qs_kind_set(akind), core_charge=zeff)
                        charge_set(imol_now) = charge_set(imol_now) + zeff
                     END IF
                  END DO
                  ! Charges of the wfc involved
                  DO istate = 1, SIZE(molecule_set(imol_now)%lmi(ispin)%states)
                     charge_set(imol_now) = charge_set(imol_now) - zwfc
                  ENDDO

                  charge_tot = charge_set(imol_now)
                  ria = twopi*MATMUL(cell%h_inv, rcc)
                  zphase = CMPLX(COS(ria), SIN(ria), KIND=dp)**charge_tot
                  ggamma = CMPLX(1.0_dp, 0.0_dp, KIND=dp)

                  ! Nuclear charges
                  IF (ispin == 1) THEN
                  DO iatom = 1, natom
                     i = first_atom + iatom - 1
                     atomic_kind => particle_set(i)%atomic_kind
                     CALL get_atomic_kind(atomic_kind, kind_number=akind)
                     CALL get_qs_kind(qs_kind_set(akind), ghost=ghost, floating=floating)
                     IF (.NOT. ghost .AND. .NOT. floating) THEN
                        CALL get_qs_kind(qs_kind_set(akind), core_charge=zeff)
                        ria = pbc(particle_set(i)%r, cell)
                        DO j = 1, 3
                           gvec = twopi*cell%h_inv(j, :)
                           theta = SUM(ria(:)*gvec(:))
                           zeta = CMPLX(COS(theta), SIN(theta), KIND=dp)**(zeff)
                           ggamma(j) = ggamma(j)*zeta
                        END DO
                     END IF
                  END DO
                  END IF

                  ! Charges of the wfc involved
                  DO istate = 1, SIZE(molecule_set(imol_now)%lmi(ispin)%states)
                     i = molecule_set(imol_now)%lmi(ispin)%states(istate)
                     ria = pbc(center(1:3, i), cell)
                     DO j = 1, 3
                        gvec = twopi*cell%h_inv(j, :)
                        theta = SUM(ria(:)*gvec(:))
                        zeta = CMPLX(COS(theta), SIN(theta), KIND=dp)**(-zwfc)
                        ggamma(j) = ggamma(j)*zeta
                     END DO
                  ENDDO

                  ggamma = ggamma*zphase
                  ci = AIMAG(LOG(ggamma))/twopi
                  dipole = MATMUL(cell%hmat, ci)
               ELSE
                  IF (ispin == 1) THEN
                     ! Nuclear charges
                     DO iatom = 1, natom
                        i = first_atom + iatom - 1
                        atomic_kind => particle_set(i)%atomic_kind
                        CALL get_atomic_kind(atomic_kind, kind_number=akind)
                        CALL get_qs_kind(qs_kind_set(akind), ghost=ghost, floating=floating)
                        IF (.NOT. ghost .AND. .NOT. floating) THEN
                           CALL get_qs_kind(qs_kind_set(akind), core_charge=zeff)
                           ria = pbc(particle_set(i)%r, cell) - rcc
                           dipole = dipole + zeff*(ria - rcc)
                           charge_set(imol_now) = charge_set(imol_now) + zeff
                        END IF
                     END DO
                  END IF
                  ! Charges of the wfc involved
                  DO istate = 1, SIZE(molecule_set(imol_now)%lmi(ispin)%states)
                     i = molecule_set(imol_now)%lmi(ispin)%states(istate)
                     ria = pbc(center(1:3, i), cell)
                     dipole = dipole - zwfc*(ria - rcc)
                     charge_set(imol_now) = charge_set(imol_now) - zwfc
                  ENDDO
               END IF
               dipole_set(:, imol_now) = dipole_set(:, imol_now) + dipole ! a.u.
            ENDDO
         ENDDO
      END DO
      CALL mp_sum(dipole_set, para_env%group)
      CALL mp_sum(charge_set, para_env%group)

      iounit = cp_print_key_unit_nr(logger, loc_print_key, "MOLECULAR_DIPOLES", &
                                    extension=".MolDip", middle_name="MOLECULAR_DIPOLES")
      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT='(A80)') &
            "# molecule nr,      charge,           dipole vector,           dipole[Debye]"
         dipole_set(:, :) = dipole_set(:, :)*debye ! Debye
         DO I = 1, SIZE(dipole_set, 2)
            WRITE (UNIT=iounit, FMT='(T8,I6,T21,5F12.6)') I, charge_set(I), dipole_set(1:3, I), &
               SQRT(DOT_PRODUCT(dipole_set(1:3, I), dipole_set(1:3, I)))
         ENDDO
         WRITE (UNIT=iounit, FMT="(T2,A,T61,E20.12)") ' DIPOLE : CheckSum  =', SUM(dipole_set)
      ENDIF
      CALL cp_print_key_finished_output(iounit, logger, loc_print_key, &
                                        "MOLECULAR_DIPOLES")

      DEALLOCATE (dipole_set, charge_set)

   END SUBROUTINE calculate_molecular_dipole
   !------------------------------------------------------------------------------

END MODULE molecular_dipoles

