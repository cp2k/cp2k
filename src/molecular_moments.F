!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Set of routines handling the localization for molecular properties
! **************************************************************************************************
MODULE molecular_moments
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_schur_product
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_constants,                 ONLY: use_mom_ref_com
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_maxloc,&
                                              mp_sum
   USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                              molecule_kind_type
   USE molecule_types,                  ONLY: molecule_type
   USE moments_utils,                   ONLY: get_reference_point
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_loc_types,                    ONLY: qs_loc_env_new_type
   USE qs_moments,                      ONLY: build_local_moment_matrix
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public ***
   PUBLIC :: calculate_molecular_moments

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'molecular_moments'

CONTAINS

! **************************************************************************************************
!> \brief Calculates electrical molecular moments using local operators r-r_ref
!>        r_ref: center of mass of the molecule
!>        Output is in atomic units
!> \param qs_env the qs_env in which the qs_env lives
!> \param qs_loc_env ...
!> \param mo_local ...
!> \param loc_print_key ...
!> \param molecule_set ...
! **************************************************************************************************
   SUBROUTINE calculate_molecular_moments(qs_env, qs_loc_env, mo_local, loc_print_key, molecule_set)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_loc_env_new_type), INTENT(IN)              :: qs_loc_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mo_local
      TYPE(section_vals_type), POINTER                   :: loc_print_key
      TYPE(molecule_type), POINTER                       :: molecule_set(:)

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_molecular_moments', &
         routineP = moduleN//':'//routineN

      INTEGER :: akind, first_atom, i, iatom, ikind, imol, imol_now, iounit, iproc, ispin, j, lx, &
         ly, lz, molkind, n, n1, n2, natom, ncol_global, nm, nmol, nmols, norder, nrow_global, ns, &
         nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: states
      INTEGER, DIMENSION(2)                              :: nstates
      LOGICAL                                            :: floating, ghost
      REAL(KIND=dp)                                      :: zeff, zwfc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: charge_set
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: moment_set
      REAL(KIND=dp), DIMENSION(3)                        :: rcc, ria
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_localized, momv, mvector, omvector
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, moments
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      logger => cp_get_default_logger()

      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins
      zwfc = 3.0_dp-REAL(nspins, KIND=dp)

      CALL section_vals_val_get(loc_print_key, "MOLECULAR_MOMENTS%ORDER", i_val=norder)
      CPASSERT(norder >= 0)
      nm = ncoset(norder)-1

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, cell=cell)
      particle_set => qs_loc_env%particle_set
      para_env => qs_loc_env%para_env
      local_molecules => qs_loc_env%local_molecules
      molkind = SIZE(local_molecules%n_el)
      nmols = SIZE(molecule_set)
      ALLOCATE (charge_set(nmols), moment_set(nm, nmols))
      charge_set = 0.0_dp
      moment_set = 0.0_dp

      IF (norder > 0) THEN
         CALL get_qs_env(qs_env, matrix_s=matrix_s)
         DO imol = 1, SIZE(molecule_set)
            molecule_kind => molecule_set(imol)%molecule_kind
            first_atom = molecule_set(imol)%first_atom
            CALL get_molecule_kind(molecule_kind=molecule_kind, natom=natom)
            ! Get reference point for this molecule
            CALL get_reference_point(rcc, qs_env=qs_env, reference=use_mom_ref_com, &
                                     ref_point=ref_point, ifirst=first_atom, &
                                     ilast=first_atom+natom-1)
            ALLOCATE (moments(nm))
            DO i = 1, nm
               ALLOCATE (moments(i)%matrix)
               CALL dbcsr_copy(moments(i)%matrix, matrix_s(1)%matrix, 'MOM MAT')
               CALL dbcsr_set(moments(i)%matrix, 0.0_dp)
            END DO
            !
            CALL build_local_moment_matrix(qs_env, moments, norder, rcc)
            !
            DO ispin = 1, nspins
               IF (ASSOCIATED(molecule_set(imol)%lmi)) THEN
                  nstates(1) = molecule_set(imol)%lmi(ispin)%nstates
               ELSE
                  nstates(1) = 0
               ENDIF
               nstates(2) = para_env%mepos
               CALL mp_maxloc(nstates, para_env%group)
               IF (nstates(1) == 0) CYCLE
               ns = nstates(1)
               iproc = nstates(2)
               ALLOCATE (states(ns))
               IF (iproc == para_env%mepos) THEN
                  states(:) = molecule_set(imol)%lmi(ispin)%states(:)
               ELSE
                  states(:) = 0
               END IF
               CALL mp_bcast(states, iproc, para_env%group)
               ! assemble local states for this molecule
               mo_localized => mo_local(ispin)%matrix
               CALL cp_fm_get_info(mo_localized, ncol_global=ncol_global, nrow_global=nrow_global)
               CALL cp_fm_struct_create(fm_struct, nrow_global=nrow_global, ncol_global=ns, &
                                        para_env=mo_localized%matrix_struct%para_env, &
                                        context=mo_localized%matrix_struct%context)
               CALL cp_fm_create(mvector, fm_struct, name="mvector")
               CALL cp_fm_create(omvector, fm_struct, name="omvector")
               CALL cp_fm_create(momv, fm_struct, name="omvector")
               CALL cp_fm_struct_release(fm_struct)
               !
               DO i = 1, ns
                  CALL cp_fm_to_fm(mo_localized, mvector, 1, states(i), i)
               END DO
               DO i = 1, nm
                  CALL cp_dbcsr_sm_fm_multiply(moments(i)%matrix, mvector, omvector, ns)
                  CALL cp_fm_schur_product(mvector, omvector, momv)
                  moment_set(i, imol) = moment_set(i, imol)-zwfc*SUM(momv%local_data)
               END DO
               !
               CALL cp_fm_release(mvector)
               CALL cp_fm_release(omvector)
               CALL cp_fm_release(momv)
               DEALLOCATE (states)
            END DO
            DO i = 1, nm
               CALL dbcsr_deallocate_matrix(moments(i)%matrix)
            END DO
            DEALLOCATE (moments)
         END DO
      END IF
      !
      DO ikind = 1, molkind ! loop over different molecules
         nmol = SIZE(local_molecules%list(ikind)%array)
         DO imol = 1, nmol ! all the molecules of the kind
            imol_now = local_molecules%list(ikind)%array(imol) ! index in the global array
            molecule_kind => molecule_set(imol_now)%molecule_kind
            first_atom = molecule_set(imol_now)%first_atom
            CALL get_molecule_kind(molecule_kind=molecule_kind, natom=natom)
            ! Get reference point for this molecule
            CALL get_reference_point(rcc, qs_env=qs_env, reference=use_mom_ref_com, &
                                     ref_point=ref_point, ifirst=first_atom, &
                                     ilast=first_atom+natom-1)
            ! charge
            DO iatom = 1, natom
               i = first_atom+iatom-1
               atomic_kind => particle_set(i)%atomic_kind
               CALL get_atomic_kind(atomic_kind, kind_number=akind)
               CALL get_qs_kind(qs_kind_set(akind), ghost=ghost, floating=floating)
               IF (.NOT. ghost .AND. .NOT. floating) THEN
                  CALL get_qs_kind(qs_kind_set(akind), core_charge=zeff)
                  charge_set(imol_now) = charge_set(imol_now)+zeff
               END IF
            END DO
            DO ispin = 1, nspins
               IF (ASSOCIATED(molecule_set(imol_now)%lmi(ispin)%states)) THEN
                  ns = SIZE(molecule_set(imol_now)%lmi(ispin)%states)
                  charge_set(imol_now) = charge_set(imol_now)-zwfc*ns
               END IF
            ENDDO
            !
            IF (norder > 0) THEN
               ! nuclear contribution
               DO i = 1, nm
                  lx = indco(1, i+1)
                  ly = indco(2, i+1)
                  lz = indco(3, i+1)
                  DO iatom = 1, natom
                     j = first_atom+iatom-1
                     atomic_kind => particle_set(j)%atomic_kind
                     CALL get_atomic_kind(atomic_kind, kind_number=akind)
                     CALL get_qs_kind(qs_kind_set(akind), ghost=ghost, floating=floating)
                     IF (.NOT. ghost .AND. .NOT. floating) THEN
                        CALL get_qs_kind(qs_kind_set(akind), core_charge=zeff)
                        ria = particle_set(j)%r-rcc
                        ria = pbc(ria, cell)
                        moment_set(i, imol_now) = moment_set(i, imol_now)+ &
                                                  zeff*ria(1)**lx*ria(2)**ly*ria(3)**lz
                     END IF
                  END DO
               END DO
            END IF
         END DO
      END DO
      CALL mp_sum(moment_set, para_env%group)
      CALL mp_sum(charge_set, para_env%group)

      iounit = cp_print_key_unit_nr(logger, loc_print_key, "MOLECULAR_MOMENTS", &
                                    extension=".MolMom", middle_name="MOLECULAR_MOMENTS")
      IF (iounit > 0) THEN
         DO i = 1, SIZE(charge_set)
            WRITE (UNIT=iounit, FMT='(A,I6,A,F12.6)') "  # molecule nr:", i, "      Charge:", charge_set(I)
            DO n = 1, norder
               n1 = ncoset(n-1)
               n2 = ncoset(n)-1
               WRITE (UNIT=iounit, FMT='(T4,A,I2,10(T16,6F12.6))') "Order:", n, moment_set(n1:n2, i)
            END DO
         ENDDO
      ENDIF
      CALL cp_print_key_finished_output(iounit, logger, loc_print_key, &
                                        "MOLECULAR_MOMENTS")

      DEALLOCATE (charge_set, moment_set)

   END SUBROUTINE calculate_molecular_moments
   !------------------------------------------------------------------------------

END MODULE molecular_moments

