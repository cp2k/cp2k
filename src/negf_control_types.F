!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input control types for NEGF based quantum transport calculations
! **************************************************************************************************

MODULE negf_control_types
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE input_constants,                 ONLY: negf_run
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi
   USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                              molecule_kind_type
   USE molecule_types,                  ONLY: get_molecule,&
                                              molecule_type
   USE negf_alloc_types,                ONLY: negf_allocatable_ivector
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: kelvin
   USE string_utilities,                ONLY: integer_to_string
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_control_types'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

   PUBLIC :: negf_control_type, negf_control_contact_type
   PUBLIC :: negf_control_create, negf_control_release, read_negf_control

! **************************************************************************************************
!> \brief Input parameters related to a single contact.
!> \author Sergey Chulkov
! **************************************************************************************************
   TYPE negf_control_contact_type
      !> atoms belonging to bulk and screening regions
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atomlist_bulk, atomlist_screening
      !> atom belonging to the primary and secondary bulk unit cells
      TYPE(negf_allocatable_ivector), ALLOCATABLE, &
         DIMENSION(:)                                    :: atomlist_cell
      !> index of the sub_force_env which should be used for bulk calculation
      INTEGER                                            :: force_env_index
      !> contact Fermi level needs to be computed prior NEGF run
      LOGICAL                                            :: compute_fermi_level
      !> when computing contact Fermi level, use the energy given in 'fermi_level' (instead of HOMO)
      !> (instead of the HOMO energy) as a starting point
      LOGICAL                                            :: refine_fermi_level
      !> Fermi level
      REAL(kind=dp)                                      :: fermi_level
      !> temperature [in a.u.]
      REAL(kind=dp)                                      :: temperature
      !> applied electric potential
      REAL(kind=dp)                                      :: v_external
   END TYPE negf_control_contact_type

! **************************************************************************************************
!> \brief Input parameters related to the NEGF run.
!> \author Sergey Chulkov
! **************************************************************************************************
   TYPE negf_control_type
      !> input options for every contact
      TYPE(negf_control_contact_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: contacts
      !> atoms belonging to the scattering region
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atomlist_S
      !> atoms belonging to the scattering region as well as atoms belonging to
      !> screening regions of all the contacts
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atomlist_S_screening
      !> do not keep contact self-energy matrices
      LOGICAL                                            :: disable_cache
      !> convergence criteria for adaptive integration methods
      REAL(kind=dp)                                      :: conv_density
      !> convergence criteria for iterative Lopez-Sancho algorithm
      REAL(kind=dp)                                      :: conv_green
      !> convergence criteria for self-consistent iterations
      REAL(kind=dp)                                      :: conv_scf
      !> accuracy in mapping atoms between different force environments
      REAL(kind=dp)                                      :: eps_geometry
      !> applied bias [in a.u.]
      REAL(kind=dp)                                      :: v_bias
      !> integration lower bound [in a.u.]
      REAL(kind=dp)                                      :: energy_lbound
      !> infinitesimal offset along the imaginary axis [in a.u.]
      REAL(kind=dp)                                      :: eta
      !> initial guess to determine the actual Fermi level of bulk contacts [in a.u.]
      REAL(kind=dp)                                      :: homo_lumo_gap
      !> number of residuals (poles of the Fermi function)
      INTEGER                                            :: delta_npoles
      !> offset along the x-axis away from the poles of the Fermi function [in units of kT]
      INTEGER                                            :: gamma_kT
      !> integration method
      INTEGER                                            :: integr_method
      !> minimal number of grid points along the closed contour
      INTEGER                                            :: integr_min_points
      !> maximal number of grid points along the closed contour
      INTEGER                                            :: integr_max_points
      !> maximal number of SCF iterations
      INTEGER                                            :: max_scf
      !> minimal number of MPI processes to be used to compute Green's function per energy point
      INTEGER                                            :: nprocs
      !> shift in Hartree potential [in a.u.]
      REAL(kind=dp)                                      :: v_shift
      !> initial offset to determine the correct shift in Hartree potential [in a.u.]
      REAL(kind=dp)                                      :: v_shift_offset
      !> maximal number of iteration to determine the shift in Hartree potential
      INTEGER                                            :: v_shift_maxiters
   END TYPE negf_control_type

   PRIVATE :: read_negf_atomlist

CONTAINS

! **************************************************************************************************
!> \brief allocate control options for Non-equilibrium Green's Function calculation
!> \param negf_control an object to create
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE negf_control_create(negf_control)
      TYPE(negf_control_type), POINTER                   :: negf_control

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_control_create'

      INTEGER                                            :: handle

      CPASSERT(.NOT. ASSOCIATED(negf_control))
      CALL timeset(routineN, handle)

      ALLOCATE (negf_control)

      CALL timestop(handle)
   END SUBROUTINE negf_control_create

! **************************************************************************************************
!> \brief release memory allocated for NEGF control options
!> \param negf_control an object to release
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE negf_control_release(negf_control)
      TYPE(negf_control_type), POINTER                   :: negf_control

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_control_release'

      INTEGER                                            :: handle, i, j

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(negf_control)) THEN
         IF (ALLOCATED(negf_control%atomlist_S)) DEALLOCATE (negf_control%atomlist_S)
         IF (ALLOCATED(negf_control%atomlist_S_screening)) DEALLOCATE (negf_control%atomlist_S_screening)

         IF (ALLOCATED(negf_control%contacts)) THEN
            DO i = SIZE(negf_control%contacts), 1, -1
               IF (ALLOCATED(negf_control%contacts(i)%atomlist_bulk)) &
                  DEALLOCATE (negf_control%contacts(i)%atomlist_bulk)

               IF (ALLOCATED(negf_control%contacts(i)%atomlist_screening)) &
                  DEALLOCATE (negf_control%contacts(i)%atomlist_screening)

               IF (ALLOCATED(negf_control%contacts(i)%atomlist_cell)) THEN
                  DO j = SIZE(negf_control%contacts(i)%atomlist_cell), 1, -1
                     IF (ALLOCATED(negf_control%contacts(i)%atomlist_cell(j)%vector)) &
                        DEALLOCATE (negf_control%contacts(i)%atomlist_cell(j)%vector)
                  END DO
                  DEALLOCATE (negf_control%contacts(i)%atomlist_cell)
               END IF
            END DO

            DEALLOCATE (negf_control%contacts)
         END IF

         DEALLOCATE (negf_control)
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_control_release

! **************************************************************************************************
!> \brief Read NEGF input parameters.
!> \param negf_control NEGF control parameters
!> \param input        root input section
!> \param subsys       subsystem environment
! **************************************************************************************************
   SUBROUTINE read_negf_control(negf_control, input, subsys)
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_subsys_type), POINTER                      :: subsys

      CHARACTER(len=*), PARAMETER                        :: routineN = 'read_negf_control'

      CHARACTER(len=default_string_length) :: contact_id_str, eta_current_str, eta_max_str, &
         npoles_current_str, npoles_min_str, temp_current_str, temp_min_str
      INTEGER                                            :: delta_npoles_min, handle, i2_rep, i_rep, &
                                                            n2_rep, n_rep, natoms_current, &
                                                            natoms_total, run_type
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      LOGICAL                                            :: do_negf, is_explicit
      REAL(kind=dp)                                      :: eta_max, temp_current, temp_min
      TYPE(section_vals_type), POINTER                   :: cell_section, contact_section, &
                                                            negf_section, region_section

      CALL timeset(routineN, handle)

      CALL section_vals_val_get(input, "GLOBAL%RUN_TYPE", i_val=run_type)
      do_negf = run_type == negf_run

      negf_section => section_vals_get_subs_vals(input, "NEGF")

      contact_section => section_vals_get_subs_vals(negf_section, "CONTACT")
      CALL section_vals_get(contact_section, n_repetition=n_rep, explicit=is_explicit)
      IF ((.NOT. is_explicit) .AND. do_negf) THEN
         CALL cp_abort(__LOCATION__, &
                       "At least one contact is needed for NEGF calculation.")
      END IF

      ALLOCATE (negf_control%contacts(n_rep))
      DO i_rep = 1, n_rep
         region_section => section_vals_get_subs_vals(contact_section, "SCREENING_REGION", i_rep_section=i_rep)
         CALL section_vals_get(region_section, explicit=is_explicit)

         IF ((.NOT. is_explicit) .AND. do_negf) THEN
            WRITE (contact_id_str, '(I11)') i_rep
            CALL cp_abort(__LOCATION__, &
                          "The screening region must be defined for the contact "//TRIM(ADJUSTL(contact_id_str))//".")
         END IF

         IF (is_explicit) THEN
            CALL read_negf_atomlist(negf_control%contacts(i_rep)%atomlist_screening, region_section, 1, subsys)
         END IF

         region_section => section_vals_get_subs_vals(contact_section, "BULK_REGION", i_rep_section=i_rep)

         CALL section_vals_get(region_section, explicit=is_explicit)

         IF ((.NOT. is_explicit) .AND. do_negf) THEN
            WRITE (contact_id_str, '(I11)') i_rep
            CALL cp_abort(__LOCATION__, &
                          "The bulk region must be defined for the contact "//TRIM(ADJUSTL(contact_id_str))//".")
         END IF

         IF (is_explicit) THEN
            CALL read_negf_atomlist(negf_control%contacts(i_rep)%atomlist_bulk, region_section, 1, subsys)
         END IF

         CALL section_vals_val_get(contact_section, "FORCE_EVAL_SECTION", &
                                   i_val=negf_control%contacts(i_rep)%force_env_index, &
                                   i_rep_section=i_rep)

         cell_section => section_vals_get_subs_vals(region_section, "CELL")
         CALL section_vals_get(cell_section, n_repetition=n2_rep, explicit=is_explicit)

         IF (((.NOT. is_explicit) .OR. n2_rep /= 2) .AND. negf_control%contacts(i_rep)%force_env_index <= 0 .AND. do_negf) THEN
            WRITE (contact_id_str, '(I11)') i_rep
            CALL cp_abort(__LOCATION__, &
                          "You must either provide indices of atoms belonging to two adjacent bulk unit cells "// &
                          "(BULK_REGION/CELL) for the contact, or the index of the FORCE_EVAL section (FORCE_EVAL_SECTION) "// &
                          "which will be used to construct Kohn-Sham matrix for the bulk contact "// &
                          TRIM(ADJUSTL(contact_id_str))//".")
         END IF

         IF (is_explicit .AND. n2_rep > 0) THEN
            ALLOCATE (negf_control%contacts(i_rep)%atomlist_cell(n2_rep))

            DO i2_rep = 1, n2_rep
               CALL read_negf_atomlist(negf_control%contacts(i_rep)%atomlist_cell(i2_rep)%vector, cell_section, i2_rep, subsys)
            END DO
         END IF

         CALL section_vals_val_get(contact_section, "REFINE_FERMI_LEVEL", &
                                   l_val=negf_control%contacts(i_rep)%refine_fermi_level, &
                                   i_rep_section=i_rep)

         CALL section_vals_val_get(contact_section, "FERMI_LEVEL", &
                                   r_val=negf_control%contacts(i_rep)%fermi_level, &
                                   i_rep_section=i_rep, explicit=is_explicit)
         negf_control%contacts(i_rep)%compute_fermi_level = (.NOT. is_explicit) .OR. &
                                                            negf_control%contacts(i_rep)%refine_fermi_level

         IF (do_negf .AND. negf_control%contacts(i_rep)%force_env_index <= 0 .AND. &
             (.NOT. (is_explicit .OR. negf_control%contacts(i_rep)%refine_fermi_level))) THEN
            WRITE (contact_id_str, '(I11)') i_rep
            CALL cp_warn(__LOCATION__, &
                         "There is no way to reasonably guess the Fermi level for the bulk contact "// &
                         TRIM(ADJUSTL(contact_id_str))//" without running a separate bulk DFT calculation first. "// &
                         "Therefore, 0.0 Hartree will be used as an initial guess. It is strongly advised to enable "// &
                         "the REFINE_FERMI_LEVEL switch and to provide an initial guess using the FERMI_LEVEL keyword. "// &
                         "Alternatively, a bulk FORCE_EVAL_SECTION can be set up.")
         END IF

         CALL section_vals_val_get(contact_section, "TEMPERATURE", &
                                   r_val=negf_control%contacts(i_rep)%temperature, &
                                   i_rep_section=i_rep)
         IF (negf_control%contacts(i_rep)%temperature <= 0.0_dp) THEN
            CALL cp_abort(__LOCATION__, "Electronic temperature must be > 0")
         END IF

         CALL section_vals_val_get(contact_section, "ELECTRIC_POTENTIAL", &
                                   r_val=negf_control%contacts(i_rep)%v_external, &
                                   i_rep_section=i_rep)
      END DO

      region_section => section_vals_get_subs_vals(negf_section, "SCATTERING_REGION")
      CALL section_vals_get(region_section, explicit=is_explicit)
      IF (is_explicit) THEN
         CALL read_negf_atomlist(negf_control%atomlist_S, region_section, 1, subsys)
      END IF

      CALL section_vals_val_get(negf_section, "DISABLE_CACHE", l_val=negf_control%disable_cache)

      CALL section_vals_val_get(negf_section, "EPS_DENSITY", r_val=negf_control%conv_density)
      CALL section_vals_val_get(negf_section, "EPS_GREEN", r_val=negf_control%conv_green)
      CALL section_vals_val_get(negf_section, "EPS_SCF", r_val=negf_control%conv_scf)

      CALL section_vals_val_get(negf_section, "EPS_GEO", r_val=negf_control%eps_geometry)

      CALL section_vals_val_get(negf_section, "ENERGY_LBOUND", r_val=negf_control%energy_lbound)
      CALL section_vals_val_get(negf_section, "ETA", r_val=negf_control%eta)
      CALL section_vals_val_get(negf_section, "HOMO_LUMO_GAP", r_val=negf_control%homo_lumo_gap)
      CALL section_vals_val_get(negf_section, "DELTA_NPOLES", i_val=negf_control%delta_npoles)
      CALL section_vals_val_get(negf_section, "GAMMA_KT", i_val=negf_control%gamma_kT)

      CALL section_vals_val_get(negf_section, "INTEGRATION_METHOD", i_val=negf_control%integr_method)
      CALL section_vals_val_get(negf_section, "INTEGRATION_MIN_POINTS", i_val=negf_control%integr_min_points)
      CALL section_vals_val_get(negf_section, "INTEGRATION_MAX_POINTS", i_val=negf_control%integr_max_points)

      IF (negf_control%integr_max_points < negf_control%integr_min_points) &
         negf_control%integr_max_points = negf_control%integr_min_points

      CALL section_vals_val_get(negf_section, "MAX_SCF", i_val=negf_control%max_scf)

      CALL section_vals_val_get(negf_section, "NPROC_POINT", i_val=negf_control%nprocs)

      CALL section_vals_val_get(negf_section, "V_SHIFT", r_val=negf_control%v_shift)
      CALL section_vals_val_get(negf_section, "V_SHIFT_OFFSET", r_val=negf_control%v_shift_offset)
      CALL section_vals_val_get(negf_section, "V_SHIFT_MAX_ITERS", i_val=negf_control%v_shift_maxiters)

      ! check consistency
      IF (negf_control%eta < 0.0_dp) THEN
         CALL cp_abort(__LOCATION__, "ETA must be >= 0")
      END IF

      IF (n_rep > 0) THEN
         delta_npoles_min = NINT(0.5_dp*(negf_control%eta/(pi*MAXVAL(negf_control%contacts(:)%temperature)) + 1.0_dp))
      ELSE
         delta_npoles_min = 1
      END IF

      IF (negf_control%delta_npoles < delta_npoles_min) THEN
         IF (n_rep > 0) THEN
            eta_max = REAL(2*negf_control%delta_npoles - 1, kind=dp)*pi*MAXVAL(negf_control%contacts(:)%temperature)
            temp_current = MAXVAL(negf_control%contacts(:)%temperature)*kelvin
            temp_min = negf_control%eta/(pi*REAL(2*negf_control%delta_npoles - 1, kind=dp))*kelvin

            WRITE (eta_current_str, '(ES11.4E2)') negf_control%eta
            WRITE (eta_max_str, '(ES11.4E2)') eta_max
            WRITE (npoles_current_str, '(I11)') negf_control%delta_npoles
            WRITE (npoles_min_str, '(I11)') delta_npoles_min
            WRITE (temp_current_str, '(F11.3)') temp_current
            WRITE (temp_min_str, '(F11.3)') temp_min

            CALL cp_abort(__LOCATION__, &
                          "Parameter DELTA_NPOLES must be at least "//TRIM(ADJUSTL(npoles_min_str))// &
                          " (instead of "//TRIM(ADJUSTL(npoles_current_str))// &
                          ") for given TEMPERATURE ("//TRIM(ADJUSTL(temp_current_str))// &
                          " K) and ETA ("//TRIM(ADJUSTL(eta_current_str))// &
                          "). Alternatively you can increase TEMPERATURE above "//TRIM(ADJUSTL(temp_min_str))// &
                          " K, or decrease ETA below "//TRIM(ADJUSTL(eta_max_str))// &
                          ". Please keep in mind that very tight ETA may result in dramatical precision loss"// &
                          " due to inversion of ill-conditioned matrices.")
         ELSE
            ! no leads have been defined, so calculation will abort anyway
            negf_control%delta_npoles = delta_npoles_min
         END IF
      END IF

      ! expand scattering region by adding atoms from contact screening regions
      n_rep = SIZE(negf_control%contacts)
      IF (ALLOCATED(negf_control%atomlist_S)) THEN
         natoms_total = SIZE(negf_control%atomlist_S)
      ELSE
         natoms_total = 0
      END IF

      DO i_rep = 1, n_rep
         IF (ALLOCATED(negf_control%contacts(i_rep)%atomlist_screening)) THEN
            IF (ALLOCATED(negf_control%contacts(i_rep)%atomlist_screening)) &
               natoms_total = natoms_total + SIZE(negf_control%contacts(i_rep)%atomlist_screening)
         END IF
      END DO

      IF (natoms_total > 0) THEN
         ALLOCATE (negf_control%atomlist_S_screening(natoms_total))
         IF (ALLOCATED(negf_control%atomlist_S)) THEN
            natoms_total = SIZE(negf_control%atomlist_S)
            negf_control%atomlist_S_screening(1:natoms_total) = negf_control%atomlist_S(1:natoms_total)
         ELSE
            natoms_total = 0
         END IF

         DO i_rep = 1, n_rep
            IF (ALLOCATED(negf_control%contacts(i_rep)%atomlist_screening)) THEN
               natoms_current = SIZE(negf_control%contacts(i_rep)%atomlist_screening)

               negf_control%atomlist_S_screening(natoms_total + 1:natoms_total + natoms_current) = &
                  negf_control%contacts(i_rep)%atomlist_screening(1:natoms_current)

               natoms_total = natoms_total + natoms_current
            END IF
         END DO

         ! sort and remove duplicated atoms
         ALLOCATE (inds(natoms_total))
         CALL sort(negf_control%atomlist_S_screening, natoms_total, inds)
         DEALLOCATE (inds)

         natoms_current = 1
         DO i_rep = natoms_current + 1, natoms_total
            IF (negf_control%atomlist_S_screening(i_rep) /= negf_control%atomlist_S_screening(natoms_current)) THEN
               natoms_current = natoms_current + 1
               negf_control%atomlist_S_screening(natoms_current) = negf_control%atomlist_S_screening(i_rep)
            END IF
         END DO

         IF (natoms_current < natoms_total) THEN
            CALL MOVE_ALLOC(negf_control%atomlist_S_screening, inds)

            ALLOCATE (negf_control%atomlist_S_screening(natoms_current))
            negf_control%atomlist_S_screening(1:natoms_current) = inds(1:natoms_current)
            DEALLOCATE (inds)
         END IF
      END IF

      IF (do_negf .AND. SIZE(negf_control%contacts) > 2) THEN
         CALL cp_abort(__LOCATION__, &
                       "General case (> 2 contacts) has not been implemented yet")
      END IF

      CALL timestop(handle)
   END SUBROUTINE read_negf_control

! **************************************************************************************************
!> \brief Read region-specific list of atoms.
!> \param atomlist        list of atoms
!> \param input_section   input section which contains 'LIST' and 'MOLNAME' keywords
!> \param i_rep_section   repetition index of the input_section
!> \param subsys          subsystem environment
! **************************************************************************************************
   SUBROUTINE read_negf_atomlist(atomlist, input_section, i_rep_section, subsys)
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(out)    :: atomlist
      TYPE(section_vals_type), POINTER                   :: input_section
      INTEGER, INTENT(in)                                :: i_rep_section
      TYPE(cp_subsys_type), POINTER                      :: subsys

      CHARACTER(len=*), PARAMETER :: routineN = 'read_negf_atomlist'

      CHARACTER(len=default_string_length)               :: index_str, natoms_str
      CHARACTER(len=default_string_length), &
         DIMENSION(:), POINTER                           :: cptr
      INTEGER :: first_atom, handle, iatom, ikind, imol, iname, irep, last_atom, natoms_current, &
         natoms_max, natoms_total, nkinds, nmols, nnames, nrep_list, nrep_molname
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      INTEGER, DIMENSION(:), POINTER                     :: iptr
      LOGICAL                                            :: is_list, is_molname
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(molecule_type), POINTER                       :: molecule
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      CALL cp_subsys_get(subsys, particle_set=particle_set, &
                         molecule_set=molecule_set, &
                         molecule_kind_set=molecule_kind_set)
      natoms_max = SIZE(particle_set)
      nkinds = SIZE(molecule_kind_set)

      CALL section_vals_val_get(input_section, "LIST", i_rep_section=i_rep_section, &
                                n_rep_val=nrep_list, explicit=is_list)
      CALL section_vals_val_get(input_section, "MOLNAME", i_rep_section=i_rep_section, &
                                n_rep_val=nrep_molname, explicit=is_molname)

      ! compute the number of atoms in the NEGF region, and check the validity of giben atomic indices
      natoms_total = 0
      IF (is_list .AND. nrep_list > 0) THEN
         DO irep = 1, nrep_list
            CALL section_vals_val_get(input_section, "LIST", i_rep_section=i_rep_section, i_rep_val=irep, i_vals=iptr)

            natoms_current = SIZE(iptr)
            DO iatom = 1, natoms_current
               IF (iptr(iatom) > natoms_max) THEN
                  CALL integer_to_string(iptr(iatom), index_str)
                  CALL integer_to_string(natoms_max, natoms_str)
                  CALL cp_abort(__LOCATION__, &
                                "NEGF: Atomic index "//TRIM(index_str)//" given in section "// &
                                TRIM(input_section%section%name)//" exceeds the maximum number of atoms ("// &
                                TRIM(natoms_str)//").")
               END IF
            END DO

            natoms_total = natoms_total + natoms_current
         END DO
      END IF

      IF (is_molname .AND. nrep_molname > 0) THEN
         DO irep = 1, nrep_molname
            CALL section_vals_val_get(input_section, "MOLNAME", i_rep_section=i_rep_section, i_rep_val=irep, c_vals=cptr)
            nnames = SIZE(cptr)

            DO iname = 1, nnames
               DO ikind = 1, nkinds
                  IF (molecule_kind_set(ikind)%name .EQ. cptr(iname)) EXIT
               END DO

               IF (ikind <= nkinds) THEN
                  molecule_kind => molecule_kind_set(ikind)
                  CALL get_molecule_kind(molecule_kind, nmolecule=nmols, molecule_list=iptr)

                  DO imol = 1, nmols
                     molecule => molecule_set(iptr(imol))
                     CALL get_molecule(molecule, first_atom=first_atom, last_atom=last_atom)
                     natoms_current = last_atom - first_atom + 1
                     natoms_total = natoms_total + natoms_current
                  END DO
               ELSE
                  CALL cp_abort(__LOCATION__, &
                                "NEGF: A molecule with the name '"//TRIM(cptr(iname))//"' mentioned in section "// &
                                TRIM(input_section%section%name)//" has not been defined. Note that names are case sensitive.")
               END IF
            END DO
         END DO
      END IF

      ! create a list of atomic indices
      IF (natoms_total > 0) THEN
         ALLOCATE (atomlist(natoms_total))

         natoms_total = 0

         IF (is_list .AND. nrep_list > 0) THEN
            DO irep = 1, nrep_list
               CALL section_vals_val_get(input_section, "LIST", i_rep_section=i_rep_section, i_rep_val=irep, i_vals=iptr)

               natoms_current = SIZE(iptr)
               atomlist(natoms_total + 1:natoms_total + natoms_current) = iptr(1:natoms_current)
               natoms_total = natoms_total + natoms_current
            END DO
         END IF

         IF (is_molname .AND. nrep_molname > 0) THEN
            DO irep = 1, nrep_molname
               CALL section_vals_val_get(input_section, "MOLNAME", i_rep_section=i_rep_section, i_rep_val=irep, c_vals=cptr)
               nnames = SIZE(cptr)

               DO iname = 1, nnames
                  DO ikind = 1, nkinds
                     IF (molecule_kind_set(ikind)%name .EQ. cptr(iname)) EXIT
                  END DO

                  IF (ikind <= nkinds) THEN
                     molecule_kind => molecule_kind_set(ikind)
                     CALL get_molecule_kind(molecule_kind, nmolecule=nmols, molecule_list=iptr)

                     DO imol = 1, nmols
                        molecule => molecule_set(iptr(imol))
                        CALL get_molecule(molecule, first_atom=first_atom, last_atom=last_atom)

                        DO natoms_current = first_atom, last_atom
                           natoms_total = natoms_total + 1
                           atomlist(natoms_total) = natoms_current
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END IF

         ! remove duplicated atoms
         ALLOCATE (inds(natoms_total))
         CALL sort(atomlist, natoms_total, inds)
         DEALLOCATE (inds)

         natoms_current = 1
         DO iatom = natoms_current + 1, natoms_total
            IF (atomlist(iatom) /= atomlist(natoms_current)) THEN
               natoms_current = natoms_current + 1
               atomlist(natoms_current) = atomlist(iatom)
            END IF
         END DO

         IF (natoms_current < natoms_total) THEN
            CALL MOVE_ALLOC(atomlist, inds)

            ALLOCATE (atomlist(natoms_current))
            atomlist(1:natoms_current) = inds(1:natoms_current)
            DEALLOCATE (inds)
         END IF
      END IF

      CALL timestop(handle)
   END SUBROUTINE read_negf_atomlist
END MODULE negf_control_types
