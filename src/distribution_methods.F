!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Distribution methods for atoms, particles, or molecules
!> \par History
!>      - 1d-distribution of molecules and particles (Sep. 2003, MK)
!>      - 2d-distribution for Quickstep updated with molecules (Oct. 2003, MK)
!> \author MK (22.08.2003)
! *****************************************************************************
MODULE distribution_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE cp_array_utils_i,                ONLY: cp_1d_i_p_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_dbcsr_interface,              ONLY: cp_distribution_get_num_images,&
                                             heap_fill,&
                                             heap_get_first,&
                                             heap_new,&
                                             heap_release,&
                                             heap_reset_first,&
                                             heap_t
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE distribution_1d_types,           ONLY: distribution_1d_create,&
                                             distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_create,&
                                             distribution_2d_type,&
                                             distribution_2d_write
  USE input_constants,                 ONLY: model_block_count,&
                                             model_block_lmax
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: pi
  USE mathlib,                         ONLY: gcd,&
                                             lcm
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             get_molecule_kind_set,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type
  USE particle_types,                  ONLY: particle_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'distribution_methods'

! *** Public subroutines ***

  PUBLIC :: distribute_molecules_1d, &
            distribute_molecules_2d

CONTAINS

! *****************************************************************************
!> \brief Distribute molecules and particles
!> \param atomic_kind_set particle (atomic) kind information
!> \param particle_set particle information
!> \param local_particles distribution of particles created by this routine
!> \param molecule_kind_set molecule kind information
!> \param molecule_set molecule information
!> \param local_molecules distribution of molecules created by this routine
!> \param force_env_section ...
!> \param prev_molecule_kind_set previous molecule kind information, used with
!>        prev_local_molecules
!> \param prev_local_molecules previous distribution of molecules, new one will 
!>        be identical if all the prev_* arguments are present and associated
!> \par History
!>      none
!> \author MK (Jun. 2003)
! *****************************************************************************
  SUBROUTINE distribute_molecules_1d(atomic_kind_set, particle_set, &
                                     local_particles, &
                                     molecule_kind_set, molecule_set, &
                                     local_molecules, force_env_section, &
                                     prev_molecule_kind_set, &
                                     prev_local_molecules)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(molecule_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: prev_molecule_kind_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: prev_local_molecules

    CHARACTER(len=*), PARAMETER :: routineN = 'distribute_molecules_1d', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, bin, bin_price, group, handle, iatom, imolecule, &
      imolecule_kind, imolecule_local, imolecule_prev_kind, iparticle_kind, &
      ipe, istat, iw, kind_a, molecule_a, mype, n, natom, nbins, nload, &
      nmolecule, nmolecule_kind, nparticle_kind, npe, nsgf, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nmolecule_local, &
                                                nparticle_local, work, &
                                                workload_count, workload_fill
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    LOGICAL                                  :: found, has_prev_subsys_info, &
                                                heap_error, is_local
    TYPE(cp_1d_i_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: local_molecule
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(heap_t)                             :: bin_heap_count, bin_heap_fill
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

! integer n_local_molecules
! integer, allocatable :: all_local_molecules(:)

    CALL timeset(routineN, handle)


    has_prev_subsys_info = .FALSE.
    IF (PRESENT(prev_local_molecules) .AND.  &
        PRESENT(prev_molecule_kind_set)) THEN
       IF (ASSOCIATED(prev_local_molecules) .AND. &
           ASSOCIATED(prev_molecule_kind_set)) THEN
         has_prev_subsys_info = .TRUE.
       ENDIF
    ENDIF

    logger => cp_get_default_logger()

    group = logger%para_env%group
    mype = logger%para_env%mepos + 1
    npe = logger%para_env%num_pe

    ALLOCATE (workload_count(npe), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "workload", npe*int_size)
    workload_count(:) = 0

    ALLOCATE (workload_fill(npe), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "workload_old", npe*int_size)
    workload_fill(:) = 0

    nmolecule_kind = SIZE(molecule_kind_set)

    ALLOCATE (nmolecule_local(nmolecule_kind), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "nmolecule_local", nmolecule_kind*int_size)
    nmolecule_local(:) = 0

    ALLOCATE (local_molecule(nmolecule_kind), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "local_molecule", nmolecule_kind*int_size)

    nparticle_kind = SIZE(atomic_kind_set)

    ALLOCATE (nparticle_local(nparticle_kind), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "nparticle_local", nparticle_kind*int_size)
    nparticle_local(:) = 0

    nbins=npe

    CALL heap_new (bin_heap_count, nbins, heap_error)
    CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error, &
         routineN, "Error creating heap.")
    CALL heap_fill (bin_heap_count, &
         (/(bin, bin=1, nbins)/), workload_count , heap_error)
    CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error, &
         routineN, "Error initially filling the heap.")

    CALL heap_new (bin_heap_fill, nbins, heap_error)
    CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error, &
         routineN, "Error creating heap.")
    CALL heap_fill (bin_heap_fill, &
         (/(bin, bin=1, nbins)/), workload_fill , heap_error)
    CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error, &
         routineN, "Error initially filling the heap.")

    DO imolecule_kind=1, nmolecule_kind

      molecule_kind => molecule_kind_set(imolecule_kind)

      NULLIFY (molecule_list)

!     *** Get the number of molecules and the number of ***
!     *** atoms in each molecule of that molecular kind ***

      CALL get_molecule_kind(molecule_kind=molecule_kind, &
                             molecule_list=molecule_list, &
                             natom=natom, &
                             nsgf=nsgf)

!     *** Consider the number of atoms or basis ***
!     *** functions which depends on the method ***

      nload = MAX(natom, nsgf)
      nmolecule = SIZE(molecule_list)

!     *** Get the number of local molecules of the current molecule kind ***

      DO imolecule=1, nmolecule
         IF (has_prev_subsys_info) THEN
            DO imolecule_prev_kind=1, SIZE(prev_molecule_kind_set)
               IF (ANY(prev_local_molecules%list(imolecule_prev_kind)%array(1:prev_local_molecules%n_el(imolecule_prev_kind)) == &
                        molecule_list(imolecule))) THEN
                  ! molecule used to be local
                  nmolecule_local(imolecule_kind) = nmolecule_local(imolecule_kind) + 1
               ENDIF
            END DO
         ELSE
            CALL heap_get_first (bin_heap_count, bin, bin_price, found, heap_error)
            CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error, &
                 routineN, "Error getting topmost heap element.")
            CALL cp_assert (found, cp_fatal_level, cp_internal_error, &
                 routineN, "No topmost heap element found.")

            ipe = bin
            CALL cp_assert (bin_price==workload_count(ipe), cp_fatal_level, cp_internal_error, &
                 routineN, "inconsistent heap")

            workload_count(ipe) = workload_count(ipe) + nload
            IF (ipe == mype) THEN
              nmolecule_local(imolecule_kind) = nmolecule_local(imolecule_kind) + 1
            END IF

            bin_price = workload_count(ipe)
            CALL heap_reset_first (bin_heap_count, bin_price, heap_error)
            CALL cp_assert (.NOT. heap_error, cp_warning_level, cp_internal_error, &
                 routineN, "Error setting price of top heap element.")
         END IF
      END DO


!     *** Distribute the molecules ***
      n = nmolecule_local(imolecule_kind)

      IF (n > 0) THEN
        ALLOCATE (local_molecule(imolecule_kind)%array(n), STAT=istat)
        IF (istat /= 0) THEN
           CALL stop_memory(routineN, moduleN, __LINE__, &
                            "local_molecule(imolecule_kind)%array", &
                            n*int_size)
        END IF
      ELSE
        NULLIFY (local_molecule(imolecule_kind)%array)
      END IF

      imolecule_local = 0
      DO imolecule=1, nmolecule
         is_local = .FALSE.
         IF (has_prev_subsys_info) THEN
            DO imolecule_prev_kind=1, SIZE(prev_molecule_kind_set)
               IF (ANY(prev_local_molecules%list(imolecule_prev_kind)%array(1:prev_local_molecules%n_el(imolecule_prev_kind)) == &
                        molecule_list(imolecule))) THEN
                  is_local = .TRUE.
               END IF
            END DO
         ELSE
            CALL heap_get_first (bin_heap_fill, bin, bin_price, found, heap_error)
            CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error, &
                 routineN, "Error getting topmost heap element.")
            CALL cp_assert (found, cp_fatal_level, cp_internal_error, &
                 routineN, "No topmost heap element found.")

            ipe = bin
            CALL cp_assert (bin_price==workload_fill(ipe), cp_fatal_level, cp_internal_error, &
                 routineN, "inconsistent heap")

            workload_fill(ipe) = workload_fill(ipe) + nload
            is_local = (ipe == mype)
         ENDIF
         IF (is_local) THEN
            imolecule_local = imolecule_local + 1
            molecule_a = molecule_list(imolecule)
            local_molecule(imolecule_kind)%array(imolecule_local) = molecule_a
            DO iatom=1, natom
              atom_a = molecule_set(molecule_a)%first_atom + iatom - 1

              CALL get_atomic_kind(atomic_kind=particle_set(atom_a)%atomic_kind, &
                                   kind_number=kind_a)
              nparticle_local(kind_a) = nparticle_local(kind_a) + 1
            END DO
         END IF
         IF (.NOT. has_prev_subsys_info) THEN
            bin_price = workload_fill(ipe)
            CALL heap_reset_first (bin_heap_fill, bin_price, heap_error)
            CALL cp_assert (.NOT. heap_error, cp_warning_level, cp_internal_error, &
                 routineN, "Error setting price of top heap element.")
         END IF
      END DO

    END DO

    CALL cp_assert (ALL(workload_fill.EQ.workload_count), cp_fatal_level, cp_internal_error, &
              routineN, "Inconsistent heaps encountered")

    CALL heap_release (bin_heap_count, heap_error)
    CALL cp_assert (.NOT. heap_error, cp_warning_level, cp_internal_error, &
         routineN, "Error releasing heap.")
    CALL heap_release (bin_heap_fill, heap_error)
    CALL cp_assert (.NOT. heap_error, cp_warning_level, cp_internal_error, &
         routineN, "Error releasing heap.")

!   *** Create the local molecule structure ***

    CALL distribution_1d_create(local_molecules, &
                                n_el=nmolecule_local, &
                                para_env=logger%para_env)

!   *** Create the local particle structure ***

    CALL distribution_1d_create(local_particles, &
                                n_el=nparticle_local, &
                                para_env=logger%para_env)

!   *** Store the generated local molecule and particle distributions ***

    nparticle_local(:) = 0

    DO imolecule_kind=1, nmolecule_kind

      IF (nmolecule_local(imolecule_kind) == 0) CYCLE

      local_molecules%list(imolecule_kind)%array(:) =&
        local_molecule(imolecule_kind)%array(:)

      molecule_kind => molecule_kind_set(imolecule_kind)

      CALL get_molecule_kind(molecule_kind=molecule_kind, &
                             natom=natom)

      DO imolecule=1, nmolecule_local(imolecule_kind)
        molecule_a = local_molecule(imolecule_kind)%array(imolecule)
        DO iatom=1, natom
          atom_a = molecule_set(molecule_a)%first_atom + iatom - 1
          CALL get_atomic_kind(atomic_kind=particle_set(atom_a)%atomic_kind, &
                               kind_number=kind_a)
          nparticle_local(kind_a) = nparticle_local(kind_a) + 1
          local_particles%list(kind_a)%array(nparticle_local(kind_a)) = atom_a
        END DO
      END DO

    END DO

!   *** Print distribution, if requested ***

    IF (BTEST(cp_print_key_should_output(logger%iter_info, &
         force_env_section, "PRINT%DISTRIBUTION1D"), cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger, force_env_section, "PRINT%DISTRIBUTION1D", &
            extension=".Log")

       iw = output_unit
       IF (output_unit<0) iw = cp_logger_get_default_unit_nr(logger, LOCAL=.TRUE.)

!     *** Print molecule distribution ***

      ALLOCATE (work(npe), STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                       "work", npe*int_size)
      work(:) = 0

      work(mype) = SUM(nmolecule_local)
      CALL mp_sum(work, group)

      IF (output_unit>0) THEN
         WRITE (UNIT=output_unit, &
              FMT="(/, T2, A, T51, A, /, (T52, I6, T73, I8))")&
              "DISTRIBUTION OF THE MOLECULES", &
              "Process    Number of molecules", &
              (ipe-1, work(ipe), ipe=1, npe)
         WRITE (UNIT=output_unit, FMT="(T55, A3, T73, I8)")&
              "Sum", SUM(work)
         CALL m_flush(output_unit)
      END IF

      CALL mp_sync(group)

      DO ipe=1, npe
         IF (ipe == mype) THEN
            WRITE (UNIT=iw, FMT="(/, T3, A)")&
                 "Process   Kind   Local molecules (global indices)"
            DO imolecule_kind=1, nmolecule_kind
               IF (imolecule_kind == 1) THEN
                  WRITE (UNIT=iw, FMT="(T4, I6, 2X, I5, (T21, 10I6))")&
                       ipe-1, imolecule_kind, &
                       (local_molecules%list(imolecule_kind)%array(imolecule), &
                       imolecule=1, nmolecule_local(imolecule_kind))
               ELSE
                  WRITE (UNIT=iw, FMT="(T12, I5, (T21, 10I6))")&
                       imolecule_kind, &
                       (local_molecules%list(imolecule_kind)%array(imolecule), &
                       imolecule=1, nmolecule_local(imolecule_kind))
               END IF
            END DO
         END IF
         CALL m_flush(iw)
         CALL mp_sync(group)
      END DO

!     *** Print particle distribution ***

      work(:) = 0

      work(mype) = SUM(nparticle_local)
      CALL mp_sum(work, group)

      IF (output_unit>0) THEN
         WRITE (UNIT=output_unit, &
              FMT="(/, T2, A, T51, A, /, (T52, I6, T73, I8))")&
              "DISTRIBUTION OF THE PARTICLES", &
              "Process    Number of particles", &
              (ipe-1, work(ipe), ipe=1, npe)
         WRITE (UNIT=output_unit, FMT="(T55, A3, T73, I8)")&
              "Sum", SUM(work)
         CALL m_flush(output_unit)
      END IF

      CALL mp_sync(group)

      DO ipe=1, npe
         IF (ipe == mype) THEN
            WRITE (UNIT=iw, FMT="(/, T3, A)")&
                 "Process   Kind   Local particles (global indices)"
            DO iparticle_kind=1, nparticle_kind
               IF (iparticle_kind == 1) THEN
                  WRITE (UNIT=iw, FMT="(T4, I6, 2X, I5, (T20, 10I6))")&
                       ipe-1, iparticle_kind, &
                       (local_particles%list(iparticle_kind)%array(iatom), &
                       iatom=1, nparticle_local(iparticle_kind))
               ELSE
                  WRITE (UNIT=iw, FMT="(T12, I5, (T20, 10I6))")&
                       iparticle_kind, &
                       (local_particles%list(iparticle_kind)%array(iatom), &
                       iatom=1, nparticle_local(iparticle_kind))
               END IF
            END DO
         END IF
         CALL m_flush(iw)
         CALL mp_sync(group)
      END DO
      DEALLOCATE (work, STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, "work")

      CALL cp_print_key_finished_output(output_unit, logger, force_env_section, &
           "PRINT%DISTRIBUTION1D")
   END IF
!   *** Release work storage ***

    DEALLOCATE (workload_count, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "workload")

    DEALLOCATE (workload_fill, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "workload_old")

    DEALLOCATE (nmolecule_local, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "nmolecule_local")

    DEALLOCATE (nparticle_local, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "nparticle_local")

    DO imolecule_kind=1, nmolecule_kind
      IF (ASSOCIATED(local_molecule(imolecule_kind)%array)) THEN
        DEALLOCATE (local_molecule(imolecule_kind)%array, STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN, moduleN, __LINE__, &
                           "local_molecule(imolecule_kind)%array")
        END IF
      END IF
    END DO
    DEALLOCATE (local_molecule, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "local_molecule")

    CALL timestop(handle)

  END SUBROUTINE distribute_molecules_1d

! *****************************************************************************
!> \brief Distributes the particle pairs creating a 2d distribution optimally
!>      suited for quickstep
!> \param cell ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param distribution_2d the distribution that will be created by this
!>                         method
!> \param blacs_env the parallel environement at the basis of the
!>                   distribution
!> \param force_env_section ...
!> \par History
!>      - local_rows & cols blocksize optimizations (Aug. 2003, MK)
!>      - cleanup of distribution_2d (Sep. 2003, fawzi)
!>      - update for molecules (Oct. 2003, MK)
!> \author fawzi (Feb. 2003)
!> \note
!>      Intermediate generation of a 2d distribution of the molecules, but
!>      only the corresponding particle (atomic) distribution is currently
!>      used. The 2d distribution of the molecules is deleted, but may easily
!>      be recovered (MK).
! *****************************************************************************
  SUBROUTINE distribute_molecules_2d(cell, atomic_kind_set, particle_set, &
                                     qs_kind_set, molecule_kind_set, molecule_set, &
                                     distribution_2d, blacs_env, force_env_section)
    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(section_vals_type), POINTER         :: force_env_section

    CHARACTER(len=*), PARAMETER :: routineN = 'distribute_molecules_2d', &
      routineP = moduleN//':'//routineN

    INTEGER :: cluster_price, cost_model, group, handle, iatom, iatom_mol, &
      iatom_one, ikind, imol, imolecule, imolecule_kind, iparticle_kind, &
      ipcol, iprow, istat, iw, kind_a, mypcol, myprow, n, natom, natom_mol, &
      nclusters, nmolecule, nmolecule_kind, nparticle_kind, npcol, nprow, &
      nsgf, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cluster_list, cluster_prices, &
                                                nparticle_local_col, &
                                                nparticle_local_row, work
    INTEGER, DIMENSION(:), POINTER           :: lmax_basis, molecule_list
    INTEGER, DIMENSION(:, :), POINTER :: cluster_col_distribution, &
      cluster_row_distribution, col_distribution, row_distribution
    LOGICAL :: basic_cluster_optimization, basic_optimization, &
      basic_spatial_optimization, failure, molecular_distribution, &
      skip_optimization
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: coords, pbc_scaled_coords
    REAL(KIND=dp), DIMENSION(3)              :: center
    TYPE(cp_1d_i_p_type), DIMENSION(:), &
      POINTER                                :: local_particle_col, &
                                                local_particle_row
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(section_vals_type), POINTER         :: distribution_section

!...

    failure = .FALSE.

    CALL timeset(routineN, handle)

    logger => cp_get_default_logger()

    distribution_section => section_vals_get_subs_vals(force_env_section, "DFT%QS%DISTRIBUTION")

    CALL section_vals_val_get(distribution_section, "2D_MOLECULAR_DISTRIBUTION", l_val=molecular_distribution)
    CALL section_vals_val_get(distribution_section, "SKIP_OPTIMIZATION", l_val=skip_optimization)
    CALL section_vals_val_get(distribution_section, "BASIC_OPTIMIZATION", l_val=basic_optimization)
    CALL section_vals_val_get(distribution_section, "BASIC_SPATIAL_OPTIMIZATION", l_val=basic_spatial_optimization)
    CALL section_vals_val_get(distribution_section, "BASIC_CLUSTER_OPTIMIZATION", l_val=basic_cluster_optimization)

    CALL section_vals_val_get(distribution_section, "COST_MODEL", i_val=cost_model)
    !

    group = blacs_env%para_env%group
    myprow = blacs_env%mepos(1) + 1
    mypcol = blacs_env%mepos(2) + 1
    nprow = blacs_env%num_pe(1)
    npcol = blacs_env%num_pe(2)

    nmolecule_kind = SIZE(molecule_kind_set)
    CALL get_molecule_kind_set(molecule_kind_set, nmolecule=nmolecule)

    nparticle_kind = SIZE(atomic_kind_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)

    !
    ! we need to generate two representations of the distribution, one as a straight array with global particles
    ! one ordered wrt to kinds and only listing the local particles
    !
    ALLOCATE (row_distribution(natom,2), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "row_distribution", natom*int_size)

    ALLOCATE (col_distribution(natom,2), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "col_distribution", natom*int_size)
    ! Initialize the distributions to -1, as the second dimension only gets set with cluster optimization
    ! but the information is needed by dbcsr
    row_distribution=-1; col_distribution=-1

    ALLOCATE (local_particle_col(nparticle_kind), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "local_particle_col", &
                                     nparticle_kind*int_size)

    ALLOCATE (local_particle_row(nparticle_kind), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "local_particle_col", &
                                     nparticle_kind*int_size)

    ALLOCATE (nparticle_local_row(nparticle_kind), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "nparticle_local_row", &
                                     nparticle_kind*int_size)

    ALLOCATE (nparticle_local_col(nparticle_kind), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "nparticle_local_col", &
                                     nparticle_kind*int_size)

    IF (basic_optimization .OR. basic_spatial_optimization .OR. basic_cluster_optimization) THEN

       IF (molecular_distribution) THEN
           nclusters=nmolecule
       ELSE
           nclusters=natom
       ENDIF

       ALLOCATE (cluster_list (nclusters), stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
            "cluster_list", nclusters*int_size)
       ALLOCATE (cluster_prices (nclusters), stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
            "cluster_prices", nclusters*int_size)
       ALLOCATE (cluster_row_distribution (nclusters,2), stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
            "cluster_prices", nclusters*int_size)
       ALLOCATE (cluster_col_distribution (nclusters,2), stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
            "cluster_prices", nclusters*int_size)
       cluster_row_distribution=-1; cluster_col_distribution=-1


       ! Fill in the clusters and their prices
       CALL section_vals_val_get(distribution_section, "COST_MODEL", i_val=cost_model)
       IF (.NOT.molecular_distribution) THEN
          DO iatom = 1, natom
             CALL cp_assert (iatom .LE. nclusters, cp_fatal_level, &
                  cp_internal_error, routineN, "Bounds error")
             CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
             cluster_list(iatom) = iatom
             SELECT CASE (cost_model)
             CASE (model_block_count)
                CALL get_qs_kind(qs_kind_set(ikind), nsgf=nsgf)
                cluster_price = nsgf
             CASE (model_block_lmax)
                CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
                CALL get_gto_basis_set(orb_basis_set, lmax=lmax_basis)
                cluster_price = MAXVAL (lmax_basis)
             CASE default
                CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
                CALL get_gto_basis_set(orb_basis_set, lmax=lmax_basis)
                cluster_price = 8 + (MAXVAL (lmax_basis)**2)
             END SELECT
             cluster_prices(iatom) = cluster_price
          ENDDO
       ELSE
          imol = 0
          DO imolecule_kind=1, nmolecule_kind
             molecule_kind => molecule_kind_set(imolecule_kind)
             CALL get_molecule_kind(molecule_kind=molecule_kind, molecule_list=molecule_list, natom=natom_mol)
             DO imolecule=1, SIZE(molecule_list)
                imol = imol + 1
                cluster_list(imol) = imol
                cluster_price = 0
                DO iatom_mol=1, natom_mol
                   iatom = molecule_set(molecule_list(imolecule))%first_atom + iatom_mol - 1
                   CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
                   SELECT CASE (cost_model)
                   CASE (model_block_count)
                      CALL get_qs_kind(qs_kind_set(ikind), nsgf=nsgf)
                      cluster_price = cluster_price + nsgf
                   CASE (model_block_lmax)
                      CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
                      CALL get_gto_basis_set(orb_basis_set, lmax=lmax_basis)
                      cluster_price = cluster_price + MAXVAL (lmax_basis)
                   CASE default
                      CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
                      CALL get_gto_basis_set(orb_basis_set, lmax=lmax_basis)
                      cluster_price = cluster_price +  8 + (MAXVAL (lmax_basis)**2)
                   END SELECT
                ENDDO
                cluster_prices(imol) = cluster_price
             ENDDO
          ENDDO
       ENDIF

       ! And distribute
       IF (basic_optimization) THEN
          CALL make_basic_distribution (cluster_list, cluster_prices, &
               nprow, cluster_row_distribution(:,1), npcol, cluster_col_distribution(:,1))
       ELSE
          IF(basic_cluster_optimization)THEN
          CALL cp_assert(.NOT. molecular_distribution,cp_failure_level,cp_assertion_failed,&
                         routineP,"clustering and molecular blocking NYI",failure)
             ALLOCATE(pbc_scaled_coords(3, natom),coords(3, natom))
             DO iatom=1, natom
                CALL real_to_scaled(pbc_scaled_coords(:, iatom), pbc(particle_set(iatom)%r(:), cell), cell)
                coords(:,iatom)=pbc(particle_set(iatom)%r(:),cell)
             ENDDO
             CALL make_cluster_distribution(coords, pbc_scaled_coords, cell, cluster_prices,&
                    nprow, cluster_row_distribution, npcol, cluster_col_distribution)
          ELSE! basic_spatial_optimization
             ALLOCATE(pbc_scaled_coords(3, nclusters))
             IF (.NOT. molecular_distribution) THEN
                ! just scaled coords
                DO iatom=1, natom
                   CALL real_to_scaled(pbc_scaled_coords(:, iatom), pbc(particle_set(iatom)%r(:), cell), cell)
                ENDDO
             ELSE
                ! use scaled coords of geometric center, folding when appropriate
                imol = 0
                DO imolecule_kind=1, nmolecule_kind
                   molecule_kind => molecule_kind_set(imolecule_kind)
                   CALL get_molecule_kind(molecule_kind=molecule_kind, molecule_list=molecule_list, natom=natom_mol)
                   DO imolecule=1, SIZE(molecule_list)
                      imol = imol + 1
                      iatom_one = molecule_set(molecule_list(imolecule))%first_atom
                      center = 0.0_dp
                      DO iatom_mol=1, natom_mol
                         iatom = molecule_set(molecule_list(imolecule))%first_atom + iatom_mol - 1
                         center = center + &
                           pbc(particle_set(iatom)%r(:)-particle_set(iatom_one)%r(:), cell)+particle_set(iatom_one)%r(:)
                      ENDDO
                      center = center / natom_mol
                      CALL real_to_scaled(pbc_scaled_coords(:, imol), pbc(center, cell), cell)
                   ENDDO
                ENDDO
             ENDIF
   
             CALL make_basic_spatial_distribution (pbc_scaled_coords, cluster_prices, &
                            nprow, cluster_row_distribution(:,1), npcol, cluster_col_distribution(:,1))
   
             DEALLOCATE(pbc_scaled_coords)
          END IF
       ENDIF

       ! And assign back
       IF (.NOT. molecular_distribution) THEN
          row_distribution = cluster_row_distribution
          col_distribution = cluster_col_distribution
       ELSE
          imol = 0
          DO imolecule_kind=1, nmolecule_kind
             molecule_kind => molecule_kind_set(imolecule_kind)
             CALL get_molecule_kind(molecule_kind=molecule_kind, molecule_list=molecule_list, natom=natom_mol)
             DO imolecule=1, SIZE(molecule_list)
                imol = imol + 1
                DO iatom_mol=1, natom_mol
                   iatom = molecule_set(molecule_list(imolecule))%first_atom + iatom_mol - 1
                   row_distribution(iatom,:) = cluster_row_distribution(imol,:)
                   col_distribution(iatom,:) = cluster_col_distribution(imol,:)
                ENDDO
             ENDDO
          ENDDO
       ENDIF

       ! cleanup
       DEALLOCATE (cluster_list, stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
            "cluster_list", nclusters*int_size)
       DEALLOCATE (cluster_prices, stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
            "cluster_prices", nclusters*int_size)
       DEALLOCATE (cluster_row_distribution, stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
            "cluster_prices", nclusters*int_size)
       DEALLOCATE (cluster_col_distribution, stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
            "cluster_prices", nclusters*int_size)

    ELSE
       ! expects nothing else
       CPASSERT(.FALSE.)
    ENDIF


    ! prepare the lists of local particles

    ! count local particles of a given kind
    nparticle_local_col=0
    nparticle_local_row=0
    DO iatom=1, natom
       CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, kind_number=kind_a)
       IF (row_distribution(iatom,1)==myprow) nparticle_local_row(kind_a)=nparticle_local_row(kind_a)+1
       IF (col_distribution(iatom,1)==mypcol) nparticle_local_col(kind_a)=nparticle_local_col(kind_a)+1
    ENDDO

    ! allocate space
    DO iparticle_kind=1, nparticle_kind
      n = nparticle_local_row(iparticle_kind)
      ALLOCATE (local_particle_row(iparticle_kind)%array(n), STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                       "local_particle_row(iparticle_kind)%array", n*int_size)

      n = nparticle_local_col(iparticle_kind)
      ALLOCATE (local_particle_col(iparticle_kind)%array(n), STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                       "local_particle_col(iparticle_kind)%array", n*int_size)
    ENDDO

    ! store
    nparticle_local_col=0
    nparticle_local_row=0
    DO iatom=1, natom
       CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, kind_number=kind_a)
       IF (row_distribution(iatom,1)==myprow) THEN
           nparticle_local_row(kind_a)=nparticle_local_row(kind_a)+1
           local_particle_row(kind_a)%array(nparticle_local_row(kind_a)) = iatom
       ENDIF
       IF (col_distribution(iatom,1)==mypcol) THEN
           nparticle_local_col(kind_a)=nparticle_local_col(kind_a)+1
           local_particle_col(kind_a)%array(nparticle_local_col(kind_a)) = iatom
       ENDIF
    ENDDO

!   *** Generate the 2d distribution structure  but take care of the zero offsets required
    row_distribution(:,1) = row_distribution(:,1) - 1
    col_distribution(:,1) = col_distribution(:,1) - 1
    CALL distribution_2d_create(distribution_2d, &
                                row_distribution_ptr=row_distribution, &
                                col_distribution_ptr=col_distribution, &
                                local_rows_ptr=local_particle_row, &
                                local_cols_ptr=local_particle_col, &
                                blacs_env=blacs_env)

    NULLIFY (local_particle_row)
    NULLIFY (local_particle_col)
    NULLIFY (row_distribution)
    NULLIFY (col_distribution)

!   *** Print distribution, if requested ***
    IF (BTEST(cp_print_key_should_output(logger%iter_info, &
         force_env_section, "PRINT%DISTRIBUTION"), cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger, force_env_section, "PRINT%DISTRIBUTION", &
            extension=".Log")

!     *** Print row distribution ***

      ALLOCATE (work(nprow), STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN, moduleN, __LINE__, &
                         "work", nprow*int_size)
      END IF
      work(:) = 0

      IF (mypcol == 1) work(myprow) = SUM(distribution_2d%n_local_rows)

      CALL mp_sum(work, group)

      IF (output_unit>0) THEN
         WRITE (UNIT=output_unit, &
              FMT="(/, T2, A, /, T15, A, /, (T16, I10, T41, I10, T71, I10))")&
              "DISTRIBUTION OF THE PARTICLES (ROWS)", &
              "Process row      Number of particles         Number of matrix rows", &
              (iprow-1, work(iprow), -1, iprow=1, nprow)
         WRITE (UNIT=output_unit, FMT="(T23, A3, T41, I10, T71, I10)")&
              "Sum", SUM(work), -1
         CALL m_flush(output_unit)
      END IF

      DEALLOCATE (work, STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN, moduleN, __LINE__, "work")
      END IF

!     *** Print column distribution ***

      ALLOCATE (work(npcol), STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN, moduleN, __LINE__, &
                         "work", npcol*int_size)
      END IF
      work(:) = 0

      IF (myprow == 1) work(mypcol) = SUM(distribution_2d%n_local_cols)

      CALL mp_sum(work, group)

      IF (output_unit>0) THEN
         WRITE (UNIT=output_unit, &
              FMT="(/, T2, A, /, T15, A, /, (T16, I10, T41, I10, T71, I10))")&
              "DISTRIBUTION OF THE PARTICLES (COLUMNS)", &
              "Process col      Number of particles      Number of matrix columns", &
              (ipcol-1, work(ipcol), -1, ipcol=1, npcol)
         WRITE (UNIT=output_unit, FMT="(T23, A3, T41, I10, T71, I10)")&
              "Sum", SUM(work), -1
         CALL m_flush(output_unit)
      END IF

      DEALLOCATE (work, STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN, moduleN, __LINE__, "work")
      END IF

      CALL cp_print_key_finished_output(output_unit, logger, force_env_section, &
           "PRINT%DISTRIBUTION")
    END IF

    IF (BTEST(cp_print_key_should_output(logger%iter_info, &
         force_env_section, "PRINT%DISTRIBUTION2D"), cp_p_file)) THEN

         iw = cp_logger_get_default_unit_nr(logger, LOCAL=.TRUE.)
         CALL distribution_2d_write(distribution_2d, &
                                    unit_nr=iw, &
                                    local=.TRUE., &
                                    long_description=.TRUE.)

    ENDIF

!   *** Release work storage ***

    DEALLOCATE (nparticle_local_row, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "nparticle_local_row")

    DEALLOCATE (nparticle_local_col, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "nparticle_local_col")

    CALL timestop(handle)

  END SUBROUTINE distribute_molecules_2d

! *****************************************************************************
!> \brief Creates a basic distribution
!> \param cluster_list ...
!> \param cluster_prices ...
!> \param nprows ...
!> \param row_distribution ...
!> \param npcols ...
!> \param col_distribution ...
!> \par History
!> - Created 2010-08-06 UB
! *****************************************************************************
  SUBROUTINE make_basic_distribution (cluster_list, cluster_prices, &
       nprows, row_distribution, npcols, col_distribution)
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: cluster_list, cluster_prices
    INTEGER, INTENT(IN)                      :: nprows
    INTEGER, DIMENSION(:), INTENT(OUT)       :: row_distribution
    INTEGER, INTENT(IN)                      :: npcols
    INTEGER, DIMENSION(:), INTENT(OUT)       :: col_distribution

    CHARACTER(len=*), PARAMETER :: routineN = 'make_basic_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER :: bin, bin_price, cluster, cluster_index, cluster_price, nbins, &
      nclusters, pcol, pgrid_gcd, prow, timing_handle
    LOGICAL                                  :: failure, found, heap_error
    TYPE(heap_t)                             :: bin_heap

!   ---------------------------------------------------------------------------

    CALL timeset (routineN, timing_handle)
    nbins = lcm (nprows, npcols)
    pgrid_gcd = gcd (nprows, npcols)
    CALL sort (cluster_prices, SIZE(cluster_list), cluster_list)
    CALL heap_new (bin_heap, nbins, heap_error)
    CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error, &
         routineN, "Error creating heap.")
    CALL heap_fill (bin_heap, &
         (/(bin, bin=0, nbins-1)/), (/(0, bin=1, nbins)/), heap_error)
    CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error, &
         routineN, "Error initially filling the heap.")
    !
    nclusters = SIZE (cluster_list)
    ! Put the most expensive cluster in the bin with the smallest
    ! price and repeat.
    DO cluster_index = nclusters, 1, -1
       cluster = cluster_list(cluster_index)
       CALL heap_get_first (bin_heap, bin, bin_price, found, heap_error)
       CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error, &
            routineN, "Error getting topmost heap element.")
       CALL cp_assert (found, cp_fatal_level, cp_internal_error, &
            routineN, "No topmost heap element found.")
       !
       failure = .FALSE.
       prow = INT (bin*pgrid_gcd / npcols)
       CALL cp_assert (prow .LT. nprows, cp_fatal_level, cp_internal_error, &
            routineN, "Invalid process row.", failure=failure)
       failure = .TRUE.
       pcol = INT (bin*pgrid_gcd / nprows)
       CALL cp_assert (pcol .LT. npcols, cp_fatal_level, cp_internal_error, &
            routineN, "Invalid process column.", failure=failure)
       row_distribution (cluster) = prow + 1
       col_distribution (cluster) = pcol + 1
       !
       cluster_price = cluster_prices(cluster_index)
       bin_price = bin_price + cluster_price
       CALL heap_reset_first (bin_heap, bin_price, heap_error)
       CALL cp_assert (.NOT. heap_error, cp_warning_level, cp_internal_error, &
            routineN, "Error setting price of top heap element.")
    ENDDO
    CALL heap_release (bin_heap, heap_error)
    CALL cp_assert (.NOT. heap_error, cp_warning_level, cp_internal_error, &
         routineN, "Error releasing heap.")
    CALL timestop (timing_handle)
  END SUBROUTINE make_basic_distribution

! *****************************************************************************
!> \brief Creates a basic spatial distribution
!>        that tries to make the corresponding blocks as homogeneous as possible
!> \param pbc_scaled_coords ...
!> \param costs ...
!> \param nprows ...
!> \param row_distribution ...
!> \param npcols ...
!> \param col_distribution ...
!> \par History
!> - Created 2010-11-11 Joost VandeVondele
! *****************************************************************************
  SUBROUTINE make_basic_spatial_distribution (pbc_scaled_coords, costs, &
       nprows, row_distribution, npcols, col_distribution)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: pbc_scaled_coords
    INTEGER, DIMENSION(:), INTENT(IN)        :: costs
    INTEGER, INTENT(IN)                      :: nprows
    INTEGER, DIMENSION(:), INTENT(OUT)       :: row_distribution
    INTEGER, INTENT(IN)                      :: npcols
    INTEGER, DIMENSION(:), INTENT(OUT)       :: col_distribution

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'make_basic_spatial_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, natoms, nbins, &
                                                pgrid_gcd
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bin_costs, distribution

    CALL timeset (routineN, handle)

    natoms=SIZE(costs)
    nbins = lcm (nprows, npcols)
    pgrid_gcd = gcd (nprows, npcols)
    ALLOCATE(bin_costs(nbins), distribution(natoms))
    bin_costs=0

    CALL spatial_recurse(pbc_scaled_coords, costs, (/(iatom, iatom=1, natoms)/) , bin_costs, distribution, 0)

    ! WRITE(*, *) "Final bin costs: ", bin_costs

    ! final row_distribution / col_distribution
    DO iatom=1, natoms
       row_distribution (iatom) = (distribution(iatom)-1)*pgrid_gcd / npcols +1
       col_distribution (iatom) = (distribution(iatom)-1)*pgrid_gcd / nprows +1
    ENDDO

    DEALLOCATE(bin_costs, distribution)

    CALL timestop(handle)

  END SUBROUTINE make_basic_spatial_distribution

! *****************************************************************************
!> \brief ...
!> \param pbc_scaled_coords ...
!> \param costs ...
!> \param indices ...
!> \param bin_costs ...
!> \param distribution ...
!> \param level ...
! *****************************************************************************
  RECURSIVE SUBROUTINE spatial_recurse(pbc_scaled_coords, costs, indices, bin_costs, distribution, level)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: pbc_scaled_coords
    INTEGER, DIMENSION(:), INTENT(IN)        :: costs, indices
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: bin_costs, distribution
    INTEGER, INTENT(IN)                      :: level

    INTEGER                                  :: iatom, ibin, natoms, nbins, &
                                                nhalf
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_costs_sorted, &
                                                atom_permutation, &
                                                bin_costs_sorted, permutation
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: coord

    natoms=SIZE(costs)
    nbins=SIZE(bin_costs)
    nhalf=(natoms+1)/2

    IF (natoms<=nbins) THEN
       ! assign the most expensive atom to the least costly bin
       ALLOCATE(bin_costs_sorted(nbins), permutation(nbins))
       bin_costs_sorted(:)=bin_costs
       CALL sort(bin_costs_sorted, nbins, permutation)
       ALLOCATE(atom_costs_sorted(natoms), atom_permutation(natoms))
       atom_costs_sorted(:)=costs
       CALL sort(atom_costs_sorted, natoms, atom_permutation)
       ibin=0
       ! WRITE(*, *) "Dealing with a new bunch of atoms "
       DO iatom=natoms, 1, -1
          ibin=ibin+1
          ! WRITE(*, *) "atom", indices(atom_permutation(iatom)), "cost", atom_costs_sorted(iatom), &
          !            "bin", permutation(ibin), "its cost", bin_costs(permutation(ibin))
          ! WRITE(100, '(A, I0, 3F12.6)') "A", permutation(ibin), pbc_scaled_coords(:, atom_permutation(iatom))
          bin_costs(permutation(ibin))=bin_costs(permutation(ibin))+atom_costs_sorted(iatom)
          distribution(indices(atom_permutation(iatom)))=permutation(ibin)
       ENDDO
       DEALLOCATE(bin_costs_sorted, permutation, atom_costs_sorted, atom_permutation)
    ELSE
       ! divide atoms in two subsets, sorting according to their coordinates, alternatively x, y, z
       ! recursively do this for both subsets
       ALLOCATE(coord(natoms), permutation(natoms))
       coord(:)=pbc_scaled_coords(MOD(level, 3)+1, :)
       CALL sort(coord, natoms, permutation)
       CALL spatial_recurse(pbc_scaled_coords(:, permutation(1:nhalf)), costs(permutation(1:nhalf)), &
                            indices(permutation(1:nhalf)), bin_costs, distribution, level+1)
       CALL spatial_recurse(pbc_scaled_coords(:, permutation(nhalf+1:)), costs(permutation(nhalf+1:)), &
                            indices(permutation(nhalf+1:)), bin_costs, distribution, level+1)
       DEALLOCATE(coord, permutation)
    ENDIF

  END SUBROUTINE spatial_recurse


! *****************************************************************************
!> \brief creates a distribution placing close by atoms into clusters and 
!>        putting them on the smae processors. Load balancing is 
!>        performed by balancing sum of the cluster costs per processor
!> \param coords coordinates of the system
!> \param scaled_coords scaled coordinates
!> \param cell the cell_type
!> \param costs costs per atomic block
!> \param nprows number of precessors per row on the 2d grid
!> \param row_distribution the resulting distribution over proc_rows of atomic blocks
!> \param npcols number of precessors per col on the 2d grid
!> \param col_distribution the resulting distribution over proc_cols of atomic blocks
! *****************************************************************************
  SUBROUTINE make_cluster_distribution(coords, scaled_coords, cell, costs, &
       nprows, row_distribution, npcols, col_distribution)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: coords, scaled_coords
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, DIMENSION(:), INTENT(IN)        :: costs
    INTEGER, INTENT(IN)                      :: nprows
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: row_distribution
    INTEGER, INTENT(IN)                      :: npcols
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: col_distribution

    CHARACTER(len=*), PARAMETER :: routineN = 'make_cluster_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, icluster, level, &
                                                natom, ncluster
    INTEGER, ALLOCATABLE, DIMENSION(:) :: atom_to_cluster, cluster_cost, &
      cluster_to_col, cluster_to_row, piv_cost, proc_cost, sorted_cost

    CALL timeset (routineN, handle)
    natom=SIZE(costs)
    ncluster=cp_distribution_get_num_images(SUM(costs),natom,nprows,npcols)
    ALLOCATE(atom_to_cluster(natom))
    ALLOCATE(cluster_cost(ncluster))
    ALLOCATE(cluster_to_row(ncluster))
    ALLOCATE(cluster_to_col(ncluster))
    ALLOCATE(sorted_cost(ncluster))
    ALLOCATE(piv_cost(ncluster))

    icluster=0
    CALL cluster_recurse(coords, scaled_coords, cell, costs, atom_to_cluster, ncluster, icluster, cluster_cost)

    sorted_cost(:)=cluster_cost(:)
    CALL sort(sorted_cost,ncluster,piv_cost)

    ALLOCATE(proc_cost(nprows))
    proc_cost=0; level=1
    CALL assign_clusters(cluster_cost,piv_cost,proc_cost,cluster_to_row,nprows)

    DEALLOCATE(proc_cost); ALLOCATE(proc_cost(npcols))
    proc_cost=0; level=1
    CALL assign_clusters(cluster_cost,piv_cost,proc_cost,cluster_to_col,npcols)

    DO i=1,natom
       row_distribution(i,1)=cluster_to_row(atom_to_cluster(i))
       row_distribution(i,2)=atom_to_cluster(i)
       col_distribution(i,1)=cluster_to_col(atom_to_cluster(i))
       col_distribution(i,2)=atom_to_cluster(i)
    END DO

    DEALLOCATE(atom_to_cluster,cluster_cost,cluster_to_row,cluster_to_col,sorted_cost,piv_cost,proc_cost)
    CALL timestop(handle)

  END SUBROUTINE make_cluster_distribution

! *****************************************************************************
!> \brief assigns the clusters to processors, tryimg to balance the cost on the nodes
!> \param cluster_cost vector with the cost of each cluster
!> \param piv_cost pivoting vector sorting the cluster_cost
!> \param proc_cost cost per processor, on input 0 everywhere
!> \param cluster_assign assgnment of clusters on proc
!> \param nproc number of processor over which clusters are distributed
! *****************************************************************************
  SUBROUTINE  assign_clusters(cluster_cost,piv_cost,proc_cost,cluster_assign, nproc)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cluster_cost, piv_cost, &
                                                proc_cost, cluster_assign
    INTEGER                                  :: nproc

    CHARACTER(len=*), PARAMETER :: routineN = 'assign_clusters', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ilevel, offset, &
                                                piv_pcost(nproc), &
                                                sort_proc_cost(nproc)

     DO ilevel=1,SIZE(cluster_cost)/nproc
        sort_proc_cost(:)=proc_cost(:)
        CALL sort(sort_proc_cost,nproc,piv_pcost)
   
        offset=(SIZE(cluster_cost)/nproc-ilevel+1)*nproc+1
        DO i=1,nproc
           cluster_assign(piv_cost(offset-i))=piv_pcost(i)
           proc_cost(piv_pcost(i))=proc_cost(piv_pcost(i))+cluster_cost(piv_cost(offset-i))
        END DO
     END DO

 END SUBROUTINE assign_clusters
   
! *****************************************************************************
!> \brief recursive routine to cluster atoms.
!>        Low level uses a modified KMEANS algorithm
!>        recursion is used to reduce cost.
!>        each level will subdivide a cluster into smaller clusters
!>        If only a single split is necessary atoms are assigned to the current cluster
!> \param coord coordinates of the system
!> \param scaled_coord scaled coordinates
!> \param cell the cell_type
!> \param costs costs per atomic block
!> \param cluster_inds the atom_to cluster mapping
!> \param ncluster number of clusters still to be created on a given recursion level
!> \param icluster the index of the current cluster to be created
!> \param fin_cluster_cost total cost of the final clusters
! *****************************************************************************
  RECURSIVE SUBROUTINE cluster_recurse(coord, scaled_coord, cell, costs, cluster_inds, ncluster, icluster, fin_cluster_cost)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: coord, scaled_coord
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, DIMENSION(:), INTENT(IN)        :: costs
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: cluster_inds
    INTEGER, INTENT(INOUT)                   :: ncluster, icluster
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: fin_cluster_cost

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_recurse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ibeg, iend, maxv(1), &
                                                min_seed, natoms, nleft, &
                                                nsplits, seed, tot_cost
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cluster_cost, inds_tmp, &
                                                nat_cluster, ncluster_new, piv
    LOGICAL                                  :: found
    REAL(KIND=dp)                            :: balance, balance_new, conv

    natoms=SIZE(coord, 2)
    ! This is a bit of an arbitrary choice, simply a try to avoid too many clusters on large systems and too few for balancing on
    ! small systems or subclusters
    IF (natoms.le.1) THEN
       nsplits=1
    ELSE
       nsplits=MIN(MIN(MAX(6, INT(60.00/LOG(REAL(natoms)))), ncluster), natoms)
    ENDIF
    IF (nsplits==1) THEN
       icluster=icluster+1
       cluster_inds=icluster
       fin_cluster_cost(icluster)=SUM(costs)
    ELSE
       ALLOCATE(cluster_cost(nsplits), ncluster_new(nsplits), inds_tmp(natoms), piv(natoms), nat_cluster(nsplits))
       ! initialise some values
       cluster_cost=0; seed=300; found=.TRUE. ; min_seed=seed
       CALL kmeans(nsplits, coord, scaled_coord, cell, cluster_inds, nat_cluster, seed, conv)
       balance=MAXVAL(REAL(nat_cluster))/MINVAL(REAL(nat_cluster))

       ! If the system is small enough try to do better in terms of balancing number of atoms per cluster
       ! by changing the seed for the initial guess
       IF(natoms.lt.1000.AND. balance.gt.1.1)THEN
          found=.FALSE.
          DO i=1, 5
             IF(balance.gt.1.1)THEN
                 CALL kmeans(nsplits, coord, scaled_coord, cell, cluster_inds, nat_cluster, seed+i*40, conv)
                 balance_new=MAXVAL(REAL(nat_cluster))/MINVAL(REAL(nat_cluster))
                 IF(balance_new.LT.balance)THEN
                    balance=balance_new
                    min_seed=seed+i*40; 
                 END IF
             ELSE
                 found=.TRUE.
                 EXIT
             END IF
          END DO
       END IF
       !If we do not match the convergence than recompute at least the best assignment
       IF(.NOT.found)CALL kmeans(nsplits, coord, scaled_coord, cell, cluster_inds, nat_cluster, min_seed, conv)

       ! compute the cost of each cluster to decide how many splits have to be performed on the next lower level
       DO i=1, natoms
          cluster_cost(cluster_inds(i))=cluster_cost(cluster_inds(i))+costs(i)
       END DO
       tot_cost=SUM(cluster_cost)
       ! compute new splitting, can be done more elegant 
       ncluster_new(:)=ncluster*cluster_cost(:)/tot_cost
       nleft=ncluster-SUM(ncluster_new)
       ! As we won't have empty clusters, we can not have 0 as new size, so we correct for this at first
       DO i=1, nsplits
           IF(ncluster_new(i)==0)THEN
                ncluster_new(i)=1
                nleft=nleft-1
           END  IF
       END DO
       ! now comes the next part that the number of clusters will not match anymore, so try to correct in a meaningful way without
       ! introducing 0 sized blocks again
       IF(nleft.NE.0)THEN
         DO i=1, ABS(nleft)
            IF(nleft<0)THEN
               maxv=MINLOC(cluster_cost/ncluster_new)
               IF(ncluster_new(maxv(1)).NE.1)THEN
                  ncluster_new(maxv)=ncluster_new(maxv)-1
               ELSE
                  maxv=MAXLOC(ncluster_new)
                  ncluster_new(maxv)=ncluster_new(maxv)-1
               END IF
            ELSE
               maxv=MAXLOC(cluster_cost/ncluster_new)
               ncluster_new(maxv)=ncluster_new(maxv)+1
            END IF
         END DO
       END IF

       !Now get the permutations to sort the atoms in the nsplits clusters for the next level of iteration
       inds_tmp(:)=cluster_inds(:)
       CALL sort(inds_tmp, natoms, piv)

       ibeg=1; iend=0
       DO i=1, nsplits
          IF(nat_cluster(i)==0)CYCLE
          iend=iend+nat_cluster(i)
          CALL cluster_recurse(coord(:, piv(ibeg:iend)), scaled_coord(:, piv(ibeg:iend)), cell, costs(piv(ibeg:iend)),&
                               inds_tmp(ibeg:iend), ncluster_new(i), icluster, fin_cluster_cost)
          ibeg=ibeg+nat_cluster(i)
       END DO
       ! copy the sorted cluster IDs on the old layout, inds_tmp gets set at the lowest level of recursion
       cluster_inds(piv(:))=inds_tmp
       DEALLOCATE(cluster_cost, ncluster_new, inds_tmp, piv, nat_cluster)

   END IF

 END SUBROUTINE cluster_recurse


! *****************************************************************************
!> \brief A modified version of the kmeans algorithm.
!>        The assignment has a penalty function in case clusters become
!>        larger than average. Like this more even sized clusters are created
!>        trading it for locality
!> \param ncent number of centers to be created
!> \param coord coordinates
!> \param scaled_coord scaled coord
!> \param cell the cell_type
!> \param cluster atom to cluster assignment
!> \param nat_cl atoms per cluster
!> \param seed seed for the RNG. Algorithm might need multiple tries to deliver best results
!> \param tot_var the total variance of the clusters around the centers
! *****************************************************************************
 SUBROUTINE kmeans(ncent,coord,scaled_coord,cell,cluster,nat_cl,seed,tot_var)
    INTEGER                                  :: ncent
    REAL(KIND=dp), DIMENSION(:, :)           :: coord, scaled_coord
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, DIMENSION(:)                    :: cluster, nat_cl
    INTEGER                                  :: seed
    REAL(KIND=dp)                            :: tot_var

    INTEGER                                  :: i, ind, itn, j, nat, oldc
    LOGICAL                                  :: changed
    REAL(KIND=dp) :: average(3,ncent,2), cent_coord(3,ncent), devi, &
      deviat(ncent), dist, dvec(3), old_var, rn, scaled_cent(3,ncent), &
      var_cl(ncent)
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dmat
    REAL(KIND=dp), DIMENSION(3, 2)           :: initial_seed
    TYPE(rng_stream_type), POINTER           :: rng_stream

   initial_seed=REAL(seed,dp); nat=SIZE(coord,2)
   NULLIFY(rng_stream)
   ALLOCATE(dmat(ncent,nat))


   CALL create_rng_stream(rng_stream=rng_stream,&
        name="kmeans uniform distribution [0,1]",&
        distribution_type=UNIFORM,seed=initial_seed)

! try to find a clever initial guess with centers being somewhat distributed
   rn=next_random_number(rng_stream)
   ind=CEILING(rn*nat)
   cent_coord(:,1)=coord(:,ind)
   DO i=2,ncent
      DO
        rn=next_random_number(rng_stream)
        ind=CEILING(rn*nat)
        cent_coord(:,i)=coord(:,ind)
        devi=HUGE(1.0_dp)
        DO j=1,i-1
           dvec=pbc(cent_coord(:,j),cent_coord(:,i),cell)
           dist=SQRT(DOT_PRODUCT(dvec,dvec))
           IF(dist.LT.devi)devi=dist
        END DO
        rn=next_random_number(rng_stream)
        IF(rn.lt.devi**2/169.0)EXIT
      END DO
   END DO

! Now start the KMEANS but penalise it in case it starts packing too many atoms into a single set
! Unfoirtunatelz as this is dependent on what happend before it cant be parallel
   cluster=0; old_var=HUGE(1.0_dp)
   DO itn=1,1000
      changed=.FALSE.; var_cl=0.0_dp; tot_var=0.0_dp; nat_cl=0; deviat=0.0_dp
!      !$OMP PARALLEL DO PRIVATE(i,j,dvec)
      DO i=1,nat
         DO j=1,ncent
            dvec=pbc(cent_coord(:,j),coord(:,i),cell)
            dmat(j,i)=DOT_PRODUCT(dvec,dvec)
         END DO
      END DO
      DO i=1,nat
         devi=HUGE(1.0_dp); oldc=cluster(i)
         DO j=1,ncent
            dist=dmat(j,i)+MAX(nat_cl(j)**2/nat*ncent,nat/ncent)
            IF(dist.lt.devi)THEN
               devi=dist; cluster(i)=j
            END IF
         END DO
         deviat(cluster(i))=deviat(cluster(i))+SQRT(devi)
         nat_cl(cluster(i))=nat_cl(cluster(i))+1
         tot_var=tot_var+devi
         IF(oldc.NE.cluster(i))changed=.TRUE.
      END DO
      ! get the update of the centers done, add a new one in case one center lost all its atoms
      ! the algorithm would survive, but its nice to really create what you demand
      IF(tot_var.ge.old_var)EXIT
      IF(changed)THEN
        ! Here misery of computing the center of geometry of the clusters in PBC.
        ! The mapping on the unit circle allows to circumvent all problems
        average=0.0_dp
        DO i=1,SIZE(coord,2)
           average(:,cluster(i),1)=average(:,cluster(i),1)+COS(scaled_coord(:,i)*2.0_dp*pi)
           average(:,cluster(i),2)=average(:,cluster(i),2)+SIN(scaled_coord(:,i)*2.0_dp*pi)
        END DO

        DO i=1,ncent
           IF(nat_cl(i)==0)THEN
             rn=next_random_number(rng_stream)
             scaled_cent(:,i)=scaled_coord(:,CEILING(rn*nat))
           ELSE
             average(:,i,1)=average(:,i,1)/REAL(nat_cl(i),dp)
             average(:,i,2)=average(:,i,2)/REAL(nat_cl(i),dp)
             scaled_cent(:,i)=(ATAN2(-average(:,i,2),-average(:,i,1))+pi)/(2.0_dp*pi) 
             CALL scaled_to_real(cent_coord(:,i),scaled_cent(:,i),cell)
           END IF
        END DO
      ELSE
        EXIT
      END IF
   END DO

   CALL delete_rng_stream(rng_stream)
        
END SUBROUTINE kmeans


END MODULE distribution_methods
