!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Distribution methods for atoms, particles, or molecules
!> \par History
!>      - 1d-distribution of molecules and particles (Sep. 2003,MK)
!>      - 2d-distribution for Quickstep updated with molecules (Oct. 2003,MK)
!> \author MK (22.08.2003)
! *****************************************************************************
MODULE distribution_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc,&
                                             real_to_scaled
  USE cp_array_i_utils,                ONLY: cp_1d_i_p_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_dbcsr_interface,              ONLY: heap_fill,&
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
  USE distribution_optimize,           ONLY: distribute_2d_monte_carlo
  USE f77_blas
  USE input_constants,                 ONLY: model_block_count,&
                                             model_block_lmax,&
                                             model_block_surface
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE mathlib,                         ONLY: gcd,&
                                             lcm
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             get_molecule_kind_set,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
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

  PUBLIC :: distribute_molecules_1d,&
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
!> \param error variable to control error logging, stopping, ...
!>               see module cp_error_handling
!> \par History
!>      none
!> \author MK (Jun. 2003)
! *****************************************************************************
  SUBROUTINE distribute_molecules_1d(atomic_kind_set,particle_set,&
                                     local_particles,&
                                     molecule_kind_set,molecule_set,&
                                     local_molecules, force_env_section, &
                                     prev_molecule_kind_set, &
                                     prev_local_molecules, &
                                     error)

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
    TYPE(cp_error_type), INTENT(INOUT)       :: error

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

    CALL timeset(routineN,handle)


    has_prev_subsys_info = .FALSE.
    IF (PRESENT(prev_local_molecules) .AND.  &
        PRESENT(prev_molecule_kind_set)) THEN
       IF (ASSOCIATED(prev_local_molecules) .AND. &
           ASSOCIATED(prev_molecule_kind_set)) THEN
         has_prev_subsys_info = .TRUE.
       ENDIF
    ENDIF

    logger => cp_error_get_logger(error)

    group = logger%para_env%group
    mype = logger%para_env%mepos + 1
    npe = logger%para_env%num_pe

    ALLOCATE (workload_count(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "workload",npe*int_size)
    workload_count(:) = 0

    ALLOCATE (workload_fill(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "workload_old",npe*int_size)
    workload_fill(:) = 0

    nmolecule_kind = SIZE(molecule_kind_set)

    ALLOCATE (nmolecule_local(nmolecule_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "nmolecule_local",nmolecule_kind*int_size)
    nmolecule_local(:) = 0

    ALLOCATE (local_molecule(nmolecule_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "local_molecule",nmolecule_kind*int_size)

    nparticle_kind = SIZE(atomic_kind_set)

    ALLOCATE (nparticle_local(nparticle_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "nparticle_local",nparticle_kind*int_size)
    nparticle_local(:) = 0

    nbins=npe

    CALL heap_new (bin_heap_count, nbins, heap_error)
    CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error,&
         routineN, "Error creating heap.", error=error)
    CALL heap_fill (bin_heap_count,&
         (/(bin,bin=1,nbins)/), workload_count , heap_error)
    CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error,&
         routineN, "Error initially filling the heap.", error=error)

    CALL heap_new (bin_heap_fill, nbins, heap_error)
    CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error,&
         routineN, "Error creating heap.", error=error)
    CALL heap_fill (bin_heap_fill,&
         (/(bin,bin=1,nbins)/), workload_fill , heap_error)
    CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error,&
         routineN, "Error initially filling the heap.", error=error)

    DO imolecule_kind=1,nmolecule_kind

      molecule_kind => molecule_kind_set(imolecule_kind)

      NULLIFY (molecule_list)

!     *** Get the number of molecules and the number of ***
!     *** atoms in each molecule of that molecular kind ***

      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             molecule_list=molecule_list,&
                             natom=natom,&
                             nsgf=nsgf)

!     *** Consider the number of atoms or basis ***
!     *** functions which depends on the method ***

      nload = MAX(natom,nsgf)
      nmolecule = SIZE(molecule_list)

!     *** Get the number of local molecules of the current molecule kind ***

      DO imolecule=1,nmolecule
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
            CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error,&
                 routineN, "Error getting topmost heap element.", error=error)
            CALL cp_assert (found, cp_fatal_level, cp_internal_error,&
                 routineN, "No topmost heap element found.", error=error)

            ipe = bin
            CALL cp_assert (bin_price==workload_count(ipe), cp_fatal_level, cp_internal_error,&
                 routineN, "inconsistent heap", error=error)

            workload_count(ipe) = workload_count(ipe) + nload
            IF (ipe == mype) THEN
              nmolecule_local(imolecule_kind) = nmolecule_local(imolecule_kind) + 1
            END IF

            bin_price = workload_count(ipe)
            CALL heap_reset_first (bin_heap_count, bin_price, heap_error)
            CALL cp_assert (.NOT. heap_error, cp_warning_level, cp_internal_error,&
                 routineN, "Error setting price of top heap element.", error=error)
         END IF
      END DO


!     *** Distribute the molecules ***
      n = nmolecule_local(imolecule_kind)

      IF (n > 0) THEN
        ALLOCATE (local_molecule(imolecule_kind)%array(n),STAT=istat)
        IF (istat /= 0) THEN
           CALL stop_memory(routineN,moduleN,__LINE__,&
                            "local_molecule(imolecule_kind)%array",&
                            n*int_size)
        END IF
      ELSE
        NULLIFY (local_molecule(imolecule_kind)%array)
      END IF

      imolecule_local = 0
      DO imolecule=1,nmolecule
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
            CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error,&
                 routineN, "Error getting topmost heap element.", error=error)
            CALL cp_assert (found, cp_fatal_level, cp_internal_error,&
                 routineN, "No topmost heap element found.", error=error)

            ipe = bin
            CALL cp_assert (bin_price==workload_fill(ipe), cp_fatal_level, cp_internal_error,&
                 routineN, "inconsistent heap", error=error)

            workload_fill(ipe) = workload_fill(ipe) + nload
            is_local = (ipe == mype)
         ENDIF
         IF (is_local) THEN
            imolecule_local = imolecule_local + 1
            molecule_a = molecule_list(imolecule)
            local_molecule(imolecule_kind)%array(imolecule_local) = molecule_a
            DO iatom=1,natom
              atom_a = molecule_set(molecule_a)%first_atom + iatom - 1

              CALL get_atomic_kind(atomic_kind=particle_set(atom_a)%atomic_kind,&
                                   kind_number=kind_a)
              nparticle_local(kind_a) = nparticle_local(kind_a) + 1
            END DO
         END IF
         IF (.NOT. has_prev_subsys_info) THEN
            bin_price = workload_fill(ipe)
            CALL heap_reset_first (bin_heap_fill, bin_price, heap_error)
            CALL cp_assert (.NOT. heap_error, cp_warning_level, cp_internal_error,&
                 routineN, "Error setting price of top heap element.", error=error)
         END IF
      END DO

    END DO

    CALL cp_assert (ALL(workload_fill.EQ.workload_count), cp_fatal_level, cp_internal_error,&
              routineN, "Inconsistent heaps encountered", error=error)

    CALL heap_release (bin_heap_count, heap_error)
    CALL cp_assert (.NOT. heap_error, cp_warning_level, cp_internal_error,&
         routineN, "Error releasing heap.", error=error)
    CALL heap_release (bin_heap_fill, heap_error)
    CALL cp_assert (.NOT. heap_error, cp_warning_level, cp_internal_error,&
         routineN, "Error releasing heap.", error=error)

!   *** Create the local molecule structure ***

    CALL distribution_1d_create(local_molecules,&
                                n_el=nmolecule_local,&
                                para_env=logger%para_env,error=error)

!   *** Create the local particle structure ***

    CALL distribution_1d_create(local_particles,&
                                n_el=nparticle_local,&
                                para_env=logger%para_env,error=error)

!   *** Store the generated local molecule and particle distributions ***

    nparticle_local(:) = 0

    DO imolecule_kind=1,nmolecule_kind

      IF (nmolecule_local(imolecule_kind) == 0) CYCLE

      local_molecules%list(imolecule_kind)%array(:) =&
        local_molecule(imolecule_kind)%array(:)

      molecule_kind => molecule_kind_set(imolecule_kind)

      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             natom=natom)

      DO imolecule=1,nmolecule_local(imolecule_kind)
        molecule_a = local_molecule(imolecule_kind)%array(imolecule)
        DO iatom=1,natom
          atom_a = molecule_set(molecule_a)%first_atom + iatom - 1
          CALL get_atomic_kind(atomic_kind=particle_set(atom_a)%atomic_kind,&
                               kind_number=kind_a)
          nparticle_local(kind_a) = nparticle_local(kind_a) + 1
          local_particles%list(kind_a)%array(nparticle_local(kind_a)) = atom_a
        END DO
      END DO

    END DO

!   *** Print distribution, if requested ***

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         force_env_section,"PRINT%DISTRIBUTION1D",error=error),cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger,force_env_section,"PRINT%DISTRIBUTION1D",&
            extension=".Log",error=error)

       iw = output_unit
       IF (output_unit<0) iw = cp_logger_get_default_unit_nr(logger,LOCAL=.TRUE.)

!     *** Print molecule distribution ***

      ALLOCATE (work(npe),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "work",npe*int_size)
      work(:) = 0

      work(mype) = SUM(nmolecule_local)
      CALL mp_sum(work,group)

      IF (output_unit>0) THEN
         WRITE (UNIT=output_unit,&
              FMT="(/,T2,A,T51,A,/,(T52,I6,T73,I8))")&
              "DISTRIBUTION OF THE MOLECULES",&
              "Process    Number of molecules",&
              (ipe-1,work(ipe),ipe=1,npe)
         WRITE (UNIT=output_unit,FMT="(T55,A3,T73,I8)")&
              "Sum",SUM(work)
         CALL m_flush(output_unit)
      END IF

      CALL mp_sync(group)

      DO ipe=1,npe
         IF (ipe == mype) THEN
            WRITE (UNIT=iw,FMT="(/,T3,A)")&
                 "Process   Kind   Local molecules (global indices)"
            DO imolecule_kind=1,nmolecule_kind
               IF (imolecule_kind == 1) THEN
                  WRITE (UNIT=iw,FMT="(T4,I6,2X,I5,(T21,10I6))")&
                       ipe-1,imolecule_kind,&
                       (local_molecules%list(imolecule_kind)%array(imolecule),&
                       imolecule=1,nmolecule_local(imolecule_kind))
               ELSE
                  WRITE (UNIT=iw,FMT="(T12,I5,(T21,10I6))")&
                       imolecule_kind,&
                       (local_molecules%list(imolecule_kind)%array(imolecule),&
                       imolecule=1,nmolecule_local(imolecule_kind))
               END IF
            END DO
         END IF
         CALL m_flush(iw)
         CALL mp_sync(group)
      END DO

!     *** Print particle distribution ***

      work(:) = 0

      work(mype) = SUM(nparticle_local)
      CALL mp_sum(work,group)

      IF (output_unit>0) THEN
         WRITE (UNIT=output_unit,&
              FMT="(/,T2,A,T51,A,/,(T52,I6,T73,I8))")&
              "DISTRIBUTION OF THE PARTICLES",&
              "Process    Number of particles",&
              (ipe-1,work(ipe),ipe=1,npe)
         WRITE (UNIT=output_unit,FMT="(T55,A3,T73,I8)")&
              "Sum",SUM(work)
         CALL m_flush(output_unit)
      END IF

      CALL mp_sync(group)

      DO ipe=1,npe
         IF (ipe == mype) THEN
            WRITE (UNIT=iw,FMT="(/,T3,A)")&
                 "Process   Kind   Local particles (global indices)"
            DO iparticle_kind=1,nparticle_kind
               IF (iparticle_kind == 1) THEN
                  WRITE (UNIT=iw,FMT="(T4,I6,2X,I5,(T20,10I6))")&
                       ipe-1,iparticle_kind,&
                       (local_particles%list(iparticle_kind)%array(iatom),&
                       iatom=1,nparticle_local(iparticle_kind))
               ELSE
                  WRITE (UNIT=iw,FMT="(T12,I5,(T20,10I6))")&
                       iparticle_kind,&
                       (local_particles%list(iparticle_kind)%array(iatom),&
                       iatom=1,nparticle_local(iparticle_kind))
               END IF
            END DO
         END IF
         CALL m_flush(iw)
         CALL mp_sync(group)
      END DO
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"work")

      CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
           "PRINT%DISTRIBUTION1D",error=error)
   END IF
!   *** Release work storage ***

    DEALLOCATE (workload_count,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "workload")

    DEALLOCATE (workload_fill,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "workload_old")

    DEALLOCATE (nmolecule_local,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "nmolecule_local")

    DEALLOCATE (nparticle_local,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "nparticle_local")

    DO imolecule_kind=1,nmolecule_kind
      IF (ASSOCIATED(local_molecule(imolecule_kind)%array)) THEN
        DEALLOCATE (local_molecule(imolecule_kind)%array,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,&
                           "local_molecule(imolecule_kind)%array")
        END IF
      END IF
    END DO
    DEALLOCATE (local_molecule,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "local_molecule")

    CALL timestop(handle)

  END SUBROUTINE distribute_molecules_1d

! *****************************************************************************
!> \brief Distributes the particle pairs creating a 2d distribution optimally
!>      suited for quickstep
!> \param cell ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param distribution_2d the distribution that will be created by this
!>                         method
!> \param blacs_env the parallel environement at the basis of the
!>                   distribution
!> \param force_env_section ...
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \par History
!>      - local_rows & cols blocksize optimizations (Aug. 2003,MK)
!>      - cleanup of distribution_2d (Sep. 2003,fawzi)
!>      - update for molecules (Oct. 2003,MK)
!> \author fawzi (Feb. 2003)
!> \note
!>      Intermediate generation of a 2d distribution of the molecules, but
!>      only the corresponding particle (atomic) distribution is currently
!>      used. The 2d distribution of the molecules is deleted, but may easily
!>      be recovered (MK).
! *****************************************************************************
  SUBROUTINE distribute_molecules_2d(cell, atomic_kind_set,particle_set,&
                                     qs_kind_set, molecule_kind_set,molecule_set,&
                                     distribution_2d,blacs_env,force_env_section,&
                                     error)
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
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribute_molecules_2d', &
      routineP = moduleN//':'//routineN

    INTEGER :: cluster_price, cost_model, group, handle, iatom, iatom_mol, &
      iatom_one, ikind, imol, imolecule, imolecule_kind, iparticle_kind, &
      ipcol, iprow, istat, iw, kind_a, mypcol, myprow, n, natom, natom_mol, &
      nclusters, nele, nmolecule, nmolecule_kind, nparticle_kind, npcol, &
      nprow, nsgf, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:) :: cluster_list, cluster_prices, &
      col_ele, nparticle_local_col, nparticle_local_row, row_ele, work
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: cost_ele_pair
    INTEGER, DIMENSION(:), POINTER :: cluster_col_distribution, &
      cluster_row_distribution, col_distribution, lmax_basis, molecule_list, &
      row_distribution
    LOGICAL :: basic_optimization, basic_spatial_optimization, failure, &
      molecular_distribution, skip_optimization
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pbc_scaled_coords
    REAL(KIND=dp), DIMENSION(3)              :: center
    TYPE(cp_1d_i_p_type), DIMENSION(:), &
      POINTER                                :: local_particle_col, &
                                                local_particle_row
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(section_vals_type), POINTER         :: distribution_section, &
                                                mc_section

!...

    failure = .FALSE.

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)

    distribution_section => section_vals_get_subs_vals(force_env_section,"DFT%QS%DISTRIBUTION",error=error)
    mc_section => section_vals_get_subs_vals(force_env_section,"DFT%QS%DISTRIBUTION%2D_MC",error=error)

    CALL section_vals_val_get(distribution_section,"2D_MOLECULAR_DISTRIBUTION",l_val=molecular_distribution,error=error)
    CALL section_vals_val_get(distribution_section,"SKIP_OPTIMIZATION",l_val=skip_optimization,error=error)
    CALL section_vals_val_get(distribution_section,"BASIC_OPTIMIZATION", l_val=basic_optimization,error=error)
    CALL section_vals_val_get(distribution_section,"BASIC_SPATIAL_OPTIMIZATION", l_val=basic_spatial_optimization,error=error)

    CALL section_vals_val_get(distribution_section,"COST_MODEL",i_val=cost_model,error=error)
    !

    group = blacs_env%para_env%group
    myprow = blacs_env%mepos(1) + 1
    mypcol = blacs_env%mepos(2) + 1
    nprow = blacs_env%num_pe(1)
    npcol = blacs_env%num_pe(2)

    nmolecule_kind = SIZE(molecule_kind_set)
    CALL get_molecule_kind_set(molecule_kind_set,nmolecule=nmolecule)

    nparticle_kind = SIZE(atomic_kind_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,natom=natom)

    !
    ! we need to generate two representations of the distribution, one as a straight array with global particles
    ! one ordered wrt to kinds and only listing the local particles
    !
    ALLOCATE (row_distribution(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "row_distribution",natom*int_size)

    ALLOCATE (col_distribution(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "col_distribution",natom*int_size)

    ALLOCATE (local_particle_row(nparticle_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "local_particle_row",&
                                     nparticle_kind*int_size)

    ALLOCATE (local_particle_col(nparticle_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "local_particle_col",&
                                     nparticle_kind*int_size)

    ALLOCATE (nparticle_local_row(nparticle_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "nparticle_local_row",&
                                     nparticle_kind*int_size)

    ALLOCATE (nparticle_local_col(nparticle_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "nparticle_local_col",&
                                     nparticle_kind*int_size)

    IF (basic_optimization .OR. basic_spatial_optimization) THEN


       IF (molecular_distribution) THEN
           nclusters=nmolecule
       ELSE
           nclusters=natom
       ENDIF

       ALLOCATE (cluster_list (nclusters), stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "cluster_list", nclusters*int_size)
       ALLOCATE (cluster_prices (nclusters), stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "cluster_prices", nclusters*int_size)
       ALLOCATE (cluster_row_distribution (nclusters), stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "cluster_prices", nclusters*int_size)
       ALLOCATE (cluster_col_distribution (nclusters), stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "cluster_prices", nclusters*int_size)


       ! Fill in the clusters and their prices
       CALL section_vals_val_get(distribution_section,"COST_MODEL",i_val=cost_model,error=error)
       IF (.NOT.molecular_distribution) THEN
          DO iatom = 1, natom
             CALL cp_assert (iatom .LE. nclusters, cp_fatal_level,&
                  cp_internal_error, routineN, "Bounds error", error=error)
             CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
             cluster_list(iatom) = iatom
             SELECT CASE (cost_model)
             CASE (model_block_count)
                CALL get_qs_kind(qs_kind_set(ikind),nsgf=nsgf)
                cluster_price = nsgf
             CASE (model_block_lmax)
                CALL get_qs_kind(qs_kind_set(ikind),orb_basis_set=orb_basis_set)
                CALL get_gto_basis_set(orb_basis_set, lmax=lmax_basis)
                cluster_price = MAXVAL (lmax_basis)
             CASE default
                CALL get_qs_kind(qs_kind_set(ikind),orb_basis_set=orb_basis_set)
                CALL get_gto_basis_set(orb_basis_set, lmax=lmax_basis)
                cluster_price = 8 + (MAXVAL (lmax_basis)**2)
             END SELECT
             cluster_prices(iatom) = cluster_price
          ENDDO
       ELSE
          imol = 0
          DO imolecule_kind=1,nmolecule_kind
             molecule_kind => molecule_kind_set(imolecule_kind)
             CALL get_molecule_kind(molecule_kind=molecule_kind, molecule_list=molecule_list, natom=natom_mol)
             DO imolecule=1,SIZE(molecule_list)
                imol = imol + 1
                cluster_list(imol) = imol
                cluster_price = 0
                DO iatom_mol=1,natom_mol
                   iatom = molecule_set(molecule_list(imolecule))%first_atom + iatom_mol - 1
                   CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
                   SELECT CASE (cost_model)
                   CASE (model_block_count)
                      CALL get_qs_kind(qs_kind_set(ikind),nsgf=nsgf)
                      cluster_price = cluster_price + nsgf
                   CASE (model_block_lmax)
                      CALL get_qs_kind(qs_kind_set(ikind), orb_basis_set=orb_basis_set)
                      CALL get_gto_basis_set(orb_basis_set, lmax=lmax_basis)
                      cluster_price = cluster_price + MAXVAL (lmax_basis)
                   CASE default
                      CALL get_qs_kind(qs_kind_set(ikind), orb_basis_set=orb_basis_set)
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
          CALL make_basic_distribution (cluster_list, cluster_prices,&
               nprow, cluster_row_distribution, npcol, cluster_col_distribution, error=error)
       ELSE ! basic_spatial_optimization
          ALLOCATE(pbc_scaled_coords(3,nclusters))
          IF (.NOT. molecular_distribution) THEN
             ! just scaled coords
             DO iatom=1, natom
                CALL real_to_scaled(pbc_scaled_coords(:,iatom),pbc(particle_set(iatom)%r(:),cell),cell)
             ENDDO
          ELSE
             ! use scaled coords of geometric center, folding when appropriate
             imol = 0
             DO imolecule_kind=1,nmolecule_kind
                molecule_kind => molecule_kind_set(imolecule_kind)
                CALL get_molecule_kind(molecule_kind=molecule_kind, molecule_list=molecule_list, natom=natom_mol)
                DO imolecule=1,SIZE(molecule_list)
                   imol = imol + 1
                   iatom_one = molecule_set(molecule_list(imolecule))%first_atom
                   center = 0.0_dp
                   DO iatom_mol=1,natom_mol
                      iatom = molecule_set(molecule_list(imolecule))%first_atom + iatom_mol - 1
                      center = center + &
                        pbc(particle_set(iatom)%r(:)-particle_set(iatom_one)%r(:),cell)+particle_set(iatom_one)%r(:)
                   ENDDO
                   center = center / natom_mol
                   CALL real_to_scaled(pbc_scaled_coords(:,imol),pbc(center,cell),cell)
                ENDDO
             ENDDO
          ENDIF

          CALL make_basic_spatial_distribution (pbc_scaled_coords, cluster_prices,&
                         nprow, cluster_row_distribution, npcol, cluster_col_distribution, error)

          DEALLOCATE(pbc_scaled_coords)
       ENDIF

       ! And assign back
       IF (.NOT. molecular_distribution) THEN
          row_distribution = cluster_row_distribution
          col_distribution = cluster_col_distribution
       ELSE
          imol = 0
          DO imolecule_kind=1,nmolecule_kind
             molecule_kind => molecule_kind_set(imolecule_kind)
             CALL get_molecule_kind(molecule_kind=molecule_kind, molecule_list=molecule_list, natom=natom_mol)
             DO imolecule=1,SIZE(molecule_list)
                imol = imol + 1
                DO iatom_mol=1,natom_mol
                   iatom = molecule_set(molecule_list(imolecule))%first_atom + iatom_mol - 1
                   row_distribution(iatom) = cluster_row_distribution(imol)
                   col_distribution(iatom) = cluster_col_distribution(imol)
                ENDDO
             ENDDO
          ENDDO
       ENDIF

       ! cleanup
       DEALLOCATE (cluster_list, stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "cluster_list", nclusters*int_size)
       DEALLOCATE (cluster_prices, stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "cluster_prices", nclusters*int_size)
       DEALLOCATE (cluster_row_distribution, stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "cluster_prices", nclusters*int_size)
       DEALLOCATE (cluster_col_distribution, stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "cluster_prices", nclusters*int_size)

    ELSE
       !
       ! this could be a lot of memory, even in QS case
       ! one solution is to turn this into a function, but this would presumably
       ! slow down the distribute_2d_monte_carlo
       ! making this a function might be needed if the cost of a block depends on its
       ! assigned row / col.
       !
       IF (molecular_distribution) THEN
           nele=nmolecule
       ELSE
           nele=natom
       ENDIF

       ALLOCATE (row_ele(nele),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"row_ele",nele*int_size)
       ALLOCATE (col_ele(nele),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"col_ele",nele*int_size)
       ALLOCATE(cost_ele_pair(nele,nele),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"cost_ele_pair",natom*natom*int_size)

       !
       ! have the cost model provide the costs for each atomic block
       !
       output_unit = cp_print_key_unit_nr(logger,force_env_section,"PRINT%DISTRIBUTION",&
                                                            extension=".Log",error=error)
       CALL distribute_2d_cost_model(cost_ele_pair,distribution_section,particle_set,&
                                     qs_kind_set,molecule_kind_set,molecule_set,output_unit, error=error)
       CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
                                             "PRINT%DISTRIBUTION",error=error)

       !
       ! optimize the distribution
       !
       CALL distribute_2d_monte_carlo(cost_ele_pair,nprow,npcol,nele,row_ele,col_ele, &
                                      mc_section,blacs_env%para_env,error)

       DEALLOCATE(cost_ele_pair)

       ! prepare to add the distribution to the data, i.e. get the atomic distribution in all cases

       IF (.NOT. molecular_distribution) THEN
          row_distribution = row_ele
          col_distribution = col_ele
       ELSE
          imol = 0
          DO imolecule_kind=1,nmolecule_kind
             molecule_kind => molecule_kind_set(imolecule_kind)
             CALL get_molecule_kind(molecule_kind=molecule_kind, molecule_list=molecule_list, natom=natom_mol)
             DO imolecule=1,SIZE(molecule_list)
                imol = imol + 1
                DO iatom_mol=1,natom_mol
                   iatom = molecule_set(molecule_list(imolecule))%first_atom + iatom_mol - 1
                   row_distribution(iatom) = row_ele(imol)
                   col_distribution(iatom) = col_ele(imol)
                ENDDO
             ENDDO
          ENDDO
       ENDIF

       DEALLOCATE(row_ele)
       DEALLOCATE(col_ele)
    ENDIF


    ! prepare the lists of local particles

    ! count local particles of a given kind
    nparticle_local_col=0
    nparticle_local_row=0
    DO iatom=1,natom
       CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,kind_number=kind_a)
       IF (row_distribution(iatom)==myprow) nparticle_local_row(kind_a)=nparticle_local_row(kind_a)+1
       IF (col_distribution(iatom)==mypcol) nparticle_local_col(kind_a)=nparticle_local_col(kind_a)+1
    ENDDO

    ! allocate space
    DO iparticle_kind=1,nparticle_kind
      n = nparticle_local_row(iparticle_kind)
      ALLOCATE (local_particle_row(iparticle_kind)%array(n),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "local_particle_row(iparticle_kind)%array",n*int_size)

      n = nparticle_local_col(iparticle_kind)
      ALLOCATE (local_particle_col(iparticle_kind)%array(n),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "local_particle_col(iparticle_kind)%array",n*int_size)
    ENDDO

    ! store
    nparticle_local_col=0
    nparticle_local_row=0
    DO iatom=1,natom
       CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,kind_number=kind_a)
       IF (row_distribution(iatom)==myprow) THEN
           nparticle_local_row(kind_a)=nparticle_local_row(kind_a)+1
           local_particle_row(kind_a)%array(nparticle_local_row(kind_a)) = iatom
       ENDIF
       IF (col_distribution(iatom)==mypcol) THEN
           nparticle_local_col(kind_a)=nparticle_local_col(kind_a)+1
           local_particle_col(kind_a)%array(nparticle_local_col(kind_a)) = iatom
       ENDIF
    ENDDO

!   *** Generate the 2d distribution structure  but take care of the zero offsets required
    row_distribution = row_distribution - 1
    col_distribution = col_distribution - 1
    CALL distribution_2d_create(distribution_2d,&
                                row_distribution_ptr=row_distribution,&
                                col_distribution_ptr=col_distribution,&
                                local_rows_ptr=local_particle_row,&
                                local_cols_ptr=local_particle_col,&
                                blacs_env=blacs_env,&
                                error=error)

    NULLIFY (local_particle_row)
    NULLIFY (local_particle_col)
    NULLIFY (row_distribution)
    NULLIFY (col_distribution)

!   *** Print distribution, if requested ***
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         force_env_section,"PRINT%DISTRIBUTION",error=error),cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger,force_env_section,"PRINT%DISTRIBUTION",&
            extension=".Log",error=error)

!     *** Print row distribution ***

      ALLOCATE (work(nprow),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "work",nprow*int_size)
      END IF
      work(:) = 0

      IF (mypcol == 1) work(myprow) = SUM(distribution_2d%n_local_rows)

      CALL mp_sum(work,group)

      IF (output_unit>0) THEN
         WRITE (UNIT=output_unit,&
              FMT="(/,T2,A,/,T15,A,/,(T16,I10,T41,I10,T71,I10))")&
              "DISTRIBUTION OF THE PARTICLES (ROWS)",&
              "Process row      Number of particles         Number of matrix rows",&
              (iprow-1,work(iprow),-1,iprow=1,nprow)
         WRITE (UNIT=output_unit,FMT="(T23,A3,T41,I10,T71,I10)")&
              "Sum",SUM(work),-1
         CALL m_flush(output_unit)
      END IF

      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"work")
      END IF

!     *** Print column distribution ***

      ALLOCATE (work(npcol),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "work",npcol*int_size)
      END IF
      work(:) = 0

      IF (myprow == 1) work(mypcol) = SUM(distribution_2d%n_local_cols)

      CALL mp_sum(work,group)

      IF (output_unit>0) THEN
         WRITE (UNIT=output_unit,&
              FMT="(/,T2,A,/,T15,A,/,(T16,I10,T41,I10,T71,I10))")&
              "DISTRIBUTION OF THE PARTICLES (COLUMNS)",&
              "Process col      Number of particles      Number of matrix columns",&
              (ipcol-1,work(ipcol),-1,ipcol=1,npcol)
         WRITE (UNIT=output_unit,FMT="(T23,A3,T41,I10,T71,I10)")&
              "Sum",SUM(work),-1
         CALL m_flush(output_unit)
      END IF

      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"work")
      END IF

      CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
           "PRINT%DISTRIBUTION",error=error)
    END IF

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         force_env_section,"PRINT%DISTRIBUTION2D",error=error),cp_p_file)) THEN

         iw = cp_logger_get_default_unit_nr(logger,LOCAL=.TRUE.)
         CALL distribution_2d_write(distribution_2d,&
                                    unit_nr=iw,&
                                    local=.TRUE.,&
                                    long_description=.TRUE.,&
                                    error=error)

    ENDIF

!   *** Release work storage ***

    DEALLOCATE (nparticle_local_row,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "nparticle_local_row")

    DEALLOCATE (nparticle_local_col,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "nparticle_local_col")

    CALL timestop(handle)

  END SUBROUTINE distribute_molecules_2d

! *****************************************************************************
!> \brief this is the crucial part in order to get a good distribution. The cost
!>      model must model properly the cost of each cost_ele_pair
!> \param cost_ele_pair ...
!> \param distribution_section ...
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param output_unit ...
!> \param error ...
!> \par History
!>      04.2007 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE distribute_2d_cost_model(cost_ele_pair,distribution_section,particle_set,&
                                      qs_kind_set,molecule_kind_set,molecule_set,output_unit, error)

    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: cost_ele_pair
    TYPE(section_vals_type), POINTER         :: distribution_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribute_2d_cost_model', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nscale = 1000

    CHARACTER(LEN=20)                        :: cost_model_used
    INTEGER :: cost_block_ab, cost_model, handle, iatom_a, iatom_b, &
      iatom_mol_a, iatom_mol_b, ikind, imol_a, imol_b, imolecule_a, &
      imolecule_b, imolecule_kind_a, imolecule_kind_b, iparticle, natom, &
      natom_mol_a, natom_mol_b, nmolecule_kind
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: lmax, nsgf
    INTEGER, DIMENSION(:), POINTER           :: lmax_basis, molecule_list_a, &
                                                molecule_list_b
    LOGICAL                                  :: failure, include_ab, &
                                                molecular_distribution, &
                                                symmetric
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_a, &
                                                molecule_kind_b

! use kilointegrals

    CALL timeset(routineN,handle)

    failure = .FALSE.
    cost_ele_pair = 0

    natom=SIZE(particle_set)

    CALL section_vals_val_get(distribution_section,"2D_MOLECULAR_DISTRIBUTION",l_val=molecular_distribution,error=error)
    CALL section_vals_val_get(distribution_section,"COST_MODEL",i_val=cost_model,error=error)
    CALL section_vals_val_get(distribution_section,"SYMMETRIC",l_val=symmetric,error=error)

    IF (output_unit>0) THEN
        ! Cost Model (CM) info
        WRITE(output_unit,'(T2,A)') ""
        WRITE(output_unit,'(T2,A)') "2D_CM| distribution_2d cost model info"
        WRITE(output_unit,'(T2,A,T78,A)') "2D_CM| molecular_distribution ",MERGE("YES"," NO",molecular_distribution)
        WRITE(output_unit,'(T2,A,T78,A)') "2D_CM| account for symmetry ",MERGE("YES"," NO",symmetric)

        cost_model_used="UNKNOWN"
        SELECT CASE(cost_model)
        CASE(model_block_count)
             cost_model_used="BLOCK COUNT"
        CASE(model_block_surface)
             cost_model_used="BLOCK SURFACE"
        CASE(model_block_lmax)
             cost_model_used="BLOCK LMAX"
        CASE DEFAULT
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        END SELECT

        WRITE(output_unit,'(T2,A,T61,A)') "2D_CM| cost model ",ADJUSTR(cost_model_used)
    ENDIF

    SELECT CASE(cost_model)
    CASE(model_block_count)
       ! do nothing
    CASE(model_block_surface)
       ALLOCATE(nsgf(natom))
       DO iparticle=1,natom
          CALL get_atomic_kind(particle_set(iparticle)%atomic_kind, kind_number=ikind)
          CALL get_qs_kind(qs_kind_set(ikind), nsgf=nsgf(iparticle))
       ENDDO
    CASE(model_block_lmax)
        ALLOCATE(lmax(natom))
        DO iparticle=1,natom
           CALL get_atomic_kind(particle_set(iparticle)%atomic_kind, kind_number=ikind)
           CALL get_qs_kind(qs_kind_set(ikind), orb_basis_set=orb_basis_set)
           CALL get_gto_basis_set(orb_basis_set,lmax=lmax_basis)
           lmax(iparticle)=MAXVAL(lmax_basis)
        ENDDO
    END SELECT

    !
    ! loop over the molecules and atoms in the molecules so that we know both molecule index and atom index
    !
    nmolecule_kind = SIZE(molecule_kind_set)
    imol_a = 0
    DO imolecule_kind_a=1,nmolecule_kind
       molecule_kind_a => molecule_kind_set(imolecule_kind_a)
       CALL get_molecule_kind(molecule_kind=molecule_kind_a, molecule_list=molecule_list_a, natom=natom_mol_a)
       DO imolecule_a=1,SIZE(molecule_list_a)
          imol_a = imol_a + 1
          DO iatom_mol_a=1,natom_mol_a
             iatom_a = molecule_set(molecule_list_a(imolecule_a))%first_atom + iatom_mol_a - 1

             imol_b = 0
             DO imolecule_kind_b=1,nmolecule_kind
                molecule_kind_b => molecule_kind_set(imolecule_kind_b)
                CALL get_molecule_kind(molecule_kind=molecule_kind_b, molecule_list=molecule_list_b, natom=natom_mol_b)
                DO imolecule_b=1,SIZE(molecule_list_b)
                   imol_b = imol_b + 1
                   DO iatom_mol_b=1,natom_mol_b
                      iatom_b = molecule_set(molecule_list_b(imolecule_b))%first_atom + iatom_mol_b - 1

                      ! is this block actually include in the S matrix
                      IF (symmetric) THEN
                        IF (iatom_a > iatom_b) THEN
                          include_ab = (MODULO(iatom_a + iatom_b,2) /= 0)
                        ELSE
                          include_ab = (MODULO(iatom_a + iatom_b,2) == 0)
                        END IF
                      ELSE
                        include_ab = .TRUE.
                      END IF

                      ! here we estimate the computational cost of the block ab based on the model of choice
                      IF (include_ab) THEN
                         SELECT CASE(cost_model)
                         CASE(model_block_count)
                             cost_block_ab = 1
                         CASE(model_block_surface)
                             cost_block_ab = nsgf(iatom_a)*nsgf(iatom_b)
                         CASE(model_block_lmax)
                             cost_block_ab = lmax(iatom_a)+lmax(iatom_b)+1
                         CASE DEFAULT
                             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                         END SELECT
                      ELSE
                         cost_block_ab=0
                      ENDIF

                      IF (molecular_distribution) THEN
                          cost_ele_pair(imol_a,imol_b) = cost_ele_pair(imol_a,imol_b) + cost_block_ab
                      ELSE
                          cost_ele_pair(iatom_a,iatom_b) = cost_block_ab
                      ENDIF

                   ENDDO ! iatom_mol_b
                ENDDO ! imolecule_b
             ENDDO ! imolecule_kind_b

          ENDDO ! iatom_mol_a
       ENDDO ! imolecule_a
    ENDDO ! imolecule_kind_a

    CALL timestop(handle)

  END SUBROUTINE distribute_2d_cost_model


! *****************************************************************************
!> \brief Creates a basic distribution
!> \param cluster_list ...
!> \param cluster_prices ...
!> \param nprows ...
!> \param row_distribution ...
!> \param npcols ...
!> \param col_distribution ...
!> \param error ...
!> \par History
!> - Created 2010-08-06 UB
! *****************************************************************************
  SUBROUTINE make_basic_distribution (cluster_list, cluster_prices,&
       nprows, row_distribution, npcols, col_distribution, error)
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: cluster_list, cluster_prices
    INTEGER, INTENT(IN)                      :: nprows
    INTEGER, DIMENSION(:), INTENT(OUT)       :: row_distribution
    INTEGER, INTENT(IN)                      :: npcols
    INTEGER, DIMENSION(:), INTENT(OUT)       :: col_distribution
    TYPE(cp_error_type), INTENT(INOUT)       :: error

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
    CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error,&
         routineN, "Error creating heap.", error=error)
    CALL heap_fill (bin_heap,&
         (/(bin,bin=0,nbins-1)/), (/(0, bin=1,nbins)/), heap_error)
    CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error,&
         routineN, "Error initially filling the heap.", error=error)
    !
    nclusters = SIZE (cluster_list)
    ! Put the most expensive cluster in the bin with the smallest
    ! price and repeat.
    DO cluster_index = nclusters, 1, -1
       cluster = cluster_list(cluster_index)
       CALL heap_get_first (bin_heap, bin, bin_price, found, heap_error)
       CALL cp_assert (.NOT. heap_error, cp_fatal_level, cp_internal_error,&
            routineN, "Error getting topmost heap element.", error=error)
       CALL cp_assert (found, cp_fatal_level, cp_internal_error,&
            routineN, "No topmost heap element found.", error=error)
       !
       failure = .FALSE.
       prow = INT (bin*pgrid_gcd / npcols)
       CALL cp_assert (prow .LT. nprows, cp_fatal_level, cp_internal_error,&
            routineN, "Invalid process row.", failure=failure, error=error)
       failure = .TRUE.
       pcol = INT (bin*pgrid_gcd / nprows)
       CALL cp_assert (pcol .LT. npcols, cp_fatal_level, cp_internal_error,&
            routineN, "Invalid process column.", failure=failure, error=error)
       row_distribution (cluster) = prow + 1
       col_distribution (cluster) = pcol + 1
       !
       cluster_price = cluster_prices(cluster_index)
       bin_price = bin_price + cluster_price
       CALL heap_reset_first (bin_heap, bin_price, heap_error)
       CALL cp_assert (.NOT. heap_error, cp_warning_level, cp_internal_error,&
            routineN, "Error setting price of top heap element.", error=error)
    ENDDO
    CALL heap_release (bin_heap, heap_error)
    CALL cp_assert (.NOT. heap_error, cp_warning_level, cp_internal_error,&
         routineN, "Error releasing heap.", error=error)
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
!> \param error ...
!> \par History
!> - Created 2010-11-11 Joost VandeVondele
! *****************************************************************************
  SUBROUTINE make_basic_spatial_distribution (pbc_scaled_coords, costs,&
       nprows, row_distribution, npcols, col_distribution, error)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: pbc_scaled_coords
    INTEGER, DIMENSION(:), INTENT(IN)        :: costs
    INTEGER, INTENT(IN)                      :: nprows
    INTEGER, DIMENSION(:), INTENT(OUT)       :: row_distribution
    INTEGER, INTENT(IN)                      :: npcols
    INTEGER, DIMENSION(:), INTENT(OUT)       :: col_distribution
    TYPE(cp_error_type), INTENT(INOUT)       :: error

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
    ALLOCATE(bin_costs(nbins),distribution(natoms))
    bin_costs=0

    CALL spatial_recurse(pbc_scaled_coords, costs, (/(iatom,iatom=1,natoms)/) , bin_costs, distribution,0)

    ! WRITE(*,*) "Final bin costs: ",bin_costs

    ! final row_distribution / col_distribution
    DO iatom=1,natoms
       row_distribution (iatom) = (distribution(iatom)-1)*pgrid_gcd / npcols +1
       col_distribution (iatom) = (distribution(iatom)-1)*pgrid_gcd / nprows +1
    ENDDO

    DEALLOCATE(bin_costs,distribution)

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
  RECURSIVE SUBROUTINE spatial_recurse(pbc_scaled_coords, costs, indices, bin_costs, distribution,level)
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
       ALLOCATE(bin_costs_sorted(nbins),permutation(nbins))
       bin_costs_sorted(:)=bin_costs
       CALL sort(bin_costs_sorted,nbins,permutation)
       ALLOCATE(atom_costs_sorted(natoms),atom_permutation(natoms))
       atom_costs_sorted(:)=costs
       CALL sort(atom_costs_sorted,natoms,atom_permutation)
       ibin=0
       ! WRITE(*,*) "Dealing with a new bunch of atoms "
       DO iatom=natoms,1,-1
          ibin=ibin+1
          ! WRITE(*,*) "atom",indices(atom_permutation(iatom)),"cost",atom_costs_sorted(iatom),&
          !            "bin",permutation(ibin),"its cost",bin_costs(permutation(ibin))
          ! WRITE(100,'(A,I0,3F12.6)') "A",permutation(ibin),pbc_scaled_coords(:,atom_permutation(iatom))
          bin_costs(permutation(ibin))=bin_costs(permutation(ibin))+atom_costs_sorted(iatom)
          distribution(indices(atom_permutation(iatom)))=permutation(ibin)
       ENDDO
       DEALLOCATE(bin_costs_sorted,permutation,atom_costs_sorted,atom_permutation)
    ELSE
       ! divide atoms in two subsets, sorting according to their coordinates, alternatively x,y,z
       ! recursively do this for both subsets
       ALLOCATE(coord(natoms),permutation(natoms))
       coord(:)=pbc_scaled_coords(MOD(level,3)+1,:)
       CALL sort(coord,natoms,permutation)
       CALL spatial_recurse(pbc_scaled_coords(:,permutation(1:nhalf)), costs(permutation(1:nhalf)),&
                            indices(permutation(1:nhalf)), bin_costs, distribution, level+1)
       CALL spatial_recurse(pbc_scaled_coords(:,permutation(nhalf+1:)), costs(permutation(nhalf+1:)),&
                            indices(permutation(nhalf+1:)), bin_costs, distribution, level+1)
       DEALLOCATE(coord,permutation)
    ENDIF

  END SUBROUTINE spatial_recurse

END MODULE distribution_methods
