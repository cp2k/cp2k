!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****s* cp2k/distribution_methods [1.0] *
!!
!!   NAME
!!     distribution_methods
!!
!!   FUNCTION
!!     Distribution methods for atoms, particles, or molecules
!!
!!   AUTHOR
!!     MK (22.08.2003)
!!
!!   MODIFICATION HISTORY
!!     - 1d-distribution of molecules and particles (Sep. 2003,MK)
!!     - 2d-distribution for Quickstep updated with molecules (Oct. 2003,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE distribution_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_array_i_utils,                ONLY: cp_1d_i_p_type
  USE cp_error_handling,               ONLY: cp_error_get_logger,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_logger_get_default_unit_nr,&
                                             cp_logger_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_create,&
                                             distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_create,&
                                             distribution_2d_type,&
                                             distribution_2d_write
  USE global_types,                    ONLY: LOW
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "distribution_methods"

! *** Public subroutines ***

  PUBLIC :: distribute_molecules_1d,&
            distribute_molecules_2d

! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****f* distribution_methods/distribute_molecules_1d [1.0] *
!!
!!   NAME
!!     distribute_molecules_1d
!!
!!   FUNCTION
!!     Distribute molecules and particles
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - particle_kind_set: particle (atomic) kind information
!!     - particle_set: particle information
!!     - local_particles: distribution of particles created by this routine
!!     - molecule_kind_set: molecule kind information
!!     - molecule_set: molecule information
!!     - local_molecules: distribution of molecules created by this routine
!!     - error: variable to control error logging, stopping, ...
!!              see module cp_error_handling
!!
!!   AUTHOR
!!     MK (Jun. 2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************

  SUBROUTINE distribute_molecules_1d(particle_kind_set,particle_set,&
                                     local_particles,&
                                     molecule_kind_set,molecule_set,&
                                     local_molecules, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: particle_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routine_name = "distribute_molecules_1d"

    INTEGER :: atom_a, group, iatom, imolecule, imolecule_kind, &
      imolecule_local, iparticle_kind, ipe, istat, jpe, kind_a, molecule_a, &
      mype, n, natom, nload, nmolecule, nmolecule_kind, nparticle_kind, npe, &
      nsgf, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nmolecule_local, &
                                                nparticle_local, work, &
                                                workload, workload_old
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    LOGICAL                                  :: ionode
    TYPE(cp_1d_i_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: local_molecule
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

!   ---------------------------------------------------------------------------

    logger => cp_error_get_logger(error)

    group = logger%para_env%group
    mype = logger%para_env%mepos + 1
    npe = logger%para_env%num_pe
    ionode= (logger%para_env%mepos == logger%para_env%source)

    IF (logger%print_keys%level > LOW) THEN
       output_unit = cp_logger_get_default_unit_nr(logger,LOCAL=.TRUE.)
    ELSE
       IF (ionode) output_unit = cp_logger_get_default_unit_nr(logger)
    ENDIF

    ALLOCATE (workload(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "workload",npe*int_size)
    workload(:) = 0

    ALLOCATE (workload_old(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "workload_old",npe*int_size)
    workload_old(:) = 0

    nmolecule_kind = SIZE(molecule_kind_set)

    ALLOCATE (nmolecule_local(nmolecule_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nmolecule_local",nmolecule_kind*int_size)
    nmolecule_local(:) = 0

    ALLOCATE (local_molecule(nmolecule_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "local_molecule",nmolecule_kind*int_size)

    nparticle_kind = SIZE(particle_kind_set)

    ALLOCATE (nparticle_local(nparticle_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nparticle_local",nparticle_kind*int_size)
    nparticle_local(:) = 0

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

!     *** Save current estimated work load ***

      workload_old(:) = workload(:)

!     *** Get the number of local molecules of the current molecule kind ***

      DO imolecule=1,nmolecule
        ipe = 1
        DO jpe=1,npe
          IF (workload(jpe) < workload(ipe)) ipe = jpe
        END DO
        workload(ipe) = workload(ipe) + nload
        IF (ipe == mype) THEN
          nmolecule_local(imolecule_kind) = nmolecule_local(imolecule_kind) + 1
        END IF
      END DO

!     *** Reset work load ***

      workload(:) = workload_old(:)

!     *** Distribute the molecules ***

      n = nmolecule_local(imolecule_kind)

      IF (n > 0) THEN
        ALLOCATE (local_molecule(imolecule_kind)%array(n),STAT=istat)
        IF (istat /= 0) THEN
           CALL stop_memory(routine_name,module_name,__LINE__,&
                            "local_molecule(imolecule_kind)%array",&
                            n*int_size)
        END IF
      ELSE
        NULLIFY (local_molecule(imolecule_kind)%array)
      END IF

      imolecule_local = 0

      DO imolecule=1,nmolecule
        ipe = 1
        DO jpe=1,npe
          IF (workload(jpe) < workload(ipe)) ipe = jpe
        END DO
        workload(ipe) = workload(ipe) + nload
        IF (ipe == mype) THEN
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
      END DO

    END DO


!   *** Create the local molecule structure ***

    CALL distribution_1d_create(local_molecules,&
                                n_el=nmolecule_local,&
                                para_env=logger%para_env)

!   *** Create the local particle structure ***

    CALL distribution_1d_create(local_particles,&
                                n_el=nparticle_local,&
                                para_env=logger%para_env)

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
        !  ALLOCATE (particle_set(atom_a)%lpi,STAT=istat)
        !  IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
        !                                   "particle_set(atom_a)%lpi",0)
        !  NULLIFY (particle_set(atom_a)%lpi%exclusion)
        END DO
      END DO

    END DO

!   *** Print distribution, if requested ***

    IF (logger%print_keys%distribution) THEN

!     *** Print molecule distribution ***

      ALLOCATE (work(npe),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "work",npe*int_size)
      work(:) = 0

      work(mype) = SUM(nmolecule_local)
      CALL mp_sum(work,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,&
               FMT="(/,/,T2,A,/,/,T3,A,/,/,(T4,I6,T24,I8))")&
          "DISTRIBUTION OF THE MOLECULES",&
          "Process   Number of molecules",&
          (ipe-1,work(ipe),ipe=1,npe)
        WRITE (UNIT=output_unit,FMT="(/,T7,A3,T24,I8)")&
          "Sum",SUM(work)
         CALL m_flush(output_unit)
      END IF

      CALL mp_sync(group)

      IF (logger%print_keys%level > LOW) THEN
        DO ipe=1,npe
          IF (ipe == mype) THEN
            WRITE (UNIT=output_unit,FMT="(/,T3,A)")&
                "Process   Kind   Local molecules (global indices)"
            DO imolecule_kind=1,nmolecule_kind
              IF (imolecule_kind == 1) THEN
                WRITE (UNIT=output_unit,FMT="(/,T4,I6,2X,I5,(T20,10I6))")&
                  ipe-1,imolecule_kind,&
                  (local_molecules%list(imolecule_kind)%array(imolecule),&
                   imolecule=1,nmolecule_local(imolecule_kind))
              ELSE
                WRITE (UNIT=output_unit,FMT="(T12,I5,(T20,10I6))")&
                  imolecule_kind,&
                  (local_molecules%list(imolecule_kind)%array(imolecule),&
                   imolecule=1,nmolecule_local(imolecule_kind))
              END IF
            END DO
          END IF
          CALL m_flush(output_unit)
          CALL mp_sync(group)
        END DO
      END IF

!     *** Print particle distribution ***

      work(:) = 0

      work(mype) = SUM(nparticle_local)
      CALL mp_sum(work,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,&
               FMT="(/,/,T2,A,/,/,T3,A,/,/,(T4,I6,T24,I8))")&
          "DISTRIBUTION OF THE PARTICLES",&
          "Process   Number of particles",&
          (ipe-1,work(ipe),ipe=1,npe)
        WRITE (UNIT=output_unit,FMT="(/,T7,A3,T24,I8)")&
          "Sum",SUM(work)
        CALL m_flush(output_unit)
      END IF

      CALL mp_sync(group)

      IF (logger%print_keys%level > LOW) THEN
        DO ipe=1,npe
          IF (ipe == mype) THEN
            WRITE (UNIT=output_unit,FMT="(/,T3,A)")&
                "Process   Kind   Local particles (global indices)"
            DO iparticle_kind=1,nparticle_kind
              IF (iparticle_kind == 1) THEN
                WRITE (UNIT=output_unit,FMT="(/,T4,I6,2X,I5,(T20,10I6))")&
                  ipe-1,iparticle_kind,&
                  (local_particles%list(iparticle_kind)%array(iatom),&
                   iatom=1,nparticle_local(iparticle_kind))
              ELSE
                WRITE (UNIT=output_unit,FMT="(T12,I5,(T20,10I6))")&
                  iparticle_kind,&
                  (local_particles%list(iparticle_kind)%array(iatom),&
                   iatom=1,nparticle_local(iparticle_kind))
              END IF
            END DO
          END IF
          CALL m_flush(output_unit)
          CALL mp_sync(group)
        END DO
      END IF

      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"work")

    END IF

!   *** Release work storage ***

    DEALLOCATE (workload,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "workload")

    DEALLOCATE (workload_old,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "workload_old")

    DEALLOCATE (nmolecule_local,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nmolecule_local")

    DEALLOCATE (nparticle_local,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nparticle_local")

    DO imolecule_kind=1,nmolecule_kind
      IF (ASSOCIATED(local_molecule(imolecule_kind)%array)) THEN
        DEALLOCATE (local_molecule(imolecule_kind)%array,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine_name,module_name,__LINE__,&
                           "local_molecule(imolecule_kind)%array")
        END IF
      END IF
    END DO
    DEALLOCATE (local_molecule,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "local_molecule")

  END SUBROUTINE distribute_molecules_1d

! *****************************************************************************
!!****f* distribution_methods/distribute_molecules_2d [1.0] *
!!
!!   NAME
!!     distribute_molecules_2d
!!
!!   FUNCTION
!!     Distributes the particle pairs creating a 2d distribution optimally
!!     suited for quickstep
!!
!!   NOTES
!!     Intermediate generation of a 2d distribution of the molecules, but
!!     only the corresponding particle (atomic) distribution is currently
!!     used. The 2d distribution of the molecules is deleted, but may easily
!!     be recovered (MK).
!!
!!   ARGUMENTS
!!     - distribution_2d: the distribution that will be created by this
!!                        method
!!     - atomic_kind_set: the atomic kind set to distribute
!!     - blacs_env: the parallel environement at the basis of the
!!                  distribution
!!     - error: variable to control error logging, stopping,...
!!              see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi (Feb. 2003)
!!
!!   MODIFICATION HISTORY
!!     - local_rows & cols blocksize optimizations (Aug. 2003,MK)
!!     - cleanup of distribution_2d (Sep. 2003,fawzi)
!!     - update for molecules (Oct. 2003,MK)
!!
!!*** **********************************************************************

  SUBROUTINE distribute_molecules_2d(particle_kind_set,particle_set,&
                                     molecule_kind_set,molecule_set,&
                                     distribution_2d,blacs_env,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: particle_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routine_name = "distribute_molecules_2d"

    INTEGER :: atom_a, group, iatom, imolecule, imolecule_kind, &
      imolecule_local, iparticle_kind, ipcol, iprow, istat, jpcol, jprow, &
      kind_a, molecule_a, mypcol, myprow, n, natom, nload, nmolecule, &
      nmolecule_kind, nparticle_kind, npcol, nprow, nsgf, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:) :: nmolecule_local_col, &
      nmolecule_local_row, nparticle_local_col, nparticle_local_row, work, &
      workload_col, workload_col_old, workload_row, workload_row_old
    INTEGER, DIMENSION(:), POINTER           :: col_distribution, &
                                                molecule_list, &
                                                row_distribution
    LOGICAL                                  :: ionode
    TYPE(cp_1d_i_p_type), DIMENSION(:), &
      POINTER                                :: local_molecule_col, &
                                                local_molecule_row, &
                                                local_particle_col, &
                                                local_particle_row
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

!   ---------------------------------------------------------------------------

    logger => cp_error_get_logger(error)

    group = blacs_env%para_env%group
    ionode= (logger%para_env%mepos == logger%para_env%source)
    IF (ionode) THEN
        output_unit = cp_logger_get_default_unit_nr(logger)
    ENDIF

    myprow = blacs_env%mepos(1) + 1
    mypcol = blacs_env%mepos(2) + 1
    nprow = blacs_env%num_pe(1)
    npcol = blacs_env%num_pe(2)

    ALLOCATE (workload_row(nprow),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "workload_row",nprow*int_size)
    workload_row(:) = 0

    ALLOCATE (workload_row_old(nprow),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "workload_row_old",nprow*int_size)

    ALLOCATE (workload_col(npcol),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "workload_col",npcol*int_size)
    workload_col(:) = 0

    ALLOCATE (workload_col_old(npcol),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "workload_col_old",npcol*int_size)

    nmolecule_kind = SIZE(molecule_kind_set)

    ALLOCATE (local_molecule_row(nmolecule_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "local_molecule_row",&
                                     nmolecule_kind*int_size)

    ALLOCATE (local_molecule_col(nmolecule_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "local_molecule_col",&
                                     nmolecule_kind*int_size)

    ALLOCATE (nmolecule_local_row(nmolecule_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nmolecule_local_row",&
                                     nmolecule_kind*int_size)
    nmolecule_local_row(:) = 0

    ALLOCATE (nmolecule_local_col(nmolecule_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nmolecule_local_col",&
                                     nmolecule_kind*int_size)
    nmolecule_local_col(:) = 0

    nparticle_kind = SIZE(particle_kind_set)

    ALLOCATE (local_particle_row(nparticle_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "local_particle_row",&
                                     nparticle_kind*int_size)

    ALLOCATE (local_particle_col(nparticle_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "local_particle_col",&
                                     nparticle_kind*int_size)

    ALLOCATE (nparticle_local_row(nparticle_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nparticle_local_row",&
                                     nparticle_kind*int_size)
    nparticle_local_row(:) = 0

    ALLOCATE (nparticle_local_col(nparticle_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nparticle_local_col",&
                                     nparticle_kind*int_size)
    nparticle_local_col(:) = 0

    CALL get_atomic_kind_set(atomic_kind_set=particle_kind_set,natom=natom)

    ALLOCATE (row_distribution(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "row_distribution",natom*int_size)
    row_distribution(:) = 0

    ALLOCATE (col_distribution(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "col_distribution",natom*int_size)
    col_distribution(:) = 0

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

!     *** Save current estimated work load ***

      workload_row_old(:) = workload_row(:)
      workload_col_old(:) = workload_col(:)

!     *** Get the number of local molecules for each process row ***

      DO imolecule=1,nmolecule
        iprow = 1
        DO jprow=1,nprow
          IF (workload_row(jprow) < workload_row(iprow)) iprow = jprow
        END DO
        workload_row(iprow) = workload_row(iprow) + nload
        IF (iprow == myprow) THEN
          nmolecule_local_row(imolecule_kind) =&
            nmolecule_local_row(imolecule_kind) + 1
        END IF
      END DO

!     *** Get the number of local molecules for each process column ***

      DO imolecule=1,nmolecule
        ipcol = 1
        DO jpcol=1,npcol
          IF (workload_col(jpcol) < workload_col(ipcol)) ipcol = jpcol
        END DO
        workload_col(ipcol) = workload_col(ipcol) + nload
        IF (ipcol == mypcol) THEN
          nmolecule_local_col(imolecule_kind) =&
            nmolecule_local_col(imolecule_kind) + 1
        END IF
      END DO

!     *** Reset work load ***

      workload_row(:) = workload_row_old(:)
      workload_col(:) = workload_col_old(:)

!     *** Distribute the molecules ***

!     *** Block rows ***

      n = nmolecule_local_row(imolecule_kind)

      !JVDV IF (n > 0) THEN
        ALLOCATE (local_molecule_row(imolecule_kind)%array(n),STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine_name,module_name,__LINE__,&
                           "local_molecule_row(imolecule_kind)%array",&
                           n*int_size)
        END IF
      !JVDV ELSE
      !JVDV NULLIFY (local_molecule_row(imolecule_kind)%array)
      !JVDV END IF

      imolecule_local = 0

      DO imolecule=1,nmolecule
        iprow = 1
        DO jprow=1,nprow
          IF (workload_row(jprow) < workload_row(iprow)) iprow = jprow
        END DO
        workload_row(iprow) = workload_row(iprow) + nload
        molecule_a = molecule_list(imolecule)
        DO iatom=1,natom
          atom_a = molecule_set(molecule_a)%first_atom + iatom - 1
          row_distribution(atom_a) = iprow - 1 ! generate blacs rows
          IF (iprow == myprow) THEN
            IF (iatom == 1) THEN
              imolecule_local = imolecule_local + 1
              local_molecule_row(imolecule_kind)%array(imolecule_local) = molecule_a
            END IF
            CALL get_atomic_kind(atomic_kind=particle_set(atom_a)%atomic_kind,&
                                 kind_number=kind_a)
            nparticle_local_row(kind_a) = nparticle_local_row(kind_a) + 1
          END IF
        END DO
      END DO

!     *** Block columns ***

      n = nmolecule_local_col(imolecule_kind)

      !JVDV IF (n > 0) THEN
        ALLOCATE (local_molecule_col(imolecule_kind)%array(n),STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine_name,module_name,__LINE__,&
                           "local_molecule_col(imolecule_kind)%array",&
                           n*int_size)
        END IF
      !JVDV ELSE
      !JVDV   NULLIFY (local_molecule_col(imolecule_kind)%array)
      !JVDV END IF

      imolecule_local = 0

      DO imolecule=1,nmolecule
        ipcol = 1
        DO jpcol=1,npcol
          IF (workload_col(jpcol) < workload_col(ipcol)) ipcol = jpcol
        END DO
        workload_col(ipcol) = workload_col(ipcol) + nload
        molecule_a = molecule_list(imolecule)
        DO iatom=1,natom
          atom_a = molecule_set(molecule_a)%first_atom + iatom - 1
          col_distribution(atom_a) = ipcol - 1 ! generate blacs cols
          IF (ipcol == mypcol) THEN
            IF (iatom == 1) THEN
              imolecule_local = imolecule_local + 1
              local_molecule_col(imolecule_kind)%array(imolecule_local) = molecule_a
            END IF
            CALL get_atomic_kind(atomic_kind=particle_set(atom_a)%atomic_kind,&
                                 kind_number=kind_a)
            nparticle_local_col(kind_a) = nparticle_local_col(kind_a) + 1
          END IF
        END DO
      END DO

    END DO

!   *** Build the particle distribution w.r.t. the molecule distribution ***

    DO iparticle_kind=1,nparticle_kind

      n = nparticle_local_row(iparticle_kind)

      !JVDV IF (n > 0) THEN
        ALLOCATE (local_particle_row(iparticle_kind)%array(n),STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine_name,module_name,__LINE__,&
                           "local_particle_row(iparticle_kind)%array",&
                           n*int_size)
        END IF
      !JVDV ELSE
      !JVDV   NULLIFY (local_particle_row(iparticle_kind)%array)
      !JVDV END IF

      n = nparticle_local_col(iparticle_kind)

      !JVDV IF (n > 0) THEN
        ALLOCATE (local_particle_col(iparticle_kind)%array(n),STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine_name,module_name,__LINE__,&
                           "local_particle_col(iparticle_kind)%array",&
                           n*int_size)
        END IF
      !JVDV ELSE
      !JVDV   NULLIFY (local_particle_col(iparticle_kind)%array)
      !JVDV END IF

    END DO

    nparticle_local_row(:) = 0
    nparticle_local_col(:) = 0

    DO imolecule_kind=1,nmolecule_kind

      molecule_kind => molecule_kind_set(imolecule_kind)

      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             natom=natom)

      DO imolecule=1,nmolecule_local_row(imolecule_kind)
        molecule_a = local_molecule_row(imolecule_kind)%array(imolecule)
        DO iatom=1,natom
          atom_a = molecule_set(molecule_a)%first_atom + iatom - 1
          CALL get_atomic_kind(atomic_kind=particle_set(atom_a)%atomic_kind,&
                               kind_number=kind_a)
          nparticle_local_row(kind_a) = nparticle_local_row(kind_a) + 1
          local_particle_row(kind_a)%array(nparticle_local_row(kind_a)) = atom_a
        END DO
      END DO

      DO imolecule=1,nmolecule_local_col(imolecule_kind)
        molecule_a = local_molecule_col(imolecule_kind)%array(imolecule)
        DO iatom=1,natom
          atom_a = molecule_set(molecule_a)%first_atom + iatom - 1
          CALL get_atomic_kind(atomic_kind=particle_set(atom_a)%atomic_kind,&
                               kind_number=kind_a)
          nparticle_local_col(kind_a) = nparticle_local_col(kind_a) + 1
          local_particle_col(kind_a)%array(nparticle_local_col(kind_a)) = atom_a
        END DO
      END DO

    END DO

!   *** Generate the 2d distribution structure ***

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

    IF (logger%print_keys%distribution) THEN

!     *** Print row distribution ***

      ALLOCATE (work(nprow),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,&
                         "work",nprow*int_size)
      END IF
      work(:) = 0

      IF (mypcol == 1) work(myprow) = SUM(distribution_2d%n_local_rows)

      CALL mp_sum(work,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,&
               FMT="(/,/,T2,A,/,/,T3,A,/,/,(T4,I10,T26,I10,T53,I10))")&
          "DISTRIBUTION OF THE PARTICLES (ROWS)",&
          "Process row   Number of particles      Number of matrix rows",&
          (iprow-1,work(iprow),workload_row(iprow),iprow=1,nprow)
        WRITE (UNIT=output_unit,FMT="(/,T11,A3,T26,I10,T53,I10)")&
          "Sum",SUM(work),SUM(workload_row)
        CALL m_flush(output_unit)
      END IF

      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF

!     *** Print column distribution ***

      ALLOCATE (work(npcol),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,&
                         "work",npcol*int_size)
      END IF
      work(:) = 0

      IF (myprow == 1) work(mypcol) = SUM(distribution_2d%n_local_cols)

      CALL mp_sum(work,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,&
               FMT="(/,/,T2,A,/,/,T3,A,/,/,(T4,I10,T26,I10,T53,I10))")&
          "DISTRIBUTION OF THE PARTICLES (COLUMNS)",&
          "Process col   Number of particles   Number of matrix columns",&
          (ipcol-1,work(ipcol),workload_col(ipcol),ipcol=1,npcol)
        WRITE (UNIT=output_unit,FMT="(/,T11,A3,T26,I10,T53,I10)")&
          "Sum",SUM(work),SUM(workload_col)
        CALL m_flush(output_unit)
      END IF


      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF

      IF (logger%print_keys%level > LOW) THEN
        output_unit = cp_logger_get_default_unit_nr(logger,local=.TRUE.)
        CALL distribution_2d_write(distribution_2d,&
                                   unit_nr=output_unit,&
                                   local=.TRUE.,&
                                   long_description=.TRUE.,&
                                   error=error)
      END IF

    END IF

!   *** Release work storage ***

    DEALLOCATE (workload_row,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                    "workload_row")

    DEALLOCATE (workload_row_old,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                    "workload_row_old")

    DEALLOCATE (workload_col,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                    "workload_col")

    DEALLOCATE (workload_col_old,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                    "workload_col_old")

    DO imolecule_kind=1,nmolecule_kind
      IF (ASSOCIATED(local_molecule_row(imolecule_kind)%array)) THEN
        DEALLOCATE (local_molecule_row(imolecule_kind)%array,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine_name,module_name,__LINE__,&
                           "local_molecule_row(imolecule_kind)%array")
        END IF
      END IF
    END DO
    DEALLOCATE (local_molecule_row,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "local_molecule_row")

    DO imolecule_kind=1,nmolecule_kind
      IF (ASSOCIATED(local_molecule_col(imolecule_kind)%array)) THEN
        DEALLOCATE (local_molecule_col(imolecule_kind)%array,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine_name,module_name,__LINE__,&
                           "local_molecule_col(imolecule_kind)%array")
        END IF
      END IF
    END DO
    DEALLOCATE (local_molecule_col,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "local_molecule_col")

    DEALLOCATE (nmolecule_local_row,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nmolecule_local_row")

    DEALLOCATE (nmolecule_local_col,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nmolecule_local_col")

    DEALLOCATE (nparticle_local_row,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nparticle_local_row")

    DEALLOCATE (nparticle_local_col,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nparticle_local_col")

  END SUBROUTINE distribute_molecules_2d

! *****************************************************************************

END MODULE distribution_methods
