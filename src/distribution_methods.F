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
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE distribution_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_create,&
                                             distribution_1d_type
  USE global_types,                    ONLY: LOW,&
                                             global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             wp => dp,&
                                             wp_size => dp_size
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "distribution_methods"

! *** Public subroutines ***

  PUBLIC :: distribute_molecules

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE distribute_molecules(particle_kind_set,particle_set,&
                                  local_particles,&
                                  molecule_kind_set,molecule_set,&
                                  local_molecules,globenv)

    TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: particle_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), POINTER :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), POINTER      :: molecule_set
    TYPE(distribution_1d_type), POINTER             :: local_molecules,&
                                                       local_particles
    TYPE(global_environment_type), INTENT(IN)       :: globenv

!   *** Local types ***

    TYPE local_list_type
      INTEGER, DIMENSION(:), POINTER :: list
    END TYPE local_list_type

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "distribute_molecules"

!   *** Local variables ***

    TYPE(molecule_kind_type), POINTER :: molecule_kind

    INTEGER :: atom_a,group,i,iatom,imolecule,imolecule_kind,imolecule_local,&
               iparticle_kind,ipe,istat,j,kind_a,molecule_a,mype,n,natom,&
               nmolecule,nmolecule_kind,nparticle_kind,npe,output_unit
    LOGICAL :: ionode

    TYPE(local_list_type), DIMENSION(:), ALLOCATABLE :: local_molecule,&
                                                        local_particle
    INTEGER, DIMENSION(:), ALLOCATABLE               :: natom_local,&
                                                        natom_local_old,&
                                                        nmolecule_local,&
                                                        nparticle_local,&
                                                        work
    INTEGER, DIMENSION(:), POINTER                   :: molecule_list

!   ---------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    mype = globenv%mepos + 1
    npe = globenv%num_pe
    output_unit = globenv%scr

    ALLOCATE (natom_local(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "natom_local",npe*int_size)
    natom_local(:) = 0

    ALLOCATE (natom_local_old(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "natom_local_old",npe*int_size)
    natom_local_old(:) = 0

    nmolecule_kind = SIZE(molecule_kind_set)

    ALLOCATE (nmolecule_local(nmolecule_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nmolecule_local",nmolecule_kind*int_size)
    nmolecule_local(:) = 0

    nparticle_kind = SIZE(particle_kind_set)

    ALLOCATE (nparticle_local(nparticle_kind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nparticle_local",nparticle_kind*int_size)
    nparticle_local(:) = 0

    DO imolecule_kind=1,nmolecule_kind

      molecule_kind => molecule_kind_set(imolecule_kind)

!     *** Get the number of molecules and the number of ***
!     *** atoms in each molecule of that molecular kind ***

      NULLIFY (molecule_list)

      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             molecule_list=molecule_list,&
                             natom=natom)

      nmolecule = SIZE(molecule_list)
PRINT*,imolecule_kind,nmolecule
!     *** Save current estimated work load ***

      natom_local_old(:) = natom_local(:)

!     *** Get the number of local molecules of the current molecule kind ***

      DO imolecule=1,nmolecule
        i = 1
        DO j=1,npe
          IF (natom_local(j) < natom_local(i)) i = j
        END DO
        natom_local(i) = natom_local(i) + natom
        IF (i == mype) THEN
          nmolecule_local(imolecule_kind) = nmolecule_local(imolecule_kind) + 1
        END IF
      END DO

!     *** Reset work load ***

      natom_local(:) = natom_local_old(:)

!     *** Distribute the molecules ***

      n = nmolecule_local(imolecule_kind)

      IF (n > 0) THEN
        ALLOCATE (local_molecule(imolecule_kind)%list(n),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "local_molecule(imolecule_kind)%list",&
                                         n*int_size)
      ELSE
        NULLIFY (local_molecule(imolecule_kind)%list)
      END IF

      imolecule_local = 0

      DO imolecule=1,nmolecule
        i = 1
        DO j=1,npe
          IF (natom_local(j) < natom_local(i)) i = j
        END DO
        natom_local(i) = natom_local(i) + natom
        IF (i == mype) THEN
          imolecule_local = imolecule_local + 1
          molecule_a = molecule_list(imolecule)
          local_molecule(imolecule_kind)%list(imolecule_local) = molecule_a
          DO iatom=1,natom
            atom_a = molecule_set(molecule_a)%first_atom + iatom - 1
            kind_a = particle_set(atom_a)%kind
            nparticle_local(kind_a) = nparticle_local(kind_a) + 1
          END DO
        END IF
      END DO

    END DO

!   *** Create the local molecule structure ***

    CALL distribution_1d_create(local_molecules,&
                      n_el=nmolecule_local,&
                      para_env=globenv%para_env)

!   *** Create the local particle structure ***

    DO iparticle_kind=1,nparticle_kind
      n = nparticle_local(iparticle_kind)
      IF (n > 0) THEN
        ALLOCATE (local_particle(iparticle_kind)%list(n),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "local_particle(iparticle_kind)%list",&
                                         n*int_size)
      ELSE
        NULLIFY (local_particle(iparticle_kind)%list)
      END IF
    END DO

    CALL distribution_1d_create(local_particles,&
                      n_el=nparticle_local,&
                      para_env=globenv%para_env)

!   *** Store the generated local molecule and particle distributions ***

    nparticle_local(:) = 0

    DO imolecule_kind=1,nmolecule_kind

      local_molecules%list(imolecule_kind)%array(:) =&
        local_molecule(imolecule_kind)%list(:)

      molecule_kind => molecule_kind_set(imolecule_kind)

      NULLIFY (molecule_list)

      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             molecule_list=molecule_list,&
                             natom=natom)

      nmolecule = SIZE(molecule_list)

      DO imolecule=1,nmolecule
        molecule_a = molecule_list(imolecule)
        DO iatom=1,natom
          atom_a = molecule_set(molecule_a)%first_atom + iatom - 1
          kind_a = particle_set(atom_a)%kind
          nparticle_local(kind_a) = nparticle_local(kind_a) + 1
          local_particles%list(kind_a)%array(nparticle_local(kind_a)) = atom_a
        END DO
      END DO

    END DO

!   *** Print distribution, if requested ***

    IF (globenv%print%distribution) THEN

!     *** Print molecule distribution ***

      ALLOCATE (work(npe),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "work",npe*int_size)

      work(mype) = SUM(nmolecule_local)

      CALL mp_sum(work,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,&
               FMT="(/,/,T2,A,/,/,T3,A,/,/,(T4,I10,T26,I10))")&
          "DISTRIBUTION OF THE MOLECULES",&
          "Process   Number of molecules",&
          (ipe-1,work(ipe),ipe=1,npe)
        WRITE (UNIT=output_unit,FMT="(/,T11,A3,T26,I10)")&
          "Sum",SUM(work)
      END IF

      CALL m_flush(output_unit)
      CALL mp_sync(group)

      IF (globenv%print%level > LOW) THEN
        DO ipe=1,npe
          IF (ipe == mype) THEN
            IF (ipe == 1) THEN
              WRITE (UNIT=output_unit,FMT="(/,T3,A)")&
                "Process   Kind   Local molecules (global indices)"
            END IF
            DO imolecule_kind=1,nmolecule_kind
              IF (imolecule_kind == 1) THEN
                WRITE (UNIT=output_unit,FMT="(/,T4,I10,2X,I5,(T24,10I6))")&
                  ipe-1,imolecule_kind,&
                  (local_molecule(imolecule_kind)%list(imolecule),&
                   imolecule=1,nmolecule_local(imolecule_kind))
              ELSE
                WRITE (UNIT=output_unit,FMT="(T16,I5,(T24,10I6))")&
                  imolecule_kind,&
                  (local_molecule(imolecule_kind)%list(imolecule),&
                   imolecule=1,nmolecule_local(imolecule_kind))
              END IF
            END DO
          END IF
          CALL m_flush(output_unit)
          CALL mp_sync(group)
        END DO
      END IF

!     *** Print particle distribution ***

      work(mype) = SUM(nparticle_local)

      CALL mp_sum(work,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,&
               FMT="(/,/,T2,A,/,/,T3,A,/,/,(T4,I10,T26,I10))")&
          "DISTRIBUTION OF THE PARTICLES",&
          "Process   Number of particles",&
          (ipe-1,work(ipe),ipe=1,npe)
        WRITE (UNIT=output_unit,FMT="(/,T11,A3,T26,I10)")&
          "Sum",SUM(work)
      END IF

      CALL m_flush(output_unit)
      CALL mp_sync(group)

      IF (globenv%print%level > LOW) THEN
        DO ipe=1,npe
          IF (ipe == mype) THEN
            IF (ipe == 1) THEN
              WRITE (UNIT=output_unit,FMT="(/,T3,A)")&
                "Process   Kind   Local particles (global indices)"
            END IF
            DO iparticle_kind=1,nparticle_kind
              IF (iparticle_kind == 1) THEN
                WRITE (UNIT=output_unit,FMT="(/,T4,I10,2X,I5,(T24,10I6))")&
                  ipe-1,iparticle_kind,&
                  (local_particle(iparticle_kind)%list(iatom),&
                   iatom=1,nparticle_local(iparticle_kind))
              ELSE
                WRITE (UNIT=output_unit,FMT="(T16,I5,(T24,10I6))")&
                  iparticle_kind,&
                  (local_particle(iparticle_kind)%list(iatom),&
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

    DEALLOCATE (natom_local,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "natom_local")

    DEALLOCATE (natom_local_old,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "natom_local_old")

    DEALLOCATE (nmolecule_local,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nmolecule_local")

    DEALLOCATE (nparticle_local,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nparticle_local")

    DO imolecule_kind=1,nmolecule_kind
      IF (ASSOCIATED(local_molecule(imolecule_kind)%list)) THEN
        DEALLOCATE (local_molecule(imolecule_kind)%list,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "local_molecule(imolecule_kind)%list")
      END IF
    END DO

    DEALLOCATE (local_molecule,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "local_molecule")

    DO iparticle_kind=1,nparticle_kind
      IF (ASSOCIATED(local_particle(iparticle_kind)%list)) THEN
        DEALLOCATE (local_particle(iparticle_kind)%list,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "local_particle(iparticle_kind)%list")
      END IF
    END DO

    DEALLOCATE (local_particle,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "local_particle")

  END SUBROUTINE distribute_molecules

! *****************************************************************************

END MODULE distribution_methods
