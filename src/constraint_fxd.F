!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      none
! *****************************************************************************
MODULE constraint_fxd
  
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: get_atomic_kind_set
  USE colvar_types,                    ONLY: colvar_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_constants,                 ONLY: use_perd_x,&
                                             use_perd_xy,&
                                             use_perd_xyz,&
                                             use_perd_xz,&
                                             use_perd_y,&
                                             use_perd_yz,&
                                             use_perd_z
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: fixd_constraint_type,&
                                             get_molecule_kind,&
                                             local_fixd_constraint_type,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: local_g3x3_constraint_type,&
                                             local_g4x6_constraint_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fix_atom_control,&
            check_fixed_atom_cns_g3x3,&
            check_fixed_atom_cns_g4x6,&
            check_fixed_atom_cns_colv,&
            create_local_fixd_list,&
            release_local_fixd_list

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'constraint_fxd'

CONTAINS

! *****************************************************************************
!> \brief allows for fix atom constraints
!> \par History
!>      - optionally apply fix atom constraint to random forces (Langevin)
!>        (04.10.206,MK)
! *****************************************************************************
  SUBROUTINE fix_atom_control( force_env, error, w)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL :: w

    CHARACTER(len=*), PARAMETER :: routineN = 'fix_atom_control', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ifixd, ii, ikind, iparticle, iparticle_local, &
      my_atm_fixed, natom, ncore, nfixed_atoms, nkind, nparticle, &
      nparticle_local, nshell, shell_index, stat
    LOGICAL                                  :: failure, shell_present
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: force
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_fixd_constraint_type), &
      POINTER                                :: lfixd_list(:)
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set

    failure = .FALSE.

    CALL timeset(routineN,handle)

    NULLIFY (atomic_kinds)
    NULLIFY (core_particles)
    NULLIFY (particles)
    NULLIFY (shell_particles)
    shell_present = .FALSE.

    IF (.NOT.failure) THEN
       NULLIFY(lfixd_list)
       CALL force_env_get(force_env=force_env,&
                          subsys=subsys,&
                          error=error)
       CALL cp_subsys_get(subsys=subsys,&
                          atomic_kinds=atomic_kinds,&
                          core_particles=core_particles,&
                          local_particles=local_particles,&
                          molecule_kinds_new=molecule_kinds,&
                          natom=natom,&
                          ncore=ncore,&
                          nshell=nshell,&
                          particles=particles,&
                          shell_particles=shell_particles,&
                          error=error)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,&
                                shell_present=shell_present)

       particle_set => particles%els
       CPPostcondition((SIZE(particle_set) == natom),cp_failure_level,routineP,error,failure)
       IF (shell_present) THEN
          core_particle_set => core_particles%els
          CPPostcondition((SIZE(core_particle_set) == ncore),cp_failure_level,routineP,error,failure)
          shell_particle_set => shell_particles%els
          CPPostcondition((SIZE(shell_particle_set) == nshell),cp_failure_level,routineP,error,failure)
       END IF
       nparticle = natom + nshell
       molecule_kind_set => molecule_kinds%els

       nkind = molecule_kinds%n_els
       my_atm_fixed = 0
       DO ikind=1,nkind
          molecule_kind => molecule_kind_set(ikind)
          CALL get_molecule_kind(molecule_kind,nfixd=nfixed_atoms)
          my_atm_fixed = my_atm_fixed + nfixed_atoms
       END DO

       IF (my_atm_fixed /= 0) THEN
          IF (.NOT.PRESENT(w)) THEN
             ! Allocate scratch array
             ALLOCATE (force(3,nparticle),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             force(:,:) = 0.0_dp
             DO i=1,SIZE(local_particles%n_el)
                nparticle_local = local_particles%n_el(i)
                DO iparticle_local=1,nparticle_local
                   iparticle = local_particles%list(i)%array(iparticle_local)
                   shell_index = particle_set(iparticle)%shell_index
                   IF (shell_index == 0) THEN
                      force(:,iparticle) = particle_set(iparticle)%f(:)
                   ELSE
                      force(:,iparticle) = core_particle_set(shell_index)%f(:)
                      force(:,natom+shell_index) = shell_particle_set(shell_index)%f(:)
                   END IF
                END DO
             END DO
          END IF

          ! Create the list of locally fixed atoms
          CALL create_local_fixd_list(lfixd_list,nkind,molecule_kind_set,local_particles,error)

          ! Apply fixed atom constraint
          DO ifixd=1,SIZE(lfixd_list)
             ikind = lfixd_list(ifixd)%ikind
             ii = lfixd_list(ifixd)%ifixd_index
             molecule_kind => molecule_kind_set(ikind)
             CALL get_molecule_kind(molecule_kind,fixd_list=fixd_list)
             IF (.NOT.fixd_list(ii)%restraint%active) THEN
                iparticle = fixd_list(ii)%fixd
                shell_index = particle_set(iparticle)%shell_index
                ! Select constraint type
                IF (PRESENT(w)) THEN
                   SELECT CASE(fixd_list(ii)%itype)
                   CASE (use_perd_x)
                      w(1,iparticle)     = 0.0_dp
                   CASE (use_perd_y)
                      w(2,iparticle)     = 0.0_dp
                   CASE (use_perd_z)
                      w(3,iparticle)     = 0.0_dp
                   CASE (use_perd_xy)
                      w(1,iparticle)     = 0.0_dp
                      w(2,iparticle)     = 0.0_dp
                   CASE (use_perd_xz)
                      w(1,iparticle)     = 0.0_dp
                      w(3,iparticle)     = 0.0_dp
                   CASE (use_perd_yz)
                      w(2,iparticle)     = 0.0_dp
                      w(3,iparticle)     = 0.0_dp
                   CASE (use_perd_xyz)
                      w(:,iparticle)     = 0.0_dp
                   END SELECT
                ELSE
                   SELECT CASE(fixd_list(ii)%itype)
                   CASE (use_perd_x)
                      force(1,iparticle) = 0.0_dp
                      IF (shell_index /= 0) THEN
                         force(1,natom+shell_index) = 0.0_dp
                      END IF
                   CASE (use_perd_y)
                      force(2,iparticle) = 0.0_dp
                      IF (shell_index /= 0) THEN
                         force(2,natom+shell_index) = 0.0_dp
                      END IF
                   CASE (use_perd_z)
                      force(3,iparticle) = 0.0_dp
                      IF (shell_index /= 0) THEN
                         force(3,natom+shell_index) = 0.0_dp
                      END IF
                   CASE (use_perd_xy)
                      force(1,iparticle) = 0.0_dp
                      force(2,iparticle) = 0.0_dp
                      IF (shell_index /= 0) THEN
                         force(1,natom+shell_index) = 0.0_dp
                         force(2,natom+shell_index) = 0.0_dp
                      END IF
                   CASE (use_perd_xz)
                      force(1,iparticle) = 0.0_dp
                      force(3,iparticle) = 0.0_dp
                      IF (shell_index /= 0) THEN
                         force(1,natom+shell_index) = 0.0_dp
                         force(3,natom+shell_index) = 0.0_dp
                      END IF
                   CASE (use_perd_yz)
                      force(2,iparticle) = 0.0_dp
                      force(3,iparticle) = 0.0_dp
                      IF (shell_index /= 0) THEN
                         force(2,natom+shell_index) = 0.0_dp
                         force(3,natom+shell_index) = 0.0_dp
                      END IF
                   CASE (use_perd_xyz)
                      force(:,iparticle) = 0.0_dp
                      IF (shell_index /= 0) THEN
                         force(:,natom+shell_index) = 0.0_dp
                      END IF
                   END SELECT
                END IF
             END IF
          END DO
          CALL release_local_fixd_list(lfixd_list,error=error)

          IF (.NOT.PRESENT(w)) THEN
             CALL mp_sum(force,force_env%para_env%group)
             DO iparticle=1,natom
                shell_index = particle_set(iparticle)%shell_index
                IF (shell_index == 0) THEN
                   particle_set(iparticle)%f(:) = force(:,iparticle)
                ELSE
                   core_particle_set(shell_index)%f(:) = force(:,iparticle)
                   shell_particle_set(shell_index)%f(:) = force(:,natom+shell_index)
                END IF
             END DO
             DEALLOCATE (force,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE fix_atom_control

! *****************************************************************************
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE check_fixed_atom_cns_g3x3(imass1, imass2, imass3,&
       index_a, index_b, index_c, fixd_list, lg3x3)
    REAL(KIND=dp), INTENT(INOUT)             :: imass1, imass2, imass3
    INTEGER, INTENT(IN)                      :: index_a, index_b, index_c
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_g3x3_constraint_type)         :: lg3x3

    INTEGER                                  :: i

    IF (lg3x3%init) THEN
       imass1 = lg3x3%imass1
       imass2 = lg3x3%imass2
       imass3 = lg3x3%imass3
    ELSE
       IF (ASSOCIATED(fixd_list)) THEN
          IF (SIZE(fixd_list)>0) THEN
             DO i = 1, SIZE(fixd_list)
                IF (fixd_list(i)%fixd==index_a) THEN
                   IF (fixd_list(i)%itype/=use_perd_xyz) CYCLE
                   IF (.NOT.fixd_list(i)%restraint%active) imass1 = 0.0_dp
                   EXIT
                END IF
             END DO
             DO i = 1, SIZE(fixd_list)
                IF (fixd_list(i)%fixd==index_b) THEN
                   IF (fixd_list(i)%itype/=use_perd_xyz) CYCLE
                   IF (.NOT.fixd_list(i)%restraint%active) imass2 = 0.0_dp
                   EXIT
                END IF
             END DO
             DO i = 1, SIZE(fixd_list)
                IF (fixd_list(i)%fixd==index_c) THEN
                   IF (fixd_list(i)%itype/=use_perd_xyz) CYCLE
                   IF (.NOT.fixd_list(i)%restraint%active) imass3 = 0.0_dp
                   EXIT
                END IF
             END DO
          END IF
       END IF
       lg3x3 % imass1 = imass1
       lg3x3 % imass2 = imass2
       lg3x3 % imass3 = imass3
       lg3x3 % init = .TRUE.
    END IF
  END SUBROUTINE check_fixed_atom_cns_g3x3

! *****************************************************************************
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE check_fixed_atom_cns_g4x6(imass1, imass2, imass3, imass4,&
       index_a, index_b, index_c, index_d, fixd_list, lg4x6)
    REAL(KIND=dp), INTENT(INOUT)             :: imass1, imass2, imass3, imass4
    INTEGER, INTENT(IN)                      :: index_a, index_b, index_c, &
                                                index_d
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_g4x6_constraint_type)         :: lg4x6

    INTEGER                                  :: i

    IF (lg4x6%init) THEN
       imass1 = lg4x6%imass1
       imass2 = lg4x6%imass2
       imass3 = lg4x6%imass3
       imass4 = lg4x6%imass4
    ELSE
       IF (ASSOCIATED(fixd_list)) THEN
          IF (SIZE(fixd_list)>0) THEN
             DO i = 1, SIZE(fixd_list)
                IF (fixd_list(i)%fixd==index_a) THEN
                   IF (fixd_list(i)%itype/=use_perd_xyz) CYCLE
                   IF (.NOT.fixd_list(i)%restraint%active) imass1 = 0.0_dp
                   EXIT
                END IF
             END DO
             DO i = 1, SIZE(fixd_list)
                IF (fixd_list(i)%fixd==index_b) THEN
                   IF (fixd_list(i)%itype/=use_perd_xyz) CYCLE
                   IF (.NOT.fixd_list(i)%restraint%active) imass2 = 0.0_dp
                   EXIT
                END IF
             END DO
             DO i = 1, SIZE(fixd_list)
                IF (fixd_list(i)%fixd==index_c) THEN
                   IF (fixd_list(i)%itype/=use_perd_xyz) CYCLE
                   IF (.NOT.fixd_list(i)%restraint%active) imass3 = 0.0_dp
                   EXIT
                END IF
             END DO
             DO i = 1, SIZE(fixd_list)
                IF (fixd_list(i)%fixd==index_d) THEN
                   IF (fixd_list(i)%itype/=use_perd_xyz) CYCLE
                   IF (.NOT.fixd_list(i)%restraint%active) imass4 = 0.0_dp
                   EXIT
                END IF
             END DO
          END IF
       END IF
       lg4x6 % imass1 = imass1
       lg4x6 % imass2 = imass2
       lg4x6 % imass3 = imass3
       lg4x6 % imass4 = imass4
       lg4x6 % init = .TRUE.
    END IF
  END SUBROUTINE check_fixed_atom_cns_g4x6

! *****************************************************************************
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE check_fixed_atom_cns_colv(fixd_list, colvar)
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(colvar_type), POINTER               :: colvar

    INTEGER                                  :: i, j, k

    IF (ASSOCIATED(fixd_list)) THEN
       IF (ASSOCIATED(fixd_list)) THEN
          IF (SIZE(fixd_list)>0) THEN
             DO i = 1, SIZE(colvar%i_atom)
                j = colvar%i_atom(i)
                DO k = 1, SIZE(fixd_list)
                   IF (fixd_list(k)%fixd==j) THEN
                      IF (fixd_list(k)%itype/=use_perd_xyz) CYCLE
                      IF (.NOT.fixd_list(k)%restraint%active)&
                           colvar%dsdr(:,i) = 0.0_dp
                      EXIT
                   END IF
                END DO
             END DO
          END IF
       END IF
    END IF

  END SUBROUTINE check_fixed_atom_cns_colv

! *****************************************************************************
!> \brief setup a list of local atoms on which to apply constraints/restraints
!> \author Teodoro Laino [tlaino] - 11.2008
! *****************************************************************************
  SUBROUTINE create_local_fixd_list(lfixd_list, nkind, molecule_kind_set, &
       local_particles, error)
    TYPE(local_fixd_constraint_type), &
      POINTER                                :: lfixd_list(:)
    INTEGER, INTENT(IN)                      :: nkind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_local_fixd_list', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ikind, iparticle, iparticle_local, isize, jsize, &
      ncnst, nparticle_local, nparticle_local_all, nsize, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: fixed_atom_all, &
                                                kind_index_all, &
                                                local_particle_all, work0, &
                                                work1, work2
    LOGICAL                                  :: failure
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    CALL timeset(routineN,handle)
    failure = .FALSE.
    CPPostcondition(.NOT.ASSOCIATED(lfixd_list),cp_failure_level,routineP,error,failure)
    nsize = 0
    DO ikind = 1, nkind
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind ( molecule_kind, fixd_list=fixd_list )
       IF (ASSOCIATED(fixd_list)) THEN
          nsize = nsize + SIZE(fixd_list)
       END IF
    END DO
    IF (nsize/=0) THEN
       ALLOCATE(fixed_atom_all(nsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(work0(nsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(work1(nsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(kind_index_all(nsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       nsize = 0
       DO ikind = 1, nkind
          molecule_kind => molecule_kind_set(ikind)
          CALL get_molecule_kind(molecule_kind, fixd_list=fixd_list)
          IF (ASSOCIATED(fixd_list)) THEN
             DO i = 1, SIZE(fixd_list)
                nsize = nsize + 1
                work0(nsize)          = i
                kind_index_all(nsize) = ikind
                fixed_atom_all(nsize) = fixd_list(i)%fixd
             END DO
          END IF
       END DO
       ! Sort the number of all atoms to be constrained/restrained
       CALL sort(fixed_atom_all, nsize, work1)

       ! Sort the local particles
       nparticle_local_all = 0
       DO i = 1, SIZE(local_particles%n_el)
          nparticle_local_all   = nparticle_local_all + local_particles%n_el(i)
       END DO
       ALLOCATE(local_particle_all(nparticle_local_all),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(work2(nparticle_local_all),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       nparticle_local_all = 0
       DO i = 1, SIZE(local_particles%n_el)
          nparticle_local   = local_particles%n_el(i)
          DO iparticle_local=1,nparticle_local
             nparticle_local_all = nparticle_local_all + 1
             iparticle = local_particles%list(i)%array(iparticle_local)
             local_particle_all(nparticle_local_all) = iparticle
          END DO
       END DO
       CALL sort(local_particle_all, nparticle_local_all, work2)

       ! Count the amount of local constraints/restraints
       ncnst = 0
       jsize = 1
       Loop_count: DO isize = 1, nparticle_local_all
          DO WHILE (local_particle_all(isize)>fixed_atom_all(jsize))
             jsize = jsize + 1
             IF (jsize>nsize) THEN
                jsize = nsize
                EXIT Loop_count
             END IF
          END DO
          IF  (local_particle_all(isize)==fixed_atom_all(jsize)) ncnst = ncnst + 1
       END DO Loop_count

       ! Allocate local fixed atom array
       ALLOCATE(lfixd_list(ncnst),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ! Fill array with constraints infos
       ncnst = 0
       jsize = 1
       Loop_fill: DO isize = 1, nparticle_local_all
          DO WHILE (local_particle_all(isize)>fixed_atom_all(jsize))
             jsize = jsize + 1
             IF (jsize>nsize) THEN
                jsize = nsize
                EXIT Loop_fill
             END IF
          END DO
          IF  (local_particle_all(isize)==fixed_atom_all(jsize)) THEN
             ncnst = ncnst + 1
             lfixd_list(ncnst)%ifixd_index = work0(work1(jsize))
             lfixd_list(ncnst)%ikind       = kind_index_all(work1(jsize))
          END IF
       END DO Loop_fill

       ! Deallocate working arrays
       DEALLOCATE(local_particle_all,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(work2,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(fixed_atom_all, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(work1, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(kind_index_all, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
       ! Allocate local fixed atom array with dimension 0
       ALLOCATE(lfixd_list(0),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE create_local_fixd_list

! *****************************************************************************
!> \brief destroy the list of local atoms on which to apply constraints/restraints
!>      Teodoro Laino [tlaino] - 11.2008
! *****************************************************************************
  SUBROUTINE release_local_fixd_list(lfixd_list, error)
    TYPE(local_fixd_constraint_type), &
      POINTER                                :: lfixd_list(:)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'release_local_fixd_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(lfixd_list),cp_failure_level,routineP,error,failure)
    DEALLOCATE(lfixd_list, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE release_local_fixd_list

END MODULE constraint_fxd
