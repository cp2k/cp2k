!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****** cp2k/constraint_fxd [1.0] *
!!
!!   NAME
!!     constraint_fxd
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
MODULE constraint_fxd
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
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: fixd_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: local_g3x3_constraint_type,&
                                             local_g4x6_constraint_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             update_particle_set
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fix_atom_control,&
            check_fixed_atom_cns_g3x3,&
            check_fixed_atom_cns_g4x6,&
            check_fixed_atom_cns_colv

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'constraint_fxd'

CONTAINS

!!****** constraint_fxd/fix_atom_control [1.0] *
!!
!!   NAME
!!     fix_atom_control
!!
!!   FUNCTION
!!     allows for fix atom constraints
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     - optionally apply fix atom constraint to random forces (Langevin)
!!       (04.10.206,MK)
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE fix_atom_control( force_env, error, w)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL :: w

    CHARACTER(len=*), PARAMETER :: routineN = 'fix_atom_control', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ii, ikind, iparticle, iparticle_local, &
      my_atm_fixed, nfixed_atoms, nkind, nparticle_local, stat
    LOGICAL                                  :: failure, zero_force
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: force
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.
    CALL timeset(routineN,'I','',handle)
    IF (.NOT. failure) THEN
       CALL force_env_get(force_env=force_env, subsys=subsys,error=error)
       CALL cp_subsys_get(subsys=subsys, particles=particles, local_particles=local_particles,&
            molecule_kinds_new=molecule_kinds ,error=error)

       nkind             =  molecule_kinds % n_els
       molecule_kind_set => molecule_kinds % els
       particle_set      => particles%els
       my_atm_fixed      =  0
       DO ikind = 1, nkind
          molecule_kind => molecule_kind_set(ikind)
          CALL get_molecule_kind ( molecule_kind, nfixd=nfixed_atoms )
          my_atm_fixed = my_atm_fixed + nfixed_atoms
       END DO
       IF (my_atm_fixed /=0) THEN
          IF (.NOT.PRESENT(w)) THEN
             ! Allocate scratch array
             ALLOCATE(force(3,SIZE(particles%els)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             force = 0.0_dp
             DO i = 1, SIZE(local_particles%n_el)
                nparticle_local   = local_particles%n_el(i)
                DO iparticle_local=1,nparticle_local
                   iparticle = local_particles%list(i)%array(iparticle_local)
                   force(:,iparticle)  = particle_set(iparticle)%f(:)
                END DO
             END DO
          END IF
          
          Mol_kind:  DO ikind = 1, nkind
             molecule_kind => molecule_kind_set(ikind)
             CALL get_molecule_kind ( molecule_kind, nfixd=nfixed_atoms, fixd_list=fixd_list )
             
             IF (nfixed_atoms ==0) CYCLE Mol_Kind
             DO i = 1, SIZE(local_particles%n_el)
                nparticle_local   = local_particles%n_el(i)
                DO iparticle_local=1,nparticle_local
                   iparticle = local_particles%list(i)%array(iparticle_local)
                   DO ii = 1, SIZE(fixd_list)
                      IF (fixd_list(ii)%fixd == iparticle) THEN
                         zero_force = ((fixd_list(ii)%fixd == iparticle).AND.&
                                       (.NOT.fixd_list(ii)%restraint%active))
                         IF (zero_force) THEN
                            ! Constraint
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
                               CASE (use_perd_y)
                                  force(2,iparticle) = 0.0_dp
                               CASE (use_perd_z)
                                  force(3,iparticle) = 0.0_dp
                               CASE (use_perd_xy)
                                  force(1,iparticle) = 0.0_dp
                                  force(2,iparticle) = 0.0_dp
                               CASE (use_perd_xz)
                                  force(1,iparticle) = 0.0_dp
                                  force(3,iparticle) = 0.0_dp
                               CASE (use_perd_yz)
                                  force(2,iparticle) = 0.0_dp
                                  force(3,iparticle) = 0.0_dp
                               CASE (use_perd_xyz)
                                  force(:,iparticle) = 0.0_dp
                               END SELECT
                            END IF
                            EXIT
                         END IF
                      END IF
                   END DO
                END DO
             END DO
          END DO Mol_kind
          IF (.NOT.PRESENT(w)) THEN
             CALL update_particle_set(particle_set,force_env%para_env%group,for=force)
             DEALLOCATE(force,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
       END IF
    END IF
    CALL timestop(0.0_dp,handle)
    
  END SUBROUTINE fix_atom_control

!!****** constraint_fxd/check_fixed_atom_cns_g3x3 [1.0] *
!!
!!   NAME
!!     check_fixed_atom_cns_g3x3
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
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

!!****** constraint_fxd/check_fixed_atom_cns_g4x6 [1.0] *
!!
!!   NAME
!!     check_fixed_atom_cns_g4x6
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
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

!!****** constraint_fxd/check_fixed_atom_cns_colv [1.0] *
!!
!!   NAME
!!     check_fixed_atom_cns_colv
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
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

END MODULE constraint_fxd
