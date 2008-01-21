!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****** cp2k/force_fields_all [1.0] *
!!
!!   NAME
!!     force_fields_all
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     Splitting and cleaning the original force_field_pack - May 2007
!!     Teodoro Laino - Zurich University
!!
!!   SOURCE
!******************************************************************************

MODULE force_fields_all

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set,&
                                             set_atomic_kind
  USE atoms_input,                     ONLY: read_shell_coord_input
  USE cell_types,                      ONLY: cell_type
  USE ewald_environment_types,         ONLY: ewald_environment_type
  USE external_potential_types,        ONLY: elp_potential_type,&
                                             get_potential,&
                                             set_potential
  USE f77_blas
  USE force_field_types,               ONLY: &
       allocate_bend_kind_set, allocate_bond_kind_set, &
       allocate_impr_kind_set, allocate_onfo_kind_set, &
       allocate_torsion_kind_set, allocate_ub_kind_set, bend_kind_type, &
       bond_kind_type, charmm_info_type, force_field_type, gromos_info_type, &
       impr_kind_type, input_info_type, onfo_kind_type, torsion_kind_type, &
       ub_kind_type
  USE input_constants,                 ONLY: do_ff_charmm,&
                                             do_ff_type_g87,&
                                             do_ff_type_g96,&
                                             do_ff_undef,&
                                             do_qmmm_none
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: &
       bend_type, bond_type, get_molecule_kind, impr_type, &
       molecule_kind_type, onfo_type, set_molecule_kind, shell_type, &
       torsion_type, ub_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE pair_potential,                  ONLY: get_nonbond_storage,&
                                             spline_nonbond_control
  USE pair_potential_types,            ONLY: &
       ea_type, lj_charmm_type, lj_type, nn_type, nosh_nosh, nosh_sh, &
       pair_potential_lj_create, pair_potential_pp_create, &
       pair_potential_pp_release, pair_potential_pp_type, &
       pair_potential_single_add, pair_potential_single_clean, &
       pair_potential_single_copy, pair_potential_single_type, sh_sh, &
       tersoff_type
  USE particle_types,                  ONLY: allocate_particle_set,&
                                             particle_type
  USE physcon,                         ONLY: bohr
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE shell_potential_types,           ONLY: shell_create,&
                                             shell_kind_type,&
                                             shell_release
  USE splines_types,                   ONLY: spline_data_p_release,&
                                             spline_data_p_retain,&
                                             spline_data_p_type,&
                                             spline_environment_type,&
                                             spline_factor_copy
  USE string_utilities,                ONLY: compress,&
                                             integer_to_string,&
                                             uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_fields_all'

  PRIVATE
  PUBLIC :: force_field_unique_bond,&
            force_field_unique_bend,&
            force_field_unique_ub,&
            force_field_unique_tors,&
            force_field_unique_impr,&
            force_field_unique_onfo,&
            force_field_pack_bond,&
            force_field_pack_bend,&
            force_field_pack_ub,&
            force_field_pack_tors,&
            force_field_pack_impr,&
            force_field_pack_ch,&
            force_field_pack_shell,&
            force_field_pack_sbe,&
            force_field_pack_onfo,&
            force_field_pack_sne,&
            force_field_pack_ppse,&
            force_field_pack_pp

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** force_fields_all/force_field_unique_bond [1.0] *
!!
!!   NAME
!!     force_field_unique_bond
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Determine the number of unique bond kind and allocate bond_kind_set
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_unique_bond (particle_set, &
       molecule_kind_set, molecule_set, ff_type, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_unique_bond', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_a2, &
                                                name_atm_b, name_atm_b2
    INTEGER                                  :: atm_a, atm_b, counter, first, &
                                                handle2, i, j, k, last, &
                                                natom, nbond, stat
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    INTEGER, POINTER                         :: map_bond_kind(:)
    LOGICAL                                  :: failure, found
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(bond_kind_type), DIMENSION(:), &
      POINTER                                :: bond_kind_set
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    failure = .FALSE.
    CALL timeset(routineN,'I','',handle2)
    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            nbond=nbond,bond_list=bond_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       IF(nbond>0) THEN
          ALLOCATE(map_bond_kind(nbond),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          counter=0
         IF((ff_type%ff_type==do_ff_type_g96).OR.(ff_type%ff_type==do_ff_type_g87)) THEN
          DO j=1,nbond
            map_bond_kind(j)=j
          END DO
          counter=nbond
         ELSE
          DO j=1,nbond
             atm_a = bond_list(j)%a
             atomic_kind => particle_set(atm_a+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_a)
             atm_b = bond_list(j)%b
             atomic_kind => particle_set(atm_b+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_b)
             found = .FALSE.
             DO k=1,j-1
                atm_a = bond_list(k)%a
                atomic_kind => particle_set(atm_a+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_a2)
                atm_b = bond_list(k)%b
                atomic_kind => particle_set(atm_b+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_b2)
                IF( (((name_atm_a)==(name_atm_a2)) .AND. &
                     ((name_atm_b)==(name_atm_b2))) .OR. &
                     (((name_atm_a)==(name_atm_b2)) .AND. &
                     ((name_atm_b)==(name_atm_a2))) ) THEN
                   found = .TRUE.
                   map_bond_kind(j) = map_bond_kind(k)
                   EXIT
                END IF
             END DO
             IF(.NOT.found) THEN
                counter=counter+1
                map_bond_kind(j) = counter
             END IF
          END DO
         END IF
          NULLIFY(bond_kind_set)
          CALL allocate_bond_kind_set(bond_kind_set,counter,error)
          DO j=1,nbond
             bond_list(j)%bond_kind => bond_kind_set(map_bond_kind(j))
          END DO
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               bond_kind_set=bond_kind_set,bond_list=bond_list)
          DEALLOCATE(map_bond_kind,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_unique_bond

!******************************************************************************
!!****** force_fields_all/force_field_unique_bend [1.0] *
!!
!!   NAME
!!    force_field_unique_bend 
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!      Determine the number of unique bend kind and allocate bend_kind_set
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_unique_bend (particle_set, &
       molecule_kind_set, molecule_set, ff_type, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_unique_bend', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_a2, &
                                                name_atm_b, name_atm_b2, &
                                                name_atm_c, name_atm_c2
    INTEGER                                  :: atm_a, atm_b, atm_c, counter, &
                                                first, handle2, i, j, k, &
                                                last, natom, nbend, stat
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    INTEGER, POINTER                         :: map_bend_kind(:)
    LOGICAL                                  :: failure, found
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(bend_kind_type), DIMENSION(:), &
      POINTER                                :: bend_kind_set
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    failure = .FALSE.

    CALL timeset(routineN,'I','',handle2)
    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            nbend=nbend,bend_list=bend_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       IF(nbend>0) THEN
          ALLOCATE(map_bend_kind(nbend),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          counter=0
         IF((ff_type%ff_type==do_ff_type_g96).OR.(ff_type%ff_type==do_ff_type_g87)) THEN
          DO j=1,nbend
            map_bend_kind(j)=j
          END DO
          counter=nbend
         ELSE
          DO j=1,nbend
             atm_a = bend_list(j)%a
             atomic_kind => particle_set(atm_a+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_a)
             atm_b = bend_list(j)%b
             atomic_kind => particle_set(atm_b+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_b)
             atm_c = bend_list(j)%c
             atomic_kind => particle_set(atm_c+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_c)
             found = .FALSE.
             DO k=1,j-1
                atm_a = bend_list(k)%a
                atomic_kind => particle_set(atm_a+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_a2)
                atm_b = bend_list(k)%b
                atomic_kind => particle_set(atm_b+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_b2)
                atm_c = bend_list(k)%c
                atomic_kind => particle_set(atm_c+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_c2)
                IF( (((name_atm_a)==(name_atm_a2)) .AND. &
                     ((name_atm_b)==(name_atm_b2)) .AND. &
                     ((name_atm_c)==(name_atm_c2))) .OR. &
                     (((name_atm_a)==(name_atm_c2)) .AND. &
                     ((name_atm_b)==(name_atm_b2)) .AND. &
                     ((name_atm_c)==(name_atm_a2))) ) THEN
                   found = .TRUE.
                   map_bend_kind(j) = map_bend_kind(k)
                   EXIT
                END IF
             END DO
             IF(.NOT.found) THEN
                counter=counter+1
                map_bend_kind(j) = counter
             END IF
          END DO
         END IF
          NULLIFY(bend_kind_set)
          CALL allocate_bend_kind_set(bend_kind_set,counter,error)
          DO j=1,nbend
             bend_list(j)%bend_kind => bend_kind_set(map_bend_kind(j))
          END DO
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               bend_kind_set=bend_kind_set,bend_list=bend_list)
          DEALLOCATE(map_bend_kind,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_unique_bend

!******************************************************************************
!!****** force_fields_all/force_field_unique_ub [1.0] *
!!
!!   NAME
!!     force_field_unique_ub
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Determine the number of unique Urey-Bradley kind and allocate ub_kind_set
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_unique_ub(particle_set, &
       molecule_kind_set, molecule_set,  error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_unique_ub', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_a2, &
                                                name_atm_b, name_atm_b2, &
                                                name_atm_c, name_atm_c2
    INTEGER                                  :: atm_a, atm_b, atm_c, counter, &
                                                first, handle2, i, j, k, &
                                                last, natom, nub, stat
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    INTEGER, POINTER                         :: map_ub_kind(:)
    LOGICAL                                  :: failure, found
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(ub_kind_type), DIMENSION(:), &
      POINTER                                :: ub_kind_set
    TYPE(ub_type), DIMENSION(:), POINTER     :: ub_list

    failure = .FALSE.

    CALL timeset(routineN,'I','',handle2)
    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            nub=nub,ub_list=ub_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       IF(nub>0) THEN
          ALLOCATE(map_ub_kind(nub),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          counter=0
          DO j=1,nub
             atm_a = ub_list(j)%a
             atomic_kind => particle_set(atm_a+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_a)
             atm_b = ub_list(j)%b
             atomic_kind => particle_set(atm_b+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_b)
             atm_c = ub_list(j)%c
             atomic_kind => particle_set(atm_c+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_c)
             found = .FALSE.
             DO k=1,j-1
                atm_a = ub_list(k)%a
                atomic_kind => particle_set(atm_a+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_a2)
                atm_b = ub_list(k)%b
                atomic_kind => particle_set(atm_b+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_b2)
                atm_c = ub_list(k)%c
                atomic_kind => particle_set(atm_c+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_c2)
                IF( (((name_atm_a)==(name_atm_a2)) .AND. &
                     ((name_atm_b)==(name_atm_b2)) .AND. &
                     ((name_atm_c)==(name_atm_c2))) .OR. &
                     (((name_atm_a)==(name_atm_c2)) .AND. &
                     ((name_atm_b)==(name_atm_b2)) .AND. &
                     ((name_atm_c)==(name_atm_a2))) ) THEN
                   found = .TRUE.
                   map_ub_kind(j) = map_ub_kind(k)
                   EXIT
                END IF
             END DO
             IF(.NOT.found) THEN
                counter=counter+1
                map_ub_kind(j) = counter
             END IF
          END DO
          CALL allocate_ub_kind_set(ub_kind_set,counter,error)
          DO j=1,nub
             ub_list(j)%ub_kind => ub_kind_set(map_ub_kind(j))
          END DO
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               ub_kind_set=ub_kind_set,ub_list=ub_list)
          DEALLOCATE(map_ub_kind,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_unique_ub

!******************************************************************************
!!****** force_fields_all/force_field_unique_tors [1.0] *
!!
!!   NAME
!!     force_field_unique_tors
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!      Determine the number of unique torsion kind and allocate torsion_kind_set
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_unique_tors(particle_set, &
       molecule_kind_set, molecule_set, ff_type, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_unique_tors', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length) :: name_atm_a, name_atm_a2, &
      name_atm_b, name_atm_b2, name_atm_c, name_atm_c2, name_atm_d, &
      name_atm_d2
    INTEGER                                  :: atm_a, atm_b, atm_c, atm_d, &
                                                counter, first, handle2, i, &
                                                j, k, last, natom, ntorsion, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    INTEGER, POINTER                         :: map_torsion_kind(:)
    LOGICAL                                  :: failure, found
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(torsion_kind_type), DIMENSION(:), &
      POINTER                                :: torsion_kind_set
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: torsion_list

    failure = .FALSE.
    CALL timeset(routineN,'I','',handle2)
    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            ntorsion=ntorsion,torsion_list=torsion_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       IF(ntorsion>0) THEN
          ALLOCATE(map_torsion_kind(ntorsion),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          counter=0
         IF((ff_type%ff_type==do_ff_type_g96).OR.(ff_type%ff_type==do_ff_type_g87)) THEN
          DO j=1,ntorsion
            map_torsion_kind(j)=j
          END DO
          counter=ntorsion
         ELSE
          DO j=1,ntorsion
             atm_a = torsion_list(j)%a
             atomic_kind => particle_set(atm_a+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_a)
             atm_b = torsion_list(j)%b
             atomic_kind => particle_set(atm_b+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_b)
             atm_c = torsion_list(j)%c
             atomic_kind => particle_set(atm_c+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_c)
             atm_d = torsion_list(j)%d
             atomic_kind => particle_set(atm_d+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_d)
             found = .FALSE.
             DO k=1,j-1
                atm_a = torsion_list(k)%a
                atomic_kind => particle_set(atm_a+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_a2)
                atm_b = torsion_list(k)%b
                atomic_kind => particle_set(atm_b+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_b2)
                atm_c = torsion_list(k)%c
                atomic_kind => particle_set(atm_c+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_c2)
                atm_d = torsion_list(k)%d
                atomic_kind => particle_set(atm_d+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_d2)
                IF( (((name_atm_a)==(name_atm_a2)) .AND. &
                     ((name_atm_b)==(name_atm_b2)) .AND. &
                     ((name_atm_c)==(name_atm_c2)) .AND. &
                     ((name_atm_d)==(name_atm_d2))) .OR. &
                     (((name_atm_a)==(name_atm_d2)) .AND. &
                     ((name_atm_b)==(name_atm_c2)) .AND. &
                     ((name_atm_c)==(name_atm_b2)) .AND. &
                     ((name_atm_d)==(name_atm_a2))) ) THEN
                   found = .TRUE.
                   map_torsion_kind(j) = map_torsion_kind(k)
                   EXIT
                END IF
             END DO
             IF(.NOT.found) THEN
                counter=counter+1
                map_torsion_kind(j) = counter
             END IF
          END DO
         END IF
          NULLIFY(torsion_kind_set)
          CALL allocate_torsion_kind_set(torsion_kind_set,counter,error)
          DO j=1,ntorsion
             torsion_list(j)%torsion_kind => torsion_kind_set(map_torsion_kind(j))
          END DO
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               torsion_kind_set=torsion_kind_set,torsion_list=torsion_list)
          DEALLOCATE(map_torsion_kind,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_unique_tors

!******************************************************************************
!!****** force_fields_all/force_field_unique_impr [1.0] *
!!
!!   NAME
!!     force_field_unique_impr
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Determine the number of unique impr kind and allocate impr_kind_set
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_unique_impr (particle_set, &
       molecule_kind_set, molecule_set, ff_type, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_unique_impr', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length) :: name_atm_a, name_atm_a2, &
      name_atm_b, name_atm_b2, name_atm_c, name_atm_c2, name_atm_d, &
      name_atm_d2
    INTEGER                                  :: atm_a, atm_b, atm_c, atm_d, &
                                                counter, first, handle2, i, &
                                                j, k, last, natom, nimpr, stat
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    INTEGER, POINTER                         :: map_impr_kind(:)
    LOGICAL                                  :: failure, found
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(impr_kind_type), DIMENSION(:), &
      POINTER                                :: impr_kind_set
    TYPE(impr_type), DIMENSION(:), POINTER   :: impr_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    failure = .FALSE.

    CALL timeset(routineN,'I','',handle2)
    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            nimpr=nimpr,impr_list=impr_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       IF(nimpr>0) THEN
          ALLOCATE(map_impr_kind(nimpr),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          counter=0
         IF((ff_type%ff_type==do_ff_type_g96).OR.(ff_type%ff_type==do_ff_type_g87)) THEN
          DO j=1,nimpr
            map_impr_kind(j)=j
          END DO
          counter=nimpr
         ELSE
          DO j=1,nimpr
             atm_a = impr_list(j)%a
             atomic_kind => particle_set(atm_a+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_a)
             atm_b = impr_list(j)%b
             atomic_kind => particle_set(atm_b+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_b)
             atm_c = impr_list(j)%c
             atomic_kind => particle_set(atm_c+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_c)
             atm_d = impr_list(j)%d
             atomic_kind => particle_set(atm_d+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_d)
             found = .FALSE.
             DO k=1,j-1
                atm_a = impr_list(k)%a
                atomic_kind => particle_set(atm_a+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_a2)
                atm_b = impr_list(k)%b
                atomic_kind => particle_set(atm_b+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_b2)
                atm_c = impr_list(k)%c
                atomic_kind => particle_set(atm_c+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_c2)
                atm_d = impr_list(k)%d
                atomic_kind => particle_set(atm_d+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_d2)
                IF( (((name_atm_a)==(name_atm_a2)) .AND. &
                     ((name_atm_b)==(name_atm_b2)) .AND. &
                     ((name_atm_c)==(name_atm_c2)) .AND. &
                     ((name_atm_d)==(name_atm_d2))) .OR. &
                     (((name_atm_a)==(name_atm_d2)) .AND. &
                     ((name_atm_b)==(name_atm_c2)) .AND. &
                     ((name_atm_c)==(name_atm_b2)) .AND. &
                     ((name_atm_d)==(name_atm_a2))) ) THEN
                   found = .TRUE.
                   map_impr_kind(j) = map_impr_kind(k)
                   EXIT
                END IF
             END DO
             IF(.NOT.found) THEN
                counter=counter+1
                map_impr_kind(j) = counter
             END IF
          END DO
         END IF
          NULLIFY(impr_kind_set)
          CALL allocate_impr_kind_set(impr_kind_set,counter,error)
          DO j=1,nimpr
             impr_list(j)%impr_kind => impr_kind_set(map_impr_kind(j))
          END DO
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               impr_kind_set=impr_kind_set,impr_list=impr_list)
          DEALLOCATE(map_impr_kind,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_unique_impr

!******************************************************************************
!!****** force_fields_all/force_field_unique_onfo [1.0] *
!!
!!   NAME
!!     force_field_unique_onfo
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Determine the number of unique onfo kind and allocate onfo_kind_set
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_unique_onfo  (particle_set, &
       molecule_kind_set, molecule_set, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_unique_onfo', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_a2, &
                                                name_atm_b, name_atm_b2
    INTEGER                                  :: atm_a, atm_b, counter, first, &
                                                handle2, i, j, k, last, &
                                                natom, nonfo, stat
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    INTEGER, POINTER                         :: map_onfo_kind(:)
    LOGICAL                                  :: failure, found
    REAL(KIND=dp)                            :: qi1, qi2, qi3, qi4
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(onfo_kind_type), DIMENSION(:), &
      POINTER                                :: onfo_kind_set
    TYPE(onfo_type), DIMENSION(:), POINTER   :: onfo_list

    failure = .FALSE.

    CALL timeset(routineN,'I','',handle2)
    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            nonfo=nonfo,onfo_list=onfo_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       IF(nonfo>0) THEN
          ALLOCATE(map_onfo_kind(nonfo),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          counter=0
          DO j=1,nonfo
             atm_a = onfo_list(j)%a
             atomic_kind => particle_set(atm_a+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_a,qeff=qi1)
             atm_b = onfo_list(j)%b
             atomic_kind => particle_set(atm_b+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_b,qeff=qi2)
             found = .FALSE.
             DO k=1,j-1
                atm_a = onfo_list(k)%a
                atomic_kind => particle_set(atm_a+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_a2,qeff=qi3)
                atm_b = onfo_list(k)%b
                atomic_kind => particle_set(atm_b+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_b2,qeff=qi4)
                IF( (((name_atm_a)==(name_atm_a2)) .AND. &
                     ((name_atm_b)==(name_atm_b2))) .OR. &
                     (((name_atm_a)==(name_atm_b2)) .AND. &
                     ((name_atm_b)==(name_atm_a2))) ) THEN
                   IF(((qi1==qi3).AND.(qi2==qi4)).OR.((qi1==qi4).AND.(qi2==qi3)))THEN
                      found = .TRUE.
                      map_onfo_kind(j) = map_onfo_kind(k)
                      EXIT
                   END IF
                END IF
             END DO
             IF(.NOT.found) THEN
                counter=counter+1
                map_onfo_kind(j) = counter
             END IF
          END DO
          CALL allocate_onfo_kind_set(onfo_kind_set,counter,error)
          DO j=1,nonfo
             onfo_list(j)%onfo_kind => onfo_kind_set(map_onfo_kind(j))
          END DO
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               onfo_kind_set=onfo_kind_set,onfo_list=onfo_list)
          DEALLOCATE(map_onfo_kind,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_unique_onfo

!******************************************************************************
!!****** force_fields_all/force_field_pack_bond [1.0] *
!!
!!   NAME
!!     force_field_pack_bond
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Pack in bonds information needed for the force_field
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_pack_bond (particle_set, molecule_kind_set, molecule_set, &
      fatal, Ainfo, chm_info, inp_info, gro_info, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    LOGICAL                                  :: fatal
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: Ainfo
    TYPE(charmm_info_type), POINTER          :: chm_info
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(gromos_info_type), POINTER          :: gro_info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_bond', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_b
    INTEGER                                  :: atm_a, atm_b, first, handle2, &
                                                i, itype, j, k, last, natom, &
                                                nbond
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    LOGICAL                                  :: failure, found, only_qm
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    failure = .FALSE.
    CALL timeset(routineN,'I','',handle2)

    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            nbond=nbond,bond_list=bond_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       DO j=1,nbond
          atm_a = bond_list(j)%a
          atomic_kind => particle_set(atm_a+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_a)
          atm_b = bond_list(j)%b
          atomic_kind => particle_set(atm_b+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_b)
          found = .FALSE.
          only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b)
          CALL uppercase(name_atm_a)
          CALL uppercase(name_atm_b)

          ! loop over params from gromos
          IF(ASSOCIATED(gro_info%bond_k)) THEN
             k=SIZE(gro_info%bond_k)
             itype = bond_list(j)%itype
             IF(itype<=k) THEN
                bond_list(j)%bond_kind%k(1) = gro_info%bond_k(itype)
                bond_list(j)%bond_kind%r0   = gro_info%bond_r0(itype)
             ELSE
                itype=itype-k
                bond_list(j)%bond_kind%k(1) = gro_info%solvent_k(itype)
                bond_list(j)%bond_kind%r0   = gro_info%solvent_r0(itype)
             END IF
             bond_list(j)%bond_kind%id_type =  bond_list(j)%id_type
             found = .TRUE.
          END IF

          ! loop over params from charmm
          IF(ASSOCIATED(chm_info%bond_a)) THEN
             DO k=1,SIZE(chm_info%bond_a)
                IF( (((chm_info%bond_a(k))==(name_atm_a)) .AND. &
                     ((chm_info%bond_b(k))==(name_atm_b))) .OR. &
                     (((chm_info%bond_a(k))==(name_atm_b)) .AND. &
                     ((chm_info%bond_b(k))==(name_atm_a))) ) THEN
                   bond_list(j)%bond_kind%id_type = do_ff_charmm
                   bond_list(j)%bond_kind%k(1) = chm_info%bond_k(k)
                   bond_list(j)%bond_kind%r0   = chm_info%bond_r0(k)
                   CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                        "Multiple Bond declarations: "//TRIM(name_atm_a)//" and "//&
                        TRIM(name_atm_b)//" overwriting! "//&
CPSourceFileRef,&
                        only_ionode=.TRUE.)
                   found = .TRUE.
                   EXIT
                END IF
             END DO
          END IF

          ! always have the input param last to overwrite all the other ones
          IF(ASSOCIATED(inp_info%bond_a)) THEN
             DO k=1,SIZE(inp_info%bond_a)
                IF( (((inp_info%bond_a(k))==(name_atm_a)) .AND. &
                     ((inp_info%bond_b(k))==(name_atm_b))) .OR. &
                     (((inp_info%bond_a(k))==(name_atm_b)) .AND. &
                     ((inp_info%bond_b(k))==(name_atm_a))) ) THEN
                   bond_list(j)%bond_kind%id_type = inp_info%bond_kind(k)
                   bond_list(j)%bond_kind%k(:)    = inp_info%bond_k(:,k)
                   bond_list(j)%bond_kind%r0      = inp_info%bond_r0(k)
                   CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                        "Multiple Bond declarations: "//TRIM(name_atm_a)//" and "//&
                        TRIM(name_atm_b)//" overwriting! "//&
CPSourceFileRef,&
                        only_ionode=.TRUE.)
                   found = .TRUE.
                   EXIT
                END IF
             END DO
          END IF

          IF(.NOT.found) CALL store_FF_missing_par(atm1=TRIM(name_atm_a),&
                                                   atm2=TRIM(name_atm_b),&
                                                   fatal=fatal,&
                                                   type_name="Bond",&
                                                   array=Ainfo,&
                                                   error=error)
          ! QM/MM modifications
          IF (only_qm) THEN
             bond_list(j)%id_type           = do_ff_undef
             bond_list(j)%bond_kind%id_type = do_ff_undef
          END IF
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            bond_list=bond_list)
    END DO
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_pack_bond

!******************************************************************************
!!****** force_fields_all/force_field_pack_bend [1.0] *
!!
!!   NAME
!!     force_field_pack_bend
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Pack in bends information needed for the force_field
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_pack_bend (particle_set, molecule_kind_set, molecule_set, &
       fatal, Ainfo, chm_info, inp_info, gro_info, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    LOGICAL                                  :: fatal
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: Ainfo
    TYPE(charmm_info_type), POINTER          :: chm_info
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(gromos_info_type), POINTER          :: gro_info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_bend', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_b, &
                                                name_atm_c
    INTEGER                                  :: atm_a, atm_b, atm_c, first, &
                                                handle2, i, itype, j, k, &
                                                last, natom, nbend
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    LOGICAL                                  :: failure, found, only_qm
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset(routineN,'I','',handle2)
    failure = .FALSE.

    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            nbend=nbend,bend_list=bend_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       DO j=1,nbend
          atm_a = bend_list(j)%a
          atomic_kind => particle_set(atm_a+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_a)
          atm_b = bend_list(j)%b
          atomic_kind => particle_set(atm_b+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_b)
          atm_c = bend_list(j)%c
          atomic_kind => particle_set(atm_c+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_c)
          found = .FALSE.
          only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b, id3=name_atm_c)
          CALL uppercase(name_atm_a)
          CALL uppercase(name_atm_b)
          CALL uppercase(name_atm_c)

          ! loop over params from GROMOS
          IF(ASSOCIATED(gro_info%bend_k)) THEN
             k=SIZE(gro_info%bend_k)
             itype = bend_list(j)%itype
             IF(itype>0) THEN
               bend_list(j)%bend_kind%k      = gro_info%bend_k(itype)
               bend_list(j)%bend_kind%theta0 = gro_info%bend_theta0(itype)
             ELSE
               bend_list(j)%bend_kind%k      = gro_info%bend_k(itype/k)
               bend_list(j)%bend_kind%theta0 = gro_info%bend_theta0(itype/k)
             END IF
             bend_list(j)%bend_kind%id_type = bend_list(j)%id_type
             found = .TRUE.
          END IF

          ! loop over params from CHARMM
          IF(ASSOCIATED(chm_info%bend_a)) THEN
             DO k=1,SIZE(chm_info%bend_a)
                IF( (((chm_info%bend_a(k))==(name_atm_a)) .AND. &
                     ((chm_info%bend_b(k))==(name_atm_b)) .AND. &
                     ((chm_info%bend_c(k))==(name_atm_c))) .OR. &
                     (((chm_info%bend_a(k))==(name_atm_c)) .AND. &
                     ((chm_info%bend_b(k))==(name_atm_b)) .AND. &
                     ((chm_info%bend_c(k))==(name_atm_a))) ) THEN
                   bend_list(j)%bend_kind%id_type = do_ff_charmm
                   bend_list(j)%bend_kind%k       = chm_info%bend_k(k)
                   bend_list(j)%bend_kind%theta0  = chm_info%bend_theta0(k)
                   CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                        "Multiple Bend declarations: "//TRIM(name_atm_a)//&
                        ","//TRIM(name_atm_b)//" and "//TRIM(name_atm_c)//" overwriting! "//&
CPSourceFileRef,&
                        only_ionode=.TRUE.)
                   found = .TRUE.
                   EXIT
                END IF
             END DO
          END IF

          ! always have the input param last to overwrite all the other ones
          IF(ASSOCIATED(inp_info%bend_a)) THEN
             DO k=1,SIZE(inp_info%bend_a)
                IF( (((inp_info%bend_a(k))==(name_atm_a)) .AND. &
                     ((inp_info%bend_b(k))==(name_atm_b)) .AND. &
                     ((inp_info%bend_c(k))==(name_atm_c))) .OR. &
                     (((inp_info%bend_a(k))==(name_atm_c)) .AND. &
                     ((inp_info%bend_b(k))==(name_atm_b)) .AND. &
                     ((inp_info%bend_c(k))==(name_atm_a))) ) THEN
                   bend_list(j)%bend_kind%id_type = inp_info%bend_kind(k)
                   bend_list(j)%bend_kind%k       = inp_info%bend_k(k)
                   bend_list(j)%bend_kind%theta0  = inp_info%bend_theta0(k)
                   CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                        "Multiple Bend declarations: "//TRIM(name_atm_a)//&
                        ","//TRIM(name_atm_b)//" and "//TRIM(name_atm_c)//" overwriting! "//&
CPSourceFileRef,&
                        only_ionode=.TRUE.)
                   found = .TRUE.
                   EXIT
                END IF
             END DO
          END IF

          IF(.NOT.found) CALL store_FF_missing_par(atm1=TRIM(name_atm_a),&
                                                   atm2=TRIM(name_atm_b),&
                                                   atm3=TRIM(name_atm_c),&
                                                   fatal=fatal,&
                                                   type_name="Angle",&
                                                   array=Ainfo,&
                                                   error=error)
          ! QM/MM modifications
          IF (only_qm) THEN
             bend_list(j)%id_type           = do_ff_undef
             bend_list(j)%bend_kind%id_type = do_ff_undef
          END IF
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            bend_list=bend_list)
    END DO
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_pack_bend

!******************************************************************************
!!****** force_fields_all/force_field_pack_ub [1.0] *
!!
!!   NAME
!!     force_field_pack_ub
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Pack in Urey-Bradley information needed for the force_field
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_pack_ub (particle_set, molecule_kind_set, molecule_set, &
       Ainfo, chm_info, inp_info, iw, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: Ainfo
    TYPE(charmm_info_type), POINTER          :: chm_info
    TYPE(input_info_type), POINTER           :: inp_info
    INTEGER                                  :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_ub', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_b, &
                                                name_atm_c
    INTEGER                                  :: atm_a, atm_b, atm_c, first, &
                                                handle2, i, j, k, last, &
                                                natom, nub
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    LOGICAL                                  :: failure, found, only_qm
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(ub_type), DIMENSION(:), POINTER     :: ub_list

    failure = .FALSE.

    CALL timeset(routineN,'I','',handle2)
    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            nub=nub,ub_list=ub_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       DO j=1,nub
          atm_a = ub_list(j)%a
          atomic_kind => particle_set(atm_a+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_a)
          atm_b = ub_list(j)%b
          atomic_kind => particle_set(atm_b+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_b)
          atm_c = ub_list(j)%c
          atomic_kind => particle_set(atm_c+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_c)
          found = .FALSE.
          only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b, id3=name_atm_c)
          CALL uppercase(name_atm_a)
          CALL uppercase(name_atm_b)
          CALL uppercase(name_atm_c)

          ! loop over params from GROMOS
          ! ikuo - None that I know...

          ! loop over params from charmm
          IF(ASSOCIATED(chm_info%ub_a)) THEN
             DO k=1,SIZE(chm_info%ub_a)
                IF( (((chm_info%ub_a(k))==(name_atm_a)) .AND. &
                     ((chm_info%ub_b(k))==(name_atm_b)) .AND. &
                     ((chm_info%ub_c(k))==(name_atm_c))) .OR. &
                     (((chm_info%ub_a(k))==(name_atm_c)) .AND. &
                     ((chm_info%ub_b(k))==(name_atm_b)) .AND. &
                     ((chm_info%ub_c(k))==(name_atm_a))) ) THEN
                   ub_list(j)%ub_kind%id_type = do_ff_charmm
                   ub_list(j)%ub_kind%k(1)    = chm_info%ub_k(k)
                   ub_list(j)%ub_kind%r0      = chm_info%ub_r0(k)
                   IF(iw>0) WRITE(iw,*) "    Found UB ",TRIM(name_atm_a)," ",&
                        TRIM(name_atm_b)," ",TRIM(name_atm_c)
                   CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                        "Multiple UB declarations: "//TRIM(name_atm_a)//&
                        ","//TRIM(name_atm_b)//" and "//TRIM(name_atm_c)//" overwriting! "//&
CPSourceFileRef,&
                        only_ionode=.TRUE.)
                   found = .TRUE.
                   EXIT
                END IF
             END DO
          END IF

          ! always have the input param last to overwrite all the other ones
          IF(ASSOCIATED(inp_info%ub_a)) THEN
             DO k=1,SIZE(inp_info%ub_a)
                IF( (((inp_info%ub_a(k))==(name_atm_a)) .AND. &
                     ((inp_info%ub_b(k))==(name_atm_b)) .AND. &
                     ((inp_info%ub_c(k))==(name_atm_c))) .OR. &
                    (((inp_info%ub_a(k))==(name_atm_c)) .AND. &
                     ((inp_info%ub_b(k))==(name_atm_b)) .AND. &
                     ((inp_info%ub_c(k))==(name_atm_a))) ) THEN
                   ub_list(j)%ub_kind%id_type = inp_info%ub_kind(k)
                   ub_list(j)%ub_kind%k(:)    = inp_info%ub_k(:,k)
                   ub_list(j)%ub_kind%r0      = inp_info%ub_r0(k)
                   CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                        "Multiple UB declarations: "//TRIM(name_atm_a)//&
                        ","//TRIM(name_atm_b)//" and "//TRIM(name_atm_c)//" overwriting! "//&
CPSourceFileRef,&
                        only_ionode=.TRUE.)
                   found = .TRUE.
                   EXIT
                END IF
             END DO
          END IF

          IF(.NOT.found) THEN
             CALL store_FF_missing_par(atm1=TRIM(name_atm_a),&
                                       atm2=TRIM(name_atm_b),&
                                       atm3=TRIM(name_atm_c),&
                                       type_name="Urey-Bradley",&
                                       array=Ainfo,&
                                       error=error)
             ub_list(j)%id_type         = do_ff_undef
             ub_list(j)%ub_kind%id_type = do_ff_undef
             ub_list(j)%ub_kind%k       = 0.0_dp
             ub_list(j)%ub_kind%r0      = 0.0_dp
          END IF

          ! QM/MM modifications
          IF (only_qm) THEN
             ub_list(j)%id_type         = do_ff_undef
             ub_list(j)%ub_kind%id_type = do_ff_undef
          END IF
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            ub_list=ub_list)
    END DO
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_pack_ub

!******************************************************************************
!!****** force_fields_all/force_field_pack_tors [1.0] *
!!
!!   NAME
!!     force_field_pack_tors
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Pack in torsion information needed for the force_field
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_pack_tors (particle_set, molecule_kind_set, molecule_set, &
       Ainfo, chm_info, inp_info, gro_info, iw, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: Ainfo
    TYPE(charmm_info_type), POINTER          :: chm_info
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(gromos_info_type), POINTER          :: gro_info
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_tors', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: ldum, name_atm_a, name_atm_b, &
                                                name_atm_c, name_atm_d
    INTEGER                                  :: atm_a, atm_b, atm_c, atm_d, &
                                                first, handle2, i, imul, &
                                                itype, j, k, last, natom, &
                                                ntorsion
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    LOGICAL                                  :: failure, found, only_qm
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: torsion_list

    CALL timeset(routineN,'I','',handle2)
    failure = .FALSE.

    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            ntorsion=ntorsion,torsion_list=torsion_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       DO j=1,ntorsion
          IF(torsion_list(j)%torsion_kind%id_type == do_ff_undef) THEN
             atm_a = torsion_list(j)%a
             atomic_kind => particle_set(atm_a+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_a)
             atm_b = torsion_list(j)%b
             atomic_kind => particle_set(atm_b+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_b)
             atm_c = torsion_list(j)%c
             atomic_kind => particle_set(atm_c+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_c)
             atm_d = torsion_list(j)%d
             atomic_kind => particle_set(atm_d+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_d)
             found = .FALSE.
             only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b, id3=name_atm_c, id4=name_atm_d)
             CALL uppercase(name_atm_a)
             CALL uppercase(name_atm_b)
             CALL uppercase(name_atm_c)
             CALL uppercase(name_atm_d)

             ! loop over params from GROMOS
             IF(ASSOCIATED(gro_info%torsion_k)) THEN
                k=SIZE(gro_info%torsion_k)
                itype = torsion_list(j)%itype
                IF(itype>0) THEN
                   CALL reallocate(torsion_list(j)%torsion_kind%k,1,1)
                   CALL reallocate(torsion_list(j)%torsion_kind%m,1,1)
                   CALL reallocate(torsion_list(j)%torsion_kind%cosphi0,1,1)
                   torsion_list(j)%torsion_kind%nmul   = 1
                   torsion_list(j)%torsion_kind%m(1)   = gro_info%torsion_m(itype)
                   torsion_list(j)%torsion_kind%k(1)   = gro_info%torsion_k(itype)
                   torsion_list(j)%torsion_kind%cosphi0(1) = gro_info%torsion_cosphi0(itype)
                ELSE
                   CALL reallocate(torsion_list(j)%torsion_kind%k,1,1)
                   CALL reallocate(torsion_list(j)%torsion_kind%m,1,1)
                   CALL reallocate(torsion_list(j)%torsion_kind%cosphi0,1,1)
                   torsion_list(j)%torsion_kind%nmul   = 1
                   torsion_list(j)%torsion_kind%m(1)   = gro_info%torsion_m(itype/k)
                   torsion_list(j)%torsion_kind%k(1)   = gro_info%torsion_k(itype/k)
                   torsion_list(j)%torsion_kind%cosphi0(1) = gro_info%torsion_cosphi0(itype/k)
                END IF
                torsion_list(j)%torsion_kind%id_type = torsion_list(j)%id_type
                found = .TRUE.
                imul = torsion_list(j)%torsion_kind%nmul
             END IF

             ! loop over params from CHARMM
             IF(ASSOCIATED(chm_info%torsion_a)) THEN
                DO k=1,SIZE(chm_info%torsion_a)
                   IF( (((chm_info%torsion_a(k))==(name_atm_a)) .AND. &
                        ((chm_info%torsion_b(k))==(name_atm_b)) .AND. &
                        ((chm_info%torsion_c(k))==(name_atm_c)) .AND. &
                        ((chm_info%torsion_d(k))==(name_atm_d))) .OR. &
                       (((chm_info%torsion_a(k))==(name_atm_d)) .AND. &
                        ((chm_info%torsion_b(k))==(name_atm_c)) .AND. &
                        ((chm_info%torsion_c(k))==(name_atm_b)) .AND. &
                        ((chm_info%torsion_d(k))==(name_atm_a))) ) THEN
                      imul = torsion_list(j)%torsion_kind%nmul + 1
                      CALL reallocate(torsion_list(j)%torsion_kind%k,1,imul)
                      CALL reallocate(torsion_list(j)%torsion_kind%m,1,imul)
                      CALL reallocate(torsion_list(j)%torsion_kind%cosphi0,1,imul)
                      torsion_list(j)%torsion_kind%id_type = do_ff_charmm
                      torsion_list(j)%torsion_kind%k(imul) = chm_info%torsion_k(k)
                      torsion_list(j)%torsion_kind%m(imul) = chm_info%torsion_m(k)
                      torsion_list(j)%torsion_kind%cosphi0(imul) = chm_info%torsion_cosphi0(k)
                      torsion_list(j)%torsion_kind%nmul = imul
                      found = .TRUE.
                   END IF
                END DO

                IF(.NOT.found) THEN
                   DO k=1,SIZE(chm_info%torsion_a)
                      IF( (((chm_info%torsion_a(k))==("X")) .AND. &
                           ((chm_info%torsion_b(k))==(name_atm_b)) .AND. &
                           ((chm_info%torsion_c(k))==(name_atm_c)) .AND. &
                           ((chm_info%torsion_d(k))==("X"))) .OR. &
                           (((chm_info%torsion_a(k))==("X")) .AND. &
                           ((chm_info%torsion_b(k))==(name_atm_c)) .AND. &
                           ((chm_info%torsion_c(k))==(name_atm_b)) .AND. &
                           ((chm_info%torsion_d(k))==("X"))) ) THEN
                         imul = torsion_list(j)%torsion_kind%nmul + 1
                         CALL reallocate(torsion_list(j)%torsion_kind%k,1,imul)
                         CALL reallocate(torsion_list(j)%torsion_kind%m,1,imul)
                         CALL reallocate(torsion_list(j)%torsion_kind%cosphi0,1,imul)
                         torsion_list(j)%torsion_kind%id_type = do_ff_charmm
                         torsion_list(j)%torsion_kind%k(imul) = chm_info%torsion_k(k)
                         torsion_list(j)%torsion_kind%m(imul) = chm_info%torsion_m(k)
                         torsion_list(j)%torsion_kind%cosphi0(imul) = chm_info%torsion_cosphi0(k)
                         torsion_list(j)%torsion_kind%nmul = imul
                         found = .TRUE.
                      END IF
                   END DO
                END IF
             END IF
             ! always have the input param last to overwrite all the other ones
             IF(ASSOCIATED(inp_info%torsion_a)) THEN
                DO k=1,SIZE(inp_info%torsion_a)
                   IF( (((inp_info%torsion_a(k))==(name_atm_a)) .AND. &
                        ((inp_info%torsion_b(k))==(name_atm_b)) .AND. &
                        ((inp_info%torsion_c(k))==(name_atm_c)) .AND. &
                        ((inp_info%torsion_d(k))==(name_atm_d))) .OR. &
                        (((inp_info%torsion_a(k))==(name_atm_d)) .AND. &
                        ((inp_info%torsion_b(k))==(name_atm_c)) .AND. &
                        ((inp_info%torsion_c(k))==(name_atm_b)) .AND. &
                        ((inp_info%torsion_d(k))==(name_atm_a))) ) THEN
                      imul = torsion_list(j)%torsion_kind%nmul + 1
                      CALL reallocate(torsion_list(j)%torsion_kind%k,1,imul)
                      CALL reallocate(torsion_list(j)%torsion_kind%m,1,imul)
                      CALL reallocate(torsion_list(j)%torsion_kind%cosphi0,1,imul)
                      torsion_list(j)%torsion_kind%id_type = inp_info%torsion_kind(k)
                      torsion_list(j)%torsion_kind%k(imul) = inp_info%torsion_k(k)
                      torsion_list(j)%torsion_kind%m(imul) = inp_info%torsion_m(k)
                      torsion_list(j)%torsion_kind%cosphi0(imul) = COS(inp_info%torsion_cosphi0(k))
                      torsion_list(j)%torsion_kind%nmul    = imul
                      found = .TRUE.
                   END IF
                END DO
             END IF

             IF(.NOT.found) THEN
                CALL store_FF_missing_par(atm1=TRIM(name_atm_a),&
                                          atm2=TRIM(name_atm_b),&
                                          atm3=TRIM(name_atm_c),&
                                          atm4=TRIM(name_atm_c),&
                                          type_name="Torsion",&
                                          array=Ainfo,&
                                          error=error)
                torsion_list(j)%torsion_kind%id_type = do_ff_undef
                torsion_list(j)%id_type              = do_ff_undef
             ELSE
                ldum = cp_to_string(imul)
                IF ((imul /= 1).AND.(iw>0))&
                     WRITE(iw,'(/,2("UTIL_INFO| ",A,/))')&
                     "Multiple Torsion declarations: "//TRIM(name_atm_a)//&
                     ","//TRIM(name_atm_b)//","//TRIM(name_atm_c)//" and "//TRIM(name_atm_d),&
                     "Number of defined torsions: "//TRIM(ldum)//" ."
             END IF
             !
             ! QM/MM modifications
             !
             IF (only_qm) THEN
                IF (iw>0) WRITE(iw,*)"    Torsion PARAM between QM atoms ",j," : ",&
                     TRIM(name_atm_a)," ",&
                     TRIM(name_atm_b)," ",&
                     TRIM(name_atm_c)," ",&
                     TRIM(name_atm_d)," ",&
                     torsion_list(j)%a,&
                     torsion_list(j)%b,&
                     torsion_list(j)%c,&
                     torsion_list(j)%d
                torsion_list(j)%torsion_kind%id_type = do_ff_undef
                torsion_list(j)%id_type              = do_ff_undef
             END IF
          END IF
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            torsion_list=torsion_list)
    END DO
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_pack_tors

!******************************************************************************
!!****** force_fields_all/force_field_pack_impr [1.0] *
!!
!!   NAME
!!     force_field_pack_impr
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Pack in impropers information needed for the force_field
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_pack_impr (particle_set, molecule_kind_set, molecule_set, &
       Ainfo, chm_info, inp_info, gro_info, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: Ainfo
    TYPE(charmm_info_type), POINTER          :: chm_info
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(gromos_info_type), POINTER          :: gro_info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_impr', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_b, &
                                                name_atm_c, name_atm_d
    INTEGER                                  :: atm_a, atm_b, atm_c, atm_d, &
                                                first, handle2, i, itype, j, &
                                                k, last, natom, nimpr
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    LOGICAL                                  :: failure, found, only_qm
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(impr_type), DIMENSION(:), POINTER   :: impr_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset(routineN,'I','',handle2)
    failure = .FALSE.

    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            nimpr=nimpr,impr_list=impr_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       DO j=1,nimpr
          atm_a = impr_list(j)%a
          atomic_kind => particle_set(atm_a+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_a)
          atm_b = impr_list(j)%b
          atomic_kind => particle_set(atm_b+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_b)
          atm_c = impr_list(j)%c
          atomic_kind => particle_set(atm_c+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_c)
          atm_d = impr_list(j)%d
          atomic_kind => particle_set(atm_d+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_d)
          found = .FALSE.
          only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b, id3=name_atm_c, id4=name_atm_d)
          CALL uppercase(name_atm_a)
          CALL uppercase(name_atm_b)
          CALL uppercase(name_atm_c)
          CALL uppercase(name_atm_d)

          ! loop over params from GROMOS
          IF(ASSOCIATED(gro_info%impr_k)) THEN
            k=SIZE(gro_info%impr_k)
            itype = impr_list(j)%itype
            IF(itype>0) THEN
              impr_list(j)%impr_kind%k    = gro_info%impr_k(itype)
              impr_list(j)%impr_kind%phi0 = gro_info%impr_phi0(itype)
            ELSE
              impr_list(j)%impr_kind%k    = gro_info%impr_k(itype)
              impr_list(j)%impr_kind%phi0 = gro_info%impr_phi0(itype)
            END IF
            found = .TRUE.
            impr_list(j)%impr_kind%id_type = impr_list(j)%id_type
          END IF

          ! loop over params from charmm
          IF(ASSOCIATED(chm_info%impr_a)) THEN
             DO k=1,SIZE(chm_info%impr_a)
                IF( (((chm_info%impr_a(k))==(name_atm_a)) .AND. &
                     ((chm_info%impr_b(k))==(name_atm_b)) .AND. &
                     ((chm_info%impr_c(k))==(name_atm_c)) .AND. &
                     ((chm_info%impr_d(k))==(name_atm_d))) .OR. &
                     (((chm_info%impr_a(k))==(name_atm_d)) .AND. &
                     ((chm_info%impr_b(k))==(name_atm_c)) .AND. &
                     ((chm_info%impr_c(k))==(name_atm_b)) .AND. &
                     ((chm_info%impr_d(k))==(name_atm_a))) ) THEN
                   impr_list(j)%impr_kind%id_type = do_ff_charmm
                   impr_list(j)%impr_kind%k    = chm_info%impr_k(k)
                   impr_list(j)%impr_kind%phi0 = chm_info%impr_phi0(k)
                   CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                        "Multiple Improper declarations: "//TRIM(name_atm_a)//&
                        ","//TRIM(name_atm_b)//","//TRIM(name_atm_c)//" and "//TRIM(name_atm_d)//&
                        " overwriting! "//&
CPSourceFileRef,&
                        only_ionode=.TRUE.)
                   found = .TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.found) THEN
                DO k=1,SIZE(chm_info%impr_a)
                   IF( (((chm_info%impr_a(k))==(name_atm_a)) .AND. &
                        ((chm_info%impr_b(k))==("X")) .AND. &
                        ((chm_info%impr_c(k))==("X")) .AND. &
                        ((chm_info%impr_d(k))==(name_atm_d))) .OR. &
                        (((chm_info%impr_a(k))==(name_atm_d)) .AND. &
                        ((chm_info%impr_b(k))==("X")) .AND. &
                        ((chm_info%impr_c(k))==("X")) .AND. &
                        ((chm_info%impr_d(k))==(name_atm_a))) ) THEN
                      impr_list(j)%impr_kind%id_type = do_ff_charmm
                      impr_list(j)%impr_kind%k    = chm_info%impr_k(k)
                      impr_list(j)%impr_kind%phi0 = chm_info%impr_phi0(k)
                      CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                           "Multiple Improper declarations: "//TRIM(name_atm_a)//&
                           ","//TRIM(name_atm_b)//","//TRIM(name_atm_c)//" and "//TRIM(name_atm_d)//&
                           " overwriting! "//&
CPSourceFileRef,&
                           only_ionode=.TRUE.)
                      found = .TRUE.
                      EXIT
                   END IF
                END DO
             END IF
          END IF

          ! always have the input param last to overwrite all the other ones
          IF(ASSOCIATED(inp_info%impr_a)) THEN
             DO k=1,SIZE(inp_info%impr_a)
                IF( (((inp_info%impr_a(k))==(name_atm_a)) .AND. &
                     ((inp_info%impr_b(k))==(name_atm_b)) .AND. &
                     ((inp_info%impr_c(k))==(name_atm_c)) .AND. &
                     ((inp_info%impr_d(k))==(name_atm_d))) .OR. &
                     (((inp_info%impr_a(k))==(name_atm_d)) .AND. &
                     ((inp_info%impr_b(k))==(name_atm_c)) .AND. &
                     ((inp_info%impr_c(k))==(name_atm_b)) .AND. &
                     ((inp_info%impr_d(k))==(name_atm_a))) ) THEN
                   impr_list(j)%impr_kind%id_type = inp_info%impr_kind(k)
                   impr_list(j)%impr_kind%k       = inp_info%impr_k(k)
                   impr_list(j)%impr_kind%phi0    = inp_info%impr_phi0(k)
                   CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                        "Multiple Improper declarations: "//TRIM(name_atm_a)//&
                        ","//TRIM(name_atm_b)//","//TRIM(name_atm_c)//" and "//TRIM(name_atm_d)//&
                        " overwriting! "//&
CPSourceFileRef,&
                        only_ionode=.TRUE.)
                   found = .TRUE.
                   EXIT
                END IF
             END DO
          END IF

          IF(.NOT.found) THEN
             CALL store_FF_missing_par(atm1=TRIM(name_atm_a),&
                                       atm2=TRIM(name_atm_b),&
                                       atm3=TRIM(name_atm_c),&
                                       atm4=TRIM(name_atm_c),&
                                       type_name="Improper",&
                                       array=Ainfo,&
                                       error=error)
             impr_list(j)%impr_kind%k       = 0.0_dp
             impr_list(j)%impr_kind%phi0    = 0.0_dp
             impr_list(j)%impr_kind%id_type = do_ff_undef
             impr_list(j)%id_type           = do_ff_undef
          END IF
          !
          ! QM/MM modifications
          !
          IF (only_qm) THEN
             impr_list(j)%impr_kind%id_type = do_ff_undef
             impr_list(j)%id_type           = do_ff_undef
          END IF

       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,impr_list=impr_list)
    END DO
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_pack_impr

!******************************************************************************
!!****** force_fields_all/force_field_pack_ch [1.0] *
!!
!!   NAME
!!     force_field_pack_ch
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Set up atomic_kind_set()%elp_potentail%qeff (PART 2) and shell potential
!!     parameters
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_pack_ch  ( atomic_kind_set, qmmm_env, fatal, iw,&
       Ainfo, my_qmmm, chm_info, inp_info, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qmmm_env_mm_type), POINTER          :: qmmm_env
    LOGICAL                                  :: fatal
    INTEGER                                  :: iw
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: Ainfo
    LOGICAL                                  :: my_qmmm
    TYPE(charmm_info_type), POINTER          :: chm_info
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_ch', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atmname, name_atm_a
    INTEGER                                  :: handle2, i, ilink, j, stat
    INTEGER, DIMENSION(:), POINTER           :: my_atom_list
    LOGICAL                                  :: failure, found, found_p, &
                                                is_link_atom, only_manybody, &
                                                only_qm
    REAL(KIND=dp)                            :: scale_factor
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charge
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(elp_potential_type), POINTER        :: elp_potential

    failure = .FALSE.

    CALL timeset(routineN,'I','',handle2)
    ALLOCATE(charge(SIZE(atomic_kind_set)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i=1,SIZE(atomic_kind_set)
       atomic_kind => atomic_kind_set(i)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            elp_potential=elp_potential,&
            atom_list=my_atom_list,&
            name=atmname)
       CALL get_potential(potential=elp_potential,&
            qeff=charge(i))

       found = .FALSE.
       only_qm = qmmm_ff_precond_only_qm(id1=atmname,is_link=is_link_atom)
       CALL uppercase(atmname)
       IF(charge(i)/=-HUGE(0.0_dp)) found = .TRUE.
       ! loop over params from charmm
       IF(ASSOCIATED(chm_info%charge_atm)) THEN
          DO j=1,SIZE(chm_info%charge_atm)
             IF((chm_info%charge_atm(j))==atmname) THEN
                charge(i) = chm_info%charge(j)
                CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                     "Multiple charge declarations: "//TRIM(atmname)//&
                     " overwriting! "//&
CPSourceFileRef,&
                     only_ionode=.TRUE.)
                found = .TRUE.
             END IF
          END DO
       END IF
       ! always have the input param last to overwrite all the other ones
       IF(ASSOCIATED(inp_info%charge_atm)) THEN
          DO j=1,SIZE(inp_info%charge_atm)
             IF (iw>0) WRITE(iw,*)"Charge Checking ::",TRIM(inp_info%charge_atm(j)),atmname
             IF((inp_info%charge_atm(j))==atmname) THEN
                charge(i) = inp_info%charge(j)
                CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                     "Multiple charge declarations: "//TRIM(atmname)//&
                     " overwriting! "//&
CPSourceFileRef,&
                     only_ionode=.TRUE.)
                found = .TRUE.
             END IF
          END DO
       END IF
       IF (ASSOCIATED(inp_info%nonbonded)) THEN
          IF (ASSOCIATED(inp_info%nonbonded%pot)) THEN
             ! Let's find the nonbonded potential where this atom is involved
             only_manybody = .TRUE.
             found_p       = .FALSE.
             DO j = 1, SIZE(inp_info%nonbonded%pot)
                IF ( atmname==inp_info%nonbonded%pot(j)%pot %at1.OR.&
                     atmname==inp_info%nonbonded%pot(j)%pot %at2) THEN
                   SELECT CASE(inp_info%nonbonded%pot(j)%pot%type(1))
                   CASE (ea_type,tersoff_type)
                      ! Charge is zero for EAM and TERSOFF type potential
                      ! Do nothing..
                   CASE DEFAULT
                      only_manybody = .FALSE.
                      EXIT
                   END SELECT
                   found_p = .TRUE.
                END IF
             END DO
             IF (only_manybody.AND.found_p) THEN
                charge(i) = 0.0_dp
                found     = .TRUE.
             END IF
          END IF
       END IF
       IF (.NOT.found) CALL store_FF_missing_par(atm1=TRIM(name_atm_a),&
                                                 fatal=fatal,&
                                                 type_name="Charge",&
                                                 array=Ainfo,&
                                                 error=error)
       !
       ! QM/MM modifications
       !
       IF (only_qm.AND.my_qmmm) THEN
          IF (qmmm_env%qmmm_coupl_type /= do_qmmm_none) THEN
             scale_factor = 0.0_dp
             IF (is_link_atom) THEN
                !
                ! Find the scaling factor...
                !
                DO ilink = 1, SIZE(qmmm_env%mm_link_atoms)
                   IF (ANY(my_atom_list == qmmm_env%mm_link_atoms(ilink))) EXIT
                END DO
                CPPostcondition(ilink <= SIZE(qmmm_env%mm_link_atoms),cp_failure_level,routineP,error,failure)
                scale_factor = qmmm_env%fist_scale_charge_link(ilink)
             END IF
             charge(i) = charge(i) * scale_factor
          END IF
       END IF

       CALL set_potential(potential=elp_potential,&
            qeff=charge(i))
       CALL set_atomic_kind(atomic_kind=atomic_kind,&
            elp_potential=elp_potential)

    END DO
    DEALLOCATE (charge,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_pack_ch

!******************************************************************************
!!****** force_fields_all/force_field_pack_shell [1.0] *
!!
!!   NAME
!!     force_field_pack_shell
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Set up shell potential parameters
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_pack_shell  (particle_set, atomic_kind_set,&
       molecule_kind_set, molecule_set, root_section, subsys_section,&
       shell_particle_set, core_particle_set, cell, iw, inp_info, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(section_vals_type), POINTER         :: root_section, subsys_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: shell_particle_set, &
                                                core_particle_set
    TYPE(cell_type), POINTER                 :: cell
    INTEGER                                  :: iw
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_shell', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atmname
    INTEGER                                  :: counter, first, first_shell, &
                                                handle2, i, j, last, &
                                                last_shell, n, natom, nmol, &
                                                nshell_tot, stat
    INTEGER, DIMENSION(:), POINTER           :: molecule_list, shell_list_tmp
    LOGICAL :: failure, found_shell, is_a_shell, is_link_atom, null_massfrac, &
      only_qm, save_mem, shell_adiabatic
    REAL(KIND=dp)                            :: atmmass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(section_vals_type), POINTER         :: global_section
    TYPE(shell_kind_type), POINTER           :: shell
    TYPE(shell_type), DIMENSION(:), POINTER  :: shell_list

    failure = .FALSE.

    CALL timeset(routineN,'I','',handle2)
    nshell_tot = 0
    n = 0
    first_shell = 1
    null_massfrac =.FALSE.

    NULLIFY(global_section)
    global_section => section_vals_get_subs_vals(root_section,"GLOBAL",error=error)
    CALL section_vals_val_get(global_section,"SAVE_MEM",l_val=save_mem,error=error)

    DO i=1,SIZE(atomic_kind_set)
       atomic_kind => atomic_kind_set(i)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            name=atmname)

       found_shell = .FALSE.
       only_qm = qmmm_ff_precond_only_qm(id1=atmname,is_link=is_link_atom)
       CALL uppercase(atmname)

       ! The shell potential can be defined only from input
       IF(ASSOCIATED(inp_info%shell_list)) THEN
          DO j=1,SIZE(inp_info%shell_list)
             IF (iw>0) WRITE(iw,*)"Shell Checking ::",TRIM(inp_info%shell_list(j)%atm_name),atmname

             IF((inp_info%shell_list(j)%atm_name)==atmname) THEN
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     shell=shell, mass=atmmass, natom=natom)
                IF(.NOT. ASSOCIATED(shell)) THEN
                   CALL shell_create(shell,error)
                END IF
                nshell_tot = nshell_tot + natom
                shell%charge_core=inp_info%shell_list(j)%shell%charge_core
                shell%charge_shell=inp_info%shell_list(j)%shell%charge_shell
                shell%massfrac=inp_info%shell_list(j)%shell%massfrac
                IF(shell%massfrac <EPSILON( 1.0_dp )) null_massfrac = .TRUE.
                shell%k_spring=inp_info%shell_list(j)%shell%k_spring
                shell%max_dist=inp_info%shell_list(j)%shell%max_dist
                shell%shell_cutoff=inp_info%shell_list(j)%shell%shell_cutoff
                
                shell%mass_shell=shell%massfrac*atmmass
                shell%mass_core=atmmass-shell%mass_shell
                CALL cp_assert(.NOT.found_shell,cp_warning_level,cp_assertion_failed,routineP,&
                     "Multiple Shell declarations: "//TRIM(atmname)//&
                     " overwriting! "//&
CPSourceFileRef,&
                     only_ionode=.TRUE.)
                found_shell = .TRUE.
                CALL set_atomic_kind(atomic_kind=atomic_kind,&
                     shell=shell, shell_active=.TRUE.,error=error)
                CALL shell_release(shell, error)
             END IF
          END DO  ! j shell kind
       END IF   ! associated shell_list
    END DO  ! i atomic kind

    IF (iw>0) WRITE(iw,*)"Total number of particles with a shell :: ", nshell_tot
    ! If shell-model is present: Create particle_set of shells (coord. vel. force)
    NULLIFY(shell_particle_set,core_particle_set)
    IF(nshell_tot>0) THEN
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, shell_adiabatic=shell_adiabatic)
       IF(shell_adiabatic .AND. null_massfrac) THEN
          CALL stop_program(TRIM(routineN),"Shell-model adiabatic: at least one shell_kind has mass=0 Abort..")
       END IF
       CALL allocate_particle_set(shell_particle_set,nshell_tot,error)
       CALL allocate_particle_set(core_particle_set,nshell_tot,error)
       counter = 0
       DO i = 1,SIZE(particle_set)
          NULLIFY(atomic_kind,shell)
          atomic_kind => particle_set(i)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind, shell_active=is_a_shell)
          IF(is_a_shell) THEN
             counter = counter+1
             particle_set(i)%shell_index = counter
             shell_particle_set(counter)%shell_index=counter
             shell_particle_set(counter)%atomic_kind => particle_set(i)%atomic_kind
             shell_particle_set(counter)%r(1:3) =  particle_set(i)%r(1:3)
             shell_particle_set(counter)%atom_index =   i
             core_particle_set(counter)%shell_index=counter
             core_particle_set(counter)%atomic_kind => particle_set(i)%atomic_kind
             core_particle_set(counter)%r(1:3) =  particle_set(i)%r(1:3)
             core_particle_set(counter)%atom_index =   i
          ELSE
             particle_set(i)%shell_index = 0
          END IF
       END DO
       CPPostcondition(counter==nshell_tot,cp_failure_level,routineP,error,failure)

      ! Initialize the shells coordinates
       IF(shell_adiabatic) THEN
          CALL read_shell_coord_input(particle_set, shell_particle_set, cell,&
               root_section, subsys_section, core_particle_set, save_mem=save_mem,&
               error=error)
       ELSE
          CALL read_shell_coord_input(particle_set, shell_particle_set, cell,&
               root_section, subsys_section,save_mem=save_mem, error=error)
       END IF

       ! Determine the number of shells per molecule kind
       n = 0
       DO i=1,SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(i)
          CALL get_molecule_kind(molecule_kind=molecule_kind,molecule_list=molecule_list,&
               natom=natom, nmolecule=nmol)
          molecule=>molecule_set(molecule_list(1))
          CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
          ALLOCATE(shell_list_tmp(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          counter = 0
          DO j = first, last
             atomic_kind => particle_set(j)%atomic_kind
             first_shell = MIN(first_shell, MAX(1,particle_set(j)%shell_index))
             CALL get_atomic_kind(atomic_kind=atomic_kind,shell_active=is_a_shell)
             IF(is_a_shell) THEN
                counter = counter+1
                shell_list_tmp(counter) = j - first + 1
             END IF
          END DO  ! j atom in molecule_kind i, molecule 1 of the molecule_list
          IF(counter/=0) THEN
             ! Setup of fist_shell and last_shell for all molecules..
             DO j = 1, SIZE(molecule_list)
                last_shell = first_shell + counter -1
                molecule=>molecule_set(molecule_list(j))
                molecule%first_shell = first_shell
                molecule%last_shell  = last_shell
                first_shell = last_shell + 1
             END DO
             ! Setup of shell_list
             CALL get_molecule_kind(molecule_kind=molecule_kind,shell_list=shell_list)
             IF(ASSOCIATED(shell_list)) THEN
                DEALLOCATE(shell_list,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             ALLOCATE(shell_list(counter), STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO j = 1, counter
                shell_list(j)%a = shell_list_tmp(j)
                atomic_kind => particle_set(shell_list_tmp(j)+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,name=atmname, shell=shell)
                CALL uppercase(atmname)
                shell_list(j)%name = atmname 
                shell_list(j)%shell_kind => shell
             END DO
             CALL set_molecule_kind(molecule_kind=molecule_kind, nshell=counter,shell_list=shell_list)
          END IF
          DEALLOCATE(shell_list_tmp,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          n = n + nmol * counter 
       END DO  ! i molecule kind
    END IF  ! nshell_tot>0
    CPPostcondition(first_shell-1==nshell_tot,cp_failure_level,routineP,error,failure)
    CPPostcondition(n==nshell_tot,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_pack_shell

!******************************************************************************
!!****** force_fields_all/force_field_pack_sbe [1.0] *
!!
!!   NAME
!!     force_field_pack_sbe
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!      Set spline_bond_env, use for ONFO interaction
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_pack_sbe  (atomic_kind_set, ff_type, qmmm_env,&
       iw, Ainfo,  chm_info, inp_info, gro_info, potparm_bond, ewald_env, &
       spline_bond_env, iw2, iw3, iw4, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(qmmm_env_mm_type), POINTER          :: qmmm_env
    INTEGER                                  :: iw
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: Ainfo
    TYPE(charmm_info_type), POINTER          :: chm_info
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(gromos_info_type), POINTER          :: gro_info
    TYPE(pair_potential_pp_type), POINTER    :: potparm_bond
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(spline_environment_type), POINTER   :: spline_bond_env
    INTEGER                                  :: iw2, iw3, iw4
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_sbe', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_a_local, &
                                                name_atm_b, name_atm_b_local
    INTEGER                                  :: handle2, i, ii, ikind, j, jj, &
                                                jkind, k, match_names, n
    LOGICAL                                  :: allocate_spline_bond_env, &
                                                failure, found, found_a, &
                                                found_b, only_qm, use_qmmm_ff
    REAL(KIND=dp)                            :: epsilon0, epsilon_a, &
                                                epsilon_b, rmin, rmin2_a, &
                                                rmin2_b
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(spline_data_p_type), DIMENSION(:), &
      POINTER                                :: spl_p

    failure = .FALSE.
    use_qmmm_ff = qmmm_env%use_qmmm_ff
    allocate_spline_bond_env = .FALSE.
    NULLIFY(pot)
    CALL timeset(routineN,'I','',handle2)
    CALL pair_potential_pp_create (potparm_bond, SIZE(atomic_kind_set), error)
    DO i=1,SIZE(atomic_kind_set)
       atomic_kind => atomic_kind_set(i)
       CALL get_atomic_kind(atomic_kind=atomic_kind,name=name_atm_a_local)
       DO j=i,SIZE(atomic_kind_set)
          atomic_kind => atomic_kind_set(j)
          CALL get_atomic_kind(atomic_kind=atomic_kind,name=name_atm_b_local)
          found = .FALSE.
          found_a = .FALSE.
          found_b = .FALSE.
          name_atm_a = name_atm_a_local
          name_atm_b = name_atm_b_local
          only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b)
          CALL uppercase(name_atm_a)
          CALL uppercase(name_atm_b)
          pot => potparm_bond%pot(i,j)%pot

          IF(ASSOCIATED(gro_info%nonbond_a_14)) THEN
             ii = 0 
             jj = 0
             DO k=1,SIZE(gro_info%nonbond_a_14)
                IF(TRIM(name_atm_a)==TRIM(gro_info%nonbond_a_14(k))) THEN
                   ii = k
                   found_a = .TRUE.
                   EXIT
                END IF
             END DO
             DO k=1,SIZE(gro_info%nonbond_a_14)
                IF(TRIM(name_atm_b)==TRIM(gro_info%nonbond_a_14(k))) THEN
                   jj = k
                   found_b   = .TRUE.
                   EXIT
                END IF
             END DO
             IF(ii/=0 .AND. jj/=0) THEN
                CALL  pair_potential_lj_create(pot%set(1)%lj, error)
                pot%type       = lj_type
                pot%at1        = name_atm_a
                pot%at2        = name_atm_b
                pot%set(1)%lj%epsilon = 1.0_dp
                pot%set(1)%lj%sigma6  = gro_info%nonbond_c6_14(ii,jj)
                pot%set(1)%lj%sigma12 = gro_info%nonbond_c12_14(ii,jj)
                pot%rcutsq     = (10.0_dp*bohr)**2
                CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                     "Multiple Lennard-Jones declarations: "//TRIM(name_atm_a)//&
                     " and "//TRIM(name_atm_b)//" overwriting! "//&
CPSourceFileRef,&
                     only_ionode=.TRUE.)
                found = .TRUE.
             END IF
          END IF

          ! loop over params from CHARMM
          ii = 0
          jj = 0
          IF(ASSOCIATED(chm_info%nonbond_a_14)) THEN
             DO k=1,SIZE(chm_info%nonbond_a_14)
                IF((name_atm_a)==(chm_info%nonbond_a_14(k))) THEN
                   ii = k
                   rmin2_a = chm_info%nonbond_rmin2_14(k)
                   epsilon_a = chm_info%nonbond_eps_14(k)
                   found_a = .TRUE.
                END IF
             END DO
             DO k=1,SIZE(chm_info%nonbond_a_14)
                IF((name_atm_b)==(chm_info%nonbond_a_14(k))) THEN
                   jj = k
                   rmin2_b   = chm_info%nonbond_rmin2_14(k)
                   epsilon_b = chm_info%nonbond_eps_14(k)
                   found_b   = .TRUE.
                END IF
             END DO
          END IF
          IF(ASSOCIATED(chm_info%nonbond_a)) THEN
             IF(.NOT.found_a) THEN
                DO k=1,SIZE(chm_info%nonbond_a)
                   IF((name_atm_a)==(chm_info%nonbond_a(k))) THEN
                      ii = k
                      rmin2_a   = chm_info%nonbond_rmin2(k)
                      epsilon_a = chm_info%nonbond_eps(k)
                   END IF
                END DO
             END IF
             IF(.NOT.found_b) THEN
                DO k=1,SIZE(chm_info%nonbond_a)
                   IF((name_atm_b)==(chm_info%nonbond_a(k))) THEN
                      jj = k
                      rmin2_b   = chm_info%nonbond_rmin2(k)
                      epsilon_b = chm_info%nonbond_eps(k)
                   END IF
                END DO
             END IF
          END IF
          IF(ii/=0 .AND. jj/=0) THEN
             rmin = rmin2_a + rmin2_b
             ! ABS to allow for mixing the two different sign conventions for epsilon
             epsilon0 = SQRT(ABS(epsilon_a*epsilon_b))
             CALL  pair_potential_lj_create(pot%set(1)%lj, error)
             pot%type       = lj_charmm_type
             pot%at1        = name_atm_a
             pot%at2        = name_atm_b
             pot%set(1)%lj%epsilon = epsilon0
             pot%set(1)%lj%sigma6  = 0.5_dp*rmin**6
             pot%set(1)%lj%sigma12 = 0.25_dp*rmin**12
             pot%rcutsq     = (10.0_dp*bohr)**2
             CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                  "Multiple Lennard-Jones declarations: "//TRIM(name_atm_a)//&
                  " and "//TRIM(name_atm_b)//" overwriting! "//&
CPSourceFileRef,&
                  only_ionode=.TRUE.)
             found = .TRUE.
          END IF

          ! always have the input param last to overwrite all the other ones
          IF(ASSOCIATED(inp_info%nonbonded14)) THEN
             DO k=1,SIZE(inp_info%nonbonded14%pot)
                IF(iw>0) WRITE(iw,*) "    TESTING ",TRIM(name_atm_a),TRIM(name_atm_b),&
                     " with ",TRIM(inp_info%nonbonded14%pot(k)%pot%at1),&
                     TRIM(inp_info%nonbonded14%pot(k)%pot%at2)
                IF((((name_atm_a)==(inp_info%nonbonded14%pot(k)%pot%at1))  .AND. &
                     ((name_atm_b)==(inp_info%nonbonded14%pot(k)%pot%at2))) .OR.  &
                     (((name_atm_b)==(inp_info%nonbonded14%pot(k)%pot%at1))  .AND. &
                     ((name_atm_a)==(inp_info%nonbonded14%pot(k)%pot%at2))) ) THEN
                   IF (ff_type%multiple_potential) THEN
                      CALL pair_potential_single_add(inp_info%nonbonded14%pot(k)%pot,pot,error)
                      CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                           "Multiple ONFO declaration: "//TRIM(name_atm_a)//&
                           " and "//TRIM(name_atm_b)//" ADDING! "//&
CPSourceFileRef,&
                           only_ionode=.TRUE.)
                      potparm_bond%pot(i,j)%pot => pot
                      potparm_bond%pot(j,i)%pot => pot
                   ELSE
                      CALL pair_potential_single_copy(inp_info%nonbonded14%pot(k)%pot,pot,error)
                      CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                           "Multiple ONFO declarations: "//TRIM(name_atm_a)//&
                           " and "//TRIM(name_atm_b)//" OVERWRITING! "//&
CPSourceFileRef,&
                           only_ionode=.TRUE.)
                   END IF
                   IF(iw>0) WRITE(iw,*) "    FOUND ",TRIM(name_atm_a)," ",TRIM(name_atm_b)
                   found = .TRUE.
                END IF
             END DO
          END IF
          ! at the very end we offer the possibility to overwrite the parameters for QM/MM
          ! nonbonded interactions
          IF (use_qmmm_ff) THEN
             match_names = 0
             IF ((name_atm_a) == (name_atm_a_local)) match_names = match_names + 1
             IF ((name_atm_b) == (name_atm_b_local)) match_names = match_names + 1
             IF (match_names == 1) THEN
                IF (ASSOCIATED(qmmm_env%inp_info%nonbonded14)) THEN
                   DO k=1,SIZE(qmmm_env%inp_info%nonbonded14%pot)
                      IF(iw>0) WRITE(iw,*) "    TESTING ",TRIM(name_atm_a),TRIM(name_atm_b),&
                           " with ",TRIM(qmmm_env%inp_info%nonbonded14%pot(k)%pot%at1),&
                           TRIM(qmmm_env%inp_info%nonbonded14%pot(k)%pot%at2)
                      IF(( ((name_atm_a) ==(qmmm_env%inp_info%nonbonded14%pot(k)%pot%at1))  .AND. &
                           ((name_atm_b) ==(qmmm_env%inp_info%nonbonded14%pot(k)%pot%at2))) .OR.  &
                           (((name_atm_b)==(qmmm_env%inp_info%nonbonded14%pot(k)%pot%at1))  .AND. &
                           ((name_atm_a) ==(qmmm_env%inp_info%nonbonded14%pot(k)%pot%at2))) ) THEN
                         IF (qmmm_env%multiple_potential) THEN
                            CALL pair_potential_single_add(qmmm_env%inp_info%nonbonded14%pot(k)%pot,pot,error)
                            CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                                 "Multiple ONFO declaration: "//TRIM(name_atm_a)//&
                                 " and "//TRIM(name_atm_b)//" ADDING QM/MM forcefield specifications! "//&
CPSourceFileRef,&
                                 only_ionode=.TRUE.)
                            potparm_bond%pot(i,j)%pot => pot
                            potparm_bond%pot(j,i)%pot => pot
                         ELSE
                            CALL pair_potential_single_copy(qmmm_env%inp_info%nonbonded14%pot(k)%pot,pot,error)
                            CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                                 "Multiple ONFO declaration: "//TRIM(name_atm_a)//&
                                 " and "//TRIM(name_atm_b)//" OVERWRITING QM/MM forcefield specifications! "//&
CPSourceFileRef,&
                                 only_ionode=.TRUE.)
                         END IF
                         IF(iw>0) WRITE(iw,*) "    FOUND ",TRIM(name_atm_a),&
                              " ",TRIM(name_atm_b)
                         found = .TRUE.
                      END IF
                   END DO
                END IF
             END IF
          END IF

          IF(.NOT.found) THEN
             CALL store_FF_missing_par(atm1=TRIM(name_atm_a),&
                                       atm2=TRIM(name_atm_b),&
                                       type_name="Spline_Bond_Env",&
                                       array=Ainfo,&
                                       error=error)
             CALL pair_potential_single_clean( pot ,error=error)
             pot%type       = nn_type
             pot%at1        = name_atm_a
             pot%at2        = name_atm_b
          ELSE
             allocate_spline_bond_env = .TRUE.
          END IF
          ! If defined global RCUT let's use it
          IF(ff_type%rcut_nb>0.0_dp) THEN
             pot%rcutsq     = ff_type%rcut_nb*ff_type%rcut_nb
          END IF
          IF (only_qm) THEN
             CALL pair_potential_single_clean( pot ,error=error)
          END IF
       END DO  ! atom kind j
    END DO  ! atom kind i

    ! Allocate Spline Bond Env
    CALL get_nonbond_storage(spline_bond_env, potparm_bond, atomic_kind_set, &
          ei_scale14=ff_type%ei_scale14, vdw_scale14=ff_type%vdw_scale14,&
          do_14=.TRUE., ewald_env=ewald_env, scale_cutoff=ff_type%scale_cutoff,&
          error=error)
    CALL spline_nonbond_control ( spline_bond_env, potparm_bond, &
         atomic_kind_set , ewald_env, ei_scale14=ff_type%ei_scale14, &
         vdw_scale14=ff_type%vdw_scale14, eps_spline=ff_type%eps_spline, &
         max_energy=ff_type%max_energy, rlow_nb=ff_type%rlow_nb, &
         emax_spline=ff_type%emax_spline, iw=iw2, iw2=iw3, iw3=iw4, &
         do_14=.TRUE., scale_cutoff=ff_type%scale_cutoff, error=error)

    DO ikind = 1, SIZE ( potparm_bond%pot, 1 )
       DO jkind = ikind, SIZE ( potparm_bond%pot, 2 )
          n = spline_bond_env % spltab ( ikind, jkind )
          spl_p => spline_bond_env%spl_pp(n)%spl_p
          CALL spline_data_p_retain(spl_p, error)
          CALL spline_data_p_release(potparm_bond%pot(ikind,jkind)%pot%pair_spline_data,error)
          potparm_bond%pot(ikind,jkind)%pot%pair_spline_data => spl_p
       END DO
    END DO
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_pack_sbe

!******************************************************************************
!!****** force_fields_all/force_field_pack_onfo [1.0] *
!!
!!   NAME
!!     force_field_pack_onfo 
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Pack in ONFO information needed for the force_field
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_pack_onfo  (particle_set, molecule_kind_set, molecule_set,&
       iw, potparm_bond, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    INTEGER                                  :: iw
    TYPE(pair_potential_pp_type), POINTER    :: potparm_bond
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_onfo', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_b
    INTEGER                                  :: atm_a, atm_b, first, handle2, &
                                                i, ikind, j, jkind, last, &
                                                natom, nonfo
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    LOGICAL                                  :: failure, only_qm
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(onfo_type), DIMENSION(:), POINTER   :: onfo_list

    CALL timeset(routineN,'I','',handle2)
    failure = .FALSE.

    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            nonfo=nonfo,&
            onfo_list=onfo_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)

       DO j=1,nonfo
          atm_a = onfo_list(j)%a
          atomic_kind => particle_set(atm_a+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_a, kind_number = ikind )
          atm_b = onfo_list(j)%b
          atomic_kind => particle_set(atm_b+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_b, kind_number = jkind )
          only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b)
          CALL uppercase(name_atm_a)
          CALL uppercase(name_atm_b)

          CALL spline_data_p_retain (potparm_bond%pot(ikind,jkind)%pot%pair_spline_data,error=error)
          CALL spline_data_p_release(onfo_list(j)%onfo_kind%spl_p, error=error)
          onfo_list(j)%onfo_kind%spl_p => potparm_bond%pot(ikind,jkind)%pot%pair_spline_data
          CALL spline_factor_copy(potparm_bond%pot(ikind,jkind)%pot%spl_f,onfo_list(j)%onfo_kind%spl_f,error)

          IF(iw>0) WRITE(iw,*) "    ONE-FOUR ",j," : ",&
               TRIM(name_atm_a)," ",&
               TRIM(name_atm_b)
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            onfo_list=onfo_list)
    END DO
    CALL pair_potential_pp_release ( potparm_bond, error )
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_pack_onfo

!******************************************************************************
!!****** force_fields_all/force_field_pack_sne [1.0] *
!!
!!   NAME
!!     force_field_pack_sne
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Pack in spline nonbond env the information needed for the force_field
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_pack_sne (atomic_kind_set, ff_type, qmmm_env,&
       fatal, iw, Ainfo, chm_info, inp_info, gro_info, potparm_nonbond, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(qmmm_env_mm_type), POINTER          :: qmmm_env
    LOGICAL                                  :: fatal
    INTEGER                                  :: iw
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: Ainfo
    TYPE(charmm_info_type), POINTER          :: chm_info
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(gromos_info_type), POINTER          :: gro_info
    TYPE(pair_potential_pp_type), POINTER    :: potparm_nonbond
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_sne', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_a_local, &
                                                name_atm_b, name_atm_b_local
    INTEGER                                  :: handle2, i, ii, j, jj, k, &
                                                match_names
    LOGICAL                                  :: failure, found, is_a_shell, &
                                                is_b_shell, only_qm, &
                                                use_qmmm_ff
    REAL(KIND=dp)                            :: epsilon0, rmin
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot

    CALL timeset(routineN,'I','',handle2)
    failure = .FALSE.
    use_qmmm_ff = qmmm_env%use_qmmm_ff
    NULLIFY(pot)
    
    CALL pair_potential_pp_create ( potparm_nonbond, SIZE(atomic_kind_set),error )
    DO i=1,SIZE(atomic_kind_set)
       atomic_kind => atomic_kind_set(i)
       CALL get_atomic_kind(atomic_kind=atomic_kind,name=name_atm_a_local,&
           shell_active=is_a_shell)
       DO j=i,SIZE(atomic_kind_set)
          atomic_kind => atomic_kind_set(j)
          CALL get_atomic_kind(atomic_kind=atomic_kind,name=name_atm_b_local,&
               shell_active=is_b_shell)
          found = .FALSE.
          name_atm_a = name_atm_a_local
          name_atm_b = name_atm_b_local
          only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b)
          CALL uppercase(name_atm_a)
          CALL uppercase(name_atm_b)
          pot => potparm_nonbond%pot(i,j)%pot

          ! loop over params from GROMOS
          IF(ASSOCIATED(gro_info%nonbond_a)) THEN
             ii = 0 
             jj = 0
             DO k=1,SIZE(gro_info%nonbond_a)
                IF(TRIM(name_atm_a)==TRIM(gro_info%nonbond_a(k))) THEN
                   ii = k
                   EXIT
                END IF
             END DO
             DO k=1,SIZE(gro_info%nonbond_a)
                IF(TRIM(name_atm_b)==TRIM(gro_info%nonbond_a(k))) THEN
                   jj = k
                   EXIT
                END IF
             END DO

             IF(ii/=0 .AND. jj/=0) THEN
                CALL  pair_potential_lj_create(pot%set(1)%lj, error)
                pot%type       = lj_type
                pot%at1        = name_atm_a
                pot%at2        = name_atm_b
                pot%set(1)%lj%epsilon = 1.0_dp
                pot%set(1)%lj%sigma6  = gro_info%nonbond_c6(ii,jj)
                pot%set(1)%lj%sigma12 = gro_info%nonbond_c12(ii,jj)
                pot%rcutsq     = (10.0_dp*bohr)**2
                CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                     "Multiple Lennard-Jones declarations: "//TRIM(name_atm_a)//&
                     " and "//TRIM(name_atm_b)//" overwriting! "//&
CPSourceFileRef,&
                     only_ionode=.TRUE.)
                found = .TRUE.
             END IF
          END IF

          ! loop over params from CHARMM
          IF(ASSOCIATED(chm_info%nonbond_a)) THEN
             ii = 0
             jj = 0
             DO k=1,SIZE(chm_info%nonbond_a)
                IF((name_atm_a)==(chm_info%nonbond_a(k))) THEN
                   ii = k
                END IF
             END DO
             DO k=1,SIZE(chm_info%nonbond_a)
                IF((name_atm_b)==(chm_info%nonbond_a(k))) THEN
                   jj = k
                END IF
             END DO

             IF(ii/=0 .AND. jj/=0) THEN
                rmin = chm_info%nonbond_rmin2(ii)+chm_info%nonbond_rmin2(jj)
                epsilon0 = SQRT(chm_info%nonbond_eps(ii)*&
                     chm_info%nonbond_eps(jj))
                CALL  pair_potential_lj_create(pot%set(1)%lj, error)
                pot%type       = lj_charmm_type
                pot%at1        = name_atm_a
                pot%at2        = name_atm_b
                pot%set(1)%lj%epsilon = epsilon0
                pot%set(1)%lj%sigma6  = 0.5_dp*rmin**6
                pot%set(1)%lj%sigma12 = 0.25_dp*rmin**12
                pot%rcutsq     = (10.0_dp*bohr)**2
                CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                     "Multiple Lennard-Jones declarations: "//TRIM(name_atm_a)//&
                     " and "//TRIM(name_atm_b)//" overwriting! "//&
CPSourceFileRef,&
                     only_ionode=.TRUE.)
                found = .TRUE.
             END IF
          END IF
          
          ! always have the input param last to overwrite all the other ones
          IF(ASSOCIATED(inp_info%nonbonded)) THEN
             DO k=1,SIZE(inp_info%nonbonded%pot)
                IF(iw>0) WRITE(iw,*) "    TESTING ",TRIM(name_atm_a),TRIM(name_atm_b),&
                     " with ",TRIM(inp_info%nonbonded%pot(k)%pot%at1),&
                     TRIM(inp_info%nonbonded%pot(k)%pot%at2)
                IF(( ((name_atm_a)==(inp_info%nonbonded%pot(k)%pot%at1))  .AND. &
                     ((name_atm_b)==(inp_info%nonbonded%pot(k)%pot%at2))) .OR.  &
                    (((name_atm_b)==(inp_info%nonbonded%pot(k)%pot%at1))  .AND. &
                     ((name_atm_a)==(inp_info%nonbonded%pot(k)%pot%at2))) ) THEN
                   IF (ff_type%multiple_potential) THEN
                      CALL pair_potential_single_add(inp_info%nonbonded%pot(k)%pot,pot,error)
                      CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                           "Multiple NONBONDED declaration: "//TRIM(name_atm_a)//&
                           " and "//TRIM(name_atm_b)//" ADDING! "//&
CPSourceFileRef,&
                           only_ionode=.TRUE.)
                      potparm_nonbond%pot(i,j)%pot => pot
                      potparm_nonbond%pot(j,i)%pot => pot
                   ELSE
                      CALL pair_potential_single_copy(inp_info%nonbonded%pot(k)%pot,pot,error)
                      CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                           "Multiple NONBONDED declaration: "//TRIM(name_atm_a)//&
                           " and "//TRIM(name_atm_b)//" OVERWRITING! "//&
CPSourceFileRef,&
                           only_ionode=.TRUE.)
                   END IF
                   IF(iw>0) WRITE(iw,*) "    FOUND ",TRIM(name_atm_a)," ",TRIM(name_atm_b)
                   found = .TRUE.
                END IF
             END DO
          END IF
          ! at the very end we offer the possibility to overwrite the parameters for QM/MM
          ! nonbonded interactions
          IF (use_qmmm_ff) THEN
             match_names = 0
             IF ((name_atm_a) == (name_atm_a_local)) match_names = match_names + 1
             IF ((name_atm_b) == (name_atm_b_local)) match_names = match_names + 1
             IF (match_names == 1) THEN
                IF (ASSOCIATED(qmmm_env%inp_info%nonbonded)) THEN
                   DO k=1,SIZE(qmmm_env%inp_info%nonbonded%pot)
                      IF(iw>0) WRITE(iw,*) "    TESTING ",TRIM(name_atm_a),TRIM(name_atm_b),&
                           " with ",TRIM(qmmm_env%inp_info%nonbonded%pot(k)%pot%at1),&
                                    TRIM(qmmm_env%inp_info%nonbonded%pot(k)%pot%at2)
                      IF(( ((name_atm_a) ==(qmmm_env%inp_info%nonbonded%pot(k)%pot%at1))  .AND. &
                           ((name_atm_b) ==(qmmm_env%inp_info%nonbonded%pot(k)%pot%at2))) .OR.  &
                           (((name_atm_b)==(qmmm_env%inp_info%nonbonded%pot(k)%pot%at1))  .AND. &
                           ((name_atm_a) ==(qmmm_env%inp_info%nonbonded%pot(k)%pot%at2))) ) THEN
                         IF (qmmm_env%multiple_potential) THEN
                            CALL pair_potential_single_add(qmmm_env%inp_info%nonbonded%pot(k)%pot,pot,error)
                            CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                                 "Multiple NONBONDED declaration: "//TRIM(name_atm_a)//&
                                 " and "//TRIM(name_atm_b)//" ADDING QM/MM forcefield specifications! "//&
CPSourceFileRef,&
                                 only_ionode=.TRUE.)
                            potparm_nonbond%pot(i,j)%pot => pot
                            potparm_nonbond%pot(j,i)%pot => pot
                         ELSE
                            CALL pair_potential_single_copy(qmmm_env%inp_info%nonbonded%pot(k)%pot,pot,error)
                            CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                                 "Multiple NONBONDED declaration: "//TRIM(name_atm_a)//&
                                 " and "//TRIM(name_atm_b)//" OVERWRITING QM/MM forcefield specifications! "//&
CPSourceFileRef,&
                                 only_ionode=.TRUE.)
                         END IF
                         IF(iw>0) WRITE(iw,*) "    FOUND ",TRIM(name_atm_a)," ",TRIM(name_atm_b)
                         found = .TRUE.
                      END IF
                   END DO
                END IF
             END IF
          END IF
          IF(.NOT.found) THEN
             CALL store_FF_missing_par(atm1=TRIM(name_atm_a),&
                                       atm2=TRIM(name_atm_b),&
                                       type_name="Spline_Non_Bond_Env",&
                                       fatal=fatal,&
                                       array=Ainfo,&
                                       error=error)
          END IF
          ! If defined global RCUT let's use it
          IF(ff_type%rcut_nb>0.0_dp) THEN
             pot%rcutsq     = ff_type%rcut_nb*ff_type%rcut_nb
          END IF
          IF((is_a_shell .AND. .NOT. is_b_shell) .OR. (is_b_shell .AND. .NOT. is_a_shell)) THEN
             pot%shell_type  = nosh_sh
          ELSE IF (is_a_shell .AND. is_b_shell) THEN
             pot%shell_type  = sh_sh
          ELSE
             pot%shell_type  = nosh_nosh
          END IF
          IF (only_qm) THEN
             CALL pair_potential_single_clean( pot ,error=error)
          END IF
       END DO  ! jkind
    END DO   ! ikind
    CALL timestop(0.0_dp,handle2)
  END SUBROUTINE force_field_pack_sne

!******************************************************************************
!!****** force_fields_all/force_field_pack_ppse [1.0] *
!!
!!   NAME
!!     force_field_pack_ppse
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     create the pair potential spline environment
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_pack_ppse (atomic_kind_set, ewald_env, ff_type, &
       iw2, iw3, iw4, spline_nonbond_env, spline_c_coul_env, spline_cs_coul_env,&
       spline_sc_coul_env, potparm_nonbond, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    INTEGER                                  :: iw2, iw3, iw4
    TYPE(spline_environment_type), POINTER   :: spline_nonbond_env, &
                                                spline_c_coul_env, &
                                                spline_cs_coul_env, &
                                                spline_sc_coul_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm_nonbond
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_ppse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle2
    LOGICAL                                  :: failure, shell_present

    CALL timeset(routineN,'I','',handle2)
    failure = .FALSE.
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, shell_present=shell_present)
    IF(shell_present) THEN
      CALL get_nonbond_storage(spline_nonbond_env, potparm_nonbond, atomic_kind_set, &
           spline_c_coul_env=spline_c_coul_env, spline_cs_coul_env=spline_cs_coul_env,&
           spline_sc_coul_env=spline_sc_coul_env, do_14=.FALSE., ewald_env=ewald_env, &
           scale_cutoff=ff_type%scale_cutoff, error=error)
      CALL spline_nonbond_control(spline_nonbond_env, potparm_nonbond, &
           atomic_kind_set,ewald_env,eps_spline=ff_type%eps_spline,&
           max_energy=ff_type%max_energy, rlow_nb=ff_type%rlow_nb,&
           emax_spline=ff_type%emax_spline, iw=iw2, iw2=iw3,iw3=iw4,&
           scale_cutoff=ff_type%scale_cutoff, spline_c_coul=spline_c_coul_env, &
           spline_cs_coul=spline_cs_coul_env, spline_sc_coul=spline_sc_coul_env, error=error)
    ELSE
      CALL get_nonbond_storage(spline_nonbond_env, potparm_nonbond, atomic_kind_set, &
           do_14=.FALSE., ewald_env=ewald_env, scale_cutoff=ff_type%scale_cutoff,error=error)
      NULLIFY (spline_c_coul_env,spline_cs_coul_env,spline_sc_coul_env) 
      CALL spline_nonbond_control(spline_nonbond_env, potparm_nonbond, &
           atomic_kind_set,ewald_env,eps_spline=ff_type%eps_spline,&
           max_energy=ff_type%max_energy, rlow_nb=ff_type%rlow_nb,&
           emax_spline=ff_type%emax_spline, iw=iw2, iw2=iw3,iw3=iw4,&
           scale_cutoff=ff_type%scale_cutoff, error=error)
    END IF
    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_pack_ppse

!******************************************************************************
!!****** force_fields_all/force_field_pack_pp [1.0] *
!!
!!   NAME
!!     force_field_pack_pp
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Setup the pair potential spline environment
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE  force_field_pack_pp (atomic_kind_set, spline_nonbond_env, spline_c_coul_env,&
       spline_cs_coul_env, spline_sc_coul_env, potparm_nonbond, shell_present, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(spline_environment_type), POINTER   :: spline_nonbond_env, &
                                                spline_c_coul_env, &
                                                spline_cs_coul_env, &
                                                spline_sc_coul_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm_nonbond
    LOGICAL                                  :: shell_present
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_pp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle2, ikind, jkind, n, &
                                                nkinds
    LOGICAL                                  :: failure
    TYPE(spline_data_p_type), DIMENSION(:), &
      POINTER                                :: spl_p

    CALL timeset(routineN,'I','',handle2)
    failure = .FALSE.

    nkinds = SIZE(atomic_kind_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, shell_present=shell_present)

    DO ikind = 1, SIZE ( potparm_nonbond%pot, 1 )
       DO jkind = ikind, SIZE ( potparm_nonbond%pot, 2)
          !---------------------------------------
          ! 1st: The pair-potential contribution, 
          !      if shell-model, it is ion-shell or shell-shell 
          !---------------------------------------
          n = spline_nonbond_env % spltab ( ikind, jkind )
          spl_p => spline_nonbond_env%spl_pp(n)%spl_p
          CALL spline_data_p_retain ( spl_p, error )
          CALL spline_data_p_release ( potparm_nonbond%pot(ikind,jkind)%pot%pair_spline_data, error )
          potparm_nonbond%pot(ikind,jkind)%pot%pair_spline_data => spl_p

          IF(shell_present) THEN
             !---------------------------------------
             ! 2nd: The ion-core or core-core contribution (only shell-model)
             !---------------------------------------
             n = spline_c_coul_env % spltab ( ikind, jkind )
             IF (n>0) THEN
                spl_p => spline_c_coul_env%spl_pp(n)%spl_p
                ! do i,j
                IF ( ASSOCIATED(spl_p) ) THEN
                   CALL spline_data_p_retain( spl_p, error )
                   CALL spline_data_p_release(potparm_nonbond%pot( ikind,jkind)%pot%c_coul_spline_data, error )
                   potparm_nonbond%pot ( ikind, jkind )%pot % c_coul_spline_data => spl_p
                ENDIF
             END IF
             !---------------------------------------
             ! 3rd: The core-shell contribution (only shell-model)
             !---------------------------------------
             n = spline_cs_coul_env % spltab ( ikind, jkind )
             IF (n>0) THEN
                spl_p => spline_cs_coul_env%spl_pp(n)%spl_p
                ! do i,j
                IF ( ASSOCIATED(spl_p) ) THEN
                   CALL spline_data_p_retain(spl_p, error )
                   CALL spline_data_p_release(potparm_nonbond%pot(ikind,jkind)%pot%cs_coul_spline_data, error )
                   potparm_nonbond%pot ( ikind, jkind )%pot % cs_coul_spline_data => spl_p
                ENDIF
             END IF
             !---------------------------------------
             ! 4th: The core-shell contribution (only shell-model)
             !---------------------------------------
             n = spline_sc_coul_env % spltab ( ikind, jkind )
             IF (n>0) THEN
                spl_p => spline_sc_coul_env%spl_pp(n)%spl_p
                ! do i,j
                IF ( ASSOCIATED(spl_p) ) THEN
                   CALL spline_data_p_retain (spl_p, error )
                   CALL spline_data_p_release(potparm_nonbond%pot(ikind,jkind)%pot%sc_coul_spline_data, error )
                   potparm_nonbond%pot(ikind,jkind)%pot%sc_coul_spline_data => spl_p
                ENDIF
             END IF
          END IF  ! shell_present
       END DO
    END DO

    CALL timestop(0.0_dp,handle2)

  END SUBROUTINE force_field_pack_pp

!******************************************************************************
!!****** force_fields_all/store_FF_missing_par [1.0] *
!!
!!   NAME
!!     store_FF_missing_par
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Store informations on possible missing ForceFields parameters
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************
  SUBROUTINE store_FF_missing_par(atm1,atm2,atm3,atm4,type_name,fatal,array,error)
    CHARACTER(LEN=*), INTENT(IN)             :: atm1
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: atm2, atm3, atm4
    CHARACTER(LEN=*), INTENT(IN)             :: type_name
    LOGICAL, INTENT(INOUT), OPTIONAL         :: fatal
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: array
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'store_FF_missing_par', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: sfmt
    CHARACTER(LEN=4)                         :: my_atm1, my_atm2, my_atm3, &
                                                my_atm4
    CHARACTER(LEN=default_path_length)       :: my_format
    INTEGER                                  :: fmt, i, nsize
    LOGICAL                                  :: failure, found

    failure = .FALSE.
    nsize   = 0
    fmt     = 1
    my_format = '(T2,"FORCEFIELD| Missing ","'//TRIM(type_name)//&
                '",T40,"(",A4,")")'
    IF (PRESENT(atm2)) fmt = fmt+1
    IF (PRESENT(atm3)) fmt = fmt+1
    IF (PRESENT(atm4)) fmt = fmt+1
    CALL integer_to_string(fmt-1,sfmt)
    IF (fmt>1) &
    my_format = '(T2,"FORCEFIELD| Missing ","'//TRIM(type_name)//&
                '",T40,"(",A4,'//TRIM(sfmt)//'(",",A4),")")'
    IF (.NOT.failure) THEN
       IF (PRESENT(fatal)) fatal = .TRUE.
       ! Check for previous already stored equal force fields
       IF (ASSOCIATED(array)) nsize = SIZE(array)
       found = .FALSE.
       IF (nsize>=1) THEN
          DO i = 1, nsize
             SELECT CASE (type_name)
             CASE("Bond")
                IF (INDEX(array(i)(21:39),"Bond")==0)     CYCLE
                my_atm1 = array(i)(41:44)
                my_atm2 = array(i)(46:49)
                CALL compress(my_atm1,.TRUE.)
                CALL compress(my_atm2,.TRUE.)
    IF  (((atm1==my_atm1).AND.(atm2==my_atm2)).OR.&
         ((atm1==my_atm2).AND.(atm2==my_atm1))) found = .TRUE.
             CASE("Angle")
                IF (INDEX(array(i)(21:39),"Angle")==0)    CYCLE
                my_atm1 = array(i)(41:44)
                my_atm2 = array(i)(46:49)
                my_atm3 = array(i)(51:54)
                CALL compress(my_atm1,.TRUE.)
                CALL compress(my_atm2,.TRUE.)
                CALL compress(my_atm3,.TRUE.)
    IF  (((atm1==my_atm1).AND.(atm2==my_atm2).AND.(atm3==my_atm3)).OR.&
         ((atm1==my_atm3).AND.(atm2==my_atm2).AND.(atm3==my_atm1)))&
         found = .TRUE.
             CASE("Urey-Bradley")
                IF (INDEX(array(i)(21:39),"Urey-Bradley")==0)    CYCLE
                my_atm1 = array(i)(41:44)
                my_atm2 = array(i)(46:49)
                my_atm3 = array(i)(51:54)
                CALL compress(my_atm1,.TRUE.)
                CALL compress(my_atm2,.TRUE.)
                CALL compress(my_atm3,.TRUE.)
    IF  (((atm1==my_atm1).AND.(atm2==my_atm2).AND.(atm3==my_atm3)).OR.&
         ((atm1==my_atm3).AND.(atm2==my_atm2).AND.(atm3==my_atm1)))&
         found = .TRUE.
             CASE("Torsion")
                IF (INDEX(array(i)(21:39),"Torsion")==0)  CYCLE
                my_atm1 = array(i)(41:44)
                my_atm2 = array(i)(46:49)
                my_atm3 = array(i)(51:54)
                my_atm4 = array(i)(56:59)
                CALL compress(my_atm1,.TRUE.)
                CALL compress(my_atm2,.TRUE.)
                CALL compress(my_atm3,.TRUE.)
                CALL compress(my_atm4,.TRUE.)
    IF  (((atm1==my_atm1).AND.(atm2==my_atm2).AND.(atm3==my_atm3).AND.(atm4==my_atm4)).OR.&
         ((atm1==my_atm4).AND.(atm2==my_atm3).AND.(atm3==my_atm2).AND.(atm4==my_atm1)))&
         found = .TRUE.
             CASE("Improper")
                IF (INDEX(array(i)(21:39),"Improper")==0) CYCLE
                my_atm1 = array(i)(41:44)
                my_atm2 = array(i)(46:49)
                my_atm3 = array(i)(51:54)
                my_atm4 = array(i)(56:59)
                CALL compress(my_atm1,.TRUE.)
                CALL compress(my_atm2,.TRUE.)
                CALL compress(my_atm3,.TRUE.)
                CALL compress(my_atm4,.TRUE.)
    IF  (((atm1==my_atm1).AND.(atm2==my_atm2).AND.(atm3==my_atm3).AND.(atm4==my_atm4)).OR.&
         ((atm1==my_atm1).AND.(atm2==my_atm3).AND.(atm3==my_atm2).AND.(atm4==my_atm4)).OR.&
         ((atm1==my_atm1).AND.(atm2==my_atm3).AND.(atm3==my_atm4).AND.(atm4==my_atm3)).OR.&
         ((atm1==my_atm1).AND.(atm2==my_atm4).AND.(atm3==my_atm3).AND.(atm4==my_atm2)).OR.&
         ((atm1==my_atm1).AND.(atm2==my_atm4).AND.(atm3==my_atm2).AND.(atm4==my_atm3)).OR.&
         ((atm1==my_atm1).AND.(atm2==my_atm2).AND.(atm3==my_atm4).AND.(atm4==my_atm3)))&
         found = .TRUE.
             CASE("Charge")
                IF (INDEX(array(i)(21:39),"Charge")==0)   CYCLE
                my_atm1 = array(i)(41:44)
                CALL compress(my_atm1,.TRUE.)
    IF  (atm1==my_atm1) found = .TRUE.
             CASE("Spline_Bond_Env","Spline_Non_Bond_Env")
                IF (INDEX(array(i)(21:39),"Spline_")==0)  CYCLE
                fmt = 0
                my_atm1 = array(i)(41:44)
                my_atm2 = array(i)(46:49)
                CALL compress(my_atm1,.TRUE.)
                CALL compress(my_atm2,.TRUE.)
    IF  (((atm1==my_atm1).AND.(atm2==my_atm2)).OR.&
         ((atm1==my_atm2).AND.(atm2==my_atm1))) found = .TRUE.
             CASE DEFAULT
                ! Should never reach this point
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
             IF (found) EXIT
          END DO
       ENDIF
       IF (.NOT.found) THEN
          nsize = nsize + 1
          CALL reallocate(array,1,nsize)
          SELECT CASE(fmt)
          CASE(1)
             WRITE(array(nsize),FMT=TRIM(my_format))atm1
          CASE(2)
             WRITE(array(nsize),FMT=TRIM(my_format))atm1,atm2
          CASE(3)
             WRITE(array(nsize),FMT=TRIM(my_format))atm1,atm2,atm3
          CASE(4)
             WRITE(array(nsize),FMT=TRIM(my_format))atm1,atm2,atm3,atm4
          END SELECT
       END IF
    END IF

  END SUBROUTINE store_FF_missing_par

END MODULE force_fields_all

