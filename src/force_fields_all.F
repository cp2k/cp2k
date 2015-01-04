!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Splitting and cleaning the original force_field_pack - May 2007
!>      Teodoro Laino - Zurich University
!> \author CJM
! *****************************************************************************
MODULE force_fields_all

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set,&
                                             set_atomic_kind
  USE atoms_input,                     ONLY: read_shell_coord_input
  USE cell_types,                      ONLY: cell_type
  USE cp_linked_list_val,              ONLY: cp_sll_val_next,&
                                             cp_sll_val_type
  USE damping_dipole_types,            ONLY: damping_p_create,&
                                             damping_p_type,&
                                             tang_toennies
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_env_set,&
                                             ewald_environment_type
  USE external_potential_types,        ONLY: fist_potential_type,&
                                             get_potential,&
                                             set_potential
  USE force_field_kind_types,          ONLY: &
       allocate_bend_kind_set, allocate_bond_kind_set, &
       allocate_impr_kind_set, allocate_opbend_kind_set, &
       allocate_torsion_kind_set, allocate_ub_kind_set, bend_kind_type, &
       bond_kind_type, do_ff_amber, do_ff_charmm, do_ff_g87, do_ff_g96, &
       do_ff_undef, impr_kind_type, opbend_kind_type, torsion_kind_type, &
       ub_kind_type
  USE force_field_types,               ONLY: amber_info_type,&
                                             charmm_info_type,&
                                             force_field_type,&
                                             gromos_info_type,&
                                             input_info_type
  USE input_constants,                 ONLY: do_qmmm_none
  USE input_cp2k_binary_restarts,      ONLY: read_binary_cs_coordinates
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_list_get,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE input_val_types,                 ONLY: val_get,&
                                             val_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: sqrthalf
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: &
       bend_type, bond_type, get_molecule_kind, impr_type, &
       molecule_kind_type, opbend_type, set_molecule_kind, shell_type, &
       torsion_type, ub_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE pair_potential,                  ONLY: get_nonbond_storage,&
                                             spline_nonbond_control
  USE pair_potential_coulomb,          ONLY: potential_coulomb
  USE pair_potential_types,            ONLY: &
       ea_type, lj_charmm_type, lj_type, nn_type, nosh_nosh, nosh_sh, &
       pair_potential_lj_create, pair_potential_pp_create, &
       pair_potential_pp_type, pair_potential_single_add, &
       pair_potential_single_clean, pair_potential_single_copy, &
       pair_potential_single_type, quip_type, sh_sh, siepmann_type, &
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
                                             spline_env_release,&
                                             spline_environment_type
  USE string_utilities,                ONLY: compress,&
                                             integer_to_string,&
                                             uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_fields_all'

  PRIVATE
  PUBLIC :: force_field_unique_bond,&
            force_field_unique_bend,&
            force_field_unique_ub,&
            force_field_unique_tors,&
            force_field_unique_impr,&
            force_field_unique_opbend,&
            force_field_pack_bond,&
            force_field_pack_bend,&
            force_field_pack_ub,&
            force_field_pack_tors,&
            force_field_pack_impr,&
            force_field_pack_opbend,&
            force_field_pack_charge,&
            force_field_pack_charges,&
            force_field_pack_radius,&
            force_field_pack_pol, &
            force_field_pack_shell,&
            force_field_pack_nonbond14,&
            force_field_pack_nonbond,&
            force_field_pack_splines,&
            force_field_pack_eicut,&
            force_field_pack_damp

CONTAINS

! *****************************************************************************
!> \brief Determine the number of unique bond kind and allocate bond_kind_set
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param ff_type ...
!> \param error ...
! *****************************************************************************
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
    CALL timeset(routineN,handle2)
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
         IF((ff_type%ff_type==do_ff_g96).OR.(ff_type%ff_type==do_ff_g87)) THEN
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
    CALL timestop(handle2)

  END SUBROUTINE force_field_unique_bond

! *****************************************************************************
!> \brief Determine the number of unique bend kind and allocate bend_kind_set
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param ff_type ...
!> \param error ...
! *****************************************************************************
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

    CALL timeset(routineN,handle2)
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
         IF((ff_type%ff_type==do_ff_g96).OR.(ff_type%ff_type==do_ff_g87)) THEN
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
    CALL timestop(handle2)

  END SUBROUTINE force_field_unique_bend

! *****************************************************************************
!> \brief Determine the number of unique Urey-Bradley kind and allocate ub_kind_set
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param error ...
! *****************************************************************************
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

    CALL timeset(routineN,handle2)
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
    CALL timestop(handle2)

  END SUBROUTINE force_field_unique_ub

! *****************************************************************************
!> \brief Determine the number of unique torsion kind and allocate torsion_kind_set
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param ff_type ...
!> \param error ...
! *****************************************************************************
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
    CALL timeset(routineN,handle2)
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
         IF((ff_type%ff_type==do_ff_g96).OR.(ff_type%ff_type==do_ff_g87)) THEN
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
    CALL timestop(handle2)

  END SUBROUTINE force_field_unique_tors

! *****************************************************************************
!> \brief Determine the number of unique impr kind and allocate impr_kind_set
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param ff_type ...
!> \param error ...
! *****************************************************************************
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

    CALL timeset(routineN,handle2)
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
         IF((ff_type%ff_type==do_ff_g96).OR.(ff_type%ff_type==do_ff_g87)) THEN
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
                     (((name_atm_a)==(name_atm_a2)) .AND. &
                     ((name_atm_b)==(name_atm_b2)) .AND. &
                     ((name_atm_c)==(name_atm_d2)) .AND. &
                     ((name_atm_d)==(name_atm_c2))) ) THEN
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
    CALL timestop(handle2)

  END SUBROUTINE force_field_unique_impr

! *****************************************************************************
!> \brief Determine the number of unique opbend kind and allocate opbend_kind_set
!>        based on the present impropers. With each improper, there also
!>        corresponds a opbend
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param ff_type ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE force_field_unique_opbend (particle_set, &
       molecule_kind_set, molecule_set, ff_type, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_unique_opbend', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length) :: name_atm_a, name_atm_a2, &
      name_atm_b, name_atm_b2, name_atm_c, name_atm_c2, name_atm_d, &
      name_atm_d2
    INTEGER                                  :: atm_a, atm_b, atm_c, atm_d, &
                                                counter, first, handle2, i, &
                                                j, k, last, natom, nopbend, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    INTEGER, POINTER                         :: map_opbend_kind(:)
    LOGICAL                                  :: failure, found
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(opbend_kind_type), DIMENSION(:), &
      POINTER                                :: opbend_kind_set
    TYPE(opbend_type), DIMENSION(:), POINTER :: opbend_list

    failure = .FALSE.

    CALL timeset(routineN,handle2)
    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            nopbend=nopbend,opbend_list=opbend_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       IF(nopbend>0) THEN
          ALLOCATE(map_opbend_kind(nopbend),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          counter=0
         IF((ff_type%ff_type==do_ff_g96).OR.(ff_type%ff_type==do_ff_g87)) THEN
          DO j=1,nopbend
            map_opbend_kind(j)=j
          END DO
          counter=nopbend
         ELSE
          DO j=1,nopbend
             atm_a = opbend_list(j)%a
             atomic_kind => particle_set(atm_a+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_a)
             atm_b = opbend_list(j)%b
             atomic_kind => particle_set(atm_b+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_b)
             atm_c = opbend_list(j)%c
             atomic_kind => particle_set(atm_c+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_c)
             atm_d = opbend_list(j)%d
             atomic_kind => particle_set(atm_d+first-1)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name_atm_d)
             found = .FALSE.
             DO k=1,j-1
                atm_a = opbend_list(k)%a
                atomic_kind => particle_set(atm_a+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_a2)
                atm_b = opbend_list(k)%b
                atomic_kind => particle_set(atm_b+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_b2)
                atm_c = opbend_list(k)%c
                atomic_kind => particle_set(atm_c+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_c2)
                atm_d = opbend_list(k)%d
                atomic_kind => particle_set(atm_d+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name_atm_d2)
                IF( (((name_atm_a)==(name_atm_a2)) .AND. &
                     ((name_atm_b)==(name_atm_b2)) .AND. &
                     ((name_atm_c)==(name_atm_c2)) .AND. &
                     ((name_atm_d)==(name_atm_d2))) .OR. &
                     (((name_atm_a)==(name_atm_a2)) .AND. &
                     ((name_atm_b)==(name_atm_c2)) .AND. &
                     ((name_atm_c)==(name_atm_b2)) .AND. &
                     ((name_atm_d)==(name_atm_d2))) ) THEN
                   found = .TRUE.
                   map_opbend_kind(j) = map_opbend_kind(k)
                   EXIT
                END IF
             END DO
             IF(.NOT.found) THEN
                counter=counter+1
                map_opbend_kind(j) = counter
             END IF
          END DO
         END IF
         NULLIFY(opbend_kind_set)
         CALL allocate_opbend_kind_set(opbend_kind_set,counter,error)
         DO j=1,nopbend
            opbend_list(j)%opbend_kind => opbend_kind_set(map_opbend_kind(j))
         END DO
         CALL set_molecule_kind(molecule_kind=molecule_kind,&
              opbend_kind_set=opbend_kind_set,opbend_list=opbend_list)
         DEALLOCATE(map_opbend_kind,STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO
    CALL timestop(handle2)

  END SUBROUTINE force_field_unique_opbend

! *****************************************************************************
!> \brief Pack in bonds information needed for the force_field
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param fatal ...
!> \param Ainfo ...
!> \param chm_info ...
!> \param inp_info ...
!> \param gro_info ...
!> \param amb_info ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE force_field_pack_bond (particle_set, molecule_kind_set, molecule_set, &
      fatal, Ainfo, chm_info, inp_info, gro_info, amb_info, error)

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
    TYPE(amber_info_type), POINTER           :: amb_info
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
    CALL timeset(routineN,handle2)

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

          ! loop over params from GROMOS
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
             bond_list(j)%bond_kind%id_type =  gro_info%ff_gromos_type
             bond_list(j)%id_type           =  gro_info%ff_gromos_type
             found = .TRUE.
          END IF

          ! loop over params from CHARMM
          IF(ASSOCIATED(chm_info%bond_a)) THEN
             DO k=1,SIZE(chm_info%bond_a)
                IF( (((chm_info%bond_a(k))==(name_atm_a)) .AND. &
                     ((chm_info%bond_b(k))==(name_atm_b))) .OR. &
                     (((chm_info%bond_a(k))==(name_atm_b)) .AND. &
                     ((chm_info%bond_b(k))==(name_atm_a))) ) THEN
                   bond_list(j)%bond_kind%id_type = do_ff_charmm
                   bond_list(j)%bond_kind%k(1) = chm_info%bond_k(k)
                   bond_list(j)%bond_kind%r0   = chm_info%bond_r0(k)
                   CALL issue_duplications(found,routineP,"Bond", name_atm_a, name_atm_b,&
                        error=error)
                   found = .TRUE.
                   EXIT
                END IF
             END DO
          END IF

          ! loop over params from AMBER
          IF(ASSOCIATED(amb_info%bond_a)) THEN
             DO k=1,SIZE(amb_info%bond_a)
                IF( (((amb_info%bond_a(k))==(name_atm_a)) .AND. &
                     ((amb_info%bond_b(k))==(name_atm_b))) .OR. &
                     (((amb_info%bond_a(k))==(name_atm_b)) .AND. &
                     ((amb_info%bond_b(k))==(name_atm_a))) ) THEN
                   bond_list(j)%bond_kind%id_type = do_ff_amber
                   bond_list(j)%bond_kind%k(1) = amb_info%bond_k(k)
                   bond_list(j)%bond_kind%r0   = amb_info%bond_r0(k)
                   CALL issue_duplications(found,routineP,"Bond", name_atm_a, name_atm_b,&
                        error=error)
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
                   bond_list(j)%bond_kind%cs      = inp_info%bond_cs(k)
                   CALL issue_duplications(found,routineP,"Bond", name_atm_a, name_atm_b,&
                        error=error)
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
    CALL timestop(handle2)

  END SUBROUTINE force_field_pack_bond

! *****************************************************************************
!> \brief Pack in bends information needed for the force_field
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param fatal ...
!> \param Ainfo ...
!> \param chm_info ...
!> \param inp_info ...
!> \param gro_info ...
!> \param amb_info ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE force_field_pack_bend (particle_set, molecule_kind_set, molecule_set, &
       fatal, Ainfo, chm_info, inp_info, gro_info, amb_info, error)

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
    TYPE(amber_info_type), POINTER           :: amb_info
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

    CALL timeset(routineN,handle2)
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
             bend_list(j)%bend_kind%id_type = gro_info%ff_gromos_type
             bend_list(j)%id_type           = gro_info%ff_gromos_type
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
                   CALL issue_duplications(found,routineP, "Bend", name_atm_a, name_atm_b,&
                        name_atm_c, error=error)
                   found = .TRUE.
                   EXIT
                END IF
             END DO
          END IF

          ! loop over params from AMBER
          IF(ASSOCIATED(amb_info%bend_a)) THEN
             DO k=1,SIZE(amb_info%bend_a)
                IF( (((amb_info%bend_a(k))==(name_atm_a)) .AND. &
                     ((amb_info%bend_b(k))==(name_atm_b)) .AND. &
                     ((amb_info%bend_c(k))==(name_atm_c))) .OR. &
                     (((amb_info%bend_a(k))==(name_atm_c)) .AND. &
                     ((amb_info%bend_b(k))==(name_atm_b)) .AND. &
                     ((amb_info%bend_c(k))==(name_atm_a))) ) THEN
                   bend_list(j)%bend_kind%id_type = do_ff_amber
                   bend_list(j)%bend_kind%k       = amb_info%bend_k(k)
                   bend_list(j)%bend_kind%theta0  = amb_info%bend_theta0(k)
                   CALL issue_duplications(found,routineP, "Bend", name_atm_a, name_atm_b,&
                        name_atm_c, error=error)
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
                   bend_list(j)%bend_kind%cb      = inp_info%bend_cb(k)
                   bend_list(j)%bend_kind%r012    = inp_info%bend_r012(k)
                   bend_list(j)%bend_kind%r032    = inp_info%bend_r032(k)
                   bend_list(j)%bend_kind%kbs12   = inp_info%bend_kbs12(k)
                   bend_list(j)%bend_kind%kbs32   = inp_info%bend_kbs32(k)
                   bend_list(j)%bend_kind%kss     = inp_info%bend_kss(k)
                   CALL issue_duplications(found,routineP, "Bend", name_atm_a, name_atm_b,&
                        name_atm_c, error=error)
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
    CALL timestop(handle2)

  END SUBROUTINE force_field_pack_bend

! *****************************************************************************
!> \brief Pack in Urey-Bradley information needed for the force_field
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param Ainfo ...
!> \param chm_info ...
!> \param inp_info ...
!> \param iw ...
!> \param error ...
! *****************************************************************************
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

    CALL timeset(routineN,handle2)
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

          ! loop over params from CHARMM
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
                   CALL issue_duplications(found,routineP, "Urey-Bradley", name_atm_a,&
                        name_atm_b, name_atm_c, error=error)
                   found = .TRUE.
                   EXIT
                END IF
             END DO
          END IF

          ! loop over params from AMBER
          ! teo - None that I know...

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
                   CALL issue_duplications(found,routineP, "Urey-Bradley", name_atm_a,&
                        name_atm_b, name_atm_c, error=error)
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
    CALL timestop(handle2)

  END SUBROUTINE force_field_pack_ub

! *****************************************************************************
!> \brief Pack in torsion information needed for the force_field
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param Ainfo ...
!> \param chm_info ...
!> \param inp_info ...
!> \param gro_info ...
!> \param amb_info ...
!> \param iw ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE force_field_pack_tors (particle_set, molecule_kind_set, molecule_set, &
       Ainfo, chm_info, inp_info, gro_info, amb_info, iw, error)

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
    TYPE(amber_info_type), POINTER           :: amb_info
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

    CALL timeset(routineN,handle2)
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
                   CALL reallocate(torsion_list(j)%torsion_kind%phi0,1,1)
                   torsion_list(j)%torsion_kind%nmul   = 1
                   torsion_list(j)%torsion_kind%m(1)   = gro_info%torsion_m(itype)
                   torsion_list(j)%torsion_kind%k(1)   = gro_info%torsion_k(itype)
                   torsion_list(j)%torsion_kind%phi0(1) = gro_info%torsion_phi0(itype)
                ELSE
                   CALL reallocate(torsion_list(j)%torsion_kind%k,1,1)
                   CALL reallocate(torsion_list(j)%torsion_kind%m,1,1)
                   CALL reallocate(torsion_list(j)%torsion_kind%phi0,1,1)
                   torsion_list(j)%torsion_kind%nmul   = 1
                   torsion_list(j)%torsion_kind%m(1)   = gro_info%torsion_m(itype/k)
                   torsion_list(j)%torsion_kind%k(1)   = gro_info%torsion_k(itype/k)
                   torsion_list(j)%torsion_kind%phi0(1) = gro_info%torsion_phi0(itype/k)
                END IF
                torsion_list(j)%torsion_kind%id_type = gro_info%ff_gromos_type
                torsion_list(j)%id_type              = gro_info%ff_gromos_type
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
                      CALL reallocate(torsion_list(j)%torsion_kind%phi0,1,imul)
                      torsion_list(j)%torsion_kind%id_type = do_ff_charmm
                      torsion_list(j)%torsion_kind%k(imul) = chm_info%torsion_k(k)
                      torsion_list(j)%torsion_kind%m(imul) = chm_info%torsion_m(k)
                      torsion_list(j)%torsion_kind%phi0(imul) = chm_info%torsion_phi0(k)
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
                         CALL reallocate(torsion_list(j)%torsion_kind%phi0,1,imul)
                         torsion_list(j)%torsion_kind%id_type = do_ff_charmm
                         torsion_list(j)%torsion_kind%k(imul) = chm_info%torsion_k(k)
                         torsion_list(j)%torsion_kind%m(imul) = chm_info%torsion_m(k)
                         torsion_list(j)%torsion_kind%phi0(imul) = chm_info%torsion_phi0(k)
                         torsion_list(j)%torsion_kind%nmul = imul
                         found = .TRUE.
                      END IF
                   END DO
                END IF
             END IF

             ! loop over params from AMBER
             IF(ASSOCIATED(amb_info%torsion_a)) THEN
                DO k=1,SIZE(amb_info%torsion_a)
                   IF( (((amb_info%torsion_a(k))==(name_atm_a)) .AND. &
                        ((amb_info%torsion_b(k))==(name_atm_b)) .AND. &
                        ((amb_info%torsion_c(k))==(name_atm_c)) .AND. &
                        ((amb_info%torsion_d(k))==(name_atm_d))) .OR. &
                       (((amb_info%torsion_a(k))==(name_atm_d)) .AND. &
                        ((amb_info%torsion_b(k))==(name_atm_c)) .AND. &
                        ((amb_info%torsion_c(k))==(name_atm_b)) .AND. &
                        ((amb_info%torsion_d(k))==(name_atm_a))) ) THEN
                      imul = torsion_list(j)%torsion_kind%nmul + 1
                      CALL reallocate(torsion_list(j)%torsion_kind%k,1,imul)
                      CALL reallocate(torsion_list(j)%torsion_kind%m,1,imul)
                      CALL reallocate(torsion_list(j)%torsion_kind%phi0,1,imul)
                      torsion_list(j)%torsion_kind%id_type = do_ff_amber
                      torsion_list(j)%torsion_kind%k(imul) = amb_info%torsion_k(k)
                      torsion_list(j)%torsion_kind%m(imul) = amb_info%torsion_m(k)
                      torsion_list(j)%torsion_kind%phi0(imul) = amb_info%torsion_phi0(k)
                      torsion_list(j)%torsion_kind%nmul = imul
                      found = .TRUE.
                   END IF
                END DO

                IF(.NOT.found) THEN
                   DO k=1,SIZE(amb_info%torsion_a)
                      IF( (((amb_info%torsion_a(k))==("X")) .AND. &
                           ((amb_info%torsion_b(k))==(name_atm_b)) .AND. &
                           ((amb_info%torsion_c(k))==(name_atm_c)) .AND. &
                           ((amb_info%torsion_d(k))==("X"))) .OR. &
                           (((amb_info%torsion_a(k))==("X")) .AND. &
                           ((amb_info%torsion_b(k))==(name_atm_c)) .AND. &
                           ((amb_info%torsion_c(k))==(name_atm_b)) .AND. &
                           ((amb_info%torsion_d(k))==("X"))) ) THEN
                         imul = torsion_list(j)%torsion_kind%nmul + 1
                         CALL reallocate(torsion_list(j)%torsion_kind%k,1,imul)
                         CALL reallocate(torsion_list(j)%torsion_kind%m,1,imul)
                         CALL reallocate(torsion_list(j)%torsion_kind%phi0,1,imul)
                         torsion_list(j)%torsion_kind%id_type = do_ff_amber
                         torsion_list(j)%torsion_kind%k(imul) = amb_info%torsion_k(k)
                         torsion_list(j)%torsion_kind%m(imul) = amb_info%torsion_m(k)
                         torsion_list(j)%torsion_kind%phi0(imul) = amb_info%torsion_phi0(k)
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
                      CALL reallocate(torsion_list(j)%torsion_kind%phi0,1,imul)
                      torsion_list(j)%torsion_kind%id_type = inp_info%torsion_kind(k)
                      torsion_list(j)%torsion_kind%k(imul) = inp_info%torsion_k(k)
                      torsion_list(j)%torsion_kind%m(imul) = inp_info%torsion_m(k)
                      torsion_list(j)%torsion_kind%phi0(imul) = inp_info%torsion_phi0(k)
                      torsion_list(j)%torsion_kind%nmul    = imul
                      found = .TRUE.
                   END IF
                END DO
             END IF

             IF(.NOT.found) THEN
                CALL store_FF_missing_par(atm1=TRIM(name_atm_a),&
                                          atm2=TRIM(name_atm_b),&
                                          atm3=TRIM(name_atm_c),&
                                          atm4=TRIM(name_atm_d),&
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
    CALL timestop(handle2)

  END SUBROUTINE force_field_pack_tors

! *****************************************************************************
!> \brief Pack in impropers information needed for the force_field
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param Ainfo ...
!> \param chm_info ...
!> \param inp_info ...
!> \param gro_info ...
!> \param amb_info ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE force_field_pack_impr (particle_set, molecule_kind_set, molecule_set, &
       Ainfo, chm_info, inp_info, gro_info, amb_info, error)

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
    TYPE(amber_info_type), POINTER           :: amb_info
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

    CALL timeset(routineN,handle2)
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
            impr_list(j)%impr_kind%id_type = gro_info%ff_gromos_type
            impr_list(j)%id_type           = gro_info%ff_gromos_type
          END IF

          ! loop over params from CHARMM
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
                   CALL issue_duplications(found,routineP, "Impropers", name_atm_a, name_atm_b,&
                        name_atm_c, name_atm_d, error=error)
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
                      CALL issue_duplications(found,routineP, "Impropers", name_atm_a, name_atm_b,&
                           name_atm_c,name_atm_d, error=error)
                      found = .TRUE.
                      EXIT
                   END IF
                END DO
             END IF
          END IF

          ! loop over params from AMBER not needed since impropers in AMBER
          ! are treated like standard torsions

          ! always have the input param last to overwrite all the other ones
          IF(ASSOCIATED(inp_info%impr_a)) THEN
             DO k=1,SIZE(inp_info%impr_a)
                IF( ((inp_info%impr_a(k))==(name_atm_a)) .AND. &
                    ((inp_info%impr_b(k))==(name_atm_b)) .AND. &
                    ((((inp_info%impr_c(k))==(name_atm_c)) .AND. &
                     ((inp_info%impr_d(k))==(name_atm_d))) .OR. &
                    (((inp_info%impr_c(k))==(name_atm_d)) .AND. &
                     ((inp_info%impr_d(k))==(name_atm_c)))) ) THEN
                   impr_list(j)%impr_kind%id_type = inp_info%impr_kind(k)
                   impr_list(j)%impr_kind%k       = inp_info%impr_k(k)
                   IF ( ((inp_info%impr_c(k))==(name_atm_c)) .AND. &
                      ((inp_info%impr_d(k))==(name_atm_d)) ) THEN
                     impr_list(j)%impr_kind%phi0    = inp_info%impr_phi0(k)
                   ELSE
                     impr_list(j)%impr_kind%phi0    = -inp_info%impr_phi0(k)
                     ! alternative solutions:
                     !  - swap impr_list(j)%c with impr_list(j)%d and
                     !    name_atom_c with name_atom_d and
                     !    atm_c with atm_d
                     !  - introduce impr_list(j)%impr_kind%sign. if one, the
                     !    sign of phi is not changed in mol_force.f90. if minus
                     !    one, the sign of phi is changed in mol_force.f90
                     ! similar problems with parameters from charmm pot file
                     ! above
                   END IF
                   CALL issue_duplications(found,routineP, "Impropers", name_atm_a, name_atm_b,&
                        name_atm_c,name_atm_d, error=error)
                   found = .TRUE.
                   EXIT
                END IF
             END DO
          END IF

          IF(.NOT.found) THEN
             CALL store_FF_missing_par(atm1=TRIM(name_atm_a),&
                                       atm2=TRIM(name_atm_b),&
                                       atm3=TRIM(name_atm_c),&
                                       atm4=TRIM(name_atm_d),&
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
    CALL timestop(handle2)

  END SUBROUTINE force_field_pack_impr

! *****************************************************************************
!> \brief Pack in opbend information needed for the force_field.
!>        No loop over params for charmm, amber and gromos since these force
!>        fields have no opbends
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param Ainfo ...
!> \param inp_info ...
!> \param error ...
!> \author Louis Vanduyfhuys
! *****************************************************************************
  SUBROUTINE force_field_pack_opbend (particle_set, molecule_kind_set, molecule_set, &
       Ainfo, inp_info, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: Ainfo
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_opbend', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_b, &
                                                name_atm_c, name_atm_d
    INTEGER                                  :: atm_a, atm_b, atm_c, atm_d, &
                                                first, handle2, i, j, k, &
                                                last, natom, nopbend
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    LOGICAL                                  :: failure, found, only_qm
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(opbend_type), DIMENSION(:), POINTER :: opbend_list

    CALL timeset(routineN,handle2)
    failure = .FALSE.

    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            nopbend=nopbend,opbend_list=opbend_list)
       molecule=>molecule_set(molecule_list(1))

       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       DO j=1,nopbend
          atm_a = opbend_list(j)%a
          atomic_kind => particle_set(atm_a+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_a)
          atm_b = opbend_list(j)%b
          atomic_kind => particle_set(atm_b+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_b)
          atm_c = opbend_list(j)%c
          atomic_kind => particle_set(atm_c+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_c)
          atm_d = opbend_list(j)%d
          atomic_kind => particle_set(atm_d+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name_atm_d)
          found = .FALSE.
          only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b, id3=name_atm_c, id4=name_atm_d)
          CALL uppercase(name_atm_a)
          CALL uppercase(name_atm_b)
          CALL uppercase(name_atm_c)
          CALL uppercase(name_atm_d)

          ! always have the input param last to overwrite all the other ones
          IF(ASSOCIATED(inp_info%opbend_a)) THEN
             DO k=1,SIZE(inp_info%opbend_a)
                IF( ((inp_info%opbend_a(k))==(name_atm_a)) .AND. &
                    ((inp_info%opbend_d(k))==(name_atm_d)) .AND. &
                    ((((inp_info%opbend_c(k))==(name_atm_c)) .AND. &
                     ((inp_info%opbend_b(k))==(name_atm_b))) .OR. &
                    (((inp_info%opbend_c(k))==(name_atm_b)) .AND. &
                     ((inp_info%opbend_b(k))==(name_atm_c)))) ) THEN
                   opbend_list(j)%opbend_kind%id_type = inp_info%opbend_kind(k)
                   opbend_list(j)%opbend_kind%k       = inp_info%opbend_k(k)
                   IF ( ((inp_info%opbend_c(k))==(name_atm_c)) .AND. &
                      ((inp_info%opbend_b(k))==(name_atm_b)) ) THEN
                     opbend_list(j)%opbend_kind%phi0    = inp_info%opbend_phi0(k)
                   ELSE
                     opbend_list(j)%opbend_kind%phi0    = -inp_info%opbend_phi0(k)
                     ! alternative solutions:
                     !  - swap opbend_list(j)%c with opbend_list(j)%b and
                     !    name_atom_c with name_atom_b and
                     !    atm_c with atm_b
                     !  - introduce opbend_list(j)%opbend_kind%sign. if one, the
                     !    sign of phi is not changed in mol_force.f90. if minus
                     !    one, the sign of phi is changed in mol_force.f90

                   END IF
                   CALL issue_duplications(found,routineP, "Out of plane bend", name_atm_a, name_atm_b,&
                        name_atm_c,name_atm_d, error=error)
                   found = .TRUE.
                   EXIT
                END IF
             END DO
          END IF

          IF(.NOT.found) THEN
             CALL store_FF_missing_par(atm1=TRIM(name_atm_a),&
                                       atm2=TRIM(name_atm_b),&
                                       atm3=TRIM(name_atm_c),&
                                       atm4=TRIM(name_atm_d),&
                                       type_name="Out of plane bend",&
                                       array=Ainfo,&
                                       error=error)
             opbend_list(j)%opbend_kind%k       = 0.0_dp
             opbend_list(j)%opbend_kind%phi0    = 0.0_dp
             opbend_list(j)%opbend_kind%id_type = do_ff_undef
             opbend_list(j)%id_type             = do_ff_undef
          END IF
          !
          ! QM/MM modifications
          !
          IF (only_qm) THEN
             opbend_list(j)%opbend_kind%id_type = do_ff_undef
             opbend_list(j)%id_type             = do_ff_undef
          END IF

       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,opbend_list=opbend_list)
    END DO
    CALL timestop(handle2)

  END SUBROUTINE force_field_pack_opbend

! *****************************************************************************
!> \brief Set up array of full charges
!> \param charges ...
!> \param charges_section ...
!> \param particle_set ...
!> \param my_qmmm ...
!> \param qmmm_env ...
!> \param inp_info ...
!> \param iw4 ...
!> \param error ...
!> \date 12.2010
!> \author Teodoro Laino (teodoro.laino@gmail.com)
! *****************************************************************************
  SUBROUTINE force_field_pack_charges(charges, charges_section, particle_set, &
             my_qmmm, qmmm_env, inp_info, iw4, error)

    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    TYPE(section_vals_type), POINTER         :: charges_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    LOGICAL                                  :: my_qmmm
    TYPE(qmmm_env_mm_type), POINTER          :: qmmm_env
    TYPE(input_info_type), POINTER           :: inp_info
    INTEGER                                  :: iw4
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_charges', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atmname
    INTEGER                                  :: handle, iatom, ilink, j, &
                                                nval, stat
    LOGICAL                                  :: failure, found_p, &
                                                is_link_atom, is_ok, &
                                                only_manybody, only_qm
    REAL(KIND=dp)                            :: charge, charge_tot, rval, &
                                                scale_factor
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_sll_val_type), POINTER           :: list
    TYPE(fist_potential_type), POINTER       :: fist_potential
    TYPE(val_type), POINTER                  :: val

    failure = .FALSE.

    CALL timeset(routineN,handle)
    charge_tot = 0.0_dp
    NULLIFY(list)

    ! Not implemented for core-shell
    IF(ASSOCIATED(inp_info%shell_list)) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
                 message="Array of charges not implemented for CORE-SHELL model!!",&
                 error=error, error_level=cp_failure_level)
    END IF

    ! Allocate array to particle_set size
    CPPostcondition(.NOT.(ASSOCIATED(charges)),cp_failure_level,routineP,error,failure)
    ALLOCATE(charges(SIZE(particle_set)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Fill with input values
    CALL section_vals_val_get(charges_section,"_DEFAULT_KEYWORD_",n_rep_val=nval,error=error)
    CPPostcondition(nval==SIZE(charges),cp_failure_level,routineP,error,failure)
    CALL section_vals_list_get(charges_section,"_DEFAULT_KEYWORD_",list=list,error=error)
    DO iatom=1,nval
       ! we use only the first default_string_length characters of each line
       is_ok=cp_sll_val_next(list,val,error=error)
       CALL val_get(val,r_val=rval,error=error)
       ! assign values
       charges(iatom)=rval

       ! Perform a post-processing
       atomic_kind => particle_set(iatom)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            fist_potential=fist_potential,&
                            name=atmname)
       CALL get_potential(potential=fist_potential, qeff=charge)

       only_qm  = qmmm_ff_precond_only_qm(id1=atmname,is_link=is_link_atom)
       CALL uppercase(atmname)
       CALL cp_assert(charge==-HUGE(0.0_dp),cp_warning_level,cp_assertion_failed,routineP,&
            "The charge for atom index ("//cp_to_string(iatom)//") and atom name ("//&
            TRIM(atmname)//") was already defined. The charge associated to this kind"//&
            " will be set to an uninitialized value and only the atom specific charge will be used! "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
       charge=-HUGE(0.0_dp)

       ! Check if the potential really requires the charge definition..
       IF (ASSOCIATED(inp_info%nonbonded)) THEN
          IF (ASSOCIATED(inp_info%nonbonded%pot)) THEN
             ! Let's find the nonbonded potential where this atom is involved
             only_manybody = .TRUE.
             found_p       = .FALSE.
             DO j = 1, SIZE(inp_info%nonbonded%pot)
                IF ( atmname==inp_info%nonbonded%pot(j)%pot %at1.OR.&
                     atmname==inp_info%nonbonded%pot(j)%pot %at2) THEN
                   SELECT CASE(inp_info%nonbonded%pot(j)%pot%type(1))
                   CASE (ea_type,tersoff_type,siepmann_type)
                      ! Charge is zero for EAM, TERSOFF and SIEPMANN  type potential
                      ! Do nothing..
                   CASE DEFAULT
                      only_manybody = .FALSE.
                      EXIT
                   END SELECT
                   found_p = .TRUE.
                END IF
             END DO
             IF (only_manybody.AND.found_p) THEN
                charges(iatom)= 0.0_dp
             END IF
          END IF
       END IF
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
                   IF (iatom == qmmm_env%mm_link_atoms(ilink)) EXIT
                END DO
                CPPostcondition(ilink <= SIZE(qmmm_env%mm_link_atoms),cp_failure_level,routineP,error,failure)
                scale_factor = qmmm_env%fist_scale_charge_link(ilink)
             END IF
             charges(iatom) = charges(iatom) * scale_factor
          END IF
       END IF
    END DO
    ! Sum up total charges for IO
    charge_tot = SUM(charges)
    ! Print Total Charge of the system
    IF (iw4>0) THEN
       WRITE(iw4,FMT='(T2,"CHARGE_INFO| Total Charge of the Classical System: ",T69,F12.6)')charge_tot
    END IF
    CALL timestop(handle)
  END SUBROUTINE force_field_pack_charges

! *****************************************************************************
!> \brief Set up atomic_kind_set()%fist_potentail%[qeff]
!>      and shell potential parameters
!> \param atomic_kind_set ...
!> \param qmmm_env ...
!> \param fatal ...
!> \param iw ...
!> \param iw4 ...
!> \param Ainfo ...
!> \param my_qmmm ...
!> \param inp_info ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE force_field_pack_charge(atomic_kind_set, qmmm_env, fatal, iw, iw4,&
       Ainfo, my_qmmm, inp_info, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qmmm_env_mm_type), POINTER          :: qmmm_env
    LOGICAL                                  :: fatal
    INTEGER                                  :: iw, iw4
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: Ainfo
    LOGICAL                                  :: my_qmmm
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_charge', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atmname
    INTEGER                                  :: handle, i, ilink, j
    INTEGER, DIMENSION(:), POINTER           :: my_atom_list
    LOGICAL                                  :: failure, found, found_p, &
                                                is_link_atom, is_shell, &
                                                only_manybody, only_qm
    REAL(KIND=dp)                            :: charge, charge_tot, &
                                                cs_charge, scale_factor
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(fist_potential_type), POINTER       :: fist_potential

    failure = .FALSE.

    CALL timeset(routineN,handle)
    charge_tot = 0.0_dp
    DO i=1,SIZE(atomic_kind_set)
       atomic_kind => atomic_kind_set(i)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            fist_potential=fist_potential,&
                            atom_list=my_atom_list,&
                            name=atmname)
       CALL get_potential(potential=fist_potential, qeff=charge)

       is_shell = .FALSE.
       found    = .FALSE.
       only_qm  = qmmm_ff_precond_only_qm(id1=atmname,is_link=is_link_atom)
       CALL uppercase(atmname)
       IF(charge/=-HUGE(0.0_dp)) found = .TRUE.

       ! Always have the input param last to overwrite all the other ones
       IF(ASSOCIATED(inp_info%charge_atm)) THEN
          DO j=1,SIZE(inp_info%charge_atm)
             IF (iw>0) WRITE(iw,*)"Charge Checking ::",TRIM(inp_info%charge_atm(j)),atmname
             IF((inp_info%charge_atm(j))==atmname) THEN
                charge = inp_info%charge(j)
                CALL issue_duplications(found,routineP, "Charge", atmname, error=error)
                found = .TRUE.
             END IF
          END DO
       END IF
       ! Check if the ATOM type has a core-shell associated.. In this case
       ! print a warning: the CHARGE will not be used if defined.. or we can
       ! even skip its definition..
       IF(ASSOCIATED(inp_info%shell_list)) THEN
          DO j=1,SIZE(inp_info%shell_list)
             IF((inp_info%shell_list(j)%atm_name)==atmname) THEN
                is_shell  = .TRUE.
                cs_charge = inp_info%shell_list(j)%shell%charge_core+&
                            inp_info%shell_list(j)%shell%charge_shell
                charge = 0.0_dp
                IF (found) THEN
                   CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routineP,&
                        "CORE-SHELL model defined for KIND ("//TRIM(atmname)//")"//&
                        " ignoring charge definition! "//&
CPSourceFileRef,&
                        only_ionode=.TRUE.)
                ELSE
                   found = .TRUE.
                END IF
             END IF
          END DO
       END IF
       ! Check if the potential really requires the charge definition..
       IF (ASSOCIATED(inp_info%nonbonded)) THEN
          IF (ASSOCIATED(inp_info%nonbonded%pot)) THEN
             ! Let's find the nonbonded potential where this atom is involved
             only_manybody = .TRUE.
             found_p       = .FALSE.
             DO j = 1, SIZE(inp_info%nonbonded%pot)
                IF ( atmname==inp_info%nonbonded%pot(j)%pot %at1.OR.&
                     atmname==inp_info%nonbonded%pot(j)%pot %at2) THEN
                   SELECT CASE(inp_info%nonbonded%pot(j)%pot%type(1))
                   CASE (ea_type,tersoff_type,siepmann_type,quip_type)
                      ! Charge is zero for EAM, TERSOFF and SIEPMANN type potential
                      ! Do nothing..
                   CASE DEFAULT
                      only_manybody = .FALSE.
                      EXIT
                   END SELECT
                   found_p = .TRUE.
                END IF
             END DO
             IF (only_manybody.AND.found_p) THEN
                charge = 0.0_dp
                found     = .TRUE.
             END IF
          END IF
       END IF
       IF (.NOT.found) THEN
         ! Set the charge to zero anyway in case the user decides to ignore
         ! missing critical parameters.
         charge = 0.0_dp
         CALL store_FF_missing_par(atm1=TRIM(atmname),&
                                   fatal=fatal,&
                                   type_name="Charge",&
                                   array=Ainfo,&
                                   error=error)
       END IF
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
             charge = charge * scale_factor
          END IF
       END IF

       CALL set_potential(potential=fist_potential, qeff=charge)
       ! Sum up total charges for IO
       IF (found) THEN
          IF (is_shell) THEN
             charge_tot = charge_tot + atomic_kind%natom*cs_charge
          ELSE
             charge_tot = charge_tot + atomic_kind%natom*charge
          END IF
       END IF
    END DO
    ! Print Total Charge of the system
    IF (iw4>0) THEN
       WRITE(iw4,FMT='(T2,"CHARGE_INFO| Total Charge of the Classical System: ",T69,F12.6)')charge_tot
    END IF
    CALL timestop(handle)
  END SUBROUTINE force_field_pack_charge

! *****************************************************************************
!> \brief Set up the radius of the electrostatic multipole in Fist
!> \param atomic_kind_set ...
!> \param iw ...
!> \param subsys_section ...
!> \param error ...
!> \author Toon.Verstraelen@gmail.com
! *****************************************************************************
  SUBROUTINE force_field_pack_radius(atomic_kind_set, iw, subsys_section, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, INTENT(IN)                      :: iw
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_radius', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: inp_kind_name, kind_name
    INTEGER                                  :: handle, i, i_rep, n_rep
    LOGICAL                                  :: found
    REAL(KIND=dp)                            :: mm_radius
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(fist_potential_type), POINTER       :: fist_potential
    TYPE(section_vals_type), POINTER         :: kind_section

    CALL timeset(routineN,handle)

    kind_section => section_vals_get_subs_vals(subsys_section,"KIND",error=error)
    CALL section_vals_get(kind_section,n_repetition=n_rep,error=error)

    DO i = 1, SIZE(atomic_kind_set)
       atomic_kind => atomic_kind_set(i)
       CALL get_atomic_kind(atomic_kind=atomic_kind, &
            fist_potential=fist_potential, name=kind_name)
       CALL uppercase(kind_name)
       found = .FALSE.

       ! Try to find a matching KIND section in the SUBSYS section and read the
       ! MM_RADIUS field if it is present. In case the kind section is never
       ! encountered, the mm_radius remains zero.
       mm_radius = 0.0_dp
       DO i_rep = 1, n_rep
          CALL section_vals_val_get(kind_section, "_SECTION_PARAMETERS_",&
            c_val=inp_kind_name, i_rep_section=i_rep, error=error)
          CALL uppercase(inp_kind_name)
          IF (iw>0) WRITE(iw,*) "Matching kinds for MM_RADIUS :: '", &
            TRIM(kind_name), "' with '", TRIM(inp_kind_name), "'"
          IF (TRIM(kind_name)==TRIM(inp_kind_name)) THEN
             CALL section_vals_val_get(kind_section, i_rep_section=i_rep,&
               keyword_name="MM_RADIUS", r_val=mm_radius, error=error)
             CALL issue_duplications(found, routineP, "MM_RADIUS", kind_name, error=error)
             found = .TRUE.
          END IF
       END DO

       CALL set_potential(potential=fist_potential, mm_radius=mm_radius)
    END DO
    CALL timestop(handle)
  END SUBROUTINE force_field_pack_radius

! *****************************************************************************
!> \brief Set up the polarizable FF parameters
!> \param atomic_kind_set ...
!> \param iw ...
!> \param inp_info ...
!> \param error ...
!> \author Toon.Verstraelen@gmail.com
! *****************************************************************************
  SUBROUTINE force_field_pack_pol(atomic_kind_set, iw, inp_info, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, INTENT(IN)                      :: iw
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_pol', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: kind_name
    INTEGER                                  :: handle, i, j
    LOGICAL                                  :: found
    REAL(KIND=dp)                            :: apol, cpol
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(fist_potential_type), POINTER       :: fist_potential

    CALL timeset(routineN,handle)

    DO i=1,SIZE(atomic_kind_set)
       atomic_kind => atomic_kind_set(i)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            fist_potential=fist_potential,&
                            name=kind_name)
       CALL get_potential(potential=fist_potential, apol=apol, cpol=cpol)
       CALL uppercase(kind_name)
       found = .FALSE.

       ! Always have the input param last to overwrite all the other ones
       IF(ASSOCIATED(inp_info%apol_atm)) THEN
          DO j=1,SIZE(inp_info%apol_atm)
             IF (iw>0) WRITE(iw,*) "Matching kinds for APOL :: '", &
                TRIM(kind_name), "' with '", TRIM(inp_info%apol_atm(j)), "'"
             IF((inp_info%apol_atm(j))==kind_name) THEN
                apol = inp_info%apol(j)
                CALL issue_duplications(found,routineP, "APOL", kind_name, error=error)
                found = .TRUE.
             END IF
          END DO
       END IF

       IF(ASSOCIATED(inp_info%cpol_atm)) THEN
          DO j=1,SIZE(inp_info%cpol_atm)
             IF (iw>0) WRITE(iw,*) "Matching kinds for CPOL :: '", &
                TRIM(kind_name), "' with '", TRIM(inp_info%cpol_atm(j)), "'"
             IF((inp_info%cpol_atm(j))==kind_name) THEN
                cpol = inp_info%cpol(j)
                CALL issue_duplications(found,routineP, "CPOL", kind_name, error=error)
                found = .TRUE.
             END IF
          END DO
       END IF

       CALL set_potential(potential=fist_potential, apol=apol, cpol=cpol)
    END DO
    CALL timestop(handle)
  END SUBROUTINE force_field_pack_pol

! *****************************************************************************
!> \brief Set up damping parameters
!> \param atomic_kind_set ...
!> \param iw ...
!> \param inp_info ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE force_field_pack_damp (atomic_kind_set,&
                                    iw,inp_info,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER                                  :: iw
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_damp', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: atm_name1, atm_name2, &
                                                my_atm_name1, my_atm_name2
    INTEGER                                  :: handle2, i, j, k, nkinds
    LOGICAL                                  :: found
    TYPE(atomic_kind_type), POINTER          :: atomic_kind, atomic_kind2
    TYPE(damping_p_type), POINTER            :: damping

     CALL timeset(routineN,handle2)
     NULLIFY(damping)
     nkinds=SIZE(atomic_kind_set)

#ifdef DEBUG_DAMPING
     WRITE(*,*) "DAMPING COUCOU"
!     CALL timestop(handle2)
!     RETURN
#endif
     DO j=1,SIZE(atomic_kind_set)
       atomic_kind=>atomic_kind_set(j)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            name=atm_name1)
       CALL UPPERCASE(atm_name1)
#ifdef DEBUG_DAMPING
       WRITE(*,*) "DAMPING in",j,TRIM(atm_name1)
#endif
       IF (ASSOCIATED(inp_info%damping_list)) THEN
         DO i=1,SIZE(inp_info%damping_list)
           my_atm_name1=inp_info%damping_list(i)%atm_name1
           my_atm_name2=inp_info%damping_list(i)%atm_name2

           IF (iw>0) WRITE(iw,*) "Damping Checking ::",TRIM(my_atm_name1),&
                                                       TRIM(atm_name1)
#ifdef DEBUG_DAMPING
           WRITE(*,*) "DAMPING DEBUG, checking ::",TRIM(my_atm_name1),&
                                                   TRIM(atm_name1)
#endif
           IF (my_atm_name1==atm_name1) THEN
              IF (.NOT.ASSOCIATED(damping)) THEN
                 CALL damping_p_create(damping,nkinds,error)
              END IF

              found=.FALSE.
              DO k=1,SIZE(atomic_kind_set)
                 atomic_kind2=>atomic_kind_set(k)
                 CALL get_atomic_kind(atomic_kind=atomic_kind2,&
                                      name=atm_name2)
                 CALL UPPERCASE(atm_name2)

                 IF (my_atm_name2==atm_name2) THEN
#ifdef DEBUG_DAMPING
                   WRITE(*,*) "DAMPING FOUND",TRIM(my_atm_name1),&
                                              TRIM(atm_name1),&
                                              TRIM(my_atm_name2),&
                                              TRIM(atm_name2)
#endif
                    IF(damping%damp(k)%bij/=HUGE(0.0_dp)) found=.TRUE.
                    CALL issue_duplications(found,routineP, "Damping",&
                       atm_name1,error=error)
                    found=.TRUE.

                    SELECT CASE(TRIM(inp_info%damping_list(i)%dtype))
                    CASE('TANG-TOENNIES')
                        damping%damp(k)%itype=tang_toennies
                    CASE DEFAULT
                        CALL stop_program(routineN,moduleN,__LINE__,&
                                          "Unknown damping type.")
                    END SELECT
                    damping%damp(k)%order=inp_info%damping_list(i)%order
                    damping%damp(k)%bij=inp_info%damping_list(i)%bij
                    damping%damp(k)%cij=inp_info%damping_list(i)%cij
                 ENDIF

              END DO
              CALL cp_assert(found,cp_warning_level,cp_assertion_failed,&
                    routineP,&
                    "Atom "//TRIM(my_atm_name2)//&
                    " in damping parameters for atom "//TRIM(my_atm_name1)//&
                    " not found! "//&
CPSourceFileRef,&
                   only_ionode=.TRUE.)
           ENDIF
         END DO

       ENDIF

       CALL set_atomic_kind(atomic_kind=atomic_kind,damping=damping)
       NULLIFY(damping)
     END DO

     CALL timestop(handle2)

  END SUBROUTINE force_field_pack_damp

! *****************************************************************************
!> \brief Set up shell potential parameters
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param root_section ...
!> \param subsys_section ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param cell ...
!> \param iw ...
!> \param inp_info ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE force_field_pack_shell(particle_set, atomic_kind_set,&
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
    LOGICAL :: core_coord_read, failure, found_shell, is_a_shell, &
      is_link_atom, null_massfrac, only_qm, save_mem, shell_adiabatic, &
      shell_coord_read
    REAL(KIND=dp)                            :: atmmass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(section_vals_type), POINTER         :: global_section
    TYPE(shell_kind_type), POINTER           :: shell
    TYPE(shell_type), DIMENSION(:), POINTER  :: shell_list

    failure = .FALSE.

    CALL timeset(routineN,handle2)
    nshell_tot = 0
    n = 0
    first_shell = 1
    null_massfrac =.FALSE.
    core_coord_read = .FALSE.
    shell_coord_read = .FALSE.

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
             IF (iw > 0) WRITE(iw,*)"Shell Checking ::",TRIM(inp_info%shell_list(j)%atm_name),atmname
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
                IF (shell%massfrac < EPSILON(1.0_dp)) null_massfrac = .TRUE.
                shell%k2_spring = inp_info%shell_list(j)%shell%k2_spring
                shell%k4_spring = inp_info%shell_list(j)%shell%k4_spring
                shell%max_dist=inp_info%shell_list(j)%shell%max_dist
                shell%shell_cutoff=inp_info%shell_list(j)%shell%shell_cutoff
                shell%mass_shell=shell%massfrac*atmmass
                shell%mass_core=atmmass-shell%mass_shell
                CALL issue_duplications(found_shell,routineP, "Shell", atmname, error=error)
                found_shell = .TRUE.
                CALL set_atomic_kind(atomic_kind=atomic_kind,&
                     shell=shell, shell_active=.TRUE.,error=error)
                CALL shell_release(shell, error)
             END IF
          END DO ! j shell kind
       END IF ! associated shell_list
    END DO ! i atomic kind

    IF (iw > 0) WRITE(iw,*) "Total number of particles with a shell :: ", nshell_tot
    ! If shell-model is present: Create particle_set of shells (coord. vel. force)
    NULLIFY (shell_particle_set)
    NULLIFY (core_particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,shell_adiabatic=shell_adiabatic)
    IF (nshell_tot > 0) THEN
       IF (shell_adiabatic .AND. null_massfrac) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Shell-model adiabatic: at least one shell_kind has mass zero")
       END IF
       CALL allocate_particle_set(shell_particle_set,nshell_tot,error)
       CALL allocate_particle_set(core_particle_set,nshell_tot,error)
       counter = 0
       ! Initialise the shell (and core) coordinates with the particle (atomic) coordinates,
       ! count the shell and set pointers
       DO i=1,SIZE(particle_set)
          NULLIFY (atomic_kind)
          NULLIFY (shell)
          atomic_kind => particle_set(i)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,shell_active=is_a_shell)
          IF (is_a_shell) THEN
             counter = counter + 1
             particle_set(i)%shell_index = counter
             shell_particle_set(counter)%shell_index = counter
             shell_particle_set(counter)%atomic_kind => particle_set(i)%atomic_kind
             shell_particle_set(counter)%r(1:3) = particle_set(i)%r(1:3)
             shell_particle_set(counter)%atom_index = i
             core_particle_set(counter)%shell_index = counter
             core_particle_set(counter)%atomic_kind => particle_set(i)%atomic_kind
             core_particle_set(counter)%r(1:3) = particle_set(i)%r(1:3)
             core_particle_set(counter)%atom_index = i
          ELSE
             particle_set(i)%shell_index = 0
          END IF
       END DO
       CPPostcondition(counter==nshell_tot,cp_failure_level,routineP,error,failure)
    END IF

    ! Read the shell (and core) coordinates from the restart file, if available
    CALL read_binary_cs_coordinates("SHELL",shell_particle_set,root_section,&
                                    subsys_section,shell_coord_read,error)
    CALL read_binary_cs_coordinates("CORE",core_particle_set,root_section,&
                                    subsys_section,core_coord_read,error)

    IF (nshell_tot > 0) THEN
       ! Read the shell (and core) coordinates from the input, if no coordinates were found
       ! in the restart file
       IF (shell_adiabatic) THEN
          IF (.NOT.(core_coord_read.AND.shell_coord_read)) THEN
             CALL read_shell_coord_input(particle_set,shell_particle_set,cell,&
                                         subsys_section,core_particle_set,&
                                         save_mem=save_mem,error=error)
          END IF
       ELSE
          IF (.NOT.shell_coord_read) THEN
             CALL read_shell_coord_input(particle_set,shell_particle_set,cell,&
                                         subsys_section,save_mem=save_mem,error=error)
          END IF
       END IF
       ! Determine the number of shells per molecule kind
       n = 0
       DO i=1,SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(i)
          CALL get_molecule_kind(molecule_kind=molecule_kind,molecule_list=molecule_list,&
                                 natom=natom, nmolecule=nmol)
          molecule=>molecule_set(molecule_list(1))
          CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
          ALLOCATE (shell_list_tmp(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          counter = 0
          DO j = first, last
             atomic_kind => particle_set(j)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,shell_active=is_a_shell)
             IF (is_a_shell) THEN
                counter = counter+1
                shell_list_tmp(counter) = j - first + 1
                first_shell = MIN(first_shell, MAX(1,particle_set(j)%shell_index))
             END IF
          END DO  ! j atom in molecule_kind i, molecule 1 of the molecule_list
          IF (counter /= 0) THEN
             ! Setup of fist_shell and last_shell for all molecules..
             DO j=1,SIZE(molecule_list)
                last_shell = first_shell + counter -1
                molecule=>molecule_set(molecule_list(j))
                molecule%first_shell = first_shell
                molecule%last_shell  = last_shell
                first_shell = last_shell + 1
             END DO
             ! Setup of shell_list
             CALL get_molecule_kind(molecule_kind=molecule_kind,shell_list=shell_list)
             IF (ASSOCIATED(shell_list)) THEN
                DEALLOCATE(shell_list,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             ALLOCATE (shell_list(counter), STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO j=1,counter
                shell_list(j)%a = shell_list_tmp(j)
                atomic_kind => particle_set(shell_list_tmp(j)+first-1)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,name=atmname, shell=shell)
                CALL uppercase(atmname)
                shell_list(j)%name = atmname
                shell_list(j)%shell_kind => shell
             END DO
             CALL set_molecule_kind(molecule_kind=molecule_kind, nshell=counter,shell_list=shell_list)
          END IF
          DEALLOCATE (shell_list_tmp,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          n = n + nmol*counter
       END DO ! i molecule kind
    END IF
    CPPostcondition(first_shell-1==nshell_tot,cp_failure_level,routineP,error,failure)
    CPPostcondition(n==nshell_tot,cp_failure_level,routineP,error,failure)
    CALL timestop(handle2)

  END SUBROUTINE force_field_pack_shell

! *****************************************************************************
!> \brief Assign input and potential info to potparm_nonbond14
!> \param atomic_kind_set ...
!> \param ff_type ...
!> \param qmmm_env ...
!> \param iw ...
!> \param Ainfo ...
!> \param chm_info ...
!> \param inp_info ...
!> \param gro_info ...
!> \param amb_info ...
!> \param potparm_nonbond14 ...
!> \param ewald_env ...
!> \param iw2 ...
!> \param iw3 ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE force_field_pack_nonbond14(atomic_kind_set, ff_type, qmmm_env, iw, &
       Ainfo, chm_info, inp_info, gro_info, amb_info, potparm_nonbond14, ewald_env, &
       iw2, iw3, error)

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
    TYPE(amber_info_type), POINTER           :: amb_info
    TYPE(pair_potential_pp_type), POINTER    :: potparm_nonbond14
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    INTEGER                                  :: iw2, iw3
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_nonbond14', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_a_local, &
                                                name_atm_b, name_atm_b_local
    INTEGER                                  :: handle2, i, ii, j, jj, k, &
                                                match_names
    LOGICAL                                  :: failure, found, found_a, &
                                                found_b, only_qm, use_qmmm_ff
    REAL(KIND=dp)                            :: epsilon0, epsilon_a, &
                                                epsilon_b, ewald_rcut, rmin, &
                                                rmin2_a, rmin2_b
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot

    failure = .FALSE.
    use_qmmm_ff = qmmm_env%use_qmmm_ff
    NULLIFY(pot)
    CALL ewald_env_get(ewald_env, rcut=ewald_rcut, error=error)
    CALL timeset(routineN,handle2)
    CALL pair_potential_pp_create (potparm_nonbond14, SIZE(atomic_kind_set), error)
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
          pot => potparm_nonbond14%pot(i,j)%pot

          ! loop over params from GROMOS
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
                CALL issue_duplications(found,routineP, "Lennard-Jones", name_atm_a, name_atm_b, error=error)
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
             CALL issue_duplications(found,routineP, "Lennard-Jones", name_atm_a, name_atm_b,&
                  error=error)
             found = .TRUE.
          END IF

          ! loop over params from AMBER
          IF(ASSOCIATED(amb_info%nonbond_a)) THEN
             ii = 0
             jj = 0
             IF(.NOT.found_a) THEN
                DO k=1,SIZE(amb_info%nonbond_a)
                   IF((name_atm_a)==(amb_info%nonbond_a(k))) THEN
                      ii = k
                      rmin2_a   = amb_info%nonbond_rmin2(k)
                      epsilon_a = amb_info%nonbond_eps(k)
                   END IF
                END DO
             END IF
             IF(.NOT.found_b) THEN
                DO k=1,SIZE(amb_info%nonbond_a)
                   IF((name_atm_b)==(amb_info%nonbond_a(k))) THEN
                      jj = k
                      rmin2_b   = amb_info%nonbond_rmin2(k)
                      epsilon_b = amb_info%nonbond_eps(k)
                   END IF
                END DO
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
                CALL issue_duplications(found,routineP, "Lennard-Jones", name_atm_a,&
                     name_atm_b, error=error)
                found = .TRUE.
             END IF
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
                      potparm_nonbond14%pot(i,j)%pot => pot
                      potparm_nonbond14%pot(j,i)%pot => pot
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
                            potparm_nonbond14%pot(i,j)%pot => pot
                            potparm_nonbond14%pot(j,i)%pot => pot
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
             CALL pair_potential_single_clean(pot, error=error)
             pot%type       = nn_type
             pot%at1        = name_atm_a
             pot%at2        = name_atm_b
          END IF
          ! If defined global RCUT let's use it
          IF(ff_type%rcut_nb>0.0_dp) THEN
             pot%rcutsq = ff_type%rcut_nb*ff_type%rcut_nb
          END IF
          ! Cutoff is defined always as the maximum between the FF and Ewald
          pot%rcutsq = MAX(pot%rcutsq, ewald_rcut*ewald_rcut)
          IF (only_qm) THEN
             CALL pair_potential_single_clean(pot, error=error)
          END IF
       END DO  ! atom kind j
    END DO  ! atom kind i
    CALL timestop(handle2)

  END SUBROUTINE force_field_pack_nonbond14

! *****************************************************************************
!> \brief Assign input and potential info to potparm_nonbond
!> \param atomic_kind_set ...
!> \param ff_type ...
!> \param qmmm_env ...
!> \param fatal ...
!> \param iw ...
!> \param Ainfo ...
!> \param chm_info ...
!> \param inp_info ...
!> \param gro_info ...
!> \param amb_info ...
!> \param potparm_nonbond ...
!> \param ewald_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE force_field_pack_nonbond(atomic_kind_set, ff_type, qmmm_env, fatal, &
       iw, Ainfo, chm_info, inp_info, gro_info, amb_info, potparm_nonbond, &
       ewald_env, error)

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
    TYPE(amber_info_type), POINTER           :: amb_info
    TYPE(pair_potential_pp_type), POINTER    :: potparm_nonbond
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_nonbond', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_atm_a, name_atm_a_local, &
                                                name_atm_b, name_atm_b_local
    INTEGER                                  :: handle2, i, ii, j, jj, k, &
                                                match_names
    LOGICAL                                  :: failure, found, is_a_shell, &
                                                is_b_shell, only_qm, &
                                                use_qmmm_ff
    REAL(KIND=dp)                            :: epsilon0, ewald_rcut, rmin
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot

    CALL timeset(routineN,handle2)
    failure = .FALSE.
    use_qmmm_ff = qmmm_env%use_qmmm_ff
    NULLIFY(pot)
    CALL ewald_env_get(ewald_env, rcut=ewald_rcut, error=error)
    CALL pair_potential_pp_create ( potparm_nonbond, SIZE(atomic_kind_set),error )
    DO i=1,SIZE(atomic_kind_set)
       atomic_kind => atomic_kind_set(i)
       CALL get_atomic_kind(atomic_kind=atomic_kind, name=name_atm_a_local, &
            shell_active=is_a_shell)
       DO j=i,SIZE(atomic_kind_set)
          atomic_kind => atomic_kind_set(j)
          CALL get_atomic_kind(atomic_kind=atomic_kind, name=name_atm_b_local, &
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
                CALL issue_duplications(found,routineP, "Lennard-Jones", name_atm_a, name_atm_b, error=error)
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
                CALL issue_duplications(found,routineP, "Lennard-Jones", name_atm_a, name_atm_b, error=error)
                found = .TRUE.
             END IF
          END IF

          ! loop over params from AMBER
          IF(ASSOCIATED(amb_info%nonbond_a)) THEN
             ii = 0
             jj = 0
             DO k=1,SIZE(amb_info%nonbond_a)
                IF((name_atm_a)==(amb_info%nonbond_a(k))) THEN
                   ii = k
                END IF
             END DO
             DO k=1,SIZE(amb_info%nonbond_a)
                IF((name_atm_b)==(amb_info%nonbond_a(k))) THEN
                   jj = k
                END IF
             END DO

             IF(ii/=0 .AND. jj/=0) THEN
                rmin = amb_info%nonbond_rmin2(ii)+amb_info%nonbond_rmin2(jj)
                epsilon0 = SQRT(amb_info%nonbond_eps(ii)*amb_info%nonbond_eps(jj))
                CALL  pair_potential_lj_create(pot%set(1)%lj, error)
                pot%type       = lj_charmm_type
                pot%at1        = name_atm_a
                pot%at2        = name_atm_b
                pot%set(1)%lj%epsilon = epsilon0
                pot%set(1)%lj%sigma6  = 0.5_dp*rmin**6
                pot%set(1)%lj%sigma12 = 0.25_dp*rmin**12
                pot%rcutsq     = (10.0_dp*bohr)**2
                CALL issue_duplications(found,routineP, "Lennard-Jones", name_atm_a, name_atm_b, error=error)
                found = .TRUE.
             END IF
          END IF

          ! always have the input param last to overwrite all the other ones
          IF(ASSOCIATED(inp_info%nonbonded)) THEN
             DO k=1,SIZE(inp_info%nonbonded%pot)
                IF ((TRIM(inp_info%nonbonded%pot(k)%pot%at1)=="*").OR.&
                    (TRIM(inp_info%nonbonded%pot(k)%pot%at2)=="*")) CYCLE

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
             ! Check for wildcards for one of the two types (if not associated yet)
             IF (.not.found) THEN
                DO k=1,SIZE(inp_info%nonbonded%pot)
                   IF ((TRIM(inp_info%nonbonded%pot(k)%pot%at1)=="*").EQV.&
                       (TRIM(inp_info%nonbonded%pot(k)%pot%at2)=="*")) CYCLE

                   IF(iw>0) WRITE(iw,*) "    TESTING ",TRIM(name_atm_a),TRIM(name_atm_b),&
                        " with ",TRIM(inp_info%nonbonded%pot(k)%pot%at1),&
                        TRIM(inp_info%nonbonded%pot(k)%pot%at2)

                   IF( (name_atm_a==inp_info%nonbonded%pot(k)%pot%at1) .OR.&
                       (name_atm_b==inp_info%nonbonded%pot(k)%pot%at2) .OR.&
                       (name_atm_b==inp_info%nonbonded%pot(k)%pot%at1) .OR.&
                       (name_atm_a==inp_info%nonbonded%pot(k)%pot%at2) ) THEN
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
                      IF(iw>0) WRITE(iw,*) "    FOUND (one WILDCARD)",TRIM(name_atm_a)," ",TRIM(name_atm_b)
                      found = .TRUE.
                   END IF
                END DO
             END IF
             ! Check for wildcards for both types (if not associated yet)
             IF (.not.found) THEN
                DO k=1,SIZE(inp_info%nonbonded%pot)
                   IF ((TRIM(inp_info%nonbonded%pot(k)%pot%at1)/="*").OR.&
                       (TRIM(inp_info%nonbonded%pot(k)%pot%at2)/="*")) CYCLE

                   IF(iw>0) WRITE(iw,*) "    TESTING ",TRIM(name_atm_a),TRIM(name_atm_b),&
                        " with ",TRIM(inp_info%nonbonded%pot(k)%pot%at1),&
                        TRIM(inp_info%nonbonded%pot(k)%pot%at2)

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
                   IF(iw>0) WRITE(iw,*) "    FOUND (both WILDCARDS)",TRIM(name_atm_a)," ",TRIM(name_atm_b)
                   found = .TRUE.
                END DO
             END IF
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
             pot%rcutsq = ff_type%rcut_nb*ff_type%rcut_nb
          END IF
          ! Cutoff is defined always as the maximum between the FF and Ewald
          pot%rcutsq = MAX(pot%rcutsq, ewald_rcut*ewald_rcut)
          ! Set the shell type
          IF((is_a_shell .AND. .NOT. is_b_shell) .OR. (is_b_shell .AND. .NOT. is_a_shell)) THEN
             pot%shell_type  = nosh_sh
          ELSE IF (is_a_shell .AND. is_b_shell) THEN
             pot%shell_type  = sh_sh
          ELSE
             pot%shell_type  = nosh_nosh
          END IF
          IF (only_qm) THEN
             CALL pair_potential_single_clean(pot, error=error)
          END IF
       END DO  ! jkind
    END DO   ! ikind
    CALL timestop(handle2)
  END SUBROUTINE force_field_pack_nonbond

! *****************************************************************************
!> \brief create the pair potential spline environment
!> \param atomic_kind_set ...
!> \param ff_type ...
!> \param iw2 ...
!> \param iw3 ...
!> \param iw4 ...
!> \param potparm ...
!> \param do_zbl ...
!> \param nonbonded_type ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE force_field_pack_splines(atomic_kind_set, ff_type, iw2, iw3, iw4, &
       potparm, do_zbl, nonbonded_type, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    INTEGER                                  :: iw2, iw3, iw4
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    LOGICAL, INTENT(IN)                      :: do_zbl
    CHARACTER(LEN=*), INTENT(IN)             :: nonbonded_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_splines', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle2, ikind, jkind, n
    LOGICAL                                  :: failure
    TYPE(spline_data_p_type), DIMENSION(:), &
      POINTER                                :: spl_p
    TYPE(spline_environment_type), POINTER   :: spline_env

    CALL timeset(routineN,handle2)
    failure = .FALSE.
    ! Figure out which nonbonded interactions happen to be indentical, and
    ! prepare storage for these, avoiding duplicates.
    NULLIFY(spline_env)
    CALL get_nonbond_storage(spline_env, potparm, atomic_kind_set, &
         do_zbl, shift_cutoff=ff_type%shift_cutoff, error=error)
    ! Effectively compute the spline data.
    CALL spline_nonbond_control(spline_env, potparm, &
         atomic_kind_set, eps_spline=ff_type%eps_spline, &
         max_energy=ff_type%max_energy, rlow_nb=ff_type%rlow_nb, &
         emax_spline=ff_type%emax_spline, npoints=ff_type%npoints, iw=iw2, iw2=iw3, iw3=iw4, &
         do_zbl=do_zbl, shift_cutoff=ff_type%shift_cutoff, &
         nonbonded_type=nonbonded_type, error=error)
    ! Let the pointers on potparm point to the splines generated in
    ! spline_nonbond_control.
    DO ikind = 1, SIZE ( potparm%pot, 1 )
       DO jkind = ikind, SIZE ( potparm%pot, 2)
          n = spline_env % spltab ( ikind, jkind )
          spl_p => spline_env%spl_pp(n)%spl_p
          CALL spline_data_p_retain ( spl_p, error )
          CALL spline_data_p_release ( potparm%pot(ikind,jkind)%pot%pair_spline_data, error )
          potparm%pot(ikind,jkind)%pot%pair_spline_data => spl_p
       END DO
    END DO
    CALL spline_env_release(spline_env,error)
    CALL timestop(handle2)

  END SUBROUTINE force_field_pack_splines

! *****************************************************************************
!> \brief Compute the electrostatic interaction cutoffs
!> \param atomic_kind_set ...
!> \param ff_type ...
!> \param potparm_nonbond ...
!> \param ewald_env ...
!> \param error ...
!> \author Toon.Verstraelen@gmail.com
! *****************************************************************************
  SUBROUTINE force_field_pack_eicut(atomic_kind_set, ff_type, &
      potparm_nonbond, ewald_env, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(force_field_type), INTENT(IN)       :: ff_type
    TYPE(pair_potential_pp_type), POINTER    :: potparm_nonbond
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack_eicut', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ewald_type, handle, i1, i2, &
                                                nkinds, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, beta, mm_radius1, &
                                                mm_radius2, rcut2, &
                                                rcut2_ewald, tmp
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: interaction_cutoffs
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    CALL timeset(routineN,handle)
    failure = .FALSE.

    tmp = 0.0_dp
    nkinds = SIZE(atomic_kind_set)
    ! allocate the array with interaction cutoffs for the electrostatics, used
    ! to make the electrostatic interaction continuous at ewald_env%rcut
    ALLOCATE(interaction_cutoffs(3,nkinds,nkinds), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    interaction_cutoffs = 0.0_dp

    ! compute the interaction cutoff if SHIFT_CUTOFF is active
    IF (ff_type%shift_cutoff) THEN
       CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type, &
            rcut=rcut2_ewald)
       rcut2_ewald = rcut2_ewald*rcut2_ewald
       DO i1 = 1, nkinds
          atomic_kind => atomic_kind_set(i1)
          CALL get_atomic_kind(atomic_kind=atomic_kind, mm_radius=mm_radius1)
          DO i2 = 1, nkinds
             rcut2=rcut2_ewald
             IF (ASSOCIATED(potparm_nonbond)) THEN
                rcut2 = MAX(potparm_nonbond%pot(i1,i2)%pot%rcutsq, rcut2_ewald)
             END IF
             IF (rcut2 > 0) THEN
                atomic_kind => atomic_kind_set(i2)
                CALL get_atomic_kind(atomic_kind=atomic_kind, mm_radius=mm_radius2)
                ! cutoff for core-core
                interaction_cutoffs(1,i1,i2) = potential_coulomb(rcut2, tmp, &
                  1.0_dp, ewald_type, alpha, 0.0_dp, 0.0_dp)
                ! cutoff for core-shell, core-ion, shell-core or ion-core
                IF (mm_radius1 > 0.0_dp) THEN
                   beta = sqrthalf/mm_radius1
                ELSE
                   beta = 0.0_dp
                END IF
                interaction_cutoffs(2,i1,i2) = potential_coulomb(rcut2, tmp, &
                  1.0_dp, ewald_type, alpha, beta, 0.0_dp)
                ! cutoff for shell-shell or ion-ion
                IF (mm_radius1 + mm_radius2 > 0.0_dp) THEN
                   beta = sqrthalf/SQRT(mm_radius1*mm_radius1+mm_radius2*mm_radius2)
                ELSE
                   beta = 0.0_dp
                END IF
                interaction_cutoffs(3,i1,i2) = potential_coulomb(rcut2, tmp, &
                  1.0_dp, ewald_type, alpha, beta, 0.0_dp)
             END IF
          END DO
       END DO
    END IF

    CALL ewald_env_set(ewald_env, interaction_cutoffs=interaction_cutoffs)

    CALL timestop(handle)
  END SUBROUTINE force_field_pack_eicut

! *****************************************************************************
!> \brief Issues on screen a warning when repetitions are present in the
!>        definition of the forcefield
!> \param found ...
!> \param routinePext ...
!> \param tag_label ...
!> \param name_atm_a ...
!> \param name_atm_b ...
!> \param name_atm_c ...
!> \param name_atm_d ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  SUBROUTINE issue_duplications(found, routinePext, tag_label, name_atm_a, name_atm_b,&
       name_atm_c, name_atm_d, error)

    LOGICAL, INTENT(IN)                      :: found
    CHARACTER(LEN=*), INTENT(IN)             :: routinePext, tag_label, &
                                                name_atm_a
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name_atm_b, name_atm_c, &
                                                name_atm_d
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'issue_duplications', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: item
    LOGICAL                                  :: failure

    failure = .FALSE.
    item="( "//TRIM(name_atm_a)
    IF (PRESENT(name_atm_b)) THEN
       item = TRIM(item)//" , "//TRIM(name_atm_b)
    END IF
    IF (PRESENT(name_atm_c)) THEN
       item = TRIM(item)//" , "//TRIM(name_atm_c)
    END IF
    IF (PRESENT(name_atm_d)) THEN
       item = TRIM(item)//" , "//TRIM(name_atm_d)
    END IF
    item = TRIM(item)//" )"
    CALL cp_assert(.NOT.found,cp_warning_level,cp_assertion_failed,routinePext,&
         "Multiple "//TRIM(tag_label)//" declarations: "//TRIM(item)//" overwriting! "//&
CPSourceFileRef,&
         only_ionode=.TRUE.)

  END SUBROUTINE issue_duplications

! *****************************************************************************
!> \brief Store informations on possible missing ForceFields parameters
!> \param atm1 ...
!> \param atm2 ...
!> \param atm3 ...
!> \param atm4 ...
!> \param type_name ...
!> \param fatal ...
!> \param array ...
!> \param error ...
! *****************************************************************************
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

            CASE("Out of plane bend")
                IF (INDEX(array(i)(21:39),"Out of plane bend")==0) CYCLE
                my_atm1 = array(i)(41:44)
                my_atm2 = array(i)(46:49)
                my_atm3 = array(i)(51:54)
                my_atm4 = array(i)(56:59)
                CALL compress(my_atm1,.TRUE.)
                CALL compress(my_atm2,.TRUE.)
                CALL compress(my_atm3,.TRUE.)
                CALL compress(my_atm4,.TRUE.)
    IF  (((atm1==my_atm1).AND.(atm2==my_atm2).AND.(atm3==my_atm3).AND.(atm4==my_atm4)).OR.&
         ((atm1==my_atm1).AND.(atm2==my_atm3).AND.(atm3==my_atm2).AND.(atm4==my_atm4)))&
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

