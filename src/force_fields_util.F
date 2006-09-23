!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****** cp2k/force_fields_util [1.0] *
!!
!!   NAME
!!     force_fields_util
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     Subroutine input_torsions changed (DG) 05-Dec-2000
!!     Output formats changed (DG) 05-Dec-2000
!!     JGH (26-01-2002) : force field parameters stored in tables, not in
!!       matrices. Input changed to have parameters labeled by the position
!!       and not atom pairs (triples etc)
!!     Teo (11.2005) : Moved all information on force field  pair_potential to
!!                     a much lighter memory structure
!!
!!   SOURCE
!******************************************************************************

MODULE force_fields_util

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             set_atomic_kind
  USE checkpoint_handler,              ONLY: Pwarning
  USE colvar_types,                    ONLY: dist_colvar_id
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE ewald_environment_types,         ONLY: ewald_environment_type
  USE external_potential_types,        ONLY: elp_potential_type,&
                                             get_potential,&
                                             set_potential
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_create,&
                                             fist_nonbond_env_type
  USE force_field_types,               ONLY: &
       allocate_bend_kind_set, allocate_bond_kind_set, &
       allocate_impr_kind_set, allocate_onfo_kind_set, &
       allocate_torsion_kind_set, allocate_ub_kind_set, bend_kind_type, &
       bond_kind_type, charmm_info_type, deallocate_bend_kind_set, &
       deallocate_bond_kind_set, force_field_type, gromos_info_type, &
       impr_kind_dealloc_ref, impr_kind_type, input_info_type, &
       onfo_kind_type, torsion_kind_dealloc_ref, torsion_kind_type, &
       ub_kind_dealloc_ref, ub_kind_type
  USE input_constants,                 ONLY: do_qmmm_none
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE manybody_potential,              ONLY: spline_manybody_control
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: &
       atom_type, bend_type, bond_type, colvar_constraint_type, &
       g3x3_constraint_type, g4x6_constraint_type, get_molecule_kind, &
       impr_type, molecule_kind_type, onfo_type, set_molecule_kind, &
       torsion_type, ub_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE pair_potential,                  ONLY: get_nonbond_storage,&
                                             spline_nonbond_control
  USE pair_potential_types,            ONLY: &
       ea_type, lj_charmm_type, lj_type, nn_type, pair_potential_lj_create, &
       pair_potential_pp_create, pair_potential_pp_release, &
       pair_potential_pp_type, pair_potential_single_clean, &
       pair_potential_single_copy
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: bohr
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE splines,                         ONLY: spline_data_release,&
                                             spline_data_retain,&
                                             spline_data_type,&
                                             spline_env_create,&
                                             spline_env_release,&
                                             spline_environment_type
  USE string_utilities,                ONLY: compress,&
                                             integer_to_string,&
                                             uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_fields_util'

  PRIVATE
  PUBLIC :: force_field_pack,&
            force_field_qeff_output,&
            clean_intra_force_kind

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** force_fields_util/force_field_pack [1.0] *
!!
!!   NAME
!!     force_field_pack
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Pack in all the information needed for the force_field
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

  SUBROUTINE force_field_pack (particle_set,atomic_kind_set,&
       molecule_kind_set,molecule_set,&
       ewald_env,fist_nonbond_env,ff_type,root_section,&
       qmmm, qmmm_env, mm_section, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(section_vals_type), POINTER         :: root_section
    LOGICAL, INTENT(IN), OPTIONAL            :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(INOUT)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length) :: atmname, ldum, name_atm_a, &
      name_atm_a2, name_atm_a_local, name_atm_b, name_atm_b2, &
      name_atm_b_local, name_atm_c, name_atm_c2, name_atm_d, name_atm_d2
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: Ainfo
    INTEGER :: atm_a, atm_b, atm_c, atm_d, counter, first, handle, handle2, &
      i, ii, ikind, ilink, imul, itype, iw, iw2, iw3, iw4, j, jj, jkind, k, &
      last, match_names, n, natom, nbend, nbond, nimpr, nkinds, nonfo, &
      ntorsion, nub, output_unit, stat
    INTEGER, DIMENSION(:), POINTER           :: molecule_list, my_atom_list
    INTEGER, POINTER :: map_bend_kind(:), map_bond_kind(:), map_impr_kind(:), &
      map_onfo_kind(:), map_torsion_kind(:), map_ub_kind(:)
    LOGICAL :: allocate_spline_bond_env, failure, fatal, found, found_a, &
      found_b, ionode, is_link_atom, my_qmmm, only_qm, use_qmmm_ff
    REAL(KIND=dp)                            :: epsilon, epsilon_a, &
                                                epsilon_b, qi1, qi2, qi3, &
                                                qi4, rmin, rmin2_a, rmin2_b, &
                                                scale_factor
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charge
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(bend_kind_type), DIMENSION(:), &
      POINTER                                :: bend_kind_set
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list
    TYPE(bond_kind_type), DIMENSION(:), &
      POINTER                                :: bond_kind_set
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(charmm_info_type), POINTER          :: chm_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(gromos_info_type), POINTER          :: gro_info
    TYPE(impr_kind_type), DIMENSION(:), &
      POINTER                                :: impr_kind_set
    TYPE(impr_type), DIMENSION(:), POINTER   :: impr_list
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(onfo_kind_type), DIMENSION(:), &
      POINTER                                :: onfo_kind_set
    TYPE(onfo_type), DIMENSION(:), POINTER   :: onfo_list
    TYPE(pair_potential_pp_type), POINTER    :: potparm_bond, potparm_nonbond
    TYPE(spline_data_type), POINTER          :: spline_data
    TYPE(spline_environment_type), POINTER   :: spline_bond_env, &
                                                spline_manybody_env, &
                                                spline_nonbond_env
    TYPE(torsion_kind_type), DIMENSION(:), &
      POINTER                                :: torsion_kind_set
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: torsion_list
    TYPE(ub_kind_type), DIMENSION(:), &
      POINTER                                :: ub_kind_set
    TYPE(ub_type), DIMENSION(:), POINTER     :: ub_list

    CALL timeset(routineN,'I','',handle)
    failure = .FALSE.
    fatal   = .FALSE.
    NULLIFY(logger, Ainfo)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
      output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
      output_unit = -1
    END IF    
    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)
    iw2= cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO/SPLINE_INFO",&
         extension=".mmLog",error=error)
    iw3= cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO/SPLINE_DATA",&
         extension=".mmLog",error=error)
    iw4= cp_print_key_unit_nr(logger,mm_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".mmLog",error=error)
    NULLIFY(atomic_kind,molecule_kind,bond_list,bend_list,&
         torsion_list,onfo_list,map_bond_kind,map_bend_kind,map_ub_kind,&
         map_torsion_kind,bond_kind_set,bend_kind_set,ub_kind_set,&
         onfo_kind_set,impr_kind_set,torsion_kind_set,my_atom_list,&
         potparm_bond, potparm_nonbond)
    only_qm = .FALSE.
    my_qmmm = .FALSE.
    use_qmmm_ff = qmmm_env%use_qmmm_ff
    allocate_spline_bond_env = .FALSE.
    IF (PRESENT(qmmm).AND.PRESENT(qmmm_env)) my_qmmm = qmmm
    inp_info => ff_type%inp_info
    chm_info => ff_type%chm_info
    gro_info => ff_type%gro_info
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Determine the number of unique bond kind and allocate bond_kind_set
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_bond1','I','',handle2)
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
         IF((ff_type%ff_type=="G96").OR.(ff_type%ff_type=="G87")) THEN
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
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Determine the number of unique bend kind and allocate bend_kind_set
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_bend1',handle2)
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
         IF((ff_type%ff_type=="G96").OR.(ff_type%ff_type=="G87")) THEN
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
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Determine the number of unique Urey-Bradley kind and allocate ub_kind_set
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_UB1',handle2)
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
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Determine the number of unique torsion kind and allocate torsion_kind_set
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_torsion1',handle2)
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
         IF((ff_type%ff_type=="G96").OR.(ff_type%ff_type=="G87")) THEN
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
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Determine the number of unique impr kind and allocate impr_kind_set
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_impr1',handle2)
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
         IF((ff_type%ff_type=="G96").OR.(ff_type%ff_type=="G87")) THEN
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
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Determine the number of unique onfo kind and allocate onfo_kind_set
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_onfo1',handle2)
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
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 7. BONDS
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_bonds2',handle2)
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
             READ(bond_list(j)%name(6:80),*) itype
             IF(itype<=k) THEN
               bond_list(j)%bond_kind%k(1) = gro_info%bond_k(itype)
               bond_list(j)%bond_kind%r0   = gro_info%bond_r0(itype)
               IF(ff_type%ff_type=="G87") bond_list(j)%bond_kind%name = "G87"
               IF(ff_type%ff_type=="G96") bond_list(j)%bond_kind%name = "G96"
               IF(ff_type%ff_type=="G87") bond_list(j)%name = "G87"
               IF(ff_type%ff_type=="G96") bond_list(j)%name = "G96"
               found = .TRUE.
             ELSE
               !k=SIZE(gro_info%solvent_k)
               !bond_list(j)%bond_kind%k(1) = gro_info%solvent_k(itype/k)
               !bond_list(j)%bond_kind%r0   = gro_info%solvent_r0(itype/k)
               itype=itype-k
               bond_list(j)%bond_kind%k(1) = gro_info%solvent_k(itype)
               bond_list(j)%bond_kind%r0   = gro_info%solvent_r0(itype)
               IF(ff_type%ff_type=="G87") bond_list(j)%bond_kind%name = "G87"
               IF(ff_type%ff_type=="G96") bond_list(j)%bond_kind%name = "G96"
               IF(ff_type%ff_type=="G87") bond_list(j)%name = "G87"
               IF(ff_type%ff_type=="G96") bond_list(j)%name = "G96"
               found = .TRUE.
             END IF
          END IF

          ! loop over params from charmm
          IF(ASSOCIATED(chm_info%bond_a)) THEN
             DO k=1,SIZE(chm_info%bond_a)
                IF( (((chm_info%bond_a(k))==(name_atm_a)) .AND. &
                     ((chm_info%bond_b(k))==(name_atm_b))) .OR. &
                     (((chm_info%bond_a(k))==(name_atm_b)) .AND. &
                     ((chm_info%bond_b(k))==(name_atm_a))) ) THEN
                   bond_list(j)%bond_kind%name = "CHARMM"
                   bond_list(j)%bond_kind%k(1) = chm_info%bond_k(k)
                   bond_list(j)%bond_kind%r0   = chm_info%bond_r0(k)
                   IF (found) CALL Pwarning("Multiple Bond declarations: "//TRIM(name_atm_a)//&
                        " and "//TRIM(name_atm_b)//" overwriting!",root_section,error)
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
                   bond_list(j)%bond_kind%name = inp_info%bond_kind(k)
                   bond_list(j)%bond_kind%k(:) = inp_info%bond_k(:,k)
                   bond_list(j)%bond_kind%r0   = inp_info%bond_r0(k)
                   IF (found) CALL Pwarning("Multiple Bond declarations: "//TRIM(name_atm_a)//&
                        " and "//TRIM(name_atm_b)//" overwriting!",root_section,error)
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
          !
          ! QM/MM modifications
          !
          IF (only_qm) THEN
             bond_list(j)%name           = "UNSET"
             bond_list(j)%bond_kind%name = "UNSET"
          END IF
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            bond_list=bond_list)
    END DO
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 8. BENDS
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_bends2',handle2)
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
             READ(bend_list(j)%name(6:80),*) itype
             IF(itype>0) THEN
               bend_list(j)%bend_kind%k      = gro_info%bend_k(itype)
               bend_list(j)%bend_kind%theta0 = gro_info%bend_theta0(itype)
               IF(ff_type%ff_type=="G87") bend_list(j)%bend_kind%name   = "G87"
               IF(ff_type%ff_type=="G96") bend_list(j)%bend_kind%name   = "G96"
               IF(ff_type%ff_type=="G87") bend_list(j)%name   = "G87"
               IF(ff_type%ff_type=="G96") bend_list(j)%name   = "G96"
               found = .TRUE.
             ELSE
               bend_list(j)%bend_kind%k      = gro_info%bend_k(itype/k)
               bend_list(j)%bend_kind%theta0 = gro_info%bend_theta0(itype/k)
               IF(ff_type%ff_type=="G87") bend_list(j)%bend_kind%name   = "G87"
               IF(ff_type%ff_type=="G96") bend_list(j)%bend_kind%name   = "G96"
               IF(ff_type%ff_type=="G87") bend_list(j)%name   = "G87"
               IF(ff_type%ff_type=="G96") bend_list(j)%name   = "G96"
               found = .TRUE.
             END IF
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
                   bend_list(j)%bend_kind%name   = "CHARMM"
                   bend_list(j)%bend_kind%k      = chm_info%bend_k(k)
                   bend_list(j)%bend_kind%theta0 = chm_info%bend_theta0(k)
                   IF (found) CALL Pwarning("Multiple Bend declarations: "//TRIM(name_atm_a)//&
                        ","//TRIM(name_atm_b)//" and "//TRIM(name_atm_c)//" overwriting!",root_section,error)
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
                   bend_list(j)%bend_kind%name   = inp_info%bend_kind(k)
                   bend_list(j)%bend_kind%k      = inp_info%bend_k(k)
                   bend_list(j)%bend_kind%theta0 = inp_info%bend_theta0(k)
                   IF (found) CALL Pwarning("Multiple Bend declarations: "//TRIM(name_atm_a)//&
                        ","//TRIM(name_atm_b)//" and "//TRIM(name_atm_c)//" overwriting!",root_section,error)
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
          !
          ! QM/MM modifications
          !
          IF (only_qm) THEN
             bend_list(j)%name           = "UNSET"
             bend_list(j)%bend_kind%name = "UNSET"
          END IF
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            bend_list=bend_list)
    END DO
    ! Give information and abort if any bond or angle parameter is missing..
    CALL release_FF_missing_par(fatal,AInfo,output_unit,iw,error)
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 9. Urey-Bradley
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_UB2',handle2)
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
                   ub_list(j)%ub_kind%name   = "CHARMM"
                   ub_list(j)%ub_kind%k      = chm_info%ub_k(k)
                   ub_list(j)%ub_kind%r0     = chm_info%ub_r0(k)
                   IF(iw>0) WRITE(iw,*) "    Found UB ",TRIM(name_atm_a)," ",&
                        TRIM(name_atm_b)," ",TRIM(name_atm_c)
                   IF (found) CALL Pwarning("Multiple UB declarations: "//TRIM(name_atm_a)//&
                        ","//TRIM(name_atm_b)//" and "//TRIM(name_atm_c)//" overwriting!",root_section,error)
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
                   ub_list(j)%ub_kind%name   = inp_info%ub_kind(k)
                   ub_list(j)%ub_kind%k      = inp_info%ub_k(k)
                   ub_list(j)%ub_kind%r0     = inp_info%ub_r0(k)
                   IF (found) CALL Pwarning("Multiple UB declarations: "//TRIM(name_atm_a)//&
                        ","//TRIM(name_atm_b)//" and "//TRIM(name_atm_c)//" overwriting!",root_section,error)
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
             ub_list(j)%name           = "UNSET"
             ub_list(j)%ub_kind%name   = "UNSET"
             ub_list(j)%ub_kind%k      = 0.0_dp
             ub_list(j)%ub_kind%r0     = 0.0_dp
          END IF
          !
          ! QM/MM modifications
          !
          IF (only_qm) THEN
             ub_list(j)%name         = "UNSET"
             ub_list(j)%ub_kind%name = "UNSET"
          END IF
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            ub_list=ub_list)
    END DO
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 10. TORSION
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_torsion2',handle2)
    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            ntorsion=ntorsion,torsion_list=torsion_list)
       molecule=>molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
       DO j=1,ntorsion
          IF(torsion_list(j)%torsion_kind%name == "UNDEF") THEN
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
               READ(torsion_list(j)%name(6:80),*) itype
               IF(itype>0) THEN
                 CALL reallocate(torsion_list(j)%torsion_kind%k,1,1)
                 CALL reallocate(torsion_list(j)%torsion_kind%m,1,1)
                 CALL reallocate(torsion_list(j)%torsion_kind%cosphi0,1,1)
                 torsion_list(j)%torsion_kind%nmul   = 1
                 torsion_list(j)%torsion_kind%m(1)   = gro_info%torsion_m(itype)
                 torsion_list(j)%torsion_kind%k(1)   = gro_info%torsion_k(itype)
                 torsion_list(j)%torsion_kind%cosphi0(1) = gro_info%torsion_cosphi0(itype)
                 IF(ff_type%ff_type=="G87") torsion_list(j)%torsion_kind%name   = "G87"
                 IF(ff_type%ff_type=="G96") torsion_list(j)%torsion_kind%name   = "G96"
                 IF(ff_type%ff_type=="G87") torsion_list(j)%name   = "G87"
                 IF(ff_type%ff_type=="G96") torsion_list(j)%name   = "G96"
                 found = .TRUE.
               ELSE
                 CALL reallocate(torsion_list(j)%torsion_kind%k,1,1)
                 CALL reallocate(torsion_list(j)%torsion_kind%m,1,1)
                 CALL reallocate(torsion_list(j)%torsion_kind%cosphi0,1,1)
                 torsion_list(j)%torsion_kind%nmul   = 1
                 torsion_list(j)%torsion_kind%m(1)   = gro_info%torsion_m(itype/k)
                 torsion_list(j)%torsion_kind%k(1)   = gro_info%torsion_k(itype/k)
                 torsion_list(j)%torsion_kind%cosphi0(1) = gro_info%torsion_cosphi0(itype/k)
                 IF(ff_type%ff_type=="G87") torsion_list(j)%torsion_kind%name   = "G87"
                 IF(ff_type%ff_type=="G96") torsion_list(j)%torsion_kind%name   = "G96"
                 IF(ff_type%ff_type=="G87") torsion_list(j)%name   = "G87"
                 IF(ff_type%ff_type=="G96") torsion_list(j)%name   = "G96"
                 found = .TRUE.
               END IF
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
                      torsion_list(j)%torsion_kind%name    = "CHARMM"
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
                         torsion_list(j)%torsion_kind%name    = "CHARMM"
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
                      torsion_list(j)%torsion_kind%name    = inp_info%torsion_kind(k)
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
                torsion_list(j)%torsion_kind%name = "UNSET"
                torsion_list(j)%name              = "UNSET"
             ELSE
                ldum = cp_to_string(imul)
                IF ((imul /= 1).AND.(iw>0))&
                     WRITE(iw,'(/,2("UTIL_INFO| ",A,/))')&
                     "Multiple Torsion declarations: "//TRIM(name_atm_a)//&
                     ","//TRIM(name_atm_b)//","//TRIM(name_atm_c)//" and "//TRIM(name_atm_d),&
                     "Present torsion multiplicity "//TRIM(ldum)//" ."
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
                torsion_list(j)%torsion_kind%name = "UNSET"
                torsion_list(j)%name              = "UNSET"
             END IF
          END IF
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            torsion_list=torsion_list)
    END DO
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 11. IMPROPERS
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_impr2',handle2)
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
            READ(impr_list(j)%name(6:80),*) itype
            IF(itype>0) THEN
              impr_list(j)%impr_kind%k    = gro_info%impr_k(itype)
              impr_list(j)%impr_kind%phi0 = gro_info%impr_phi0(itype)
              IF(ff_type%ff_type=="G87") impr_list(j)%impr_kind%name = "G87"
              IF(ff_type%ff_type=="G96") impr_list(j)%impr_kind%name = "G96"
              IF(ff_type%ff_type=="G87") impr_list(j)%name = "G87"
              IF(ff_type%ff_type=="G96") impr_list(j)%name = "G96"
              found = .TRUE.
            ELSE
              impr_list(j)%impr_kind%k    = gro_info%impr_k(itype)
              impr_list(j)%impr_kind%phi0 = gro_info%impr_phi0(itype)
              IF(ff_type%ff_type=="G87") impr_list(j)%impr_kind%name = "G87"
              IF(ff_type%ff_type=="G96") impr_list(j)%impr_kind%name = "G96"
              IF(ff_type%ff_type=="G87") impr_list(j)%name = "G87"
              IF(ff_type%ff_type=="G96") impr_list(j)%name = "G96"
              found = .TRUE.
            END IF
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
                   impr_list(j)%impr_kind%name = "CHARMM"
                   impr_list(j)%impr_kind%k    = chm_info%impr_k(k)
                   impr_list(j)%impr_kind%phi0 = chm_info%impr_phi0(k)
                   IF (found) CALL Pwarning("Multiple Improper declarations: "//TRIM(name_atm_a)//&
                        ","//TRIM(name_atm_b)//","//TRIM(name_atm_c)//" and "//TRIM(name_atm_d)//&
                        " overwriting!",root_section,error)
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
                      impr_list(j)%impr_kind%name = "CHARMM"
                      impr_list(j)%impr_kind%k    = chm_info%impr_k(k)
                      impr_list(j)%impr_kind%phi0 = chm_info%impr_phi0(k)
                      IF (found) CALL Pwarning("Multiple Improper declarations: "//TRIM(name_atm_a)//&
                           ","//TRIM(name_atm_b)//","//TRIM(name_atm_c)//" and "//TRIM(name_atm_d)//&
                           " overwriting!",root_section,error)
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
                   impr_list(j)%impr_kind%name = inp_info%impr_kind(k)
                   impr_list(j)%impr_kind%k    = inp_info%impr_k(k)
                   impr_list(j)%impr_kind%phi0 = inp_info%impr_phi0(k)
                   IF (found) CALL Pwarning("Multiple Improper declarations: "//TRIM(name_atm_a)//&
                        ","//TRIM(name_atm_b)//","//TRIM(name_atm_c)//" and "//TRIM(name_atm_d)//&
                        " overwriting!",root_section,error)
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
             impr_list(j)%impr_kind%k      = 0.0_dp
             impr_list(j)%impr_kind%phi0   = 0.0_dp
             impr_list(j)%impr_kind%name   = "UNSET"
             impr_list(j)%name             = "UNSET"
          END IF
          !
          ! QM/MM modifications
          !
          IF (only_qm) THEN
             impr_list(j)%impr_kind%name = "UNSET"
             impr_list(j)%name           = "UNSET"
          END IF

       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            impr_list=impr_list)
    END DO
    ! Give information only if any Urey-Bradley, Torsion or improper is missing
    ! continue calculation..
    CALL release_FF_missing_par(fatal,AInfo,output_unit,iw,error)
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    !12. Set up atomic_kind_set()%elp_potentail%qeff (PART 2)
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_qeff',handle2)
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
       IF(ASSOCIATED(ff_type%chm_info%charge_atm)) THEN
          DO j=1,SIZE(ff_type%chm_info%charge_atm)
             IF((ff_type%chm_info%charge_atm(j))==atmname) THEN
                charge(i) = ff_type%chm_info%charge(j)
                IF (found) CALL Pwarning("Multiple charge declarations: "//TRIM(atmname)//&
                     " overwriting!",root_section,error)
                found = .TRUE.
             END IF
          END DO
       END IF
       ! always have the input param last to overwrite all the other ones
       IF(ASSOCIATED(ff_type%inp_info%charge_atm)) THEN
          DO j=1,SIZE(ff_type%inp_info%charge_atm)
             IF (iw>0) WRITE(iw,*)"Charge Checking ::",TRIM(ff_type%inp_info%charge_atm(j)),atmname
             IF((ff_type%inp_info%charge_atm(j))==atmname) THEN
                charge(i) = ff_type%inp_info%charge(j)
                IF (found) CALL Pwarning("Multiple charge declarations: "//TRIM(atmname)//&
                     " overwriting!",root_section,error)
                found = .TRUE.
             END IF
          END DO
       END IF
       IF (ASSOCIATED(ff_type%inp_info%nonbonded)) THEN
          IF (ASSOCIATED(ff_type%inp_info%nonbonded%pot)) THEN
             ! Let's find the nonbonded potential where this atom is involved
             DO j = 1, SIZE(ff_type%inp_info%nonbonded%pot)
                IF (atmname==ff_type%inp_info%nonbonded%pot(j)%pot %at1.OR.&
                     atmname==ff_type%inp_info%nonbonded%pot(j)%pot %at2) THEN
                   SELECT CASE(ff_type%inp_info%nonbonded%pot(j)%pot%type)
                   CASE (ea_type)
                      ! Nothing to complain.. we don't need charges
                   CASE default
                      IF (.NOT.found) CALL store_FF_missing_par(atm1=TRIM(name_atm_a),&
                                                                fatal=fatal,&
                                                                type_name="Charge",&
                                                                array=Ainfo,&
                                                                error=error)
                   END SELECT
                END IF
             END DO
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
    ! Give information only if charge is missing and abort..
    CALL release_FF_missing_par(fatal,AInfo,output_unit,iw,error)
    DEALLOCATE (charge,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 13. Set spline_bond_env, use for ONFO interaction
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_onfo2',handle2)
    CALL pair_potential_pp_create ( potparm_bond,  &
         SIZE ( atomic_kind_set ), &
         error )
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

          ! loop over params from GROMOS
          IF(ASSOCIATED(gro_info%nonbond_a_14)) THEN
             ii = 0
             jj = 0
             DO k=1,SIZE(gro_info%nonbond_a_14)
                IF((name_atm_a)==(gro_info%nonbond_a_14(k))) THEN
                   ii = k
                   found_a = .TRUE.
                   EXIT
                END IF
             END DO
             DO k=1,SIZE(gro_info%nonbond_a_14)
                IF((name_atm_b)==(gro_info%nonbond_a_14(k))) THEN
                   jj = k
                   found_b   = .TRUE.
                   EXIT
                END IF
             END DO
             IF(ii/=0 .AND. jj/=0) THEN
                CALL  pair_potential_lj_create(potparm_bond%pot(i,j)%pot, error)
                potparm_bond%pot(i,j)%pot%type       = lj_type
                potparm_bond%pot(i,j)%pot%at1        = name_atm_a
                potparm_bond%pot(i,j)%pot%at2        = name_atm_b
                potparm_bond%pot(i,j)%pot%lj%epsilon = 1.0_dp
                potparm_bond%pot(i,j)%pot%lj%sigma6  = &
                  SQRT(gro_info%nonbond_c6_14(ii)*gro_info%nonbond_c6_14(jj))
                potparm_bond%pot(i,j)%pot%lj%sigma12 = &
                  SQRT(gro_info%nonbond_c12_14(ii)*gro_info%nonbond_c12_14(jj))
                potparm_bond%pot(i,j)%pot%rcutsq     = (10.0_dp*bohr)**2
                IF (found) CALL Pwarning("Multiple Lennard-Jones declarations: "//TRIM(name_atm_a)//&
                     " and "//TRIM(name_atm_b)//" overwriting!",root_section,error)
                found = .TRUE.
             END IF
          END IF

          ! loop over params from CHARMM
          IF(ASSOCIATED(chm_info%nonbond_a_14)) THEN
             ii = 0
             jj = 0
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
             IF(ii/=0 .AND. jj/=0) THEN
                rmin = rmin2_a + rmin2_b
                epsilon = SQRT(epsilon_a*epsilon_b)
                CALL  pair_potential_lj_create(potparm_bond%pot(i,j)%pot, error)
                potparm_bond%pot(i,j)%pot%type       = lj_charmm_type
                potparm_bond%pot(i,j)%pot%at1        = name_atm_a
                potparm_bond%pot(i,j)%pot%at2        = name_atm_b
                potparm_bond%pot(i,j)%pot%lj%epsilon = epsilon
                potparm_bond%pot(i,j)%pot%lj%sigma6  = 0.5_dp*rmin**6
                potparm_bond%pot(i,j)%pot%lj%sigma12 = 0.25_dp*rmin**12
                potparm_bond%pot(i,j)%pot%rcutsq     = (10.0_dp*bohr)**2
                IF (found) CALL Pwarning("Multiple Lennard-Jones declarations: "//TRIM(name_atm_a)//&
                     " and "//TRIM(name_atm_b)//" overwriting!",root_section,error)
                found = .TRUE.
             END IF
          END IF

          ! always have the input param last to overwrite all the other ones
          IF(ASSOCIATED(inp_info%nonbonded14)) THEN
             DO k=1,SIZE(inp_info%nonbonded14%pot)
                IF(iw>0) WRITE(iw,*) "    TESTING ",TRIM(name_atm_a),TRIM(name_atm_b),&
                     " with ",TRIM(inp_info%nonbonded14%pot(k)%pot%at1),TRIM(inp_info%nonbonded14%pot(k)%pot%at2)
                IF((((name_atm_a)==(inp_info%nonbonded14%pot(k)%pot%at1))  .AND. &
                     ((name_atm_b)==(inp_info%nonbonded14%pot(k)%pot%at2))) .OR.  &
                     (((name_atm_b)==(inp_info%nonbonded14%pot(k)%pot%at1))  .AND. &
                     ((name_atm_a)==(inp_info%nonbonded14%pot(k)%pot%at2))) ) THEN
                   CALL pair_potential_single_copy(inp_info%nonbonded14%pot(k)%pot,potparm_bond%pot(i,j)%pot,error)
                   IF (found) CALL Pwarning("Multiple Lennard-Jones declarations: "//TRIM(name_atm_a)//&
                        " and "//TRIM(name_atm_b)//" overwriting!",root_section,error)
                   IF(iw>0) WRITE(iw,*) "    FOUND ",TRIM(name_atm_a)," ",TRIM(name_atm_b)
                   found = .TRUE.
                   EXIT
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
                         CALL pair_potential_single_copy(qmmm_env%inp_info%nonbonded14%pot(k)%pot,&
                              potparm_bond%pot(i,j)%pot,error)
                         IF (found) CALL Pwarning("Multiple Lennard-Jones declarations: "//TRIM(name_atm_a)//&
                              " and "//TRIM(name_atm_b)//" overwriting with QM/MM forcefield specifications",&
                              root_section,error)
                         IF(iw>0) WRITE(iw,*) "    FOUND ",TRIM(name_atm_a),&
                              " ",TRIM(name_atm_b)
                         found = .TRUE.
                         EXIT
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
             CALL pair_potential_single_clean( potparm_bond%pot(i,j)%pot ,error=error)
             potparm_bond%pot(i,j)%pot%type       = nn_type
             potparm_bond%pot(i,j)%pot%at1        = name_atm_a
             potparm_bond%pot(i,j)%pot%at2        = name_atm_b
          ELSE
             allocate_spline_bond_env = .TRUE.
          END IF
          ! If defined global RCUT let's use it
          IF(ff_type%rcut_nb>0.0_dp) THEN
             potparm_bond%pot(i,j)%pot%rcutsq     = ff_type%rcut_nb*ff_type%rcut_nb
          END IF
          IF (only_qm) THEN
             CALL pair_potential_single_clean( potparm_bond%pot(i,j)%pot ,error=error)
          END IF
       END DO
    END DO
    ! Give information if any 1-4 is missing.. continue calculation..
    CALL release_FF_missing_par(fatal,AInfo,output_unit,iw,error)
    NULLIFY(spline_bond_env)
    CALL get_nonbond_storage(spline_bond_env, potparm_bond, atomic_kind_set, &
          ff_type%ei_scale14,ff_type%vdw_scale14, error=error)
    CALL spline_nonbond_control ( spline_bond_env, potparm_bond, &
         atomic_kind_set , ewald_env, ff_type%ei_scale14,ff_type%vdw_scale14,&
         ff_type%eps_spline,max_energy=ff_type%max_energy, rlow_nb=ff_type%rlow_nb,&
         emax_spline=ff_type%emax_spline, iw=iw2, iw2=iw3, iw3=iw4,&
         do_14=.TRUE., scale_cutoff=ff_type%scale_cutoff, error=error)
    DO ikind = 1, SIZE ( potparm_bond%pot, 1 )
       DO jkind = ikind, SIZE ( potparm_bond%pot, 2 )
          n = spline_bond_env % spltab ( ikind, jkind )
          spline_data => spline_bond_env % spline_data_p ( n ) % spline_data
          !
          CALL spline_data_retain ( spline_data, error )
          CALL spline_data_release ( potparm_bond %pot( ikind, jkind )%pot % pair_spline_data, error )
          potparm_bond %pot( ikind, jkind )%pot % pair_spline_data => spline_data
       END DO
    END DO
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 14. ONFO
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_onfo3',handle2)
    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list,&
            natom=natom,&
            nonfo=nonfo,onfo_list=onfo_list)
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

          CALL spline_data_retain (potparm_bond%pot(ikind, jkind )%pot % pair_spline_data,error=error)
          CALL spline_data_release(onfo_list(j)%onfo_kind%spline_data, error=error)
          onfo_list(j)%onfo_kind%spline_data => &
               potparm_bond %pot( ikind, jkind )%pot % pair_spline_data

          IF(iw>0) WRITE(iw,*) "    ONE-FOUR ",j," : ",&
               TRIM(name_atm_a)," ",&
               TRIM(name_atm_b)
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            onfo_list=onfo_list)
    END DO
    CALL pair_potential_pp_release ( potparm_bond, error )
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 15. Set spline_nonbond_env
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_nonbond',handle2)
    CALL pair_potential_pp_create ( potparm_nonbond,  &
         SIZE ( atomic_kind_set ), &
         error )
    DO i=1,SIZE(atomic_kind_set)
       atomic_kind => atomic_kind_set(i)
       CALL get_atomic_kind(atomic_kind=atomic_kind,name=name_atm_a_local)
       DO j=i,SIZE(atomic_kind_set)
          atomic_kind => atomic_kind_set(j)
          CALL get_atomic_kind(atomic_kind=atomic_kind,name=name_atm_b_local)
          found = .FALSE.
          name_atm_a = name_atm_a_local
          name_atm_b = name_atm_b_local
          only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b)
          CALL uppercase(name_atm_a)
          CALL uppercase(name_atm_b)

          ! loop over params from GROMOS
          IF(ASSOCIATED(gro_info%nonbond_a)) THEN
             ii = 0
             jj = 0
             DO k=1,SIZE(gro_info%nonbond_a)
                IF((name_atm_a)==(gro_info%nonbond_a(k))) THEN
                   ii = k
                   EXIT
                END IF
             END DO
             DO k=1,SIZE(gro_info%nonbond_a)
                IF((name_atm_b)==(gro_info%nonbond_a(k))) THEN
                   jj = k
                   EXIT
                END IF
             END DO

             IF(ii/=0 .AND. jj/=0) THEN
                CALL  pair_potential_lj_create(potparm_nonbond%pot(i,j)%pot, error)
                potparm_nonbond%pot(i,j)%pot%type       = lj_type
                potparm_nonbond%pot(i,j)%pot%at1        = name_atm_a
                potparm_nonbond%pot(i,j)%pot%at2        = name_atm_b
                potparm_nonbond%pot(i,j)%pot%lj%epsilon = 1.0_dp
                potparm_nonbond%pot(i,j)%pot%lj%sigma6  = &
                  SQRT(gro_info%nonbond_c6(ii)*gro_info%nonbond_c6(jj))
                potparm_nonbond%pot(i,j)%pot%lj%sigma12 = &
                  SQRT(gro_info%nonbond_c12(ii)*gro_info%nonbond_c12(jj))
                potparm_nonbond%pot(i,j)%pot%rcutsq     = (10.0_dp*bohr)**2
                IF (found) CALL Pwarning("Multiple Lennard-Jones declarations: "//TRIM(name_atm_a)//&
                     " and "//TRIM(name_atm_b)//" overwriting!",&
                     root_section,error)
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
                epsilon = SQRT(chm_info%nonbond_eps(ii)*&
                     chm_info%nonbond_eps(jj))
                CALL  pair_potential_lj_create(potparm_nonbond%pot(i,j)%pot, error)
                potparm_nonbond%pot(i,j)%pot%type       = lj_charmm_type
                potparm_nonbond%pot(i,j)%pot%at1        = name_atm_a
                potparm_nonbond%pot(i,j)%pot%at2        = name_atm_b
                potparm_nonbond%pot(i,j)%pot%lj%epsilon = epsilon
                potparm_nonbond%pot(i,j)%pot%lj%sigma6  = 0.5_dp*rmin**6
                potparm_nonbond%pot(i,j)%pot%lj%sigma12 = 0.25_dp*rmin**12
                potparm_nonbond%pot(i,j)%pot%rcutsq     = (10.0_dp*bohr)**2
                IF (found) CALL Pwarning("Multiple Lennard-Jones declarations: "//TRIM(name_atm_a)//&
                     " and "//TRIM(name_atm_b)//" overwriting!",&
                     root_section,error)
                found = .TRUE.
             END IF
          END IF


          ! always have the input param last to overwrite all the other ones
          IF(ASSOCIATED(inp_info%nonbonded)) THEN
             DO k=1,SIZE(inp_info%nonbonded%pot)
                IF(iw>0) WRITE(iw,*) "    TESTING ",TRIM(name_atm_a),TRIM(name_atm_b),&
                     " with ",TRIM(inp_info%nonbonded%pot(k)%pot%at1),&
                     TRIM(inp_info%nonbonded%pot(k)%pot%at2)
                IF((((name_atm_a)==(inp_info%nonbonded%pot(k)%pot%at1))  .AND. &
                     ((name_atm_b)==(inp_info%nonbonded%pot(k)%pot%at2))) .OR.  &
                     (((name_atm_b)==(inp_info%nonbonded%pot(k)%pot%at1))  .AND. &
                     ((name_atm_a)==(inp_info%nonbonded%pot(k)%pot%at2))) ) THEN
                   CALL pair_potential_single_copy(inp_info%nonbonded%pot(k)%pot,potparm_nonbond%pot(i,j)%pot,error)
                   IF (found) CALL Pwarning("Multiple Lennard-Jones declarations: "//TRIM(name_atm_a)//&
                        " and "//TRIM(name_atm_b)//" overwriting!",&
                        root_section,error)
                   IF(iw>0) WRITE(iw,*) "    FOUND ",TRIM(name_atm_a),&
                        " ",TRIM(name_atm_b)
                   found = .TRUE.
                   EXIT
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
                         CALL pair_potential_single_copy(qmmm_env%inp_info%nonbonded%pot(k)%pot,&
                              potparm_nonbond%pot(i,j)%pot,error)
                         IF (found) CALL Pwarning("Multiple Lennard-Jones declarations: "//TRIM(name_atm_a)//&
                              " and "//TRIM(name_atm_b)//" overwriting with QM/MM forcefield specifications",&
                              root_section,error)
                         IF(iw>0) WRITE(iw,*) "    FOUND ",TRIM(name_atm_a),&
                              " ",TRIM(name_atm_b)
                         found = .TRUE.
                         EXIT
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
             potparm_nonbond%pot(i,j)%pot%rcutsq     = ff_type%rcut_nb*ff_type%rcut_nb
          END IF
          IF (only_qm) THEN
             CALL pair_potential_single_clean( potparm_nonbond%pot(i,j)%pot ,error=error)
          END IF
       END DO
    END DO
    ! Give information and abort if any pair potential spline is missing..
    CALL release_FF_missing_par(fatal,AInfo,output_unit,iw,error)
    !---------------------------------------
    ! create the pair potential spline environment
    !---------------------------------------
    NULLIFY(spline_nonbond_env)
    CALL get_nonbond_storage(spline_nonbond_env, potparm_nonbond, atomic_kind_set, &
         error=error)
    CALL spline_nonbond_control(spline_nonbond_env, potparm_nonbond, &
         atomic_kind_set,ewald_env,eps_spline=ff_type%eps_spline,&
         max_energy=ff_type%max_energy, rlow_nb=ff_type%rlow_nb,&
         emax_spline=ff_type%emax_spline, iw=iw2, iw2=iw3,iw3=iw4,&
         scale_cutoff=ff_type%scale_cutoff, error=error)
    !---------------------------------------
    ! create the manybody spline environment
    !---------------------------------------
    nkinds = SIZE(atomic_kind_set)
    CALL spline_env_create(spline_manybody_env, nkinds, error=error)
    CALL spline_manybody_control(spline_manybody_env, potparm_nonbond,&
         eps_spline=ff_type%eps_spline,max_energy=ff_type%max_energy,&
         rlow_nb=ff_type%rlow_nb,&
         emax_spline=ff_type%emax_spline, iw=iw2, iw2=iw3, error=error)

    DO ikind = 1, SIZE ( potparm_nonbond%pot, 1 )
       DO jkind = ikind, SIZE ( potparm_nonbond%pot, 2)
          !---------------------------------------
          ! 1st: The pair-potential contribution
          !---------------------------------------
          n = spline_nonbond_env % spltab ( ikind, jkind )
          spline_data => spline_nonbond_env % spline_data_p ( n ) % spline_data
          CALL spline_data_retain ( spline_data, error )
          CALL spline_data_release (  &
               potparm_nonbond %pot( ikind, jkind )%pot % pair_spline_data, error )
          potparm_nonbond%pot ( ikind, jkind )%pot % pair_spline_data => spline_data
          !---------------------------------------
          ! 2nd: The many-body contribution
          !---------------------------------------
          n = spline_manybody_env % spltab ( ikind, jkind )
          IF (n>0) THEN
             spline_data => spline_manybody_env % spline_data_p ( n ) % spline_data
             ! do i,j
             IF ( ASSOCIATED ( spline_data ) ) THEN
                CALL spline_data_retain ( spline_data, error )
                CALL spline_data_release (  &
                     potparm_nonbond%pot ( ikind, jkind )%pot % rho_spline_data, error )
                potparm_nonbond%pot ( ikind, jkind )%pot % rho_spline_data => spline_data
             ENDIF
          END IF
       END DO
    END DO
    CALL fist_nonbond_env_create (  fist_nonbond_env, potparm_nonbond, error )
    CALL timestop(handle2)

    ! Release all spline_environments
    CALL spline_env_release (spline_manybody_env,error)
    CALL spline_env_release (spline_nonbond_env,error)
    CALL spline_env_release (spline_bond_env,error)

    CALL cp_print_key_finished_output(iw4,logger,mm_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)
    CALL cp_print_key_finished_output(iw3,logger,mm_section,&
         "PRINT%FF_INFO/SPLINE_DATA",error=error)
    CALL cp_print_key_finished_output(iw2,logger,mm_section,&
         "PRINT%FF_INFO/SPLINE_INFO",error=error)
    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%FF_INFO",error=error)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE force_field_pack

!******************************************************************************
!!****** force_fields_util/release_FF_missing_par [1.0] *
!!
!!   NAME
!!     release_FF_missing_par
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
  SUBROUTINE release_FF_missing_par(fatal, array, output_unit, iw, error)
    LOGICAL, INTENT(INOUT), OPTIONAL         :: fatal
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: array
    INTEGER, INTENT(IN)                      :: output_unit, iw
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_FF_missing_par', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(array)) THEN
       IF (output_unit>0) THEN
          WRITE(output_unit,*)
          WRITE(output_unit,'(T2,"FORCEFIELD|",A)')&
               " WARNING: Non Critical ForceField Parameter missing!   Calculation proceeding.",&
               " Non Critical Parameters are:   Urey-Bradley,   Torsions,  Impropers  and  1-4",&
               " All missing parameters will give a zero contribution to the potential energy!"
          IF (fatal.OR.iw>0) THEN
             WRITE(output_unit,*)
             DO i = 1,SIZE(array)
                WRITE(output_unit,'(A)')array(i)
             END DO
          END IF
          IF (.NOT.fatal.AND.iw<0) THEN
             WRITE(output_unit,'(T2,"FORCEFIELD|",A)')&
                  " Activate the print key FF_INFO to have a list of missing parameters."
             WRITE(output_unit,*)
          END IF
       END IF
       DEALLOCATE(array,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (fatal) CALL stop_program(TRIM(routineN),"Missing critical ForceField parameters! Abort..")
  END SUBROUTINE release_FF_missing_par
!******************************************************************************
!!****** force_fields_util/store_FF_missing_par [1.0] *
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
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

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
    IF (PRESENT(atm2)) fmt = fmt+1
    IF (PRESENT(atm3)) fmt = fmt+1
    IF (PRESENT(atm4)) fmt = fmt+1
    CALL integer_to_string(fmt-1,sfmt)
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

!******************************************************************************
!!****** force_fields_util/force_field_qeff_output [1.0] *
!!
!!   NAME
!!     force_field_qeff_output
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Compute the total qeff charges for each molecule kind and total system
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************
  SUBROUTINE force_field_qeff_output (particle_set,atomic_kind_set,&
       molecule_kind_set,molecule_set,mm_section,error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_qeff_output', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atmname, molname
    INTEGER                                  :: handle, iatom, imol, iw
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: qeff, qeff_mol, qeff_sum
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    CALL timeset(routineN,'I','',handle)
    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)

    qeff     = 0.0_dp
    qeff_mol = 0.0_dp
    qeff_sum = 0.0_dp

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Sum of qeff for each molecule_kind
    !-----------------------------------------------------------------------------
    DO imol=1,SIZE(molecule_kind_set)
       qeff_mol=0.0_dp
       molecule_kind => molecule_kind_set(imol)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            name=molname,atom_list=atom_list)
       DO iatom=1,SIZE(atom_list)
          atomic_kind => atom_list(iatom)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=atmname,qeff=qeff)
          IF(iw>0) WRITE(iw,*) "      atom ",iatom," ",TRIM(atmname)," charge = ",qeff
          qeff_mol = qeff_mol + qeff
       END DO
       IF(iw>0) WRITE(iw,*) "    Mol Kind ",TRIM(molname)," charge = ",qeff_mol
    END DO

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Sum of qeff for particle_set
    !-----------------------------------------------------------------------------
    DO iatom=1,SIZE(particle_set)
       atomic_kind => particle_set(iatom)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            name=atmname,qeff=qeff)
       IF(iw>0) WRITE(iw,*) "      atom ",iatom," ",TRIM(atmname),&
            " charge = ",qeff
       qeff_sum = qeff_sum + qeff
    END DO
    IF(iw>0) WRITE(iw,*) "    Total system charge = ",qeff_sum

    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%FF_INFO",error=error)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE force_field_qeff_output

!******************************************************************************
!!****** force_fields_util/clean_intra_forc_kind [1.0] *
!!
!!   NAME
!!     clean_intra_forc_kind
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Removes UNSET force field types
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************
  SUBROUTINE clean_intra_force_kind (molecule_kind_set,mm_section,error)

    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'clean_intra_force_kind', &
      routineP = moduleN//':'//routineN

    INTEGER :: atm2_a, atm2_b, atm2_c, atm_a, atm_b, atm_c, atm_d, central, &
      counter, handle, i, ibend, ibond, icolv, ig3x3, ig4x6, iimpr, ikind, &
      itorsion, iub, iw, j, k, nbend, nbond, newkind, ng3x3, ng4x6, nimpr, &
      ntorsion, nub, stat
    INTEGER, POINTER                         :: bad1(:), bad2(:)
    LOGICAL                                  :: failure, unsetme, valid_kind
    TYPE(bend_kind_type), DIMENSION(:), &
      POINTER                                :: bend_kind_set, &
                                                new_bend_kind_set
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list, new_bend_list
    TYPE(bond_kind_type), DIMENSION(:), &
      POINTER                                :: bond_kind_set, &
                                                new_bond_kind_set
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list, new_bond_list
    TYPE(colvar_constraint_type), &
      DIMENSION(:), POINTER                  :: colv_list
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(g3x3_constraint_type), &
      DIMENSION(:), POINTER                  :: g3x3_list
    TYPE(g4x6_constraint_type), &
      DIMENSION(:), POINTER                  :: g4x6_list
    TYPE(impr_kind_type), DIMENSION(:), &
      POINTER                                :: impr_kind_set, &
                                                new_impr_kind_set
    TYPE(impr_type), DIMENSION(:), POINTER   :: impr_list, new_impr_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(torsion_kind_type), DIMENSION(:), &
      POINTER                                :: new_torsion_kind_set, &
                                                torsion_kind_set
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: new_torsion_list, torsion_list
    TYPE(ub_kind_type), DIMENSION(:), &
      POINTER                                :: new_ub_kind_set, ub_kind_set
    TYPE(ub_type), DIMENSION(:), POINTER     :: new_ub_list, ub_list

    CALL timeset(routineN,'I','',handle)
    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Lets Tag the unwanted bonds due to the use of distance constraint
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            colv_list=colv_list,&
            nbond=nbond,&
            bond_kind_set=bond_kind_set,&
            bond_list=bond_list)
       IF (ASSOCIATED(colv_list)) THEN
          DO icolv=1,SIZE(colv_list)
             IF (colv_list(icolv)%type_id == dist_colvar_id) THEN
                atm_a = colv_list(icolv)%i_atoms(1)
                atm_b = colv_list(icolv)%i_atoms(2)
                DO ibond=1,nbond
                   unsetme = .FALSE.
                   atm2_a = bond_list(ibond)%a
                   atm2_b = bond_list(ibond)%b
                   IF(atm2_a==atm_a .AND. atm2_b==atm_b) unsetme=.TRUE.
                   IF(atm2_a==atm_b .AND. atm2_b==atm_a) unsetme=.TRUE.
                   IF(unsetme) bond_list(ibond)%name = "UNSET"
                END DO
             END IF
          END DO
       END IF
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Lets Tag the unwanted bends due to the use of distance constraint
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            colv_list=colv_list,&
            nbend=nbend,&
            bend_kind_set=bend_kind_set,&
            bend_list=bend_list)
       IF (ASSOCIATED(colv_list)) THEN
          DO ibend=1,nbend
             counter=0
             unsetme = .FALSE.
             atm_a = bend_list(ibend)%a
             atm_b = bend_list(ibend)%b
             atm_c = bend_list(ibend)%c
             central = 0
             DO icolv=1,SIZE(colv_list)
                IF (colv_list(icolv)%type_id == dist_colvar_id) THEN
                   atm2_a = colv_list(icolv)%i_atoms(1)
                   atm2_b = colv_list(icolv)%i_atoms(2)
                   ! Check to see if we have a central bend atom in dist list
                   IF ( central == 0 ) THEN
                      IF ( atm2_a == atm_b ) central = atm2_a
                      IF ( atm2_b == atm_b ) central = atm2_b
                   ENDIF
                   ! If so, check the other two atoms for a match
                   IF ( atm2_a == central ) THEN
                      IF(atm2_b==atm_a.OR.atm2_b==atm_c) counter=counter+1
                      ! If so, check the other two atoms for a match
                   ELSEIF ( atm2_b == central ) THEN
                      IF(atm2_a==atm_a.OR.atm2_a==atm_c) counter=counter+1
                      ! If not, check the other two atoms for a match
                   ELSE
                      IF(atm2_a==atm_a.OR.atm2_a==atm_c) counter=counter+1
                      IF(atm2_b==atm_a.OR.atm2_b==atm_c) counter=counter+1
                   ENDIF
                END IF
             END DO
             ! must have central atom to be involved in bend
             IF ( central /= 0 ) THEN
                IF(counter>=4) unsetme=.TRUE.
                IF(unsetme) bend_list(ibend)%name = "UNSET"
             ENDIF
          END DO
       END IF
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Lets Tag the unwanted bonds due to the use of 3x3
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            ng3x3=ng3x3,&
            g3x3_list=g3x3_list,&
            nbond=nbond,&
            bond_kind_set=bond_kind_set,&
            bond_list=bond_list)
       DO ig3x3=1,ng3x3
          atm_a = g3x3_list(ig3x3)%a
          atm_b = g3x3_list(ig3x3)%b
          atm_c = g3x3_list(ig3x3)%c
          DO ibond=1,nbond
             unsetme = .FALSE.
             atm2_a = bond_list(ibond)%a
             atm2_b = bond_list(ibond)%b
             IF(atm2_a==atm_a .AND. atm2_b==atm_b) unsetme=.TRUE.
             IF(atm2_a==atm_a .AND. atm2_b==atm_c) unsetme=.TRUE.
             IF(atm2_a==atm_c .AND. atm2_b==atm_c) unsetme=.TRUE.
             IF(unsetme) bond_list(ibond)%name = "UNSET"
          END DO
       END DO
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Lets Tag the unwanted bends due to the use of 3x3
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            ng3x3=ng3x3,&
            g3x3_list=g3x3_list,&
            nbend=nbend,&
            bend_kind_set=bend_kind_set,&
            bend_list=bend_list)
       DO ig3x3=1,ng3x3
          atm_a = g3x3_list(ig3x3)%a
          atm_b = g3x3_list(ig3x3)%b
          atm_c = g3x3_list(ig3x3)%c
          DO ibend=1,nbend
             unsetme = .FALSE.
             atm2_a = bend_list(ibend)%a
             atm2_b = bend_list(ibend)%b
             atm2_c = bend_list(ibend)%c
             IF(atm2_a==atm_a .AND. atm2_b==atm_b .AND. atm2_c==atm_c) unsetme=.TRUE.
             IF(atm2_a==atm_a .AND. atm2_b==atm_c .AND. atm2_c==atm_b) unsetme=.TRUE.
             IF(atm2_a==atm_b .AND. atm2_b==atm_a .AND. atm2_c==atm_c) unsetme=.TRUE.
             IF(atm2_a==atm_b .AND. atm2_b==atm_c .AND. atm2_c==atm_a) unsetme=.TRUE.
             IF(unsetme) bend_list(ibend)%name = "UNSET"
          END DO
       END DO
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Lets Tag the unwanted bonds due to the use of 4x6
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            ng4x6=ng4x6,&
            g4x6_list=g4x6_list,&
            nbond=nbond,&
            bond_kind_set=bond_kind_set,&
            bond_list=bond_list)
       DO ig4x6=1,ng4x6
          atm_a = g4x6_list(ig4x6)%a
          atm_b = g4x6_list(ig4x6)%b
          atm_c = g4x6_list(ig4x6)%c
          atm_d = g4x6_list(ig4x6)%d
          DO ibond=1,nbond
             unsetme = .FALSE.
             atm2_a = bond_list(ibond)%a
             atm2_b = bond_list(ibond)%b
             IF(atm2_a==atm_a .AND. atm2_b==atm_b) unsetme=.TRUE.
             IF(atm2_a==atm_a .AND. atm2_b==atm_c) unsetme=.TRUE.
             IF(atm2_a==atm_a .AND. atm2_b==atm_d) unsetme=.TRUE.
             IF(unsetme) bond_list(ibond)%name = "UNSET"
          END DO
       END DO
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Lets Tag the unwanted bends due to the use of 4x6
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            ng4x6=ng4x6,&
            g4x6_list=g4x6_list,&
            nbend=nbend,&
            bend_kind_set=bend_kind_set,&
            bend_list=bend_list)
       DO ig4x6=1,ng4x6
          atm_a = g4x6_list(ig4x6)%a
          atm_b = g4x6_list(ig4x6)%b
          atm_c = g4x6_list(ig4x6)%c
          atm_d = g4x6_list(ig4x6)%d
          DO ibend=1,nbend
             unsetme = .FALSE.
             atm2_a = bend_list(ibend)%a
             atm2_b = bend_list(ibend)%b
             atm2_c = bend_list(ibend)%c
             IF(atm2_a==atm_b .AND. atm2_b==atm_a .AND. atm2_c==atm_c) unsetme=.TRUE.
             IF(atm2_a==atm_b .AND. atm2_b==atm_a .AND. atm2_c==atm_d) unsetme=.TRUE.
             IF(atm2_a==atm_c .AND. atm2_b==atm_a .AND. atm2_c==atm_d) unsetme=.TRUE.
             IF(unsetme) bend_list(ibend)%name = "UNSET"
          END DO
       END DO
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 7. Count the number of UNSET bond kinds there are
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            nbond=nbond,&
            bond_kind_set=bond_kind_set,&
            bond_list=bond_list)
       IF(nbond>0) THEN
          IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") Old BOND Count: ",&
               SIZE(bond_list),SIZE(bond_kind_set)
          IF (iw>0) WRITE(iw,'(2I6)')(bond_list(ibond)%a,bond_list(ibond)%b,ibond=1,SIZE(bond_list))
          NULLIFY(bad1,bad2)
          ALLOCATE(bad1(SIZE(bond_kind_set)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          bad1(:)=0
          DO ibond=1,SIZE(bond_kind_set)
             unsetme=.FALSE.
             IF(bond_kind_set(ibond)%name=="UNSET") unsetme = .TRUE.
             valid_kind=.FALSE.
             DO i=1,SIZE(bond_list)
                IF(bond_list(i)%name/="UNSET".AND.&
                     bond_list(i)%bond_kind%kind_number==ibond) THEN
                   valid_kind=.TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.valid_kind) unsetme = .TRUE.
             IF(unsetme) bad1(ibond) = 1
          END DO
          IF(SUM(bad1)/=0) THEN
             counter = SIZE(bond_kind_set)-SUM(bad1)
             CALL allocate_bond_kind_set(new_bond_kind_set,counter,error)
             counter=0
             DO ibond=1,SIZE(bond_kind_set)
                IF(bad1(ibond)==0) THEN
                   counter=counter+1
                   new_bond_kind_set(counter)= bond_kind_set(ibond)
                END IF
             END DO
             counter=0
             ALLOCATE(bad2(SIZE(bond_list)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             bad2(:)=0
             DO ibond=1,SIZE(bond_list)
                unsetme = .FALSE.
                IF(bond_list(ibond)%bond_kind%name=="UNSET") unsetme = .TRUE.
                IF(bond_list(ibond)%name=="UNSET")           unsetme = .TRUE.
                IF(unsetme) bad2(ibond) = 1
             END DO
             IF(SUM(bad2)/=0) THEN
                counter = SIZE(bond_list)-SUM(bad2)
                ALLOCATE(new_bond_list(counter),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                counter=0
                DO ibond=1,SIZE(bond_list)
                   IF(bad2(ibond)==0) THEN
                      counter=counter+1
                      new_bond_list(counter) = bond_list(ibond)
                      newkind = bond_list(ibond)%bond_kind%kind_number
                      newkind = newkind - SUM(bad1(1:newkind))
                      new_bond_list(counter)%bond_kind => new_bond_kind_set(newkind)
                   END IF
                END DO
                CALL set_molecule_kind(molecule_kind=molecule_kind,&
                     nbond=SIZE(new_bond_list),&
                     bond_kind_set=new_bond_kind_set,&
                     bond_list=new_bond_list)
                DO ibond=1,SIZE(new_bond_kind_set)
                   new_bond_kind_set(ibond)%kind_number=ibond
                END DO
             END IF
             DEALLOCATE(bad2,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CALL deallocate_bond_kind_set(bond_kind_set,error)
             DEALLOCATE(bond_list,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") New BOND Count: ",&
                  SIZE(new_bond_list),SIZE(new_bond_kind_set)
             IF (iw>0) WRITE(iw,'(2I6)')(new_bond_list(ibond)%a,new_bond_list(ibond)%b,&
                  ibond=1,SIZE(new_bond_list))
          END IF
          DEALLOCATE(bad1,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 8. Count the number of UNSET bend kinds there are
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            nbend=nbend,&
            bend_kind_set=bend_kind_set,&
            bend_list=bend_list)
       IF(nbend>0) THEN
          IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") Old BEND Count: ",&
               SIZE(bend_list),SIZE(bend_kind_set)
          IF (iw>0) WRITE(iw,'(3I6)')(bend_list(ibend)%a,bend_list(ibend)%b,&
               bend_list(ibend)%c,ibend=1,SIZE(bend_list))
          NULLIFY(bad1,bad2)
          ALLOCATE(bad1(SIZE(bend_kind_set)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          bad1(:)=0
          DO ibend=1,SIZE(bend_kind_set)
             unsetme=.FALSE.
             IF(bend_kind_set(ibend)%name=="UNSET") unsetme = .TRUE.
             valid_kind=.FALSE.
             DO i=1,SIZE(bend_list)
                IF(bend_list(i)%name/="UNSET".AND.&
                     bend_list(i)%bend_kind%kind_number==ibend) THEN
                   valid_kind=.TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.valid_kind) unsetme = .TRUE.
             IF(unsetme) bad1(ibend) = 1
          END DO
          IF(SUM(bad1)/=0) THEN
             counter = SIZE(bend_kind_set)-SUM(bad1)
             CALL allocate_bend_kind_set(new_bend_kind_set,counter,error)
             counter=0
             DO ibend=1,SIZE(bend_kind_set)
                IF(bad1(ibend)==0) THEN
                   counter=counter+1
                   new_bend_kind_set(counter)= bend_kind_set(ibend)
                END IF
             END DO
             counter=0
             ALLOCATE(bad2(SIZE(bend_list)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             bad2(:)=0
             DO ibend=1,SIZE(bend_list)
                unsetme = .FALSE.
                IF(bend_list(ibend)%bend_kind%name=="UNSET") unsetme = .TRUE.
                IF(bend_list(ibend)%name=="UNSET")           unsetme = .TRUE.
                IF(unsetme) bad2(ibend) = 1
             END DO
             IF(SUM(bad2)/=0) THEN
                counter = SIZE(bend_list)-SUM(bad2)
                ALLOCATE(new_bend_list(counter),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                counter=0
                DO ibend=1,SIZE(bend_list)
                   IF(bad2(ibend)==0) THEN
                      counter=counter+1
                      new_bend_list(counter) = bend_list(ibend)
                      newkind = bend_list(ibend)%bend_kind%kind_number
                      newkind = newkind - SUM(bad1(1:newkind))
                      new_bend_list(counter)%bend_kind => new_bend_kind_set(newkind)
                   END IF
                END DO
                CALL set_molecule_kind(molecule_kind=molecule_kind,&
                     nbend=SIZE(new_bend_list),&
                     bend_kind_set=new_bend_kind_set,&
                     bend_list=new_bend_list)
                DO ibend=1,SIZE(new_bend_kind_set)
                   new_bend_kind_set(ibend)%kind_number=ibend
                END DO
             END IF
             DEALLOCATE(bad2,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CALL deallocate_bend_kind_set(bend_kind_set,error)
             DEALLOCATE(bend_list,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") New BEND Count: ",&
                  SIZE(new_bend_list),SIZE(new_bend_kind_set)
             IF (iw>0) WRITE(iw,'(3I6)')(new_bend_list(ibend)%a,new_bend_list(ibend)%b,&
                  new_bend_list(ibend)%c,ibend=1,SIZE(new_bend_list))
          END IF
          DEALLOCATE(bad1,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 9. Count the number of UNSET Urey-Bradley kinds there are
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            nub=nub,&
            ub_kind_set=ub_kind_set,&
            ub_list=ub_list)
       IF(nub>0) THEN
          IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") Old UB Count: ",&
               SIZE(ub_list),SIZE(ub_kind_set)
          IF (iw>0) WRITE(iw,'(3I6)')(ub_list(iub)%a,ub_list(iub)%b,&
               ub_list(iub)%c,iub=1,SIZE(ub_list))
          NULLIFY(bad1,bad2)
          ALLOCATE(bad1(SIZE(ub_kind_set)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          bad1(:)=0
          DO iub=1,SIZE(ub_kind_set)
             unsetme=.FALSE.
             IF(ub_kind_set(iub)%name=="UNSET") unsetme = .TRUE.
             valid_kind=.FALSE.
             DO i=1,SIZE(ub_list)
                IF(ub_list(i)%name/="UNSET".AND.&
                     ub_list(i)%ub_kind%kind_number==iub) THEN
                   valid_kind=.TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.valid_kind) unsetme = .TRUE.
             IF(unsetme) bad1(iub) = 1
          END DO
          IF(SUM(bad1)/=0) THEN
             counter = SIZE(ub_kind_set)-SUM(bad1)
             CALL allocate_ub_kind_set(new_ub_kind_set,counter,error)
             counter=0
             DO iub=1,SIZE(ub_kind_set)
                IF(bad1(iub)==0) THEN
                   counter=counter+1
                   new_ub_kind_set(counter)= ub_kind_set(iub)
                END IF
             END DO
             counter=0
             ALLOCATE(bad2(SIZE(ub_list)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             bad2(:)=0
             DO iub=1,SIZE(ub_list)
                unsetme = .FALSE.
                IF(ub_list(iub)%ub_kind%name=="UNSET") unsetme = .TRUE.
                IF(ub_list(iub)%name=="UNSET")           unsetme = .TRUE.
                IF(unsetme) bad2(iub) = 1
             END DO
             IF(SUM(bad2)/=0) THEN
                counter = SIZE(ub_list)-SUM(bad2)
                ALLOCATE(new_ub_list(counter),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                counter=0
                DO iub=1,SIZE(ub_list)
                   IF(bad2(iub)==0) THEN
                      counter=counter+1
                      new_ub_list(counter) = ub_list(iub)
                      newkind = ub_list(iub)%ub_kind%kind_number
                      newkind = newkind - SUM(bad1(1:newkind))
                      new_ub_list(counter)%ub_kind => new_ub_kind_set(newkind)
                   END IF
                END DO
                CALL set_molecule_kind(molecule_kind=molecule_kind,&
                     nub=SIZE(new_ub_list),&
                     ub_kind_set=new_ub_kind_set,&
                     ub_list=new_ub_list)
                DO iub=1,SIZE(new_ub_kind_set)
                   new_ub_kind_set(iub)%kind_number=iub
                END DO
             END IF
             DEALLOCATE(bad2,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CALL ub_kind_dealloc_ref(ub_kind_set,error=error)
             DEALLOCATE(ub_list,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") New UB Count: ",&
                  SIZE(new_ub_list),SIZE(new_ub_kind_set)
             IF (iw>0) WRITE(iw,'(3I6)')(new_ub_list(iub)%a,new_ub_list(iub)%b,&
                  new_ub_list(iub)%c,iub=1,SIZE(new_ub_list))
          END IF
          DEALLOCATE(bad1,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 10. Count the number of UNSET torsion kinds there are
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            ntorsion=ntorsion,&
            torsion_kind_set=torsion_kind_set,&
            torsion_list=torsion_list)
       IF(ntorsion>0) THEN
          IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") Old TORSION Count: ",&
               SIZE(torsion_list),SIZE(torsion_kind_set)
          IF (iw>0) WRITE(iw,'(4I6)')(torsion_list(itorsion)%a,torsion_list(itorsion)%b,&
               torsion_list(itorsion)%c,torsion_list(itorsion)%d,itorsion=1,SIZE(torsion_list))
          NULLIFY(bad1,bad2)
          ALLOCATE(bad1(SIZE(torsion_kind_set)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          bad1(:)=0
          DO itorsion=1,SIZE(torsion_kind_set)
             unsetme=.FALSE.
             IF(torsion_kind_set(itorsion)%name=="UNSET") unsetme = .TRUE.
             valid_kind=.FALSE.
             DO i=1,SIZE(torsion_list)
                IF(torsion_list(i)%name/="UNSET".AND.&
                     torsion_list(i)%torsion_kind%kind_number==itorsion) THEN
                   valid_kind=.TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.valid_kind) unsetme = .TRUE.
             IF(unsetme) bad1(itorsion) = 1
          END DO
          IF(SUM(bad1)/=0) THEN
             counter = SIZE(torsion_kind_set)-SUM(bad1)
             CALL allocate_torsion_kind_set(new_torsion_kind_set,counter,error)
             counter=0
             DO itorsion=1,SIZE(torsion_kind_set)
                IF(bad1(itorsion)==0) THEN
                   counter=counter+1
                   new_torsion_kind_set(counter)= torsion_kind_set(itorsion)
                   i = SIZE(torsion_kind_set(itorsion)%m)
                   j = SIZE(torsion_kind_set(itorsion)%k)
                   k = SIZE(torsion_kind_set(itorsion)%cosphi0)
                   ALLOCATE(new_torsion_kind_set(counter)%m(i),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   ALLOCATE(new_torsion_kind_set(counter)%k(i),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   ALLOCATE(new_torsion_kind_set(counter)%cosphi0(i),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   new_torsion_kind_set(counter)%m= torsion_kind_set(itorsion)%m
                   new_torsion_kind_set(counter)%k= torsion_kind_set(itorsion)%k
                   new_torsion_kind_set(counter)%cosphi0= torsion_kind_set(itorsion)%cosphi0
                END IF
             END DO
             counter=0
             ALLOCATE(bad2(SIZE(torsion_list)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             bad2(:)=0
             DO itorsion=1,SIZE(torsion_list)
                unsetme = .FALSE.
                IF(torsion_list(itorsion)%torsion_kind%name=="UNSET") unsetme = .TRUE.
                IF(torsion_list(itorsion)%name=="UNSET")           unsetme = .TRUE.
                IF(unsetme) bad2(itorsion) = 1
             END DO
             IF(SUM(bad2)/=0) THEN
                counter = SIZE(torsion_list)-SUM(bad2)
                ALLOCATE(new_torsion_list(counter),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                counter=0
                DO itorsion=1,SIZE(torsion_list)
                   IF(bad2(itorsion)==0) THEN
                      counter=counter+1
                      new_torsion_list(counter) = torsion_list(itorsion)
                      newkind = torsion_list(itorsion)%torsion_kind%kind_number
                      newkind = newkind - SUM(bad1(1:newkind))
                      new_torsion_list(counter)%torsion_kind => new_torsion_kind_set(newkind)
                   END IF
                END DO
                CALL set_molecule_kind(molecule_kind=molecule_kind,&
                     ntorsion=SIZE(new_torsion_list),&
                     torsion_kind_set=new_torsion_kind_set,&
                     torsion_list=new_torsion_list)
                DO itorsion=1,SIZE(new_torsion_kind_set)
                   new_torsion_kind_set(itorsion)%kind_number=itorsion
                END DO
             END IF
             DEALLOCATE(bad2,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO itorsion=1,SIZE(torsion_kind_set)
                CALL torsion_kind_dealloc_ref(torsion_kind_set(itorsion),error=error)
             END DO
             DEALLOCATE(torsion_kind_set,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(torsion_list,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") New TORSION Count: ",&
                  SIZE(new_torsion_list),SIZE(new_torsion_kind_set)
             IF (iw>0) WRITE(iw,'(4I6)')(new_torsion_list(itorsion)%a,new_torsion_list(itorsion)%b,&
                  new_torsion_list(itorsion)%c,new_torsion_list(itorsion)%d,itorsion=1,&
                  SIZE(new_torsion_list))
          END IF
          DEALLOCATE(bad1,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 11. Count the number of UNSET improper kinds there are
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            nimpr=nimpr,&
            impr_kind_set=impr_kind_set,&
            impr_list=impr_list)
       IF(nimpr>0) THEN
          IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") Old IMPROPER Count: ",&
               SIZE(impr_list),SIZE(impr_kind_set)
          NULLIFY(bad1,bad2)
          ALLOCATE(bad1(SIZE(impr_kind_set)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          bad1(:)=0
          DO iimpr=1,SIZE(impr_kind_set)
             unsetme=.FALSE.
             IF(impr_kind_set(iimpr)%name=="UNSET") unsetme = .TRUE.
             valid_kind=.FALSE.
             DO i=1,SIZE(impr_list)
                IF(impr_list(i)%name/="UNSET".AND.&
                     impr_list(i)%impr_kind%kind_number==iimpr) THEN
                   valid_kind=.TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.valid_kind) unsetme = .TRUE.
             IF(unsetme) bad1(iimpr) = 1
          END DO
          IF(SUM(bad1)/=0) THEN
             counter = SIZE(impr_kind_set)-SUM(bad1)
             CALL allocate_impr_kind_set(new_impr_kind_set,counter,error)
             counter=0
             DO iimpr=1,SIZE(impr_kind_set)
                IF(bad1(iimpr)==0) THEN
                   counter=counter+1
                   new_impr_kind_set(counter)= impr_kind_set(iimpr)
                END IF
             END DO
             counter=0
             ALLOCATE(bad2(SIZE(impr_list)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             bad2(:)=0
             DO iimpr=1,SIZE(impr_list)
                unsetme = .FALSE.
                IF(impr_list(iimpr)%impr_kind%name=="UNSET") unsetme = .TRUE.
                IF(impr_list(iimpr)%name=="UNSET")           unsetme = .TRUE.
                IF(unsetme) bad2(iimpr) = 1
             END DO
             IF(SUM(bad2)/=0) THEN
                counter = SIZE(impr_list)-SUM(bad2)
                ALLOCATE(new_impr_list(counter),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                counter=0
                DO iimpr=1,SIZE(impr_list)
                   IF(bad2(iimpr)==0) THEN
                      counter=counter+1
                      new_impr_list(counter) = impr_list(iimpr)
                      newkind = impr_list(iimpr)%impr_kind%kind_number
                      newkind = newkind - SUM(bad1(1:newkind))
                      new_impr_list(counter)%impr_kind => new_impr_kind_set(newkind)
                   END IF
                END DO
                CALL set_molecule_kind(molecule_kind=molecule_kind,&
                     nimpr=SIZE(new_impr_list),&
                     impr_kind_set=new_impr_kind_set,&
                     impr_list=new_impr_list)
                DO iimpr=1,SIZE(new_impr_kind_set)
                   new_impr_kind_set(iimpr)%kind_number=iimpr
                END DO
             END IF
             DEALLOCATE(bad2,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO iimpr=1,SIZE(impr_kind_set)
                CALL impr_kind_dealloc_ref(impr_kind_set(iimpr),error=error)
             END DO
             DEALLOCATE(impr_kind_set,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(impr_list,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") New IMPROPER Count: ",&
                  SIZE(new_impr_list),SIZE(new_impr_kind_set)
          END IF
          DEALLOCATE(bad1,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! 12. Count the number of UNSET NONBOND14 kinds there are
    !-                NEED TO REMOVE EXTRAS HERE   - IKUO
    !---------------------------------------------------------------------------

    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%FF_INFO",error=error)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE clean_intra_force_kind

END MODULE force_fields_util
!******************************************************************************
