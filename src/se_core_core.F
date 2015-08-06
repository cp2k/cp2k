!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Split and build its own idependent core_core SE interaction module
!> \author Teodoro Laino [tlaino] - 05.2009
!> \par History
!>      Teodoro Laino (05.2009) [tlaino] - create
! *****************************************************************************
MODULE se_core_core
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE atprop_types,                    ONLY: atprop_array_init,&
                                             atprop_type
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             semi_empirical_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                             ewald_pw_type
  USE input_constants,                 ONLY: do_method_am1,&
                                             do_method_mndo,&
                                             do_method_mndod,&
                                             do_method_pdg,&
                                             do_method_pm3,&
                                             do_method_pm6,&
                                             do_method_pnnl,&
                                             do_method_rm1
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE semi_empirical_int_arrays,       ONLY: rij_threshold
  USE semi_empirical_integrals,        ONLY: corecore,&
                                             dcorecore
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             se_int_control_type,&
                                             se_taper_type,&
                                             semi_empirical_p_type,&
                                             semi_empirical_type,&
                                             setup_se_int_control_type
  USE semi_empirical_utils,            ONLY: finalize_se_taper,&
                                             get_se_type,&
                                             initialize_se_taper
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'se_core_core'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module       = .FALSE.

  PUBLIC :: se_core_core_interaction

CONTAINS

! *****************************************************************************
!> \brief Evaluates the core-core interactions for NDDO methods
!> \param qs_env ...
!> \param para_env ...
!> \param calculate_forces ...
!> \param error ...
!> \date 04.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE se_core_core_interaction(qs_env, para_env, calculate_forces, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'se_core_core_interaction', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, handle, iab, &
                                                iatom, ikind, itype, jatom, &
                                                jkind, natom, nkind, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: anag, atener, defined, &
                                                failure, use_virial
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: se_defined
    REAL(KIND=dp)                            :: delta, dr1, dr3inv(3), enuc, &
                                                enucij, enuclear, r2inv, &
                                                r3inv, rinv
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atprop_type), POINTER               :: atprop
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_se
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(se_int_control_type)                :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_p_type), &
      DIMENSION(:), POINTER                  :: se_kind_param
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b
    TYPE(virial_type), POINTER               :: virial

    failure  = .FALSE.
    enuclear = 0.0_dp
    NULLIFY(dft_control,cell,force,particle_set,se_control,se_taper,atomic_kind_set,&
            virial, atprop)

    CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, se_taper=se_taper,&
         virial=virial,atprop=atprop, energy=energy, error=error)

    CALL initialize_se_taper(se_taper,coulomb=.TRUE.,error=error)
    ! Parameters
    se_control => dft_control%qs_control%se_control
    anag       =  se_control%analytical_gradients
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
    CALL setup_se_int_control_type(se_int_control,do_ewald_r3=se_control%do_ewald_r3,&
         do_ewald_gks=se_control%do_ewald_gks, integral_screening=se_control%integral_screening,&
         shortrange=(se_control%do_ewald.OR.se_control%do_ewald_gks),&
         max_multipole=se_control%max_multipole, pc_coulomb_int=.FALSE.)

    ! atomic energy decomposition
    atener = atprop%energy
    IF (atener) THEN
      CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,error=error)
      natom = SIZE (particle_set)
      CALL atprop_array_init(atprop%atecc,natom,error)
    END IF

    ! Retrieve some information if GKS ewald scheme is used
    IF(se_control%do_ewald_gks) THEN
       CALL get_qs_env(qs_env=qs_env,ewald_env=ewald_env,ewald_pw=ewald_pw,error=error)
       CALL ewald_env_get (ewald_env, alpha=se_int_control%ewald_gks%alpha, error=error)
       CALL ewald_pw_get (ewald_pw, pw_big_pool=se_int_control%ewald_gks%pw_pool, &
            dg=se_int_control%ewald_gks%dg)
       ! Virial not implemented
       CPPrecondition(.NOT.use_virial,cp_failure_level,routineP,error,failure)
    END IF

    CALL get_qs_env(qs_env=qs_env,sab_se=sab_se,atomic_kind_set=atomic_kind_set,&
         qs_kind_set=qs_kind_set, error=error)

    nkind = SIZE(atomic_kind_set)
    ! Possibly compute forces
    IF(calculate_forces) THEN
       CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,force=force,error=error)
       natom = SIZE (particle_set)
       ALLOCATE (atom_of_kind(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       delta = se_control%delta
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)
    END IF

    itype    = get_se_type(dft_control%qs_control%method_id)

    ALLOCATE (se_kind_param(nkind),se_defined(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
       CALL get_qs_kind(qs_kind_set(ikind), se_parameter=se_kind_a, error=error)
       se_kind_param(ikind)%se_param => se_kind_a
       CALL get_se_param(se_kind_a,defined=defined)
       se_defined(ikind)=defined
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_se)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,r=rij)
       IF (.NOT.se_defined(ikind)) CYCLE
       IF (.NOT.se_defined(jkind)) CYCLE
       se_kind_a => se_kind_param(ikind)%se_param
       se_kind_b => se_kind_param(jkind)%se_param
       iab = ikind + nkind*(jkind - 1)
       dr1 = DOT_PRODUCT(rij,rij)
       enucij = 0._dp
       IF ( dr1 > rij_threshold ) THEN
          SELECT CASE (dft_control%qs_control%method_id)
          CASE (do_method_mndo,do_method_am1,do_method_pm3,do_method_pm6,do_method_pdg,&
               do_method_rm1,do_method_mndod, do_method_pnnl)

             ! Core-Core energy term
             CALL corecore (se_kind_a,se_kind_b,rij,enuc=enuc,itype=itype,anag=anag,&
                  se_int_control=se_int_control, se_taper=se_taper, error=error)
             enucij = enucij + enuc
             ! Residual integral (1/R^3) correction
             IF (se_int_control%do_ewald_r3) THEN
                r2inv = 1.0_dp/dr1
                rinv  = SQRT(r2inv)
                r3inv = rinv**3
                ! Core-Core term
                enucij = enucij + se_kind_a%expns3_int(jkind)%expns3%core_core*r3inv
             END IF

             ! Core-Core Derivatives
             IF(calculate_forces) THEN
                atom_a = atom_of_kind(iatom)
                atom_b = atom_of_kind(jatom)

                CALL dcorecore (se_kind_a,se_kind_b,rij,denuc=force_ab,itype=itype,delta=delta,&
                     anag=anag,se_int_control=se_int_control,se_taper=se_taper,error=error)

                ! Residual integral (1/R^3) correction
                IF (se_int_control%do_ewald_r3) THEN
                   dr3inv   = -3.0_dp*rij*r3inv*r2inv
                   ! Derivatives of core-core terms
                   force_ab = force_ab + se_kind_a%expns3_int(jkind)%expns3%core_core*dr3inv
                END IF
                IF (use_virial) THEN
                   CALL virial_pair_force ( virial%pv_virial, -1.0_dp, force_ab, rij, error)
                END IF

                ! Sum up force components
                force(ikind)%all_potential(1,atom_a) = force(ikind)%all_potential(1,atom_a) - force_ab(1)
                force(jkind)%all_potential(1,atom_b) = force(jkind)%all_potential(1,atom_b) + force_ab(1)

                force(ikind)%all_potential(2,atom_a) = force(ikind)%all_potential(2,atom_a) - force_ab(2)
                force(jkind)%all_potential(2,atom_b) = force(jkind)%all_potential(2,atom_b) + force_ab(2)

                force(ikind)%all_potential(3,atom_a) = force(ikind)%all_potential(3,atom_a) - force_ab(3)
                force(jkind)%all_potential(3,atom_b) = force(jkind)%all_potential(3,atom_b) + force_ab(3)
             END IF
          CASE DEFAULT
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
       ELSE
          IF ( se_int_control%do_ewald_gks ) THEN
             ! Core-Core energy term (self term in periodic systems)
             CALL corecore (se_kind_a,se_kind_b,rij,enuc=enuc,itype=itype,anag=anag,&
                  se_int_control=se_int_control, se_taper=se_taper, error=error)
             enucij = enucij + 0.5_dp*enuc
          END IF
       END IF
       IF (atener) THEN
         atprop%atecc(iatom) = atprop%atecc(iatom) + 0.5_dp*enucij
         atprop%atecc(jatom) = atprop%atecc(jatom) + 0.5_dp*enucij
       END IF
       enuclear = enuclear + enucij
    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    DEALLOCATE (se_kind_param,se_defined,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (calculate_forces) THEN
       DEALLOCATE(atom_of_kind,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL mp_sum(enuclear,para_env%group)
    energy%core_overlap  = enuclear
    energy%core_overlap0 = enuclear

    CALL finalize_se_taper(se_taper,error=error)
    CALL timestop(handle)
  END SUBROUTINE se_core_core_interaction

END MODULE se_core_core

