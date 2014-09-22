!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief SCPTB short range repulsive core interaction
!> \author Juerg Hutter - 12.2011
!> \par History
!>      JHU (12.2011) - create
! *****************************************************************************
MODULE scptb_core_interactions
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE atprop_types,                    ONLY: atprop_array_init,&
                                             atprop_type
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
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
  USE scptb_types,                     ONLY: get_scptb_parameter,&
                                             scptb_parameter_p_type,&
                                             scptb_parameter_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scptb_core_interactions'

  PUBLIC :: scptb_core_interaction

CONTAINS

! *****************************************************************************
!> \brief Evaluates the repulsive core interactions for SCPTB methods
!> \param qs_env ...
!> \param calculate_forces ...
!> \param error ...
!> \date 12.2011 [jhu]
!> \author Juerg Hutter [jhu] - University of Zurich
! *****************************************************************************
  SUBROUTINE scptb_core_interaction(qs_env, calculate_forces, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scptb_core_interaction', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, handle, &
                                                iatom, ikind, jatom, jkind, &
                                                natom, nkind, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: atener, defined, failure, &
                                                use_virial
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: scptb_defined
    REAL(KIND=dp)                            :: dr1, enuc, enucij, enuclear
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atprop_type), POINTER               :: atprop
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_core
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(scptb_parameter_p_type), &
      DIMENSION(:), POINTER                  :: scptb_kind_param
    TYPE(scptb_parameter_type), POINTER      :: scptb_kind_a, scptb_kind_b
    TYPE(virial_type), POINTER               :: virial

    failure  = .FALSE.
    enuclear = 0.0_dp
    NULLIFY(atprop)
    CALL timeset(routineN,handle)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CALL get_qs_env(qs_env=qs_env,para_env=para_env,dft_control=dft_control,&
                    cell=cell,virial=virial,energy=energy,atprop=atprop,error=error)

    ! Parameters
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    ! atomic energy decomposition
    atener = atprop%energy
    IF (atener) THEN
      CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,error=error)
      natom = SIZE (particle_set)
      CALL atprop_array_init(atprop%atecc,natom,error)
    END IF

    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env=qs_env,sab_core=sab_core, atomic_kind_set=atomic_kind_set,&
                      qs_kind_set=qs_kind_set,error=error)

       nkind = SIZE(atomic_kind_set)
       ! Possibly compute forces
       IF(calculate_forces) THEN
          CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,force=force,error=error)
          natom = SIZE (particle_set)
          ALLOCATE (atom_of_kind(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL get_atomic_kind_set(atomic_kind_set,atom_of_kind=atom_of_kind)
       END IF

       ALLOCATE (scptb_kind_param(nkind),scptb_defined(nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO ikind=1,nkind
          CALL get_qs_kind(qs_kind_set(ikind),scptb_parameter=scptb_kind_a)
          scptb_kind_param(ikind)%scptb_param => scptb_kind_a
          CALL get_scptb_parameter(scptb_kind_a,defined=defined)
          scptb_defined(ikind)=defined
       END DO
       CALL neighbor_list_iterator_create(nl_iterator,sab_core)
       DO WHILE (neighbor_list_iterate(nl_iterator)==0)
          CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,r=rij)
          IF (.NOT.scptb_defined(ikind)) CYCLE
          IF (.NOT.scptb_defined(jkind)) CYCLE
          scptb_kind_a => scptb_kind_param(ikind)%scptb_param
          scptb_kind_b => scptb_kind_param(jkind)%scptb_param
          dr1 = DOT_PRODUCT(rij,rij)
          enucij = 0._dp
          IF ( dr1 > 0.00001_dp ) THEN
             ! Core-Core energy and derivatives
             IF(calculate_forces) THEN
                CALL corecore (scptb_kind_a,scptb_kind_b,rij,enuc=enuc,denuc=force_ab,error=error)

                atom_a = atom_of_kind(iatom)
                atom_b = atom_of_kind(jatom)

                ! Sum up force components
                force(ikind)%all_potential(1,atom_a) = force(ikind)%all_potential(1,atom_a) - force_ab(1)
                force(jkind)%all_potential(1,atom_b) = force(jkind)%all_potential(1,atom_b) + force_ab(1)

                force(ikind)%all_potential(2,atom_a) = force(ikind)%all_potential(2,atom_a) - force_ab(2)
                force(jkind)%all_potential(2,atom_b) = force(jkind)%all_potential(2,atom_b) + force_ab(2)

                force(ikind)%all_potential(3,atom_a) = force(ikind)%all_potential(3,atom_a) - force_ab(3)
                force(jkind)%all_potential(3,atom_b) = force(jkind)%all_potential(3,atom_b) + force_ab(3)

                IF (use_virial) THEN
                   CALL virial_pair_force ( virial%pv_virial, -1._dp, force_ab, rij, error)
                END IF

             ELSE
                CALL corecore (scptb_kind_a,scptb_kind_b,rij,enuc=enuc,error=error)
             END IF
             enucij = enucij + enuc
          END IF
          IF (atener) THEN
            atprop%atecc(iatom) = atprop%atecc(iatom) + 0.5_dp*enucij
            atprop%atecc(jatom) = atprop%atecc(jatom) + 0.5_dp*enucij
          END IF
          enuclear = enuclear + enucij
       END DO
       CALL neighbor_list_iterator_release(nl_iterator)

       DEALLOCATE (scptb_kind_param,scptb_defined,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       IF (calculate_forces) THEN
          DEALLOCATE(atom_of_kind,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       CALL mp_sum(enuclear,para_env%group)
       energy%repulsive  = enuclear
    END IF

    CALL timestop(handle)

  END SUBROUTINE scptb_core_interaction

! *****************************************************************************
!> \brief ...
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param enuc ...
!> \param denuc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE corecore (sepi,sepj,rij,enuc,denuc,error)
    TYPE(scptb_parameter_type), POINTER      :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), INTENT(OUT)                    :: enuc
    REAL(dp), DIMENSION(3), INTENT(OUT), &
      OPTIONAL                               :: denuc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'corecore', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: a, c, d, de, expa, x

!
!  (c+d*x*x) * exp[-(a*x)]
!

    a = SQRT(sepi%crep(1)*sepj%crep(1))
    c = 0.5_dp*(sepi%crep(2)+sepj%crep(2))
    d = 0.5_dp*(sepi%crep(3)+sepj%crep(3))

    x = SQRT(rij(1)**2 + rij(2)**2 + rij(3)**2)
    expa = EXP(-(a*x))

    enuc = (c+d*x**2)*expa

    IF(PRESENT(denuc)) THEN
       !de = (d-(c+d*x)*(a+2._dp*b*x)) * expa
       de = (-a*c + (-a*x+2.0_dp)*(d*x)) * expa
       denuc = de * rij/x
    END IF

  END SUBROUTINE corecore

END MODULE scptb_core_interactions

