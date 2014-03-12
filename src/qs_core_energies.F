!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the energies concerning the core charge distribution
!> \par History
!>      - Full refactoring of calculate_ecore and calculate_ecore_overlap (jhu)
!> \author Matthias Krack (27.04.2001)
! *****************************************************************************
MODULE qs_core_energies

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE atprop_types,                    ONLY: atprop_array_init,&
                                             atprop_type
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_trace
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE erf_fn,                          ONLY: erfc
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: oorootpi,&
                                             twopi
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_core_energies'

  PUBLIC :: calculate_ecore,&
            calculate_ecore_overlap,&
            calculate_ecore_self

CONTAINS

! *****************************************************************************
!> \brief  Calculate the core Hamiltonian energy which includes the kinetic
!>          and the potential energy of the electrons. It is assumed, that
!>          the core Hamiltonian matrix h and the density matrix p have the
!>          same sparse matrix structure (same atomic blocks and block
!>          ordering)
!> \author  MK
!> \date    03.05.2001
!> \par History
!>         - simplified taking advantage of new non-redundant matrix
!>           structure (27.06.2003,MK)
!>         - simplified using DBCSR trace function (21.07.2010, jhu)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE calculate_ecore(h,p,ecore,error)

    TYPE(cp_dbcsr_type), POINTER             :: h, p
    REAL(KIND=dp), INTENT(OUT)               :: ecore
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_ecore', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    ecore = 0.0_dp
    CALL cp_dbcsr_trace(h,p,ecore,error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_ecore

! *****************************************************************************
!> \brief   Calculate the overlap energy of the core charge distribution.
!> \author  MK
!> \date    30.04.2001
!> \par History
!>       - Force calculation added (03.06.2002,MK)
!>       - Parallelized using a list of local atoms for rows and
!>         columns (19.07.2003,MK)
!>       - Use precomputed neighborlists (sab_core) and nl iterator (28.07.2010,jhu)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE calculate_ecore_overlap(qs_env,para_env,calculate_forces,molecular, &
                                     E_overlap_core,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    LOGICAL, INTENT(IN), OPTIONAL            :: molecular
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: E_overlap_core
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_ecore_overlap', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, group, &
                                                handle, iatom, ikind, jatom, &
                                                jkind, natom, nkind, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(3)                    :: periodic
    LOGICAL                                  :: atenergy, failure, &
                                                only_molecule, use_virial
    REAL(KIND=dp)                            :: aab, dab, eab, ecore_overlap, &
                                                f, fab, rab2, rootaab, zab
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: alpha, radius, zeff
    REAL(KIND=dp), DIMENSION(3)              :: deab, rab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(atprop_type), POINTER               :: atprop
    TYPE(cell_type), POINTER                 :: cell
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_core
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.

    CALL timeset(routineN,handle)

    NULLIFY (atomic_kind)
    NULLIFY (atomic_kind_set)
    NULLIFY (cell)
    NULLIFY (energy)
    NULLIFY (atprop)
    NULLIFY (force)
    NULLIFY (particle_set)

    group = para_env%group

    only_molecule = .FALSE.
    IF (PRESENT(molecular)) only_molecule = molecular

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    particle_set=particle_set,&
                    energy=energy,&
                    force=force,&
                    sab_core=sab_core,&
                    atprop=atprop,&
                    virial = virial,error=error)

    CALL get_cell(cell=cell,periodic=periodic)

    ! Allocate work storage
    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    ALLOCATE (alpha(nkind),radius(nkind),zeff(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    alpha(:) = 0.0_dp
    radius(:) = 0.0_dp
    zeff(:) = 0.0_dp

    IF (calculate_forces) THEN
       ALLOCATE (atom_of_kind(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)
    END IF

    atenergy = .FALSE.
    IF (ASSOCIATED(atprop)) THEN
       IF (atprop%energy) THEN
          atenergy = .TRUE.
          CALL atprop_array_init(atprop%atecc,natom,error)
       END IF
    END IF

    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           alpha_core_charge=alpha(ikind),&
                           core_charge_radius=radius(ikind),&
                           zeff=zeff(ikind))
    END DO

    ecore_overlap = 0.0_dp

    CALL neighbor_list_iterator_create(nl_iterator,sab_core)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,r=rab)
       zab = zeff(ikind)*zeff(jkind)
       aab = alpha(ikind)*alpha(jkind)/(alpha(ikind) + alpha(jkind))
       rootaab = SQRT(aab)
       fab = 2.0_dp*oorootpi*zab*rootaab
       rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       IF (rab2 > 1.e-8_dp) THEN
          IF (ikind==jkind .AND. iatom==jatom) THEN
            f=0.5_dp
          ELSE
            f=1.0_dp
          END IF
          dab = SQRT(rab2)
          eab = zab*erfc(rootaab*dab)/dab
          ecore_overlap = ecore_overlap + f*eab
          IF (atenergy) THEN
             atprop%atecc(iatom) = atprop%atecc(iatom) + 0.5_dp*f*eab
             atprop%atecc(jatom) = atprop%atecc(jatom) + 0.5_dp*f*eab
          END IF
          IF (calculate_forces) THEN
             deab(:) = rab(:)*f*(eab + fab*EXP(-aab*rab2))/rab2
             atom_a = atom_of_kind(iatom)
             atom_b = atom_of_kind(jatom)
             force(ikind)%core_overlap(:,atom_a) = force(ikind)%core_overlap(:,atom_a) + deab(:)
             force(jkind)%core_overlap(:,atom_b) = force(jkind)%core_overlap(:,atom_b) - deab(:)
             IF (use_virial) THEN
                CALL virial_pair_force ( virial%pv_virial, 1._dp, deab, rab, error)
             END IF
          END IF
       END IF
    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    DEALLOCATE (alpha,radius,zeff,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (calculate_forces) THEN
       DEALLOCATE (atom_of_kind,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL mp_sum(ecore_overlap,group)

    energy%core_overlap = ecore_overlap

    IF (PRESENT(E_overlap_core)) THEN
      E_overlap_core = energy%core_overlap
    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_ecore_overlap

! *****************************************************************************
!> \brief   Calculate the self energy of the core charge distribution.
!> \author  MK
!> \date    27.04.2001
!> \version 1.0
! *****************************************************************************
  SUBROUTINE calculate_ecore_self(qs_env,E_self_core,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: E_self_core
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_ecore_self', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, ikind, &
                                                iparticle_local, natom, &
                                                nparticle_local
    REAL(KIND=dp)                            :: alpha_core_charge, &
                                                ecore_self, es, zeff
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,energy=energy,error=error)

    ecore_self = 0.0_dp

    DO ikind=1,SIZE(atomic_kind_set)
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom, zeff=zeff,&
                           alpha_core_charge=alpha_core_charge)
      ecore_self = ecore_self - REAL(natom,dp)*zeff**2*SQRT(alpha_core_charge)
    END DO

    energy%core_self = ecore_self/SQRT(twopi)
    IF (PRESENT(E_self_core)) THEN
      E_self_core = energy%core_self
    END IF

    IF (ASSOCIATED(qs_env%atprop)) THEN
      IF (qs_env%atprop%energy) THEN
        ! atomic energy
        CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,&
                        local_particles=local_particles,error=error)
        natom = SIZE(particle_set)
        CALL atprop_array_init(qs_env%atprop%ateself,natom,error)

        DO ikind=1,SIZE(atomic_kind_set)
           atomic_kind => atomic_kind_set(ikind)
           nparticle_local = local_particles%n_el(ikind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                zeff=zeff,alpha_core_charge=alpha_core_charge)
           es = zeff**2*SQRT(alpha_core_charge)/SQRT(twopi)
           DO iparticle_local = 1, nparticle_local
              iatom = local_particles%list(ikind)%array(iparticle_local)
              qs_env%atprop%ateself(iatom) = qs_env%atprop%ateself(iatom) - es
           END DO
        END DO
      END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_ecore_self

END MODULE qs_core_energies
