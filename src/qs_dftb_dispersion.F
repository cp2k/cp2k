!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of dispersion in DFTB
!> \author JGH
! *****************************************************************************
MODULE qs_dftb_dispersion

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_control_types,                ONLY: dft_control_type,&
                                             dftb_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type,&
                                             qs_dftb_pairpot_type
  USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb_dispersion'

  PUBLIC :: calculate_dftb_dispersion

CONTAINS

! *****************************************************************************
  SUBROUTINE calculate_dftb_dispersion(qs_env,para_env,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dftb_dispersion', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iab, iatom, ikind, ilist, inode, &
      istat, jatom, jkind, natom, nkind, nlist, nnode
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: defined, failure, use_virial
    REAL(KIND=dp)                            :: a, b, c, devdw, dij, dr, &
                                                evdw, rc, rc_a, rc_b, x0ij, &
                                                xij, xp
    REAL(KIND=dp), DIMENSION(3)              :: fdij, rij
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dftb_control_type), POINTER         :: dftb_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_vdw
    TYPE(neighbor_list_type), POINTER        :: sab_vdw_neighbor_list, &
                                                sab_vdw_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_vdw_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_kind_a, dftb_kind_b
    TYPE(qs_dftb_pairpot_type), &
      DIMENSION(:, :), POINTER               :: dftb_potential
    TYPE(qs_dftb_pairpot_type), POINTER      :: dftb_param_ij
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(virial_type), POINTER               :: virial

!   ----------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    NULLIFY (atomic_kind_set,sab_vdw)

    CALL get_qs_env(qs_env=qs_env,&
                    energy=energy,&
                    atomic_kind_set=atomic_kind_set,&
                    virial=virial,&
                    dft_control=dft_control,error=error)

    energy%dispersion = 0._dp

    dftb_control => dft_control%qs_control%dftb_control

    IF ( dftb_control%dispersion ) THEN

      NULLIFY (dftb_potential)
      CALL get_qs_env(qs_env=qs_env,&
                      dftb_potential=dftb_potential,error=error)
      IF(calculate_forces) THEN
        NULLIFY (force,particle_set)
        CALL get_qs_env(qs_env=qs_env,&
                        particle_set=particle_set,&
                        force=force,error=error)
        natom = SIZE (particle_set)
        ALLOCATE (atom_of_kind(natom),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                 atom_of_kind=atom_of_kind)
        use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
      END IF

      evdw = 0._dp

      CALL get_qs_env(qs_env=qs_env,sab_vdw=sab_vdw,error=error)

      nkind = SIZE(atomic_kind_set)

      DO ikind=1,nkind
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
                              natom=natom,&
                              dftb_parameter=dftb_kind_a)
       CALL get_dftb_atom_param(dftb_kind_a,&
              defined=defined,rcdisp=rc_a)

       IF (.NOT.defined) CYCLE

       DO jkind=1,nkind
         atomic_kind => atomic_kind_set(jkind)
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
                              dftb_parameter=dftb_kind_b)
         CALL get_dftb_atom_param(dftb_kind_b,&
                defined=defined,rcdisp=rc_b)

         rc = rc_a + rc_b

         IF (.NOT.defined) CYCLE

           ! retrieve information on potential
           dftb_param_ij => dftb_potential(ikind,jkind)
           ! vdW parameter
           xij = dftb_param_ij%xij
           dij = dftb_param_ij%dij
           x0ij = dftb_param_ij%x0ij
           a = dftb_param_ij%a
           b = dftb_param_ij%b
           c = dftb_param_ij%c

           iab = ikind + nkind*(jkind - 1)
           IF (.NOT.ASSOCIATED(sab_vdw(iab)%neighbor_list_set)) CYCLE
           CALL get_neighbor_list_set(neighbor_list_set=&
                                      sab_vdw(iab)%neighbor_list_set,&
                                      nlist=nlist)
           NULLIFY(sab_vdw_neighbor_list)

           DO ilist=1,nlist
             IF ( .NOT. ASSOCIATED(sab_vdw_neighbor_list) ) THEN
               sab_vdw_neighbor_list => &
                          first_list(sab_vdw(iab)%neighbor_list_set)
             ELSE
               sab_vdw_neighbor_list => next(sab_vdw_neighbor_list)
             END IF
             sab_vdw_neighbor_list_local => sab_vdw_neighbor_list
             CALL get_neighbor_list(neighbor_list=sab_vdw_neighbor_list_local,&
                                    atom=iatom,nnode=nnode)

             sab_vdw_neighbor_node => first_node(sab_vdw_neighbor_list_local)

             DO inode=1,nnode
               CALL get_neighbor_node(neighbor_node=sab_vdw_neighbor_node,&
                                      neighbor=jatom,r=rij)

               ! vdW potential
               dr = SQRT(SUM(rij(:)**2))
               IF (dr <= rc) THEN
                 IF ( dr > x0ij ) THEN
                   ! This is the standard London contribution.
                   ! UFF1 - Eq. 20 (long-range)
                   xp = xij/dr
                   evdw = evdw + dij*(-2._dp*xp**6 + xp**12 ) 
                   IF(calculate_forces .AND. (dr > 0.001_dp)) THEN
                     devdw = dij*12._dp*(xp**6 - xp**12)/dr
                     atom_a = atom_of_kind(iatom)
                     atom_b = atom_of_kind(jatom)
                     fdij(:) = devdw * rij(:)/dr
                     force(ikind)%dispersion(:,atom_a) =&
                         force(ikind)%dispersion(:,atom_a) - fdij(:)
                     force(jkind)%dispersion(:,atom_b) =&
                         force(jkind)%dispersion(:,atom_b) + fdij(:)
                   END IF
                 ELSE
                   ! Shorter distance. 
                   ! London contribution should converge to a finite value.
                   ! Using a parabola of the form (y = A - Bx**5 -Cx**10). 
                   ! Analytic parameters by forcing energy, first and second
                   ! derivatives to be continuous.
                   evdw = evdw + A-B*dr**5-C*dr**10
                   IF(calculate_forces .AND. (dr > 0.001_dp)) THEN
                     atom_a = atom_of_kind(iatom)
                     atom_b = atom_of_kind(jatom)
                     devdw = -5*B*dr**4 -10*C*dr**9
                     fdij(:) = devdw * rij(:)/dr
                     force(ikind)%dispersion(:,atom_a) =&
                         force(ikind)%dispersion(:,atom_a) - fdij(:)
                     force(jkind)%dispersion(:,atom_b) =&
                         force(jkind)%dispersion(:,atom_b) + fdij(:)
                   END IF
                 END IF
                 IF(calculate_forces .AND. (dr > 0.001_dp) .AND. use_virial) THEN
                   CALL virial_pair_force ( virial%pv_virial, -1._dp, &
                                            fdij, rij, error)
                 END IF
               END IF

               sab_vdw_neighbor_node => next(sab_vdw_neighbor_node)

             END DO ! inode => jatom(atom B)

           END DO ! ilist => iatom(atom A)

        END DO ! jkind

     END DO ! ikind

     IF(calculate_forces) THEN
       DEALLOCATE(atom_of_kind,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     END IF

     ! set dispersion energy
     CALL mp_sum(evdw,para_env%group)
     energy%dispersion = evdw

    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_dftb_dispersion

END MODULE qs_dftb_dispersion

