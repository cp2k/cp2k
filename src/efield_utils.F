!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief all routins needed for a nonperiodic  electric field
! *****************************************************************************


MODULE efield_utils
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type,&
                                             efield_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: constant_env,&
                                             gaussian_env,&
                                             ramp_env
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_types,                        ONLY: pw_p_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'efield_utils'

! *** Public subroutines ***

  PUBLIC :: efield_potential,&
            calculate_ecore_efield

CONTAINS

! *****************************************************************************
!> \brief computes the time dependend potential on the grid
!> \param qs_env ...
!> \param v_efield_rspace ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE efield_potential(qs_env,v_efield_rspace,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type)                          :: v_efield_rspace
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'efield_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k
    INTEGER, DIMENSION(2, 3)                 :: bo_global, bo_local
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dvol, efield_ener, field(3)
    REAL(kind=dp), DIMENSION(3)              :: dr, grid_p
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho

    NULLIFY(dft_control, para_env)
    failure=.FALSE.
    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,&
                    energy=energy,&
                    rho=rho,&
                    dft_control=dft_control,&
                    para_env=para_env,&
                    error=error)

    v_efield_rspace%pw%cr3d=0.0_dp

    bo_local=v_efield_rspace%pw%pw_grid%bounds_local
    bo_global=v_efield_rspace%pw%pw_grid%bounds

    dvol=v_efield_rspace%pw%pw_grid%dvol
    dr=v_efield_rspace%pw%pw_grid%dr

    CALL make_field(dft_control,field,qs_env%sim_step,qs_env%sim_time,error)

    DO k=bo_local(1,3),bo_local(2,3)
       DO j=bo_local(1,2),bo_local(2,2)
          DO i=bo_local(1,1),bo_local(2,1)
             grid_p(1)=(i-bo_global(1,1))*dr(1)
             grid_p(2)=(j-bo_global(1,2))*dr(2)
             grid_p(3)=(k-bo_global(1,3))*dr(3)
             v_efield_rspace%pw%cr3d(i,j,k)= v_efield_rspace%pw%cr3d(i,j,k)+DOT_PRODUCT(field(:),grid_p(:))
          END DO
       END DO
    END DO
    efield_ener=0.0_dp
    DO i=1,dft_control%nspins
       efield_ener=efield_ener+accurate_sum(v_efield_rspace%pw%cr3d*rho%rho_r(i)%pw%cr3d)*dvol
    END DO
    CALL mp_sum(efield_ener,para_env%group)
    energy%efield=efield_ener
    CALL timestop(handle)


  END SUBROUTINE efield_potential

! *****************************************************************************
!> \brief computes the amplitude of the efield within a given envelop
!> \param dft_control ...
!> \param field ...
!> \param sim_step ...
!> \param sim_time ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE make_field(dft_control,field,sim_step,sim_time,error)
    TYPE(dft_control_type)                   :: dft_control
    REAL(dp)                                 :: field(3)
    INTEGER                                  :: sim_step
    REAL(KIND=dp)                            :: sim_time
    TYPE(cp_error_type)                      :: error

    INTEGER                                  :: i, nfield
    REAL(dp)                                 :: c, env, nu, pol(3), strength
    TYPE(efield_type), POINTER               :: efield

    c=137.03599962875_dp
    field=0._dp
    nfield=SIZE(dft_control%efield_fields)
    DO i=1,nfield
       efield=>dft_control%efield_fields(i)%efield
       nu=c/(efield%wavelength)
       strength=SQRT(efield%strength/(3.50944_dp*10.0_dp**16))
       IF(DOT_PRODUCT(efield%polarisation,efield%polarisation)==0)THEN
          pol(:)=1.0_dp/3.0_dp
       ELSE
          pol(:)=efield%polarisation(:)/(SQRT(DOT_PRODUCT(efield%polarisation,efield%polarisation)))
       END IF
       IF(efield%envelop_id==constant_env)THEN
          IF(sim_step.GE.efield%envelop_i_vars(1).AND.&
               (sim_step.LE.efield%envelop_i_vars(2).OR.efield%envelop_i_vars(2).LT.0))THEN
             field=field+strength*COS(sim_time*nu*2.0_dp*pi+&
                  efield%phase_offset*pi)*pol(:)
          END IF
       END IF
       IF(efield%envelop_id==ramp_env)THEN
! check this part, it seems to be wrong
          IF(sim_step.GE.efield%envelop_i_vars(1).AND.sim_step.GE.efield%envelop_i_vars(2)) &
               strength=strength*(sim_step-efield%envelop_i_vars(1))/(efield%envelop_i_vars(2)-efield%envelop_i_vars(1))
          IF(sim_step.GE.efield%envelop_i_vars(3).AND.sim_step.GE.efield%envelop_i_vars(4)) &
               strength=strength*(efield%envelop_i_vars(4)-sim_step)/(efield%envelop_i_vars(4)-efield%envelop_i_vars(3))
          IF(sim_step.GT.efield%envelop_i_vars(4).AND.efield%envelop_i_vars(4).gt.0)strength=0.0_dp
          field=field+strength*COS(sim_time*nu*2.0_dp*pi+&
               efield%phase_offset*pi)*pol(:)
       END IF
       IF(efield%envelop_id==gaussian_env)THEN
          env=EXP(-0.5_dp*((sim_time-efield%envelop_r_vars(1))/efield%envelop_r_vars(2))**2.0_dp)
          field=field+strength*env*COS(sim_time*nu*2.0_dp*pi+&
               efield%phase_offset*pi)*pol(:)
       END IF
    END DO

  END SUBROUTINE make_field

! *****************************************************************************
!> \brief computes the force and the enrgy due to a efield on the cores
!> \param qs_env ...
!> \param calculate_forces ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calculate_ecore_efield(qs_env,calculate_forces,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, OPTIONAL                        :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_ecore_efield', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, ikind, &
                                                natom, nkind
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: failure, my_force
    REAL(KIND=dp)                            :: efield_ener, zeff
    REAL(KIND=dp), DIMENSION(3)              :: field, r
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    NULLIFY(dft_control)
    CALL timeset(routineN,handle)
    CALL get_qs_env(qs_env, dft_control=dft_control, error=error)
    IF(dft_control%apply_efield_field)THEN

       failure=.FALSE.

       my_force=.FALSE.
       IF(PRESENT(calculate_forces))my_force=calculate_forces

       CALL get_qs_env(qs_env=qs_env,&
            atomic_kind_set=atomic_kind_set,&
            qs_kind_set=qs_kind_set,&
            energy=energy,&
            particle_set=particle_set,&
            cell=cell,&
            error=error)
       efield_ener=0.0_dp
       nkind = SIZE(atomic_kind_set)
       CALL make_field(dft_control,field,qs_env%sim_step,qs_env%sim_time,error)

       DO ikind=1,SIZE(atomic_kind_set)
          CALL get_atomic_kind(atomic_kind_set(ikind),atom_list=list,natom=natom)
          CALL get_qs_kind(qs_kind_set(ikind),zeff=zeff)

          natom = SIZE(list)
          DO iatom=1,natom

             atom_a = list(iatom)
             r(:) = pbc(particle_set(atom_a)%r(:),cell)
             efield_ener=efield_ener-zeff*DOT_PRODUCT(r,field)
              IF(my_force)THEN
                CALL get_qs_env(qs_env=qs_env,force=force,error=error)
                force(ikind)%efield(:,iatom)=force(ikind)%efield(:,iatom)+field*zeff
             END IF
          END DO

       END DO
       energy%efield_core=efield_ener
    END IF
    CALL timestop(handle)
  END SUBROUTINE calculate_ecore_efield
END MODULE efield_utils
