!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Util force_env module
!> \author Teodoro Laino [tlaino] - 02.2011
! *****************************************************************************
MODULE force_env_utils

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE cell_types,                      ONLY: cell_type
  USE constraint,                      ONLY: rattle_control,&
                                             shake_control
  USE constraint_util,                 ONLY: getold
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: det_3x3,&
                                             jacobi
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_types_new,              ONLY: global_constraint_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: update_particle_set
  USE physcon,                         ONLY: pascal
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_env_utils'

  PUBLIC :: force_env_shake,&
            force_env_rattle,&
            rescale_forces,&
            write_stress_tensor,&
            write_forces

CONTAINS

! *****************************************************************************
!> \brief perform shake (enforcing of constraints)
!> \param force_env the force env to shake
!> \param dt the dt for shake (if you are not interested in the velocities
!>        it can be any positive number)
!> \param shake_tol the tolerance for shake
!> \param log_unit if >0 then some information on the shake is printed,
!>        defaults to -1
!> \param lagrange_mult ...
!> \param dump_lm ...
!> \param pos ...
!> \param vel ...
!> \param compold ...
!> \param reset ...
!> \author fawzi
! *****************************************************************************
  SUBROUTINE force_env_shake(force_env,dt,shake_tol,log_unit,lagrange_mult,dump_lm,&
       pos,vel,compold,reset)

    TYPE(force_env_type), POINTER            :: force_env
    REAL(kind=dp), INTENT(IN), OPTIONAL      :: dt
    REAL(kind=dp), INTENT(IN)                :: shake_tol
    INTEGER, INTENT(in), OPTIONAL            :: log_unit, lagrange_mult
    LOGICAL, INTENT(IN), OPTIONAL            :: dump_lm
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL, TARGET        :: pos, vel
    LOGICAL, INTENT(IN), OPTIONAL            :: compold, reset

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_shake', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iparticle, iparticle_kind, iparticle_local, j, &
      my_lagrange_mult, my_log_unit, nparticle_kind, nparticle_local
    LOGICAL                                  :: failure, has_pos, has_vel, &
                                                my_dump_lm
    REAL(KIND=dp)                            :: mydt
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: my_pos, my_vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.
    CALL timeset(routineN,handle)
    CPASSERT(ASSOCIATED(force_env))
    CPASSERT(force_env%ref_count>0)
    my_log_unit=-1
    IF (PRESENT(log_unit)) my_log_unit=log_unit
    my_lagrange_mult=-1
    IF (PRESENT(lagrange_mult)) my_lagrange_mult=lagrange_mult
    my_dump_lm = .FALSE.
    IF (PRESENT(dump_lm)) my_dump_lm = dump_lm
    NULLIFY(subsys,cell,molecules,molecule_kinds,local_molecules,particles,&
         my_pos,my_vel,gci)
    IF (PRESENT(pos)) my_pos => pos
    IF (PRESENT(vel)) my_vel => vel
    mydt = 0.1_dp
    IF (PRESENT(dt)) mydt = dt
    CALL force_env_get(force_env,subsys=subsys,cell=cell)
    CALL cp_subsys_get(subsys, &
         atomic_kinds=atomic_kinds,&
         local_molecules_new=local_molecules,&
         local_particles=local_particles,&
         molecules_new=molecules,&
         molecule_kinds_new=molecule_kinds,&
         particles=particles,&
         gci=gci)
    nparticle_kind = atomic_kinds%n_els
    IF (PRESENT(compold)) THEN
       IF (compold) THEN
          CALL getold( gci, local_molecules, molecules%els, molecule_kinds%els,&
               particles%els, cell)
       END IF
    END IF
    has_pos=.FALSE.
    IF (.NOT.ASSOCIATED(my_pos)) THEN
       has_pos=.TRUE.
       ALLOCATE(my_pos(3,particles%n_els))
       my_pos = 0.0_dp
       DO iparticle_kind=1,nparticle_kind
          nparticle_local = local_particles%n_el(iparticle_kind)
          DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             my_pos (:,iparticle) = particles%els(iparticle)%r(:)
          END DO
       END DO
    END IF
    has_vel=.FALSE.
    IF (.NOT.ASSOCIATED(my_vel)) THEN
       has_vel=.TRUE.
       ALLOCATE(my_vel(3,particles%n_els))
       my_vel = 0.0_dp
       DO iparticle_kind=1,nparticle_kind
          nparticle_local = local_particles%n_el(iparticle_kind)
          DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             my_vel (:,iparticle) = particles%els(iparticle)%v(:)
          END DO
       END DO
    END IF

    CALL shake_control( gci=gci, local_molecules=local_molecules,&
         molecule_set=molecules%els, molecule_kind_set=molecule_kinds%els, &
         particle_set=particles%els, pos=my_pos, vel=my_vel, dt=mydt,&
         shake_tol=shake_tol, log_unit=my_log_unit, lagrange_mult=my_lagrange_mult,&
         dump_lm= my_dump_lm, cell=cell,group=force_env%para_env%group,&
         local_particles=local_particles)

    ! Possibly reset the lagrange multipliers
    IF (PRESENT(reset)) THEN
       IF (reset) THEN
          ! Reset Intramolecular constraints
          DO i = 1,SIZE(molecules%els)
             IF (ASSOCIATED(molecules%els(i)%lci%lcolv)) THEN
                DO j = 1, SIZE(molecules%els(i)%lci%lcolv)
                   ! Reset langrange multiplier
                   molecules%els(i)%lci%lcolv(j)%lambda = 0.0_dp
                END DO
             END IF
             IF (ASSOCIATED(molecules%els(i)%lci%lg3x3)) THEN
                DO j = 1, SIZE(molecules%els(i)%lci%lg3x3)
                   ! Reset langrange multiplier
                   molecules%els(i)%lci%lg3x3(j)%lambda = 0.0_dp
                END DO
             END IF
             IF (ASSOCIATED(molecules%els(i)%lci%lg4x6)) THEN
                DO j = 1, SIZE(molecules%els(i)%lci%lg4x6)
                   ! Reset langrange multiplier
                   molecules%els(i)%lci%lg4x6(j)%lambda = 0.0_dp
                END DO
             END IF
          END DO
          ! Reset Intermolecular constraints
          IF (ASSOCIATED(gci)) THEN
             IF (ASSOCIATED(gci%lcolv)) THEN
                DO j = 1, SIZE(gci%lcolv)
                   ! Reset langrange multiplier
                   gci%lcolv(j)%lambda = 0.0_dp
                END DO
             END IF
             IF (ASSOCIATED(gci%lg3x3)) THEN
                DO j = 1, SIZE(gci%lg3x3)
                   ! Reset langrange multiplier
                   gci%lg3x3(j)%lambda = 0.0_dp
                END DO
             END IF
             IF (ASSOCIATED(gci%lg4x6)) THEN
                DO j = 1, SIZE(gci%lg4x6)
                   ! Reset langrange multiplier
                   gci%lg4x6(j)%lambda = 0.0_dp
                END DO
             END IF
          END IF
       END IF
    END IF

    IF (has_pos) THEN
       CALL update_particle_set ( particles%els, force_env%para_env%group, pos=my_pos)
       DEALLOCATE(my_pos)
    END IF
    IF (has_vel) THEN
       CALL update_particle_set ( particles%els, force_env%para_env%group, vel=my_vel)
       DEALLOCATE(my_vel)
    END IF
    CALL timestop(handle)
  END SUBROUTINE force_env_shake

! *****************************************************************************
!> \brief perform rattle (enforcing of constraints on velocities)
!>      This routine can be easily adapted to performe rattle on whatever
!>      other vector different from forces..
!> \param force_env the force env to shake
!> \param dt the dt for shake (if you are not interested in the velocities
!>        it can be any positive number)
!> \param shake_tol the tolerance for shake
!> \param log_unit if >0 then some information on the shake is printed,
!>        defaults to -1
!> \param lagrange_mult ...
!> \param dump_lm ...
!> \param vel ...
!> \param reset ...
!> \author tlaino
! *****************************************************************************
  SUBROUTINE force_env_rattle(force_env,dt,shake_tol,log_unit,lagrange_mult,dump_lm,&
       vel,reset)

    TYPE(force_env_type), POINTER            :: force_env
    REAL(kind=dp), INTENT(in), OPTIONAL      :: dt
    REAL(kind=dp), INTENT(in)                :: shake_tol
    INTEGER, INTENT(in), OPTIONAL            :: log_unit, lagrange_mult
    LOGICAL, INTENT(IN), OPTIONAL            :: dump_lm
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL, TARGET        :: vel
    LOGICAL, INTENT(IN), OPTIONAL            :: reset

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_rattle', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iparticle, iparticle_kind, iparticle_local, j, &
      my_lagrange_mult, my_log_unit, nparticle_kind, nparticle_local
    LOGICAL                                  :: failure, has_vel, my_dump_lm
    REAL(KIND=dp)                            :: mydt
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: my_vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.

    CALL timeset(routineN,handle)
    CPASSERT(ASSOCIATED(force_env))
    CPASSERT(force_env%ref_count>0)
    my_log_unit=-1
    IF (PRESENT(log_unit)) my_log_unit=log_unit
    my_lagrange_mult=-1
    IF (PRESENT(lagrange_mult)) my_lagrange_mult=lagrange_mult
    my_dump_lm = .FALSE.
    IF (PRESENT(dump_lm)) my_dump_lm = dump_lm
    NULLIFY(subsys,cell,molecules,molecule_kinds,local_molecules,particles,&
         my_vel)
    IF (PRESENT(vel)) my_vel => vel
    mydt = 0.1_dp
    IF (PRESENT(dt)) mydt = dt
    CALL force_env_get(force_env,subsys=subsys,cell=cell)
    CALL cp_subsys_get(subsys, &
         atomic_kinds=atomic_kinds,&
         local_molecules_new=local_molecules,&
         local_particles=local_particles,&
         molecules_new=molecules,&
         molecule_kinds_new=molecule_kinds,&
         particles=particles,&
         gci=gci)
    nparticle_kind = atomic_kinds%n_els
    has_vel=.FALSE.
    IF (.NOT.ASSOCIATED(my_vel)) THEN
       has_vel=.TRUE.
       ALLOCATE(my_vel(3,particles%n_els))
       my_vel = 0.0_dp
       DO iparticle_kind=1,nparticle_kind
          nparticle_local = local_particles%n_el(iparticle_kind)
          DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             my_vel (:,iparticle) = particles%els(iparticle)%v(:)
          END DO
       END DO
    END IF

    CALL rattle_control( gci=gci, local_molecules=local_molecules,&
         molecule_set=molecules%els, molecule_kind_set=molecule_kinds%els, &
         particle_set=particles%els, vel=my_vel, dt=mydt,&
         rattle_tol=shake_tol, log_unit=my_log_unit, lagrange_mult=my_lagrange_mult,&
         dump_lm=my_dump_lm, cell=cell, group=force_env%para_env%group,&
         local_particles=local_particles)

    ! Possibly reset the lagrange multipliers
    IF (PRESENT(reset)) THEN
       IF (reset) THEN
          ! Reset Intramolecular constraints
          DO i = 1,SIZE(molecules%els)
             IF (ASSOCIATED(molecules%els(i)%lci%lcolv)) THEN
                DO j = 1, SIZE(molecules%els(i)%lci%lcolv)
                   ! Reset langrange multiplier
                   molecules%els(i)%lci%lcolv(j)%lambda = 0.0_dp
                END DO
             END IF
             IF (ASSOCIATED(molecules%els(i)%lci%lg3x3)) THEN
                DO j = 1, SIZE(molecules%els(i)%lci%lg3x3)
                   ! Reset langrange multiplier
                   molecules%els(i)%lci%lg3x3(j)%lambda = 0.0_dp
                END DO
             END IF
             IF (ASSOCIATED(molecules%els(i)%lci%lg4x6)) THEN
                DO j = 1, SIZE(molecules%els(i)%lci%lg4x6)
                   ! Reset langrange multiplier
                   molecules%els(i)%lci%lg4x6(j)%lambda = 0.0_dp
                END DO
             END IF
          END DO
          ! Reset Intermolecular constraints
          IF (ASSOCIATED(gci)) THEN
             IF (ASSOCIATED(gci%lcolv)) THEN
                DO j = 1, SIZE(gci%lcolv)
                   ! Reset langrange multiplier
                   gci%lcolv(j)%lambda = 0.0_dp
                END DO
             END IF
             IF (ASSOCIATED(gci%lg3x3)) THEN
                DO j = 1, SIZE(gci%lg3x3)
                   ! Reset langrange multiplier
                   gci%lg3x3(j)%lambda = 0.0_dp
                END DO
             END IF
             IF (ASSOCIATED(gci%lg4x6)) THEN
                DO j = 1, SIZE(gci%lg4x6)
                   ! Reset langrange multiplier
                   gci%lg4x6(j)%lambda = 0.0_dp
                END DO
             END IF
          END IF
       END IF
    END IF

    IF (has_vel) THEN
       CALL update_particle_set ( particles%els, force_env%para_env%group, vel=my_vel)
    END IF
    DEALLOCATE(my_vel)
    CALL timestop(handle)
  END SUBROUTINE force_env_rattle

! *****************************************************************************
!> \brief Rescale forces if requested
!> \param force_env the force env to shake
!> \author tlaino
! *****************************************************************************
  SUBROUTINE rescale_forces (force_env)
    TYPE(force_env_type), POINTER            :: force_env

    CHARACTER(len=*), PARAMETER :: routineN = 'rescale_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iparticle
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: force(3), max_value, mod_force
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: rescale_force_section

    failure=.FALSE.
    CALL timeset(routineN,handle)
    CPASSERT(ASSOCIATED(force_env))
    CPASSERT(force_env%ref_count>0)
    rescale_force_section => section_vals_get_subs_vals(force_env%force_env_section,"RESCALE_FORCES")
    CALL section_vals_get(rescale_force_section, explicit=explicit)
    IF (.NOT.failure.AND.explicit) THEN
       CALL section_vals_val_get(rescale_force_section,"MAX_FORCE",r_val=max_value)
       CALL force_env_get(force_env,subsys=subsys)
       CALL cp_subsys_get(subsys,particles=particles)
       DO iparticle = 1, SIZE(particles%els)
          force = particles%els(iparticle)%f(:)
          mod_force = SQRT(DOT_PRODUCT(force,force))
          IF ((mod_force > max_value).AND.(mod_force /= 0.0_dp)) THEN
             force = force / mod_force * max_value
             particles%els(iparticle)%f(:) = force
          END IF
       END DO
    END IF
    CALL timestop(handle)
  END SUBROUTINE rescale_forces

! *****************************************************************************
!> \brief Variable precision output of the stress tensor
!>
!> \param pv_virial ...
!> \param output_unit ...
!> \param cell ...
!> \param ndigits ...
!> \param numerical ...
!> \author MK (26.08.2010)
! *****************************************************************************
  SUBROUTINE write_stress_tensor(pv_virial,output_unit,cell,ndigits,numerical)

    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: pv_virial
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: ndigits
    LOGICAL, INTENT(IN)                      :: numerical

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_stress_tensor', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=15)                        :: fmtstr3
    CHARACTER(LEN=16)                        :: fmtstr4
    CHARACTER(LEN=22)                        :: fmtstr2
    CHARACTER(LEN=27)                        :: fmtstr5
    CHARACTER(LEN=31)                        :: fmtstr1
    INTEGER                                  :: n
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: eigval
    REAL(KIND=dp), DIMENSION(3, 3)           :: eigvec, stress_tensor

    failure = .FALSE.

    IF (output_unit > 0) THEN
      CPASSERT(ASSOCIATED(cell))
      stress_tensor(:,:) = pv_virial(:,:)/cell%deth*pascal*1.0E-9_dp
      n = MIN(MAX(1,ndigits),20)
      fmtstr1 = "(/,T2,A,/,/,T13,A1,2(  X,A1))"
      WRITE (UNIT=fmtstr1(22:23),FMT="(I2)") n + 7
      fmtstr2 = "(T3,A,T5,3(1X,F  .  ))"
      WRITE (UNIT=fmtstr2(16:17),FMT="(I2)") n + 7
      WRITE (UNIT=fmtstr2(19:20),FMT="(I2)") n
      fmtstr3 = "(/,T3,A,F  .  )"
      WRITE (UNIT=fmtstr3(10:11),FMT="(I2)") n + 8
      WRITE (UNIT=fmtstr3(13:14),FMT="(I2)") n
      IF (numerical) THEN
        WRITE (UNIT=output_unit,FMT=fmtstr1)&
          "NUMERICAL STRESS TENSOR [GPa]","X","Y","Z"
      ELSE
        WRITE (UNIT=output_unit,FMT=fmtstr1)&
          "STRESS TENSOR [GPa]","X","Y","Z"
      END IF
      WRITE (UNIT=output_unit,FMT=fmtstr2) "X",stress_tensor(1,1:3)
      WRITE (UNIT=output_unit,FMT=fmtstr2) "Y",stress_tensor(2,1:3)
      WRITE (UNIT=output_unit,FMT=fmtstr2) "Z",stress_tensor(3,1:3)
      fmtstr4 = "(/,T3,A,ES  .  )"
      WRITE (UNIT=fmtstr4(11:12),FMT="(I2)") n + 8
      WRITE (UNIT=fmtstr4(14:15),FMT="(I2)") n
      WRITE (UNIT=output_unit,FMT=fmtstr4)&
        "1/3 Trace(stress tensor): ",(stress_tensor(1,1) +&
                                      stress_tensor(2,2) +&
                                      stress_tensor(3,3))/3.0_dp,&
        "Det(stress tensor)      : ",det_3x3(stress_tensor(:,1),&
                                             stress_tensor(:,2),&
                                             stress_tensor(:,3))
      eigval(:) = 0.0_dp
      eigvec(:,:) = 0.0_dp
      CALL jacobi(stress_tensor,eigval,eigvec)
      fmtstr5 = "(/,/,T2,A,/,/,T5,3F  .  ,/)"
      WRITE (UNIT=fmtstr5(20:21),FMT="(I2)") n + 8
      WRITE (UNIT=fmtstr5(23:24),FMT="(I2)") n
      WRITE (UNIT=output_unit,FMT=fmtstr5)&
        "EIGENVECTORS AND EIGENVALUES OF THE STRESS TENSOR",&
        eigval(1:3)
      WRITE (UNIT=output_unit,FMT=fmtstr2) " ",eigvec(1,1:3)
      WRITE (UNIT=output_unit,FMT=fmtstr2) " ",eigvec(2,1:3)
      WRITE (UNIT=output_unit,FMT=fmtstr2) " ",eigvec(3,1:3)
    END IF

  END SUBROUTINE write_stress_tensor

! *****************************************************************************
!> \brief Write forces
!>
!> \param particles ...
!> \param output_unit ...
!> \param label ...
!> \param ndigits ...
!> \param total_force ...
!> \param grand_total_force ...
!> \author MK (06.09.2010)
! *****************************************************************************
  SUBROUTINE write_forces(particles,output_unit,label,ndigits,total_force,&
                          grand_total_force)

    TYPE(particle_list_type), POINTER        :: particles
    INTEGER, INTENT(IN)                      :: output_unit
    CHARACTER(LEN=*), INTENT(IN)             :: label
    INTEGER, INTENT(IN)                      :: ndigits
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: total_force
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT), OPTIONAL                :: grand_total_force

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_forces', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=23)                        :: fmtstr3
    CHARACTER(LEN=36)                        :: fmtstr2
    CHARACTER(LEN=46)                        :: fmtstr1
    INTEGER                                  :: i, ikind, iparticle, n
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (output_unit > 0) THEN
      CPASSERT(ASSOCIATED(particles))
      n = MIN(MAX(1,ndigits),20)
      fmtstr1 = "(/,T2,A,/,/,T2,A,T11,A,T18,A,T35,A1,2(  X,A1))"
      WRITE (UNIT=fmtstr1(39:40),FMT="(I2)") n + 6
      fmtstr2 = "(T2,I6,1X,I6,T21,A,T28,3(1X,F  .  ))"
      WRITE (UNIT=fmtstr2(33:34),FMT="(I2)") n
      WRITE (UNIT=fmtstr2(30:31),FMT="(I2)") n + 6
      fmtstr3 = "(T2,A,T28,4(1X,F  .  ))"
      WRITE (UNIT=fmtstr3(20:21),FMT="(I2)") n
      WRITE (UNIT=fmtstr3(17:18),FMT="(I2)") n + 6
      WRITE (UNIT=output_unit,FMT=fmtstr1)&
        label//" FORCES in [a.u.]","# Atom","Kind","Element","X","Y","Z"
      total_force(1:3) = 0.0_dp
      DO iparticle=1,particles%n_els
        ikind = particles%els(iparticle)%atomic_kind%kind_number
        IF (particles%els(iparticle)%atom_index /= 0) THEN
          i = particles%els(iparticle)%atom_index
        ELSE
          i = iparticle
        END IF
        WRITE (UNIT=output_unit,FMT=fmtstr2)&
          i,ikind,particles%els(iparticle)%atomic_kind%element_symbol,&
          particles%els(iparticle)%f(1:3)
        total_force(1:3) = total_force(1:3) + particles%els(iparticle)%f(1:3)
      END DO
      WRITE (UNIT=output_unit,FMT=fmtstr3)&
        "SUM OF "//label//" FORCES",total_force(1:3),SQRT(SUM(total_force(:)**2))
    END IF

    IF (PRESENT(grand_total_force)) THEN
      grand_total_force(1:3) = grand_total_force(1:3) + total_force(1:3)
      WRITE (UNIT=output_unit,FMT="(A)") ""
      WRITE (UNIT=output_unit,FMT=fmtstr3)&
        "GRAND TOTAL FORCE",grand_total_force(1:3),SQRT(SUM(grand_total_force(:)**2))
    END IF

  END SUBROUTINE write_forces

END MODULE force_env_utils
