!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Treats the electrostatic for multipoles (up to quadrupoles)
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! *****************************************************************************
MODULE ewalds_multipole
  USE bibliography,                    ONLY: Aguado2003,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_type
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE dg_rho0_types,                   ONLY: dg_rho0_type
  USE dg_types,                        ONLY: dg_get,&
                                             dg_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE erf_fn,                          ONLY: erfc
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                             ewald_pw_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE input_constants,                 ONLY: do_ewald_none
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE mathlib,                         ONLY: matvec_3x3
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE structure_factor_types,          ONLY: structure_factor_type
  USE structure_factors,               ONLY: structure_factor_allocate,&
                                             structure_factor_deallocate,&
                                             structure_factor_evaluate
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewalds_multipole'

  PRIVATE
  PUBLIC :: ewald_multipole_evaluate

CONTAINS

! *****************************************************************************
!> \brief computes the potential and the force for a lattice sum of multipoles
!>      (up to quadrupole)
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, &
       particle_set, local_particles, energy_local, energy_glob, e_neut, e_self,&
       task, do_forces, do_stress, do_efield, charges, dipoles, quadrupoles,&
       forces_local, forces_glob, pv_local, pv_glob, efield1, efield2, error )
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), POINTER             :: particle_set(:)
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(OUT)               :: energy_local
    REAL(KIND=dp), INTENT(INOUT)             :: energy_glob
    REAL(KIND=dp), INTENT(OUT)               :: e_neut, e_self
    LOGICAL, DIMENSION(3), INTENT(IN)        :: task
    LOGICAL, INTENT(IN)                      :: do_forces, do_stress, &
                                                do_efield
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: charges
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN), OPTIONAL                   :: quadrupoles
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: forces_local, forces_glob
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: pv_local, pv_glob
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: efield1, efield2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_evaluate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j
    LOGICAL                                  :: do_task(3), failure
    LOGICAL, DIMENSION(3, 3)                 :: my_task

    failure = .FALSE.
    CALL timeset(routineN,'I','',handle)
    CPPostcondition(ASSOCIATED(nonbond_env),cp_failure_level,routineP,error,failure)

    ! Setup the tasks (needed to skip useless parts in the real-space term)
    do_task = task
    DO i = 1, 3
       IF (do_task(i)) THEN
          SELECT CASE(i)
          CASE(1)
             do_task(1) = ANY(charges/=0.0_dp)
          CASE(2)
             do_task(2) = ANY(dipoles/=0.0_dp)
          CASE(3)
             do_task(3) = ANY(quadrupoles/=0.0_dp)
          END SELECT
       END IF
    END DO
    DO i = 1,3
       DO j =i,3
          my_task(j,i) = do_task(i).AND.do_task(j)
          my_task(i,j) = my_task(j,i)
       END DO
    END DO

    ! Compute the Real Space (Short-Range) part of the Ewald sum
    CALL ewald_multipole_SR (nonbond_env, ewald_env, cell, energy_glob, my_task,&
         do_forces, do_efield, do_stress, charges, dipoles, quadrupoles, forces_glob,&
         pv_glob, efield1, efield2, error)
    
    ! Self-Interactions corrections
    CALL ewald_multipole_self (ewald_env, local_particles, e_self, e_neut, my_task,&
         do_efield, charges, dipoles, quadrupoles, efield1, efield2, error)

    ! Compute the Reciprocal Space (Long-Range) part of the Ewald sum
    CALL ewald_multipole_LR(ewald_env, ewald_pw, cell, particle_set, &
         local_particles, energy_local, my_task, do_forces, do_efield, do_stress,&
         charges, dipoles, quadrupoles, forces_local, pv_local, efield1, efield2,&
         error)

    CALL timestop(0.0_dp, handle)
  END SUBROUTINE ewald_multipole_evaluate
  
! *****************************************************************************
!> \brief computes the potential and the force for a lattice sum of multipoles
!>      up to quadrupole - Short Range (Real Space) Term
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE ewald_multipole_SR (nonbond_env, ewald_env, cell, energy, task,&
       do_forces, do_efield, do_stress, charges, dipoles, quadrupoles, forces,&
       pv, efield1, efield2, error)
    TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    LOGICAL, DIMENSION(3, 3), INTENT(IN)     :: task
    LOGICAL, INTENT(IN)                      :: do_forces, do_efield, &
                                                do_stress
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: charges
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN), OPTIONAL                   :: quadrupoles
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: forces, pv, efield1, efield2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_SR', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a, atom_a, atom_b, b, c, d, &
                                                e, handle, i, iend, igrp, &
                                                ikind, ilist, ipair, istart, &
                                                jkind, k, nkinds, npairs
    INTEGER, DIMENSION(:, :), POINTER        :: list
    LOGICAL                                  :: do_efield1, do_efield2, &
                                                eval_forces, failure
    REAL(KIND=dp) :: alpha, ch_a, ch_b, eps0, fac, invsqpi, ir, irab2, &
      loc_energy, ptens11, ptens12, ptens13, ptens21, ptens22, ptens23, &
      ptens31, ptens32, ptens33, r, rab2, rab2_max, rcut, tij, tmp, tmp1, &
      tmp11, tmp12, tmp13, tmp2, tmp21, tmp22, tmp23, tmp31, tmp32, tmp33, &
      tmp_ab, tmp_ba
    REAL(KIND=dp), DIMENSION(0:5)            :: f
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, cvi, dp_a, dp_b, &
                                                ef_a, ef_b, fr, rab, tij_a
    REAL(KIND=dp), DIMENSION(3, 3)           :: ef2_a, ef2_b, qp_a, qp_b, &
                                                tij_ab
    REAL(KIND=dp), DIMENSION(3, 3, 3)        :: tij_abc
    REAL(KIND=dp), DIMENSION(3, 3, 3, 3)     :: tij_abcd
    REAL(KIND=dp), DIMENSION(3, 3, 3, 3, 3)  :: tij_abcde
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update, &
                                                r_last_update_pbc

    failure = .FALSE.
    CALL timeset ( routineN, handle )
    NULLIFY(nonbonded,r_last_update, r_last_update_pbc)
    ! Get nonbond_env info
    CALL fist_nonbond_env_get (nonbond_env, nonbonded=nonbonded, natom_types = nkinds,&
         r_last_update=r_last_update,r_last_update_pbc=r_last_update_pbc, error=error)
    CALL ewald_env_get (ewald_env, alpha=alpha, eps0=eps0, rcut=rcut, error=error)
    rab2_max = rcut**2
    invsqpi= 1.0_dp/SQRT(pi)
    do_efield1 = do_efield.AND.PRESENT(efield1)
    do_efield2 = do_efield.AND.PRESENT(efield2)
    IF (do_efield1) efield1 = 0.0_dp
    IF (do_efield2) efield2 = 0.0_dp
    IF (do_stress) THEN 
       ptens11 = 0.0_dp ; ptens12 = 0.0_dp ; ptens13 = 0.0_dp
       ptens21 = 0.0_dp ; ptens22 = 0.0_dp ; ptens23 = 0.0_dp
       ptens31 = 0.0_dp ; ptens32 = 0.0_dp ; ptens33 = 0.0_dp
    END IF
    ! Starting the force loop
    Lists: DO ilist=1,nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
       npairs=neighbor_kind_pair%npairs
       IF (npairs ==0) CYCLE
       list  => neighbor_kind_pair%list
       cvi   =  neighbor_kind_pair%cell_vector
       CALL matvec_3x3(cell_v, cell%hmat, cvi)
       Kind_Group_Loop: DO igrp = 1, neighbor_kind_pair%ngrp_kind
          istart  = neighbor_kind_pair%grp_kind_start(igrp)
          iend    = neighbor_kind_pair%grp_kind_end(igrp)
          ikind   = neighbor_kind_pair%ij_kind(1,igrp)
          jkind   = neighbor_kind_pair%ij_kind(2,igrp)
          Pairs: DO ipair = istart, iend
             atom_a = list(1,ipair)
             atom_b = list(2,ipair)
             rab    = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
             rab    = rab + cell_v
             rab2   = rab(1)**2 + rab(2)**2 + rab(3)**2
             IF (rab2 <= rab2_max) THEN
                ! Compute the Short Range constribution according the task
                ! these initializations can be deleted after debug
                f         = HUGE(0.0_dp)
                tij       = HUGE(0.0_dp)
                tij_a     = HUGE(0.0_dp)
                tij_ab    = HUGE(0.0_dp)
                tij_abc   = HUGE(0.0_dp)
                tij_abcd  = HUGE(0.0_dp)
                tij_abcde = HUGE(0.0_dp)
                
                r = SQRT(rab2)
                irab2 = 1.0_dp/rab2
                ir    = 1.0_dp/r
                ! Compute the radial function
                f(0) = erfc(alpha*r)*ir
                tmp  = EXP(-alpha**2*rab2) * invsqpi
                DO i = 1, 5
                   f(i) = irab2*(f(i-1)+ tmp*((2.0_dp*alpha**2)**i)/(alpha*REAL(2*i-1,KIND=dp)))
                END DO
                ! Compute the Tensor components
                eval_forces = .FALSE.
                IF (task(1,1)) THEN
                   tij         = f(0)
                   eval_forces = do_forces
                END IF
                IF (task(1,2).OR.eval_forces) THEN
                   eval_forces = .FALSE.
                   tij_a    = - rab * f(1)
                   IF (task(1,2)) eval_forces = do_forces
                END IF
                IF ((task(2,2).OR.ANY(task(3,:))).OR.eval_forces) THEN
                   eval_forces = .FALSE.
                   DO b = 1,3
                      DO a = b,3
                         tmp = rab(a)*rab(b)
                         tij_ab(a,b) = 3.0_dp*tmp*f(2)
                         IF (a==b) tij_ab(a,b) = tij_ab(a,b) - rab2*f(1)
                         ! Symmetrize the tensor
                         tij_ab(b,a) = tij_ab(a,b)
                      END DO
                   END DO
                   IF (task(2,2).OR.ANY(task(3,:))) eval_forces = do_forces
                END IF
                IF (ANY(task(3,2:3)).OR.eval_forces) THEN
                   eval_forces = .FALSE.
                   DO c = 1, 3
                      DO b = 1, 3
                         DO a = 1, 3
                            tmp = rab(a)*rab(b)*rab(c)
                            tij_abc(a,b,c) = - 15.0_dp*tmp*f(3)
                            tmp = 3.0_dp*f(2)
                            IF (b==c) tij_abc(a,b,c) = tij_abc(a,b,c) + tmp*rab(a)
                            IF (a==c) tij_abc(a,b,c) = tij_abc(a,b,c) + tmp*rab(b)
                            IF (a==b) tij_abc(a,b,c) = tij_abc(a,b,c) + tmp*rab(c)
                         END DO
                      END DO
                   END DO
                   IF (ANY(task(3,2:3))) eval_forces = do_forces
                END IF
                IF (task(3,3).OR.eval_forces) THEN
                   eval_forces = .FALSE.
                   DO d = 1, 3
                      DO c = 1, 3
                         DO b = 1, 3
                            DO a = 1, 3
                               tmp = rab(a)*rab(b)*rab(c)*rab(d)
                               tij_abcd(a,b,c,d) = 105.0_dp*tmp*f(4)
                               tmp1 = 15.0_dp*f(3)
                               tmp2 =  3.0_dp*f(2)
                               IF (c==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(a)*rab(b)
                               IF (b==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(a)*rab(c)
                               IF (a==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(b)*rab(c)
                               IF (a==b) THEN
                                  tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(c)*rab(d)
                                  IF (c==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) + tmp2
                               END IF
                               IF (a==c) THEN 
                                  tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(b)*rab(d)
                                  IF (b==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) + tmp2
                               END IF
                               IF (b==c) THEN
                                  tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(a)*rab(d)
                                  IF (a==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) + tmp2
                               END IF
                            END DO
                         END DO
                      END DO
                   END DO
                   IF (task(3,3)) eval_forces = do_forces
                END IF
                IF (eval_forces) THEN
                   eval_forces = .FALSE.
                   DO e = 1, 3
                      DO d = 1, 3
                         DO c = 1, 3
                            DO b = 1, 3
                               DO a = 1, 3
                                  tmp = rab(a)*rab(b)*rab(c)*rab(d)*rab(e)
                                  tij_abcde(a,b,c,d,e) = -945.0_dp*tmp*f(5)
                                  tmp1 = 105.0_dp*f(4)
                                  tmp2 =  15.0_dp*f(3)
                                  IF (d==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(b)*rab(c)
                                  IF (c==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(b)*rab(d)
                                  IF (b==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(c)*rab(d)
                                  IF (a==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(b)*rab(c)*rab(d)
                                  IF (a==b) THEN
                                     tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(c)*rab(d)*rab(e)
                                     IF (d==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(c)
                                     IF (c==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(d)
                                     IF (c==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(e)
                                  END IF
                                  IF (a==c) THEN
                                     tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(b)*rab(d)*rab(e)
                                     IF (b==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(e)
                                     IF (b==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(d)
                                     IF (d==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(b)
                                  END IF
                                  IF (a==d) THEN
                                     tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(b)*rab(c)*rab(e)
                                     IF (b==c) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(e)
                                     IF (b==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(c)
                                     IF (c==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(b)
                                  END IF
                                  IF (b==c) THEN
                                     tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(d)*rab(e)
                                     IF (a==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(d)
                                     IF (d==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(a)
                                  END IF
                                  IF (b==d) THEN
                                     tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(c)*rab(e)
                                     IF (a==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(c)
                                     IF (c==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(a)
                                  END IF
                                  IF (c==d) THEN
                                     tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(b)*rab(e)
                                     IF (a==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(b)
                                     IF (b==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(a)
                                  END IF
                               END DO
                            END DO
                         END DO
                      END DO
                   END DO
                END IF
                loc_energy = 0.0_dp
                ch_a = charges(atom_a)
                ch_b = charges(atom_b)
                dp_a = dipoles(:,atom_a)
                dp_b = dipoles(:,atom_b)
                qp_a = quadrupoles(:,:,atom_a)
                qp_b = quadrupoles(:,:,atom_b)
                IF (task(1,1)) THEN
                   ! Charge - Charge
                   loc_energy = loc_energy + ch_a*ch_b*tij
                   ! Forces
                   IF (do_forces) THEN
                      fr   = 0.0_dp
                      fr(1) = fr(1) - ch_a * tij_a(1) * ch_b
                      fr(2) = fr(2) - ch_a * tij_a(2) * ch_b
                      fr(3) = fr(3) - ch_a * tij_a(3) * ch_b
                   END IF
                   ! Electric fields
                   IF (do_efield) THEN
                      IF (do_efield1) THEN
                         ef_a = 0.0_dp
                         ef_b = 0.0_dp
                         ef_a(1) = ef_a(1) - tij_a(1) * ch_b
                         ef_a(2) = ef_a(2) - tij_a(1) * ch_b
                         ef_a(3) = ef_a(3) - tij_a(1) * ch_b

                         ef_b(1) = ef_b(1) + tij_a(1) * ch_a
                         ef_b(2) = ef_b(2) + tij_a(1) * ch_a
                         ef_b(3) = ef_b(3) + tij_a(1) * ch_a
                      END IF
                      IF (do_efield2) THEN
                         ef2_a = 0.0_dp
                         ef2_b = 0.0_dp
                         ef2_a(1,1) = ef2_a(1,1) - tij_ab(1,1) * ch_b
                         ef2_a(1,2) = ef2_a(1,2) - tij_ab(1,2) * ch_b
                         ef2_a(1,3) = ef2_a(1,3) - tij_ab(1,3) * ch_b
                         ef2_a(2,1) = ef2_a(2,1) - tij_ab(2,1) * ch_b
                         ef2_a(2,2) = ef2_a(2,2) - tij_ab(2,2) * ch_b
                         ef2_a(2,3) = ef2_a(2,3) - tij_ab(2,3) * ch_b
                         ef2_a(3,1) = ef2_a(3,1) - tij_ab(3,1) * ch_b
                         ef2_a(3,2) = ef2_a(3,2) - tij_ab(3,2) * ch_b
                         ef2_a(3,3) = ef2_a(3,3) - tij_ab(3,3) * ch_b

                         ef2_b(1,1) = ef2_b(1,1) + tij_ab(1,1) * ch_a
                         ef2_b(1,2) = ef2_b(1,2) + tij_ab(1,2) * ch_a
                         ef2_b(1,3) = ef2_b(1,3) + tij_ab(1,3) * ch_a
                         ef2_b(2,1) = ef2_b(2,1) + tij_ab(2,1) * ch_a
                         ef2_b(2,2) = ef2_b(2,2) + tij_ab(2,2) * ch_a
                         ef2_b(2,3) = ef2_b(2,3) + tij_ab(2,3) * ch_a
                         ef2_b(3,1) = ef2_b(3,1) + tij_ab(3,1) * ch_a
                         ef2_b(3,2) = ef2_b(3,2) + tij_ab(3,2) * ch_a
                         ef2_b(3,3) = ef2_b(3,3) + tij_ab(3,3) * ch_a
                      END IF
                   END IF
                END IF
                IF (task(2,2)) THEN
                   ! Dipole - Dipole
                   tmp= - (dp_a(1)*(tij_ab(1,1)*dp_b(1)+&
                                    tij_ab(1,2)*dp_b(2)+&
                                    tij_ab(1,3)*dp_b(3))+&
                           dp_a(2)*(tij_ab(2,1)*dp_b(1)+&
                                    tij_ab(2,2)*dp_b(2)+&
                                    tij_ab(2,3)*dp_b(3))+&
                           dp_a(3)*(tij_ab(3,1)*dp_b(1)+&
                                    tij_ab(3,2)*dp_b(2)+&
                                    tij_ab(3,3)*dp_b(3)))
                   loc_energy = loc_energy + tmp
                   ! Forces
                   IF (do_forces) THEN
                      DO k = 1, 3
                         fr(k) = fr(k) +  dp_a(k)*(tij_abc(k,1,1)*dp_b(1)+&
                                                   tij_abc(k,1,2)*dp_b(2)+&
                                                   tij_abc(k,1,3)*dp_b(3))&
                                       +  dp_a(k)*(tij_abc(k,2,1)*dp_b(1)+&
                                                   tij_abc(k,2,2)*dp_b(2)+&
                                                   tij_abc(k,2,3)*dp_b(3))&
                                       +  dp_a(k)*(tij_abc(k,3,1)*dp_b(1)+&
                                                   tij_abc(k,3,2)*dp_b(2)+&
                                                   tij_abc(k,3,3)*dp_b(3))
                      END DO
                   END IF
                   ! Electric field
                   IF (do_efield) THEN
                      IF (do_efield1) THEN
                         ef_a(1) = ef_a(1) + (tij_ab(1,1)*dp_b(1)+&
                                              tij_ab(1,2)*dp_b(2)+&
                                              tij_ab(1,3)*dp_b(3)) 
                         ef_a(2) = ef_a(2) + (tij_ab(2,1)*dp_b(1)+&
                                              tij_ab(2,2)*dp_b(2)+&
                                              tij_ab(2,3)*dp_b(3))
                         ef_a(3) = ef_a(3) + (tij_ab(3,1)*dp_b(1)+&
                                              tij_ab(3,2)*dp_b(2)+&
                                              tij_ab(3,3)*dp_b(3))

                         ef_b(1) = ef_b(1) - (tij_ab(1,1)*dp_a(1)+&
                                              tij_ab(1,2)*dp_a(2)+&
                                              tij_ab(1,3)*dp_a(3))
                         ef_b(2) = ef_b(2) - (tij_ab(2,1)*dp_a(1)+&
                                              tij_ab(2,2)*dp_a(2)+&
                                              tij_ab(2,3)*dp_a(3))
                         ef_b(3) = ef_b(3) - (tij_ab(3,1)*dp_a(1)+&
                                              tij_ab(3,2)*dp_a(2)+&
                                              tij_ab(3,3)*dp_a(3))
                      END IF
                      IF (do_efield2) THEN
                         ef2_a(1,1) = ef2_a(1,1) + (tij_abc(1,1,1)*dp_b(1)+&
                                                    tij_abc(1,1,2)*dp_b(2)+&
                                                    tij_abc(1,1,3)*dp_b(3)) 
                         ef2_a(1,2) = ef2_a(1,2) + (tij_abc(1,2,1)*dp_b(1)+&
                                                    tij_abc(1,2,2)*dp_b(2)+&
                                                    tij_abc(1,2,3)*dp_b(3)) 
                         ef2_a(1,3) = ef2_a(1,3) + (tij_abc(1,3,1)*dp_b(1)+&
                                                    tij_abc(1,3,2)*dp_b(2)+&
                                                    tij_abc(1,3,3)*dp_b(3)) 
                         ef2_a(2,1) = ef2_a(2,1) + (tij_abc(2,1,1)*dp_b(1)+&
                                                    tij_abc(2,1,2)*dp_b(2)+&
                                                    tij_abc(2,1,3)*dp_b(3)) 
                         ef2_a(2,2) = ef2_a(2,2) + (tij_abc(2,2,1)*dp_b(1)+&
                                                    tij_abc(2,2,2)*dp_b(2)+&
                                                    tij_abc(2,2,3)*dp_b(3)) 
                         ef2_a(2,3) = ef2_a(2,3) + (tij_abc(2,3,1)*dp_b(1)+&
                                                    tij_abc(2,3,2)*dp_b(2)+&
                                                    tij_abc(2,3,3)*dp_b(3)) 
                         ef2_a(3,1) = ef2_a(3,1) + (tij_abc(3,1,1)*dp_b(1)+&
                                                    tij_abc(3,1,2)*dp_b(2)+&
                                                    tij_abc(3,1,3)*dp_b(3)) 
                         ef2_a(3,2) = ef2_a(3,2) + (tij_abc(3,2,1)*dp_b(1)+&
                                                    tij_abc(3,2,2)*dp_b(2)+&
                                                    tij_abc(3,2,3)*dp_b(3)) 
                         ef2_a(3,3) = ef2_a(3,3) + (tij_abc(3,3,1)*dp_b(1)+&
                                                    tij_abc(3,3,2)*dp_b(2)+&
                                                    tij_abc(3,3,3)*dp_b(3))

                         ef2_b(1,1) = ef2_b(1,1) + (tij_abc(1,1,1)*dp_a(1)+&
                                                    tij_abc(1,1,2)*dp_a(2)+&
                                                    tij_abc(1,1,3)*dp_a(3)) 
                         ef2_b(1,2) = ef2_b(1,2) + (tij_abc(1,2,1)*dp_a(1)+&
                                                    tij_abc(1,2,2)*dp_a(2)+&
                                                    tij_abc(1,2,3)*dp_a(3)) 
                         ef2_b(1,3) = ef2_b(1,3) + (tij_abc(1,3,1)*dp_a(1)+&
                                                    tij_abc(1,3,2)*dp_a(2)+&
                                                    tij_abc(1,3,3)*dp_a(3)) 
                         ef2_b(2,1) = ef2_b(2,1) + (tij_abc(2,1,1)*dp_a(1)+&
                                                    tij_abc(2,1,2)*dp_a(2)+&
                                                    tij_abc(2,1,3)*dp_a(3)) 
                         ef2_b(2,2) = ef2_b(2,2) + (tij_abc(2,2,1)*dp_a(1)+&
                                                    tij_abc(2,2,2)*dp_a(2)+&
                                                    tij_abc(2,2,3)*dp_a(3)) 
                         ef2_b(2,3) = ef2_b(2,3) + (tij_abc(2,3,1)*dp_a(1)+&
                                                    tij_abc(2,3,2)*dp_a(2)+&
                                                    tij_abc(2,3,3)*dp_a(3)) 
                         ef2_b(3,1) = ef2_b(3,1) + (tij_abc(3,1,1)*dp_a(1)+&
                                                    tij_abc(3,1,2)*dp_a(2)+&
                                                    tij_abc(3,1,3)*dp_a(3)) 
                         ef2_b(3,2) = ef2_b(3,2) + (tij_abc(3,2,1)*dp_a(1)+&
                                                    tij_abc(3,2,2)*dp_a(2)+&
                                                    tij_abc(3,2,3)*dp_a(3)) 
                         ef2_b(3,3) = ef2_b(3,3) + (tij_abc(3,3,1)*dp_a(1)+&
                                                    tij_abc(3,3,2)*dp_a(2)+&
                                                    tij_abc(3,3,3)*dp_a(3))
                      END IF
                   END IF
                END IF
                IF (task(2,1)) THEN
                   ! Dipole - Charge
                   tmp= - ch_a*(tij_a(1)*dp_b(1)+&
                                tij_a(2)*dp_b(2)+&
                                tij_a(3)*dp_b(3))&
                        + ch_b*(tij_a(1)*dp_a(1)+&
                                tij_a(2)*dp_a(2)+&
                                tij_a(3)*dp_a(3))
                   loc_energy = loc_energy + tmp
                   IF (do_forces) THEN
                      DO k = 1, 3
                         fr(k) = fr(k) +  ch_a *(tij_ab(k,1)*dp_b(1)+&
                                                 tij_ab(k,2)*dp_b(2)+&
                                                 tij_ab(k,3)*dp_b(3))&
                                       -  ch_b *(tij_ab(k,1)*dp_b(1)+&
                                                 tij_ab(k,2)*dp_b(2)+&
                                                 tij_ab(k,3)*dp_b(3))
                      END DO                      
                   END IF
                END IF
                IF (task(3,3)) THEN
                   ! Quadrupole - Quadrupole
                   fac  = 1.0_dp/9.0_dp
                   tmp11 = qp_a(1,1)*(tij_abcd(1,1,1,1)*qp_b(1,1)+&
                                      tij_abcd(1,1,1,2)*qp_b(1,2)+&
                                      tij_abcd(1,1,1,3)*qp_b(1,3)+&
                                      tij_abcd(1,1,2,1)*qp_b(2,1)+&
                                      tij_abcd(1,1,2,2)*qp_b(2,2)+&
                                      tij_abcd(1,1,2,3)*qp_b(2,3)+&
                                      tij_abcd(1,1,3,1)*qp_b(3,1)+&
                                      tij_abcd(1,1,3,2)*qp_b(3,2)+&
                                      tij_abcd(1,1,3,3)*qp_b(3,3)) 
                   tmp12 = qp_a(1,2)*(tij_abcd(1,2,1,1)*qp_b(1,1)+&
                                      tij_abcd(1,2,1,2)*qp_b(1,2)+&
                                      tij_abcd(1,2,1,3)*qp_b(1,3)+&
                                      tij_abcd(1,2,2,1)*qp_b(2,1)+&
                                      tij_abcd(1,2,2,2)*qp_b(2,2)+&
                                      tij_abcd(1,2,2,3)*qp_b(2,3)+&
                                      tij_abcd(1,2,3,1)*qp_b(3,1)+&
                                      tij_abcd(1,2,3,2)*qp_b(3,2)+&
                                      tij_abcd(1,2,3,3)*qp_b(3,3))
                   tmp13 = qp_a(1,3)*(tij_abcd(1,3,1,1)*qp_b(1,1)+&
                                      tij_abcd(1,3,1,2)*qp_b(1,2)+&
                                      tij_abcd(1,3,1,3)*qp_b(1,3)+&
                                      tij_abcd(1,3,2,1)*qp_b(2,1)+&
                                      tij_abcd(1,3,2,2)*qp_b(2,2)+&
                                      tij_abcd(1,3,2,3)*qp_b(2,3)+&
                                      tij_abcd(1,3,3,1)*qp_b(3,1)+&
                                      tij_abcd(1,3,3,2)*qp_b(3,2)+&
                                      tij_abcd(1,3,3,3)*qp_b(3,3))
                   tmp21 = qp_a(2,1)*(tij_abcd(2,1,1,1)*qp_b(1,1)+&
                                      tij_abcd(2,1,1,2)*qp_b(1,2)+&
                                      tij_abcd(2,1,1,3)*qp_b(1,3)+&
                                      tij_abcd(2,1,2,1)*qp_b(2,1)+&
                                      tij_abcd(2,1,2,2)*qp_b(2,2)+&
                                      tij_abcd(2,1,2,3)*qp_b(2,3)+&
                                      tij_abcd(2,1,3,1)*qp_b(3,1)+&
                                      tij_abcd(2,1,3,2)*qp_b(3,2)+&
                                      tij_abcd(2,1,3,3)*qp_b(3,3))
                   tmp22 = qp_a(2,2)*(tij_abcd(2,2,1,1)*qp_b(1,1)+&
                                      tij_abcd(2,2,1,2)*qp_b(1,2)+&
                                      tij_abcd(2,2,1,3)*qp_b(1,3)+&
                                      tij_abcd(2,2,2,1)*qp_b(2,1)+&
                                      tij_abcd(2,2,2,2)*qp_b(2,2)+&
                                      tij_abcd(2,2,2,3)*qp_b(2,3)+&
                                      tij_abcd(2,2,3,1)*qp_b(3,1)+&
                                      tij_abcd(2,2,3,2)*qp_b(3,2)+&
                                      tij_abcd(2,2,3,3)*qp_b(3,3))
                   tmp23 = qp_a(2,3)*(tij_abcd(2,3,1,1)*qp_b(1,1)+&
                                      tij_abcd(2,3,1,2)*qp_b(1,2)+&
                                      tij_abcd(2,3,1,3)*qp_b(1,3)+&
                                      tij_abcd(2,3,2,1)*qp_b(2,1)+&
                                      tij_abcd(2,3,2,2)*qp_b(2,2)+&
                                      tij_abcd(2,3,2,3)*qp_b(2,3)+&
                                      tij_abcd(2,3,3,1)*qp_b(3,1)+&
                                      tij_abcd(2,3,3,2)*qp_b(3,2)+&
                                      tij_abcd(2,3,3,3)*qp_b(3,3))
                   tmp31 = qp_a(3,1)*(tij_abcd(3,1,1,1)*qp_b(1,1)+&
                                      tij_abcd(3,1,1,2)*qp_b(1,2)+&
                                      tij_abcd(3,1,1,3)*qp_b(1,3)+&
                                      tij_abcd(3,1,2,1)*qp_b(2,1)+&
                                      tij_abcd(3,1,2,2)*qp_b(2,2)+&
                                      tij_abcd(3,1,2,3)*qp_b(2,3)+&
                                      tij_abcd(3,1,3,1)*qp_b(3,1)+&
                                      tij_abcd(3,1,3,2)*qp_b(3,2)+&
                                      tij_abcd(3,1,3,3)*qp_b(3,3))
                   tmp32 = qp_a(3,2)*(tij_abcd(3,2,1,1)*qp_b(1,1)+&
                                      tij_abcd(3,2,1,2)*qp_b(1,2)+&
                                      tij_abcd(3,2,1,3)*qp_b(1,3)+&
                                      tij_abcd(3,2,2,1)*qp_b(2,1)+&
                                      tij_abcd(3,2,2,2)*qp_b(2,2)+&
                                      tij_abcd(3,2,2,3)*qp_b(2,3)+&
                                      tij_abcd(3,2,3,1)*qp_b(3,1)+&
                                      tij_abcd(3,2,3,2)*qp_b(3,2)+&
                                      tij_abcd(3,2,3,3)*qp_b(3,3))
                   tmp33 = qp_a(3,3)*(tij_abcd(3,3,1,1)*qp_b(1,1)+&
                                      tij_abcd(3,3,1,2)*qp_b(1,2)+&
                                      tij_abcd(3,3,1,3)*qp_b(1,3)+&
                                      tij_abcd(3,3,2,1)*qp_b(2,1)+&
                                      tij_abcd(3,3,2,2)*qp_b(2,2)+&
                                      tij_abcd(3,3,2,3)*qp_b(2,3)+&
                                      tij_abcd(3,3,3,1)*qp_b(3,1)+&
                                      tij_abcd(3,3,3,2)*qp_b(3,2)+&
                                      tij_abcd(3,3,3,3)*qp_b(3,3))

                   tmp= tmp11 + tmp12 + tmp13 + &
                        tmp21 + tmp22 + tmp23 + &
                        tmp31 + tmp32 + tmp33

                   loc_energy = loc_energy + fac*tmp
                   ! Forces
                   IF (do_forces) THEN
                      DO k = 1, 3
                         tmp11 = qp_a(1,1)*(tij_abcde(k,1,1,1,1)*qp_b(1,1)+&  
                                            tij_abcde(k,1,1,1,2)*qp_b(1,2)+&  
                                            tij_abcde(k,1,1,1,3)*qp_b(1,3)+&  
                                            tij_abcde(k,1,1,2,1)*qp_b(2,1)+&  
                                            tij_abcde(k,1,1,2,2)*qp_b(2,2)+&  
                                            tij_abcde(k,1,1,2,3)*qp_b(2,3)+&  
                                            tij_abcde(k,1,1,3,1)*qp_b(3,1)+&  
                                            tij_abcde(k,1,1,3,2)*qp_b(3,2)+&  
                                            tij_abcde(k,1,1,3,3)*qp_b(3,3))
                         tmp12 = qp_a(1,2)*(tij_abcde(k,1,2,1,1)*qp_b(1,1)+&  
                                            tij_abcde(k,1,2,1,2)*qp_b(1,2)+&  
                                            tij_abcde(k,1,2,1,3)*qp_b(1,3)+&  
                                            tij_abcde(k,1,2,2,1)*qp_b(2,1)+&  
                                            tij_abcde(k,1,2,2,2)*qp_b(2,2)+&  
                                            tij_abcde(k,1,2,2,3)*qp_b(2,3)+&  
                                            tij_abcde(k,1,2,3,1)*qp_b(3,1)+&  
                                            tij_abcde(k,1,2,3,2)*qp_b(3,2)+&  
                                            tij_abcde(k,1,2,3,3)*qp_b(3,3))
                         tmp13 = qp_a(1,3)*(tij_abcde(k,1,3,1,1)*qp_b(1,1)+&  
                                            tij_abcde(k,1,3,1,2)*qp_b(1,2)+&  
                                            tij_abcde(k,1,3,1,3)*qp_b(1,3)+&  
                                            tij_abcde(k,1,3,2,1)*qp_b(2,1)+&  
                                            tij_abcde(k,1,3,2,2)*qp_b(2,2)+&  
                                            tij_abcde(k,1,3,2,3)*qp_b(2,3)+&  
                                            tij_abcde(k,1,3,3,1)*qp_b(3,1)+&  
                                            tij_abcde(k,1,3,3,2)*qp_b(3,2)+&  
                                            tij_abcde(k,1,3,3,3)*qp_b(3,3))
                         tmp21 = qp_a(2,1)*(tij_abcde(k,2,1,1,1)*qp_b(1,1)+&  
                                            tij_abcde(k,2,1,1,2)*qp_b(1,2)+&  
                                            tij_abcde(k,2,1,1,3)*qp_b(1,3)+&  
                                            tij_abcde(k,2,1,2,1)*qp_b(2,1)+&  
                                            tij_abcde(k,2,1,2,2)*qp_b(2,2)+&  
                                            tij_abcde(k,2,1,2,3)*qp_b(2,3)+&  
                                            tij_abcde(k,2,1,3,1)*qp_b(3,1)+&  
                                            tij_abcde(k,2,1,3,2)*qp_b(3,2)+&  
                                            tij_abcde(k,2,1,3,3)*qp_b(3,3))
                         tmp22 = qp_a(2,2)*(tij_abcde(k,2,2,1,1)*qp_b(1,1)+&  
                                            tij_abcde(k,2,2,1,2)*qp_b(1,2)+&  
                                            tij_abcde(k,2,2,1,3)*qp_b(1,3)+&  
                                            tij_abcde(k,2,2,2,1)*qp_b(2,1)+&  
                                            tij_abcde(k,2,2,2,2)*qp_b(2,2)+&  
                                            tij_abcde(k,2,2,2,3)*qp_b(2,3)+&  
                                            tij_abcde(k,2,2,3,1)*qp_b(3,1)+&  
                                            tij_abcde(k,2,2,3,2)*qp_b(3,2)+&  
                                            tij_abcde(k,2,2,3,3)*qp_b(3,3))
                         tmp23 = qp_a(2,3)*(tij_abcde(k,2,3,1,1)*qp_b(1,1)+&  
                                            tij_abcde(k,2,3,1,2)*qp_b(1,2)+&  
                                            tij_abcde(k,2,3,1,3)*qp_b(1,3)+&  
                                            tij_abcde(k,2,3,2,1)*qp_b(2,1)+&  
                                            tij_abcde(k,2,3,2,2)*qp_b(2,2)+&  
                                            tij_abcde(k,2,3,2,3)*qp_b(2,3)+&  
                                            tij_abcde(k,2,3,3,1)*qp_b(3,1)+&  
                                            tij_abcde(k,2,3,3,2)*qp_b(3,2)+&  
                                            tij_abcde(k,2,3,3,3)*qp_b(3,3))
                         tmp31 = qp_a(3,1)*(tij_abcde(k,3,1,1,1)*qp_b(1,1)+&  
                                            tij_abcde(k,3,1,1,2)*qp_b(1,2)+&  
                                            tij_abcde(k,3,1,1,3)*qp_b(1,3)+&  
                                            tij_abcde(k,3,1,2,1)*qp_b(2,1)+&  
                                            tij_abcde(k,3,1,2,2)*qp_b(2,2)+&  
                                            tij_abcde(k,3,1,2,3)*qp_b(2,3)+&  
                                            tij_abcde(k,3,1,3,1)*qp_b(3,1)+&  
                                            tij_abcde(k,3,1,3,2)*qp_b(3,2)+&  
                                            tij_abcde(k,3,1,3,3)*qp_b(3,3))
                         tmp32 = qp_a(3,2)*(tij_abcde(k,3,2,1,1)*qp_b(1,1)+&  
                                            tij_abcde(k,3,2,1,2)*qp_b(1,2)+&  
                                            tij_abcde(k,3,2,1,3)*qp_b(1,3)+&  
                                            tij_abcde(k,3,2,2,1)*qp_b(2,1)+&  
                                            tij_abcde(k,3,2,2,2)*qp_b(2,2)+&  
                                            tij_abcde(k,3,2,2,3)*qp_b(2,3)+&  
                                            tij_abcde(k,3,2,3,1)*qp_b(3,1)+&  
                                            tij_abcde(k,3,2,3,2)*qp_b(3,2)+&  
                                            tij_abcde(k,3,2,3,3)*qp_b(3,3))
                         tmp33 = qp_a(3,3)*(tij_abcde(k,3,3,1,1)*qp_b(1,1)+&  
                                            tij_abcde(k,3,3,1,2)*qp_b(1,2)+&  
                                            tij_abcde(k,3,3,1,3)*qp_b(1,3)+&  
                                            tij_abcde(k,3,3,2,1)*qp_b(2,1)+&  
                                            tij_abcde(k,3,3,2,2)*qp_b(2,2)+&  
                                            tij_abcde(k,3,3,2,3)*qp_b(2,3)+&  
                                            tij_abcde(k,3,3,3,1)*qp_b(3,1)+&  
                                            tij_abcde(k,3,3,3,2)*qp_b(3,2)+&  
                                            tij_abcde(k,3,3,3,3)*qp_b(3,3))
                         fr(k) = fr(k) - fac * ( tmp11 + tmp12 + tmp13 +&
                                                 tmp21 + tmp22 + tmp23 +&
                                                 tmp31 + tmp32 + tmp33  )
                      END DO
                   END IF
                   ! Electric field
                   IF (do_efield) THEN
                      IF (do_efield1) THEN
                         ef_a(1) = ef_a(1) - (tij_abc(1,1,1)*qp_b(1,1)+&
                                              tij_abc(1,1,2)*qp_b(1,2)+&
                                              tij_abc(1,1,3)*qp_b(1,3)+&
                                              tij_abc(1,2,1)*qp_b(2,1)+&
                                              tij_abc(1,2,2)*qp_b(2,2)+&
                                              tij_abc(1,2,3)*qp_b(2,3)+&
                                              tij_abc(1,3,1)*qp_b(3,1)+&
                                              tij_abc(1,3,2)*qp_b(3,2)+&
                                              tij_abc(1,3,3)*qp_b(3,3)) 
                         ef_a(2) = ef_a(2) - (tij_abc(2,1,1)*qp_b(1,1)+&
                                              tij_abc(2,1,2)*qp_b(1,2)+&
                                              tij_abc(2,1,3)*qp_b(1,3)+&
                                              tij_abc(2,2,1)*qp_b(2,1)+&
                                              tij_abc(2,2,2)*qp_b(2,2)+&
                                              tij_abc(2,2,3)*qp_b(2,3)+&
                                              tij_abc(2,3,1)*qp_b(3,1)+&
                                              tij_abc(2,3,2)*qp_b(3,2)+&
                                              tij_abc(2,3,3)*qp_b(3,3)) 
                         ef_a(3) = ef_a(3) - (tij_abc(3,1,1)*qp_b(1,1)+&
                                              tij_abc(3,1,2)*qp_b(1,2)+&
                                              tij_abc(3,1,3)*qp_b(1,3)+&
                                              tij_abc(3,2,1)*qp_b(2,1)+&
                                              tij_abc(3,2,2)*qp_b(2,2)+&
                                              tij_abc(3,2,3)*qp_b(2,3)+&
                                              tij_abc(3,3,1)*qp_b(3,1)+&
                                              tij_abc(3,3,2)*qp_b(3,2)+&
                                              tij_abc(3,3,3)*qp_b(3,3)) 

                         ef_b(1) = ef_b(1) + (tij_abc(1,1,1)*qp_a(1,1)+&
                                              tij_abc(1,1,2)*qp_a(1,2)+&
                                              tij_abc(1,1,3)*qp_a(1,3)+&
                                              tij_abc(1,2,1)*qp_a(2,1)+&
                                              tij_abc(1,2,2)*qp_a(2,2)+&
                                              tij_abc(1,2,3)*qp_a(2,3)+&
                                              tij_abc(1,3,1)*qp_a(3,1)+&
                                              tij_abc(1,3,2)*qp_a(3,2)+&
                                              tij_abc(1,3,3)*qp_a(3,3)) 
                         ef_b(2) = ef_b(2) + (tij_abc(2,1,1)*qp_a(1,1)+&
                                              tij_abc(2,1,2)*qp_a(1,2)+&
                                              tij_abc(2,1,3)*qp_a(1,3)+&
                                              tij_abc(2,2,1)*qp_a(2,1)+&
                                              tij_abc(2,2,2)*qp_a(2,2)+&
                                              tij_abc(2,2,3)*qp_a(2,3)+&
                                              tij_abc(2,3,1)*qp_a(3,1)+&
                                              tij_abc(2,3,2)*qp_a(3,2)+&
                                              tij_abc(2,3,3)*qp_a(3,3)) 
                         ef_b(3) = ef_b(3) + (tij_abc(3,1,1)*qp_a(1,1)+&
                                              tij_abc(3,1,2)*qp_a(1,2)+&
                                              tij_abc(3,1,3)*qp_a(1,3)+&
                                              tij_abc(3,2,1)*qp_a(2,1)+&
                                              tij_abc(3,2,2)*qp_a(2,2)+&
                                              tij_abc(3,2,3)*qp_a(2,3)+&
                                              tij_abc(3,3,1)*qp_a(3,1)+&
                                              tij_abc(3,3,2)*qp_a(3,2)+&
                                              tij_abc(3,3,3)*qp_a(3,3))
                      END IF
                      IF (do_efield2) THEN
                         ef2_a(1,1) = ef2_a(1,1) - fac * (tij_abcd(1,1,1,1)*qp_b(1,1)+&
                                                          tij_abcd(1,1,1,2)*qp_b(1,2)+&
                                                          tij_abcd(1,1,1,3)*qp_b(1,3)+&
                                                          tij_abcd(1,1,2,1)*qp_b(2,1)+&
                                                          tij_abcd(1,1,2,2)*qp_b(2,2)+&
                                                          tij_abcd(1,1,2,3)*qp_b(2,3)+&
                                                          tij_abcd(1,1,3,1)*qp_b(3,1)+&
                                                          tij_abcd(1,1,3,2)*qp_b(3,2)+&
                                                          tij_abcd(1,1,3,3)*qp_b(3,3)) 
                         ef2_a(1,2) = ef2_a(1,2) - fac * (tij_abcd(1,2,1,1)*qp_b(1,1)+&
                                                          tij_abcd(1,2,1,2)*qp_b(1,2)+&
                                                          tij_abcd(1,2,1,3)*qp_b(1,3)+&
                                                          tij_abcd(1,2,2,1)*qp_b(2,1)+&
                                                          tij_abcd(1,2,2,2)*qp_b(2,2)+&
                                                          tij_abcd(1,2,2,3)*qp_b(2,3)+&
                                                          tij_abcd(1,2,3,1)*qp_b(3,1)+&
                                                          tij_abcd(1,2,3,2)*qp_b(3,2)+&
                                                          tij_abcd(1,2,3,3)*qp_b(3,3)) 
                         ef2_a(1,3) = ef2_a(1,3) - fac * (tij_abcd(1,3,1,1)*qp_b(1,1)+&
                                                          tij_abcd(1,3,1,2)*qp_b(1,2)+&
                                                          tij_abcd(1,3,1,3)*qp_b(1,3)+&
                                                          tij_abcd(1,3,2,1)*qp_b(2,1)+&
                                                          tij_abcd(1,3,2,2)*qp_b(2,2)+&
                                                          tij_abcd(1,3,2,3)*qp_b(2,3)+&
                                                          tij_abcd(1,3,3,1)*qp_b(3,1)+&
                                                          tij_abcd(1,3,3,2)*qp_b(3,2)+&
                                                          tij_abcd(1,3,3,3)*qp_b(3,3)) 
                         ef2_a(2,1) = ef2_a(2,1) - fac * (tij_abcd(2,1,1,1)*qp_b(1,1)+&
                                                          tij_abcd(2,1,1,2)*qp_b(1,2)+&
                                                          tij_abcd(2,1,1,3)*qp_b(1,3)+&
                                                          tij_abcd(2,1,2,1)*qp_b(2,1)+&
                                                          tij_abcd(2,1,2,2)*qp_b(2,2)+&
                                                          tij_abcd(2,1,2,3)*qp_b(2,3)+&
                                                          tij_abcd(2,1,3,1)*qp_b(3,1)+&
                                                          tij_abcd(2,1,3,2)*qp_b(3,2)+&
                                                          tij_abcd(2,1,3,3)*qp_b(3,3)) 
                         ef2_a(2,2) = ef2_a(2,2) - fac * (tij_abcd(2,2,1,1)*qp_b(1,1)+&
                                                          tij_abcd(2,2,1,2)*qp_b(1,2)+&
                                                          tij_abcd(2,2,1,3)*qp_b(1,3)+&
                                                          tij_abcd(2,2,2,1)*qp_b(2,1)+&
                                                          tij_abcd(2,2,2,2)*qp_b(2,2)+&
                                                          tij_abcd(2,2,2,3)*qp_b(2,3)+&
                                                          tij_abcd(2,2,3,1)*qp_b(3,1)+&
                                                          tij_abcd(2,2,3,2)*qp_b(3,2)+&
                                                          tij_abcd(2,2,3,3)*qp_b(3,3)) 
                         ef2_a(2,3) = ef2_a(2,3) - fac * (tij_abcd(2,3,1,1)*qp_b(1,1)+&
                                                          tij_abcd(2,3,1,2)*qp_b(1,2)+&
                                                          tij_abcd(2,3,1,3)*qp_b(1,3)+&
                                                          tij_abcd(2,3,2,1)*qp_b(2,1)+&
                                                          tij_abcd(2,3,2,2)*qp_b(2,2)+&
                                                          tij_abcd(2,3,2,3)*qp_b(2,3)+&
                                                          tij_abcd(2,3,3,1)*qp_b(3,1)+&
                                                          tij_abcd(2,3,3,2)*qp_b(3,2)+&
                                                          tij_abcd(2,3,3,3)*qp_b(3,3)) 
                         ef2_a(3,1) = ef2_a(3,1) - fac * (tij_abcd(3,1,1,1)*qp_b(1,1)+&
                                                          tij_abcd(3,1,1,2)*qp_b(1,2)+&
                                                          tij_abcd(3,1,1,3)*qp_b(1,3)+&
                                                          tij_abcd(3,1,2,1)*qp_b(2,1)+&
                                                          tij_abcd(3,1,2,2)*qp_b(2,2)+&
                                                          tij_abcd(3,1,2,3)*qp_b(2,3)+&
                                                          tij_abcd(3,1,3,1)*qp_b(3,1)+&
                                                          tij_abcd(3,1,3,2)*qp_b(3,2)+&
                                                          tij_abcd(3,1,3,3)*qp_b(3,3)) 
                         ef2_a(3,2) = ef2_a(3,2) - fac * (tij_abcd(3,2,1,1)*qp_b(1,1)+&
                                                          tij_abcd(3,2,1,2)*qp_b(1,2)+&
                                                          tij_abcd(3,2,1,3)*qp_b(1,3)+&
                                                          tij_abcd(3,2,2,1)*qp_b(2,1)+&
                                                          tij_abcd(3,2,2,2)*qp_b(2,2)+&
                                                          tij_abcd(3,2,2,3)*qp_b(2,3)+&
                                                          tij_abcd(3,2,3,1)*qp_b(3,1)+&
                                                          tij_abcd(3,2,3,2)*qp_b(3,2)+&
                                                          tij_abcd(3,2,3,3)*qp_b(3,3)) 
                         ef2_a(3,3) = ef2_a(3,3) - fac * (tij_abcd(3,3,1,1)*qp_b(1,1)+&
                                                          tij_abcd(3,3,1,2)*qp_b(1,2)+&
                                                          tij_abcd(3,3,1,3)*qp_b(1,3)+&
                                                          tij_abcd(3,3,2,1)*qp_b(2,1)+&
                                                          tij_abcd(3,3,2,2)*qp_b(2,2)+&
                                                          tij_abcd(3,3,2,3)*qp_b(2,3)+&
                                                          tij_abcd(3,3,3,1)*qp_b(3,1)+&
                                                          tij_abcd(3,3,3,2)*qp_b(3,2)+&
                                                          tij_abcd(3,3,3,3)*qp_b(3,3)) 

                         ef2_b(1,1) = ef2_b(1,1) + fac *(tij_abcd(1,1,1,1)*qp_a(1,1)+&
                                                         tij_abcd(1,1,1,2)*qp_a(1,2)+&
                                                         tij_abcd(1,1,1,3)*qp_a(1,3)+&
                                                         tij_abcd(1,1,2,1)*qp_a(2,1)+&
                                                         tij_abcd(1,1,2,2)*qp_a(2,2)+&
                                                         tij_abcd(1,1,2,3)*qp_a(2,3)+&
                                                         tij_abcd(1,1,3,1)*qp_a(3,1)+&
                                                         tij_abcd(1,1,3,2)*qp_a(3,2)+&
                                                         tij_abcd(1,1,3,3)*qp_a(3,3)) 
                         ef2_b(1,2) = ef2_b(1,2) + fac *(tij_abcd(1,2,1,1)*qp_a(1,1)+&
                                                         tij_abcd(1,2,1,2)*qp_a(1,2)+&
                                                         tij_abcd(1,2,1,3)*qp_a(1,3)+&
                                                         tij_abcd(1,2,2,1)*qp_a(2,1)+&
                                                         tij_abcd(1,2,2,2)*qp_a(2,2)+&
                                                         tij_abcd(1,2,2,3)*qp_a(2,3)+&
                                                         tij_abcd(1,2,3,1)*qp_a(3,1)+&
                                                         tij_abcd(1,2,3,2)*qp_a(3,2)+&
                                                         tij_abcd(1,2,3,3)*qp_a(3,3)) 
                         ef2_b(1,3) = ef2_b(1,3) + fac *(tij_abcd(1,3,1,1)*qp_a(1,1)+&
                                                         tij_abcd(1,3,1,2)*qp_a(1,2)+&
                                                         tij_abcd(1,3,1,3)*qp_a(1,3)+&
                                                         tij_abcd(1,3,2,1)*qp_a(2,1)+&
                                                         tij_abcd(1,3,2,2)*qp_a(2,2)+&
                                                         tij_abcd(1,3,2,3)*qp_a(2,3)+&
                                                         tij_abcd(1,3,3,1)*qp_a(3,1)+&
                                                         tij_abcd(1,3,3,2)*qp_a(3,2)+&
                                                         tij_abcd(1,3,3,3)*qp_a(3,3)) 
                         ef2_b(2,1) = ef2_b(2,1) + fac *(tij_abcd(2,1,1,1)*qp_a(1,1)+&
                                                         tij_abcd(2,1,1,2)*qp_a(1,2)+&
                                                         tij_abcd(2,1,1,3)*qp_a(1,3)+&
                                                         tij_abcd(2,1,2,1)*qp_a(2,1)+&
                                                         tij_abcd(2,1,2,2)*qp_a(2,2)+&
                                                         tij_abcd(2,1,2,3)*qp_a(2,3)+&
                                                         tij_abcd(2,1,3,1)*qp_a(3,1)+&
                                                         tij_abcd(2,1,3,2)*qp_a(3,2)+&
                                                         tij_abcd(2,1,3,3)*qp_a(3,3)) 
                         ef2_b(2,2) = ef2_b(2,2) + fac *(tij_abcd(2,2,1,1)*qp_a(1,1)+&
                                                         tij_abcd(2,2,1,2)*qp_a(1,2)+&
                                                         tij_abcd(2,2,1,3)*qp_a(1,3)+&
                                                         tij_abcd(2,2,2,1)*qp_a(2,1)+&
                                                         tij_abcd(2,2,2,2)*qp_a(2,2)+&
                                                         tij_abcd(2,2,2,3)*qp_a(2,3)+&
                                                         tij_abcd(2,2,3,1)*qp_a(3,1)+&
                                                         tij_abcd(2,2,3,2)*qp_a(3,2)+&
                                                         tij_abcd(2,2,3,3)*qp_a(3,3)) 
                         ef2_b(2,3) = ef2_b(2,3) + fac *(tij_abcd(2,3,1,1)*qp_a(1,1)+&
                                                         tij_abcd(2,3,1,2)*qp_a(1,2)+&
                                                         tij_abcd(2,3,1,3)*qp_a(1,3)+&
                                                         tij_abcd(2,3,2,1)*qp_a(2,1)+&
                                                         tij_abcd(2,3,2,2)*qp_a(2,2)+&
                                                         tij_abcd(2,3,2,3)*qp_a(2,3)+&
                                                         tij_abcd(2,3,3,1)*qp_a(3,1)+&
                                                         tij_abcd(2,3,3,2)*qp_a(3,2)+&
                                                         tij_abcd(2,3,3,3)*qp_a(3,3)) 
                         ef2_b(3,1) = ef2_b(3,1) + fac *(tij_abcd(3,1,1,1)*qp_a(1,1)+&
                                                         tij_abcd(3,1,1,2)*qp_a(1,2)+&
                                                         tij_abcd(3,1,1,3)*qp_a(1,3)+&
                                                         tij_abcd(3,1,2,1)*qp_a(2,1)+&
                                                         tij_abcd(3,1,2,2)*qp_a(2,2)+&
                                                         tij_abcd(3,1,2,3)*qp_a(2,3)+&
                                                         tij_abcd(3,1,3,1)*qp_a(3,1)+&
                                                         tij_abcd(3,1,3,2)*qp_a(3,2)+&
                                                         tij_abcd(3,1,3,3)*qp_a(3,3)) 
                         ef2_b(3,2) = ef2_b(3,2) + fac *(tij_abcd(3,2,1,1)*qp_a(1,1)+&
                                                         tij_abcd(3,2,1,2)*qp_a(1,2)+&
                                                         tij_abcd(3,2,1,3)*qp_a(1,3)+&
                                                         tij_abcd(3,2,2,1)*qp_a(2,1)+&
                                                         tij_abcd(3,2,2,2)*qp_a(2,2)+&
                                                         tij_abcd(3,2,2,3)*qp_a(2,3)+&
                                                         tij_abcd(3,2,3,1)*qp_a(3,1)+&
                                                         tij_abcd(3,2,3,2)*qp_a(3,2)+&
                                                         tij_abcd(3,2,3,3)*qp_a(3,3)) 
                         ef2_b(3,3) = ef2_b(3,3) + fac *(tij_abcd(3,3,1,1)*qp_a(1,1)+&
                                                         tij_abcd(3,3,1,2)*qp_a(1,2)+&
                                                         tij_abcd(3,3,1,3)*qp_a(1,3)+&
                                                         tij_abcd(3,3,2,1)*qp_a(2,1)+&
                                                         tij_abcd(3,3,2,2)*qp_a(2,2)+&
                                                         tij_abcd(3,3,2,3)*qp_a(2,3)+&
                                                         tij_abcd(3,3,3,1)*qp_a(3,1)+&
                                                         tij_abcd(3,3,3,2)*qp_a(3,2)+&
                                                         tij_abcd(3,3,3,3)*qp_a(3,3)) 
                      END IF
                   END IF
                END IF
                IF (task(3,2)) THEN
                   ! Quadrupole - Dipole
                   fac = 1.0_dp/3.0_dp
                   ! Dipole A - Quadrupole B
                   tmp_ab = dp_a(1)*(tij_abc(1,1,1)*qp_b(1,1)+&
                                     tij_abc(1,1,2)*qp_b(1,2)+&
                                     tij_abc(1,1,3)*qp_b(1,3)+&
                                     tij_abc(1,2,1)*qp_b(2,1)+&
                                     tij_abc(1,2,2)*qp_b(2,2)+&
                                     tij_abc(1,2,3)*qp_b(2,3)+&
                                     tij_abc(1,3,1)*qp_b(3,1)+&
                                     tij_abc(1,3,2)*qp_b(3,2)+&
                                     tij_abc(1,3,3)*qp_b(3,3))+&
                            dp_a(2)*(tij_abc(2,1,1)*qp_b(1,1)+&
                                     tij_abc(2,1,2)*qp_b(1,2)+&
                                     tij_abc(2,1,3)*qp_b(1,3)+&
                                     tij_abc(2,2,1)*qp_b(2,1)+&
                                     tij_abc(2,2,2)*qp_b(2,2)+&
                                     tij_abc(2,2,3)*qp_b(2,3)+&
                                     tij_abc(2,3,1)*qp_b(3,1)+&
                                     tij_abc(2,3,2)*qp_b(3,2)+&
                                     tij_abc(2,3,3)*qp_b(3,3))+&
                            dp_a(3)*(tij_abc(3,1,1)*qp_b(1,1)+&
                                     tij_abc(3,1,2)*qp_b(1,2)+&
                                     tij_abc(3,1,3)*qp_b(1,3)+&
                                     tij_abc(3,2,1)*qp_b(2,1)+&
                                     tij_abc(3,2,2)*qp_b(2,2)+&
                                     tij_abc(3,2,3)*qp_b(2,3)+&
                                     tij_abc(3,3,1)*qp_b(3,1)+&
                                     tij_abc(3,3,2)*qp_b(3,2)+&
                                     tij_abc(3,3,3)*qp_b(3,3))
                   ! Dipole B - Quadrupole A
                   tmp_ba = dp_b(1)*(tij_abc(1,1,1)*qp_a(1,1)+&
                                     tij_abc(1,1,2)*qp_a(1,2)+&
                                     tij_abc(1,1,3)*qp_a(1,3)+&
                                     tij_abc(1,2,1)*qp_a(2,1)+&
                                     tij_abc(1,2,2)*qp_a(2,2)+&
                                     tij_abc(1,2,3)*qp_a(2,3)+&
                                     tij_abc(1,3,1)*qp_a(3,1)+&
                                     tij_abc(1,3,2)*qp_a(3,2)+&
                                     tij_abc(1,3,3)*qp_a(3,3))+&
                            dp_b(2)*(tij_abc(2,1,1)*qp_a(1,1)+&
                                     tij_abc(2,1,2)*qp_a(1,2)+&
                                     tij_abc(2,1,3)*qp_a(1,3)+&
                                     tij_abc(2,2,1)*qp_a(2,1)+&
                                     tij_abc(2,2,2)*qp_a(2,2)+&
                                     tij_abc(2,2,3)*qp_a(2,3)+&
                                     tij_abc(2,3,1)*qp_a(3,1)+&
                                     tij_abc(2,3,2)*qp_a(3,2)+&
                                     tij_abc(2,3,3)*qp_a(3,3))+&
                            dp_b(3)*(tij_abc(3,1,1)*qp_a(1,1)+&
                                     tij_abc(3,1,2)*qp_a(1,2)+&
                                     tij_abc(3,1,3)*qp_a(1,3)+&
                                     tij_abc(3,2,1)*qp_a(2,1)+&
                                     tij_abc(3,2,2)*qp_a(2,2)+&
                                     tij_abc(3,2,3)*qp_a(2,3)+&
                                     tij_abc(3,3,1)*qp_a(3,1)+&
                                     tij_abc(3,3,2)*qp_a(3,2)+&
                                     tij_abc(3,3,3)*qp_a(3,3))
                   tmp= fac * (tmp_ab - tmp_ba)
                   loc_energy = loc_energy + tmp
                   IF (do_forces) THEN
                      DO k = 1, 3
                         ! Dipole A - Quadrupole B
                         tmp_ab = dp_a(1)*(tij_abcd(k,1,1,1)*qp_b(1,1)+& 
                                           tij_abcd(k,1,1,2)*qp_b(1,2)+& 
                                           tij_abcd(k,1,1,3)*qp_b(1,3)+& 
                                           tij_abcd(k,1,2,1)*qp_b(2,1)+& 
                                           tij_abcd(k,1,2,2)*qp_b(2,2)+& 
                                           tij_abcd(k,1,2,3)*qp_b(2,3)+& 
                                           tij_abcd(k,1,3,1)*qp_b(3,1)+& 
                                           tij_abcd(k,1,3,2)*qp_b(3,2)+& 
                                           tij_abcd(k,1,3,3)*qp_b(3,3))+&
                                  dp_a(2)*(tij_abcd(k,2,1,1)*qp_b(1,1)+& 
                                           tij_abcd(k,2,1,2)*qp_b(1,2)+& 
                                           tij_abcd(k,2,1,3)*qp_b(1,3)+& 
                                           tij_abcd(k,2,2,1)*qp_b(2,1)+& 
                                           tij_abcd(k,2,2,2)*qp_b(2,2)+& 
                                           tij_abcd(k,2,2,3)*qp_b(2,3)+& 
                                           tij_abcd(k,2,3,1)*qp_b(3,1)+& 
                                           tij_abcd(k,2,3,2)*qp_b(3,2)+& 
                                           tij_abcd(k,2,3,3)*qp_b(3,3))+&
                                  dp_a(3)*(tij_abcd(k,3,1,1)*qp_b(1,1)+& 
                                           tij_abcd(k,3,1,2)*qp_b(1,2)+& 
                                           tij_abcd(k,3,1,3)*qp_b(1,3)+& 
                                           tij_abcd(k,3,2,1)*qp_b(2,1)+& 
                                           tij_abcd(k,3,2,2)*qp_b(2,2)+& 
                                           tij_abcd(k,3,2,3)*qp_b(2,3)+& 
                                           tij_abcd(k,3,3,1)*qp_b(3,1)+& 
                                           tij_abcd(k,3,3,2)*qp_b(3,2)+& 
                                           tij_abcd(k,3,3,3)*qp_b(3,3))

                         ! Dipole B - Quadrupole A   
                         tmp_ba = dp_b(1)*(tij_abcd(k,1,1,1)*qp_a(1,1)+& 
                                           tij_abcd(k,1,1,2)*qp_a(1,2)+& 
                                           tij_abcd(k,1,1,3)*qp_a(1,3)+& 
                                           tij_abcd(k,1,2,1)*qp_a(2,1)+& 
                                           tij_abcd(k,1,2,2)*qp_a(2,2)+& 
                                           tij_abcd(k,1,2,3)*qp_a(2,3)+& 
                                           tij_abcd(k,1,3,1)*qp_a(3,1)+& 
                                           tij_abcd(k,1,3,2)*qp_a(3,2)+& 
                                           tij_abcd(k,1,3,3)*qp_a(3,3))+&
                                  dp_b(2)*(tij_abcd(k,2,1,1)*qp_a(1,1)+& 
                                           tij_abcd(k,2,1,2)*qp_a(1,2)+& 
                                           tij_abcd(k,2,1,3)*qp_a(1,3)+& 
                                           tij_abcd(k,2,2,1)*qp_a(2,1)+& 
                                           tij_abcd(k,2,2,2)*qp_a(2,2)+& 
                                           tij_abcd(k,2,2,3)*qp_a(2,3)+& 
                                           tij_abcd(k,2,3,1)*qp_a(3,1)+& 
                                           tij_abcd(k,2,3,2)*qp_a(3,2)+& 
                                           tij_abcd(k,2,3,3)*qp_a(3,3))+&
                                  dp_b(3)*(tij_abcd(k,3,1,1)*qp_a(1,1)+& 
                                           tij_abcd(k,3,1,2)*qp_a(1,2)+& 
                                           tij_abcd(k,3,1,3)*qp_a(1,3)+& 
                                           tij_abcd(k,3,2,1)*qp_a(2,1)+& 
                                           tij_abcd(k,3,2,2)*qp_a(2,2)+& 
                                           tij_abcd(k,3,2,3)*qp_a(2,3)+& 
                                           tij_abcd(k,3,3,1)*qp_a(3,1)+& 
                                           tij_abcd(k,3,3,2)*qp_a(3,2)+& 
                                           tij_abcd(k,3,3,3)*qp_a(3,3))

                         fr(k) = fr(k) - fac * (tmp_ab - tmp_ba)
                      END DO
                   END IF
                END IF
                IF (task(3,1)) THEN
                   ! Quadrupole - Charge
                   fac = 1.0_dp/3.0_dp
                                ! Quadrupole B - Charge A 
                   tmp= ch_a * (tij_ab(1,1)*qp_b(1,1)+&
                                tij_ab(1,2)*qp_b(1,2)+&
                                tij_ab(1,3)*qp_b(1,3)+&
                                tij_ab(2,1)*qp_b(2,1)+&
                                tij_ab(2,2)*qp_b(2,2)+&
                                tij_ab(2,3)*qp_b(2,3)+&
                                tij_ab(3,1)*qp_b(3,1)+&
                                tij_ab(3,2)*qp_b(3,2)+&
                                tij_ab(3,3)*qp_b(3,3))+&
                                ! Quadrupole A - Charge B
                        ch_b * (tij_ab(1,1)*qp_a(1,1)+&
                                tij_ab(1,2)*qp_a(1,2)+&
                                tij_ab(1,3)*qp_a(1,3)+&
                                tij_ab(2,1)*qp_a(2,1)+&
                                tij_ab(2,2)*qp_a(2,2)+&
                                tij_ab(2,3)*qp_a(2,3)+&
                                tij_ab(3,1)*qp_a(3,1)+&
                                tij_ab(3,2)*qp_a(3,2)+&
                                tij_ab(3,3)*qp_a(3,3))     
                   loc_energy = loc_energy + fac*tmp
                   IF (do_forces) THEN
                      DO k = 1, 3
                                                       ! Quadrupole B - Charge A 
                         fr(k) = fr(k) - fac *(ch_a * (tij_abc(k,1,1)*qp_b(1,1)+&    
                                                       tij_abc(k,1,2)*qp_b(1,2)+&    
                                                       tij_abc(k,1,3)*qp_b(1,3)+&    
                                                       tij_abc(k,2,1)*qp_b(2,1)+&    
                                                       tij_abc(k,2,2)*qp_b(2,2)+&    
                                                       tij_abc(k,2,3)*qp_b(2,3)+&    
                                                       tij_abc(k,3,1)*qp_b(3,1)+&    
                                                       tij_abc(k,3,2)*qp_b(3,2)+&    
                                                       tij_abc(k,3,3)*qp_b(3,3))+&    
                                                       ! Quadrupole A - Charge B  
                                               ch_b * (tij_abc(k,1,1)*qp_a(1,1)+&    
                                                       tij_abc(k,1,2)*qp_a(1,2)+&    
                                                       tij_abc(k,1,3)*qp_a(1,3)+&    
                                                       tij_abc(k,2,1)*qp_a(2,1)+&    
                                                       tij_abc(k,2,2)*qp_a(2,2)+&    
                                                       tij_abc(k,2,3)*qp_a(2,3)+&    
                                                       tij_abc(k,3,1)*qp_a(3,1)+&    
                                                       tij_abc(k,3,2)*qp_a(3,2)+&    
                                                       tij_abc(k,3,3)*qp_a(3,3)))
                      END DO
                   END IF
                END IF
                energy = energy + loc_energy
                IF (do_forces) THEN
                   forces(1,atom_a) = forces(1,atom_a) - fr(1)
                   forces(2,atom_a) = forces(2,atom_a) - fr(2)
                   forces(3,atom_a) = forces(3,atom_a) - fr(3)
                   forces(1,atom_b) = forces(1,atom_b) + fr(1)
                   forces(2,atom_b) = forces(2,atom_b) + fr(2)
                   forces(3,atom_b) = forces(3,atom_b) + fr(3)
                END IF
                IF (do_efield) THEN
                   IF (do_efield1) THEN
                      efield1(1,atom_a) = efield1(1,atom_a) + ef_a(1)
                      efield1(2,atom_a) = efield1(2,atom_a) + ef_a(2)
                      efield1(3,atom_a) = efield1(3,atom_a) + ef_a(3)

                      efield1(1,atom_b) = efield1(1,atom_b) + ef_b(1)
                      efield1(2,atom_b) = efield1(2,atom_b) + ef_b(2)
                      efield1(3,atom_b) = efield1(3,atom_b) + ef_b(3)
                   END IF
                   IF (do_efield2) THEN
                      efield2(1,atom_a) = efield2(1,atom_a) + ef2_a(1,1)
                      efield2(2,atom_a) = efield2(2,atom_a) + ef2_a(1,2)
                      efield2(3,atom_a) = efield2(3,atom_a) + ef2_a(1,3)
                      efield2(4,atom_a) = efield2(4,atom_a) + ef2_a(2,1)
                      efield2(5,atom_a) = efield2(5,atom_a) + ef2_a(2,2)
                      efield2(6,atom_a) = efield2(6,atom_a) + ef2_a(2,3)
                      efield2(7,atom_a) = efield2(7,atom_a) + ef2_a(3,1)
                      efield2(8,atom_a) = efield2(8,atom_a) + ef2_a(3,2)
                      efield2(9,atom_a) = efield2(9,atom_a) + ef2_a(3,3)

                      efield2(1,atom_b) = efield2(1,atom_b) + ef2_b(1,1)
                      efield2(2,atom_b) = efield2(2,atom_b) + ef2_b(1,2)
                      efield2(3,atom_b) = efield2(3,atom_b) + ef2_b(1,3)
                      efield2(4,atom_b) = efield2(4,atom_b) + ef2_b(2,1)
                      efield2(5,atom_b) = efield2(5,atom_b) + ef2_b(2,2)
                      efield2(6,atom_b) = efield2(6,atom_b) + ef2_b(2,3)
                      efield2(7,atom_b) = efield2(7,atom_b) + ef2_b(3,1)
                      efield2(8,atom_b) = efield2(8,atom_b) + ef2_b(3,2)
                      efield2(9,atom_b) = efield2(9,atom_b) + ef2_b(3,3)
                   END IF
                END IF
                IF (do_stress) THEN
                   ptens11 = ptens11 + rab(1) * fr(1)
                   ptens21 = ptens21 + rab(2) * fr(1)
                   ptens31 = ptens31 + rab(3) * fr(1)
                   ptens12 = ptens12 + rab(1) * fr(2)
                   ptens22 = ptens22 + rab(2) * fr(2)
                   ptens32 = ptens32 + rab(3) * fr(2)
                   ptens13 = ptens13 + rab(1) * fr(3)
                   ptens23 = ptens23 + rab(2) * fr(3)
                   ptens33 = ptens33 + rab(3) * fr(3)
                END IF
             END IF
          END DO Pairs
       END DO Kind_Group_Loop
    END DO Lists
    IF (do_stress) THEN
       pv(1,1) = pv(1,1) + ptens11
       pv(1,2) = pv(1,2) + ptens12
       pv(1,3) = pv(1,3) + ptens13
       pv(2,1) = pv(2,1) + ptens21
       pv(2,2) = pv(2,2) + ptens22
       pv(2,3) = pv(2,3) + ptens23
       pv(3,1) = pv(3,1) + ptens31
       pv(3,2) = pv(3,2) + ptens32
       pv(3,3) = pv(3,3) + ptens33
    END IF

    CALL timestop ( handle )
  END SUBROUTINE ewald_multipole_SR

! *****************************************************************************
!> \brief computes the potential and the force for a lattice sum of multipoles
!>      up to quadrupole - Long Range (Reciprocal Space) Term
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE ewald_multipole_LR(ewald_env, ewald_pw, cell, particle_set, &
       local_particles, energy, task, do_forces, do_efield, do_stress, &
       charges, dipoles, quadrupoles, forces, pv, efield1, efield2, error)
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    LOGICAL, DIMENSION(3, 3), INTENT(IN)     :: task
    LOGICAL, INTENT(IN)                      :: do_forces, do_efield, &
                                                do_stress
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: charges
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN), OPTIONAL                   :: quadrupoles
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: forces, pv, efield1, efield2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_LR', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: atm_factor, cnjg_fac, &
                                                summe_tmp
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: summe_ef
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: summe_f1, summe_f2, summe_st
    INTEGER :: gpt, group, handle, i, icc, iparticle, iparticle_kind, &
      iparticle_local, j, k, lp, mp, nnodes, node, np, nparticle_kind, &
      nparticle_local, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    LOGICAL                                  :: do_efield1, do_efield2, &
                                                failure
    REAL(KIND=dp)                            :: alpha, denom, e_igdotr, eps0, &
                                                factor, flops, four_alpha_sq, &
                                                gauss, pref, tmp
    REAL(KIND=dp), DIMENSION(3)              :: vec
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: rho0
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(dg_type), POINTER                   :: dg
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(structure_factor_type)              :: exp_igr

    CALL timeset ( routineN, 'E', 'Mflops', handle )
    CALL cite_reference(Aguado2003)
    flops = 0.0_dp
    do_efield1 = do_efield.AND.PRESENT(efield1)
    do_efield2 = do_efield.AND.PRESENT(efield2)

    ! Gathering data from the ewald environment
    CALL ewald_env_get (ewald_env, alpha=alpha, eps0 = eps0, group = group ,error=error)
    CALL ewald_pw_get (ewald_pw, pw_big_pool=pw_pool, dg=dg)
    CALL dg_get (dg, dg_rho0=dg_rho0)
    rho0    => dg_rho0%density%pw%cr3d
    pw_grid => pw_pool%pw_grid
    bds     => pw_grid%bounds

    ! Allocation of working arrays
    nparticle_kind = SIZE(local_particles%n_el)
    nnodes = 0
    DO iparticle_kind = 1, nparticle_kind
       nnodes = nnodes + local_particles%n_el(iparticle_kind)
    ENDDO
    CALL structure_factor_allocate ( pw_grid%bounds, nnodes, exp_igr)

    ALLOCATE (summe_ef(1:pw_grid % ngpts_cut), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    summe_ef = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
    ! Stress Tensor
    IF (do_stress) THEN
       ALLOCATE (summe_st(9,1:pw_grid % ngpts_cut), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       summe_st = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
    END IF
    ! Electric Fields 
    IF (do_efield1) THEN
       ALLOCATE (summe_f1(3,1:pw_grid % ngpts_cut), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       summe_f1 = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
    END IF
    IF (do_efield2) THEN
       ALLOCATE (summe_f2(3,1:pw_grid % ngpts_cut), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       summe_f2 = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
    END IF

    ! defining four_alpha_sq
    four_alpha_sq = 4.0_dp * alpha ** 2
     ! zero node count
    node = 0
    DO iparticle_kind = 1, nparticle_kind
       nparticle_local = local_particles%n_el(iparticle_kind)
       DO iparticle_local=1,nparticle_local
          node = node + 1
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          CALL matvec_3x3 (vec, cell%h_inv, particle_set(iparticle)%r )
          CALL structure_factor_evaluate ( vec, pw_grid%npts, exp_igr%lb, &
               exp_igr%ex(:,node), exp_igr%ey(:,node), exp_igr%ez(:,node) )
       END DO
    END DO

    ! Looping over the positive g-vectors
    DO gpt = 1, pw_grid % ngpts_cut
       lp = pw_grid % mapl % pos ( pw_grid % g_hat ( 1, gpt ) )
       mp = pw_grid % mapm % pos ( pw_grid % g_hat ( 2, gpt ) )
       np = pw_grid % mapn % pos ( pw_grid % g_hat ( 3, gpt ) )

       lp = lp + bds (1,1)
       mp = mp + bds (1,2)
       np = np + bds (1,3)

       ! Initializing sum to be used in the energy and force
       node = 0
       DO iparticle_kind = 1, nparticle_kind
          nparticle_local = local_particles%n_el(iparticle_kind)
          DO iparticle_local=1,nparticle_local
             node = node + 1
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             CALL get_atom_factor(atm_factor, pw_grid, gpt, iparticle, task, charges,&
                  dipoles, quadrupoles, error)
             summe_tmp     = exp_igr%ex(lp,node)*exp_igr%ey(mp,node)*exp_igr%ez(np,node)
             summe_ef(gpt) = summe_ef(gpt) + atm_factor *summe_tmp
          END DO
       END DO

       ! Precompute pseudo-density for stress tensor calculation
       IF (do_stress) THEN
          icc = 0
          DO i = 1,3
             DO j =1,3
                icc = icc + 1
                DO iparticle_kind = 1, nparticle_kind
                   nparticle_local = local_particles%n_el(iparticle_kind)
                   DO iparticle_local=1,nparticle_local
                      node = node + 1
                      iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                      atm_factor = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
                      IF (task(2,2)) THEN
                         ! Dipole - Dipole
                         tmp = 0.5_dp*(dipoles(i,iparticle)*pw_grid%g(j,gpt)+&
                                       dipoles(j,iparticle)*pw_grid%g(i,gpt))
                         atm_factor = atm_factor - 2.0_dp * tmp * CMPLX(0.0_dp, -1.0_dp, KIND=dp)
                      END IF
                      IF (task(3,3)) THEN
                         ! Quadrupole - Quadrupole
                         tmp = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
                         DO k = 1,3
                            tmp = tmp + quadrupoles(j,k,iparticle)*pw_grid%g(i,gpt)*pw_grid%g(k,gpt)+&
                                        quadrupoles(i,k,iparticle)*pw_grid%g(j,gpt)*pw_grid%g(k,gpt)
                         END DO
                         tmp        = 0.5_dp * tmp
                         atm_factor = atm_factor - 4.0_dp/3.0_dp * tmp
                      END IF
                      summe_tmp = exp_igr%ex(lp,node)*exp_igr%ey(mp,node)*exp_igr%ez(np,node)
                      summe_st(icc,gpt) = summe_st(icc,gpt) + atm_factor *summe_tmp
                   END DO
                END DO
             END DO
          END DO
       END IF
       ! Precompute pseudo-density for stress tensor calculation
       IF (do_efield) THEN
          IF (do_efield1) THEN
             icc = 0
             DO i = 1,3
                icc = icc + 1
                DO iparticle_kind = 1, nparticle_kind
                   nparticle_local = local_particles%n_el(iparticle_kind)
                   DO iparticle_local=1,nparticle_local
                      node = node + 1
                      iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                      IF (task(2,2)) THEN
                         ! Dipole - Dipole
                         atm_factor = -  pw_grid%g(i,gpt) * CMPLX(0.0_dp, -1.0_dp, KIND=dp)
                      END IF
                      summe_tmp = exp_igr%ex(lp,node)*exp_igr%ey(mp,node)*exp_igr%ez(np,node)
                      summe_f1(icc,gpt) = summe_f1(icc,gpt) + atm_factor*summe_tmp
                   END DO
                END DO
             END DO
          END IF
          IF (do_efield2) THEN
             icc = 0
             DO i = 1,3
                DO j =1,3
                   icc = icc + 1
                   DO iparticle_kind = 1, nparticle_kind
                      nparticle_local = local_particles%n_el(iparticle_kind)
                      DO iparticle_local=1,nparticle_local
                         node = node + 1
                         iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                         IF (task(3,3)) THEN
                            ! Quadrupole - Quadrupole
                            atm_factor = - 1.0_dp/3.0_dp * pw_grid%g(j,gpt)*pw_grid%g(i,gpt)
                         END IF
                         summe_tmp = exp_igr%ex(lp,node)*exp_igr%ey(mp,node)*exp_igr%ez(np,node)
                         summe_f2(icc,gpt) = summe_f2(icc,gpt) + atm_factor*summe_tmp
                      END DO
                   END DO
                END DO
             END DO
          END IF
       END IF
    END DO
    flops = flops + 13.0_dp*(pw_grid%ngpts_cut*10.0_dp&
                  +          pw_grid%ngpts_cut*node*20.0_dp)

    CALL mp_sum ( summe_ef, group )
    IF (do_stress)  CALL mp_sum ( summe_st, group )
    IF (do_efield1) CALL mp_sum ( summe_f1, group )
    IF (do_efield2) CALL mp_sum ( summe_f2, group )

    ! looping over the positive g-vectors
    DO gpt = 1, pw_grid % ngpts_cut
       ! computing the potential energy
       lp = pw_grid%mapl%pos(pw_grid%g_hat(1,gpt))
       mp = pw_grid%mapm%pos(pw_grid%g_hat(2,gpt))
       np = pw_grid%mapn%pos(pw_grid%g_hat(3,gpt))

       lp = lp + bds(1,1)
       mp = mp + bds(1,2)
       np = np + bds(1,3)

       IF ( pw_grid % gsq ( gpt ) <= 1.0E-10_dp ) CYCLE

       gauss  = (rho0(lp,mp,np) * pw_grid%vol)**2 / pw_grid%gsq(gpt)
       factor = gauss * REAL(summe_ef(gpt) * CONJG(summe_ef(gpt)))
       energy = energy + factor

       ! Compute the force
       IF (do_forces) THEN
          node = 0
          DO iparticle_kind = 1, nparticle_kind
             nparticle_local = local_particles%n_el(iparticle_kind)
             DO iparticle_local=1,nparticle_local
                node = node + 1
                iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                CALL get_atom_factor(atm_factor, pw_grid, gpt, iparticle, task, charges,&
                     dipoles, quadrupoles, error)
                cnjg_fac = CONJG(exp_igr%ex(lp,node)*exp_igr%ey(mp,node)*exp_igr%ez(np,node))
                e_igdotr = AIMAG ( summe_ef(gpt) * cnjg_fac)
                forces(:,node) = forces(:,node) + atm_factor * gauss * e_igdotr * pw_grid%g(:,gpt)

                ! Electric field
                IF (do_efield) THEN
                   IF (do_efield1) THEN
                      efield1(1,iparticle) = efield1(1,iparticle) - 2.0_dp * gauss * REAL(summe_f1(1,gpt) * cnjg_fac)
                      efield1(2,iparticle) = efield1(2,iparticle) - 2.0_dp * gauss * REAL(summe_f1(2,gpt) * cnjg_fac)
                      efield1(3,iparticle) = efield1(3,iparticle) - 2.0_dp * gauss * REAL(summe_f1(3,gpt) * cnjg_fac)
                   END IF
                   IF (do_efield2) THEN
                      efield2(1,iparticle) = efield2(1,iparticle) - 2.0_dp * gauss * REAL(summe_f2(1,gpt) * cnjg_fac)
                      efield2(2,iparticle) = efield2(2,iparticle) - 2.0_dp * gauss * REAL(summe_f2(2,gpt) * cnjg_fac)
                      efield2(3,iparticle) = efield2(3,iparticle) - 2.0_dp * gauss * REAL(summe_f2(3,gpt) * cnjg_fac)
                      efield2(1,iparticle) = efield2(1,iparticle) - 2.0_dp * gauss * REAL(summe_f2(4,gpt) * cnjg_fac)
                      efield2(2,iparticle) = efield2(2,iparticle) - 2.0_dp * gauss * REAL(summe_f2(5,gpt) * cnjg_fac)
                      efield2(3,iparticle) = efield2(3,iparticle) - 2.0_dp * gauss * REAL(summe_f2(6,gpt) * cnjg_fac)
                      efield2(1,iparticle) = efield2(1,iparticle) - 2.0_dp * gauss * REAL(summe_f2(7,gpt) * cnjg_fac)
                      efield2(2,iparticle) = efield2(2,iparticle) - 2.0_dp * gauss * REAL(summe_f2(8,gpt) * cnjg_fac)
                      efield2(3,iparticle) = efield2(3,iparticle) - 2.0_dp * gauss * REAL(summe_f2(9,gpt) * cnjg_fac)
                   END IF
                END IF

             END DO
          END DO
       END IF

       ! Compute the virial P*V
       IF (do_stress) THEN
          ! The Stress Tensor can be decomposed in two main components.
          ! The first one is just a normal ewald component for reciprocal space
          denom = 1.0_dp / four_alpha_sq + 1.0_dp / pw_grid % gsq ( gpt )
          pv(1,1) = pv(1,1) + factor*(1.0_dp - 2.0_dp*pw_grid%g(1,gpt)*pw_grid%g(1,gpt)*denom)
          pv(1,2) = pv(1,2) - factor*(         2.0_dp*pw_grid%g(1,gpt)*pw_grid%g(2,gpt)*denom)
          pv(1,3) = pv(1,3) - factor*(         2.0_dp*pw_grid%g(1,gpt)*pw_grid%g(3,gpt)*denom)
          pv(2,1) = pv(2,1) - factor*(         2.0_dp*pw_grid%g(2,gpt)*pw_grid%g(1,gpt)*denom)
          pv(2,2) = pv(2,2) + factor*(1.0_dp - 2.0_dp*pw_grid%g(2,gpt)*pw_grid%g(2,gpt)*denom)
          pv(2,3) = pv(2,3) - factor*(         2.0_dp*pw_grid%g(2,gpt)*pw_grid%g(3,gpt)*denom)
          pv(3,1) = pv(3,1) - factor*(         2.0_dp*pw_grid%g(3,gpt)*pw_grid%g(1,gpt)*denom)
          pv(3,2) = pv(3,2) - factor*(         2.0_dp*pw_grid%g(3,gpt)*pw_grid%g(2,gpt)*denom)
          pv(3,3) = pv(3,3) + factor*(1.0_dp - 2.0_dp*pw_grid%g(3,gpt)*pw_grid%g(3,gpt)*denom)
          ! The second one can be written in the following way
          pv(1,1) = pv(1,1) + gauss * REAL(summe_st(1,gpt) * CONJG(summe_ef(gpt)))
          pv(1,2) = pv(1,2) + gauss * REAL(summe_st(2,gpt) * CONJG(summe_ef(gpt)))
          pv(1,3) = pv(1,3) + gauss * REAL(summe_st(3,gpt) * CONJG(summe_ef(gpt)))
          pv(2,1) = pv(2,1) + gauss * REAL(summe_st(4,gpt) * CONJG(summe_ef(gpt)))
          pv(2,2) = pv(2,2) + gauss * REAL(summe_st(5,gpt) * CONJG(summe_ef(gpt)))
          pv(2,3) = pv(2,3) + gauss * REAL(summe_st(6,gpt) * CONJG(summe_ef(gpt)))
          pv(3,1) = pv(3,1) + gauss * REAL(summe_st(7,gpt) * CONJG(summe_ef(gpt)))
          pv(3,2) = pv(3,2) + gauss * REAL(summe_st(8,gpt) * CONJG(summe_ef(gpt)))
          pv(3,3) = pv(3,3) + gauss * REAL(summe_st(9,gpt) * CONJG(summe_ef(gpt)))
       END IF
       
    END DO

    flops = flops + pw_grid % ngpts_cut * 55.0_dp &
                  + pw_grid % ngpts_cut * node * 23.0_dp

    pref   = 1.0_dp / eps0 / pw_grid % vol
    energy = energy * pref

    CALL structure_factor_deallocate (exp_igr)
    DEALLOCATE (summe_ef, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (do_forces) THEN
       forces = 2.0_dp * forces * pref
    END IF
    IF (do_stress) THEN
       pv = pv * pref
       DEALLOCATE (summe_st, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (do_efield1) THEN
       efield1 = efield1 * pref
       DEALLOCATE (summe_f1, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (do_efield2) THEN
       efield2 = efield2 * pref
       DEALLOCATE (summe_f2, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    flops = flops * 1.0E-6_dp
    CALL timestop ( flops, handle )

  END SUBROUTINE ewald_multipole_LR

! *****************************************************************************
!> \brief Computes the atom factor including charge, dipole and quadrupole
!> \par History
!>      none
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich 
! *****************************************************************************
  SUBROUTINE get_atom_factor(atm_factor, pw_grid, gpt, iparticle, task, charges,&
       dipoles, quadrupoles, error)
    COMPLEX(KIND=dp), INTENT(OUT)            :: atm_factor
    TYPE(pw_grid_type), POINTER              :: pw_grid
    INTEGER, INTENT(IN)                      :: gpt
    INTEGER                                  :: iparticle
    LOGICAL, DIMENSION(3, 3), INTENT(IN)     :: task
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: charges
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN), OPTIONAL                   :: quadrupoles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_atom_factor', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: tmp
    INTEGER                                  :: i, j

    atm_factor = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
    IF (task(1,1)) THEN
       ! Charge - Charge
       atm_factor = atm_factor + charges(iparticle)
    END IF
    IF (task(2,2)) THEN
       ! Dipole - Dipole
       tmp = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
       DO i = 1,3
          tmp = tmp + dipoles(i,iparticle)*pw_grid%g(i,gpt)
       END DO
       atm_factor = atm_factor - tmp * CMPLX(0.0_dp, -1.0_dp, KIND=dp)
    END IF
    IF (task(3,3)) THEN
       ! Quadrupole - Quadrupole
       tmp = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
       DO i = 1,3
          DO j = 1,3
             tmp = tmp + quadrupoles(j,i,iparticle)*pw_grid%g(j,gpt)*pw_grid%g(i,gpt)
          END DO
       END DO
       atm_factor = atm_factor - 1.0_dp/3.0_dp * tmp
    END IF

  END SUBROUTINE get_atom_factor

! *****************************************************************************
!> \brief Computes the self interaction from g-space and the neutralizing background
!>     when using multipoles
!> \author Teodoro Laino [tlaino] - University of Zurich - 12.2007
! *****************************************************************************
  SUBROUTINE ewald_multipole_self (ewald_env, local_particles, e_self, e_neut,&
       task, do_efield, charges, dipoles, quadrupoles, efield1, efield2, error)

    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(OUT)               :: e_self, e_neut
    LOGICAL, DIMENSION(3, 3), INTENT(IN)     :: task
    LOGICAL, INTENT(IN)                      :: do_efield
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: charges
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN), OPTIONAL                   :: quadrupoles
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: efield1, efield2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_self', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: f23 = 2.0_dp/3.0_dp, &
                                                f415 = 4.0_dp/15.0_dp

    INTEGER                                  :: ewald_type, group, i, &
                                                iparticle, iparticle_kind, &
                                                iparticle_local, j, &
                                                nparticle_local
    LOGICAL                                  :: do_efield1, do_efield2
    REAL(KIND=dp)                            :: alpha, ch_qu_self, &
                                                dipole_self, eps0, &
                                                invfoureps0, q, q_self, &
                                                q_sum, qu_qu_self

    CALL ewald_env_get ( ewald_env, ewald_type=ewald_type, eps0=eps0, alpha=alpha,&
         group=group, error=error)

    do_efield1 = do_efield.AND.PRESENT(efield1)
    do_efield2 = do_efield.AND.PRESENT(efield2)
    q_self      = 0.0_dp
    q_sum       = 0.0_dp
    dipole_self = 0.0_dp
    ch_qu_self  = 0.0_dp
    qu_qu_self  = 0.0_dp
    invfoureps0 = 1.0_dp / (4.0_dp *  pi * eps0)
    DO iparticle_kind=1,SIZE(local_particles%n_el)
       nparticle_local = local_particles%n_el(iparticle_kind)
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          IF (task(1,1)) THEN
             ! Charge - Charge
             q = charges(iparticle)
             q_self = q_self + q*q
             q_sum  = q_sum  + q
          END IF
          IF (task(2,2)) THEN
             ! Dipole - Dipole
             DO i = 1,3
                dipole_self = dipole_self + dipoles(i,iparticle)**2
             END DO
             IF (do_efield1) THEN
                efield1(1,iparticle) = efield1(1,iparticle) - (2.0_dp/SQRT(pi))*f23*alpha**3*dipoles(1,iparticle)
                efield1(2,iparticle) = efield1(2,iparticle) - (2.0_dp/SQRT(pi))*f23*alpha**3*dipoles(2,iparticle)
                efield1(3,iparticle) = efield1(3,iparticle) - (2.0_dp/SQRT(pi))*f23*alpha**3*dipoles(3,iparticle)
             END IF
          END IF
          IF (task(1,3)) THEN
             DO i = 1,3
                ch_qu_self = ch_qu_self + quadrupoles(i,i,iparticle) * charges(iparticle)
             END DO
             IF (do_efield1) THEN
                efield1(1,iparticle) = efield1(1,iparticle) + (2.0_dp/SQRT(pi))*f23*alpha**3*charges(iparticle)
                efield1(5,iparticle) = efield1(5,iparticle) + (2.0_dp/SQRT(pi))*f23*alpha**3*charges(iparticle)
                efield1(9,iparticle) = efield1(9,iparticle) + (2.0_dp/SQRT(pi))*f23*alpha**3*charges(iparticle)
             END IF
          END IF
          IF (task(3,3)) THEN
             ! Quadrupole - Quadrupole
             DO i = 1,3 
                DO j = 1,3
                   qu_qu_self = qu_qu_self + quadrupoles(j,i,iparticle)**2
                END DO
             END DO
             IF (do_efield2) THEN
                efield2(1,iparticle) = efield2(1,iparticle) - (4.0_dp/SQRT(pi))*f415*alpha**5*quadrupoles(1,1,iparticle)
                efield2(2,iparticle) = efield2(2,iparticle) - (4.0_dp/SQRT(pi))*f415*alpha**5*quadrupoles(2,1,iparticle)
                efield2(3,iparticle) = efield2(3,iparticle) - (4.0_dp/SQRT(pi))*f415*alpha**5*quadrupoles(3,1,iparticle)
                efield2(4,iparticle) = efield2(4,iparticle) - (4.0_dp/SQRT(pi))*f415*alpha**5*quadrupoles(1,2,iparticle)
                efield2(5,iparticle) = efield2(5,iparticle) - (4.0_dp/SQRT(pi))*f415*alpha**5*quadrupoles(2,2,iparticle)
                efield2(6,iparticle) = efield2(6,iparticle) - (4.0_dp/SQRT(pi))*f415*alpha**5*quadrupoles(3,2,iparticle)
                efield2(7,iparticle) = efield2(7,iparticle) - (4.0_dp/SQRT(pi))*f415*alpha**5*quadrupoles(1,3,iparticle)
                efield2(8,iparticle) = efield2(8,iparticle) - (4.0_dp/SQRT(pi))*f415*alpha**5*quadrupoles(2,3,iparticle)
                efield2(9,iparticle) = efield2(9,iparticle) - (4.0_dp/SQRT(pi))*f415*alpha**5*quadrupoles(3,3,iparticle)
             END IF
          END IF
       END DO
    END DO

    CALL mp_sum (q_self, group)
    CALL mp_sum (q_sum, group)
    CALL mp_sum (dipole_self, group)
    CALL mp_sum (ch_qu_self, group)
    CALL mp_sum (qu_qu_self, group)

    e_neut = 0.0_dp
    e_self = 0.0_dp
    IF ( ewald_type /= do_ewald_none) THEN
       e_self = -(q_self+f23*(dipole_self-f23*ch_qu_self+f415*qu_qu_self*alpha**2)*alpha**2)&
            *alpha*invfoureps0/SQRT(pi)
       e_neut = -q_sum**2 * pi * invfoureps0/(2.0_dp * alpha ** 2)
    END IF

  END SUBROUTINE ewald_multipole_self

! *****************************************************************************
!> \author Teodoro Laino [tlaino] - University of Zurich - 12.2007
! *****************************************************************************
  SUBROUTINE ewald_multipole_print ( iw, cell, e_gspace, e_self, e_neut, e_bonded )

    INTEGER, INTENT(IN)                      :: iw
    TYPE(cell_type), INTENT(IN)              :: cell
    REAL(KIND=dp), INTENT(IN)                :: e_gspace, e_self, e_neut, &
                                                e_bonded

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_print', &
      routineP = moduleN//':'//routineN

    TYPE(unit_convert_type), POINTER         :: units

    CALL get_cp2k_units ( units )

    WRITE ( iw, '( A, A )' ) ' *********************************', &
         '**********************************************'
    WRITE ( iw, '( A, A, T35, A, T66, E15.7 )' ) ' INITIAL GSPACE ENERGY', &
         units % e_label, '= ', e_gspace
    WRITE ( iw, '( A, A, T35, A, T66, E15.7 )' ) ' SELF ENERGY CORRECTION', &
         units % e_label, '= ', e_self
    WRITE ( iw, '( A, A, T35, A, T66, E15.7 )' ) ' NEUT. BACKGROUND', units % e_label, &
         '= ', e_neut / cell % deth
    WRITE ( iw, '( A, A, T35, A, T66, E15.7 )' ) ' BONDED CORRECTION', units % e_label, &
         '= ', e_bonded
    WRITE ( iw, '( A, A )' ) ' *********************************', &
         '**********************************************'

  END SUBROUTINE ewald_multipole_print

END MODULE ewalds_multipole
