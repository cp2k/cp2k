!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Treats the electrostatic for multipoles (up to quadrupoles)
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! *****************************************************************************
MODULE ewalds_multipole
  USE atomic_kind_types,               ONLY: atomic_kind_type
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
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             oorootpi,&
                                             pi
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
  PRIVATE
#include "ewalds_multipole_debug.h"

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  LOGICAL, PRIVATE, PARAMETER :: debug_r_space    =.FALSE.
  LOGICAL, PRIVATE, PARAMETER :: debug_g_space    =.FALSE.
  LOGICAL, PRIVATE, PARAMETER :: debug_e_field    =.FALSE.
  LOGICAL, PRIVATE, PARAMETER :: debug_e_field_en =.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewalds_multipole'

  PUBLIC :: ewald_multipole_evaluate

CONTAINS

! *****************************************************************************
!> \brief  Computes the potential and the force for a lattice sum of multipoles
!>      (up to quadrupole)
!> \par    Note 
!>         atomic_kind_set and mm_section are between the arguments only
!>         for debug purpose (therefore optional) and can be avoided when this
!>         function is called in other part of the program
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! *****************************************************************************
  RECURSIVE SUBROUTINE ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env,&
       cell, particle_set, local_particles, energy_local, energy_glob, e_neut, e_self,&
       task, do_forces, do_stress, do_efield, charges, dipoles, quadrupoles,&
       forces_local, forces_glob, pv_local, pv_glob, efield0, efield1, efield2,&
       iw, do_debug, atomic_kind_set, mm_section, error )
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), POINTER             :: particle_set(:)
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(INOUT)             :: energy_local, energy_glob
    REAL(KIND=dp), INTENT(OUT)               :: e_neut, e_self
    LOGICAL, DIMENSION(3), INTENT(IN)        :: task
    LOGICAL, INTENT(IN)                      :: do_forces, do_stress, &
                                                do_efield
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: quadrupoles
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: forces_local, forces_glob, &
                                                pv_local, pv_glob
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(OUT), OPTIONAL                  :: efield0
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: efield1, efield2
    INTEGER, INTENT(IN)                      :: iw
    LOGICAL, INTENT(IN)                      :: do_debug
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: mm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_evaluate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, size1, size2, &
                                                stat
    LOGICAL                                  :: check_debug, check_efield, &
                                                check_forces, do_task(3), &
                                                failure
    LOGICAL, DIMENSION(3, 3)                 :: my_task
    REAL(KIND=dp)                            :: e_rspace
    REAL(KIND=dp), DIMENSION(:), POINTER     :: efield0_lr, efield0_sr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: efield1_lr, efield1_sr, &
                                                efield2_lr, efield2_sr

    failure = .FALSE.
    CALL timeset(routineN,'I','',handle)
    CPPostcondition(ASSOCIATED(nonbond_env),cp_failure_level,routineP,error,failure)
    check_debug  = (debug_this_module.OR.debug_r_space.OR.debug_g_space.OR.debug_e_field.OR.debug_e_field_en)&
         .EQV.debug_this_module
    CPPostcondition(check_debug,cp_failure_level,routineP,error,failure)
    check_forces = do_forces.EQV.(PRESENT(forces_local).AND.PRESENT(forces_glob))
    CPPostcondition(check_forces,cp_failure_level,routineP,error,failure)
    check_efield = do_efield.EQV.(PRESENT(efield0).OR.PRESENT(efield1).OR.PRESENT(efield2))
    CPPostcondition(check_efield,cp_failure_level,routineP,error,failure)
    ! Debugging this module
    IF (debug_this_module.AND.do_debug) THEN
       ! Debug specifically real space part
       IF (debug_r_space) THEN
          CALL debug_ewald_multipoles(ewald_env, ewald_pw, nonbond_env, cell, &
               particle_set, local_particles, iw, debug_r_space, error)
          STOP "Debug Multipole Requested:  Real Part!"
       END IF
       ! Debug electric fields and gradients as pure derivatives
       IF (debug_e_field) THEN
          CPPostcondition(PRESENT(atomic_kind_set),cp_failure_level,routineP,error,failure)
          CPPostcondition(PRESENT(mm_section),cp_failure_level,routineP,error,failure)
          CALL debug_ewald_multipoles_fields(ewald_env, ewald_pw, nonbond_env,&
               cell, particle_set, local_particles, charges, dipoles, quadrupoles,&
               task, iw, atomic_kind_set, mm_section, error)
          STOP "Debug Multipole Requested:  POT+EFIELDS+GRAD!"
       END IF
       ! Debug the potential, electric fields and electric fields gradient in oder
       ! to retrieve the correct energy
       IF (debug_e_field_en) THEN
          CALL debug_ewald_multipoles_fields2(ewald_env, ewald_pw, nonbond_env,&
               cell, particle_set, local_particles, charges, dipoles, quadrupoles,&
               task, iw, error)
          STOP "Debug Multipole Requested:  POT+EFIELDS+GRAD to give the correct energy!!" 
       END IF
    END IF

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

    ! Allocate arrays for the evaluation of the potential, fields and electrostatic field gradients
    NULLIFY(efield0_sr, efield0_lr, efield1_sr, efield1_lr, efield2_sr, efield2_lr)
    IF (do_efield) THEN
       IF (PRESENT(efield0)) THEN
          size1 = SIZE(efield0)
          ALLOCATE (efield0_sr(size1), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (efield0_lr(size1), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          efield0_sr = 0.0_dp
          efield0_lr = 0.0_dp
       END IF
       IF (PRESENT(efield1)) THEN
          size1 = SIZE(efield1,1)
          size2 = SIZE(efield1,2)
          ALLOCATE (efield1_sr(size1,size2), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (efield1_lr(size1,size2), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          efield1_sr = 0.0_dp
          efield1_lr = 0.0_dp
       END IF
       IF (PRESENT(efield2)) THEN
          size1 = SIZE(efield2,1)
          size2 = SIZE(efield2,2)
          ALLOCATE (efield2_sr(size1,size2), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (efield2_lr(size1,size2), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          efield2_sr = 0.0_dp
          efield2_lr = 0.0_dp
       END IF
    END IF

    e_rspace = 0.0_dp
    IF (.NOT.debug_g_space) THEN
       ! Compute the Real Space (Short-Range) part of the Ewald sum
       CALL ewald_multipole_SR (nonbond_env, ewald_env, cell, e_rspace, my_task,&
            do_forces, do_efield, do_stress, charges, dipoles, quadrupoles, forces_glob,&
            pv_glob, efield0_sr, efield1_sr, efield2_sr, error)
       energy_glob = energy_glob + e_rspace
    END IF

    e_neut = 0.0_dp
    e_self = 0.0_dp
    IF (.NOT.debug_r_space) THEN
       ! Compute the Reciprocal Space (Long-Range) part of the Ewald sum
       CALL ewald_multipole_LR(ewald_env, ewald_pw, cell, particle_set, &
            local_particles, energy_local, my_task, do_forces, do_efield, do_stress,&
            charges, dipoles, quadrupoles, forces_local, pv_local, efield0_lr, efield1_lr,&
            efield2_lr, error)

       ! Self-Interactions corrections
       CALL ewald_multipole_self (ewald_env, cell, local_particles, e_self, e_neut, my_task,&
            do_efield, charges, dipoles, quadrupoles, efield0_lr, efield1_lr, efield2_lr, error)
    END IF

    ! Print some info about energetics
    CALL ewald_multipole_print (iw, energy_local, energy_glob, e_rspace, e_self, e_neut)

    ! Gather the components of the potential, fields and electrostatic field gradients
    IF (do_efield) THEN
       IF (PRESENT(efield0)) THEN
          efield0 = efield0_sr + efield0_lr
          DEALLOCATE (efield0_sr, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (efield0_lr, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(efield1)) THEN
          efield1 = efield1_sr + efield1_lr
          DEALLOCATE (efield1_sr, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (efield1_lr, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(efield2)) THEN
          efield2 = efield2_sr + efield2_lr
          DEALLOCATE (efield2_sr, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (efield2_lr, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    CALL timestop(0.0_dp, handle)
  END SUBROUTINE ewald_multipole_evaluate
  
! *****************************************************************************
!> \brief computes the potential and the force for a lattice sum of multipoles
!>      up to quadrupole - Short Range (Real Space) Term
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE ewald_multipole_SR (nonbond_env, ewald_env, cell, energy, task,&
       do_forces, do_efield, do_stress, charges, dipoles, quadrupoles, forces,&
       pv, efield0, efield1, efield2, error)
    TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    LOGICAL, DIMENSION(3, 3), INTENT(IN)     :: task
    LOGICAL, INTENT(IN)                      :: do_forces, do_efield, &
                                                do_stress
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: quadrupoles
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: forces, pv
    REAL(KIND=dp), DIMENSION(:), POINTER     :: efield0
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: efield1, efield2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_SR', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a, atom_a, atom_b, b, c, d, &
                                                e, handle, i, iend, igrp, &
                                                ikind, ilist, ipair, istart, &
                                                jkind, k, nkinds, npairs
    INTEGER, DIMENSION(:, :), POINTER        :: list
    LOGICAL                                  :: do_efield0, do_efield1, &
                                                do_efield2, failure, &
                                                force_eval
    REAL(KIND=dp) :: alpha, ch_i, ch_j, ef0_i, ef0_j, eloc, fac, fac_ij, ir, &
      irab2, ptens11, ptens12, ptens13, ptens21, ptens22, ptens23, ptens31, &
      ptens32, ptens33, r, rab2, rab2_max, rcut, tij, tmp, tmp1, tmp11, &
      tmp12, tmp13, tmp2, tmp21, tmp22, tmp23, tmp31, tmp32, tmp33, tmp_ij, &
      tmp_ji
    REAL(KIND=dp), DIMENSION(0:5)            :: f
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, cvi, dp_i, dp_j, &
                                                ef1_i, ef1_j, fr, rab, tij_a
    REAL(KIND=dp), DIMENSION(3, 3)           :: ef2_i, ef2_j, qp_i, qp_j, &
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
    CALL cite_reference(Aguado2003)
    NULLIFY(nonbonded,r_last_update, r_last_update_pbc)
    do_efield0 = do_efield.AND.ASSOCIATED(efield0)
    do_efield1 = do_efield.AND.ASSOCIATED(efield1)
    do_efield2 = do_efield.AND.ASSOCIATED(efield2)
    IF (do_stress) THEN 
       ptens11 = 0.0_dp ; ptens12 = 0.0_dp ; ptens13 = 0.0_dp
       ptens21 = 0.0_dp ; ptens22 = 0.0_dp ; ptens23 = 0.0_dp
       ptens31 = 0.0_dp ; ptens32 = 0.0_dp ; ptens33 = 0.0_dp
    END IF
    ! Get nonbond_env info
    CALL fist_nonbond_env_get (nonbond_env, nonbonded=nonbonded, natom_types = nkinds,&
         r_last_update=r_last_update,r_last_update_pbc=r_last_update_pbc, error=error)
    CALL ewald_env_get (ewald_env, alpha=alpha, rcut=rcut, error=error)
    rab2_max = rcut**2
    IF (debug_r_space) THEN
       rab2_max = HUGE(0.0_dp)
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
             fac_ij = 1.0_dp
             atom_a = list(1,ipair)
             atom_b = list(2,ipair)
             IF (atom_a==atom_b) fac_ij = 0.5_dp
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
                
                r     = SQRT(rab2)
                irab2 = 1.0_dp/rab2
                ir    = 1.0_dp/r
                ! Compute the radial function
                IF (debug_this_module.AND.debug_r_space.AND.(.NOT.debug_g_space)) THEN
                   f(0)  = ir
                   tmp   = 0.0_dp
                ELSE
                   f(0)  = erfc(alpha*r)*ir
                   tmp   = EXP(-alpha**2*rab2)*oorootpi
                END IF
                fac = 1.0_dp
                DO i = 1, 5
                   fac  = fac*REAL(2*i-1,KIND=dp)
                   f(i) = irab2*(f(i-1)+ tmp*((2.0_dp*alpha**2)**i)/(fac*alpha))
                END DO
                ! Compute the Tensor components
                force_eval = .FALSE.
                IF (task(1,1)) THEN
                   tij         = f(0)*fac_ij
                                                 force_eval = do_forces .OR.do_efield1
                END IF
                IF (task(2,2))                   force_eval = force_eval.OR.do_efield0  
                IF (task(1,2).OR.force_eval) THEN
                   force_eval = .FALSE.
                   tij_a    = - rab*f(1)*fac_ij
                   IF (task(1,2))                force_eval = do_forces
                END IF
                IF (task(1,1))                   force_eval = force_eval.OR.do_efield2
                IF (task(3,3))                   force_eval = force_eval.OR.do_efield0
                IF (task(2,2).OR.task(3,1).OR.force_eval) THEN
                   force_eval = .FALSE.
                   DO b = 1,3
                      DO a = 1,3
                         tmp = rab(a)*rab(b)*fac_ij
                         tij_ab(a,b) = 3.0_dp*tmp*f(2)
                         IF (a==b) tij_ab(a,b) = tij_ab(a,b) - f(1)*fac_ij
                      END DO
                   END DO
                   IF (task(2,2).OR.task(3,1))   force_eval = do_forces
                END IF
                IF (task(2,2))                   force_eval = force_eval.OR.do_efield2
                IF (task(3,3))                   force_eval = force_eval.OR.do_efield1
                IF (task(3,2).OR.force_eval) THEN
                   force_eval = .FALSE.
                   DO c = 1, 3
                      DO b = 1, 3
                         DO a = 1, 3
                            tmp = rab(a)*rab(b)*rab(c)*fac_ij
                            tij_abc(a,b,c) = - 15.0_dp*tmp*f(3)
                            tmp = 3.0_dp*f(2)*fac_ij
                            IF (a==b) tij_abc(a,b,c) = tij_abc(a,b,c) + tmp*rab(c)
                            IF (a==c) tij_abc(a,b,c) = tij_abc(a,b,c) + tmp*rab(b)
                            IF (b==c) tij_abc(a,b,c) = tij_abc(a,b,c) + tmp*rab(a)
                         END DO
                      END DO
                   END DO
                   IF (task(3,2))                force_eval = do_forces
                END IF
                IF (task(3,3).OR.force_eval) THEN
                   force_eval = .FALSE.
                   DO d = 1, 3
                      DO c = 1, 3
                         DO b = 1, 3
                            DO a = 1, 3
                               tmp = rab(a)*rab(b)*rab(c)*rab(d)*fac_ij
                               tij_abcd(a,b,c,d) = 105.0_dp*tmp*f(4)
                               tmp1 = 15.0_dp*f(3)*fac_ij
                               tmp2 =  3.0_dp*f(2)*fac_ij
                               IF (a==b) THEN
                                            tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(c)*rab(d)
                                  IF (c==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) + tmp2
                               END IF
                               IF (a==c) THEN 
                                            tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(b)*rab(d)
                                  IF (b==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) + tmp2
                               END IF
                               IF (a==d)    tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(b)*rab(c)
                               IF (b==c) THEN
                                            tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(a)*rab(d)
                                  IF (a==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) + tmp2
                               END IF
                               IF (b==d)    tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(a)*rab(c)
                               IF (c==d)    tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(a)*rab(b)
                            END DO
                         END DO
                      END DO
                   END DO
                   IF (task(3,3))                force_eval = do_forces
                END IF
                IF (force_eval) THEN
                   force_eval = .FALSE.
                   DO e = 1, 3
                      DO d = 1, 3
                         DO c = 1, 3
                            DO b = 1, 3
                               DO a = 1, 3
                                  tmp = rab(a)*rab(b)*rab(c)*rab(d)*rab(e)*fac_ij
                                  tij_abcde(a,b,c,d,e) = -945.0_dp*tmp*f(5)
                                  tmp1 = 105.0_dp*f(4)*fac_ij
                                  tmp2 =  15.0_dp*f(3)*fac_ij
                                  IF (a==b) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(c)*rab(d)*rab(e)
                                     IF (c==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(e)
                                     IF (c==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(d)
                                     IF (d==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(c)
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
                                  IF (a==e) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(b)*rab(c)*rab(d)
                                     IF (b==c) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(d)
                                     IF (b==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(c)
                                     IF (c==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(b)
                                  END IF
                                  IF (b==c) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(d)*rab(e)
                                     IF (d==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(a)
                                  END IF
                                  IF (b==d) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(c)*rab(e)
                                     IF (c==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(a)
                                  END IF
                                  IF (b==e) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(c)*rab(d)
                                     IF (c==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(a)
                                  END IF
                                  IF (c==d)    tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(b)*rab(e)
                                  IF (c==e)    tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(b)*rab(d)
                                  IF (d==e)    tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(b)*rab(c)
                               END DO
                            END DO
                         END DO
                      END DO
                   END DO
                END IF
                eloc  = 0.0_dp
                fr    = 0.0_dp
                ef0_i = 0.0_dp
                ef0_j = 0.0_dp
                ef1_j = 0.0_dp
                ef1_i = 0.0_dp
                ef2_j = 0.0_dp
                ef2_i = 0.0_dp
                ch_j  = HUGE(0.0_dp)
                ch_i  = HUGE(0.0_dp)
                dp_j  = HUGE(0.0_dp)
                dp_i  = HUGE(0.0_dp)
                qp_j  = HUGE(0.0_dp)
                qp_i  = HUGE(0.0_dp)
                IF (ANY(task(1,:))) THEN
                   ch_j  = charges(atom_a)
                   ch_i  = charges(atom_b)
                END IF
                IF (ANY(task(2,:))) THEN
                   dp_j  = dipoles(:,atom_a)
                   dp_i  = dipoles(:,atom_b)
                END IF
                IF (ANY(task(3,:))) THEN
                   qp_j  = quadrupoles(:,:,atom_a)
                   qp_i  = quadrupoles(:,:,atom_b)
                END IF
                IF (task(1,1)) THEN
                   ! Charge - Charge
                   eloc = eloc + ch_i*tij*ch_j
                   ! Forces on particle i (locally b)
                   IF (do_forces) THEN
                      fr(1) = fr(1) - ch_j * tij_a(1) * ch_i
                      fr(2) = fr(2) - ch_j * tij_a(2) * ch_i
                      fr(3) = fr(3) - ch_j * tij_a(3) * ch_i
                   END IF
                   ! Electric fields
                   IF (do_efield) THEN
                      ! Potential
                      IF (do_efield0) THEN
                         ef0_i = ef0_i + tij * ch_j 

                         ef0_j = ef0_j + tij * ch_i
                      END IF
                      ! Electric field
                      IF (do_efield1) THEN
                         ef1_i(1) = ef1_i(1) - tij_a(1) * ch_j
                         ef1_i(2) = ef1_i(2) - tij_a(2) * ch_j
                         ef1_i(3) = ef1_i(3) - tij_a(3) * ch_j

                         ef1_j(1) = ef1_j(1) + tij_a(1) * ch_i
                         ef1_j(2) = ef1_j(2) + tij_a(2) * ch_i
                         ef1_j(3) = ef1_j(3) + tij_a(3) * ch_i
                      END IF
                      ! Electric field gradient
                      IF (do_efield2) THEN
                         ef2_i(1,1) = ef2_i(1,1) - tij_ab(1,1) * ch_j
                         ef2_i(2,1) = ef2_i(2,1) - tij_ab(2,1) * ch_j
                         ef2_i(3,1) = ef2_i(3,1) - tij_ab(3,1) * ch_j
                         ef2_i(1,2) = ef2_i(1,2) - tij_ab(1,2) * ch_j
                         ef2_i(2,2) = ef2_i(2,2) - tij_ab(2,2) * ch_j
                         ef2_i(3,2) = ef2_i(3,2) - tij_ab(3,2) * ch_j
                         ef2_i(1,3) = ef2_i(1,3) - tij_ab(1,3) * ch_j
                         ef2_i(2,3) = ef2_i(2,3) - tij_ab(2,3) * ch_j
                         ef2_i(3,3) = ef2_i(3,3) - tij_ab(3,3) * ch_j

                         ef2_j(1,1) = ef2_j(1,1) - tij_ab(1,1) * ch_i
                         ef2_j(2,1) = ef2_j(2,1) - tij_ab(2,1) * ch_i
                         ef2_j(3,1) = ef2_j(3,1) - tij_ab(3,1) * ch_i
                         ef2_j(1,2) = ef2_j(1,2) - tij_ab(1,2) * ch_i
                         ef2_j(2,2) = ef2_j(2,2) - tij_ab(2,2) * ch_i
                         ef2_j(3,2) = ef2_j(3,2) - tij_ab(3,2) * ch_i
                         ef2_j(1,3) = ef2_j(1,3) - tij_ab(1,3) * ch_i
                         ef2_j(2,3) = ef2_j(2,3) - tij_ab(2,3) * ch_i
                         ef2_j(3,3) = ef2_j(3,3) - tij_ab(3,3) * ch_i
                      END IF
                   END IF
                END IF
                IF (task(2,2)) THEN
                   ! Dipole - Dipole
                   tmp= - (dp_i(1)*(tij_ab(1,1)*dp_j(1)+&
                                    tij_ab(2,1)*dp_j(2)+&
                                    tij_ab(3,1)*dp_j(3))+&
                           dp_i(2)*(tij_ab(1,2)*dp_j(1)+&
                                    tij_ab(2,2)*dp_j(2)+&
                                    tij_ab(3,2)*dp_j(3))+&
                           dp_i(3)*(tij_ab(1,3)*dp_j(1)+&
                                    tij_ab(2,3)*dp_j(2)+&
                                    tij_ab(3,3)*dp_j(3)))
                   eloc = eloc + tmp
                   ! Forces on particle i (locally b)
                   IF (do_forces) THEN
                      DO k = 1, 3
                         fr(k) = fr(k) +  dp_i(1)*(tij_abc(1,1,k)*dp_j(1)+&
                                                   tij_abc(2,1,k)*dp_j(2)+&
                                                   tij_abc(3,1,k)*dp_j(3))&
                                       +  dp_i(2)*(tij_abc(1,2,k)*dp_j(1)+&
                                                   tij_abc(2,2,k)*dp_j(2)+&
                                                   tij_abc(3,2,k)*dp_j(3))&
                                       +  dp_i(3)*(tij_abc(1,3,k)*dp_j(1)+&
                                                   tij_abc(2,3,k)*dp_j(2)+&
                                                   tij_abc(3,3,k)*dp_j(3))
                      END DO
                   END IF
                   ! Electric fields
                   IF (do_efield) THEN
                      ! Potential
                      IF (do_efield0) THEN
                         ef0_i = ef0_i - (tij_a(1)*dp_j(1)+&
                                          tij_a(2)*dp_j(2)+&
                                          tij_a(3)*dp_j(3))

                         ef0_j = ef0_j + (tij_a(1)*dp_i(1)+&
                                          tij_a(2)*dp_i(2)+&
                                          tij_a(3)*dp_i(3))
                      END IF
                      ! Electric field
                      IF (do_efield1) THEN
                         ef1_i(1) = ef1_i(1) + (tij_ab(1,1)*dp_j(1)+&
                                                tij_ab(2,1)*dp_j(2)+&
                                                tij_ab(3,1)*dp_j(3))
                         ef1_i(2) = ef1_i(2) + (tij_ab(1,2)*dp_j(1)+&
                                                tij_ab(2,2)*dp_j(2)+&
                                                tij_ab(3,2)*dp_j(3))
                         ef1_i(3) = ef1_i(3) + (tij_ab(1,3)*dp_j(1)+&
                                                tij_ab(2,3)*dp_j(2)+&
                                                tij_ab(3,3)*dp_j(3))

                         ef1_j(1) = ef1_j(1) + (tij_ab(1,1)*dp_i(1)+&
                                                tij_ab(2,1)*dp_i(2)+&
                                                tij_ab(3,1)*dp_i(3)) 
                         ef1_j(2) = ef1_j(2) + (tij_ab(1,2)*dp_i(1)+&
                                                tij_ab(2,2)*dp_i(2)+&
                                                tij_ab(3,2)*dp_i(3))
                         ef1_j(3) = ef1_j(3) + (tij_ab(1,3)*dp_i(1)+&
                                                tij_ab(2,3)*dp_i(2)+&
                                                tij_ab(3,3)*dp_i(3))
                      END IF
                      ! Electric field gradient
                      IF (do_efield2) THEN
                         ef2_i(1,1) = ef2_i(1,1) + (tij_abc(1,1,1)*dp_j(1)+&
                                                    tij_abc(2,1,1)*dp_j(2)+&
                                                    tij_abc(3,1,1)*dp_j(3)) 
                         ef2_i(1,2) = ef2_i(1,2) + (tij_abc(1,1,2)*dp_j(1)+&
                                                    tij_abc(2,1,2)*dp_j(2)+&
                                                    tij_abc(3,1,2)*dp_j(3)) 
                         ef2_i(1,3) = ef2_i(1,3) + (tij_abc(1,1,3)*dp_j(1)+&
                                                    tij_abc(2,1,3)*dp_j(2)+&
                                                    tij_abc(3,1,3)*dp_j(3)) 
                         ef2_i(2,1) = ef2_i(2,1) + (tij_abc(1,2,1)*dp_j(1)+&
                                                    tij_abc(2,2,1)*dp_j(2)+&
                                                    tij_abc(3,2,1)*dp_j(3)) 
                         ef2_i(2,2) = ef2_i(2,2) + (tij_abc(1,2,2)*dp_j(1)+&
                                                    tij_abc(2,2,2)*dp_j(2)+&
                                                    tij_abc(3,2,2)*dp_j(3)) 
                         ef2_i(2,3) = ef2_i(2,3) + (tij_abc(1,2,3)*dp_j(1)+&
                                                    tij_abc(2,2,3)*dp_j(2)+&
                                                    tij_abc(3,2,3)*dp_j(3)) 
                         ef2_i(3,1) = ef2_i(3,1) + (tij_abc(1,3,1)*dp_j(1)+&
                                                    tij_abc(2,3,1)*dp_j(2)+&
                                                    tij_abc(3,3,1)*dp_j(3)) 
                         ef2_i(3,2) = ef2_i(3,2) + (tij_abc(1,3,2)*dp_j(1)+&
                                                    tij_abc(2,3,2)*dp_j(2)+&
                                                    tij_abc(3,3,2)*dp_j(3)) 
                         ef2_i(3,3) = ef2_i(3,3) + (tij_abc(1,3,3)*dp_j(1)+&
                                                    tij_abc(2,3,3)*dp_j(2)+&
                                                    tij_abc(3,3,3)*dp_j(3))

                         ef2_j(1,1) = ef2_j(1,1) - (tij_abc(1,1,1)*dp_i(1)+&
                                                    tij_abc(2,1,1)*dp_i(2)+&
                                                    tij_abc(3,1,1)*dp_i(3)) 
                         ef2_j(1,2) = ef2_j(1,2) - (tij_abc(1,1,2)*dp_i(1)+&
                                                    tij_abc(2,1,2)*dp_i(2)+&
                                                    tij_abc(3,1,2)*dp_i(3)) 
                         ef2_j(1,3) = ef2_j(1,3) - (tij_abc(1,1,3)*dp_i(1)+&
                                                    tij_abc(2,1,3)*dp_i(2)+&
                                                    tij_abc(3,1,3)*dp_i(3)) 
                         ef2_j(2,1) = ef2_j(2,1) - (tij_abc(1,2,1)*dp_i(1)+&
                                                    tij_abc(2,2,1)*dp_i(2)+&
                                                    tij_abc(3,2,1)*dp_i(3)) 
                         ef2_j(2,2) = ef2_j(2,2) - (tij_abc(1,2,2)*dp_i(1)+&
                                                    tij_abc(2,2,2)*dp_i(2)+&
                                                    tij_abc(3,2,2)*dp_i(3)) 
                         ef2_j(2,3) = ef2_j(2,3) - (tij_abc(1,2,3)*dp_i(1)+&
                                                    tij_abc(2,2,3)*dp_i(2)+&
                                                    tij_abc(3,2,3)*dp_i(3)) 
                         ef2_j(3,1) = ef2_j(3,1) - (tij_abc(1,3,1)*dp_i(1)+&
                                                    tij_abc(2,3,1)*dp_i(2)+&
                                                    tij_abc(3,3,1)*dp_i(3)) 
                         ef2_j(3,2) = ef2_j(3,2) - (tij_abc(1,3,2)*dp_i(1)+&
                                                    tij_abc(2,3,2)*dp_i(2)+&
                                                    tij_abc(3,3,2)*dp_i(3)) 
                         ef2_j(3,3) = ef2_j(3,3) - (tij_abc(1,3,3)*dp_i(1)+&
                                                    tij_abc(2,3,3)*dp_i(2)+&
                                                    tij_abc(3,3,3)*dp_i(3))
                      END IF
                   END IF
                END IF
                IF (task(2,1)) THEN
                   ! Dipole - Charge
                   tmp=   ch_j*(tij_a(1)*dp_i(1)+&
                                tij_a(2)*dp_i(2)+&
                                tij_a(3)*dp_i(3))&
                        - ch_i*(tij_a(1)*dp_j(1)+&
                                tij_a(2)*dp_j(2)+&
                                tij_a(3)*dp_j(3))
                   eloc = eloc + tmp
                   ! Forces on particle i (locally b)
                   IF (do_forces) THEN
                      DO k = 1, 3
                         fr(k) = fr(k) -  ch_j *(tij_ab(1,k)*dp_i(1)+&
                                                 tij_ab(2,k)*dp_i(2)+&
                                                 tij_ab(3,k)*dp_i(3))&
                                       +  ch_i *(tij_ab(1,k)*dp_j(1)+&
                                                 tij_ab(2,k)*dp_j(2)+&
                                                 tij_ab(3,k)*dp_j(3))
                      END DO                      
                   END IF
                END IF
                IF (task(3,3)) THEN
                   ! Quadrupole - Quadrupole
                   fac  = 1.0_dp/9.0_dp
                   tmp11 = qp_i(1,1)*(tij_abcd(1,1,1,1)*qp_j(1,1)+&
                                      tij_abcd(2,1,1,1)*qp_j(2,1)+&
                                      tij_abcd(3,1,1,1)*qp_j(3,1)+&
                                      tij_abcd(1,2,1,1)*qp_j(1,2)+&
                                      tij_abcd(2,2,1,1)*qp_j(2,2)+&
                                      tij_abcd(3,2,1,1)*qp_j(3,2)+&
                                      tij_abcd(1,3,1,1)*qp_j(1,3)+&
                                      tij_abcd(2,3,1,1)*qp_j(2,3)+&
                                      tij_abcd(3,3,1,1)*qp_j(3,3)) 
                   tmp21 = qp_i(2,1)*(tij_abcd(1,1,1,2)*qp_j(1,1)+&
                                      tij_abcd(2,1,1,2)*qp_j(2,1)+&
                                      tij_abcd(3,1,1,2)*qp_j(3,1)+&
                                      tij_abcd(1,2,1,2)*qp_j(1,2)+&
                                      tij_abcd(2,2,1,2)*qp_j(2,2)+&
                                      tij_abcd(3,2,1,2)*qp_j(3,2)+&
                                      tij_abcd(1,3,1,2)*qp_j(1,3)+&
                                      tij_abcd(2,3,1,2)*qp_j(2,3)+&
                                      tij_abcd(3,3,1,2)*qp_j(3,3))
                   tmp31 = qp_i(3,1)*(tij_abcd(1,1,1,3)*qp_j(1,1)+&
                                      tij_abcd(2,1,1,3)*qp_j(2,1)+&
                                      tij_abcd(3,1,1,3)*qp_j(3,1)+&
                                      tij_abcd(1,2,1,3)*qp_j(1,2)+&
                                      tij_abcd(2,2,1,3)*qp_j(2,2)+&
                                      tij_abcd(3,2,1,3)*qp_j(3,2)+&
                                      tij_abcd(1,3,1,3)*qp_j(1,3)+&
                                      tij_abcd(2,3,1,3)*qp_j(2,3)+&
                                      tij_abcd(3,3,1,3)*qp_j(3,3))
                   tmp22 = qp_i(2,2)*(tij_abcd(1,1,2,2)*qp_j(1,1)+&
                                      tij_abcd(2,1,2,2)*qp_j(2,1)+&
                                      tij_abcd(3,1,2,2)*qp_j(3,1)+&
                                      tij_abcd(1,2,2,2)*qp_j(1,2)+&
                                      tij_abcd(2,2,2,2)*qp_j(2,2)+&
                                      tij_abcd(3,2,2,2)*qp_j(3,2)+&
                                      tij_abcd(1,3,2,2)*qp_j(1,3)+&
                                      tij_abcd(2,3,2,2)*qp_j(2,3)+&
                                      tij_abcd(3,3,2,2)*qp_j(3,3))
                   tmp32 = qp_i(3,2)*(tij_abcd(1,1,2,3)*qp_j(1,1)+&
                                      tij_abcd(2,1,2,3)*qp_j(2,1)+&
                                      tij_abcd(3,1,2,3)*qp_j(3,1)+&
                                      tij_abcd(1,2,2,3)*qp_j(1,2)+&
                                      tij_abcd(2,2,2,3)*qp_j(2,2)+&
                                      tij_abcd(3,2,2,3)*qp_j(3,2)+&
                                      tij_abcd(1,3,2,3)*qp_j(1,3)+&
                                      tij_abcd(2,3,2,3)*qp_j(2,3)+&
                                      tij_abcd(3,3,2,3)*qp_j(3,3))
                   tmp33 = qp_i(3,3)*(tij_abcd(1,1,3,3)*qp_j(1,1)+&
                                      tij_abcd(2,1,3,3)*qp_j(2,1)+&
                                      tij_abcd(3,1,3,3)*qp_j(3,1)+&
                                      tij_abcd(1,2,3,3)*qp_j(1,2)+&
                                      tij_abcd(2,2,3,3)*qp_j(2,2)+&
                                      tij_abcd(3,2,3,3)*qp_j(3,2)+&
                                      tij_abcd(1,3,3,3)*qp_j(1,3)+&
                                      tij_abcd(2,3,3,3)*qp_j(2,3)+&
                                      tij_abcd(3,3,3,3)*qp_j(3,3))
                   tmp12 = tmp21
                   tmp13 = tmp31
                   tmp23 = tmp32
                   tmp   = tmp11 + tmp12 + tmp13 + &
                           tmp21 + tmp22 + tmp23 + &
                           tmp31 + tmp32 + tmp33

                   eloc = eloc + fac*tmp
                   ! Forces on particle i (locally b)
                   IF (do_forces) THEN
                      DO k = 1, 3
                         tmp11 = qp_i(1,1)*(tij_abcde(1,1,1,1,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,1,1,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,1,1,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,1,1,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,1,1,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,1,1,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,1,1,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,1,1,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,1,1,k)*qp_j(3,3))
                         tmp21 = qp_i(2,1)*(tij_abcde(1,1,2,1,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,2,1,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,2,1,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,2,1,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,2,1,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,2,1,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,2,1,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,2,1,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,2,1,k)*qp_j(3,3))
                         tmp31 = qp_i(3,1)*(tij_abcde(1,1,3,1,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,3,1,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,3,1,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,3,1,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,3,1,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,3,1,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,3,1,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,3,1,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,3,1,k)*qp_j(3,3))
                         tmp22 = qp_i(2,2)*(tij_abcde(1,1,2,2,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,2,2,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,2,2,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,2,2,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,2,2,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,2,2,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,2,2,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,2,2,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,2,2,k)*qp_j(3,3))
                         tmp32 = qp_i(3,2)*(tij_abcde(1,1,3,2,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,3,2,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,3,2,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,3,2,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,3,2,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,3,2,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,3,2,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,3,2,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,3,2,k)*qp_j(3,3))
                         tmp33 = qp_i(3,3)*(tij_abcde(1,1,3,3,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,3,3,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,3,3,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,3,3,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,3,3,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,3,3,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,3,3,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,3,3,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,3,3,k)*qp_j(3,3))
                         tmp12 = tmp21
                         tmp13 = tmp31
                         tmp23 = tmp32
                         fr(k) = fr(k) - fac * ( tmp11 + tmp12 + tmp13 +&
                                                 tmp21 + tmp22 + tmp23 +&
                                                 tmp31 + tmp32 + tmp33  )
                      END DO
                   END IF
                   ! Electric field
                   IF (do_efield) THEN
                      fac = 1.0_dp/3.0_dp
                      ! Potential
                      IF (do_efield0) THEN
                         ef0_i = ef0_i + fac*(tij_ab(1,1)*qp_j(1,1)+&
                                              tij_ab(2,1)*qp_j(2,1)+&
                                              tij_ab(3,1)*qp_j(3,1)+&
                                              tij_ab(1,2)*qp_j(1,2)+&
                                              tij_ab(2,2)*qp_j(2,2)+&
                                              tij_ab(3,2)*qp_j(3,2)+&
                                              tij_ab(1,3)*qp_j(1,3)+&
                                              tij_ab(2,3)*qp_j(2,3)+&
                                              tij_ab(3,3)*qp_j(3,3)) 

                         ef0_j = ef0_j + fac*(tij_ab(1,1)*qp_i(1,1)+&
                                              tij_ab(2,1)*qp_i(2,1)+&
                                              tij_ab(3,1)*qp_i(3,1)+&
                                              tij_ab(1,2)*qp_i(1,2)+&
                                              tij_ab(2,2)*qp_i(2,2)+&
                                              tij_ab(3,2)*qp_i(3,2)+&
                                              tij_ab(1,3)*qp_i(1,3)+&
                                              tij_ab(2,3)*qp_i(2,3)+&
                                              tij_ab(3,3)*qp_i(3,3))
                      END IF
                      ! Electric field
                      IF (do_efield1) THEN
                         ef1_i(1) = ef1_i(1) - fac*(tij_abc(1,1,1)*qp_j(1,1)+&
                                                    tij_abc(2,1,1)*qp_j(2,1)+&
                                                    tij_abc(3,1,1)*qp_j(3,1)+&
                                                    tij_abc(1,2,1)*qp_j(1,2)+&
                                                    tij_abc(2,2,1)*qp_j(2,2)+&
                                                    tij_abc(3,2,1)*qp_j(3,2)+&
                                                    tij_abc(1,3,1)*qp_j(1,3)+&
                                                    tij_abc(2,3,1)*qp_j(2,3)+&
                                                    tij_abc(3,3,1)*qp_j(3,3)) 
                         ef1_i(2) = ef1_i(2) - fac*(tij_abc(1,1,2)*qp_j(1,1)+&
                                                    tij_abc(2,1,2)*qp_j(2,1)+&
                                                    tij_abc(3,1,2)*qp_j(3,1)+&
                                                    tij_abc(1,2,2)*qp_j(1,2)+&
                                                    tij_abc(2,2,2)*qp_j(2,2)+&
                                                    tij_abc(3,2,2)*qp_j(3,2)+&
                                                    tij_abc(1,3,2)*qp_j(1,3)+&
                                                    tij_abc(2,3,2)*qp_j(2,3)+&
                                                    tij_abc(3,3,2)*qp_j(3,3)) 
                         ef1_i(3) = ef1_i(3) - fac*(tij_abc(1,1,3)*qp_j(1,1)+&
                                                    tij_abc(2,1,3)*qp_j(2,1)+&
                                                    tij_abc(3,1,3)*qp_j(3,1)+&
                                                    tij_abc(1,2,3)*qp_j(1,2)+&
                                                    tij_abc(2,2,3)*qp_j(2,2)+&
                                                    tij_abc(3,2,3)*qp_j(3,2)+&
                                                    tij_abc(1,3,3)*qp_j(1,3)+&
                                                    tij_abc(2,3,3)*qp_j(2,3)+&
                                                    tij_abc(3,3,3)*qp_j(3,3))

                         ef1_j(1) = ef1_j(1) + fac*(tij_abc(1,1,1)*qp_i(1,1)+&
                                                    tij_abc(2,1,1)*qp_i(2,1)+&
                                                    tij_abc(3,1,1)*qp_i(3,1)+&
                                                    tij_abc(1,2,1)*qp_i(1,2)+&
                                                    tij_abc(2,2,1)*qp_i(2,2)+&
                                                    tij_abc(3,2,1)*qp_i(3,2)+&
                                                    tij_abc(1,3,1)*qp_i(1,3)+&
                                                    tij_abc(2,3,1)*qp_i(2,3)+&
                                                    tij_abc(3,3,1)*qp_i(3,3)) 
                         ef1_j(2) = ef1_j(2) + fac*(tij_abc(1,1,2)*qp_i(1,1)+&
                                                    tij_abc(2,1,2)*qp_i(2,1)+&
                                                    tij_abc(3,1,2)*qp_i(3,1)+&
                                                    tij_abc(1,2,2)*qp_i(1,2)+&
                                                    tij_abc(2,2,2)*qp_i(2,2)+&
                                                    tij_abc(3,2,2)*qp_i(3,2)+&
                                                    tij_abc(1,3,2)*qp_i(1,3)+&
                                                    tij_abc(2,3,2)*qp_i(2,3)+&
                                                    tij_abc(3,3,2)*qp_i(3,3)) 
                         ef1_j(3) = ef1_j(3) + fac*(tij_abc(1,1,3)*qp_i(1,1)+&
                                                    tij_abc(2,1,3)*qp_i(2,1)+&
                                                    tij_abc(3,1,3)*qp_i(3,1)+&
                                                    tij_abc(1,2,3)*qp_i(1,2)+&
                                                    tij_abc(2,2,3)*qp_i(2,2)+&
                                                    tij_abc(3,2,3)*qp_i(3,2)+&
                                                    tij_abc(1,3,3)*qp_i(1,3)+&
                                                    tij_abc(2,3,3)*qp_i(2,3)+&
                                                    tij_abc(3,3,3)*qp_i(3,3)) 
                      END IF
                      ! Electric field gradient
                      IF (do_efield2) THEN
                         tmp11 =   fac *(tij_abcd(1,1,1,1)*qp_j(1,1)+&
                                         tij_abcd(2,1,1,1)*qp_j(2,1)+&
                                         tij_abcd(3,1,1,1)*qp_j(3,1)+&
                                         tij_abcd(1,2,1,1)*qp_j(1,2)+&
                                         tij_abcd(2,2,1,1)*qp_j(2,2)+&
                                         tij_abcd(3,2,1,1)*qp_j(3,2)+&
                                         tij_abcd(1,3,1,1)*qp_j(1,3)+&
                                         tij_abcd(2,3,1,1)*qp_j(2,3)+&
                                         tij_abcd(3,3,1,1)*qp_j(3,3)) 
                         tmp12 =   fac *(tij_abcd(1,1,1,2)*qp_j(1,1)+&
                                         tij_abcd(2,1,1,2)*qp_j(2,1)+&
                                         tij_abcd(3,1,1,2)*qp_j(3,1)+&
                                         tij_abcd(1,2,1,2)*qp_j(1,2)+&
                                         tij_abcd(2,2,1,2)*qp_j(2,2)+&
                                         tij_abcd(3,2,1,2)*qp_j(3,2)+&
                                         tij_abcd(1,3,1,2)*qp_j(1,3)+&
                                         tij_abcd(2,3,1,2)*qp_j(2,3)+&
                                         tij_abcd(3,3,1,2)*qp_j(3,3)) 
                         tmp13 =   fac *(tij_abcd(1,1,1,3)*qp_j(1,1)+&
                                         tij_abcd(2,1,1,3)*qp_j(2,1)+&
                                         tij_abcd(3,1,1,3)*qp_j(3,1)+&
                                         tij_abcd(1,2,1,3)*qp_j(1,2)+&
                                         tij_abcd(2,2,1,3)*qp_j(2,2)+&
                                         tij_abcd(3,2,1,3)*qp_j(3,2)+&
                                         tij_abcd(1,3,1,3)*qp_j(1,3)+&
                                         tij_abcd(2,3,1,3)*qp_j(2,3)+&
                                         tij_abcd(3,3,1,3)*qp_j(3,3)) 
                         tmp22 =   fac *(tij_abcd(1,1,2,2)*qp_j(1,1)+&
                                         tij_abcd(2,1,2,2)*qp_j(2,1)+&
                                         tij_abcd(3,1,2,2)*qp_j(3,1)+&
                                         tij_abcd(1,2,2,2)*qp_j(1,2)+&
                                         tij_abcd(2,2,2,2)*qp_j(2,2)+&
                                         tij_abcd(3,2,2,2)*qp_j(3,2)+&
                                         tij_abcd(1,3,2,2)*qp_j(1,3)+&
                                         tij_abcd(2,3,2,2)*qp_j(2,3)+&
                                         tij_abcd(3,3,2,2)*qp_j(3,3)) 
                         tmp23 =   fac *(tij_abcd(1,1,2,3)*qp_j(1,1)+&
                                         tij_abcd(2,1,2,3)*qp_j(2,1)+&
                                         tij_abcd(3,1,2,3)*qp_j(3,1)+&
                                         tij_abcd(1,2,2,3)*qp_j(1,2)+&
                                         tij_abcd(2,2,2,3)*qp_j(2,2)+&
                                         tij_abcd(3,2,2,3)*qp_j(3,2)+&
                                         tij_abcd(1,3,2,3)*qp_j(1,3)+&
                                         tij_abcd(2,3,2,3)*qp_j(2,3)+&
                                         tij_abcd(3,3,2,3)*qp_j(3,3)) 
                         tmp33 =   fac *(tij_abcd(1,1,3,3)*qp_j(1,1)+&
                                         tij_abcd(2,1,3,3)*qp_j(2,1)+&
                                         tij_abcd(3,1,3,3)*qp_j(3,1)+&
                                         tij_abcd(1,2,3,3)*qp_j(1,2)+&
                                         tij_abcd(2,2,3,3)*qp_j(2,2)+&
                                         tij_abcd(3,2,3,3)*qp_j(3,2)+&
                                         tij_abcd(1,3,3,3)*qp_j(1,3)+&
                                         tij_abcd(2,3,3,3)*qp_j(2,3)+&
                                         tij_abcd(3,3,3,3)*qp_j(3,3)) 

                         ef2_i(1,1) = ef2_i(1,1) - tmp11
                         ef2_i(1,2) = ef2_i(1,2) - tmp12
                         ef2_i(1,3) = ef2_i(1,3) - tmp13
                         ef2_i(2,1) = ef2_i(2,1) - tmp12
                         ef2_i(2,2) = ef2_i(2,2) - tmp22
                         ef2_i(2,3) = ef2_i(2,3) - tmp23
                         ef2_i(3,1) = ef2_i(3,1) - tmp13
                         ef2_i(3,2) = ef2_i(3,2) - tmp23
                         ef2_i(3,3) = ef2_i(3,3) - tmp33

                         tmp11 =   fac *(tij_abcd(1,1,1,1)*qp_i(1,1)+&
                                         tij_abcd(2,1,1,1)*qp_i(2,1)+&
                                         tij_abcd(3,1,1,1)*qp_i(3,1)+&
                                         tij_abcd(1,2,1,1)*qp_i(1,2)+&
                                         tij_abcd(2,2,1,1)*qp_i(2,2)+&
                                         tij_abcd(3,2,1,1)*qp_i(3,2)+&
                                         tij_abcd(1,3,1,1)*qp_i(1,3)+&
                                         tij_abcd(2,3,1,1)*qp_i(2,3)+&
                                         tij_abcd(3,3,1,1)*qp_i(3,3)) 
                         tmp12 =   fac *(tij_abcd(1,1,1,2)*qp_i(1,1)+&
                                         tij_abcd(2,1,1,2)*qp_i(2,1)+&
                                         tij_abcd(3,1,1,2)*qp_i(3,1)+&
                                         tij_abcd(1,2,1,2)*qp_i(1,2)+&
                                         tij_abcd(2,2,1,2)*qp_i(2,2)+&
                                         tij_abcd(3,2,1,2)*qp_i(3,2)+&
                                         tij_abcd(1,3,1,2)*qp_i(1,3)+&
                                         tij_abcd(2,3,1,2)*qp_i(2,3)+&
                                         tij_abcd(3,3,1,2)*qp_i(3,3)) 
                         tmp13 =   fac *(tij_abcd(1,1,1,3)*qp_i(1,1)+&
                                         tij_abcd(2,1,1,3)*qp_i(2,1)+&
                                         tij_abcd(3,1,1,3)*qp_i(3,1)+&
                                         tij_abcd(1,2,1,3)*qp_i(1,2)+&
                                         tij_abcd(2,2,1,3)*qp_i(2,2)+&
                                         tij_abcd(3,2,1,3)*qp_i(3,2)+&
                                         tij_abcd(1,3,1,3)*qp_i(1,3)+&
                                         tij_abcd(2,3,1,3)*qp_i(2,3)+&
                                         tij_abcd(3,3,1,3)*qp_i(3,3)) 
                         tmp22 =   fac *(tij_abcd(1,1,2,2)*qp_i(1,1)+&
                                         tij_abcd(2,1,2,2)*qp_i(2,1)+&
                                         tij_abcd(3,1,2,2)*qp_i(3,1)+&
                                         tij_abcd(1,2,2,2)*qp_i(1,2)+&
                                         tij_abcd(2,2,2,2)*qp_i(2,2)+&
                                         tij_abcd(3,2,2,2)*qp_i(3,2)+&
                                         tij_abcd(1,3,2,2)*qp_i(1,3)+&
                                         tij_abcd(2,3,2,2)*qp_i(2,3)+&
                                         tij_abcd(3,3,2,2)*qp_i(3,3)) 
                         tmp23 =   fac *(tij_abcd(1,1,2,3)*qp_i(1,1)+&
                                         tij_abcd(2,1,2,3)*qp_i(2,1)+&
                                         tij_abcd(3,1,2,3)*qp_i(3,1)+&
                                         tij_abcd(1,2,2,3)*qp_i(1,2)+&
                                         tij_abcd(2,2,2,3)*qp_i(2,2)+&
                                         tij_abcd(3,2,2,3)*qp_i(3,2)+&
                                         tij_abcd(1,3,2,3)*qp_i(1,3)+&
                                         tij_abcd(2,3,2,3)*qp_i(2,3)+&
                                         tij_abcd(3,3,2,3)*qp_i(3,3)) 
                         tmp33 =   fac *(tij_abcd(1,1,3,3)*qp_i(1,1)+&
                                         tij_abcd(2,1,3,3)*qp_i(2,1)+&
                                         tij_abcd(3,1,3,3)*qp_i(3,1)+&
                                         tij_abcd(1,2,3,3)*qp_i(1,2)+&
                                         tij_abcd(2,2,3,3)*qp_i(2,2)+&
                                         tij_abcd(3,2,3,3)*qp_i(3,2)+&
                                         tij_abcd(1,3,3,3)*qp_i(1,3)+&
                                         tij_abcd(2,3,3,3)*qp_i(2,3)+&
                                         tij_abcd(3,3,3,3)*qp_i(3,3)) 

                         ef2_j(1,1) = ef2_j(1,1) - tmp11
                         ef2_j(1,2) = ef2_j(1,2) - tmp12
                         ef2_j(1,3) = ef2_j(1,3) - tmp13
                         ef2_j(2,1) = ef2_j(2,1) - tmp12
                         ef2_j(2,2) = ef2_j(2,2) - tmp22
                         ef2_j(2,3) = ef2_j(2,3) - tmp23
                         ef2_j(3,1) = ef2_j(3,1) - tmp13
                         ef2_j(3,2) = ef2_j(3,2) - tmp23
                         ef2_j(3,3) = ef2_j(3,3) - tmp33
                      END IF
                   END IF
                END IF
                IF (task(3,2)) THEN
                   ! Quadrupole - Dipole
                   fac = 1.0_dp/3.0_dp
                   ! Dipole i (locally B) - Quadrupole j (locally A)
                   tmp_ij = dp_i(1)*(tij_abc(1,1,1)*qp_j(1,1)+&
                                     tij_abc(2,1,1)*qp_j(2,1)+&
                                     tij_abc(3,1,1)*qp_j(3,1)+&
                                     tij_abc(1,2,1)*qp_j(1,2)+&
                                     tij_abc(2,2,1)*qp_j(2,2)+&
                                     tij_abc(3,2,1)*qp_j(3,2)+&
                                     tij_abc(1,3,1)*qp_j(1,3)+&
                                     tij_abc(2,3,1)*qp_j(2,3)+&
                                     tij_abc(3,3,1)*qp_j(3,3))+&
                            dp_i(2)*(tij_abc(1,1,2)*qp_j(1,1)+&
                                     tij_abc(2,1,2)*qp_j(2,1)+&
                                     tij_abc(3,1,2)*qp_j(3,1)+&
                                     tij_abc(1,2,2)*qp_j(1,2)+&
                                     tij_abc(2,2,2)*qp_j(2,2)+&
                                     tij_abc(3,2,2)*qp_j(3,2)+&
                                     tij_abc(1,3,2)*qp_j(1,3)+&
                                     tij_abc(2,3,2)*qp_j(2,3)+&
                                     tij_abc(3,3,2)*qp_j(3,3))+&
                            dp_i(3)*(tij_abc(1,1,3)*qp_j(1,1)+&
                                     tij_abc(2,1,3)*qp_j(2,1)+&
                                     tij_abc(3,1,3)*qp_j(3,1)+&
                                     tij_abc(1,2,3)*qp_j(1,2)+&
                                     tij_abc(2,2,3)*qp_j(2,2)+&
                                     tij_abc(3,2,3)*qp_j(3,2)+&
                                     tij_abc(1,3,3)*qp_j(1,3)+&
                                     tij_abc(2,3,3)*qp_j(2,3)+&
                                     tij_abc(3,3,3)*qp_j(3,3))

                   ! Dipole j (locally A) - Quadrupole i (locally B)
                   tmp_ji = dp_j(1)*(tij_abc(1,1,1)*qp_i(1,1)+&
                                     tij_abc(2,1,1)*qp_i(2,1)+&
                                     tij_abc(3,1,1)*qp_i(3,1)+&
                                     tij_abc(1,2,1)*qp_i(1,2)+&
                                     tij_abc(2,2,1)*qp_i(2,2)+&
                                     tij_abc(3,2,1)*qp_i(3,2)+&
                                     tij_abc(1,3,1)*qp_i(1,3)+&
                                     tij_abc(2,3,1)*qp_i(2,3)+&
                                     tij_abc(3,3,1)*qp_i(3,3))+&
                            dp_j(2)*(tij_abc(1,1,2)*qp_i(1,1)+&
                                     tij_abc(2,1,2)*qp_i(2,1)+&
                                     tij_abc(3,1,2)*qp_i(3,1)+&
                                     tij_abc(1,2,2)*qp_i(1,2)+&
                                     tij_abc(2,2,2)*qp_i(2,2)+&
                                     tij_abc(3,2,2)*qp_i(3,2)+&
                                     tij_abc(1,3,2)*qp_i(1,3)+&
                                     tij_abc(2,3,2)*qp_i(2,3)+&
                                     tij_abc(3,3,2)*qp_i(3,3))+&
                            dp_j(3)*(tij_abc(1,1,3)*qp_i(1,1)+&
                                     tij_abc(2,1,3)*qp_i(2,1)+&
                                     tij_abc(3,1,3)*qp_i(3,1)+&
                                     tij_abc(1,2,3)*qp_i(1,2)+&
                                     tij_abc(2,2,3)*qp_i(2,2)+&
                                     tij_abc(3,2,3)*qp_i(3,2)+&
                                     tij_abc(1,3,3)*qp_i(1,3)+&
                                     tij_abc(2,3,3)*qp_i(2,3)+&
                                     tij_abc(3,3,3)*qp_i(3,3))

                   tmp= fac * (tmp_ij - tmp_ji)
                   eloc = eloc + tmp
                   IF (do_forces) THEN
                      DO k = 1, 3
                         ! Dipole i (locally B) - Quadrupole j (locally A)
                         tmp_ij = dp_i(1)*(tij_abcd(1,1,1,k)*qp_j(1,1)+& 
                                           tij_abcd(2,1,1,k)*qp_j(2,1)+& 
                                           tij_abcd(3,1,1,k)*qp_j(3,1)+& 
                                           tij_abcd(1,2,1,k)*qp_j(1,2)+& 
                                           tij_abcd(2,2,1,k)*qp_j(2,2)+& 
                                           tij_abcd(3,2,1,k)*qp_j(3,2)+& 
                                           tij_abcd(1,3,1,k)*qp_j(1,3)+& 
                                           tij_abcd(2,3,1,k)*qp_j(2,3)+& 
                                           tij_abcd(3,3,1,k)*qp_j(3,3))+&
                                  dp_i(2)*(tij_abcd(1,1,2,k)*qp_j(1,1)+& 
                                           tij_abcd(2,1,2,k)*qp_j(2,1)+& 
                                           tij_abcd(3,1,2,k)*qp_j(3,1)+& 
                                           tij_abcd(1,2,2,k)*qp_j(1,2)+& 
                                           tij_abcd(2,2,2,k)*qp_j(2,2)+& 
                                           tij_abcd(3,2,2,k)*qp_j(3,2)+& 
                                           tij_abcd(1,3,2,k)*qp_j(1,3)+& 
                                           tij_abcd(2,3,2,k)*qp_j(2,3)+& 
                                           tij_abcd(3,3,2,k)*qp_j(3,3))+&
                                  dp_i(3)*(tij_abcd(1,1,3,k)*qp_j(1,1)+& 
                                           tij_abcd(2,1,3,k)*qp_j(2,1)+& 
                                           tij_abcd(3,1,3,k)*qp_j(3,1)+& 
                                           tij_abcd(1,2,3,k)*qp_j(1,2)+& 
                                           tij_abcd(2,2,3,k)*qp_j(2,2)+& 
                                           tij_abcd(3,2,3,k)*qp_j(3,2)+& 
                                           tij_abcd(1,3,3,k)*qp_j(1,3)+& 
                                           tij_abcd(2,3,3,k)*qp_j(2,3)+& 
                                           tij_abcd(3,3,3,k)*qp_j(3,3))

                         ! Dipole j (locally A) - Quadrupole i (locally B)
                         tmp_ji = dp_j(1)*(tij_abcd(1,1,1,k)*qp_i(1,1)+& 
                                           tij_abcd(2,1,1,k)*qp_i(2,1)+& 
                                           tij_abcd(3,1,1,k)*qp_i(3,1)+& 
                                           tij_abcd(1,2,1,k)*qp_i(1,2)+& 
                                           tij_abcd(2,2,1,k)*qp_i(2,2)+& 
                                           tij_abcd(3,2,1,k)*qp_i(3,2)+& 
                                           tij_abcd(1,3,1,k)*qp_i(1,3)+& 
                                           tij_abcd(2,3,1,k)*qp_i(2,3)+& 
                                           tij_abcd(3,3,1,k)*qp_i(3,3))+&
                                  dp_j(2)*(tij_abcd(1,1,2,k)*qp_i(1,1)+& 
                                           tij_abcd(2,1,2,k)*qp_i(2,1)+& 
                                           tij_abcd(3,1,2,k)*qp_i(3,1)+& 
                                           tij_abcd(1,2,2,k)*qp_i(1,2)+& 
                                           tij_abcd(2,2,2,k)*qp_i(2,2)+& 
                                           tij_abcd(3,2,2,k)*qp_i(3,2)+& 
                                           tij_abcd(1,3,2,k)*qp_i(1,3)+& 
                                           tij_abcd(2,3,2,k)*qp_i(2,3)+& 
                                           tij_abcd(3,3,2,k)*qp_i(3,3))+&
                                  dp_j(3)*(tij_abcd(1,1,3,k)*qp_i(1,1)+& 
                                           tij_abcd(2,1,3,k)*qp_i(2,1)+& 
                                           tij_abcd(3,1,3,k)*qp_i(3,1)+& 
                                           tij_abcd(1,2,3,k)*qp_i(1,2)+& 
                                           tij_abcd(2,2,3,k)*qp_i(2,2)+& 
                                           tij_abcd(3,2,3,k)*qp_i(3,2)+& 
                                           tij_abcd(1,3,3,k)*qp_i(1,3)+& 
                                           tij_abcd(2,3,3,k)*qp_i(2,3)+& 
                                           tij_abcd(3,3,3,k)*qp_i(3,3))

                         fr(k) = fr(k) - fac * (tmp_ij - tmp_ji)
                      END DO
                   END IF
                END IF
                IF (task(3,1)) THEN
                   ! Quadrupole - Charge
                   fac = 1.0_dp/3.0_dp

                   ! Quadrupole j (locally A) - Charge j (locally B)
                   tmp_ij = ch_i * (tij_ab(1,1)*qp_j(1,1)+&
                                    tij_ab(2,1)*qp_j(2,1)+&
                                    tij_ab(3,1)*qp_j(3,1)+&
                                    tij_ab(1,2)*qp_j(1,2)+&
                                    tij_ab(2,2)*qp_j(2,2)+&
                                    tij_ab(3,2)*qp_j(3,2)+&
                                    tij_ab(1,3)*qp_j(1,3)+&
                                    tij_ab(2,3)*qp_j(2,3)+&
                                    tij_ab(3,3)*qp_j(3,3)) 

                   ! Quadrupole i (locally B) - Charge j (locally A) 
                   tmp_ji = ch_j * (tij_ab(1,1)*qp_i(1,1)+&
                                    tij_ab(2,1)*qp_i(2,1)+&
                                    tij_ab(3,1)*qp_i(3,1)+&
                                    tij_ab(1,2)*qp_i(1,2)+&
                                    tij_ab(2,2)*qp_i(2,2)+&
                                    tij_ab(3,2)*qp_i(3,2)+&
                                    tij_ab(1,3)*qp_i(1,3)+&
                                    tij_ab(2,3)*qp_i(2,3)+&
                                    tij_ab(3,3)*qp_i(3,3))  

                   eloc = eloc + fac*(tmp_ij+tmp_ji)
                   IF (do_forces) THEN
                      DO k = 1, 3
                         ! Quadrupole j (locally A) - Charge i (locally B)  
                         tmp_ij = ch_i * (tij_abc(1,1,k)*qp_j(1,1)+&    
                                          tij_abc(2,1,k)*qp_j(2,1)+&    
                                          tij_abc(3,1,k)*qp_j(3,1)+&    
                                          tij_abc(1,2,k)*qp_j(1,2)+&    
                                          tij_abc(2,2,k)*qp_j(2,2)+&    
                                          tij_abc(3,2,k)*qp_j(3,2)+&    
                                          tij_abc(1,3,k)*qp_j(1,3)+&    
                                          tij_abc(2,3,k)*qp_j(2,3)+&    
                                          tij_abc(3,3,k)*qp_j(3,3))

                         ! Quadrupole i (locally B) - Charge j (locally A) 
                         tmp_ji = ch_j * (tij_abc(1,1,k)*qp_i(1,1)+&    
                                          tij_abc(2,1,k)*qp_i(2,1)+&    
                                          tij_abc(3,1,k)*qp_i(3,1)+&    
                                          tij_abc(1,2,k)*qp_i(1,2)+&    
                                          tij_abc(2,2,k)*qp_i(2,2)+&    
                                          tij_abc(3,2,k)*qp_i(3,2)+&    
                                          tij_abc(1,3,k)*qp_i(1,3)+&    
                                          tij_abc(2,3,k)*qp_i(2,3)+&    
                                          tij_abc(3,3,k)*qp_i(3,3))

                         fr(k) = fr(k) - fac *(tmp_ij + tmp_ji)
                      END DO
                   END IF
                END IF
                energy = energy + eloc
                IF (do_forces) THEN
                   forces(1,atom_a) = forces(1,atom_a) - fr(1)
                   forces(2,atom_a) = forces(2,atom_a) - fr(2)
                   forces(3,atom_a) = forces(3,atom_a) - fr(3)
                   forces(1,atom_b) = forces(1,atom_b) + fr(1)
                   forces(2,atom_b) = forces(2,atom_b) + fr(2)
                   forces(3,atom_b) = forces(3,atom_b) + fr(3)
                END IF
                ! Electric fields
                IF (do_efield) THEN
                   ! Potential
                   IF (do_efield0) THEN
                      efield0(  atom_a) = efield0(  atom_a) + ef0_j 

                      efield0(  atom_b) = efield0(  atom_b) + ef0_i 
                   END IF
                   ! Electric field
                   IF (do_efield1) THEN
                      efield1(1,atom_a) = efield1(1,atom_a) + ef1_j(1)
                      efield1(2,atom_a) = efield1(2,atom_a) + ef1_j(2)
                      efield1(3,atom_a) = efield1(3,atom_a) + ef1_j(3)

                      efield1(1,atom_b) = efield1(1,atom_b) + ef1_i(1)
                      efield1(2,atom_b) = efield1(2,atom_b) + ef1_i(2)
                      efield1(3,atom_b) = efield1(3,atom_b) + ef1_i(3)
                   END IF
                   ! Electric field gradient
                   IF (do_efield2) THEN
                      efield2(1,atom_a) = efield2(1,atom_a) + ef2_j(1,1)
                      efield2(2,atom_a) = efield2(2,atom_a) + ef2_j(1,2)
                      efield2(3,atom_a) = efield2(3,atom_a) + ef2_j(1,3)
                      efield2(4,atom_a) = efield2(4,atom_a) + ef2_j(2,1)
                      efield2(5,atom_a) = efield2(5,atom_a) + ef2_j(2,2)
                      efield2(6,atom_a) = efield2(6,atom_a) + ef2_j(2,3)
                      efield2(7,atom_a) = efield2(7,atom_a) + ef2_j(3,1)
                      efield2(8,atom_a) = efield2(8,atom_a) + ef2_j(3,2)
                      efield2(9,atom_a) = efield2(9,atom_a) + ef2_j(3,3)

                      efield2(1,atom_b) = efield2(1,atom_b) + ef2_i(1,1)
                      efield2(2,atom_b) = efield2(2,atom_b) + ef2_i(1,2)
                      efield2(3,atom_b) = efield2(3,atom_b) + ef2_i(1,3)
                      efield2(4,atom_b) = efield2(4,atom_b) + ef2_i(2,1)
                      efield2(5,atom_b) = efield2(5,atom_b) + ef2_i(2,2)
                      efield2(6,atom_b) = efield2(6,atom_b) + ef2_i(2,3)
                      efield2(7,atom_b) = efield2(7,atom_b) + ef2_i(3,1)
                      efield2(8,atom_b) = efield2(8,atom_b) + ef2_i(3,2)
                      efield2(9,atom_b) = efield2(9,atom_b) + ef2_i(3,3)
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
       charges, dipoles, quadrupoles, forces, pv, efield0, efield1, efield2,&
       error)
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    LOGICAL, DIMENSION(3, 3), INTENT(IN)     :: task
    LOGICAL, INTENT(IN)                      :: do_forces, do_efield, &
                                                do_stress
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: quadrupoles
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: forces, pv
    REAL(KIND=dp), DIMENSION(:), POINTER     :: efield0
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: efield1, efield2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_LR', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: atm_factor, atm_factor_st(3), &
                                                cnjg_fac, fac, summe_tmp
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: summe_ef
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: summe_st
    INTEGER :: gpt, group, handle, i, iparticle, iparticle_kind, &
      iparticle_local, j, lp, mp, nnodes, node, np, nparticle_kind, &
      nparticle_local, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    LOGICAL                                  :: do_efield0, do_efield1, &
                                                do_efield2, failure
    REAL(KIND=dp)                            :: alpha, denom, dipole_t(3), &
                                                factor, flops, four_alpha_sq, &
                                                gauss, pref, q_t, tmp, trq_t
    REAL(KIND=dp), DIMENSION(3)              :: tmp_v, vec
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
    do_efield0 = do_efield.AND.ASSOCIATED(efield0)
    do_efield1 = do_efield.AND.ASSOCIATED(efield1)
    do_efield2 = do_efield.AND.ASSOCIATED(efield2)

    ! Gathering data from the ewald environment
    CALL ewald_env_get (ewald_env, alpha=alpha, group=group, error=error)
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
    CALL structure_factor_allocate(pw_grid%bounds, nnodes, exp_igr)

    ALLOCATE (summe_ef(1:pw_grid%ngpts_cut), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    summe_ef = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
    ! Stress Tensor
    IF (do_stress) THEN
       ALLOCATE (summe_st(3,1:pw_grid%ngpts_cut), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       summe_st = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
    END IF

    ! Defining four_alpha_sq
    four_alpha_sq = 4.0_dp * alpha **2
    dipole_t      = 0.0_dp
    q_t           = 0.0_dp
    trq_t         = 0.0_dp
    ! Zero node count
    node = 0
    DO iparticle_kind = 1, nparticle_kind
       nparticle_local = local_particles%n_el(iparticle_kind)
       DO iparticle_local=1,nparticle_local
          node = node + 1
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          CALL matvec_3x3 (vec, cell%h_inv, particle_set(iparticle)%r)
          CALL structure_factor_evaluate (vec, pw_grid%npts, exp_igr%lb, &
               exp_igr%ex(:,node), exp_igr%ey(:,node), exp_igr%ez(:,node))

          ! Computing the total charge, dipole and quadrupole trace (if any)
          IF (ANY(task(1,:))) THEN
             q_t      = q_t      + charges(  iparticle)
          END IF
          IF (ANY(task(2,:))) THEN
             dipole_t = dipole_t + dipoles(:,iparticle)
          END IF
          IF (ANY(task(3,:))) THEN
             trq_t    = trq_t    + quadrupoles(1,1,iparticle)+&
                                   quadrupoles(2,2,iparticle)+&
                                   quadrupoles(3,3,iparticle)
          END IF
       END DO
    END DO

    ! Looping over the positive g-vectors
    DO gpt = 1, pw_grid%ngpts_cut
       lp = pw_grid%mapl%pos(pw_grid%g_hat(1, gpt))
       mp = pw_grid%mapm%pos(pw_grid%g_hat(2, gpt))
       np = pw_grid%mapn%pos(pw_grid%g_hat(3, gpt))

       lp = lp + bds(1,1)
       mp = mp + bds(1,2)
       np = np + bds(1,3)

       ! Initializing sum to be used in the energy and force
       node = 0
       DO iparticle_kind = 1, nparticle_kind
          nparticle_local = local_particles%n_el(iparticle_kind)
          DO iparticle_local=1,nparticle_local
             node = node + 1
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             ! Density for energy and forces
             CALL get_atom_factor(atm_factor, pw_grid, gpt, iparticle, task, charges,&
                  dipoles, quadrupoles, error)
             summe_tmp     = exp_igr%ex(lp,node)*exp_igr%ey(mp,node)*exp_igr%ez(np,node)
             summe_ef(gpt) = summe_ef(gpt) + atm_factor*summe_tmp

             ! Precompute pseudo-density for stress tensor calculation
             IF (do_stress) THEN
                CALL get_atom_factor_stress(atm_factor_st, pw_grid, gpt, iparticle, task,&
                     dipoles, quadrupoles, error)
                summe_st(1:3,gpt) = summe_st(1:3,gpt) + atm_factor_st(1:3) *summe_tmp
             END IF
          END DO
       END DO
    END DO
    flops = flops + 13.0_dp*(pw_grid%ngpts_cut*10.0_dp&
                  +          pw_grid%ngpts_cut*node*20.0_dp)

                    CALL mp_sum (     q_t, group)
                    CALL mp_sum (dipole_t, group)
                    CALL mp_sum (   trq_t, group)
                    CALL mp_sum (summe_ef, group)
    IF (do_stress)  CALL mp_sum (summe_st, group)

    ! Looping over the positive g-vectors
    DO gpt = 1, pw_grid%ngpts_cut
       ! computing the potential energy
       lp = pw_grid%mapl%pos(pw_grid%g_hat(1,gpt))
       mp = pw_grid%mapm%pos(pw_grid%g_hat(2,gpt))
       np = pw_grid%mapn%pos(pw_grid%g_hat(3,gpt))

       lp = lp + bds(1,1)
       mp = mp + bds(1,2)
       np = np + bds(1,3)

       IF (pw_grid%gsq(gpt) == 0.0_dp) THEN
          ! G=0 vector for dipole-dipole and charge-quadrupole
          energy = energy + (1.0_dp/6.0_dp)*DOT_PRODUCT(dipole_t,dipole_t)&
                          - (1.0_dp/9.0_dp)*q_t*trq_t
          IF (debug_e_field_en.OR.(.NOT.debug_this_module)) THEN
             ! This term is important and may give problems if one is debugging
             ! VS finite differences since it comes from a residual integral in
             ! the complex plane (cannot be reproduced with finite differences)
             node = 0
             DO iparticle_kind = 1, nparticle_kind
                nparticle_local = local_particles%n_el(iparticle_kind)
                DO iparticle_local=1,nparticle_local
                   node = node + 1
                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                   ! Corrections for G=0 to potential, field
                   IF (do_efield) THEN
                      ! Potential
                      IF (do_efield0) THEN
                         efield0(    iparticle) = efield0(    iparticle) - (1.0_dp/(18.0_dp))*trq_t
                      END IF
                      ! Electrostatic field
                      IF (do_efield1) THEN
                         efield1(1:3,iparticle) = efield1(1:3,iparticle) - (1.0_dp/6.0_dp)*(dipole_t(1:3)-dipoles(:,iparticle))
                      END IF
                      ! Electrostatic field gradients
                      IF (do_efield2) THEN
                         efield2(1,iparticle) = efield2(1,iparticle) - (1.0_dp/(18.0_dp))*(q_t-charges(  iparticle))
                         efield2(5,iparticle) = efield2(5,iparticle) - (1.0_dp/(18.0_dp))*(q_t-charges(  iparticle))
                         efield2(9,iparticle) = efield2(9,iparticle) - (1.0_dp/(18.0_dp))*(q_t-charges(  iparticle))
                      END IF
                   END IF
                END DO
             END DO
          END IF
          CYCLE
       END IF
       gauss  = (rho0(lp,mp,np) * pw_grid%vol)**2 / pw_grid%gsq(gpt)
       factor = gauss * REAL(summe_ef(gpt) * CONJG(summe_ef(gpt)),KIND=dp)
       energy = energy + factor

       node = 0
       DO iparticle_kind = 1, nparticle_kind
          nparticle_local = local_particles%n_el(iparticle_kind)
          DO iparticle_local=1,nparticle_local
             node = node + 1
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             fac       = exp_igr%ex(lp,node)*exp_igr%ey(mp,node)*exp_igr%ez(np,node)
             cnjg_fac  = CONJG(fac)

             ! Forces
             IF (do_forces) THEN
                CALL get_atom_factor(atm_factor, pw_grid, gpt, iparticle, task, charges,&
                     dipoles, quadrupoles, error)
                tmp      = gauss * AIMAG(summe_ef(gpt) * (cnjg_fac * CONJG(atm_factor)))
                forces(1,node) = forces(1,node) + tmp * pw_grid%g(1,gpt)
                forces(2,node) = forces(2,node) + tmp * pw_grid%g(2,gpt)
                forces(3,node) = forces(3,node) + tmp * pw_grid%g(3,gpt)
             END IF

             ! Electric field
             IF (do_efield) THEN
                ! Potential
                IF (do_efield0) THEN
                   efield0(  iparticle) = efield0(  iparticle) + gauss *  REAL(fac*CONJG(summe_ef(gpt)),KIND=dp) 
                END IF
                ! Electric field
                IF (do_efield1) THEN
                   tmp = 0.0_dp
                   IF (ANY(task(2,:))) THEN
                      tmp = DOT_PRODUCT(dipoles(:,iparticle),pw_grid%g(:,gpt))
                   END IF
                   tmp = (AIMAG(fac*CONJG(summe_ef(gpt)))-tmp)*gauss
                   efield1(1,iparticle) = efield1(1,iparticle) - tmp * pw_grid%g(1,gpt)
                   efield1(2,iparticle) = efield1(2,iparticle) - tmp * pw_grid%g(2,gpt)
                   efield1(3,iparticle) = efield1(3,iparticle) - tmp * pw_grid%g(3,gpt)
                END IF
                ! Electric field gradient
                IF (do_efield2) THEN
                   tmp = 0.0_dp
                   IF (ANY(task(1,:))) THEN
                      tmp = -charges(iparticle)
                   END IF
                   IF (ANY(task(3,:))) THEN
                      DO i = 1,3
                         DO j = 1,3
                            tmp = tmp + 1.0_dp/3.0_dp*quadrupoles(j,i,iparticle)*pw_grid%g(j,gpt)*pw_grid%g(i,gpt)
                         END DO
                      END DO
                   END IF
                   tmp_v(1) = (tmp + REAL(fac*CONJG(summe_ef(gpt)),KIND=dp))*pw_grid%g(1,gpt)*gauss
                   tmp_v(2) = (tmp + REAL(fac*CONJG(summe_ef(gpt)),KIND=dp))*pw_grid%g(2,gpt)*gauss
                   tmp_v(3) = (tmp + REAL(fac*CONJG(summe_ef(gpt)),KIND=dp))*pw_grid%g(3,gpt)*gauss

                   efield2(1,iparticle) = efield2(1,iparticle) + tmp_v(1) * pw_grid%g(1,gpt)
                   efield2(2,iparticle) = efield2(2,iparticle) + tmp_v(1) * pw_grid%g(2,gpt)
                   efield2(3,iparticle) = efield2(3,iparticle) + tmp_v(1) * pw_grid%g(3,gpt)
                   efield2(4,iparticle) = efield2(4,iparticle) + tmp_v(2) * pw_grid%g(1,gpt)
                   efield2(5,iparticle) = efield2(5,iparticle) + tmp_v(2) * pw_grid%g(2,gpt)
                   efield2(6,iparticle) = efield2(6,iparticle) + tmp_v(2) * pw_grid%g(3,gpt)
                   efield2(7,iparticle) = efield2(7,iparticle) + tmp_v(3) * pw_grid%g(1,gpt)
                   efield2(8,iparticle) = efield2(8,iparticle) + tmp_v(3) * pw_grid%g(2,gpt)
                   efield2(9,iparticle) = efield2(9,iparticle) + tmp_v(3) * pw_grid%g(3,gpt)
                END IF
             END IF
          END DO
       END DO

       ! Compute the virial P*V
       IF (do_stress) THEN
          ! The Stress Tensor can be decomposed in two main components.
          ! The first one is just a normal ewald component for reciprocal space
          denom = 1.0_dp/four_alpha_sq + 1.0_dp/pw_grid%gsq(gpt)
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
          pv(1,1) = pv(1,1) + 2.0_dp * gauss * pw_grid%g(1,gpt) * REAL(summe_st(1,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv(1,2) = pv(1,2) + 2.0_dp * gauss * pw_grid%g(1,gpt) * REAL(summe_st(2,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv(1,3) = pv(1,3) + 2.0_dp * gauss * pw_grid%g(1,gpt) * REAL(summe_st(3,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv(2,1) = pv(2,1) + 2.0_dp * gauss * pw_grid%g(2,gpt) * REAL(summe_st(1,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv(2,2) = pv(2,2) + 2.0_dp * gauss * pw_grid%g(2,gpt) * REAL(summe_st(2,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv(2,3) = pv(2,3) + 2.0_dp * gauss * pw_grid%g(2,gpt) * REAL(summe_st(3,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv(3,1) = pv(3,1) + 2.0_dp * gauss * pw_grid%g(3,gpt) * REAL(summe_st(1,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv(3,2) = pv(3,2) + 2.0_dp * gauss * pw_grid%g(3,gpt) * REAL(summe_st(2,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv(3,3) = pv(3,3) + 2.0_dp * gauss * pw_grid%g(3,gpt) * REAL(summe_st(3,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
       END IF
    END DO

    flops = flops + pw_grid%ngpts_cut * 55.0_dp &
                  + pw_grid%ngpts_cut * node * 23.0_dp

    pref   = fourpi/pw_grid%vol
    energy = energy * pref

    CALL structure_factor_deallocate (exp_igr)
    DEALLOCATE (summe_ef, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (do_stress) THEN
       pv = pv * pref
       DEALLOCATE (summe_st, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (do_forces) THEN
       forces  = 2.0_dp * forces  * pref
    END IF
    IF (do_efield0) THEN
       efield0 = 2.0_dp * efield0 * pref
    END IF
    IF (do_efield1) THEN
       efield1 = 2.0_dp * efield1 * pref
    END IF
    IF (do_efield2) THEN
       efield2 = 2.0_dp * efield2 * pref
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
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: quadrupoles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_atom_factor', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: tmp
    INTEGER                                  :: i, j

    atm_factor = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
    IF (task(1,1)) THEN
       ! Charge
       atm_factor = atm_factor + charges(iparticle)
    END IF
    IF (task(2,2)) THEN
       ! Dipole
       tmp = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
       DO i = 1,3
          tmp = tmp + dipoles(i,iparticle)*pw_grid%g(i,gpt)
       END DO
       atm_factor = atm_factor + tmp * CMPLX(0.0_dp, -1.0_dp, KIND=dp)
    END IF
    IF (task(3,3)) THEN
       ! Quadrupole
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
!> \brief Computes the atom factor including charge, dipole and quadrupole
!> \par History
!>      none
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich 
! *****************************************************************************
  SUBROUTINE get_atom_factor_stress(atm_factor, pw_grid, gpt, iparticle, task,&
       dipoles, quadrupoles, error)
    COMPLEX(KIND=dp), INTENT(OUT)            :: atm_factor(3)
    TYPE(pw_grid_type), POINTER              :: pw_grid
    INTEGER, INTENT(IN)                      :: gpt
    INTEGER                                  :: iparticle
    LOGICAL, DIMENSION(3, 3), INTENT(IN)     :: task
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: quadrupoles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_atom_factor_stress', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    atm_factor = CMPLX (0.0_dp, 0.0_dp,KIND=dp)
    IF (ANY(task(2,:))) THEN
       ! Dipole
       atm_factor = dipoles(:,iparticle) * CMPLX(0.0_dp, -1.0_dp, KIND=dp)
    END IF
    IF (ANY(task(3,:))) THEN
       ! Quadrupole
       DO i = 1,3
          atm_factor(1) = atm_factor(1) - 1.0_dp/3.0_dp *&
              (quadrupoles(1,i,iparticle)*pw_grid%g(i,gpt)+&
               quadrupoles(i,1,iparticle)*pw_grid%g(i,gpt))
          atm_factor(2) = atm_factor(2) - 1.0_dp/3.0_dp *&
              (quadrupoles(2,i,iparticle)*pw_grid%g(i,gpt)+&
               quadrupoles(i,2,iparticle)*pw_grid%g(i,gpt))
          atm_factor(3) = atm_factor(3) - 1.0_dp/3.0_dp *&
              (quadrupoles(3,i,iparticle)*pw_grid%g(i,gpt)+&
               quadrupoles(i,3,iparticle)*pw_grid%g(i,gpt))
       END DO
    END IF

  END SUBROUTINE get_atom_factor_stress

! *****************************************************************************
!> \brief Computes the self interaction from g-space and the neutralizing background
!>     when using multipoles
!> \author Teodoro Laino [tlaino] - University of Zurich - 12.2007
! *****************************************************************************
  SUBROUTINE ewald_multipole_self (ewald_env, cell, local_particles, e_self, e_neut,&
       task, do_efield, charges, dipoles, quadrupoles, efield0, efield1, efield2,&
       error)
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cell_type), INTENT(IN)              :: cell
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(OUT)               :: e_self, e_neut
    LOGICAL, DIMENSION(3, 3), INTENT(IN)     :: task
    LOGICAL, INTENT(IN)                      :: do_efield
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: quadrupoles
    REAL(KIND=dp), DIMENSION(:), POINTER     :: efield0
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: efield1, efield2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_self', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: f23 = 2.0_dp/3.0_dp, &
                                                f415 = 4.0_dp/15.0_dp

    INTEGER                                  :: ewald_type, group, i, &
                                                iparticle, iparticle_kind, &
                                                iparticle_local, j, &
                                                nparticle_local
    LOGICAL                                  :: do_efield0, do_efield1, &
                                                do_efield2
    REAL(KIND=dp) :: alpha, ch_qu_self, ch_qu_self_tmp, dipole_self, fac1, &
      fac2, fac3, fac4, q, q_self, q_sum, qu_qu_self, qu_trace

    CALL ewald_env_get ( ewald_env, ewald_type=ewald_type, alpha=alpha,&
         group=group, error=error)

    do_efield0 = do_efield.AND.ASSOCIATED(efield0)
    do_efield1 = do_efield.AND.ASSOCIATED(efield1)
    do_efield2 = do_efield.AND.ASSOCIATED(efield2)
    q_self      = 0.0_dp
    q_sum       = 0.0_dp
    dipole_self = 0.0_dp
    ch_qu_self  = 0.0_dp
    qu_qu_self  = 0.0_dp
    fac1        = 2.0_dp*alpha*oorootpi
    fac2        = (f23**2)*(alpha**3)*oorootpi
    fac3        = (2.0_dp*oorootpi)*f23*alpha**3
    fac4        = (4.0_dp*oorootpi)*f415*alpha**5
    DO iparticle_kind=1,SIZE(local_particles%n_el)
       nparticle_local = local_particles%n_el(iparticle_kind)
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          IF (ANY(task(1,:))) THEN
             ! Charge - Charge
             q = charges(iparticle)
             q_self = q_self + q*q
             q_sum  = q_sum  + q
             ! Potential energy
             IF (do_efield0) THEN
                efield0(  iparticle) = efield0(  iparticle) - q*fac1
             END IF

             IF (task(1,3)) THEN
                ! Charge - Quadrupole
                ch_qu_self_tmp = 0.0_dp
                qu_trace = 0.0_dp
                DO i = 1,3
                   ch_qu_self_tmp = ch_qu_self_tmp + quadrupoles(i,i,iparticle) * q
                   qu_trace       = qu_trace       + quadrupoles(i,i,iparticle)
                END DO
                ch_qu_self = ch_qu_self + ch_qu_self_tmp
                ! Potential energy
                IF (do_efield0) THEN
                   efield0(  iparticle) = efield0(  iparticle) + fac2 * qu_trace
                END IF
             END IF
          END IF
          IF (ANY(task(2,:))) THEN
             ! Dipole - Dipole
             DO i = 1,3
                dipole_self = dipole_self + dipoles(i,iparticle)**2
             END DO
          END IF
          IF (ANY(task(3,:))) THEN
             ! Quadrupole - Quadrupole
             DO i = 1,3 
                DO j = 1,3
                   qu_qu_self = qu_qu_self + quadrupoles(j,i,iparticle)**2
                END DO
             END DO
          END IF
       END DO
    END DO

    CALL mp_sum (q_self, group)
    CALL mp_sum (q_sum, group)
    CALL mp_sum (dipole_self, group)
    CALL mp_sum (ch_qu_self, group)
    CALL mp_sum (qu_qu_self, group)

    e_self = -(q_self+f23*(dipole_self-f23*ch_qu_self+f415*qu_qu_self*alpha**2)*alpha**2)*alpha*oorootpi
    fac1   = pi/(2.0_dp * alpha ** 2 * cell % deth)
    e_neut = -q_sum**2 * fac1

    ! Correcting Potential for the neutralizing background charge
    DO iparticle_kind=1,SIZE(local_particles%n_el)
       nparticle_local = local_particles%n_el(iparticle_kind)
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          IF (ANY(task(1,:))) THEN
             q = charges(iparticle)
             ! Potential energy
             IF (do_efield0) THEN
                efield0(  iparticle) = efield0(  iparticle) - 2.0_dp*ABS(q_sum)*fac1
             END IF
          END IF
       END DO
    END DO
  END SUBROUTINE ewald_multipole_self

! *****************************************************************************
!> \author Teodoro Laino [tlaino] - University of Zurich - 12.2007
! *****************************************************************************
  SUBROUTINE ewald_multipole_print(iw, e_gspace, e_nonbond, e_rspace, e_self, e_neut)

    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), INTENT(IN)                :: e_gspace, e_nonbond, &
                                                e_rspace, e_self, e_neut

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_print', &
      routineP = moduleN//':'//routineN

    TYPE(unit_convert_type), POINTER         :: units

    CALL get_cp2k_units ( units )

    IF (iw>0) THEN
       WRITE ( iw, '( A, A )' ) ' *********************************', &
            '**********************************************'
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' INITIAL GSPACE ENERGY', &
            units%e_label, '= ', e_gspace
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' INITIAL RSPACE ENERGY', &
            units%e_label, '= ', e_rspace
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' INITIAL NONBONDED ENERGY', &
            units%e_label, '= ', e_nonbond-e_rspace
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' SELF ENERGY CORRECTION', &
            units%e_label, '= ', e_self
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' NEUTRALIZ. BCKGR. ENERGY', &
            units%e_label, '= ', e_neut
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' TOTAL ELECTROSTATIC EN.', &
            units%e_label, '= ', e_rspace+e_gspace+e_self+e_neut
       WRITE ( iw, '( A, A )' ) ' *********************************', &
            '**********************************************'
    END IF
  END SUBROUTINE ewald_multipole_print

END MODULE ewalds_multipole
