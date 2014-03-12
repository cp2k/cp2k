!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Treats the electrostatic for multipoles (up to quadrupoles)
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
!> inclusion of optional electric field damping for the polarizable atoms
!> Rodolphe Vuilleumier and Mathieu Salanne - 12.2009
! *****************************************************************************
MODULE ewalds_multipole
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE bibliography,                    ONLY: Aguado2003,&
                                             Laino2008,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE damping_dipole_types,            ONLY: damping_type,&
                                             no_damping,&
                                             tang_toennies
  USE dg_rho0_types,                   ONLY: dg_rho0_type
  USE dg_types,                        ONLY: dg_get,&
                                             dg_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE erf_fn,                          ONLY: erf,&
                                             erfc
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                             ewald_pw_type
  USE f77_blas
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             oorootpi,&
                                             pi,&
                                             sqrthalf
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
#include "./common/cp_common_uses.f90"

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
!> \par    Note
!>         When a gaussian multipole is used instead of point multipole, i.e.
!>         when radii(i)>0, the electrostatic fields (efield0, efield1, efield2)
!>         become derivatives of the electrostatic potential energy towards
!>         these gaussian multipoles.
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! *****************************************************************************
  RECURSIVE SUBROUTINE ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env,&
       cell, particle_set, local_particles, energy_local, energy_glob, e_neut, e_self,&
       task, do_correction_bonded, do_forces, do_stress, do_efield, radii, charges, dipoles,&
       quadrupoles, forces_local, forces_glob, pv_local, pv_glob, efield0, efield1,&
       efield2, iw, do_debug, atomic_kind_set, mm_section, error )
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), POINTER             :: particle_set(:)
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(INOUT)             :: energy_local, energy_glob
    REAL(KIND=dp), INTENT(OUT)               :: e_neut, e_self
    LOGICAL, DIMENSION(3), INTENT(IN)        :: task
    LOGICAL, INTENT(IN)                      :: do_correction_bonded, &
                                                do_forces, do_stress, &
                                                do_efield
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: radii, charges
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

    INTEGER                                  :: group, handle, i, j, size1, &
                                                size2, stat
    LOGICAL                                  :: check_debug, check_efield, &
                                                check_forces, do_task(3), &
                                                failure
    LOGICAL, DIMENSION(3, 3)                 :: my_task
    REAL(KIND=dp)                            :: e_bonded, e_bonded_t, &
                                                e_rspace, e_rspace_t, &
                                                energy_glob_t
    REAL(KIND=dp), DIMENSION(:), POINTER     :: efield0_lr, efield0_sr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: efield1_lr, efield1_sr, &
                                                efield2_lr, efield2_sr

    failure = .FALSE.
    CALL cite_reference(Aguado2003)
    CALL cite_reference(Laino2008)
    CALL timeset(routineN,handle)
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
               cell, particle_set, local_particles, radii, charges, dipoles,&
               quadrupoles, task, iw, atomic_kind_set, mm_section, error)
          STOP "Debug Multipole Requested:  POT+EFIELDS+GRAD!"
       END IF
       ! Debug the potential, electric fields and electric fields gradient in oder
       ! to retrieve the correct energy
       IF (debug_e_field_en) THEN
          CALL debug_ewald_multipoles_fields2(ewald_env, ewald_pw, nonbond_env,&
               cell, particle_set, local_particles, radii, charges, dipoles,&
               quadrupoles, task, iw, error)
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
    e_bonded = 0.0_dp
    IF ((.NOT.debug_g_space) .AND. (nonbond_env%do_nonbonded)) THEN
       ! Compute the Real Space (Short-Range) part of the Ewald sum.
       ! This contribution is only added when the nonbonded flag in the input
       ! is set, because these contributions depend. the neighborlists.
       CALL ewald_multipole_SR (nonbond_env, ewald_env, atomic_kind_set,&
            particle_set, cell, e_rspace, my_task,&
            do_forces, do_efield, do_stress, radii, charges, dipoles, quadrupoles,&
            forces_glob, pv_glob, efield0_sr, efield1_sr, efield2_sr, error)
       energy_glob = energy_glob + e_rspace

       IF (do_correction_bonded) THEN
          ! The corrections for bonded interactions are stored in the Real Space
          ! (Short-Range) part of the fields array.
          CALL ewald_multipole_bonded(nonbond_env, particle_set, ewald_env, &
               cell, e_bonded, my_task, do_forces, do_efield, do_stress, &
               charges, dipoles, quadrupoles, forces_glob, pv_glob, &
               efield0_sr, efield1_sr, efield2_sr, error)
          energy_glob = energy_glob + e_bonded
       END IF
    END IF

    e_neut       = 0.0_dp
    e_self       = 0.0_dp
    energy_local = 0.0_dp
    IF (.NOT.debug_r_space) THEN
       ! Compute the Reciprocal Space (Long-Range) part of the Ewald sum
       CALL ewald_multipole_LR(ewald_env, ewald_pw, cell, particle_set, &
            local_particles, energy_local, my_task, do_forces, do_efield, do_stress,&
            charges, dipoles, quadrupoles, forces_local, pv_local, efield0_lr, efield1_lr,&
            efield2_lr, error)

       ! Self-Interactions corrections
       CALL ewald_multipole_self (ewald_env, cell, local_particles, e_self, &
            e_neut, my_task, do_efield, radii, charges, dipoles, quadrupoles, &
            efield0_lr, efield1_lr, efield2_lr, error)
    END IF

    ! Sumup energy contributions for possible IO
    CALL ewald_env_get (ewald_env, group=group, error=error)
    energy_glob_t = energy_glob
    e_rspace_t    = e_rspace
    e_bonded_t    = e_bonded
    CALL mp_sum(energy_glob_t, group)
    CALL mp_sum(e_rspace_t, group)
    CALL mp_sum(e_bonded_t, group)
    ! Print some info about energetics
    CALL ewald_multipole_print (iw, energy_local, e_rspace_t, e_bonded_t, e_self, e_neut)

    ! Gather the components of the potential, fields and electrostatic field gradients
    IF (do_efield) THEN
       IF (PRESENT(efield0)) THEN
          efield0 = efield0_sr + efield0_lr
          CALL mp_sum(efield0, group)
          DEALLOCATE (efield0_sr, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (efield0_lr, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(efield1)) THEN
          efield1 = efield1_sr + efield1_lr
          CALL mp_sum(efield1, group)
          DEALLOCATE (efield1_sr, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (efield1_lr, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(efield2)) THEN
          efield2 = efield2_sr + efield2_lr
          CALL mp_sum(efield2, group)
          DEALLOCATE (efield2_sr, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (efield2_lr, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE ewald_multipole_evaluate

! *****************************************************************************
!> \brief computes the potential and the force for a lattice sum of multipoles
!>      up to quadrupole - Short Range (Real Space) Term
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE ewald_multipole_SR (nonbond_env, ewald_env, atomic_kind_set,&
       particle_set, cell, energy, task,&
       do_forces, do_efield, do_stress, radii, charges, dipoles, quadrupoles,&
       forces, pv, efield0, efield1, efield2, error)
    TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(particle_type), POINTER             :: particle_set(:)
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    LOGICAL, DIMENSION(3, 3), INTENT(IN)     :: task
    LOGICAL, INTENT(IN)                      :: do_forces, do_efield, &
                                                do_stress
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: radii, charges
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

    INTEGER :: a, atom_a, atom_b, b, c, d, e, handle, i, iend, igrp, ikind, &
      ilist, ipair, istart, itype_ij, itype_ji, jkind, k, kind_a, kind_b, kk, &
      nkdamp_ij, nkdamp_ji, nkinds, npairs
    INTEGER, DIMENSION(:, :), POINTER        :: list
    LOGICAL                                  :: do_efield0, do_efield1, &
                                                do_efield2, failure, &
                                                force_eval
    REAL(KIND=dp) :: alpha, beta, ch_i, ch_j, dampa_ij, dampa_ji, dampaexpi, &
      dampaexpj, dampfac_ij, dampfac_ji, dampfuncdiffi, dampfuncdiffj, &
      dampfunci, dampfuncj, dampsumfi, dampsumfj, ef0_i, ef0_j, eloc, fac, &
      fac_ij, factorial, ir, irab2, ptens11, ptens12, ptens13, ptens21, &
      ptens22, ptens23, ptens31, ptens32, ptens33, r, rab2, rab2_max, radius, &
      rcut, tij, tmp, tmp1, tmp11, tmp12, tmp13, tmp2, tmp21, tmp22, tmp23, &
      tmp31, tmp32, tmp33, tmp_ij, tmp_ji, xf
    REAL(KIND=dp), DIMENSION(0:5)            :: f
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, cvi, damptij_a, &
                                                damptji_a, dp_i, dp_j, ef1_i, &
                                                ef1_j, fr, rab, tij_a
    REAL(KIND=dp), DIMENSION(3, 3)           :: damptij_ab, damptji_ab, &
                                                ef2_i, ef2_j, qp_i, qp_j, &
                                                tij_ab
    REAL(KIND=dp), DIMENSION(3, 3, 3)        :: tij_abc
    REAL(KIND=dp), DIMENSION(3, 3, 3, 3)     :: tij_abcd
    REAL(KIND=dp), DIMENSION(3, 3, 3, 3, 3)  :: tij_abcde
    TYPE(damping_type)                       :: damping_ij, damping_ji
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update, &
                                                r_last_update_pbc

    failure = .FALSE.
    CALL timeset ( routineN, handle )
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

          itype_ij=no_damping
          nkdamp_ij=1
          dampa_ij=1.0_dp
          dampfac_ij=0.0_dp

          itype_ji=no_damping
          nkdamp_ji=1
          dampa_ji=1.0_dp
          dampfac_ji=0.0_dp
          IF (PRESENT(atomic_kind_set)) THEN
             IF (ASSOCIATED(atomic_kind_set(jkind)%damping)) THEN
                damping_ij=atomic_kind_set(jkind)%damping%damp(ikind)
                itype_ij=damping_ij%itype
                nkdamp_ij=damping_ij%order
                dampa_ij=damping_ij%bij
                dampfac_ij=damping_ij%cij
             END IF

             IF (ASSOCIATED(atomic_kind_set(ikind)%damping)) THEN
                damping_ji=atomic_kind_set(ikind)%damping%damp(jkind)
                itype_ji=damping_ji%itype
                nkdamp_ji=damping_ji%order
                dampa_ji=damping_ji%bij
                dampfac_ji=damping_ji%cij
             END IF
          END IF

          Pairs: DO ipair = istart, iend
             IF (ipair <= neighbor_kind_pair%nscale) THEN
                ! scale the electrostatic interaction if needed
                ! (most often scaled to zero)
                fac_ij = neighbor_kind_pair%ei_scale(ipair)
                IF (fac_ij<=0) CYCLE
             ELSE
                fac_ij = 1.0_dp
             END IF
             atom_a = list(1,ipair)
             atom_b = list(2,ipair)
             kind_a=particle_set(atom_a)%atomic_kind%kind_number
             kind_b=particle_set(atom_b)%atomic_kind%kind_number
             IF (atom_a==atom_b) fac_ij = 0.5_dp
             rab    = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
             rab    = rab + cell_v
             rab2   = rab(1)**2 + rab(2)**2 + rab(3)**2
             IF (rab2 <= rab2_max) THEN
                IF (PRESENT(radii)) THEN
                   radius = SQRT(radii(atom_a)*radii(atom_a) + radii(atom_b)*radii(atom_b))
                ELSE
                   radius = 0.0_dp
                END IF
                IF (radius > 0.0_dp) THEN
                   beta = sqrthalf/radius
#define SCREENED_COULOMB_GAUSS
#define STORE_ENERGY
#define STORE_FORCES
#include "ewalds_multipole_sr.f90"
#undef STORE_FORCES
#undef STORE_ENERGY
#undef SCREENED_COULOMB_GAUSS
               ELSE
#define SCREENED_COULOMB_ERFC
#define STORE_ENERGY
#define STORE_FORCES
#define DAMPING
#include "ewalds_multipole_sr.f90"
#undef DAMPING
#undef STORE_FORCES
#undef STORE_ENERGY
#undef SCREENED_COULOMB_ERFC
               END IF
             END IF
          END DO Pairs
       END DO Kind_Group_Loop
    END DO Lists
    IF (do_stress) THEN
       pv(1,1) = pv(1,1) +  ptens11
       pv(1,2) = pv(1,2) + (ptens12+ptens21)*0.5_dp
       pv(1,3) = pv(1,3) + (ptens13+ptens31)*0.5_dp
       pv(2,1) = pv(1,2)
       pv(2,2) = pv(2,2) +  ptens22
       pv(2,3) = pv(2,3) + (ptens23+ptens32)*0.5_dp
       pv(3,1) = pv(1,3)
       pv(3,2) = pv(2,3)
       pv(3,3) = pv(3,3) +  ptens33
    END IF

    CALL timestop ( handle )
  END SUBROUTINE ewald_multipole_SR

! *****************************************************************************
!> \brief computes the bonded correction for the potential and the force for a
!>        lattice sum of multipoles up to quadrupole
!> \author Teodoro Laino [tlaino] - 05.2009
! *****************************************************************************
  SUBROUTINE ewald_multipole_bonded (nonbond_env, particle_set, ewald_env, &
       cell, energy, task, do_forces, do_efield, do_stress, charges, &
       dipoles, quadrupoles, forces, pv, efield0, efield1, efield2, error)

    TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
    TYPE(particle_type), POINTER             :: particle_set( : )
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

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_bonded', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a, atom_a, atom_b, b, c, d, &
                                                e, handle, i, iend, igrp, &
                                                ilist, ipair, istart, k, &
                                                nscale
    INTEGER, DIMENSION(:, :), POINTER        :: list
    LOGICAL                                  :: do_efield0, do_efield1, &
                                                do_efield2, force_eval
    REAL(KIND=dp) :: alpha, ch_i, ch_j, ef0_i, ef0_j, eloc, fac, fac_ij, ir, &
      irab2, ptens11, ptens12, ptens13, ptens21, ptens22, ptens23, ptens31, &
      ptens32, ptens33, r, rab2, tij, tmp, tmp1, tmp11, tmp12, tmp13, tmp2, &
      tmp21, tmp22, tmp23, tmp31, tmp32, tmp33, tmp_ij, tmp_ji
    REAL(KIND=dp), DIMENSION(0:5)            :: f
    REAL(KIND=dp), DIMENSION(3)              :: dp_i, dp_j, ef1_i, ef1_j, fr, &
                                                rab, tij_a
    REAL(KIND=dp), DIMENSION(3, 3)           :: ef2_i, ef2_j, qp_i, qp_j, &
                                                tij_ab
    REAL(KIND=dp), DIMENSION(3, 3, 3)        :: tij_abc
    REAL(KIND=dp), DIMENSION(3, 3, 3, 3)     :: tij_abcd
    REAL(KIND=dp), DIMENSION(3, 3, 3, 3, 3)  :: tij_abcde
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair

    CALL timeset ( routineN, handle )
    do_efield0 = do_efield.AND.ASSOCIATED(efield0)
    do_efield1 = do_efield.AND.ASSOCIATED(efield1)
    do_efield2 = do_efield.AND.ASSOCIATED(efield2)
    IF (do_stress) THEN
       ptens11 = 0.0_dp ; ptens12 = 0.0_dp ; ptens13 = 0.0_dp
       ptens21 = 0.0_dp ; ptens22 = 0.0_dp ; ptens23 = 0.0_dp
       ptens31 = 0.0_dp ; ptens32 = 0.0_dp ; ptens33 = 0.0_dp
    END IF
    CALL ewald_env_get (ewald_env, alpha=alpha, error=error)
    CALL fist_nonbond_env_get(nonbond_env, nonbonded=nonbonded, error=error)

    ! Starting the force loop
    Lists: DO ilist=1,nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
       nscale = neighbor_kind_pair%nscale
       IF (nscale==0) CYCLE
       list => neighbor_kind_pair%list
       Kind_Group_Loop: DO igrp = 1, neighbor_kind_pair%ngrp_kind
          istart = neighbor_kind_pair%grp_kind_start(igrp)
          IF (istart>nscale) CYCLE
          iend = MIN(neighbor_kind_pair%grp_kind_end(igrp), nscale)
          Pairs: DO ipair = istart, iend
             ! only use pairs that are (partially) excluded for electrostatics
             fac_ij = -1.0_dp + neighbor_kind_pair%ei_scale(ipair)
             IF (fac_ij>=0) CYCLE

             atom_a = list(1,ipair)
             atom_b = list(2,ipair)

             rab = particle_set(atom_b)%r - particle_set(atom_a)%r
             rab = pbc(rab, cell)
             rab2 = rab(1)**2 + rab(2)**2 + rab(3)**2
#define SCREENED_COULOMB_ERF
#define STORE_ENERGY
#define STORE_FORCES
#include "ewalds_multipole_sr.f90"
#undef STORE_FORCES
#undef STORE_ENERGY
#undef SCREENED_COULOMB_ERF
          END DO Pairs
       END DO Kind_Group_Loop
    END DO Lists
    IF (do_stress) THEN
       pv(1,1) = pv(1,1) +  ptens11
       pv(1,2) = pv(1,2) + (ptens12+ptens21)*0.5_dp
       pv(1,3) = pv(1,3) + (ptens13+ptens31)*0.5_dp
       pv(2,1) = pv(1,2)
       pv(2,2) = pv(2,2) +  ptens22
       pv(2,3) = pv(2,3) + (ptens23+ptens32)*0.5_dp
       pv(3,1) = pv(1,3)
       pv(3,2) = pv(2,3)
       pv(3,3) = pv(3,3) +  ptens33
    END IF

    CALL timestop ( handle )
  END SUBROUTINE ewald_multipole_bonded

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
    INTEGER :: gpt, group, handle, iparticle, iparticle_kind, &
      iparticle_local, lp, mp, nnodes, node, np, nparticle_kind, &
      nparticle_local, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    LOGICAL                                  :: do_efield0, do_efield1, &
                                                do_efield2, failure
    REAL(KIND=dp)                            :: alpha, denom, dipole_t(3), &
                                                f0, factor, four_alpha_sq, &
                                                gauss, pref, q_t, tmp, trq_t
    REAL(KIND=dp), DIMENSION(3)              :: tmp_v, vec
    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_tmp
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: rho0
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(dg_type), POINTER                   :: dg
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(structure_factor_type)              :: exp_igr

    CALL timeset(routineN,handle)
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
       pv_tmp = 0.0_dp
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
          ! Stress tensor
          IF (do_stress) THEN
             pv_tmp(1,1) = pv_tmp(1,1) + (1.0_dp/6.0_dp)*DOT_PRODUCT(dipole_t,dipole_t)
             pv_tmp(2,2) = pv_tmp(2,2) + (1.0_dp/6.0_dp)*DOT_PRODUCT(dipole_t,dipole_t)
             pv_tmp(3,3) = pv_tmp(3,3) + (1.0_dp/6.0_dp)*DOT_PRODUCT(dipole_t,dipole_t)
          END IF
          ! Corrections for G=0 to potential, field and field gradient
          IF (do_efield.AND.(debug_e_field_en.OR.(.NOT.debug_this_module))) THEN
             ! This term is important and may give problems if one is debugging
             ! VS finite differences since it comes from a residual integral in
             ! the complex plane (cannot be reproduced with finite differences)
             node = 0
             DO iparticle_kind = 1, nparticle_kind
                nparticle_local = local_particles%n_el(iparticle_kind)
                DO iparticle_local=1,nparticle_local
                   node = node + 1
                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)

                   ! Potential
                   IF (do_efield0) THEN
                      efield0(    iparticle) = efield0(    iparticle)
                   END IF
                   ! Electrostatic field
                   IF (do_efield1) THEN
                      efield1(1:3,iparticle) = efield1(1:3,iparticle) - (1.0_dp/6.0_dp)*dipole_t(1:3)
                   END IF
                   ! Electrostatic field gradients
                   IF (do_efield2) THEN
                      efield2(1,iparticle) = efield2(1,iparticle) - (1.0_dp/(18.0_dp))*q_t
                      efield2(5,iparticle) = efield2(5,iparticle) - (1.0_dp/(18.0_dp))*q_t
                      efield2(9,iparticle) = efield2(9,iparticle) - (1.0_dp/(18.0_dp))*q_t
                   END IF
                END DO
             END DO
          END IF
          CYCLE
       END IF
       gauss  = (rho0(lp,mp,np) * pw_grid%vol)**2 / pw_grid%gsq(gpt)
       factor = gauss * REAL(summe_ef(gpt) * CONJG(summe_ef(gpt)),KIND=dp)
       energy = energy + factor

       IF (do_forces.OR.do_efield) THEN
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
                      tmp = AIMAG(fac*CONJG(summe_ef(gpt)))*gauss
                      efield1(1,iparticle) = efield1(1,iparticle) - tmp * pw_grid%g(1,gpt)
                      efield1(2,iparticle) = efield1(2,iparticle) - tmp * pw_grid%g(2,gpt)
                      efield1(3,iparticle) = efield1(3,iparticle) - tmp * pw_grid%g(3,gpt)
                   END IF
                   ! Electric field gradient
                   IF (do_efield2) THEN
                      tmp_v(1) = REAL(fac*CONJG(summe_ef(gpt)),KIND=dp)*pw_grid%g(1,gpt)*gauss
                      tmp_v(2) = REAL(fac*CONJG(summe_ef(gpt)),KIND=dp)*pw_grid%g(2,gpt)*gauss
                      tmp_v(3) = REAL(fac*CONJG(summe_ef(gpt)),KIND=dp)*pw_grid%g(3,gpt)*gauss

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
       END IF

       ! Compute the virial P*V
       IF (do_stress) THEN
          ! The Stress Tensor can be decomposed in two main components.
          ! The first one is just a normal ewald component for reciprocal space
          denom = 1.0_dp/four_alpha_sq + 1.0_dp/pw_grid%gsq(gpt)
          pv_tmp(1,1) = pv_tmp(1,1) + factor*(1.0_dp - 2.0_dp*pw_grid%g(1,gpt)*pw_grid%g(1,gpt)*denom)
          pv_tmp(1,2) = pv_tmp(1,2) - factor*(         2.0_dp*pw_grid%g(1,gpt)*pw_grid%g(2,gpt)*denom)
          pv_tmp(1,3) = pv_tmp(1,3) - factor*(         2.0_dp*pw_grid%g(1,gpt)*pw_grid%g(3,gpt)*denom)
          pv_tmp(2,1) = pv_tmp(2,1) - factor*(         2.0_dp*pw_grid%g(2,gpt)*pw_grid%g(1,gpt)*denom)
          pv_tmp(2,2) = pv_tmp(2,2) + factor*(1.0_dp - 2.0_dp*pw_grid%g(2,gpt)*pw_grid%g(2,gpt)*denom)
          pv_tmp(2,3) = pv_tmp(2,3) - factor*(         2.0_dp*pw_grid%g(2,gpt)*pw_grid%g(3,gpt)*denom)
          pv_tmp(3,1) = pv_tmp(3,1) - factor*(         2.0_dp*pw_grid%g(3,gpt)*pw_grid%g(1,gpt)*denom)
          pv_tmp(3,2) = pv_tmp(3,2) - factor*(         2.0_dp*pw_grid%g(3,gpt)*pw_grid%g(2,gpt)*denom)
          pv_tmp(3,3) = pv_tmp(3,3) + factor*(1.0_dp - 2.0_dp*pw_grid%g(3,gpt)*pw_grid%g(3,gpt)*denom)
          ! The second one can be written in the following way
          f0 = 2.0_dp * gauss
          pv_tmp(1,1) = pv_tmp(1,1) + f0 * pw_grid%g(1,gpt) * REAL(summe_st(1,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv_tmp(1,2) = pv_tmp(1,2) + f0 * pw_grid%g(1,gpt) * REAL(summe_st(2,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv_tmp(1,3) = pv_tmp(1,3) + f0 * pw_grid%g(1,gpt) * REAL(summe_st(3,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv_tmp(2,1) = pv_tmp(2,1) + f0 * pw_grid%g(2,gpt) * REAL(summe_st(1,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv_tmp(2,2) = pv_tmp(2,2) + f0 * pw_grid%g(2,gpt) * REAL(summe_st(2,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv_tmp(2,3) = pv_tmp(2,3) + f0 * pw_grid%g(2,gpt) * REAL(summe_st(3,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv_tmp(3,1) = pv_tmp(3,1) + f0 * pw_grid%g(3,gpt) * REAL(summe_st(1,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv_tmp(3,2) = pv_tmp(3,2) + f0 * pw_grid%g(3,gpt) * REAL(summe_st(2,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
          pv_tmp(3,3) = pv_tmp(3,3) + f0 * pw_grid%g(3,gpt) * REAL(summe_st(3,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
       END IF
    END DO
    pref   = fourpi/pw_grid%vol
    energy = energy * pref

    CALL structure_factor_deallocate (exp_igr)
    DEALLOCATE (summe_ef, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (do_stress) THEN
       pv_tmp = pv_tmp * pref
       ! Symmetrize the tensor
       pv(1,1) = pv(1,1) +  pv_tmp(1,1)
       pv(1,2) = pv(1,2) + (pv_tmp(1,2) + pv_tmp(2,1))*0.5_dp
       pv(1,3) = pv(1,3) + (pv_tmp(1,3) + pv_tmp(3,1))*0.5_dp
       pv(2,1) = pv(1,2)
       pv(2,2) = pv(2,2) +  pv_tmp(2,2)
       pv(2,3) = pv(2,3) + (pv_tmp(2,3) + pv_tmp(3,2))*0.5_dp
       pv(3,1) = pv(1,3)
       pv(3,2) = pv(2,3)
       pv(3,3) = pv(3,3) +  pv_tmp(3,3)
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
    CALL timestop(handle)

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
  SUBROUTINE ewald_multipole_self (ewald_env, cell, local_particles, e_self, &
       e_neut, task, do_efield, radii, charges, dipoles, quadrupoles, efield0, &
       efield1, efield2, error)
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cell_type), INTENT(IN)              :: cell
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(OUT)               :: e_self, e_neut
    LOGICAL, DIMENSION(3, 3), INTENT(IN)     :: task
    LOGICAL, INTENT(IN)                      :: do_efield
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: radii, charges
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
                                                do_efield2, lradii
    REAL(KIND=dp) :: alpha, ch_qu_self, ch_qu_self_tmp, dipole_self, fac1, &
      fac2, fac3, fac4, q, q_neutg, q_self, q_sum, qu_qu_self, radius

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
    fac2        = 6.0_dp*(f23**2)*(alpha**3)*oorootpi
    fac3        = (2.0_dp*oorootpi)*f23*alpha**3
    fac4        = (4.0_dp*oorootpi)*f415*alpha**5
    lradii      = PRESENT(radii)
    radius      = 0.0_dp
    q_neutg     = 0.0_dp
    DO iparticle_kind=1,SIZE(local_particles%n_el)
       nparticle_local = local_particles%n_el(iparticle_kind)
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          IF (ANY(task(1,:))) THEN
             ! Charge - Charge
             q = charges(iparticle)
             IF (lradii) radius = radii(iparticle)
             IF (radius>0) THEN
                q_neutg = q_neutg + 2.0_dp*q*radius**2
             END IF
             q_self = q_self + q*q
             q_sum  = q_sum  + q
             ! Potential
             IF (do_efield0) THEN
                efield0(  iparticle) = efield0(  iparticle) - q*fac1
             END IF

             IF (task(1,3)) THEN
                ! Charge - Quadrupole
                ch_qu_self_tmp = 0.0_dp
                DO i = 1,3
                   ch_qu_self_tmp = ch_qu_self_tmp + quadrupoles(i,i,iparticle) * q
                END DO
                ch_qu_self = ch_qu_self + ch_qu_self_tmp
                ! Electric Field Gradient
                IF (do_efield2) THEN
                   efield2(1,iparticle) = efield2(1,iparticle) + fac2*q
                   efield2(5,iparticle) = efield2(5,iparticle) + fac2*q
                   efield2(9,iparticle) = efield2(9,iparticle) + fac2*q
                END IF
             END IF
          END IF
          IF (ANY(task(2,:))) THEN
             ! Dipole - Dipole
             DO i = 1,3
                dipole_self = dipole_self + dipoles(i,iparticle)**2
             END DO
             ! Electric Field
             IF (do_efield1) THEN
                efield1(1,iparticle) = efield1(1,iparticle) + fac3 * dipoles(1,iparticle)
                efield1(2,iparticle) = efield1(2,iparticle) + fac3 * dipoles(2,iparticle)
                efield1(3,iparticle) = efield1(3,iparticle) + fac3 * dipoles(3,iparticle)
             END IF
          END IF
          IF (ANY(task(3,:))) THEN
             ! Quadrupole - Quadrupole
             DO i = 1,3
                DO j = 1,3
                   qu_qu_self = qu_qu_self + quadrupoles(j,i,iparticle)**2
                END DO
             END DO
             ! Electric Field Gradient
             IF (do_efield2) THEN
                efield2(1,iparticle) = efield2(1,iparticle) + fac4 * quadrupoles(1,1,iparticle)
                efield2(2,iparticle) = efield2(2,iparticle) + fac4 * quadrupoles(2,1,iparticle)
                efield2(3,iparticle) = efield2(3,iparticle) + fac4 * quadrupoles(3,1,iparticle)
                efield2(4,iparticle) = efield2(4,iparticle) + fac4 * quadrupoles(1,2,iparticle)
                efield2(5,iparticle) = efield2(5,iparticle) + fac4 * quadrupoles(2,2,iparticle)
                efield2(6,iparticle) = efield2(6,iparticle) + fac4 * quadrupoles(3,2,iparticle)
                efield2(7,iparticle) = efield2(7,iparticle) + fac4 * quadrupoles(1,3,iparticle)
                efield2(8,iparticle) = efield2(8,iparticle) + fac4 * quadrupoles(2,3,iparticle)
                efield2(9,iparticle) = efield2(9,iparticle) + fac4 * quadrupoles(3,3,iparticle)
             END IF
          END IF
       END DO
    END DO

    CALL mp_sum (q_neutg, group)
    CALL mp_sum (q_self, group)
    CALL mp_sum (q_sum, group)
    CALL mp_sum (dipole_self, group)
    CALL mp_sum (ch_qu_self, group)
    CALL mp_sum (qu_qu_self, group)

    e_self = -(q_self+f23*(dipole_self-f23*ch_qu_self+f415*qu_qu_self*alpha**2)*alpha**2)*alpha*oorootpi
    fac1 = pi/(2.0_dp*cell%deth)
    e_neut = -q_sum*fac1*(q_sum/alpha**2 - q_neutg)

    ! Correcting Potential for the neutralizing background charge
    DO iparticle_kind=1,SIZE(local_particles%n_el)
       nparticle_local = local_particles%n_el(iparticle_kind)
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          IF (ANY(task(1,:))) THEN
             ! Potential energy
             IF (do_efield0) THEN
                efield0(  iparticle) = efield0(  iparticle) - q_sum*2.0_dp*fac1/alpha**2
                IF (lradii) radius = radii(iparticle)
                IF (radius>0) THEN
                   q = charges(iparticle)
                   efield0(  iparticle) = efield0(  iparticle) + fac1*radius**2*(q_sum+q)
                END IF
             END IF
          END IF
       END DO
    END DO
  END SUBROUTINE ewald_multipole_self

! *****************************************************************************
!> \author Teodoro Laino [tlaino] - University of Zurich - 12.2007
! *****************************************************************************
  SUBROUTINE ewald_multipole_print(iw, e_gspace, e_rspace, e_bonded, e_self, e_neut)

    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), INTENT(IN)                :: e_gspace, e_rspace, e_bonded, &
                                                e_self, e_neut

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_print', &
      routineP = moduleN//':'//routineN

    IF (iw>0) THEN
       WRITE ( iw, '( A, A )' ) ' *********************************', &
            '**********************************************'
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' INITIAL GSPACE ENERGY', &
            '[hartree]', '= ', e_gspace
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' INITIAL RSPACE ENERGY', &
            '[hartree]', '= ', e_rspace
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' BONDED CORRECTION', &
            '[hartree]', '= ', e_bonded
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' SELF ENERGY CORRECTION', &
            '[hartree]', '= ', e_self
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' NEUTRALIZ. BCKGR. ENERGY', &
            '[hartree]', '= ', e_neut
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' TOTAL ELECTROSTATIC EN.', &
            '[hartree]', '= ', e_rspace+e_bonded+e_gspace+e_self+e_neut
       WRITE ( iw, '( A, A )' ) ' *********************************', &
            '**********************************************'
    END IF
  END SUBROUTINE ewald_multipole_print

END MODULE ewalds_multipole
