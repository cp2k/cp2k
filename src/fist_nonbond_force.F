!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH (11 May 2001) : cleaning up of support structures
!>      CJM & HAF (27 July 2001): fixed bug with handling of cutoff larger than
!>                                half the boxsize.
!>      07.02.2005: getting rid of scaled_to_real calls in force loop (MK)
!> \author CJM
! *****************************************************************************
MODULE fist_nonbond_force

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE f77_blas
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: oorootpi
  USE mathlib,                         ONLY: matvec_3x3
  USE message_passing,                 ONLY: mp_sum
  USE pair_potential_types,            ONLY: nosh_nosh,&
                                             nosh_sh,&
                                             pair_potential_pp_type,&
                                             pair_potential_single_type,&
                                             sh_sh
  USE particle_types,                  ONLY: particle_type
  USE shell_potential_types,           ONLY: shell_kind_type
  USE splines_methods,                 ONLY: potential_s
  USE splines_types,                   ONLY: spline_data_p_type,&
                                             spline_factor_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_nonbond_force'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.

  PUBLIC :: force_nonbond,&
            bonded_correct_gaussian

CONTAINS

! *****************************************************************************
!> \brief Calculates the force and the potential of the minimum image, and
!>      the pressure tensor
! *****************************************************************************
  SUBROUTINE force_nonbond ( fist_nonbond_env, particle_set, cell, pot_nonbond,&
       f_nonbond, pv_nonbond, fshell_nonbond, fcore_nonbond, &
       atomic_kind_set, use_virial, error )

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(OUT)               :: pot_nonbond
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_nonbond, pv_nonbond
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: fshell_nonbond, fcore_nonbond
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set(:)
    LOGICAL, INTENT(IN)                      :: use_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_nonbond', &
      routineP = moduleN//':'//routineN

    INTEGER :: a_kind, atom_a, atom_b, b_kind, handle, i, iend, igrp, ikind, &
      ilist, ipair, istart, jkind, nkinds, npairs, shell_a, shell_b, &
      shell_type
    INTEGER, DIMENSION(:, :), POINTER        :: list
    LOGICAL                                  :: all_terms, do_shell, failure, &
                                                shell_present, swap_ab
    REAL(KIND=dp) :: energy, fac, fscalar, ptens11, ptens12, ptens13, &
      ptens21, ptens22, ptens23, ptens31, ptens32, ptens33, rab2, rab2_com, &
      rab2_max
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, cvi, fr, rab, rab_cc, &
                                                rab_com, rab_cs, rab_sc, &
                                                rab_ss
    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_com
    REAL(KIND=dp), DIMENSION(3, 4)           :: rab_list
    REAL(KIND=dp), DIMENSION(4)              :: rab2_list
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ij_kind_full_fac
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update, &
                                                r_last_update_pbc, &
                                                rcore_last_update_pbc, &
                                                rshell_last_update_pbc
    TYPE(spline_data_p_type), DIMENSION(:), &
      POINTER                                :: c_spline_data, &
                                                cs_spline_data, &
                                                sc_spline_data, spline_data
    TYPE(spline_factor_type), POINTER        :: spl_c, spl_cs, spl_f, spl_sc

    CALL timeset ( routineN, handle )
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    NULLIFY(pot, c_spline_data, cs_spline_data, sc_spline_data, rshell_last_update_pbc)
    NULLIFY(spl_f, spl_c, spl_cs, spl_sc, ij_kind_full_fac)
    CALL fist_nonbond_env_get ( fist_nonbond_env, nonbonded = nonbonded, potparm = potparm,&
         r_last_update=r_last_update, r_last_update_pbc=r_last_update_pbc,natom_types = nkinds,&
         rshell_last_update_pbc=rshell_last_update_pbc,rcore_last_update_pbc=rcore_last_update_pbc,&
         ij_kind_full_fac=ij_kind_full_fac, error=error)

    ! Initializing the potential energy, pressure tensor and force
    pot_nonbond        = 0.0_dp
    f_nonbond(:,:)     = 0.0_dp
    
    IF (use_virial) THEN
       pv_nonbond(:,:) = 0.0_dp
       ptens11 = 0.0_dp ; ptens12 = 0.0_dp ; ptens13 = 0.0_dp
       ptens21 = 0.0_dp ; ptens22 = 0.0_dp ; ptens23 = 0.0_dp
       ptens31 = 0.0_dp ; ptens32 = 0.0_dp ; ptens33 = 0.0_dp
    END IF
    shell_present = .FALSE.
    IF(PRESENT(fshell_nonbond)) THEN
       CPPostcondition(PRESENT(fcore_nonbond),cp_failure_level,routineP,error,failure)
       fshell_nonbond = 0.0_dp
       fcore_nonbond = 0.0_dp
       shell_present = .TRUE.
    END IF
    ! Starting the force loop
! incorrect omp (i.e. threads can write to the variable an other thread might be reading (e.g. pot_nonbonded))
! !$OMP parallel do default(none) &
! !$OMP shared(nkinds) &
! !$OMP private(ikind,jkind,id_kind,istart,iend) &
! !$OMP private(spline_data) &
! !$OMP shared(potparm,use_virial) &
! !$OMP private(ab) &
! !$OMP shared(nonbonded) &
! !$OMP private(rab,rab2_max,fac) &
! !$OMP private(ilist,nlist,igrp,atom_a,atom_b,npairs,cvi,list,inode,nnode) &
! !$OMP shared(cell) &
! !$OMP shared(pot_nonbond) &
! !$OMP private(f_eam,pot,fr,fscalar,energy) &
! !$OMP shared(f_nonbond,pv_nonbond,r_last_update,ij_kind_full_fac) &
! !$OMP private(rab2) &
! !$OMP private(neighbor_kind_pair,cell_v) &
! !$OMP shared(particle_set,r_last_update_pbc) &
! !$OMP reduction(+:ptens11,ptens12,ptens13,ptens21,ptens22,ptens23, &
! !$OMP             ptens31,ptens32,ptens33) &
! !$OMP private(ipair) 
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
          fac     = ij_kind_full_fac(ikind, jkind)
          pot     => potparm%pot(ikind,jkind)%pot
          ! In case we have manybody potentials only they will be handled in the manybody routine...
          IF (pot%only_mb) CYCLE
          ! Setup spline_data set
          spl_f       => pot%spl_f
          spline_data => pot%pair_spline_data
          shell_type  =  pot%shell_type 
          IF(shell_type/=nosh_nosh) THEN
             c_spline_data => pot%c_coul_spline_data
             spl_c         => pot%spl_c
             CPPostcondition(shell_present,cp_failure_level,routineP,error,failure)
             do_shell = .TRUE.
          ELSE
             c_spline_data => pot%pair_spline_data
             spl_c         => pot%spl_f
             do_shell = .FALSE.
          END IF
          IF(shell_type==sh_sh) THEN
             cs_spline_data => pot%cs_coul_spline_data
             sc_spline_data => pot%sc_coul_spline_data
             spl_cs         => pot%spl_cs
             spl_sc         => pot%spl_sc
          END IF
          rab2_max = pot%rcutsq
          Pairs: DO ipair = istart, iend
             atom_a = list(1,ipair)
             atom_b = list(2,ipair)
             IF (do_shell) THEN
                all_terms = .TRUE.
                a_kind  = particle_set(atom_a)%atomic_kind%kind_number
                b_kind  = particle_set(atom_b)%atomic_kind%kind_number
                swap_ab = (b_kind>a_kind)
                IF (swap_ab) THEN
                   cs_spline_data => pot%sc_coul_spline_data
                   sc_spline_data => pot%cs_coul_spline_data
                   spl_cs         => pot%spl_sc
                   spl_sc         => pot%spl_cs
                END IF
                IF(shell_type==sh_sh) THEN
                   shell_a = particle_set(atom_a)%shell_index
                   shell_b = particle_set(atom_b)%shell_index
                   rab_cc =  rcore_last_update_pbc(shell_b)%r -  rcore_last_update_pbc(shell_a)%r
                   rab_cs = rshell_last_update_pbc(shell_b)%r -  rcore_last_update_pbc(shell_a)%r
                   rab_sc =  rcore_last_update_pbc(shell_b)%r - rshell_last_update_pbc(shell_a)%r
                   rab_ss = rshell_last_update_pbc(shell_b)%r - rshell_last_update_pbc(shell_a)%r
                   rab_list(1:3,1) = rab_cc(1:3)+cell_v(1:3)
                   rab_list(1:3,2) = rab_cs(1:3)+cell_v(1:3)
                   rab_list(1:3,3) = rab_sc(1:3)+cell_v(1:3)
                   rab_list(1:3,4) = rab_ss(1:3)+cell_v(1:3)
                ELSEIF(shell_type==nosh_sh .AND. particle_set(atom_a)%shell_index/=0) THEN
                   shell_a = particle_set(atom_a)%shell_index
                   shell_b = 0
                   rab_cc = r_last_update_pbc(atom_b)%r -  rcore_last_update_pbc(shell_a)%r
                   rab_sc = 0.0_dp
                   rab_cs = 0.0_dp
                   rab_ss = r_last_update_pbc(atom_b)%r - rshell_last_update_pbc(shell_a)%r
                   rab_list(1:3,1) = rab_cc(1:3)+cell_v(1:3)
                   rab_list(1:3,2) = 0.0_dp
                   rab_list(1:3,3) = 0.0_dp
                   rab_list(1:3,4) = rab_ss(1:3)+cell_v(1:3)
                ELSEIF(shell_type==nosh_sh .AND. particle_set(atom_b)%shell_index/=0) THEN
                   shell_b = particle_set(atom_b)%shell_index
                   shell_a = 0
                   rab_cc =  rcore_last_update_pbc(shell_b)%r - r_last_update_pbc(atom_a)%r
                   rab_sc = 0.0_dp
                   rab_cs = 0.0_dp
                   rab_ss = rshell_last_update_pbc(shell_b)%r - r_last_update_pbc(atom_a)%r
                   rab_list(1:3,1) = rab_cc(1:3)+cell_v(1:3)
                   rab_list(1:3,2) = 0.0_dp
                   rab_list(1:3,3) = 0.0_dp
                   rab_list(1:3,4) = rab_ss(1:3)+cell_v(1:3)
                END IF
                ! Compute the term only if all the pairs (cc,cs,sc,ss) are within the cut-off
                Check_terms: DO i = 1,4
                   rab2_list(i) = rab_list(1,i)**2+rab_list(2,i)**2+rab_list(3,i)**2
                   IF (rab2_list(i) >= rab2_max) THEN
                      all_terms = .FALSE.
                      EXIT Check_terms
                   END IF
                END DO Check_terms
                rab_com  = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
             ELSE
                rab_cc   = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
                rab_com  = rab_cc
             END IF
             rab_com  = rab_com + cell_v
             rab2_com = rab_com(1)**2 + rab_com(2)**2 + rab_com(3)**2
             IF (do_shell) THEN
                IF (rab2_com <= rab2_max .AND. all_terms) THEN
                   ! SHELL MODEL
                   ! core-core or core-ion/ion-core: Coulomb only
                   rab = rab_list(:,1)
                   rab2 = rab2_list(1)
                   energy =  potential_s (c_spline_data,rab2,fscalar, spl_c, logger)
                   pot_nonbond = pot_nonbond + energy*fac
                   fscalar = fscalar*fac
                   fr(1) = fscalar*rab(1)
                   fr(2) = fscalar*rab(2)
                   fr(3) = fscalar*rab(3)
                   IF(shell_a /= 0 ) THEN
                      fcore_nonbond(1,shell_a) = fcore_nonbond(1,shell_a) - fr(1)
                      fcore_nonbond(2,shell_a) = fcore_nonbond(2,shell_a) - fr(2)
                      fcore_nonbond(3,shell_a) = fcore_nonbond(3,shell_a) - fr(3)
                   ELSE
                      f_nonbond(1,atom_a) = f_nonbond(1,atom_a) - fr(1)
                      f_nonbond(2,atom_a) = f_nonbond(2,atom_a) - fr(2)
                      f_nonbond(3,atom_a) = f_nonbond(3,atom_a) - fr(3)
                      CPPostcondition(shell_b/=0,cp_failure_level,routineP,error,failure)
                   END IF
                   
                   IF(shell_b /= 0 ) THEN
                      fcore_nonbond(1,shell_b) = fcore_nonbond(1,shell_b) + fr(1)
                      fcore_nonbond(2,shell_b) = fcore_nonbond(2,shell_b) + fr(2)
                      fcore_nonbond(3,shell_b) = fcore_nonbond(3,shell_b) + fr(3)
                   ELSE
                      f_nonbond(1,atom_b) = f_nonbond(1,atom_b) + fr(1)
                      f_nonbond(2,atom_b) = f_nonbond(2,atom_b) + fr(2)
                      f_nonbond(3,atom_b) = f_nonbond(3,atom_b) + fr(3)
                      CPPostcondition(shell_a/=0,cp_failure_level,routineP,error,failure)
                   END IF
                   IF (use_virial) THEN
                      pv_com(1,1) = rab(1) * fr(1)
                      pv_com(1,2) = rab(1) * fr(2)
                      pv_com(1,3) = rab(1) * fr(3)
                      pv_com(2,1) = rab(2) * fr(1)
                      pv_com(2,2) = rab(2) * fr(2)
                      pv_com(2,3) = rab(2) * fr(3)
                      pv_com(3,1) = rab(3) * fr(1)
                      pv_com(3,2) = rab(3) * fr(2)
                      pv_com(3,3) = rab(3) * fr(3)
                   END IF
          
                   IF(shell_type==sh_sh) THEN
                      ! shell-shell : VDW + Coulomb
                      rab  = rab_list(:,4)
                      rab2 = rab2_list(4)
                      energy =  potential_s (spline_data,rab2,fscalar,spl_f,logger)
                      pot_nonbond = pot_nonbond + energy*fac
                      fscalar = fscalar * fac
                      fr(1) = fscalar*rab(1)
                      fr(2) = fscalar*rab(2)
                      fr(3) = fscalar*rab(3)
                      fshell_nonbond(1,shell_a) = fshell_nonbond(1,shell_a) - fr(1)
                      fshell_nonbond(2,shell_a) = fshell_nonbond(2,shell_a) - fr(2)
                      fshell_nonbond(3,shell_a) = fshell_nonbond(3,shell_a) - fr(3)
                      fshell_nonbond(1,shell_b) = fshell_nonbond(1,shell_b) + fr(1)
                      fshell_nonbond(2,shell_b) = fshell_nonbond(2,shell_b) + fr(2)
                      fshell_nonbond(3,shell_b) = fshell_nonbond(3,shell_b) + fr(3)
                      IF (use_virial) THEN
                         pv_com(1,1) = pv_com(1,1) + rab(1) * fr(1)
                         pv_com(1,2) = pv_com(1,2) + rab(1) * fr(2)
                         pv_com(1,3) = pv_com(1,3) + rab(1) * fr(3)
                         pv_com(2,1) = pv_com(2,1) + rab(2) * fr(1)
                         pv_com(2,2) = pv_com(2,2) + rab(2) * fr(2)
                         pv_com(2,3) = pv_com(2,3) + rab(2) * fr(3)
                         pv_com(3,1) = pv_com(3,1) + rab(3) * fr(1)
                         pv_com(3,2) = pv_com(3,2) + rab(3) * fr(2)
                         pv_com(3,3) = pv_com(3,3) + rab(3) * fr(3)
                      END IF
          
                      ! core-shell : Coulomb only
                      rab  = rab_list(:,2)
                      rab2 = rab2_list(2)
                      energy =  potential_s (cs_spline_data,rab2,fscalar,spl_cs,logger)
                      pot_nonbond = pot_nonbond + energy*fac
                      fscalar = fscalar*fac
                      fr(1) = fscalar*rab(1)
                      fr(2) = fscalar*rab(2)
                      fr(3) = fscalar*rab(3)
                      fcore_nonbond(1,shell_a) = fcore_nonbond(1,shell_a) - fr(1)
                      fcore_nonbond(2,shell_a) = fcore_nonbond(2,shell_a) - fr(2)
                      fcore_nonbond(3,shell_a) = fcore_nonbond(3,shell_a) - fr(3)
                      fshell_nonbond(1,shell_b) = fshell_nonbond(1,shell_b) + fr(1)
                      fshell_nonbond(2,shell_b) = fshell_nonbond(2,shell_b) + fr(2)
                      fshell_nonbond(3,shell_b) = fshell_nonbond(3,shell_b) + fr(3)
                      IF (use_virial) THEN
                         pv_com(1,1) = pv_com(1,1) + rab(1) * fr(1)
                         pv_com(1,2) = pv_com(1,2) + rab(1) * fr(2)
                         pv_com(1,3) = pv_com(1,3) + rab(1) * fr(3)
                         pv_com(2,1) = pv_com(2,1) + rab(2) * fr(1)
                         pv_com(2,2) = pv_com(2,2) + rab(2) * fr(2)
                         pv_com(2,3) = pv_com(2,3) + rab(2) * fr(3)
                         pv_com(3,1) = pv_com(3,1) + rab(3) * fr(1)
                         pv_com(3,2) = pv_com(3,2) + rab(3) * fr(2)
                         pv_com(3,3) = pv_com(3,3) + rab(3) * fr(3)
                      END IF
                      
                      ! shell-core : Coulomb only
                      rab  = rab_list(:,3)
                      rab2 = rab2_list(3)
                      energy =  potential_s (sc_spline_data,rab2,fscalar,spl_sc,logger)
                      pot_nonbond = pot_nonbond + energy*fac
                      fscalar = fscalar * fac
                      fr(1) = fscalar*rab(1)
                      fr(2) = fscalar*rab(2)
                      fr(3) = fscalar*rab(3)
                      fshell_nonbond(1,shell_a) = fshell_nonbond(1,shell_a) - fr(1)
                      fshell_nonbond(2,shell_a) = fshell_nonbond(2,shell_a) - fr(2)
                      fshell_nonbond(3,shell_a) = fshell_nonbond(3,shell_a) - fr(3)
                      fcore_nonbond(1,shell_b) = fcore_nonbond(1,shell_b) + fr(1)
                      fcore_nonbond(2,shell_b) = fcore_nonbond(2,shell_b) + fr(2)
                      fcore_nonbond(3,shell_b) = fcore_nonbond(3,shell_b) + fr(3)
                      IF (use_virial) THEN
                         pv_com(1,1) = pv_com(1,1) + rab(1) * fr(1)
                         pv_com(1,2) = pv_com(1,2) + rab(1) * fr(2)
                         pv_com(1,3) = pv_com(1,3) + rab(1) * fr(3)
                         pv_com(2,1) = pv_com(2,1) + rab(2) * fr(1)
                         pv_com(2,2) = pv_com(2,2) + rab(2) * fr(2)
                         pv_com(2,3) = pv_com(2,3) + rab(2) * fr(3)
                         pv_com(3,1) = pv_com(3,1) + rab(3) * fr(1)
                         pv_com(3,2) = pv_com(3,2) + rab(3) * fr(2)
                         pv_com(3,3) = pv_com(3,3) + rab(3) * fr(3)
                      END IF
                   ELSEIF(shell_type==nosh_sh) THEN
                      ! shell-ion or ion-shell : VDW + Coulomb
                      rab = rab_list(:,4)
                      rab2 = rab2_list(4)
                      energy =  potential_s (spline_data,rab2,fscalar,spl_f,logger)
                      pot_nonbond = pot_nonbond + energy*fac
                      fscalar = fscalar * fac
                      fr(1) = fscalar*rab(1)
                      fr(2) = fscalar*rab(2)
                      fr(3) = fscalar*rab(3)
                      IF(shell_a/=0) THEN
                         fshell_nonbond(1,shell_a) = fshell_nonbond(1,shell_a) - fr(1)
                         fshell_nonbond(2,shell_a) = fshell_nonbond(2,shell_a) - fr(2)
                         fshell_nonbond(3,shell_a) = fshell_nonbond(3,shell_a) - fr(3)
                         f_nonbond(1,atom_b) = f_nonbond(1,atom_b) + fr(1)
                         f_nonbond(2,atom_b) = f_nonbond(2,atom_b) + fr(2)
                         f_nonbond(3,atom_b) = f_nonbond(3,atom_b) + fr(3)
                      ELSE
                         f_nonbond(1,atom_a) = f_nonbond(1,atom_a) - fr(1)
                         f_nonbond(2,atom_a) = f_nonbond(2,atom_a) - fr(2)
                         f_nonbond(3,atom_a) = f_nonbond(3,atom_a) - fr(3)
                         fshell_nonbond(1,shell_b) = fshell_nonbond(1,shell_b) + fr(1)
                         fshell_nonbond(2,shell_b) = fshell_nonbond(2,shell_b) + fr(2)
                         fshell_nonbond(3,shell_b) = fshell_nonbond(3,shell_b) + fr(3)
                      END IF
                      IF (use_virial) THEN
                         pv_com(1,1) = pv_com(1,1) + rab(1) * fr(1)
                         pv_com(1,2) = pv_com(1,2) + rab(1) * fr(2)
                         pv_com(1,3) = pv_com(1,3) + rab(1) * fr(3)
                         pv_com(2,1) = pv_com(2,1) + rab(2) * fr(1)
                         pv_com(2,2) = pv_com(2,2) + rab(2) * fr(2)
                         pv_com(2,3) = pv_com(2,3) + rab(2) * fr(3)
                         pv_com(3,1) = pv_com(3,1) + rab(3) * fr(1)
                         pv_com(3,2) = pv_com(3,2) + rab(3) * fr(2)
                         pv_com(3,3) = pv_com(3,3) + rab(3) * fr(3)
                      END IF
                   END IF
                   IF (use_virial) THEN
                      ptens11 = ptens11 + pv_com(1,1)
                      ptens21 = ptens21 + pv_com(2,1)
                      ptens31 = ptens31 + pv_com(3,1)
                      ptens12 = ptens12 + pv_com(1,2)
                      ptens22 = ptens22 + pv_com(2,2)
                      ptens32 = ptens32 + pv_com(3,2)
                      ptens13 = ptens13 + pv_com(1,3)
                      ptens23 = ptens23 + pv_com(2,3)
                      ptens33 = ptens33 + pv_com(3,3)
                   END IF
                   ! Let's restore the proper order to save time...
                   IF (swap_ab) THEN
                      cs_spline_data => pot%cs_coul_spline_data
                      sc_spline_data => pot%sc_coul_spline_data
                      spl_cs         => pot%spl_cs
                      spl_sc         => pot%spl_sc
                   END IF
                ELSE
                   IF (swap_ab) THEN
                      cs_spline_data => pot%cs_coul_spline_data
                      sc_spline_data => pot%sc_coul_spline_data
                      spl_cs         => pot%spl_cs
                      spl_sc         => pot%spl_sc
                   END IF
                END IF
             ELSE
                IF (rab2_com <= rab2_max) THEN
                   ! NO SHELL MODEL...
                   ! Ion-Ion : no shell model
                   rab  = rab_com
                   rab2 = rab2_com
                   energy =  potential_s(c_spline_data,rab2,fscalar,spl_c,logger)
                   pot_nonbond = pot_nonbond + energy*fac
                   fscalar = fscalar*fac
                   fr(1) = fscalar*rab(1)
                   fr(2) = fscalar*rab(2)
                   fr(3) = fscalar*rab(3)
                   f_nonbond(1,atom_a) = f_nonbond(1,atom_a) - fr(1)
                   f_nonbond(2,atom_a) = f_nonbond(2,atom_a) - fr(2)
                   f_nonbond(3,atom_a) = f_nonbond(3,atom_a) - fr(3)
                   f_nonbond(1,atom_b) = f_nonbond(1,atom_b) + fr(1)
                   f_nonbond(2,atom_b) = f_nonbond(2,atom_b) + fr(2)
                   f_nonbond(3,atom_b) = f_nonbond(3,atom_b) + fr(3)
                   IF (use_virial) THEN
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
             END IF
          END DO Pairs
       END DO Kind_Group_Loop
    END DO Lists
    IF (use_virial) THEN
       pv_nonbond(1,1) = pv_nonbond(1,1) + ptens11
       pv_nonbond(1,2) = pv_nonbond(1,2) + ptens12
       pv_nonbond(1,3) = pv_nonbond(1,3) + ptens13
       pv_nonbond(2,1) = pv_nonbond(2,1) + ptens21
       pv_nonbond(2,2) = pv_nonbond(2,2) + ptens22
       pv_nonbond(2,3) = pv_nonbond(2,3) + ptens23
       pv_nonbond(3,1) = pv_nonbond(3,1) + ptens31
       pv_nonbond(3,2) = pv_nonbond(3,2) + ptens32
       pv_nonbond(3,3) = pv_nonbond(3,3) + ptens33
    END IF

    CALL timestop ( handle )
  END SUBROUTINE force_nonbond

! *****************************************************************************
!> \brief corrects electrostatics for bonded terms
!> \par History
!>      Splitted routines to clean and to fix a bug with the tensor whose 
!>      original definition was not correct for PBC.. [Teodoro Laino -06/2007]
! *****************************************************************************
  SUBROUTINE bonded_correct_gaussian ( atomic_kind_set, local_particles, &
             particle_set, ewald_env, v_bonded_corr, pv_bc, shell_particle_set, &
             core_particle_set, cell, use_virial, error )

    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    REAL(KIND=dp), INTENT(OUT)               :: v_bonded_corr
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: pv_bc
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set( : ), &
                                                core_particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(IN)                      :: use_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'bonded_correct_gaussian', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, handle, i, iatom, &
                                                iexl, ikind, j, &
                                                natoms_per_kind, nkind, &
                                                shell_i_index, shell_j_index
    LOGICAL                                  :: do_shell, i_is_shell, &
                                                j_is_shell, shell_adiabatic
    REAL(KIND=dp)                            :: alpha, const, qci, qcj, qsi, &
                                                qsj
    REAL(KIND=dp), DIMENSION(3)              :: fij_com, rci, rcj, rsi, rsj

    CALL timeset(routineN,handle)
    ! Initializing values
    IF (use_virial) pv_bc = 0.0_dp
    v_bonded_corr = 0.0_dp
    ! Defining the constants
    CALL ewald_env_get ( ewald_env, alpha = alpha, group = group ,error=error)
    const = 2.0_dp * alpha * oorootpi

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             shell_adiabatic=shell_adiabatic)
    nkind = SIZE ( atomic_kind_set )
    DO ikind = 1, nkind
       natoms_per_kind = local_particles % n_el ( ikind )
       DO iatom = 1, natoms_per_kind
          i = local_particles % list ( ikind ) % array ( iatom )
          CALL get_atom_info(i, particle_set, shell_particle_set, core_particle_set,&
               shell_i_index, i_is_shell, rsi, rci, qci, qsi)

          DO iexl = 1, SIZE(particle_set(i)%list_exclude_ei)
             fij_com = 0.0_dp
             j = particle_set(i)%list_exclude_ei(iexl)
             IF ( j>=i ) CYCLE

             CALL get_atom_info(j, particle_set, shell_particle_set, core_particle_set,&
               shell_j_index, j_is_shell, rsj, rcj, qcj, qsj)
             do_shell = .FALSE.
             IF (i_is_shell.OR.j_is_shell) do_shell = .TRUE.
             !
             ! This term excludes always ION-ION interactions
             IF(i_is_shell) THEN
                CALL bonded_correct_gaussian_low(rci, rcj, cell, v_bonded_corr, &
                     core_particle_set, core_particle_set, shell_i_index, shell_j_index,&
                     .TRUE., alpha, qci, qcj, const, pv_bc, use_virial)
             ELSE
                CALL bonded_correct_gaussian_low(rci, rcj, cell, v_bonded_corr, particle_set,&
                     particle_set, i, j, .TRUE., alpha, qci, qcj, const, pv_bc,&
                     use_virial)
             END IF

             ! Shell-Model
             IF (do_shell) THEN
                ! shell-model: exclude shell_i-shell_j interactions
                IF(i_is_shell .AND. j_is_shell) THEN
                   CALL bonded_correct_gaussian_low(rsi, rsj, cell, v_bonded_corr,&
                        shell_particle_set, shell_particle_set, shell_i_index, shell_j_index,&
                        shell_adiabatic, alpha, qsi, qsj, const, pv_bc, use_virial)
                END IF
                ! shell-model: exclude shell_i-core_j interactions
                IF(i_is_shell ) THEN
                   IF (j_is_shell) THEN
                      CALL bonded_correct_gaussian_low(rsi, rcj, cell, v_bonded_corr,&
                           shell_particle_set, core_particle_set, shell_i_index, shell_j_index,&
                           shell_adiabatic, alpha, qsi, qcj, const, pv_bc, use_virial)
                   ELSE
                      CALL bonded_correct_gaussian_low(rsi, rcj, cell, v_bonded_corr,&
                           shell_particle_set, particle_set, shell_i_index, j, shell_adiabatic,&
                           alpha, qsi, qcj, const, pv_bc,use_virial)
                   END IF
                END IF
                ! shell-model: exclude shell_j-core_i interactions
                IF(j_is_shell ) THEN
                   IF (i_is_shell) THEN
                      CALL bonded_correct_gaussian_low(rci, rsj, cell, v_bonded_corr,&
                           core_particle_set, shell_particle_set, shell_i_index, shell_j_index,&
                           shell_adiabatic, alpha, qci, qsj, const, pv_bc, use_virial)
                   ELSE
                      CALL bonded_correct_gaussian_low(rci, rsj, cell, v_bonded_corr,&
                           particle_set, shell_particle_set, i, shell_j_index, shell_adiabatic,&
                           alpha, qci, qsj, const, pv_bc, use_virial)
                   END IF
                END IF
             END IF
          ENDDO

          ! Always Exclude  shell_i-core_i interaction
          IF(i_is_shell ) THEN
             CALL bonded_correct_gaussian_low_sh(rci, rsi, cell, v_bonded_corr,&
                  core_particle_set, shell_particle_set, shell_i_index, shell_adiabatic,&
                  alpha, qci, qsi, const, pv_bc, use_virial)
          END IF

       END DO
    ENDDO

    CALL mp_sum ( v_bonded_corr, group )
    CALL timestop(handle)

  END SUBROUTINE bonded_correct_gaussian

! *****************************************************************************
!> \par History
!>      Splitted routines to clean and to fix a bug with the tensor whose 
!>      original definition was not correct for PBC..
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE bonded_correct_gaussian_low(r1, r2, cell, v_bonded_corr,&
       particle_set1, particle_set2, i, j, shell_adiabatic, alpha, q1, q2, &
       const, pv_bc, use_virial)
    REAL(KIND=dp), DIMENSION(3)              :: r1, r2
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(INOUT)             :: v_bonded_corr
    TYPE(particle_type), POINTER             :: particle_set1(:), &
                                                particle_set2(:)
    INTEGER, INTENT(IN)                      :: i, j
    LOGICAL, INTENT(IN)                      :: shell_adiabatic
    REAL(KIND=dp), INTENT(IN)                :: alpha, q1, q2, const
    REAL(KIND=dp), INTENT(INOUT)             :: pv_bc(3,3)
    LOGICAL, INTENT(IN)                      :: use_virial

    CHARACTER(LEN=*), PARAMETER :: routineN = 'bonded_correct_gaussian_low', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER :: ac1 = 0.254829592_dp, ac2 = -0.284496736_dp, &
      ac3 = 1.421413741_dp, ac4 = -1.453152027_dp, ac5 = 1.061405429_dp, &
      pc = 0.3275911_dp

    REAL(KIND=dp)                            :: arg, dij, e_arg_arg, errf, &
                                                fscalar, idij, rijsq, tc
    REAL(KIND=dp), DIMENSION(3)              :: fij_com, rij

    rij = r1 - r2
    rij = pbc(rij, cell)
    rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
    idij = 1.0_dp / SQRT ( rijsq )
    dij = rijsq * idij
    arg = alpha * dij
    e_arg_arg = EXP ( -arg ** 2 )
    tc = 1.0_dp / ( 1.0_dp + pc * arg )

    ! Defining errf=1-erfc
    errf = 1.0_dp &
         - ((((ac5*tc+ac4)*tc+ac3)*tc+ac2)*tc+ac1) * tc * e_arg_arg

    ! Getting the potential
    v_bonded_corr = v_bonded_corr - q1*q2*idij*errf

    ! Subtracting the force from the total force
    fscalar = q1 * q2 * idij**2 * (idij*errf-const*e_arg_arg)

    particle_set1(i)%f(1) = particle_set1(i)%f(1) - fscalar*rij(1)
    particle_set1(i)%f(2) = particle_set1(i)%f(2) - fscalar*rij(2)
    particle_set1(i)%f(3) = particle_set1(i)%f(3) - fscalar*rij(3)

    particle_set2(j)%f(1) = particle_set2(j)%f(1) + fscalar*rij(1)
    particle_set2(j)%f(2) = particle_set2(j)%f(2) + fscalar*rij(2)
    particle_set2(j)%f(3) = particle_set2(j)%f(3) + fscalar*rij(3)

    IF (use_virial.AND.(shell_adiabatic)) THEN
       fij_com = fscalar*rij
       pv_bc(1,1) = pv_bc(1,1) - fij_com(1) * rij(1)
       pv_bc(1,2) = pv_bc(1,2) - fij_com(1) * rij(2)
       pv_bc(1,3) = pv_bc(1,3) - fij_com(1) * rij(3)
       pv_bc(2,1) = pv_bc(2,1) - fij_com(2) * rij(1)
       pv_bc(2,2) = pv_bc(2,2) - fij_com(2) * rij(2)
       pv_bc(2,3) = pv_bc(2,3) - fij_com(2) * rij(3)
       pv_bc(3,1) = pv_bc(3,1) - fij_com(3) * rij(1)
       pv_bc(3,2) = pv_bc(3,2) - fij_com(3) * rij(2)
       pv_bc(3,3) = pv_bc(3,3) - fij_com(3) * rij(3)
    END IF
  END SUBROUTINE bonded_correct_gaussian_low

! *****************************************************************************
!> \brief specific for shell models cleans the interaction core-shell on the same
!>      atom
!> \par History
!>      Splitted routines to clean and to fix a bug with the tensor whose 
!>      original definition was not correct for PBC..
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE bonded_correct_gaussian_low_sh(r1, r2, cell, v_bonded_corr, &
       core_particle_set, shell_particle_set, i, shell_adiabatic, alpha, q1, q2,&
       const, pv_bc, use_virial)
    REAL(KIND=dp), DIMENSION(3)              :: r1, r2
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(INOUT)             :: v_bonded_corr
    TYPE(particle_type), POINTER             :: core_particle_set(:), &
                                                shell_particle_set(:)
    INTEGER, INTENT(IN)                      :: i
    LOGICAL, INTENT(IN)                      :: shell_adiabatic
    REAL(KIND=dp), INTENT(IN)                :: alpha, q1, q2, const
    REAL(KIND=dp), INTENT(INOUT)             :: pv_bc(3,3)
    LOGICAL, INTENT(IN)                      :: use_virial

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'bonded_correct_gaussian_low_sh', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER :: ac1 = 0.254829592_dp, ac2 = -0.284496736_dp, &
      ac3 = 1.421413741_dp, ac4 = -1.453152027_dp, ac5 = 1.061405429_dp, &
      pc = 0.3275911_dp

    REAL(KIND=dp)                            :: arg, dij, e_arg_arg, efac, &
                                                errf, ffac, fscalar, idij, &
                                                rijsq, tc, tc2
    REAL(KIND=dp), DIMENSION(3)              :: fij_com, rij

    rij = r1 - r2
    rij = pbc(rij, cell)
    rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
    dij = SQRT(rijsq)
    ! Two possible limiting cases according the value of dij
    arg = alpha * dij
    ! and this is a magic number.. it is related to the order expansion
    ! and to the value of the polynomial coefficients
    IF(arg > 0.355_dp) THEN
       idij = 1.0_dp / dij
       e_arg_arg = EXP ( -arg ** 2 )
       tc = 1.0_dp / ( 1.0_dp + pc * arg )
       ! defining errf=1-erfc
       errf = 1.0_dp &
            - ((((ac5*tc+ac4)*tc+ac3)*tc+ac2)*tc+ac1) * tc * e_arg_arg
       efac = idij*errf
       ffac = idij**2 * (idij*errf - const*e_arg_arg)
    ELSE
       tc  = arg*arg
       tc2 = tc*tc
       efac = 1.0_dp - tc/3.0_dp + tc2/10.0_dp - tc*tc2/42.0_dp + tc2*tc2/216.0_dp
       efac = efac * const
       ffac = - 2.0_dp/3.0_dp + 2.0_dp/5.0_dp*tc - 1.0_dp/7.0_dp*tc2 + 1.0_dp/27.0_dp*tc2*tc
       ffac = const*alpha**2*ffac
    ENDIF
    ! getting the potential
    v_bonded_corr = v_bonded_corr - q1*q2*efac
    
    ! subtracting the force from the total force
    fscalar = q1 * q2 * ffac
    
    core_particle_set(i)%f(1)=core_particle_set(i)%f(1)-fscalar*rij(1)
    core_particle_set(i)%f(2)=core_particle_set(i)%f(2)-fscalar*rij(2)
    core_particle_set(i)%f(3)=core_particle_set(i)%f(3)-fscalar*rij(3)
    
    shell_particle_set(i)%f(1)=shell_particle_set(i)%f(1)+fscalar*rij(1)
    shell_particle_set(i)%f(2)=shell_particle_set(i)%f(2)+fscalar*rij(2)
    shell_particle_set(i)%f(3)=shell_particle_set(i)%f(3)+fscalar*rij(3)
    
    IF((shell_adiabatic).AND.use_virial) THEN
       fij_com = fscalar*rij
       pv_bc(1,1) = pv_bc(1,1) - fij_com(1) * rij(1)
       pv_bc(1,2) = pv_bc(1,2) - fij_com(1) * rij(2)
       pv_bc(1,3) = pv_bc(1,3) - fij_com(1) * rij(3)
       pv_bc(2,1) = pv_bc(2,1) - fij_com(2) * rij(1)
       pv_bc(2,2) = pv_bc(2,2) - fij_com(2) * rij(2)
       pv_bc(2,3) = pv_bc(2,3) - fij_com(2) * rij(3)
       pv_bc(3,1) = pv_bc(3,1) - fij_com(3) * rij(1)
       pv_bc(3,2) = pv_bc(3,2) - fij_com(3) * rij(2)
       pv_bc(3,3) = pv_bc(3,3) - fij_com(3) * rij(3)
    END IF

  END SUBROUTINE bonded_correct_gaussian_low_sh

! *****************************************************************************
!> \brief gives back information on atoms.. specific for bonded correct gaussian
!> \par History
!>      Splitted routines to clean and to fix a bug with the tensor whose 
!>      original definition was not correct for PBC..
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE get_atom_info(i, particle_set, shell_particle_set, core_particle_set,&
       shell_i_index, i_is_shell, rsi, rci, qci, qsi, zi)

    INTEGER, INTENT(IN)                      :: i
    TYPE(particle_type), POINTER             :: particle_set(:)
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set(:), &
                                                core_particle_set(:)
    INTEGER, INTENT(OUT)                     :: shell_i_index
    LOGICAL, INTENT(OUT)                     :: i_is_shell
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: rsi, rci
    REAL(KIND=dp), INTENT(OUT)               :: qci, qsi
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: zi

    CHARACTER(len=*), PARAMETER :: routineN = 'get_atom_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: zeta
    REAL(KIND=dp)                            :: qi
    REAL(KIND=dp), DIMENSION(3)              :: ri
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(shell_kind_type), POINTER           :: shell_i

    NULLIFY(atomic_kind, shell_i)
    atomic_kind => particle_set ( i ) % atomic_kind
    ri (:) = particle_set ( i ) % r (:)
    CALL get_atomic_kind ( atomic_kind,qeff = qi,&
                           shell_active=i_is_shell,&
                           shell=shell_i, z = zeta )
    IF(i_is_shell) THEN
       shell_i_index = particle_set(i)%shell_index
       rsi = shell_particle_set(shell_i_index)%r
       rci = core_particle_set(shell_i_index)%r
       qci = shell_i%charge_core
       qsi = shell_i%charge_shell
       IF(PRESENT(zi)) zi  = REAL(zeta)
    ELSE
       qsi = 0.0_dp
       qci = qi
       rci = ri
       IF(PRESENT(zi)) zi  = 0.0_dp
    END IF

  END SUBROUTINE get_atom_info

END MODULE fist_nonbond_force
