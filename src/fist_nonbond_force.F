!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_nonbond_force [1.0] *
!!
!!   NAME
!!     fist_nonbond_force
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (11 May 2001) : cleaning up of support structures
!!     CJM & HAF (27 July 2001): fixed bug with handling of cutoff larger than
!!                               half the boxsize.
!!     07.02.2005: getting rid of scaled_to_real calls in force loop (MK)
!!
!!   SOURCE
!******************************************************************************

MODULE fist_nonbond_force

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE mathlib,                         ONLY: matvec_3x3
  USE message_passing,                 ONLY: mp_sum
  USE pair_potential_types,            ONLY: ea_type,&
                                             nosh_nosh,&
                                             nosh_sh,&
                                             pair_potential_pp_type,&
                                             pair_potential_single_type,&
                                             sh_sh,&
                                             tersoff_type
  USE particle_types,                  ONLY: particle_type
  USE shell_potential_types,           ONLY: shell_kind_type
  USE splines,                         ONLY: potential_s,&
                                             spline_data_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_nonbond_force'

  PUBLIC :: force_nonbond,&
       bonded_correct_gaussian

  !!***
  !******************************************************************************
CONTAINS
  !******************************************************************************

  SUBROUTINE force_nonbond ( fist_nonbond_env, particle_set, cell,  &
       pot_nonbond, f_nonbond, ptens_nonbond, &
       fshell_nonbond, fcore_nonbond,  error )

    ! Calculates the force and the potential of the minimum image, and
    ! the pressure tensor
    !

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(OUT)               :: pot_nonbond
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_nonbond, ptens_nonbond
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: fshell_nonbond, fcore_nonbond
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_nonbond', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, handle, &
                                                ikind, ilist, ipair, jkind, &
                                                nkinds, npairs, shell_a, &
                                                shell_b, shell_type
    INTEGER, DIMENSION(:, :), POINTER        :: list
    LOGICAL                                  :: all_terms, do_shell, failure, &
                                                shell_present
    REAL(KIND=dp) :: cell_v(3), energy, fr(3), fr_com(3), fscalar, ptens11, &
      ptens12, ptens13, ptens21, ptens22, ptens23, ptens31, ptens32, ptens33, &
      rab(3), rab2, rab2_com, rab2_max, rab_cc(3), rab_com(3), rab_cs(3), &
      rab_sc(3), rab_ss(3)
    REAL(KIND=dp), DIMENSION(3)              :: cvi
    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_com
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update, &
                                                r_last_update_pbc, &
                                                rcore_last_update_pbc, &
                                                rshell_last_update_pbc
    TYPE(spline_data_type), POINTER          :: c_spline_data, &
                                                cs_spline_data, &
                                                sc_spline_data, spline_data

    CALL timeset ( routineN, handle )

    NULLIFY(c_spline_data,cs_spline_data,sc_spline_data,rshell_last_update_pbc)
    CALL fist_nonbond_env_get ( fist_nonbond_env, nonbonded = nonbonded, &
         potparm = potparm, r_last_update=r_last_update, &
         r_last_update_pbc=r_last_update_pbc,&
         rshell_last_update_pbc=rshell_last_update_pbc,&
         rcore_last_update_pbc=rcore_last_update_pbc,&
         natom_types = nkinds ,error=error)

    ! initializing the potential energy, pressure tensor and force
    pot_nonbond = 0.0_dp
    f_nonbond(:,:) = 0.0_dp
    ptens_nonbond(:,:) = 0.0_dp
    IF(PRESENT(fshell_nonbond)) fshell_nonbond = 0.0_dp
    IF(PRESENT(fcore_nonbond)) fcore_nonbond = 0.0_dp
    shell_present = .FALSE.
    IF(PRESENT(fshell_nonbond)) shell_present = .TRUE.
    ptens11 = 0.0_dp ; ptens12 = 0.0_dp ; ptens13 = 0.0_dp
    ptens21 = 0.0_dp ; ptens22 = 0.0_dp ; ptens23 = 0.0_dp
    ptens31 = 0.0_dp ; ptens32 = 0.0_dp ; ptens33 = 0.0_dp
    ! starting the force loop
    nkinds = SIZE ( potparm%pot, 1 )
! incorrect omp (i.e. threads can write to the variable an other thread might be reading (e.g. pot_nonbonded))
! !$OMP parallel do default(none) &
! !$OMP shared(nkinds) &
! !$OMP private(ikind,jkind) &
! !$OMP private(spline_data) &
! !$OMP shared(potparm) &
! !$OMP private(ab) &
! !$OMP shared(nonbonded) &
! !$OMP private(rab,rab2_max) &
! !$OMP private(ilist,nlist,atom_a,atom_b,npairs,cvi,list,inode,nnode) &
! !$OMP shared(cell) &
! !$OMP shared(pot_nonbond) &
! !$OMP private(f_eam,pot,fr,fscalar,energy) &
! !$OMP shared(f_nonbond,ptens_nonbond,r_last_update) &
! !$OMP private(rab2) &
! !$OMP private(neighbor_kind_pair,cell_v) &
! !$OMP shared(particle_set,r_last_update_pbc) &
! !$OMP reduction(+:ptens11,ptens12,ptens13,ptens21,ptens22,ptens23, &
! !$OMP             ptens31,ptens32,ptens33) &
! !$OMP private(ipair) 

    DO ilist=1,nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
       npairs=neighbor_kind_pair%npairs
       IF (npairs ==0) CYCLE
       ikind = neighbor_kind_pair%ikind
       jkind = neighbor_kind_pair%jkind
       list  =>neighbor_kind_pair%list
       cvi   = neighbor_kind_pair%cell_vector
       spline_data => potparm % pot( ikind, jkind ) %pot % pair_spline_data
       IF(potparm %pot( ikind, jkind ) %pot % TYPE == ea_type) CYCLE
       IF(potparm %pot( ikind, jkind ) %pot % TYPE == tersoff_type) CYCLE
       shell_type=potparm %pot( ikind, jkind ) %pot % shell_type 
       IF(shell_type/=nosh_nosh) THEN
          c_spline_data => potparm % pot( ikind, jkind ) %pot % c_coul_spline_data
          CPPostcondition(PRESENT(fshell_nonbond),cp_failure_level,routineP,error,failure)
          do_shell = .TRUE.
       ELSE
          c_spline_data => potparm % pot( ikind, jkind ) %pot % pair_spline_data
          do_shell = .FALSE.
       END IF
       IF(shell_type==sh_sh) THEN
          cs_spline_data => potparm % pot( ikind, jkind ) %pot % cs_coul_spline_data
          sc_spline_data => potparm % pot( ikind, jkind ) %pot % sc_coul_spline_data
       END IF
       rab2_max = potparm%pot(ikind,jkind)%pot%rcutsq
       CALL matvec_3x3(cell_v, cell%hmat, cvi)
       DO ipair = 1, npairs
          atom_a = list(1,ipair)
          atom_b = list(2,ipair)
          all_terms = .TRUE.
!          done_sh = 0
          IF (do_shell) THEN
             IF(shell_type==sh_sh) THEN
                shell_a = particle_set(atom_a)%shell_index
                shell_b = particle_set(atom_b)%shell_index
                rab_cc = rcore_last_update_pbc(shell_b)%r - rcore_last_update_pbc(shell_a)%r
                rab_cs = rshell_last_update_pbc(shell_b)%r - rcore_last_update_pbc(shell_a)%r
                rab_sc = rcore_last_update_pbc(shell_b)%r - rshell_last_update_pbc(shell_a)%r
                rab_ss = rshell_last_update_pbc(shell_b)%r - rshell_last_update_pbc(shell_a)%r
             ELSEIF(shell_type==nosh_sh .AND. particle_set(atom_a)%shell_index/=0) THEN
                shell_a = particle_set(atom_a)%shell_index
                shell_b = 0
                rab_cc = r_last_update_pbc(atom_b)%r - rcore_last_update_pbc(shell_a)%r
                rab_sc = 0.0_dp
                rab_cs = 0.0_dp
                rab_ss = r_last_update_pbc(atom_b)%r - rshell_last_update_pbc(shell_a)%r
             ELSEIF(shell_type==nosh_sh .AND. particle_set(atom_b)%shell_index/=0) THEN
                shell_b = particle_set(atom_b)%shell_index
                shell_a = 0
                rab_cc = rcore_last_update_pbc(shell_b)%r - r_last_update_pbc(atom_a)%r
                rab_sc = 0.0_dp
                rab_cs = 0.0_dp
                rab_ss = rshell_last_update_pbc(shell_b)%r - r_last_update_pbc(atom_a)%r
             END IF
          ELSE
             shell_a=0
             shell_b=0
             rab_cc = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
             rab_sc = 0.0_dp
             rab_cs = 0.0_dp
             rab_ss = 0.0_dp
          END IF

          fr_com(:) = 0.0_dp
          pv_com = 0.0_dp
          rab_com = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
          rab_com = rab_com +cell_v
          rab2_com =  rab_com(1)*rab_com(1) + rab_com(2)*rab_com(2) + rab_com(3)*rab_com(3)
           fr = 0.0_dp
           rab = rab_cc
           rab = rab+cell_v
           rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
           IF (rab2 <= rab2_max) THEN
              ! if shell-model core-core or core-ion/ion-core: Coulomb only
              ! else ion-ion : VDW + Coulomb
              energy =  potential_s (c_spline_data,rab2,fscalar)
              pot_nonbond = pot_nonbond + energy
              fr(1) = fscalar*rab(1)
              fr(2) = fscalar*rab(2)
              fr(3) = fscalar*rab(3)
           END IF
           
           IF(shell_a /= 0 ) THEN
              fcore_nonbond(1,shell_a) = fcore_nonbond(1,shell_a) - fr(1)
              fcore_nonbond(2,shell_a) = fcore_nonbond(2,shell_a) - fr(2)
              fcore_nonbond(3,shell_a) = fcore_nonbond(3,shell_a) - fr(3)
           ELSE
              f_nonbond(1,atom_a) = f_nonbond(1,atom_a) - fr(1)
              f_nonbond(2,atom_a) = f_nonbond(2,atom_a) - fr(2)
              f_nonbond(3,atom_a) = f_nonbond(3,atom_a) - fr(3)
           END IF

           IF(shell_b /= 0 ) THEN
              fcore_nonbond(1,shell_b) = fcore_nonbond(1,shell_b) + fr(1)
              fcore_nonbond(2,shell_b) = fcore_nonbond(2,shell_b) + fr(2)
              fcore_nonbond(3,shell_b) = fcore_nonbond(3,shell_b) + fr(3)
           ELSE
              f_nonbond(1,atom_b) = f_nonbond(1,atom_b) + fr(1)
              f_nonbond(2,atom_b) = f_nonbond(2,atom_b) + fr(2)
              f_nonbond(3,atom_b) = f_nonbond(3,atom_b) + fr(3)
           END IF

           fr_com(1) = fr_com(1) + fr(1) 
           fr_com(2) = fr_com(2) + fr(2) 
           fr_com(3) = fr_com(3) + fr(3) 

           pv_com(1,1) = pv_com(1,1) + rab(1) * fr(1) 
           pv_com(1,2) = pv_com(1,2) + rab(1) * fr(2) 
           pv_com(1,3) = pv_com(1,3) + rab(1) * fr(3) 
           pv_com(2,1) = pv_com(2,1) + rab(2) * fr(1) 
           pv_com(2,2) = pv_com(2,2) + rab(2) * fr(2) 
           pv_com(2,3) = pv_com(2,3) + rab(2) * fr(3) 
           pv_com(3,1) = pv_com(3,1) + rab(3) * fr(1) 
           pv_com(3,2) = pv_com(3,2) + rab(3) * fr(2) 
           pv_com(3,3) = pv_com(3,3) + rab(3) * fr(3) 

           IF (do_shell) THEN
              IF(shell_type==sh_sh) THEN
                 ! shell-shell : VDW + Coulomb
                 rab = rab_ss 
                 rab = rab+cell_v
                 rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                 IF (rab2 <= rab2_max) THEN
                    energy =  potential_s (spline_data,rab2,fscalar)
                    pot_nonbond = pot_nonbond + energy
                    fr(1) = fscalar*rab(1)
                    fr(2) = fscalar*rab(2)
                    fr(3) = fscalar*rab(3)
                    fshell_nonbond(1,shell_a) = fshell_nonbond(1,shell_a) - fr(1)
                    fshell_nonbond(2,shell_a) = fshell_nonbond(2,shell_a) - fr(2)
                    fshell_nonbond(3,shell_a) = fshell_nonbond(3,shell_a) - fr(3)
                    fshell_nonbond(1,shell_b) = fshell_nonbond(1,shell_b) + fr(1)
                    fshell_nonbond(2,shell_b) = fshell_nonbond(2,shell_b) + fr(2)
                    fshell_nonbond(3,shell_b) = fshell_nonbond(3,shell_b) + fr(3)
                    
                    fr_com(1) = fr_com(1) + fr(1) 
                    fr_com(2) = fr_com(2) + fr(2) 
                    fr_com(3) = fr_com(3) + fr(3) 
                    
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
                 rab = rab_cs 
                 rab = rab+cell_v
                 rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                 IF (rab2 <= rab2_max) THEN
                    energy =  potential_s (cs_spline_data,rab2,fscalar)
                    pot_nonbond = pot_nonbond + energy
                    fr(1) = fscalar*rab(1)
                    fr(2) = fscalar*rab(2)
                    fr(3) = fscalar*rab(3)
                    fcore_nonbond(1,shell_a) = fcore_nonbond(1,shell_a) - fr(1)
                    fcore_nonbond(2,shell_a) = fcore_nonbond(2,shell_a) - fr(2)
                    fcore_nonbond(3,shell_a) = fcore_nonbond(3,shell_a) - fr(3)
                    fshell_nonbond(1,shell_b) = fshell_nonbond(1,shell_b) + fr(1)
                    fshell_nonbond(2,shell_b) = fshell_nonbond(2,shell_b) + fr(2)
                    fshell_nonbond(3,shell_b) = fshell_nonbond(3,shell_b) + fr(3)
                    
                    fr_com(1) = fr_com(1) + fr(1) 
                    fr_com(2) = fr_com(2) + fr(2) 
                    fr_com(3) = fr_com(3) + fr(3) 
                    
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
                 rab = rab_sc
                 rab = rab+cell_v
                 rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                 IF (rab2 <= rab2_max) THEN
                    energy =  potential_s (sc_spline_data,rab2,fscalar)
                    pot_nonbond = pot_nonbond + energy
                    fr(1) = fscalar*rab(1)
                    fr(2) = fscalar*rab(2)
                    fr(3) = fscalar*rab(3)
                    fshell_nonbond(1,shell_a) = fshell_nonbond(1,shell_a) - fr(1)
                    fshell_nonbond(2,shell_a) = fshell_nonbond(2,shell_a) - fr(2)
                    fshell_nonbond(3,shell_a) = fshell_nonbond(3,shell_a) - fr(3)
                    fcore_nonbond(1,shell_b) = fcore_nonbond(1,shell_b) + fr(1)
                    fcore_nonbond(2,shell_b) = fcore_nonbond(2,shell_b) + fr(2)
                    fcore_nonbond(3,shell_b) = fcore_nonbond(3,shell_b) + fr(3)
                    
                    fr_com(1) = fr_com(1) + fr(1) 
                    fr_com(2) = fr_com(2) + fr(2) 
                    fr_com(3) = fr_com(3) + fr(3) 
                    
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
                 rab = rab_ss
                 rab = rab+cell_v
                 rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                 IF (rab2 <= rab2_max) THEN
                    energy =  potential_s (spline_data,rab2,fscalar)
                    pot_nonbond = pot_nonbond + energy
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
                    
                    fr_com(1) = fr_com(1) + fr(1) 
                    fr_com(2) = fr_com(2) + fr(2) 
                    fr_com(3) = fr_com(3) + fr(3) 
                    
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
           END IF
           ptens11 = ptens11 + pv_com(1,1)
           ptens21 = ptens21 + pv_com(2,1)
           ptens31 = ptens31 + pv_com(3,1)
           ptens12 = ptens12 + pv_com(1,2)
           ptens22 = ptens22 + pv_com(2,2)
           ptens32 = ptens32 + pv_com(3,2)
           ptens13 = ptens13 + pv_com(1,3)
           ptens23 = ptens23 + pv_com(2,3)
           ptens33 = ptens33 + pv_com(3,3)
        END DO
     END DO
     ptens_nonbond(1,1) = ptens_nonbond(1,1) + ptens11
     ptens_nonbond(1,2) = ptens_nonbond(1,2) + ptens12
     ptens_nonbond(1,3) = ptens_nonbond(1,3) + ptens13
     ptens_nonbond(2,1) = ptens_nonbond(2,1) + ptens21
     ptens_nonbond(2,2) = ptens_nonbond(2,2) + ptens22
     ptens_nonbond(2,3) = ptens_nonbond(2,3) + ptens23
     ptens_nonbond(3,1) = ptens_nonbond(3,1) + ptens31
     ptens_nonbond(3,2) = ptens_nonbond(3,2) + ptens32
     ptens_nonbond(3,3) = ptens_nonbond(3,3) + ptens33
     
     CALL timestop ( handle )
   END SUBROUTINE force_nonbond
   
  !******************************************************************************

  SUBROUTINE bonded_correct_gaussian ( atomic_kind_set, local_particles, &
             particle_set, ewald_env, v_bonded_corr, pv_bc, shell_particle_set, &
             core_particle_set, error )


    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    REAL(KIND=dp), INTENT(OUT)               :: v_bonded_corr
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: pv_bc
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set( : ), &
                                                core_particle_set( : )
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'bonded_correct_gaussian'
    REAL(KIND=dp), PARAMETER :: ac1 = 0.254829592_dp, ac2 = -0.284496736_dp, &
      ac3 = 1.421413741_dp, ac4 = -1.453152027_dp, ac5 = 1.061405429_dp, &
      pc = 0.3275911_dp

    INTEGER                                  :: group, handle, i, iatom, &
                                                iexl, ikind, j, &
                                                natoms_per_kind, nkind, &
                                                shell_i_index, shell_j_index
    LOGICAL                                  :: do_shell, i_is_shell, &
                                                j_is_shell, shell_adiabatic
    REAL(KIND=dp) :: alpha, arg, const, dij, e_arg_arg, efac, eps0, errf, &
      ffac, flops, fscalar, idij, ifourpieps0, qci, qcj, qi, qj, qsi, qsj, &
      rijsq, tc, tc2
    REAL(KIND=dp), DIMENSION(3)              :: fij_com, rci, rcj, ri, rij, &
                                                rj, rsi, rsj
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(shell_kind_type), POINTER           :: shell_i, shell_j

! parameters for the numerical erf in Abramowitz and Stegun Eq. 7.1.26 pg 299
!------------------------------------------------------------------------------

    CALL timeset ( routineN, 'I', 'Mflops', handle )
    flops = 0.0_dp

    NULLIFY(atomic_kind, shell_i, shell_j)
    ! defining the constants
    pv_bc = 0.0_dp
    v_bonded_corr = 0.0_dp

    CALL ewald_env_get ( ewald_env, eps0 = eps0, alpha = alpha, group = group ,error=error)
    const = 2.0_dp * alpha / SQRT( pi )
    ifourpieps0 = 1.0_dp / ( 4.0_dp * pi * eps0 )

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, shell_adiabatic=shell_adiabatic)
    nkind = SIZE ( atomic_kind_set )
    DO ikind = 1, nkind
       natoms_per_kind = local_particles % n_el ( ikind )
       DO iatom = 1, natoms_per_kind
          i = local_particles % list ( ikind ) % array ( iatom )
          atomic_kind => particle_set ( i ) % atomic_kind
          ri ( : ) = particle_set ( i ) % r ( : )
          CALL get_atomic_kind ( atomic_kind,qeff = qi,&
               shell_active=i_is_shell,&
               shell=shell_i )
          IF(i_is_shell) THEN
             shell_i_index = particle_set(i)%shell_index
             rsi = shell_particle_set(shell_i_index)%r
             rci = core_particle_set(shell_i_index)%r
             qci = shell_i%charge_core
             qsi = shell_i%charge_shell
             qi = qci
          ELSE
             qsi = 0.0_dp
             qci = qi
             rci = ri
          END IF
          DO iexl = 1, SIZE(particle_set(i)%list_exclude_ei)
             j = particle_set(i)%list_exclude_ei(iexl)
             IF ( j==i ) CYCLE
             fij_com = 0.0_dp
             atomic_kind => particle_set ( j ) % atomic_kind
             rj ( : ) = particle_set ( j ) % r ( : )
             CALL get_atomic_kind ( atomic_kind,qeff = qj,&
                  shell_active=j_is_shell, &
                  shell=shell_j )
             IF(j_is_shell) THEN
               shell_j_index = particle_set(j)%shell_index
               rsj = shell_particle_set(shell_j_index)%r
               rcj = core_particle_set(shell_j_index)%r
               qcj = shell_j%charge_core
               qsj = shell_j%charge_shell
               qj = qcj
             ELSE
               qsj = 0.0_dp
               qcj = qj
               rcj = rj
             END IF
             do_shell = .FALSE.
             IF (i_is_shell.OR.j_is_shell) do_shell = .TRUE.
             ! shell-model: exclude core-core, core-ion, or ion-core interactions
             ! if not shell-model: exclude ion-ion interactions
             rij(1) = rci(1) - rcj(1)
             rij(2) = rci(2) - rcj(2)
             rij(3) = rci(3) - rcj(3)
             rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
             idij = 1.0_dp / SQRT ( rijsq )
             dij = rijsq * idij
             arg = alpha * dij
             e_arg_arg = EXP ( -arg ** 2 )
             tc = 1.0_dp / ( 1.0_dp + pc * arg )

             ! defining errf=1-erfc
             errf = 1.0_dp &
                  - ((((ac5*tc+ac4)*tc+ac3)*tc+ac2)*tc+ac1) * tc * e_arg_arg

             ! getting the potential
             v_bonded_corr = v_bonded_corr - qci*qcj*idij*errf

             ! subtracting the force from the total force
             fscalar = ifourpieps0 * qci * qcj * idij ** 2 &
                  * ( idij * errf - const * e_arg_arg  )

             IF(i_is_shell) THEN
               core_particle_set ( shell_i_index ) % f(1) = &
                    core_particle_set ( shell_i_index ) % f(1) - fscalar*rij(1)
               core_particle_set ( shell_i_index ) % f(2) = &
                    core_particle_set ( shell_i_index ) % f(2) - fscalar*rij(2)
               core_particle_set ( shell_i_index ) % f(3) = &
                    core_particle_set ( shell_i_index ) % f(3) - fscalar*rij(3)
             ELSE
               particle_set ( i ) % f(1) = particle_set ( i ) % f(1) - fscalar*rij(1)
               particle_set ( i ) % f(2) = particle_set ( i ) % f(2) - fscalar*rij(2)
               particle_set ( i ) % f(3) = particle_set ( i ) % f(3) - fscalar*rij(3)
             END IF

             fij_com = fscalar*rij
             ! shell-model: exclude shell-shell interactions
             IF (do_shell) THEN
                IF(i_is_shell .AND. j_is_shell) THEN
                   rij(1) = rsi(1) - rsj(1)
                   rij(2) = rsi(2) - rsj(2)
                   rij(3) = rsi(3) - rsj(3)
                   rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
                   idij = 1.0_dp / SQRT ( rijsq )
                   dij = rijsq * idij
                   arg = alpha * dij
                   e_arg_arg = EXP ( -arg ** 2 )
                   tc = 1.0_dp / ( 1.0_dp + pc * arg )
                   
                   ! defining errf=1-erfc
                   errf = 1.0_dp &
                        - ((((ac5*tc+ac4)*tc+ac3)*tc+ac2)*tc+ac1) * tc * e_arg_arg
                   
                   ! getting the potential
                   v_bonded_corr = v_bonded_corr - qsi*qsj*idij*errf
                   
                   ! subtracting the force from the total force
                   fscalar = ifourpieps0 * qsi * qsj * idij ** 2 &
                        * ( idij * errf - const * e_arg_arg  )
                   
                   shell_particle_set ( shell_i_index) % f(1) = &
                        shell_particle_set ( shell_i_index ) % f(1) - fscalar*rij(1)
                   shell_particle_set ( shell_i_index ) % f(2) = &
                        shell_particle_set ( shell_i_index ) % f(2) - fscalar*rij(2)
                   shell_particle_set ( shell_i_index ) % f(3) = &
                        shell_particle_set ( shell_i_index ) % f(3) - fscalar*rij(3)
                   
                   IF(shell_adiabatic) THEN
                      fij_com = fij_com + fscalar*rij
                   END IF
                   
                END IF
                ! shell-model: exclude shell_i-core interactions
                IF(i_is_shell ) THEN
                   rij(1) = rsi(1) - rcj(1)
                   rij(2) = rsi(2) - rcj(2)
                   rij(3) = rsi(3) - rcj(3)
                   rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
                   idij = 1.0_dp / SQRT ( rijsq )
                   dij = rijsq * idij
                   arg = alpha * dij
                   e_arg_arg = EXP ( -arg ** 2 )
                   tc = 1.0_dp / ( 1.0_dp + pc * arg )
                   
                   ! defining errf=1-erfc
                   errf = 1.0_dp &
                        - ((((ac5*tc+ac4)*tc+ac3)*tc+ac2)*tc+ac1) * tc * e_arg_arg
                   
                   ! getting the potential
                   v_bonded_corr = v_bonded_corr - qsi*qcj*idij*errf
                   
                   ! subtracting the force from the total force
                   fscalar = ifourpieps0 * qsi * qcj * idij ** 2 &
                        * ( idij * errf - const * e_arg_arg  )
                   
                   shell_particle_set ( shell_i_index) % f(1) = &
                        shell_particle_set ( shell_i_index ) % f(1) - fscalar*rij(1)
                   shell_particle_set ( shell_i_index ) % f(2) = &
                        shell_particle_set ( shell_i_index ) % f(2) - fscalar*rij(2)
                   shell_particle_set ( shell_i_index ) % f(3) = &
                        shell_particle_set ( shell_i_index ) % f(3) - fscalar*rij(3)
                   
                   IF(shell_adiabatic) THEN
                      fij_com = fij_com + fscalar*rij
                   END IF
                   
                END IF
                IF(j_is_shell ) THEN
                   rij(1) = rci(1) - rsj(1)
                   rij(2) = rci(2) - rsj(2)
                   rij(3) = rci(3) - rsj(3)
                   rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
                   idij = 1.0_dp / SQRT ( rijsq )
                   dij = rijsq * idij
                   arg = alpha * dij
                   e_arg_arg = EXP ( -arg ** 2 )
                   tc = 1.0_dp / ( 1.0_dp + pc * arg )
                   
                   ! defining errf=1-erfc
                   errf = 1.0_dp &
                        - ((((ac5*tc+ac4)*tc+ac3)*tc+ac2)*tc+ac1) * tc * e_arg_arg
                   
                   ! getting the potential
                   v_bonded_corr = v_bonded_corr - qci*qsj*idij*errf
                   
                   ! subtracting the force from the total force
                   fscalar = ifourpieps0 * qci * qsj * idij ** 2 &
                        * ( idij * errf - const * e_arg_arg  )
                   
                   IF(i_is_shell) THEN
                      core_particle_set ( shell_i_index) % f(1) = &
                           core_particle_set ( shell_i_index ) % f(1) - fscalar*rij(1)
                      core_particle_set ( i) % f(2) = &
                           core_particle_set ( shell_i_index ) % f(2) - fscalar*rij(2)
                      core_particle_set (i) % f(3) = &
                           core_particle_set ( shell_i_index ) % f(3) - fscalar*rij(3)
                   ELSE
                      particle_set ( i) % f(1) = &
                           particle_set ( i ) % f(1) - fscalar*rij(1)
                      particle_set ( i) % f(2) = &
                           particle_set ( i ) % f(2) - fscalar*rij(2)
                      particle_set (i) % f(3) = &
                           particle_set ( i ) % f(3) - fscalar*rij(3)
                   END IF
                   
                   IF(shell_adiabatic) THEN
                      fij_com = fij_com + fscalar*rij
                   END IF
                   
                END IF
             END IF
             pv_bc(1,1) = pv_bc(1,1) - fij_com(1) * particle_set ( i ) % r(1)
             pv_bc(1,2) = pv_bc(1,2) - fij_com(1) * particle_set ( i ) % r(2)
             pv_bc(1,3) = pv_bc(1,3) - fij_com(1) * particle_set ( i ) % r(3)
             pv_bc(2,1) = pv_bc(2,1) - fij_com(2) * particle_set ( i ) % r(1)
             pv_bc(2,2) = pv_bc(2,2) - fij_com(2) * particle_set ( i ) % r(2)
             pv_bc(2,3) = pv_bc(2,3) - fij_com(2) * particle_set ( i ) % r(3)
             pv_bc(3,1) = pv_bc(3,1) - fij_com(3) * particle_set ( i ) % r(1)
             pv_bc(3,2) = pv_bc(3,2) - fij_com(3) * particle_set ( i ) % r(2)
             pv_bc(3,3) = pv_bc(3,3) - fij_com(3) * particle_set ( i ) % r(3)
             flops = flops + 62.0_dp
          ENDDO

          ! Exclusion of shell-core interaction
          IF(i_is_shell ) THEN
             rij(1) = rci(1) - rsi(1)
             rij(2) = rci(2) - rsi(2)
             rij(3) = rci(3) - rsi(3)
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
             ! th factor 2 comes to compensate for the factor 0.5 
             ! introduced later to avoid double counting
             v_bonded_corr = v_bonded_corr - 2.0_dp * qci*qsi*efac

             ! subtracting the force from the total force
             fscalar = ifourpieps0 * qci * qsi * ffac

             core_particle_set ( shell_i_index) % f(1) = &
                          core_particle_set ( shell_i_index ) % f(1) - fscalar*rij(1)
             core_particle_set ( shell_i_index) % f(2) = &
                          core_particle_set ( shell_i_index ) % f(2) - fscalar*rij(2)
             core_particle_set (shell_i_index) % f(3) = &
                          core_particle_set ( shell_i_index ) % f(3) - fscalar*rij(3)

             shell_particle_set (shell_i_index) % f(1) = &
                          shell_particle_set ( shell_i_index ) % f(1) + fscalar*rij(1)
             shell_particle_set ( shell_i_index) % f(2) = &
                          shell_particle_set ( shell_i_index ) % f(2) + fscalar*rij(2)
             shell_particle_set (shell_i_index) % f(3) = &
                          shell_particle_set ( shell_i_index ) % f(3) + fscalar*rij(3)
          END IF

       END DO
    ENDDO
    ! the factor of 1/2 comes from double counting in the exclusion list
    v_bonded_corr = v_bonded_corr * ifourpieps0 * 0.5_dp
    CALL mp_sum ( v_bonded_corr, group )

    flops = flops * 1.0E-6_dp
    CALL timestop ( flops, handle )

  END SUBROUTINE bonded_correct_gaussian

  !******************************************************************************

END MODULE fist_nonbond_force
