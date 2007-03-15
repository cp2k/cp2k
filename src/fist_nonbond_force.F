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
  USE cell_types,                      ONLY: cell_type,pbc
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mathconstants,                   ONLY: pi
  USE mathlib,                         ONLY: matvec_3x3
  USE message_passing,                 ONLY: mp_sum
  USE pair_potential_types,            ONLY: ea_type,&
                                             tersoff_type,&
                                             pair_potential_pp_type,&
                                             pair_potential_single_type,&
                                             nosh_sh,&
                                             nosh_nosh,&
                                             sh_sh
  USE particle_types,                  ONLY: particle_type
  USE shell_potential_types,           ONLY: shell_kind_type
  USE splines,                         ONLY: potential_s,&
                                             spline_data_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_nonbond_force'

  PUBLIC :: force_nonbond,&
       bonded_correct_gaussian,&
       density_nonbond,&
       energy_manybody

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
      INTENT(OUT)                            :: f_nonbond, ptens_nonbond
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: fshell_nonbond
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: fcore_nonbond
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_nonbond', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ab, atom_a, atom_b, handle, &
                                                ikind, ilist, inode, ipair, &
                                                jkind, nkinds, nlist, nnode, &
                                                npairs, shell_a, shell_b, &
                                                shell_type
    INTEGER, DIMENSION(:, :), POINTER        :: list
    LOGICAL                                  :: failure, is_eam, is_shell, shell_present, is_tersoff
    REAL(KIND=dp) :: cell_v(3), energy, f_eam, fr(3), fr_com(3), fscalar, ptens11, &
      ptens12, ptens13, ptens21, ptens22, ptens23, ptens31, ptens32, ptens33, &
      rab(3), rab2, rab2_com, rab2_max, rab_com(3), rab_cc(3), rab_cs(3), rab_sc(3), rab_ss(3)
    REAL(KIND=dp), DIMENSION(3)              :: cvi
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update, &
                                                r_last_update_pbc, &
                                                rcore_last_update_pbc,&
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
! !$OMP private(is_eam,ab) &
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
       is_eam=(potparm %pot( ikind, jkind ) %pot % TYPE == ea_type)
       is_tersoff=(potparm %pot( ikind, jkind ) %pot % TYPE == tersoff_type)
       shell_type=potparm %pot( ikind, jkind ) %pot % shell_type 
       IF(shell_type/=nosh_nosh) THEN
         c_spline_data => potparm % pot( ikind, jkind ) %pot % c_coul_spline_data
         CPPostcondition(PRESENT(fshell_nonbond),cp_failure_level,routineP,error,failure)
       ELSE
         c_spline_data => potparm % pot( ikind, jkind ) %pot % pair_spline_data
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
          ! which atom has a shell, if any
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
          ELSE
            shell_a=0
            shell_b=0
            rab_cc = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
            rab_sc = 0.0_dp
            rab_cs = 0.0_dp
            rab_ss = 0.0_dp
          END IF

          fr_com(:) = 0.0_dp
          rab_com = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
          rab_com = rab_com +cell_v
          rab2_com =  rab_com(1)*rab_com(1) + rab_com(2)*rab_com(2) + rab_com(3)*rab_com(3)

          IF (rab2_com <= rab2_max) THEN
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


             IF (is_eam) THEN
                pot => potparm %pot( ikind, jkind )%pot
                CALL get_force_eam ( rab2, pot, particle_set, atom_a, atom_b, &
                     f_eam )
                fr(1) =  fr(1) - f_eam * rab ( 1 )
                fr(2) =  fr(2) - f_eam * rab ( 2 )
                fr(3) =  fr(3) - f_eam * rab ( 3 )
             ENDIF
             IF (is_tersoff) THEN
                pot => potparm %pot( ikind, jkind )%pot
                CALL get_force_tersoff ( pot, particle_set, cell, &
                                         atom_a, atom_b, f_nonbond, ptens_nonbond )
             ENDIF

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
              END IF

             END IF
              ptens11 = ptens11 + rab_com(1)*fr_com(1)
              ptens21 = ptens21 + rab_com(2)*fr_com(1)
              ptens31 = ptens31 + rab_com(3)*fr_com(1)
              ptens12 = ptens12 + rab_com(1)*fr_com(2)
              ptens22 = ptens22 + rab_com(2)*fr_com(2)
              ptens32 = ptens32 + rab_com(3)*fr_com(2)
              ptens13 = ptens13 + rab_com(1)*fr_com(3)
              ptens23 = ptens23 + rab_com(2)*fr_com(3)
              ptens33 = ptens33 + rab_com(3)*fr_com(3)
          END IF
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

  SUBROUTINE density_nonbond ( fist_nonbond_env, particle_set, cell, para_env, error )
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'density_nonbond'

    INTEGER                                  :: atom_a, atom_b, handle, &
                                                ikind, ilist, ipair, &
                                                iparticle, istat, jkind, &
                                                nkinds, nparticle
    LOGICAL                                  :: do_eam
    REAL(KIND=dp)                            :: flops, rab2, rab2_max
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rho
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update, &
                                                r_last_update_pbc

    CALL timeset ( routineN, 'I', 'Mflops', handle )
    flops = 0.0_dp
    do_eam = .FALSE.
    CALL fist_nonbond_env_get ( fist_nonbond_env, nonbonded = nonbonded, &
         potparm = potparm, r_last_update=r_last_update,&
         r_last_update_pbc=r_last_update_pbc,error=error)
    nkinds = SIZE ( potparm%pot, 1 )
    DO ikind = 1, nkinds
       DO jkind = ikind, nkinds
          IF (potparm%pot( ikind, jkind )%pot % TYPE== ea_type) do_eam = .TRUE.
       END DO
    END DO
    ! Only if EAM potential are present
    IF (do_eam) THEN
       nparticle = SIZE ( particle_set )
       ALLOCATE ( rho ( nparticle ), STAT = istat )
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "rho",nparticle*dp_size)
       rho ( : ) = 0._dp
       ! Starting the force loop
       DO ilist=1,nonbonded%nlists
          neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
          IF (neighbor_kind_pair%npairs ==0) CYCLE
          ikind = neighbor_kind_pair%ikind
          jkind = neighbor_kind_pair%jkind
          SELECT CASE ( potparm%pot( ikind, jkind )%pot % TYPE )
          CASE ( ea_type)
             rab2_max = potparm%pot(ikind,jkind)%pot%rcutsq
             CALL matvec_3x3(cell_v, cell%hmat,REAL(neighbor_kind_pair%cell_vector,KIND=dp))
             DO ipair = 1, neighbor_kind_pair%npairs
                atom_a = neighbor_kind_pair%list(1,ipair)
                atom_b = neighbor_kind_pair%list(2,ipair)
                rab = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
                rab = rab+cell_v
                rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                IF (rab2 <= rab2_max) THEN
                   pot => potparm %pot( ikind, jkind )%pot
                   CALL get_rho_eam ( pot, rab2, atom_a, atom_b, rho  )
                END IF
             END DO
          END SELECT
       END DO
       CALL mp_sum( rho, para_env%group )
       DO iparticle=1,nparticle
          particle_set(iparticle)%rho = rho(iparticle)
       END DO

       DEALLOCATE ( rho, STAT = istat )
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "deallocate rho")
    END IF
    CALL timestop ( flops, handle )

  END SUBROUTINE density_nonbond
  !******************************************************************************
  SUBROUTINE energy_manybody ( fist_nonbond_env, atomic_kind_set, &
       local_particles, particle_set, cell, pot_manybody, para_env, error )

    ! computes the embedding contribution to the energy

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), INTENT(INOUT)                  :: pot_manybody
    TYPE(cp_para_env_type), POINTER,OPTIONAL :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'energy_manybody'

    INTEGER                                  :: ikind, index, index1,index2, &
                                                iparticle, iparticle_local, &
                                                jparticle, kparticle, &
                                                istat, nparticle, &
                                                nparticle_local
    REAL(dp)                                 :: embed, pot_loc, qr
    REAL(dp), POINTER                        :: fembed( : )
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot



    REAL(KIND=dp)                            :: rij(3),drij
    REAL(dp)                                 :: a_ij,b_ij,b_ji
    REAL(dp)                                 :: f_C, f_R, f_A
    REAL(dp)                                 :: RpD

    CALL fist_nonbond_env_get ( fist_nonbond_env, potparm = potparm ,error=error)
    DO ikind = 1, SIZE ( atomic_kind_set )
       pot => potparm %pot ( ikind, ikind ) % pot
       SELECT CASE ( pot % type )
       CASE DEFAULT
          CYCLE
       CASE ( ea_type )

          nparticle = SIZE ( particle_set )
          ALLOCATE ( fembed ( nparticle ), STAT = istat )
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
               "fembed",nparticle*dp_size)
          fembed ( : ) = 0._dp
          ! computation of embedding function and energy
          nparticle_local = local_particles%n_el(ikind)
          DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(ikind)%array(iparticle_local)
             index = INT ( particle_set ( iparticle ) % rho / pot % eam % drhoar) + 1
             IF ( index > 4999 ) index = 5000
             qr = particle_set ( iparticle ) % rho  -  pot % eam % rhoval ( index )

             embed = pot % eam % frho ( index ) + &
                  qr * pot % eam % frhop ( index )
             fembed ( iparticle ) = pot % eam % frhop ( index ) + &
                  qr * ( pot % eam % frhop ( index + 1 ) -    &
                  pot % eam % frhop ( index ) ) / pot % eam %  drhoar

             pot_manybody = pot_manybody + embed
          END DO
          ! communicate data
          CALL mp_sum( fembed, para_env%group )
          DO iparticle=1,nparticle
             particle_set(iparticle)%f_embed = fembed(iparticle)
          END DO

          DEALLOCATE ( fembed, STAT = istat )
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
               "deallocate fembed")

       CASE ( tersoff_type )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          RpD = pot%tersoff%bigR + pot%tersoff%bigD
          nparticle = SIZE ( particle_set )
          nparticle_local = local_particles%n_el(ikind)
          DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(ikind)%array(iparticle_local)
             DO jparticle=iparticle+1,nparticle
                pot_loc = 0.0_dp
                rij(:) = pbc(particle_set(jparticle)%r(:),particle_set(iparticle)%r(:),cell)
                drij = SQRT(DOT_PRODUCT(rij,rij))
                IF(drij>RpD) CYCLE
                !a_ij   = ter_a_ij(pot,drij)
                a_ij   = 1.0_dp
                b_ij   = 0.5_dp * ter_b_ij(pot,particle_set,cell,iparticle,jparticle)
                b_ji   = 0.5_dp * ter_b_ij(pot,particle_set,cell,jparticle,iparticle)
                f_C    = ter_f_C(pot,drij)
                f_A    = ter_f_A(pot,drij)
                f_R    = ter_f_R(pot,drij)
                pot_loc = f_C*(a_ij*f_R + (b_ij+b_ji)*f_A)
                pot_manybody = pot_manybody + pot_loc
             END DO
          END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END SELECT
    END DO
  END SUBROUTINE energy_manybody


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_f_C(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    REAL(KIND=dp)                                       :: ter_f_C
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: bigR,bigD
    REAL(KIND=dp)                                       :: RpD,RmD

    bigR = pot%tersoff%bigR
    bigD = pot%tersoff%bigD
    RmD = pot%tersoff%bigR-pot%tersoff%bigD
    RpD = pot%tersoff%bigR+pot%tersoff%bigD
    ter_f_C = 0.0_dp
    IF(r < RmD) ter_f_C = 1.0_dp
    IF(r > RpD) ter_f_C = 0.0_dp
    IF( (r < RpD) .AND. (r > RmD) ) THEN
      ter_f_C = 0.5_dp *(1.0_dp - SIN(0.5_dp*PI*(r-bigR)/(bigD)))
    END IF
  END FUNCTION ter_f_C

  FUNCTION ter_f_C_d(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    REAL(KIND=dp)                                       :: ter_f_C_d
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: bigR,bigD
    REAL(KIND=dp)                                       :: RpD,RmD

    bigR = pot%tersoff%bigR
    bigD = pot%tersoff%bigD
    RmD = pot%tersoff%bigR-pot%tersoff%bigD
    RpD = pot%tersoff%bigR+pot%tersoff%bigD
    ter_f_C_d = 0.0_dp
    IF(r < RmD) ter_f_C_d = 0.0_dp
    IF(r > RpD) ter_f_C_d = 0.0_dp
    IF( (r < RpD) .AND. (r > RmD) ) THEN
      ter_f_C_d = (0.25_dp*PI/bigD) * COS(0.5_dp*PI*(r-bigR)/(bigD)) / r
    END IF
  END FUNCTION ter_f_C_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_f_R(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    REAL(KIND=dp)                                       :: ter_f_R
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h

    A = pot%tersoff%A
    lambda1 = pot%tersoff%lambda1
    ter_f_R = 0.0_dp
    ter_f_R = A * EXP(-lambda1*r)
  END FUNCTION ter_f_R

  FUNCTION ter_f_R_d(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    REAL(KIND=dp)                                       :: ter_f_R_d
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: f_R

    A = pot%tersoff%A
    lambda1 = pot%tersoff%lambda1
    f_R = A * EXP(-lambda1*r)
    ter_f_R_d = 0.0_dp
    ter_f_R_d = lambda1 * f_R / r
  END FUNCTION ter_f_R_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_f_A(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    REAL(KIND=dp)                                       :: ter_f_A
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h

    B = pot%tersoff%B
    lambda2 = pot%tersoff%lambda2
    ter_f_A = 0.0_dp
    ter_f_A = - B * EXP(-lambda2*r)
  END FUNCTION ter_f_A

  FUNCTION ter_f_A_d(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    REAL(KIND=dp)                                       :: ter_f_A_d
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h

    B = pot%tersoff%B
    lambda2 = pot%tersoff%lambda2
    ter_f_A_d = 0.0_dp
    ter_f_A_d = - B *lambda2 * EXP(-lambda2*r) / r
  END FUNCTION ter_f_A_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_a_ij(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    REAL(KIND=dp)                                       :: ter_a_ij
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h

    n = pot%tersoff%n
    alpha = pot%tersoff%alpha
    ter_a_ij = 0.0_dp
    !Note alpha = 0.0_dp for the parameters in the paper so using simplified term
    !ter_a_ij = (1.0_dp+(alpha*ter_n_ij(pot,iparticle,jparticle,r))**n)**(-0.5_dp/n)
    ter_a_ij = 1.0_dp
  END FUNCTION ter_a_ij
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_b_ij(pot,particle_set,cell,iparticle,jparticle)
    TYPE(pair_potential_single_type), POINTER           :: pot
    TYPE(particle_type), INTENT(IN)                     :: particle_set( : )
    TYPE(cell_type), POINTER                            :: cell
    INTEGER, INTENT(IN)                                 :: iparticle, jparticle
    REAL(KIND=dp)                                       :: ter_b_ij
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: zeta_ij

    n = pot%tersoff%n
    beta = pot%tersoff%beta
    ter_b_ij = 0.0_dp
    zeta_ij = ter_zeta_ij(pot,particle_set,cell,iparticle,jparticle)
    ter_b_ij = (1.0_dp+(beta*zeta_ij)**n)**(-0.5_dp/n)
  END FUNCTION ter_b_ij

  FUNCTION ter_b_ij_d(pot,particle_set,cell,iparticle,jparticle)
    TYPE(pair_potential_single_type), POINTER           :: pot
    TYPE(particle_type), INTENT(IN)                     :: particle_set( : )
    TYPE(cell_type), POINTER                            :: cell
    INTEGER, INTENT(IN)                                 :: iparticle, jparticle
    REAL(KIND=dp)                                       :: ter_b_ij_d
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: beta_n,zeta_ij
    REAL(KIND=dp)                                       :: zeta_ij_nm1,zeta_ij_n

    n           = pot%tersoff%n
    beta        = pot%tersoff%beta
    beta_n      = beta**n
    zeta_ij     = ter_zeta_ij(pot,particle_set,cell,iparticle,jparticle)
    zeta_ij_nm1 = zeta_ij**(n-1.0_dp)
    zeta_ij_n   = zeta_ij**(n)

    ter_b_ij_d = 0.0_dp
    ter_b_ij_d = -0.5_dp * beta_n * zeta_ij_nm1 * &
                   ((1.0_dp + beta_n * zeta_ij_n)**((-0.5_dp/n) - 1.0_dp))
  END FUNCTION ter_b_ij_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_zeta_ij(pot,particle_set,cell,iparticle,jparticle)
    TYPE(pair_potential_single_type), POINTER           :: pot
    TYPE(particle_type), INTENT(IN)                     :: particle_set( : )
    TYPE(cell_type), POINTER                            :: cell
    INTEGER, INTENT(IN)                                 :: iparticle, jparticle
    REAL(KIND=dp)                                       :: ter_zeta_ij
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: RpD
    INTEGER                                             :: kparticle, nparticle
    REAL(KIND=dp)                                       :: rij(3),rik(3),drij,drik
    REAL(KIND=dp)                                       :: costheta,sintheta
    REAL(KIND=dp)                                       :: f_C,gterm,expterm

    ter_zeta_ij = 0.0_dp

    n = pot%tersoff%n
    lambda3 = pot%tersoff%lambda3
    RpD = pot%tersoff%bigR + pot%tersoff%bigD
    rij(:) = pbc(particle_set(jparticle)%r(:),particle_set(iparticle)%r(:),cell)
    drij = SQRT(DOT_PRODUCT(rij,rij))
    ter_zeta_ij = 0.0_dp
    nparticle = SIZE ( particle_set )
    DO kparticle=1,nparticle
      IF((iparticle==kparticle).OR.(jparticle==kparticle)) CYCLE
      rik(:) = pbc(particle_set(kparticle)%r(:),particle_set(iparticle)%r(:),cell)
      drik = SQRT(DOT_PRODUCT(rik,rik))
      IF(drik>RpD) CYCLE
      costheta = DOT_PRODUCT(rij,rik) / (drij*drik)
      IF(costheta<-1.0_dp) costheta=-1.0_dp
      IF(costheta>+1.0_dp) costheta=+1.0_dp
      f_C = ter_f_C(pot,drik)
      gterm = ter_g(pot,costheta)
      expterm = exp((lambda3*(drij-drik))**3)
      ter_zeta_ij = ter_zeta_ij + f_C*gterm*expterm
    END DO
  END FUNCTION ter_zeta_ij

  SUBROUTINE ter_zeta_ij_d(pot,particle_set,cell,iparticle,jparticle,f_nonbond,ptens_nonbond,prefactor)
    TYPE(pair_potential_single_type), POINTER           :: pot
    TYPE(particle_type), INTENT(IN)                     :: particle_set( : )
    TYPE(cell_type), POINTER                            :: cell
    INTEGER, INTENT(IN)                                 :: iparticle, jparticle
    REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)       :: f_nonbond,ptens_nonbond
    REAL(KIND=dp), INTENT(IN)                           :: prefactor
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: RpD
    INTEGER                                             :: kparticle, nparticle
    REAL(KIND=dp), DIMENSION(3)                         :: rij,rik,rij_hat,rik_hat
    REAL(KIND=dp)                                       :: costheta
    REAL(KIND=dp)                                       :: drij,drik
    REAL(KIND=dp), DIMENSION(3)                         :: dri,drj,drk
    REAL(KIND=dp), DIMENSION(3)                         :: dcosdri,dcosdrj,dcosdrk
    REAL(KIND=dp)                                       :: f_C,f_C_d
    REAL(KIND=dp)                                       :: gterm,gterm_d
    REAL(KIND=dp)                                       :: expterm,expterm_d
    REAL(KIND=dp)                                       :: crapi,crapj,crapk

    n = pot%tersoff%n
    lambda3 = pot%tersoff%lambda3
    RpD = pot%tersoff%bigR + pot%tersoff%bigD

    rij(:) = pbc(particle_set(jparticle)%r(:),particle_set(iparticle)%r(:),cell)
    drij = SQRT(DOT_PRODUCT(rij,rij))
    rij_hat(:) = rij(:) / drij

    nparticle = SIZE ( particle_set )
    DO kparticle=1,nparticle
      IF((iparticle==kparticle).OR.(jparticle==kparticle)) CYCLE
      rik(:) = pbc(particle_set(kparticle)%r(:),particle_set(iparticle)%r(:),cell)
      drik = SQRT(DOT_PRODUCT(rik,rik))
      rik_hat(:) = rik(:) / drik

      IF(drik>RpD) CYCLE

      costheta = DOT_PRODUCT(rij,rik) / (drij*drik)
      IF(costheta<-1.0_dp) costheta=-1.0_dp
      IF(costheta>+1.0_dp) costheta=+1.0_dp

      CALL ter_costheta_d(particle_set,cell,iparticle,jparticle,kparticle,&
                          dcosdri,dcosdrj,dcosdrk)

      f_C       = ter_f_C(pot,drik)
      f_C_d     = ter_f_C_d(pot,drik)
      gterm     = ter_g(pot,costheta)
      gterm_d   = ter_g_d(pot,costheta) !still need d(costheta)/dR term
      expterm   = exp((lambda3*(drij-drik))**3)
      expterm_d = (3.0_dp)*(lambda3**3)*((drij-drik)**2)*expterm

      dri =  f_C_d * gterm   * expterm   * (rik) &
           + f_C   * gterm_d * expterm   * (dcosdri) &
           + f_C   * gterm   * expterm_d * (-rij_hat+rik_hat)

      !No f_C_d component for Rj
      drj =  f_C   * gterm_d * expterm   * (dcosdrj) &
           + f_C   * gterm   * expterm_d * (rij_hat)

      drk =  f_C_d * gterm   * expterm   * (-rik) &
           + f_C   * gterm_d * expterm   * (dcosdrk) &
           + f_C   * gterm   * expterm_d * (-rik_hat)

      f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + prefactor*dri(1)
      f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + prefactor*dri(2)
      f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + prefactor*dri(3)
  
      f_nonbond(1,jparticle) = f_nonbond(1,jparticle) + prefactor*drj(1)
      f_nonbond(2,jparticle) = f_nonbond(2,jparticle) + prefactor*drj(2)
      f_nonbond(3,jparticle) = f_nonbond(3,jparticle) + prefactor*drj(3)

      f_nonbond(1,kparticle) = f_nonbond(1,kparticle) + prefactor*drk(1)
      f_nonbond(2,kparticle) = f_nonbond(2,kparticle) + prefactor*drk(2)
      f_nonbond(3,kparticle) = f_nonbond(3,kparticle) + prefactor*drk(3)

      ptens_nonbond(1,1) = ptens_nonbond(1,1) - prefactor*(rij(1)*drj(1) + rik(1)*drk(1))
      ptens_nonbond(1,2) = ptens_nonbond(1,2) - prefactor*(rij(1)*drj(2) + rik(1)*drk(2))
      ptens_nonbond(1,3) = ptens_nonbond(1,3) - prefactor*(rij(1)*drj(3) + rik(1)*drk(3))

      ptens_nonbond(2,1) = ptens_nonbond(2,1) - prefactor*(rij(2)*drj(1) + rik(2)*drk(1))
      ptens_nonbond(2,2) = ptens_nonbond(2,2) - prefactor*(rij(2)*drj(2) + rik(2)*drk(2))
      ptens_nonbond(2,3) = ptens_nonbond(2,3) - prefactor*(rij(2)*drj(3) + rik(2)*drk(3))

      ptens_nonbond(3,1) = ptens_nonbond(3,1) - prefactor*(rij(3)*drj(1) + rik(3)*drk(1))
      ptens_nonbond(3,2) = ptens_nonbond(3,2) - prefactor*(rij(3)*drj(2) + rik(3)*drk(2))
      ptens_nonbond(3,3) = ptens_nonbond(3,3) - prefactor*(rij(3)*drj(3) + rik(3)*drk(3))

    END DO
  END SUBROUTINE ter_zeta_ij_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_g(pot,costheta)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: costheta
    REAL(KIND=dp)                                       :: ter_g
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: c2,d2

    c = pot%tersoff%c
    d = pot%tersoff%d
    h = pot%tersoff%h
    c2 = c*c
    d2 = d*d
    ter_g = 0.0_dp
    ter_g = 1.0_dp + (c2/d2) - (c2)/(d2+(h-costheta)**2)
  END FUNCTION ter_g

  FUNCTION ter_g_d(pot,costheta)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: costheta
    REAL(KIND=dp)                                       :: ter_g_d
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: c2,d2,hc
    REAL(KIND=dp)                                       :: sintheta

    c = pot%tersoff%c
    d = pot%tersoff%d
    h = pot%tersoff%h
    c2 = c*c
    d2 = d*d
    hc = h-costheta

    sintheta = SQRT(1.0 - costheta**2)

    ter_g_d = 0.0_dp

    ! Still need d(theta)/dR
    !ter_g_d = 2.0_dp * c**2 * sintheta * (h - costheta) / &
    !          ( d**2 + (h-costheta)**2 )**2

    ! Still need d(costheta)/dR
    ter_g_d = (-2.0_dp * c2 * hc)/(d2 + hc**2)**2
  END FUNCTION ter_g_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ter_costheta_d(particle_set,cell,iparticle,jparticle,kparticle,dri,drj,drk)
    TYPE(pair_potential_single_type), POINTER           :: pot
    TYPE(particle_type), INTENT(IN)                     :: particle_set( : )
    TYPE(cell_type), POINTER                            :: cell
    INTEGER, INTENT(IN)                                 :: iparticle,jparticle,kparticle
    REAL(KIND=dp), INTENT(OUT)                          :: dri(3),drj(3),drk(3)
    REAL(KIND=dp)                                       :: rij(3),drij,rij_hat(3)
    REAL(KIND=dp)                                       :: rik(3),drik,rik_hat(3)
    REAL(KIND=dp)                                       :: costheta

    !Note that the center atom is at atom i not the usual atom j

    rij(:) = pbc(particle_set(jparticle)%r(:),particle_set(iparticle)%r(:),cell)
    drij = SQRT(DOT_PRODUCT(rij,rij))
    rij_hat(:) = rij(:)/drij

    rik(:) = pbc(particle_set(kparticle)%r(:),particle_set(iparticle)%r(:),cell)
    drik = SQRT(DOT_PRODUCT(rik,rik))
    rik_hat(:) = rik(:)/drik

    costheta = DOT_PRODUCT(rij,rik) / (drij*drik)
    IF(costheta<-1.0_dp) costheta=-1.0_dp
    IF(costheta>+1.0_dp) costheta=+1.0_dp

    drj(:) = (1.0_dp/(drij)) * ( rik_hat(:) - costheta*rij_hat(:) )
    drk(:) = (1.0_dp/(drik)) * ( rij_hat(:) - costheta*rik_hat(:) )
    dri(:) = - ( drj(:) + drk(:) )
  END SUBROUTINE ter_costheta_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !******************************************************************************
  SUBROUTINE get_rho_eam ( potparm, rab2, atom_a, atom_b, rho )
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    REAL(dp), INTENT(IN)                     :: rab2
    INTEGER, INTENT(IN)                      :: atom_a, atom_b
    REAL(dp), INTENT(INOUT)                  :: rho( : )

    INTEGER                                  :: index
    REAL(dp)                                 :: qq, rab, rhoi, rhoj

! Computation of rhotot

    rab = SQRT ( rab2 )
    index = INT ( rab / potparm % eam % drar ) + 1

    IF ( index > 5000 ) THEN
       index = 5000
    ELSEIF ( index < 1 ) THEN
       index = 1
    ENDIF

    qq = rab - potparm % eam % rval ( index )

    rhoi = potparm % eam % rho ( index ) + qq * potparm % eam % rhop ( index )
    rhoj = potparm % eam % rho ( index ) + qq * potparm % eam % rhop ( index )

    rho ( atom_a ) = rho ( atom_a ) + rhoi
    rho ( atom_b ) = rho ( atom_b ) + rhoj
  END SUBROUTINE get_rho_eam
  !******************************************************************************
  SUBROUTINE get_force_eam ( rab2, potparm, particle_set, atom_a, atom_b, f_eam )
    REAL(dp), INTENT(IN)                     :: rab2
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(particle_type), INTENT(IN)          :: particle_set( : )
    INTEGER, INTENT(IN)                      :: atom_a, atom_b
    REAL(dp), INTENT(OUT)                    :: f_eam

    INTEGER                                  :: index
    REAL(dp)                                 :: denspi, denspj, fcp, qq, rab

! compute scaler part of the force

    f_eam = 0._dp
    rab = SQRT ( rab2 )
    index = INT ( rab / potparm % eam % drar ) + 1

    IF ( index > 5000 ) THEN
       index = 5000
    ELSEIF ( index < 1 ) THEN
       index = 1
    ENDIF
    qq = rab - potparm % eam % rval ( index )

    denspi = potparm % eam % rhop ( index ) +  &
         qq * ( potparm % eam % rhop ( index + 1 ) - &
         potparm % eam % rhop ( index ) ) / potparm % eam % drar
    denspj = potparm % eam % rhop ( index ) + &
         qq *( potparm % eam % rhop ( index + 1 ) - &
         potparm % eam % rhop ( index ) ) / potparm % eam % drar

    fcp = denspj * particle_set ( atom_a ) % f_embed + &
         denspi * particle_set ( atom_b ) % f_embed
    f_eam =  fcp / rab
  END SUBROUTINE get_force_eam

  !******************************************************************************

  SUBROUTINE get_force_tersoff(pot,particle_set,cell,iparticle,jparticle,f_nonbond,&
                               ptens_nonbond)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(particle_type), INTENT(IN)          :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: iparticle, jparticle
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_nonbond,ptens_nonbond
    REAL(KIND=dp)                            :: rij(3),drij,rij_hat(3)
    REAL(KIND=dp)                            :: RpD,RmD
    REAL(KIND=dp)                            :: f_R1, f_R2
    REAL(KIND=dp)                            :: f_A1, f_A2
    REAL(KIND=dp)                            :: b_ij,b_ij_d
    REAL(KIND=dp)                            :: b_ji,b_ji_d
    REAL(KIND=dp)                            :: f_A,f_A_d
    REAL(KIND=dp)                            :: f_C,f_C_d
    REAL(KIND=dp)                            :: f_R,f_R_d
    REAL(KIND=dp)                            :: prefactor



    rij(:) = pbc(particle_set(jparticle)%r(:),particle_set(iparticle)%r(:),cell)
    drij = SQRT(DOT_PRODUCT(rij,rij))
    rij_hat(:) = rij(:) / drij 

    RpD = pot%tersoff%bigR + pot%tersoff%bigD
    IF(drij>RpD) RETURN

    b_ij  = 0.5_dp * ter_b_ij(pot,particle_set,cell,iparticle,jparticle)
    b_ij_d= 0.5_dp * ter_b_ij_d(pot,particle_set,cell,iparticle,jparticle)
    b_ji  = 0.5_dp * ter_b_ij(pot,particle_set,cell,jparticle,iparticle)
    b_ji_d= 0.5_dp * ter_b_ij_d(pot,particle_set,cell,jparticle,iparticle)
    f_A   = ter_f_A(pot,drij)
    f_A_d = ter_f_A_d(pot,drij)
    f_C   = ter_f_C(pot,drij)
    f_C_d = ter_f_C_d(pot,drij)
    f_R   = ter_f_R(pot,drij)
    f_R_d = ter_f_R_d(pot,drij)

    ! Lets do the easy one first, the repulsive term
    ! Note a_ij = 1.0_dp so just going to ignore it...
    f_R1 = f_C_d * f_R
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_R1*rij(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_R1*rij(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_R1*rij(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_R1*rij(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_R1*rij(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_R1*rij(3)

    ptens_nonbond(1,1) = ptens_nonbond(1,1) + f_R1*rij(1)*rij(1)
    ptens_nonbond(1,2) = ptens_nonbond(1,2) + f_R1*rij(1)*rij(2)
    ptens_nonbond(1,3) = ptens_nonbond(1,3) + f_R1*rij(1)*rij(3)
    ptens_nonbond(2,1) = ptens_nonbond(2,1) + f_R1*rij(2)*rij(1)
    ptens_nonbond(2,2) = ptens_nonbond(2,2) + f_R1*rij(2)*rij(2)
    ptens_nonbond(2,3) = ptens_nonbond(2,3) + f_R1*rij(2)*rij(3)
    ptens_nonbond(3,1) = ptens_nonbond(3,1) + f_R1*rij(3)*rij(1)
    ptens_nonbond(3,2) = ptens_nonbond(3,2) + f_R1*rij(3)*rij(2)
    ptens_nonbond(3,3) = ptens_nonbond(3,3) + f_R1*rij(3)*rij(3)

    f_R2 = f_C * f_R_d
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_R2*rij(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_R2*rij(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_R2*rij(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_R2*rij(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_R2*rij(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_R2*rij(3)

    ptens_nonbond(1,1) = ptens_nonbond(1,1) + f_R2*rij(1)*rij(1)
    ptens_nonbond(1,2) = ptens_nonbond(1,2) + f_R2*rij(1)*rij(2)
    ptens_nonbond(1,3) = ptens_nonbond(1,3) + f_R2*rij(1)*rij(3)
    ptens_nonbond(2,1) = ptens_nonbond(2,1) + f_R2*rij(2)*rij(1)
    ptens_nonbond(2,2) = ptens_nonbond(2,2) + f_R2*rij(2)*rij(2)
    ptens_nonbond(2,3) = ptens_nonbond(2,3) + f_R2*rij(2)*rij(3)
    ptens_nonbond(3,1) = ptens_nonbond(3,1) + f_R2*rij(3)*rij(1)
    ptens_nonbond(3,2) = ptens_nonbond(3,2) + f_R2*rij(3)*rij(2)
    ptens_nonbond(3,3) = ptens_nonbond(3,3) + f_R2*rij(3)*rij(3)

    ! Lets do the f_A1 piece derivative of F_C
    f_A1 = f_C_d * (b_ij+b_ji) * f_A
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_A1*rij(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_A1*rij(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_A1*rij(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_A1*rij(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_A1*rij(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_A1*rij(3)

    ptens_nonbond(1,1) = ptens_nonbond(1,1) + f_A1*rij(1)*rij(1)
    ptens_nonbond(1,2) = ptens_nonbond(1,2) + f_A1*rij(1)*rij(2)
    ptens_nonbond(1,3) = ptens_nonbond(1,3) + f_A1*rij(1)*rij(3)
    ptens_nonbond(2,1) = ptens_nonbond(2,1) + f_A1*rij(2)*rij(1)
    ptens_nonbond(2,2) = ptens_nonbond(2,2) + f_A1*rij(2)*rij(2)
    ptens_nonbond(2,3) = ptens_nonbond(2,3) + f_A1*rij(2)*rij(3)
    ptens_nonbond(3,1) = ptens_nonbond(3,1) + f_A1*rij(3)*rij(1)
    ptens_nonbond(3,2) = ptens_nonbond(3,2) + f_A1*rij(3)*rij(2)
    ptens_nonbond(3,3) = ptens_nonbond(3,3) + f_A1*rij(3)*rij(3)

    ! Lets do the f_A2 piece derivative of F_A
    f_A2 = f_C * (b_ij+b_ji) * f_A_d
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_A2*rij(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_A2*rij(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_A2*rij(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_A2*rij(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_A2*rij(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_A2*rij(3)

    ptens_nonbond(1,1) = ptens_nonbond(1,1) + f_A2*rij(1)*rij(1)
    ptens_nonbond(1,2) = ptens_nonbond(1,2) + f_A2*rij(1)*rij(2)
    ptens_nonbond(1,3) = ptens_nonbond(1,3) + f_A2*rij(1)*rij(3)
    ptens_nonbond(2,1) = ptens_nonbond(2,1) + f_A2*rij(2)*rij(1)
    ptens_nonbond(2,2) = ptens_nonbond(2,2) + f_A2*rij(2)*rij(2)
    ptens_nonbond(2,3) = ptens_nonbond(2,3) + f_A2*rij(2)*rij(3)
    ptens_nonbond(3,1) = ptens_nonbond(3,1) + f_A2*rij(3)*rij(1)
    ptens_nonbond(3,2) = ptens_nonbond(3,2) + f_A2*rij(3)*rij(2)
    ptens_nonbond(3,3) = ptens_nonbond(3,3) + f_A2*rij(3)*rij(3)

    ! Lets do the f_A3 piece derivative of b_ij
    prefactor = f_C * b_ij_d * f_A    ! Note need to do d(Zeta_ij)/dR
    CALL ter_zeta_ij_d(pot,particle_set,cell,iparticle,jparticle,f_nonbond,ptens_nonbond,prefactor)
    prefactor = f_C * b_ji_d * f_A    ! Note need to do d(Zeta_ij)/dR
    CALL ter_zeta_ij_d(pot,particle_set,cell,jparticle,iparticle,f_nonbond,ptens_nonbond,prefactor)

  END SUBROUTINE get_force_tersoff

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
    TYPE(particle_type), POINTER, OPTIONAL   :: shell_particle_set( : )
    TYPE(particle_type), POINTER, OPTIONAL   :: core_particle_set( : )
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'bonded_correct_gaussian'
    REAL(KIND=dp), PARAMETER :: ac1 = 0.254829592_dp, ac2 = -0.284496736_dp, &
      ac3 = 1.421413741_dp, ac4 = -1.453152027_dp, ac5 = 1.061405429_dp, &
      pc = 0.3275911_dp

    INTEGER                                  :: group, handle, i, iatom, &
                                                iexl, ikind, j, &
                                                natoms_per_kind, nkind,&
                                                shell_i_index, shell_j_index
    LOGICAL                                  :: i_is_shell, j_is_shell, shell_adiabatic
    REAL(KIND=dp)                            :: alpha, arg, const, dij, &
                                                e_arg_arg, eps0, errf, flops, &
                                                fscalar, idij, ifourpieps0, &
                                                qci, qcj, qi, qj, qsi, qsj, rijsq, tc
    REAL(KIND=dp), DIMENSION(3)              :: fij_com, rci, rcj, ri, rij, rj, rsi, rsj
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

    !
    !
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

! Exclusion of shell-core interaction  (???????????)
          IF(i_is_shell .AND. .FALSE.) THEN
             rij(1) = rci(1) - rsi(1)
             rij(2) = rci(2) - rsi(2)
             rij(3) = rci(3) - rsi(3)
             rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
             IF(rijsq > EPSILON(1.0_dp)) THEN
               idij = 1.0_dp / SQRT ( rijsq )
               dij = rijsq * idij
               arg = alpha * dij
               e_arg_arg = EXP ( -arg ** 2 )
               tc = 1.0_dp / ( 1.0_dp + pc * arg )

               ! defining errf=1-erfc
               errf = 1.0_dp &
                  - ((((ac5*tc+ac4)*tc+ac3)*tc+ac2)*tc+ac1) * tc * e_arg_arg

               ! getting the potential
               v_bonded_corr = v_bonded_corr - qci*qsi*idij*errf

               ! subtracting the force from the total force
               fscalar = ifourpieps0 * qci * qsi * idij ** 2 &
                    * ( idij * errf - const * e_arg_arg  )

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

          END IF
       ENDDO
    ENDDO

    ! the factor of 1/2 comes from double counting in the exclusion list
    v_bonded_corr = v_bonded_corr * ifourpieps0 * 0.5_dp
    CALL mp_sum ( v_bonded_corr, group )

    flops = flops * 1.0E-6_dp
    CALL timestop ( flops, handle )

  END SUBROUTINE bonded_correct_gaussian

  !******************************************************************************

END MODULE fist_nonbond_force
