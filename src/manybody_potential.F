!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Efficient tersoff implementation and general "lifting" of manybody_potential module
!>      12.2007 [tlaino] - Splitting manybody module : In this module we should only
!>                         keep the main routines for computing energy and forces of
!>                         manybody potentials. Each potential should have his own module!
!> \author CJM, I-Feng W. Kuo, Teodoro Laino
! *****************************************************************************
MODULE manybody_potential

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: eam_type,&
                                             fist_nonbond_env_get,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE manybody_eam,                    ONLY: get_force_eam
  USE manybody_quip,                   ONLY: quip_add_force_virial,&
                                             quip_energy_store_force_virial
  USE manybody_siepmann,               ONLY: destroy_siepmann_arrays,&
                                             print_nr_oh_siepmann,&
                                             setup_siepmann_arrays,&
                                             siepmann_energy,&
                                             siepmann_forces_v2,&
                                             siepmann_forces_v3
  USE manybody_tersoff,                ONLY: destroy_tersoff_arrays,&
                                             setup_tersoff_arrays,&
                                             tersoff_energy,&
                                             tersoff_forces
  USE mathlib,                         ONLY: matvec_3x3
  USE message_passing,                 ONLY: mp_sum
  USE pair_potential_types,            ONLY: &
       ea_type, eam_pot_type, pair_potential_pp_type, &
       pair_potential_single_type, quip_type, siepmann_pot_type, &
       siepmann_type, tersoff_pot_type, tersoff_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: energy_manybody
  PUBLIC :: force_nonbond_manybody
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_potential'

CONTAINS

! *****************************************************************************
!> \brief computes the embedding contribution to the energy
!> \param fist_nonbond_env ...
!> \param atomic_kind_set ...
!> \param local_particles ...
!> \param particle_set ...
!> \param cell ...
!> \param pot_manybody ...
!> \param para_env ...
!> \param mm_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      tlaino [2007] - New algorithm for tersoff potential
!> \author CJM, I-Feng W. Kuo, Teodoro Laino
! *****************************************************************************
  SUBROUTINE energy_manybody ( fist_nonbond_env, atomic_kind_set, local_particles,&
       particle_set, cell, pot_manybody, para_env, mm_section, error )

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), INTENT(INOUT)                  :: pot_manybody
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'energy_manybody', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, i, iend, ifirst, igrp, ikind, ilast, &
      ilist, indexa, ipair, iparticle, iparticle_local, istart, iunique, &
      jkind, junique, mpair, nkinds, nloc_size, npairs, nparticle, &
      nparticle_local, nr_oh, nunique, stat
    INTEGER, DIMENSION(:), POINTER           :: glob_loc_list_a, work_list
    INTEGER, DIMENSION(:, :), POINTER        :: glob_loc_list, list, sort_list
    LOGICAL                                  :: any_quip, any_siepmann, &
                                                any_tersoff, failure
    REAL(KIND=dp)                            :: drij, embed, pot_loc, &
                                                pot_quip, qr, rab2_max, rij(3)
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, cvi
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: glob_cell_v
    REAL(KIND=dp), POINTER                   :: fembed( : )
    TYPE(eam_pot_type), POINTER              :: eam
    TYPE(eam_type), DIMENSION(:), POINTER    :: eam_data
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    TYPE(siepmann_pot_type), POINTER         :: siepmann
    TYPE(tersoff_pot_type), POINTER          :: tersoff

    NULLIFY(eam, siepmann, tersoff)
    any_tersoff = .FALSE.
    any_siepmann = .FALSE.
    any_quip = .FALSE.
    failure = .FALSE.
    CALL timeset ( routineN, handle )
    CALL fist_nonbond_env_get ( fist_nonbond_env, r_last_update_pbc=r_last_update_pbc,&
         potparm = potparm , eam_data=eam_data, error=error)
    ! EAM requires a single loop
    DO ikind = 1, SIZE ( atomic_kind_set )
       pot => potparm %pot ( ikind, ikind ) % pot
       DO i = 1, SIZE(pot%type)
         IF(pot%type(i) /= ea_type) CYCLE
         eam => pot%set(i)%eam
         nparticle = SIZE ( particle_set )
         ALLOCATE(fembed(nparticle), stat=stat )
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         fembed(:) = 0._dp
         CPPostcondition(ASSOCIATED(eam_data),cp_failure_level,routineP,error,failure)
         ! computation of embedding function and energy
         nparticle_local = local_particles%n_el(ikind)
         DO iparticle_local=1,nparticle_local
            iparticle = local_particles%list(ikind)%array(iparticle_local)
            indexa = INT ( eam_data(iparticle)%rho / eam%drhoar) + 1
            IF (indexa > eam%npoints-1) indexa = eam%npoints-1
            qr = eam_data(iparticle)%rho  -  eam%rhoval(indexa)

            embed = eam%frho(indexa) +  qr*eam%frhop(indexa)
            fembed (iparticle) = eam%frhop(indexa) + qr*(eam%frhop(indexa+1)-eam%frhop(indexa))/eam%drhoar

            pot_manybody = pot_manybody + embed
         END DO
         ! communicate data
         CALL mp_sum(fembed, para_env%group)
         DO iparticle=1,nparticle
            IF (particle_set(iparticle)%atomic_kind%kind_number==ikind) THEN
               eam_data(iparticle)%f_embed = fembed(iparticle)
            END IF
         END DO

         DEALLOCATE ( fembed, stat=stat )
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO
    END DO
    ! Other manybody potential
    DO ikind = 1, SIZE ( atomic_kind_set )
       DO jkind =  ikind, SIZE ( atomic_kind_set )
          pot => potparm %pot ( ikind, jkind ) % pot
          any_tersoff = any_tersoff .OR. ANY(pot%type==tersoff_type)
          any_quip = any_quip .OR. ANY(pot%type==quip_type)
          any_siepmann = any_siepmann .OR. ANY(pot%type==siepmann_type)
       END DO
    END DO
    CALL fist_nonbond_env_get(fist_nonbond_env,nonbonded=nonbonded,natom_types=nkinds,error=error)
    ! QUIP
    IF (any_quip) THEN
      CALL quip_energy_store_force_virial(particle_set, cell, atomic_kind_set, potparm, &
         fist_nonbond_env, pot_quip, para_env, error=error)
      pot_manybody = pot_manybody + pot_quip
    ENDIF
    ! TERSOFF
    IF (any_tersoff) THEN 
       NULLIFY(glob_loc_list, glob_cell_v, glob_loc_list_a)
       CALL setup_tersoff_arrays(nonbonded, potparm, glob_loc_list, glob_cell_v, glob_loc_list_a, cell, error)
       DO ilist=1,nonbonded%nlists
          neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
          npairs=neighbor_kind_pair%npairs
          IF (npairs ==0) CYCLE
          Kind_Group_Loop: DO igrp = 1, neighbor_kind_pair%ngrp_kind
             istart =  neighbor_kind_pair%grp_kind_start(igrp)
             iend   =  neighbor_kind_pair%grp_kind_end(igrp)
             ikind  =  neighbor_kind_pair%ij_kind(1,igrp)
             jkind  =  neighbor_kind_pair%ij_kind(2,igrp)
             list   => neighbor_kind_pair%list
             cvi    =  neighbor_kind_pair%cell_vector
             pot    => potparm %pot ( ikind, jkind ) % pot
             DO i = 1, SIZE(pot%type)
                IF(pot%type(i) /= tersoff_type) CYCLE
                rab2_max = pot%set(i)%tersoff%rcutsq
                CALL matvec_3x3(cell_v, cell%hmat,cvi)
                pot     => potparm%pot( ikind, jkind )%pot
                tersoff => pot%set(i)%tersoff
                npairs  = iend-istart+1
                IF (npairs /=0) THEN
                   ALLOCATE(sort_list(2,npairs),work_list (npairs),stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   sort_list = list(:,istart:iend)
                   ! Sort the list of neighbors, this increases the efficiency for single
                   ! potential contributions
                   CALL sort(sort_list(1,:),npairs,work_list)
                   DO ipair = 1, npairs
                      work_list(ipair)=sort_list(2,work_list(ipair))
                   END DO
                   sort_list(2,:) = work_list
                   ! find number of unique elements of array index 1
                   nunique = 1
                   DO ipair = 1, npairs-1
                      IF (sort_list(1,ipair+1)/=sort_list(1,ipair)) nunique = nunique + 1
                   END DO
                   ipair   = 1
                   junique = sort_list(1,ipair)
                   ifirst  = 1
                   DO iunique = 1, nunique
                      atom_a = junique
                      IF (glob_loc_list_a(ifirst)>atom_a) CYCLE
                      DO mpair = ifirst, SIZE(glob_loc_list_a)
                         IF (glob_loc_list_a(mpair)==atom_a) EXIT
                      END DO
                      ifirst = mpair
                      DO mpair = ifirst, SIZE(glob_loc_list_a)
                         IF (glob_loc_list_a(mpair)/=atom_a) EXIT
                      END DO
                      ilast = mpair-1
                      nloc_size = 0
                      IF (ifirst/=0) nloc_size = ilast-ifirst+1
                      DO WHILE (ipair<=npairs)
                         IF (sort_list(1,ipair)/=junique) EXIT
                         atom_b = sort_list(2,ipair)
                         ! Energy terms
                         pot_loc = 0.0_dp
                         rij(:) = r_last_update_pbc(atom_b)%r(:)-r_last_update_pbc(atom_a)%r(:)+cell_v
                         drij = DOT_PRODUCT(rij,rij)
                         ipair = ipair + 1
                         IF(drij>rab2_max) CYCLE
                         drij = SQRT(drij)
                         CALL tersoff_energy(pot_loc, tersoff, r_last_update_pbc, atom_a, atom_b, nloc_size,&
                              glob_loc_list(:,ifirst:ilast),glob_cell_v(:,ifirst:ilast), cell_v, drij)
                         pot_manybody = pot_manybody + 0.5_dp*pot_loc
                      END DO
                      ifirst = ilast + 1
                      IF (ipair<=npairs) junique = sort_list(1,ipair)
                   END DO
                   DEALLOCATE(sort_list,work_list,stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                END IF
             END DO
          END DO Kind_Group_Loop
       END DO
       CALL destroy_tersoff_arrays(glob_loc_list, glob_cell_v, glob_loc_list_a, error)
    END IF

    !SIEPMANN POTENTIAL
    IF (any_siepmann) THEN 
       NULLIFY(glob_loc_list, glob_cell_v, glob_loc_list_a)
       CALL setup_siepmann_arrays(nonbonded, potparm, glob_loc_list, glob_cell_v, glob_loc_list_a, cell, error)
       nr_oh = 0
       DO ilist=1,nonbonded%nlists
          neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
          npairs=neighbor_kind_pair%npairs
          IF (npairs ==0) CYCLE
          Kind_Group_Loop_2: DO igrp = 1, neighbor_kind_pair%ngrp_kind
             istart =  neighbor_kind_pair%grp_kind_start(igrp)
             iend   =  neighbor_kind_pair%grp_kind_end(igrp)
             ikind  =  neighbor_kind_pair%ij_kind(1,igrp)
             jkind  =  neighbor_kind_pair%ij_kind(2,igrp)
             list   => neighbor_kind_pair%list
             cvi    =  neighbor_kind_pair%cell_vector
             pot    => potparm %pot ( ikind, jkind ) % pot
             DO i = 1, SIZE(pot%type)
                IF(pot%type(i) /= siepmann_type) CYCLE
                rab2_max = pot%set(i)%siepmann%rcutsq
                CALL matvec_3x3(cell_v, cell%hmat,cvi)
                pot     => potparm%pot( ikind, jkind )%pot
                siepmann => pot%set(i)%siepmann
                npairs  = iend-istart+1
                IF (npairs /=0) THEN
                   ALLOCATE(sort_list(2,npairs),work_list (npairs),stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   sort_list = list(:,istart:iend)
                   ! Sort the list of neighbors, this increases the efficiency for single
                   ! potential contributions
                   CALL sort(sort_list(1,:),npairs,work_list)
                   DO ipair = 1, npairs
                      work_list(ipair)=sort_list(2,work_list(ipair))
                   END DO
                   sort_list(2,:) = work_list
                   ! find number of unique elements of array index 1
                   nunique = 1
                   DO ipair = 1, npairs-1
                      IF (sort_list(1,ipair+1)/=sort_list(1,ipair)) nunique = nunique + 1
                   END DO
                   ipair   = 1
                   junique = sort_list(1,ipair)
                   ifirst  = 1
                   DO iunique = 1, nunique
                      atom_a = junique
                      IF (glob_loc_list_a(ifirst)>atom_a) CYCLE
                      DO mpair = ifirst, SIZE(glob_loc_list_a)
                         IF (glob_loc_list_a(mpair)==atom_a) EXIT
                      END DO
                      ifirst = mpair
                      DO mpair = ifirst, SIZE(glob_loc_list_a)
                         IF (glob_loc_list_a(mpair)/=atom_a) EXIT
                      END DO
                      ilast = mpair-1
                      nloc_size = 0
                      IF (ifirst/=0) nloc_size = ilast-ifirst+1
                      DO WHILE (ipair<=npairs)
                         IF (sort_list(1,ipair)/=junique) EXIT
                         atom_b = sort_list(2,ipair)
                         ! Energy terms
                         pot_loc = 0.0_dp
                         rij(:) = r_last_update_pbc(atom_b)%r(:)-r_last_update_pbc(atom_a)%r(:)+cell_v
                         drij = DOT_PRODUCT(rij,rij)
                         ipair = ipair + 1
                         IF(drij>rab2_max) CYCLE
                         drij = SQRT(drij)
                         CALL siepmann_energy(pot_loc, siepmann, r_last_update_pbc, atom_a, atom_b, nloc_size,&
                              glob_loc_list(:,ifirst:ilast),glob_cell_v(:,ifirst:ilast), cell_v, cell, drij,&
                              particle_set, nr_oh)
                         pot_manybody = pot_manybody + pot_loc 
                      END DO
                      ifirst = ilast + 1
                      IF (ipair<=npairs) junique = sort_list(1,ipair)
                   END DO
                   DEALLOCATE(sort_list,work_list,stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                END IF
             END DO
          END DO Kind_Group_Loop_2
       END DO
       CALL destroy_siepmann_arrays(glob_loc_list, glob_cell_v, glob_loc_list_a, error)
       IF(siepmann%allow_oh_formation) THEN
         CALL print_nr_oh_siepmann(nr_oh,mm_section,para_env,error)
       ENDIF
    END IF

    CALL timestop ( handle )
  END SUBROUTINE energy_manybody

! *****************************************************************************
!> \brief ...
!> \param fist_nonbond_env ...
!> \param particle_set ...
!> \param cell ...
!> \param f_nonbond ...
!> \param pv_nonbond ...
!> \param use_virial ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      Fast implementation of the tersoff potential - [tlaino] 2007
!> \author I-Feng W. Kuo, Teodoro Laino
! *****************************************************************************
  SUBROUTINE force_nonbond_manybody ( fist_nonbond_env, particle_set, cell,  &
       f_nonbond, pv_nonbond, use_virial, error )

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_nonbond, pv_nonbond
    LOGICAL, INTENT(IN)                      :: use_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_nonbond_manybody', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, i, i_a, i_b, iend, ifirst, igrp, &
      ikind, ilast, ilist, ipair, istart, iunique, jkind, junique, kind_a, &
      kind_b, mpair, nkinds, nloc_size, npairs, nunique, stat
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: eam_kinds_index
    INTEGER, DIMENSION(:), POINTER           :: glob_loc_list_a, work_list
    INTEGER, DIMENSION(:, :), POINTER        :: glob_loc_list, list, sort_list
    LOGICAL                                  :: any_quip, any_siepmann, &
                                                any_tersoff, failure
    REAL(KIND=dp) :: f_eam, fac, fr(3), ptens11, ptens12, ptens13, ptens21, &
      ptens22, ptens23, ptens31, ptens32, ptens33, rab(3), rab2, rab2_max, &
      rtmp(3)
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, cvi
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: glob_cell_v
    TYPE(eam_pot_type), POINTER              :: eam_a, eam_b
    TYPE(eam_type), DIMENSION(:), POINTER    :: eam_data
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    TYPE(siepmann_pot_type), POINTER         :: siepmann
    TYPE(tersoff_pot_type), POINTER          :: tersoff

    failure     = .FALSE.
    any_tersoff = .FALSE.
    any_quip = .FALSE.
    any_siepmann = .FALSE.
    CALL timeset ( routineN, handle )
    NULLIFY(eam_a, eam_b, tersoff, siepmann)

    CALL fist_nonbond_env_get(fist_nonbond_env,nonbonded=nonbonded,potparm=potparm,&
         natom_types=nkinds, eam_data=eam_data, r_last_update_pbc=r_last_update_pbc,error=error)

    ! Initializing the potential energy, pressure tensor and force
    IF (use_virial) THEN
       ptens11 = 0.0_dp ; ptens12 = 0.0_dp ; ptens13 = 0.0_dp
       ptens21 = 0.0_dp ; ptens22 = 0.0_dp ; ptens23 = 0.0_dp
       ptens31 = 0.0_dp ; ptens32 = 0.0_dp ; ptens33 = 0.0_dp
    END IF

    nkinds = SIZE(potparm%pot,1)
    ALLOCATE(eam_kinds_index(nkinds,nkinds),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    eam_kinds_index = -1
    DO ikind = 1, nkinds
       DO jkind = ikind, nkinds
          DO i = 1, SIZE(potparm%pot(ikind,jkind)%pot%type)
             IF (potparm%pot(ikind,jkind)%pot%type(i)==ea_type) THEN
                ! At the moment we allow only 1 EAM per each kinds pair..
                CPPostcondition(eam_kinds_index(ikind,jkind)==-1,cp_failure_level,routineP,error,failure)
                CPPostcondition(eam_kinds_index(jkind,ikind)==-1,cp_failure_level,routineP,error,failure)
                eam_kinds_index(ikind,jkind) = i
                eam_kinds_index(jkind,ikind) = i
             END IF
          END DO
       END DO
    END DO

    ! QUIP
    IF (use_virial) &
       CALL quip_add_force_virial(fist_nonbond_env, f_nonbond, pv_nonbond, error=error)

    ! starting the force loop
    DO ilist=1,nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
       npairs=neighbor_kind_pair%npairs
       IF (npairs ==0) CYCLE
       Kind_Group_Loop1: DO igrp = 1, neighbor_kind_pair%ngrp_kind
          istart  = neighbor_kind_pair%grp_kind_start(igrp)
          iend    = neighbor_kind_pair%grp_kind_end(igrp)
          ikind   = neighbor_kind_pair%ij_kind(1,igrp)
          jkind   = neighbor_kind_pair%ij_kind(2,igrp)
          list  => neighbor_kind_pair%list
          cvi   =  neighbor_kind_pair%cell_vector
          pot   => potparm %pot( ikind, jkind )%pot
          IF (pot%no_mb) CYCLE
          rab2_max = pot%rcutsq
          CALL matvec_3x3(cell_v, cell%hmat,cvi)
          any_tersoff = any_tersoff.OR.ANY(pot%type==tersoff_type)
          any_siepmann = any_siepmann.OR.ANY(pot%type==siepmann_type)
          i = eam_kinds_index(ikind,jkind)
          IF (i==-1) CYCLE
          ! EAM
          CPPostcondition(ASSOCIATED(eam_data),cp_failure_level,routineP,error,failure)
          DO ipair = istart, iend
             atom_a = list(1,ipair)
             atom_b = list(2,ipair)
             fac    = 1.0_dp
             IF (atom_a==atom_b) fac = 0.5_dp
             kind_a = particle_set(atom_a)%atomic_kind%kind_number
             kind_b = particle_set(atom_b)%atomic_kind%kind_number
             i_a    = eam_kinds_index(kind_a,kind_a)
             i_b    = eam_kinds_index(kind_b,kind_b)
             eam_a => potparm%pot(kind_a,kind_a)%pot%set(i_a)%eam
             eam_b => potparm%pot(kind_b,kind_b)%pot%set(i_b)%eam

             !set this outside the potential type in case need multiple potentials
             !Do everything necessary for EAM here
             rab  = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
             rab  = rab+cell_v
             rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
             IF (rab2 <= rab2_max) THEN
                CALL get_force_eam ( rab2, eam_a, eam_b, eam_data, atom_a, atom_b, f_eam )
                f_eam = f_eam * fac

                fr(1) = - f_eam * rab(1)
                fr(2) = - f_eam * rab(2)
                fr(3) = - f_eam * rab(3)
                f_nonbond(1,atom_a) = f_nonbond(1,atom_a) - fr(1)
                f_nonbond(2,atom_a) = f_nonbond(2,atom_a) - fr(2)
                f_nonbond(3,atom_a) = f_nonbond(3,atom_a) - fr(3)

                f_nonbond(1,atom_b) = f_nonbond(1,atom_b) + fr(1)
                f_nonbond(2,atom_b) = f_nonbond(2,atom_b) + fr(2)
                f_nonbond(3,atom_b) = f_nonbond(3,atom_b) + fr(3)
                IF (use_virial) THEN
                   ptens11 = ptens11 + rab(1)*fr(1)
                   ptens21 = ptens21 + rab(2)*fr(1)
                   ptens31 = ptens31 + rab(3)*fr(1)
                   ptens12 = ptens12 + rab(1)*fr(2)
                   ptens22 = ptens22 + rab(2)*fr(2)
                   ptens32 = ptens32 + rab(3)*fr(2)
                   ptens13 = ptens13 + rab(1)*fr(3)
                   ptens23 = ptens23 + rab(2)*fr(3)
                   ptens33 = ptens33 + rab(3)*fr(3)
                END IF
             ENDIF
          END DO
       END DO Kind_Group_Loop1
    END DO
    DEALLOCATE (eam_kinds_index, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Special way of handling the tersoff potential..
    IF (any_tersoff) THEN
       NULLIFY(glob_loc_list, glob_cell_v, glob_loc_list_a)
       CALL setup_tersoff_arrays(nonbonded, potparm, glob_loc_list, glob_cell_v, glob_loc_list_a, cell, error)
       DO ilist=1,nonbonded%nlists
          neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
          npairs=neighbor_kind_pair%npairs
          IF (npairs ==0) CYCLE
          Kind_Group_Loop2: DO igrp = 1, neighbor_kind_pair%ngrp_kind
             istart  = neighbor_kind_pair%grp_kind_start(igrp)
             iend    = neighbor_kind_pair%grp_kind_end(igrp)
             ikind   = neighbor_kind_pair%ij_kind(1,igrp)
             jkind   = neighbor_kind_pair%ij_kind(2,igrp)
             list  => neighbor_kind_pair%list
             cvi   =  neighbor_kind_pair%cell_vector
             pot   => potparm %pot( ikind, jkind )%pot

             IF (pot%no_mb) CYCLE
             rab2_max = pot%rcutsq
             CALL matvec_3x3(cell_v, cell%hmat,cvi)
             DO i = 1, SIZE(pot%type)
                ! TERSOFF
                IF (pot%type(i)==tersoff_type) THEN
                   npairs = iend-istart+1
                   tersoff => pot%set(i)%tersoff
                   ALLOCATE(sort_list(2,npairs),work_list (npairs),stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   sort_list = list(:,istart:iend)
                   ! Sort the list of neighbors, this increases the efficiency for single
                   ! potential contributions
                   CALL sort(sort_list(1,:),npairs,work_list)
                   DO ipair = 1, npairs
                      work_list(ipair)=sort_list(2,work_list(ipair))
                   END DO
                   sort_list(2,:) = work_list
                   ! find number of unique elements of array index 1
                   nunique = 1
                   DO ipair = 1, npairs-1
                      IF (sort_list(1,ipair+1)/=sort_list(1,ipair)) nunique = nunique + 1
                   END DO
                   ipair   = 1
                   junique = sort_list(1,ipair)
                   ifirst  = 1
                   DO iunique = 1, nunique
                      atom_a = junique
                      IF (glob_loc_list_a(ifirst)>atom_a) CYCLE
                      DO mpair = ifirst, SIZE(glob_loc_list_a)
                         IF (glob_loc_list_a(mpair)==atom_a) EXIT
                      END DO
                      ifirst = mpair
                      DO mpair = ifirst, SIZE(glob_loc_list_a)
                         IF (glob_loc_list_a(mpair)/=atom_a) EXIT
                      END DO
                      ilast = mpair-1
                      nloc_size = 0
                      IF (ifirst/=0) nloc_size = ilast-ifirst+1
                      DO WHILE (ipair<=npairs)
                         IF (sort_list(1,ipair)/=junique) EXIT
                         atom_b = sort_list(2,ipair)
                         ! Derivative terms
                         rtmp = r_last_update_pbc(atom_b)%r(:)-r_last_update_pbc(atom_a)%r(:) + cell_v
                         ipair = ipair + 1
                         IF (DOT_PRODUCT(rtmp,rtmp)<=tersoff%rcutsq) THEN
                             CALL tersoff_forces ( tersoff, r_last_update_pbc, cell_v,&
                                  nloc_size, glob_loc_list(:,ifirst:ilast), glob_cell_v(:,ifirst:ilast),&
                                  atom_a, atom_b, f_nonbond, pv_nonbond, use_virial, tersoff%rcutsq)
                         END IF
                      END DO
                      ifirst = ilast + 1
                      IF (ipair<=npairs) junique = sort_list(1,ipair)
                   END DO
                   DEALLOCATE(sort_list,work_list,stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                END IF
             END DO
          END DO Kind_Group_Loop2
       END DO
       CALL destroy_tersoff_arrays(glob_loc_list, glob_cell_v, glob_loc_list_a, error)
    END IF
    ! Special way of handling the siepmann potential..
    IF (any_siepmann) THEN 
       NULLIFY(glob_loc_list, glob_cell_v, glob_loc_list_a)
       CALL setup_siepmann_arrays(nonbonded, potparm, glob_loc_list, glob_cell_v, glob_loc_list_a, cell, error)
       DO ilist=1,nonbonded%nlists
          neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
          npairs=neighbor_kind_pair%npairs
          IF (npairs ==0) CYCLE
          Kind_Group_Loop3: DO igrp = 1, neighbor_kind_pair%ngrp_kind
             istart  = neighbor_kind_pair%grp_kind_start(igrp)
             iend    = neighbor_kind_pair%grp_kind_end(igrp)
             ikind   = neighbor_kind_pair%ij_kind(1,igrp)
             jkind   = neighbor_kind_pair%ij_kind(2,igrp)
             list  => neighbor_kind_pair%list
             cvi   =  neighbor_kind_pair%cell_vector
             pot   => potparm %pot( ikind, jkind )%pot

             IF (pot%no_mb) CYCLE
             rab2_max = pot%rcutsq
             CALL matvec_3x3(cell_v, cell%hmat,cvi)
             DO i = 1, SIZE(pot%type)
                ! SIEPMANN
                IF (pot%type(i)==siepmann_type) THEN
                   npairs = iend-istart+1
                   siepmann => pot%set(i)%siepmann
                   ALLOCATE(sort_list(2,npairs),work_list (npairs),stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   sort_list = list(:,istart:iend)
                   ! Sort the list of neighbors, this increases the efficiency for single
                   ! potential contributions
                   CALL sort(sort_list(1,:),npairs,work_list)
                   DO ipair = 1, npairs
                      work_list(ipair)=sort_list(2,work_list(ipair))
                   END DO
                   sort_list(2,:) = work_list
                   ! find number of unique elements of array index 1
                   nunique = 1
                   DO ipair = 1, npairs-1
                      IF (sort_list(1,ipair+1)/=sort_list(1,ipair)) nunique = nunique + 1
                   END DO
                   ipair   = 1
                   junique = sort_list(1,ipair)
                   ifirst  = 1
                   DO iunique = 1, nunique
                      atom_a = junique
                      IF (glob_loc_list_a(ifirst)>atom_a) CYCLE
                      DO mpair = ifirst, SIZE(glob_loc_list_a)
                         IF (glob_loc_list_a(mpair)==atom_a) EXIT
                      END DO
                      ifirst = mpair
                      DO mpair = ifirst, SIZE(glob_loc_list_a)
                         IF (glob_loc_list_a(mpair)/=atom_a) EXIT
                      END DO
                      ilast = mpair-1
                      nloc_size = 0
                      IF (ifirst/=0) nloc_size = ilast-ifirst+1
                      DO WHILE (ipair<=npairs)
                         IF (sort_list(1,ipair)/=junique) EXIT
                         atom_b = sort_list(2,ipair)
                         ! Derivative terms
                         rtmp = r_last_update_pbc(atom_b)%r(:)-r_last_update_pbc(atom_a)%r(:) + cell_v
                         ipair = ipair + 1
                         IF (DOT_PRODUCT(rtmp,rtmp)<=siepmann%rcutsq) THEN
                             CALL siepmann_forces_v2 ( siepmann, r_last_update_pbc, cell_v, cell,&
                                  nloc_size, glob_loc_list(:,ifirst:ilast), glob_cell_v(:,ifirst:ilast),&
                                  atom_a, atom_b, f_nonbond, pv_nonbond, use_virial, siepmann%rcutsq,&
                                  particle_set, error)
                             CALL siepmann_forces_v3( siepmann, r_last_update_pbc, cell_v,&
                                  nloc_size, glob_loc_list(:,ifirst:ilast), glob_cell_v(:,ifirst:ilast),&
                                  atom_a, atom_b, f_nonbond, pv_nonbond, use_virial, siepmann%rcutsq, &
                                  cell, particle_set, error)
                         END IF
                      END DO
                      ifirst = ilast + 1
                      IF (ipair<=npairs) junique = sort_list(1,ipair)
                   END DO
                   DEALLOCATE(sort_list,work_list,stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                END IF
             END DO
          END DO Kind_Group_Loop3
       END DO
       CALL destroy_siepmann_arrays(glob_loc_list, glob_cell_v, glob_loc_list_a, error)
    END IF
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
  END SUBROUTINE force_nonbond_manybody

END MODULE manybody_potential

