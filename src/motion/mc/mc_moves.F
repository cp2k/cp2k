!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief the various moves in Monte Carlo (MC) simulations, including
!>      change of internal conformation, translation of a molecule, rotation
!>      of a molecule, and changing the size of the simulation box
!> \par History
!>      none
!> \author Matthew J. McGrath  (10.16.2003)
! *****************************************************************************
MODULE mc_moves
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cell_types,                      ONLY: cell_clone,&
                                             cell_create,&
                                             cell_release,&
                                             cell_type,&
                                             get_cell
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_set,&
                                             cp_subsys_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_fist_force
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size
  USE mathconstants,                   ONLY: pi
  USE mc_coordinates,                  ONLY: check_for_overlap,&
                                             create_discrete_array,&
                                             generate_cbmc_swap_config,&
                                             get_center_of_mass
  USE mc_types,                        ONLY: get_mc_molecule_info,&
                                             get_mc_par,&
                                             mc_ekin_type,&
                                             mc_molecule_info_type,&
                                             mc_moves_type,&
                                             mc_simpar_type
  USE md_run,                          ONLY: qs_mol_dyn
  USE message_passing,                 ONLY: mp_bcast
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: bend_type,&
                                             bond_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type,&
                                             torsion_type
  USE parallel_rng_types,              ONLY: next_random_number,&
                                             rng_stream_type
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: angstrom
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../../common/cp_common_uses.f90"

  IMPLICIT NONE
 
  PRIVATE

  PRIVATE  :: change_bond_angle,change_bond_length,depth_first_search,&
      change_dihedral

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_moves'

  PUBLIC :: mc_conformation_change,mc_molecule_translation,&
            mc_molecule_rotation,mc_volume_move,mc_avbmc_move,&
            mc_hmc_move

CONTAINS

! *****************************************************************************
!> \brief essentially performs a depth-first search of the molecule structure
!>      to find all atoms connected to a specific atom excluding one branch...
!>      for instance, if water is labelled 1-2-3 for O-H-H, calling this
!>      routine with current_atom=1,avoid_atom=2 returns the array
!>      atom=(0,0,1)
!> \param current_atom the atom whose connections we're looking at
!> \param avoid_atom the atom whose direction the search is not supposed to go
!> \param connectivity an array telling us the neigbors of all atoms
!> \param atom the array that tells us if one can get to a given atom by
!>        starting at current_atom and not going through avoid_atom...0 is no,
!>        1 is yes
!> \author MJM
! *****************************************************************************
  RECURSIVE SUBROUTINE depth_first_search ( current_atom,avoid_atom,&
      connectivity,atom)

    INTEGER, INTENT(IN)                      :: current_atom, avoid_atom
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: connectivity
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: atom

    INTEGER                                  :: iatom

      DO iatom=1,6
         IF(connectivity(iatom,current_atom) .NE. 0) THEN
            IF(connectivity(iatom,current_atom) .NE. avoid_atom) THEN
               atom(connectivity(iatom,current_atom))=1
               CALL depth_first_search ( connectivity(iatom,current_atom),&
                 current_atom,connectivity,atom)
            ENDIF
         ELSE
            RETURN
         ENDIF
      ENDDO

  END SUBROUTINE depth_first_search

! *****************************************************************************
!> \brief performs either a bond or angle change move for a given molecule
!> \param mc_par the mc parameters for the force env
!> \param force_env the force environment used in the move
!> \param bias_env the force environment used to bias the move, if any (it may
!>            be null if lbias=.false. in mc_par)
!> \param moves the structure that keeps track of how many moves have been
!>               accepted/rejected
!> \param move_updates the structure that keeps track of how many moves have
!>               been accepted/rejected since the last time the displacements
!>               were updated
!> \param start_atom the number of the molecule's first atom, assuming the rest
!>        of the atoms follow sequentially
!> \param molecule_type the type of the molecule we're moving
!> \param box_number the box the molecule is in
!> \param bias_energy the biased energy of the system before the move
!> \param move_type dictates what kind of conformational change we do
!> \param lreject set to .true. if there is an overlap
!> \param rng_stream the random number stream that we draw from
!> \param error the cp_error type
!>
!>      This subroutine is designed to run in parallel.
!> \author MJM
! *****************************************************************************
  SUBROUTINE mc_conformation_change ( mc_par,force_env,bias_env, moves,&
                        move_updates,start_atom,molecule_type,box_number,&
                        bias_energy,move_type,lreject,&
                        rng_stream,error)

    TYPE(mc_simpar_type), POINTER            :: mc_par
    TYPE(force_env_type), POINTER            :: force_env, bias_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: start_atom, molecule_type, &
                                                box_number
    REAL(KIND=dp), INTENT(INOUT)             :: bias_energy
    CHARACTER(LEN=*), INTENT(IN)             :: move_type
    LOGICAL, INTENT(OUT)                     :: lreject
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_conformation_change', &
      routineP = moduleN//':'//routineN

    CHARACTER(default_string_length)         :: name
    CHARACTER(default_string_length), &
      DIMENSION(:), POINTER                  :: names
    INTEGER :: atom_number, end_atom, end_mol, group, handle, imol_type, &
      imolecule, ipart, istat, jbox, molecule_number, nunits_mol, source, &
      start_mol
    INTEGER, DIMENSION(:), POINTER           :: mol_type, nunits
    INTEGER, DIMENSION(:, :), POINTER        :: nchains
    LOGICAL                                  :: ionode, lbias, loverlap
    REAL(KIND=dp) :: BETA, bias_energy_new, bias_energy_old, dis_length, &
      exp_max_val, exp_min_val, rand, value, w
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r_new, r_old
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind, &
                                                molecule_kind_test
    TYPE(particle_list_type), POINTER        :: particles

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

! nullify some pointers
      NULLIFY(particles,subsys,molecule_kinds_new,molecule_kind,&
           molecule_kind_test)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,lbias=lbias,mc_molecule_info=mc_molecule_info,&
         BETA=BETA,exp_max_val=exp_max_val,&
         exp_min_val=exp_min_val,group=group,source=source,ionode=ionode)
      CALL get_mc_molecule_info(mc_molecule_info,nchains=nchains,nunits=nunits,&
           mol_type=mol_type,names=names)

! do some allocation
      nunits_mol=nunits(molecule_type)
      ALLOCATE (r_old(1:3,1:nunits_mol),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r",3*nunits_mol*dp_size)
      ALLOCATE (r_new(1:3,1:nunits_mol),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r_new",3*nunits_mol*dp_size)

! find out some bounds for mol_type
      start_mol=1
      DO jbox=1,box_number-1
         start_mol=start_mol+SUM(nchains(:,jbox))
      ENDDO
      end_mol=start_mol+SUM(nchains(:,box_number))-1

! figure out which molecule number we are
      end_atom=start_atom+nunits_mol-1
      molecule_number=0
      atom_number=1
      DO imolecule=1,SUM(nchains(:,box_number))
         IF(atom_number == start_atom) THEN
            molecule_number=imolecule
            EXIT
         ENDIF
         atom_number=atom_number+nunits(mol_type(imolecule+start_mol-1))
      ENDDO
      IF(molecule_number == 0) CALL stop_program(routineN,moduleN,__LINE__,&
           'Cannot find the molecule number')

! are we biasing this move?
      IF(lbias) THEN

! grab the coordinates
         CALL force_env_get(bias_env,subsys=subsys,error=error)
! save the energy
         bias_energy_old=bias_energy

      ELSE

! grab the coordinates
         CALL force_env_get(force_env,subsys=subsys,error=error)
      ENDIF

! now find the molecule type associated with this guy
      CALL cp_subsys_get(subsys, &
           particles=particles, molecule_kinds_new=molecule_kinds_new,&
           error=error)
      DO imol_type=1,SIZE(molecule_kinds_new%els(:))
         molecule_kind_test => molecule_kinds_new%els(imol_type)
         CALL get_molecule_kind(molecule_kind_test,name=name)
         IF(TRIM(ADJUSTL(name)) == TRIM(ADJUSTL(names(molecule_type)))) THEN
            molecule_kind => molecule_kinds_new%els(imol_type)
            EXIT
         ENDIF
      ENDDO

! save the coordinates
      DO ipart=start_atom,end_atom
         r_old(1:3,ipart-start_atom+1)=particles%els(ipart)%r(1:3)
      ENDDO

      IF(.NOT. ASSOCIATED(molecule_kind)) CALL &
           stop_program(routineN,moduleN,__LINE__,&
           'Cannot find the molecule type')
! do the move
      IF (move_type == 'bond') THEN

! record the attempt
         moves%bond%attempts=moves%bond%attempts+1
         move_updates%bond%attempts=move_updates%bond%attempts+1
         moves%bias_bond%attempts=moves%bias_bond%attempts+1
         move_updates%bias_bond%attempts=move_updates%bias_bond%attempts+1
         IF ( .NOT. lbias ) THEN
            moves%bond%qsuccesses=moves%bond%qsuccesses+1
            move_updates%bond%qsuccesses=&
                 move_updates%bond%qsuccesses+1
            moves%bias_bond%qsuccesses=moves%bias_bond%qsuccesses+1
            move_updates%bias_bond%qsuccesses=&
                 move_updates%bias_bond%qsuccesses+1
         ENDIF

! do the move
         CALL change_bond_length(r_old,r_new,mc_par,molecule_type,&
              molecule_kind,dis_length,particles,rng_stream,error=error)

      ELSEIF( move_type == 'angle') THEN

! record the attempt
         moves%angle%attempts=moves%angle%attempts+1
         move_updates%angle%attempts=move_updates%angle%attempts+1
         moves%bias_angle%attempts=moves%bias_angle%attempts+1
         move_updates%bias_angle%attempts=move_updates%bias_angle%attempts+1
         IF ( .NOT. lbias ) THEN
            moves%angle%qsuccesses=moves%angle%qsuccesses+1
            move_updates%angle%qsuccesses=&
                 move_updates%angle%qsuccesses+1
            moves%bias_angle%qsuccesses=moves%bias_angle%qsuccesses+1
            move_updates%bias_angle%qsuccesses=&
                 move_updates%bias_angle%qsuccesses+1
         ENDIF

! do the move
         CALL change_bond_angle(r_old,r_new,mc_par,molecule_type,&
              molecule_kind,particles,rng_stream,error=error)
         dis_length=1.0E0_dp
      ELSE
! record the attempt
         moves%dihedral%attempts=moves%dihedral%attempts+1
         move_updates%dihedral%attempts=move_updates%dihedral%attempts+1
         moves%bias_dihedral%attempts=moves%bias_dihedral%attempts+1
         move_updates%bias_dihedral%attempts=move_updates%bias_dihedral%attempts+1
         IF ( .NOT. lbias ) THEN
            moves%dihedral%qsuccesses=moves%dihedral%qsuccesses+1
            move_updates%dihedral%qsuccesses=&
                 move_updates%dihedral%qsuccesses+1
            moves%bias_dihedral%qsuccesses=moves%bias_dihedral%qsuccesses+1
            move_updates%bias_dihedral%qsuccesses=&
                 move_updates%bias_dihedral%qsuccesses+1
         ENDIF

! do the move
         CALL change_dihedral(r_old,r_new,mc_par,molecule_type,&
              molecule_kind,particles,rng_stream,error=error)
         dis_length=1.0E0_dp

      ENDIF

! set the coordinates
      DO ipart=start_atom,end_atom
         particles%els(ipart)%r(1:3)=r_new(1:3,ipart-start_atom+1)
      ENDDO

! check for overlap
      lreject=.FALSE.
      IF(lbias) THEN
         CALL check_for_overlap(bias_env,nchains(:,box_number),&
              nunits(:),loverlap,mol_type(start_mol:end_mol),&
              molecule_number=molecule_number)
      ELSE
         CALL check_for_overlap(force_env,nchains(:,box_number),&
              nunits(:),loverlap,mol_type(start_mol:end_mol),&
              molecule_number=molecule_number)
         IF(loverlap) lreject=.TRUE.
      ENDIF

! if we're biasing classical, check for acceptance
      IF(lbias) THEN

! here's where we bias the moves

         IF(loverlap) THEN
            w=0.0E0_dp
         ELSE
            CALL force_env_calc_energy_force(bias_env,calc_force=.FALSE.,error=error)
            CALL force_env_get(bias_env,&
               potential_energy=bias_energy_new,error=error)
! accept or reject the move based on the Metropolis rule with a
! correction factor for the change in phase space...dis_length is
! made unitless in change_bond_length
            value=-BETA*(bias_energy_new-bias_energy_old)
            IF    (value .GT. exp_max_val) THEN
               w=10.0_dp
            ELSEIF(value .LT. exp_min_val) THEN
               w=0.0_dp
            ELSE
               w=EXP(value)*dis_length**2
            ENDIF

         ENDIF

         IF ( w .GE. 1.0E0_dp ) THEN
            w=1.0E0_dp
            rand=0.0E0_dp
         ELSE
            IF(ionode) THEN
               rand=next_random_number(rng_stream,error=error)
            ENDIF
            CALL mp_bcast(rand,source,group)
         ENDIF

         IF (rand .LT. w) THEN

! accept the move
            IF (move_type == 'bond') THEN
               moves%bond%qsuccesses=moves%bond%qsuccesses+1
               move_updates%bond%successes=&
                  move_updates%bond%successes+1
               moves%bias_bond%successes=moves%bias_bond%successes+1
               move_updates%bias_bond%successes=&
                  move_updates%bias_bond%successes+1
            ELSEIF(move_type == 'angle') THEN
               moves%angle%qsuccesses=moves%angle%qsuccesses+1
               move_updates%angle%successes=&
                  move_updates%angle%successes+1
               moves%bias_angle%successes=moves%bias_angle%successes+1
               move_updates%bias_angle%successes=&
                  move_updates%bias_angle%successes+1
            ELSE
               moves%dihedral%qsuccesses=moves%dihedral%qsuccesses+1
               move_updates%dihedral%successes=&
                  move_updates%dihedral%successes+1
               moves%bias_dihedral%successes=moves%bias_dihedral%successes+1
               move_updates%bias_dihedral%successes=&
                  move_updates%bias_dihedral%successes+1
            ENDIF

            bias_energy=bias_energy+bias_energy_new-&
                                 bias_energy_old

         ELSE

! reject the move
! restore the coordinates
            CALL force_env_get(bias_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys,particles=particles, error=error)
            DO ipart=start_atom,end_atom
               particles%els(ipart)%r(1:3)=r_old(1:3,ipart-start_atom+1)
            ENDDO
            CALL cp_subsys_set(subsys,particles=particles,error=error)

         ENDIF

      ENDIF

! deallocate some stuff
      DEALLOCATE(r_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r")
      DEALLOCATE(r_new,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r_new")

! end the timing
  CALL timestop(handle)

END SUBROUTINE mc_conformation_change

! *****************************************************************************
!> \brief translates the given molecule randomly in either the x,y, or z direction
!> \param mc_par the mc parameters for the force env
!> \param force_env the force environment used in the move
!> \param bias_env the force environment used to bias the move, if any (it may
!>            be null if lbias=.false. in mc_par)
!> \param moves the structure that keeps track of how many moves have been
!>               accepted/rejected
!> \param move_updates the structure that keeps track of how many moves have
!>               been accepted/rejected since the last time the displacements
!>               were updated
!> \param start_atom the number of the molecule's first atom, assuming the rest of
!>        the atoms follow sequentially
!> \param box_number the box the molecule is in
!> \param bias_energy the biased energy of the system before the move
!> \param molecule_type the type of molecule we're moving
!> \param lreject set to .true. if there is an overlap
!> \param rng_stream the random number stream that we draw from
!> \param error the cp_error type
!>
!>      This subroutine designed to be parallel.
!> \author MJM
! *****************************************************************************
 SUBROUTINE mc_molecule_translation(  mc_par,force_env, bias_env,moves,&
                      move_updates,start_atom,box_number,&
                      bias_energy,molecule_type,&
                        lreject,rng_stream,error)

    TYPE(mc_simpar_type), POINTER            :: mc_par
    TYPE(force_env_type), POINTER            :: force_env, bias_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: start_atom, box_number
    REAL(KIND=dp), INTENT(INOUT)             :: bias_energy
    INTEGER, INTENT(IN)                      :: molecule_type
    LOGICAL, INTENT(OUT)                     :: lreject
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_molecule_translation', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_number, end_atom, end_mol, group, handle, imolecule, &
      ipart, iparticle, istat, jbox, molecule_number, move_direction, &
      nunits_mol, source, start_mol
    INTEGER, DIMENSION(:), POINTER           :: mol_type, nunits, nunits_tot
    INTEGER, DIMENSION(:, :), POINTER        :: nchains
    LOGICAL                                  :: ionode, lbias, loverlap
    REAL(dp), DIMENSION(:), POINTER          :: rmtrans
    REAL(KIND=dp)                            :: BETA, bias_energy_new, &
                                                bias_energy_old, dis_mol, &
                                                exp_max_val, exp_min_val, &
                                                rand, value, w
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r_old
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    TYPE(particle_list_type), POINTER        :: particles

!   *** Local Counters ***
! begin the timing of the subroutine

      CALL timeset(routineN,handle)

! nullify some pointers
      NULLIFY(particles,subsys)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,lbias=lbias,&
         BETA=BETA,exp_max_val=exp_max_val,&
         exp_min_val=exp_min_val,rmtrans=rmtrans,ionode=ionode,source=source,&
         group=group,mc_molecule_info=mc_molecule_info)
      CALL get_mc_molecule_info(mc_molecule_info,nunits_tot=nunits_tot,&
           nchains=nchains,nunits=nunits,mol_type=mol_type)

! find out some bounds for mol_type
      start_mol=1
      DO jbox=1,box_number-1
         start_mol=start_mol+SUM(nchains(:,jbox))
      ENDDO
      end_mol=start_mol+SUM(nchains(:,box_number))-1

! do some allocation
      ALLOCATE (r_old(1:3,1:nunits_tot(box_number)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r",3*nunits_tot(box_number)*dp_size)

! find the index of the last atom of this molecule, and the molecule number
      nunits_mol=nunits(molecule_type)
      end_atom=start_atom+nunits_mol-1
      molecule_number=0
      atom_number=1
      DO imolecule=1,SUM(nchains(:,box_number))
         IF(atom_number == start_atom) THEN
            molecule_number=imolecule
            EXIT
         ENDIF
         atom_number=atom_number+nunits(mol_type(imolecule+start_mol-1))
      ENDDO
      IF(molecule_number == 0) CALL stop_program(routineN,moduleN,__LINE__,&
           'Cannot find the molecule number')

! are we biasing this move?
      IF(lbias) THEN

! grab the coordinates
         CALL force_env_get(bias_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys,particles=particles, error=error)

! save the coordinates
         DO ipart=1,nunits_tot(box_number)
            r_old(1:3,ipart)=particles%els(ipart)%r(1:3)
        ENDDO

! save the energy
         bias_energy_old=bias_energy

      ELSE

! grab the coordinates
         CALL force_env_get(force_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys,particles=particles, error=error)
      ENDIF

! record the attempt
      moves%trans%attempts=moves%trans%attempts+1
      move_updates%trans%attempts=move_updates%trans%attempts+1
      moves%bias_trans%attempts=moves%bias_trans%attempts+1
      move_updates%bias_trans%attempts=move_updates%bias_trans%attempts+1
      IF ( .NOT. lbias ) THEN
         moves%trans%qsuccesses=moves%trans%qsuccesses+1
         move_updates%trans%qsuccesses=move_updates%trans%qsuccesses+1
         moves%bias_trans%qsuccesses=moves%bias_trans%qsuccesses+1
         move_updates%bias_trans%qsuccesses=move_updates%bias_trans%qsuccesses+1
      ENDIF

! move one molecule in the system

! call a random number to figure out which direction we're moving
      IF(ionode) rand=next_random_number(rng_stream,error=error)
      CALL mp_bcast(rand,source,group)
      ! 1,2,3 with equal prob
      move_direction=INT(3*rand)+1

! call a random number to figure out how far we're moving
      IF(ionode) rand=next_random_number(rng_stream,error=error)
      CALL mp_bcast(rand,source,group)
      dis_mol=rmtrans(molecule_type)*(rand-0.5E0_dp)*2.0E0_dp

! do the move
      DO iparticle=start_atom,end_atom
         particles%els(iparticle)%r(move_direction)=&
             particles%els(iparticle)%r(move_direction)+dis_mol
      ENDDO
      CALL cp_subsys_set(subsys,particles=particles,error=error)

! figure out if there is any overlap...need the number of the molecule
      lreject=.FALSE.
      IF(lbias) THEN
         CALL check_for_overlap(bias_env,nchains(:,box_number),&
              nunits(:),loverlap,mol_type(start_mol:end_mol),&
              molecule_number=molecule_number)
      ELSE
         CALL check_for_overlap(force_env,nchains(:,box_number),&
              nunits(:),loverlap,mol_type(start_mol:end_mol),&
              molecule_number=molecule_number)
         IF(loverlap) lreject=.TRUE.
      ENDIF

! if we're biasing with a cheaper potential, check for acceptance
      IF(lbias) THEN

! here's where we bias the moves
         IF(loverlap) THEN
            w=0.0E0_dp
         ELSE
            CALL force_env_calc_energy_force(bias_env,calc_force=.FALSE.,error=error)
            CALL force_env_get(bias_env,&
               potential_energy=bias_energy_new,error=error)
! accept or reject the move based on the Metropolis rule
            value=-BETA*(bias_energy_new-bias_energy_old)
            IF    (value .GT. exp_max_val) THEN
               w=10.0_dp
            ELSEIF(value .LT. exp_min_val) THEN
               w=0.0_dp
            ELSE
               w=EXP(value)
            ENDIF

         ENDIF

         IF ( w .GE. 1.0E0_dp ) THEN
            w=1.0E0_dp
            rand=0.0E0_dp
         ELSE
            IF(ionode) rand=next_random_number(rng_stream,error=error)
            CALL mp_bcast(rand,source,group)
         ENDIF

         IF (rand .LT. w) THEN

! accept the move
            moves%bias_trans%successes=moves%bias_trans%successes+1
            move_updates%bias_trans%successes=move_updates%bias_trans%successes+1
            moves%trans%qsuccesses=moves%trans%qsuccesses+1
            move_updates%trans%successes=&
                      move_updates%trans%successes+1
            moves%qtrans_dis=moves%qtrans_dis+ABS(dis_mol)
            bias_energy=bias_energy+bias_energy_new-&
                                 bias_energy_old

         ELSE

! reject the move
! restore the coordinates
            CALL force_env_get(bias_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys,particles=particles, error=error)
            DO ipart=1,nunits_tot(box_number)
               particles%els(ipart)%r(1:3)=r_old(1:3,ipart)
            ENDDO
            CALL cp_subsys_set(subsys,particles=particles,error=error)

         ENDIF

      ENDIF

! deallocate some stuff
      DEALLOCATE(r_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r")

! end the timing
  CALL timestop(handle)

  END SUBROUTINE mc_molecule_translation

! *****************************************************************************
!> \brief rotates the given molecule randomly around the x,y, or z axis...
!>      only works for water at the moment
!> \param mc_par the mc parameters for the force env
!> \param force_env the force environment used in the move
!> \param bias_env the force environment used to bias the move, if any (it may
!>            be null if lbias=.false. in mc_par)
!> \param moves the structure that keeps track of how many moves have been
!>               accepted/rejected
!> \param move_updates the structure that keeps track of how many moves have
!>               been accepted/rejected since the last time the displacements
!>               were updated
!> \param box_number the box the molecule is in
!> \param start_atom the number of the molecule's first atom, assuming the rest of
!>        the atoms follow sequentially
!> \param molecule_type the type of molecule we're moving
!> \param bias_energy the biased energy of the system before the move
!> \param lreject set to .true. if there is an overlap
!> \param rng_stream the random number stream that we draw from
!> \param error the cp_error type
!>
!>      This subroutine has been designed for parallel use.
!> \author MJM
! *****************************************************************************
 SUBROUTINE mc_molecule_rotation ( mc_par,force_env, bias_env,moves,&
                      move_updates,box_number,&
                      start_atom,molecule_type,bias_energy,lreject,&
                      rng_stream,error)

    TYPE(mc_simpar_type), POINTER            :: mc_par
    TYPE(force_env_type), POINTER            :: force_env, bias_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: box_number, start_atom, &
                                                molecule_type
    REAL(KIND=dp), INTENT(INOUT)             :: bias_energy
    LOGICAL, INTENT(OUT)                     :: lreject
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_molecule_rotation', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_number, dir, end_atom, end_mol, group, handle, ii, &
      imolecule, ipart, istat, iunit, jbox, molecule_number, nunits_mol, &
      source, start_mol
    INTEGER, DIMENSION(:), POINTER           :: mol_type, nunits, nunits_tot
    INTEGER, DIMENSION(:, :), POINTER        :: nchains
    LOGICAL                                  :: ionode, lbias, loverlap, lx, &
                                                ly
    REAL(dp), DIMENSION(:), POINTER          :: rmrot
    REAL(dp), DIMENSION(:, :), POINTER       :: mass
    REAL(KIND=dp) :: BETA, bias_energy_new, bias_energy_old, cosdg, dgamma, &
      exp_max_val, exp_min_val, masstot, nxcm, nycm, nzcm, rand, rx, rxnew, &
      ry, rynew, rz, rznew, sindg, value, w
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r_old
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    TYPE(particle_list_type), POINTER        :: particles

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

      NULLIFY(rmrot,subsys,particles)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,lbias=lbias,&
         BETA=BETA,exp_max_val=exp_max_val,&
         exp_min_val=exp_min_val,rmrot=rmrot,mc_molecule_info=mc_molecule_info,&
         ionode=ionode,group=group,source=source)
      CALL get_mc_molecule_info(mc_molecule_info,nunits=nunits,&
           nunits_tot=nunits_tot,nchains=nchains,mass=mass,&
           mol_type=mol_type)

! figure out some bounds for mol_type
      start_mol=1
      DO jbox=1,box_number-1
         start_mol=start_mol+SUM(nchains(:,jbox))
      ENDDO
      end_mol=start_mol+SUM(nchains(:,box_number))-1

      nunits_mol=nunits(molecule_type)

! nullify some pointers
      NULLIFY(particles,subsys)

! do some allocation
      ALLOCATE (r_old(1:3,1:nunits_tot(box_number)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r",3*nunits_tot(box_number)*dp_size)

! initialize some stuff
      lx =.FALSE.
      ly =.FALSE.

! determine what the final atom in the molecule is numbered, and which
! molecule number this is
      end_atom=start_atom+nunits_mol-1
      molecule_number=0
      atom_number=1
      DO imolecule=1,SUM(nchains(:,box_number))
         IF(atom_number == start_atom) THEN
            molecule_number=imolecule
            EXIT
         ENDIF
         atom_number=atom_number+nunits(mol_type(imolecule+start_mol-1))
      ENDDO
      IF(molecule_number == 0) CALL stop_program(routineN,moduleN,__LINE__,&
               'Cannot find the molecule number')

! are we biasing this move?
      IF(lbias) THEN

! grab the coordinates
         CALL force_env_get(bias_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys,particles=particles, error=error)

! save the coordinates
         DO ipart=1,nunits_tot(box_number)
            r_old(1:3,ipart)=particles%els(ipart)%r(1:3)
         ENDDO

! save the energy
         bias_energy_old=bias_energy

      ELSE

! grab the coordinates
         CALL force_env_get(force_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys,particles=particles, error=error)
      ENDIF

! grab the masses
      masstot=SUM(mass(1:nunits(molecule_type),molecule_type))

! record the attempt
      moves%bias_rot%attempts=moves%bias_rot%attempts+1
      move_updates%bias_rot%attempts=move_updates%bias_rot%attempts+1
      moves%rot%attempts=moves%rot%attempts+1
      move_updates%rot%attempts=move_updates%rot%attempts+1
      IF ( .NOT. lbias ) THEN
         moves%rot%qsuccesses=moves%rot%qsuccesses+1
         move_updates%rot%qsuccesses=move_updates%rot%qsuccesses+1
         moves%bias_rot%qsuccesses=moves%bias_rot%qsuccesses+1
         move_updates%bias_rot%qsuccesses=move_updates%bias_rot%qsuccesses+1
      ENDIF

! rotate one molecule in the system

! call a random number to figure out which direction we're moving
      IF(ionode) rand=next_random_number(rng_stream,error=error)
!      CALL RANDOM_NUMBER(rand)
      CALL mp_bcast(rand,source,group)
      ! 1,2,3 with equal prob
      dir=INT(3*rand)+1

      IF (dir .EQ. 1) THEN
         lx = .TRUE.
      ELSEIF (dir .EQ. 2) THEN
         ly = .TRUE.
      ENDIF

! Determine new center of mass for chain i by finding the sum
! of m*r for each unit, then dividing by the total mass of the chain
      nxcm = 0.0E0_dp
      nycm = 0.0E0_dp
      nzcm = 0.0E0_dp
      DO ii = 1, nunits_mol
         nxcm = nxcm + particles%els(start_atom-1+ii)%r(1)* mass(ii,molecule_type)
         nycm = nycm + particles%els(start_atom-1+ii)%r(2)* mass(ii,molecule_type)
         nzcm = nzcm + particles%els(start_atom-1+ii)%r(3)* mass(ii,molecule_type)
      ENDDO
      nxcm = nxcm / masstot
      nycm = nycm / masstot
      nzcm = nzcm / masstot

! call a random number to figure out how far we're moving
      IF(ionode) rand=next_random_number(rng_stream,error=error)
      CALL mp_bcast(rand,source,group)
      dgamma=rmrot(molecule_type)*(rand-0.5E0_dp)*2.0E0_dp

! *** set up the rotation matrix ***

      cosdg = COS( dgamma )
      sindg = SIN( dgamma )

      IF (lx) THEN

! ***    ROTATE UNITS OF I AROUND X-AXIS ***

         DO  iunit = start_atom,end_atom
            ry = particles%els(iunit)%r(2) - nycm
            rz = particles%els(iunit)%r(3) - nzcm
            rynew = cosdg * ry - sindg * rz
            rznew = cosdg * rz + sindg * ry

            particles%els(iunit)%r(2) = rynew + nycm
            particles%els(iunit)%r(3) = rznew + nzcm

         ENDDO
      ELSEIF (ly) THEN

! ***    ROTATE UNITS OF I AROUND y-AXIS ***

         DO  iunit = start_atom,end_atom
            rx = particles%els(iunit)%r(1) - nxcm
            rz = particles%els(iunit)%r(3) - nzcm
            rxnew = cosdg * rx + sindg * rz
            rznew = cosdg * rz - sindg * rx

            particles%els(iunit)%r(1) = rxnew + nxcm
            particles%els(iunit)%r(3) = rznew + nzcm

         ENDDO

      ELSE

! ***    ROTATE UNITS OF I AROUND z-AXIS ***

         DO  iunit = start_atom,end_atom
            rx = particles%els(iunit)%r(1) - nxcm
            ry = particles%els(iunit)%r(2) - nycm

            rxnew = cosdg * rx - sindg * ry
            rynew = cosdg * ry + sindg * rx

            particles%els(iunit)%r(1) = rxnew + nxcm
            particles%els(iunit)%r(2) = rynew + nycm

         ENDDO

      ENDIF
      CALL cp_subsys_set(subsys,particles=particles,error=error)

! check for overlap
      lreject=.FALSE.
      IF(lbias) THEN
         CALL check_for_overlap(bias_env,nchains(:,box_number),&
              nunits(:),loverlap,mol_type(start_mol:end_mol),&
              molecule_number=molecule_number)
      ELSE
         CALL check_for_overlap(force_env,nchains(:,box_number),&
              nunits(:),loverlap,mol_type(start_mol:end_mol),&
              molecule_number=molecule_number)
         IF(loverlap) lreject=.TRUE.
      ENDIF

! if we're biasing classical, check for acceptance
      IF(lbias) THEN

! here's where we bias the moves

         IF(loverlap) THEN
            w=0.0E0_dp
         ELSE
            CALL force_env_calc_energy_force(bias_env,calc_force=.FALSE.,error=error)
            CALL force_env_get(bias_env,&
            potential_energy=bias_energy_new,error=error)
! accept or reject the move based on the Metropolis rule
            value=-BETA*(bias_energy_new-bias_energy_old)
            IF    (value .GT. exp_max_val) THEN
               w=10.0_dp
            ELSEIF(value .LT. exp_min_val) THEN
               w=0.0_dp
            ELSE
               w=EXP(value)
            ENDIF

         ENDIF

         IF ( w .GE. 1.0E0_dp ) THEN
            w=1.0E0_dp
            rand=0.0E0_dp
         ELSE
            IF(ionode) rand=next_random_number(rng_stream,error=error)
            CALL mp_bcast(rand,source,group)
         ENDIF

         IF (rand .LT. w) THEN

! accept the move
            moves%bias_rot%successes=moves%bias_rot%successes+1
            move_updates%bias_rot%successes=move_updates%bias_rot%successes+1
            moves%rot%qsuccesses=moves%rot%qsuccesses+1
            move_updates%rot%successes=move_updates%rot%successes+1
            bias_energy=bias_energy+bias_energy_new-&
                                 bias_energy_old

         ELSE

! reject the move
! restore the coordinates
            CALL force_env_get(bias_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys,particles=particles, error=error)
            DO ipart=1,nunits_tot(box_number)
               particles%els(ipart)%r(1:3)=r_old(1:3,ipart)
            ENDDO
            CALL cp_subsys_set(subsys,particles=particles,error=error)

         ENDIF

      ENDIF

! deallocate some stuff
      DEALLOCATE(r_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r_old")

! end the timing
  CALL timestop(handle)

  END SUBROUTINE mc_molecule_rotation

! *****************************************************************************
!> \brief performs a Monte Carlo move that alters the volume of the simulation box
!> \param mc_par the mc parameters for the force env
!> \param force_env the force environment whose cell we're changing
!> \param moves the structure that keeps track of how many moves have been
!>               accepted/rejected
!> \param move_updates the structure that keeps track of how many moves have
!>               been accepted/rejected since the last time the displacements
!>               were updated
!> \param old_energy the energy of the last accepted move involving an
!>                    unbiased calculation
!> \param box_number the box we're changing the volume of
!> \param energy_check the running total of how much the energy has changed
!>                      since the initial configuration
!> \param r_old the coordinates of the last accepted move involving an
!>               unbiased calculation
!> \param iw the unit number that writes to the screen
!> \param discrete_array tells use which volumes we can do for the discrete
!>            case
!> \param rng_stream the random number stream that we draw from
!> \param error the cp error type
!>
!> \author MJM
!> \note     Designed for parallel use.
! *****************************************************************************
  SUBROUTINE mc_volume_move ( mc_par,force_env, moves,move_updates,&
                        old_energy,box_number,&
                        energy_check,r_old,iw,discrete_array,rng_stream,error)

    TYPE(mc_simpar_type), POINTER            :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    REAL(KIND=dp), INTENT(INOUT)             :: old_energy
    INTEGER, INTENT(IN)                      :: box_number
    REAL(KIND=dp), INTENT(INOUT)             :: energy_check
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: r_old
    INTEGER, INTENT(IN)                      :: iw
    INTEGER, DIMENSION(1:3, 1:2), &
      INTENT(INOUT)                          :: discrete_array
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_volume_move', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=200)                       :: fft_lib
    CHARACTER(LEN=40)                        :: dat_file
    INTEGER :: cl, end_atom, end_mol, group, handle, iatom, idim, imolecule, &
      iside, iside_change, istat, iunit, jbox, nunits_mol, output_unit, &
      print_level, source, start_atom, start_mol
    INTEGER, DIMENSION(:), POINTER           :: mol_type, nunits, nunits_tot
    INTEGER, DIMENSION(:, :), POINTER        :: nchains
    LOGICAL                                  :: ionode, ldiscrete, lincrease, &
                                                loverlap, ltoo_small
    REAL(dp), DIMENSION(:, :), POINTER       :: mass
    REAL(KIND=dp) :: BETA, discrete_step, energy_term, exp_max_val, &
      exp_min_val, new_energy, pressure, pressure_term, rand, rcut, rmvolume, &
      temp_var, value, vol_dis, volume_term, w
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r
    REAL(KIND=dp), DIMENSION(1:3)            :: abc, center_of_mass, &
                                                center_of_mass_new, diff, &
                                                new_cell_length, &
                                                old_cell_length
    REAL(KIND=dp), DIMENSION(1:3, 1:3)       :: hmat_test
    TYPE(cell_type), POINTER                 :: cell, cell_old, cell_test
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: oldsys
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    TYPE(particle_list_type), POINTER        :: particles_old

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,ionode=ionode,&
         BETA=BETA,exp_max_val=exp_max_val,&
         exp_min_val=exp_min_val,source=source,group=group,&
         dat_file=dat_file,rmvolume=rmvolume,pressure=pressure,cl=cl,&
         fft_lib=fft_lib,discrete_step=discrete_step,&
         ldiscrete=ldiscrete,mc_molecule_info=mc_molecule_info)
      CALL get_mc_molecule_info(mc_molecule_info,nchains=nchains,&
           nunits=nunits,nunits_tot=nunits_tot,mol_type=mol_type,&
           mass=mass)
! figure out some bounds for mol_type
      start_mol=1
      DO jbox=1,box_number-1
         start_mol=start_mol+SUM(nchains(:,jbox))
      ENDDO
      end_mol=start_mol+SUM(nchains(:,box_number))-1

      print_level = 1 ! hack, printlevel is for print_keys

! nullify some pointers
      NULLIFY(particles_old,cell_old,oldsys,cell_test,cell)

! do some allocation
      ALLOCATE (r(1:3,1:nunits_tot(box_number)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "r",3*nunits_tot(box_number)*dp_size)

! record the attempt
      moves%volume%attempts=moves%volume%attempts+1
      move_updates%volume%attempts=move_updates%volume%attempts+1

! now let's grab the cell length and particle positions
      CALL force_env_get(force_env,subsys=oldsys,cell=cell,&
         error=error)
      CALL get_cell(cell,abc=abc)
      CALL cell_create(cell_old,error=error)
      CALL cell_clone(cell,cell_old,error=error)
      CALL cp_subsys_get(oldsys,particles=particles_old, &
         error=error)

! find the old cell length
      old_cell_length(1)=abc(1)
      old_cell_length(2)=abc(2)
      old_cell_length(3)=abc(3)

! save the old coordiantes
      DO iatom=1,nunits_tot(box_number)
         r(1:3,iatom)=particles_old%els(iatom)%r(1:3)
      ENDDO

! now do the move

! call a random number to figure out how far we're moving
      IF (ionode) rand=next_random_number(rng_stream,error=error)
      CALL mp_bcast(rand,source,group)

! find the test cell lenghts for the discrete volume move
      IF(ldiscrete) THEN
         IF(rand .LT. 0.5_dp) THEN
            lincrease=.TRUE.
         ELSE
            lincrease=.FALSE.
         ENDIF

         new_cell_length(1:3)=old_cell_length(1:3)

! if we're increasing the volume, we need to find a side we can increase
         IF(lincrease) THEN
            DO
               IF(ionode) rand=next_random_number(rng_stream,error=error)
               CALL mp_bcast(rand,source,group)
               iside_change=CEILING(3.0_dp*rand)
               IF(discrete_array(iside_change,1) .EQ. 1) THEN
                  new_cell_length(iside_change)=&
                  new_cell_length(iside_change)+discrete_step
                  EXIT
               ENDIF
            ENDDO
         ELSE
            DO
               IF(ionode) rand=next_random_number(rng_stream,error=error)
               CALL mp_bcast(rand,source,group)
               iside_change=CEILING(3.0_dp*rand)
               IF(discrete_array(iside_change,2) .EQ. 1) THEN
                  new_cell_length(iside_change)=&
                     new_cell_length(iside_change)-discrete_step
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         vol_dis=(new_cell_length(1)*new_cell_length(2)*new_cell_length(3))&
         -old_cell_length(1)*old_cell_length(2)*old_cell_length(3)
      ELSE
! now for the not discrete volume move
!!!!!!!!!!!!!!!! for E_V curves
         vol_dis=rmvolume*(rand-0.5E0_dp)*2.0E0_dp
!         WRITE(output_unit,*) '************************ be sure to change back!',&
!                 old_cell_length(1),14.64_dp/angstrom
!         vol_dis=-56.423592_dp/angstrom**3
!         IF(old_cell_length(1) .LE. 14.64_dp/angstrom) THEN
!            vol_dis=0.0_dp
!            WRITE(output_unit,*) 'Found the correct box length!'
!         ENDIF

         temp_var=vol_dis+&
              old_cell_length(1)*old_cell_length(2)*&
              old_cell_length(3)

         IF(temp_var .LE. 0.0E0_dp) THEN
            loverlap=.TRUE.  ! cannot have a negative volume
         ELSE
            new_cell_length(1)=(temp_var)**(1.0E0_dp/3.0E0_dp)
            new_cell_length(2)=new_cell_length(1)
            new_cell_length(3)=new_cell_length(1)
            loverlap=.FALSE.
         ENDIF
      ENDIF
      CALL mp_bcast(loverlap,source,group)

      IF(loverlap) THEN
! deallocate some stuff
         DEALLOCATE(r,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "r")
         logger=>cp_error_get_logger(error)
         output_unit=cp_logger_get_default_io_unit(logger)
         IF(output_unit>0) WRITE(output_unit,*) &
            "Volume move rejected because we tried to make too small of box.",vol_dis
!     end the timing
         CALL timestop(handle)
         RETURN
      ENDIF

! now we need to make the new cell
      hmat_test(:,:)=0.0e0_dp
      hmat_test(1,1)=new_cell_length(1)
      hmat_test(2,2)=new_cell_length(2)
      hmat_test(3,3)=new_cell_length(3)
      CALL cell_create(cell_test,hmat=hmat_test(:,:),periodic=cell%perd,error=error)
      CALL cp_subsys_set(oldsys, cell=cell_test, error=error)

! now we need to scale the coordinates of all the molecules by the
! center of mass, using the minimum image (not all molecules are in
! the central box)

! now we need to scale the coordinates of all the molecules by the
! center of mass
      end_atom=0
      DO imolecule=1,SUM(nchains(:,box_number))
         nunits_mol=nunits(mol_type(imolecule+start_mol-1))
         start_atom=end_atom+1
         end_atom=start_atom+nunits_mol-1
! now find the center of mass
         CALL get_center_of_mass(r(:,start_atom:end_atom),nunits_mol,&
            center_of_mass(:),mass(:,mol_type(imolecule+start_mol-1)))

! scale the center of mass and determine the vector that points from the
!    old COM to the new one
         DO iside=1,3
            center_of_mass_new(iside)=center_of_mass(iside)*&
               new_cell_length(iside)/old_cell_length(iside)
         ENDDO

         DO idim=1,3
            diff(idim)=center_of_mass_new(idim)-center_of_mass(idim)
! now change the particle positions
            DO iunit=start_atom,end_atom
               particles_old%els(iunit)%r(idim)=&
                  particles_old%els(iunit)%r(idim)+diff(idim)
            ENDDO
         ENDDO
      ENDDO

! check for overlap
      CALL check_for_overlap(force_env,nchains(:,box_number),&
           nunits(:),loverlap,mol_type(start_mol:end_mol),&
           cell_length=new_cell_length)

! figure out if we have overlap problems
      CALL mp_bcast(loverlap,source,group)
      IF(loverlap) THEN
! deallocate some stuff
         DEALLOCATE(r,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "r")

         logger=>cp_error_get_logger(error)
         output_unit=cp_logger_get_default_io_unit(logger)
         IF(output_unit>0) WRITE(output_unit,*) &
            "Volume move rejected due to overlap.",vol_dis
!     end the timing
         CALL timestop(handle)
! reset the cell and particle positions
         CALL cp_subsys_set(oldsys, cell=cell_old, error=error)
         DO iatom=1,nunits_tot(box_number)
            particles_old%els(iatom)%r(1:3)=r_old(1:3,iatom)
         ENDDO
         RETURN
      ENDIF

! stop if we're trying to change a box to a boxlength smaller than rcut
      IF(ionode) THEN
         ltoo_small=.FALSE.
         IF(force_env%in_use .EQ. use_fist_force) THEN
            CALL get_mc_par(mc_par,rcut=rcut)
            IF(new_cell_length(1) .LT. 2.0_dp*rcut) ltoo_small=.TRUE.
            IF(new_cell_length(2) .LT. 2.0_dp*rcut) ltoo_small=.TRUE.
            IF(new_cell_length(3) .LT. 2.0_dp*rcut) ltoo_small=.TRUE.

            IF(ltoo_small) THEN
               WRITE(iw,*) 'new_cell_lengths ',&
               new_cell_length(1:3)/angstrom
               WRITE(iw,*) 'rcut ',rcut/angstrom
            ENDIF
         ENDIF
      ENDIF
      CALL mp_bcast(ltoo_small,source,group)
      IF(ltoo_small) &
         CALL stop_program(routineN,moduleN,__LINE__,&
              "Attempted a volume move where box size got too small.")

! now compute the energy
      CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,error=error)
      CALL force_env_get(force_env,&
         potential_energy=new_energy,error=error)

! accept or reject the move
! to prevent overflows
      energy_term=new_energy-old_energy
      volume_term=-REAL(SUM(nchains(:,box_number)),dp)/BETA*&
      LOG(new_cell_length(1)*new_cell_length(2)*new_cell_length(3)/&
         (old_cell_length(1)*old_cell_length(2)*old_cell_length(3)))
      pressure_term=pressure*vol_dis

      value=-BETA*(energy_term+volume_term+pressure_term)
      IF    (value .GT. exp_max_val) THEN
         w=10.0_dp
      ELSEIF(value .LT. exp_min_val) THEN
         w=0.0_dp
      ELSE
         w=EXP(value)
      ENDIF

!!!!!!!!!!!!!!!! for E_V curves
!         w=1.0E0_dp
!         w=0.0E0_dp

      IF ( w .GE. 1.0E0_dp ) THEN
         w=1.0E0_dp
         rand=0.0E0_dp
      ELSE
         IF(ionode) rand=next_random_number(rng_stream,error=error)
         CALL mp_bcast(rand,source,group)
      ENDIF

      IF (rand .LT. w) THEN

! accept the move
         moves%volume%successes=moves%volume%successes+1
         move_updates%volume%successes=move_updates%volume%successes+1

! update energies
         energy_check=energy_check+(new_energy-old_energy)
         old_energy=new_energy

         DO iatom=1,nunits_tot(box_number)
            r_old(1:3,iatom)=particles_old%els(iatom)%r(1:3)
         ENDDO

! update discrete_array if we're doing a discrete volume move
         IF(ldiscrete) THEN
            CALL create_discrete_array(new_cell_length(:),&
               discrete_array(:,:),discrete_step)
         ENDIF

      ELSE

! reset the cell and particle positions
         CALL cp_subsys_set(oldsys, cell=cell_old, error=error)
         DO iatom=1,nunits_tot(box_number)
            particles_old%els(iatom)%r(1:3)=r_old(1:3,iatom)
         ENDDO

      ENDIF

! deallocate some stuff
      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r")
      CALL cell_release(cell_test,error=error)
      CALL cell_release(cell_old,error=error)

! end the timing
      CALL timestop(handle)

  END SUBROUTINE mc_volume_move

! *****************************************************************************
!> \brief alters the length of a random bond for the given molecule, using
!>      a mass weighted scheme so the lightest atoms move the most
!> \param r_old the initial coordinates of all molecules in the system
!> \param r_new the new coordinates of all molecules in the system
!> \param mc_par the mc parameters for the force env
!> \param molecule_type the molecule type that we're moving
!> \param molecule_kind the structure containing the molecule information
!> \param dis_length the ratio of the new bond length to the old bond length,
!>                    used in the acceptance rule
!> \param particles the particle_list_type for all particles in the force_env..
!>             used to grab the mass of each atom
!> \param rng_stream the random number stream that we draw from
!>
!>    This subroutine is written to be parallel.
!> \param error ...
!> \author MJM
! *****************************************************************************
  SUBROUTINE change_bond_length ( r_old,r_new,mc_par,molecule_type,molecule_kind,&
      dis_length,particles,rng_stream,error)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_old
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: r_new
    TYPE(mc_simpar_type), POINTER            :: mc_par
    INTEGER, INTENT(IN)                      :: molecule_type
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    REAL(KIND=dp), INTENT(OUT)               :: dis_length
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'change_bond_length', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bond_number, group, handle, &
                                                i, iatom, ibond, ipart, &
                                                istat, natom, nbond, source
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_a, atom_b, counter
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: connection, connectivity
    INTEGER, DIMENSION(:), POINTER           :: nunits
    LOGICAL                                  :: ionode
    REAL(dp), DIMENSION(:), POINTER          :: rmbond
    REAL(KIND=dp)                            :: atom_mass, mass_a, mass_b, &
                                                new_length_a, new_length_b, &
                                                old_length, rand
    REAL(KIND=dp), DIMENSION(1:3)            :: bond_a, bond_b
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

      NULLIFY(rmbond,mc_molecule_info)

! get some stuff from mc_par
      CALL get_mc_par(mc_par,mc_molecule_info=mc_molecule_info,source=source,&
         group=group,rmbond=rmbond,ionode=ionode)
      CALL get_mc_molecule_info(mc_molecule_info,nunits=nunits)

! copy the incoming coordinates so we can change them
      DO ipart=1,nunits(molecule_type)
         r_new(1:3,ipart)=r_old(1:3,ipart)
      ENDDO

! pick which bond in the molecule at random
      IF(ionode) THEN
         rand=next_random_number(rng_stream,error=error)
      ENDIF
      CALL mp_bcast(rand,source,group)
      CALL get_molecule_kind(molecule_kind,natom=natom,nbond=nbond,&
         bond_list=bond_list)
      bond_number=CEILING(rand*REAL(nbond,dp))

      ALLOCATE(connection(1:natom,1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "connection",2*natom*dp_size)
! assume at most six bonds per atom
      ALLOCATE(connectivity(1:6,1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "connectivity",6*natom*dp_size)
      ALLOCATE(counter(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "counter",natom)
      ALLOCATE(atom_a(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_a",natom)
      ALLOCATE(atom_b(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_b",natom)
      connection(:,:)=0
      connectivity(:,:)=0
      counter(:)=0
      atom_a(:)=0
      atom_b(:)=0

! now we need to find a list of atoms that each atom in this bond is connected
! to
      DO iatom=1,natom
         DO ibond=1,nbond
            IF(bond_list(ibond)%a == iatom) THEN
               counter(iatom)=counter(iatom)+1
               connectivity(counter(iatom),iatom)=bond_list(ibond)%b
            ELSEIF(bond_list(ibond)%b == iatom)THEN
               counter(iatom)=counter(iatom)+1
               connectivity(counter(iatom),iatom)=bond_list(ibond)%a
            ENDIF
         ENDDO
      ENDDO

! now I need to do a depth first search to figure out which atoms are on atom a's
! side and which are on atom b's
      atom_a(:)=0
      atom_a(bond_list(bond_number)%a)=1
      CALL depth_first_search(bond_list(bond_number)%a,bond_list(bond_number)%b,&
           connectivity(:,:),atom_a(:))
      atom_b(:)=0
      atom_b(bond_list(bond_number)%b)=1
      CALL depth_first_search(bond_list(bond_number)%b,bond_list(bond_number)%a,&
           connectivity(:,:),atom_b(:))

! now figure out the masses of the various sides, so we can weight how far we move each
! group of atoms
      mass_a=0.0_dp
      mass_b=0.0_dp
      DO iatom=1,natom
         CALL get_atomic_kind(particles%els(iatom)%atomic_kind,&
            mass=atom_mass)
         IF(atom_a(iatom) == 1) THEN
            mass_a=mass_a+atom_mass
         ELSE
            mass_b=mass_b+atom_mass
         ENDIF
      ENDDO

! choose a displacement
      IF(ionode) rand=next_random_number(rng_stream,error=error)
      CALL mp_bcast(rand,source,group)

      dis_length=rmbond(molecule_type)*2.0E0_dp*(rand-0.5E0_dp)

! find the bond vector that atom a will be moving
      DO i=1,3
         bond_a(i)=r_new(i,bond_list(bond_number)%a)-&
            r_new(i,bond_list(bond_number)%b)
         bond_b(i)=-bond_a(i)
      ENDDO

! notice we weight by the opposite masses...therefore lighter segments
! will move further
      old_length=SQRT(DOT_PRODUCT(bond_a,bond_a))
      new_length_a=dis_length*mass_b/(mass_a+mass_b)
      new_length_b=dis_length*mass_a/(mass_a+mass_b)

      DO i=1,3
         bond_a(i)=bond_a(i)/old_length*new_length_a
         bond_b(i)=bond_b(i)/old_length*new_length_b
      ENDDO

      DO iatom=1,natom
         IF(atom_a(iatom) == 1) THEN
            r_new(1,iatom)=r_new(1,iatom)+bond_a(1)
            r_new(2,iatom)=r_new(2,iatom)+bond_a(2)
            r_new(3,iatom)=r_new(3,iatom)+bond_a(3)
         ELSE
            r_new(1,iatom)=r_new(1,iatom)+bond_b(1)
            r_new(2,iatom)=r_new(2,iatom)+bond_b(2)
            r_new(3,iatom)=r_new(3,iatom)+bond_b(3)
         ENDIF
      ENDDO

! correct the value of dis_length for the acceptance rule
      dis_length=(old_length+dis_length)/old_length

      DEALLOCATE(connection,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "connection")
      DEALLOCATE(connectivity,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "connectivity")
      DEALLOCATE(counter,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "counter")
      DEALLOCATE(atom_a,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_a")
      DEALLOCATE(atom_b,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_b")
! end the timing
      CALL timestop(handle)

  END SUBROUTINE change_bond_length

! *****************************************************************************
!> \brief Alters the magnitude of a random angle in a molecule centered on atom C
!>      (connected to atoms A and B).  Atoms A and B are moved amounts related
!>      to their masses (and masses of all connecting atoms), so that heavier
!>      segments are moved less.
!> \param r_old the initial coordinates of all molecules in the system
!> \param r_new the new coordinates of all molecules in the system
!> \param mc_par the mc parameters for the force env
!> \param molecule_type the type of molecule we're playing with
!> \param molecule_kind the structure containing the molecule information
!> \param particles the particle_list_type for all particles in the force_env...
!>             used to grab the mass of each atom
!> \param rng_stream the random number stream that we draw from
!> \param error ...
!> \author MJM
! *****************************************************************************
  SUBROUTINE change_bond_angle ( r_old,r_new,mc_par,molecule_type,molecule_kind,&
      particles,rng_stream,error)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_old
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: r_new
    TYPE(mc_simpar_type), POINTER            :: mc_par
    INTEGER, INTENT(IN)                      :: molecule_type
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'change_bond_angle', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bend_number, group, handle, &
                                                i, iatom, ibond, ipart, &
                                                istat, natom, nbend, nbond, &
                                                source
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_a, atom_c, counter
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: connection, connectivity
    INTEGER, DIMENSION(:), POINTER           :: nunits
    LOGICAL                                  :: ionode
    REAL(dp), DIMENSION(:), POINTER          :: rmangle
    REAL(KIND=dp) :: atom_mass, bis_length, dis_angle, dis_angle_a, &
      dis_angle_c, mass_a, mass_c, new_angle_a, new_angle_c, old_angle, &
      old_length_a, old_length_c, rand, temp_length
    REAL(KIND=dp), DIMENSION(1:3)            :: bisector, bond_a, bond_c, &
                                                cross_prod, cross_prod_plane, &
                                                temp
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

      NULLIFY(bend_list,bond_list,rmangle,mc_molecule_info)

! get some stuff from mc_par
      CALL get_mc_par(mc_par,rmangle=rmangle,source=source,&
         group=group,ionode=ionode,mc_molecule_info=mc_molecule_info)
      CALL get_mc_molecule_info(mc_molecule_info,nunits=nunits)

! copy the incoming coordinates so we can change them
      DO ipart=1,nunits(molecule_type)
         r_new(1:3,ipart)=r_old(1:3,ipart)
      ENDDO

! pick which bond in the molecule at random
      IF(ionode) THEN
         rand=next_random_number(rng_stream,error=error)
      ENDIF
      CALL mp_bcast(rand,source,group)
      CALL get_molecule_kind(molecule_kind,natom=natom,nbend=nbend,&
         bend_list=bend_list,bond_list=bond_list,nbond=nbond)
      bend_number=CEILING(rand*REAL(nbend,dp))

      ALLOCATE(connection(1:natom,1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "connection",2*natom*dp_size)
! assume at most six bonds per atom
      ALLOCATE(connectivity(1:6,1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "connectivity",6*natom*dp_size)
      ALLOCATE(counter(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "counter",natom)
      ALLOCATE(atom_a(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_a",natom)
      ALLOCATE(atom_c(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_c",natom)
      connection(:,:)=0
      connectivity(:,:)=0
      counter(:)=0
      atom_a(:)=0
      atom_c(:)=0

! now we need to find a list of atoms that each atom in this bond is connected
! to
      DO iatom=1,natom
         DO ibond=1,nbond
            IF(bond_list(ibond)%a == iatom) THEN
               counter(iatom)=counter(iatom)+1
               connectivity(counter(iatom),iatom)=bond_list(ibond)%b
            ELSEIF(bond_list(ibond)%b == iatom)THEN
               counter(iatom)=counter(iatom)+1
               connectivity(counter(iatom),iatom)=bond_list(ibond)%a
            ENDIF
         ENDDO
      ENDDO

! now I need to do a depth first search to figure out which atoms are on atom a's
! side and which are on atom c's
      atom_a(:)=0
      atom_a(bend_list(bend_number)%a)=1
      CALL depth_first_search(bend_list(bend_number)%a,bend_list(bend_number)%b,&
           connectivity(:,:),atom_a(:))
      atom_c(:)=0
      atom_c(bend_list(bend_number)%c)=1
      CALL depth_first_search(bend_list(bend_number)%c,bend_list(bend_number)%b,&
           connectivity(:,:),atom_c(:))

! now figure out the masses of the various sides, so we can weight how far we move each
! group of atoms
      mass_a=0.0_dp
      mass_c=0.0_dp
      DO iatom=1,natom
         CALL get_atomic_kind(particles%els(iatom)%atomic_kind,&
            mass=atom_mass)
         IF(atom_a(iatom) == 1) mass_a=mass_a+atom_mass
         IF(atom_c(iatom) == 1) mass_c=mass_c+atom_mass
      ENDDO

! choose a displacement
      IF(ionode) rand=next_random_number(rng_stream,error=error)
      CALL mp_bcast(rand,source,group)

      dis_angle=rmangle(molecule_type)*2.0E0_dp*(rand-0.5E0_dp)

! need to find the A-B-C bisector

! this going to be tough...we need to find the plane of the A-B-C bond and only shift
! that component for all atoms connected to A and C...otherwise we change other
! internal degrees of freedom

! find the bond vectors
      DO i=1,3
         bond_a(i)=r_new(i,bend_list(bend_number)%a)-&
            r_new(i,bend_list(bend_number)%b)
         bond_c(i)=r_new(i,bend_list(bend_number)%c)-&
            r_new(i,bend_list(bend_number)%b)
      ENDDO
      old_length_a=SQRT(DOT_PRODUCT(bond_a,bond_a))
      old_length_c=SQRT(DOT_PRODUCT(bond_c,bond_c))
      old_angle=ACOS(DOT_PRODUCT(bond_a,bond_c)/(old_length_a*old_length_c))

      DO i=1,3
         bisector(i)=bond_a(i)/old_length_a+& ! not yet normalized
                 bond_c(i)/old_length_c
      ENDDO
      bis_length=SQRT(DOT_PRODUCT(bisector,bisector))
      bisector(1:3)=bisector(1:3)/bis_length

! now we need to find the cross product of the B-A and B-C vectors and normalize
! it, so we have a vector that defines the bend plane
      cross_prod(1)=bond_a(2)*bond_c(3)-bond_a(3)*bond_c(2)
      cross_prod(2)=bond_a(3)*bond_c(1)-bond_a(1)*bond_c(3)
      cross_prod(3)=bond_a(1)*bond_c(2)-bond_a(2)*bond_c(1)
      cross_prod(1:3)=cross_prod(1:3)/SQRT(DOT_PRODUCT(cross_prod,cross_prod))

! we have two axis of a coordinate system...let's get the third
      cross_prod_plane(1)=cross_prod(2)*bisector(3)-cross_prod(3)*bisector(2)
      cross_prod_plane(2)=cross_prod(3)*bisector(1)-cross_prod(1)*bisector(3)
      cross_prod_plane(3)=cross_prod(1)*bisector(2)-cross_prod(2)*bisector(1)
      cross_prod_plane(1:3)=cross_prod_plane(1:3)/&
         SQRT(DOT_PRODUCT(cross_prod_plane,cross_prod_plane))

! now bisector is x, cross_prod_plane is the y vector (pointing towards c),
! and cross_prod is z
! shift the molecule so that atom b is at the origin
      DO iatom=1,natom
         r_new(1:3,iatom)=r_new(1:3,iatom)-&
            r_old(1:3,bend_list(bend_number)%b)
      ENDDO

! figure out how much we move each side, since we're mass-weighting, by the
! opposite masses, so lighter moves farther..this angle is the angle between
! the bond vector BA or BC and the bisector
      dis_angle_a=dis_angle*mass_c/(mass_a+mass_c)
      dis_angle_c=dis_angle*mass_a/(mass_a+mass_c)

! now loop through all the atoms, moving the ones that are connected to a or c
      DO iatom=1,natom
! subtract out the z component (perpendicular to the angle plane)
         temp(1:3)=r_new(1:3,iatom)-&
            DOT_PRODUCT(cross_prod(1:3),r_new(1:3,iatom))*&
            cross_prod(1:3)
         temp_length=SQRT(DOT_PRODUCT(temp,temp))

! we can now compute all three components of the new bond vector along the
! axis defined above
         IF(atom_a(iatom) == 1) THEN

! if the y-coordinate is less than zero, we need to switch the sign when we make the vector,
! as the angle computed by the dot product can't distinguish between that
            IF(DOT_PRODUCT(cross_prod_plane(1:3),r_new(1:3,iatom)) &
               .LT. 0.0_dp) THEN

! need to figure out the current iatom-B-bisector angle, so we know what the new angle is
            new_angle_a=ACOS(DOT_PRODUCT(bisector,temp(1:3))/&
               (temp_length))+dis_angle_a

            r_new(1:3,iatom)=COS(new_angle_a)*temp_length*bisector(1:3)-&
               SIN(new_angle_a)*temp_length*cross_prod_plane(1:3)+&
               DOT_PRODUCT(cross_prod(1:3),r_new(1:3,iatom))*&
               cross_prod(1:3)
            ELSE

! need to figure out the current iatom-B-bisector angle, so we know what the new angle is
            new_angle_a=ACOS(DOT_PRODUCT(bisector,temp(1:3))/&
               (temp_length))-dis_angle_a

            r_new(1:3,iatom)=COS(new_angle_a)*temp_length*bisector(1:3)+&
               SIN(new_angle_a)*temp_length*cross_prod_plane(1:3)+&
               DOT_PRODUCT(cross_prod(1:3),r_new(1:3,iatom))*&
               cross_prod(1:3)
            ENDIF

         ELSEIF(atom_c(iatom) == 1) THEN

! if the y-coordinate is less than zero, we need to switch the sign when we make the vector,
! as the angle computed by the dot product can't distinguish between that
            IF(DOT_PRODUCT(cross_prod_plane(1:3),r_new(1:3,iatom)) &
               .LT. 0.0_dp) THEN
! need to figure out the current iatom-B-bisector angle, so we know what the new angle is
            new_angle_c=ACOS(DOT_PRODUCT(bisector(1:3),temp(1:3))/&
               (temp_length))-dis_angle_c

            r_new(1:3,iatom)=COS(new_angle_c)*temp_length*bisector(1:3)-&
               SIN(new_angle_c)*temp_length*cross_prod_plane(1:3)+&
               DOT_PRODUCT(cross_prod(1:3),r_new(1:3,iatom))*&
               cross_prod(1:3)
            ELSE
            new_angle_c=ACOS(DOT_PRODUCT(bisector(1:3),temp(1:3))/&
               (temp_length))+dis_angle_c

            r_new(1:3,iatom)=COS(new_angle_c)*temp_length*bisector(1:3)+&
               SIN(new_angle_c)*temp_length*cross_prod_plane(1:3)+&
               DOT_PRODUCT(cross_prod(1:3),r_new(1:3,iatom))*&
               cross_prod(1:3)
            ENDIF
         ENDIF

      ENDDO

      DO iatom=1,natom
         r_new(1:3,iatom)=r_new(1:3,iatom)+&
            r_old(1:3,bend_list(bend_number)%b)
      ENDDO

! deallocate some stuff
      DEALLOCATE(connection,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "connection")
      DEALLOCATE(connectivity,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "connectivity")
      DEALLOCATE(counter,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "counter")
      DEALLOCATE(atom_a,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_a")
      DEALLOCATE(atom_c,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_c")

! end the timing
      CALL timestop(handle)

  END SUBROUTINE change_bond_angle

! *****************************************************************************
!> \brief Alters a dihedral (A-B-C-D) in the molecule so that all other internal
!>      degrees of freedom remain the same.  If other dihedrals are centered
!>      on B-C, they rotate as well to keep the relationship between the
!>      dihedrals the same.  Atoms A and D are moved amounts related to their
!>      masses (and masses of all connecting atoms), so that heavier segments
!>      are moved less.  All atoms except B and C are rotated around the
!>      B-C bond vector (B and C are not moved).
!> \param r_old the initial coordinates of all molecules in the system
!> \param r_new the new coordinates of all molecules in the system
!> \param mc_par the mc parameters for the force env
!> \param molecule_type the type of molecule we're playing with
!> \param molecule_kind the structure containing the molecule information
!> \param particles the particle_list_type for all particles in the force_env..
!>             used to grab the mass of each atom
!> \param rng_stream the random number stream that we draw from
!> \param error ...
!> \author MJM
! *****************************************************************************
  SUBROUTINE change_dihedral ( r_old,r_new,mc_par,molecule_type,molecule_kind,&
      particles,rng_stream,error)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_old
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: r_new
    TYPE(mc_simpar_type), POINTER            :: mc_par
    INTEGER, INTENT(IN)                      :: molecule_type
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'change_dihedral', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, handle, i, iatom, &
                                                ibond, ipart, istat, natom, &
                                                nbond, ntorsion, source, &
                                                torsion_number
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_a, atom_d, counter
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: connection, connectivity
    INTEGER, DIMENSION(:), POINTER           :: nunits
    LOGICAL                                  :: ionode
    REAL(dp), DIMENSION(:), POINTER          :: rmdihedral
    REAL(KIND=dp)                            :: atom_mass, dis_angle, &
                                                dis_angle_a, dis_angle_d, &
                                                mass_a, mass_d, old_length_a, &
                                                rand, u, v, w, x, y, z
    REAL(KIND=dp), DIMENSION(1:3)            :: bond_a, temp
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: torsion_list

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

      NULLIFY(rmdihedral,torsion_list,bond_list,mc_molecule_info)

! get some stuff from mc_par
      CALL get_mc_par(mc_par,rmdihedral=rmdihedral,&
         source=source,group=group,ionode=ionode,&
         mc_molecule_info=mc_molecule_info)
      CALL get_mc_molecule_info(mc_molecule_info,nunits=nunits)

! copy the incoming coordinates so we can change them
      DO ipart=1,nunits(molecule_type)
         r_new(1:3,ipart)=r_old(1:3,ipart)
      ENDDO

! pick which bond in the molecule at random
      IF(ionode) THEN
         rand=next_random_number(rng_stream,error=error)
!      CALL RANDOM_NUMBER(rand)
      ENDIF
      CALL mp_bcast(rand,source,group)
      CALL get_molecule_kind(molecule_kind,natom=natom,&
         bond_list=bond_list,nbond=nbond,&
         ntorsion=ntorsion,torsion_list=torsion_list)
      torsion_number=CEILING(rand*REAL(ntorsion,dp))

      ALLOCATE(connection(1:natom,1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "connection",2*natom*dp_size)
! assume at most six bonds per atom
      ALLOCATE(connectivity(1:6,1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "connectivity",6*natom*dp_size)
      ALLOCATE(counter(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "counter",natom)
      ALLOCATE(atom_a(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_a",natom)
      ALLOCATE(atom_d(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_d",natom)
      connection(:,:)=0
      connectivity(:,:)=0
      counter(:)=0
      atom_a(:)=0
      atom_d(:)=0

! now we need to find a list of atoms that each atom in this bond is connected
! to
      DO iatom=1,natom
         DO ibond=1,nbond
            IF(bond_list(ibond)%a == iatom) THEN
               counter(iatom)=counter(iatom)+1
               connectivity(counter(iatom),iatom)=bond_list(ibond)%b
            ELSEIF(bond_list(ibond)%b == iatom)THEN
               counter(iatom)=counter(iatom)+1
               connectivity(counter(iatom),iatom)=bond_list(ibond)%a
            ENDIF
         ENDDO
      ENDDO

! now I need to do a depth first search to figure out which atoms are on atom
! a's side and which are on atom d's, but remember we're moving all atoms on a's
! side of b, including atoms not in a's branch
      atom_a(:)=0
      atom_a(torsion_list(torsion_number)%a)=1
      CALL depth_first_search(torsion_list(torsion_number)%b,&
         torsion_list(torsion_number)%c,connectivity(:,:),atom_a(:))
      atom_d(:)=0
      atom_d(torsion_list(torsion_number)%d)=1
      CALL depth_first_search(torsion_list(torsion_number)%c,&
         torsion_list(torsion_number)%b,connectivity(:,:),atom_d(:))

! now figure out the masses of the various sides, so we can weight how far we
! move each group of atoms
      mass_a=0.0_dp
      mass_d=0.0_dp
      DO iatom=1,natom
         CALL get_atomic_kind(particles%els(iatom)%atomic_kind,&
            mass=atom_mass)
         IF(atom_a(iatom) == 1) mass_a=mass_a+atom_mass
         IF(atom_d(iatom) == 1) mass_d=mass_d+atom_mass
      ENDDO

! choose a displacement
      IF(ionode) rand=next_random_number(rng_stream,error=error)
      CALL mp_bcast(rand,source,group)

      dis_angle=rmdihedral(molecule_type)*2.0E0_dp*(rand-0.5E0_dp)

! find the bond vectors, B-C, so we know what to rotate around
      DO i=1,3
         bond_a(i)=r_new(i,torsion_list(torsion_number)%c)-&
            r_new(i,torsion_list(torsion_number)%b)
      ENDDO
      old_length_a=SQRT(DOT_PRODUCT(bond_a,bond_a))
      bond_a(1:3)=bond_a(1:3)/old_length_a

! figure out how much we move each side, since we're mass-weighting, by the
! opposite masses, so lighter moves farther...we take the opposite sign of d
! so we're not rotating both angles in the same direction
      dis_angle_a=dis_angle*mass_d/(mass_a+mass_d)
      dis_angle_d=-dis_angle*mass_a/(mass_a+mass_d)

      DO iatom=1,natom

         IF(atom_a(iatom) == 1) THEN
! shift the coords so b is at the origin
            r_new(1:3,iatom)=r_new(1:3,iatom)-&
               r_new(1:3,torsion_list(torsion_number)%b)

! multiply by the rotation matrix
               u=bond_a(1)
               v=bond_a(2)
               w=bond_a(3)
               x=r_new(1,iatom)
               y=r_new(2,iatom)
               z=r_new(3,iatom)
               temp(1)=(u*(u*x+v*y+w*z)+(x*(v**2+w**2)-u*(v*y+w*z))*COS(dis_angle_a)+&
                  SQRT(u**2+v**2+w**2)*(v*z-w*y)*SIN(dis_angle_a))/(u**2+v**2+w**2)
               temp(2)=(v*(u*x+v*y+w*z)+(y*(u**2+w**2)-v*(u*x+w*z))*COS(dis_angle_a)+&
                  SQRT(u**2+v**2+w**2)*(w*x-u*z)*SIN(dis_angle_a))/(u**2+v**2+w**2)
               temp(3)=(w*(u*x+v*y+w*z)+(z*(v**2+u**2)-w*(u*x+v*y))*COS(dis_angle_a)+&
                  SQRT(u**2+v**2+w**2)*(u*y-v*x)*SIN(dis_angle_a))/(u**2+v**2+w**2)

! shift back to the original position
            temp(1:3)=temp(1:3)+r_new(1:3,torsion_list(torsion_number)%b)
            r_new(1:3,iatom)=temp(1:3)

         ELSEIF(atom_d(iatom) == 1) THEN

! shift the coords so c is at the origin
            r_new(1:3,iatom)=r_new(1:3,iatom)-&
               r_new(1:3,torsion_list(torsion_number)%c)

! multiply by the rotation matrix
               u=bond_a(1)
               v=bond_a(2)
               w=bond_a(3)
               x=r_new(1,iatom)
               y=r_new(2,iatom)
               z=r_new(3,iatom)
               temp(1)=(u*(u*x+v*y+w*z)+(x*(v**2+w**2)-u*(v*y+w*z))*COS(dis_angle_d)+&
                  SQRT(u**2+v**2+w**2)*(v*z-w*y)*SIN(dis_angle_d))/(u**2+v**2+w**2)
               temp(2)=(v*(u*x+v*y+w*z)+(y*(u**2+w**2)-v*(u*x+w*z))*COS(dis_angle_d)+&
                  SQRT(u**2+v**2+w**2)*(w*x-u*z)*SIN(dis_angle_d))/(u**2+v**2+w**2)
               temp(3)=(w*(u*x+v*y+w*z)+(z*(v**2+u**2)-w*(u*x+v*y))*COS(dis_angle_d)+&
                  SQRT(u**2+v**2+w**2)*(u*y-v*x)*SIN(dis_angle_d))/(u**2+v**2+w**2)

! shift back to the original position
            temp(1:3)=temp(1:3)+r_new(1:3,torsion_list(torsion_number)%c)
            r_new(1:3,iatom)=temp(1:3)
         ENDIF
      ENDDO

! deallocate some stuff
      DEALLOCATE(connection,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "connection")
      DEALLOCATE(connectivity,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "connectivity")
      DEALLOCATE(counter,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "counter")
      DEALLOCATE(atom_a,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_a")
      DEALLOCATE(atom_d,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_d")

! end the timing
      CALL timestop(handle)

  END SUBROUTINE change_dihedral

! *****************************************************************************
!> \brief performs either a bond or angle change move for a given molecule
!> \param mc_par the mc parameters for the force env
!> \param force_env the force environment used in the move
!> \param bias_env the force environment used to bias the move, if any (it may
!>            be null if lbias=.false. in mc_par)
!> \param moves the structure that keeps track of how many moves have been
!>               accepted/rejected
!> \param energy_check the running energy difference between now and the inital
!>        energy
!> \param r_old the coordinates of force_env before the move
!> \param old_energy the energy of the force_env before the move
!> \param start_atom_swap the number of the swap molecule's first atom, assuming the rest of
!>        the atoms follow sequentially
!> \param target_atom the number of the target atom for swapping
!> \param molecule_type the molecule type for the atom we're swapping
!> \param box_number the number of the box we're doing this move in
!> \param bias_energy_old the biased energy of the system before the move
!> \param last_bias_energy the last biased energy of the system
!> \param move_type dictates if we're moving to an "in" or "out" region
!> \param rng_stream the random number stream that we draw from
!> \param error the cp_error type
!>
!> \author MJM
!> \note     Designed for parallel.
! *****************************************************************************
  SUBROUTINE mc_avbmc_move ( mc_par,force_env,bias_env, moves,&
                        energy_check,r_old,old_energy,start_atom_swap,&
                        target_atom,&
                        molecule_type,box_number,bias_energy_old,last_bias_energy,&
                        move_type,rng_stream,error)

    TYPE(mc_simpar_type), POINTER            :: mc_par
    TYPE(force_env_type), POINTER            :: force_env, bias_env
    TYPE(mc_moves_type), POINTER             :: moves
    REAL(KIND=dp), INTENT(INOUT)             :: energy_check
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: r_old
    REAL(KIND=dp), INTENT(INOUT)             :: old_energy
    INTEGER, INTENT(IN)                      :: start_atom_swap, target_atom, &
                                                molecule_type, box_number
    REAL(KIND=dp), INTENT(INOUT)             :: bias_energy_old, &
                                                last_bias_energy
    CHARACTER(LEN=*), INTENT(IN)             :: move_type
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_avbmc_move', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: end_mol, group, handle, &
                                                ipart, istat, jbox, natom, &
                                                nswapmoves, source, start_mol
    INTEGER, DIMENSION(:), POINTER           :: avbmc_atom, mol_type, nunits, &
                                                nunits_tot
    INTEGER, DIMENSION(:, :), POINTER        :: nchains
    LOGICAL                                  :: ionode, lbias, ldum, lin, &
                                                loverlap
    REAL(dp), DIMENSION(:), POINTER          :: avbmc_rmax, avbmc_rmin, pbias
    REAL(dp), DIMENSION(:, :), POINTER       :: mass
    REAL(KIND=dp) :: BETA, bias_energy_new, del_quickstep_energy, distance, &
      exp_max_val, exp_min_val, max_val, min_val, new_energy, prefactor, &
      rand, rdum, volume_in, volume_out, w, weight_new, weight_old
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r_new
    REAL(KIND=dp), DIMENSION(1:3)            :: abc, RIJ
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys, subsys_force
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_list_type), POINTER        :: particles, particles_force

      rdum=1.0_dp

! begin the timing of the subroutine
      CALL timeset(routineN,handle)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,lbias=lbias,&
         BETA=BETA,max_val=max_val, min_val=min_val, exp_max_val=exp_max_val,&
         exp_min_val=exp_min_val,avbmc_atom=avbmc_atom,&
         avbmc_rmin=avbmc_rmin,avbmc_rmax=avbmc_rmax,&
         nswapmoves=nswapmoves,ionode=ionode,source=source,&
         group=group,pbias=pbias,mc_molecule_info=mc_molecule_info)
      CALL get_mc_molecule_info(mc_molecule_info,nchains=nchains,&
           mass=mass,nunits=nunits,nunits_tot=nunits_tot,mol_type=mol_type)
! figure out some bounds for mol_type
      start_mol=1
      DO jbox=1,box_number-1
         start_mol=start_mol+SUM(nchains(:,jbox))
      ENDDO
      end_mol=start_mol+SUM(nchains(:,box_number))-1

! nullify some pointers
      NULLIFY(particles,subsys,molecule_kinds_new,molecule_kind,&
         particles_force,subsys_force)

! do some allocation
      ALLOCATE (r_new(1:3,1:nunits_tot(box_number)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r_new",3*nunits_tot(box_number)*dp_size)

! now we need to grab and save coordinates, in case we reject
! are we biasing this move?
      IF(lbias) THEN

! grab the coordinates
         CALL force_env_get(bias_env,cell=cell,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys, &
            particles=particles, molecule_kinds_new=molecule_kinds_new,&
            error=error)
         molecule_kind => molecule_kinds_new%els(1)
         CALL get_molecule_kind(molecule_kind,natom=natom)
         CALL get_cell(cell,abc=abc)

! save the energy
!         bias_energy_old=bias_energy

      ELSE

! grab the coordinates
         CALL force_env_get(force_env,cell=cell,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys, &
            particles=particles,  molecule_kinds_new=molecule_kinds_new,&
            error=error)
         molecule_kind => molecule_kinds_new%els(1)
         CALL get_molecule_kind(molecule_kind,natom=natom)
         CALL get_cell(cell,abc=abc)

      ENDIF

! let's determine if the molecule to be moved is in the "in" region or the
! "out" region of the target
      RIJ(1)=particles%els(start_atom_swap+avbmc_atom(molecule_type)-1)%r(1)-&
         particles%els(target_atom)%r(1)-abc(1)*ANINT(&
         (particles%els(start_atom_swap+avbmc_atom(molecule_type)-1)%r(1)-&
         particles%els(target_atom)%r(1))/abc(1))
      RIJ(2)=particles%els(start_atom_swap+avbmc_atom(molecule_type)-1)%r(2)-&
         particles%els(target_atom)%r(2)-abc(2)*ANINT(&
         (particles%els(start_atom_swap+avbmc_atom(molecule_type)-1)%r(2)-&
         particles%els(target_atom)%r(2))/abc(2))
      RIJ(3)=particles%els(start_atom_swap+avbmc_atom(molecule_type)-1)%r(3)-&
         particles%els(target_atom)%r(3)-abc(3)*ANINT(&
         (particles%els(start_atom_swap+avbmc_atom(molecule_type)-1)%r(3)-&
         particles%els(target_atom)%r(3))/abc(3))
      distance=SQRT(RIJ(1)**2+RIJ(2)**2+RIJ(3)**2)
      IF(distance .LE. avbmc_rmax(molecule_type) .AND. distance .GE. avbmc_rmin(molecule_type)) THEN
         lin=.TRUE.
      ELSE
         lin=.FALSE.
      ENDIF

! increment the counter of the particular move we've done
!     swapping into the "in" region of mol_target
      IF(lin) THEN
         IF(move_type == 'in' ) THEN
            moves%avbmc_inin%attempts=&
                 moves%avbmc_inin%attempts+1
         ELSE
            moves%avbmc_inout%attempts=&
                 moves%avbmc_inout%attempts+1
         ENDIF
      ELSE
         IF(move_type == 'in' ) THEN
            moves%avbmc_outin%attempts=&
                 moves%avbmc_outin%attempts+1
         ELSE
            moves%avbmc_outout%attempts=&
                 moves%avbmc_outout%attempts+1
         ENDIF
      ENDIF

      IF(lbias) THEN

         IF(move_type == 'in') THEN

! do CBMC for the old config
            CALL generate_cbmc_swap_config( bias_env, BETA,  max_val, min_val, exp_max_val,&
                 exp_min_val,nswapmoves,&
                 weight_old,start_atom_swap,nunits_tot(box_number),nunits,nunits(molecule_type),&
                 mass(:,molecule_type),ldum, rdum,&
                 bias_energy_old,ionode,.TRUE.,mol_type(start_mol:end_mol),nchains(:,box_number),&
                 source,group,rng_stream,error,&
                 avbmc_atom=avbmc_atom(molecule_type),&
                 rmin=avbmc_rmin(molecule_type),rmax=avbmc_rmax(molecule_type),move_type='out',&
                 target_atom=target_atom)

         ELSE

! do CBMC for the old config
           CALL generate_cbmc_swap_config( bias_env, BETA,  max_val, min_val, exp_max_val,&
                exp_min_val,nswapmoves,&
                weight_old,start_atom_swap,nunits_tot(box_number),nunits,nunits(molecule_type),&
                mass(:,molecule_type),ldum, rdum,&
                bias_energy_old,ionode,.TRUE.,mol_type(start_mol:end_mol),nchains(:,box_number),&
                source,group,rng_stream,error,&
                avbmc_atom=avbmc_atom(molecule_type),&
                rmin=avbmc_rmin(molecule_type),rmax=avbmc_rmax(molecule_type),move_type='in',&
                target_atom=target_atom)

         ENDIF

! generate the new config
         CALL generate_cbmc_swap_config( bias_env, BETA, max_val, min_val, exp_max_val,&
              exp_min_val, nswapmoves,&
              weight_new,start_atom_swap,nunits_tot(box_number),nunits,nunits(molecule_type),&
              mass(:,molecule_type),loverlap, bias_energy_new,&
              bias_energy_old,ionode,.FALSE.,mol_type(start_mol:end_mol),nchains(:,box_number),&
              source,group,rng_stream,error,&
              avbmc_atom=avbmc_atom(molecule_type),&
              rmin=avbmc_rmin(molecule_type),rmax=avbmc_rmax(molecule_type),move_type=move_type,&
              target_atom=target_atom)

! the energy that comes out of the above routine is the difference...we want
! the real energy for the acceptance rule...we don't do this for the
! lbias=.false. case because it doesn't appear in the acceptance rule, and
! we compensate in case of acceptance
         bias_energy_new=bias_energy_new+bias_energy_old

      ELSE

         IF(move_type == 'in') THEN

! find the weight of the old config
            CALL generate_cbmc_swap_config( force_env, BETA,  max_val, min_val, exp_max_val,&
                 exp_min_val,nswapmoves,&
                 weight_old,start_atom_swap,nunits_tot(box_number),nunits,nunits(molecule_type),&
                 mass(:,molecule_type),ldum,rdum,old_energy,&
                 ionode,.TRUE.,mol_type(start_mol:end_mol),nchains(:,box_number),source,group,rng_stream,error,&
                 avbmc_atom=avbmc_atom(molecule_type),&
                 rmin=avbmc_rmin(molecule_type),rmax=avbmc_rmax(molecule_type),move_type='out',&
                 target_atom=target_atom)

         ELSE

! find the weight of the old config
            CALL generate_cbmc_swap_config( force_env, BETA,  max_val, min_val, exp_max_val,&
                 exp_min_val,nswapmoves,&
                 weight_old,start_atom_swap,nunits_tot(box_number),nunits,nunits(molecule_type),&
                 mass(:,molecule_type),ldum,rdum,old_energy,&
                 ionode,.TRUE.,mol_type(start_mol:end_mol),nchains(:,box_number),source,group,rng_stream,error,&
                 avbmc_atom=avbmc_atom(molecule_type),&
                 rmin=avbmc_rmin(molecule_type),rmax=avbmc_rmax(molecule_type),move_type='in',&
                 target_atom=target_atom)

         ENDIF

 ! generate the new config...do this after, because it changes the force_env
        CALL generate_cbmc_swap_config( force_env, BETA,  max_val, min_val, exp_max_val,&
             exp_min_val,nswapmoves,&
             weight_new,start_atom_swap,nunits_tot(box_number),nunits,nunits(molecule_type),&
             mass(:,molecule_type),loverlap,new_energy,old_energy,&
             ionode,.FALSE.,mol_type(start_mol:end_mol),nchains(:,box_number),source,group,rng_stream,error,&
             avbmc_atom=avbmc_atom(molecule_type),&
             rmin=avbmc_rmin(molecule_type),rmax=avbmc_rmax(molecule_type),move_type=move_type,&
             target_atom=target_atom)

      ENDIF

      IF(loverlap) THEN
         DEALLOCATE(r_new,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "r_new")

! need to reset the old coordinates
         IF(lbias) THEN
            CALL force_env_get(bias_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys,particles=particles,error=error)
         ELSE
            CALL force_env_get(force_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys,particles=particles,error=error)
         ENDIF
         DO ipart=1,nunits_tot(box_number)
            particles%els(ipart)%r(1:3)=r_old(1:3,ipart)
         ENDDO

         CALL timestop(handle)

         RETURN
      ENDIF

! if we're biasing, we need to compute the new energy with the full
! potential
      IF(lbias) THEN
! need to give the force_env the coords from the bias_env
         CALL force_env_get(force_env,subsys=subsys_force,error=error)
         CALL cp_subsys_get(subsys_force,particles=particles_force,error=error)
         CALL force_env_get(bias_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys,particles=particles,error=error)
         DO ipart=1,nunits_tot(box_number)
            particles_force%els(ipart)%r(1:3)=particles%els(ipart)%r(1:3)
         ENDDO

         CALL force_env_calc_energy_force(force_env,&
            calc_force=.FALSE.,error=error)
         CALL force_env_get(force_env,&
            potential_energy=new_energy,error=error)

      ENDIF

      volume_in=4.0_dp/3.0_dp*pi*(avbmc_rmax(molecule_type)**3-avbmc_rmin(molecule_type)**3)
      volume_out=abc(1)*abc(2)*abc(3)-volume_in

      IF(lin .AND. move_type == 'in' .OR. &
         .NOT. lin .AND. move_type == 'out') THEN
! standard Metropolis rule
         prefactor=1.0_dp
      ELSEIF(.NOT. lin .AND. move_type == 'in') THEN
         prefactor=(1.0_dp-pbias(molecule_type))*volume_in/(pbias(molecule_type)*volume_out)
      ELSE
         prefactor=pbias(molecule_type)*volume_out/((1.0_dp-pbias(molecule_type))*volume_in)
      ENDIF

      IF(lbias) THEN
! AVBMC with CBMC and a biasing potential...notice that if the biasing
! potential equals the quickstep potential, this cancels out to the
! acceptance below
         del_quickstep_energy=(-BETA)*(new_energy-old_energy-&
            (bias_energy_new-bias_energy_old))

         IF     (del_quickstep_energy .GT. exp_max_val) THEN
            del_quickstep_energy=max_val
         ELSEIF (del_quickstep_energy .LT. exp_min_val) THEN
            del_quickstep_energy=0.0_dp
         ELSE
            del_quickstep_energy=EXP(del_quickstep_energy)
         ENDIF

         w=prefactor*del_quickstep_energy*weight_new/weight_old

      ELSE

! AVBMC with CBMC
         w=prefactor*weight_new/weight_old
      ENDIF

! check if the move is accepted
      IF(w .GE. 1.0E0_dp) THEN
         rand=0.0E0_dp
      ELSE
         IF(ionode) rand=next_random_number(rng_stream,error=error)
         CALL mp_bcast(rand,source,group)
      ENDIF

      IF ( rand .LT. w ) THEN

! accept the move

         IF(lin) THEN
            IF(move_type == 'in' ) THEN
               moves%avbmc_inin%successes=&
                    moves%avbmc_inin%successes+1
            ELSE
               moves%avbmc_inout%successes=&
                    moves%avbmc_inout%successes+1
            ENDIF
         ELSE
            IF(move_type == 'in' ) THEN
               moves%avbmc_outin%successes=&
                    moves%avbmc_outin%successes+1
            ELSE
               moves%avbmc_outout%successes=&
                    moves%avbmc_outout%successes+1
            ENDIF
         ENDIF

! we need to compensate for the fact that we take the difference in
! generate_cbmc_config to keep the exponetials small
         IF(.NOT. lbias) THEN
            new_energy=new_energy+old_energy
         ENDIF

! update energies
         energy_check=energy_check+(new_energy-old_energy)
         old_energy=new_energy

! if we're biasing the update the biasing energy
         IF (lbias) THEN
! need to do this outside of the routine
            last_bias_energy=bias_energy_new
            bias_energy_old=bias_energy_new
         ENDIF

! update coordinates
         CALL force_env_get(force_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys,particles=particles,error=error)
         DO ipart=1,nunits_tot(box_number)
            r_old(1:3,ipart)=particles%els(ipart)%r(1:3)
         ENDDO
      ELSE
! reject the move...need to restore the old coordinates
         IF(lbias) THEN
            CALL force_env_get(bias_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys,particles=particles, error=error)
            DO ipart=1,nunits_tot(box_number)
               particles%els(ipart)%r(1:3)=r_old(1:3,ipart)
            ENDDO
            CALL cp_subsys_set(subsys,particles=particles,error=error)
         ENDIF
         CALL force_env_get(force_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys,particles=particles, error=error)
         DO ipart=1,nunits_tot(box_number)
            particles%els(ipart)%r(1:3)=r_old(1:3,ipart)
         ENDDO
         CALL cp_subsys_set(subsys,particles=particles,error=error)

      ENDIF

! deallocate some stuff
      DEALLOCATE(r_new,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r_new")
! end the timing
  CALL timestop(handle)

  END SUBROUTINE mc_avbmc_move

! *****************************************************************************
!> \brief performs a hybrid Monte Carlo move that runs a short MD sequence
!> \param mc_par the mc parameters for the force env
!> \param force_env the force environment whose cell we're changing
!> \param globenv ...
!> \param moves the structure that keeps track of how many moves have been
!>               accepted/rejected
!> \param move_updates the structure that keeps track of how many moves have
!>               been accepted/rejected since the last time the displacements
!>               were updated
!> \param old_energy the energy of the last accepted move involving an
!>                    unbiased calculation
!> \param box_number the box we're changing the volume of
!> \param energy_check the running total of how much the energy has changed
!>                      since the initial configuration
!> \param r_old the coordinates of the last accepted move involving an
!>               unbiased calculation
!> \param rng_stream the random number stream that we draw from
!> \param error the cp error type
!>
!> \author MJM
!> \note     Designed for parallel use.
! *****************************************************************************
  SUBROUTINE mc_hmc_move ( mc_par,force_env, globenv, moves,move_updates,&
                        old_energy,box_number,&
                        energy_check,r_old,rng_stream,error)

    TYPE(mc_simpar_type), POINTER            :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    REAL(KIND=dp), INTENT(INOUT)             :: old_energy
    INTEGER, INTENT(IN)                      :: box_number
    REAL(KIND=dp), INTENT(INOUT)             :: energy_check
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: r_old
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_hmc_move', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, handle, iatom, istat, &
                                                source
    INTEGER, DIMENSION(:), POINTER           :: nunits_tot
    LOGICAL                                  :: ionode
    REAL(KIND=dp)                            :: BETA, energy_term, &
                                                exp_max_val, exp_min_val, &
                                                new_energy, rand, value, w
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r
    TYPE(cp_subsys_type), POINTER            :: oldsys
    TYPE(mc_ekin_type), POINTER              :: hmc_ekin
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    TYPE(particle_list_type), POINTER        :: particles_old

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,ionode=ionode,&
         BETA=BETA,exp_max_val=exp_max_val,&
         exp_min_val=exp_min_val,source=source,group=group,&
         mc_molecule_info=mc_molecule_info)
      CALL get_mc_molecule_info(mc_molecule_info,nunits_tot=nunits_tot)

! nullify some pointers
      NULLIFY(particles_old,oldsys,hmc_ekin)

! do some allocation
      ALLOCATE (r(1:3,1:nunits_tot(box_number)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "r",3*nunits_tot(box_number)*dp_size)
      ALLOCATE (hmc_ekin,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "hmc_ekin")

! record the attempt
      moves%hmc%attempts=moves%hmc%attempts+1
      move_updates%hmc%attempts=move_updates%hmc%attempts+1

! now let's grab the particle positions
      CALL force_env_get(force_env,subsys=oldsys,&
         error=error)
      CALL cp_subsys_get(oldsys,particles=particles_old, &
         error=error)

! save the old coordinates
      DO iatom=1,nunits_tot(box_number)
         r(1:3,iatom)=particles_old%els(iatom)%r(1:3)
      ENDDO

! now run the MD simulation
      CALL qs_mol_dyn(force_env,globenv,error=error,hmc_e_initial=hmc_ekin%initial_ekin,hmc_e_final=hmc_ekin%final_ekin)


! get the energy
      CALL force_env_get(force_env,&
         potential_energy=new_energy,error=error)

! accept or reject the move
! to prevent overflows
      energy_term=new_energy+hmc_ekin%final_ekin-old_energy-hmc_ekin%initial_ekin

      value=-BETA*(energy_term)
      IF    (value .GT. exp_max_val) THEN
         w=10.0_dp
      ELSEIF(value .LT. exp_min_val) THEN
         w=0.0_dp
      ELSE
         w=EXP(value)
      ENDIF

      IF ( w .GE. 1.0E0_dp ) THEN
         w=1.0E0_dp
         rand=0.0E0_dp
      ELSE
         IF(ionode) rand=next_random_number(rng_stream,error=error)
         CALL mp_bcast(rand,source,group)
      ENDIF

      IF (rand .LT. w) THEN

! accept the move
         moves%hmc%successes=moves%hmc%successes+1
         move_updates%hmc%successes=move_updates%hmc%successes+1

! update energies
         energy_check=energy_check+(new_energy-old_energy)
         old_energy=new_energy

         DO iatom=1,nunits_tot(box_number)
            r_old(1:3,iatom)=particles_old%els(iatom)%r(1:3)
         ENDDO

      ELSE

! reset the cell and particle positions
         DO iatom=1,nunits_tot(box_number)
            particles_old%els(iatom)%r(1:3)=r_old(1:3,iatom)
         ENDDO

      ENDIF

! deallocate some stuff
      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r")
      DEALLOCATE(hmc_ekin,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "hmc_ekin")

! end the timing
      CALL timestop(handle)

    END SUBROUTINE mc_hmc_move
END MODULE mc_moves

