!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief different move types are applied
!> \par History
!>      11.2012 created [Mandes Schönherr] 
!> \author Mandes 11/2012
! *****************************************************************************

MODULE tmc_moves
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE mathlib,                         ONLY: dihedral_angle,&
                                             rotate_vector
  USE parallel_rng_types,              ONLY: get_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type,&
                                             set_rng_stream
  USE physcon,                         ONLY: boltzmann,&
                                             joule
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_calculations,                ONLY: center_of_mass,&
                                             geometrical_center,&
                                             get_scaled_cell,&
                                             nearest_distance
  USE tmc_move_types,                  ONLY: &
       mv_type_MD, mv_type_atom_swap, mv_type_atom_trans, &
       mv_type_gausian_adapt, mv_type_mol_rot, mv_type_mol_trans, &
       mv_type_proton_reorder, mv_type_volume_move, tmc_move_type
  USE tmc_tree_types,                  ONLY: status_frozen,&
                                             status_ok,&
                                             tree_type
  USE tmc_types,                       ONLY: tmc_atom_type,&
                                             tmc_param_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_moves'

  PUBLIC :: change_pos
  PUBLIC :: elements_in_new_subbox

  INTEGER, PARAMETER :: not_selected    = 0
  INTEGER, PARAMETER :: proton_donor    = -1
  INTEGER, PARAMETER :: proton_acceptor = 1

CONTAINS
! *****************************************************************************
!> \brief applying the preselected move type
!> \param tmc_params TMC parameters with dimensions ...
!> \param rng_stream random number stream
!> \param elem configuration to change
!> \param mv_conf temperature index for determinig the move size
!> \param new_subbox flag if new sub box should be crated
!> \param move_rejected return flag if during configurational change 
!>        configuration should still be accepted (not if e.g. atom/molecule 
!>        leave the sub box
!> \param last_acc_elem last accepted configuration for Gaussian Adaptation
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE change_pos(tmc_params, move_types, rng_stream, elem, mv_conf, &
                        new_subbox, move_rejected, last_acc_elem, error)
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(tmc_move_type), POINTER             :: move_types
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(tree_type), POINTER                 :: elem
    INTEGER                                  :: mv_conf
    LOGICAL                                  :: new_subbox, move_rejected
    TYPE(tree_type), OPTIONAL, POINTER       :: last_acc_elem
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'change_pos', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: act_nr_elem_mv, counter, d, &
                                                i, ind, ind_e, m, nr_molec, &
                                                nr_sub_box_elem
    INTEGER, DIMENSION(:), POINTER           :: mol_in_sb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rnd
    REAL(KIND=dp), DIMENSION(:), POINTER     :: direction, elem_center

    failure = .FALSE.
    NULLIFY(direction, elem_center, mol_in_sb)

    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(move_types),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rng_stream),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)

    move_rejected = .FALSE.

    CALL set_rng_stream(rng_stream=rng_stream, bg=elem%rng_seed(:,:,1), &
           cg=elem%rng_seed(:,:,2), ig=elem%rng_seed(:,:,3), error=error)

    IF(new_subbox) THEN
       IF(ALL(tmc_params%sub_box_size.GT.0.0_dp)) THEN
         CALL elements_in_new_subbox(tmc_params=tmc_params, &
                rng_stream=rng_stream, elem=elem, &
                nr_of_sub_box_elements=nr_sub_box_elem, error=error)
       ELSE
         elem%elem_stat(:) = status_ok
       END IF
    END IF

    ! at least one atom should be in the sub box
    CPPrecondition(ANY(elem%elem_stat(:).EQ.status_ok),cp_warning_level,routineP,error,failure)
    IF(failure) THEN
      act_nr_elem_mv = -1
    ELSE IF(tmc_params%nr_elem_mv.EQ.0) THEN
      ! move all elements (could be all atoms or all molecules)
      act_nr_elem_mv = 0 
    ELSE
      act_nr_elem_mv = tmc_params%nr_elem_mv
    END IF
    !-- select the type of move (looked up in list, using the move type index)
    !-- for each move type exist single moves of certain number of elements 
    !-- or move of all elements
    !-- one element is a position or velocity of an atom. 
    !-- Always all dimension are changed.
    SELECT CASE(elem%move_type)
    CASE(mv_type_gausian_adapt)
       ! just for Gaussian Adaptation
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                 "gaussian adaptation is not imlemented yet.",&
                 error,failure)
!TODO       CALL new_pos_gauss_adapt(acc=ASSOCIATED(elem%parent%acc, elem), &
!                    pos=elem%pos, covari=elem%frc, pot=elem%potential, &
!                    step_size=elem%ekin, pos_aver=elem%vel, temp=elem%ekin_before_md, &
!                    rng_seed=elem%rng_seed, rng_seed_last_acc=last_acc_elem%rng_seed)
    !-- atom translation
    CASE(mv_type_atom_trans)
      IF(act_nr_elem_mv.EQ.0) &
        act_nr_elem_mv = SIZE(elem%pos)/tmc_params%dim_per_elem
      ALLOCATE(elem_center(tmc_params%dim_per_elem))
      i=1
      move_elements_loop: DO
        ! select atom
        IF(tmc_params%nr_elem_mv.EQ.0) THEN
          ind = (i-1) * (tmc_params%dim_per_elem) +1
        ELSE
          rnd = next_random_number(rng_stream, error=error) ! next random number
          ind = tmc_params%dim_per_elem*&
                INT(rnd*(SIZE(elem%pos)/tmc_params%dim_per_elem))+1
        END IF
        ! apply move
        IF(elem%elem_stat(ind).EQ.status_ok)THEN
          ! displace atom
          DO d=0, tmc_params%dim_per_elem-1
            rnd = next_random_number(rng_stream, error=error) ! next random number
            elem%pos(ind+d) = elem%pos(ind+d) + (rnd-0.5)*2.0*&
              move_types%mv_size(mv_type_atom_trans,mv_conf)
          END DO
          ! check if new position is in subbox
          elem_center = elem%pos(ind:ind+tmc_params%dim_per_elem-1)
          IF(.NOT.check_pos_in_subbox(pos=elem_center, &
                   subbox_center=elem%subbox_center, &
                   box_scale=elem%box_scale, tmc_params=tmc_params, &
                   error=error)) THEN
            move_rejected = .TRUE.
            EXIT move_elements_loop
          END IF
        ELSE
          ! element was not in sub box, search new one instead
          IF(tmc_params%nr_elem_mv.GT.0) i = i-1
        END IF
        i=i+1
        IF(i.GT.act_nr_elem_mv) EXIT move_elements_loop
      END DO move_elements_loop
      DEALLOCATE(elem_center)

    !-- molecule translation
    CASE(mv_type_mol_trans)
      nr_molec = MAXVAL(elem%mol(:))
      ! if all particles should be displaced, set the amount of molecules
      IF(act_nr_elem_mv.EQ.0) &
        act_nr_elem_mv = nr_molec
      ALLOCATE(mol_in_sb(nr_molec))
      ALLOCATE(elem_center(tmc_params%dim_per_elem))
      mol_in_sb(:) = status_frozen
      ! check if any molecule is in sub_box
      DO m=1, nr_molec
        CALL get_mol_indeces(tmc_params=tmc_params, mol_arr=elem%mol, mol=m, &
                             start_ind=ind, end_ind=ind_e, error=error)
        CALL geometrical_center(pos=elem%pos(ind:ind_e+tmc_params%dim_per_elem-1), &
                                center=elem_center, error=error)
        IF(check_pos_in_subbox(pos=elem_center, &
             subbox_center=elem%subbox_center, &
             box_scale=elem%box_scale, tmc_params=tmc_params, &
             error=error)) &
          mol_in_sb(m) = status_ok
      END DO
      ! displace the selected amount of molecules
      IF(ANY(mol_in_sb(:).EQ.status_ok)) THEN
        ALLOCATE(direction(tmc_params%dim_per_elem))
        counter = 1
        move_molecule_loop: DO 
          ! select molecule
          IF(tmc_params%nr_elem_mv.EQ.0) THEN
            m = counter
          ELSE
            rnd = next_random_number(rng_stream, error=error) ! next random number
            m = INT(rnd*nr_molec)+1
          END IF
          CALL get_mol_indeces(tmc_params=tmc_params, mol_arr=elem%mol, mol=m, &
                               start_ind=ind, end_ind=ind_e, error=error)
          ! when "molecule" is single atom, search a new one
          IF(ind .EQ. ind_e) CYCLE move_molecule_loop

          ! calculate displacement
          !  move only molecules, with geom. center in subbox
          IF(mol_in_sb(m) .EQ. status_ok) THEN
            ! calculate displacement
            DO d=1, tmc_params%dim_per_elem
              rnd = next_random_number(rng_stream, error=error) ! next random number
              direction(d) = (rnd-0.5)* 2.0_dp*move_types%mv_size(&
                                         mv_type_mol_trans,mv_conf)
            END DO
            ! check if displaced position is still in subbox
            elem_center(:) = elem_center(:) + direction(:)
            IF(check_pos_in_subbox(pos=elem_center,&
                 subbox_center=elem%subbox_center, &
                 box_scale=elem%box_scale, tmc_params=tmc_params, &
                 error=error)) THEN
              ! apply move
              atom_in_mol_loop:DO i=ind, ind_e+tmc_params%dim_per_elem-1, tmc_params%dim_per_elem
                dim_loop:DO d=0, tmc_params%dim_per_elem-1
                  elem%pos(i+d) = elem%pos(i+d) + direction(d+1)
                END DO dim_loop
              END DO atom_in_mol_loop
            ELSE
              ! the whole move is rejected, because one element is outside the subbox
              move_rejected = .TRUE.
              EXIT move_molecule_loop
            END IF
          ELSE 
            ! element was not in sub box, search new one instead
            IF(tmc_params%nr_elem_mv.GT.0) counter = counter-1
          END IF
          counter = counter+1
          IF(counter.GT.act_nr_elem_mv) EXIT move_molecule_loop
        END DO move_molecule_loop
        DEALLOCATE(direction)
      END IF
        DEALLOCATE(elem_center)
      DEALLOCATE(mol_in_sb)

    !-- molecule rotation
    CASE(mv_type_mol_rot)
      nr_molec = MAXVAL(elem%mol(:))
      IF(act_nr_elem_mv.EQ.0) &
        act_nr_elem_mv = nr_molec
      ALLOCATE(mol_in_sb(nr_molec))
      ALLOCATE(elem_center(tmc_params%dim_per_elem))
      mol_in_sb(:) = status_frozen
      ! check if any molecule is in sub_box
      DO m=1, nr_molec
        CALL get_mol_indeces(tmc_params=tmc_params, mol_arr=elem%mol, mol=m, &
                             start_ind=ind, end_ind=ind_e, error=error)
        CALL geometrical_center(pos=elem%pos(ind:ind_e+tmc_params%dim_per_elem-1), &
                                center=elem_center, error=error)
        IF(check_pos_in_subbox(pos=elem_center, &
             subbox_center=elem%subbox_center, &
             box_scale=elem%box_scale, tmc_params=tmc_params, &
             error=error)) &
          mol_in_sb(m) = status_ok
      END DO
      ! rotate the selected amount of molecules
      IF(ANY(mol_in_sb(:).EQ.status_ok)) THEN
        counter=1
        rot_molecule_loop:DO
          ! select molecule
          IF(tmc_params%nr_elem_mv.EQ.0) THEN
            m = counter
          ELSE
            rnd = next_random_number(rng_stream, error=error) ! next random number
            m = INT(rnd*nr_molec)+1
          END IF
          CALL get_mol_indeces(tmc_params=tmc_params, mol_arr=elem%mol, mol=m, &
                               start_ind=ind, end_ind=ind_e, error=error)
          ! when "molecule" is single atom, search a new one
          IF(ind .EQ. ind_e) CYCLE rot_molecule_loop

          ! apply move
          IF(mol_in_sb(m).EQ.status_ok)THEN
            CALL do_mol_rot(pos=elem%pos, ind_start=ind, ind_end=ind_e, &
                   max_angle=move_types%mv_size(&
                               mv_type_mol_rot,mv_conf), &
                   move_types=move_types, rng_stream=rng_stream, &
                   dim_per_elem=tmc_params%dim_per_elem, error=error)
            ! update sub box status of single atom
            DO i=ind, ind_e+tmc_params%dim_per_elem-1, tmc_params%dim_per_elem
              elem_center = elem%pos(i:i+tmc_params%dim_per_elem-1)
              IF(check_pos_in_subbox(pos=elem_center, &
                   subbox_center=elem%subbox_center, &
                   box_scale=elem%box_scale, tmc_params=tmc_params, &
                   error=error))THEN
                elem%elem_stat(i:i+tmc_params%dim_per_elem-1) = status_ok
              ELSE
                elem%elem_stat(i:i+tmc_params%dim_per_elem-1) = status_frozen
              END IF
            END DO
          ELSE
            ! element was not in sub box, search new one instead
            IF(tmc_params%nr_elem_mv.GT.0) counter = counter-1
          END IF
          counter = counter+1
          IF(counter.GT.act_nr_elem_mv) EXIT rot_molecule_loop
        END DO rot_molecule_loop
      END IF
      DEALLOCATE(elem_center)
      DEALLOCATE(mol_in_sb)

    !-- velocity changes for MD 
    !-- here all velocities are changed
    CASE(mv_type_MD)
      CPPreconditionNoFail(ASSOCIATED(tmc_params%atoms),cp_failure_level,routineP,error)
      change_all_velocities_loop: DO i=1, SIZE(elem%pos)
        !-- attention, move type size is in atomic units of velocity
        IF(elem%elem_stat(i).NE.status_frozen) THEN
          CALL vel_change(vel=elem%vel(i), &
                          atom_kind=tmc_params%atoms(INT(i/REAL(tmc_params%dim_per_elem,KIND=dp))+1),&
                          phi=move_types%mv_size(mv_type_MD,1),& ! TODO parallel tempering move sizes for vel_change
                          temp=tmc_params%Temp(mv_conf), &
                          rnd_sign_change=.TRUE., &! MD_vel_invert, &
                          rng_stream=rng_stream, error=error)
        END IF
      END DO change_all_velocities_loop

    !-- proton order and disorder
    !   a loop of molecules is build an in this loop proton acceptors become proton donators
    !   Therefor the molecules are rotated along the not involved O-H bond
    CASE(mv_type_proton_reorder)
      CALL search_and_do_proton_displace_loop(elem=elem, &
             short_loop=move_rejected, rng_stream=rng_stream, &
             tmc_params=tmc_params, error=error)

    !-- volume move
    ! the box is increased or decreased and with it the coordinates
    CASE(mv_type_volume_move)
      CALL change_volume(conf=elem, T_ind=mv_conf, move_types=move_types, &
                         rng_stream=rng_stream, tmc_params=tmc_params, &
                         mv_cen_of_mass=tmc_params%mv_cen_of_mass, error=error)

    !-- atom swap
    ! two atoms of different types are swapped
    CASE(mv_type_atom_swap)
      CALL swap_atoms(conf=elem, move_types=move_types, rng_stream=rng_stream,&
                         tmc_params=tmc_params, error=error)

    CASE DEFAULT
      CALL cp_assert(.FALSE.,&
                    cp_failure_level,cp_assertion_failed,routineP,&
                    "unknown move type "//&
                    cp_to_string(elem%move_type),error)
    END SELECT

    CALL get_rng_stream(rng_stream=rng_stream, bg=elem%rng_seed(:,:,1), &
           cg=elem%rng_seed(:,:,2), ig=elem%rng_seed(:,:,3), error=error)

  END SUBROUTINE change_pos

! *****************************************************************************
!> \brief gets the index of the first molecule element position and the size
!> \param tmc_params TMC parameters with dim_per_elem
!> \param mol_arr array with molecule information (which atom attend which mol)
!> \param mol the selected molecule number
!> \param start_ind start index of the first atom in molecule
!> \param end_ind index of the last atom in molecule
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 10.2013
! *****************************************************************************
  SUBROUTINE get_mol_indeces(tmc_params, mol_arr, mol, start_ind, end_ind, error)
    TYPE(tmc_param_type), POINTER            :: tmc_params
    INTEGER, DIMENSION(:), INTENT(IN), &
      POINTER                                :: mol_arr
    INTEGER, INTENT(IN)                      :: mol
    INTEGER, INTENT(OUT)                     :: start_ind, end_ind
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mol_indeces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure = .FALSE.

    start_ind = -1
    end_ind   = -1

    CPPrecondition(ASSOCIATED(mol_arr),cp_failure_level,routineP,error,failure)
    CPPrecondition(mol.LE.MAXVAL(mol_arr(:)),cp_failure_level,routineP,error,failure)
    ! get start index
    loop_start: DO i=1, SIZE(mol_arr)
      IF(mol_arr(i) .EQ. mol) THEN
        start_ind = i
        EXIT loop_start
      END IF
    END DO loop_start
    ! get end index
    loop_end: DO i=SIZE(mol_arr), i, -1
      IF(mol_arr(i) .EQ. mol) THEN
        end_ind = i
        EXIT loop_end
      END IF
    END DO loop_end
    ! check if all atoms inbetween attend to molecule
    CPPostcondition(ALL(mol_arr(start_ind:end_ind).EQ.mol),cp_failure_level,routineP,error,failure)
    CPPostcondition(start_ind.GT.0,cp_failure_level,routineP,error,failure)
    CPPostcondition(end_ind.GT.0,cp_failure_level,routineP,error,failure)
    ! convert to indeces mapped for the position array (multiple dim per atom)
    start_ind = (start_ind-1)*tmc_params%dim_per_elem+1
    end_ind   = (end_ind-1)*tmc_params%dim_per_elem+1
  END SUBROUTINE

! *****************************************************************************
!> \brief checks if a position is within the sub box
!>        returns true if position is inside
!> \param pos array with positions
!> \param subbox_center actual center of sub box
!> \param subbox_size size of the sub box
!> \param box_scale scaling factors for the cell
!> \param tmc_params TMC parameters with sub box size and cell
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  FUNCTION check_pos_in_subbox(pos, subbox_center, box_scale, tmc_params, error)&
                               RESULT(inside)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pos, subbox_center, box_scale
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: inside

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_pos_in_subbox', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, flag
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: pos_tmp

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(pos),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(subbox_center),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(box_scale),cp_failure_level,routineP,error,failure)
    ! if pressure is defined, no scale should be 0
    flag = .NOT.((tmc_params%pressure.GT.0.0_dp).AND.(ANY(box_scale.EQ.0.0_dp)))
    CPPrecondition(flag,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(pos).EQ.3,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(pos).EQ.SIZE(subbox_center),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    ALLOCATE(pos_tmp(SIZE(pos)))

    inside = .TRUE.
    ! return if no subbox is defined
    IF(.NOT.ANY(tmc_params%sub_box_size(:).LE.0.1_dp)) THEN
      pos_tmp = pos(:)-subbox_center(:)
      CALL get_scaled_cell(cell=tmc_params%cell, box_scale=box_scale, &
                           vec=pos_tmp, error=error)
      ! check
      IF(ANY(pos_tmp(:).GE.tmc_params%sub_box_size(:)/2.0) .OR. &
         ANY(pos_tmp(:).LE.-tmc_params%sub_box_size(:)/2.0)) THEN
        inside = .FALSE.
      END IF
    END IF
    DEALLOCATE(pos_tmp)
    ! end the timing
    CALL timestop(handle)
  END FUNCTION check_pos_in_subbox

! *****************************************************************************
!> \brief set a new random sub box center and counte the number of atoms in it
!> \param 
!> \param 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE elements_in_new_subbox(tmc_params, rng_stream, elem, &
                                    nr_of_sub_box_elements, error)
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(tree_type), POINTER                 :: elem
    INTEGER, INTENT(OUT)                     :: nr_of_sub_box_elements
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'elements_in_new_subbox', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rnd
    REAL(KIND=dp), DIMENSION(3)              :: box_size
    REAL(KIND=dp), DIMENSION(:), POINTER     :: atom_tmp, center_of_sub_box

    failure = .FALSE.
    NULLIFY(center_of_sub_box, atom_tmp)

    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    IF(ANY(tmc_params%sub_box_size(:).LE.0.1_dp)) THEN
      !CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
      !       "try to count elements in sub box without sub box.",&
      !       error,failure)
      elem%elem_stat = status_ok
      nr_of_sub_box_elements = SIZE(elem%elem_stat)
    ELSE
      ALLOCATE(center_of_sub_box(tmc_params%dim_per_elem))
      ALLOCATE(atom_tmp(tmc_params%dim_per_elem))
      nr_of_sub_box_elements = 0
      ! -- define the center of the sub box 
      CALL set_rng_stream(rng_stream=rng_stream, &
                          bg=elem%rng_seed(:,:,1), cg=elem%rng_seed(:,:,2), &
                          ig=elem%rng_seed(:,:,3), error=error)

      CALL get_cell(cell=tmc_params%cell, abc=box_size)
      DO i=1, SIZE(tmc_params%sub_box_size)
         rnd = next_random_number(rng_stream, error=error) ! next random number
         center_of_sub_box(i) = rnd * box_size(i)
      END DO
      elem%subbox_center(:) = center_of_sub_box(:)
  
      CALL get_rng_stream(rng_stream=rng_stream, &
                          bg=elem%rng_seed(:,:,1), cg=elem%rng_seed(:,:,2), &
                          ig=elem%rng_seed(:,:,3), error=error)
  
      ! check all elements if they are in subbox
      DO i=1, SIZE(elem%pos), tmc_params%dim_per_elem
         atom_tmp(:) = elem%pos(i:i+tmc_params%dim_per_elem-1)
         IF(check_pos_in_subbox(pos=atom_tmp, &
              subbox_center=center_of_sub_box, box_scale=elem%box_scale, &
              tmc_params=tmc_params, error=error)) THEN
           elem%elem_stat(i:i+tmc_params%dim_per_elem-1) = status_ok
           nr_of_sub_box_elements = nr_of_sub_box_elements + 1
         ELSE
           elem%elem_stat(i:i+tmc_params%dim_per_elem-1) = status_frozen
         END IF
      END DO
      DEALLOCATE(atom_tmp)
      DEALLOCATE(center_of_sub_box)
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE elements_in_new_subbox

! *****************************************************************************
!> \brief molecule rotation using quaternions
!> \param pos atom positions
!> \param ind_start starting index in the array
!> \param ind_end index of last atom in the array
!> \param max_angle maximal angle in each direction
!> \param rng_stream ramdon stream
!> \param dim_per_elem dimension per atom
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE do_mol_rot(pos, ind_start, ind_end, max_angle, move_types, &
                        rng_stream, dim_per_elem, error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pos
    INTEGER                                  :: ind_start, ind_end
    REAL(KIND=dp)                            :: max_angle
    TYPE(tmc_move_type), POINTER             :: move_types
    TYPE(rng_stream_type), POINTER           :: rng_stream
    INTEGER                                  :: dim_per_elem
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_mol_rot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a1, a2, a3, q0, q1, q2, q3, &
                                                rnd
    REAL(KIND=dp), DIMENSION(3, 3)           :: rot
    REAL(KIND=dp), DIMENSION(:), POINTER     :: elem_center

    failure = .FALSE.
    NULLIFY(elem_center)

    CPPrecondition(ASSOCIATED(pos),cp_failure_level,routineP,error,failure)
    CPPrecondition(dim_per_elem.EQ.3,cp_failure_level,routineP,error,failure)
    CPPrecondition(ind_start.GT.0.AND.ind_start.LT.SIZE(pos),cp_failure_level,routineP,error,failure)
    CPPrecondition(ind_end.GT.0.AND.ind_end.LT.SIZE(pos),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(move_types),cp_failure_level,routineP,error,failure)

    ! calculate rotation matrix (using quanternions)
    rnd = next_random_number(rng_stream, error=error) ! next random number
    a1 = (rnd-0.5)*2.0*max_angle !move_types%mv_size(mv_type_mol_rot,mv_conf)
    rnd = next_random_number(rng_stream, error=error) ! next random number
    a2 = (rnd-0.5)*2.0*max_angle !move_types%mv_size(mv_type_mol_rot,mv_conf)
    rnd = next_random_number(rng_stream, error=error) ! next random number
    a3 = (rnd-0.5)*2.0*max_angle !move_types%mv_size(mv_type_mol_rot,mv_conf)
    q0 = COS(a2/2)*COS((a1+a3)/2.0_dp)
    q1 = SIN(a2/2)*COS((a1-a3)/2.0_dp)
    q2 = SIN(a2/2)*SIN((a1-a3)/2.0_dp)
    q3 = COS(a2/2)*SIN((a1+a3)/2.0_dp)
    rot = RESHAPE((/q0*q0+q1*q1-q2*q2-q3*q3, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2),&
                    2*(q1*q2+q0*q3), q0*q0-q1*q1+q2*q2-q3*q3, 2*(q2*q3-q0*q1),&
                    2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1),  q0*q0-q1*q1-q2*q2+q3*q3/),(/3,3/))

    ALLOCATE(elem_center(dim_per_elem))
    ! calculate geometrical center
    CALL geometrical_center(pos=pos(ind_start:ind_end+dim_per_elem-1), &
                            center=elem_center, error=error)

    ! proceed rotation
    atom_loop:DO i=ind_start, ind_end+dim_per_elem-1, dim_per_elem
       pos(i:i+2) = MATMUL(pos(i:i+2)-elem_center(:), rot) + elem_center(:)
    END DO atom_loop
    DEALLOCATE(elem_center)
  END SUBROUTINE do_mol_rot
  
! *****************************************************************************
!> \brief velocity change should be gaussian distributed 
!>        around the old velocity with respect to kB*T/m
!> \param vel velocity of atom (one direction)
!> \param mass atom mass
!> \param phi angle for mixing old with random gaussian distributed velocity
!>        phi =90 degree -> only gaussian velocity around 0
!>        phi = 0 degree -> only old velocity (with sign change)
!> \param temp temperature for gaussian distributed velocity
!> \param rnd_sign_change if sign of old velocity should change randomly
!> \param rng_stream random number stream
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE vel_change(vel, atom_kind, phi, temp, rnd_sign_change, rng_stream, &
                        error)
    REAL(KIND=dp), INTENT(INOUT)             :: vel
    TYPE(tmc_atom_type)                      :: atom_kind
    REAL(KIND=dp), INTENT(IN)                :: phi, temp
    LOGICAL                                  :: rnd_sign_change
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'vel_change', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: d
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: delta_vel, kB, rnd1, rnd2, &
                                                rnd3, rnd_g

    failure = .FALSE.
    kB = boltzmann/joule

    CPPrecondition(ASSOCIATED(rng_stream),cp_failure_level,routineP,error,failure)

    !phi = move_types%mv_size(mv_type_MD,1) ! TODO parallel tempering move sizes for vel_change
    ! hence first producing a gaussian random number
    rnd1 = next_random_number(rng_stream, error=error)
    rnd2 = next_random_number(rng_stream, error=error)

    rnd_g = SQRT(-2.0_dp*LOG(rnd1))*COS(2.0_dp*PI*rnd2)
    !we can also produce a second one in the same step:
    !rnd_g2 = SQRT(-2.0_dp*LOG(rnd1))*SIN(2.0_dp*PI*rnd2)

    ! adapting the variance with respect to kB*T/m
    delta_vel = SQRT(kB*temp/atom_kind%mass)*rnd_g
    ! check if TODO random velocity sign change
      ! using detailed balance, velocity sign changes are necessary, 
      ! which are done randomly and
      ! can be switched of using MD_vel_invert
      ! without still the balance condition should be fullfilled

    rnd3 = next_random_number(rng_stream, error=error)
    IF(rnd3.GE.0.5.AND.rnd_sign_change) THEN
       d = -1
    ELSE
       d = 1
    END IF
    vel = SIN(phi)*delta_vel + COS(phi)*vel*d*1.0_dp
  END SUBROUTINE vel_change

! *****************************************************************************
!> \brief proton order and disorder (specialized move for ice Ih)
!>        a loop of molecules is build an 
!>        in this loop proton acceptors become proton donators
!>        Therefor the molecules are rotated along the not involved O-H bond
!> \param elem sub tree element with actual positions
!> \param short_loop return if the a loop shorter than 6 molecules is found
!>        (should not be in ice structure)
!> \param rng_stream random number stream
!> \param tmc_params TMC parameters with numbers of dimensions per element
!>        number of atoms per molecule
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE search_and_do_proton_displace_loop(elem, short_loop, rng_stream, &
                                                tmc_params, error)
    TYPE(tree_type), POINTER                 :: elem
    LOGICAL                                  :: short_loop
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'search_and_do_proton_displace_loop', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1000)                      :: tmp_chr
    INTEGER                                  :: counter, donor_acceptor, &
                                                handle, k, mol, nr_mol
    INTEGER, DIMENSION(:), POINTER           :: mol_arr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rnd

    failure = .FALSE.
    NULLIFY(mol_arr)

    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rng_stream),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    short_loop = .FALSE.
    counter = 0
    nr_mol  = MAXVAL(elem%mol(:))
    ! ind_arr: one array element for each molecule
    ALLOCATE(mol_arr(nr_mol))
    mol_arr(:) = -1
    donor_acceptor = not_selected
    ! select randomly if neighboring molecule is donor / acceptor
    IF(next_random_number(rng_stream, error=error).LT.0.5_dp) THEN
      donor_acceptor = proton_acceptor
    ELSE
      donor_acceptor = proton_donor
    END IF

    ! first step build loop
    !  select randomly one atom
    rnd = next_random_number(rng_stream, error=error)
    ! the randomly selected first atom
    mol = INT(rnd*nr_mol) +1
    counter = counter +1
    mol_arr(counter) = mol

    ! do until the loop is closed 
    !  (until path connects back to any spot of the path)
    chain_completition_loop: DO
       counter = counter+1
       ! find nearest neighbor
       !  (with same state, in the chain, proton donator or proton accptor)
       CALL find_nearest_proton_acceptor_donator(elem=elem, mol=mol, &
              donor_acceptor=donor_acceptor, tmc_params=tmc_params, &
              rng_stream=rng_stream, error=error)
       IF(ANY(mol_arr(:).EQ.mol))&
         EXIT chain_completition_loop
       mol_arr(counter) = mol
    END DO chain_completition_loop
    counter = counter -1 ! last searched element is equal to one other in list

    ! just take the loop of molecules out of the chain
    DO k = 1,counter
      IF(mol_arr(k).EQ.mol)&
        EXIT
    END DO
    mol_arr(1:counter-k+1) = mol_arr(k:counter)
    counter = counter-k+1

    ! check if loop is minimum size of 6 molecules
    CALL cp_assert(counter.GE.6,cp_warning_level,cp_assertion_failed,&
                 routineP,"short proton loop with"//cp_to_string(counter)//&
                 "molecules.", error,failure)
    IF(failure) THEN
      tmp_chr=""
      WRITE(tmp_chr, *) mol_arr(1:counter)
      CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
                 routineP,"selected molecules:"//TRIM(tmp_chr),&
                 error,failure)
      short_loop = .TRUE.
    END IF
    ! rotate the molecule along the not involved O-H bond 
    !   (about the angle in of the neighboring chain elements)
    CALL rotate_molecules_in_chain(tmc_params=tmc_params, elem=elem, &
           mol_arr_in=mol_arr(1:counter), donor_acceptor=donor_acceptor, &
           error=error)
    DEALLOCATE(mol_arr)

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE search_and_do_proton_displace_loop

! *****************************************************************************
!> \brief searches the next (first atom of) neighboring molecule 
!>        which is proton donor / acceptor
!> \param elem sub tree element with actual positions 
!> \param mol (in_out) actual regarded molecule, which neighbor is searched for
!> \param donor_acceptor type of searched neighbor 
!>        (proton donor or proton acceptor)
!> \param tmc_params TMC parameters with numbers of dimensions per element
!>        number of atoms per molecule
!> \param rng_stream random number stream
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE find_nearest_proton_acceptor_donator(elem, mol, donor_acceptor, &
                                               tmc_params, rng_stream, error)
    TYPE(tree_type), POINTER                 :: elem
    INTEGER                                  :: mol, donor_acceptor
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'find_nearest_proton_acceptor_donator', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ind, ind_e, ind_n, &
                                                mol_tmp, nr_mol
    INTEGER, DIMENSION(2)                    :: neighbor_mol
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dist_tmp, rnd
    REAL(KIND=dp), DIMENSION(:), POINTER     :: distH1, distH2, distO

    failure = .FALSE.
    NULLIFY(distO, distH1, distH2)
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    nr_mol = MAXVAL(elem%mol)
    ALLOCATE(distO(nr_mol))
    ALLOCATE(distH1(nr_mol))
    ALLOCATE(distH2(nr_mol))
    !-- initialize the distances to huge values
    ! distance of nearest proton of certain molecule to preselected O
    distO(:)  = HUGE(distO(1))
    ! distance of (first) proton of preselected molecule to certain molecule
    distH1(:) = HUGE(distH1(1))
    ! distance of (second) proton of preselected molecule to certain molecule
    distH2(:) = HUGE(distH2(1))

    ! get the indices of the old O atom (assuming the first atom of the molecule the first atom)
    CALL get_mol_indeces(tmc_params=tmc_params, mol_arr=elem%mol, mol=mol, &
                               start_ind=ind, end_ind=ind_e, error=error)

    ! calculate distances to all molecules
    list_distances: DO mol_tmp=1, nr_mol
      IF(mol_tmp.EQ.mol) CYCLE list_distances
      ! index of the molecule (the O atom)
      ! assume the first atom of the molecule the first atom
      CALL get_mol_indeces(tmc_params=tmc_params, mol_arr=elem%mol, &
             mol=mol_tmp, start_ind=ind_n, end_ind=ind_e, error=error)
      ! check if selected molecule is water respectively consists of 3 atoms
      IF(MOD(ind_e-ind_n,3).GT.0) THEN
        CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
                 routineP,"selected a molecule with more than 3 atoms, "//&
                 "the proton reordering does not support, skip molecule",&
                 error,failure)
        CYCLE list_distances
      END IF
      IF(donor_acceptor.EQ.proton_acceptor) THEN
        IF(check_donor_acceptor(elem=elem, i_orig=ind, i_neighbor=ind_n, &
             tmc_params=tmc_params, error=error).EQ.proton_acceptor) THEN
          !distance of fist proton to certain O
          distH1(mol_tmp) = nearest_distance(&
                          x1=elem%pos(ind+tmc_params%dim_per_elem:&
                                      ind+2*tmc_params%dim_per_elem-1), &
                          x2=elem%pos(ind_n:ind_n+tmc_params%dim_per_elem-1),&
                          cell=tmc_params%cell, box_scale=elem%box_scale, error=error)
          !distance of second proton to certain O
          distH2(mol_tmp) = nearest_distance(&
                          x1=elem%pos(ind+2*tmc_params%dim_per_elem:&
                                      ind+3*tmc_params%dim_per_elem-1), &
                          x2=elem%pos(ind_n:ind_n+tmc_params%dim_per_elem-1), &
                          cell=tmc_params%cell, box_scale=elem%box_scale, error=error)
        END IF
      END IF
      !check for neighboring proton donors
      IF(donor_acceptor.EQ.proton_donor) THEN
        IF(check_donor_acceptor(elem=elem, i_orig=ind, i_neighbor=ind_n, &
             tmc_params=tmc_params, error=error).EQ.proton_donor) THEN
          !distance of selected O to all first protons of other melecules
          distO(mol_tmp) = nearest_distance(&
                         x1=elem%pos(ind:ind+tmc_params%dim_per_elem-1), &
                         x2=elem%pos(ind_n+tmc_params%dim_per_elem:&
                                     ind_n+2*tmc_params%dim_per_elem-1), &
                         cell=tmc_params%cell, box_scale=elem%box_scale, error=error)
          dist_tmp   = nearest_distance(&
                         x1=elem%pos(ind:ind+tmc_params%dim_per_elem-1), &
                         x2=elem%pos(ind_n+2*tmc_params%dim_per_elem:&
                                     ind_n+3*tmc_params%dim_per_elem-1), &
                         cell=tmc_params%cell, box_scale=elem%box_scale, error=error)
          IF(dist_tmp.LT.distO(mol_tmp)) distO(mol_tmp) = dist_tmp
         END IF
      END IF
    END DO list_distances

    mol_tmp = 1
    ! select the nearest neighbors
    !check for neighboring proton acceptors
    IF(donor_acceptor.EQ.proton_acceptor) THEN
      neighbor_mol(mol_tmp) = MINLOC(distH1(:),1)
      neighbor_mol(mol_tmp+1) = MINLOC(distH2(:),1)
      ! if both smallest distances points to the shortest molecule search also the second next shortest distance
      IF(neighbor_mol(mol_tmp).EQ.neighbor_mol(mol_tmp+1))THEN
        distH1(neighbor_mol(mol_tmp)) = HUGE(distH1(1))
        distH2(neighbor_mol(mol_tmp+1)) = HUGE(distH2(1))
        IF(MINVAL(distH1(:),1).LT.MINVAL(distH2(:),1)) THEN
           neighbor_mol(mol_tmp) = MINLOC(distH1(:),1)
        ELSE
           neighbor_mol(mol_tmp+1) = MINLOC(distH2(:),1)
        END IF
      END IF
      mol_tmp = mol_tmp +2
    END IF

    !check for neighboring proton donors
    IF(donor_acceptor.EQ.proton_donor) THEN
      neighbor_mol(mol_tmp) = MINLOC(distO(:),1)
      distO(neighbor_mol(mol_tmp)) = HUGE(distO(1))
      neighbor_mol(mol_tmp+1) = MINLOC(distO(:),1)
    END IF

    ! select randomly the next neighboring molecule
    rnd = next_random_number(rng_stream, error=error)
    ! the randomly selected atom: return value!
    mol_tmp = neighbor_mol(INT(rnd*SIZE(neighbor_mol(:)))+1)
    mol = mol_tmp

    DEALLOCATE(distO)
    DEALLOCATE(distH1)
    DEALLOCATE(distH2)

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE find_nearest_proton_acceptor_donator

! *****************************************************************************
!> \brief checks if neighbor of the selected/orig element 
!>        is a proron donator or acceptor
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  FUNCTION check_donor_acceptor(elem, i_orig, i_neighbor, tmc_params, error) &
                                RESULT(donor_acceptor)
    TYPE(tree_type), POINTER                 :: elem
    INTEGER                                  :: i_orig, i_neighbor
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: donor_acceptor

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_donor_acceptor', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(4)              :: distances

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(i_orig.GE.1.AND.i_orig.LE.SIZE(elem%pos),cp_failure_level,routineP,error,failure)
    CPPrecondition(i_neighbor.GE.1.AND.i_neighbor.LE.SIZE(elem%pos),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    ! 1. proton of orig with neighbor O
    distances(1) = nearest_distance(&
                     x1=elem%pos(i_neighbor:i_neighbor+tmc_params%dim_per_elem-1), &
                     x2=elem%pos(i_orig+tmc_params%dim_per_elem:&
                                 i_orig+2*tmc_params%dim_per_elem-1), &
                     cell=tmc_params%cell, box_scale=elem%box_scale, error=error)
    ! 2. proton of orig with neighbor O
    distances(2) = nearest_distance(&
                     x1=elem%pos(i_neighbor:i_neighbor+tmc_params%dim_per_elem-1), &
                     x2=elem%pos(i_orig+2*tmc_params%dim_per_elem:&
                                 i_orig+3*tmc_params%dim_per_elem-1), &
                     cell=tmc_params%cell, box_scale=elem%box_scale, error=error)
    ! 1. proton of neighbor with orig O
    distances(3) = nearest_distance(&
                     x1=elem%pos(i_orig:i_orig+tmc_params%dim_per_elem-1), &
                     x2=elem%pos(i_neighbor+tmc_params%dim_per_elem:&
                                 i_neighbor+2*tmc_params%dim_per_elem-1), &
                     cell=tmc_params%cell, box_scale=elem%box_scale, error=error)
    ! 2. proton of neigbor with orig O
    distances(4) = nearest_distance(&
                     x1=elem%pos(i_orig:i_orig+tmc_params%dim_per_elem-1), &
                     x2=elem%pos(i_neighbor+2*tmc_params%dim_per_elem:&
                                 i_neighbor+3*tmc_params%dim_per_elem-1), &
                     cell=tmc_params%cell, box_scale=elem%box_scale, error=error)

    IF(MINLOC(distances(:),1).LE.2) THEN
       donor_acceptor = proton_acceptor
    ELSE
       donor_acceptor = proton_donor
    END IF
  END FUNCTION check_donor_acceptor

! *****************************************************************************
!> \brief rotates all the molecules in the chain
!>        the protons were flipped from the donor to the acceptor
!> \param tmc_params TMC environment parameters
!> \param elem sub tree element the pos of the molecules in chain should be 
!>        changed by rotating
!> \param mol_arr_in array of indeces of molecules, should be rotated
!> \param donor_acceptor gives the direction of rotation
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE rotate_molecules_in_chain(tmc_params, elem, mol_arr_in, &
                                       donor_acceptor, error)
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(tree_type), POINTER                 :: elem
    INTEGER, DIMENSION(:)                    :: mol_arr_in
    INTEGER                                  :: donor_acceptor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rotate_molecules_in_chain', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: H_offset, handle, i, ind
    INTEGER, DIMENSION(:), POINTER           :: ind_arr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dihe_angle, dist_near, tmp
    REAL(KIND=dp), DIMENSION(3)              :: rot_axis, tmp_1, tmp_2, &
                                                vec_1O, vec_2H_f, vec_2H_m, &
                                                vec_2O, vec_3O, vec_4O, &
                                                vec_rotated
    TYPE(cell_type), POINTER                 :: tmp_cell

    failure = .FALSE.
    NULLIFY(ind_arr, tmp_cell)

    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    ALLOCATE(ind_arr(0:SIZE(mol_arr_in)+1))
    DO i=1, SIZE(mol_arr_in)
      CALL get_mol_indeces(tmc_params=tmc_params, mol_arr=elem%mol, &
                           mol=mol_arr_in(i), &
                           start_ind=ind_arr(i), end_ind=ind, error=error)
    END DO
    ind_arr(0) = ind_arr(SIZE(ind_arr)-2)
    ind_arr(SIZE(ind_arr)-1) = ind_arr(1)

    ! get the scaled cell
    ALLOCATE(tmp_cell)
    CALL get_scaled_cell(cell=tmc_params%cell, box_scale=elem%box_scale, &
           scaled_cell=tmp_cell, error=error)

    ! rotate single molecules
    DO i=1, SIZE(ind_arr)-2
      ! the 3 O atoms
      vec_1O(:) = elem%pos(ind_arr(i-1):ind_arr(i-1)+tmc_params%dim_per_elem-1)
      vec_2O(:) = elem%pos(ind_arr( i ):ind_arr( i )+tmc_params%dim_per_elem-1)
      vec_3O(:) = elem%pos(ind_arr(i+1):ind_arr(i+1)+tmc_params%dim_per_elem-1)
      ! the H atoms 
      ! distinguished between the one fixed (rotation axis with 2 O) 
      ! and the moved one
      ! if true the first H atom is between the O atoms
      IF(nearest_distance(&
           x1=elem%pos(ind_arr(i+donor_acceptor):&
                       ind_arr(i+donor_acceptor)+tmc_params%dim_per_elem-1),&
           x2=elem%pos(ind_arr(i)+tmc_params%dim_per_elem:&
                       ind_arr(i)+2*tmc_params%dim_per_elem-1), &
           cell=tmc_params%cell, box_scale=elem%box_scale, &
           error=error) .LT. &
         nearest_distance(&
           x1=elem%pos(ind_arr(i+donor_acceptor):&
                       ind_arr(i+donor_acceptor)+tmc_params%dim_per_elem-1),&
           x2=elem%pos(ind_arr(i)+2*tmc_params%dim_per_elem:&
                       ind_arr(i)+3*tmc_params%dim_per_elem-1), &
           cell=tmc_params%cell, box_scale=elem%box_scale, &
           error=error)) THEN
        vec_2H_m = elem%pos(ind_arr(i)+tmc_params%dim_per_elem:&
                            ind_arr(i)+2*tmc_params%dim_per_elem-1)
        vec_2H_f = elem%pos(ind_arr(i)+2*tmc_params%dim_per_elem:&
                            ind_arr(i)+3*tmc_params%dim_per_elem-1)
        H_offset = 1
      ELSE
        vec_2H_f = elem%pos(ind_arr(i)+tmc_params%dim_per_elem:&
                            ind_arr(i)+2*tmc_params%dim_per_elem-1)
        vec_2H_m = elem%pos(ind_arr(i)+2*tmc_params%dim_per_elem:&
                            ind_arr(i)+3*tmc_params%dim_per_elem-1)
        H_offset = 2
      END IF

      IF (.TRUE.) THEN !TODO find a better switch for the pauling model
        
        ! do rotation (NOT pauling model)
        tmp_1 = pbc(vec_2O-vec_1O, tmp_cell)
        tmp_2 = pbc(vec_3O-vec_2H_f, tmp_cell)

        dihe_angle = donor_acceptor*dihedral_angle(tmp_1,vec_2H_f-vec_2O,tmp_2)
        DO ind=ind_arr(i), ind_arr(i)+tmc_params%dim_per_elem*3-1, tmc_params%dim_per_elem
          ! set rotation vector
          !vec_rotated = rotate_vector(vec_2H_m-vec_2O, dihe_angle, vec_2H_f-vec_2O)
          vec_rotated = rotate_vector(elem%pos(ind:&
                            ind+tmc_params%dim_per_elem-1)-vec_2O, &
                          dihe_angle, vec_2H_f-vec_2O)

          ! set new position
          !elem%pos(ind_arr(i)+H_offset*dim_per_elem:ind_arr(i)+(H_offset+1)*dim_per_elem-1) = vec_2O+vec_rotated
          elem%pos(ind:ind+tmc_params%dim_per_elem-1) = vec_2O+vec_rotated
        END DO
      ELSE
        ! using the pauling model 
        !  (see Aragones and Vega: Dielectric constant of ices...)
        ! the rotation axis is defined using the 4th not involved O 
        !  (next to the not involved H)
        ! O atom next to not involved proton for axis calculation
        dist_near = HUGE(dist_near)
        search_O_loop: DO ind=1, SIZE(elem%pos), &
                       tmc_params%dim_per_elem*3
          IF(ind.EQ.ind_arr(i)) CYCLE search_O_loop
          tmp = nearest_distance(x1=vec_2H_f, &
                  x2=elem%pos(ind:ind+tmc_params%dim_per_elem-1), &
                  cell=tmc_params%cell, box_scale=elem%box_scale, &
                  error=error)
          IF(dist_near .GT. tmp) THEN
             dist_near = tmp
             vec_4O = elem%pos(ind:ind+tmc_params%dim_per_elem-1)
          END IF
        END DO search_O_loop
        rot_axis = pbc(-vec_2O(:)+vec_4O(:), tmp_cell)
        tmp_1 = pbc(vec_2O-vec_1O, tmp_cell)
        tmp_2 = pbc(vec_3O-vec_4O, tmp_cell)
        dihe_angle = donor_acceptor*dihedral_angle(tmp_1, rot_axis, tmp_2)
        vec_rotated = rotate_vector(vec_2H_m-vec_2O , dihe_angle, rot_axis)
        ! set new position
        elem%pos(ind_arr(i)+H_offset*tmc_params%dim_per_elem:&
                 ind_arr(i)+(H_offset+1)*tmc_params%dim_per_elem-1) &
          = vec_2O+vec_rotated
        vec_rotated = rotate_vector(vec_2H_f-vec_2O , dihe_angle, rot_axis)
        IF(H_offset.EQ.1) THEN
          H_offset = 2
        ELSE
          H_offset = 1
        END IF
        elem%pos(ind_arr(i)+H_offset*tmc_params%dim_per_elem:&
                 ind_arr(i)+(H_offset+1)*tmc_params%dim_per_elem-1) &
          = vec_2O+vec_rotated
      END IF
    END DO
    DEALLOCATE(tmp_cell)
    DEALLOCATE(ind_arr)
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE rotate_molecules_in_chain

! *****************************************************************************
!> \brief volume move, the box size is increased or decreased, 
!>        using the mv_size a the factor.
!>        the coordinated are scaled moleculewise 
!>        (the is moved like the center of mass is moves)
!> \param conf configuration to change with positions
!> \param T_ind temperature index, to select the correct temperature 
!>        for move size
!> \param rng_stream random number generator stream
!> \param tmc_params TMC parameters with e.g. dimensions of atoms and molecules
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE change_volume(conf, T_ind, move_types, rng_stream, tmc_params, &
                           mv_cen_of_mass, error)
    TYPE(tree_type), POINTER                 :: conf
    INTEGER                                  :: T_ind
    TYPE(tmc_move_type), POINTER             :: move_types
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(tmc_param_type), POINTER            :: tmc_params
    LOGICAL                                  :: mv_cen_of_mass
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'change_volume', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, dir, handle, ind, &
                                                ind_e, mol
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rnd, vol
    REAL(KIND=dp), DIMENSION(3)              :: box_length_new, &
                                                box_length_orig, box_scale_old
    REAL(KIND=dp), DIMENSION(:), POINTER     :: disp, scaling

    failure = .FALSE.
    NULLIFY(scaling, disp)

    CPPrecondition(ASSOCIATED(conf),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(move_types),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rng_stream),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(T_ind.GT.0.AND.T_ind.LE.tmc_params%nr_temp,cp_failure_level,routineP,error,failure)
    CPPrecondition(tmc_params%dim_per_elem.EQ.3,cp_failure_level,routineP,error,failure)
    CPPrecondition(tmc_params%cell%orthorhombic,cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    ALLOCATE(scaling(tmc_params%dim_per_elem))
    ALLOCATE(disp(tmc_params%dim_per_elem))

    box_scale_old(:) = conf%box_scale
    ! get the cell vector length of the configuration (before move)
    CALL get_scaled_cell(cell=tmc_params%cell, box_scale=conf%box_scale, &
                         abc=box_length_new, error=error)

    IF(.FALSE.)THEN
      ! the volume move in volume space (dV)
      IF(tmc_params%v_isotropic) THEN
        CALL get_scaled_cell(cell=tmc_params%cell, box_scale=conf%box_scale, &
        abc=box_length_new, vol=vol, error=error)
        rnd = next_random_number(rng_stream, error=error)
        vol=vol+(rnd-0.5_dp)*2.0_dp*move_types%mv_size(mv_type_volume_move, T_ind)
        box_length_new(:) = vol**(1/REAL(3,KIND=dp))
      ELSE
        CALL get_scaled_cell(cell=tmc_params%cell, box_scale=conf%box_scale, &
        abc=box_length_new, vol=vol, error=error)
        rnd = next_random_number(rng_stream, error=error)
        vol=vol+(rnd-0.5_dp)*2.0_dp*move_types%mv_size(mv_type_volume_move, T_ind)
        rnd = next_random_number(rng_stream, error=error)
        dir = 1+INT(rnd*3)
        box_length_new(dir) = 1.0_dp
        box_length_new(dir) = vol/PRODUCT(box_length_new(:))
      END IF
    ELSE
      ! the volume move in box length space (dL)
      ! increase / decrease box lenght in this direction
      ! l_n = l_o +- rnd * mv_size
      IF(tmc_params%v_isotropic) THEN
        rnd = next_random_number(rng_stream, error=error)
        box_length_new(:) = box_length_new(:) + &
          (rnd-0.5_dp)*2.0_dp*&
          move_types%mv_size(mv_type_volume_move, T_ind)
      ELSE
        ! select a random direction
        rnd = next_random_number(rng_stream, error=error)
        dir = 1+INT(rnd*3)
        rnd = next_random_number(rng_stream, error=error)
        box_length_new(dir) = box_length_new(dir) + &
          (rnd-0.5_dp)*2.0_dp*&
          move_types%mv_size(mv_type_volume_move, T_ind)
      END IF
    END IF

    ! get the original box length 
    scaling(:) = 1.0_dp
    CALL get_scaled_cell(cell=tmc_params%cell, &
                         box_scale=scaling,&
                         abc=box_length_orig, error=error)
    ! get the new box scale
    conf%box_scale(:) = box_length_new(:)/box_length_orig(:)
    ! molecule scaling
    scaling(:) = conf%box_scale(:)/box_scale_old(:)


    IF(mv_cen_of_mass .EQV. .FALSE.) THEN
       ! homogene scaling of atomic coordinates
       DO atom=1, SIZE(conf%pos), tmc_params%dim_per_elem
          conf%pos(atom:atom+tmc_params%dim_per_elem-1) = &
            conf%pos(atom:atom+tmc_params%dim_per_elem-1)*scaling(:)
       END DO
    ELSE
       DO mol=1, MAXVAL(conf%mol(:))
         ! move the molecule related to the molecule center of mass
         ! get center of mass
         CPPreconditionNoFail(ASSOCIATED(tmc_params%atoms),cp_failure_level,routineP,error)

         CALL get_mol_indeces(tmc_params=tmc_params, mol_arr=conf%mol, mol=mol, &
                              start_ind=ind, end_ind=ind_e, error=error)
         CALL center_of_mass(&
                    pos=conf%pos(ind:ind_e+tmc_params%dim_per_elem-1),&
                    atoms=tmc_params%atoms(INT(ind/REAL(tmc_params%dim_per_elem,KIND=dp))+1:&
                            INT(ind_e/REAL(tmc_params%dim_per_elem,KIND=dp))+1),&
                    center=disp, error=error)
         ! calculate the center of mass DISPLACEMENT
         disp(:) = disp(:) * (scaling(:)-1.0_dp)
         ! displace all atoms of the molecule
         DO atom=ind, ind_e+tmc_params%dim_per_elem-1, tmc_params%dim_per_elem
                conf%pos(atom:atom+tmc_params%dim_per_elem-1) = &
                  conf%pos(atom:atom+tmc_params%dim_per_elem-1) + disp(:)
             END DO
       END DO
    END IF

    DEALLOCATE(scaling)
    DEALLOCATE(disp)

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE change_volume

! *****************************************************************************
!> \brief volume move, two atoms of different types are swapped, both selected 
!>        randomly
!> \param conf configuration to change with positions
!> \param rng_stream random number generator stream
!> \param tmc_params TMC parameters with e.g. dimensions of atoms and molecules
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE swap_atoms(conf, move_types, rng_stream, tmc_params, error)
    TYPE(tree_type), POINTER                 :: conf
    TYPE(tmc_move_type), POINTER             :: move_types
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'swap_atoms', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_1, a_2, ind_1, ind_2
    LOGICAL                                  :: failure, found
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: pos_tmp

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(conf),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(move_types),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rng_stream),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params%atoms),cp_failure_level,routineP,error,failure)

    ! loop until two different atoms are found
    atom_search_loop: DO 
      ! select one atom randomly
      a_1 = INT(SIZE(conf%pos)/REAL(tmc_params%dim_per_elem,KIND=dp)* &
                next_random_number(rng_stream, error=error))+1
      ! select the second atom randomly
      a_2 = INT(SIZE(conf%pos)/REAL(tmc_params%dim_per_elem,KIND=dp)* &
                next_random_number(rng_stream, error=error))+1
      ! check if they have different kinds
      IF(tmc_params%atoms(a_1)%name .NE. tmc_params%atoms(a_2)%name) THEN
        ! if present, check if atoms have different type related to the specified table
        IF(ASSOCIATED(move_types%atom_lists)) THEN
          DO ind_1=1, SIZE(move_types%atom_lists)
            IF(ANY(move_types%atom_lists(ind_1)%atoms(:) .EQ. &
                   tmc_params%atoms(a_1)%name) .AND. &
               ANY(move_types%atom_lists(ind_1)%atoms(:) .EQ. &  
                   tmc_params%atoms(a_2)%name)) THEN
              found = .TRUE.
              EXIT atom_search_loop
            END IF
          END DO
        ELSE
          found = .TRUE.
          EXIT atom_search_loop
        END IF
      END IF
    END DO atom_search_loop
    IF(found) THEN
      ! perform coordinate exchange
      ALLOCATE(pos_tmp(tmc_params%dim_per_elem))
      ind_1 = (a_1-1)*tmc_params%dim_per_elem+1
      pos_tmp(:) = conf%pos(ind_1:ind_1+tmc_params%dim_per_elem-1)
      ind_2 = (a_2-1)*tmc_params%dim_per_elem+1
      conf%pos(ind_1:ind_1+tmc_params%dim_per_elem-1) = &
        conf%pos(ind_2:ind_2+tmc_params%dim_per_elem-1)
      conf%pos(ind_2:ind_2+tmc_params%dim_per_elem-1) = pos_tmp(:)
      DEALLOCATE(pos_tmp)
    END IF
  END SUBROUTINE swap_atoms

END MODULE tmc_moves
