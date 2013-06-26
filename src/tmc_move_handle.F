!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief acceptance retio handling of the different Monte Carlo Moves types 
!>        For each move type and each temperature average acceptence is 
!>        determined. 
!> \par History
!>      11.2012 created [Mandes Schönherr]
!> \author Mandes
! *****************************************************************************

MODULE tmc_move_handle
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE physcon,                         ONLY: au2a => angstrom
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_move_types,                  ONLY: &
       move_types_create, move_types_release, mv_type_MD, mv_type_NMC_moves, &
       mv_type_atom_trans, mv_type_gausian_adapt, mv_type_mol_rot, &
       mv_type_mol_trans, mv_type_none, mv_type_proton_reorder, &
       mv_type_swap_conf, mv_type_volume_move, nr_nmc_mv_types, tmc_move_type
  USE tmc_stati,                       ONLY: task_type_MC,&
                                             task_type_gaussian_adaptation,&
                                             task_type_ideal_gas
  USE tmc_tree_types,                  ONLY: global_tree_type,&
                                             status_accepted_result,&
                                             status_rejected_result,&
                                             tree_type
  USE tmc_types,                       ONLY: tmc_param_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_move_handle'

  PUBLIC :: finalize_mv_types, init_move_types, print_move_types
  PUBLIC :: select_random_move_type
  PUBLIC :: prob_update, add_mv_prob
  PUBLIC :: clear_move_probs
  PUBLIC :: select_random_NMC_move_type

CONTAINS

! *****************************************************************************
!> \brief initialization of the different moves, with sizes and probabilities
!> \param move_types pointer to the structure which should be deallocated
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE init_move_types(tmc_params, error)
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_move_types', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ind
    LOGICAL                                  :: failure, flag
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: probs

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tmc_params%moves_init_sizes),cp_failure_level,routineP,error,failure)
    flag = (SIZE(tmc_params%moves_init_sizes%mv_types).EQ.&
            SIZE(tmc_params%moves_init_sizes%mv_sizes))
    CPPrecondition( flag , cp_failure_level,routineP,error,failure)
    flag = (SIZE(tmc_params%moves_init_sizes%mv_probs).EQ.1) .OR. &
            SIZE(tmc_params%moves_init_sizes%mv_probs).EQ.&
            SIZE(tmc_params%moves_init_sizes%mv_types)
    CPPrecondition( flag, cp_failure_level,routineP,error,failure)

    ALLOCATE(probs(SIZE(tmc_params%moves_init_sizes%mv_types)))
    IF((SIZE(tmc_params%moves_init_sizes%mv_probs).EQ.1).AND.&
       (tmc_params%moves_init_sizes%mv_probs(1).LT.0.0_dp))THEN
      probs(:) = 1.0_dp
    ELSE
      probs(:) = tmc_params%moves_init_sizes%mv_probs
    END IF

    CALL move_types_create(tmc_params%move_types, tmc_params%nr_temp, error)

    tmc_params%move_types%acc_prob(0,:) = &
      tmc_params%moves_init_sizes%start_acc_prob
    ! mv_weights are used within NMC, special moves like NMC, 
    !  swap and volume moves are treated seperately.
    ! the probability to do a swap move is defined by the move size of this type.
    DO ind=1, SIZE(tmc_params%moves_init_sizes%mv_types)
      SELECT CASE(tmc_params%moves_init_sizes%mv_types(ind))
      ! molecule rotation
      CASE(mv_type_mol_rot)
         tmc_params%move_types%mv_weight(tmc_params%moves_init_sizes%mv_types(ind)) = probs(ind)/SUM(probs(:))
         tmc_params%move_types%mv_size(tmc_params%moves_init_sizes%mv_types(ind),1) = tmc_params%moves_init_sizes%mv_sizes(ind)
         ! convert units
         SELECT CASE(tmc_params%task_type)
         CASE(task_type_MC, task_type_ideal_gas)
            tmc_params%move_types%mv_size(tmc_params%moves_init_sizes%mv_types(ind),1) = &
                          tmc_params%move_types%mv_size(tmc_params%moves_init_sizes%mv_types(ind),1) *PI/180.0_dp
         CASE DEFAULT
           CALL cp_assert(.FALSE.,&
                    cp_failure_level,cp_assertion_failed,routineP,&
                    "move type mol_rot is not defined for this task type",&
                    error)
         END SELECT
         tmc_params%move_types%acc_prob(tmc_params%moves_init_sizes%mv_types(ind),:) = &
                                             tmc_params%moves_init_sizes%start_acc_prob
      ! atom / molecule translation
      CASE(mv_type_atom_trans, mv_type_mol_trans)
         tmc_params%move_types%mv_weight(tmc_params%moves_init_sizes%mv_types(ind)) = probs(ind)/SUM(probs(:))
         tmc_params%move_types%mv_size(tmc_params%moves_init_sizes%mv_types(ind),1) = tmc_params%moves_init_sizes%mv_sizes(ind)
         ! convert units
         SELECT CASE(tmc_params%task_type)
         CASE(task_type_MC, task_type_ideal_gas)
            tmc_params%move_types%mv_size(tmc_params%moves_init_sizes%mv_types(ind),1) = &
                                          tmc_params%move_types%mv_size(tmc_params%moves_init_sizes%mv_types(ind),1)/au2a
         CASE(task_type_gaussian_adaptation)
            !nothing to do (no unit conversion)
         CASE DEFAULT
           CALL cp_assert(.FALSE.,&
                    cp_failure_level,cp_assertion_failed,routineP,&
                    "move type atom / mol trans is not defined for this task type",&
                    error)
         END SELECT
         tmc_params%move_types%acc_prob(tmc_params%moves_init_sizes%mv_types(ind),:) = &
                                             tmc_params%moves_init_sizes%start_acc_prob
      ! Hybrid MC (MD)
      CASE(mv_type_MD)
         tmc_params%move_types%mv_weight(tmc_params%moves_init_sizes%mv_types(ind)) = probs(ind)/SUM(probs(:))
         tmc_params%move_types%mv_size(tmc_params%moves_init_sizes%mv_types(ind),1) = &
                          tmc_params%moves_init_sizes%mv_sizes(ind) *Pi/180.0_dp !input in degree, calculating in rad
         tmc_params%move_types%acc_prob(tmc_params%moves_init_sizes%mv_types(ind),:) = &
                          tmc_params%moves_init_sizes%start_acc_prob
         tmc_params%print_forces = .TRUE.
      ! proton reordering
      CASE(mv_type_proton_reorder)
         tmc_params%move_types%mv_weight(tmc_params%moves_init_sizes%mv_types(ind)) = &
                          probs(ind)/SUM(probs(:))
         tmc_params%move_types%acc_prob(tmc_params%moves_init_sizes%mv_types(ind),:) = &
                          tmc_params%moves_init_sizes%start_acc_prob
      ! volume moves
      CASE(mv_type_volume_move)
         CALL cp_assert(tmc_params%pressure .GE. 0.0_dp,&
                    cp_warning_level,cp_assertion_failed,routineP,&
                    "no valid pressure defined, but volume move defined. Volume move disabled.",&
                    error)
         IF(tmc_params%pressure .GE. 0.0_dp)THEN
            tmc_params%print_cell = .TRUE. ! print the cell sizes by default
            tmc_params%move_types%mv_size(tmc_params%moves_init_sizes%mv_types(ind),1) = &
                          tmc_params%moves_init_sizes%mv_sizes(ind)/au2a
            tmc_params%move_types%mv_weight(tmc_params%moves_init_sizes%mv_types(ind)) = &
                          probs(ind)/SUM(probs(:))
            tmc_params%move_types%acc_prob(tmc_params%moves_init_sizes%mv_types(ind),:) = &
                          tmc_params%moves_init_sizes%start_acc_prob
            tmc_params%print_cell = .TRUE.
         END IF
      ! parallel tempering swap move
      CASE(mv_type_swap_conf)
         IF(tmc_params%nr_temp.GT.1) THEN
            tmc_params%move_types%mv_weight(tmc_params%moves_init_sizes%mv_types(ind)) = &
                          probs(ind)/SUM(probs(:))
             IF(probs(ind).GT.0.0_dp) &
               tmc_params%move_types%mv_size(tmc_params%moves_init_sizes%mv_types(ind),1) = &
                          probs(ind)
         END IF
         tmc_params%move_types%acc_prob(tmc_params%moves_init_sizes%mv_types(ind),:) = &
                          tmc_params%moves_init_sizes%start_acc_prob
      ! gaussian adaptation
      CASE(mv_type_gausian_adapt)
         tmc_params%move_types%mv_weight(tmc_params%moves_init_sizes%mv_types(ind)) = probs(ind)/SUM(probs(:))
         tmc_params%move_types%mv_size(tmc_params%moves_init_sizes%mv_types(ind),1) = tmc_params%moves_init_sizes%mv_sizes(ind)
         tmc_params%move_types%acc_prob(tmc_params%moves_init_sizes%mv_types(ind),:) = 0.5_dp
      CASE(mv_type_none, mv_type_NMC_moves)
      CASE DEFAULT
        CALL cp_assert(.FALSE.,&
                    cp_failure_level,cp_assertion_failed,routineP,&
                    "entered unknown move type"//cp_to_string(tmc_params%moves_init_sizes%mv_types(ind)),&
                    error)
      END SELECT
    END DO

    ! set NMC start acceptance probability if necessary
    IF(tmc_params%NMC_inp_file.NE."") THEN
      tmc_params%move_types%acc_prob(mv_type_NMC_moves,:) = &
        tmc_params%moves_init_sizes%start_acc_prob
    END IF

    ! distribute to all temperatures
    DO ind=2, tmc_params%nr_temp, 1
       tmc_params%move_types%mv_size(:,ind) = tmc_params%move_types%mv_size(:,1)
    END DO
    DEALLOCATE(probs)
  END SUBROUTINE init_move_types

! *****************************************************************************
!> \brief deallocating the module variables
!> \param move_types pointer to the structure which should be deallocated
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  ! deallocating the module variables
  SUBROUTINE finalize_mv_types(move_types, error)
    TYPE(tmc_move_type), POINTER             :: move_types
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'finalize_mv_types', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(move_types),cp_failure_level,routineP,error,failure)
    CALL move_types_release(move_types, error)
  END SUBROUTINE finalize_mv_types

! *****************************************************************************
!> \brief routine pronts out the probabilities and sized for each type and 
!>        temperature the output is divided into two parts the init, 
!>        which is printed out at the beginning of the programm and 
!>        .NOT.init which are the probabilites and counter printed out every 
!>        print cycle
!> \param move_types pointer to the structure which should be deallocated
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE print_move_types(init, file_io, tmc_params, error)
    LOGICAL                                  :: init
    INTEGER                                  :: file_io
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_move_types', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: c_t
    CHARACTER(LEN=50)                        :: FMT_c, FMT_i, FMT_r
    CHARACTER(LEN=500)                       :: c_a = "", c_b = "", c_c = "", &
                                                c_d = "", c_e = "", c_tit = ""
    INTEGER                                  :: column_size = 10, temper, typ
    LOGICAL                                  :: failure, &
                                                subbox_out = .FALSE., &
                                                type_title = .FALSE.

    failure = .FALSE.
    CPPrecondition(file_io.GT.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params%move_types),cp_failure_level,routineP,error,failure)

    FLUSH(file_io)

    IF(.NOT.init.AND.&
       tmc_params%move_types%mv_weight(mv_type_NMC_moves).GT.0 .AND. &
       ANY(tmc_params%sub_box_size.GT.0.0_dp)) subbox_out = .TRUE.

    ! set the format for each typ to add one column
    WRITE(FMT_c,'("(A,1X,A", I0, ")")') column_size
    WRITE(FMT_i,'("(A,1X,I", I0, ")")') column_size
    WRITE(FMT_r,'("(A,1X,F", I0, ".3)")') column_size
    !IF(init) &
    type_title=.TRUE.

    temp_loop: DO temper=1, tmc_params%nr_temp
       c_tit=""; c_a=""; c_b=""; c_c=""
       IF(init .AND. temper.GT.1) EXIT temp_loop
       WRITE(c_t, "(F10.2)") tmc_params%Temp(temper)
       typ_loop: DO typ=0,SIZE(tmc_params%move_types%mv_weight)
          ! total average
          IF(typ.EQ.0) THEN
             ! line start
             IF(type_title) WRITE(c_tit,TRIM(FMT_c))" type  temperature  |" 
             IF(init) WRITE(c_b,TRIM(FMT_c))        "   I       I        |"
             IF(init) WRITE(c_c,TRIM(FMT_c))        "   V       V        |"
             IF(.NOT.init) WRITE(c_a,TRIM(FMT_c))   "probs  T="//c_t//" |"
             IF(.NOT.init) WRITE(c_b,TRIM(FMT_c))   "counts T="//c_t//" |"
             IF(.NOT.init) WRITE(c_c,TRIM(FMT_c))   "nr_acc T="//c_t//" |"
             IF(subbox_out) THEN
                           WRITE(c_d,TRIM(FMT_c))   "sb_acc T="//c_t//" |"
                           WRITE(c_e,TRIM(FMT_c))   "sb_cou T="//c_t//" |"
             END IF
             ! overall column
             IF(type_title) &
                WRITE(c_tit,TRIM(FMT_c))TRIM(c_tit),   " trajec"
             IF(init) WRITE(c_b,TRIM(FMT_c))TRIM(c_b), "  weight->"
             IF(init) WRITE(c_c,TRIM(FMT_c))TRIM(c_c), "  size  ->"
             IF(.NOT.init) WRITE(c_a,TRIM(FMT_r)) TRIM(c_a), &
                             tmc_params%move_types%acc_prob(typ,temper)
             IF(.NOT.init) WRITE(c_b,TRIM(FMT_i)) TRIM(c_b), &
                             tmc_params%move_types%mv_count(typ,temper)
             IF(.NOT.init) WRITE(c_c,TRIM(FMT_i)) TRIM(c_c), &
                             tmc_params%move_types%acc_count(typ,temper)
             IF(subbox_out) THEN
                           WRITE(c_d,TRIM(FMT_c)) TRIM(c_d),  "."
                           WRITE(c_e,TRIM(FMT_c)) TRIM(c_e),  "."
             END IF
          ELSE
          ! certain move types
             IF(tmc_params%move_types%mv_weight(typ).GT.0 .OR. &
                ! NMC 
                (typ.EQ.mv_type_NMC_moves .AND. tmc_params%NMC_inp_file.NE.""&
                  .AND. tmc_params%nr_NMC_steps.GE.1) .OR. &
                ! swap
                (typ.EQ.mv_type_swap_conf .AND. &
                  ANY(tmc_params%move_types%mv_size(typ,:).GT.0.0_dp)) .OR. &
                ! volume move
                (typ.EQ.mv_type_volume_move .AND. &
                  tmc_params%pressure.GT.0.0_dp) ) THEN

                ! INIT: the weights in the initialisation output
                IF(init) WRITE(c_b,TRIM(FMT_r)) TRIM(c_b), tmc_params%move_types%mv_weight(typ)

                ! acc probabilities
                IF(typ.EQ.mv_type_swap_conf.AND.&
                         temper.EQ.tmc_params%nr_temp)THEN
                  IF(.NOT.init) WRITE(c_a,TRIM(FMT_c)) TRIM(c_a), "---"
                ELSE
                  IF(.NOT.init) WRITE(c_a,TRIM(FMT_r)) TRIM(c_a), tmc_params%move_types%acc_prob(typ,temper)
                END IF
                IF(.NOT.init) WRITE(c_b,TRIM(FMT_i)) TRIM(c_b), tmc_params%move_types%mv_count(typ,temper)
                IF(.NOT.init) WRITE(c_c,TRIM(FMT_i)) TRIM(c_c), tmc_params%move_types%acc_count(typ,temper)
                ! sub box
                IF(subbox_out)THEN
                   !IF(tmc_params%move_types%mv_weight(typ).GT.0) THEN
                   IF(typ.LE.nr_nmc_mv_types) THEN
                      WRITE(c_d,TRIM(FMT_r))TRIM(c_d),tmc_params%move_types%subbox_acc_count(typ,temper)/&
                                                      REAL(tmc_params%move_types%subbox_count(typ,temper),KIND=dp)
                      WRITE(c_e,TRIM(FMT_i))TRIM(c_e),tmc_params%move_types%subbox_count(typ,temper)
                   ELSE
                      WRITE(c_d,TRIM(FMT_c))TRIM(c_d),"-"
                      WRITE(c_e,TRIM(FMT_c))TRIM(c_e),"-"
                   END IF
                END IF
                
                SELECT CASE(typ)
                CASE(mv_type_atom_trans)
                   IF(type_title) WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tit), "atom trans."
                   IF(init) WRITE(c_c,TRIM(FMT_r)) TRIM(c_c), tmc_params%move_types%mv_size(typ,temper)*au2a
                CASE(mv_type_mol_trans)
                   IF(type_title) WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tit), "mol trans"
                   IF(init) WRITE(c_c,TRIM(FMT_r)) TRIM(c_c), tmc_params%move_types%mv_size(typ,temper)*au2a
                CASE(mv_type_mol_rot)
                   IF(type_title) WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tit), "mol rot"
                   IF(init) WRITE(c_c,TRIM(FMT_r)) TRIM(c_c), tmc_params%move_types%mv_size(typ,temper)/(PI/180.0_dp)
                CASE(mv_type_MD)
CALL cp_assert(.TRUE.,&
cp_warning_level,cp_assertion_failed,routineP,&
"md_time_step and nr md_steps not implemented...",&
error)
!                   IF(type_title) WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tit), "HybridMC"
!                   IF(init) WRITE(c_c,TRIM(FMT_c)) TRIM(c_c), "s.above"
!                   IF(init) THEN
!                      WRITE(file_io,*)"   move type: molecular dynamics with file ",NMC_inp_file
!                      WRITE(file_io,*)"                                 with time step [fs] ",md_time_step*au2fs
!                      WRITE(file_io,*)"                                 with number of steps ",md_steps
!                      WRITE(file_io,*)"                                 with velocity changes consists of old vel and ",&
!                         sin(tmc_params%move_types%mv_size(typ,1))*100.0_dp,"% random Gaussian with variance to temperature,"
!                   END IF
                CASE(mv_type_proton_reorder)
                   IF(type_title) WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tit), "H-Reorder"
                   IF(init) WRITE(c_c,TRIM(FMT_c)) TRIM(c_c), "XXX"
                CASE(mv_type_swap_conf)
                   IF(type_title) WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tit), "(swap)"
                   IF(init) WRITE(c_c,TRIM(FMT_r)) TRIM(c_c), tmc_params%move_types%mv_size(mv_type_swap_conf,1)
                CASE(mv_type_NMC_moves)
                   IF(type_title) WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tit), "NMC"
                   IF(init) WRITE(c_c,TRIM(FMT_i)) TRIM(c_c), tmc_params%nr_NMC_steps
                   
                CASE(mv_type_volume_move)
                   IF(type_title) WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tit), "volume"
                   IF(init) WRITE(c_c,TRIM(FMT_r)) TRIM(c_c), tmc_params%move_types%mv_size(typ,temper)*au2a
                CASE(mv_type_gausian_adapt)
                   IF(type_title) WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tit), "gauss adap"
                   IF(init) WRITE(c_c,TRIM(FMT_r)) TRIM(c_c), tmc_params%move_types%mv_size(typ,temper)
                CASE DEFAULT
                   CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                                  "unknown move type "//cp_to_string(typ)//"with weight"//&
                                  cp_to_string(tmc_params%move_types%mv_weight(typ)),&
                                  error)
                END SELECT
             END IF
          END IF
       END DO typ_loop
       IF(init) WRITE (UNIT=file_io,FMT="(/,T2,A)") REPEAT("-",79)
       IF(type_title.AND.temper.LE.1) WRITE(file_io,*) TRIM(c_tit)
       IF(.NOT.init)  WRITE(file_io,*) TRIM(c_a)
       WRITE(file_io,*) TRIM(c_b)
       WRITE(file_io,*) TRIM(c_c)
       IF(subbox_out) WRITE(file_io,*) TRIM(c_d)
       IF(subbox_out) WRITE(file_io,*) TRIM(c_e)
       IF(init) WRITE (UNIT=file_io,FMT="(/,T2,A)") REPEAT("-",79)
    END DO temp_loop
  END SUBROUTINE print_move_types

! *****************************************************************************
!> \brief adaptation of acceptance probability of every kind of change/move 
!>        and the overall acc prob, 
!>        using the acceptance and rejectance information
!> \param move_types structure for storing sizes and probabilities of moves
!> \param pt_el global tree element
!> \param elem sub tree element
!> \param acc input if the element is accepted
!> \param subbox logical if move was with respect to the sub box
!> \param prob_opt if the average probability should be adapted
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE prob_update(move_types, pt_el, elem, acc, subbox, prob_opt, error)
    TYPE(tmc_move_type), POINTER             :: move_types
    TYPE(global_tree_type), OPTIONAL, &
      POINTER                                :: pt_el
    TYPE(tree_type), OPTIONAL, POINTER       :: elem
    LOGICAL, INTENT(IN), OPTIONAL            :: acc, subbox
    LOGICAL, INTENT(IN)                      :: prob_opt
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'prob_update', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: change_res, change_sb_type, &
                                                change_type, conf_moved, &
                                                handle, mv_type
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(move_types),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.(PRESENT(pt_el).AND.PRESENT(subbox)),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    mv_type    = -1
    conf_moved = -1

    change_type = 0
    change_res = 0
    change_sb_type = 0
    ! updating probability of the trajectory 
    IF(PRESENT(pt_el)) THEN
      CPPrecondition(ASSOCIATED(pt_el),cp_failure_level,routineP,error,failure)
      conf_moved = pt_el%mv_conf
      SELECT CASE(pt_el%stat)
      CASE(status_accepted_result)
        change_res = 1
        !-- swaped move is not noted in subtree elements
        IF(pt_el%swaped)THEN
          mv_type = mv_type_swap_conf
          change_type = 1 
        END IF
      CASE(status_rejected_result)
        change_res = -1
        !-- swaped move is not noted in subtree elements
        IF(pt_el%swaped)THEN
          mv_type = mv_type_swap_conf
          change_type = -1 
        END IF
      CASE DEFAULT
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
              routineP,"global elem"//cp_to_string(pt_el%nr)//&
              "has unknown status"//cp_to_string(pt_el%stat),&
              error)
      END SELECT
    END IF 

    IF(PRESENT(elem)) THEN
      CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
      !conf_moved = elem%sub_tree_nr
      conf_moved = elem%temp_created
      mv_type = elem%move_type
      ! for NMC prob update the acceptance is needed
      CPPrecondition(PRESENT(acc),cp_failure_level,routineP,error,failure)
      IF(PRESENT(subbox)) THEN
        ! only update subbox acceptance
        IF(acc) &
          move_types%subbox_acc_count(mv_type, conf_moved) = move_types%subbox_acc_count(mv_type, conf_moved) +1 
        move_types%subbox_count(mv_type, conf_moved) = move_types%subbox_count(mv_type, conf_moved) +1
        ! No more to do
        change_type = 0
        change_res  = 0
        conf_moved  = 0
        ! RETURN
      ELSE
        ! update move type acceptance
        IF(acc) THEN
          change_type = 1
        ELSE
          change_type = -1
        END IF
      END IF
    END IF

    !-- INcrease or DEcrease accaptance rate
    ! MOVE types
    IF(change_type.GT.0) THEN
      move_types%acc_count(mv_type, conf_moved) = move_types%acc_count(mv_type, conf_moved) +1
    END IF

    ! RESULTs
    IF(change_res.GT.0) THEN
      move_types%acc_count(0, conf_moved) = move_types%acc_count(0, conf_moved) +1
    END IF

    IF(conf_moved.GT.0) move_types%mv_count(0, conf_moved) = move_types%mv_count(0, conf_moved) + ABS(change_res)
    IF(mv_type.GE.0 .AND.conf_moved.GT.0) &
      move_types%mv_count(mv_type, conf_moved) = move_types%mv_count(mv_type, conf_moved) + ABS(change_type)

    IF(prob_opt) THEN
      WHERE(move_types%mv_count.GT.0) &
        move_types%acc_prob(:,:) = move_types%acc_count(:,:)/REAL(move_types%mv_count(:,:),KIND=dp)
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE prob_update

! *****************************************************************************
!> \brief add the actual moves to the average probabilities
!> \param move_types structure with move counters and probabilities
!> \param mv_counter move counter for actual performed moves of certain types
!> \param acc_counter counters of acceptance for these moves
!> \param subbox_counter, subbox_acc_counter same for sub box moves
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! ***************************************************************************** 
  SUBROUTINE add_mv_prob(move_types, prob_opt, mv_counter, acc_counter, &
                         subbox_counter, subbox_acc_counter, error)
    TYPE(tmc_move_type), POINTER             :: move_types
    LOGICAL                                  :: prob_opt
    INTEGER, DIMENSION(:, :), OPTIONAL       :: mv_counter, acc_counter, &
                                                subbox_counter, &
                                                subbox_acc_counter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'add_mv_prob', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(move_types),cp_failure_level,routineP,error,failure)
    CPPrecondition(PRESENT(mv_counter).OR.PRESENT(subbox_counter),cp_failure_level,routineP,error,failure)
 
    IF(PRESENT(mv_counter)) THEN
       CPPrecondition(PRESENT(acc_counter),cp_failure_level,routineP,error,failure)
       move_types%mv_count(:,:) = move_types%mv_count(:,:) + mv_counter(:,:)
       move_types%acc_count(:,:)= move_types%acc_count(:,:) + acc_counter(:,:) 
       IF(prob_opt) THEN
         WHERE(move_types%mv_count.GT.0) &
           move_types%acc_prob(:,:) = move_types%acc_count(:,:)/REAL(move_types%mv_count(:,:),KIND=dp)
       END IF
    END IF
 
    IF(PRESENT(subbox_counter)) THEN
       CPPrecondition(PRESENT(subbox_acc_counter),cp_failure_level,routineP,error,failure)
       move_types%subbox_count(:,:)     = move_types%subbox_count(:,:) + subbox_counter(:,:)
       move_types%subbox_acc_count(:,:) = move_types%subbox_acc_count(:,:) + subbox_acc_counter(:,:)
    END IF
  END SUBROUTINE add_mv_prob
  
! *****************************************************************************
!> \brief clear the statistics of accepting/rejection moves
!>        because worker statistics will be add seperatelly on masters counters
!> \param move_types counters for acceptance/rejection
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE clear_move_probs(move_types, error)
    TYPE(tmc_move_type), POINTER             :: move_types
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'clear_move_probs', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    CPPrecondition(ASSOCIATED(move_types),cp_failure_level,routineP,error,failure)

    move_types%acc_prob(:,:) = 0.5_dp
    move_types%acc_count(:,:) = 0
    move_types%mv_count(:,:) = 0
    move_types%subbox_acc_count(:,:) = 0.0_dp
    move_types%subbox_count(:,:) = 0
  END SUBROUTINE clear_move_probs

! *****************************************************************************
!> \brief selects a move type related to the weighings and the entered rnd nr
!> \param move_types structure for storing sizes and probabilities of moves
!> \param rnd random number
!> \param mv_type (result) move type
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  ! function returns a possible move type
  FUNCTION select_random_move_type(move_types, rnd, error) RESULT(mv_type)
    TYPE(tmc_move_type), POINTER             :: move_types
    REAL(KIND=dp)                            :: rnd
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: mv_type

    CHARACTER(LEN=*), PARAMETER :: routineN = 'select_random_move_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rnd_mv, total_moves

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(move_types),cp_failure_level,routineP,error,failure)
    CPPrecondition(rnd.GE.0.0_dp.AND.rnd.LT.1.0_dp,cp_failure_level,routineP,error,failure)

    CALL timeset(routineN,handle)

    total_moves = SUM(move_types%mv_weight(1:nr_nmc_mv_types))
    rnd_mv = total_moves*rnd
    search_loop: DO i=1, SIZE(move_types%mv_weight(1:nr_nmc_mv_types))
       IF(SUM(move_types%mv_weight(1:i)).GE.rnd_mv) THEN
          mv_type = i
          EXIT search_loop
       END IF
    END DO search_loop
    CALL timestop(handle)
  END FUNCTION select_random_move_type

! *****************************************************************************
!> \brief function returns a possible move type
!> \param move_types weights for each move type
!> \param rnd random number for selecting the move type
!> \param mv_type result: the move type
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  FUNCTION select_random_NMC_move_type(move_types, rnd, error) RESULT(mv_type)
    TYPE(tmc_move_type), POINTER             :: move_types
    REAL(KIND=dp)                            :: rnd
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: mv_type

    CHARACTER(LEN=*), PARAMETER :: routineN = 'select_random_NMC_move_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rnd_mv, total_moves

    CPPrecondition(ASSOCIATED(move_types),cp_failure_level,routineP,error,failure)
 
    total_moves = SUM(move_types%mv_weight(1:nr_nmc_mv_types))
    rnd_mv = total_moves*rnd
    search_loop: DO i=1, SIZE(move_types%mv_weight(1:nr_nmc_mv_types))
      IF(SUM(move_types%mv_weight(1:i)).GE.rnd_mv) THEN
        mv_type = i
        EXIT search_loop
      END IF
    END DO search_loop
  END FUNCTION select_random_NMC_move_type
   
END MODULE tmc_move_handle
