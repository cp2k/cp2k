!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief acceptance ratio handling of the different Monte Carlo Moves types 
!>        For each move type and each temperature average acceptence is 
!>        determined. 
!>        For each move is a weight (mv_weight) defined, which defines the 
!>        probability to perform the move.
!>        We distinguish between moves performed on the exact potential 
!>        (move on the master, energy on the energy worker) and
!>        NMC moves, which are performed on the worker using the approximate 
!>        potential. The energies are calculated as usual on the energy worker 
!>        with the exact potential.
!>        The move probabilities to perform a NMC is stored in the NMC move.
!>        The probilities of the single move types (performed with the 
!>        approximate potential) are only compared within the NMC move
!> \par History
!>      11.2012 created [Mandes Schönherr]
!> \author Mandes
! *****************************************************************************

MODULE tmc_move_handle
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE physcon,                         ONLY: au2a => angstrom
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_move_types,                  ONLY: &
       move_types_create, move_types_release, mv_type_MD, mv_type_NMC_moves, &
       mv_type_atom_swap, mv_type_atom_trans, mv_type_gausian_adapt, &
       mv_type_mol_rot, mv_type_mol_trans, mv_type_proton_reorder, &
       mv_type_swap_conf, mv_type_volume_move, tmc_move_type
  USE tmc_stati,                       ONLY: task_type_MC,&
                                             task_type_gaussian_adaptation,&
                                             task_type_ideal_gas
  USE tmc_tree_types,                  ONLY: global_tree_type,&
                                             status_accepted_result,&
                                             status_rejected_result,&
                                             tree_type
  USE tmc_types,                       ONLY: tmc_param_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_move_handle'

  PUBLIC :: finalize_mv_types, print_move_types, read_init_move_types
  PUBLIC :: check_moves
  PUBLIC :: select_random_move_type
  PUBLIC :: prob_update, add_mv_prob
  PUBLIC :: clear_move_probs

CONTAINS

! *****************************************************************************
!> \brief initialization of the different moves, with sizes and probabilities
!> \param move_types pointer to the structure which should be deallocated
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 10.2013
! *****************************************************************************
  SUBROUTINE read_init_move_types(tmc_params, tmc_section, error)
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(section_vals_type), POINTER         :: tmc_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_init_move_types', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: inp_kind_name
    INTEGER                                  :: i, i_rep, i_rep_val, i_tmp, &
                                                ind, n_items, n_NMC_items, &
                                                n_rep_val, nmc_steps
    LOGICAL                                  :: explicit, flag
    REAL(KIND=dp) :: delta_x, init_acc_prob, mv_prob, mv_prob_sum, &
      nmc_init_acc_prob, nmc_prob, nmc_prob_sum, prob_ex
    TYPE(section_vals_type), POINTER         :: move_type_section, nmc_section
    TYPE(tmc_move_type), POINTER             :: move_types

    NULLIFY(move_types, move_type_section, nmc_section)

    n_items         = 0
    n_NMC_items     = 0
    delta_x         = 0.0_dp
    nmc_prob        = 0.0_dp
    mv_prob         = 0.0_dp
    nmc_prob        = 0.0_dp
    mv_prob_sum     = 0.0_dp
    nmc_prob_sum    = 0.0_dp
    prob_ex         = 0.0_dp
    init_acc_prob   = 0.0_dp

    ! the move types on exact potential
    move_type_section => section_vals_get_subs_vals(tmc_section,"MOVE_TYPE",error=error)
    CALL section_vals_get(move_type_section, explicit=explicit, error=error)
    IF(explicit) THEN
      CALL section_vals_get(move_type_section,n_repetition=n_items,error=error)
      mv_prob_sum = 0.0_dp
      DO i_rep = 1, n_items
        CALL section_vals_val_get(move_type_section,"PROB",i_rep_section=i_rep,&
                                  r_val=mv_prob,error=error)
        mv_prob_sum = mv_prob_sum + mv_prob
      END DO
    END IF

    ! get the NMC prameters
    nmc_section => section_vals_get_subs_vals(tmc_section,"NMC_MOVES",error=error)
    CALL section_vals_get(nmc_section, explicit=explicit, error=error)
    IF(explicit) THEN
      ! check the approx potential file, already read
      CALL cp_assert(tmc_params%NMC_inp_file.NE."",&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "Please specify a valid approximate potential.",&
                     error, failure=flag, only_ionode=.TRUE.)

      CALL section_vals_val_get(nmc_section,"NR_NMC_STEPS",&
                                i_val=nmc_steps,error=error)
      CALL cp_assert(nmc_steps.GT.0,&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "Please specify a valid amount of NMC steps (NR_NMC_STEPS {INTEGER}).",&
                     error, failure=flag, only_ionode=.TRUE.)

      CALL section_vals_val_get(nmc_section,"PROB", r_val=nmc_prob,error=error)

      CALL section_vals_val_get(move_type_section,"INIT_ACC_PROB",&
                                r_val=nmc_init_acc_prob,error=error)
      CALL cp_assert(nmc_init_acc_prob.GT.0.0_dp,&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "Please select a valid initial acceptance probability (>0.0) "//&
                     "for the move type "//inp_kind_name,&
                     error)

      move_type_section => section_vals_get_subs_vals(nmc_section,"MOVE_TYPE",error=error)
      CALL section_vals_get(move_type_section,n_repetition=n_NMC_items,error=error)

      ! get the NMC move probability sum
      nmc_prob_sum = 0.0_dp
      DO i_rep = 1, n_NMC_items
        CALL section_vals_val_get(move_type_section,"PROB",i_rep_section=i_rep,&
                                  r_val=mv_prob,error=error)
        nmc_prob_sum = nmc_prob_sum + mv_prob
      END DO
    END IF


    ! get the total weight/amount of move probabilities
    mv_prob_sum = mv_prob_sum + nmc_prob

    IF(n_items+n_NMC_items.GT.0) THEN
      ! initilaize the move array with related sizes, probs, etc.
      CALL move_types_create(tmc_params%move_types, tmc_params%nr_temp, error)

      CALL cp_assert(mv_prob_sum.GT.0.0,&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "The probabilities to perform the moves are "//&
                     "in total less equal 0", error)

      
      ! get the sizes, probs, etc. for each move type and convert units
      DO i_tmp = 1, n_items+n_NMC_items
        ! select the corect section
        IF(i_tmp .GT. n_items) THEN
          i_rep = i_tmp - n_items
          IF(i_rep .EQ. 1) THEN
            ! set the NMC stuff (approx potential)
            tmc_params%move_types%mv_weight(mv_type_NMC_moves) = &
              nmc_prob/REAL(mv_prob_sum, KIND=dp)
            tmc_params%move_types%mv_size(mv_type_NMC_moves,:) = nmc_steps
            tmc_params%move_types%acc_prob(mv_type_NMC_moves,:) = nmc_init_acc_prob
          
            move_type_section => section_vals_get_subs_vals(tmc_section,"NMC_MOVES%MOVE_TYPE",error=error)
            mv_prob_sum = nmc_prob_sum
            ! allocate the NMC move types
            CALL move_types_create(tmc_params%nmc_move_types, tmc_params%nr_temp, error)
            move_types => tmc_params%nmc_move_types
          END IF
        ELSE
          ! the moves on exact potential
          move_type_section => section_vals_get_subs_vals(tmc_section,"MOVE_TYPE",error=error)
          i_rep = i_tmp
          move_types => tmc_params%move_types
        END IF

        CALL section_vals_val_get(move_type_section, "_SECTION_PARAMETERS_",&
                                  c_val=inp_kind_name, i_rep_section=i_rep, &
                                  error=error)
        CALL uppercase(inp_kind_name)
        CALL section_vals_val_get(move_type_section,"SIZE",i_rep_section=i_rep,&
                                  r_val=delta_x,error=error)
        ! move sizes are checked afterwards, because not all moves require a valid move size
        CALL section_vals_val_get(move_type_section,"PROB",i_rep_section=i_rep,&
                                  r_val=mv_prob,error=error)
        CALL cp_assert(mv_prob.GE.0.0_dp,&
                      cp_failure_level,cp_assertion_failed,routineP,&
                      "Please select a valid move probability (>0.0) "//&
                      "for the move type "//inp_kind_name,&
                      error)
        CALL section_vals_val_get(move_type_section,"INIT_ACC_PROB",i_rep_section=i_rep,&
                                  r_val=init_acc_prob,error=error)
        CALL cp_assert(init_acc_prob.GE.0.0_dp,&
                      cp_failure_level,cp_assertion_failed,routineP,&
                      "Please select a valid initial acceptance probability (>0.0) "//&
                      "for the move type "//inp_kind_name,&
                      error)
        ! set the related index and perform unit conversion of move sizes
        SELECT CASE(inp_kind_name)
        ! atom / molecule translation
        CASE("ATOM_TRANS", "MOL_TRANS")
          SELECT CASE(inp_kind_name)
          CASE("ATOM_TRANS")
            ind = mv_type_atom_trans
          CASE("MOL_TRANS")
            ind = mv_type_mol_trans
          CASE DEFAULT
            CALL cp_assert(.FALSE.,&
                      cp_failure_level,cp_assertion_failed,routineP,&
                      "move type is not defined in the translation types",&
                      error)
          END SELECT
          ! convert units
          SELECT CASE(tmc_params%task_type)
          CASE(task_type_MC, task_type_ideal_gas)
            delta_x = delta_x /au2a
          CASE(task_type_gaussian_adaptation)
             !nothing to do (no unit conversion)
          CASE DEFAULT
            CALL cp_assert(.FALSE.,&
                   cp_failure_level,cp_assertion_failed,routineP,&
                   "move type atom / mol trans is not defined for this TMC run type",&
                   error)
          END SELECT
        ! molecule rotation
        CASE("MOL_ROT")
          ind = mv_type_mol_rot
          ! convert units
          SELECT CASE(tmc_params%task_type)
          CASE(task_type_MC, task_type_ideal_gas)
            delta_x = delta_x *PI/180.0_dp
          CASE DEFAULT
            CALL cp_assert(.FALSE.,&
                   cp_failure_level,cp_assertion_failed,routineP,&
                   "move type MOL_ROT is not defined for this TMC run type",&
                    error)
          END SELECT
        ! proton reordering
        CASE("PROT_REORDER")
          ind = mv_type_proton_reorder
          ! the move size is not necessary
          delta_x = 0.0_dp
        ! Hybrid MC (MD)
        CASE("HYBRID_MC")
          ind = mv_type_MD
          delta_x = delta_x *Pi/180.0_dp !input in degree, calculating in rad
          tmc_params%print_forces = .TRUE.
        ! parallel tempering swap move
        CASE("PT_SWAP")
          ind = mv_type_swap_conf
          ! the move size is not necessary
          delta_x = 0.0_dp
          IF(tmc_params%nr_temp.LE.1) THEN
            ! no configurational swapping if only one temperature
            mv_prob = 0.0_dp
            CALL cp_assert(.FALSE.,&
                   cp_warning_level,cp_assertion_failed,routineP,&
                   "Configurational swap disabled, because "//&
                   "Parallel Tempering requires more than one temperature.",&
                   error)
           END IF
        ! volume moves
        CASE("VOL_MOVE")
          ind = mv_type_volume_move
          ! check the selected pressure
          IF(tmc_params%pressure .GE. 0.0_dp)THEN
            delta_x = delta_x /au2a
            tmc_params%print_cell = .TRUE. ! print the cell sizes by default
          ELSE
            CALL cp_assert(.FALSE.,&
                 cp_warning_level,cp_assertion_failed,routineP,&
                 "no valid pressure defined, but volume move defined. "//&
                 "Consequently, the volume move is disabled.",&
                 error)
            mv_prob = 0.0_dp
          END IF
        ! parallel tempering swap move
        CASE("ATOM_SWAP")
          ind = mv_type_atom_swap
          ! the move size is not necessary
          delta_x = 0.0_dp
          ! select the types of atoms swapped
          CALL section_vals_val_get(move_type_section,"ATOMS",i_rep_section=i_rep,&
                 n_rep_val=n_rep_val, error=error)
          IF(n_rep_val.GT.0) THEN
            ALLOCATE(move_types%atom_lists(n_rep_val))
            DO i= 1, n_rep_val
              CALL section_vals_val_get(move_type_section,"ATOMS",&
                     i_rep_section=i_rep, i_rep_val=i, &
                     c_vals=move_types%atom_lists(i)%atoms, &
                     error=error)
              CALL cp_assert(SIZE(move_types%atom_lists(i)%atoms).GT.1,&
                      cp_failure_level,cp_assertion_failed,routineP,&
                      "ATOM_SWAP requires minimum two atom kinds selected. ",&
                      error)
            END DO
          END IF
        ! gaussian adaptation
        CASE("GAUSS_ADAPT")
          ind = mv_type_gausian_adapt
          init_acc_prob = 0.5_dp
        CASE DEFAULT
          CALL cp_assert(.FALSE.,&
                      cp_failure_level,cp_assertion_failed,routineP,&
                      "A unknown move type is selected: "//inp_kind_name,&
                      error)
        END SELECT
        ! check for valid move sizes
        CALL cp_assert(delta_x.GE.0.0_dp,&
                      cp_failure_level,cp_assertion_failed,routineP,&
                      "Please select a valid move size (>0.0) "//&
                      "for the move type "//inp_kind_name,&
                      error)
        ! check if not already set
        IF(move_types%mv_weight(ind).GT.0.0) THEN
          CALL cp_assert(.FALSE.,&
                         cp_failure_level,cp_assertion_failed,routineP,&
                         "TMC: Each move type can be set only once. ",&
                         error)
        END IF

        ! set the move size
        move_types%mv_size(ind,:) = delta_x
        ! set the probability to perform move
        move_types%mv_weight(ind) = mv_prob/mv_prob_sum
        ! set the initial acceptance probability
        move_types%acc_prob(ind,:) = init_acc_prob
      END DO
    ELSE
      CALL cp_assert(.FALSE.,&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "No move type selected, please select at least one.", &
                     error)
    END IF
    mv_prob_sum = SUM(tmc_params%move_types%mv_weight(:))
    flag = .TRUE.
    CPPostcondition(ABS(mv_prob_sum-1.0_dp).LT.0.01,cp_failure_level,routineP,error,failure=flag)
    IF(ASSOCIATED(tmc_params%nmc_move_types)) THEN
      mv_prob_sum = SUM(tmc_params%nmc_move_types%mv_weight(:))
      CPPostcondition(mv_prob_sum.EQ.1.0,cp_failure_level,routineP,error,failure=flag)
    END IF
  END SUBROUTINE read_init_move_types

! *****************************************************************************
!> \brief checks if the moves are possible
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 10.2013
! *****************************************************************************
  SUBROUTINE check_moves(tmc_params, move_types, mol_array, error)
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(tmc_move_type), POINTER             :: move_types
    INTEGER, DIMENSION(:), POINTER           :: mol_array
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_moves', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_j, list_i, ref_k
    LOGICAL                                  :: found

    CPPreconditionNoFail(ASSOCIATED(tmc_params),cp_failure_level,routineP,error)
    CPPreconditionNoFail(ASSOCIATED(move_types),cp_failure_level,routineP,error)

    ! molecule moves need molecule info
    IF(move_types%mv_weight(mv_type_mol_trans).GT.0.0_dp .OR. &
       move_types%mv_weight(mv_type_mol_rot).GT.0.0_dp ) THEN
      ! if there is no molecule information available, 
      !   molecules moves can not be performed
      IF(mol_array(SIZE(mol_array)).EQ.SIZE(mol_array)) &
        CALL cp_assert(.FALSE., &
                       cp_failure_level,cp_assertion_failed,&
                       routineP, "molecule move: there is no molecule "//&
                       "information available. Please specify molecules when "//&
                       "using molecule moves.", &
                       error)
    END IF

    ! for the atom swap move
    IF(move_types%mv_weight(mv_type_atom_swap).GT.0.0_dp) THEN
      ! check if the selected atom swaps are possible
      IF(ASSOCIATED(move_types%atom_lists)) THEN
        DO list_i=1, SIZE(move_types%atom_lists(:))
          DO atom_j=1, SIZE(move_types%atom_lists(list_i)%atoms(:))
            ! check if atoms exists
            found = .FALSE.
            ref_loop: DO ref_k=1, SIZE(tmc_params%atoms(:))
              IF(move_types%atom_lists(list_i)%atoms(atom_j).EQ.&
                 tmc_params%atoms(ref_k)%name) THEN
                found = .TRUE.
                EXIT ref_loop
              END IF
            END DO ref_loop
            CALL cp_assert(found, cp_failure_level,cp_assertion_failed,&
                           routineP, "ATOM_SWAP: The selected atom type ("//&
                           TRIM(move_types%atom_lists(list_i)%atoms(atom_j))//&
                           ") is not contained in the system. ",&
                           error)
            ! check if not be swapped with the same atom type
            IF(ANY(move_types%atom_lists(list_i)%atoms(atom_j).EQ.&
                   move_types%atom_lists(list_i)%atoms(atom_j+1:))) THEN
              CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,&
                     routineP, "ATOM_SWAP can not swap two atoms of same kind ("//&
                     TRIM(move_types%atom_lists(list_i)%atoms(atom_j))//&
                     ")", error)
            END IF
          END DO
        END DO
      ELSE
        ! check if there exisit different atoms
        found = .FALSE.
        IF(SIZE(tmc_params%atoms(:)).GT.1) THEN
          ref_lop: DO ref_k=2, SIZE(tmc_params%atoms(:))
            IF(tmc_params%atoms(1)%name.NE. tmc_params%atoms(ref_k)%name) THEN
              found = .TRUE.
              EXIT ref_lop
            END IF
          END DO ref_lop
        END IF
        CALL cp_assert(found, cp_failure_level,cp_assertion_failed,&
               routineP, "The system contains only a single atom type,"//&
                         " atom_swap is not possible.", error)
      END IF
    END IF
  END SUBROUTINE check_moves

! *****************************************************************************
!> \brief deallocating the module variables
!> \param move_types pointer to the structure which should be deallocated
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  ! deallocating the module variables
  SUBROUTINE finalize_mv_types(tmc_params, error)
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'finalize_mv_types', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CALL move_types_release(tmc_params%move_types, error)
    IF(ASSOCIATED(tmc_params%nmc_move_types)) &
      CALL move_types_release(tmc_params%nmc_move_types, error)
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
    CHARACTER(LEN=500)                       :: c_a, c_b, c_c, c_d, c_e, &
                                                c_tit, c_tmp
    INTEGER                                  :: column_size, move, &
                                                nr_nmc_moves, temper, typ
    LOGICAL                                  :: failure, subbox_out, &
                                                type_title
    TYPE(tmc_move_type), POINTER             :: move_types

    NULLIFY(move_types)
    

    c_a=""; c_b=""; c_c=""
    c_d=""; c_e=""; c_tit=""
    column_size = 10
    subbox_out = .FALSE.
    type_title = .FALSE.
    failure    = .FALSE.
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

    nr_nmc_moves = 0
    IF(ASSOCIATED(tmc_params%nmc_move_types)) THEN
      nr_nmc_moves = SIZE(tmc_params%nmc_move_types%mv_weight(1:))
    END IF

    temp_loop: DO temper=1, tmc_params%nr_temp
       c_tit=""; c_a=""; c_b=""; c_c=""
       IF(init .AND. temper.GT.1) EXIT temp_loop
       WRITE(c_t, "(F10.2)") tmc_params%Temp(temper)
       typ_loop: DO move=0,SIZE(tmc_params%move_types%mv_weight)+nr_nmc_moves
         ! the NMC moves
         IF(move .LE. SIZE(tmc_params%move_types%mv_weight)) THEN
           typ = move
           move_types => tmc_params%move_types
         ELSE
           typ = move - SIZE(tmc_params%move_types%mv_weight)
           move_types => tmc_params%nmc_move_types
         END IF
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
             IF(type_title) THEN
               c_tmp = TRIM(c_tit)
               WRITE(c_tit,TRIM(FMT_c))TRIM(c_tmp),   " trajec"
             END IF
             IF(init)THEN
               c_tmp = TRIM(c_b)
               WRITE(c_b,TRIM(FMT_c))TRIM(c_tmp), "  weight->"
             END IF
             IF(init) THEN
               c_tmp = TRIM(c_c)
               WRITE(c_c,TRIM(FMT_c))TRIM(c_tmp), "  size  ->"
             END IF
             IF(.NOT.init) THEN
               c_tmp = TRIM(c_a)
               WRITE(c_a,TRIM(FMT_r)) TRIM(c_tmp), &
                 move_types%acc_prob(typ,temper)
             END IF
             IF(.NOT.init) THEN
               c_tmp = TRIM(c_b)
               WRITE(c_b,TRIM(FMT_i)) TRIM(c_tmp), &
                 move_types%mv_count(typ,temper)
             END IF
             IF(.NOT.init) THEN
               c_tmp = TRIM(c_c)
               WRITE(c_c,TRIM(FMT_i)) TRIM(c_tmp), &
                 move_types%acc_count(typ,temper)
             END IF
             IF(subbox_out) THEN
               c_tmp = TRIM(c_d)
               WRITE(c_d,TRIM(FMT_c)) TRIM(c_tmp),  "."
               c_tmp = TRIM(c_e)
               WRITE(c_e,TRIM(FMT_c)) TRIM(c_tmp),  "."
             END IF
          ELSE
          ! certain move types
             IF(move_types%mv_weight(typ).GT.0.0_dp)THEN 
                ! INIT: the weights in the initialisation output
                IF(init) THEN
                  c_tmp = TRIM(c_b)
                  WRITE(c_b,TRIM(FMT_r)) TRIM(c_tmp), move_types%mv_weight(typ)
                END IF
                ! acc probabilities
                IF(typ.EQ.mv_type_swap_conf.AND.&
                   temper.EQ.tmc_params%nr_temp)THEN
                  IF(.NOT.init) THEN
                    c_tmp = TRIM(c_a)
                    WRITE(c_a,TRIM(FMT_c)) TRIM(c_tmp), "---"
                  END IF
                ELSE
                  IF(.NOT.init) THEN
                    c_tmp = TRIM(c_a)
                    WRITE(c_a,TRIM(FMT_r)) TRIM(c_tmp), move_types%acc_prob(typ,temper)
                  END IF
                END IF
                IF(.NOT.init) THEN
                  c_tmp = TRIM(c_b)
                  WRITE(c_b,TRIM(FMT_i)) TRIM(c_tmp), move_types%mv_count(typ,temper)
                END IF
                IF(.NOT.init) THEN
                  c_tmp = TRIM(c_c)
                  WRITE(c_c,TRIM(FMT_i)) TRIM(c_tmp), move_types%acc_count(typ,temper)
                END IF
                ! sub box
                IF(subbox_out)THEN
                   IF(move.GT.SIZE(tmc_params%move_types%mv_weight)) THEN
                      c_tmp = TRIM(c_d)
                      WRITE(c_d,TRIM(FMT_r))TRIM(c_tmp),&
                        move_types%subbox_acc_count(typ,temper)/&
                        REAL(move_types%subbox_count(typ,temper),KIND=dp)
                      c_tmp = TRIM(c_e)
                      WRITE(c_e,TRIM(FMT_i))TRIM(c_tmp),&
                        move_types%subbox_count(typ,temper)
                   ELSE
                      c_tmp = TRIM(c_d)
                      WRITE(c_d,TRIM(FMT_c))TRIM(c_tmp),"-"
                      c_tmp = TRIM(c_e)
                      WRITE(c_e,TRIM(FMT_c))TRIM(c_tmp),"-"
                   END IF
                END IF
                
                SELECT CASE(typ)
                CASE(mv_type_atom_trans)
                   IF(type_title) THEN
                     c_tmp = TRIM(c_tit)
                     WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tmp), "atom trans."
                   END IF
                   IF(init) THEN
                     c_tmp = TRIM(c_c)
                     WRITE(c_c,TRIM(FMT_r)) TRIM(c_tmp), &
                       move_types%mv_size(typ,temper)*au2a
                   END IF
                CASE(mv_type_mol_trans)
                   IF(type_title) THEN
                     c_tmp = TRIM(c_tit)
                     WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tmp), "mol trans"
                   END IF
                   IF(init) THEN
                     c_tmp = TRIM(c_c)
                     WRITE(c_c,TRIM(FMT_r)) TRIM(c_tmp), &
                       move_types%mv_size(typ,temper)*au2a
                   END IF
                CASE(mv_type_mol_rot)
                   IF(type_title) THEN
                     c_tmp = TRIM(c_tit)
                     WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tmp), "mol rot"
                   END IF
                   IF(init) THEN
                     c_tmp = TRIM(c_c)
                     WRITE(c_c,TRIM(FMT_r)) TRIM(c_tmp), &
                       move_types%mv_size(typ,temper)/(PI/180.0_dp)
                   END IF
                CASE(mv_type_MD)
CALL cp_assert(.FALSE.,&
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
!                         sin(move_types%mv_size(typ,1))*100.0_dp,"% random Gaussian with variance to temperature,"
!                   END IF
                CASE(mv_type_proton_reorder)
                   IF(type_title) THEN
                     c_tmp = TRIM(c_tit)
                     WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tmp), "H-Reorder"
                   END IF
                   IF(init) THEN
                     c_tmp = TRIM(c_c)
                     WRITE(c_c,TRIM(FMT_c)) TRIM(c_tmp), "XXX"
                   END IF
                CASE(mv_type_swap_conf)
                   IF(type_title) THEN
                     c_tmp = TRIM(c_tit)
                     WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tmp), "PT(swap)"
                   END IF
                   IF(init) THEN
                     c_tmp = TRIM(c_c)
                     WRITE(c_c,TRIM(FMT_c)) TRIM(c_tmp), "XXX" !move_types%mv_size(mv_type_swap_conf,1)
                   END IF
                CASE(mv_type_NMC_moves)
                   IF(type_title) THEN
                     c_tmp = TRIM(c_tit)
                     WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tmp), "NMC:"
                   END IF
                   IF(init) THEN
                     c_tmp = TRIM(c_c)
                     WRITE(c_c,TRIM(FMT_i)) TRIM(c_tmp), &
                       INT(move_types%mv_size(typ,temper))
                   END IF
                CASE(mv_type_volume_move)
                   IF(type_title) THEN
                     c_tmp = TRIM(c_tit)
                     WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tmp), "volume"
                   END IF
                   IF(init) THEN
                     c_tmp = TRIM(c_c)
                     WRITE(c_c,TRIM(FMT_r)) TRIM(c_tmp), &
                       move_types%mv_size(typ,temper)*au2a
                   END IF
                CASE(mv_type_atom_swap)
                   IF(type_title) THEN
                     c_tmp = TRIM(c_tit)
                     WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tmp), "atom swap"
                   END IF
                   IF(init) THEN
                     c_tmp = TRIM(c_c)
                     WRITE(c_c,TRIM(FMT_c)) TRIM(c_tmp), "XXX"
                   END IF
                CASE(mv_type_gausian_adapt)
                   IF(type_title) THEN
                     c_tmp = TRIM(c_tit)
                     WRITE(c_tit,TRIM(FMT_c)) TRIM(c_tmp), "gauss adap"
                   END IF
                   IF(init) THEN
                     c_tmp = TRIM(c_c)
                     WRITE(c_c,TRIM(FMT_r)) TRIM(c_tmp), &
                       move_types%mv_size(typ,temper)
                   END IF
                CASE DEFAULT
                   CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                                  "unknown move type "//cp_to_string(typ)//" with weight"//&
                                  cp_to_string(move_types%mv_weight(typ)),&
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
  ! function returns a possible move type without the PT swap moves 
  !   (are selected in global tree, this routine is for sub tree elements)
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

    total_moves = SUM(move_types%mv_weight(2:))
    rnd_mv = total_moves*rnd
    mv_type = 0
    search_loop: DO i=2, SIZE(move_types%mv_weight(:))
      IF(SUM(move_types%mv_weight(2:i)).GE.rnd_mv) THEN
        mv_type = i
        EXIT search_loop
      END IF
    END DO search_loop

    CALL timestop(handle)
  END FUNCTION select_random_move_type
   
END MODULE tmc_move_handle
