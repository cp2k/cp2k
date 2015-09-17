!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief sets variables for the qmmm pool of pw_types
!> \author Teodoro Laino
! *****************************************************************************
MODULE qmmm_pw_grid
  USE input_constants,                 ONLY: do_par_atom,&
                                             do_qmmm_gauss,&
                                             do_qmmm_swave
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE message_passing,                 ONLY: mp_comm_dup,&
                                             mp_environ
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             PW_MODE_DISTRIBUTED,&
                                             PW_MODE_LOCAL,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_create,&
                                             pw_grid_release
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_dealloc
  USE qmmm_types_low,                  ONLY: qmmm_env_qm_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: qmmm_pw_grid_init
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_pw_grid'
  INTEGER :: qmmm_grid_tag = 0

CONTAINS

! *****************************************************************************
!> \brief Initialize the qmmm pool of pw_type.
!>      Then Main difference w.r.t. QS pw_type pools is that this pool
!>      has [0,L] as boundaries.
!> \param qmmm_env ...
!> \param pw_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_pw_grid_init(qmmm_env, pw_env, error)
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_pw_grid_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: auxbas_grid, Ilevel, pw_mode, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: Maxdr, Mindr
    TYPE(pw_grid_type), POINTER              :: el_struct
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: pool

    NULLIFY(el_struct)
    failure = .FALSE.
    Maxdr   = TINY(0.0_dp)
    Mindr   = HUGE(0.0_dp)
    IF ((qmmm_env%qmmm_coupl_type == do_qmmm_gauss).OR.(qmmm_env%qmmm_coupl_type == do_qmmm_swave)) THEN
       CALL pw_env_get(pw_env=pw_env,&
                       pw_pools=pw_pools,&
                       auxbas_grid=auxbas_grid,&
                       error=error)
       !
       IF (ASSOCIATED(qmmm_env%aug_pools)) THEN
          CALL pw_pools_dealloc(qmmm_env%aug_pools,error=error)
       END IF
       ALLOCATE(qmmm_env%aug_pools(SIZE(pw_pools)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !
       DO Ilevel = 1, SIZE(pw_pools)
          NULLIFY(pool, qmmm_env%aug_pools(Ilevel)%pool)
          pool      => pw_pools(Ilevel)%pool
          NULLIFY(el_struct)
          pw_mode = PW_MODE_DISTRIBUTED
          ! Parallelization scheme
          IF (qmmm_env%par_scheme==do_par_atom) THEN
             pw_mode = PW_MODE_LOCAL
          END IF

          CALL pw_grid_create_copy_no_pbc(pool%pw_grid,el_struct, &
               pw_mode=pw_mode, error=error)
          CALL pw_pool_create(qmmm_env%aug_pools(Ilevel)%pool,&
               pw_grid=el_struct,error=error)

          Maxdr = MAX(Maxdr,el_struct%dr)
          Mindr = MIN(Mindr,el_struct%dr)
          IF (ALL(Maxdr.EQ.el_struct%dr)) qmmm_env%gridlevel_info%coarser_grid = Ilevel
          IF (ALL(Mindr.EQ.el_struct%dr)) qmmm_env%gridlevel_info%auxbas_grid  = Ilevel

          CALL pw_grid_release(el_struct,error=error)

       END DO
    END IF

  END SUBROUTINE qmmm_pw_grid_init

! *****************************************************************************
!> \brief creates a copy of pw_grid_in in which the pbc have been removed
!>      (by adding a point for the upper boundary)
!> \param pw_grid_in the pw grid to duplicate
!> \param pw_grid_out the output pw_grid_type
!> \param pw_mode ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2004 created [tlaino]
!>      04.2005 completly rewritten the duplicate routine, fixed parallel
!>              behaviour, narrowed scope to copy to non pbc and renamed
!>              accordingly [fawzi]
!>      06.2007 moved to new module [jgh]
!> \author Fawzi, Teo
! *****************************************************************************
  SUBROUTINE pw_grid_create_copy_no_pbc(pw_grid_in, pw_grid_out, pw_mode, error)
    TYPE(pw_grid_type), POINTER              :: pw_grid_in, pw_grid_out
    INTEGER, INTENT(IN), OPTIONAL            :: pw_mode
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_grid_create_copy_no_pbc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: pw_mode_loc, stat
    INTEGER, DIMENSION(:), POINTER           :: pos_of_x
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(pw_grid_in%ngpts_cut>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(pw_grid_out),cp_failure_level,routineP,error,failure)
    pw_mode_loc = pw_grid_in%para%mode
    IF (PRESENT(pw_mode)) pw_mode_loc = pw_mode
    CALL pw_grid_create(pw_grid_out, pw_grid_in%para%group,error=error)
    qmmm_grid_tag = qmmm_grid_tag + 1
    pw_grid_out %id_nr = qmmm_grid_tag
    pw_grid_out % ref_count  = 1
    pw_grid_out % reference  = 0
    pw_grid_out%bounds=pw_grid_in%bounds
    pw_grid_out%bounds(2,:)=pw_grid_out%bounds(2,:) + 1
    IF (pw_mode_loc==PW_MODE_DISTRIBUTED) THEN
       pw_grid_out%bounds_local=pw_grid_in%bounds_local
       IF ( pw_grid_in%bounds_local(2,1)==pw_grid_in%bounds(2,1).AND.&
            pw_grid_in%bounds_local(1,1)<=pw_grid_in%bounds(2,1)) THEN
          pw_grid_out%bounds_local(2,1)=pw_grid_out%bounds_local(2,1)+1
       END IF
       pw_grid_out%bounds_local(2,2)=pw_grid_out%bounds_local(2,2)+1
       pw_grid_out%bounds_local(2,3)=pw_grid_out%bounds_local(2,3)+1
    ELSE
       pw_grid_out%bounds_local = pw_grid_out%bounds
    END IF
    pw_grid_out%npts            =  pw_grid_in%npts  + 1
    pw_grid_out%ngpts           =  PRODUCT ( INT(pw_grid_out% npts, KIND=int_8) )
    pw_grid_out%ngpts_cut=0
    pw_grid_out%npts_local=pw_grid_out%bounds_local(2,:)-pw_grid_out%bounds_local(1,:)+1
    pw_grid_out%ngpts_local=PRODUCT(pw_grid_out%npts_local)
    pw_grid_out%ngpts_cut_local=0
    pw_grid_out%dr              =  pw_grid_in%dr
    pw_grid_out%dh              =  pw_grid_in%dh
    pw_grid_out%dh_inv          =  pw_grid_in%dh_inv
    pw_grid_out%orthorhombic    =  pw_grid_in%orthorhombic
    pw_grid_out%dvol            =  pw_grid_in%dvol
    pw_grid_out%vol             =  pw_grid_in%vol*REAL(pw_grid_out%ngpts,dp)&
         /REAL(pw_grid_in%ngpts,dp) !FM do not modify?
    pw_grid_out%cutoff          =  pw_grid_in%cutoff
    NULLIFY(pw_grid_out%mapl%pos, pw_grid_out%mapl%neg,&
         pw_grid_out%mapm%pos,pw_grid_out%mapm%neg,&
         pw_grid_out%mapn%pos,pw_grid_out%mapn%neg)

    !para
    CALL mp_environ ( pw_grid_out % para % group_size, &
                      pw_grid_out % para % my_pos, &
                      pw_grid_out % para % group )
    pw_grid_out % para % group_head_id = pw_grid_in % para % group_head_id
    pw_grid_out % para % group_head = &
         ( pw_grid_out % para % group_head_id == pw_grid_out % para % my_pos )
    pw_grid_out % para % mode = pw_mode_loc
    NULLIFY(pw_grid_out % para % yzp,&
            pw_grid_out % para % yzq,&
            pw_grid_out % para % nyzray,&
            pw_grid_out % para % bo)
    ALLOCATE(pos_of_x(pw_grid_out%bounds(1,1):pw_grid_out%bounds(2,1)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    pos_of_x(:pw_grid_out%bounds(2,1)-1)=pw_grid_in%para%pos_of_x
    pos_of_x(pw_grid_out%bounds(2,1))=pos_of_x(pw_grid_out%bounds(2,1)-1)
    pw_grid_out%para%pos_of_x => pos_of_x
    pw_grid_out % para % rs_dims = pw_grid_in % para % rs_dims
    IF (PRODUCT(pw_grid_in % para % rs_dims)/=0) THEN
       CALL mp_comm_dup  ( pw_grid_in  % para % rs_group, &
                           pw_grid_out % para % rs_group )
    END IF
    pw_grid_out % para % rs_pos = pw_grid_in % para % rs_pos
    pw_grid_out % para % rs_mpo = pw_grid_in % para % rs_mpo

    NULLIFY(pw_grid_out%g,pw_grid_out%gsq,pw_grid_out%g_hat)
    CPPrecondition(pw_grid_in%grid_span==FULLSPACE,cp_failure_level,routineP,error,failure)
    pw_grid_out%grid_span=pw_grid_in%grid_span
    pw_grid_out%have_g0=.FALSE.
    pw_grid_out%first_gne0=HUGE(0)
    NULLIFY(pw_grid_out%gidx)
    pw_grid_out%spherical                 = .FALSE.
    pw_grid_out % para % ray_distribution = .FALSE.
    pw_grid_out % para % blocked          = .FALSE.
  END SUBROUTINE pw_grid_create_copy_no_pbc
END MODULE qmmm_pw_grid
