!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_proj_methods [1.0] *
!!
!!   NAME
!!     pao_proj_methods
!!
!!   FUNCTION
!!     projection and injection between min bas and full bas
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     2.2002 created [fawzi]
!!     4.2002 uses projection object and qs_env [fawzi]
!!     6.2002 added atomic orthogonalization
!!
!!   SOURCE
!*****************************************************************************
module pao_proj_methods
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level
  use cp_error_handling, only: cp_error_type, cp_assert, cp_error_message,&
       cp_assertion_failed, cp_debug, cp_precondition_failed, cp_error_init,&
       cp_dealloc_ref, cp_error_propagate_error
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use pao_types, only: pao_env_type, pao_glob_angles_type, pao_env_get
  use cp_matrix_utils, only: cp_matrix_block_iterator, cp_next, cp_init,&
       cp_dealloc_ref, cp_get
  use cp_sparse_matrix, only: cp_sparse_matrix_p_type, cp_sparse_matrix_type,&
       cp_sp_create
  use pao_l_angles_methods, only: pao_local_angles_type, cp_get, cp_init, &
       cp_dealloc_ref, cp_set
  use pao_qs_env_methods, only: pao_qs_env_get
  use sparse_matrix_types, only: real_matrix_type, real_matrix_p_type,&
       get_block_node, add_block_node
  use qs_environment_types, only: qs_environment_type, get_qs_env
  use message_passing, only: mp_sum
  use cp_b_matrix_structure, only: cp_b_matrix_struct_type
  implicit none
  private

  logical, private, parameter :: debug_this_module=.false.
  character(len=*), private, parameter :: moduleN='pao_proj_methods'

  public :: cp_get

! projection - injection
  public :: cp_m2f, cp_f2m, cp_m2f_with_B_M_pi

  public :: pao_g_ang_get_min_m, sm_matrix_transf_to_m,&
       sm_matrix_transf_to_f, pao_m_m2f_with_b_m_pi,&
       sm_matrix_transf_p_to_m

  interface cp_get
     module procedure pao_g_ang_get_min_m
  end interface

  ! transfer from full basis to min basis
  interface cp_f2m
     module procedure sm_matrix_transf_to_m
  end interface

  !transfer from min basis to full basis
  interface cp_m2f
     module procedure sm_matrix_transf_to_f
  end interface
  
  ! tranfser to full basis with B matrix pi
  interface cp_m2f_with_B_M_pi
     module procedure pao_m_m2f_with_b_m_pi
  end interface
!!***
!****************************************************************************

contains

!!****f* pao_types/pao_m_m2f_with_b_m_pi [1.0] *
!!
!!   NAME
!!     pao_m_m2f_with_b_m_pi
!!
!!   SYNOPSIS
!!     Subroutine pao_m_m2f_with_b_m_pi(source_m, target_m, glob_angles,&
!!         qs_env, global_env, error)
!!       Type(real_matrix_type), Pointer:: source_m
!!       Type(real_matrix_type), Pointer:: target_m
!!       Type(pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_m_m2f_with_b_m_pi
!!
!!   FUNCTION
!!     Transform the given matrix to the full basis with B= N U i and pi=i^*
!!     Puts N U i source_m pi  into the target_m
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     source_m: matrix in the min basis
!!     target_m: the result matrix in the full basis
!!     error: variable to control error logging, stopping,...
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine pao_m_m2f_with_b_m_pi(source_m, target_m, glob_angles,&
       qs_env, global_env, error)
    type(real_matrix_type), pointer :: source_m ! make target?
    type(real_matrix_type), pointer :: target_m ! make target?
    type(pao_glob_angles_type), intent(in) :: glob_angles
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_types:pao_m_m2f_with_b_m_pi'
    type(cp_matrix_block_iterator) :: iterator
    integer :: brow,bcol, stat, max_bsize,i,j
    integer, dimension(:), pointer :: l_index
    real(kind=wp), dimension(:,:), pointer :: block_val, u_row
    real(kind=wp),allocatable, dimension(:,:) :: t1_m
    type(pao_local_angles_type) :: l_angles_col, l_angles_row
    type(cp_error_type) :: iError
    type(pao_env_type), pointer :: pao_env
    call cp_error_init(iError,template_error=error)
    nullify(block_val,l_index,u_row,pao_env)
    failure=.false.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       call cp_init(l_angles_col,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       call cp_init(l_angles_row,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       call pao_qs_env_get(qs_env,global_env=global_env, &
            pao_max_full_bas=max_bsize)
       allocate(t1_m(max_bsize,max_bsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       call cp_init(iterator,source_m,error=iError)
       call cp_error_propagate_error(iError,routineP,&
            CPSourceFileRef,&
            error,failure)
       if (.not.failure) then
          blocksDo: do
             if ((.not.cp_next(iterator,error=error)).or.failure) exit
             call cp_get(iterator,block_row=brow,block_col=bcol,&
                  block_val=block_val)
! get u
             call cp_set(l_angles_col, atom_nr=bcol,error=error)
             call cp_get(l_angles_col, proj_indexes=l_index)
             call cp_set(l_angles_row, atom_nr=brow,error=error)
             call cp_get(l_angles_row,NUi_injection=u_row)

             ! t1_m = block_val pi
             t1_m=0.0_wp
             CPPrecondition(associated(l_index),cp_failure_level,routineP,error,failure)
             if (.not.failure) then
                call cp_assert(size(l_index)==size(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                call cp_assert(size(block_val,1)==size(u_row,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             end if
             if (.not.failure) then ! t1_m = block_val pi
                do j=1,size(block_val,2)
                   do i=1,size(block_val,1)
                      t1_m(i,l_index(j))=block_val(i,j)
                   end do
                end do
             end if

! get full basis block
             call get_block_node(target_m, block_col=bcol,&
                  block_row=brow, block=block_val)
             if (.not.associated(block_val)) then
                call add_block_node(target_m, block_col=bcol,&
                     block_row=brow, block=block_val)
                CPPostcondition(associated(block_val),cp_failure_level,routineP,error,failure)
             end if
             if (.not.failure) then
                call cp_assert(size(block_val,1)==size(u_row,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             end if
! block_val = NUi t1_m
             if (.not.failure) then
                call dgemm('N','N',size(u_row,1),size(block_val,2),&
                     size(u_row,2),1.0_wp,&
                     u_row,size(u_row,1),&
                     t1_m,size(t1_m,1),&
                     0.0_wp,block_val,size(block_val,1))
             end if

          end do blocksDo
       end if
       call cp_dealloc_ref(iterator,error=error)
       call cp_dealloc_ref(l_angles_col,error=error)
       call cp_dealloc_ref(l_angles_row,error=error)
       deallocate(t1_m,stat=stat) ! just to be sure...
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine pao_m_m2f_with_b_m_pi
!***************************************************************************


!!****f* pao_types/sm_matrix_transf_to_f [1.0] *
!!
!!   NAME
!!     sm_matrix_transf_to_f
!!
!!   SYNOPSIS
!!     Subroutine sm_matrix_transf_to_f(source_m, target_m, glob_angles,&
!!         qs_env, global_env, error)
!!       Type(real_matrix_type), Pointer:: source_m
!!       Type(real_matrix_type), Pointer:: target_m
!!       Type(pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_matrix_transf_to_f
!!
!!   FUNCTION
!!     Transform the given matrix to the full basis.
!!     Puts N^(-1) U i source_m pi U^* N^(-1*)  into the target_m
!!
!!   NOTES
!!     could be optimized by caching N^(-1) (i.e. atomic S^(1/2))
!!
!!   INPUTS
!!     source_m: matrix in the min basis
!!     target_m: the result matrix in the full basis
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine sm_matrix_transf_to_f(source_m, target_m, glob_angles,&
       qs_env, global_env, error)
    type(real_matrix_type), pointer :: source_m ! make target?
    type(real_matrix_type), pointer :: target_m ! make target ?
    type(pao_glob_angles_type), intent(in) :: glob_angles
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='sm_matrix_transf_to_f',&
         routineP=moduleN//':'//routineN
    type(cp_matrix_block_iterator) :: iterator
    integer :: brow,bcol, stat, max_bsize,i,j
    real(kind=wp), dimension(:,:), pointer :: block_val, u_col, u_row
    real(kind=wp),allocatable, dimension(:,:) :: t1_m
    type(pao_local_angles_type) :: l_angles_col, l_angles_row
    type(cp_error_type) :: iError
    type(pao_env_type), pointer :: pao_env
    call cp_error_init(iError,template_error=error)
    nullify(block_val,u_col, u_row,pao_env)
    failure=.false.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       call cp_init(l_angles_col,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       call cp_init(l_angles_row,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       call pao_qs_env_get(qs_env,global_env=global_env,&
            pao_max_full_bas=max_bsize)
       allocate(t1_m(max_bsize,max_bsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       call cp_init(iterator,source_m,error=iError)
       call cp_error_propagate_error(iError,routineP,&
            CPSourceFileRef,&
            error,failure)
       if (.not.failure) then
          blocksDo: do
             if ((.not.cp_next(iterator,error=error)).or.failure) exit
             call cp_get(iterator,block_row=brow,block_col=bcol,&
                  block_val=block_val)
! get u
             call cp_set(l_angles_col, atom_nr=bcol,error=error)
             call cp_get(l_angles_col,NUi_injection=u_col)
             call cp_set(l_angles_row, atom_nr=brow,error=error)
             call cp_get(l_angles_row,NUi_injection=u_row)

             if (.not.failure) then
                call cp_assert(size(u_col,2)==size(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                call cp_assert(size(u_row,2)==size(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             end if

! t1_m = block_val (N U i)^*
             if (.not.failure) then
                call dgemm('N','T',size(block_val,1),size(u_col,1),&
                     size(u_col,2),1.0_wp,&
                     block_val,size(block_val,1),&
                     u_col,size(u_col,1),&
                     0.0_wp,t1_m,size(t1_m,1))

                ! get full basis block
                call get_block_node(target_m, block_col=bcol,&
                     block_row=brow, block=block_val)
                if (.not.associated(block_val)) then
                   call add_block_node(target_m, block_col=bcol,&
                        block_row=brow, block=block_val)
                   CPPostcondition(associated(block_val),cp_failure_level,routineP,error,failure)
                end if
             end if

             if (.not.failure) then
                call cp_assert(size(u_col,1)==size(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                call cp_assert(size(u_row,1)==size(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             end if

! block_val = NUi t1_m
             if (.not.failure) then
                call dgemm('N','N',size(u_row,1),size(u_col,1),&
                     size(u_row,2),1.0_wp,&
                     u_row,size(u_row,1),&
                     t1_m,size(t1_m,1),&
                     0.0_wp,block_val,size(block_val,1))
             end if

          end do blocksDo
       end if
       call cp_dealloc_ref(iterator,error=error)
       call cp_dealloc_ref(l_angles_col,error=error)
       call cp_dealloc_ref(l_angles_row,error=error)
       deallocate(t1_m,stat=stat) ! just to be sure...
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine sm_matrix_transf_to_f
!***************************************************************************

!!****f* pao_types/sm_matrix_transf_to_m [1.0] *
!!
!!   NAME
!!     sm_matrix_transf_to_m
!!
!!   SYNOPSIS
!!     Subroutine sm_matrix_transf_to_m(source_m, target_m, glob_angles,&
!!         qs_env, global_env, error)
!!       Type(real_matrix_type), Pointer:: source_m
!!       Type(real_matrix_type), Pointer:: target_m
!!       Type(pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_matrix_transf_to_m
!!
!!   FUNCTION
!!     Transform the given matrix to the minimal basis.
!!     Puts pi U^* N source_m N U i into the target_m
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     source_m: matrix in the full basis
!!     target_m: the result matrix in the small basis
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine sm_matrix_transf_to_m(source_m, target_m, glob_angles,&
       qs_env, global_env, error)
    type(real_matrix_type), pointer :: source_m ! make target?
    type(real_matrix_type), pointer :: target_m ! make target?
    type(pao_glob_angles_type), intent(in) :: glob_angles
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='sm_matrix_transf_to_m',&
         routineP=moduleN//':'//routineN
    type(cp_matrix_block_iterator) :: iterator
    integer :: brow,bcol, stat, max_bsize,i,j
    real(kind=wp), dimension(:,:), pointer :: block_val, u_col, u_row
    real(kind=wp),allocatable, dimension(:,:) :: t1_m
    type(pao_local_angles_type) :: l_angles_col, l_angles_row
    type(cp_error_type) :: iError
    type(pao_env_type), pointer :: pao_env
    call cp_error_init(iError,template_error=error)
    nullify(block_val,block_val, u_col, u_row,pao_env)
    failure=.false.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       call cp_init(l_angles_col,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       call cp_init(l_angles_row,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       call pao_qs_env_get(qs_env,global_env=global_env,&
            pao_max_full_bas=max_bsize)
       allocate(t1_m(max_bsize,max_bsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       call cp_init(iterator,source_m,error=iError)
       call cp_error_propagate_error(iError,routineP,&
            CPSourceFileRef,&
            error,failure)
       if (.not.failure) then
          blocksDo: do
             if ((.not.cp_next(iterator,error=error)).or.failure) exit
             call cp_get(iterator,block_row=brow,block_col=bcol,&
                  block_val=block_val)
! get u
             call cp_set(l_angles_col, atom_nr=bcol,error=error)
             call cp_get(l_angles_col,NUi_injection=u_col,error=error)
             call cp_set(l_angles_row, atom_nr=brow,error=error)
             call cp_get(l_angles_row,NUi_injection=u_row,error=error)

             if (.not.failure) then
                call cp_assert(size(u_row,1)==size(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                call cp_assert(size(u_col,1)==size(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             end if

! t1_m = source NUi
             call dgemm('N','N',size(block_val,1),size(u_col,2),&
                  size(u_col,1),1.0_wp,&
                  block_val,size(block_val,1),&
                  u_col,size(u_col,1),&
                  0.0_wp,t1_m,size(t1_m,1))
! get small basis block
             call get_block_node(target_m, block_col=bcol,&
                  block_row=brow, block=block_val)
             if (.not.associated(block_val)) then
                call add_block_node(target_m, block_col=bcol,&
                     block_row=brow, block=block_val)
                CPPostcondition(associated(block_val),cp_failure_level,routineP,error,failure)
             end if

             if (.not.failure) then
                call cp_assert(size(u_row,2)==size(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                call cp_assert(size(u_col,2)==size(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             end if
             ! target_m = (NUi)^* t1_m
             if (.not.failure) then 
                call dgemm('T','N',size(u_row,2),size(u_col,2),&
                     size(u_row,1),1.0_wp,&
                     u_row,size(u_row,1),&
                     t1_m,size(t1_m,1),&
                     0.0_wp,block_val,size(block_val,1))
             end if

          end do blocksDo
       end if
       call cp_dealloc_ref(iterator,error=error)
       call cp_dealloc_ref(l_angles_col,error=error)
       call cp_dealloc_ref(l_angles_row,error=error)
       deallocate(t1_m,stat=stat) ! just to be sure...
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine sm_matrix_transf_to_m
!***************************************************************************

!!****f* pao_types/sm_matrix_transf_p_to_m [1.0] *
!!
!!   NAME
!!     sm_matrix_transf_p_to_m
!!
!!   SYNOPSIS
!!     Subroutine sm_matrix_transf_p_to_m(source_m, target_m, glob_angles,&
!!         qs_env, global_env, error)
!!       Type(real_matrix_type), Pointer:: source_m
!!       Type(real_matrix_type), Pointer:: target_m
!!       Type(pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_matrix_transf_p_to_m
!!
!!   FUNCTION
!!     Transform the given matrix to the minimal basis.
!!     Puts pi U^* N^(-1) source_m N^(-1) U i into the target_m
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     source_m: matrix in the full basis
!!     target_m: the result matrix in the small basis
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine sm_matrix_transf_p_to_m(source_m, target_m, glob_angles,&
       qs_env, global_env, error)
    type(real_matrix_type), pointer :: source_m ! make target?
    type(real_matrix_type), pointer :: target_m ! make target?
    type(pao_glob_angles_type), intent(in) :: glob_angles
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='sm_matrix_transf_p_to_m',&
         routineP=moduleN//':'//routineN
    type(cp_matrix_block_iterator) :: iterator
    integer :: brow,bcol, stat, max_bsize,i,j
    real(kind=wp), dimension(:,:), pointer :: block_val, u_col, u_row,&
         s_col,s_row
    real(kind=wp),allocatable, dimension(:,:) :: t1_m,t2_m
    type(pao_local_angles_type) :: l_angles_col, l_angles_row
    type(real_matrix_p_type), dimension(:), pointer :: s
    type(cp_error_type) :: iError
    type(pao_env_type), pointer :: pao_env
    call cp_error_init(iError,template_error=error)
    nullify(block_val,block_val, u_col, u_row,pao_env)
    failure=.false.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       call cp_init(l_angles_col,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       call cp_init(l_angles_row,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       call get_qs_env(qs_env,s=s)
       call pao_qs_env_get(qs_env,global_env=global_env,&
            pao_max_full_bas=max_bsize)
       allocate(t1_m(max_bsize,max_bsize),t2_m(max_bsize,max_bsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       call cp_init(iterator,source_m,error=iError)
       call cp_error_propagate_error(iError,routineP,&
            CPSourceFileRef,&
            error,failure)
       if (.not.failure) then
          blocksDo: do
             if ((.not.cp_next(iterator,error=error)).or.failure) exit
             call cp_get(iterator,block_row=brow,block_col=bcol,&
                  block_val=block_val)
! get u
             call cp_set(l_angles_col, atom_nr=bcol,error=error)
             call cp_get(l_angles_col,NUi_injection=u_col,error=error)
             call cp_set(l_angles_row, atom_nr=brow,error=error)
             call cp_get(l_angles_row,NUi_injection=u_row,error=error)

             if (.not.failure) then
                call cp_assert(size(u_row,1)==size(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                call cp_assert(size(u_col,1)==size(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             end if

             call get_block_node(s(1)%matrix,block_row=bcol,block_col=bcol,&
                  block=s_col)
             call get_block_node(s(1)%matrix,block_row=brow,block_col=brow,&
                  block=s_row)
             if (.not.failure) then
                call cp_assert(size(s_col,1)==size(block_val,2).and.&
                     size(s_col,1)==size(s_col,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                call cp_assert(size(s_row,2)==size(block_val,1).and.&
                     size(s_row,1)==size(s_row,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             end if
             ! t1_m = block_val S
             if (.not.failure) then
                call dgemm('N','N',size(block_val,1),size(block_val,2),&
                     size(block_val,2),1.0_wp,&
                     block_val,size(block_val,1),&
                     s_col,size(s_col,1),&
                     0.0_wp,t1_m,size(t1_m,1))
             end if
             ! t2_m = S t1_m (could be partially faster if done later)
             if (.not.failure) then
                call dgemm('N','N',size(block_val,1),size(block_val,2),&
                     size(block_val,1),1.0_wp,&
                     s_row,size(s_row,1),&
                     t1_m,size(t1_m,1),&
                     0.0_wp,t2_m,size(t2_m,1))
             end if

! t1_m = t2_m NUi
             call dgemm('N','N',size(block_val,1),size(u_col,2),&
                  size(u_col,1),1.0_wp,&
                  t2_m,size(t2_m,1),&
                  u_col,size(u_col,1),&
                  0.0_wp,t1_m,size(t1_m,1))
! get small basis block
             call get_block_node(target_m, block_col=bcol,&
                  block_row=brow, block=block_val)
             if (.not.associated(block_val)) then
                call add_block_node(target_m, block_col=bcol,&
                     block_row=brow, block=block_val)
                CPPostcondition(associated(block_val),cp_failure_level,routineP,error,failure)
             end if

             if (.not.failure) then
                call cp_assert(size(u_row,2)==size(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                call cp_assert(size(u_col,2)==size(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             end if
             ! target_m = (NUi)^* t1_m
             if (.not.failure) then 
                call dgemm('T','N',size(u_row,2),size(u_col,2),&
                     size(u_row,1),1.0_wp,&
                     u_row,size(u_row,1),&
                     t1_m,size(t1_m,1),&
                     0.0_wp,block_val,size(block_val,1))
             end if

          end do blocksDo
       end if
       call cp_dealloc_ref(iterator,error=error)
       call cp_dealloc_ref(l_angles_col,error=error)
       call cp_dealloc_ref(l_angles_row,error=error)
       deallocate(t1_m,stat=stat) ! just to be sure...
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine sm_matrix_transf_p_to_m
!***************************************************************************

!!****f* pao_proj_methods/pao_g_ang_get_min_m [1.0] *
!!
!!   NAME
!!     pao_g_ang_get_min_m
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_get_min_m(glob_angles, qs_env, global_env,&
!!         min_overlap_m, min_hamiltonian_m, error)
!!       Type(pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_sparse_matrix_p_type), Dimension(:), Pointer, Optional::&
!!         min_hamiltonian_m
!!       Type(cp_sparse_matrix_p_type), Dimension(:), Pointer::&
!!         min_overlap_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_get_min_m
!!
!!   FUNCTION
!!     returns the minimal matrices, calculating them if necessary
!!
!!   NOTES
!!     very ugly: to avoid circular dependencies this is not in the get
!!     of the pao glob angles. to do: fix it!!!
!!
!!   INPUTS
!!     glob_angles: 
!!     min_hamiltonian_m: the hamiltonian in the small basis. 
!!            If not yet calculated calculates it (allocating it if necessary)
!!     min_overlap_m: the overlap matrix in the small basis
!!            If not yet calculated calculates it (allocating it if necessary)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine pao_g_ang_get_min_m(glob_angles,qs_env,global_env,&
       min_overlap_m,  min_hamiltonian_m, error)
    type(pao_glob_angles_type), intent(inout) :: glob_angles
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    type(cp_sparse_matrix_p_type), dimension(:), pointer, optional :: &
         min_hamiltonian_m
    type(cp_sparse_matrix_p_type), dimension(:), pointer ::  min_overlap_m
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: handle, stat
    character(len=*), parameter :: routineN='pao_g_ang_get_min_m',&
         routineP=moduleN//':'//routineN
    type(cp_error_type) :: iError
    type(real_matrix_p_type), dimension(:), pointer :: full_h, full_s
    type(cp_b_matrix_struct_type), pointer :: min_bas_m_struct, &
         full_bas_m_struct
    type(pao_env_type), pointer :: pao_env
    call cp_error_init(iError, template_error=error)
    failure=.false.
    nullify(full_h, full_s, pao_env)

    !call timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       call get_qs_env(qs_env,k=full_h,s=full_s,pao_env=pao_env)
       CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       call pao_env_get(pao_env,min_bas_m_struct=min_bas_m_struct,&
            full_bas_diag_m_struct=full_bas_m_struct,error=error)
       !       if (present(min_overlap_m)) then
       if (.not.glob_angles%m_s_m_valid) then
          if (.not.associated(glob_angles%min_overlap_m)) then
             allocate(glob_angles%min_overlap_m(1),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             call cp_sp_create(glob_angles%min_overlap_m(1)%matrix,&
                  matrix_structure=min_bas_m_struct,&
                  global_env=global_env,error=iError)
          end if
          call sm_matrix_transf_to_m(source_m=full_s(1)%matrix, &
               target_m=glob_angles%min_overlap_m(1)%matrix%matrix,&
               glob_angles=glob_angles,qs_env=qs_env,&
               global_env=global_env,error=error)
          glob_angles%m_s_m_valid=.true.
       end if
       CPAssert(associated(glob_angles%min_overlap_m),cp_failure_level,routineP,error,failure)
       if (glob_angles%m_s_m_valid) then
          min_overlap_m => glob_angles%min_overlap_m
       else
          nullify(min_overlap_m)
       end if
       !       end if
       if (present(min_hamiltonian_m)) then
          if (.not.glob_angles%m_h_m_valid) then
             if (.not.associated(glob_angles%min_hamiltonian_m)) then
                allocate(glob_angles%min_hamiltonian_m(1),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                call cp_sp_create(glob_angles%min_hamiltonian_m(1)%matrix, &
                     matrix_structure=min_bas_m_struct,&
                     global_env=global_env,error=iError)
                call cp_error_propagate_error(iError,routineP,&
                     CPSourceFileRef,&
                     error,failure)
             end if
             if (.not.failure) then
                call sm_matrix_transf_to_m(source_m=full_h(1)%matrix,&
                     target_m=glob_angles%min_hamiltonian_m(1)%matrix%matrix,&
                     glob_angles=glob_angles,qs_env=qs_env,&
                     global_env=global_env,error=error)
                glob_angles%m_h_m_valid=.true.
             end if
          end if
          CPAssert(associated(glob_angles%min_overlap_m),cp_failure_level,routineP,error,failure)
          if (glob_angles%m_h_m_valid) then
             min_hamiltonian_m => glob_angles%min_hamiltonian_m
          else
             nullify(min_hamiltonian_m)
          end if
       end if
    end if
    !call timestop(0.0_wp,handle)
  end subroutine pao_g_ang_get_min_m
!***************************************************************************

end module pao_proj_methods
