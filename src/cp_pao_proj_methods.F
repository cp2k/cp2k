!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_pao_proj_methods [1.0] *
!!
!!   NAME
!!     cp_pao_proj_methods
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
!!     none
!!
!!   SOURCE
!*****************************************************************************
module cp_pao_proj_methods
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use cp_pao_types
  use cp_matrix_utils 
  use cp_pao_utils
  use cp_pao_matrix_methods
  use cp_pao_l_angles_methods
  use sparse_matrix_types
  use message_passing, only: mp_sum
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_pao_proj_methods'

  public :: cp_get

! core procedures
!  public :: cp_init, cp_dealloc_ref, cp_set, cp_get, cp_create, cp_dealloc

! core functions
!  public :: cp_validate, cp_valid

! projection - injection
  public :: cp_pao_l_atom_projection, cp_pao_l_atom_injection, &
       cp_pao_r_atom_projection, cp_pao_r_atom_injection, &
       cp_pao_atom_m_to_f_matrix, cp_pao_atom_f_to_m_matrix, &
       cp_m2f, cp_f2m, cp_m2f_with_B_M_pi

  public :: cp_pao_g_ang_get_min_m, cp_pao_matrix_transf_to_m,&
       cp_pao_matrix_transf_to_f, cp_pao_m_m2f_with_b_m_pi

  interface cp_get
     module procedure cp_pao_g_ang_get_min_m
  end interface

  ! transfer from full basis to min basis
  interface cp_f2m
     module procedure cp_pao_matrix_transf_to_m
  end interface

  !transfer from min basis to full basis
  interface cp_m2f
     module procedure cp_pao_matrix_transf_to_f
  end interface
  
  ! tranfser to full basis with B matrix pi
  interface cp_m2f_with_B_M_pi
     module procedure cp_pao_m_m2f_with_b_m_pi
  end interface
!!***
!****************************************************************************

contains

! these routines WILL CHANGE, they are extremly ugly and inefficient
! I will come back to them as soo as the programs works
! to do

!!****f* cp_pao_projection/cp_pao_l_atom_projection [1.0] *
!!
!!   NAME
!!     cp_pao_l_atom_projection
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_l_atom_projection(pao_env, atom_nr, source_m,&
!!         target_m, error)
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Integer, Intent (IN):: atom_nr
!!       Real(Kind=wp), Dimension(:,:), Intent (IN):: source_m
!!       Real(Kind=wp), Dimension(:,:), Intent (OUT):: target_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_l_atom_projection
!!
!!   FUNCTION
!!     Left multiplication with the projection in the small basis
!!     in one atom block.
!!     Puts pi*source_m into target_m, where pi is the projection in the
!!     minimal basis.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao env where the projection is defined
!!     atom_nr: the actual atom block
!!     source_m: the source matrix
!!     target_m: the target matrix
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
  subroutine cp_pao_l_atom_projection(pao_env, atom_nr, source_m,&
       target_m, error)
    type(cp_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    real(kind=wp), dimension(:,:), intent(in) :: source_m
    real(kind=wp), dimension(:,:), intent(out) :: target_m
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_types:cp_pao_l_atom_projection'
    integer :: size_min_bas, size_full_bas, i
    failure=.false.

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       size_full_bas=cp_get_full_bas_size(pao_env,atom_nr=atom_nr,error=error)
       size_min_bas=cp_get_min_bas_size(pao_env,atom_nr=atom_nr,error=error)
       CPPrecondition(size(source_m,1)==size_full_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(target_m,1)==size_min_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(source_m,2)==size(target_m,2),cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       do i=1,size_min_bas
          target_m(i,:)=source_m(cp_get_p_index(pao_env,atom_nr=atom_nr,i=i,error=error),:)
! reorganization of the operations could probabily make it faster
       end do
    end if
  end subroutine cp_pao_l_atom_projection
!***************************************************************************

!!****f* cp_pao_projection/cp_pao_l_atom_injection [1.0] *
!!
!!   NAME
!!     cp_pao_l_atom_injection
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_l_atom_injection(pao_env, atom_nr, source_m,&
!!         target_m, error)
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Integer, Intent (IN):: atom_nr
!!       Real(Kind=wp), Dimension(:,:), Intent (IN):: source_m
!!       Real(Kind=wp), Dimension(:,:), Intent (OUT):: target_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_l_atom_injection
!!
!!   FUNCTION
!!     Left multiplication with the injection in the full basis
!!     in one atom block.
!!     Puts i*source_m into target_m, where i is the injection in the
!!     full basis.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao env where the projection is defined
!!     atom_nr: the actual atom block
!!     source_m: the source matrix
!!     target_m: the target matrix
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
  subroutine cp_pao_l_atom_injection(pao_env, atom_nr, source_m,&
       target_m, error)
    type(cp_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    real(kind=wp), dimension(:,:), intent(in) :: source_m
    real(kind=wp), dimension(:,:), intent(out) :: target_m
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_types:cp_pao_l_atom_injection'
    integer :: size_min_bas, size_full_bas, i
    failure=.false.

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       size_full_bas=cp_get_full_bas_size(pao_env,atom_nr=atom_nr,error=error)
       size_min_bas=cp_get_min_bas_size(pao_env,atom_nr=atom_nr,error=error)
       CPPrecondition(size(source_m,1)==size_min_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(target_m,1)==size_full_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(source_m,2)==size(target_m,2),cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       target_m=0.0_wp
       do i=1,size_min_bas
          target_m(cp_get_p_index(pao_env,atom_nr=atom_nr,i=i,error=error),:)=&
               source_m(i,:)
! reorganization of the operations could probabily make it faster
! and do not set target_m to 0 if it gets overwritten
       end do
    end if
  end subroutine cp_pao_l_atom_injection
!***************************************************************************

!!****f* cp_pao_projection/cp_pao_r_atom_projection [1.0] *
!!
!!   NAME
!!     cp_pao_r_atom_projection
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_r_atom_projection(pao_env, atom_nr, source_m,&
!!         target_m, error)
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Integer, Intent (IN):: atom_nr
!!       Real(Kind=wp), Dimension(:,:), Intent (IN):: source_m
!!       Real(Kind=wp), Dimension(:,:), Intent (OUT):: target_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_r_atom_projection
!!
!!   FUNCTION
!!     Right multiplication with the projection in the small basis
!!     in one atom block.
!!     Puts source_m*pi into target_m, where pi is the projection in the
!!     minimal basis.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao env where the projection is defined
!!     atom_nr: the actual atom block
!!     source_m: the source matrix
!!     target_m: the target matrix
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
  subroutine cp_pao_r_atom_projection(pao_env, atom_nr, source_m,&
       target_m, error)
    type(cp_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    real(kind=wp), dimension(:,:), intent(in) :: source_m
    real(kind=wp), dimension(:,:), intent(out) :: target_m
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_types:cp_pao_r_atom_projection'
    integer :: size_min_bas, size_full_bas, i
    failure=.false.

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       size_full_bas=cp_get_full_bas_size(pao_env,atom_nr=atom_nr,error=error)
       size_min_bas=cp_get_min_bas_size(pao_env,atom_nr=atom_nr,error=error)
       CPPrecondition(size(source_m,2)==size_min_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(target_m,2)==size_full_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(source_m,1)==size(target_m,1),cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       target_m=0.0_wp
       do i=1,size_min_bas
          target_m(:,cp_get_p_index(pao_env,atom_nr=atom_nr,i=i,error=error))=&
               source_m(:,i)
! reorganization of the operations could probabily make it faster
       end do
    end if
  end subroutine cp_pao_r_atom_projection
!***************************************************************************

!!****f* cp_pao_projection/cp_pao_r_atom_injection [1.0] *
!!
!!   NAME
!!     cp_pao_r_atom_injection
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_r_atom_injection(pao_env, atom_nr, source_m,&
!!         target_m, error)
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Integer, Intent (IN):: atom_nr
!!       Real(Kind=wp), Dimension(:,:), Intent (IN):: source_m
!!       Real(Kind=wp), Dimension(:,:), Intent (OUT):: target_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_r_atom_injection
!!
!!   FUNCTION
!!     Right multiplication with the injection in the full basis
!!     in one atom block.
!!     Puts source_m*i into target, where i is the injection in the
!!     full basis.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao env where the projection is defined
!!     atom_nr: the actual atom block
!!     source_m: the source matrix
!!     target_m: the target matrix
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
  subroutine cp_pao_r_atom_injection(pao_env, atom_nr, source_m,&
       target_m, error)
    type(cp_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    real(kind=wp), dimension(:,:), intent(in) :: source_m
    real(kind=wp), dimension(:,:), intent(out) :: target_m
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_types:cp_pao_r_atom_injection'
    integer :: size_min_bas, size_full_bas, i
    type(cp_error_type) :: iError
    failure=.false.
    call cp_init(iError,template_error=error)

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       size_full_bas=cp_get_full_bas_size(pao_env,atom_nr=atom_nr,error=iError)
       CPAssert(cp_error_get_level(iError)==0,cp_error_get_level(iError),routineP,error,failure)
       size_min_bas=cp_get_min_bas_size(pao_env,atom_nr=atom_nr,error=iError)
       CPAssert(cp_error_get_level(iError)==0,cp_error_get_level(iError),routineP,error,failure)
       CPPrecondition(size(source_m,2)==size_full_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(target_m,2)==size_min_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(source_m,1)==size(target_m,1),cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       do i=1,size_min_bas
          target_m(:,i) = source_m(:,&
               cp_get_p_index(pao_env,atom_nr=atom_nr,i=i,error=error))
! reorganization of the operations could probabily make it faster
       end do
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine cp_pao_r_atom_injection
!***************************************************************************

!!****f* cp_pao_types/cp_pao_atom_m_to_f_matrix [1.0] *
!!
!!   NAME
!!     cp_pao_atom_m_to_f_matrix
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_atom_m_to_f_matrix(pao_env, l_atom_nr, r_atom_nr,&
!!         source_m, target_m, error)
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Integer, Intent (IN):: r_atom_nr, l_atom_nr
!!       Real(Kind=wp), Dimension(:,:), Intent (IN):: source_m
!!       Real(Kind=wp), Dimension(:,:), Intent (OUT):: target_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_atom_m_to_f_matrix
!!
!!   FUNCTION
!!     inject an atomic matrix block into the full basis:
!!     puts i*source_m*pi into target_m, where i the injection, and 
!!     pi(=i^*) is the projection in the small basis.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao env where the projection is defined
!!     l_atom_nr: the left atom block
!!            (pao_env%nr_min_bas(l_atom_nr)= rows of source_m )
!!     r_atom_nr: the right atom block
!!            (pao_env%nr_min_bas(r_atom_nr)= colss of source_m )
!!     source_m: the source matrix
!!     target_m: the target matrix
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
  subroutine cp_pao_atom_m_to_f_matrix(pao_env,l_atom_nr,r_atom_nr,&
       source_m, target_m, error)
    type(cp_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: r_atom_nr, l_atom_nr
    real(kind=wp), dimension(:,:), intent(in) :: source_m
    real(kind=wp), dimension(:,:), intent(out) :: target_m
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_types:cp_pao_atom_m_to_f_matrix'
    integer :: size_l_min_bas, size_l_full_bas, size_r_min_bas, size_r_full_bas,&
         i,j
    integer, dimension(:), pointer :: l_index, r_index
    type(cp_error_type) :: iError
    failure=.false.
    nullify(l_index,r_index)
    call cp_init(iError,template_error=error)

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       size_l_full_bas=cp_get_full_bas_size(pao_env,atom_nr=l_atom_nr,&
            error=iError)
       CPAssert(cp_error_get_level(iError)==0,cp_error_get_level(iError),routineP,error,failure)
       size_l_min_bas=cp_get_min_bas_size(pao_env,atom_nr=l_atom_nr,&
            error=iError)
       CPAssert(cp_error_get_level(iError)==0,cp_error_get_level(iError),routineP,error,failure)
       size_r_full_bas=cp_get_full_bas_size(pao_env,atom_nr=r_atom_nr,error=iError)
       CPAssert(cp_error_get_level(iError)==0,cp_error_get_level(iError),routineP,error,failure)
       size_r_min_bas=cp_get_min_bas_size(pao_env,atom_nr=r_atom_nr,error=iError)
       CPAssert(cp_error_get_level(iError)==0,cp_error_get_level(iError),routineP,error,failure)
       CPPrecondition(size(source_m,1)==size_l_min_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(source_m,2)==size_r_min_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(target_m,1)==size_l_full_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(target_m,2)==size_r_full_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(pao_env%projection),cp_failure_level,routineP,error,failure)
! check projection more?
    end if
    if (.not. failure) then
       target_m=0.0_wp
       l_index => pao_env%projection(l_atom_nr)%array
       r_index => pao_env%projection(r_atom_nr)%array
       do i=1,size_l_min_bas
          do j=1,size_r_min_bas
             target_m(l_index(i),r_index(j))=source_m(i,j)
          end do
       end do
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine cp_pao_atom_m_to_f_matrix
!***************************************************************************

!!****f* cp_pao_types/cp_pao_atom_f_to_m_matrix [1.0] *
!!
!!   NAME
!!     cp_pao_atom_f_to_m_matrix
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_atom_f_to_m_matrix(pao_env, l_atom_nr, r_atom_nr,&
!!         source_m, target_m, error)
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Integer, Intent (IN):: r_atom_nr, l_atom_nr
!!       Real(Kind=wp), Dimension(:,:), Intent (IN):: source_m
!!       Real(Kind=wp), Dimension(:,:), Intent (OUT):: target_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_atom_f_to_m_matrix
!!
!!   FUNCTION
!!     project an atomic matrix block into the minimal basis:
!!     puts pi*source_m*i into target_m, where i the injection, and 
!!     pi(=i^*) is the projection in the small basis.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao env where the projection is defined
!!     l_atom_nr: the left atom block 
!!            (pao_env%nr_full_bas(l_atom_nr)= rows of source_m )
!!     r_atom_nr: the right atom block
!!            (pao_env%nr_full_bas(r_atom_nr)= colss of source_m )
!!     source_m: the source matrix
!!     target_m: the target matrix
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
  subroutine cp_pao_atom_f_to_m_matrix(pao_env,l_atom_nr,r_atom_nr,&
       source_m, target_m, error)
    type(cp_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: r_atom_nr, l_atom_nr
    real(kind=wp), dimension(:,:), intent(in) :: source_m
    real(kind=wp), dimension(:,:), intent(out) :: target_m
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_types:cp_pao_atom_f_to_m_matrix'
    integer :: size_l_min_bas, size_l_full_bas, size_r_min_bas, size_r_full_bas,&
         i,j
    integer, dimension(:), pointer :: l_index, r_index
    type(cp_error_type) :: iError
    failure=.false.
    nullify(l_index,r_index)
    call cp_init(iError,template_error=error)

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       size_l_full_bas=cp_get_full_bas_size(pao_env,atom_nr=l_atom_nr,&
            error=iError)
       CPAssert(cp_error_get_level(iError)==0,cp_error_get_level(iError),routineP,error,failure)
       size_l_min_bas=cp_get_min_bas_size(pao_env,atom_nr=l_atom_nr,&
            error=iError)
       CPAssert(cp_error_get_level(iError)==0,cp_error_get_level(iError),routineP,error,failure)
       size_r_full_bas=cp_get_full_bas_size(pao_env,atom_nr=r_atom_nr,error=iError)
       CPAssert(cp_error_get_level(iError)==0,cp_error_get_level(iError),routineP,error,failure)
       size_r_min_bas=cp_get_min_bas_size(pao_env,atom_nr=r_atom_nr,error=iError)
       CPAssert(cp_error_get_level(iError)==0,cp_error_get_level(iError),routineP,error,failure)
       CPPrecondition(size(source_m,1)==size_l_min_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(source_m,2)==size_r_min_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(target_m,1)==size_l_full_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(size(target_m,2)==size_r_full_bas,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(pao_env%projection),cp_failure_level,routineP,error,failure)
! check projection more?
    end if
    if (.not. failure) then
       l_index => pao_env%projection(l_atom_nr)%array
       r_index => pao_env%projection(r_atom_nr)%array
       do i=1,size_l_min_bas
          do j=1,size_r_min_bas
             target_m(i,j)=source_m(l_index(i),r_index(j))
          end do
       end do
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine cp_pao_atom_f_to_m_matrix
!***************************************************************************

!!****f* cp_pao_types/cp_pao_m_m2f_with_b_m_pi [1.0] *
!!
!!   NAME
!!     cp_pao_m_m2f_with_b_m_pi
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_m_m2f_with_b_m_pi(source_m, target_m, glob_angles,&
!!         error)
!!       Type(cp_pao_matrix_type), Intent (IN):: source_m
!!       Type(cp_pao_matrix_type), Intent (INOUT):: target_m
!!       Type(cp_pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_m_m2f_with_b_m_pi
!!
!!   FUNCTION
!!     Transform the given matrix to the full basis with B=i U and pi=i^*
!!     Puts U i source_m pi  into the target_m
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
  subroutine cp_pao_m_m2f_with_b_m_pi(source_m, target_m, glob_angles, error)
    type(cp_pao_matrix_type), intent(in) :: source_m
    type(cp_pao_matrix_type), intent(inout) :: target_m
    type(cp_pao_glob_angles_type), intent(in) :: glob_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_types:cp_pao_m_m2f_with_b_m_pi'
    type(cp_matrix_block_iterator) :: iterator
    integer :: brow,bcol, stat, max_bsize,i,j
    integer, dimension(:), pointer :: l_index,r_index
    real(kind=wp), dimension(:,:), pointer :: block_val, u_col, u_row
    real(kind=wp),allocatable, dimension(:,:) :: t1_m
    type(cp_pao_local_angles_type) :: l_angles
    type(cp_error_type) :: iError
    type(cp_pao_env_type), pointer :: pao_env
    call cp_init(iError,template_error=error)
    nullify(block_val)
    failure=.false.

    CPPrecondition(source_m%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(target_m%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       CPPrecondition(source_m%min_bas_cols,cp_failure_level,routineP,error,failure)
       CPPrecondition(source_m%min_bas_rows,cp_failure_level,routineP,error,failure)
       CPPrecondition(.not.target_m%min_bas_cols,cp_failure_level,routineP,error,failure)
       CPPrecondition(.not.target_m%min_bas_rows,cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       pao_env => cp_get_pao_env(source_m)
       CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       CPPrecondition(cp_valid(pao_env),cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       call cp_init(l_angles,atom_nr=1, &
            glob_angles=glob_angles, error=error)
       max_bsize=maxval(pao_env%nr_full_bas)
       allocate(t1_m(max_bsize,max_bsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       call cp_init(iterator,source_m%matrix,error=iError)
       call cp_error_propagate_error(iError,routineP,&
            CPSourceFileRef,&
            error,failure)
       if (.not.failure) then
          blocksDo: do
             if ((.not.cp_next(iterator,error=error)).or.failure) exit
             call cp_get(iterator,block_row=brow,block_col=bcol,&
                  block_val=block_val)
! get u
             call cp_set(l_angles, atom_nr=bcol,error=error)
             call cp_get(l_angles,unitary_t=u_col)
             call cp_set(l_angles, atom_nr=brow,error=error)
             call cp_get(l_angles,unitary_t=u_row)

             ! t1_m = i block_val pi
             t1_m=0.0_wp
             l_index => pao_env%projection(bcol)%array
             r_index => pao_env%projection(brow)%array
             CPPrecondition(associated(l_index),cp_failure_level,routineP,error,failure)
             CPPrecondition(associated(r_index),cp_failure_level,routineP,error,failure)
             if (.not.failure) then
                call cp_assert(size(l_index)==size(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                call cp_assert(size(r_index)==size(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             end if
             if (.not.failure) then ! put result in the small basis block
                do i=1,size(block_val,1)
                   do j=1,size(block_val,2)
                      t1_m(l_index(i),r_index(j))=block_val(i,j)
                   end do
                end do
             end if


! size(u_row,1)==size(u_row,2)==size(target_m,1)
! size(u_col,1)==size(u_col,2)==size(target_m,2)

! get full basis block
             call get_block_node(target_m%matrix, block_col=bcol,&
                  block_row=bcol, block=block_val)
             if (.not.associated(block_val)) then
                call add_block_node(target_m%matrix, block_col=bcol,&
                     block_row=bcol, block=block_val)
                CPPostcondition(associated(block_val),cp_failure_level,routineP,error,failure)
             end if
! block_val = u t1_m
             if (.not.failure) then
                call dgemm('N','N',size(u_row,1),size(u_row,2),&
                     size(u_col,2),1.0_wp,&
                     u_row,size(u_row,1),&
                     block_val,size(block_val,1),&
                     0.0_wp,t1_m,size(t1_m,1))
             end if

          end do blocksDo
       end if
       call cp_dealloc_ref(iterator,error=error)
       call cp_dealloc_ref(l_angles,error=error)
       deallocate(t1_m,stat=stat) ! just to be sure...
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine cp_pao_m_m2f_with_b_m_pi
!***************************************************************************


!!****f* cp_pao_types/cp_pao_matrix_transf_to_f [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_transf_to_f
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_matrix_transf_to_f(source_m, target_m, glob_angles,&
!!         error)
!!       Type(cp_pao_matrix_type), Intent (IN):: source_m
!!       Type(cp_pao_matrix_type), Intent (INOUT):: target_m
!!       Type(cp_pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_matrix_transf_to_f
!!
!!   FUNCTION
!!     Transform the given matrix to the full basis.
!!     Puts U i source_m pi U^*  into the target_m
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
  subroutine cp_pao_matrix_transf_to_f(source_m, target_m, glob_angles, error)
    type(cp_pao_matrix_type), intent(in) :: source_m
    type(cp_pao_matrix_type), intent(inout) :: target_m
    type(cp_pao_glob_angles_type), intent(in) :: glob_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_types:cp_pao_matrix_transf_to_f'
    type(cp_matrix_block_iterator) :: iterator
    integer :: brow,bcol, stat, max_bsize,i,j
    integer, dimension(:), pointer :: l_index,r_index
    real(kind=wp), dimension(:,:), pointer :: block_val, u_col, u_row
    real(kind=wp),allocatable, dimension(:,:) :: t1_m,t2_m
    type(cp_pao_local_angles_type) :: l_angles
    type(cp_error_type) :: iError
    type(cp_pao_env_type), pointer :: pao_env
    call cp_init(iError,template_error=error)
    nullify(block_val)
    failure=.false.

    CPPrecondition(source_m%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(target_m%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       CPPrecondition(source_m%min_bas_cols,cp_failure_level,routineP,error,failure)
       CPPrecondition(source_m%min_bas_rows,cp_failure_level,routineP,error,failure)
       CPPrecondition(.not.target_m%min_bas_cols,cp_failure_level,routineP,error,failure)
       CPPrecondition(.not.target_m%min_bas_rows,cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       pao_env => cp_get_pao_env(source_m)
       CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       CPPrecondition(cp_valid(pao_env),cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       call cp_init(l_angles,atom_nr=1, &
            glob_angles=glob_angles, error=error)
       max_bsize=maxval(pao_env%nr_full_bas)
       allocate(t1_m(max_bsize,max_bsize), t2_m(max_bsize,max_bsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       call cp_init(iterator,source_m%matrix,error=iError)
       call cp_error_propagate_error(iError,routineP,&
            CPSourceFileRef,&
            error,failure)
       if (.not.failure) then
          blocksDo: do
             if ((.not.cp_next(iterator,error=error)).or.failure) exit
             call cp_get(iterator,block_row=brow,block_col=bcol,&
                  block_val=block_val)
! get u
             call cp_set(l_angles, atom_nr=bcol,error=error)
             call cp_get(l_angles,unitary_t=u_col)
             call cp_set(l_angles, atom_nr=brow,error=error)
             call cp_get(l_angles,unitary_t=u_row)

             ! t2_m = i block_val pi
             t2_m=0.0_wp
             l_index => pao_env%projection(bcol)%array
             r_index => pao_env%projection(brow)%array
             CPPrecondition(associated(l_index),cp_failure_level,routineP,error,failure)
             CPPrecondition(associated(r_index),cp_failure_level,routineP,error,failure)
             if (.not.failure) then
                call cp_assert(size(l_index)==size(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                call cp_assert(size(r_index)==size(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             end if
             if (.not.failure) then ! put result in the small basis block
                do i=1,size(block_val,1)
                   do j=1,size(block_val,2)
                      t2_m(l_index(i),r_index(j))=block_val(i,j)
                   end do
                end do
             end if


! size(u_row,1)==size(u_row,2)==size(target_m,1)
! size(u_col,1)==size(u_col,2)==size(target_m,2)

! t1_m = t2_m u^*
             call dgemm('N','T',size(u_row,1),size(u_col,1),&
                  size(u_col,2),1.0_wp,&
                  t2_m,size(t2_m,1),&
                  u_col,size(u_col,1),&
                  0.0_wp,t1_m,size(t1_m,1))
! get full basis block
             call get_block_node(target_m%matrix, block_col=bcol,&
                  block_row=bcol, block=block_val)
             if (.not.associated(block_val)) then
                call add_block_node(target_m%matrix, block_col=bcol,&
                     block_row=bcol, block=block_val)
                CPPostcondition(associated(block_val),cp_failure_level,routineP,error,failure)
             end if
! block_val = u t1_m
             if (.not.failure) then
                call dgemm('N','N',size(u_row,1),size(u_row,2),&
                     size(u_col,2),1.0_wp,&
                     u_row,size(u_row,1),&
                     block_val,size(block_val,1),&
                     0.0_wp,t2_m,size(t2_m,1))
             end if

          end do blocksDo
       end if
       call cp_dealloc_ref(iterator,error=error)
       call cp_dealloc_ref(l_angles,error=error)
       deallocate(t1_m,t2_m,stat=stat) ! just to be sure...
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine cp_pao_matrix_transf_to_f
!***************************************************************************

!!****f* cp_pao_types/cp_pao_matrix_transf_to_m [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_transf_to_m
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_matrix_transf_to_m(source_m, target_m, glob_angles,&
!!         error)
!!       Type(cp_pao_matrix_type), Intent (IN):: source_m
!!       Type(cp_pao_matrix_type), Intent (INOUT):: target_m
!!       Type(cp_pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_matrix_transf_to_m
!!
!!   FUNCTION
!!     Transform the given matrix to the minimal basis.
!!     Puts pi U^* source_m U i into the target_m
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
  subroutine cp_pao_matrix_transf_to_m(source_m, target_m, glob_angles, error)
    type(cp_pao_matrix_type), intent(in) :: source_m
    type(cp_pao_matrix_type), intent(inout) :: target_m
    type(cp_pao_glob_angles_type), intent(in) :: glob_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_types:cp_pao_matrix_transf_to_m'
    type(cp_matrix_block_iterator) :: iterator
    integer :: brow,bcol, stat, max_bsize,i,j
    integer, dimension(:), pointer :: l_index,r_index
    real(kind=wp), dimension(:,:), pointer :: block_val, u_col, u_row
    real(kind=wp),allocatable, dimension(:,:) :: t1_m,t2_m
    type(cp_pao_local_angles_type) :: l_angles
    type(cp_error_type) :: iError
    type(cp_pao_env_type), pointer :: pao_env
    call cp_init(iError,template_error=error)
    nullify(block_val)
    failure=.false.

    CPPrecondition(source_m%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(target_m%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       CPPrecondition(.not.source_m%min_bas_cols,cp_failure_level,routineP,error,failure)
       CPPrecondition(.not.source_m%min_bas_rows,cp_failure_level,routineP,error,failure)
       CPPrecondition(target_m%min_bas_cols,cp_failure_level,routineP,error,failure)
       CPPrecondition(target_m%min_bas_rows,cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       pao_env => cp_get_pao_env(source_m)
       CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       CPPrecondition(associated(pao_env%nr_full_bas),cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(pao_env%projection),cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       call cp_init(l_angles,atom_nr=1, &
            glob_angles=glob_angles, error=error)
       max_bsize=maxval(pao_env%nr_full_bas)
       allocate(t1_m(max_bsize,max_bsize), t2_m(max_bsize,max_bsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       call cp_init(iterator,source_m%matrix,error=iError)
       call cp_error_propagate_error(iError,routineP,&
            CPSourceFileRef,&
            error,failure)
       if (.not.failure) then
          blocksDo: do
             if ((.not.cp_next(iterator,error=error)).or.failure) exit
             call cp_get(iterator,block_row=brow,block_col=bcol,&
                  block_val=block_val)
! get u
             call cp_set(l_angles, atom_nr=bcol,error=error)
             call cp_get(l_angles,unitary_t=u_col)
             call cp_set(l_angles, atom_nr=brow,error=error)
             call cp_get(l_angles,unitary_t=u_row)

! in the following the sizes are put to make sense, i.e. not using
! size(u_row,1)==size(u_row,2)==size(source_m,1)
! size(u_col,1)==size(u_col,2)==size(source_m,2)

! t1_m = source u
             call dgemm('N','N',size(block_val,1),size(block_val,2),&
                  size(u_col,2),1.0_wp,&
                  block_val,size(block_val,1),&
                  u_col,size(u_col,1),&
                  0.0_wp,t1_m,size(t1_m,1))
! t2_m = u^* t1_m
             call dgemm('T','N',size(u_row,1),size(u_row,2),&
                  size(u_col,2),1.0_wp,&
                  u_row,size(u_row,1),&
                  t1_m,size(t1_m,1),&
                  0.0_wp,t2_m,size(t2_m,1))
! get small basis block
             call get_block_node(target_m%matrix, block_col=bcol,&
                  block_row=bcol, block=block_val)
             if (.not.associated(block_val)) then
                call add_block_node(target_m%matrix, block_col=bcol,&
                     block_row=bcol, block=block_val)
                CPPostcondition(associated(block_val),cp_failure_level,routineP,error,failure)
             end if

             l_index => pao_env%projection(bcol)%array
             r_index => pao_env%projection(brow)%array
             CPPrecondition(associated(l_index),cp_failure_level,routineP,error,failure)
             CPPrecondition(associated(r_index),cp_failure_level,routineP,error,failure)
             if (.not.failure) then
                call cp_assert(size(l_index)==size(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                call cp_assert(size(r_index)==size(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             end if
             if (.not.failure) then ! put result in the small basis block
                do i=1,size(block_val,1)
                   do j=1,size(block_val,2)
                      block_val(i,j)=t2_m(l_index(i),r_index(j))
                   end do
                end do
             end if

          end do blocksDo
       end if
       call cp_dealloc_ref(iterator,error=error)
       call cp_dealloc_ref(l_angles,error=error)
       deallocate(t1_m,t2_m,stat=stat) ! just to be sure...
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine cp_pao_matrix_transf_to_m
!***************************************************************************

!!****f* cp_pao_proj_methods/cp_pao_g_ang_get_min_m [1.0] *
!!
!!   NAME
!!     cp_pao_g_ang_get_min_m
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_g_ang_get_min_m(glob_angles, min_overlap_m,&
!!         min_hamiltonian_m, error)
!!       Type(cp_pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Type(cp_pao_matrix_type), Pointer, Optional:: min_hamiltonian_m
!!       Type(cp_pao_matrix_type), Pointer:: min_overlap_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_g_ang_get_min_m
!!
!!   FUNCTION
!!     returns the minimal matrices
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
subroutine cp_pao_g_ang_get_min_m(glob_angles,&
     min_overlap_m,  min_hamiltonian_m, error)
  type(cp_pao_glob_angles_type), intent(inout) :: glob_angles
  type(cp_pao_matrix_type), pointer, optional :: min_hamiltonian_m
  type(cp_pao_matrix_type), pointer ::  min_overlap_m
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: handle
  character(len=*), parameter :: routineN='cp_pao_g_ang_get_min_m',&
       routineP=moduleN//':'//routineN
  type(cp_error_type) :: iError
  call cp_init(iError, template_error=error)
  failure=.false.
  
  !call timeset(routineN//','//moduleN,'I',"",handle)
  CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
!       if (present(min_overlap_m)) then
          if (.not.glob_angles%m_s_m_valid) then
             if (.not.associated(glob_angles%min_overlap_m)) then
                call cp_create(glob_angles%min_overlap_m, &
                     min_bas_cols=.true., min_bas_rows=.true.,&
                     pao_env=glob_angles%pao_env)
             end if
             call cp_pao_matrix_transf_to_m(source_m=&
                  glob_angles%pao_env%full_overlap_m, &
                  target_m=glob_angles%min_overlap_m,&
                  glob_angles=glob_angles,error=error)
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
                call cp_create(glob_angles%min_hamiltonian_m, &
                     min_bas_cols=.true., min_bas_rows=.true.,&
                     pao_env=glob_angles%pao_env,error=iError)
                call cp_error_propagate_error(iError,routineP,&
                     CPSourceFileRef,&
                     error,failure)
             end if
             if (.not.failure) then
                call cp_pao_matrix_transf_to_m(source_m=&
                     glob_angles%pao_env%full_hamiltonian_m, &
                     target_m=glob_angles%min_hamiltonian_m,&
                     glob_angles=glob_angles,error=error)
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
end subroutine cp_pao_g_ang_get_min_m
!***************************************************************************

end module cp_pao_proj_methods
