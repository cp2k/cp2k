!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_kpp1_env_types [1.0] *
!!
!!   NAME
!!     qs_kpp1_env_types
!!
!!   FUNCTION
!!     basis types for the calculation of the perturbation of density theory.
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!   SOURCE
! ****************************************************************************
MODULE qs_kpp1_env_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
  USE cp_b_matrix_structure,           ONLY: cp_b_matrix_struct_type
  USE cp_block_matrix,                 ONLY: cp_block_matrix_p_type,&
                                             cp_block_matrix_type,&
                                             cp_sp_array_create,&
                                             cp_sp_array_dealloc,&
                                             cp_sp_create,&
                                             cp_sp_multiplies_blacs,&
                                             cp_sp_set_to
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_fm_struct,                    ONLY: cp_fm_struct_get,&
                                             cp_fm_struct_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level,&
                                             cp_logger_generate_filename
  USE cp_matrix_utils,                 ONLY: cp_sm_output
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_control_types,                       ONLY: dft_control_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             dp
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  use qs_parser, only: open_file, close_file
  USE realspace_grid_types,            ONLY: rs_pw_to_cube
  use cp_array_i_utils, only: cp_2d_i_output, cp_2d_i_write
  use cp_array_r_utils, only: cp_2d_r_output, cp_2d_r_write
  use cp_array_logical_utils, only: cp_2d_logical_output,&
       cp_1d_logical_output, cp_2d_logical_write, cp_1d_logical_write
  use string_utilities, only: compress
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             deallocate_matrix_set,&
                                             get_matrix_info,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,             ONLY: xc_derivative_set_type, &
                                             xc_dset_create, &
                                             xc_dset_release
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: qs_kpp1_env_type, qs_kpp1_env_p_type
  PUBLIC :: kpp1_release, kpp1_retain, kpp1_write

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_kpp1_env_types'

!!***
!****************************************************************************

!!****s* qs_p_build_kernel/qs_kpp1_env_type [1.0] *
!!
!!   NAME
!!     qs_kpp1_env_type
!!
!!   FUNCTION
!!     environement that keeps the informations and temporary
!!     val to build the kpp1 kernel matrix
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - ref_count: reference count (how may objects are sharing this one)
!!     - v_rspace: potential in r space. This is used to do an update only
!!       of what has changed. Useful if P1 converges to some density
!!       (you spare in the grid-ao conversion).
!!       With LSD contains the various orbitals.
!!     - v_ao: the potential in the ao basis (used togheter with v_rspace
!!       to update only what changed
!!     - rebuild_each: how often a full rebuild should be performed
!!     - updates_since_rebuild: how many updates have been done since the
!!       last rebuild (if -1, forcces the rebuild)
!!     - id_nr: identification number, unique for each kpp1 env
!!     - print_count: counter to create unique filename
!!     - iter: number of iterations
!!     - drho_r(idir,ispin): the derivative of rho wrt. x,y,z in the real space
!!     - deriv_xc(ii,ipot): the second derivative of the xc potential at psi0 
!!       (qs_env%c), if grad pot is true it should already be divised
!!       by the gradient
!!     - spin_pot(1:2,ipot): information about wrt. to which spins the
!!        corresponding component of deriv_xc was derived (see
!!        xc_create_2nd_deriv_info)
!!     - grad_pot(1:2,ipot): if the derivative spin_pot was wrt. to
!!       the gradient (see xc_create_2nd_deriv_info)
!!     - ndiag_term(ipot): it the term is an off diagonal term (see 
!!       xc_create_2nd_deriv_info)
!!
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_kpp1_env_type
     INTEGER :: ref_count, rebuild_each, updates_since_rebuild, id_nr,&
          print_count, iter
     TYPE(coeff_type), DIMENSION(:), POINTER :: v_rspace
     TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: v_ao
     TYPE(pw_p_type), DIMENSION(:,:), pointer :: drho_r
     !     REAL ( KIND = dp ) , DIMENSION (:,:), pointer :: deriv_xc
     type(xc_derivative_set_type), pointer :: derivative_set
     INTEGER, DIMENSION(:,:), POINTER :: spin_pot
     LOGICAL, DIMENSION(:,:), POINTER :: grad_pot
     LOGICAL, DIMENSION(:), POINTER :: ndiag_term
  END TYPE qs_kpp1_env_type
!!***
  !****************************************************************************

!!****s* qs_p_build_kernel/qs_kpp1_env_p_type [1.0] *
!!
!!   NAME
!!     qs_kpp1_env_p_type
!!
!!   FUNCTION
!!     just to build array of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     kpp1_env: the pointer to the kpp1_env
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_kpp1_env_p_type
     TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
  END TYPE qs_kpp1_env_p_type
!!***
  !****************************************************************************

CONTAINS

!!****f* qs_kpp1_env_types/kpp1_release [1.0] *
!!
!!   NAME
!!     kpp1_release
!!
!!   SYNOPSIS
!!     Subroutine kpp1_release(kpp1_env, error)
!!       Type(qs_kpp1_env_type), Pointer:: kpp1_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine kpp1_release
!!
!!   FUNCTION
!!     releases a kpp1_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kpp1_env: the environement to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE kpp1_release(kpp1_env,error)
    TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_release',&
         routineP=moduleN//':'//routineN
    INTEGER :: ispin, stat
    failure=.FALSE.

    IF (ASSOCIATED(kpp1_env)) THEN
       CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
       kpp1_env%ref_count=kpp1_env%ref_count-1
       IF (kpp1_env%ref_count<1) THEN
          IF (ASSOCIATED(kpp1_env%v_rspace)) THEN
             DO ispin=1,SIZE(kpp1_env%v_rspace)
                CALL coeff_deallocate(kpp1_env%v_rspace(ispin))
             END DO
             DEALLOCATE(kpp1_env%v_rspace,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(kpp1_env%v_ao)) THEN
             DO ispin=1,SIZE(kpp1_env%v_ao)
                IF (ASSOCIATED(kpp1_env%v_ao(ispin)%matrix)) THEN
                   CALL deallocate_matrix(kpp1_env%v_ao(ispin)%matrix)
                END IF
             END DO
             DEALLOCATE(kpp1_env%v_ao, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          if (associated(kpp1_env%drho_r)) then
             deallocate(kpp1_env%drho_r, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          end if
          if (associated(kpp1_env%derivative_set)) then
             call xc_dset_release(kpp1_env%derivative_set, error=error)
             nullify(kpp1_env%derivative_set)
          end if
          if (associated(kpp1_env%spin_pot)) then
             deallocate(kpp1_env%spin_pot, stat=stat)
          end if
          if (associated(kpp1_env%grad_pot)) then
             deallocate(kpp1_env%grad_pot, stat=stat)
          end if
          if (associated(kpp1_env%ndiag_term)) then
             deallocate(kpp1_env%ndiag_term, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          end if
          deallocate(kpp1_env, stat=stat)
       END IF
    END IF
    NULLIFY(kpp1_env)
  END SUBROUTINE kpp1_release
!***************************************************************************

!!****f* qs_kpp1_env_types/kpp1_retain [1.0] *
!!
!!   NAME
!!     kpp1_retain
!!
!!   SYNOPSIS
!!     Subroutine kpp1_retain(kpp1_env, error)
!!       Type(qs_kpp1_env_type), Pointer:: kpp1_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine kpp1_retain
!!
!!   FUNCTION
!!     retains a kpp1_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kpp1_env: the environement to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE kpp1_retain(kpp1_env,error)
    TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_retain',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
       kpp1_env%ref_count=kpp1_env%ref_count+1
    END IF
  END SUBROUTINE kpp1_retain
!***************************************************************************

!!****f* qs_kpp1_env_types/kpp1_write [1.0] *
!!
!!   NAME
!!     kpp1_write
!!
!!   FUNCTION
!!     writes the content of the kpp1_env to the given output unit
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - kpp1_env: the kpp1_env to output
!!     - unit_nr: the unit where to output
!!     - long_description: if a long description should be written
!!       (defaults to false)
!!     - local: if the unit is a local unit or a global unit 
!!       (defaults to false, i.e. global)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE kpp1_write(kpp1_env, unit_nr, long_description, local,&
       error)
    TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
    INTEGER, INTENT(in) :: unit_nr
    LOGICAL, OPTIONAL, INTENT(in) :: long_description, local
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, my_local, my_long_description, should_w
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_write',&
         routineP=moduleN//':'//routineN
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(cp_para_env_type), POINTER :: para_env
    character(len=30) :: base_filename
    character(len=60) :: filename
    integer :: ispin, idir, file_unit,my_unit_nr

    failure=.FALSE.; my_local=.FALSE.; my_long_description=.FALSE.

    IF (PRESENT(local)) my_local=local
    IF (PRESENT(long_description)) my_long_description=long_description
    logger=>cp_error_get_logger(error) 
    para_env=>logger%para_env
    should_w=my_local .OR. para_env%mepos==para_env%source
    my_unit_nr=unit_nr

    IF (ASSOCIATED(kpp1_env)) THEN
       kpp1_env%print_count=kpp1_env%print_count+1

       IF (should_w) THEN
          base_filename="kpp1-"//trim(adjustl(cp_to_string(kpp1_env%id_nr)))&
               //"-"//trim(adjustl(cp_to_string(kpp1_env%print_count)))
          call compress(base_filename,full=.TRUE.)

          if (my_long_description) then
             WRITE (unit=unit_nr,fmt="(a,a,a)") " wrote kpp1 into to files ",&
                  trim(base_filename),"*"
             call cp_logger_generate_filename(logger=logger,res=filename,&
                  root=trim(base_filename), postfix=".desc",&
                  local=my_local)
             CALL open_file(file_name=trim(filename),file_status="UNKNOWN",&
                  file_action="WRITE",unit_number=my_unit_nr)
          end if

          WRITE (unit=my_unit_nr,&
               fmt="(' <kpp1_env>:{ id_nr=',i10,' ref_count=',i10,',')")&
               kpp1_env%id_nr, kpp1_env%ref_count

          if (associated(kpp1_env%v_rspace)) then
             if (my_long_description) then
                do ispin=1,size(kpp1_env%v_rspace)
                   call cp_logger_generate_filename(logger=logger,res=filename,&
                        root=trim(base_filename)//"%v_rspace-"//&
                        trim(adjustl(cp_to_string(ispin))), postfix=".cube",&
                        local=my_local)
                   CALL open_file(file_name=trim(filename),file_status="UNKNOWN",&
                        file_action="WRITE",unit_number=file_unit)
                   CALL rs_pw_to_cube ( kpp1_env%v_rspace(ispin)%pw,&
                        iunit=file_unit,&
                        ionode=logger%para_env%mepos==logger%para_env%source,&
                        title="v_rspace")
                   CALL close_file(file_unit)
                end do
                WRITE (unit=my_unit_nr,fmt="(a,a,a)") " v_rspace=*written to ",&
                     trim(filename),"*"
             else
                WRITE (unit=my_unit_nr,fmt="(a)") " v_rspace=*associated*"
             end if
          else
             WRITE (unit=my_unit_nr,fmt="(a)") " v_rspace=*null*"
          end if

          if (associated(kpp1_env%v_ao)) then
             if (my_long_description) then
                do ispin=1,size(kpp1_env%v_ao)
                   call cp_sm_output(logger=logger,&
                        outputName=trim(base_filename)//"%v_ao-"//&
                        trim(adjustl(cp_to_string(ispin))),&
                        fromWhere=routineP,iter=kpp1_env%print_count,&
                        matrix=kpp1_env%v_ao(ispin)%matrix,&
                        para_env=logger%para_env,&
                        error=error)
                end do
                WRITE (unit=my_unit_nr,fmt="(a,a,a)") " v_ao=*written to ",&
                     filename,"*"
             else
                write (unit=my_unit_nr,fmt="(a)") " v_ao=*associated*"
             end if
          else
             write (unit=my_unit_nr,fmt="(a)") " v_ao=*null*"
          end if

          if (associated(kpp1_env%drho_r)) then
             if (my_long_description) then
                do ispin=1,size(kpp1_env%drho_r,2)
                   do idir=1,size(kpp1_env%drho_r,1)
                      call cp_logger_generate_filename(logger=logger,&
                           res=filename,&
                           root=trim(base_filename)//"%drho_r-"//&
                           trim(adjustl(cp_to_string(ispin)))//&
                           trim(adjustl(cp_to_string(idir))),&
                           postfix=".cube",&
                           local=my_local)
                      CALL open_file(file_name=trim(filename),&
                           file_status="UNKNOWN",&
                           file_action="WRITE",unit_number=file_unit)
                      CALL rs_pw_to_cube ( kpp1_env%v_rspace(ispin)%pw,&
                           iunit=file_unit,&
                           ionode=logger%para_env%mepos==logger%para_env%source,&
                           title="drho_r")
                      CALL close_file(file_unit)
                   end do
                end do
                WRITE (unit=my_unit_nr,fmt="(a,a,a)") " v_ao=*written to ",&
                     filename,"*"
             else
                write (unit=my_unit_nr,fmt="(a)") " drho_r=*associated*"
             end if
          else
             write (unit=my_unit_nr,fmt="(a)") " drho_r=*null*"
          end if

          if (associated(kpp1_env%derivative_set)) then
             if (my_long_description) then
                call cp_2d_r_output(logger=logger,&
                     outputName=base_filename//"%derivative_set%data", &
                     fromWhere=routineP, iter=kpp1_env%print_count,&
                     array=kpp1_env%derivative_set%data,&
                     local=my_local, error=error)
             else
                write (unit=my_unit_nr,fmt="(a)") " deriv_xc=*associated*,"
             end if
          else
             write (unit=my_unit_nr,fmt="(a)") " deriv_xc=*null*,"
          end if

          if (associated(kpp1_env%spin_pot)) then
             if (my_long_description) then
                write (unit=my_unit_nr,fmt="(a)") " spin_pot=("
                call cp_2d_i_write(array=kpp1_env%spin_pot,&
                     unit_nr=my_unit_nr, error=error)
                write (unit=my_unit_nr,fmt="(a)") " ),"
             else
                write (unit=my_unit_nr,fmt="(a)") " spin_pot=*associated*,"
             end if
          else
             write (unit=my_unit_nr,fmt="(a)") " spin_pot=*null*,"
          end if

          if (associated(kpp1_env%grad_pot)) then
             if (my_long_description) then
                write (unit=my_unit_nr,fmt="(a)") " grad_pot=("
                call cp_2d_logical_write(array=kpp1_env%grad_pot,&
                     unit_nr=my_unit_nr, error=error)
                write (unit=my_unit_nr,fmt="(a)") " ),"
             else
                write (unit=my_unit_nr,fmt="(a)") " grad_pot=*associated*,"
             end if
          else
             write (unit=my_unit_nr,fmt="(a)") " grad_pot=*null*"
          end if

          if (associated(kpp1_env%ndiag_term)) then
             write (unit=my_unit_nr,fmt="(a)") " grad_pot=("
             if (my_long_description) then
                write (unit=my_unit_nr,fmt="(a)") " ndiag_term=("
                call cp_1d_logical_write(array=kpp1_env%ndiag_term,&
                     unit_nr=my_unit_nr, error=error)
                write (unit=my_unit_nr,fmt="(a)") " ),"
             else
                write (unit=my_unit_nr,fmt="(a)") " ndiag_term=*associated*"
             end if
          else
             write (unit=my_unit_nr,fmt="(a)") " ndiag_term=*null*"
          end if

          if (my_long_description) then
             call close_file(my_unit_nr)
          end if
       end IF
    ELSE
       IF (should_w) THEN
          WRITE (unit=my_unit_nr,fmt="(a)") " <kpp1_env>:*null*"
       END IF
    END IF
  END SUBROUTINE kpp1_write
!***************************************************************************

END MODULE qs_kpp1_env_types
