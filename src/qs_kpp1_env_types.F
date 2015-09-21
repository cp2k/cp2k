!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief basis types for the calculation of the perturbation of density theory.
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_kpp1_env_types
  USE cp_array_utils_i,                ONLY: cp_2d_i_write
  USE cp_array_utils_logical,          ONLY: cp_1d_logical_write,&
                                             cp_2d_logical_write
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_p_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
  USE kinds,                           ONLY: default_path_length
  USE particle_list_types,             ONLY: particle_list_type
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_release
  USE string_utilities,                ONLY: compress
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_release
  USE xc_rho_set_types,                ONLY: xc_rho_set_release,&
                                             xc_rho_set_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: qs_kpp1_env_type
  PUBLIC :: kpp1_release

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kpp1_env_types'

! *****************************************************************************
!> \brief environement that keeps the informations and temporary
!>      val to build the kpp1 kernel matrix
!> \param ref_count reference count (how may objects are sharing this one)
!> \param v_rspace potential in r space. This is used to do an update only
!>        of what has changed. Useful if P1 converges to some density
!>        (you spare in the grid-ao conversion).
!>        With LSD contains the various orbitals.
!> \param v_ao the potential in the ao basis (used togheter with v_rspace
!>        to update only what changed
!> \param id_nr identification number, unique for each kpp1 env
!> \param print_count counter to create unique filename
!> \param iter number of iterations
!> \param drho_r (idir,ispin): the derivative of rho wrt. x,y,z in the real space
!> \param deriv_xc (ii,ipot): the second derivative of the xc potential at psi0
!>        (qs_env%c), if grad pot is true it should already be divised
!>        by the gradient
!> \param spin_pot (1:2,ipot): information about wrt. to which spins the
!>         corresponding component of deriv_xc was derived (see
!>         xc_create_2nd_deriv_info)
!> \param grad_pot (1:2,ipot): if the derivative spin_pot was wrt. to
!>        the gradient (see xc_create_2nd_deriv_info)
!> \param ndiag_term (ipot): it the term is an off diagonal term (see
!>        xc_create_2nd_deriv_info)
! *****************************************************************************
  TYPE qs_kpp1_env_type
     INTEGER :: ref_count, id_nr, print_count, iter
     TYPE(pw_p_type), DIMENSION(:), POINTER :: v_rspace
     TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: v_ao
     TYPE(pw_p_type), DIMENSION(:,:), POINTER :: drho_r
     TYPE(xc_derivative_set_type), POINTER :: deriv_set
     TYPE(xc_rho_set_type), POINTER        :: rho_set
     INTEGER, DIMENSION(:,:), POINTER :: spin_pot
     LOGICAL, DIMENSION(:,:), POINTER :: grad_pot
     LOGICAL, DIMENSION(:), POINTER :: ndiag_term
  END TYPE qs_kpp1_env_type

! *****************************************************************************
!> \brief just to build array of pointers
!> \param kpp1_env: the pointer to the kpp1_env
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_kpp1_env_p_type
     TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
  END TYPE qs_kpp1_env_p_type

CONTAINS

! *****************************************************************************
!> \brief releases a kpp1_env (see doc/ReferenceCounting.html)
!> \param kpp1_env the environement to release
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE kpp1_release(kpp1_env)
    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(kpp1_env)) THEN
       CPASSERT(kpp1_env%ref_count>0)
       kpp1_env%ref_count=kpp1_env%ref_count-1
       IF (kpp1_env%ref_count<1) THEN
          IF (ASSOCIATED(kpp1_env%v_rspace)) THEN
             DO ispin=1,SIZE(kpp1_env%v_rspace)
                CALL pw_release(kpp1_env%v_rspace(ispin)%pw)
             END DO
             DEALLOCATE(kpp1_env%v_rspace)
          END IF
          IF (ASSOCIATED(kpp1_env%v_ao)) THEN
             DO ispin=1,SIZE(kpp1_env%v_ao)
                IF (ASSOCIATED(kpp1_env%v_ao(ispin)%matrix)) THEN
                   CALL cp_dbcsr_deallocate_matrix(kpp1_env%v_ao(ispin)%matrix)
                END IF
             END DO
             DEALLOCATE(kpp1_env%v_ao)
          END IF
          IF (ASSOCIATED(kpp1_env%drho_r)) THEN
             DEALLOCATE(kpp1_env%drho_r)
          END IF
          IF (ASSOCIATED(kpp1_env%deriv_set)) THEN
             CALL xc_dset_release(kpp1_env%deriv_set)
             NULLIFY(kpp1_env%deriv_set)
          END IF
          IF (ASSOCIATED(kpp1_env%rho_set)) THEN
             CALL xc_rho_set_release(kpp1_env%rho_set)
             NULLIFY(kpp1_env%rho_set)
          END IF
          IF (ASSOCIATED(kpp1_env%spin_pot)) THEN
             DEALLOCATE(kpp1_env%spin_pot)
          END IF
          IF (ASSOCIATED(kpp1_env%grad_pot)) THEN
             DEALLOCATE(kpp1_env%grad_pot)
          END IF
          IF (ASSOCIATED(kpp1_env%ndiag_term)) THEN
             DEALLOCATE(kpp1_env%ndiag_term)
          END IF
          DEALLOCATE(kpp1_env)
       END IF
    END IF
    NULLIFY(kpp1_env)
  END SUBROUTINE kpp1_release

! *****************************************************************************
!> \brief retains a kpp1_env (see doc/ReferenceCounting.html)
!> \param kpp1_env the environement to retain
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE kpp1_retain(kpp1_env)
    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(ASSOCIATED(kpp1_env))
    CPASSERT(kpp1_env%ref_count>0)
    kpp1_env%ref_count=kpp1_env%ref_count+1
  END SUBROUTINE kpp1_retain

! *****************************************************************************
!> \brief writes the content of the kpp1_env to the given output unit
!> \param kpp1_env the kpp1_env to output
!> \param particles ...
!> \param unit_nr the unit where to output
!> \param long_description if a long description should be written
!>        (defaults to false)
!> \param local if the unit is a local unit or a global unit
!>        (defaults to false, i.e. global)
!> \par History
!>      12.2002 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE kpp1_write(kpp1_env, particles, unit_nr, long_description, local)
    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env
    TYPE(particle_list_type), POINTER        :: particles
    INTEGER, INTENT(in)                      :: unit_nr
    LOGICAL, INTENT(in), OPTIONAL            :: long_description, local

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_write', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: base_filename, filename
    INTEGER                                  :: file_unit, idir, ispin, &
                                                my_unit_nr
    LOGICAL                                  :: failure, my_local, &
                                                my_long_description, should_w
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env

    failure=.FALSE.; my_local=.FALSE.; my_long_description=.FALSE.

    IF (PRESENT(local)) my_local=local
    IF (PRESENT(long_description)) my_long_description=long_description
    logger=>cp_get_default_logger()
    para_env=>logger%para_env
    should_w=my_local .OR. para_env%mepos==para_env%source
    my_unit_nr=unit_nr

    IF (ASSOCIATED(kpp1_env)) THEN
       kpp1_env%print_count=kpp1_env%print_count+1

       IF (should_w) THEN
          base_filename="kpp1-"//TRIM(ADJUSTL(cp_to_string(kpp1_env%id_nr)))&
               //"-"//TRIM(ADJUSTL(cp_to_string(kpp1_env%print_count)))
          CALL compress(base_filename,full=.TRUE.)

          IF (my_long_description) THEN
             WRITE (unit=unit_nr,fmt="(a,a,a)") " wrote kpp1 into to files ",&
                  TRIM(base_filename),"*"
             CALL cp_logger_generate_filename(logger=logger,res=filename,&
                  root=TRIM(base_filename), postfix=".desc",&
                  local=my_local)
             CALL open_file(file_name=TRIM(filename),file_status="UNKNOWN",&
                  file_action="WRITE",unit_number=my_unit_nr)
          END IF

          WRITE (unit=my_unit_nr,&
               fmt="(' <kpp1_env>:{ id_nr=',i10,' ref_count=',i10,',')")&
               kpp1_env%id_nr, kpp1_env%ref_count

          IF (ASSOCIATED(kpp1_env%v_rspace)) THEN
             IF (my_long_description) THEN
                DO ispin=1,SIZE(kpp1_env%v_rspace)
                   CALL cp_logger_generate_filename(logger=logger,res=filename,&
                        root=TRIM(base_filename)//"%v_rspace-"//&
                        TRIM(ADJUSTL(cp_to_string(ispin))), postfix=".cube",&
                        local=my_local)
                   CALL open_file(file_name=TRIM(filename),file_status="UNKNOWN",&
                        file_action="WRITE",unit_number=file_unit)
                   CALL cp_pw_to_cube ( kpp1_env%v_rspace(ispin)%pw,&
                        unit_nr=file_unit,&
                        particles=particles,&
                        title="v_rspace")
                   CALL close_file(file_unit)
                END DO
                WRITE (unit=my_unit_nr,fmt="(a,a,a)") " v_rspace=*written to ",&
                     TRIM(filename),"*"
             ELSE
                WRITE (unit=my_unit_nr,fmt="(a)") " v_rspace=*associated*"
             END IF
          ELSE
             WRITE (unit=my_unit_nr,fmt="(a)") " v_rspace=*null*"
          END IF

          IF (ASSOCIATED(kpp1_env%v_ao)) THEN
!FM             IF (my_long_description) THEN
!FM                DO ispin=1,SIZE(kpp1_env%v_ao)
!FM                   CALL cp_sm_write(logger=logger,&
!FM                        outputName=TRIM(base_filename)//"%v_ao-"//&
!FM                        TRIM(ADJUSTL(cp_to_string(ispin))),&
!FM                        fromWhere=routineP,iter=kpp1_env%print_count,&
!FM                        matrix=kpp1_env%v_ao(ispin)%matrix,&
!FM                        para_env=logger%para_env,&
!FM                        
!FM                END DO
!FM                WRITE (unit=my_unit_nr,fmt="(a,a,a)") " v_ao=*written to ",&
!FM                     filename,"*"
!FM             ELSE
                WRITE (unit=my_unit_nr,fmt="(a)") " v_ao=*associated*"
!FM             END IF
          ELSE
             WRITE (unit=my_unit_nr,fmt="(a)") " v_ao=*null*"
          END IF

          IF (ASSOCIATED(kpp1_env%drho_r)) THEN
             IF (my_long_description) THEN
                DO ispin=1,SIZE(kpp1_env%drho_r,2)
                   DO idir=1,SIZE(kpp1_env%drho_r,1)
                      CALL cp_logger_generate_filename(logger=logger,&
                           res=filename,&
                           root=TRIM(base_filename)//"%drho_r-"//&
                           TRIM(ADJUSTL(cp_to_string(ispin)))//&
                           TRIM(ADJUSTL(cp_to_string(idir))),&
                           postfix=".cube",&
                           local=my_local)
                      CALL open_file(file_name=TRIM(filename),&
                           file_status="UNKNOWN",&
                           file_action="WRITE",unit_number=file_unit)
                      CALL cp_pw_to_cube ( kpp1_env%v_rspace(ispin)%pw,&
                           unit_nr=file_unit,&
                           particles=particles,&
                           title="drho_r")
                      CALL close_file(file_unit)
                   END DO
                END DO
                WRITE (unit=my_unit_nr,fmt="(a,a,a)") " v_ao=*written to ",&
                     filename,"*"
             ELSE
                WRITE (unit=my_unit_nr,fmt="(a)") " drho_r=*associated*"
             END IF
          ELSE
             WRITE (unit=my_unit_nr,fmt="(a)") " drho_r=*null*"
          END IF

          IF (ASSOCIATED(kpp1_env%deriv_set)) THEN
!FM             IF (my_long_description) THEN
!FM                CALL cp_2d_r_output(logger=logger,&
!FM                     outputName=base_filename//"%deriv_set%data", &
!FM                     fromWhere=routineP, iter=kpp1_env%print_count,&
!FM                     array=kpp1_env%deriv_set%data,&
!FM                     local=my_local)
!FM             ELSE
                WRITE (unit=my_unit_nr,fmt="(a)") " deriv_xc=*associated*,"
!FM             END IF
          ELSE
             WRITE (unit=my_unit_nr,fmt="(a)") " deriv_xc=*null*,"
          END IF

          IF (ASSOCIATED(kpp1_env%spin_pot)) THEN
             IF (my_long_description) THEN
                WRITE (unit=my_unit_nr,fmt="(a)") " spin_pot=("
                CALL cp_2d_i_write(array=kpp1_env%spin_pot,&
                     unit_nr=my_unit_nr)
                WRITE (unit=my_unit_nr,fmt="(a)") " ),"
             ELSE
                WRITE (unit=my_unit_nr,fmt="(a)") " spin_pot=*associated*,"
             END IF
          ELSE
             WRITE (unit=my_unit_nr,fmt="(a)") " spin_pot=*null*,"
          END IF

          IF (ASSOCIATED(kpp1_env%grad_pot)) THEN
             IF (my_long_description) THEN
                WRITE (unit=my_unit_nr,fmt="(a)") " grad_pot=("
                CALL cp_2d_logical_write(array=kpp1_env%grad_pot,&
                     unit_nr=my_unit_nr)
                WRITE (unit=my_unit_nr,fmt="(a)") " ),"
             ELSE
                WRITE (unit=my_unit_nr,fmt="(a)") " grad_pot=*associated*,"
             END IF
          ELSE
             WRITE (unit=my_unit_nr,fmt="(a)") " grad_pot=*null*"
          END IF

          IF (ASSOCIATED(kpp1_env%ndiag_term)) THEN
             WRITE (unit=my_unit_nr,fmt="(a)") " grad_pot=("
             IF (my_long_description) THEN
                WRITE (unit=my_unit_nr,fmt="(a)") " ndiag_term=("
                CALL cp_1d_logical_write(array=kpp1_env%ndiag_term,&
                     unit_nr=my_unit_nr)
                WRITE (unit=my_unit_nr,fmt="(a)") " ),"
             ELSE
                WRITE (unit=my_unit_nr,fmt="(a)") " ndiag_term=*associated*"
             END IF
          ELSE
             WRITE (unit=my_unit_nr,fmt="(a)") " ndiag_term=*null*"
          END IF

          IF (my_long_description) THEN
             CALL close_file(my_unit_nr)
          END IF
       END IF
    ELSE
       IF (should_w) THEN
          WRITE (unit=my_unit_nr,fmt="(a)") " <kpp1_env>:*null*"
       END IF
    END IF
  END SUBROUTINE kpp1_write

END MODULE qs_kpp1_env_types
