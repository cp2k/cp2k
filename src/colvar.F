#include "cp_prep_globals.h"
MODULE colvar
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/colvar [1.0] *
!!
!!   NAME
!!     colvar
!!
!!   FUNCTION
!!     collective variables for metadynamics
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Alessandro Laio,Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     04.2004 created
!!
!!   SOURCE
!****************************************************************************
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE metadynamics_types,              ONLY: colvar_type,&
                                             dist_colvar_id,&
                                             no_colvar_id,&
                                             coord_colvar_id
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_parser,                       ONLY: read_object,&
                                             search_object,&
                                             stop_parser,&
                                             test_object
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE physcon, ONLY: massunit
  USE convert_units, only: convert_to_cp2k_units

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='colvar'

  PUBLIC :: colvar_read
  PUBLIC :: colvar_eval_f

!***
!****************************************************************************
CONTAINS

!!****f* metadynamics/colvar_read *
!!
!!   NAME
!!     colvar_read
!!
!!   FUNCTION
!!     reads a colvar with the parser. the input should have this form:
!!     &colvar dist
!!       lambda 1.0
!!       mass   1.0
!!       atoms 1 3
!!     &end
!!
!!   NOTES
!!     the parser should already be initialized
!!
!!   INPUTS
!!     - colvar: the place where to store what will be read
!!     - globenv: inputfile,...
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!   MODIFICATION HISTORY
!!     04.2004 created
!!
!!*** **********************************************************************
  SUBROUTINE colvar_read(colvar, error)
    TYPE(colvar_type), INTENT(inout)         :: colvar
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'colvar_read', &
         routinep = modulen//':'//routinen

    CHARACTER(len=3)                         :: test_result
    CHARACTER(len=default_string_length)     :: colvar_type_name, keyword
    LOGICAL                                  :: failure, RECOGNIZED_KEYWORD, &
         section_found
    INTEGER :: n_opt,n_opt_cv,n_opt_var,ntemp,i,stat

    failure=.FALSE.
    n_opt=0
    n_opt_cv=0
    n_opt_var=0

    CALL search_object("&COLVAR",ignore_case=.TRUE.,&
         found=section_found)
    CALL cp_assert(section_found,cp_failure_level,cp_assertion_failed,routinep,&
         "an expected colvar was not found",error,failure)
    IF (.NOT. failure) THEN
       CALL read_object(colvar_type_name,lower_to_upper=.TRUE.,&
            newline=.FALSE., string_length=default_string_length)
       SELECT CASE(colvar_type_name)
       CASE("DIST","DISTANCE")
          colvar%type_id=dist_colvar_id
          colvar%dist_param%i_subsys=1
          colvar%dist_param%j_subsys=1
          n_opt_var=1
       CASE("COORD","COORDINATION")
          colvar%type_id=coord_colvar_id
          colvar%coord_param%i_subsys=1
          colvar%coord_param%j_subsys=1
          n_opt_var=3
       CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routinep,&
               "Unknown collective variable type '"//TRIM(colvar_type_name)//"'",&
               error,failure)
          colvar%type_id=no_colvar_id
       END SELECT
    END IF
    IF (.NOT.failure) THEN
       DO
          test_result = test_object(newline=.TRUE.)
          IF (test_result == "EOS") THEN
             EXIT
          ELSE IF (test_result == "EOF") THEN
             CALL stop_parser(routinen,modulen,__LINE__,"EOF")
          END IF

          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE(keyword)
          CASE ("LAMBDA")
             CALL read_object(colvar%lambda)
             n_opt=n_opt+1
          CASE ("MASS","MASS_S")
             CALL read_object(colvar%mass)
!             CALL convert_to_cp2k_units('u',mass=colvar%mass)
             colvar%mass=colvar%mass*massunit
             n_opt=n_opt+1
          CASE ("SCALE")
             CALL read_object(colvar%delta_s)
             n_opt=n_opt+1
          CASE default
             ! colvar type specific keyword
             recognized_keyword=.TRUE.

             SELECT CASE(colvar%type_id)
             CASE(dist_colvar_id)
                SELECT CASE(keyword)
                CASE ("ATOMS")
                   CALL read_object(colvar%dist_param%i_at)
                   CALL read_object(colvar%dist_param%j_at)
                   n_opt_cv=n_opt_cv+1
                CASE ("SUBSYS")
                   CALL read_object(colvar%dist_param%i_subsys)
                   CALL read_object(colvar%dist_param%j_subsys)
                CASE default
                   recognized_keyword=.FALSE.
                END SELECT
             CASE(coord_colvar_id)
                SELECT CASE(keyword)
                CASE ("ATOMS_FROM")
                   CALL read_object(ntemp)
                   colvar%coord_param%n_atoms_from=ntemp
                   ALLOCATE(colvar%coord_param%i_at_from(ntemp),stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
                   do i=1,ntemp
                      CALL read_object(colvar%coord_param%i_at_from(i))
                   enddo
                   n_opt_cv=n_opt_cv+1
                CASE ("SUBSYS_FROM")
                   CALL read_object(colvar%coord_param%i_subsys)
                CASE ("ATOMS_TO")
                   CALL read_object(ntemp)
                   colvar%coord_param%n_atoms_to=ntemp
                   ALLOCATE(colvar%coord_param%i_at_to(ntemp),stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
                   do i=1,ntemp
                      CALL read_object(colvar%coord_param%i_at_to(i))
                   enddo
                   n_opt_cv=n_opt_cv+1
                CASE ("SUBSYS_TO")
                   CALL read_object(colvar%coord_param%j_subsys)
                CASE ("R0","R_0")
                   CALL read_object(colvar%coord_param%r_0)
                   n_opt_cv=n_opt_cv+1
                CASE ("NN","EXPON_NUMERATOR")
                   CALL read_object(colvar%coord_param%nncrd)
                CASE ("ND","EXPON_DENOMINATOR")
                   CALL read_object(colvar%coord_param%ndcrd)
                CASE default
                   recognized_keyword=.FALSE.
                END SELECT
             CASE default
                CPAssert(.FALSE.,cp_failure_level,routinep,error,failure)
             END SELECT

             CALL cp_assert(recognized_keyword,cp_failure_level,&
                  cp_assertion_failed,routinep,"unknown keyword '"//TRIM(keyword)//&
                  "' in section colvar", error,failure)

          END SELECT
       END DO
    END IF
    CALL cp_assert(n_opt==3,cp_failure_level,cp_assertion_failed,routinep,&
         "mass, k or scaling factor missing from cv"//cp_to_string(colvar%index),error,failure)
    CALL cp_assert(n_opt_cv==n_opt_var,cp_failure_level,cp_assertion_failed,routinep,&
         "missing non optional parameter from CV "//cp_to_string(colvar%index), &
      error,failure)
  END SUBROUTINE colvar_read
!***************************************************************************

!!****f* colvar/colvar_eval_f [1.0] *
!!
!!   NAME
!!     colvar_eval_f
!!
!!   FUNCTION
!!     evaluates the force given and due to the given colvar
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - colvar: the collective variable to evaluate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Alessandro Laio and fawzi
!!
!!*** **********************************************************************
SUBROUTINE colvar_eval_f(colvar,force_env,n_step,error)
    TYPE(colvar_type), INTENT(inout)         :: colvar
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, INTENT(in)                      :: n_step
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_eval_f', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys

  failure=.FALSE.
  CPPrecondition(n_step>0,cp_failure_level,routineP,error,failure)
  NULLIFY(subsys,cell)
  CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)
  IF (.NOT. failure) THEN
     SELECT CASE(colvar%type_id)
     CASE (dist_colvar_id)
        CALL dist_colvar(colvar,subsys,cell,&
             n_step=n_step,error=error)
     CASE (coord_colvar_id)
        CALL coord_colvar(colvar,subsys,cell,&
             n_step=n_step,error=error)
     CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     END SELECT
  END IF
END SUBROUTINE colvar_eval_f
!***************************************************************************

!!****f* colvar/dist_colvar [1.0] *
!!
!!   NAME
!!     dist_colvar
!!
!!   FUNCTION
!!     evaluates the force due (and on) the distance collective variable
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - subsys: the atom positions,...
!!     - cell: cell information (to enforce pbc)
!!     - n_step: the number of steps done so far (at the first one
!!       some initialization is performed)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Alessandro Laio, Fawzi Mohamed
!!
!!*** **********************************************************************
  SUBROUTINE dist_colvar(colvar,subsys,cell,n_step,error)
    TYPE(colvar_type)                        :: colvar
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(cell_type), POINTER                 :: cell
    INTEGER                                  :: n_step
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dist_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dfunc, r12, ss(3), xij(3), &
                                                xpi(3), xpj(3)
    TYPE(particle_list_type), POINTER        :: particles_i, particles_j

    failure=.FALSE.
    NULLIFY(particles_i,particles_j)
    CPPrecondition(colvar%type_id==dist_colvar_id,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys(colvar%dist_param%i_subsys)%subsys,&
         particles=particles_i,error=error)
    CALL cp_subsys_get(subsys(colvar%dist_param%j_subsys)%subsys,&
         particles=particles_j, error=error)

    i=colvar%dist_param%i_at
    j=colvar%dist_param%j_at
    xpi=particles_i%els(i)%r(1:3)
    xpj=particles_j%els(j)%r(1:3)
    ss=MATMUL(cell%h_inv,xpi-xpj)
    ss=ss-NINT(ss)
    xij=MATMUL(cell%hmat,ss)
    r12=SQRT(xij(1)**2+xij(2)**2+xij(3)**2)
    colvar%ss=r12
    IF(n_step==1)THEN
     colvar%ss0=colvar%ss
    ENDIF

    colvar%cv_energy=0.5_dp*colvar%lambda*(colvar%ss-colvar%ss0)**2
    colvar%ff_s=colvar%lambda* (colvar%ss-colvar%ss0)
    dfunc=colvar%ff_s/ r12
    particles_i%els(i)%f(1:3)=particles_i%els(i)%f(1:3)-dfunc*xij
    particles_j%els(j)%f(1:3)=particles_j%els(j)%f(1:3)+dfunc*xij

  END SUBROUTINE dist_colvar
!***************************************************************************


!!****f* colvar/coord_colvar [1.0] *
!!
!!   NAME
!!     coord_colvar
!!
!!   FUNCTION
!!     evaluates the force due (and on) the coordination collective variable
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - subsys: the atom positions,...
!!     - cell: cell information (to enforce pbc)
!!     - n_step: the number of steps done so far (at the first one
!!       some initialization is performed)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Alessandro Laio, Fawzi Mohamed
!!
!!*** **********************************************************************
  SUBROUTINE coord_colvar(colvar,subsys,cell,n_step,error)
    INTEGER :: n_step
    TYPE(colvar_type) :: colvar
    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys
    TYPE(cell_type), POINTER :: cell
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='coord_colvar',&
         routineP=moduleN//':'//routineN
    TYPE(particle_list_type), POINTER :: particles_i,particles_j
    
!  locals
    integer :: i,ii,jj,n_atoms_to,n_atoms_from
    REAL(dp) :: xij(3),ss(3),r12,dfunc
    REAL(dp) , POINTER, DIMENSION (:,:) :: xpi,xpj,fpi,fpj
    real(dp) :: r_0,rdist,num,invden,func,ncoord
    integer :: nncrd,ndcrd,stat
    
    n_atoms_to=colvar%coord_param%n_atoms_to
    n_atoms_from=colvar%coord_param%n_atoms_from
    nncrd=colvar%coord_param%nncrd
    ndcrd=colvar%coord_param%ndcrd
    r_0=colvar%coord_param%r_0

    allocate(xpi(3,n_atoms_from),xpj(3,n_atoms_to), &
             fpi(3,n_atoms_from),fpj(3,n_atoms_to),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    failure=.FALSE.
    NULLIFY(particles_i,particles_j)
    CPPrecondition(colvar%type_id==coord_colvar_id,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys(colvar%coord_param%i_subsys)%subsys,&
         particles=particles_i,error=error)
    CALL cp_subsys_get(subsys(colvar%coord_param%j_subsys)%subsys,&
         particles=particles_j, error=error)

    do ii=1,n_atoms_from
      i=colvar%coord_param%i_at_from(ii)
      xpi(:,ii)=particles_i%els(i)%r(1:3)
    enddo
    do ii=1,n_atoms_to
      i=colvar%coord_param%i_at_to(ii)
      xpj(:,ii)=particles_i%els(i)%r(1:3)
    enddo
    ncoord=0.0_dp
    fpi=0.0_dp
    fpj=0.0_dp
    do ii=1,n_atoms_from
      do jj=1,n_atoms_to
        ss=matmul(cell%h_inv,xpi(:,ii)-xpj(:,jj))
        ss=ss-NINT(ss)
        xij=matmul(cell%hmat,ss)
        r12=sqrt(xij(1)**2+xij(2)**2+xij(3)**2)

        rdist = r12/r_0
        num=(1.0_dp-rdist**nncrd)
        invden=1.0_dp/(1.0_dp-rdist**ndcrd)
        func=num*invden
        dfunc= (- nncrd*rdist **(nncrd-1) *invden  &
          + num*(invden)**2 * ndcrd *rdist **(ndcrd-1))/(r12*r_0)
                                                                               
        ncoord=ncoord+func
        fpi(:,ii)=fpi(:,ii)+dfunc*xij
        fpj(:,jj)=fpj(:,jj)-dfunc*xij
      enddo
    enddo
    colvar%ss=ncoord
    if(n_step==1)then
     colvar%ss0=colvar%ss
    endif

    colvar%cv_energy=0.5_dp*colvar%lambda*(colvar%ss-colvar%ss0)**2
    colvar%ff_s=colvar%lambda* (colvar%ss-colvar%ss0)

    do ii=1,n_atoms_from
      i=colvar%coord_param%i_at_from(ii)
      particles_i%els(i)%f=particles_i%els(i)%f+fpi(:,ii)     
    enddo
    do ii=1,n_atoms_to
      i=colvar%coord_param%i_at_to(ii)
      particles_i%els(i)%f=particles_i%els(i)%f+fpj(:,ii)     
    enddo
    deallocate(xpi,xpj,fpi,fpj,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
  END SUBROUTINE coord_colvar
!***************************************************************************

END MODULE colvar
