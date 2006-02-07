!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/replica_methods *
!!
!!   NAME
!!     replica_methods
!!
!!   FUNCTION
!!     methods to setup replicas of the same system differing only by atom
!!     positions and velocities (as used in path integral or nudged elastic
!!     band for example)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     09.2005 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE replica_methods
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE replica_types,                   ONLY: replica_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'replica_methods'

  INTEGER, SAVE, PRIVATE :: last_rep_env_id=0

!!***
!****************************************************************************
CONTAINS

!!****f* replica_methods/rep_and_force_env_create *
!!
!!   NAME
!!     rep_env_create
!!
!!   FUNCTION
!!     creates a replica environment together with its force environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rep_env: the replica environment that will be created
!!     - para_env: the parallel enviroment that will contain the replicas
!!     - input: the input used to initialize the force environment
!!     - nrep: the number of replicas to calculate
!!     - prep: the number of processors for each replica
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE rep_env_create(rep_env, para_env, input, nrep, prep, error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    INTEGER                                  :: nrep, prep
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: comm_cart, comm_f, new_env_id
    INTEGER, DIMENSION(2)                    :: dims, pos
    LOGICAL                                  :: failure
    LOGICAL, DIMENSION(2)                    :: rdim
    TYPE(cp_para_cart_type), POINTER         :: cart
    TYPE(cp_para_env), POINTER               :: para_env_f

  failure=.FALSE.
  
  CPPrecondition(.NOT.ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(.NOT.ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  NULLIFY(cart,para_env_f)
  IF (.NOT. failure) THEN
     new_env_id=-1
     my_prep=MIN(prep,para_env%num_pe)
     dims(1)=my_prep
     dims(2)=para_env%num_pe/my_prep
     rdim(1)=.TRUE.
     rdim(2)=.FALSE.
     unit_nr=cp_logger_get_default_unit_nr(logger,local=.FALSE.)
     IF (unit_nr>0 .AND. dims(1)*dims(2)/=para_env%num_pe) THEN
        WRITE(unit_nr,"(/,' WARNING',/,' number of processors not divisible in replica creation,')")
        WRITE(unit_nr,"(' wasting ',i4,' processors')") para_env%num_pe-dims(1)*dims(2)
     END IF
     CALL mp_cart_create ( comm_old, ndims=2, dims=dims, pos=pos, comm_cart=comm_cart)
     IF (comm_cart/=MPI_COMM_NULL) THEN
        CALL cp_cart_create(cart,comm_cart,ndims=2,owns_group=.TRUE.,error=error)
        CALL mp_cart_sub( comm=comm_cart, rdim=rdim , sub_comm=comm_f)
        CALL cp_para_env_create(para_env_f,comm_f,owns_group=.TRUE.,error=error)
        input_file_path=
        CALL create_force_env(new_env_id=new_env_id,input_file_path,&
             output_file_path,mpi_comm=comm_f,ierr)
        
        ALLOCATE(rep_env,stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
        
     END IF
  END IF
  IF (ASSOCIATED(rep_env)) THEN
     last_rep_env_id=last_rep_env_id+1
     rep_env%id_nr=last_rep_env_id
     rep_env%ref_count=1
     rep_env%f_env_id=new_env_id
     CALL get_natom(new_env_id,natom,error=error)
     ALLOCATE(rep_env%replica_owner(nrep),stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

     DO i=1,nrep
        rep_env%replica_owner(i)=MODULO(i,cart
     ALLOCATE(rep_env%x(natom,nrep),rep_env%y(natom,nrep),&
          rep_env%z(natom,nrep), stat=stat)
     
     REAL(kind=dp), DIMENSION(:,:,:), POINTER :: x,v,f
     TYPE(cp_para_cart_type), POINTER :: para_cart
     TYPE(cp_para_env_type), POINTER :: para_env
END SUBROUTINE rep_env_create
!***************************************************************************

END MODULE replica_methods
