MODULE colvar_types
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_para_env,                     ONLY: cp_para_env_release
  USE kinds,                           ONLY: default_path_length,&
                                             dp
#include "cp_common_uses.h"
  IMPLICIT NONE

  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleN='colvar_types'

  PUBLIC :: colvar_type, coord_colvar_type,dist_colvar_type,colvar_init,&
        colvar_dealloc_ref,torsion_colvar_type
  INTEGER, PARAMETER, PUBLIC :: no_colvar_id=-2,dist_colvar_id=1,&
       coord_colvar_id=2,torsion_colvar_id=3

!!****s* colvar/dist_colvar_param [1.0] *
!!
!!   NAME
!!     dist_colvar_param
!!
!!   FUNCTION
!!     parameters for the distance collective variable
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - i_at,j_at: indexes of the two atoms between which you calculate
!!       the distance
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!   SOURCE
  !***************************************************************************
  TYPE dist_colvar_type
     INTEGER :: i_at,j_at
  END TYPE dist_colvar_type
  TYPE coord_colvar_type
     INTEGER :: n_atoms_to,n_atoms_from,nncrd,ndcrd
     INTEGER, POINTER, DIMENSION (:) :: i_at_from,i_at_to
     REAL(dp) :: r_0
  END TYPE coord_colvar_type
  TYPE torsion_colvar_type
     INTEGER, DIMENSION (4) :: i_at_tors
  END TYPE torsion_colvar_type
!!***
  !****************************************************************************

!!****s* colvar_types/dist_colvar_param [1.0] *
!!
!!   NAME
!!     dist_colvar_param
!!
!!   FUNCTION
!!     parameters for a collective variable
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - index: index of this the colvar
!!     - type_id: the type of this colvar
!!     - mass_s mass of the collective variable
!!     - lambda: lambda of the gaussian in the direction of this colvar
!!     - dist_param: parameters for the distance collective variable
!!       (meaningful only if type_id==dist_colvar_id)
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!   SOURCE
  !***************************************************************************
  TYPE colvar_type
     INTEGER :: index,type_id
     REAL(kind=dp) :: mass,lambda,ss, ss0, ff_s,vv, vvp, &
                      cv_energy,delta_s,ff_hills
     LOGICAL :: do_wall_p,do_wall_m
     REAL(kind=dp) :: wall_p,wall_m
     REAL(kind=dp), DIMENSION (:,:), POINTER :: dsdr
     INTEGER, DIMENSION (:), POINTER :: i_atom
     INTEGER :: n_atom_s,i_subsys
     TYPE(dist_colvar_type) :: dist_param
     TYPE(coord_colvar_type) :: coord_param
     TYPE(torsion_colvar_type) :: torsion_param
  END TYPE colvar_type
!!***
  !****************************************************************************

CONTAINS

!!****f* colvar_types/colvar_init *
!!
!!   NAME
!!     colvar_init
!!
!!   FUNCTION
!!     initializes a colvar_param type
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - colvar: the colvat to initialize
!!     - index: the index of the colvar
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!*** **********************************************************************
SUBROUTINE colvar_init(colvar,index,error)
    TYPE(colvar_type), INTENT(inout)         :: colvar
    INTEGER, INTENT(in)                      :: index
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'colvar_init', &
      routinep = modulen//':'//routinen

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  colvar%index=index
  colvar%type_id=no_colvar_id
  colvar%mass=-HUGE(0.0_dp)
  colvar%lambda=-HUGE(0.0_dp)
  colvar%dist_param%i_at=-HUGE(0)
  colvar%dist_param%j_at=-HUGE(0)
  colvar%torsion_param%i_at_tors(1)=-HUGE(0)
  colvar%torsion_param%i_at_tors(2)=-HUGE(0)
  colvar%torsion_param%i_at_tors(3)=-HUGE(0)
  colvar%torsion_param%i_at_tors(4)=-HUGE(0)
  colvar%i_subsys=-HUGE(0)
  colvar%coord_param%n_atoms_from=-HUGE(0)
  colvar%coord_param%n_atoms_to=-HUGE(0)
  colvar%coord_param%nncrd=6
  colvar%coord_param%ndcrd=12
  colvar%coord_param%R_0=3.0_dp
  NULLIFY(colvar%coord_param%i_at_from,colvar%coord_param%i_at_to)

  colvar%ss=0.0_dp
  colvar%ss0=0.0_dp
  colvar%ff_s=0.0_dp
  colvar%vv=0.0_dp
  colvar%vvp=0.0_dp
  colvar%wall_p=0.0_dp
  colvar%wall_m=0.0_dp
END SUBROUTINE colvar_init
!***************************************************************************

!!****f* colvar_types/colvar_dealloc_ref *
!!
!!   NAME
!!     colvar_dealloc_ref
!!
!!   FUNCTION
!!     releases the memory that might have been allocated by the colvar
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - colvar: the colvar to deallocate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!*** **********************************************************************
SUBROUTINE colvar_dealloc_ref(colvar,error)
    TYPE(colvar_type), INTENT(inout)         :: colvar
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'colvar_dealloc_ref', &
      routinep = modulen//':'//routinen

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  SELECT CASE(colvar%type_id)
  CASE (dist_colvar_id)
     DEALLOCATE(colvar%dsdr, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
     DEALLOCATE(colvar%i_atom, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)

     colvar%dist_param%i_at=-HUGE(0)
     colvar%dist_param%j_at=-HUGE(0)
  CASE (coord_colvar_id)
     DEALLOCATE(colvar%dsdr, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
     DEALLOCATE(colvar%i_atom, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)

     colvar%coord_param%n_atoms_from=-HUGE(0)
     colvar%coord_param%n_atoms_to=-HUGE(0)
     DEALLOCATE(colvar%coord_param%i_at_from,&
          colvar%coord_param%i_at_to,stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  CASE (torsion_colvar_id)
     DEALLOCATE(colvar%dsdr, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
     DEALLOCATE(colvar%i_atom, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)

     colvar%torsion_param%i_at_tors(1)=-HUGE(0)
     colvar%torsion_param%i_at_tors(2)=-HUGE(0)
     colvar%torsion_param%i_at_tors(3)=-HUGE(0)
     colvar%torsion_param%i_at_tors(4)=-HUGE(0)
  CASE(no_colvar_id)
  CASE default
     CPAssert(.FALSE.,cp_failure_level,routinep,error,failure)
  END SELECT
  colvar%index=-HUGE(0)
  colvar%type_id=no_colvar_id
  colvar%mass=-HUGE(0.0_dp)
  colvar%lambda=-HUGE(0.0_dp)
  colvar%delta_s=-HUGE(0.0_dp)
END SUBROUTINE colvar_dealloc_ref

END MODULE colvar_types
