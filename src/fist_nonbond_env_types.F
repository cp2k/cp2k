!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****s* cp2k/fist_nonbond_env_types [1.0] *
!!
!!   NAME
!!     fist_nonbond_env_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     HAF
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE fist_nonbond_env_types
  USE cell_types,                      ONLY: cell_release,&
                                             cell_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_deallocate,&
                                             fist_neighbor_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE pair_potential_types,            ONLY: pair_potential_pp_release,&
                                             pair_potential_pp_type,&
                                             tersoff_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_nonbond_env_types'
  INTEGER, PRIVATE, SAVE :: last_fist_nonbond_env_id_nr=0
  PUBLIC :: fist_nonbond_env_type, fist_nonbond_env_set, &
            fist_nonbond_env_get, fist_nonbond_env_create, &
            fist_nonbond_env_retain, fist_nonbond_env_release, pos_type

  TYPE pos_type
    REAL (KIND=dp) :: r(3)
  END TYPE

  TYPE fist_nonbond_env_type
    INTEGER :: ref_count, id_nr, natom_types
    INTEGER :: counter, last_update, num_update
    LOGICAL :: do_nonbonded14, do_nonbonded
    CHARACTER ( len = default_string_length )   :: unit_type
    REAL (KIND=dp)                              :: lup, aup
    REAL (KIND=dp)                              :: long_range_correction
    REAL (KIND=dp), DIMENSION ( :, : ), POINTER :: rlist_cutsq, rlist_lowsq, ij_kind_full_fac
    TYPE(fist_neighbor_type), POINTER           :: nonbonded
    TYPE(pair_potential_pp_type), POINTER       :: potparm
    TYPE(cell_type), POINTER                    :: cell_last_update
    TYPE(pos_type), DIMENSION(:), POINTER       :: r_last_update, r_last_update_pbc
    TYPE(pos_type), DIMENSION(:), POINTER       :: rshell_last_update_pbc
    TYPE(pos_type), DIMENSION(:), POINTER       :: rcore_last_update_pbc
  END TYPE fist_nonbond_env_type

!!***
!******************************************************************************
CONTAINS

!!****f* fist_nonbond_env_types/fist_nonbond_env_get [1.0] *
!!
!!   NAME
!!     fist_nonbond_env_get
!!
!!   SYNOPSIS
!!     Subroutine fist_nonbond_env_set (fist_nonbond_env )
!!       Type(fist_nonbond_env_type), Pointer:: fist_nonbond_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine fist_nonbond_env_create
!!
!!   FUNCTION
!!     sets a fist_nonbond_env
!!
!!   ARGUMENTS
!!     - fist_nonbond_env: the object to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE fist_nonbond_env_get (fist_nonbond_env, potparm, nonbonded,&
       rlist_cutsq, rlist_lowsq, aup, lup, r_last_update, r_last_update_pbc,&
       rshell_last_update_pbc, rcore_last_update_pbc,cell_last_update,num_update,&
       last_update, counter, natom_types, long_range_correction, ij_kind_full_fac,&
       error )

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(pair_potential_pp_type), OPTIONAL, &
      POINTER                                :: potparm
    TYPE(fist_neighbor_type), OPTIONAL, &
      POINTER                                :: nonbonded
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: rlist_cutsq, rlist_lowsq
    REAL(KIND=dp), OPTIONAL                  :: aup, lup
    TYPE(pos_type), DIMENSION(:), OPTIONAL, &
      POINTER                                :: r_last_update, &
                                                r_last_update_pbc, &
                                                rshell_last_update_pbc, &
                                                rcore_last_update_pbc
    TYPE(cell_type), OPTIONAL, POINTER       :: cell_last_update
    INTEGER, OPTIONAL                        :: num_update, last_update, &
                                                counter, natom_types
    REAL(KIND=dp), OPTIONAL                  :: long_range_correction
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: ij_kind_full_fac
    TYPE(cp_error_type), INTENT(inout)       :: error

    IF ( PRESENT ( potparm ) ) potparm => fist_nonbond_env % potparm
    IF ( PRESENT ( rlist_cutsq ) ) rlist_cutsq => fist_nonbond_env % rlist_cutsq
    IF ( PRESENT ( rlist_lowsq ) ) rlist_lowsq => fist_nonbond_env % rlist_lowsq
    IF ( PRESENT ( ij_kind_full_fac ) ) ij_kind_full_fac => fist_nonbond_env % ij_kind_full_fac
    IF ( PRESENT ( nonbonded ) ) nonbonded => fist_nonbond_env % nonbonded
    IF ( PRESENT ( r_last_update ) ) &
         r_last_update => fist_nonbond_env % r_last_update
    IF ( PRESENT ( r_last_update_pbc ) ) &
         r_last_update_pbc => fist_nonbond_env % r_last_update_pbc
    IF ( PRESENT ( rshell_last_update_pbc ) ) &
         rshell_last_update_pbc => fist_nonbond_env % rshell_last_update_pbc
    IF ( PRESENT ( rcore_last_update_pbc ) ) &
         rcore_last_update_pbc => fist_nonbond_env % rcore_last_update_pbc
    IF ( PRESENT ( cell_last_update ) ) &
         cell_last_update => fist_nonbond_env % cell_last_update
    IF ( PRESENT ( lup ) ) lup = fist_nonbond_env % lup 
    IF ( PRESENT ( aup ) ) aup = fist_nonbond_env % aup
    IF ( PRESENT ( natom_types ) ) natom_types = fist_nonbond_env % natom_types
    IF ( PRESENT ( counter ) ) counter = fist_nonbond_env % counter
    IF ( PRESENT ( last_update ) ) last_update = fist_nonbond_env % last_update
    IF ( PRESENT ( num_update ) ) num_update = fist_nonbond_env % num_update
    IF ( PRESENT ( long_range_correction ) ) &
         long_range_correction = fist_nonbond_env % long_range_correction 
  END SUBROUTINE fist_nonbond_env_get

!!****f* fist_nonbond_env_types/fist_nonbond_env_set [1.0] *
!!
!!   NAME
!!     fist_nonbond_env_set
!!
!!   SYNOPSIS
!!     Subroutine fist_nonbond_env_set (fist_nonbond_env )
!!       Type(fist_nonbond_env_type), Pointer:: fist_nonbond_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine fist_nonbond_env_create
!!
!!   FUNCTION
!!     sets a fist_nonbond_env
!!
!!   ARGUMENTS
!!     - fist_nonbond_env: the object to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE fist_nonbond_env_set ( fist_nonbond_env,  potparm, rlist_cutsq,&
       rlist_lowsq, nonbonded, aup, lup, r_last_update, r_last_update_pbc,&
       rshell_last_update_pbc, rcore_last_update_pbc, cell_last_update, num_update,&
       last_update, counter, natom_types, long_range_correction, error )

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(pair_potential_pp_type), OPTIONAL, &
      POINTER                                :: potparm
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: rlist_cutsq, rlist_lowsq
    TYPE(fist_neighbor_type), OPTIONAL, &
      POINTER                                :: nonbonded
    REAL(KIND=dp), OPTIONAL                  :: aup, lup
    TYPE(pos_type), DIMENSION(:), OPTIONAL, &
      POINTER                                :: r_last_update, &
                                                r_last_update_pbc, &
                                                rshell_last_update_pbc, &
                                                rcore_last_update_pbc
    TYPE(cell_type), OPTIONAL, POINTER       :: cell_last_update
    INTEGER, OPTIONAL                        :: num_update, last_update, &
                                                counter, natom_types
    REAL(KIND=dp), OPTIONAL                  :: long_range_correction
    TYPE(cp_error_type), INTENT(inout)       :: error

    IF ( PRESENT ( potparm ) ) THEN
       fist_nonbond_env % potparm => potparm
    ENDIF
    IF ( PRESENT ( rlist_cutsq ) ) fist_nonbond_env % rlist_cutsq => rlist_cutsq
    IF ( PRESENT ( rlist_lowsq ) ) fist_nonbond_env % rlist_lowsq => rlist_lowsq
    IF ( PRESENT ( nonbonded ) ) fist_nonbond_env % nonbonded => nonbonded
    IF ( PRESENT ( r_last_update ) ) &
         fist_nonbond_env % r_last_update => r_last_update
    IF ( PRESENT ( r_last_update_pbc ) ) &
         fist_nonbond_env % r_last_update_pbc => r_last_update_pbc
    IF ( PRESENT ( rshell_last_update_pbc ) ) &
         fist_nonbond_env % rshell_last_update_pbc => rshell_last_update_pbc
    IF ( PRESENT ( rcore_last_update_pbc ) ) &
         fist_nonbond_env % rcore_last_update_pbc => rcore_last_update_pbc
    IF ( PRESENT ( cell_last_update ) ) &
         fist_nonbond_env % cell_last_update => cell_last_update
    IF ( PRESENT ( lup ) ) fist_nonbond_env % lup = lup
    IF ( PRESENT ( aup ) ) fist_nonbond_env % aup = aup
    IF ( PRESENT ( natom_types ) ) fist_nonbond_env % natom_types = natom_types
    IF ( PRESENT ( counter ) ) fist_nonbond_env % counter = counter
    IF ( PRESENT ( last_update ) ) fist_nonbond_env % last_update = last_update
    IF ( PRESENT ( num_update ) ) fist_nonbond_env % num_update = num_update
    IF ( PRESENT ( long_range_correction ) ) &
         fist_nonbond_env % long_range_correction = long_range_correction
  END SUBROUTINE fist_nonbond_env_set

!!****f* fist_nonbond_env_types/fist_nonbond_env_create [1.0] *
!!
!!   NAME
!!     fist_nonbond_env_create
!!
!!   SYNOPSIS
!!     Subroutine fist_nonbond_env_creat (fist_nonbond_env )
!!       Type(fist_nonbond_env_type), Pointer:: fist_nonbond_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine fist_nonbond_env_create
!!
!!   FUNCTION
!!     allocates and intitializes a fist_nonbond_env
!!
!!   ARGUMENTS
!!     - fist_nonbond_env: the object to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE fist_nonbond_env_create(fist_nonbond_env, potparm, nonbonded14,&
       nonbonded, mm_section, ewald_env, error)
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    LOGICAL, INTENT(IN)                      :: nonbonded14, nonbonded
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_nonbond_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(fist_nonbond_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL init_fist_nonbond_env( fist_nonbond_env, potparm, nonbonded14,&
            nonbonded, mm_section, ewald_env, error)
    END IF
  END SUBROUTINE fist_nonbond_env_create

!!****f* fist_nonbond_env_types/init_fist_nonbond_env [1.0] *
!!
!!   NAME
!!     init_fist_nonbond_env
!!
!!
!!   FUNCTION
!!     Purpose: Initialise the FIST nonbond environment.
!!
!!   ARGUMENTS
!!     - fist_nonbond_env: the object to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!    
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  SUBROUTINE init_fist_nonbond_env ( fist_nonbond_env, potparm, nonbonded14,&
     nonbonded, mm_section, ewald_env, error)

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    LOGICAL, INTENT(IN)                      :: nonbonded14, nonbonded
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_fist_nonbond_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idim, jdim, natom_types, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ewald_rcut, fac, rcut, tmp, &
                                                verlet_skin

    NULLIFY (fist_nonbond_env%nonbonded)
    NULLIFY (fist_nonbond_env%r_last_update)
    NULLIFY (fist_nonbond_env%r_last_update_pbc)
    NULLIFY (fist_nonbond_env%rshell_last_update_pbc)
    NULLIFY (fist_nonbond_env%rcore_last_update_pbc)
    NULLIFY (fist_nonbond_env%cell_last_update)
    NULLIFY (fist_nonbond_env%rlist_cutsq)
    NULLIFY (fist_nonbond_env%rlist_lowsq)
    NULLIFY (fist_nonbond_env%ij_kind_full_fac)
    fist_nonbond_env % unit_type = "ANGSTROM"
    fist_nonbond_env % do_nonbonded14 = nonbonded14
    fist_nonbond_env % do_nonbonded = nonbonded
    fist_nonbond_env % lup = 0
    fist_nonbond_env % aup = 0
    fist_nonbond_env % counter = 0
    fist_nonbond_env % last_update = 0
    fist_nonbond_env % num_update = 0
    fist_nonbond_env % long_range_correction = 0
    fist_nonbond_env%ref_count=1
    fist_nonbond_env%id_nr=last_fist_nonbond_env_id_nr
    last_fist_nonbond_env_id_nr=last_fist_nonbond_env_id_nr+1
    IF (nonbonded) THEN
       CALL ewald_env_get(ewald_env, rcut=ewald_rcut, error=error)
       CALL section_vals_val_get(mm_section,"FORCEFIELD%NEIGHBOR_LISTS%VERLET_SKIN",&
            r_val=verlet_skin,error=error)
       natom_types = SIZE ( potparm%pot, 1 )
       ALLOCATE ( fist_nonbond_env%rlist_cutsq(natom_types, natom_types), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( fist_nonbond_env%rlist_lowsq(natom_types, natom_types), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( fist_nonbond_env%ij_kind_full_fac(natom_types, natom_types), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       fist_nonbond_env%ij_kind_full_fac = 1.0_dp
       DO idim = 1, natom_types
          DO jdim = idim, natom_types
             fac  = potparm%pot(idim,jdim)%pot%spl_f%rscale(2)
             rcut = SQRT(potparm%pot(idim,jdim)%pot%rcutsq) + verlet_skin
             rcut = MAX(rcut, ewald_rcut)
             fist_nonbond_env % rlist_cutsq(idim,jdim) = rcut
             fist_nonbond_env % rlist_cutsq(jdim,idim) = rcut
             rcut = fac / (potparm%pot(idim,jdim)%pot%pair_spline_data(1)%spline_data%xn)
             tmp  = rcut*(1.06_dp)**2 ! 1.06_dp in order to have 1/2 Emax_spline
             fist_nonbond_env%rlist_lowsq(idim,jdim) = tmp
             fist_nonbond_env%rlist_lowsq(jdim,idim) = tmp
             ! In case of manybody potential the NL will be full
             IF (ANY(potparm%pot(idim,jdim)%pot%type==tersoff_type)) THEN
                fist_nonbond_env%ij_kind_full_fac(idim,jdim) = 0.5_dp
                fist_nonbond_env%ij_kind_full_fac(idim,jdim) = 0.5_dp
             END IF
          END DO
       END DO
       fist_nonbond_env%potparm => potparm
       fist_nonbond_env%natom_types = natom_types
    ELSE
       NULLIFY(fist_nonbond_env % potparm)
    END IF
  END SUBROUTINE init_fist_nonbond_env

!!****f* fist_nonbond_env_types/fist_nonbond_env_retain [1.0] *
!!
!!   NAME
!!     fist_nonbond_env_retain
!!
!!   SYNOPSIS
!!     Subroutine fist_nonbond_env_retain(fist_nonbond_env, error)
!!       Type(fist_nonbond_env_type), Pointer:: fist_nonbond_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine fist_nonbond_env_retain
!!
!!   FUNCTION
!!     retains the given fist_nonbond_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - fist_nonbond_env: the object to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE fist_nonbond_env_retain(fist_nonbond_env,error)
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_nonbond_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(fist_nonbond_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(fist_nonbond_env%ref_count>0,cp_failure_level,routineP,error,failure)
       fist_nonbond_env%ref_count=fist_nonbond_env%ref_count+1
    END IF
  END SUBROUTINE fist_nonbond_env_retain

!!****f* fist_nonbond_env_types/fist_nonbond_env_release [1.0] *
!!
!!   NAME
!!     fist_nonbond_env_release
!!
!!   SYNOPSIS
!!     Subroutine fist_nonbond_env_release(fist_nonbond_env, error)
!!       Type(fist_nonbond_env_type), Pointer:: fist_nonbond_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine fist_nonbond_env_release
!!
!!   FUNCTION
!!     releases the given fist_nonbond_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - fist_nonbond_env: the object to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE fist_nonbond_env_release(fist_nonbond_env,error)
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_nonbond_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(fist_nonbond_env)) THEN
       CPPrecondition(fist_nonbond_env%ref_count>0,cp_failure_level,routineP,error,failure)
       fist_nonbond_env%ref_count=fist_nonbond_env%ref_count-1
       IF (fist_nonbond_env%ref_count<1) THEN
          IF (ASSOCIATED(fist_nonbond_env % nonbonded)) THEN
             CALL fist_neighbor_deallocate (fist_nonbond_env % nonbonded, error )
          END IF
          ! Release potparm
          CALL pair_potential_pp_release ( fist_nonbond_env % potparm, error )
          IF (ASSOCIATED(fist_nonbond_env % r_last_update)) THEN
             DEALLOCATE (fist_nonbond_env % r_last_update,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env % r_last_update_pbc)) THEN
             DEALLOCATE (fist_nonbond_env % r_last_update_pbc,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env % rshell_last_update_pbc)) THEN
             DEALLOCATE (fist_nonbond_env % rshell_last_update_pbc,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env % rcore_last_update_pbc)) THEN
             DEALLOCATE (fist_nonbond_env % rcore_last_update_pbc,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env % cell_last_update)) THEN
             CALL cell_release(fist_nonbond_env % cell_last_update,error)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env % ij_kind_full_fac)) THEN
             DEALLOCATE (fist_nonbond_env % ij_kind_full_fac,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env % rlist_cutsq)) THEN
             DEALLOCATE (fist_nonbond_env % rlist_cutsq,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env % rlist_lowsq)) THEN
             DEALLOCATE (fist_nonbond_env % rlist_lowsq,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          DEALLOCATE(fist_nonbond_env, STAT=stat)
          CPPrecondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE fist_nonbond_env_release

END MODULE fist_nonbond_env_types
!***************************************************************************
