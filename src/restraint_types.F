!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****h* cp2k/restraint_types *
!!
!!   NAME
!!     restraint_types
!!
!!   FUNCTION
!!     types used in restraint
!!
!!   NOTES
!!
!!   AUTHOR
!!     ikuo
!!
!!   MODIFICATION HISTORY
!!     07.2006 [ikuo]
!!
!!   SOURCE
!****************************************************************************
MODULE restraint_types
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: restraint_env_type 
  PUBLIC :: restraint_env_create
  PUBLIC :: restraint_env_release
  PUBLIC :: restraint_env_retain
  PUBLIC :: restraint_env_read
  PUBLIC :: restraint_env_write

  PUBLIC :: fp_type
  PUBLIC :: harmonic_type
  PUBLIC :: dist_type
  PUBLIC :: angle_type
  PUBLIC :: torsion_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'restraint_types'

  TYPE fp_type
    ! input related objects
    LOGICAL                          :: in_use

    INTEGER                          :: central_atom
    INTEGER, DIMENSION(:), POINTER   :: inner_atoms, outer_atoms
    REAL(KIND=dp)                    :: inner_radius, outer_radius
    REAL(KIND=dp)                    :: strength, smooth_width
    LOGICAL                          :: bias
    REAL(KIND=dp)                    :: temperature
    TYPE(section_vals_type), POINTER :: print_section
    ! computed during runs
    INTEGER                          :: i1,i2,o1,o2 
    REAL(KIND=dp)                    :: ri1,ri2,ro1,ro2
    REAL(KIND=dp)                    :: weight, comb_weight, bias_weight
    REAL(KIND=dp)                    :: energy, bias_energy, restraint_energy
  END TYPE fp_type

  TYPE harmonic_type
    LOGICAL                          :: in_use
    REAL(KIND=dp)                    :: weight, comb_weight, bias_weight
    REAL(KIND=dp)                    :: energy, bias_energy, restraint_energy
    REAL(KIND=dp),DIMENSION(:,:),POINTER &
                                     :: r_orig
    REAL(KIND=dp),DIMENSION(:),POINTER &
                                     :: k
  END TYPE harmonic_type

  TYPE dist_type
    LOGICAL                          :: in_use
    REAL(KIND=dp)                    :: weight, comb_weight, bias_weight
    REAL(KIND=dp)                    :: energy, bias_energy, restraint_energy
    INTEGER                          :: nrestraint
    INTEGER, DIMENSION(:), POINTER   :: atom_a, atom_b
    REAL(KIND=dp),DIMENSION(:),POINTER &
                                     :: r0
    REAL(KIND=dp),DIMENSION(:),POINTER &
                                     :: k
  END TYPE dist_type

  TYPE angle_type
    LOGICAL                          :: in_use
    REAL(KIND=dp)                    :: weight, comb_weight, bias_weight
    REAL(KIND=dp)                    :: energy, bias_energy, restraint_energy
    INTEGER                          :: nrestraint
    INTEGER, DIMENSION(:), POINTER   :: atom_a, atom_b, atom_c
    REAL(KIND=dp),DIMENSION(:),POINTER &
                                     :: theta0
    REAL(KIND=dp),DIMENSION(:),POINTER &
                                     :: k
  END TYPE angle_type

  TYPE torsion_type
    LOGICAL                          :: in_use
    REAL(KIND=dp)                    :: weight, comb_weight, bias_weight
    REAL(KIND=dp)                    :: energy, bias_energy, restraint_energy
    INTEGER                          :: nrestraint
    INTEGER, DIMENSION(:), POINTER   :: atom_a, atom_b, atom_c, atom_d
    REAL(KIND=dp),DIMENSION(:),POINTER &
                                     :: phi0
    REAL(KIND=dp),DIMENSION(:),POINTER &
                                     :: k
    INTEGER(KIND=dp),DIMENSION(:),POINTER &
                                     :: imul
  END TYPE torsion_type

  TYPE restraint_env_type
    LOGICAL                        :: in_use
    INTEGER                        :: ref_count
    TYPE (fp_type), POINTER        :: fp_env
    TYPE (harmonic_type), POINTER  :: harmonic_env
    TYPE (dist_type), POINTER      :: dist_env
    TYPE (angle_type), POINTER     :: angle_env
    TYPE (torsion_type), POINTER   :: torsion_env
    REAL(KIND=dp)                  :: energy
  END TYPE restraint_env_type

CONTAINS

!!****f* restraint_types/restraint_env_create/retain/release *
!!
!!   NAME
!!      restraint_env_create/retain/release
!!   FUNCTION
!!      create retain release the restraint environment
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE restraint_env_create(restraint_env,error)
    TYPE(restraint_env_type), POINTER        :: restraint_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'restraint_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.
 
    CPPrecondition(.NOT.ASSOCIATED(restraint_env),cp_failure_level,routineP,error,failure)
    ALLOCATE(restraint_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    restraint_env%ref_count=1
    restraint_env%in_use=.FALSE.
    restraint_env%energy=0.0_dp
    NULLIFY(restraint_env%fp_env)
    NULLIFY(restraint_env%harmonic_env)
    NULLIFY(restraint_env%dist_env)
    NULLIFY(restraint_env%angle_env)
    NULLIFY(restraint_env%torsion_env)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deal with fp
    ALLOCATE(restraint_env%fp_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    restraint_env%fp_env%in_use=.FALSE.
    restraint_env%fp_env%bias=.FALSE.
    NULLIFY(restraint_env%fp_env%inner_atoms)
    NULLIFY(restraint_env%fp_env%outer_atoms)
    NULLIFY(restraint_env%fp_env%print_section)
    restraint_env%fp_env%strength=0.0_dp
    restraint_env%fp_env%smooth_width=0.0_dp
    restraint_env%fp_env%temperature=0.0_dp
    restraint_env%fp_env%ri1=0.0_dp
    restraint_env%fp_env%ri2=0.0_dp
    restraint_env%fp_env%ro1=0.0_dp
    restraint_env%fp_env%ro2=0.0_dp
    restraint_env%fp_env%weight=0.0_dp
    restraint_env%fp_env%comb_weight=0.0_dp
    restraint_env%fp_env%bias_weight=0.0_dp
    restraint_env%fp_env%energy=0.0_dp
    restraint_env%fp_env%bias_energy=0.0_dp
    restraint_env%fp_env%restraint_energy=0.0_dp
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deal with harmonic
    ALLOCATE(restraint_env%harmonic_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    restraint_env%harmonic_env%in_use=.FALSE.
    NULLIFY(restraint_env%harmonic_env%r_orig)
    NULLIFY(restraint_env%harmonic_env%k)
    restraint_env%harmonic_env%energy=0.0_dp
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deal with dist
    ALLOCATE(restraint_env%dist_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    restraint_env%dist_env%in_use=.FALSE.
    restraint_env%dist_env%energy=0.0_dp
    NULLIFY(restraint_env%dist_env%atom_a)
    NULLIFY(restraint_env%dist_env%atom_b)
    NULLIFY(restraint_env%dist_env%r0)
    NULLIFY(restraint_env%dist_env%k)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deal with angle
    ALLOCATE(restraint_env%angle_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    restraint_env%angle_env%in_use=.FALSE.
    restraint_env%angle_env%energy=0.0_dp
    NULLIFY(restraint_env%angle_env%atom_a)
    NULLIFY(restraint_env%angle_env%atom_b)
    NULLIFY(restraint_env%angle_env%atom_c)
    NULLIFY(restraint_env%angle_env%theta0)
    NULLIFY(restraint_env%angle_env%k)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deal with torsion
    ALLOCATE(restraint_env%torsion_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    restraint_env%torsion_env%in_use=.FALSE.
    restraint_env%torsion_env%energy=0.0_dp
    NULLIFY(restraint_env%torsion_env%atom_a)
    NULLIFY(restraint_env%torsion_env%atom_b)
    NULLIFY(restraint_env%torsion_env%atom_c)
    NULLIFY(restraint_env%torsion_env%atom_d)
    NULLIFY(restraint_env%torsion_env%phi0)
    NULLIFY(restraint_env%torsion_env%k)
    NULLIFY(restraint_env%torsion_env%imul)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
    CALL timestop(handle)
  
  END SUBROUTINE restraint_env_create
  
  SUBROUTINE restraint_env_release(restraint_env,error)
    TYPE(restraint_env_type), POINTER        :: restraint_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'restraint_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

   CALL timeset(routineN,handle)
   failure=.FALSE.

   IF (ASSOCIATED(restraint_env)) THEN

     CPPrecondition(restraint_env%ref_count>0,cp_failure_level,routineP,error,failure)
     restraint_env%ref_count=restraint_env%ref_count-1

     IF (restraint_env%ref_count==0) THEN

       ! Deal with the fp_env
       IF (ASSOCIATED(restraint_env%fp_env)) THEN
         IF (ASSOCIATED(restraint_env%fp_env%inner_atoms)) DEALLOCATE(restraint_env%fp_env%inner_atoms)
         IF (ASSOCIATED(restraint_env%fp_env%outer_atoms)) DEALLOCATE(restraint_env%fp_env%outer_atoms)
         IF (ASSOCIATED(restraint_env%fp_env%print_section)) CALL section_vals_release(restraint_env%fp_env%print_section)
         restraint_env%fp_env%in_use=.FALSE.
         DEALLOCATE(restraint_env%fp_env)
       ENDIF

       ! Deal with the harmonic_env
       IF (ASSOCIATED(restraint_env%harmonic_env)) THEN
         IF (ASSOCIATED(restraint_env%harmonic_env%r_orig)) DEALLOCATE(restraint_env%harmonic_env%r_orig)
         IF (ASSOCIATED(restraint_env%harmonic_env%k)) DEALLOCATE(restraint_env%harmonic_env%k)
         restraint_env%harmonic_env%in_use=.FALSE.
         DEALLOCATE(restraint_env%harmonic_env)
       ENDIF

       ! Deal with the dist_env
       IF (ASSOCIATED(restraint_env%dist_env)) THEN
         IF (ASSOCIATED(restraint_env%dist_env%r0)) DEALLOCATE(restraint_env%dist_env%r0)
         IF (ASSOCIATED(restraint_env%dist_env%k)) DEALLOCATE(restraint_env%dist_env%k)
         IF (ASSOCIATED(restraint_env%dist_env%atom_a)) DEALLOCATE(restraint_env%dist_env%atom_a)
         IF (ASSOCIATED(restraint_env%dist_env%atom_b)) DEALLOCATE(restraint_env%dist_env%atom_b)
         restraint_env%dist_env%in_use=.FALSE.
         DEALLOCATE(restraint_env%dist_env)
       ENDIF

       ! Deal with the angle_env
       IF (ASSOCIATED(restraint_env%angle_env)) THEN
         IF (ASSOCIATED(restraint_env%angle_env%theta0)) DEALLOCATE(restraint_env%angle_env%theta0)
         IF (ASSOCIATED(restraint_env%angle_env%k)) DEALLOCATE(restraint_env%angle_env%k)
         IF (ASSOCIATED(restraint_env%angle_env%atom_a)) DEALLOCATE(restraint_env%angle_env%atom_a)
         IF (ASSOCIATED(restraint_env%angle_env%atom_b)) DEALLOCATE(restraint_env%angle_env%atom_b)
         IF (ASSOCIATED(restraint_env%angle_env%atom_c)) DEALLOCATE(restraint_env%angle_env%atom_c)
         restraint_env%angle_env%in_use=.FALSE.
         DEALLOCATE(restraint_env%angle_env)
       ENDIF

       ! Deal with the torsion_env
       IF (ASSOCIATED(restraint_env%torsion_env)) THEN
         IF (ASSOCIATED(restraint_env%torsion_env%phi0)) DEALLOCATE(restraint_env%torsion_env%phi0)
         IF (ASSOCIATED(restraint_env%torsion_env%k)) DEALLOCATE(restraint_env%torsion_env%k)
         IF (ASSOCIATED(restraint_env%torsion_env%imul)) DEALLOCATE(restraint_env%torsion_env%imul)
         IF (ASSOCIATED(restraint_env%torsion_env%atom_a)) DEALLOCATE(restraint_env%torsion_env%atom_a)
         IF (ASSOCIATED(restraint_env%torsion_env%atom_b)) DEALLOCATE(restraint_env%torsion_env%atom_b)
         IF (ASSOCIATED(restraint_env%torsion_env%atom_c)) DEALLOCATE(restraint_env%torsion_env%atom_c)
         IF (ASSOCIATED(restraint_env%torsion_env%atom_d)) DEALLOCATE(restraint_env%torsion_env%atom_d)
         restraint_env%torsion_env%in_use=.FALSE.
         DEALLOCATE(restraint_env%torsion_env)
       ENDIF

       restraint_env%in_use=.FALSE.
       DEALLOCATE(restraint_env)

     ENDIF

   ENDIF
   CALL timestop(handle)

  END SUBROUTINE restraint_env_release
  
  SUBROUTINE restraint_env_retain(restraint_env,error)
    TYPE(restraint_env_type), POINTER                   :: restraint_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'restraint_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

      failure=.FALSE.
      CPPrecondition(ASSOCIATED(restraint_env),cp_failure_level,routineP,error,failure)
      restraint_env%ref_count=restraint_env%ref_count+1
  
  END SUBROUTINE restraint_env_retain

!!****f* restraint_types/restraint_env_read *
!!
!!   NAME
!!     restraint_env_read
!!
!!   FUNCTION
!!     reads the corresponding input section and stores it in the restraint_env
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     07.2006 ikuo created 
!!
!!*** **********************************************************************
  SUBROUTINE restraint_env_read(restraint_env,restraint_section,error)
    TYPE(restraint_env_type), POINTER        :: restraint_env
    TYPE(section_vals_type), POINTER         :: restraint_section,subsection
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'restraint_env_read', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(restraint_env%fp_env),cp_failure_level,routineP,error,failure)

    CALL section_vals_get(restraint_section,explicit=restraint_env%in_use,error=error)

    subsection => section_vals_get_subs_vals(restraint_section,&
       "FLEXIBLE_PARTITIONING",error=error)
    CALL fp_env_read(restraint_env%fp_env,subsection,error)

    subsection => section_vals_get_subs_vals(restraint_section,&
       "HARMONIC",error=error)
    CALL harmonic_env_read(restraint_env%harmonic_env,subsection,error)

    subsection => section_vals_get_subs_vals(restraint_section,&
       "DIST",error=error)
    CALL dist_env_read(restraint_env%dist_env,subsection,error)

    subsection => section_vals_get_subs_vals(restraint_section,&
       "ANGLE",error=error)
    CALL angle_env_read(restraint_env%angle_env,subsection,error)

    subsection => section_vals_get_subs_vals(restraint_section,&
       "TORSION",error=error)
    CALL torsion_env_read(restraint_env%torsion_env,subsection,error)

    CALL timestop(handle)
  END SUBROUTINE restraint_env_read

!!****f* restraint_types/fp_env_read *
!!
!!   NAME
!!     fp_env_read
!!
!!   FUNCTION
!!     reads the corresponding input section and stores it in the fp_env
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     04.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE fp_env_read(fp_env,fp_section,error)
    TYPE(fp_type), POINTER                   :: fp_env
    TYPE(section_vals_type), POINTER         :: fp_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fp_env_read', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(fp_env),cp_failure_level,routineP,error,failure)
    CALL section_vals_get(fp_section,explicit=fp_env%in_use,error=error)
    IF (fp_env%in_use) THEN
       CALL section_vals_val_get(fp_section,"CENTRAL_ATOM",i_val=fp_env%central_atom,error=error)

       CALL section_vals_val_get(fp_section,"INNER_ATOMS",i_vals=tmplist,error=error)
       ALLOCATE(fp_env%inner_atoms(SIZE(tmplist,1)))
       fp_env%inner_atoms=tmplist

       CALL section_vals_val_get(fp_section,"OUTER_ATOMS",i_vals=tmplist,error=error)
       ALLOCATE(fp_env%outer_atoms(SIZE(tmplist,1)))
       fp_env%outer_atoms=tmplist

       CALL section_vals_val_get(fp_section,"INNER_RADIUS",r_val=fp_env%inner_radius,error=error)
       CALL section_vals_val_get(fp_section,"OUTER_RADIUS",r_val=fp_env%outer_radius,error=error)
       CALL section_vals_val_get(fp_section,"STRENGTH",r_val=fp_env%strength,error=error)
       CALL section_vals_val_get(fp_section,"SMOOTH_WIDTH",r_val=fp_env%smooth_width,error=error)
       CALL section_vals_val_get(fp_section,"BIAS",l_val=fp_env%bias,error=error)
       CALL section_vals_val_get(fp_section,"TEMPERATURE",r_val=fp_env%temperature,error=error)

       fp_env%print_section=>section_vals_get_subs_vals(fp_section,"WEIGHTS",error=error)
       CALL section_vals_retain(fp_env%print_section)
    ENDIF
    CALL timestop(handle)

  END SUBROUTINE fp_env_read



!!****f* restraint_types/harmonic_env_read *
!!
!!   NAME
!!     harmonic_env_read
!!
!!   FUNCTION
!!     reads the corresponding input section and stores it in the harmonic_env
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     07.2006 created [ikuo]
!!
!!*** **********************************************************************
  SUBROUTINE harmonic_env_read(harmonic_env,harmonic_section,error)
    TYPE(harmonic_type), POINTER             :: harmonic_env
    TYPE(section_vals_type), POINTER         :: harmonic_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harmonic_env_read', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(harmonic_env),cp_failure_level,routineP,error,failure)
    CALL section_vals_get(harmonic_section,explicit=harmonic_env%in_use,error=error)
    IF (harmonic_env%in_use) THEN



    ENDIF
    CALL timestop(handle)

  END SUBROUTINE harmonic_env_read






!!****f* restraint_types/dist_env_read *
!!
!!   NAME
!!     dist_env_read
!!
!!   FUNCTION
!!     reads the corresponding input section and stores it in the dist_env
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     07.2006 created [ikuo]
!!
!!*** **********************************************************************
  SUBROUTINE dist_env_read(dist_env,dist_section,error)
    TYPE(dist_type), POINTER                 :: dist_env
    TYPE(section_vals_type), POINTER         :: dist_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dist_env_read', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(:), POINTER           :: tmplist
    INTEGER                                  :: handle,stat
    LOGICAL                                  :: failure
    INTEGER                                  :: ires,nres
    REAL(KIND=dp),POINTER,DIMENSION(:) :: r0,k

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(dist_env),cp_failure_level,routineP,error,failure)
    CALL section_vals_get(dist_section,explicit=dist_env%in_use,error=error)
    IF (dist_env%in_use) THEN
        
      CALL section_vals_get(dist_section,n_repetition=nres, error=error)
      dist_env%nrestraint = nres
    
      ALLOCATE(dist_env%atom_a(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(dist_env%atom_b(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(dist_env%r0(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(dist_env%k(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     
      DO ires=1,nres
        CALL section_vals_val_get(dist_section,"ATOMS",i_rep_section=ires,&
                                  i_vals=tmplist,error=error)
        dist_env%atom_a(ires) = tmplist(1)
        dist_env%atom_b(ires) = tmplist(2)

        CALL section_vals_val_get(dist_section,"K",i_rep_section=ires,&
                                  r_vals=k,error=error)
        CALL section_vals_val_get(dist_section,"R0",i_rep_section=ires,&
                                  r_vals=r0,error=error)
        dist_env%k(ires) = k(1) * 2.0_dp
        dist_env%r0(ires) = r0(1)
      END DO
     
    ENDIF
    CALL timestop(handle)

  END SUBROUTINE dist_env_read




!!****f* restraint_types/angle_env_read *
!!
!!   NAME
!!     angle_env_read
!!
!!   FUNCTION
!!     reads the corresponding input section and stores it in the angle_env
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     07.2006 created [ikuo]
!!
!!*** **********************************************************************
  SUBROUTINE angle_env_read(angle_env,angle_section,error)
    TYPE(angle_type), POINTER             :: angle_env
    TYPE(section_vals_type), POINTER         :: angle_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'angle_env_read', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle,stat
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: failure
    INTEGER                                  :: ires,nres
    REAL(KIND=dp),POINTER,DIMENSION(:) :: theta0,k

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(angle_env),cp_failure_level,routineP,error,failure)
    CALL section_vals_get(angle_section,explicit=angle_env%in_use,error=error)
    IF (angle_env%in_use) THEN

      CALL section_vals_get(angle_section,n_repetition=nres, error=error)
      angle_env%nrestraint = nres
    
      ALLOCATE(angle_env%atom_a(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(angle_env%atom_b(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(angle_env%atom_c(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(angle_env%theta0(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(angle_env%k(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     
      DO ires=1,nres
        CALL section_vals_val_get(angle_section,"ATOMS",i_rep_section=ires,&
                                  i_vals=tmplist,error=error)
        angle_env%atom_a(ires) = tmplist(1)
        angle_env%atom_b(ires) = tmplist(2)
        angle_env%atom_c(ires) = tmplist(3)

        CALL section_vals_val_get(angle_section,"K",i_rep_section=ires,&
                                  r_vals=k,error=error)
        CALL section_vals_val_get(angle_section,"THETA0",i_rep_section=ires,&
                                  r_vals=theta0,error=error)
        angle_env%k(ires) = k(1) * 2.0_dp
        angle_env%theta0(ires) = theta0(1)
      END DO

    ENDIF
    CALL timestop(handle)

  END SUBROUTINE angle_env_read




!!****f* restraint_types/torsion_env_read *
!!
!!   NAME
!!     torsion_env_read
!!
!!   FUNCTION
!!     reads the corresponding input section and stores it in the torsion_env
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     07.2006 created [ikuo]
!!
!!*** **********************************************************************
  SUBROUTINE torsion_env_read(torsion_env,torsion_section,error)
    TYPE(torsion_type), POINTER             :: torsion_env
    TYPE(section_vals_type), POINTER         :: torsion_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'torsion_env_read', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle,stat
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: failure
    INTEGER                                  :: ires,nres
    REAL(KIND=dp),POINTER,DIMENSION(:) :: phi0,k
    INTEGER,POINTER,DIMENSION(:) :: imul

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(torsion_env),cp_failure_level,routineP,error,failure)
    CALL section_vals_get(torsion_section,explicit=torsion_env%in_use,error=error)
    IF (torsion_env%in_use) THEN

      CALL section_vals_get(torsion_section,n_repetition=nres, error=error)
      torsion_env%nrestraint = nres
    
      ALLOCATE(torsion_env%atom_a(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(torsion_env%atom_b(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(torsion_env%atom_c(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(torsion_env%atom_d(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(torsion_env%phi0(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(torsion_env%k(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(torsion_env%imul(nres),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     
      DO ires=1,nres
        CALL section_vals_val_get(torsion_section,"ATOMS",i_rep_section=ires,&
                                  i_vals=tmplist,error=error)
        torsion_env%atom_a(ires) = tmplist(1)
        torsion_env%atom_b(ires) = tmplist(2)
        torsion_env%atom_c(ires) = tmplist(3)
        torsion_env%atom_d(ires) = tmplist(4)

        CALL section_vals_val_get(torsion_section,"K",i_rep_section=ires,&
                                  r_vals=k,error=error)
        CALL section_vals_val_get(torsion_section,"PHI0",i_rep_section=ires,&
                                  r_vals=phi0,error=error)
        CALL section_vals_val_get(torsion_section,"M",i_rep_section=ires,&
                                  i_vals=imul,error=error)
        torsion_env%k(ires) = k(1)
        torsion_env%phi0(ires) = phi0(1)
        torsion_env%imul(ires) = imul(1)
      END DO

    ENDIF
    CALL timestop(handle)

  END SUBROUTINE torsion_env_read








!!****f* restraint_types/restraint_env_write *
!!
!!   NAME
!!     restraint_env_write
!!
!!   FUNCTION
!!     writes information concerning the restraint_env to the output
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     07.2006 created [ikuo]
!!
!!*** **********************************************************************
  SUBROUTINE restraint_env_write(restraint_env,restraint_section,error)
    TYPE(restraint_env_type), POINTER        :: restraint_env
    TYPE(section_vals_type), POINTER         :: restraint_section
    TYPE(section_vals_type), POINTER         :: subsection
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'restraint_env_write', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    subsection => section_vals_get_subs_vals(restraint_section,&
       "FLEXIBLE_PARTITIONING",error=error)
    CALL fp_env_write(restraint_env%fp_env,subsection,error)

    !CALL harmonic_env_write(restraint_env%fp_env,restraint_section,error)
    !CALL dist_env_write(restraint_env%fp_env,restraint_section,error)
    !CALL angle_env_write(restraint_env%fp_env,restraint_section,error)
    !CALL torsion_env_write(restraint_env%fp_env,restraint_section,error)

    CALL timestop(handle)

  END SUBROUTINE restraint_env_write

!!****f* restraint_types/fp_env_write *
!!
!!   NAME
!!     fp_env_write
!!
!!   FUNCTION
!!     writes information concerning the fp_env to the output
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     04.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE fp_env_write(fp_env,fp_section,error)
    TYPE(fp_type), POINTER                   :: fp_env
    TYPE(section_vals_type), POINTER         :: fp_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fp_env_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, output_unit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: PRINT

    CALL timeset(routineN,handle)
    failure=.FALSE.
    logger => cp_error_get_logger(error)
    CPPrecondition(ASSOCIATED(fp_env),cp_failure_level,routineP,error,failure)

    IF (fp_env%in_use) THEN
       PRINT=>section_vals_get_subs_vals(fp_section,"CONTROL",error=error)
       output_unit=cp_print_key_unit_nr(logger,PRINT,"",extension=".Log",error=error)
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T79,A)")&
             "FP| Flexible partitioning is ","ON"
          WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
             "FP| Central atom ", fp_env%central_atom
          WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
             "FP| number of inner atoms", SIZE(fp_env%inner_atoms,1)
          WRITE (UNIT=output_unit,FMT="(1(T2,8I8))") fp_env%inner_atoms
          WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
             "FP| number of outer atoms", SIZE(fp_env%outer_atoms,1)
          WRITE (UNIT=output_unit,FMT="(1(T2,8I8))") fp_env%outer_atoms
          WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.10)")&
             "FP| inner radius [a.u.] ", fp_env%inner_radius
          WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.10)")&
             "FP| outer radius [a.u.] ", fp_env%outer_radius
          WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.10)")&
             "FP| reflecting restraint strength ", fp_env%strength
          IF (fp_env%bias) THEN
             WRITE (UNIT=output_unit,FMT="(T2,A,T79,A)")&
                     "FP| Flexible partitioning bias is " ,"ON"
             WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.10)")&
                     "FP| bias temperature [kT a.u.]", fp_env%temperature
             WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.10)")&
                     "FP| smooth width [a.u.] ", fp_env%smooth_width
          ELSE
             WRITE (UNIT=output_unit,FMT="(T2,A,T78,A)")&
                     "FP| Flexible partitioning bias is" ,"OFF"
          ENDIF
       ENDIF
       CALL cp_print_key_finished_output(output_unit,logger,PRINT,"",error=error)
    ENDIF
    CALL timestop(handle)
  
  END SUBROUTINE fp_env_write

END MODULE restraint_types
