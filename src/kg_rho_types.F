!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/kg_rho_types [1.0] *
!!
!!   NAME
!!     kg_rho_types
!!
!!   FUNCTION
!!     rho_types    
!!
!!   AUTHOR
!!     gloria,30.09.2002, previously in kg_force 
!!
!!   MODIFICATION HISTORY
!!     JGH (22-Feb-03) PW grid options added
!!     gt 16-nov-03 moved initialization in this new module
!!
!!   SOURCE
!******************************************************************************

MODULE kg_rho_types
  USE coefficient_types,               ONLY: coeff_allocate,&
                                             coeff_deallocate,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cube_utils,                      ONLY: cube_info_type,&
                                             destroy_cube_info
  USE dft_types,                       ONLY: dft_control_type
  USE gaussian_gridlevels,             ONLY: destroy_gaussian_gridlevel,&
                                             gridlevel_info_type  
  USE kinds,                           ONLY: dbl,dp_size
  USE l_utils,                         ONLY: destroy_l_info,&
                                             l_info_type
  USE mathconstants,                   ONLY: zero 
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_max,&
                                             mp_range,&
                                             mp_sum,&
                                             mp_sync
  USE pw_grid_types,                   ONLY: pw_grid_type  
  USE pw_pool_types,                   ONLY: pw_pool_release,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='kg_rho_types'
  PUBLIC :: kg_rho_type, kg_rho_release, &
            kg_rho_retain, allocate_kg_densities, deallocate_kg_densities

! definition of all variables (densities, potentials, grids....)
! needed to perform a force calculation 
  TYPE kg_rho_type
     ! definition of all grid types used (working grids and multigrids)
     TYPE(pw_grid_type), DIMENSION(:), POINTER :: pw_grid
     ! number of the grid used as reference grid for the density
     INTEGER :: auxbas_grid
     INTEGER :: ref_count
     LOGICAL :: kgpol       
     REAL(dbl) :: total_rho_gspace, total_rho_core_rspace, &
                  total_rho_rspace  
     TYPE(coeff_type)                             :: rho_gspace, &
                                                     rho_p_gspace, &
                                                     v_gspace, &
                                                     v_rspace, &
                                                     vp_rspace, &
                                                     v_xc_rspace, &
                                                     v1_xc_rspace, &
                                                     v2_xc_rspace
     ! multigrids in qs-style
     TYPE(coeff_type), DIMENSION(:), POINTER        :: mgrid_gspace, &
                                                      mgrid_rspace
     ! rho_elec_rspace contains the frozen electron density (for cube output)
     TYPE(coeff_type)                :: rho_elec_rspace
     ! rho_pol_rspace contains the polarization density (for cube output)
     TYPE(coeff_type)                :: rho_p_rspace
     ! hartree contains the total hartree potential (for cube ouput)
     TYPE(coeff_type)          :: v_hartree_rspace
     ! the core density (in KG, density of the Zeff nuclear charges)
     TYPE(coeff_type)          :: rho_core
     TYPE(gridlevel_info_type),POINTER :: gridlevel_info
     TYPE(l_info_type),        POINTER         :: l_info
     TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
     TYPE(pw_pool_type), POINTER :: pw_pool
  END TYPE kg_rho_type
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
!!***
!******************************************************************************

  SUBROUTINE allocate_kg_densities(kg_rho) 

! Arguments
    TYPE(kg_rho_type), POINTER  :: kg_rho

! Locals
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_kg_densities (MODULE kg_rho_types)"
    TYPE(pw_grid_type), POINTER, DIMENSION (:) :: pw_grid 
    INTEGER :: auxbas, igrid_level, istat, ngrid_level
    INTEGER :: handle

!---------------------------------------------------------------------------

    CALL timeset(routine,"I","",handle)

    NULLIFY(kg_rho%mgrid_gspace, kg_rho%mgrid_rspace)    

    auxbas = kg_rho % auxbas_grid
    pw_grid => kg_rho % pw_grid
    ngrid_level = SIZE ( pw_grid )

! Allocate the rho's on the largest reference grid

    CALL coeff_allocate(kg_rho%rho_elec_rspace,pw_grid(auxbas),REALDATA3D)
    kg_rho%rho_elec_rspace%pw%in_space = REALSPACE
    CALL coeff_zero (kg_rho%rho_elec_rspace)

    CALL coeff_allocate(kg_rho%rho_gspace,pw_grid(auxbas),COMPLEXDATA1D)
    kg_rho%rho_gspace%pw%in_space = RECIPROCALSPACE
    CALL coeff_zero (kg_rho%rho_gspace)

    CALL coeff_allocate(kg_rho%v_xc_rspace,pw_grid(auxbas),REALDATA3D)
    kg_rho%v_xc_rspace%pw%in_space = REALSPACE
    CALL coeff_zero (kg_rho%v_xc_rspace)

    CALL coeff_allocate(kg_rho%v_rspace,pw_grid(auxbas),REALDATA3D)
    kg_rho%v_rspace%pw%in_space = REALSPACE
    CALL coeff_zero (kg_rho%v_rspace)

    CALL coeff_allocate(kg_rho%v_gspace,pw_grid(auxbas),COMPLEXDATA1D)
    kg_rho%v_gspace%pw%in_space = RECIPROCALSPACE
    CALL coeff_zero (kg_rho%v_gspace)

    CALL coeff_allocate(kg_rho%v_hartree_rspace,pw_grid(auxbas),REALDATA3D)
    kg_rho%v_hartree_rspace%pw%in_space = REALSPACE
    CALL coeff_zero (kg_rho%v_hartree_rspace)

    CALL coeff_allocate(kg_rho%rho_core,pw_grid(auxbas),COMPLEXDATA1D)
    kg_rho%rho_core%pw%in_space = RECIPROCALSPACE
    CALL coeff_zero (kg_rho%rho_core)

    IF ( kg_rho%kgpol ) THEN

      CALL coeff_allocate(kg_rho%rho_p_gspace,pw_grid(auxbas),COMPLEXDATA1D)
      kg_rho%rho_p_gspace%pw%in_space = RECIPROCALSPACE
      CALL coeff_zero (kg_rho%rho_p_gspace)

      CALL coeff_allocate(kg_rho%rho_p_rspace,pw_grid(auxbas),REALDATA3D)
      kg_rho%rho_p_rspace%pw%in_space = REALSPACE
      CALL coeff_zero (kg_rho%rho_p_rspace)

      CALL coeff_allocate(kg_rho%v1_xc_rspace,pw_grid(auxbas),REALDATA3D)
      kg_rho%v1_xc_rspace%pw%in_space = REALSPACE
      CALL coeff_zero (kg_rho%v1_xc_rspace)

      CALL coeff_allocate(kg_rho%v2_xc_rspace,pw_grid(auxbas),REALDATA3D)
      kg_rho%v2_xc_rspace%pw%in_space = REALSPACE
      CALL coeff_zero (kg_rho%v2_xc_rspace)

      CALL coeff_allocate(kg_rho%vp_rspace,pw_grid(auxbas),REALDATA3D)
      kg_rho%vp_rspace%pw%in_space = REALSPACE
      CALL coeff_zero (kg_rho%vp_rspace)

    END IF

! allocate the multi-grid rho's
    ALLOCATE (kg_rho%mgrid_rspace(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
    ALLOCATE (kg_rho%mgrid_gspace(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)

! setup the multigrids
    DO igrid_level=1,ngrid_level
       CALL coeff_allocate(kg_rho%mgrid_rspace(igrid_level),pw_grid(igrid_level),&
            REALDATA3D)
       kg_rho%mgrid_rspace(igrid_level)%pw%in_space = REALSPACE
       CALL coeff_zero (kg_rho % mgrid_rspace(igrid_level))
       CALL coeff_allocate(kg_rho%mgrid_gspace(igrid_level),pw_grid(igrid_level),&
            COMPLEXDATA1D)
       kg_rho%mgrid_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE
       CALL coeff_zero (kg_rho % mgrid_gspace(igrid_level))
    END DO

    CALL timestop(0.0_dbl,handle)

END SUBROUTINE allocate_kg_densities

! *************************************************************************
  SUBROUTINE deallocate_kg_densities(kg_rho)

! Arguments
    TYPE(kg_rho_type), POINTER     :: kg_rho

! Locals
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_kg_densities (MODULE kg_rho)"
    TYPE(pw_grid_type), POINTER, DIMENSION (:) :: pw_grid 
    INTEGER :: auxbas, igrid_level, istat, ngrid_level
    INTEGER :: handle

!---------------------------------------------------------------------------

    CALL timeset("deallocate_density","I","",handle)

    auxbas = kg_rho % auxbas_grid
    pw_grid => kg_rho % pw_grid
    ngrid_level = SIZE ( pw_grid )

! DeAllocate the rho's on the largest reference grid

    CALL coeff_deallocate(kg_rho%rho_elec_rspace)

    CALL coeff_deallocate(kg_rho%rho_gspace)

    CALL coeff_deallocate(kg_rho%v_xc_rspace)

    CALL coeff_deallocate(kg_rho%v_rspace)

    CALL coeff_deallocate(kg_rho%v_gspace)

    CALL coeff_deallocate(kg_rho%v_hartree_rspace)

    CALL coeff_deallocate(kg_rho%rho_core)

    IF (kg_rho % kgpol) THEN

      CALL coeff_deallocate(kg_rho%rho_p_rspace)

      CALL coeff_deallocate(kg_rho%rho_p_gspace)

      CALL coeff_deallocate(kg_rho%v1_xc_rspace)

      CALL coeff_deallocate(kg_rho%v2_xc_rspace)

      CALL coeff_deallocate(kg_rho%vp_rspace)

   END IF


! deallocate the multi-grid rho's

    DO igrid_level=1,ngrid_level
       CALL coeff_deallocate(kg_rho%mgrid_rspace(igrid_level))
       CALL coeff_deallocate(kg_rho%mgrid_gspace(igrid_level))
    END DO

    DEALLOCATE (kg_rho%mgrid_rspace, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
    DEALLOCATE (kg_rho%mgrid_gspace,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)

    CALL timestop(0.0_dbl,handle)

END SUBROUTINE deallocate_kg_densities
!!*** **********************************************************************
SUBROUTINE kg_rho_retain(kg_rho,error)
  TYPE(kg_rho_type), POINTER :: kg_rho
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
 
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='kg_rho_retain',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
 
  CPPrecondition(ASSOCIATED(kg_rho),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(kg_rho%ref_count>0,cp_failure_level,routineP,error,failure)
     kg_rho%ref_count=kg_rho%ref_count+1
  END IF
END SUBROUTINE kg_rho_retain
!***************************************************************************

!!****f* kg_force/release_kg_rho [1.0] *
!!
!!   NAME
!!     release_kg_rho
!!
!!   FUNCTION
!!     releases the memory used by the kg_rho and kg_rho 
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE kg_rho_release(kg_rho,error)               

  IMPLICIT NONE

!arguments
    TYPE(kg_rho_type), POINTER       :: kg_rho
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)            :: error

!Locals
    CHARACTER(len=*), PARAMETER :: routineN='kg_rho_release',&
       routineP=moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: igrid_level,stat,handle

!   ---------------------------------------------------------------------------
 
   CALL timeset(routineN,"I","",handle)
   failure=.FALSE.  
   IF (ASSOCIATED(kg_rho)) THEN
     CPPrecondition(kg_rho%ref_count>0,cp_failure_level,routineP,error,failure)
     kg_rho%ref_count=kg_rho%ref_count-1
     IF (kg_rho%ref_count<1) THEN
       IF (ASSOCIATED(kg_rho%gridlevel_info)) THEN
         CALL destroy_gaussian_gridlevel(kg_rho%gridlevel_info,error=error)
         DEALLOCATE(kg_rho%gridlevel_info,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(kg_rho%l_info)) THEN
         CALL destroy_l_info(kg_rho%l_info)
         DEALLOCATE(kg_rho%l_info, stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(kg_rho%cube_info)) THEN
         DO igrid_level=1,SIZE(kg_rho%cube_info)
           CALL destroy_cube_info(kg_rho%cube_info(igrid_level))
         END DO
        DEALLOCATE (kg_rho%cube_info,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF   
      CALL pw_pool_release (kg_rho%pw_pool)
      CALL deallocate_kg_densities(kg_rho)
    END IF
  END IF
  NULLIFY(kg_rho)
  CALL timestop(0.0_dbl,handle)
  END SUBROUTINE kg_rho_release

!***************************************************************************

END MODULE kg_rho_types

!******************************************************************************
