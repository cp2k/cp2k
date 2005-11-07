!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/pair_potential_types [1.0] *
!!
!!   NAME
!!     pair_potential_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     Teodoro Laino [Teo] 11.2005 : Reorganizing the structures to optimize
!!                                   memory management
!!
!!   SOURCE
!******************************************************************************

MODULE pair_potential_types

  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE splines,                         ONLY: spline_data_type, &
                                             spline_data_release,&
                                             spline_data_copy
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pair_potential_types'

  PRIVATE
  PUBLIC :: pair_potential_type,&           
            pair_potential_reallocate

  PUBLIC :: pair_potential_single_release,&
            pair_potential_single_create,&
            pair_potential_single_copy,&
            pair_potential_single_clean,&
            pair_potential_single_type

  PUBLIC :: pair_potential_pp_create,&
            pair_potential_pp_release,&
            pair_potential_pp_type

  PUBLIC :: pair_potential_p_type,&
            pair_potential_p_release

  PUBLIC :: lj_pot_type,&
            ft_pot_type,&
            williams_pot_type,&
            goodwin_pot_type,&
            ipbv_pot_type,&
            eam_pot_type

  PUBLIC :: pair_potential_lj_create

  TYPE ipbv_pot_type
     REAL (KIND=dp), DIMENSION (2:15) :: a
     REAL (KIND=dp) :: rcore
     REAL (KIND=dp) :: m
     REAL (KIND=dp) :: b
  END TYPE ipbv_pot_type

  TYPE lj_pot_type 
     REAL (KIND=dp) :: epsilon
     REAL (KIND=dp) :: sigma6
     REAL (KIND=dp) :: sigma12
  END TYPE Lj_pot_type

  TYPE ft_pot_type 
     REAL (KIND=dp) :: A
     REAL (KIND=dp) :: B
     REAL (KIND=dp) :: C
     REAL (KIND=dp) :: D
  END TYPE ft_pot_type

  TYPE williams_pot_type
     REAL (KIND=dp) :: a
     REAL (KIND=dp) :: b
     REAL (KIND=dp) :: c
  END TYPE williams_pot_type

  TYPE goodwin_pot_type
     REAL (KIND=dp) :: vr0
     REAL (KIND=dp) :: m, mc
     REAL (KIND=dp) :: d, dc
  END TYPE goodwin_pot_type

  TYPE eam_pot_type
     CHARACTER ( len = default_string_length ) :: eam_file_name
     REAL ( KIND = dp ) drar, drhoar, acutal
     REAL ( KIND = dp ), DIMENSION ( 5000 ) :: rho, phi, frho, rhoval, rval
     REAL ( KIND = dp ), DIMENSION ( 5000 ) :: rhop, phip, frhop
  END TYPE eam_pot_type

  TYPE pair_potential_single_type
     REAL (KIND=dp) :: rcutsq
     REAL (KIND=dp) :: energy_cutoff
     REAL (KIND=dp) :: e_cutoff_coul
     CHARACTER ( LEN = default_string_length ) :: TYPE
     CHARACTER ( LEN = default_string_length ) :: at1
     CHARACTER ( LEN = default_string_length ) :: at2
     TYPE ( ipbv_pot_type ),     POINTER :: ipbv
     TYPE ( lj_pot_type ),       POINTER :: lj
     TYPE ( ft_pot_type ),       POINTER :: ft
     TYPE ( williams_pot_type ), POINTER :: willis
     TYPE ( goodwin_pot_type ),  POINTER :: goodwin
     TYPE ( eam_pot_type ),      POINTER :: eam
     TYPE ( spline_data_type ),  POINTER :: pair_spline_data
     TYPE ( spline_data_type ),  POINTER :: rho_spline_data     ! for EAM type potentials
  END TYPE pair_potential_single_type

  TYPE pair_potential_type
     TYPE(pair_potential_single_type),          POINTER :: pot
  END TYPE pair_potential_type

  TYPE pair_potential_p_type
     TYPE(pair_potential_type), DIMENSION(:),   POINTER :: pot
  END TYPE pair_potential_p_type

  TYPE pair_potential_pp_type
     TYPE(pair_potential_type), DIMENSION(:,:), POINTER :: pot
  END TYPE pair_potential_pp_type

CONTAINS

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_single_create [1.0] *
  !!
  !!   NAME
  !!     pair_potential_single_create
  !!
  !!   FUNCTION
  !!     Creates the potential parameter type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_single_create(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(pair_potential_single_type), POINTER    :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN='pair_potential_single_create',&
         routineP=moduleN//':'//routineN
    INTEGER :: i, stat
    LOGICAL :: failure 

    failure = .FALSE.
    CPPostcondition(.NOT.ASSOCIATED(potparm),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       ALLOCATE(potparm, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY( potparm%ipbv,&
                potparm%lj,&
                potparm%ft,&
                potparm%willis,&
                potparm%goodwin,&
                potparm%eam,&
                potparm%pair_spline_data,&
                potparm%rho_spline_data )
       CALL pair_potential_single_clean(potparm, error)
    END IF
  END SUBROUTINE pair_potential_single_create

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_single_clean [1.0] *
  !!
  !!   NAME
  !!     pair_potential_single_clean
  !!
  !!   FUNCTION
  !!     Cleans the potential parameter type
  !!
  !!   AUTHOR
  !!     unknown
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_single_clean(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(pair_potential_single_type), POINTER    :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN='pair_potential_single_clean',&
         routineP=moduleN//':'//routineN
    INTEGER :: i

    potparm % TYPE          = 'NULL'
    potparm % at1           = 'NULL'
    potparm % at2           = 'NULL'
    potparm % rcutsq        = 0.0_dp
    potparm % energy_cutoff = 0.0_dp
    potparm % e_cutoff_coul = 0.0_dp
    IF (ASSOCIATED( potparm % pair_spline_data ))&
         CALL spline_data_release ( potparm % pair_spline_data )
    IF (ASSOCIATED( potparm % rho_spline_data )) &
         CALL spline_data_release ( potparm % rho_spline_data )

    CALL pair_potential_lj_clean (potparm, error=error)
    CALL pair_potential_williams_clean (potparm, error=error)
    CALL pair_potential_goodwin_clean (potparm, error=error)
    CALL pair_potential_eam_clean (potparm, error=error)
    CALL pair_potential_bmhft_clean (potparm, error=error)
    CALL pair_potential_ipbv_clean (potparm, error=error)
    
  END SUBROUTINE pair_potential_single_clean

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_single_copy [1.0] *
  !!
  !!   NAME
  !!     pair_potential_single_copy
  !!
  !!   FUNCTION
  !!     Copy two potential parameter type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_single_copy(potparm_source, potparm_dest, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(pair_potential_single_type), POINTER    :: potparm_source, potparm_dest
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN='pair_potential_single_copy',&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm_dest)) THEN
          CALL pair_potential_single_create(potparm_dest,error)
       END IF
       potparm_dest % TYPE          = potparm_source % TYPE
       potparm_dest % at1           = potparm_source % at1
       potparm_dest % at2           = potparm_source % at2
       potparm_dest % rcutsq        = potparm_source % rcutsq       
       potparm_dest % energy_cutoff = potparm_source % energy_cutoff
       potparm_dest % e_cutoff_coul = potparm_source % e_cutoff_coul
       
       IF (ASSOCIATED( potparm_source % pair_spline_data )) THEN
          CALL spline_data_copy ( potparm_source % pair_spline_data, potparm_dest % pair_spline_data, error=error  )
       END IF
       IF (ASSOCIATED( potparm_source % rho_spline_data )) THEN
          CALL spline_data_copy ( potparm_source % rho_spline_data, potparm_dest % rho_spline_data, error=error  )
       END IF
       
       CALL pair_potential_lj_copy (potparm_source, potparm_dest, error=error)
       CALL pair_potential_williams_copy (potparm_source, potparm_dest, error=error)
       CALL pair_potential_goodwin_copy (potparm_source, potparm_dest, error=error)
       CALL pair_potential_eam_copy (potparm_source, potparm_dest, error=error)
       CALL pair_potential_bmhft_copy (potparm_source, potparm_dest, error=error)
       CALL pair_potential_ipbv_copy (potparm_source, potparm_dest, error=error)    
    END IF
  END SUBROUTINE pair_potential_single_copy

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_copy [1.0] *
  !!
  !!   NAME
  !!     pair_potential_copy
  !!
  !!   FUNCTION
  !!     Copy two potential parameter type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_copy(potparm_source, potparm_dest, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(pair_potential_type), POINTER           :: potparm_source, potparm_dest
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN='pair_potential_copy',&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure

    failure = .FALSE.
    
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(potparm_dest),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm_dest%pot)) CALL pair_potential_single_create(potparm_dest%pot, error)
       CALL pair_potential_single_copy(potparm_source%pot, potparm_dest%pot, error=error)
    END IF
  END SUBROUTINE pair_potential_copy

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_single_release [1.0] *
  !!
  !!   NAME
  !!     pair_potential_single_release
  !!
  !!   FUNCTION
  !!     Release Data-structure that constains potential parameters of a single pair
  !!
  !!   AUTHOR
  !!     Teodoro Laino [Teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!     
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_single_release ( potparm, error )
    IMPLICIT NONE
    ! Arguments
    TYPE(pair_potential_single_type), POINTER :: potparm 
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN='pair_potential_single_release',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat, i, j
    LOGICAL :: failure 

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(ASSOCIATED(potparm),cp_failure_level,routineP,error,failure)
       CALL spline_data_release ( potparm% pair_spline_data ,error=error)
       CALL spline_data_release ( potparm% rho_spline_data, error=error )
       CALL pair_potential_ipbv_release     ( potparm, error=error )
       CALL pair_potential_lj_release       ( potparm, error=error )
       CALL pair_potential_bmhft_release    ( potparm, error=error )
       CALL pair_potential_williams_release ( potparm, error=error )
       CALL pair_potential_goodwin_release  ( potparm, error=error )
       CALL pair_potential_eam_release      ( potparm, error=error )
       DEALLOCATE(potparm, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    NULLIFY ( potparm )
  END SUBROUTINE pair_potential_single_release

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_pp_create [1.0] *
  !!
  !!   NAME
  !!     pair_potential_pp_create
  !!
  !!   FUNCTION
  !!     Data-structure that constains potential parameters
  !!
  !!   AUTHOR
  !!     unknown
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_pp_create ( potparm, nkinds, error )
    IMPLICIT NONE
    ! Arguments
    TYPE(pair_potential_pp_type), POINTER :: potparm 
    INTEGER, INTENT ( IN ) :: nkinds
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pair_potential_pp_create',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat, i, j

    failure=.FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(potparm),cp_failure_level,routineP,error,failure)
       ALLOCATE ( potparm, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( potparm%pot( nkinds, nkinds ), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, nkinds
          DO j = 1, nkinds             
             NULLIFY ( potparm%pot(i,j)%pot )
          END DO
       END DO
       ! Use no-redundancy in the potential definition
       DO i = 1, nkinds
          DO j = i, nkinds
             CALL pair_potential_single_create(potparm%pot( i, j)%pot, error=error)
             potparm%pot( j, i )%pot => potparm%pot( i, j )%pot
          END DO
       END DO
    END IF
  END SUBROUTINE pair_potential_pp_create

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_pp_release [1.0] *
  !!
  !!   NAME
  !!     pair_potential_pp_release
  !!
  !!   FUNCTION
  !!     Release Data-structure that constains potential parameters
  !!
  !!   AUTHOR
  !!     unknown
  !!
  !!   MODIFICATION HISTORY
  !!     Teodoro Laino [Teo] 11.2005 : Reorganizing the structures to optimize 
  !!                                   memory management
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_pp_release ( potparm, error )
    IMPLICIT NONE
    ! Arguments
    TYPE(pair_potential_pp_type), POINTER :: potparm 
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN='pair_potential_pp_release',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat, i, j
    LOGICAL :: failure

    failure = .FALSE.

    IF ( ASSOCIATED ( potparm ) ) THEN
       IF (ASSOCIATED (potparm%pot)) THEN
          DO i = 1, SIZE ( potparm%pot, 1 )
             DO j = i, SIZE ( potparm%pot, 2 )
                CALL pair_potential_single_release(potparm%pot( i, j )%pot, error)
                NULLIFY(potparm%pot( j, i )%pot)
             END DO
          END DO
          DEALLOCATE(potparm%pot, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE ( potparm, stat = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    NULLIFY ( potparm )
  END SUBROUTINE pair_potential_pp_release

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_p_create [1.0] *
  !!
  !!   NAME
  !!     pair_potential_p_create
  !!
  !!   FUNCTION
  !!     Data-structure that constains potential parameters
  !!
  !!   AUTHOR
  !!     unknown
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_p_create ( potparm, ndim, ub, lb, error )
    IMPLICIT NONE
    ! Arguments
    TYPE(pair_potential_p_type), POINTER :: potparm 
    INTEGER, INTENT ( IN ), OPTIONAL     :: ndim, lb, ub
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pair_potential_p_create',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat, i, loc_lb, loc_ub

    failure=.FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(potparm),cp_failure_level,routineP,error,failure)
       ALLOCATE ( potparm, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (PRESENT(ndim)) THEN
          loc_lb = 1
          loc_ub = ndim
          ALLOCATE ( potparm%pot( loc_lb : loc_ub ), stat=stat)
          IF (PRESENT(lb).OR.PRESENT(ub)) THEN
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
       ELSE IF (PRESENT(lb).AND.PRESENT(ub)) THEN
          loc_lb = lb
          loc_ub = ub
          ALLOCATE ( potparm%pot( loc_lb : loc_ub ), stat=stat)
          IF (PRESENT(ndim)) THEN
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
       ELSE
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = loc_lb, loc_ub
          NULLIFY ( potparm%pot(i)%pot)
          CALL pair_potential_single_create(potparm%pot(i)%pot, error=error)
       END DO
    END IF
  END SUBROUTINE pair_potential_p_create

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_p_release [1.0] *
  !!
  !!   NAME
  !!     pair_potential_p_release
  !!
  !!   FUNCTION
  !!     Release Data-structure that constains potential parameters
  !!
  !!   AUTHOR
  !!     unknown
  !!
  !!   MODIFICATION HISTORY
  !!     Teodoro Laino [Teo] 11.2005 : Reorganizing the structures to optimize 
  !!                                   memory management
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_p_release ( potparm, error )
    IMPLICIT NONE
    ! Arguments
    TYPE(pair_potential_p_type), POINTER :: potparm 
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN='pair_potential_p_release',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat, i
    logical :: failure 

    failure = .FALSE.

    IF ( ASSOCIATED ( potparm ) ) THEN
       IF (ASSOCIATED (potparm%pot)) THEN
          DO i = 1, SIZE ( potparm%pot )
             CALL pair_potential_single_release(potparm%pot( i )%pot, error)
          END DO
          DEALLOCATE(potparm%pot, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE ( potparm, stat = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    NULLIFY ( potparm )
  END SUBROUTINE pair_potential_p_release

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_p_copy [1.0] *
  !!
  !!   NAME
  !!     pair_potential_p_copy
  !!
  !!   FUNCTION
  !!     Copy structures between two pair_potential_p_type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [Teo] 11.2005 
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_p_copy(source, dest, istart, iend, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_p_type), POINTER             :: source, dest
    INTEGER, INTENT(IN), OPTIONAL   :: istart, iend
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_p_copy",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: i, l_start, l_end

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(source),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(dest),cp_failure_level,routineP,error,failure)
    l_start = LBOUND(source%pot,1) 
    l_end   = UBOUND(source%pot,1)
    IF (PRESENT(istart)) l_start = istart
    IF (PRESENT(iend))   l_end   = iend
    IF (.NOT.failure) THEN
       DO i = l_start, l_end
          IF (.NOT.ASSOCIATED(source%pot(i)%pot)) &
               CALL pair_potential_single_create(source%pot(i)%pot, error=error)
          CALL pair_potential_single_copy(source%pot(i)%pot,dest%pot(i)%pot,error=error)
       END DO
    END IF
  END SUBROUTINE pair_potential_p_copy

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_reallocate [1.0] *
  !!
  !!   NAME
  !!     pair_potential_reallocate
  !!
  !!   FUNCTION
  !!     Cleans the potential parameter type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [Teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_reallocate(p,lb1_new,ub1_new,lj,williams,goodwin,eam,bmhft,ipbv,error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_p_type), POINTER             :: p
    INTEGER, INTENT(IN)                               :: lb1_new, ub1_new
    LOGICAL, OPTIONAL, INTENT(IN)                     :: lj,williams,goodwin,eam,bmhft,ipbv
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_reallocate",&
         routineP=moduleN//':'//routineN

    INTEGER                            :: stat, lb1, lb1_old, &
         new_size, old_size, ub1, &
         ub1_old, i, ipot
    TYPE (pair_potential_p_type), POINTER             :: work
    LOGICAL :: failure, llj, lwilliams, lgoodwin, leam, lbmhft, lipbv, check

    NULLIFY(work)
    failure = .FALSE.
    ipot = 0
    llj       = .FALSE. ; IF (PRESENT(lj))        llj       = lj
    lwilliams = .FALSE. ; IF (PRESENT(williams))  lwilliams = williams
    lgoodwin  = .FALSE. ; IF (PRESENT(goodwin))   lgoodwin  = goodwin
    leam      = .FALSE. ; IF (PRESENT(eam))       leam      = eam
    lbmhft    = .FALSE. ; IF (PRESENT(bmhft))     lbmhft    = bmhft
    lipbv     = .FALSE. ; IF (PRESENT(ipbv))      lipbv     = ipbv

    IF (llj) THEN 
       ipot  = 1 ! LENNARD-JONES POTENTIAL
       check = .NOT.(lwilliams.OR.lgoodwin.OR.leam.OR.lbmhft.OR.lipbv)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lwilliams) THEN
       ipot  = 2 ! WILLIAMS POTENTIAL 
       check = .NOT.(llj.OR.lgoodwin.OR.leam.OR.lbmhft.OR.lipbv)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lgoodwin) THEN
       ipot  = 3 ! GOODWIN POTENTIAL
       check = .NOT.(llj.OR.lwilliams.OR.leam.OR.lbmhft.OR.lipbv)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (leam) THEN
       ipot  = 4 ! EAM POTENTIAL  
       check = .NOT.(llj.OR.lgoodwin.OR.lwilliams.OR.lbmhft.OR.lipbv)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lbmhft) THEN
       ipot  = 5 ! TOSI-FUMI POTENTIAL
       check = .NOT.(llj.OR.lgoodwin.OR.leam.OR.lwilliams.OR.lipbv)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lipbv) THEN
       ipot  = 6 ! IPBV POTENTIAL
       check = .NOT.(llj.OR.lgoodwin.OR.leam.OR.lwilliams.OR.lbmhft)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       lb1_old = 0
       ub1_old = 0
       IF (ASSOCIATED(p)) THEN
          lb1_old = LBOUND(p%pot,1)
          ub1_old = UBOUND(p%pot,1)
          CALL pair_potential_p_create(work, lb=lb1_old, ub=ub1_old, error=error)
          CALL pair_potential_p_copy(p, work, error=error)
          CALL pair_potential_p_release(p, error=error)
       END IF
       
       CALL pair_potential_p_create(p, lb=lb1_new, ub=ub1_new, error=error)
       IF (ASSOCIATED(work)) THEN
          CALL pair_potential_p_copy(work, p, istart=lb1_old, iend=ub1_old, error=error)
       END IF
       DO i = ub1_old+1, ub1_new
          p % pot ( i ) % pot % TYPE = 'NULL'
          p % pot ( i ) % pot % at1  = 'NULL'
          p % pot ( i ) % pot % at2  = 'NULL'
          SELECT CASE (ipot)
          CASE (1)
             CALL pair_potential_lj_create(p%pot(i)%pot, error=error)
          CASE (2)
             CALL pair_potential_williams_create(p%pot(i)%pot, error=error)
          CASE (3)
             CALL pair_potential_goodwin_create(p%pot(i)%pot, error=error)
          CASE (4)
             CALL pair_potential_eam_create(p%pot(i)%pot, error=error)
          CASE (5)
             CALL pair_potential_bmhft_create(p%pot(i)%pot, error=error)
          CASE (6)
             CALL pair_potential_ipbv_create(p%pot(i)%pot, error=error)
          END SELECT
          NULLIFY ( p%pot(i)%pot%pair_spline_data )
          NULLIFY ( p%pot(i)%pot%rho_spline_data )
       END DO

       IF (ASSOCIATED(work)) CALL pair_potential_p_release(work, error=error)
    END IF
  END SUBROUTINE pair_potential_reallocate

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_lj_create [1.0] *
  !!
  !!   NAME
  !!     pair_potential_lj_create
  !!
  !!   FUNCTION
  !!     Creates the LJ potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_lj_create(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_lj_create",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(potparm%lj),cp_failure_level,routineP,error,failure)
       ALLOCATE(potparm%lj, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_lj_clean(potparm, error=error)
    END IF
  END SUBROUTINE pair_potential_lj_create

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_lj_copy [1.0] *
  !!
  !!   NAME
  !!     pair_potential_lj_copy
  !!
  !!   FUNCTION
  !!     Copy two LJ potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_lj_copy(potparm_source, potparm_dest, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm_source, potparm_dest
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_lj_copy",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(potparm_dest),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm_source%lj)) RETURN 
       IF (ASSOCIATED(potparm_dest%lj)) CALL pair_potential_lj_release(potparm_dest, error=error)
       CALL pair_potential_lj_create(potparm_dest, error=error)
       potparm_dest% lj % epsilon = potparm_source% lj % epsilon 
       potparm_dest% lj % sigma6  = potparm_source% lj % sigma6 
       potparm_dest% lj % sigma12 = potparm_source% lj % sigma12      
    END IF
  END SUBROUTINE pair_potential_lj_copy

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_lj_clean [1.0] *
  !!
  !!   NAME
  !!     pair_potential_lj_clean
  !!
  !!   FUNCTION
  !!     Creates the LJ potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_lj_clean(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_lj_clean",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm%lj)) RETURN
       potparm% lj % epsilon = 0.0_dp
       potparm% lj % sigma6  = 0.0_dp
       potparm% lj % sigma12 = 0.0_dp
    END IF
  END SUBROUTINE pair_potential_lj_clean

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_lj_release [1.0] *
  !!
  !!   NAME
  !!     pair_potential_lj_release
  !!
  !!   FUNCTION
  !!     Destroys the LJ potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_lj_release(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_lj_release",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(potparm%lj)) THEN
          DEALLOCATE(potparm%lj, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(potparm%lj)
    END IF
  END SUBROUTINE pair_potential_lj_release

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_williams_create [1.0] *
  !!
  !!   NAME
  !!     pair_potential_williams_create
  !!
  !!   FUNCTION
  !!     Creates the WILLIAMS potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_williams_create(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_williams_create",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(potparm%willis),cp_failure_level,routineP,error,failure)
       ALLOCATE(potparm%willis, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_williams_clean(potparm, error=error)
    END IF    
  END SUBROUTINE pair_potential_williams_create

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_williams_copy [1.0] *
  !!
  !!   NAME
  !!     pair_potential_williams_copy
  !!
  !!   FUNCTION
  !!     Copy two WILLIAMS potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_williams_copy(potparm_source, potparm_dest, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm_source, potparm_dest
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_williams_copy",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(potparm_dest),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm_source%willis)) RETURN 
       IF (ASSOCIATED(potparm_dest%willis)) CALL pair_potential_williams_release(potparm_dest, error=error)
       CALL pair_potential_williams_create(potparm_dest, error=error)
       potparm_dest% willis % a = potparm_source% willis % a
       potparm_dest% willis % b = potparm_source% willis % b
       potparm_dest% willis % c = potparm_source% willis % c
    END IF
  END SUBROUTINE pair_potential_williams_copy

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_williams_clean [1.0] *
  !!
  !!   NAME
  !!     pair_potential_williams_clean
  !!
  !!   FUNCTION
  !!     Creates the WILLIAMS potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_williams_clean(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_williams_clean",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(potparm%willis)) RETURN
       potparm% willis % a = 0.0_dp
       potparm% willis % b = 0.0_dp
       potparm% willis % c = 0.0_dp
    END IF    
  END SUBROUTINE pair_potential_williams_clean

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_williams_release [1.0] *
  !!
  !!   NAME
  !!     pair_potential_williams_release
  !!
  !!   FUNCTION
  !!     Destroys the WILLIAMS potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_williams_release(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_williams_release",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(potparm%willis)) THEN
          DEALLOCATE(potparm%willis, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(potparm%willis)
    END IF
  END SUBROUTINE pair_potential_williams_release

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_goodwin_create [1.0] *
  !!
  !!   NAME
  !!     pair_potential_goodwin_create
  !!
  !!   FUNCTION
  !!     Creates the GOODWIN potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_goodwin_create(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_goodwin_create",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(potparm%goodwin),cp_failure_level,routineP,error,failure)
       ALLOCATE(potparm%goodwin, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_goodwin_clean(potparm, error=error)
    END IF        
  END SUBROUTINE pair_potential_goodwin_create

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_goodwin_copy [1.0] *
  !!
  !!   NAME
  !!     pair_potential_goodwin_copy
  !!
  !!   FUNCTION
  !!     Copy two GOODWIN potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_goodwin_copy(potparm_source, potparm_dest, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm_source, potparm_dest
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_goodwin_copy",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(potparm_dest),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm_source%goodwin)) RETURN 
       IF (ASSOCIATED(potparm_dest%goodwin)) CALL pair_potential_goodwin_release(potparm_dest, error=error)
       CALL pair_potential_goodwin_create(potparm_dest, error=error)
       potparm_dest% goodwin % vr0 = potparm_source% goodwin % vr0
       potparm_dest% goodwin % d   = potparm_source% goodwin % d  
       potparm_dest% goodwin % dc  = potparm_source% goodwin % dc 
       potparm_dest% goodwin % m   = potparm_source% goodwin % m  
       potparm_dest% goodwin % mc  = potparm_source% goodwin % mc 
    END IF
  END SUBROUTINE pair_potential_goodwin_copy

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_goodwin_clean [1.0] *
  !!
  !!   NAME
  !!     pair_potential_goodwin_clean
  !!
  !!   FUNCTION
  !!     Creates the GOODWIN potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_goodwin_clean(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_goodwin_clean",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(potparm%goodwin)) RETURN
       potparm% goodwin % vr0 = 0.0_dp
       potparm% goodwin % d   = 0.0_dp
       potparm% goodwin % dc  = 0.0_dp
       potparm% goodwin % m   = 0.0_dp
       potparm% goodwin % mc  = 0.0_dp
    END IF        
  END SUBROUTINE pair_potential_goodwin_clean

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_goodwin_release [1.0] *
  !!
  !!   NAME
  !!     pair_potential_goodwin_release
  !!
  !!   FUNCTION
  !!     Destroys the GOODWIN potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_goodwin_release(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_goodwin_release",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(potparm%goodwin)) THEN
          DEALLOCATE(potparm%goodwin, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(potparm%goodwin)
    END IF
  END SUBROUTINE pair_potential_goodwin_release

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_eam_create [1.0] *
  !!
  !!   NAME
  !!     pair_potential_eam_create
  !!
  !!   FUNCTION
  !!     Creates the EAM potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_eam_create(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_eam_create",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(potparm%eam),cp_failure_level,routineP,error,failure)
       ALLOCATE(potparm%eam, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_eam_clean(potparm, error=error)
    END IF            
  END SUBROUTINE pair_potential_eam_create

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_eam_copy [1.0] *
  !!
  !!   NAME
  !!     pair_potential_eam_copy
  !!
  !!   FUNCTION
  !!     Copy two EAM potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_eam_copy(potparm_source, potparm_dest, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm_source, potparm_dest
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_eam_copy",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(potparm_dest),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm_source%eam)) RETURN 
       IF (ASSOCIATED(potparm_dest%eam)) CALL pair_potential_eam_release(potparm_dest, error=error)
       CALL pair_potential_eam_create(potparm_dest, error=error)
       potparm_dest% eam % eam_file_name = potparm_source% eam % eam_file_name
       potparm_dest% eam % drar          = potparm_source% eam % drar         
       potparm_dest% eam % drhoar        = potparm_source% eam % drhoar       
       potparm_dest% eam % acutal        = potparm_source% eam % acutal
       potparm_dest% eam % rho           = potparm_source% eam % rho    
       potparm_dest% eam % phi           = potparm_source% eam % phi    
       potparm_dest% eam % frho          = potparm_source% eam % frho   
       potparm_dest% eam % rhoval        = potparm_source% eam % rhoval 
       potparm_dest% eam % rval          = potparm_source% eam % rval   
       potparm_dest% eam % rhop          = potparm_source% eam % rhop   
       potparm_dest% eam % phip          = potparm_source% eam % phip   
       potparm_dest% eam % frhop         = potparm_source% eam % frhop  
    END IF
  END SUBROUTINE pair_potential_eam_copy

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_eam_clean [1.0] *
  !!
  !!   NAME
  !!     pair_potential_eam_clean
  !!
  !!   FUNCTION
  !!     Creates the EAM potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_eam_clean(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_eam_clean",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(potparm%eam)) RETURN
       potparm% eam % eam_file_name = 'NULL'
       potparm% eam % drar          = 0.0_dp
       potparm% eam % drhoar        = 0.0_dp
       potparm% eam % acutal        = 0.0_dp
       ! Arrays
       potparm% eam % rho           = 0.0_dp
       potparm% eam % phi           = 0.0_dp
       potparm% eam % frho          = 0.0_dp
       potparm% eam % rhoval        = 0.0_dp
       potparm% eam % rval          = 0.0_dp
       potparm% eam % rhop          = 0.0_dp
       potparm% eam % phip          = 0.0_dp
       potparm% eam % frhop         = 0.0_dp
    END IF            
  END SUBROUTINE pair_potential_eam_clean

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_eam_release [1.0] *
  !!
  !!   NAME
  !!     pair_potential_eam_release
  !!
  !!   FUNCTION
  !!     Destroys the EAM potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_eam_release(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_eam_release",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(potparm%eam)) THEN
          DEALLOCATE(potparm%eam, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(potparm%eam)
    END IF
  END SUBROUTINE pair_potential_eam_release

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_bmhft_create [1.0] *
  !!
  !!   NAME
  !!     pair_potential_bmhft_create
  !!
  !!   FUNCTION
  !!     Creates the BMHFT (TOSI-FUMI) potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_bmhft_create(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_bmhft_create",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(potparm%ft),cp_failure_level,routineP,error,failure)
       ALLOCATE(potparm%ft, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_bmhft_clean(potparm, error=error)
    END IF        
  END SUBROUTINE pair_potential_bmhft_create

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_bmhft_copy [1.0] *
  !!
  !!   NAME
  !!     pair_potential_bmhft_copy
  !!
  !!   FUNCTION
  !!     Copy two BMHFT (TOSI-FUMI) potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_bmhft_copy(potparm_source, potparm_dest, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm_source, potparm_dest
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_bmhft_copy",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(potparm_dest),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm_source%ft)) RETURN 
       IF (ASSOCIATED(potparm_dest%ft)) CALL pair_potential_bmhft_release(potparm_dest, error=error)
       CALL pair_potential_bmhft_create(potparm_dest, error=error)
       potparm_dest% ft %  A     = potparm_source% ft %  A
       potparm_dest% ft %  B     = potparm_source% ft %  B       
       potparm_dest% ft %  C     = potparm_source% ft %  C
       potparm_dest% ft %  D     = potparm_source% ft %  D
    END IF
  END SUBROUTINE pair_potential_bmhft_copy

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_bmhft_clean [1.0] *
  !!
  !!   NAME
  !!     pair_potential_bmhft_clean
  !!
  !!   FUNCTION
  !!     Creates the BMHFT (TOSI-FUMI) potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_bmhft_clean(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_bmhft_clean",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(potparm%ft)) RETURN
       potparm%ft % A = 0.0_dp 
       potparm%ft % B = 0.0_dp 
       potparm%ft % C = 0.0_dp 
       potparm%ft % D = 0.0_dp 
    END IF        
  END SUBROUTINE pair_potential_bmhft_clean

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_bmhft_release [1.0] *
  !!
  !!   NAME
  !!     pair_potential_bmhft_release
  !!
  !!   FUNCTION
  !!     Destroys the BMHFT potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_bmhft_release(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_bmhft_release",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(potparm%ft)) THEN
          DEALLOCATE(potparm%ft, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(potparm%ft)
    END IF
  END SUBROUTINE pair_potential_bmhft_release

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_ipbv_create [1.0] *
  !!
  !!   NAME
  !!     pair_potential_ipbv_create
  !!
  !!   FUNCTION
  !!     Creates the IPBV potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_ipbv_create(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_ipbv_create",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(potparm%ipbv),cp_failure_level,routineP,error,failure)
       ALLOCATE(potparm%ipbv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_ipbv_clean(potparm, error=error)
    END IF            
  END SUBROUTINE pair_potential_ipbv_create

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_ipbv_copy [1.0] *
  !!
  !!   NAME
  !!     pair_potential_ipbv_copy
  !!
  !!   FUNCTION
  !!     Copy two IPBV potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_ipbv_copy(potparm_source, potparm_dest, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm_source, potparm_dest
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_ipbv_copy",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(potparm_dest),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm_source%ipbv)) RETURN 
       IF (ASSOCIATED(potparm_dest%ipbv)) CALL pair_potential_ipbv_release(potparm_dest, error=error)
       CALL pair_potential_ipbv_create(potparm_dest, error=error)
       potparm_dest% ipbv %  a        = potparm_source% ipbv % a
       potparm_dest% ipbv %  rcore    = potparm_source% ipbv % rcore
       potparm_dest% ipbv %  b        = potparm_source% ipbv % b
       potparm_dest% ipbv %  m        = potparm_source% ipbv % m
    END IF
  END SUBROUTINE pair_potential_ipbv_copy

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_ipbv_clean [1.0] *
  !!
  !!   NAME
  !!     pair_potential_ipbv_clean
  !!
  !!   FUNCTION
  !!     Creates the IPBV potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_ipbv_clean(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_ipbv_clean",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(potparm%ipbv)) RETURN
       potparm%ipbv % a        = 0.0_dp
       potparm%ipbv % rcore    = 0.0_dp
       potparm%ipbv % b        = 0.0_dp
       potparm%ipbv % m        = 0.0_dp
    END IF            
  END SUBROUTINE pair_potential_ipbv_clean

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_ipbv_release [1.0] *
  !!
  !!   NAME
  !!     pair_potential_ipbv_release
  !!
  !!   FUNCTION
  !!     Destroys the IPBV potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_ipbv_release(potparm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE (pair_potential_single_type), POINTER :: potparm
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER        :: routineN = "pair_potential_ipbv_release",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(potparm%ipbv)) THEN
          DEALLOCATE(potparm%ipbv, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(potparm%ipbv)
    END IF
  END SUBROUTINE pair_potential_ipbv_release

END MODULE pair_potential_types


