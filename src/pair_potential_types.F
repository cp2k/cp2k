!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
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

  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE splines,                         ONLY: spline_data_copy,&
                                             spline_data_release,&
                                             spline_data_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pair_potential_types'

  PRIVATE
  ! when adding a new nonbonedd potential please update also the list_pot
  ! used for the linear scaling screening of potential calculation
  INTEGER, PUBLIC, PARAMETER :: nn_type= 0,&
                                lj_type= 1,&
                                lj_charmm_type= 2,&
                                ft_type= 3,&
                                wl_type= 4,&
                                gw_type= 5,&
                                ip_type= 6,&
                                ea_type= 7,&
                                b4_type= 8,&
                                bm_type= 9,&
                                gp_type=10

  INTEGER, PUBLIC, PARAMETER, DIMENSION(11) :: list_pot = (/nn_type,&
                                                           lj_type,&
                                                           lj_charmm_type,&
                                                           ft_type,&
                                                           wl_type,&
                                                           gw_type,&
                                                           ip_type,&
                                                           ea_type,&
                                                           b4_type,&
                                                           bm_type,&
                                                           gp_type/)

  INTEGER, PUBLIC, PARAMETER :: nosh_nosh=0,&
                                nosh_sh=1,&
                                sh_sh=2
  INTEGER, PUBLIC, PARAMETER, DIMENSION(3) :: list_sh_type =(/nosh_nosh,nosh_sh,sh_sh/)

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
            eam_pot_type,&
            buck4ran_pot_type,&
            buckmorse_pot_type,&
            gp_pot_type

  PUBLIC :: pair_potential_lj_create
  PUBLIC :: compare_pot

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
     REAL ( KIND = dp )                     :: drar, drhoar, acutal
     REAL ( KIND = dp ), DIMENSION ( 5000 ) :: rho, phi, frho, rhoval, rval
     REAL ( KIND = dp ), DIMENSION ( 5000 ) :: rhop, phip, frhop
  END TYPE eam_pot_type

  TYPE buck4ran_pot_type
     REAL (KIND=dp) :: a
     REAL (KIND=dp) :: b
     REAL (KIND=dp) :: c
     REAL (KIND=dp) :: r1
     REAL (KIND=dp) :: r2
     REAL (KIND=dp) :: r3
     INTEGER :: npoly1, npoly2
     REAL (KIND=dp), DIMENSION(0:10) :: poly1
     REAL (KIND=dp), DIMENSION(0:10) :: poly2
  END TYPE buck4ran_pot_type

  TYPE buckmorse_pot_type
     REAL (KIND=dp) :: f0
     REAL (KIND=dp) :: a1
     REAL (KIND=dp) :: a2
     REAL (KIND=dp) :: b1
     REAL (KIND=dp) :: b2
     REAL (KIND=dp) :: c
     REAL (KIND=dp) :: d
     REAL (KIND=dp) :: r0
     REAL (KIND=dp) :: beta
  END TYPE buckmorse_pot_type

  TYPE gp_pot_type
     INTEGER                                   :: myid
     CHARACTER ( LEN = default_path_length )   :: potential
     CHARACTER ( LEN = default_string_length ), POINTER, DIMENSION(:) :: parameters
     CHARACTER ( LEN = default_string_length ) :: variables
     REAL(KIND=dp), DIMENSION(:), POINTER      :: values
  END TYPE gp_pot_type

  TYPE pair_potential_single_type
     REAL (KIND=dp) :: rcutsq
     REAL (KIND=dp) :: energy_cutoff
     REAL (KIND=dp) :: e_cutoff_coul
     LOGICAL        :: undef
     INTEGER        :: TYPE
     INTEGER        :: shell_type
     CHARACTER ( LEN = default_string_length ) :: at1
     CHARACTER ( LEN = default_string_length ) :: at2
     TYPE ( ipbv_pot_type ),       POINTER :: ipbv
     TYPE ( gp_pot_type),          POINTER :: gp
     TYPE ( lj_pot_type ),         POINTER :: lj
     TYPE ( ft_pot_type ),         POINTER :: ft
     TYPE ( williams_pot_type ),   POINTER :: willis
     TYPE ( goodwin_pot_type ),    POINTER :: goodwin
     TYPE ( eam_pot_type ),        POINTER :: eam
     TYPE ( buck4ran_pot_type ),   POINTER :: buck4r
     TYPE ( buckmorse_pot_type ),  POINTER :: buckmo
     TYPE ( spline_data_type ),    POINTER :: pair_spline_data
     TYPE ( spline_data_type ),    POINTER :: rho_spline_data     ! for EAM type potentials
     TYPE ( spline_data_type ),    POINTER :: c_coul_spline_data      ! for shell-model (core-ion or core-core)
     TYPE ( spline_data_type ),    POINTER :: cs_coul_spline_data     ! for shell-model (core-shell) 
     TYPE ( spline_data_type ),    POINTER :: sc_coul_spline_data     ! for shell-model (shell-core) 
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
  !!****s* pair_potential_types/compare_pot [1.0] *
  !!
  !!   NAME
  !!     compare_pot
  !!
  !!   FUNCTION
  !!     compare two different potentials
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 05.2006
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE compare_pot(pot1, pot2, compare, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot1, pot2
    LOGICAL, INTENT(OUT)                     :: compare
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compare_pot', &
      routineP = moduleN//':'//routineN

    compare=.FALSE.
    IF (pot1%type/=pot2%type) RETURN

    SELECT CASE(pot1%type)
    CASE (lj_type,lj_charmm_type)
       IF  ((pot1%lj%epsilon==pot2%lj%epsilon).AND.&
            (pot1%lj%sigma6 ==pot2%lj%sigma6 ).AND.&
            (pot1%lj%sigma12==pot2%lj%sigma12)) compare=.TRUE.
    CASE (wl_type)
       IF  ((pot1%willis%a==pot2%willis%a).AND.&
            (pot1%willis%b==pot2%willis%b).AND.&
            (pot1%willis%c==pot2%willis%c)) compare=.TRUE.
    CASE (gw_type)
       IF  ((pot1%goodwin%vr0==pot2%goodwin%vr0).AND.&
            (pot1%goodwin%m  ==pot2%goodwin%m  ).AND.&
            (pot1%goodwin%mc ==pot2%goodwin%mc ).AND.&
            (pot1%goodwin%d  ==pot2%goodwin%d  ).AND.&
            (pot1%goodwin%dc ==pot2%goodwin%dc )) compare=.TRUE.
    CASE (ea_type)
       IF  ((pot1%eam%drar  ==pot2%eam%drar  ).AND.&
            (pot1%eam%drhoar==pot2%eam%drhoar).AND.&
            (pot1%eam%acutal==pot2%eam%acutal).AND.&
            (SUM(ABS(pot1%eam%rho   -pot2%eam%rho   ))==0.0_dp).AND.&
            (SUM(ABS(pot1%eam%phi   -pot2%eam%phi   ))==0.0_dp).AND.&
            (SUM(ABS(pot1%eam%frho  -pot2%eam%frho  ))==0.0_dp).AND.&
            (SUM(ABS(pot1%eam%rhoval-pot2%eam%rhoval))==0.0_dp).AND.&
            (SUM(ABS(pot1%eam%rval  -pot2%eam%rval  ))==0.0_dp).AND.&
            (SUM(ABS(pot1%eam%rhop  -pot2%eam%rhop  ))==0.0_dp).AND.&
            (SUM(ABS(pot1%eam%phip  -pot2%eam%phip  ))==0.0_dp).AND.&
            (SUM(ABS(pot1%eam%frhop -pot2%eam%frhop ))==0.0_dp)) compare=.TRUE.
    CASE (ft_type)
       IF  ((pot1%ft%A==pot2%ft%A).AND.&
            (pot1%ft%B==pot2%ft%B).AND.&
            (pot1%ft%C==pot2%ft%C).AND.&
            (pot1%ft%D==pot2%ft%D)) compare=.TRUE.
    CASE (ip_type)
       IF  ((SUM(ABS(pot1%ipbv%a-pot2%ipbv%a))==0.0_dp).AND.&
            (pot1%ipbv%rcore==pot2%ipbv%rcore).AND.&
            (pot1%ipbv%m    ==pot2%ipbv%m    ).AND.&
            (pot1%ipbv%b    ==pot2%ipbv%b    )) compare=.TRUE.
    END SELECT

  END SUBROUTINE compare_pot
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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pair_potential_single_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(.NOT.ASSOCIATED(potparm),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       ALLOCATE(potparm, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY( potparm%ipbv,&
                potparm%lj,&
                potparm%gp,&
                potparm%ft,&
                potparm%willis,&
                potparm%goodwin,&
                potparm%eam,&
                potparm%buck4r,&
                potparm%buckmo,&
                potparm%pair_spline_data,&
                potparm%rho_spline_data,&
                potparm%c_coul_spline_data,&
                potparm%cs_coul_spline_data,&
                potparm%sc_coul_spline_data)
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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pair_potential_single_clean', &
      routineP = moduleN//':'//routineN

    potparm % type          = nn_type
    potparm % shell_type    = nosh_nosh
    potparm % undef         = .TRUE.
    potparm % at1           = 'NULL'
    potparm % at2           = 'NULL'
    potparm % rcutsq        = 0.0_dp
    potparm % energy_cutoff = 0.0_dp
    potparm % e_cutoff_coul = 0.0_dp
    IF (ASSOCIATED( potparm % pair_spline_data ))&
         CALL spline_data_release ( potparm % pair_spline_data ,error=error)
    IF (ASSOCIATED( potparm % rho_spline_data )) &
         CALL spline_data_release ( potparm % rho_spline_data ,error=error)
    IF (ASSOCIATED( potparm % c_coul_spline_data )) &
         CALL spline_data_release ( potparm % c_coul_spline_data ,error=error)
    IF (ASSOCIATED( potparm % cs_coul_spline_data )) &
         CALL spline_data_release ( potparm % cs_coul_spline_data ,error=error)
    IF (ASSOCIATED( potparm % sc_coul_spline_data )) &
         CALL spline_data_release ( potparm % sc_coul_spline_data ,error=error)

    CALL pair_potential_lj_clean (potparm, error=error)
    CALL pair_potential_williams_clean (potparm, error=error)
    CALL pair_potential_goodwin_clean (potparm, error=error)
    CALL pair_potential_eam_clean (potparm, error=error)
    CALL pair_potential_buck4r_clean (potparm, error=error)
    CALL pair_potential_buckmo_clean (potparm, error=error)
    CALL pair_potential_bmhft_clean (potparm, error=error)
    CALL pair_potential_ipbv_clean (potparm, error=error)
    CALL pair_potential_gp_clean (potparm, error=error)

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm_source, potparm_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pair_potential_single_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm_dest)) THEN
          CALL pair_potential_single_create(potparm_dest,error)
       END IF
       potparm_dest % TYPE          = potparm_source % TYPE
       potparm_dest % shell_type    = potparm_source % shell_type
       potparm_dest % undef         = potparm_source % undef
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
       IF (ASSOCIATED( potparm_source % c_coul_spline_data )) THEN
          CALL spline_data_copy ( potparm_source % c_coul_spline_data, potparm_dest % c_coul_spline_data, error=error  )
       END IF
       IF (ASSOCIATED( potparm_source % cs_coul_spline_data )) THEN
          CALL spline_data_copy ( potparm_source % cs_coul_spline_data, potparm_dest % cs_coul_spline_data, error=error  )
       END IF
       IF (ASSOCIATED( potparm_source % sc_coul_spline_data )) THEN
          CALL spline_data_copy ( potparm_source % sc_coul_spline_data, potparm_dest % sc_coul_spline_data, error=error  )
       END IF

       CALL pair_potential_lj_copy (potparm_source, potparm_dest, error=error)
       CALL pair_potential_williams_copy (potparm_source, potparm_dest, error=error)
       CALL pair_potential_goodwin_copy (potparm_source, potparm_dest, error=error)
       CALL pair_potential_eam_copy (potparm_source, potparm_dest, error=error)
       CALL pair_potential_bmhft_copy (potparm_source, potparm_dest, error=error)
       CALL pair_potential_ipbv_copy (potparm_source, potparm_dest, error=error)
       CALL pair_potential_buck4r_copy (potparm_source, potparm_dest, error=error)
       CALL pair_potential_buckmo_copy (potparm_source, potparm_dest, error=error)
       CALL pair_potential_gp_copy (potparm_source, potparm_dest, error=error)
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
    TYPE(pair_potential_type), POINTER       :: potparm_source, potparm_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pair_potential_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'pair_potential_single_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(ASSOCIATED(potparm),cp_failure_level,routineP,error,failure)
       CALL spline_data_release ( potparm% pair_spline_data ,error=error)
       CALL spline_data_release ( potparm% rho_spline_data, error=error )
       CALL spline_data_release ( potparm% c_coul_spline_data, error=error )
       CALL spline_data_release ( potparm% cs_coul_spline_data, error=error )
       CALL spline_data_release ( potparm% sc_coul_spline_data, error=error )
       CALL pair_potential_ipbv_release     ( potparm, error=error )
       CALL pair_potential_lj_release       ( potparm, error=error )
       CALL pair_potential_bmhft_release    ( potparm, error=error )
       CALL pair_potential_williams_release ( potparm, error=error )
       CALL pair_potential_goodwin_release  ( potparm, error=error )
       CALL pair_potential_eam_release      ( potparm, error=error )
       CALL pair_potential_buck4r_release   ( potparm, error=error )
       CALL pair_potential_buckmo_release   ( potparm, error=error )
       CALL pair_potential_gp_release       ( potparm, error=error )
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
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    INTEGER, INTENT(IN)                      :: nkinds
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pair_potential_pp_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pair_potential_pp_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_p_type), POINTER     :: potparm
    INTEGER, INTENT(IN), OPTIONAL            :: ndim, ub, lb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pair_potential_p_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, loc_lb, loc_ub, stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_p_type), POINTER     :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pair_potential_p_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_p_type), POINTER     :: source, dest
    INTEGER, INTENT(IN), OPTIONAL            :: istart, iend
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_p_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, l_end, l_start
    LOGICAL                                  :: failure

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
  SUBROUTINE pair_potential_reallocate(p,lb1_new,ub1_new,lj,lj_charmm,williams,goodwin,&
       eam,bmhft,ipbv,buck4r,buckmo,gp,error)
    TYPE(pair_potential_p_type), POINTER     :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new
    LOGICAL, INTENT(IN), OPTIONAL            :: lj, lj_charmm, williams, &
                                                goodwin, eam, bmhft, ipbv, &
                                                buck4r, buckmo, gp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_reallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ipot, lb1_old, ub1_old
    LOGICAL                                  :: check, failure, lbmhft, &
                                                lbuck4r, lbuckmo, leam, &
                                                lgoodwin, lgp, lipbv, llj, &
                                                llj_charmm, lwilliams
    TYPE(pair_potential_p_type), POINTER     :: work

    NULLIFY(work)
    failure = .FALSE.
    ipot = 0
    llj       = .FALSE. ; IF (PRESENT(lj))        llj       = lj
    llj_charmm= .FALSE. ; IF (PRESENT(lj_charmm)) llj_charmm= lj_charmm
    lwilliams = .FALSE. ; IF (PRESENT(williams))  lwilliams = williams
    lgoodwin  = .FALSE. ; IF (PRESENT(goodwin))   lgoodwin  = goodwin
    leam      = .FALSE. ; IF (PRESENT(eam))       leam      = eam
    lbmhft    = .FALSE. ; IF (PRESENT(bmhft))     lbmhft    = bmhft
    lipbv     = .FALSE. ; IF (PRESENT(ipbv))      lipbv     = ipbv
    lbuck4r   = .FALSE. ; IF (PRESENT(buck4r))    lbuck4r   = buck4r
    lbuckmo   = .FALSE. ; IF (PRESENT(buckmo))    lbuckmo   = buckmo
    lgp       = .FALSE. ; IF (PRESENT(gp))        lgp       = gp

    IF (llj) THEN
       ipot  = lj_type
       check = .NOT.(llj_charmm.OR.lwilliams.OR.lgoodwin.OR.leam.OR.lbmhft&
                     .OR.lipbv.OR.lbuck4r.OR.lbuckmo.OR.lgp)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (llj_charmm) THEN
       ipot  = lj_charmm_type
       check = .NOT.(llj.OR.lwilliams.OR.lgoodwin.OR.leam.OR.lbmhft.OR.lipbv&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lwilliams) THEN
       ipot  = wl_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lbmhft.OR.lipbv&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lgoodwin) THEN
       ipot  = gw_type
       check = .NOT.(llj.OR.llj_charmm.OR.lwilliams.OR.leam.OR.lbmhft.OR.lipbv&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (leam) THEN
       ipot  = ea_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.lwilliams.OR.lbmhft.OR.lipbv&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lbmhft) THEN
       ipot  = ft_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lwilliams.OR.lipbv&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lipbv) THEN
       ipot  = ip_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lwilliams.OR.lbmhft&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lbuck4r) THEN
       ipot  = b4_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lwilliams.OR.lbmhft&
                     .OR.lipbv.OR.lbuckmo.OR.lgp)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lbuckmo) THEN
       ipot  = bm_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lwilliams.OR.lbmhft&
                    .OR.lipbv.OR.lbuck4r.OR.lgp)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lgp) THEN
       ipot  = gp_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lwilliams.OR.lbmhft&
                    .OR.lipbv.OR.lbuck4r.OR.lbuckmo)
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
          p % pot ( i ) % pot % type  = nn_type
          p % pot ( i ) % pot % shell_type  = nosh_nosh
          p % pot ( i ) % pot % undef = .TRUE.
          p % pot ( i ) % pot % at1   = 'NULL'
          p % pot ( i ) % pot % at2   = 'NULL'
          SELECT CASE (ipot)
          CASE (lj_type,lj_charmm_type)
             CALL pair_potential_lj_create(p%pot(i)%pot, error=error)
          CASE (wl_type)
             CALL pair_potential_williams_create(p%pot(i)%pot, error=error)
          CASE (gw_type)
             CALL pair_potential_goodwin_create(p%pot(i)%pot, error=error)
          CASE (ea_type)
             CALL pair_potential_eam_create(p%pot(i)%pot, error=error)
          CASE (ft_type)
             CALL pair_potential_bmhft_create(p%pot(i)%pot, error=error)
          CASE (ip_type)
             CALL pair_potential_ipbv_create(p%pot(i)%pot, error=error)
          CASE (b4_type)
             CALL pair_potential_buck4r_create(p%pot(i)%pot, error=error)
          CASE (bm_type)
             CALL pair_potential_buckmo_create(p%pot(i)%pot, error=error)
          CASE (gp_type)
             CALL pair_potential_gp_create(p%pot(i)%pot, error=error)
          END SELECT
          NULLIFY ( p%pot(i)%pot%pair_spline_data )
          NULLIFY ( p%pot(i)%pot%rho_spline_data )
          NULLIFY ( p%pot(i)%pot%c_coul_spline_data )
          NULLIFY ( p%pot(i)%pot%sc_coul_spline_data )
          NULLIFY ( p%pot(i)%pot%cs_coul_spline_data )
       END DO

       IF (ASSOCIATED(work)) CALL pair_potential_p_release(work, error=error)
    END IF
  END SUBROUTINE pair_potential_reallocate

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_gp_create [1.0] *
  !!
  !!   NAME
  !!     pair_potential_gp_create
  !!
  !!   FUNCTION
  !!     Creates the generic potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_gp_create(potparm, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_gp_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(potparm%gp),cp_failure_level,routineP,error,failure)
       ALLOCATE(potparm%gp, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_gp_clean(potparm, error=error)
    END IF
  END SUBROUTINE pair_potential_gp_create

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_gp_copy [1.0] *
  !!
  !!   NAME
  !!     pair_potential_gp_copy
  !!
  !!   FUNCTION
  !!     Copy two generic potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_gp_copy(potparm_source, potparm_dest, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm_source, potparm_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_gp_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idim, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(potparm_dest),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm_source%gp)) RETURN
       IF (ASSOCIATED(potparm_dest%gp)) CALL pair_potential_gp_release(potparm_dest, error=error)
       CALL pair_potential_gp_create(potparm_dest, error=error)
       potparm_dest% gp % myid      = potparm_source% gp % myid
       potparm_dest% gp % potential = potparm_source% gp % potential
       potparm_dest% gp % variables = potparm_source% gp % variables
       IF (ASSOCIATED(potparm_source% gp % parameters)) THEN
          idim = SIZE(potparm_source% gp % parameters)
          ALLOCATE(potparm_dest% gp % parameters(idim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          potparm_dest% gp % parameters = potparm_source% gp % parameters
       END IF
       IF (ASSOCIATED(potparm_source% gp % values)) THEN
          idim = SIZE(potparm_source% gp % values)
          ALLOCATE(potparm_dest% gp % values(idim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          potparm_dest% gp % values = potparm_source% gp % values
       END IF
    END IF
  END SUBROUTINE pair_potential_gp_copy

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_gp_clean [1.0] *
  !!
  !!   NAME
  !!     pair_potential_gp_clean
  !!
  !!   FUNCTION
  !!     Cleans the generic potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_gp_clean(potparm, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_gp_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm % gp)) RETURN
       potparm % gp % myid      = 0
       potparm % gp % potential = ""
       potparm % gp % variables = ""
       NULLIFY (potparm % gp % parameters)
       NULLIFY (potparm % gp % values)
    END IF
  END SUBROUTINE pair_potential_gp_clean

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_gp_release [1.0] *
  !!
  !!   NAME
  !!     pair_potential_gp_release
  !!
  !!   FUNCTION
  !!     Destroys the generic potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_gp_release(potparm, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_gp_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(potparm%gp)) THEN
          IF (ASSOCIATED(potparm%gp%parameters)) THEN
             DEALLOCATE(potparm%gp%parameters, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(potparm%gp%values)) THEN
             DEALLOCATE(potparm%gp%values, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          DEALLOCATE(potparm%gp, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(potparm%gp)
    END IF
  END SUBROUTINE pair_potential_gp_release

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_lj_create [1.0] *
  !!
  !!   NAME
  !!     pair_potential_lj_create
  !!
  !!   FUNCTION
  !!     Cleans the LJ potential type
  !!
  !!   AUTHOR
  !!     Teodoro Laino [teo] 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_lj_create(potparm, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_lj_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm_source, potparm_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_lj_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_lj_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_lj_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_williams_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm_source, potparm_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_williams_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_williams_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_williams_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_goodwin_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm_source, potparm_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_goodwin_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_goodwin_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_goodwin_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_eam_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm_source, potparm_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_eam_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_eam_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_eam_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_bmhft_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm_source, potparm_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_bmhft_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_bmhft_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_bmhft_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_ipbv_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm_source, potparm_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_ipbv_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_ipbv_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_ipbv_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(potparm%ipbv)) THEN
          DEALLOCATE(potparm%ipbv, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(potparm%ipbv)
    END IF
  END SUBROUTINE pair_potential_ipbv_release

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_buck4r_create [1.0] *
  !!
  !!   NAME
  !!     pair_potential_buck4r_create
  !!
  !!   FUNCTION
  !!     Creates the Buckingham 4 ranges  potential type
  !!
  !!   AUTHOR
  !!     MI 10.2006
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_buck4r_create(potparm, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_buck4r_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(potparm%buck4r),cp_failure_level,routineP,error,failure)
       ALLOCATE(potparm%buck4r, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_buck4r_clean(potparm, error=error)
    END IF
  END SUBROUTINE pair_potential_buck4r_create

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_buck4r_copy [1.0] *
  !!
  !!   NAME
  !!     pair_potential_buck4r_copy
  !!
  !!   FUNCTION
  !!     Copy two Buckingham 4 ranges  potential type
  !!
  !!   AUTHOR
  !!     MI 10.2006
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_buck4r_copy(potparm_source, potparm_dest, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm_source, potparm_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_buck4r_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(potparm_dest),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm_source%buck4r)) RETURN
       IF (ASSOCIATED(potparm_dest%buck4r)) CALL pair_potential_buck4r_release(potparm_dest, error=error)
       CALL pair_potential_buck4r_create(potparm_dest, error=error)
       potparm_dest% buck4r % a = potparm_source% buck4r % a
       potparm_dest% buck4r % b = potparm_source% buck4r % b
       potparm_dest% buck4r % c = potparm_source% buck4r % c
       potparm_dest% buck4r % r1 = potparm_source% buck4r % r1
       potparm_dest% buck4r % r2 = potparm_source% buck4r % r2
       potparm_dest% buck4r % r3 = potparm_source% buck4r % r3
       potparm_dest% buck4r % poly1 = potparm_source% buck4r % poly1
       potparm_dest% buck4r % poly2 = potparm_source% buck4r % poly2
       potparm_dest% buck4r % npoly1 = potparm_source% buck4r % npoly1
       potparm_dest% buck4r % npoly2 = potparm_source% buck4r % npoly2
    END IF
  END SUBROUTINE pair_potential_buck4r_copy

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_buck4r_clean [1.0] *
  !!
  !!   NAME
  !!     pair_potential_buck4r_clean
  !!
  !!   FUNCTION
  !!     Creates the Buckingham 4 ranges  potential type
  !!
  !!   AUTHOR
  !!     MI 10.2006
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_buck4r_clean(potparm, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_buck4r_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(potparm%buck4r)) RETURN
       potparm% buck4r % a = 0.0_dp
       potparm% buck4r % b = 0.0_dp
       potparm% buck4r % c = 0.0_dp
       potparm% buck4r % r1 = 0.0_dp
       potparm% buck4r % r2 = 0.0_dp
       potparm% buck4r % r3 = 0.0_dp
       potparm% buck4r % poly1 = 0.0_dp
       potparm% buck4r % npoly1 = 0
       potparm% buck4r % poly2 = 0.0_dp
       potparm% buck4r % npoly2 = 0
    END IF
  END SUBROUTINE pair_potential_buck4r_clean

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_buck4r_release [1.0] *
  !!
  !!   NAME
  !!     pair_potential_buck4r_release
  !!
  !!   FUNCTION
  !!     Destroys the Buckingham 4 ranges potential type
  !!
  !!   AUTHOR
  !!     MI  10.2006
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_buck4r_release(potparm, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_buck4r_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(potparm%buck4r)) THEN
          DEALLOCATE(potparm%buck4r, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(potparm%willis)
    END IF
  END SUBROUTINE pair_potential_buck4r_release

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_buckmo_create [1.0] *
  !!
  !!   NAME
  !!     pair_potential_buckmo_create
  !!
  !!   FUNCTION
  !!     Creates the Buckingham plus Morse potential type
  !!
  !!   AUTHOR
  !!     MI 10.2006
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_buckmo_create(potparm, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_buckmo_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(potparm%buckmo),cp_failure_level,routineP,error,failure)
       ALLOCATE(potparm%buckmo, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_buckmo_clean(potparm, error=error)
    END IF
  END SUBROUTINE pair_potential_buckmo_create

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_buckmo_copy [1.0] *
  !!
  !!   NAME
  !!     pair_potential_buckmo_copy
  !!
  !!   FUNCTION
  !!     Copy two Buckingham plus Morse  potential type
  !!
  !!   AUTHOR
  !!     MI 10.2006
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_buckmo_copy(potparm_source, potparm_dest, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm_source, potparm_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_buckmo_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(potparm_dest),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm_source%buckmo)) RETURN
       IF (ASSOCIATED(potparm_dest%buckmo)) CALL pair_potential_buckmo_release(potparm_dest, error=error)
       CALL pair_potential_buckmo_create(potparm_dest, error=error)
       potparm_dest% buckmo % f0 = potparm_source% buckmo % f0
       potparm_dest% buckmo % a1 = potparm_source% buckmo % a1
       potparm_dest% buckmo % a2 = potparm_source% buckmo % a2
       potparm_dest% buckmo % b1 = potparm_source% buckmo % b1
       potparm_dest% buckmo % b2 = potparm_source% buckmo % b2
       potparm_dest% buckmo % c = potparm_source% buckmo % c
       potparm_dest% buckmo % d  = potparm_source% buckmo % d
       potparm_dest% buckmo % r0 = potparm_source% buckmo % r0
       potparm_dest% buckmo % beta = potparm_source% buckmo % beta
    END IF
  END SUBROUTINE pair_potential_buckmo_copy

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_buckmo_clean [1.0] *
  !!
  !!   NAME
  !!     pair_potential_buckmo_clean
  !!
  !!   FUNCTION
  !!     Creates the Buckingham plus Morse  potential type
  !!
  !!   AUTHOR
  !!     MI 10.2006
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_buckmo_clean(potparm, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_buckmo_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(potparm%buckmo)) RETURN
       potparm% buckmo % f0 = 0.0_dp
       potparm% buckmo % a1 = 0.0_dp
       potparm% buckmo % a2 = 0.0_dp
       potparm% buckmo % b1 = 0.0_dp
       potparm% buckmo % b2 = 0.0_dp
       potparm% buckmo % c = 0.0_dp
       potparm% buckmo % d  = 0.0_dp
       potparm% buckmo % r0 = 0.0_dp
       potparm% buckmo % beta = 0.0_dp
    END IF
  END SUBROUTINE pair_potential_buckmo_clean

  !******************************************************************************
  !!****
  !!****s* pair_potential_types/pair_potential_buckmo_release [1.0] *
  !!
  !!   NAME
  !!     pair_potential_buckmo_release
  !!
  !!   FUNCTION
  !!     Destroys the Buckingham plus Morse potential type
  !!
  !!   AUTHOR
  !!     MI  10.2006
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pair_potential_buckmo_release(potparm, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_buckmo_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(potparm%buckmo)) THEN
          DEALLOCATE(potparm%buckmo, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(potparm%willis)
    END IF
  END SUBROUTINE pair_potential_buckmo_release

  !******************************************************************************

END MODULE pair_potential_types


