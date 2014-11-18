!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Teodoro Laino [Teo] 11.2005 : Reorganizing the structures to optimize
!>                                    memory management
!> \author CJM
! *****************************************************************************
MODULE pair_potential_types

  
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE splines_types,                   ONLY: spline_data_p_copy,&
                                             spline_data_p_release,&
                                             spline_data_p_type,&
                                             spline_factor_copy,&
                                             spline_factor_release,&
                                             spline_factor_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pair_potential_types'

  PRIVATE
  ! when adding a new nonbonedd potential please update also the list_pot
  ! used for the linear scaling screening of potential calculation
  INTEGER, PUBLIC, PARAMETER :: multi_type=-1,&
                                nn_type= 0,&
                                lj_type= 1,&
                                lj_charmm_type= 2,&
                                ft_type= 3,&
                                wl_type= 4,&
                                gw_type= 5,&
                                ip_type= 6,&
                                ea_type= 7,&
                                b4_type= 8,&
                                bm_type= 9,&
                                gp_type=10,&
                                tersoff_type=11,&
                                ftd_type=12,&
                                siepmann_type=13, &
                                quip_type= 14

  INTEGER, PUBLIC, PARAMETER, DIMENSION(15) :: list_pot = (/nn_type,&
                                                           lj_type,&
                                                           lj_charmm_type,&
                                                           ft_type,&
                                                           wl_type,&
                                                           gw_type,&
                                                           ip_type,&
                                                           ea_type,&
                                                           b4_type,&
                                                           bm_type,&
                                                           gp_type,&
                                                           tersoff_type,&
                                                           ftd_type,&
                                                           siepmann_type, &
                                                           quip_type /)
  ! Shell model
  INTEGER, PUBLIC, PARAMETER :: nosh_nosh=0,&
                                nosh_sh=1,&
                                sh_sh=2

  INTEGER, PUBLIC, PARAMETER, DIMENSION(3) :: list_sh_type =(/nosh_nosh,nosh_sh,sh_sh/)

  ! Single Spline generation info
  REAL(KIND=dp), PARAMETER, PUBLIC         :: not_initialized=-HUGE(0.0_dp)
  INTEGER, PARAMETER, DIMENSION(2), PUBLIC :: do_potential_single_allocation=(/lj_type,lj_charmm_type/)
  INTEGER, PARAMETER, DIMENSION(2), PUBLIC :: no_potential_single_allocation=(/-HUGE(0),-HUGE(0)/)
  INTEGER, DIMENSION(2), PUBLIC            :: potential_single_allocation

  PUBLIC :: pair_potential_reallocate

  PUBLIC :: pair_potential_single_copy,&
            pair_potential_single_add,&
            pair_potential_single_clean,&
            pair_potential_single_type

  PUBLIC :: pair_potential_pp_create,&
            pair_potential_pp_release,&
            pair_potential_pp_type

  PUBLIC :: pair_potential_p_type,&
            pair_potential_p_release

  PUBLIC :: ft_pot_type,&
            ipbv_pot_type,&
            eam_pot_type,&
            quip_pot_type, &
            tersoff_pot_type,&
            siepmann_pot_type,&
            ftd_pot_type

  PUBLIC :: pair_potential_lj_create
  PUBLIC :: compare_pot

! *****************************************************************************
  TYPE ipbv_pot_type
     REAL (KIND=dp), DIMENSION (2:15) :: a
     REAL (KIND=dp) :: rcore
     REAL (KIND=dp) :: m
     REAL (KIND=dp) :: b
  END TYPE ipbv_pot_type

! *****************************************************************************
  TYPE lj_pot_type
     REAL (KIND=dp) :: epsilon
     REAL (KIND=dp) :: sigma6
     REAL (KIND=dp) :: sigma12
  END TYPE Lj_pot_type

! *****************************************************************************
  TYPE ft_pot_type
     REAL (KIND=dp) :: A
     REAL (KIND=dp) :: B
     REAL (KIND=dp) :: C
     REAL (KIND=dp) :: D
  END TYPE ft_pot_type

! *****************************************************************************
  TYPE ftd_pot_type
     REAL (KIND=dp) :: A
     REAL (KIND=dp) :: B
     REAL (KIND=dp) :: C
     REAL (KIND=dp) :: D
     REAL (KIND=dp) :: BD
  END TYPE ftd_pot_type

! *****************************************************************************
  TYPE williams_pot_type
     REAL (KIND=dp) :: a
     REAL (KIND=dp) :: b
     REAL (KIND=dp) :: c
  END TYPE williams_pot_type

! *****************************************************************************
  TYPE goodwin_pot_type
     REAL (KIND=dp) :: vr0
     REAL (KIND=dp) :: m, mc
     REAL (KIND=dp) :: d, dc
  END TYPE goodwin_pot_type

! *****************************************************************************
  TYPE eam_pot_type
     CHARACTER (LEN=default_path_length)     :: eam_file_name
     INTEGER                                 :: npoints
     REAL (KIND=dp)                          :: drar, drhoar, acutal
     REAL (KIND=dp),POINTER, DIMENSION(:)    :: rho,  phi,  frho,  rhoval, rval
     REAL (KIND=dp),POINTER, DIMENSION(:)    :: rhop, phip, frhop
  END TYPE eam_pot_type

! *****************************************************************************
  TYPE quip_pot_type
     CHARACTER (LEN=default_path_length)     :: quip_file_name
     CHARACTER (LEN=1024)                    :: init_args
     CHARACTER (LEN=1024)                    :: calc_args
  END TYPE quip_pot_type

! *****************************************************************************
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

! *****************************************************************************
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

! *****************************************************************************
  TYPE gp_pot_type
     INTEGER                                   :: myid
     CHARACTER (LEN=default_path_length)       :: potential
     CHARACTER (LEN=default_string_length), &
       POINTER, DIMENSION(:)                   :: parameters, units
     CHARACTER (LEN=default_string_length)     :: variables
     REAL(KIND=dp), DIMENSION(:), POINTER      :: values
  END TYPE gp_pot_type

! *****************************************************************************
  TYPE tersoff_pot_type
     ! Get this stuff from the PRB V38, N14 9902 (1988) by Tersoff
     REAL (KIND=dp) :: A
     REAL (KIND=dp) :: B
     REAL (KIND=dp) :: lambda1
     REAL (KIND=dp) :: lambda2
     REAL (KIND=dp) :: alpha
     REAL (KIND=dp) :: beta
     REAL (KIND=dp) :: n
     REAL (KIND=dp) :: c
     REAL (KIND=dp) :: d
     REAL (KIND=dp) :: h
     REAL (KIND=dp) :: lambda3
     REAL (KIND=dp) :: bigR          ! Used to be R = Rij + D
     REAL (KIND=dp) :: bigD          ! Used to be D = Rij - D
     REAL (KIND=dp) :: rcutsq        ! Always set to (bigR+bigD)^2
  END TYPE tersoff_pot_type

! *****************************************************************************
  TYPE siepmann_pot_type
     REAL (KIND=dp) :: B
     REAL (KIND=dp) :: D
     REAL (KIND=dp) :: E
     REAL (KIND=dp) :: F
     REAL (KIND=dp) :: alpha
     REAL (KIND=dp) :: beta
     REAL (KIND=dp) :: rcutsq
  END TYPE siepmann_pot_type

! *****************************************************************************
  TYPE pot_set_type
     REAL ( KIND=dp )                      :: rmin, rmax
     TYPE ( ipbv_pot_type ),       POINTER :: ipbv
     TYPE ( gp_pot_type),          POINTER :: gp
     TYPE ( lj_pot_type ),         POINTER :: lj
     TYPE ( ft_pot_type ),         POINTER :: ft
     TYPE ( williams_pot_type ),   POINTER :: willis
     TYPE ( goodwin_pot_type ),    POINTER :: goodwin
     TYPE ( eam_pot_type ),        POINTER :: eam
     TYPE ( quip_pot_type ),       POINTER :: quip
     TYPE ( buck4ran_pot_type ),   POINTER :: buck4r
     TYPE ( buckmorse_pot_type ),  POINTER :: buckmo
     TYPE ( tersoff_pot_type ),    POINTER :: tersoff
     TYPE ( siepmann_pot_type ),   POINTER :: siepmann
     TYPE ( ftd_pot_type ),        POINTER :: ftd
  END TYPE pot_set_type

! *****************************************************************************
  TYPE pair_potential_single_type
     REAL (KIND=dp) :: rcutsq
     REAL (KIND=dp) :: e_fac
     REAL (KIND=dp) :: e_fcc
     REAL (KIND=dp) :: e_fcs
     REAL (KIND=dp) :: e_fsc
     REAL (KIND=dp) :: z1
     REAL (KIND=dp) :: z2
     REAL (KIND=dp) :: zbl_poly(0:5)
     REAL (KIND=dp) :: zbl_rcut(2)
     LOGICAL        :: undef, &  ! non-bonding interaction not defined
                       no_mb, &  ! no many-body potential
                       no_pp     ! no pair (=two-body) potential
     INTEGER        :: shell_type
     CHARACTER ( LEN = default_string_length )   :: at1
     CHARACTER ( LEN = default_string_length )   :: at2
     INTEGER, POINTER, DIMENSION(:)              :: TYPE
     TYPE (pot_set_type), POINTER, DIMENSION(:) :: set
     TYPE (spline_data_p_type), POINTER, DIMENSION(:) :: pair_spline_data
     TYPE (spline_factor_type), POINTER :: spl_f
  END TYPE pair_potential_single_type

! *****************************************************************************
  TYPE pair_potential_type
     TYPE(pair_potential_single_type),          POINTER :: pot
  END TYPE pair_potential_type

! *****************************************************************************
  TYPE pair_potential_p_type
     TYPE(pair_potential_type), DIMENSION(:),   POINTER :: pot
  END TYPE pair_potential_p_type

! *****************************************************************************
  TYPE pair_potential_pp_type
     TYPE(pair_potential_type), DIMENSION(:,:), POINTER :: pot
  END TYPE pair_potential_pp_type

CONTAINS

! *****************************************************************************
!> \brief compare two different potentials
!> \param pot1 ...
!> \param pot2 ...
!> \param compare ...
!> \param error ...
!> \author Teodoro Laino [teo] 05.2006
! *****************************************************************************
  SUBROUTINE compare_pot(pot1, pot2, compare, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot1, pot2
    LOGICAL, INTENT(OUT)                     :: compare
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compare_pot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, mycompare

    failure=.FALSE.
    compare=.FALSE.
    ! Preliminary checks

    CPPostcondition(ASSOCIATED(pot1%type),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(pot2%type),cp_failure_level,routineP,error,failure)
    IF (SIZE(pot1%type)/=SIZE(pot2%type)) RETURN
    IF (ANY(pot1%type/=pot2%type))        RETURN

    ! Checking the real values of parameters
    CPPostcondition(ASSOCIATED(pot1%set),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(pot2%set),cp_failure_level,routineP,error,failure)
    DO i = 1, SIZE(pot1%type)
       mycompare = .FALSE.
       SELECT CASE(pot1%type(i))
       CASE (lj_type,lj_charmm_type)
          IF  ((pot1%set(i)%lj%epsilon==pot2%set(i)%lj%epsilon).AND.&
               (pot1%set(i)%lj%sigma6 ==pot2%set(i)%lj%sigma6 ).AND.&
               (pot1%set(i)%lj%sigma12==pot2%set(i)%lj%sigma12)) mycompare=.TRUE.
       CASE (wl_type)
          IF  ((pot1%set(i)%willis%a==pot2%set(i)%willis%a).AND.&
               (pot1%set(i)%willis%b==pot2%set(i)%willis%b).AND.&
               (pot1%set(i)%willis%c==pot2%set(i)%willis%c)) mycompare=.TRUE.
       CASE (gw_type)
          IF  ((pot1%set(i)%goodwin%vr0==pot2%set(i)%goodwin%vr0).AND.&
               (pot1%set(i)%goodwin%m  ==pot2%set(i)%goodwin%m  ).AND.&
               (pot1%set(i)%goodwin%mc ==pot2%set(i)%goodwin%mc ).AND.&
               (pot1%set(i)%goodwin%d  ==pot2%set(i)%goodwin%d  ).AND.&
               (pot1%set(i)%goodwin%dc ==pot2%set(i)%goodwin%dc )) mycompare=.TRUE.
       CASE (ea_type)
          ! Compare only if EAM have the same number of points
          IF (pot1%set(i)%eam%npoints==pot2%set(i)%eam%npoints) THEN
             IF  ((pot1%set(i)%eam%drar   ==pot2%set(i)%eam%drar  ).AND.&
                  (pot1%set(i)%eam%drhoar ==pot2%set(i)%eam%drhoar).AND.&
                  (pot1%set(i)%eam%acutal ==pot2%set(i)%eam%acutal).AND.&
                  (SUM(ABS(pot1%set(i)%eam%rho   -pot2%set(i)%eam%rho   ))==0.0_dp).AND.&
                  (SUM(ABS(pot1%set(i)%eam%phi   -pot2%set(i)%eam%phi   ))==0.0_dp).AND.&
                  (SUM(ABS(pot1%set(i)%eam%frho  -pot2%set(i)%eam%frho  ))==0.0_dp).AND.&
                  (SUM(ABS(pot1%set(i)%eam%rhoval-pot2%set(i)%eam%rhoval))==0.0_dp).AND.&
                  (SUM(ABS(pot1%set(i)%eam%rval  -pot2%set(i)%eam%rval  ))==0.0_dp).AND.&
                  (SUM(ABS(pot1%set(i)%eam%rhop  -pot2%set(i)%eam%rhop  ))==0.0_dp).AND.&
                  (SUM(ABS(pot1%set(i)%eam%phip  -pot2%set(i)%eam%phip  ))==0.0_dp).AND.&
                  (SUM(ABS(pot1%set(i)%eam%frhop -pot2%set(i)%eam%frhop ))==0.0_dp)) mycompare=.TRUE.
          END IF
       CASE (quip_type)
          IF  ((pot1%set(i)%quip%quip_file_name ==pot2%set(i)%quip%quip_file_name  ).AND.&
               (pot1%set(i)%quip%init_args ==pot2%set(i)%quip%init_args).AND.&
               (pot1%set(i)%quip%calc_args ==pot2%set(i)%quip%calc_args)) mycompare=.TRUE.
       CASE (ft_type)
          IF  ((pot1%set(i)%ft%A==pot2%set(i)%ft%A).AND.&
               (pot1%set(i)%ft%B==pot2%set(i)%ft%B).AND.&
               (pot1%set(i)%ft%C==pot2%set(i)%ft%C).AND.&
               (pot1%set(i)%ft%D==pot2%set(i)%ft%D)) mycompare=.TRUE.
       CASE (ftd_type)
          IF  ((pot1%set(i)%ftd%A==pot2%set(i)%ftd%A).AND.&
               (pot1%set(i)%ftd%B==pot2%set(i)%ftd%B).AND.&
               (pot1%set(i)%ftd%C==pot2%set(i)%ftd%C).AND.&
               (pot1%set(i)%ftd%D==pot2%set(i)%ftd%D).AND.&
               (pot1%set(i)%ftd%BD==pot2%set(i)%ftd%BD)) mycompare=.TRUE.
       CASE (ip_type)
          IF  ((SUM(ABS(pot1%set(i)%ipbv%a-pot2%set(i)%ipbv%a))==0.0_dp).AND.&
               (pot1%set(i)%ipbv%rcore==pot2%set(i)%ipbv%rcore).AND.&
               (pot1%set(i)%ipbv%m    ==pot2%set(i)%ipbv%m    ).AND.&
               (pot1%set(i)%ipbv%b    ==pot2%set(i)%ipbv%b    )) mycompare=.TRUE.
       CASE (tersoff_type)
          IF  ((pot1%set(i)%tersoff%A       == pot2%set(i)%tersoff%A       ).AND.&
               (pot1%set(i)%tersoff%B       == pot2%set(i)%tersoff%B       ).AND.&
               (pot1%set(i)%tersoff%lambda1 == pot2%set(i)%tersoff%lambda1 ).AND.&
               (pot1%set(i)%tersoff%lambda2 == pot2%set(i)%tersoff%lambda2 ).AND.&
               (pot1%set(i)%tersoff%alpha   == pot2%set(i)%tersoff%alpha   ).AND.&
               (pot1%set(i)%tersoff%beta    == pot2%set(i)%tersoff%beta    ).AND.&
               (pot1%set(i)%tersoff%n       == pot2%set(i)%tersoff%n       ).AND.&
               (pot1%set(i)%tersoff%c       == pot2%set(i)%tersoff%c       ).AND.&
               (pot1%set(i)%tersoff%d       == pot2%set(i)%tersoff%d       ).AND.&
               (pot1%set(i)%tersoff%h       == pot2%set(i)%tersoff%h       ).AND.&
               (pot1%set(i)%tersoff%lambda3 == pot2%set(i)%tersoff%lambda3 ).AND.&
               (pot1%set(i)%tersoff%rcutsq  == pot2%set(i)%tersoff%rcutsq  ).AND.&
               (pot1%set(i)%tersoff%bigR    == pot2%set(i)%tersoff%bigR    ).AND.&
               (pot1%set(i)%tersoff%bigD    == pot2%set(i)%tersoff%bigD    )) mycompare=.TRUE.
       CASE (siepmann_type)
          IF  ((pot1%set(i)%siepmann%B       == pot2%set(i)%siepmann%B       ).AND.&
               (pot1%set(i)%siepmann%D       == pot2%set(i)%siepmann%D       ).AND.&
               (pot1%set(i)%siepmann%E       == pot2%set(i)%siepmann%E       ).AND.&
               (pot1%set(i)%siepmann%F       == pot2%set(i)%siepmann%F       ).AND.&
               (pot1%set(i)%siepmann%alpha   == pot2%set(i)%siepmann%alpha   ).AND.&
               (pot1%set(i)%siepmann%beta    == pot2%set(i)%siepmann%beta    ).AND.&
               (pot1%set(i)%siepmann%rcutsq  == pot2%set(i)%siepmann%rcutsq  )) mycompare=.TRUE.
       END SELECT
       mycompare = mycompare.AND.&
                  (pot1%set(i)%rmin == pot2%set(i)%rmin ).AND.(pot1%set(i)%rmax == pot2%set(i)%rmax )
       IF ((mycompare).AND.(i==1)) compare = .TRUE.
       compare = compare .AND. mycompare
    END DO

  END SUBROUTINE compare_pot

! *****************************************************************************
!> \brief Creates the potential parameter type
!> \param potparm ...
!> \param nset ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_single_create(potparm, nset, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    INTEGER, INTENT(IN), OPTIONAL            :: nset
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pair_potential_single_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, lnset, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(.NOT.ASSOCIATED(potparm),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       ALLOCATE(potparm, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       lnset = 1
       IF (PRESENT(nset)) lnset = nset
       ! Standard allocation to size 1
       ALLOCATE(potparm%type(lnset), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(potparm%set(lnset), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY( potparm%spl_f,&
                potparm%pair_spline_data)
       DO i = 1, lnset
          potparm%set(i)%rmin = not_initialized
          potparm%set(i)%rmax = not_initialized
          NULLIFY( potparm%set(i)%ipbv,&
                   potparm%set(i)%lj,&
                   potparm%set(i)%gp,&
                   potparm%set(i)%ft,&
                   potparm%set(i)%willis,&
                   potparm%set(i)%goodwin,&
                   potparm%set(i)%eam,&
                   potparm%set(i)%quip,&
                   potparm%set(i)%buck4r,&
                   potparm%set(i)%buckmo,&
                   potparm%set(i)%tersoff,&
                   potparm%set(i)%siepmann,&
                   potparm%set(i)%ftd)
       END DO
       CALL pair_potential_single_clean(potparm, error)
    END IF
  END SUBROUTINE pair_potential_single_create

! *****************************************************************************
!> \brief Cleans the potential parameter type
!> \param potparm ...
!> \param error ...
!> \author unknown
! *****************************************************************************
  SUBROUTINE pair_potential_single_clean(potparm, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pair_potential_single_clean', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure = .FALSE.
    potparm % type          = nn_type
    potparm % shell_type    = nosh_nosh
    potparm % undef         = .TRUE.
    potparm % no_pp         = .FALSE.
    potparm % no_mb         = .FALSE.
    potparm % at1           = 'NULL'
    potparm % at2           = 'NULL'
    potparm % rcutsq        = 0.0_dp
    IF (ASSOCIATED( potparm % pair_spline_data ))&
         CALL spline_data_p_release ( potparm % pair_spline_data ,error=error)
    IF (ASSOCIATED( potparm % spl_f ))&
         CALL spline_factor_release ( potparm % spl_f ,error=error)

    DO i = 1, SIZE(potparm%type)
       potparm%set(i)%rmin = not_initialized
       potparm%set(i)%rmax = not_initialized
       CALL pair_potential_lj_clean (potparm%set(i)%lj, error=error)
       CALL pair_potential_williams_clean (potparm%set(i)%willis, error=error)
       CALL pair_potential_goodwin_clean (potparm%set(i)%goodwin, error=error)
       CALL pair_potential_eam_clean (potparm%set(i)%eam, error=error)
       CALL pair_potential_quip_clean (potparm%set(i)%quip, error=error)
       CALL pair_potential_buck4r_clean (potparm%set(i)%buck4r, error=error)
       CALL pair_potential_buckmo_clean (potparm%set(i)%buckmo, error=error)
       CALL pair_potential_bmhft_clean (potparm%set(i)%ft, error=error)
       CALL pair_potential_bmhftd_clean (potparm%set(i)%ftd, error=error)
       CALL pair_potential_ipbv_clean (potparm%set(i)%ipbv, error=error)
       CALL pair_potential_gp_clean (potparm%set(i)%gp, error=error)
       CALL pair_potential_tersoff_clean (potparm%set(i)%tersoff, error=error)
       CALL pair_potential_siepmann_clean (potparm%set(i)%siepmann, error=error)
    END DO
  END SUBROUTINE pair_potential_single_clean

! *****************************************************************************
!> \brief Copy two potential parameter type
!> \param potparm_source ...
!> \param potparm_dest ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_single_copy(potparm_source, potparm_dest, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm_source, potparm_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pair_potential_single_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(potparm_dest)) THEN
          CALL pair_potential_single_create(potparm_dest,SIZE(potparm_source%type),error)
       ELSE
          CALL pair_potential_single_clean(potparm_dest, error)
       END IF
       potparm_dest % type          = potparm_source % type
       potparm_dest % shell_type    = potparm_source % shell_type
       potparm_dest % undef         = potparm_source % undef
       potparm_dest % no_mb         = potparm_source % no_mb
       potparm_dest % no_pp         = potparm_source % no_pp
       potparm_dest % at1           = potparm_source % at1
       potparm_dest % at2           = potparm_source % at2
       potparm_dest % rcutsq        = potparm_source % rcutsq
       IF (ASSOCIATED( potparm_source % pair_spline_data )) THEN
          CALL spline_data_p_copy ( potparm_source % pair_spline_data, potparm_dest % pair_spline_data, error=error  )
       END IF

       IF (ASSOCIATED( potparm_source % spl_f )) THEN
          CALL spline_factor_copy ( potparm_source % spl_f,  potparm_dest % spl_f, error=error  )
       END IF

       DO i = 1, SIZE(potparm_source%type)
          potparm_dest%set(i)%rmin = potparm_source%set(i)%rmin
          potparm_dest%set(i)%rmax = potparm_source%set(i)%rmax
          CALL pair_potential_lj_copy (potparm_source%set(i)%lj, potparm_dest%set(i)%lj, error=error)
          CALL pair_potential_williams_copy (potparm_source%set(i)%willis, potparm_dest%set(i)%willis, error=error)
          CALL pair_potential_goodwin_copy (potparm_source%set(i)%goodwin, potparm_dest%set(i)%goodwin, error=error)
          CALL pair_potential_eam_copy (potparm_source%set(i)%eam, potparm_dest%set(i)%eam, error=error)
          CALL pair_potential_quip_copy (potparm_source%set(i)%quip, potparm_dest%set(i)%quip, error=error)
          CALL pair_potential_bmhft_copy (potparm_source%set(i)%ft, potparm_dest%set(i)%ft, error=error)
          CALL pair_potential_bmhftd_copy (potparm_source%set(i)%ftd, potparm_dest%set(i)%ftd, error=error)
          CALL pair_potential_ipbv_copy (potparm_source%set(i)%ipbv, potparm_dest%set(i)%ipbv, error=error)
          CALL pair_potential_buck4r_copy (potparm_source%set(i)%buck4r, potparm_dest%set(i)%buck4r, error=error)
          CALL pair_potential_buckmo_copy (potparm_source%set(i)%buckmo, potparm_dest%set(i)%buckmo, error=error)
          CALL pair_potential_gp_copy (potparm_source%set(i)%gp, potparm_dest%set(i)%gp, error=error)
          CALL pair_potential_tersoff_copy (potparm_source%set(i)%tersoff, potparm_dest%set(i)%tersoff, error=error)
          CALL pair_potential_siepmann_copy (potparm_source%set(i)%siepmann, potparm_dest%set(i)%siepmann, error=error)
       END DO
    END IF
  END SUBROUTINE pair_potential_single_copy

! *****************************************************************************
!> \brief Add potential parameter type to an existing potential parameter type
!>      Used in case of multiple_potential definition
!> \param potparm_source ...
!> \param potparm_dest ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_single_add(potparm_source, potparm_dest, error)
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm_source, potparm_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pair_potential_single_add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, size_dest, size_source
    LOGICAL                                  :: allocate_new, check, failure
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm_tmp

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(potparm_source),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       ! At this level we expect all splines types
       ! be not allocated.. No sense add splines at this level.. in case fail!
       check = (.NOT.ASSOCIATED( potparm_source % pair_spline_data )).AND.&
               (.NOT.ASSOCIATED( potparm_source % spl_f ))
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       check = (.NOT.ASSOCIATED( potparm_dest % pair_spline_data )).AND.&
               (.NOT.ASSOCIATED( potparm_dest % spl_f ))
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       ! Increase the size of the destination potparm (in case) and copy the new data
       size_source = SIZE(potparm_source%type)
       allocate_new = .NOT.ASSOCIATED(potparm_dest)
       IF (.NOT.allocate_new) THEN
          size_dest = SIZE(potparm_dest%type)
          IF (size_dest==1) THEN
             check = (ASSOCIATED(potparm_dest%set(1)%lj)).OR.&
                     (ASSOCIATED(potparm_dest%set(1)%willis)).OR.&
                     (ASSOCIATED(potparm_dest%set(1)%goodwin)).OR.&
                     (ASSOCIATED(potparm_dest%set(1)%eam)).OR.&
                     (ASSOCIATED(potparm_dest%set(1)%quip)).OR.&
                     (ASSOCIATED(potparm_dest%set(1)%ft)).OR.&
                     (ASSOCIATED(potparm_dest%set(1)%ftd)).OR.&
                     (ASSOCIATED(potparm_dest%set(1)%ipbv)).OR.&
                     (ASSOCIATED(potparm_dest%set(1)%buck4r)).OR.&
                     (ASSOCIATED(potparm_dest%set(1)%buckmo)).OR.&
                     (ASSOCIATED(potparm_dest%set(1)%gp)).OR.&
                     (ASSOCIATED(potparm_dest%set(1)%tersoff)).OR.&
                     (ASSOCIATED(potparm_dest%set(1)%siepmann))
             IF (.NOT.check) THEN
                allocate_new = .TRUE.
                CALL pair_potential_single_release(potparm_dest, error)
             END IF
          END IF
       END IF
       IF (allocate_new) THEN
          size_dest = 0
          CALL pair_potential_single_create(potparm_dest,size_source,error)
          potparm_dest % shell_type    = potparm_source % shell_type
          potparm_dest % undef         = potparm_source % undef
          potparm_dest % no_mb         = potparm_source % no_mb
          potparm_dest % no_pp         = potparm_source % no_pp
          potparm_dest % at1           = potparm_source % at1
          potparm_dest % at2           = potparm_source % at2
          potparm_dest % rcutsq        = potparm_source % rcutsq
       ELSE
          size_dest = SIZE(potparm_dest%type)
          NULLIFY(potparm_tmp)
          CALL pair_potential_single_copy(potparm_dest,potparm_tmp,error)
          CALL pair_potential_single_release(potparm_dest, error)
          CALL pair_potential_single_create(potparm_dest,size_dest+size_source,error)
          ! Copy back original informations..
          potparm_dest % shell_type = potparm_tmp % shell_type
          potparm_dest % undef      = potparm_tmp % undef
          potparm_dest % no_mb      = potparm_tmp % no_mb
          potparm_dest % no_pp      = potparm_tmp % no_pp
          potparm_dest % at1        = potparm_tmp % at1
          potparm_dest % at2        = potparm_tmp % at2
          potparm_dest % rcutsq     = potparm_tmp % rcutsq
          DO i = 1, size_dest
             potparm_dest%type(i) = potparm_tmp%type(i)
             potparm_dest%set(i)%rmin = potparm_tmp%set(i)%rmin
             potparm_dest%set(i)%rmax = potparm_tmp%set(i)%rmax
             CALL pair_potential_lj_copy (potparm_tmp%set(i)%lj, potparm_dest%set(i)%lj, error=error)
             CALL pair_potential_williams_copy (potparm_tmp%set(i)%willis, potparm_dest%set(i)%willis, error=error)
             CALL pair_potential_goodwin_copy (potparm_tmp%set(i)%goodwin, potparm_dest%set(i)%goodwin, error=error)
             CALL pair_potential_eam_copy (potparm_tmp%set(i)%eam, potparm_dest%set(i)%eam, error=error)
             CALL pair_potential_quip_copy (potparm_tmp%set(i)%quip, potparm_dest%set(i)%quip, error=error)
             CALL pair_potential_bmhft_copy (potparm_tmp%set(i)%ft, potparm_dest%set(i)%ft, error=error)
             CALL pair_potential_bmhftd_copy (potparm_tmp%set(i)%ftd, potparm_dest%set(i)%ftd, error=error)
             CALL pair_potential_ipbv_copy (potparm_tmp%set(i)%ipbv, potparm_dest%set(i)%ipbv, error=error)
             CALL pair_potential_buck4r_copy (potparm_tmp%set(i)%buck4r, potparm_dest%set(i)%buck4r, error=error)
             CALL pair_potential_buckmo_copy (potparm_tmp%set(i)%buckmo, potparm_dest%set(i)%buckmo, error=error)
             CALL pair_potential_gp_copy (potparm_tmp%set(i)%gp, potparm_dest%set(i)%gp, error=error)
             CALL pair_potential_tersoff_copy (potparm_tmp%set(i)%tersoff, potparm_dest%set(i)%tersoff, error=error)
             CALL pair_potential_siepmann_copy (potparm_tmp%set(i)%siepmann, potparm_dest%set(i)%siepmann, error=error)
          END DO
          CALL pair_potential_single_release ( potparm_tmp, error )
       END IF
       ! Further check with main option with source and dest (already filled with few informations)
       check= (potparm_dest % shell_type == potparm_source % shell_type).AND.&
              (potparm_dest % undef    .EQV.potparm_source % undef     ).AND.&
              (potparm_dest % no_mb    .EQV.potparm_source % no_mb     ).AND.&
              (potparm_dest % no_pp    .EQV.potparm_source % no_pp     ).AND.&
              (potparm_dest % at1        == potparm_source % at1       ).AND.&
              (potparm_dest % at2        == potparm_source % at2       ).AND.&
              (potparm_dest % rcutsq     == potparm_source % rcutsq    )
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       ! Now copy the new pair_potential type
       DO i = size_dest+1, size_dest+size_source
          j = i-size_dest
          potparm_dest%type(i) = potparm_source%type(j)
          potparm_dest%set(i)%rmin = potparm_source%set(j)%rmin
          potparm_dest%set(i)%rmax = potparm_source%set(j)%rmax
          CALL pair_potential_lj_copy (potparm_source%set(j)%lj, potparm_dest%set(i)%lj, error=error)
          CALL pair_potential_williams_copy (potparm_source%set(j)%willis, potparm_dest%set(i)%willis, error=error)
          CALL pair_potential_goodwin_copy (potparm_source%set(j)%goodwin, potparm_dest%set(i)%goodwin, error=error)
          CALL pair_potential_eam_copy (potparm_source%set(j)%eam, potparm_dest%set(i)%eam, error=error)
          CALL pair_potential_quip_copy (potparm_source%set(j)%quip, potparm_dest%set(i)%quip, error=error)
          CALL pair_potential_bmhft_copy (potparm_source%set(j)%ft, potparm_dest%set(i)%ft, error=error)
          CALL pair_potential_bmhftd_copy (potparm_source%set(j)%ftd, potparm_dest%set(i)%ftd, error=error)
          CALL pair_potential_ipbv_copy (potparm_source%set(j)%ipbv, potparm_dest%set(i)%ipbv, error=error)
          CALL pair_potential_buck4r_copy (potparm_source%set(j)%buck4r, potparm_dest%set(i)%buck4r, error=error)
          CALL pair_potential_buckmo_copy (potparm_source%set(j)%buckmo, potparm_dest%set(i)%buckmo, error=error)
          CALL pair_potential_gp_copy (potparm_source%set(j)%gp, potparm_dest%set(i)%gp, error=error)
          CALL pair_potential_tersoff_copy (potparm_source%set(j)%tersoff, potparm_dest%set(i)%tersoff, error=error)
          CALL pair_potential_siepmann_copy (potparm_source%set(j)%siepmann, potparm_dest%set(i)%siepmann, error=error)
       END DO
    END IF
  END SUBROUTINE pair_potential_single_add

! *****************************************************************************
!> \brief Release Data-structure that constains potential parameters of a single pair
!> \param potparm ...
!> \param error ...
!> \author Teodoro Laino [Teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_single_release ( potparm, error )
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'pair_potential_single_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(ASSOCIATED(potparm),cp_failure_level,routineP,error,failure)
       CALL spline_data_p_release( potparm% pair_spline_data ,error=error)
       CALL spline_factor_release( potparm% spl_f, error=error )
       DO i = 1, SIZE(potparm%type)
          CALL pair_potential_ipbv_release     ( potparm%set(i)%ipbv, error=error )
          CALL pair_potential_lj_release       ( potparm%set(i)%lj, error=error )
          CALL pair_potential_bmhft_release    ( potparm%set(i)%ft, error=error )
          CALL pair_potential_bmhftd_release   ( potparm%set(i)%ftd, error=error )
          CALL pair_potential_williams_release ( potparm%set(i)%willis, error=error )
          CALL pair_potential_goodwin_release  ( potparm%set(i)%goodwin, error=error )
          CALL pair_potential_eam_release      ( potparm%set(i)%eam, error=error )
          CALL pair_potential_quip_release      ( potparm%set(i)%quip, error=error )
          CALL pair_potential_buck4r_release   ( potparm%set(i)%buck4r, error=error )
          CALL pair_potential_buckmo_release   ( potparm%set(i)%buckmo, error=error )
          CALL pair_potential_gp_release       ( potparm%set(i)%gp, error=error )
          CALL pair_potential_tersoff_release  ( potparm%set(i)%tersoff, error=error )
          CALL pair_potential_siepmann_release ( potparm%set(i)%siepmann, error=error )
       END DO
       DEALLOCATE(potparm%type, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(potparm%set, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(potparm, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    NULLIFY ( potparm )
  END SUBROUTINE pair_potential_single_release

! *****************************************************************************
!> \brief Data-structure that constains potential parameters
!> \param potparm ...
!> \param nkinds ...
!> \param error ...
!> \author unknown
! *****************************************************************************
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

! *****************************************************************************
!> \brief Release Data-structure that constains potential parameters
!> \param potparm ...
!> \param error ...
!> \par History
!>      Teodoro Laino [Teo] 11.2005 : Reorganizing the structures to optimize
!>                                    memory management
!> \author unknown
! *****************************************************************************
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

! *****************************************************************************
!> \brief Data-structure that constains potential parameters
!> \param potparm ...
!> \param ndim ...
!> \param ub ...
!> \param lb ...
!> \param error ...
!> \author unknown
! *****************************************************************************
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
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (PRESENT(lb).OR.PRESENT(ub)) THEN
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
       ELSE IF (PRESENT(lb).AND.PRESENT(ub)) THEN
          loc_lb = lb
          loc_ub = ub
          ALLOCATE ( potparm%pot( loc_lb : loc_ub ), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (PRESENT(ndim)) THEN
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
       ELSE
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
       DO i = loc_lb, loc_ub
          NULLIFY ( potparm%pot(i)%pot)
          CALL pair_potential_single_create(potparm%pot(i)%pot, error=error)
       END DO
    END IF
  END SUBROUTINE pair_potential_p_create

! *****************************************************************************
!> \brief Release Data-structure that constains potential parameters
!> \param potparm ...
!> \param error ...
!> \par History
!>      Teodoro Laino [Teo] 11.2005 : Reorganizing the structures to optimize
!>                                    memory management
!> \author unknown
! *****************************************************************************
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

! *****************************************************************************
!> \brief Copy structures between two pair_potential_p_type
!> \param source ...
!> \param dest ...
!> \param istart ...
!> \param iend ...
!> \param error ...
!> \author Teodoro Laino [Teo] 11.2005
! *****************************************************************************
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

! *****************************************************************************
!> \brief Cleans the potential parameter type
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lj ...
!> \param lj_charmm ...
!> \param williams ...
!> \param goodwin ...
!> \param eam ...
!> \param quip ...
!> \param bmhft ...
!> \param bmhftd ...
!> \param ipbv ...
!> \param buck4r ...
!> \param buckmo ...
!> \param gp ...
!> \param tersoff ...
!> \param siepmann ...
!> \param error ...
!> \author Teodoro Laino [Teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_reallocate(p,lb1_new,ub1_new,lj,lj_charmm,williams,&
                           goodwin,eam,quip,bmhft,bmhftd,ipbv,buck4r,buckmo,gp,tersoff,&
                           siepmann,error)
    TYPE(pair_potential_p_type), POINTER     :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new
    LOGICAL, INTENT(IN), OPTIONAL :: lj, lj_charmm, williams, goodwin, eam, &
      quip, bmhft, bmhftd, ipbv, buck4r, buckmo, gp, tersoff, siepmann
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_reallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ipot, lb1_old, std_dim, &
                                                ub1_old
    LOGICAL :: check, failure, lbmhft, lbmhftd, lbuck4r, lbuckmo, leam, &
      lgoodwin, lgp, lipbv, llj, llj_charmm, lquip, lsiepmann, ltersoff, &
      lwilliams
    TYPE(pair_potential_p_type), POINTER     :: work

    NULLIFY(work)
    failure = .FALSE.
    ipot = 0
    llj       = .FALSE. ; IF (PRESENT(lj))        llj       = lj
    llj_charmm= .FALSE. ; IF (PRESENT(lj_charmm)) llj_charmm= lj_charmm
    lwilliams = .FALSE. ; IF (PRESENT(williams))  lwilliams = williams
    lgoodwin  = .FALSE. ; IF (PRESENT(goodwin))   lgoodwin  = goodwin
    leam      = .FALSE. ; IF (PRESENT(eam))       leam      = eam
    lquip      = .FALSE. ; IF (PRESENT(quip))       lquip      = quip
    lbmhft    = .FALSE. ; IF (PRESENT(bmhft))     lbmhft    = bmhft
    lbmhftd   = .FALSE. ; IF (PRESENT(bmhftd))    lbmhftd   = bmhftd
    lipbv     = .FALSE. ; IF (PRESENT(ipbv))      lipbv     = ipbv
    lbuck4r   = .FALSE. ; IF (PRESENT(buck4r))    lbuck4r   = buck4r
    lbuckmo   = .FALSE. ; IF (PRESENT(buckmo))    lbuckmo   = buckmo
    lgp       = .FALSE. ; IF (PRESENT(gp))        lgp       = gp
    ltersoff  = .FALSE. ; IF (PRESENT(tersoff))   ltersoff  = tersoff
    lsiepmann = .FALSE. ; IF (PRESENT(siepmann))  lsiepmann = siepmann

    IF (llj) THEN
       ipot  = lj_type
       check = .NOT.(llj_charmm.OR.lwilliams.OR.lgoodwin.OR.leam.OR.lquip.OR.lbmhft&
                     .OR.lbmhftd.OR.lipbv.OR.lbuck4r.OR.lbuckmo.OR.lgp.OR.ltersoff&
                     .OR.lsiepmann)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (llj_charmm) THEN
       ipot  = lj_charmm_type
       check = .NOT.(llj.OR.lwilliams.OR.lgoodwin.OR.leam.OR.lquip.OR.lbmhft.OR.lbmhftd.OR.lipbv&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp.OR.ltersoff.OR.lsiepmann)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lwilliams) THEN
       ipot  = wl_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lquip.OR.lbmhft.OR.lbmhftd.OR.lipbv&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp.OR.ltersoff.OR.lsiepmann)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lgoodwin) THEN
       ipot  = gw_type
       check = .NOT.(llj.OR.llj_charmm.OR.lwilliams.OR.leam.OR.lquip.OR.lbmhft.OR.lbmhftd.OR.lipbv&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp.OR.ltersoff.OR.lsiepmann)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (leam) THEN
       ipot  = ea_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.lwilliams.OR.lquip.OR.lbmhft.OR.lbmhftd.OR.lipbv&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp.OR.ltersoff.OR.lsiepmann)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lquip) THEN
       ipot  = quip_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.lwilliams.OR.leam.OR.lbmhft.OR.lbmhftd.OR.lipbv&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp.OR.ltersoff.OR.lsiepmann)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lbmhft) THEN
       ipot  = ft_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lquip.OR.lwilliams.OR.lbmhftd.OR.lipbv&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp.OR.ltersoff.OR.lsiepmann)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lbmhftd) THEN
       ipot  = ftd_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lquip.OR.lwilliams.OR.lbmhft.OR.lipbv&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp.OR.ltersoff.OR.lsiepmann)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lipbv) THEN
       ipot  = ip_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lquip.OR.lwilliams.OR.lbmhft.OR.lbmhftd&
                     .OR.lbuck4r.OR.lbuckmo.OR.lgp.OR.ltersoff.OR.lsiepmann)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lbuck4r) THEN
       ipot  = b4_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lquip.OR.lwilliams.OR.lbmhft.OR.lbmhftd&
                     .OR.lipbv.OR.lbuckmo.OR.lgp.OR.ltersoff.OR.lsiepmann)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lbuckmo) THEN
       ipot  = bm_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lquip.OR.lwilliams.OR.lbmhft.OR.lbmhftd&
                    .OR.lipbv.OR.lbuck4r.OR.lgp.OR.ltersoff.OR.lsiepmann)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (ltersoff) THEN
       ipot  = tersoff_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lquip.OR.lwilliams.OR.lbmhft.OR.lbmhftd&
                    .OR.lipbv.OR.lbuck4r.OR.lbuckmo.OR.lgp.OR.lsiepmann)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lsiepmann) THEN
       ipot  = siepmann_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lquip.OR.lwilliams.OR.lbmhft.OR.lbmhftd&
                    .OR.lipbv.OR.lbuck4r.OR.lbuckmo.OR.lgp.OR.ltersoff)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    IF (lgp) THEN
       ipot  = gp_type
       check = .NOT.(llj.OR.llj_charmm.OR.lgoodwin.OR.leam.OR.lwilliams.OR.lbmhft.OR.lbmhftd&
                    .OR.lipbv.OR.lbuck4r.OR.lbuckmo.OR.ltersoff.OR.lsiepmann)
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
       std_dim = 1
       DO i = ub1_old+1, ub1_new
          check = (SIZE(p%pot(i)%pot%type)==std_dim).AND.(SIZE(p%pot(i)%pot%type)==std_dim)
          CPPostcondition(check,cp_failure_level,routineP,error,failure)
          p % pot ( i ) % pot % type  = nn_type
          p % pot ( i ) % pot % shell_type  = nosh_nosh
          p % pot ( i ) % pot % undef = .TRUE.
          p % pot ( i ) % pot % no_mb = .FALSE.
          p % pot ( i ) % pot % no_pp = .FALSE.
          p % pot ( i ) % pot % at1   = 'NULL'
          p % pot ( i ) % pot % at2   = 'NULL'
          p % pot ( i ) % pot % set(std_dim) % rmin = not_initialized
          p % pot ( i ) % pot % set(std_dim) % rmax = not_initialized
          SELECT CASE (ipot)
          CASE (lj_type,lj_charmm_type)
             CALL pair_potential_lj_create(p%pot(i)%pot%set(std_dim)%lj, error=error)
          CASE (wl_type)
             CALL pair_potential_williams_create(p%pot(i)%pot%set(std_dim)%willis, error=error)
          CASE (gw_type)
             CALL pair_potential_goodwin_create(p%pot(i)%pot%set(std_dim)%goodwin, error=error)
          CASE (ea_type)
             CALL pair_potential_eam_create(p%pot(i)%pot%set(std_dim)%eam, error=error)
          CASE (quip_type)
             CALL pair_potential_quip_create(p%pot(i)%pot%set(std_dim)%quip, error=error)
          CASE (ft_type)
             CALL pair_potential_bmhft_create(p%pot(i)%pot%set(std_dim)%ft, error=error)
          CASE (ftd_type)
             CALL pair_potential_bmhftd_create(p%pot(i)%pot%set(std_dim)%ftd, error=error)
          CASE (ip_type)
             CALL pair_potential_ipbv_create(p%pot(i)%pot%set(std_dim)%ipbv, error=error)
          CASE (b4_type)
             CALL pair_potential_buck4r_create(p%pot(i)%pot%set(std_dim)%buck4r, error=error)
          CASE (bm_type)
             CALL pair_potential_buckmo_create(p%pot(i)%pot%set(std_dim)%buckmo, error=error)
          CASE (gp_type)
             CALL pair_potential_gp_create(p%pot(i)%pot%set(std_dim)%gp, error=error)
          CASE (tersoff_type)
             CALL pair_potential_tersoff_create(p%pot(i)%pot%set(std_dim)%tersoff, error=error)
          CASE (siepmann_type)
             CALL pair_potential_siepmann_create(p%pot(i)%pot%set(std_dim)%siepmann, error=error)
          END SELECT
          NULLIFY ( p%pot(i)%pot%spl_f  )
          NULLIFY ( p%pot(i)%pot%pair_spline_data )
       END DO

       IF (ASSOCIATED(work)) CALL pair_potential_p_release(work, error=error)
    END IF
  END SUBROUTINE pair_potential_reallocate

! *****************************************************************************
!> \brief Creates the generic potential type
!> \param gp ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_gp_create(gp, error)
    TYPE(gp_pot_type), POINTER               :: gp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_gp_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(gp),cp_failure_level,routineP,error,failure)
       ALLOCATE(gp, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY (gp%parameters)
       NULLIFY (gp%values)
       CALL pair_potential_gp_clean(gp, error=error)
    END IF
  END SUBROUTINE pair_potential_gp_create

! *****************************************************************************
!> \brief Copy two generic potential type
!> \param gp_source ...
!> \param gp_dest ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_gp_copy(gp_source, gp_dest, error)
    TYPE(gp_pot_type), POINTER               :: gp_source, gp_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_gp_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idim, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(gp_source)) RETURN
       IF (ASSOCIATED(gp_dest)) CALL pair_potential_gp_release(gp_dest, error=error)
       CALL pair_potential_gp_create(gp_dest, error=error)
       gp_dest%myid      = gp_source%myid
       gp_dest%potential = gp_source%potential
       gp_dest%variables = gp_source%variables
       IF (ASSOCIATED(gp_source%parameters)) THEN
          idim = SIZE(gp_source%parameters)
          ALLOCATE(gp_dest%parameters(idim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          gp_dest%parameters = gp_source%parameters
       END IF
       IF (ASSOCIATED(gp_source%values)) THEN
          idim = SIZE(gp_source%values)
          ALLOCATE(gp_dest%values(idim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          gp_dest%values = gp_source%values
       END IF
    END IF
  END SUBROUTINE pair_potential_gp_copy

! *****************************************************************************
!> \brief Cleans the generic potential type
!> \param gp ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_gp_clean(gp, error)
    TYPE(gp_pot_type), POINTER               :: gp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_gp_clean', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(gp)) RETURN
       gp%myid      = 0
       gp%potential = ""
       gp%variables = ""
       IF (ASSOCIATED(gp%values)) THEN
          DEALLOCATE(gp%values,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(gp%parameters)) THEN
          DEALLOCATE(gp%parameters,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE pair_potential_gp_clean

! *****************************************************************************
!> \brief Destroys the generic potential type
!> \param gp ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_gp_release(gp, error)
    TYPE(gp_pot_type), POINTER               :: gp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_gp_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(gp)) THEN
          IF (ASSOCIATED(gp%parameters)) THEN
             DEALLOCATE(gp%parameters, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(gp%values)) THEN
             DEALLOCATE(gp%values, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          DEALLOCATE(gp, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(gp)
    END IF
  END SUBROUTINE pair_potential_gp_release

! *****************************************************************************
!> \brief Cleans the LJ potential type
!> \param lj ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_lj_create(lj, error)
    TYPE(lj_pot_type), POINTER               :: lj
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_lj_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(lj),cp_failure_level,routineP,error,failure)
       ALLOCATE(lj, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_lj_clean(lj, error=error)
    END IF
  END SUBROUTINE pair_potential_lj_create

! *****************************************************************************
!> \brief Copy two LJ potential type
!> \param lj_source ...
!> \param lj_dest ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_lj_copy(lj_source, lj_dest, error)
    TYPE(lj_pot_type), POINTER               :: lj_source, lj_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_lj_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(lj_source)) RETURN
       IF (ASSOCIATED(lj_dest)) CALL pair_potential_lj_release(lj_dest, error=error)
       CALL pair_potential_lj_create(lj_dest, error=error)
       lj_dest%epsilon = lj_source%epsilon
       lj_dest%sigma6  = lj_source%sigma6
       lj_dest%sigma12 = lj_source%sigma12
    END IF
  END SUBROUTINE pair_potential_lj_copy

! *****************************************************************************
!> \brief Creates the LJ potential type
!> \param lj ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_lj_clean(lj, error)
    TYPE(lj_pot_type), POINTER               :: lj
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_lj_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(lj)) RETURN
       lj%epsilon = 0.0_dp
       lj%sigma6  = 0.0_dp
       lj%sigma12 = 0.0_dp
    END IF
  END SUBROUTINE pair_potential_lj_clean

! *****************************************************************************
!> \brief Destroys the LJ potential type
!> \param lj ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_lj_release(lj, error)
    TYPE(lj_pot_type), POINTER               :: lj
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_lj_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(lj)) THEN
          DEALLOCATE(lj, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(lj)
    END IF
  END SUBROUTINE pair_potential_lj_release

! *****************************************************************************
!> \brief Creates the WILLIAMS potential type
!> \param willis ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_williams_create(willis, error)
    TYPE(williams_pot_type), POINTER         :: willis
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_williams_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(willis),cp_failure_level,routineP,error,failure)
       ALLOCATE(willis, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_williams_clean(willis, error=error)
    END IF
  END SUBROUTINE pair_potential_williams_create

! *****************************************************************************
!> \brief Copy two WILLIAMS potential type
!> \param willis_source ...
!> \param willis_dest ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_williams_copy(willis_source, willis_dest, error)
    TYPE(williams_pot_type), POINTER         :: willis_source, willis_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_williams_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(willis_source)) RETURN
       IF (ASSOCIATED(willis_dest)) CALL pair_potential_williams_release(willis_dest, error=error)
       CALL pair_potential_williams_create(willis_dest, error=error)
       willis_dest%a = willis_source%a
       willis_dest%b = willis_source%b
       willis_dest%c = willis_source%c
    END IF
  END SUBROUTINE pair_potential_williams_copy

! *****************************************************************************
!> \brief Creates the WILLIAMS potential type
!> \param willis ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_williams_clean(willis, error)
    TYPE(williams_pot_type), POINTER         :: willis
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_williams_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(willis)) RETURN
       willis%a = 0.0_dp
       willis%b = 0.0_dp
       willis%c = 0.0_dp
    END IF
  END SUBROUTINE pair_potential_williams_clean

! *****************************************************************************
!> \brief Destroys the WILLIAMS potential type
!> \param willis ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_williams_release(willis, error)
    TYPE(williams_pot_type), POINTER         :: willis
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_williams_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(willis)) THEN
          DEALLOCATE(willis, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(willis)
    END IF
  END SUBROUTINE pair_potential_williams_release

! *****************************************************************************
!> \brief Creates the GOODWIN potential type
!> \param goodwin ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_goodwin_create(goodwin, error)
    TYPE(goodwin_pot_type), POINTER          :: goodwin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_goodwin_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(goodwin),cp_failure_level,routineP,error,failure)
       ALLOCATE(goodwin, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_goodwin_clean(goodwin, error=error)
    END IF
  END SUBROUTINE pair_potential_goodwin_create

! *****************************************************************************
!> \brief Copy two GOODWIN potential type
!> \param goodwin_source ...
!> \param goodwin_dest ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_goodwin_copy(goodwin_source, goodwin_dest, error)
    TYPE(goodwin_pot_type), POINTER          :: goodwin_source, goodwin_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_goodwin_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(goodwin_source)) RETURN
       IF (ASSOCIATED(goodwin_dest)) CALL pair_potential_goodwin_release(goodwin_dest, error=error)
       CALL pair_potential_goodwin_create(goodwin_dest, error=error)
       goodwin_dest%vr0 = goodwin_source%vr0
       goodwin_dest%d   = goodwin_source%d
       goodwin_dest%dc  = goodwin_source%dc
       goodwin_dest%m   = goodwin_source%m
       goodwin_dest%mc  = goodwin_source%mc
    END IF
  END SUBROUTINE pair_potential_goodwin_copy

! *****************************************************************************
!> \brief Creates the GOODWIN potential type
!> \param goodwin ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_goodwin_clean(goodwin, error)
    TYPE(goodwin_pot_type), POINTER          :: goodwin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_goodwin_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(goodwin)) RETURN
       goodwin%vr0 = 0.0_dp
       goodwin%d   = 0.0_dp
       goodwin%dc  = 0.0_dp
       goodwin%m   = 0.0_dp
       goodwin%mc  = 0.0_dp
    END IF
  END SUBROUTINE pair_potential_goodwin_clean

! *****************************************************************************
!> \brief Destroys the GOODWIN potential type
!> \param goodwin ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_goodwin_release(goodwin, error)
    TYPE(goodwin_pot_type), POINTER          :: goodwin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_goodwin_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(goodwin)) THEN
          DEALLOCATE(goodwin, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(goodwin)
    END IF
  END SUBROUTINE pair_potential_goodwin_release

! *****************************************************************************
!> \brief Creates the EAM potential type
!> \param eam ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_eam_create(eam, error)
    TYPE(eam_pot_type), POINTER              :: eam
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_eam_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(eam),cp_failure_level,routineP,error,failure)
       ALLOCATE(eam, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY(eam%rho, eam%phi, eam%frho, eam%rhoval, eam%rval,&
               eam%rhop, eam%phip, eam%frhop)
       CALL pair_potential_eam_clean(eam, error=error)
    END IF
  END SUBROUTINE pair_potential_eam_create

! *****************************************************************************
!> \brief Copy two EAM potential type
!> \param eam_source ...
!> \param eam_dest ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_eam_copy(eam_source, eam_dest, error)
    TYPE(eam_pot_type), POINTER              :: eam_source, eam_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_eam_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(eam_source)) RETURN
       IF (ASSOCIATED(eam_dest)) CALL pair_potential_eam_release(eam_dest, error=error)
       CALL pair_potential_eam_create(eam_dest, error=error)
       eam_dest%eam_file_name = eam_source%eam_file_name
       eam_dest%drar          = eam_source%drar
       eam_dest%drhoar        = eam_source%drhoar
       eam_dest%acutal        = eam_source%acutal
       eam_dest%npoints       = eam_source%npoints
       ! Allocate arrays with the proper size
       CALL reallocate(eam_dest%rho,   1, eam_dest%npoints)
       CALL reallocate(eam_dest%rhop,  1, eam_dest%npoints)
       CALL reallocate(eam_dest%phi,   1, eam_dest%npoints)
       CALL reallocate(eam_dest%phip,  1, eam_dest%npoints)
       CALL reallocate(eam_dest%frho,  1, eam_dest%npoints)
       CALL reallocate(eam_dest%frhop, 1, eam_dest%npoints)
       CALL reallocate(eam_dest%rval,  1, eam_dest%npoints)
       CALL reallocate(eam_dest%rhoval,1, eam_dest%npoints)
       eam_dest%rho           = eam_source%rho
       eam_dest%phi           = eam_source%phi
       eam_dest%frho          = eam_source%frho
       eam_dest%rhoval        = eam_source%rhoval
       eam_dest%rval          = eam_source%rval
       eam_dest%rhop          = eam_source%rhop
       eam_dest%phip          = eam_source%phip
       eam_dest%frhop         = eam_source%frhop
    END IF
  END SUBROUTINE pair_potential_eam_copy

! *****************************************************************************
!> \brief Creates the EAM potential type
!> \param eam ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_eam_clean(eam, error)
    TYPE(eam_pot_type), POINTER              :: eam
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_eam_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(eam)) RETURN
       eam%eam_file_name = 'NULL'
       eam%drar          = 0.0_dp
       eam%drhoar        = 0.0_dp
       eam%acutal        = 0.0_dp
       eam%npoints       = 0
       CALL reallocate(eam%rho,   1, eam%npoints)
       CALL reallocate(eam%rhop,  1, eam%npoints)
       CALL reallocate(eam%phi,   1, eam%npoints)
       CALL reallocate(eam%phip,  1, eam%npoints)
       CALL reallocate(eam%frho,  1, eam%npoints)
       CALL reallocate(eam%frhop, 1, eam%npoints)
       CALL reallocate(eam%rval,  1, eam%npoints)
       CALL reallocate(eam%rhoval,1, eam%npoints)
    END IF
  END SUBROUTINE pair_potential_eam_clean

! *****************************************************************************
!> \brief Destroys the EAM potential type
!> \param eam ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_eam_release(eam, error)
    TYPE(eam_pot_type), POINTER              :: eam
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_eam_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(eam)) THEN
          IF (ASSOCIATED(eam%rho)) THEN
             DEALLOCATE(eam%rho, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(eam%rhop)) THEN
             DEALLOCATE(eam%rhop, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(eam%phi)) THEN
             DEALLOCATE(eam%phi, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(eam%phip)) THEN
             DEALLOCATE(eam%phip, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(eam%frho)) THEN
             DEALLOCATE(eam%frho, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(eam%frhop)) THEN
             DEALLOCATE(eam%frhop, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(eam%rval)) THEN
             DEALLOCATE(eam%rval, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(eam%rhoval)) THEN
             DEALLOCATE(eam%rhoval, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          DEALLOCATE(eam, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE pair_potential_eam_release

! *****************************************************************************
!> \brief Creates the QUIP potential type
!> \param quip ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_quip_create(quip, error)
    TYPE(quip_pot_type), POINTER             :: quip
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_quip_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(quip),cp_failure_level,routineP,error,failure)
       ALLOCATE(quip, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       quip%quip_file_name=""
       quip%init_args=""
       quip%calc_args=""
       CALL pair_potential_quip_clean(quip, error=error)
    END IF
  END SUBROUTINE pair_potential_quip_create

! *****************************************************************************
!> \brief Copy two QUIP potential type
!> \param quip_source ...
!> \param quip_dest ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_quip_copy(quip_source, quip_dest, error)
    TYPE(quip_pot_type), POINTER             :: quip_source, quip_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_quip_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(quip_source)) RETURN
       IF (ASSOCIATED(quip_dest)) CALL pair_potential_quip_release(quip_dest, error=error)
       CALL pair_potential_quip_create(quip_dest, error=error)
       quip_dest%quip_file_name = quip_source%quip_file_name
       quip_dest%init_args = quip_source%init_args
       quip_dest%calc_args = quip_source%calc_args
    END IF
  END SUBROUTINE pair_potential_quip_copy

! *****************************************************************************
!> \brief Creates the QUIP potential type
!> \param quip ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_quip_clean(quip, error)
    TYPE(quip_pot_type), POINTER             :: quip
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_quip_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(quip)) RETURN
       quip%quip_file_name = 'NULL'
       quip%init_args = ''
       quip%calc_args = ''
    END IF
  END SUBROUTINE pair_potential_quip_clean

! *****************************************************************************
!> \brief Destroys the QUIP potential type
!> \param quip ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_quip_release(quip, error)
    TYPE(quip_pot_type), POINTER             :: quip
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_quip_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(quip)) THEN
          DEALLOCATE(quip, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE pair_potential_quip_release

! *****************************************************************************
!> \brief Creates the BMHFT (TOSI-FUMI) potential type
!> \param ft ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_bmhft_create(ft, error)
    TYPE(ft_pot_type), POINTER               :: ft
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_bmhft_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(ft),cp_failure_level,routineP,error,failure)
       ALLOCATE(ft, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_bmhft_clean(ft, error=error)
    END IF
  END SUBROUTINE pair_potential_bmhft_create

! *****************************************************************************
!> \brief Copy two BMHFT (TOSI-FUMI) potential type
!> \param ft_source ...
!> \param ft_dest ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_bmhft_copy(ft_source, ft_dest, error)
    TYPE(ft_pot_type), POINTER               :: ft_source, ft_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_bmhft_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(ft_source)) RETURN
       IF (ASSOCIATED(ft_dest)) CALL pair_potential_bmhft_release(ft_dest, error=error)
       CALL pair_potential_bmhft_create(ft_dest, error=error)
       ft_dest%A = ft_source%A
       ft_dest%B = ft_source%B
       ft_dest%C = ft_source%C
       ft_dest%D = ft_source%D
    END IF
  END SUBROUTINE pair_potential_bmhft_copy

! *****************************************************************************
!> \brief Creates the BMHFT (TOSI-FUMI) potential type
!> \param ft ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_bmhft_clean(ft, error)
    TYPE(ft_pot_type), POINTER               :: ft
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_bmhft_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(ft)) RETURN
       ft%A = 0.0_dp
       ft%B = 0.0_dp
       ft%C = 0.0_dp
       ft%D = 0.0_dp
    END IF
  END SUBROUTINE pair_potential_bmhft_clean

! *****************************************************************************
!> \brief Destroys the BMHFT potential type
!> \param ft ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_bmhft_release(ft, error)
    TYPE(ft_pot_type), POINTER               :: ft
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_bmhft_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(ft)) THEN
          DEALLOCATE(ft, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(ft)
    END IF
  END SUBROUTINE pair_potential_bmhft_release

! *****************************************************************************
!> \brief Creates the BMHFTD (damped TOSI-FUMI) potential type
!> \param ftd ...
!> \param error ...
!> \author Mathieu Salanne 05.2010
! *****************************************************************************
  SUBROUTINE pair_potential_bmhftd_create(ftd, error)
    TYPE(ftd_pot_type), POINTER              :: ftd
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_bmhftd_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(ftd),cp_failure_level,routineP,error,failure)
       ALLOCATE(ftd, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_bmhftd_clean(ftd, error=error)
    END IF
  END SUBROUTINE pair_potential_bmhftd_create

! *****************************************************************************
!> \brief Copy two BMHFTD (Damped TOSI-FUMI) potential type
!> \param ftd_source ...
!> \param ftd_dest ...
!> \param error ...
!> \author Mathieu Salanne 05.2010
! *****************************************************************************
  SUBROUTINE pair_potential_bmhftd_copy(ftd_source, ftd_dest, error)
    TYPE(ftd_pot_type), POINTER              :: ftd_source, ftd_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_bmhftd_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(ftd_source)) RETURN
       IF (ASSOCIATED(ftd_dest)) CALL pair_potential_bmhftd_release(ftd_dest, error=error)
       CALL pair_potential_bmhftd_create(ftd_dest, error=error)
       ftd_dest%A = ftd_source%A
       ftd_dest%B = ftd_source%B
       ftd_dest%C = ftd_source%C
       ftd_dest%D = ftd_source%D
       ftd_dest%BD = ftd_source%BD
    END IF
  END SUBROUTINE pair_potential_bmhftd_copy

! *****************************************************************************
!> \brief Cleans the BMHFTD (damped TOSI-FUMI) potential type
!> \param ftd ...
!> \param error ...
!> \author Mathieu Salanne
! *****************************************************************************
  SUBROUTINE pair_potential_bmhftd_clean(ftd, error)
    TYPE(ftd_pot_type), POINTER              :: ftd
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_bmhftd_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(ftd)) RETURN
       ftd%A = 0.0_dp
       ftd%B = 0.0_dp
       ftd%C = 0.0_dp
       ftd%D = 0.0_dp
       ftd%BD = 0.0_dp
    END IF
  END SUBROUTINE pair_potential_bmhftd_clean

! *****************************************************************************
!> \brief Destroys the BMHFTD potential type
!> \param ftd ...
!> \param error ...
!> \author Mathieu Salanne 05.2010
! *****************************************************************************
  SUBROUTINE pair_potential_bmhftd_release(ftd, error)
    TYPE(ftd_pot_type), POINTER              :: ftd
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_bmhftd_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(ftd)) THEN
          DEALLOCATE(ftd, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(ftd)
    END IF
  END SUBROUTINE pair_potential_bmhftd_release

! *****************************************************************************
!> \brief Creates the IPBV potential type
!> \param ipbv ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_ipbv_create(ipbv, error)
    TYPE(ipbv_pot_type), POINTER             :: ipbv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_ipbv_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(ipbv),cp_failure_level,routineP,error,failure)
       ALLOCATE(ipbv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_ipbv_clean(ipbv, error=error)
    END IF
  END SUBROUTINE pair_potential_ipbv_create

! *****************************************************************************
!> \brief Copy two IPBV potential type
!> \param ipbv_source ...
!> \param ipbv_dest ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_ipbv_copy(ipbv_source, ipbv_dest, error)
    TYPE(ipbv_pot_type), POINTER             :: ipbv_source, ipbv_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_ipbv_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(ipbv_source)) RETURN
       IF (ASSOCIATED(ipbv_dest)) CALL pair_potential_ipbv_release(ipbv_dest, error=error)
       CALL pair_potential_ipbv_create(ipbv_dest, error=error)
       ipbv_dest%a     = ipbv_source%a
       ipbv_dest%rcore = ipbv_source%rcore
       ipbv_dest%b     = ipbv_source%b
       ipbv_dest%m     = ipbv_source%m
    END IF
  END SUBROUTINE pair_potential_ipbv_copy

! *****************************************************************************
!> \brief Creates the IPBV potential type
!> \param ipbv ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_ipbv_clean(ipbv, error)
    TYPE(ipbv_pot_type), POINTER             :: ipbv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_ipbv_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(ipbv)) RETURN
       ipbv%a     = 0.0_dp
       ipbv%rcore = 0.0_dp
       ipbv%b     = 0.0_dp
       ipbv%m     = 0.0_dp
    END IF
  END SUBROUTINE pair_potential_ipbv_clean

! *****************************************************************************
!> \brief Destroys the IPBV potential type
!> \param ipbv ...
!> \param error ...
!> \author Teodoro Laino [teo] 11.2005
! *****************************************************************************
  SUBROUTINE pair_potential_ipbv_release(ipbv, error)
    TYPE(ipbv_pot_type), POINTER             :: ipbv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_ipbv_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(ipbv)) THEN
          DEALLOCATE(ipbv, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(ipbv)
    END IF
  END SUBROUTINE pair_potential_ipbv_release

! *****************************************************************************
!> \brief Creates the Buckingham 4 ranges  potential type
!> \param buck4r ...
!> \param error ...
!> \author MI 10.2006
! *****************************************************************************
  SUBROUTINE pair_potential_buck4r_create(buck4r, error)
    TYPE(buck4ran_pot_type), POINTER         :: buck4r
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_buck4r_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(buck4r),cp_failure_level,routineP,error,failure)
       ALLOCATE(buck4r, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_buck4r_clean(buck4r, error=error)
    END IF
  END SUBROUTINE pair_potential_buck4r_create

! *****************************************************************************
!> \brief Copy two Buckingham 4 ranges  potential type
!> \param buck4r_source ...
!> \param buck4r_dest ...
!> \param error ...
!> \author MI 10.2006
! *****************************************************************************
  SUBROUTINE pair_potential_buck4r_copy(buck4r_source, buck4r_dest, error)
    TYPE(buck4ran_pot_type), POINTER         :: buck4r_source, buck4r_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_buck4r_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(buck4r_source)) RETURN
       IF (ASSOCIATED(buck4r_dest)) CALL pair_potential_buck4r_release(buck4r_dest, error=error)
       CALL pair_potential_buck4r_create(buck4r_dest, error=error)
       buck4r_dest%a      = buck4r_source%a
       buck4r_dest%b      = buck4r_source%b
       buck4r_dest%c      = buck4r_source%c
       buck4r_dest%r1     = buck4r_source%r1
       buck4r_dest%r2     = buck4r_source%r2
       buck4r_dest%r3     = buck4r_source%r3
       buck4r_dest%poly1  = buck4r_source%poly1
       buck4r_dest%poly2  = buck4r_source%poly2
       buck4r_dest%npoly1 = buck4r_source%npoly1
       buck4r_dest%npoly2 = buck4r_source%npoly2
    END IF
  END SUBROUTINE pair_potential_buck4r_copy

! *****************************************************************************
!> \brief Creates the Buckingham 4 ranges  potential type
!> \param buck4r ...
!> \param error ...
!> \author MI 10.2006
! *****************************************************************************
  SUBROUTINE pair_potential_buck4r_clean(buck4r, error)
    TYPE(buck4ran_pot_type), POINTER         :: buck4r
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_buck4r_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(buck4r)) RETURN
       buck4r%a  = 0.0_dp
       buck4r%b  = 0.0_dp
       buck4r%c  = 0.0_dp
       buck4r%r1 = 0.0_dp
       buck4r%r2 = 0.0_dp
       buck4r%r3 = 0.0_dp
       buck4r%poly1  = 0.0_dp
       buck4r%npoly1 = 0
       buck4r%poly2  = 0.0_dp
       buck4r%npoly2 = 0
    END IF
  END SUBROUTINE pair_potential_buck4r_clean

! *****************************************************************************
!> \brief Destroys the Buckingham 4 ranges potential type
!> \param buck4r ...
!> \param error ...
!> \author MI  10.2006
! *****************************************************************************
  SUBROUTINE pair_potential_buck4r_release(buck4r, error)
    TYPE(buck4ran_pot_type), POINTER         :: buck4r
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_buck4r_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(buck4r)) THEN
          DEALLOCATE(buck4r, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(buck4r)
    END IF
  END SUBROUTINE pair_potential_buck4r_release

! *****************************************************************************
!> \brief Creates the Buckingham plus Morse potential type
!> \param buckmo ...
!> \param error ...
!> \author MI 10.2006
! *****************************************************************************
  SUBROUTINE pair_potential_buckmo_create(buckmo, error)
    TYPE(buckmorse_pot_type), POINTER        :: buckmo
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_buckmo_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(buckmo),cp_failure_level,routineP,error,failure)
       ALLOCATE(buckmo, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_buckmo_clean(buckmo, error=error)
    END IF
  END SUBROUTINE pair_potential_buckmo_create

! *****************************************************************************
!> \brief Copy two Buckingham plus Morse  potential type
!> \param buckmo_source ...
!> \param buckmo_dest ...
!> \param error ...
!> \author MI 10.2006
! *****************************************************************************
  SUBROUTINE pair_potential_buckmo_copy(buckmo_source, buckmo_dest, error)
    TYPE(buckmorse_pot_type), POINTER        :: buckmo_source, buckmo_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_buckmo_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(buckmo_source)) RETURN
       IF (ASSOCIATED(buckmo_dest)) CALL pair_potential_buckmo_release(buckmo_dest, error=error)
       CALL pair_potential_buckmo_create(buckmo_dest, error=error)
       buckmo_dest%f0   = buckmo_source%f0
       buckmo_dest%a1   = buckmo_source%a1
       buckmo_dest%a2   = buckmo_source%a2
       buckmo_dest%b1   = buckmo_source%b1
       buckmo_dest%b2   = buckmo_source%b2
       buckmo_dest%c    = buckmo_source%c
       buckmo_dest%d    = buckmo_source%d
       buckmo_dest%r0   = buckmo_source%r0
       buckmo_dest%beta = buckmo_source%beta
    END IF
  END SUBROUTINE pair_potential_buckmo_copy

! *****************************************************************************
!> \brief Creates the Buckingham plus Morse  potential type
!> \param buckmo ...
!> \param error ...
!> \author MI 10.2006
! *****************************************************************************
  SUBROUTINE pair_potential_buckmo_clean(buckmo, error)
    TYPE(buckmorse_pot_type), POINTER        :: buckmo
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_buckmo_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(buckmo)) RETURN
       buckmo%f0 = 0.0_dp
       buckmo%a1 = 0.0_dp
       buckmo%a2 = 0.0_dp
       buckmo%b1 = 0.0_dp
       buckmo%b2 = 0.0_dp
       buckmo%c = 0.0_dp
       buckmo%d  = 0.0_dp
       buckmo%r0 = 0.0_dp
       buckmo%beta = 0.0_dp
    END IF
  END SUBROUTINE pair_potential_buckmo_clean

! *****************************************************************************
!> \brief Destroys the Buckingham plus Morse potential type
!> \param buckmo ...
!> \param error ...
!> \author MI  10.2006
! *****************************************************************************
  SUBROUTINE pair_potential_buckmo_release(buckmo, error)
    TYPE(buckmorse_pot_type), POINTER        :: buckmo
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_buckmo_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(buckmo)) THEN
          DEALLOCATE(buckmo, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(buckmo)
    END IF
  END SUBROUTINE pair_potential_buckmo_release

! *****************************************************************************
!> \brief Creates the Tersoff potential type
!>      (Tersoff, J. PRB 39(8), 5566, 1989)
!> \param tersoff ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pair_potential_tersoff_create(tersoff, error)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_tersoff_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(tersoff),cp_failure_level,routineP,error,failure)
       ALLOCATE(tersoff, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_tersoff_clean(tersoff, error=error)
    END IF
  END SUBROUTINE pair_potential_tersoff_create

! *****************************************************************************
!> \brief Copy two Tersoff potential type
!>      (Tersoff, J. PRB 39(8), 5566, 1989)
!> \param tersoff_source ...
!> \param tersoff_dest ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pair_potential_tersoff_copy(tersoff_source, tersoff_dest, error)
    TYPE(tersoff_pot_type), POINTER          :: tersoff_source, tersoff_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_tersoff_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(tersoff_source)) RETURN
       IF (ASSOCIATED(tersoff_dest)) CALL pair_potential_tersoff_release(tersoff_dest, error=error)
       CALL pair_potential_tersoff_create(tersoff_dest, error=error)
       tersoff_dest%A       = tersoff_source%A
       tersoff_dest%B       = tersoff_source%B
       tersoff_dest%lambda1 = tersoff_source%lambda1
       tersoff_dest%lambda2 = tersoff_source%lambda2
       tersoff_dest%alpha   = tersoff_source%alpha
       tersoff_dest%beta    = tersoff_source%beta
       tersoff_dest%n       = tersoff_source%n
       tersoff_dest%c       = tersoff_source%c
       tersoff_dest%d       = tersoff_source%d
       tersoff_dest%h       = tersoff_source%h
       tersoff_dest%lambda3 = tersoff_source%lambda3
       tersoff_dest%bigR    = tersoff_source%bigR
       tersoff_dest%bigD    = tersoff_source%bigD
       tersoff_dest%rcutsq  = tersoff_source%rcutsq
    END IF
  END SUBROUTINE pair_potential_tersoff_copy

! *****************************************************************************
!> \brief Creates the Tersoff potential type
!>      (Tersoff, J. PRB 39(8), 5566, 1989)
!> \param tersoff ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pair_potential_tersoff_clean(tersoff, error)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_tersoff_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(tersoff)) RETURN
       tersoff%A       = 0.0_dp
       tersoff%B       = 0.0_dp
       tersoff%lambda1 = 0.0_dp
       tersoff%lambda2 = 0.0_dp
       tersoff%alpha   = 0.0_dp
       tersoff%beta    = 0.0_dp
       tersoff%n       = 0.0_dp
       tersoff%c       = 0.0_dp
       tersoff%d       = 0.0_dp
       tersoff%h       = 0.0_dp
       tersoff%lambda3 = 0.0_dp
       tersoff%bigR    = 0.0_dp
       tersoff%bigD    = 0.0_dp
       tersoff%rcutsq  = 0.0_dp
    END IF
  END SUBROUTINE pair_potential_tersoff_clean

! *****************************************************************************
!> \brief Destroys the Tersoff
!>      (Tersoff, J. PRB 39(8), 5566, 1989)
!> \param tersoff ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pair_potential_tersoff_release(tersoff, error)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_tersoff_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(tersoff)) THEN
          DEALLOCATE(tersoff, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(tersoff)
    END IF
  END SUBROUTINE pair_potential_tersoff_release

! ***************************************************************************** 
!> \brief Creates the Siepmann-Sprik potential type
!>      (Siepmann and Sprik, J. Chem. Phys. 102(1) 511, 1995)  
!> \param siepmann ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pair_potential_siepmann_create(siepmann, error)
    TYPE(siepmann_pot_type), POINTER         :: siepmann
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_siepmann_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(siepmann),cp_failure_level,routineP,error,failure)
       ALLOCATE(siepmann, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pair_potential_siepmann_clean(siepmann, error=error)
    END IF
  END SUBROUTINE pair_potential_siepmann_create
! *****************************************************************************
!> \brief Copy two Siepmann potential type
!>      (Siepmann and Sprik, J. Chem. Phys. 102(1) 511, 1995)  
!> \param siepmann_source ...
!> \param siepmann_dest ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pair_potential_siepmann_copy(siepmann_source, siepmann_dest, error)
    TYPE(siepmann_pot_type), POINTER         :: siepmann_source, siepmann_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pair_potential_siepmann_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(siepmann_source)) RETURN
       IF (ASSOCIATED(siepmann_dest)) CALL pair_potential_siepmann_release(siepmann_dest, error=error)
       CALL pair_potential_siepmann_create(siepmann_dest, error=error)
       siepmann_dest%B       = siepmann_source%B
       siepmann_dest%D       = siepmann_source%D
       siepmann_dest%E       = siepmann_source%E
       siepmann_dest%F       = siepmann_source%F
       siepmann_dest%alpha   = siepmann_source%alpha
       siepmann_dest%beta    = siepmann_source%beta
       siepmann_dest%rcutsq  = siepmann_source%rcutsq
    END IF
  END SUBROUTINE pair_potential_siepmann_copy

! *****************************************************************************
!> \brief Creates the Siepmann-Sprik potential type
!>      (Siepmann and Sprik, J. Chem. Phys. 102(1) 511, 1995)  
!> \param siepmann ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pair_potential_siepmann_clean(siepmann, error)
    TYPE(siepmann_pot_type), POINTER         :: siepmann
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_siepmann_clean', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF(.NOT.ASSOCIATED(siepmann)) RETURN
       siepmann%B       = 0.0_dp
       siepmann%D       = 0.0_dp
       siepmann%E       = 0.0_dp
       siepmann%F       = 0.0_dp
       siepmann%alpha   = 0.0_dp
       siepmann%beta    = 0.0_dp
       siepmann%rcutsq  = 0.0_dp
    END IF
  END SUBROUTINE pair_potential_siepmann_clean

! *****************************************************************************
!> \brief Destroys the Siepmann-Sprik potential
!>      (Siepmann and Sprik, J. Chem. Phys. 102(1) 511, 1995)  
!> \param siepmann ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pair_potential_siepmann_release(siepmann, error)
    TYPE(siepmann_pot_type), POINTER         :: siepmann
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'pair_potential_siepmann_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(siepmann)) THEN
          DEALLOCATE(siepmann, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       NULLIFY(siepmann)
    END IF
  END SUBROUTINE pair_potential_siepmann_release

END MODULE pair_potential_types

