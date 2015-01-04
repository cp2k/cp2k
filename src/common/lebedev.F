!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Generation of the spherical Lebedev grids.
!>          All Lebedev grids were generated with a precision of at least
!>          33 digits (Mathematica). Thus the Lebedev abscissas and weights
!>          are accurate to 32 digits (quadruple precision).
!> \version 1.0
!>
!> \par Literature
!>          - V. I. Lebedev, Zh. Vychisl. Mat. Mat. Fiz. 15, 48 (1975)
!>          - V. I. Lebedev, Zh. Vychisl. Mat. Mat. Fiz. 16, 293 (1976)
!>          - V. I. Lebedev, Sibirsk Mat. Zh. 18, 132 (1977)
!>          - V. I. Lebedev and A. L. Skorokhodov, Russ. Acad. Sci. Dokl.
!>            Math. 45, 587 (1992)
!>          - V. I. Lebedev, Russ. Acad. Sci. Dokl. Math. 50, 283 (1995)
!> \par Notes
!>          Implemented Lebedev grids:
!>            No.    l    n(nonred.)    n(total)
!>             1     3         1             6
!>             2     5         2            14
!>             3     7         3            26
!>             4     9         3            38
!>             5    11         4            50
!>             6    15         5            86
!>             7    17         6           110
!>             8    19         7           146
!>             9    23         9           194
!>            10    29        12           302
!>            11    35        16           434
!>            12    41        20           590
!>            13    47        25           770
!>            14    53        30           974
!> \par Variables:
!>        - lebedev_grid: Lebedev grids.
!>                        l: Angular momentum quantum number l.
!>                        n: Number of grid points.
!>                        w: Grid point weights.
!>                        r: Abscissas (grid point coordinates)
!>        - max_np      : Maximum number of nonredundant grid points.
!>        - na1         : Number of point for subsystem A1.
!>        - na2         : Number of point for subsystem A2.
!>        - na3         : Number of point for subsystem A3.
!>        - nb          : Number of point for subsystem B.
!>        - nc          : Number of point for subsystem C.
!>        - nd          : Number of point for subsystem D.
!>        - nlg         : Number of implemented Lebedev grids.
!>        - nlgp        : Grid point counter.
!>        - r           : Working storage (grid point coordinates).
!>        - w           : Working storage (grid point weights).
! *****************************************************************************
MODULE lebedev

  
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.f90"

  IMPLICIT NONE

  SAVE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lebedev'

  INTEGER, PARAMETER  :: na1 = 6,&
                         na2 = 12,&
                         na3 = 8,&
                         nb = 24,&
                         nc = 24,&
                         nd = 48,&
                         nlg = 14,&
                         max_np = 36

! *****************************************************************************
  TYPE oh_grid
    INTEGER                           :: l,n
    REAL(KIND = dp), DIMENSION(:), POINTER   :: w => NULL()
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: r => NULL()
  END TYPE oh_grid

  TYPE(oh_grid), DIMENSION(nlg), TARGET :: lebedev_grid

  REAL(KIND = dp) :: one,rs2,rs3,zero
  INTEGER  :: nlgp
  LOGICAL  :: init_lebedev_grids_done = .FALSE.

  REAL(KIND = dp), DIMENSION(max_np/3) :: w
  REAL(KIND = dp), DIMENSION(max_np)   :: r

  PUBLIC :: lebedev_grid

  PUBLIC :: deallocate_lebedev_grids, get_number_of_lebedev_grid, init_lebedev_grids

CONTAINS

! *****************************************************************************
!> \brief  Get the number of the Lebedev grid, which has the requested
!>          angular momentum quantnum number l or size n.
!> \param l ...
!> \param n ...
!> \retval number_of_lebedev_grid ...
!> \date   05.05.99
!> \par Variables
!> \author Matthias Krack
!> \version 1.0
!> \note         - l : Minimal angular momentum quantum number of the requested Lebedev grid.
!> \note         - n : Minimal size of the requested Lebedev grid.
! *****************************************************************************
  FUNCTION get_number_of_lebedev_grid(l,n) RESULT(number_of_lebedev_grid)

    INTEGER, INTENT(IN), OPTIONAL            :: l, n
    INTEGER                                  :: number_of_lebedev_grid

    INTEGER                                  :: i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ll, nn
    INTEGER, DIMENSION(1)                    :: lgnum

    IF (.NOT.init_lebedev_grids_done) CALL init_lebedev_grids

    lgnum(1) = 0
    ALLOCATE(nn(SIZE(lebedev_grid,1)),ll(SIZE(lebedev_grid,1)))
    DO i =1,SIZE(nn,1)
      nn(i) = lebedev_grid(i)%n
      ll(i) = lebedev_grid(i)%l
    END DO
    IF (PRESENT(l)) THEN
!      lgnum(:) = MINLOC(lebedev_grid(:)%n,MASK=(lebedev_grid(:)%l >= l))
      lgnum(:) = MINLOC(nn(:),MASK=(ll(:) >= l))
    ELSE IF (PRESENT(n)) THEN
!      lgnum(:) = MINLOC(lebedev_grid(:)%n,MASK=(lebedev_grid(:)%n >= n))
      lgnum(:) = MINLOC(nn(:),MASK=(nn(:) >= n))
    ELSE
      STOP "*** ERROR in FUNCTION get_number_of_lebedev_grid ***"
    END IF

    IF (lgnum(1) == 0) THEN
      STOP "*** ERROR in FUNCTION get_number_of_lebedev_grid ***"
    ELSE
      number_of_lebedev_grid = lgnum(1)
    END IF
    DEALLOCATE(nn,ll)
  END FUNCTION get_number_of_lebedev_grid

! *****************************************************************************
!> \brief  Load the Lebedev grid points of a subsystem.
!> \param subsystem ...
!> \param lgnum ...
!> \param np ...
!> \date    06.05.99
!> \par Variables
!>       - subsystem: String, which specifies the selected subsystem.
!>       - lgnum    : Number of the current Lebedev grid.
!>       - np       : Number of the nonredundant grid points.
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE load_sub_grid(subsystem,lgnum,np)
    CHARACTER(*), INTENT(IN)                 :: subsystem
    INTEGER, INTENT(IN)                      :: lgnum, np

    CHARACTER(len=*), PARAMETER :: routineN = 'load_sub_grid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    REAL(KIND=dp)                            :: x, y, z

!   *** Check argument values ***

    IF ((lgnum < 1).OR.(lgnum > nlg)) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Argument value #2 (lgnum) is out of range")
    END IF

    IF ((np < 0).OR.(np > max_np)) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Argument value #3 (np) is out of range")
    END IF

    SELECT CASE (subsystem)
    CASE ("A1")
      lebedev_grid(lgnum)%w(1:na1) = w(1)
      lebedev_grid(lgnum)%r(1:3,1:na1) =&
        RESHAPE((/zero,zero, one,&
                  zero,zero,-one,&
                  zero, one,zero,&
                  zero,-one,zero,&
                   one,zero,zero,&
                  -one,zero,zero/),(/3,na1/))
      nlgp = na1
    CASE ("A2")
      lebedev_grid(lgnum)%w(nlgp+1:nlgp+na2) = w(1)
      lebedev_grid(lgnum)%r(1:3,nlgp+1:nlgp+na2) =&
        RESHAPE((/zero, rs2, rs2,&
                  zero, rs2,-rs2,&
                  zero,-rs2, rs2,&
                  zero,-rs2,-rs2,&
                   rs2,zero, rs2,&
                   rs2,zero,-rs2,&
                  -rs2,zero, rs2,&
                  -rs2,zero,-rs2,&
                   rs2, rs2,zero,&
                   rs2,-rs2,zero,&
                  -rs2, rs2,zero,&
                  -rs2,-rs2,zero/),(/3,na2/))
      nlgp = nlgp + na2
    CASE ("A3")
      lebedev_grid(lgnum)%w(nlgp+1:nlgp+na3) = w(1)
      lebedev_grid(lgnum)%r(1:3,nlgp+1:nlgp+na3) =&
        RESHAPE((/ rs3, rs3, rs3,&
                   rs3, rs3,-rs3,&
                   rs3,-rs3, rs3,&
                  -rs3, rs3, rs3,&
                   rs3,-rs3,-rs3,&
                  -rs3, rs3,-rs3,&
                  -rs3,-rs3, rs3,&
                  -rs3,-rs3,-rs3/),(/3,na3/))
      nlgp = nlgp + na3
    CASE ("B")
      DO i=1,np
        x = r(i)
        y = rs2*SQRT(one - x**2)
        lebedev_grid(lgnum)%w(nlgp+nb*(i-1)+1:nlgp+nb*i) = w(i)
        lebedev_grid(lgnum)%r(1:3,nlgp+nb*(i-1)+1:nlgp+nb*i) =&
          RESHAPE((/ x, y, y,&
                     x, y,-y,&
                     x,-y, y,&
                     x,-y,-y,&
                    -x, y, y,&
                    -x, y,-y,&
                    -x,-y, y,&
                    -x,-y,-y,&
                     y, x, y,&
                     y, x,-y,&
                    -y, x, y,&
                    -y, x,-y,&
                     y,-x, y,&
                     y,-x,-y,&
                    -y,-x, y,&
                    -y,-x,-y,&
                     y, y, x,&
                     y,-y, x,&
                    -y, y, x,&
                    -y,-y, x,&
                     y, y,-x,&
                     y,-y,-x,&
                    -y, y,-x,&
                    -y,-y,-x/),(/3,nb/))
      END DO
      nlgp = nlgp + nb*np
    CASE ("C")
      DO i=1,np
        x = r(i)
        y = SQRT(one - x**2)
        lebedev_grid(lgnum)%w(nlgp+nc*(i-1)+1:nlgp+nc*i) = w(i)
        lebedev_grid(lgnum)%r(1:3,nlgp+nc*(i-1)+1:nlgp+nc*i) =&
          RESHAPE((/ x, y,zero,&
                     x,-y,zero,&
                    -x, y,zero,&
                    -x,-y,zero,&
                     y, x,zero,&
                     y,-x,zero,&
                    -y, x,zero,&
                    -y,-x,zero,&
                     x,zero, y,&
                     x,zero,-y,&
                    -x,zero, y,&
                    -x,zero,-y,&
                     y,zero, x,&
                     y,zero,-x,&
                    -y,zero, x,&
                    -y,zero,-x,&
                    zero, x, y,&
                    zero, x,-y,&
                    zero,-x, y,&
                    zero,-x,-y,&
                    zero, y, x,&
                    zero, y,-x,&
                    zero,-y, x,&
                    zero,-y,-x/),(/3,nc/))
      END DO
      nlgp = nlgp + nc*np
    CASE ("D")
      IF (MODULO(np,3) == 0) THEN
        DO i=1,np,3
          j = (i + 2)/3
          x = r(i)
          y = r(i+1)
          z = r(i+2)
          lebedev_grid(lgnum)%w(nlgp+nd*(j-1)+1:nlgp+nd*j) = w(j)
          lebedev_grid(lgnum)%r(1:3,nlgp+nd*(j-1)+1:nlgp+nd*j) =&
            RESHAPE((/ x, y, z, x, y,-z, x,-y, z,-x, y, z,&
                       x,-y,-z,-x, y,-z,-x,-y, z,-x,-y,-z,&
                       x, z, y, x, z,-y, x,-z, y,-x, z, y,&
                       x,-z,-y,-x, z,-y,-x,-z, y,-x,-z,-y,&
                       y, x, z, y, x,-z, y,-x, z,-y, x, z,&
                       y,-x,-z,-y, x,-z,-y,-x, z,-y,-x,-z,&
                       y, z, x, y, z,-x, y,-z, x,-y, z, x,&
                       y,-z,-x,-y, z,-x,-y,-z, x,-y,-z,-x,&
                       z, x, y, z, x,-y, z,-x, y,-z, x, y,&
                       z,-x,-y,-z, x,-y,-z,-x, y,-z,-x,-y,&
                       z, y, x, z, y,-x, z,-y, x,-z, y, x,&
                       z,-y,-x,-z, y,-x,-z,-y, x,-z,-y,-x/),(/3,nd/))
        END DO
      ELSE
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Subsytem D: np is not modulo 3 (check argument #3)")
      END IF
      nlgp = nlgp + nd*np/3
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The invalid subsystem <"//TRIM(subsystem)//"> was "//&
                        "specified (check argument #1)")
    END SELECT

!   *** Reinitialize working storage ***

    w(:) = zero
    r(:) = zero

  END SUBROUTINE load_sub_grid

! *****************************************************************************
!> \brief ...
! *****************************************************************************
  SUBROUTINE deallocate_lebedev_grids()

    INTEGER                                  :: ilg

    DO ilg=1,nlg
      DEALLOCATE (lebedev_grid(ilg)%r,lebedev_grid(ilg)%w)
    END DO

  END SUBROUTINE deallocate_lebedev_grids

! *****************************************************************************
!> \brief Load the coordinates and weights of the nonredundant Lebedev grid
!>         points.
!> \date    04.05.99
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_lebedev_grids()

    CHARACTER(len=*), PARAMETER :: routineN = 'init_lebedev_grids', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ilg

    CALL timeset(routineN,handle)

    DO ilg=1,nlg
      IF(ASSOCIATED(lebedev_grid(ilg)%r)) DEALLOCATE(lebedev_grid(ilg)%r)
      IF(ASSOCIATED(lebedev_grid(ilg)%w)) DEALLOCATE(lebedev_grid(ilg)%w)
    END DO

!   *** Initialize constants ***

    one = 1.0_dp
    rs2 = one/SQRT(2.0_dp)
    rs3 = one/SQRT(3.0_dp)
    zero = 0.0_dp

!   *** Initialize working storage ***

    w(:) = zero
    r(:) = zero

!   *** Load the angular momentum quantum numbers l of the Lebedev grids ***

    lebedev_grid(1:nlg)%l = (/3,5,7,9,11,15,17,19,23,29,35,41,47,53/)

!   *** Load the total number of grid points for each Lebedev grid ***

    lebedev_grid(1:nlg)%n = (/6,14,26,38,50,86,110,146,194,302,434,590,770,&
                              974/)

!   *** Allocate storage for the Lebedev grids ***

    DO ilg=1,nlg
      ALLOCATE (lebedev_grid(ilg)%r(3,lebedev_grid(ilg)%n),&
                lebedev_grid(ilg)%w(lebedev_grid(ilg)%n))
    END DO

!   *** Load the coordinates and weights of the Lebedev grid points ***

!   *** 1. l = 3 (6 points, octahedron) ***

    w(1) = 1.0_dp/6.0_dp
    CALL load_sub_grid("A1",1,0)

!   *** 2. l = 5 (14 points, capped octahedron) ***

    w(1) = 1.0_dp/15.0_dp
    CALL load_sub_grid("A1",2,0)

    w(1) = 3.0_dp/40.0_dp
    CALL load_sub_grid("A3",2,0)

!   *** 3. l = 7 (26 points) ***

    w(1) = 1.0_dp/21.0_dp
    CALL load_sub_grid("A1",3,0)

    w(1) = 4.0_dp/105.0_dp
    CALL load_sub_grid("A2",3,0)

    w(1) = 9.0_dp/280.0_dp
    CALL load_sub_grid("A3",3,0)

!   *** 4. l = 9 (38 points) ***

    w(1) = 1.0_dp/105.0_dp
    CALL load_sub_grid("A1",4,0)

    w(1) = 9.0_dp/280.0_dp
    CALL load_sub_grid("A3",4,0)

    w(1) = 1.0_dp/35.0_dp
    r(1) = rs2*SQRT(1.0_dp + rs3)
    CALL load_sub_grid("C",4,1)

!   *** 5. l = 11 (50 points) ***

    w(1) = 4.0_dp/315.0_dp
    CALL load_sub_grid("A1",5,0)

    w(1) = 64.0_dp/2835.0_dp
    CALL load_sub_grid("A2",5,0)

    w(1) = 27.0_dp/1280.0_dp
    CALL load_sub_grid("A3",5,0)

    w(1) = 14641.0_dp/725760.0_dp
    r(1) = 3.0_dp/SQRT(11.0_dp)
    CALL load_sub_grid("B",5,1)

!   *** 6. l = 15 (86 points) ***

    w(1) = 1.15440115440115440115440115440115E-2_dp
    CALL load_sub_grid("A1",6,0)

    w(1) = 1.19439090858562823236989259736470E-2_dp
    CALL load_sub_grid("A3",6,0)

    w(1) = 1.11105557106034025109468482160140E-2_dp
    r(1) = 8.52518311701267605338736780155357E-1_dp
    w(2) = 1.18765012945371420137882805994025E-2_dp
    r(2) = 1.89063552885395482707075847005288E-1_dp
    CALL load_sub_grid("B",6,2)

    w(1) = 1.18123037469044753644792263073650E-2_dp
    r(1) = 9.27330657151172465678969739310097E-1_dp
    CALL load_sub_grid("C",6,1)

!   *** 7. l = 17 (110 points) ***

    w(1) = 3.82827049493716160382827049493716E-3_dp
    CALL load_sub_grid("A1",7,0)

    w(1) = 9.79373751248751248751248751248751E-3_dp
    CALL load_sub_grid("A3",7,0)

    w(1) = 8.21173728319111097598993405227308E-3_dp
    r(1) = 9.65124035086594105655529546158531E-1_dp
    w(2) = 9.59547133607096284945318117290259E-3_dp
    r(2) = 8.28769981252592210694031500711749E-1_dp
    w(3) = 9.94281489117810328140065828526450E-3_dp
    r(3) = 2.15957291845848832354930328945946E-1_dp
    CALL load_sub_grid("B",7,3)

    w(1) = 9.69499636166302832969499636166303E-3_dp
    r(1) = 8.78158910604066133449110592671618E-1_dp
    CALL load_sub_grid("C",7,1)

!   *** 8. l = 19 (146 points) ***

    w(1) = 5.99631368862138092907323676554446E-4_dp
    CALL load_sub_grid("A1",8,0)

    w(1) = 7.37299971862075642305743268410561E-3_dp
    CALL load_sub_grid("A2",8,0)

    w(1) = 7.21051536014448777763305996834282E-3_dp
    CALL load_sub_grid("A3",8,0)

    w(1) = 7.57439415905403372268748574713806E-3_dp
    r(1) = 9.74888643677173235480043489928723E-1_dp
    w(2) = 6.75382948631447744073541732486459E-3_dp
    r(2) = 8.07089818359582501629574989264033E-1_dp
    w(3) = 7.11635549311755538760089284953968E-3_dp
    r(3) = 2.91298882209526746288335878313128E-1_dp
    CALL load_sub_grid("B",8,3)

    w(1) = 6.99108735330326239417148508057599E-3_dp
    r(1) = 1.40355381171318328571556780746292E-1_dp
    r(2) = 4.49332832326955734884695382705101E-1_dp
    r(3) = 8.82270011260322631916665753046583E-1_dp
    CALL load_sub_grid("D",8,3)

!   *** 9. l = 23 (194 points) ***

    w(1) = 1.78234044724461115736727104869868E-3_dp
    CALL load_sub_grid("A1",9,0)

    w(1) = 5.71690594997710189299212838832099E-3_dp
    CALL load_sub_grid("A2",9,0)

    w(1) = 5.57338317884873796836784958446647E-3_dp
    CALL load_sub_grid("A3",9,0)

    w(1) = 5.51877146727361369172768460119380E-3_dp
    r(1) = 7.77493219314767127213777704028801E-1_dp
    w(2) = 5.15823771180538310324916154718792E-3_dp
    r(2) = 9.12509096867473692992168387555434E-1_dp
    w(3) = 5.60870408258799684374936673855184E-3_dp
    r(3) = 3.14196994182586079225390955831946E-1_dp
    w(4) = 4.10677702816939409072861128564582E-3_dp
    r(4) = 9.82972302707253296863729593965148E-1_dp
    CALL load_sub_grid("B",9,4)

    w(1) = 5.05184606461480847598931196006390E-3_dp
    r(1) = 9.38319218137591520905616389195671E-1_dp
    CALL load_sub_grid("C",9,1)

    w(1) = 5.53024891623309370129768269143303E-3_dp
    r(1) = 1.59041710538352952424263581362096E-1_dp
    r(2) = 5.25118572443642024905268207753317E-1_dp
    r(3) = 8.36036015482458885943746377793309E-1_dp
    CALL load_sub_grid("D",9,3)

!   *** 10. l = 29 (302 points) ***

    w(1) = 8.54591172512814813423121032618880E-4_dp
    CALL load_sub_grid("A1",10,0)

    w(1) = 3.59911928502557145886397858961119E-3_dp
    CALL load_sub_grid("A3",10,0)

    w(1) = 3.65004580767725542865433220112651E-3_dp
    r(1) = 1.29238672710514925339493976600550E-1_dp
    w(2) = 3.60482260141988171131480913104353E-3_dp
    r(2) = 3.71034178384821189405344300671593E-1_dp
    w(3) = 3.57672966174336707556208137560857E-3_dp
    r(3) = 7.43452042987555751133204566760209E-1_dp
    w(4) = 3.44978842430588331001302771048283E-3_dp
    r(4) = 8.67643624544083327473890295837421E-1_dp
    w(5) = 3.10895312241367525484587698082774E-3_dp
    r(5) = 9.49454317226443084214869821724903E-1_dp
    w(6) = 2.35210141368916437879217118337424E-3_dp
    r(6) = 9.90705621379408123821774729602922E-1_dp
    CALL load_sub_grid("B",10,6)

    w(1) = 3.60082093221646027279920634177099E-3_dp
    r(1) = 8.20326419827759303328870367500864E-1_dp
    w(2) = 2.98234496317180385195111046924520E-3_dp
    r(2) = 9.64408914879206014987053264857406E-1_dp
    CALL load_sub_grid("C",10,2)

    w(1) = 3.57154055427338708123297920312395E-3_dp
    r(1) = 2.51003475177046506904110494820777E-1_dp
    r(2) = 5.44867737258077380269062623720797E-1_dp
    r(3) = 8.00072749407395172686232883863933E-1_dp
    w(2) = 3.39231220500617018197882653945702E-3_dp
    r(4) = 1.23354853258332742165467555277499E-1_dp
    r(5) = 4.12772408316853095996383989497138E-1_dp
    r(6) = 9.02442529533000401095678759353334E-1_dp
    CALL load_sub_grid("D",10,6)

!   *** 11. l = 35 (434 points) ***

    w(1) = 5.26589796822443623921598397444095E-4_dp
    CALL load_sub_grid("A1",11,0)

    w(1) = 2.54821997200260718024899528063707E-3_dp
    CALL load_sub_grid("A2",11,0)

    w(1) = 2.51231741892730716751285677008112E-3_dp
    CALL load_sub_grid("A3",11,0)

    w(1) = 1.46249562159461384222611198788353E-3_dp
    r(1) = 9.94255912631277812357804047239411E-1_dp
    w(2) = 2.01427902091852819554717351102614E-3_dp
    r(2) = 9.67987158791472796844440561853288E-1_dp
    w(3) = 2.30269478222741577118243755624539E-3_dp
    r(3) = 9.14472801120872429366691177286487E-1_dp
    w(4) = 2.44537343731298000188170923067793E-3_dp
    r(4) = 8.31584400419232294499745594669204E-1_dp
    w(5) = 2.50172516840293614549327311336461E-3_dp
    r(5) = 7.19016501040843432905788800294849E-1_dp
    w(6) = 2.51326717459756436118956600426018E-3_dp
    r(6) = 4.07712664897769512483055658274718E-1_dp
    w(7) = 2.53040380118635500209853222454729E-3_dp
    r(7) = 2.12646824707552073524018894972799E-1_dp
    CALL load_sub_grid("B",11,7)

    w(1) = 1.91095128217953227363116198076372E-3_dp
    r(1) = 9.77642811118264871364737333213224E-1_dp
    w(2) = 2.41744237563898077608784837021936E-3_dp
    r(2) = 8.81813287779428812899447039158355E-1_dp
    CALL load_sub_grid("C",11,2)

    w(1) = 2.23660776043784866397235325901501E-3_dp
    r(1) = 9.92176963642923726861469072682696E-2_dp
    r(2) = 3.34436314534345491435541212141297E-1_dp
    r(3) = 9.37180985855372235789161206437430E-1_dp
    w(2) = 2.41693004432477530950382860647420E-3_dp
    r(4) = 2.05482369640304370826435780923978E-1_dp
    r(5) = 4.50233038258262537898138529993494E-1_dp
    r(6) = 8.68946032287241197946412945028781E-1_dp
    w(3) = 2.51223685456349510228853571762152E-3_dp
    r(7) = 1.06801826075804827808519981079378E-1_dp
    r(8) = 5.90515704892527108597639202779953E-1_dp
    r(9) = 7.99927854385728570494616568630393E-1_dp
    w(4) = 2.49664405455308596359511734587523E-3_dp
    r(10) = 3.10428403516654146818973681387686E-1_dp
    r(11) = 5.55015236107680716116395510647568E-1_dp
    r(12) = 7.71746262691590088133470564383333E-1_dp
    CALL load_sub_grid("D",11,12)

!   *** 12. l = 41 (590 points) ***

    w(1) = 3.09512129530618734224885916984167E-4_dp
    CALL load_sub_grid("A1",12,0)

    w(1) = 1.85237969859748902097779927120955E-3_dp
    CALL load_sub_grid("A3",12,0)

    w(1) = 9.76433116505105003063946798281256E-4_dp
    r(1) = 9.96278129754016372294779397439010E-1_dp
    w(2) = 1.38473723485169190049558481046633E-3_dp
    r(2) = 9.78480583762693839428807199465825E-1_dp
    w(3) = 1.61721064725441119207240564027320E-3_dp
    r(3) = 9.41414158220402542954812123169147E-1_dp
    w(4) = 1.74956465728115411858337657547711E-3_dp
    r(4) = 8.83078727934132546639472279693012E-1_dp
    w(5) = 1.81847177816276877817839645350131E-3_dp
    r(5) = 8.02836877335273789544708232847842E-1_dp
    w(6) = 1.84671595615124182496431643927999E-3_dp
    r(6) = 7.00768575373572920755474238992223E-1_dp
    w(7) = 1.85202882829621309713188574624172E-3_dp
    r(7) = 4.33373868777154278470427918169993E-1_dp
    w(8) = 1.85881258543831701688713685607976E-3_dp
    r(8) = 2.70356088359165054450644840302674E-1_dp
    w(9) = 1.87179063927774375088198703268059E-3_dp
    r(9) = 9.21904070768989460773325829762845E-2_dp
    CALL load_sub_grid("B",12,9)

    w(1) = 1.30032168588604773254689500223665E-3_dp
    r(1) = 9.85013335028001910431569760301413E-1_dp
    w(2) = 1.70515399639586401999426890257133E-3_dp
    r(2) = 9.18045287711453949033261544688478E-1_dp
    w(3) = 1.85716119677407798248150549753942E-3_dp
    r(3) = 7.91101929626901988791072435283520E-1_dp
    CALL load_sub_grid("C",12,3)

    w(1) = 1.55521360339680849695800543399046E-3_dp
    r(1) = 8.21302158193251139256011433978433E-2_dp
    r(2) = 2.77867319058624428710339115104831E-1_dp
    r(3) = 9.57102074310072578548507660672158E-1_dp
    w(2) = 1.80223912800852549116687840059028E-3_dp
    r(4) = 8.99920584207487492736600381248490E-2_dp
    r(5) = 5.03356427107511721798695016575784E-1_dp
    r(6) = 8.59379855890721226547417613471625E-1_dp
    w(3) = 1.84983056044366016523678105066982E-3_dp
    r(7) = 1.81664084036020946133369276322890E-1_dp
    r(8) = 5.98412649788537962548021778158698E-1_dp
    r(9) = 7.80320742479920327279647082439668E-1_dp
    w(4) = 1.71390450710670866804340083221112E-3_dp
    r(10) = 1.72079522565687812382959749101837E-1_dp
    r(11) = 3.79103540769556328153581995776882E-1_dp
    r(12) = 9.09213475092373607255560905376387E-1_dp
    w(5) = 1.80265893437745115894760427403985E-3_dp
    r(13) = 2.63471665593794963182134911777324E-1_dp
    r(14) = 4.74239284255198024276965888711084E-1_dp
    r(15) = 8.40047488359050429417090817376105E-1_dp
    w(6) = 1.84286647290528563239873317969286E-3_dp
    r(16) = 3.51828092773351899746468433493901E-1_dp
    r(17) = 5.61026380862206018993247924039616E-1_dp
    r(18) = 7.49310611904115932021924350883623E-1_dp
    CALL load_sub_grid("D",12,18)

!   *** 13. l = 47 (770 points) ***

    w(1) = 2.19294208818118413191132531546907E-4_dp
    CALL load_sub_grid("A1",13,0)

    w(1) = 1.43643361731907982089311290725143E-3_dp
    CALL load_sub_grid("A2",13,0)

    w(1) = 1.42194034433587736471699229343586E-3_dp
    CALL load_sub_grid("A3",13,0)

    w(1) = 6.79812351105050201620332392505941E-4_dp
    r(1) = 9.97408677652823043744886496130338E-1_dp
    w(2) = 9.91318423529491220896771270468706E-4_dp
    r(2) = 9.84799753572301179668683394920869E-1_dp
    w(3) = 1.18020783323894879580204544917243E-3_dp
    r(3) = 9.58036675983391439732458619060098E-1_dp
    w(4) = 1.29659960208092067423594734151509E-3_dp
    r(4) = 9.15317950483154831654282453968846E-1_dp
    w(5) = 1.36587142742831640181200652406724E-3_dp
    r(5) = 8.55901928697886435270454067917142E-1_dp
    w(6) = 1.40298860477532526267417750582642E-3_dp
    r(6) = 7.79621319527635171330859492071650E-1_dp
    w(7) = 1.41864556359560936142032713234525E-3_dp
    r(7) = 6.86644447264154195285405155618666E-1_dp
    w(8) = 1.42137674185166176103362636565150E-3_dp
    r(8) = 4.52311920313658441442123883818590E-1_dp
    w(9) = 1.42399647549096162940277119054908E-3_dp
    r(9) = 3.12521305001653125265256223110153E-1_dp
    w(10) = 1.43155404217856675285551659613462E-3_dp
    r(10) = 1.60155803498828974610387120985563E-1_dp
    CALL load_sub_grid("B",13,10)

    w(1) = 9.25440149986536789398468419487000E-4_dp
    r(1) = 9.89477537495598496934466234590453E-1_dp
    w(2) = 1.25023999505350931588903916766368E-3_dp
    r(2) = 9.40776878793758755393692765525837E-1_dp
    w(3) = 1.39436584332923012313782491233921E-3_dp
    r(3) = 8.45749305193653306776973233373898E-1_dp
    CALL load_sub_grid("C",13,3)

    w(1) = 1.12708909467174883469499826293069E-3_dp
    r(1) = 6.94402439334941301856132689108761E-2_dp
    r(2) = 2.35518789424232641745294975361896E-1_dp
    r(3) = 9.69385863498432080572651262482050E-1_dp
    w(2) = 1.34575376091067007284708159178371E-3_dp
    r(4) = 2.26900410952945985812144550510454E-1_dp
    r(5) = 4.10218247404573022502091360294078E-1_dp
    r(6) = 8.83310360522112719095559130461882E-1_dp
    w(3) = 1.42495728331678280442457616193820E-3_dp
    r(7) = 8.02557460777533890494192154117455E-2_dp
    r(8) = 6.21430241748160458082227723944198E-1_dp
    r(9) = 7.79348105702660973885784203335648E-1_dp
    w(4) = 1.26152334123774998512114871120604E-3_dp
    r(10) = 1.46799952789657197013820141783396E-1_dp
    r(11) = 3.24528434571739443574137218155171E-1_dp
    r(12) = 9.34414827052402165196138643907493E-1_dp
    w(5) = 1.39254710605269594012443620424449E-3_dp
    r(13) = 1.57150776982472710564392240262888E-1_dp
    r(14) = 5.22448218969662967424959711591943E-1_dp
    r(15) = 8.38064133458312497889043403162912E-1_dp
    w(6) = 1.41876167787765638610216672064623E-3_dp
    r(16) = 2.36570299315724563080851921006450E-1_dp
    r(17) = 6.01754663408955808458505175650313E-1_dp
    r(18) = 7.62840624604669826410820436887813E-1_dp
    w(7) = 1.33836668447955413932118842030907E-3_dp
    r(19) = 7.71481586676573204146393256189282E-2_dp
    r(20) = 4.34657551614116278920169057936696E-1_dp
    r(21) = 8.97285336132833396499244991766618E-1_dp
    w(8) = 1.39370086267613137196129133534721E-3_dp
    r(22) = 3.06293666621073020830889612421600E-1_dp
    r(23) = 4.90882658903761621546950183900134E-1_dp
    r(24) = 8.15609223203975422051589658409773E-1_dp
    w(9) = 1.41591475746693200265535028223575E-3_dp
    r(25) = 3.82247737952478700050265491776632E-1_dp
    r(26) = 5.64876814909950046651375249432877E-1_dp
    r(27) = 7.31300793659765701947278263710234E-1_dp
    CALL load_sub_grid("D",13,27)

!   *** 14. l = 53 (974 points) ***

    w(1) = 1.43829419052743111472634384792800E-4_dp
    CALL load_sub_grid("A1",14,0)

    w(1) = 1.12577228828700411922446003897444E-3_dp
    CALL load_sub_grid("A3",14,0)

    w(1) = 4.94802934194924095056628951683154E-4_dp
    r(1) = 9.98155345023846501191494527575073E-1_dp
    w(2) = 7.35799010912547049057334785150092E-4_dp
    r(2) = 9.88883224354685545466876075112182E-1_dp
    w(3) = 8.88913277130438430936840534873175E-4_dp
    r(3) = 9.68890220434707404463853311788372E-1_dp
    w(4) = 9.88834783892143487490379018339493E-4_dp
    r(4) = 9.36602730407163187902332491169310E-1_dp
    w(5) = 1.05329968170947064968554154238842E-3_dp
    r(5) = 8.91267942647606054508357503868536E-1_dp
    w(6) = 1.09277880701457848752501645357106E-3_dp
    r(6) = 8.32596723702351872809065001851482E-1_dp
    w(7) = 1.11438939406322716311544826926779E-3_dp
    r(7) = 7.60582905315251495676259855686005E-1_dp
    w(8) = 1.12372478805155528926212353749792E-3_dp
    r(8) = 6.75400969108414309734349145075511E-1_dp
    w(9) = 1.12523932524381359890224350348300E-3_dp
    r(9) = 4.66858905695743283152769352096062E-1_dp
    w(10) = 1.12615327181590500171763113615557E-3_dp
    r(10) = 3.44613654237438223690005816795002E-1_dp
    w(11) = 1.13028693112384079609639852411886E-3_dp
    r(11) = 2.11954151850184645694140401473565E-1_dp
    w(12) = 1.13498653436395488643714589963293E-3_dp
    r(12) = 7.16244014499556615845369803814026E-2_dp
    CALL load_sub_grid("B",14,12)

    w(1) = 6.82336792710993096637138363442284E-4_dp
    r(1) = 9.92323565431490196070394244122238E-1_dp
    w(2) = 9.45415816044709583567874030695969E-4_dp
    r(2) = 9.55781512496548488574596747552425E-1_dp
    w(3) = 1.07442997538567914061819358161104E-3_dp
    r(3) = 8.82785980701181710576015026411856E-1_dp
    w(4) = 1.12930008656913165822354345410170E-3_dp
    r(4) = 7.73778447257374736807517446423091E-1_dp
    CALL load_sub_grid("C",14,4)

    w(1) = 8.43688450090195438498026316671055E-4_dp
    r(1) = 5.97404861418134181366741615801962E-2_dp
    r(2) = 2.02912875277752280821562682821390E-1_dp
    r(3) = 9.77372722845309994986723236843860E-1_dp
    w(2) = 1.07525572044888463296011392358171E-3_dp
    r(4) = 1.37576040847363647998707229269443E-1_dp
    r(5) = 4.60262194248405393969672054729315E-1_dp
    r(6) = 8.77058461865802689033970747916718E-1_dp
    w(3) = 1.10857723686446203483870096411918E-3_dp
    r(7) = 3.39101652633628571879780724877562E-1_dp
    r(8) = 5.03067399966203571325278413235681E-1_dp
    r(9) = 7.94942299964208514029117644252761E-1_dp
    w(4) = 9.56647532378335729987022997768218E-4_dp
    r(10) = 1.27167519143981950887694563391455E-1_dp
    r(11) = 2.81760642244213432084049569237994E-1_dp
    r(12) = 9.51020169374389952612600596224493E-1_dp
    w(5) = 1.08066325071739071464463316273343E-3_dp
    r(13) = 2.69312074041351249441688496163892E-1_dp
    r(14) = 4.33156129172015739916409763707103E-1_dp
    r(15) = 8.60143461601761963548074094490711E-1_dp
    w(6) = 1.12679713119629459182839576027292E-3_dp
    r(16) = 1.41978645260191825726693154370171E-1_dp
    r(17) = 6.25616735858081415467234398457286E-1_dp
    r(18) = 7.67102186220558360145428572382623E-1_dp
    w(7) = 1.02256871535806121135514702022316E-3_dp
    r(19) = 6.70928460073825495111114896527700E-2_dp
    r(20) = 3.79839521685915669705644712365302E-1_dp
    r(21) = 9.22616110730809020138255560099254E-1_dp
    w(8) = 1.10896026771310753764243318312336E-3_dp
    r(22) = 7.05773818325617225225155364651637E-2_dp
    r(23) = 5.51750542142352049132113965320893E-1_dp
    r(24) = 8.31017552413474225750006298786621E-1_dp
    w(9) = 1.12279065343576582686387381292198E-3_dp
    r(25) = 2.78388847788215460217462722823059E-1_dp
    r(26) = 6.02961915615918685584436745852250E-1_dp
    r(27) = 7.47620610834085792324588316104014E-1_dp
    w(10) = 1.03240184711745983700951101978720E-3_dp
    r(28) = 1.97957893891740690574322079335028E-1_dp
    r(29) = 3.58960632958909578333772677316348E-1_dp
    r(30) = 9.12118378409121483660044157731511E-1_dp
    w(11) = 1.10724938228385393683388780275631E-3_dp
    r(31) = 2.08730706110327399856945819943227E-1_dp
    r(32) = 5.34866643813547649971323458627733E-1_dp
    r(33) = 8.18748536281021806431067610000242E-1_dp
    w(12) = 1.12178004851997205823535727177375E-3_dp
    r(34) = 4.05512213787283588596898056051174E-1_dp
    r(35) = 5.67499754607437348401393912269671E-1_dp
    r(36) = 7.16591845467023718833743633176626E-1_dp
    CALL load_sub_grid("D",14,36)

    init_lebedev_grids_done = .TRUE.

    CALL timestop(handle)

  END SUBROUTINE init_lebedev_grids

END MODULE lebedev
