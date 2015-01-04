!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Arrays of parameters used in the semi-empirical calculations
!> \References  Everywhere in this module TCA stands for:
!>   - TCA:   W. Thiel and  A. A. Voityuk - Teor. Chim. Acta (1992) 81:391-404
!>   - TCA77: M.J.S. Dewar and W. Thiel - Teor. Chim. Acta (1977) 46:89-104
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
MODULE semi_empirical_int_arrays

  
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_int_arrays'
  REAL(KINd=dp), PARAMETER, PUBLIC     :: rij_threshold = 0.00002_dp

  ! Mapping index for orbital ordering
  INTEGER, DIMENSION(9), PUBLIC        :: se_orbital_pointer=(/1,4,2,3,9,8,7,6,5/)
  INTEGER, DIMENSION(9), PUBLIC        :: se_map_alm        =(/1,3,4,2,8,6,5,7,9/)

  ! Arrays to treat the invertion of the reference frame for the integrals
  ! Using the same indexing convention of the integrals
  INTEGER, PARAMETER, DIMENSION(2,9), PUBLIC        :: map_x_to_z = RESHAPE((/&
       1,0,& ! s    <-> s
       4,0,& ! px   <-> pz
       3,0,& ! py   <-> py
       2,0,& ! pz   <-> px
       7,5,& ! dx2  <-> SQR3/2*dz2 +    0.5*dx2
       6,0,& ! dzx  <-> dzx
       7,5,& ! dz2  <->   -0.5*dz2 + SQR3/2*dx2
       9,0,& ! dzy  <-> dxy
       8,0 & ! dxy  <-> dzy
       /),(/2,9/))
  REAL(KIND=dp), PARAMETER, DIMENSION(2,9), PUBLIC  :: fac_x_to_z =  RESHAPE((/&
        1.0_dp,0.0_dp,&
        1.0_dp,0.0_dp,&
        1.0_dp,0.0_dp,&
        1.0_dp,0.0_dp,&
        0.8660254037844386_dp,0.5_dp,&
        1.0_dp,0.0_dp,&
       -0.5_dp,0.8660254037844386_dp,&
        1.0_dp,0.0_dp,&
        1.0_dp,0.0_dp &
        /),(/2,9/))

  ! Clm coefficients for d-orbitals:  see  Table [1] and [2] of TCA
  REAL(KIND=dp), DIMENSION (45,0:2,-2:2),  PUBLIC :: clm_d
  ! Clm coefficients for sp-orbitals: see original paper TCA77
  INTEGER,       DIMENSION (45,0:2,-2:2),  PUBLIC :: clm_sp
  ! alm coefficients: see Laino and Hutter (periodic SE)
  REAL(KIND=dp), DIMENSION (45,0:2,-2:2),  PUBLIC :: alm

  ! These values are absolutely arbitrary and are used only for a proper
  ! tag of the integrals
  INTEGER, PARAMETER, PUBLIC :: &
       CLMz=10, CLMp=11, CLMzz=12, CLMzp=13, CLMyy=14, CLMxy=15, CLMxx=16

  ! Indexes for diagonal storage of ij and kl multipoles
  INTEGER, DIMENSION(9,9), PUBLIC                 :: indexa, indexb

  ! Type of integral for 2electron 2centers integrals
  INTEGER, DIMENSION(45), PARAMETER, PUBLIC       :: int2c_type =(/&
       1, 2, 3, 2, 3, 3, 2, 3, 3, 3, 4, 5, 5, 5, 6, 4, 5, 5, 5, &
       6, 6, 4, 5, 5, 5, 6, 6, 6, 4, 5, 5, 5, 6, 6, 6, 6, 4, 5, &
       5, 5, 6, 6, 6, 6, 6/)

  ! Mappinf of shell index
  INTEGER, DIMENSION(9), PARAMETER, PUBLIC        :: l_index=(/&
       0, 1, 1, 1, 2, 2, 2, 2, 2/)

  ! Index for <ij|kl>
  INTEGER, DIMENSION (45, 45), PUBLIC             :: ijkl_ind

  ! Symmetry index for <ij|kl>
  INTEGER, DIMENSION (491), PUBLIC                :: ijkl_sym

  ! Index for integral rotations
  INTEGER, DIMENSION (3, 3), PUBLIC               :: indpp
  INTEGER, DIMENSION (5, 3), PUBLIC               :: inddp
  INTEGER, DIMENSION (5, 5), PUBLIC               :: inddd

  ! Indexes use for the construction of the one-center two-electron integrals
  INTEGER, DIMENSION (243), PUBLIC ::  int_ij = (/&
        1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5,&
        5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9,10,10,10,&
       10,10,10,11,11,11,11,11,11,12,12,12,12,12,13,13,13,13,13,14,14,14,15,&
       15,15,15,15,15,15,15,15,15,16,16,16,16,16,17,17,17,17,17,18,18,18,19,&
       19,19,19,19,20,20,20,20,20,21,21,21,21,21,21,21,21,21,21,21,21,22,22,&
       22,22,22,22,22,22,22,23,23,23,23,23,24,24,24,24,24,25,25,25,25,26,26,&
       26,26,26,26,27,27,27,27,27,28,28,28,28,28,28,28,28,28,28,29,29,29,29,&
       29,30,30,30,31,31,31,31,31,32,32,32,32,32,33,33,33,33,33,34,34,34,34,&
       35,35,35,35,35,36,36,36,36,36,36,36,36,36,36,36,36,37,37,37,37,38,38,&
       38,38,38,39,39,39,39,39,40,40,40,41,42,42,42,42,42,43,43,43,43,44,44,&
       44,44,44,45,45,45,45,45,45,45,45,45,45 /)
  INTEGER, DIMENSION (243), PUBLIC ::  int_kl = (/&
       15,21,28,36,45,12,19,23,39,11,15,21,22,26,28,36,45,13,24,32,38,34,37,&
       43,11,15,21,22,26,28,36,45,17,25,31,16,20,27,44,29,33,35,42,15,21,22,&
       28,36,45, 3, 6,11,21,26,36, 2,12,19,23,39, 4,13,24,32,38,14,17,31, 1,&
        3, 6,10,15,21,22,28,36,45, 8,16,20,27,44, 7,14,17,25,31,18,30,40, 2,&
       12,19,23,39, 8,16,20,27,44, 1, 3, 6,10,11,15,21,22,26,28,36,45, 3, 6,&
       10,15,21,22,28,36,45, 2,12,19,23,39, 4,13,24,32,38, 7,17,25,31, 3, 6,&
       11,21,26,36, 8,16,20,27,44, 1, 3, 6,10,15,21,22,28,36,45, 9,29,33,35,&
       42,18,30,40, 7,14,17,25,31, 4,13,24,32,38, 9,29,33,35,42, 5,34,37,43,&
        9,29,33,35,42, 1, 3, 6,10,11,15,21,22,26,28,36,45, 5,34,37,43, 4,13,&
       24,32,38, 2,12,19,23,39,18,30,40,41, 9,29,33,35,42, 5,34,37,43, 8,16,&
       20,27,44, 1, 3, 6,10,15,21,22,28,36,45 /)
  INTEGER, DIMENSION (243), PUBLIC :: int_onec2el = (/&
        1, 1, 1, 1, 1, 3, 3, 8, 3, 9, 6, 6,12,14,13, 7, 6,15, 8, 3, 3,11, 9,&
       14,17, 6, 7,12,18,13, 6, 6, 3, 2, 3, 9,11,10,11, 9,16,10,11, 7, 6, 4,&
        5, 6, 7, 9,17,19,32,22,40, 3,33,34,27,46,15,33,28,41,47,35,35,42, 1,&
        6, 6, 7,29,38,22,31,38,51, 9,19,32,21,32, 3,35,33,24,34,35,35,35, 3,&
       34,33,26,34,11,32,44,37,49, 1, 6, 7, 6,32,38,29,21,39,30,38,38,12,12,&
        4,22,21,19,20,21,22, 8,27,26,25,27, 8,28,25,26,27, 2,24,23,24,14,18,&
       22,39,48,45,10,21,37,36,37, 1,13,13, 5,31,30,20,29,30,31, 9,19,40,21,&
       32,35,35,35, 3,42,34,24,33, 3,41,26,33,34,16,40,44,43,50,11,44,32,39,&
       10,21,43,36,37, 1, 7, 6, 6,40,38,38,21,45,30,29,38, 9,32,19,22, 3,47,&
       27,34,33, 3,46,34,27,33,35,35,35,52,11,32,50,37,44,14,39,22,48,11,32,&
       49,37,44, 1, 6, 6, 7,51,38,22,31,38,29 /)


  PUBLIC :: init_se_intd_array

CONTAINS
! *****************************************************************************
!> \brief  Initialize all arrays used for the evaluation of the integrals
!>
!> \param error ...
!> \date   04.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE init_se_intd_array(error)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_se_intd_array', &
      routineP = moduleN//':'//routineN

    CALL setup_index_array(error)
    CALL setup_indrot_array(error)
    CALL setup_clm_array(error)
    CALL setup_ijkl_array(error)

  END SUBROUTINE init_se_intd_array

! *****************************************************************************
!> \brief Fills in array for the diagonal storage of the ij and kl multipoles term
!>
!> \param error ...
!> \date   03.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE setup_index_array(error)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_index_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: failure

    failure = .FALSE.
    DO i = 1, 9
       DO j = 1, i
          ! indexa:
          !           s   pz   px   py  dz2  dzx  dzy dx2-y2 dxy
          !      s    1    2    3    4    5    6    7    8    9
          !      pz   2   10   11   12   13   14   15   16   17
          !      px   3   11   18   19   20   21   22   23   24
          !      py   4   12   19   25   26   27   28   29   30
          !     dz2   5   13   20   26   31   32   33   34   35
          !     dzx   6   14   21   27   32   36   37   38   39
          !     dzy   7   15   22   28   33   37   40   41   42
          !  dx2-y2   8   16   23   29   34   38   41   43   44
          !     dxy   9   17   24   30   35   39   42   44   45
          indexa(i,j) = (9*(j-1))-(j*(j-1))/2+i
          indexa(j,i) = indexa(i,j)
          ! indexb:
          !           s   pz   px   py  dz2  dzx  dzy dx2-y2 dxy
          !      s    1    2    4    7   11   16   22   29   37
          !      pz   2    3    5    8   12   17   23   30   38
          !      px   4    5    6    9   13   18   24   31   39
          !      py   7    8    9   10   14   19   25   32   40
          !     dz2  11   12   13   14   15   20   26   33   41
          !     dzx  16   17   18   19   20   21   27   34   42
          !     dzy  22   23   24   25   26   27   28   35   43
          !  dx2-y2  29   30   31   32   33   34   35   36   44
          !     dxy  37   38   39   40   41   42   43   44   45
          indexb(i,j) = (i*(i-1))/2+j
          indexb(j,i) = indexb(i,j)
       END DO
    END DO
  END SUBROUTINE setup_index_array

! *****************************************************************************
!> \brief Fills in array for the rotation of the integrals
!>
!> \param error ...
!> \date   04.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE setup_indrot_array(error)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_indrot_array', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    ! Setup indexes for integral rotations
    ! INDPP
    indpp(1, 1) =  1
    indpp(2, 1) =  4
    indpp(3, 1) =  5
    indpp(1, 2) =  4
    indpp(2, 2) =  2
    indpp(3, 2) =  6
    indpp(1, 3) =  5
    indpp(2, 3) =  6
    indpp(3, 3) =  3
    ! INDDP
    inddp(1, 1) =  1
    inddp(2, 1) =  4
    inddp(3, 1) =  7
    inddp(4, 1) = 10
    inddp(5, 1) = 13
    inddp(1, 2) =  2
    inddp(2, 2) =  5
    inddp(3, 2) =  8
    inddp(4, 2) = 11
    inddp(5, 2) = 14
    inddp(1, 3) =  3
    inddp(2, 3) =  6
    inddp(3, 3) =  9
    inddp(4, 3) = 12
    inddp(5, 3) = 15
    ! INDDD
    inddd(1, 1) =  1
    inddd(2, 1) =  6
    inddd(3, 1) =  7
    inddd(4, 1) =  9
    inddd(5, 1) = 12
    inddd(1, 2) =  6
    inddd(2, 2) =  2
    inddd(3, 2) =  8
    inddd(4, 2) = 10
    inddd(5, 2) = 13
    inddd(1, 3) =  7
    inddd(2, 3) =  8
    inddd(3, 3) =  3
    inddd(4, 3) = 11
    inddd(5, 3) = 14
    inddd(1, 4) =  9
    inddd(2, 4) = 10
    inddd(3, 4) = 11
    inddd(4, 4) =  4
    inddd(5, 4) = 15
    inddd(1, 5) = 12
    inddd(2, 5) = 13
    inddd(3, 5) = 14
    inddd(4, 5) = 15
    inddd(5, 5) =  5
  END SUBROUTINE setup_indrot_array

! *****************************************************************************
!> \brief Fills in Clm coefficients (see Table [2] of TCA)
!>
!> \param error ...
!> \date   03.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE setup_clm_array(error)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_clm_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: CLM1, CLM1m
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: ALM1, ALMs15_49, ALMs15_49m, ALMs15m, ALMs20_49, &
      ALMs20_49m, ALMs35, ALMs35m, ALMs45, ALMs5_49, ALMs5_49m, CLM23, &
      CLM23m, CLM43, CLM43m, CLMs13, CLMs13m, CLMs43, CLMs43m

    failure = .FALSE.
    CLM1       =  1
    CLM1m      = -1
    CLMs13     =  SQRT(1.0_dp/3.0_dp)
    CLMs13m    = -SQRT(1.0_dp/3.0_dp)
    CLM23      =  2.0_dp/3.0_dp
    CLM23m     = -2.0_dp/3.0_dp
    CLM43      =  4.0_dp/3.0_dp
    CLM43m     = -4.0_dp/3.0_dp
    CLMs43     =  SQRT(4.0_dp/3.0_dp)
    CLMs43m    = -SQRT(4.0_dp/3.0_dp)
    ALM1       =  1.0_dp
    ALMs45     =  SQRT(4.0_dp/5.0_dp)
    ALMs35     =  SQRT(3.0_dp/5.0_dp)
    ALMs35m    = -ALMs35
    ALMs15m    = -SQRT(1.0_dp/5.0_dp)
    ALMs20_49  =  SQRT(20.0_dp/49.0_dp)
    ALMs20_49m = -ALMs20_49
    ALMs5_49   =  SQRT(5.0_dp/49.0_dp)
    ALMs5_49m  = -ALMs5_49
    ALMs15_49  =  SQRT(15.0_dp/49.0_dp)
    ALMs15_49m = -ALMs15_49

    ! Notation (1) s
    !          (2) pz = p_sigma             (5) d_sigma           = dz2  (8) d_delta               = dx2-y2
    !          (3) px = p_pi                (6) d_pi              = dzx  (9) d_{\overline{delta}}  = dxy
    !          (4) py = pi_{\overline{pi}}  (7) d_{\overline{pi}} = dzy

    clm_d = 0.0_dp
    clm_sp= 0
    alm   = 0.0_dp
    ! Let's fill all element of table 1 with resulting multipole lesser than 2
    ! Important Note: the value of the clm_sp does not reflect any phyisical
    !                 rule for decomposing/summing multipoles. It's just a
    !                 computational trick to put some order where conceptual mess
    !                 has been created..
    ! s s
    clm_d ( 1, 0, 0) = CLM1
    clm_sp( 1, 0, 0) = CLM1
    alm   ( 1, 0, 0) = ALM1
    ! s pz
    clm_d ( 2, 1, 0) = CLM1
    clm_sp( 2, 1, 0) = CLMz
    alm   ( 2, 1, 0) = ALM1
    ! s px
    clm_d ( 3, 1, 1) = CLM1
    clm_sp( 3, 1, 1) = CLMp
    alm   ( 3, 1, 1) = ALM1
    ! s py
    clm_d ( 4, 1,-1) = CLM1
    clm_sp( 4, 1,-1) = CLMp
    alm   ( 4, 1,-1) = ALM1
    ! s dz2
    clm_d ( 5, 2, 0) = CLMs43
    alm   ( 5, 2, 0) = ALM1
    ! s dzx
    clm_d ( 6, 2, 1) = CLM1
    alm   ( 6, 2, 1) = ALM1
    ! s dzy
    clm_d ( 7, 2,-1) = CLM1
    alm   ( 7, 2,-1) = ALM1
    ! s dx2-y2
    clm_d ( 8, 2, 2) = CLM1
    alm   ( 8, 2, 2) = ALM1
    ! s dxy
    clm_d ( 9, 2,-2) = CLM1
    alm   ( 9, 2,-2) = ALM1
    ! pz pz
    clm_d (10, 0, 0) = CLM1
    clm_d (10, 2, 0) = CLM43
    clm_sp(10, 0, 0) = CLM1
    clm_sp(10, 2, 0) = CLMzz
    alm   (10, 0, 0) = ALM1
    alm   (10, 2, 0) = ALMs45
    ! pz px
    clm_d (11, 2, 1) = CLM1
    clm_sp(11, 2, 1) = CLMzp
    alm   (11, 2, 1) = ALMs35
    ! pz py
    clm_d (12, 2,-1) = CLM1
    clm_sp(12, 2,-1) = CLMzp
    alm   (12, 2,-1) = ALMs35
    ! pz dz2
    clm_d (13, 1, 0) = CLMs43
    alm   (13, 1, 0) = ALMs45
    ! pz dzx
    clm_d (14, 1, 1) = CLM1
    alm   (14, 1, 1) = ALMs35
    ! pz dzy
    clm_d (15, 1,-1) = CLM1
    alm   (15, 1,-1) = ALMs35
    ! px px
    clm_d (18, 0, 0) = CLM1
    clm_d (18, 2, 0) = CLM23m
    clm_d (18, 2, 2) = CLM1
    clm_sp(18, 0, 0) = CLM1
    clm_sp(18, 2, 0) = CLMyy
    alm   (18, 0, 0) = ALM1
    alm   (18, 2, 0) = ALMs15m
    alm   (18, 2, 2) = ALMs35
    ! px py
    clm_d (19, 2,-2) = CLM1
    clm_sp(19, 2,-2) = CLMxy
    alm   (19, 2,-2) = ALMs35
    ! px dz2
    clm_d (20, 1, 1) = CLMs13m
    alm   (20, 1, 1) = ALMs15m
    ! px dzx
    clm_d (21, 1, 0) = CLM1
    alm   (21, 1, 0) = ALMs35
    ! px dx2-y2
    clm_d (23, 1, 1) = CLM1
    alm   (23, 1, 1) = ALMs35
    ! px dxy
    clm_d (24, 1,-1) = CLM1
    alm   (24, 1,-1) = ALMs35
    ! py py
    clm_d (25, 0, 0) = CLM1
    clm_d (25, 2, 0) = CLM23m
    clm_d (25, 2, 2) = CLM1m
    clm_sp(25, 0, 0) = CLM1
    clm_sp(25, 2, 0) = CLMxx
    alm   (25, 0, 0) = ALM1
    alm   (25, 2, 0) = ALMs15m
    alm   (25, 2, 2) = ALMs35m
    ! py dz2
    clm_d (26, 1,-1) = CLMs13m
    alm   (26, 1,-1) = ALMs15m
    ! py dzy
    clm_d (28, 1, 0) = CLM1
    alm   (28, 1, 0) = ALMs35
    ! py dx2-y2
    clm_d (29, 1,-1) = CLM1m
    alm   (29, 1,-1) = ALMs35m
    ! py dxy
    clm_d (30, 1, 1) = CLM1
    alm   (30, 1, 1) = ALMs35
    ! dz2 dz2
    clm_d (31, 0, 0) = CLM1
    clm_d (31, 2, 0) = CLM43
    alm   (31, 0, 0) = ALM1
    alm   (31, 2, 0) = ALMs20_49
    ! dz2 dzx
    clm_d (32, 2, 1) = CLMs13
    alm   (32, 2, 1) = ALMs5_49
    ! dz2 dzy
    clm_d (33, 2,-1) = CLMs13
    alm   (33, 2,-1) = ALMs5_49
    ! dz2 dx2-y2
    clm_d (34, 2, 2) = CLMs43m
    alm   (34, 2, 2) = ALMs20_49m
    ! dz2 dxy
    clm_d (35, 2,-2) = CLMs43m
    alm   (35, 2,-2) = ALMs20_49m
    ! dzx dzx
    clm_d (36, 0, 0) = CLM1
    clm_d (36, 2, 0) = CLM23
    clm_d (36, 2, 2) = CLM1
    alm   (36, 0, 0) = ALM1
    alm   (36, 2, 0) = ALMs5_49
    alm   (36, 2, 2) = ALMs15_49
    ! dzx dzy
    clm_d (37, 2,-2) = CLM1
    alm   (37, 2,-2) = ALMs15_49m
    ! dzx dxy-y2
    clm_d (38, 2, 1) = CLM1
    alm   (38, 2, 1) = ALMs15_49
    ! dzx dxy
    clm_d (39, 2,-1) = CLM1
    alm   (38, 2,-1) = ALMs15_49
    ! dzy dzy
    clm_d (40, 0, 0) = CLM1
    clm_d (40, 2, 0) = CLM23
    clm_d (40, 2, 2) = CLM1m
    alm   (40, 0, 0) = ALM1
    alm   (40, 2, 0) = ALMs5_49
    alm   (40, 2, 2) = ALMs5_49m
    ! dzy dx2-y2
    clm_d (41, 2,-1) = CLM1m
    alm   (41, 2,-1) = ALMs15_49m
    ! dzy dxy
    clm_d (42, 2, 1) = CLM1
    alm   (42, 2, 1) = ALMs15_49
    ! dx2-y2 dx2-y2
    clm_d (43, 0, 0) = CLM1
    clm_d (43, 2, 0) = CLM43m
    alm   (43, 0, 0) = ALM1
    alm   (43, 2, 0) = ALMs20_49m
    ! dxy dxy
    clm_d (45, 0, 0) = CLM1
    clm_d (45, 2, 0) = CLM43m
    alm   (45, 0, 0) = ALM1
    alm   (45, 2, 0) = ALMs20_49m
  END SUBROUTINE setup_clm_array

! *****************************************************************************
!> \brief Fills in the index number for the <ij|kl> integral as well as the
!>        symmetry index
!>
!> \param error ...
!> \date   03.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE setup_ijkl_array(error)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_ijkl_array', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    ! Address unique indexes (excluding those related by rotations)
    ! Indexes according:
    !           s   pz   px   py  dz2  dzx  dzy dx2-y2 dxy
    !      s    1    2    3    4    5    6    7    8    9
    !      pz   2   10   11   12   13   14   15   16   17
    !      px   3   11   18   19   20   21   22   23   24
    !      py   4   12   19   25   26   27   28   29   30
    !     dz2   5   13   20   26   31   32   33   34   35
    !     dzx   6   14   21   27   32   36   37   38   39
    !     dzy   7   15   22   28   33   37   40   41   42
    !  dx2-y2   8   16   23   29   34   38   41   43   44
    !     dxy   9   17   24   30   35   39   42   44   45
    ! ################  Zero the Arrays #######################
    ijkl_ind = 0
    ijkl_sym = 0
    ! ################      s       s   #######################
    !                       s       s   -        s       s
    ijkl_ind( 1,  1) =  1
    !                       s       s   -       pz       s
    ijkl_ind( 1,  2) =  2
    !                       s       s   -       pz      pz
    ijkl_ind( 1, 10) =  3
    !                       s       s   -       px      px
    ijkl_ind( 1, 18) =  4
    !                       s       s   -       py      py
    ijkl_ind( 1, 25) =  5
    ijkl_sym( 5)     =  4
    !                       s       s   -      dz2       s
    ijkl_ind( 1,  5) = 35
    !                       s       s   -      dz2      pz
    ijkl_ind( 1, 13) = 36
    !                       s       s   -      dz2     dz2
    ijkl_ind( 1, 31) = 37
    !                       s       s   -      dzx      px
    ijkl_ind( 1, 21) = 38
    !                       s       s   -      dzx     dzx
    ijkl_ind( 1, 36) = 39
    !                       s       s   -      dzy      py
    ijkl_ind( 1, 28) = 40
    ijkl_sym(40)     = 38
    !                       s       s   -      dzy     dzy
    ijkl_ind( 1, 40) = 41
    ijkl_sym(41)     = 39
    !                       s       s   -   dx2-y2  dx2-y2
    ijkl_ind( 1, 43) = 42
    !                       s       s   -      dxy     dxy
    ijkl_ind( 1, 45) = 43
    ijkl_sym(43)     = 42
    ! ################     pz       s   #######################
    !                      pz       s   -        s       s
    ijkl_ind( 2,  1) =  6
    !                      pz       s   -       pz       s
    ijkl_ind( 2,  2) =  7
    !                      pz       s   -       pz      pz
    ijkl_ind( 2, 10) =  8
    !                      pz       s   -       px      px
    ijkl_ind( 2, 18) =  9
    !                      pz       s   -       py      py
    ijkl_ind( 2, 25) = 10
    ijkl_sym(10)     =  9
    !                      pz       s   -      dz2       s
    ijkl_ind( 2,  5) = 44
    !                      pz       s   -      dz2      pz
    ijkl_ind( 2, 13) = 45
    !                      pz       s   -      dz2     dz2
    ijkl_ind( 2, 31) = 46
    !                      pz       s   -      dzx      px
    ijkl_ind( 2, 21) = 47
    !                      pz       s   -      dzx     dzx
    ijkl_ind( 2, 36) = 48
    !                      pz       s   -      dzy      py
    ijkl_ind( 2, 28) = 49
    ijkl_sym(49)     = 47
    !                      pz       s   -      dzy     dzy
    ijkl_ind( 2, 40) = 50
    ijkl_sym(50)     = 48
    !                      pz       s   -   dx2-y2  dx2-y2
    ijkl_ind( 2, 43) = 51
    !                      pz       s   -      dxy     dxy
    ijkl_ind( 2, 45) = 52
    ijkl_sym(52)     = 51
    ! ################     pz      pz   #######################
    !                      pz      pz   -        s       s
    ijkl_ind(10,  1) = 11
    !                      pz      pz   -       pz       s
    ijkl_ind(10,  2) = 12
    !                      pz      pz   -       pz      pz
    ijkl_ind(10, 10) = 13
    !                      pz      pz   -       px      px
    ijkl_ind(10, 18) = 14
    !                      pz      pz   -       py      py
    ijkl_ind(10, 25) = 15
    ijkl_sym(15)     = 14
    !                      pz      pz   -      dz2       s
    ijkl_ind(10,  5) = 53
    !                      pz      pz   -      dz2      pz
    ijkl_ind(10, 13) = 54
    !                      pz      pz   -      dz2     dz2
    ijkl_ind(10, 31) = 55
    !                      pz      pz   -      dzx      px
    ijkl_ind(10, 21) = 56
    !                      pz      pz   -      dzx     dzx
    ijkl_ind(10, 36) = 57
    !                      pz      pz   -      dzy      py
    ijkl_ind(10, 28) = 58
    ijkl_sym(58)     = 56
    !                      pz      pz   -      dzy     dzy
    ijkl_ind(10, 40) = 59
    ijkl_sym(59)     = 57
    !                      pz      pz   -   dx2-y2  dx2-y2
    ijkl_ind(10, 43) = 60
    !                      pz      pz   -      dxy     dxy
    ijkl_ind(10, 45) = 61
    ijkl_sym(61)     = 60
    ! ################     px       s   #######################
    !                      px       s   -       px       s
    ijkl_ind( 3,  3) = 16
    !                      px       s   -       px      pz
    ijkl_ind( 3, 11) = 17
    !                      px       s   -      dz2      px
    ijkl_ind( 3, 20) = 62
    !                      px       s   -      dzx       s
    ijkl_ind( 3,  6) = 63
    !                      px       s   -      dzx      pz
    ijkl_ind( 3, 14) = 64
    !                      px       s   -      dzx     dz2
    ijkl_ind( 3, 32) = 65
    !                      px       s   -   dx2-y2      px
    ijkl_ind( 3, 23) = 66
    !                      px       s   -   dx2-y2     dzx
    ijkl_ind( 3, 38) = 67
    !                      px       s   -      dxy      py
    ijkl_ind( 3, 30) = 68
    ijkl_sym(68)     = 66
    !                      px       s   -      dxy     dzy
    ijkl_ind( 3, 42) = 69
    ijkl_sym(69)     = 67
    ! ################     px      pz   #######################
    !                      px      pz   -       px       s
    ijkl_ind(11,  3) = 18
    !                      px      pz   -       px      pz
    ijkl_ind(11, 11) = 19
    !                      px      pz   -      dz2      px
    ijkl_ind(11, 20) = 70
    !                      px      pz   -      dzx       s
    ijkl_ind(11,  6) = 71
    !                      px      pz   -      dzx      pz
    ijkl_ind(11, 14) = 72
    !                      px      pz   -      dzx     dz2
    ijkl_ind(11, 32) = 73
    !                      px      pz   -   dx2-y2      px
    ijkl_ind(11, 23) = 74
    !                      px      pz   -   dx2-y2     dzx
    ijkl_ind(11, 38) = 75
    !                      px      pz   -      dxy      py
    ijkl_ind(11, 30) = 76
    ijkl_sym(76)     = 74
    !                      px      pz   -      dxy     dzy
    ijkl_ind(11, 42) = 77
    ijkl_sym(77)     = 75
    ! ################     px      px   #######################
    !                      px      px   -        s       s
    ijkl_ind(18,  1) = 20
    !                      px      px   -       pz       s
    ijkl_ind(18,  2) = 21
    !                      px      px   -       pz      pz
    ijkl_ind(18, 10) = 22
    !                      px      px   -       px      px
    ijkl_ind(18, 18) = 23
    !                      px      px   -       py      py
    ijkl_ind(18, 25) = 24
    !                      px      px   -      dz2       s
    ijkl_ind(18,  5) = 78
    !                      px      px   -      dz2      pz
    ijkl_ind(18, 13) = 79
    !                      px      px   -      dz2     dz2
    ijkl_ind(18, 31) = 80
    !                      px      px   -      dzx      px
    ijkl_ind(18, 21) = 81
    !                      px      px   -      dzx     dzx
    ijkl_ind(18, 36) = 82
    !                      px      px   -      dzy      py
    ijkl_ind(18, 28) = 83
    !                      px      px   -      dzy     dzy
    ijkl_ind(18, 40) = 84
    !                      px      px   -   dx2-y2       s
    ijkl_ind(18,  8) = 85
    !                      px      px   -   dx2-y2      pz
    ijkl_ind(18, 16) = 86
    !                      px      px   -   dx2-y2     dz2
    ijkl_ind(18, 34) = 87
    !                      px      px   -   dx2-y2  dx2-y2
    ijkl_ind(18, 43) = 88
    !                      px      px   -      dxy     dxy
    ijkl_ind(18, 45) = 89
    ijkl_sym(89)     = 88
    ! ################     py       s   #######################
    !                      py       s   -       py       s
    ijkl_ind( 4,  4) = 25
    ijkl_sym(25)     = 16
    !                      py       s   -       py      pz
    ijkl_ind( 4, 12) = 26
    ijkl_sym(26)     = 17
    !                      py       s   -      dz2      py
    ijkl_ind( 4, 26) = 90
    ijkl_sym(90)     = 62
    !                      py       s   -      dzy       s
    ijkl_ind( 4,  7) = 91
    ijkl_sym(91)     = 63
    !                      py       s   -      dzy      pz
    ijkl_ind( 4, 15) = 92
    ijkl_sym(92)     = 64
    !                      py       s   -      dzy     dz2
    ijkl_ind( 4, 33) = 93
    ijkl_sym(93)     = 65
    !                      py       s   -   dx2-y2      py
    ijkl_ind( 4, 29) = 94
    ijkl_sym(94)     = 66*(-1)
    !                      py       s   -   dx2-y2     dzy
    ijkl_ind( 4, 41) = 95
    ijkl_sym(95)     = 67*(-1)
    !                      py       s   -      dxy      px
    ijkl_ind( 4, 24) = 96
    ijkl_sym(96)     = 66
    !                      py       s   -      dxy     dzx
    ijkl_ind( 4, 39) = 97
    ijkl_sym(97)     = 67
    ! ################     py      pz   #######################
    !                      py      pz   -       py       s
    ijkl_ind(12,  4) = 27
    ijkl_sym(27)     = 18
    !                      py      pz   -       py      pz
    ijkl_ind(12, 12) = 28
    ijkl_sym(28)     = 19
    !                      py      pz   -      dzy       s
    ijkl_ind(12, 26) = 98
    ijkl_sym(98)     = 70
    !                      py      pz   -      dzy      pz
    ijkl_ind(12,  7) = 99
    ijkl_sym(99)     = 71
    !                      py      pz   -      dzy     dz2
    ijkl_ind(12, 15) = 100
    ijkl_sym(100)    = 72
    !                      py      pz   -      dzy     dz2
    ijkl_ind(12, 33) = 101
    ijkl_sym(101)    = 73
    !                      py      pz   -   dx2-y2      py
    ijkl_ind(12, 29) = 102
    ijkl_sym(102)    = 74*(-1)
    !                      py      pz   -   dx2-y2     dzy
    ijkl_ind(12, 41) = 103
    ijkl_sym(103)    = 75*(-1)
    !                      py      pz   -      dxy      px
    ijkl_ind(12, 24) = 104
    ijkl_sym(104)    = 74
    !                      py      pz   -      dxy     dzx
    ijkl_ind(12, 39) = 105
    ijkl_sym(105)    = 75
    ! ################     py      px   #######################
    !                      py      px   -       py      px
    ijkl_ind(19, 19) = 29
    !                      py      px   -      dzx      py
    ijkl_ind(19, 27) = 106
    ijkl_sym(106)    = 86
    !                      py      px   -      dzy      px
    ijkl_ind(19, 22) = 107
    ijkl_sym(107)    = 86
    !                      py      px   -      dzy     dzx
    ijkl_ind(19, 37) = 108
    !                      py      px   -      dxy       s
    ijkl_ind(19,  9) = 109
    ijkl_sym(109)    = 85
    !                      py      px   -      dxy      pz
    ijkl_ind(19, 17) = 110
    ijkl_sym(110)    = 86
    !                      py      px   -      dxy     dz2
    ijkl_ind(19, 35) = 111
    ijkl_sym(111)    = 87
    ! ################     py      py   #######################
    !                      py      py   -        s       s
    ijkl_ind(25,  1) = 30
    ijkl_sym(30)     = 20
    !                      py      py   -       pz       s
    ijkl_ind(25,  2) = 31
    ijkl_sym(31)     = 21
    !                      py      py   -       pz      pz
    ijkl_ind(25, 10) = 32
    ijkl_sym(32)     = 22
    !                      py      py   -       px      px
    ijkl_ind(25, 18) = 33
    ijkl_sym(33)     = 24
    !                      py      py   -       py      py
    ijkl_ind(25, 25) = 34
    ijkl_sym(34)     = 23
    !                      py      py   -      dz2       s
    ijkl_ind(25,  5) = 112
    ijkl_sym(112)    = 78
    !                      py      py   -      dz2      pz
    ijkl_ind(25, 13) = 113
    ijkl_sym(113)    = 79
    !                      py      py   -      dz2     dz2
    ijkl_ind(25, 31) = 114
    ijkl_sym(114)    = 80
    !                      py      py   -      dzx      px
    ijkl_ind(25, 21) = 115
    ijkl_sym(115)    = 83
    !                      py      py   -      dzx     dzx
    ijkl_ind(25, 36) = 116
    ijkl_sym(116)    = 84
    !                      py      py   -      dzy      py
    ijkl_ind(25, 28) = 117
    ijkl_sym(117)    = 81
    !                      py      py   -      dzy     dzy
    ijkl_ind(25, 40) = 118
    ijkl_sym(118)    = 82
    !                      py      py   -   dx2-y2       s
    ijkl_ind(25,  8) = 119
    ijkl_sym(119)    = 85*(-1)
    !                      py      py   -   dx2-y2      pz
    ijkl_ind(25, 16) = 120
    ijkl_sym(120)    = 86*(-1)
    !                      py      py   -   dx2-y2     dz2
    ijkl_ind(25, 34) = 121
    ijkl_sym(121)    = 87*(-1)
    !                      py      py   -   dx2-y2  dx2-y2
    ijkl_ind(25, 43) = 122
    ijkl_sym(122)    = 88
    !                      py      py   -      dxy     dxy
    ijkl_ind(25, 45) = 123
    ijkl_sym(123)    = 88
    ! ################    dz2     dz2   #######################
    !                     dz2       s   -        s       s
    ijkl_ind( 5,  1) = 124
    !                     dz2       s   -       pz       s
    ijkl_ind( 5,  2) = 125
    !                     dz2       s   -       pz      pz
    ijkl_ind( 5, 10) = 126
    !                     dz2       s   -       px      px
    ijkl_ind( 5, 18) = 127
    !                     dz2       s   -       py      py
    ijkl_ind( 5, 25) = 128
    ijkl_sym(128)    = 127
    !                     dz2       s   -      dz2       s
    ijkl_ind( 5,  5) = 129
    !                     dz2       s   -      dz2      pz
    ijkl_ind( 5, 13) = 130
    !                     dz2       s   -      dz2     dz2
    ijkl_ind( 5, 31) = 131
    !                     dz2       s   -      dzx      px
    ijkl_ind( 5, 21) = 132
    !                     dz2       s   -      dzx     dzx
    ijkl_ind( 5, 36) = 133
    !                     dz2       s   -      dzy      py
    ijkl_ind( 5, 28) = 134
    ijkl_sym(134)    = 132
    !                     dz2       s   -      dzy     dzy
    ijkl_ind( 5, 40) = 135
    ijkl_sym(135)    = 133
    !                     dz2       s   -   dx2-y2  dx2-y2
    ijkl_ind( 5, 43) = 136
    !                     dz2       s   -      dxy     dxy
    ijkl_ind( 5, 45) = 137
    ijkl_sym(137)    = 136
    ! ################    dz2      pz   #######################
    !                     dz2      pz   -        s       s
    ijkl_ind(13,  1) = 138
    !                     dz2      pz   -       pz       s
    ijkl_ind(13,  2) = 139
    !                     dz2      pz   -       pz      pz
    ijkl_ind(13, 10) = 140
    !                     dz2      pz   -       px      px
    ijkl_ind(13, 18) = 141
    !                     dz2      pz   -       py      py
    ijkl_ind(13, 25) = 142
    ijkl_sym(142)    = 141
    !                     dz2      pz   -      dz2       s
    ijkl_ind(13,  5) = 143
    !                     dz2      pz   -      dz2      pz
    ijkl_ind(13, 13) = 144
    !                     dz2      pz   -      dz2     dz2
    ijkl_ind(13, 31) = 145
    !                     dz2      pz   -      dzx      px
    ijkl_ind(13, 21) = 146
    !                     dz2      pz   -      dzx     dzx
    ijkl_ind(13, 36) = 147
    !                     dz2      pz   -      dzy      py
    ijkl_ind(13, 28) = 148
    ijkl_sym(148)    = 146
    !                     dz2      pz   -      dzy     dzy
    ijkl_ind(13, 40) = 149
    ijkl_sym(149)    = 147
    !                     dz2      pz   -   dx2-y2  dx2-y2
    ijkl_ind(13, 43) = 150
    !                     dz2      pz   -      dxy     dxy
    ijkl_ind(13, 45) = 151
    ijkl_sym(151)    = 150
    ! ################    dz2      px   #######################
    !                     dz2      px   -       px       s
    ijkl_ind(20,  3) = 152
    !                     dz2      px   -       px      pz
    ijkl_ind(20, 11) = 153
    !                     dz2      px   -      dz2      px
    ijkl_ind(20, 20) = 154
    !                     dz2      px   -      dzx       s
    ijkl_ind(20,  6) = 155
    !                     dz2      px   -      dzx      pz
    ijkl_ind(20, 14) = 156
    !                     dz2      px   -      dzx     dz2
    ijkl_ind(20, 32) = 157
    !                     dz2      px   -   dx2-y2      px
    ijkl_ind(20, 23) = 158
    !                     dz2      px   -   dx2-y2     dzx
    ijkl_ind(20, 38) = 159
    !                     dz2      px   -      dxy      py
    ijkl_ind(20, 30) = 160
    ijkl_sym(160)    = 158
    !                     dz2      px   -      dxy     dxy
    ijkl_ind(20, 42) = 161
    ijkl_sym(161)    = 159
    ! ################    dz2      py   #######################
    !                     dz2      py   -       py       s
    ijkl_ind(26,  4) = 162
    ijkl_sym(162)    = 152
    !                     dz2      py   -       py      pz
    ijkl_ind(26, 12) = 163
    ijkl_sym(163)    = 153
    !                     dz2      py   -      dz2      py
    ijkl_ind(26, 26) = 164
    ijkl_sym(164)    = 154
    !                     dz2      py   -      dzy       s
    ijkl_ind(26,  7) = 165
    ijkl_sym(165)    = 155
    !                     dz2      py   -      dzy      pz
    ijkl_ind(26, 15) = 166
    ijkl_sym(166)    = 156
    !                     dz2      py   -      dzy     dz2
    ijkl_ind(26, 33) = 167
    ijkl_sym(167)    = 157
    !                     dz2      py   -   dx2-y2      py
    ijkl_ind(26, 29) = 168
    ijkl_sym(168)    = 158*(-1)
    !                     dz2      py   -   dx2-y2     dzy
    ijkl_ind(26, 41) = 169
    ijkl_sym(169)    = 159*(-1)
    !                     dz2      py   -      dxy      px
    ijkl_ind(26, 24) = 170
    ijkl_sym(170)    = 158
    !                     dz2      py   -      dxy     dzx
    ijkl_ind(26, 39) = 171
    ijkl_sym(171)    = 159
    ! ################    dz2     dz2   #######################
    !                     dz2     dz2   -        s       s
    ijkl_ind(31,  1) = 172
    !                     dz2     dz2   -       pz       s
    ijkl_ind(31,  2) = 173
    !                     dz2     dz2   -       pz      pz
    ijkl_ind(31, 10) = 174
    !                     dz2     dz2   -       px      px
    ijkl_ind(31, 18) = 175
    !                     dz2     dz2   -       py      py
    ijkl_ind(31, 25) = 176
    ijkl_sym(176)    = 175
    !                     dz2     dz2   -      dz2       s
    ijkl_ind(31,  5) = 177
    !                     dz2     dz2   -      dz2      pz
    ijkl_ind(31, 13) = 178
    !                     dz2     dz2   -      dz2     dz2
    ijkl_ind(31, 31) = 179
    !                     dz2     dz2   -      dzx      px
    ijkl_ind(31, 21) = 180
    !                     dz2     dz2   -      dzx     dzx
    ijkl_ind(31, 36) = 181
    !                     dz2     dz2   -      dzy      py
    ijkl_ind(31, 28) = 182
    ijkl_sym(182)    = 180
    !                     dz2     dz2   -      dzy     dzy
    ijkl_ind(31, 40) = 183
    ijkl_sym(183)    = 181
    !                     dz2     dz2   -   dx2-y2  dx2-y2
    ijkl_ind(31, 43) = 184
    !                     dz2     dz2   -      dxy     dxy
    ijkl_ind(31, 45) = 185
    ijkl_sym(185)    = 184
    ! ################    dzx       s   #######################
    !                     dzx       s   -       px       s
    ijkl_ind( 6,  3) = 186
    !                     dzx       s   -       px      pz
    ijkl_ind( 6, 11) = 187
    !                     dzx       s   -      dz2      px
    ijkl_ind( 6, 20) = 188
    !                     dzx       s   -      dzx       s
    ijkl_ind( 6,  6) = 189
    !                     dzx       s   -      dzx      pz
    ijkl_ind( 6, 14) = 190
    !                     dzx       s   -      dzx     dz2
    ijkl_ind( 6, 32) = 191
    !                     dzx       s   -   dx2-y2      px
    ijkl_ind( 6, 23) = 192
    !                     dzx       s   -   dx2-y2  dx2-y2
    ijkl_ind( 6, 38) = 193
    !                     dzx       s   -      dxy      py
    ijkl_ind( 6, 30) = 194
    ijkl_sym(194)    = 192
    !                     dzx       s   -      dxy     dzy
    ijkl_ind( 6, 42) = 195
    ijkl_sym(195)    = 193
    ! ################    dzx      pz   #######################
    !                     dzx      pz   -       px       s
    ijkl_ind(14,  3) = 196
    !                     dzx      pz   -       px      pz
    ijkl_ind(14, 11) = 197
    !                     dzx      pz   -      dz2      px
    ijkl_ind(14, 20) = 198
    !                     dzx      pz   -      dzx       s
    ijkl_ind(14,  6) = 199
    !                     dzx      pz   -      dzx      pz
    ijkl_ind(14, 14) = 200
    !                     dzx      pz   -      dzx     dz2
    ijkl_ind(14, 32) = 201
    !                     dzx      pz   -   dx2-y2      px
    ijkl_ind(14, 23) = 202
    !                     dzx      pz   -   dx2-y2  dx2-y2
    ijkl_ind(14, 38) = 203
    !                     dzx      pz   -      dxy      py
    ijkl_ind(14, 30) = 204
    ijkl_sym(204)    = 202
    !                     dzx      pz   -      dxy     dzy
    ijkl_ind(14, 42) = 205
    ijkl_sym(205)    = 203
    ! ################    dzx      px   #######################
    !                     dzx      px   -        s       s
    ijkl_ind(21,  1) = 206
    !                     dzx      px   -       pz       s
    ijkl_ind(21,  2) = 207
    !                     dzx      px   -       pz      pz
    ijkl_ind(21, 10) = 208
    !                     dzx      px   -       px      px
    ijkl_ind(21, 18) = 209
    !                     dzx      px   -       py      py
    ijkl_ind(21, 25) = 210
    !                     dzx      px   -      dz2       s
    ijkl_ind(21,  5) = 211
    !                     dzx      px   -      dz2      pz
    ijkl_ind(21, 13) = 212
    !                     dzx      px   -      dz2     dz2
    ijkl_ind(21, 31) = 213
    !                     dzx      px   -      dzx      px
    ijkl_ind(21, 21) = 214
    !                     dzx      px   -      dzx     dzx
    ijkl_ind(21, 36) = 215
    !                     dzx      px   -      dzy      py
    ijkl_ind(21, 28) = 216
    !                     dzx      px   -      dzy     dzy
    ijkl_ind(21, 40) = 217
    !                     dzx      px   -   dx2-y2       s
    ijkl_ind(21,  8) = 218
    !                     dzx      px   -   dx2-y2      pz
    ijkl_ind(21, 16) = 219
    !                     dzx      px   -   dx2-y2     dz2
    ijkl_ind(21, 34) = 220
    !                     dzx      px   -   dx2-y2  dx2-y2
    ijkl_ind(21, 43) = 221
    !                     dzx      px   -      dxy     dxy
    ijkl_ind(21, 45) = 222
    ijkl_sym(222)    = 221
    ! ################    dzx      py   #######################
    !                     dzx      py   -       px      py
    ijkl_ind(27, 19) = 223
    !                     dzx      py   -      dzx      py
    ijkl_ind(27, 27) = 224
    ijkl_sym(224)    = 219
    !                     dzx      py   -      dzy      px
    ijkl_ind(27, 22) = 225
    ijkl_sym(225)    = 219
    !                     dzx      py   -      dzy     dzx
    ijkl_ind(27, 37) = 226
    !                     dzx      py   -      dxy       s
    ijkl_ind(27,  9) = 227
    ijkl_sym(227)    = 218
    !                     dzx      py   -      dxy      pz
    ijkl_ind(27, 17) = 228
    ijkl_sym(228)    = 219
    !                     dzx      py   -      dxy     dz2
    ijkl_ind(27, 35) = 229
    ijkl_sym(229)    = 220
    ! ################    dzx     dz2   #######################
    !                     dzx     dz2   -       px       s
    ijkl_ind(32,  3) = 230
    !                     dzx     dz2   -       px      pz
    ijkl_ind(32, 11) = 231
    !                     dzx     dz2   -      dz2      px
    ijkl_ind(32, 20) = 232
    !                     dzx     dz2   -      dzx       s
    ijkl_ind(32,  6) = 233
    !                     dzx     dz2   -      dzx      pz
    ijkl_ind(32, 14) = 234
    !                     dzx     dz2   -      dzx     dz2
    ijkl_ind(32, 32) = 235
    !                     dzx     dz2   -   dx2-y2      px
    ijkl_ind(32, 23) = 236
    !                     dzx     dz2   -   dx2-y2  dx2-y2
    ijkl_ind(32, 38) = 237
    !                     dzx     dz2   -      dxy      py
    ijkl_ind(32, 30) = 238
    ijkl_sym(238)    = 236
    !                     dzx     dz2   -      dxy     dzy
    ijkl_ind(32, 42) = 239
    ijkl_sym(239)    = 237
    ! ################    dzx     dzx   #######################
    !                     dzx     dzx   -        s       s
    ijkl_ind(36,  1) = 240
    !                     dzx     dzx   -       pz       s
    ijkl_ind(36,  2) = 241
    !                     dzx     dzx   -       pz      pz
    ijkl_ind(36, 10) = 242
    !                     dzx     dzx   -       px      px
    ijkl_ind(36, 18) = 243
    !                     dzx     dzx   -       py      py
    ijkl_ind(36, 25) = 244
    !                     dzx     dzx   -      dz2       s
    ijkl_ind(36,  5) = 245
    !                     dzx     dzx   -      dz2      pz
    ijkl_ind(36, 13) = 246
    !                     dzx     dzx   -      dz2     dz2
    ijkl_ind(36, 31) = 247
    !                     dzx     dzx   -      dzx      px
    ijkl_ind(36, 21) = 248
    !                     dzx     dzx   -      dzx     dzx
    ijkl_ind(36, 36) = 249
    !                     dzx     dzx   -      dzy      py
    ijkl_ind(36, 28) = 250
    !                     dzx     dzx   -      dzy     dzy
    ijkl_ind(36, 40) = 251
    !                     dzx     dzx   -   dx2-y2       s
    ijkl_ind(36,  8) = 252
    !                     dzx     dzx   -   dx2-y2      pz
    ijkl_ind(36, 16) = 253
    !                     dzx     dzx   -   dx2-y2     dz2
    ijkl_ind(36, 34) = 254
    !                     dzx     dzx   -   dx2-y2  dx2-y2
    ijkl_ind(36, 43) = 255
    !                     dzx     dzx   -      dxy     dxy
    ijkl_ind(36, 45) = 256
    ijkl_sym(256)    = 255
    ! ################    dzy       s   #######################
    !                     dzy       s   -       py       s
    ijkl_ind( 7,  4) = 257
    ijkl_sym(257)    = 186
    !                     dzy       s   -       py      pz
    ijkl_ind( 7, 12) = 258
    ijkl_sym(258)    = 187
    !                     dzy       s   -      dz2      py
    ijkl_ind( 7, 26) = 259
    ijkl_sym(259)    = 188
    !                     dzy       s   -      dzy       s
    ijkl_ind( 7,  7) = 260
    ijkl_sym(260)    = 189
    !                     dzy       s   -      dzy      pz
    ijkl_ind( 7, 15) = 261
    ijkl_sym(261)    = 190
    !                     dzy       s   -      dzy     dz2
    ijkl_ind( 7, 33) = 262
    ijkl_sym(262)    = 191
    !                     dzy       s   -   dx2-y2      py
    ijkl_ind( 7, 29) = 263
    ijkl_sym(263)    = 192*(-1)
    !                     dzy       s   -   dx2-y2     dzy
    ijkl_ind( 7, 41) = 264
    ijkl_sym(264)    = 193*(-1)
    !                     dzy       s   -      dxy      px
    ijkl_ind( 7, 24) = 265
    ijkl_sym(265)    = 192
    !                     dzy       s   -      dxy     dzx
    ijkl_ind( 7, 39) = 266
    ijkl_sym(266)    = 193
    ! ################    dzy      pz   #######################
    !                     dzy      pz   -       py       s
    ijkl_ind(15,  4) = 267
    ijkl_sym(267)    = 196
    !                     dzy      pz   -       py      pz
    ijkl_ind(15, 12) = 268
    ijkl_sym(268)    = 197
    !                     dzy      pz   -      dz2      py
    ijkl_ind(15, 26) = 269
    ijkl_sym(269)    = 198
    !                     dzy      pz   -      dzy       s
    ijkl_ind(15,  7) = 270
    ijkl_sym(270)    = 199
    !                     dzy      pz   -      dzy      pz
    ijkl_ind(15, 15) = 271
    ijkl_sym(271)    = 200
    !                     dzy      pz   -      dzy     dz2
    ijkl_ind(15, 33) = 272
    ijkl_sym(272)    = 201
    !                     dzy      pz   -   dx2-y2      py
    ijkl_ind(15, 29) = 273
    ijkl_sym(273)    = 202*(-1)
    !                     dzy      pz   -   dx2-y2     dzy
    ijkl_ind(15, 41) = 274
    ijkl_sym(274)    = 203*(-1)
    !                     dzy      pz   -      dxy      px
    ijkl_ind(15, 24) = 275
    ijkl_sym(275)    = 202
    !                     dzy      pz   -      dxy     dzx
    ijkl_ind(15, 39) = 276
    ijkl_sym(276)    = 203
    ! ################    dzy      px   #######################
    !                     dzy      px   -       px      py
    ijkl_ind(22, 19) = 277
    ijkl_sym(277)    = 223
    !                     dzy      px   -      dzx      py
    ijkl_ind(22, 27) = 278
    ijkl_sym(278)    = 219
    !                     dzy      px   -      dzy      px
    ijkl_ind(22, 22) = 279
    ijkl_sym(279)    = 219
    !                     dzy      px   -      dzy     dzx
    ijkl_ind(22, 37) = 280
    ijkl_sym(280)    = 226
    !                     dzy      px   -      dxy       s
    ijkl_ind(22,  9) = 281
    ijkl_sym(281)    = 218
    !                     dzy      px   -      dxy      pz
    ijkl_ind(22, 17) = 282
    ijkl_sym(282)    = 219
    !                     dzy      px   -      dxy     dz2
    ijkl_ind(22, 35) = 283
    ijkl_sym(283)    = 220
    ! ################    dzy      py   #######################
    !                     dzy      py   -        s       s
    ijkl_ind(28,  1) = 284
    ijkl_sym(284)    = 206
    !                     dzy      py   -       pz       s
    ijkl_ind(28,  2) = 285
    ijkl_sym(285)    = 207
    !                     dzy      py   -       pz      pz
    ijkl_ind(28, 10) = 286
    ijkl_sym(286)    = 208
    !                     dzy      py   -       px      px
    ijkl_ind(28, 18) = 287
    ijkl_sym(287)    = 210
    !                     dzy      py   -       py      py
    ijkl_ind(28, 25) = 288
    ijkl_sym(288)    = 209
    !                     dzy      py   -      dz2       s
    ijkl_ind(28,  5) = 289
    ijkl_sym(289)    = 211
    !                     dzy      py   -      dz2      pz
    ijkl_ind(28, 13) = 290
    ijkl_sym(290)    = 212
    !                     dzy      py   -      dz2     dz2
    ijkl_ind(28, 31) = 291
    ijkl_sym(291)    = 213
    !                     dzy      py   -      dzx      px
    ijkl_ind(28, 21) = 292
    ijkl_sym(292)    = 216
    !                     dzy      py   -      dzx     dzx
    ijkl_ind(28, 36) = 293
    ijkl_sym(293)    = 217
    !                     dzy      py   -      dzy      py
    ijkl_ind(28, 28) = 294
    ijkl_sym(294)    = 214
    !                     dzy      py   -      dzy     dzy
    ijkl_ind(28, 40) = 295
    ijkl_sym(295)    = 215
    !                     dzy      py   -   dx2-y2       s
    ijkl_ind(28,  8) = 296
    ijkl_sym(296)    = 218*(-1)
    !                     dzy      py   -   dx2-y2      pz
    ijkl_ind(28, 16) = 297
    ijkl_sym(297)    = 219*(-1)
    !                     dzy      py   -   dx2-y2     dz2
    ijkl_ind(28, 34) = 298
    ijkl_sym(298)    = 220*(-1)
    !                     dzy      py   -   dx2-y2  dx2-y2
    ijkl_ind(28, 43) = 299
    ijkl_sym(299)    = 221
    !                     dzy      py   -      dxy     dxy
    ijkl_ind(28, 45) = 300
    ijkl_sym(300)    = 221
    ! ################    dzy     dz2   #######################
    !                     dzy     dz2   -       py       s
    ijkl_ind(33,  4) = 301
    ijkl_sym(301)    = 230
    !                     dzy     dz2   -       py      pz
    ijkl_ind(33, 12) = 302
    ijkl_sym(302)    = 231
    !                     dzy     dz2   -      dz2      py
    ijkl_ind(33, 26) = 303
    ijkl_sym(303)    = 232
    !                     dzy     dz2   -      dzy       s
    ijkl_ind(33,  7) = 304
    ijkl_sym(304)    = 233
    !                     dzy     dz2   -      dzy      pz
    ijkl_ind(33, 15) = 305
    ijkl_sym(305)    = 234
    !                     dzy     dz2   -      dzy     dz2
    ijkl_ind(33, 33) = 306
    ijkl_sym(306)    = 235
    !                     dzy     dz2   -   dx2-y2      py
    ijkl_ind(33, 29) = 307
    ijkl_sym(307)    = 236*(-1)
    !                     dzy     dz2   -   dx2-y2     dzy
    ijkl_ind(33, 41) = 308
    ijkl_sym(308)    = 237*(-1)
    !                     dzy     dz2   -      dxy      px
    ijkl_ind(33, 24) = 309
    ijkl_sym(309)    = 236
    !                     dzy     dz2   -      dxy     dzx
    ijkl_ind(33, 39) = 310
    ijkl_sym(310)    = 237
    ! ################    dzy     dzx   #######################
    !                     dzy     dzx   -       px      py
    ijkl_ind(37, 19) = 311
    !                     dzy     dzx   -      dzx      py
    ijkl_ind(37, 27) = 312
    ijkl_sym(312)    = 253
    !                     dzy     dzx   -      dzy      px
    ijkl_ind(37, 22) = 313
    ijkl_sym(313)    = 253
    !                     dzy     dzx   -      dzy     dzx
    ijkl_ind(37, 37) = 314
    !                     dzy     dzx   -      dxy       s
    ijkl_ind(37,  9) = 315
    ijkl_sym(315)    = 252
    !                     dzy     dzx   -      dxy      pz
    ijkl_ind(37, 17) = 316
    ijkl_sym(316)    = 253
    !                     dzy     dzx   -      dxy     dz2
    ijkl_ind(37, 35) = 317
    ijkl_sym(317)    = 254
    ! ################    dzy     dzy   #######################
    !                     dzy     dzy   -        s       s
    ijkl_ind(40,  1) = 318
    ijkl_sym(318)    = 240
    !                     dzy     dzy   -       pz       s
    ijkl_ind(40,  2) = 319
    ijkl_sym(319)    = 241
    !                     dzy     dzy   -       pz      pz
    ijkl_ind(40, 10) = 320
    ijkl_sym(320)    = 242
    !                     dzy     dzy   -       px      px
    ijkl_ind(40, 18) = 321
    ijkl_sym(321)    = 244
    !                     dzy     dzy   -       py      py
    ijkl_ind(40, 25) = 322
    ijkl_sym(322)    = 243
    !                     dzy     dzy   -      dz2       s
    ijkl_ind(40,  5) = 323
    ijkl_sym(323)    = 245
    !                     dzy     dzy   -      dz2      pz
    ijkl_ind(40, 13) = 324
    ijkl_sym(324)    = 246
    !                     dzy     dzy   -      dz2     dz2
    ijkl_ind(40, 31) = 325
    ijkl_sym(325)    = 247
    !                     dzy     dzy   -      dzx      px
    ijkl_ind(40, 21) = 326
    ijkl_sym(326)    = 250
    !                     dzy     dzy   -      dzx     dzx
    ijkl_ind(40, 36) = 327
    ijkl_sym(327)    = 251
    !                     dzy     dzy   -      dzy      py
    ijkl_ind(40, 28) = 328
    ijkl_sym(328)    = 248
    !                     dzy     dzy   -      dzy     dzy
    ijkl_ind(40, 40) = 329
    ijkl_sym(329)    = 249
    !                     dzy     dzy   -   dx2-y2       s
    ijkl_ind(40,  8) = 330
    ijkl_sym(330)    = 252*(-1)
    !                     dzy     dzy   -   dx2-y2      pz
    ijkl_ind(40, 16) = 331
    ijkl_sym(331)    = 253*(-1)
    !                     dzy     dzy   -   dx2-y2     dz2
    ijkl_ind(40, 34) = 332
    ijkl_sym(332)    = 254*(-1)
    !                     dzy     dzy   -   dx2-y2  dx2-y2
    ijkl_ind(40, 43) = 333
    ijkl_sym(333)    = 255
    !                     dzy     dzy   -      dxy     dxy
    ijkl_ind(40, 45) = 334
    ijkl_sym(334)    = 255
    ! ################ dx2-y2       s   #######################
    !                  dx2-y2       s   -       px      px
    ijkl_ind( 8, 18) = 335
    !                  dx2-y2       s   -       py      py
    ijkl_ind( 8, 25) = 336
    ijkl_sym(336)    = 335*(-1)
    !                  dx2-y2       s   -      dzx      px
    ijkl_ind( 8, 21) = 337
    !                  dx2-y2       s   -      dzx     dzx
    ijkl_ind( 8, 36) = 338
    !                  dx2-y2       s   -      dzy      py
    ijkl_ind( 8, 28) = 339
    ijkl_sym(339)    = 337*(-1)
    !                  dx2-y2       s   -      dzy     dzy
    ijkl_ind( 8, 40) = 340
    ijkl_sym(340)    = 338*(-1)
    !                  dx2-y2       s   -   dx2-y2       s
    ijkl_ind( 8,  8) = 341
    !                  dx2-y2       s   -   dx2-y2      pz
    ijkl_ind( 8, 16) = 342
    ijkl_sym(342)    = 337
    !                  dx2-y2       s   -   dx2-y2     dz2
    ijkl_ind( 8, 34) = 343
    ! ################ dx2-y2      pz   #######################
    !                  dx2-y2      pz   -       px      px
    ijkl_ind(16, 18) = 344
    ijkl_sym(344)    = 223
    !                  dx2-y2      pz   -       py      py
    ijkl_ind(16, 25) = 345
    ijkl_sym(345)    = 223*(-1)
    !                  dx2-y2      pz   -      dzx      px
    ijkl_ind(16, 21) = 346
    ijkl_sym(346)    = 219
    !                  dx2-y2      pz   -      dzx     dzx
    ijkl_ind(16, 36) = 347
    ijkl_sym(347)    = 226
    !                  dx2-y2      pz   -      dzy      py
    ijkl_ind(16, 28) = 348
    ijkl_sym(348)    = 219*(-1)
    !                  dx2-y2      pz   -      dzy     dzy
    ijkl_ind(16, 40) = 349
    ijkl_sym(349)    = 226*(-1)
    !                  dx2-y2      pz   -   dx2-y2       s
    ijkl_ind(16,  8) = 350
    ijkl_sym(350)    = 218
    !                  dx2-y2      pz   -   dx2-y2      pz
    ijkl_ind(16, 16) = 351
    ijkl_sym(351)    = 219
    !                  dx2-y2      pz   -   dx2-y2     dz2
    ijkl_ind(16, 34) = 352
    ijkl_sym(352)    = 220
    ! ################ dx2-y2      px   #######################
    !                  dx2-y2      px   -       px       s
    ijkl_ind(23,  3) = 353
    !                  dx2-y2      px   -       px      pz
    ijkl_ind(23, 11) = 354
    !                  dx2-y2      px   -      dz2      px
    ijkl_ind(23, 20) = 355
    !                  dx2-y2      px   -      dzx       s
    ijkl_ind(23,  6) = 356
    !                  dx2-y2      px   -      dzx      pz
    ijkl_ind(23, 14) = 357
    !                  dx2-y2      px   -      dzx     dz2
    ijkl_ind(23, 32) = 358
    !                  dx2-y2      px   -   dx2-y2      px
    ijkl_ind(23, 23) = 359
    !                  dx2-y2      px   -   dx2-y2  dx2-y2
    ijkl_ind(23, 38) = 360
    !                  dx2-y2      px   -      dxy      py
    ijkl_ind(23, 30) = 361
    !                  dx2-y2      px   -      dxy     dzy
    ijkl_ind(23, 42) = 362
    ! ################ dx2-y2      py   #######################
    !                  dx2-y2      py   -       py       s
    ijkl_ind(29,  4) = 363
    ijkl_sym(363)    = 353*(-1)
    !                  dx2-y2      py   -       py      pz
    ijkl_ind(29, 12) = 364
    ijkl_sym(364)    = 354*(-1)
    !                  dx2-y2      py   -      dz2      py
    ijkl_ind(29, 26) = 365
    ijkl_sym(365)    = 355*(-1)
    !                  dx2-y2      py   -      dzy       s
    ijkl_ind(29,  7) = 366
    ijkl_sym(366)    = 356*(-1)
    !                  dx2-y2      py   -      dzy      pz
    ijkl_ind(29, 15) = 367
    ijkl_sym(367)    = 357*(-1)
    !                  dx2-y2      py   -      dzy     dz2
    ijkl_ind(29, 33) = 368
    ijkl_sym(368)    = 358*(-1)
    !                  dx2-y2      py   -   dx2-y2      py
    ijkl_ind(29, 29) = 369
    ijkl_sym(369)    = 359
    !                  dx2-y2      py   -   dx2-y2     dzy
    ijkl_ind(29, 41) = 370
    ijkl_sym(370)    = 360
    !                  dx2-y2      py   -      dxy      px
    ijkl_ind(29, 24) = 371
    ijkl_sym(371)    = 361*(-1)
    !                  dx2-y2      py   -      dxy     dzx
    ijkl_ind(29, 39) = 372
    ijkl_sym(372)    = 362*(-1)
    ! ################ dx2-y2     dz2   #######################
    !                  dx2-y2     dz2   -       px      px
    ijkl_ind(34, 18) = 373
    !                  dx2-y2     dz2   -       py      py
    ijkl_ind(34, 25) = 374
    ijkl_sym(374)    = 373*(-1)
    !                  dx2-y2     dz2   -      dzx      px
    ijkl_ind(34, 21) = 375
    !                  dx2-y2     dz2   -      dzx     dzx
    ijkl_ind(34, 36) = 376
    !                  dx2-y2     dz2   -      dzy      py
    ijkl_ind(34, 28) = 377
    ijkl_sym(377)    = 375*(-1)
    !                  dx2-y2     dz2   -      dzy     dzy
    ijkl_ind(34, 40) = 378
    ijkl_sym(378)    = 376*(-1)
    !                  dx2-y2     dz2   -   dx2-y2       s
    ijkl_ind(34,  8) = 379
    !                  dx2-y2     dz2   -   dx2-y2      pz
    ijkl_ind(34, 16) = 380
    ijkl_sym(380)    = 375
    !                  dx2-y2     dz2   -   dx2-y2     dz2
    ijkl_ind(34, 34) = 381
    ! ################ dx2-y2     dzx   #######################
    !                  dx2-y2     dzx   -       px       s
    ijkl_ind(38,  3) = 382
    !                  dx2-y2     dzx   -       px      pz
    ijkl_ind(38, 11) = 383
    !                  dx2-y2     dzx   -      dz2      px
    ijkl_ind(38, 20) = 384
    !                  dx2-y2     dzx   -      dzx       s
    ijkl_ind(38,  6) = 385
    !                  dx2-y2     dzx   -      dzx      pz
    ijkl_ind(38, 14) = 386
    !                  dx2-y2     dzx   -      dzx     dz2
    ijkl_ind(38, 32) = 387
    !                  dx2-y2     dzx   -   dx2-y2      px
    ijkl_ind(38, 23) = 388
    !                  dx2-y2     dzx   -   dx2-y2  dx2-y2
    ijkl_ind(38, 38) = 389
    !                  dx2-y2     dzx   -      dxy      py
    ijkl_ind(38, 30) = 390
    !                  dx2-y2     dzx   -      dxy     dzy
    ijkl_ind(38, 42) = 391
    ! ################ dx2-y2     dzy   #######################
    !                  dx2-y2     dzy   -       py       s
    ijkl_ind(41,  4) = 392
    ijkl_sym(392)    = 382*(-1)
    !                  dx2-y2     dzy   -       py      pz
    ijkl_ind(41, 12) = 393
    ijkl_sym(393)    = 383*(-1)
    !                  dx2-y2     dzy   -      dz2      py
    ijkl_ind(41, 26) = 394
    ijkl_sym(394)    = 384*(-1)
    !                  dx2-y2     dzy   -      dzy       s
    ijkl_ind(41,  7) = 395
    ijkl_sym(395)    = 385*(-1)
    !                  dx2-y2     dzy   -      dzy      pz
    ijkl_ind(41, 15) = 396
    ijkl_sym(396)    = 386*(-1)
    !                  dx2-y2     dzy   -      dzy     dz2
    ijkl_ind(41, 33) = 397
    ijkl_sym(397)    = 387*(-1)
    !                  dx2-y2     dzy   -   dx2-y2      py
    ijkl_ind(41, 29) = 398
    ijkl_sym(398)    = 388
    !                  dx2-y2     dzy   -   dx2-y2     dzy
    ijkl_ind(41, 41) = 399
    ijkl_sym(399)    = 389
    !                  dx2-y2     dzy   -      dxy      px
    ijkl_ind(41, 24) = 400
    ijkl_sym(400)    = 390*(-1)
    !                  dx2-y2     dzy   -      dxy     dzx
    ijkl_ind(41, 39) = 401
    ijkl_sym(401)    = 391*(-1)
    ! ################ dx2-y2  dx2-y2   #######################
    !                  dx2-y2  dx2-y2   -        s       s
    ijkl_ind(43,  1) = 402
    !                  dx2-y2  dx2-y2   -       pz       s
    ijkl_ind(43,  2) = 403
    !                  dx2-y2  dx2-y2   -       pz      pz
    ijkl_ind(43, 10) = 404
    !                  dx2-y2  dx2-y2   -       px      px
    ijkl_ind(43, 18) = 405
    !                  dx2-y2  dx2-y2   -       py      py
    ijkl_ind(43, 25) = 406
    ijkl_sym(406)    = 405
    !                  dx2-y2  dx2-y2   -      dz2       s
    ijkl_ind(43,  5) = 407
    !                  dx2-y2  dx2-y2   -      dz2      pz
    ijkl_ind(43, 13) = 408
    !                  dx2-y2  dx2-y2   -      dz2     dz2
    ijkl_ind(43, 31) = 409
    !                  dx2-y2  dx2-y2   -      dzx      px
    ijkl_ind(43, 21) = 410
    !                  dx2-y2  dx2-y2   -      dzx     dzx
    ijkl_ind(43, 36) = 411
    !                  dx2-y2  dx2-y2   -      dzy      py
    ijkl_ind(43, 28) = 412
    ijkl_sym(412)    = 410
    !                  dx2-y2  dx2-y2   -      dzy     dzy
    ijkl_ind(43, 40) = 413
    ijkl_sym(413)    = 411
    !                  dx2-y2  dx2-y2   -   dx2-y2  dx2-y2
    ijkl_ind(43, 43) = 414
    !                  dx2-y2  dx2-y2   -      dxy     dxy
    ijkl_ind(43, 45) = 415
    ! ################    dxy       s   #######################
    !                     dxy       s   -       px      py
    ijkl_ind( 9, 19) = 416
    ijkl_sym(416)    = 335
    !                     dxy       s   -      dzx      py
    ijkl_ind( 9, 27) = 417
    ijkl_sym(417)    = 337
    !                     dxy       s   -      dzy      px
    ijkl_ind( 9, 22) = 418
    ijkl_sym(418)    = 337
    !                     dxy       s   -      dzy     dzx
    ijkl_ind( 9, 37) = 419
    ijkl_sym(419)    = 338
    !                     dxy       s   -      dxy       s
    ijkl_ind( 9,  9) = 420
    ijkl_sym(420)    = 341
    !                     dxy       s   -      dxy      pz
    ijkl_ind( 9, 17) = 421
    ijkl_sym(421)    = 337
    !                     dxy       s   -      dxy     dz2
    ijkl_ind( 9, 35) = 422
    ijkl_sym(422)    = 343
    ! ################    dxy      pz   #######################
    !                     dxy      pz   -       px      py
    ijkl_ind(17, 19) = 423
    ijkl_sym(423)    = 223
    !                     dxy      pz   -      dzx      py
    ijkl_ind(17, 27) = 424
    ijkl_sym(424)    = 219
    !                     dxy      pz   -      dzy      px
    ijkl_ind(17, 22) = 425
    ijkl_sym(425)    = 219
    !                     dxy      pz   -      dzy     dzx
    ijkl_ind(17, 37) = 426
    ijkl_sym(426)    = 226
    !                     dxy      pz   -      dxy       s
    ijkl_ind(17,  9) = 427
    ijkl_sym(427)    = 218
    !                     dxy      pz   -      dxy      pz
    ijkl_ind(17, 17) = 428
    ijkl_sym(428)    = 219
    !                     dxy      pz   -      dxy     dz2
    ijkl_ind(17, 35) = 429
    ijkl_sym(429)    = 220
    ! ################    dxy      px   #######################
    !                     dxy      px   -       py       s
    ijkl_ind(24,  4) = 430
    ijkl_sym(430)    = 353
    !                     dxy      px   -       py      pz
    ijkl_ind(24, 12) = 431
    ijkl_sym(431)    = 354
    !                     dxy      px   -      dz2      py
    ijkl_ind(24, 26) = 432
    ijkl_sym(432)    = 355
    !                     dxy      px   -      dzy       s
    ijkl_ind(24,  7) = 433
    ijkl_sym(433)    = 356
    !                     dxy      px   -      dzy      pz
    ijkl_ind(24, 15) = 434
    ijkl_sym(434)    = 357
    !                     dxy      px   -      dzy     dz2
    ijkl_ind(24, 33) = 435
    ijkl_sym(435)    = 358
    !                     dxy      px   -   dx2-y2      py
    ijkl_ind(24, 29) = 436
    ijkl_sym(436)    = 361*(-1)
    !                     dxy      px   -   dx2-y2     dzy
    ijkl_ind(24, 41) = 437
    ijkl_sym(437)    = 362*(-1)
    !                     dxy      px   -      dxy      px
    ijkl_ind(24, 24) = 438
    ijkl_sym(438)    = 359
    !                     dxy      px   -      dxy     dzx
    ijkl_ind(24, 39) = 439
    ijkl_sym(439)    = 360
    ! ################    dxy      py   #######################
    !                     dxy      py   -       px       s
    ijkl_ind(30,  3) = 440
    ijkl_sym(440)    = 353
    !                     dxy      py   -       px      pz
    ijkl_ind(30, 11) = 441
    ijkl_sym(441)    = 354
    !                     dxy      py   -      dz2      px
    ijkl_ind(30, 20) = 442
    ijkl_sym(442)    = 355
    !                     dxy      py   -      dzx       s
    ijkl_ind(30,  6) = 443
    ijkl_sym(443)    = 356
    !                     dxy      py   -      dzx      pz
    ijkl_ind(30, 14) = 444
    ijkl_sym(444)    = 357
    !                     dxy      py   -      dzx     dz2
    ijkl_ind(30, 32) = 445
    ijkl_sym(445)    = 358
    !                     dxy      py   -   dx2-y2      px
    ijkl_ind(30, 23) = 446
    ijkl_sym(446)    = 361
    !                     dxy      py   -   dx2-y2  dx2-y2
    ijkl_ind(30, 38) = 447
    ijkl_sym(447)    = 362
    !                     dxy      py   -      dxy      py
    ijkl_ind(30, 30) = 448
    ijkl_sym(448)    = 359
    !                     dxy      py   -      dxy     dzy
    ijkl_ind(30, 42) = 449
    ijkl_sym(449)    = 360
    ! ################    dxy     dz2   #######################
    !                     dxy     dz2   -       px      py
    ijkl_ind(35, 19) = 450
    ijkl_sym(450)    = 373
    !                     dxy     dz2   -      dzx      py
    ijkl_ind(35, 27) = 451
    ijkl_sym(451)    = 375
    !                     dxy     dz2   -      dzy      px
    ijkl_ind(35, 22) = 452
    ijkl_sym(452)    = 375
    !                     dxy     dz2   -      dzy     dzx
    ijkl_ind(35, 37) = 453
    ijkl_sym(453)    = 376
    !                     dxy     dz2   -      dxy       s
    ijkl_ind(35,  9) = 454
    ijkl_sym(454)    = 379
    !                     dxy     dz2   -      dxy      pz
    ijkl_ind(35, 17) = 455
    ijkl_sym(455)    = 375
    !                     dxy     dz2   -      dxy     dz2
    ijkl_ind(35, 35) = 456
    ijkl_sym(456)    = 381
    ! ################    dxy     dzx   #######################
    !                     dxy     dzx   -       py       s
    ijkl_ind(39,  4) = 457
    ijkl_sym(457)    = 382
    !                     dxy     dzx   -       py      pz
    ijkl_ind(39, 12) = 458
    ijkl_sym(458)    = 383
    !                     dxy     dzx   -      dz2      py
    ijkl_ind(39, 26) = 459
    ijkl_sym(459)    = 384
    !                     dxy     dzx   -      dzy       s
    ijkl_ind(39,  7) = 460
    ijkl_sym(460)    = 385
    !                     dxy     dzx   -      dzy      pz
    ijkl_ind(39, 15) = 461
    ijkl_sym(461)    = 386
    !                     dxy     dzx   -      dzy     dz2
    ijkl_ind(39, 33) = 462
    ijkl_sym(462)    = 387
    !                     dxy     dzx   -   dx2-y2      py
    ijkl_ind(39, 29) = 463
    ijkl_sym(463)    = 390*(-1)
    !                     dxy     dzx   -   dx2-y2     dzy
    ijkl_ind(39, 41) = 464
    ijkl_sym(464)    = 391*(-1)
    !                     dxy     dzx   -      dxy      px
    ijkl_ind(39, 24) = 465
    ijkl_sym(465)    = 388
    !                     dxy     dzx   -      dxy     dzx
    ijkl_ind(39, 39) = 466
    ijkl_sym(466)    = 389
    ! ################    dxy     dzy   #######################
    !                     dxy     dzy   -       px       s
    ijkl_ind(42,  3) = 467
    ijkl_sym(467)    = 382
    !                     dxy     dzy   -       px      pz
    ijkl_ind(42, 11) = 468
    ijkl_sym(468)    = 383
    !                     dxy     dzy   -      dz2      px
    ijkl_ind(42, 20) = 469
    ijkl_sym(469)    = 384
    !                     dxy     dzy   -      dzx       s
    ijkl_ind(42,  6) = 470
    ijkl_sym(470)    = 385
    !                     dxy     dzy   -      dzx      pz
    ijkl_ind(42, 14) = 471
    ijkl_sym(471)    = 386
    !                     dxy     dzy   -      dzx     dz2
    ijkl_ind(42, 32) = 472
    ijkl_sym(472)    = 387
    !                     dxy     dzy   -   dx2-y2      px
    ijkl_ind(42, 23) = 473
    ijkl_sym(473)    = 390
    !                     dxy     dzy   -   dx2-y2  dx2-y2
    ijkl_ind(42, 38) = 474
    ijkl_sym(474)    = 391
    !                     dxy     dzy   -      dxy      py
    ijkl_ind(42, 30) = 475
    ijkl_sym(475)    = 388
    !                     dxy     dzy   -      dxy     dzy
    ijkl_ind(42, 42) = 476
    ijkl_sym(476)    = 389
    ! ################    dxy  dx2-y2   #######################
    !                     dxy  dx2-y2   -      dxy  dx2-y2
    ijkl_ind(44, 44) = 477
    ! ################    dxy     dxy   #######################
    !                     dxy     dxy   -        s       s
    ijkl_ind(45,  1) = 478
    ijkl_sym(478)    = 402
    !                     dxy     dxy   -       pz       s
    ijkl_ind(45,  2) = 479
    ijkl_sym(479)    = 403
    !                     dxy     dxy   -       pz      pz
    ijkl_ind(45, 10) = 480
    ijkl_sym(480)    = 404
    !                     dxy     dxy   -       px      px
    ijkl_ind(45, 18) = 481
    ijkl_sym(481)    = 405
    !                     dxy     dxy   -       py      py
    ijkl_ind(45, 25) = 482
    ijkl_sym(482)    = 405
    !                     dxy     dxy   -      dz2       s
    ijkl_ind(45,  5) = 483
    ijkl_sym(483)    = 407
    !                     dxy     dxy   -      dz2      pz
    ijkl_ind(45, 13) = 484
    ijkl_sym(484)    = 408
    !                     dxy     dxy   -      dz2     dz2
    ijkl_ind(45, 31) = 485
    ijkl_sym(485)    = 409
    !                     dxy     dxy   -      dzx      px
    ijkl_ind(45, 21) = 486
    ijkl_sym(486)    = 410
    !                     dxy     dxy   -      dzx     dzx
    ijkl_ind(45, 36) = 487
    ijkl_sym(487)    = 411
    !                     dxy     dxy   -      dzy      py
    ijkl_ind(45, 28) = 488
    ijkl_sym(488)    = 410
    !                     dxy     dxy   -      dzy     dzy
    ijkl_ind(45, 40) = 489
    ijkl_sym(489)    = 411
    !                     dxy     dxy   -   dx2-y2  dx2-y2
    ijkl_ind(45, 43) = 490
    ijkl_sym(490)    = 415
    !                     dxy     dxy   -      dxy     dxy
    ijkl_ind(45, 45) = 491
    ijkl_sym(491)    = 414
  END SUBROUTINE setup_ijkl_array

END MODULE semi_empirical_int_arrays
