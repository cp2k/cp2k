!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004 CP2K developers group                                  !
!-----------------------------------------------------------------------------!
!!****** cp2k/semi_empirical_parameters [1.0] *
!!
!!   NAME
!!     semi_empirical_parameters
!!
!!   FUNCTION
!!     Default parameter sets for semi empirical models
!!
!!   AUTHOR
!!     JGH (14.08.2004)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE semi_empirical_parameters

  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE periodic_table,                  ONLY: ptable, nelem
  USE physcon,                         ONLY: evolt
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "semi_empirical_parameters"

! *** Public subroutines ***

  PUBLIC :: eht_default_parameter

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE eht_default_parameter ( sep, z )

    TYPE(semi_empirical_type),POINTER    :: sep
    INTEGER                              :: z

    CHARACTER(LEN=*), PARAMETER :: routine_name = "eht_default_parameter"

    INTEGER, DIMENSION(0:3)              :: norb
    REAL(dp), DIMENSION(1,0:3)           :: alpha

! -----------------------------------------------------------------------------

    sep%typ = "EHT"
    sep%name = ptable(z)%name
    sep%z = z
    CALL valence_electrons ( z, sep%zeff, sep%natorb )
    norb=0
    IF(sep%natorb >= 1) norb(0)=1
    IF(sep%natorb >= 4) norb(1)=1
    IF(sep%natorb >= 9) norb(2)=1
    IF(sep%natorb >=16) norb(3)=1
    CALL slater_exponents ( z, norb, alpha, 1 ) 
    sep%sto_exponents(0:3) = alpha(1,0:3)
    sep%beta(0:3) = ptable(z)%eht_param(0:3)/evolt
    IF ( ABS ( SUM (sep%beta(0:3)) ) > 0._dp ) sep%defined = .TRUE.
    sep%ngauss = 6

  END SUBROUTINE eht_default_parameter

! *****************************************************************************
  
  SUBROUTINE valence_electrons ( z, zeff, natorb )
    INTEGER, INTENT ( IN )              :: z
    INTEGER, INTENT ( OUT )             :: zeff
    INTEGER, INTENT ( OUT )             :: natorb

    INTEGER, DIMENSION(0:nelem), PARAMETER  :: zval = (/ -1, &
       1,  2,  1,  2,  3,  4,  5,  6,  7,  0,  1,  2,  3,  4,  5,  6,  7,  0,&
       1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,  2,  3,  4,  5,  6,  7,  0,&
       1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,  2,  3,  4,  5,  6,  7,  0,&
       1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,  3,  4,&
       5,  6,  7,  8,  9, 10, 11,  2,  3,  4,  5,  6,  7, -1, -1, -1, -1,  4,&
      -1,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1         &
       /)
    INTEGER, DIMENSION(0:nelem), PARAMETER  :: norb = (/ -1, &
       1,  1,  4,  4,  4,  4,  4,  4,  4,  4,  0,  4,  4,  4,  4,  4,  4,  4,&
       0,  4,  9,  9,  9,  9,  9,  9,  9,  9,  9,  4,  4,  4,  4,  4,  4,  4,&
       4,  4,  9,  9,  9,  9,  9,  9,  9,  9,  9,  4,  4,  4,  4,  4,  4,  4,&
       4,  4,  9, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,  9,  9,&
       9,  9,  9,  9,  9,  9,  9,  4,  4,  4,  4,  4,  4, -1, -1, -1, -1,  9,&
      -1,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1         &
       /)

    zeff = zval(z)
    natorb = norb(z)

  END SUBROUTINE valence_electrons 

! *****************************************************************************

  SUBROUTINE slater_exponents (z,norb,sto,itype) 
    INTEGER, INTENT (IN)                 :: z
    INTEGER, DIMENSION(0:3), INTENT (IN) :: norb
    REAL(dp), DIMENSION(:,0:)            :: sto
    INTEGER, INTENT (IN)                 :: itype
!
!   Slater exponents according to
!   itype = 0 : Clementi-Raimondi if available, Slater rules for others
!   itype = 1 : Slater rules
!   itype = 2 : Clementi-Raimondi
!
!   Clementi-Raimondi Table
    REAL(dp),DIMENSION(18),PARAMETER   :: crtab1s = (/ &
     1.0000, 1.6875, 2.6906, 3.6848, 4.6795, 5.6727, 6.6651, 7.6579, 8.6501, 9.6421,&
    10.6259,11.6089,12.5910,13.5745,14.5578,15.5409,16.5239,17.5075 /)
    REAL(dp),DIMENSION(18),PARAMETER   :: crtab2s = (/ &
     0.0000, 0.0000, 0.6396, 0.9560, 1.2881, 1.6083, 1.9237, 2.2458, 2.5638, 2.8792,&
     3.2857, 3.6960, 4.1068, 5.5100, 4.9125, 5.3144, 5.7152, 6.1152 /)
    REAL(dp),DIMENSION(18),PARAMETER   :: crtab3s = (/ &
     0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,&
     0.8358, 1.1025, 1.3724, 1.6344, 1.8806, 2.1223, 2.3561, 2.5856 /)
    REAL(dp),DIMENSION(18),PARAMETER   :: crtab2p = (/ &
     0.0000, 0.0000, 0.0000, 0.0000, 1.2107, 1.5679, 1.9170, 2.2266, 2.5500, 2.8792,&
     3.4009, 3.9129, 4.4817, 4.9725, 5.4806, 5.9885, 6.4966, 7.0041 /)
    REAL(dp),DIMENSION(18),PARAMETER   :: crtab3p = (/ &
     0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,&
     0.0000, 0.0000, 1.3552, 1.4284, 1.6288, 1.8273, 2.0387, 2.2547 /)

    INTEGER, DIMENSION(1:4,1:7)   :: ne
    INTEGER, DIMENSION(1:4)       :: no
    INTEGER                       :: i,l,nj,j

    sto = 0._dp
    ne = 0
    DO l=1,4
      nj = 2*(l-1)+1
      no(l) = l-1
      DO i=l,7
        ne(l,i) = ptable(z)%e_conv(l-1) - 2*nj*(i-l)
        ne(l,i) = MAX(ne(l,i),0)
        ne(l,i) = MIN(ne(l,i),2*nj)
        IF(ne(l,i) > 0) no(l)=no(l)+1
      END DO
    END DO
    DO l=0,3
      DO i=1,norb(l)
        j=no(l+1)-norb(l)+1
        IF ( itype == 0 .OR. itype == 1 ) THEN
          sto(i,l)=srules(z,ne,j,l)
        END IF
        IF ( itype == 0 .OR. itype == 2 ) THEN
          IF ( z <= 18 .AND. z > 0 ) THEN
            IF ( l == 0 ) THEN
              IF ( j == 1 ) sto(i,l)=crtab1s(z)
              IF ( j == 2 ) sto(i,l)=crtab2s(z)
              IF ( j == 3 ) sto(i,l)=crtab3s(z)
            ELSE IF ( l == 1 ) THEN
              IF ( j == 2 ) sto(i,l)=crtab2p(z)
              IF ( j == 3 ) sto(i,l)=crtab3p(z)
            END IF
          END IF
        END IF
      END DO
    END DO
    
  END SUBROUTINE slater_exponents  

  FUNCTION srules(z,ne,n,l)
    ! Slater rules
    REAL(dp)                             :: srules
    INTEGER                              :: z
    INTEGER, DIMENSION(4,7)              :: ne
    INTEGER                              :: n,l

    INTEGER                              :: l1,l2,m,m1,m2,i
    REAL(dp)                             :: s
    REAL(dp), DIMENSION(7), PARAMETER    :: xns = (/ &
              1.0_dp,2.0_dp,3.0_dp,3.7_dp,4.0_dp,4.2_dp,4.4_dp /)

    ! Calculate the shielding
    s = 0.0_dp
    ! The complete shell
    l1=l+1
    IF(l1 == 1) l2=2
    IF(l1 == 2) l2=1
    IF(l1 == 3) l2=4
    IF(l1 == 4) l2=3
    ! Rule a) no contribution from shells further out
    ! Rule b) 0.35 (1s 0.3) from each other electron in the same shell
    IF(n == 1) THEN
      m=ne(1,1)
      s=s+0.3_dp*REAL(m-1,dp)
    ELSE
      m=ne(l1,n)+ne(l2,n)
      s=s+0.35_dp*REAL(m-1,dp)
    END IF
    ! Rule c) if (s,p) shell 0.85 from each electron with n-1, and 1.0
    !      from all electrons further in
    IF(l1+l2 == 3) THEN
      IF(n > 1) THEN
        m1=ne(1,n-1)+ne(2,n-1)+ne(3,n-1)+ne(4,n-1)
        m2=0
       DO i=1,n-2
          m2=m2+ne(1,i)+ne(2,i)+ne(3,i)+ne(4,I)
        END DO
        s=s+0.85_dp*REAL(m1,dp)+1._dp*REAL(m2,dp)
      END IF
    ELSE
      ! Rule d) if (d,f) shell 1.0 from each electron inside
      m=0
      DO i=1,n-1
        m=m+ne(1,i)+ne(2,i)+ne(3,i)+ne(4,i)
      END DO
      s=s+1._dp*REAL(m,dp)
    END IF
    ! Slater exponent is (Z-S)/NS
    srules = (REAL(z,dp) - s)/xns(n)

  END FUNCTION srules

! *****************************************************************************

END MODULE semi_empirical_parameters
