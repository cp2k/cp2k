!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/occupations [1.0] *
!!
!!   NAME
!!     occupations
!!
!!   FUNCTION
!!     operations with the occupation numbers
!!
!!   AUTHOR
!!     APS (12-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE occupations

  USE dft_types,                       ONLY: dft_control_type
  USE erf_fn,                          ONLY: erfc
  USE kinds,                           ONLY: dp
  USE kpoints,                         ONLY: kset
  USE mathconstants,                   ONLY: pi
  USE termination,                     ONLY: stop_program
  USE wf_types,                        ONLY: wf_kset_type,&
                                             wf_spinset_type,&
                                             wf_type

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: occupations_assign, occupations_fsum_spin
  PUBLIC :: occupations_set, occupations_fsm_set
  
  REAL (KIND=dp), PARAMETER :: occlimit_eps = 1.0E-10
  INTEGER, PARAMETER :: BR_FERMI = 511, BR_GAUSSIAN = 512
  INTEGER, PARAMETER :: BR_METHFESSELPAXTON = 513
  
CONTAINS

!******************************************************************************

SUBROUTINE occupations_assign ( efermi, broadening_method, mp_order, &
     eigenvalue, occupation, occupation_sum, entropy, ekt, cold_a )
  
    REAL(KIND=dp), INTENT(IN)                    :: efermi
    CHARACTER(LEN=*), INTENT(IN)             :: broadening_method
    INTEGER, INTENT(IN)                      :: mp_order
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: eigenvalue
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: occupation
    REAL(KIND=dp), INTENT(OUT), OPTIONAL         :: occupation_sum, entropy
    REAL(KIND=dp), INTENT(IN)                    :: ekt, cold_a

    INTEGER                                  :: order, st
    REAL(KIND=dp)                                :: AN, entr, entr_mp, HN, HNm1, &
                                                occ, pp, qq, S0, SN, x

!------------------------------------------------------------------------------

  IF ( PRESENT ( entropy ) ) entropy = 0.0_dp
  IF ( PRESENT ( occupation ) ) occupation ( : ) = 0.0_dp
  IF ( PRESENT ( occupation_sum ) ) occupation_sum = 0.0_dp
  
  DO st = 1, SIZE ( eigenvalue )
     x = ( eigenvalue ( st ) - efermi ) / ekt
     IF ( x <= -30.0_dp ) THEN
        occ = 1.0_dp
     ELSE IF ( x > 30.0_dp ) THEN
        occ = 0.0_dp
     ELSE
        IF ( broadening_method == "FERMI" ) THEN
           occ = 1.0_dp / ( 1.0_dp + EXP ( x ) )
           
           pp = occ
           qq = 1.0_dp - pp
           IF ( qq < occlimit_eps ) qq = 1.0_dp
           IF ( pp < occlimit_eps ) pp = 1.0_dp
           entr = - ( pp * LOG ( pp ) + qq * LOG ( qq ) )
           
        ELSE IF ( broadening_method == "GAUSSIAN" ) THEN
           occ = erfc ( x ) * 0.5_dp
           
           entr = EXP ( -x ** 2 ) / ( 2.0_dp * SQRT ( pi ) )
           
        ELSE IF ( broadening_method == "METHFESSELPAXTON" ) THEN
           S0 = 0.5_dp * erfc ( x )
           
           entr_mp = EXP ( -x ** 2 ) * 2.0_dp / ( 2.0_dp * SQRT ( pi ) )
           
           SN = 0.0_dp
           An = -1.0_dp / ( 4.0_dp * SQRT ( pi ) ) ! degree one
           Hn   = 2.0_dp * x
           Hnm1 = 1.0_dp
           DO order = 1, mp_order
              entr_mp = An * Hn
              SN = SN + An * Hn
              
              Hnm1 = 2.0_dp * x * Hn   &
                   - 2.0_dp * Hnm1 &
                   * ( 2.0_dp * REAL ( order,KIND=dp) - 1.0_dp )
              Hn   = 2.0_dp * x * Hnm1 &
                   - 2.0_dp * ( 2.0_dp * REAL ( order,KIND=dp) ) * Hn
              An = -An / ( 4.0_dp * REAL ( order + 1,KIND=dp) )
           END DO
           
           entr = 0.5_dp * entr_mp * EXP ( - x ** 2 )
           SN = SN * EXP ( -x ** 2 )
           occ = S0 + SN
           
        ELSE IF ( broadening_method == "COLD" ) THEN
           occ = 1.0_dp / SQRT ( pi ) * &
                ( cold_a * x ** 3 - x ** 2 - 1.5_dp * cold_a * x + 1.5_dp ) &
                * EXP ( - x ** 2 )
           
        ELSE
           CALL stop_program ( "fermisum", "no suitable broadening method" )
        END IF
     END IF
     
     IF ( PRESENT ( entropy ) ) entropy = entropy + entr
     IF ( PRESENT ( occupation ) ) occupation ( st ) = occ
     IF ( PRESENT ( occupation_sum ) ) occupation_sum = occupation_sum + occ
     
  END DO
  
END SUBROUTINE occupations_assign

!******************************************************************************

FUNCTION occupations_fsum_spin ( psi_s, efermi ) RESULT ( nelectrons_integrated )
  
  
    TYPE(wf_spinset_type), INTENT(IN), &
      TARGET                                 :: psi_s
    REAL(KIND=dp), INTENT(IN)                    :: efermi
    REAL(KIND=dp)                                :: nelectrons_integrated

    INTEGER                                  :: kpt, nkpts
    REAL(KIND=dp)                                :: cold_a = 0.0_dp , &
                                                ekt = 0.1_dp , occupation_sum
    TYPE(wf_kset_type), POINTER              :: psi_ks

!TMPTMPTMP
!TMPTMPTMP
!------------------------------------------------------------------------------

  nelectrons_integrated = 0
  
  nkpts = SIZE ( psi_s % psi_ks )
  DO kpt = 1, nkpts
     psi_ks => psi_s % psi_ks ( kpt )
     
     CALL occupations_assign ( efermi, broadening_method = "FERMI", &
          mp_order = -100, eigenvalue = psi_ks % eigenvalue, &
          occupation_sum = occupation_sum, ekt = ekt, cold_a = cold_a )
     
     nelectrons_integrated = nelectrons_integrated &
          + kset % wkpt ( kpt ) * occupation_sum
  END DO
  
END FUNCTION occupations_fsum_spin

!******************************************************************************

SUBROUTINE occupations_set ( psi, efermi, dft_control )
  
    TYPE(wf_type), DIMENSION(:), &
      INTENT(INOUT)                          :: psi
    REAL(KIND=dp), INTENT(INOUT)                 :: efermi
    TYPE(dft_control_type), INTENT(IN)       :: dft_control

    INTEGER                                  :: kpt, nkpts, nspins, nwfsets, &
                                                spin, wfset
    REAL(KIND=dp)                                :: cold_a = 0.0_dp , delta, &
                                                efermi_max, efermi_min, &
                                                ekt = 0.1_dp , nelcs_occupied
    TYPE(wf_kset_type), POINTER              :: psi_ks
    TYPE(wf_spinset_type), POINTER           :: psi_s

!TMPTMPTMP
!TMPTMPTMP
!------------------------------------------------------------------------------

  efermi_max = efermi + 0.2_dp
  
  DO
     nelcs_occupied = 0.0_dp
     nwfsets = SIZE ( psi )
     DO wfset = 1, nwfsets
        nspins = SIZE ( psi ( wfset ) % psi_s )
        DO spin = 1, nspins
           psi_s => psi ( wfset ) % psi_s ( spin )
           nelcs_occupied = nelcs_occupied &
                + occupations_fsum_spin ( psi_s, efermi_max )
        END DO
     END DO
     
     IF ( nelcs_occupied < dft_control % nelectrons &
          + dft_control % fixedspinmoment_amount ) THEN
        efermi_max = efermi_max + 2.0_dp * ( efermi_max - efermi )
     ELSE
        EXIT
     END IF
  END DO
  
  efermi_min = efermi - 0.2_dp
  
  DO
     nelcs_occupied = 0.0_dp
     nwfsets = SIZE ( psi )
     DO wfset = 1, nwfsets
        nspins = SIZE ( psi ( wfset ) % psi_s )
        DO spin = 1, nspins
           psi_s => psi ( wfset ) % psi_s ( spin )
           nelcs_occupied = nelcs_occupied &
                + occupations_fsum_spin ( psi_s, efermi_min )
        END DO
     END DO
     
     IF ( nelcs_occupied > dft_control % nelectrons &
          + dft_control % fixedspinmoment_amount ) THEN
        efermi_min = efermi_min + 2.0_dp * ( efermi_min - efermi )
     ELSE
        EXIT
     END IF
  END DO
  
  main_loop: DO
     nelcs_occupied = 0.0_dp
     nwfsets = SIZE ( psi )
     DO wfset = 1, nwfsets
        nspins = SIZE ( psi ( wfset ) % psi_s )
        DO spin = 1, nspins
           psi_s => psi ( wfset ) % psi_s ( spin )
           nelcs_occupied = nelcs_occupied &
                + occupations_fsum_spin ( psi_s, efermi )
        END DO
     END DO
     
     delta = nelcs_occupied - dft_control % nelectrons
     
     IF ( ABS ( delta ) < 1.0E-8_dp ) THEN
        EXIT
     ELSE
        IF ( delta > 0.0_dp ) THEN
           efermi_max = efermi
        ELSE
           efermi_min = efermi
        END IF
        efermi = 0.5_dp * ( efermi_max + efermi_min )
     END IF
     
  END DO main_loop
  
  nwfsets = SIZE ( psi )
  DO wfset = 1, nwfsets
     nspins = SIZE ( psi ( wfset ) % psi_s )
     DO spin = 1, nspins
        psi_s => psi ( wfset ) % psi_s ( spin )
        nkpts = SIZE ( psi_s % psi_ks )
        DO kpt = 1, nkpts
           psi_ks => psi_s % psi_ks ( kpt )
           
           CALL occupations_assign ( efermi, broadening_method = "FERMI", &
                mp_order = -100, eigenvalue = psi_ks % eigenvalue, &
                occupation = psi_ks % occupation, ekt = ekt, cold_a = cold_a )
        END DO
     END DO
  END DO
  
END SUBROUTINE occupations_set

!******************************************************************************

SUBROUTINE occupations_fsm_set ( psi, efermi, dft_control )
  
    TYPE(wf_type), DIMENSION(:), &
      INTENT(INOUT)                          :: psi
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)   :: efermi
    TYPE(dft_control_type), INTENT(IN)       :: dft_control

    INTEGER                                  :: kpt, nkpts, nspins, nwfsets, &
                                                spin, wfset
    REAL(KIND=dp) :: cold_a = 0.0_dp , delta, efermi_max, efermi_min, &
      ekt = 0.1_dp , nelcs_occupied, nelcs_wanted
    TYPE(wf_kset_type), POINTER              :: psi_ks
    TYPE(wf_spinset_type), POINTER           :: psi_s

!TMPTMPTMP
!TMPTMPTMP
!------------------------------------------------------------------------------

  IF ( SIZE ( psi ( 1 ) % psi_s ) /= 2 ) THEN
     CALL stop_program ( "occupations_fsm_set", &
          "fixed spin moment works only with nspins == 2" )
  END IF
  
  nwfsets = SIZE ( psi )
  wfset_loop: DO wfset = 1, nwfsets
     nspins = SIZE ( psi ( wfset ) % psi_s )
     spin_loop: DO spin = 1, nspins
        psi_s => psi ( wfset ) % psi_s ( spin )
        
        IF ( spin == 1 ) THEN
           nelcs_wanted = ( dft_control % nelectrons &
                + dft_control % fixedspinmoment_amount ) / 2.0_dp
        ELSE
           nelcs_wanted = ( dft_control % nelectrons &
                - dft_control % fixedspinmoment_amount ) / 2.0_dp
        END IF
        
        efermi_max = efermi ( spin ) + 0.2_dp
        
        DO
           nelcs_occupied = occupations_fsum_spin ( psi_s, efermi_max )
           IF ( nelcs_occupied < nelcs_wanted ) THEN
              efermi_max = efermi_max &
                   + 2.0_dp * ( efermi_max - efermi ( spin ) )
           ELSE
              EXIT
           END IF
        END DO
        
        efermi_min = efermi ( spin ) - 0.2_dp
        
        DO
           nelcs_occupied = occupations_fsum_spin ( psi_s, efermi_min )
           
           IF ( nelcs_occupied > nelcs_wanted ) THEN
              efermi_min = efermi_min &
                   + 1.5_dp * ( efermi_min - efermi ( spin ) )
           ELSE
              EXIT
           END IF
        END DO
        
        main_loop: DO
           nelcs_occupied = occupations_fsum_spin ( psi_s, efermi ( spin ) )
           
           delta = nelcs_occupied - nelcs_wanted
           
           IF ( ABS ( delta ) < 1.0E-8_dp ) THEN
              EXIT
           ELSE
              IF ( delta > 0.0_dp ) THEN
                 efermi_max = efermi ( spin )
              ELSE
                 efermi_min = efermi ( spin )
              END IF
              efermi = 0.5_dp * ( efermi_max + efermi_min )
           END IF
           
        END DO main_loop
        
        nkpts = SIZE ( psi_s % psi_ks )
        DO kpt = 1, nkpts
           psi_ks => psi_s % psi_ks ( kpt )
           
           CALL occupations_assign ( efermi ( spin ), &
                broadening_method = "FERMI", &
                mp_order = -100, eigenvalue = psi_ks % eigenvalue, &
                occupation = psi_ks % occupation, ekt = ekt, cold_a = cold_a )
        END DO
        
     END DO spin_loop
  END DO wfset_loop
  
END SUBROUTINE occupations_fsm_set

!******************************************************************************
!******************************************************************************
!******************************************************************************

!------------------------------------------------------------------------------!
      SUBROUTINE adjust_mu(bs)
!
! adjust the chemical potential with the bisection method
! eigenvalues do not have to be sorted
! routine adapted from CPMD
!
    TYPE(band_structure_type), INTENT(INOUT) :: bs

    INTEGER, PARAMETER                       :: it_max = 200

    INTEGER                                  :: it
    REAL(KIND=dp)                                :: amu, amu1, amu2, damu, &
                                                precision, rhint, rhint1, &
                                                rhint2

!

        precision = precpar*log(float(bs%fd%nel+1))
!..min and max eigenvalue
        amu1 = minval(bs%ev)
        amu2 = maxval(bs%ev)
!
        DO
          rhint1 = rhoint(bs,amu1)
          IF (rhint1<0.0_dp) EXIT
          amu1 = -2.0_dp*abs(amu1)
        END DO
        rhint2 = rhoint(bs,amu2)
        IF (rhint2<0.0_dp) CALL stop_program("ADJUST_MU", &
          "failed to find chem. potential: number of states is too small")
        damu = amu2 - amu1
        amu = amu1
!
        DO it = 1, it_max
          amu = amu1 + 0.5_dp*damu
          rhint = rhoint(bs,amu)
          IF (damu<delta) THEN
            IF (abs(rhint)<precision) THEN
              EXIT
            ELSE IF (abs(rhint1)<precision) THEN
              amu = amu1
              EXIT
            ELSE IF (abs(rhint1)<precision) THEN
              amu = amu2
              EXIT
            END IF
          END IF
          IF (rhint<0.0_dp) THEN
            amu1 = amu
            rhint1 = rhint
          ELSE
            amu2 = amu
            rhint2 = rhint
          END IF
          damu = amu2 - amu1
        END DO
        IF (damu>delta) THEN
          CALL stop_program("ADJUST_MU","failed to find chem. potential")
        END IF
        bs%fd%mu = amu
!
      END SUBROUTINE adjust_mu
!------------------------------------------------------------------------------!
      FUNCTION rhoint(bs,mu)
    TYPE(band_structure_type), INTENT(IN)    :: bs
    REAL(KIND=dp), INTENT(IN)                    :: mu
    REAL(KIND=dp)                                :: rhoint

    INTEGER                                  :: i, ik, is, n
    REAL(KIND=dp)                                :: arg, argmax, argmin, dramu, &
                                                drwe, remain, sdeg, wbigwe, &
                                                xlm

        argmin = log(toll)
        argmax = -argmin
        IF (bs%fd%spin_polarization==0) THEN
          sdeg = 2.0_dp
        ELSE
          sdeg = 1.0_dp
        END IF
        dramu = dround(mu,delta)
        rhoint = 0.0_dp
        wbigwe = 0.0_dp
        DO is = 1, size(bs%ev,3)
          DO ik = 1, size(bs%ev,2)
            IF (is==1) n = bs%fd%na
            IF (is==2) n = bs%fd%nb
            DO i = 1, n
              arg = -bs%fd%betael*(bs%ev(i,ik,is)-mu)
              IF (arg>argmax) THEN
                xlm = 1.0_dp
              ELSE IF (arg<argmin) THEN
                xlm = 0.0_dp
              ELSE
                xlm = exp(arg)/(exp(arg)+1.0_dp)
              END IF
              drwe = dround(bs%ev(i,ik,is),delta)
              IF (abs(drwe-dramu)<delta) THEN
                wbigwe = wbigwe + bs%kpt%weight(ik)
              ELSE
                rhoint = rhoint + bs%kpt%weight(ik)*xlm
              END IF
            END DO
          END DO
        END DO
        remain = float(bs%fd%nel) - rhoint
        IF (remain<0.0_dp) remain = sdeg
        IF (wbigwe>0.0_dp) THEN
          remain = remain/wbigwe
          remain = min(remain,sdeg)
          DO is = 1, size(bs%ev,3)
            DO ik = 1, size(bs%ev,2)
              IF (is==1) n = bs%fd%na
              IF (is==2) n = bs%fd%nb
              DO i = 1, n
                drwe = dround(bs%ev(i,ik,is),delta)
                IF (abs(drwe-dramu)<delta) rhoint = rhoint + &
                  bs%kpt%weight(ik)*remain
              END DO
            END DO
          END DO
        END IF
        rhoint = sdeg*rhoint - float(bs%fd%nel)
      END FUNCTION rhoint
!------------------------------------------------------------------------------!
      SUBROUTINE occupation(bs)
    TYPE(band_structure_type), INTENT(INOUT) :: bs

    INTEGER                                  :: i, ik, is, n
    REAL(KIND=dp)                                :: arg, argmax, argmin, dramu, &
                                                drwe, remain, rhotot, sdeg, &
                                                wbigwe, xlm

        argmin = log(toll)
        argmax = -argmin
        IF (bs%fd%spin_polarization==0) THEN
          sdeg = 2.0_dp
        ELSE
          sdeg = 1.0_dp
        END IF
        dramu = dround(bs%fd%mu,delta)
        wbigwe = 0.0_dp
        rhotot = 0.0_dp
        DO is = 1, size(bs%oc,3)
          DO ik = 1, size(bs%oc,2)
            IF (is==1) n = bs%fd%na
            IF (is==2) n = bs%fd%nb
            DO i = 1, n
              arg = -bs%fd%betael*(bs%ev(i,ik,is)-bs%fd%mu)
              IF (arg>argmax) THEN
                xlm = 1.0_dp
              ELSE IF (arg<argmin) THEN
                xlm = 0.0_dp
              ELSE
                xlm = exp(arg)/(exp(arg)+1.0_dp)
              END IF
              drwe = dround(bs%ev(i,ik,is),delta)
              IF (abs(drwe-dramu)<delta) THEN
                wbigwe = wbigwe + bs%kpt%weight(ik)
              ELSE
                rhotot = rhotot + bs%kpt%weight(ik)*xlm
                bs%oc(i,ik,is) = xlm*sdeg
              END IF
            END DO
          END DO
        END DO
        remain = float(bs%fd%nel) - rhotot
        IF (remain<0.0_dp) remain = sdeg
        IF (wbigwe>0.0_dp) THEN
          remain = remain/wbigwe
          remain = min(remain,sdeg)
          DO is = 1, size(bs%ev,3)
            DO ik = 1, size(bs%ev,2)
              IF (is==1) n = bs%fd%na
              IF (is==2) n = bs%fd%nb
              DO i = 1, n
                drwe = dround(bs%ev(i,ik,is),delta)
                IF (abs(drwe-dramu)<delta) bs%oc(i,ik,is) = remain
              END DO
            END DO
          END DO
        END IF
      END SUBROUTINE occupation
!------------------------------------------------------------------------------!
      FUNCTION dround(a,delta)
    REAL(KIND=dp), INTENT(IN)                    :: a, delta
    REAL(KIND=dp)                                :: dround

        dround = a - mod(a,delta)
      END FUNCTION dround
!------------------------------------------------------------------------------!
!******************************************************************************
!******************************************************************************
!******************************************************************************
!******************************************************************************
!******************************************************************************

END MODULE occupations

!******************************************************************************
