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
!!     ??? (12-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE occupations
  
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : pi
  USE termination, ONLY : stop_program
  USE wf_types, ONLY : wf_type, wf_spinset_type, wf_kset_type
  USE dft_types, ONLY : dft_control_type
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: occupations_assign, occupations_fsum_spin
  PUBLIC :: occupations_set, occupations_fsm_set
  
  REAL ( dbl ), PARAMETER :: occlimit_eps = 1.0E-10
  INTEGER, PARAMETER :: BR_FERMI = 511, BR_GAUSSIAN = 512
  INTEGER, PARAMETER :: BR_METHFESSELPAXTON = 513
  
CONTAINS

!******************************************************************************

SUBROUTINE occupations_assign ( efermi, broadening_method, mp_order, &
     eigenvalue, occupation, occupation_sum, entropy, ekt, cold_a )
  
! Arguments
  REAL ( dbl ), INTENT ( IN ) :: efermi
  REAL ( dbl ), INTENT ( IN ) :: ekt
  REAL ( dbl ), INTENT ( IN ) :: cold_a
  CHARACTER ( LEN = * ), INTENT ( IN ) :: broadening_method
  INTEGER, INTENT ( IN ) :: mp_order
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: eigenvalue
  REAL ( dbl ), DIMENSION ( : ), INTENT ( OUT ), OPTIONAL :: occupation
  REAL ( dbl ), INTENT ( OUT ), OPTIONAL :: entropy
  REAL ( dbl ), INTENT ( OUT ), OPTIONAL :: occupation_sum
  
! Locals
  INTEGER :: st, order
  REAL ( dbl ) :: x, h, AN, SN, S0, HN, HNm1
  REAL ( dbl ) :: pp, qq, entr, entr_mp, occ
  
  REAL ( dbl ), EXTERNAL :: erfc
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( entropy ) ) entropy = 0.0_dbl
  IF ( PRESENT ( occupation ) ) occupation ( : ) = 0.0_dbl
  IF ( PRESENT ( occupation_sum ) ) occupation_sum = 0.0_dbl
  
  DO st = 1, SIZE ( eigenvalue )
     x = ( eigenvalue ( st ) - efermi ) / ekt
     IF ( x <= -30.0_dbl ) THEN
        occ = 1.0_dbl
     ELSE IF ( x > 30.0_dbl ) THEN
        occ = 0.0_dbl
     ELSE
        IF ( broadening_method == "FERMI" ) THEN
           occ = 1.0_dbl / ( 1.0_dbl + EXP ( x ) )
           
           pp = occ
           qq = 1.0_dbl - pp
           IF ( qq < occlimit_eps ) qq = 1.0_dbl
           IF ( pp < occlimit_eps ) pp = 1.0_dbl
           entr = - ( pp * LOG ( pp ) + qq * LOG ( qq ) )
           
        ELSE IF ( broadening_method == "GAUSSIAN" ) THEN
           occ = erfc ( x ) * 0.5_dbl
           
           entr = EXP ( -x ** 2 ) / ( 2.0_dbl * SQRT ( pi ) )
           
        ELSE IF ( broadening_method == "METHFESSELPAXTON" ) THEN
           S0 = 0.5_dbl * erfc ( x )
           
           entr_mp = EXP ( -x ** 2 ) * 2.0_dbl / ( 2.0_dbl * SQRT ( pi ) )
           
           SN = 0.0_dbl
           An = -1.0_dbl / ( 4.0_dbl * SQRT ( pi ) ) ! degree one
           Hn   = 2.0_dbl * x
           Hnm1 = 1.0_dbl
           DO order = 1, mp_order
              entr_mp = An * Hn
              SN = SN + An * Hn
              
              Hnm1 = 2.0_dbl * x * Hn   &
                   - 2.0_dbl * Hnm1 &
                   * ( 2.0_dbl * REAL ( order, dbl ) - 1.0_dbl )
              Hn   = 2.0_dbl * x * Hnm1 &
                   - 2.0_dbl * ( 2.0_dbl * REAL ( order, dbl ) ) * Hn
              An = -An / ( 4.0_dbl * REAL ( order + 1, dbl ) )
           END DO
           
           entr = 0.5_dbl * entr_mp * EXP ( - x ** 2 )
           SN = SN * EXP ( -x ** 2 )
           occ = S0 + SN
           
        ELSE IF ( broadening_method == "COLD" ) THEN
           occ = 1.0_dbl / SQRT ( pi ) * &
                ( cold_a * x ** 3 - x ** 2 - 1.5_dbl * cold_a * x + 1.5_dbl ) &
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
  
  USE kpoints, ONLY : kset
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( wf_spinset_type ), INTENT ( IN ), TARGET :: psi_s
  REAL ( dbl ), INTENT ( IN ) :: efermi
  
  REAL ( dbl ) :: nelectrons_integrated
  
! Local
  INTEGER :: kpt, nkpts, spin, nspins
  TYPE ( wf_kset_type ), POINTER :: psi_ks
  REAL ( dbl ) :: occupation_sum
  REAL ( dbl ) :: ekt = 0.1_dbl !TMPTMPTMP
  REAL ( dbl ) :: cold_a = 0.0_dbl !TMPTMPTMP
  
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
  IMPLICIT NONE
  
! Arguments
  TYPE ( wf_type ), DIMENSION ( : ), INTENT ( INOUT ) :: psi
  TYPE ( dft_control_type ), INTENT ( IN ) :: dft_control
  REAL ( dbl ), INTENT ( INOUT ) :: efermi
  
! Locals
  INTEGER :: kpt, nkpts, spin, nspins, wfset, nwfsets
  REAL ( dbl ) :: efermi_min, efermi_max, nelcs_occupied, delta
  TYPE ( wf_kset_type ), POINTER :: psi_ks
  TYPE ( wf_spinset_type ), POINTER :: psi_s
  REAL ( dbl ) :: ekt = 0.1_dbl !TMPTMPTMP
  REAL ( dbl ) :: cold_a = 0.0_dbl !TMPTMPTMP
  
!------------------------------------------------------------------------------
  
  efermi_max = efermi + 0.2_dbl
  
  DO
     nelcs_occupied = 0.0_dbl
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
        efermi_max = efermi_max + 2.0_dbl * ( efermi_max - efermi )
     ELSE
        EXIT
     END IF
  END DO
  
  efermi_min = efermi - 0.2_dbl
  
  DO
     nelcs_occupied = 0.0_dbl
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
        efermi_min = efermi_min + 2.0_dbl * ( efermi_min - efermi )
     ELSE
        EXIT
     END IF
  END DO
  
  main_loop: DO
     nelcs_occupied = 0.0_dbl
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
     
     IF ( ABS ( delta ) < 1.0E-8_dbl ) THEN
        EXIT
     ELSE
        IF ( delta > 0.0_dbl ) THEN
           efermi_max = efermi
        ELSE
           efermi_min = efermi
        END IF
        efermi = 0.5_dbl * ( efermi_max + efermi_min )
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
  IMPLICIT NONE
  
! Arguments
  TYPE ( wf_type ), DIMENSION ( : ), INTENT ( INOUT ) :: psi
  TYPE ( dft_control_type ), INTENT ( IN ) :: dft_control
  REAL ( dbl ), DIMENSION ( : ), INTENT ( INOUT ) :: efermi
  
! Locals
  INTEGER :: kpt, nkpts, spin, nspins, wfset, nwfsets
  REAL ( dbl ) :: efermi_min, efermi_max, nelcs_occupied, delta
  REAL ( dbl ) :: nelcs_wanted
  TYPE ( wf_kset_type ), POINTER :: psi_ks
  TYPE ( wf_spinset_type ), POINTER :: psi_s
  REAL ( dbl ) :: ekt = 0.1_dbl !TMPTMPTMP
  REAL ( dbl ) :: cold_a = 0.0_dbl !TMPTMPTMP
  
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
                + dft_control % fixedspinmoment_amount ) / 2.0_dbl
        ELSE
           nelcs_wanted = ( dft_control % nelectrons &
                - dft_control % fixedspinmoment_amount ) / 2.0_dbl
        END IF
        
        efermi_max = efermi ( spin ) + 0.2_dbl
        
        DO
           nelcs_occupied = occupations_fsum_spin ( psi_s, efermi_max )
           IF ( nelcs_occupied < nelcs_wanted ) THEN
              efermi_max = efermi_max &
                   + 2.0_dbl * ( efermi_max - efermi ( spin ) )
           ELSE
              EXIT
           END IF
        END DO
        
        efermi_min = efermi ( spin ) - 0.2_dbl
        
        DO
           nelcs_occupied = occupations_fsum_spin ( psi_s, efermi_min )
           
           IF ( nelcs_occupied > nelcs_wanted ) THEN
              efermi_min = efermi_min &
                   + 1.5_dbl * ( efermi_min - efermi ( spin ) )
           ELSE
              EXIT
           END IF
        END DO
        
        main_loop: DO
           nelcs_occupied = occupations_fsum_spin ( psi_s, efermi ( spin ) )
           
           delta = nelcs_occupied - nelcs_wanted
           
           IF ( ABS ( delta ) < 1.0E-8_dbl ) THEN
              EXIT
           ELSE
              IF ( delta > 0.0_dbl ) THEN
                 efermi_max = efermi ( spin )
              ELSE
                 efermi_min = efermi ( spin )
              END IF
              efermi = 0.5_dbl * ( efermi_max + efermi_min )
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

END MODULE occupations
