!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Efficient tersoff implementation
!> \author CJM, I-Feng W. Kuo, Teodoro Laino
! *****************************************************************************
MODULE manybody_tersoff

  USE cell_types,                      ONLY: cell_type
  USE f77_blas
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: pos_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE mathlib,                         ONLY: matvec_3x3
  USE pair_potential_types,            ONLY: pair_potential_pp_type,&
                                             pair_potential_single_type,&
                                             tersoff_pot_type,&
                                             tersoff_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: setup_tersoff_arrays, destroy_tersoff_arrays,&
            tersoff_forces, tersoff_energy
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_tersoff'

CONTAINS

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE tersoff_energy(pot_loc, tersoff, r_last_update_pbc, atom_a, atom_b, nloc_size,&
       full_loc_list,loc_cell_v, cell_v, drij)

    REAL(KIND=dp), INTENT(OUT)               :: pot_loc
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    INTEGER, INTENT(IN)                      :: atom_a, atom_b, nloc_size
    INTEGER, DIMENSION(2, 1:nloc_size)       :: full_loc_list
    REAL(KIND=dp), DIMENSION(3, 1:nloc_size) :: loc_cell_v
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    REAL(KIND=dp)                            :: drij

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tersoff_energy', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: b_ij, f_A, f_C, f_R

    b_ij   = ter_b_ij(tersoff,r_last_update_pbc,atom_a,atom_b,nloc_size,&
         full_loc_list, loc_cell_v, cell_v, tersoff%rcutsq)
    f_C    = ter_f_C (tersoff,drij)
    f_A    = ter_f_A (tersoff,drij)
    f_R    = ter_f_R (tersoff,drij)
    pot_loc = f_C*(f_R + b_ij*f_A)

  END SUBROUTINE tersoff_energy

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author I-Feng W. Kuo
! *****************************************************************************
  FUNCTION ter_f_C(tersoff,r)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_f_C

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_C', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: bigD, bigR, RmD, RpD

    bigR = tersoff%bigR
    bigD = tersoff%bigD
    RmD  = tersoff%bigR-tersoff%bigD
    RpD  = tersoff%bigR+tersoff%bigD
    ter_f_C = 0.0_dp
    IF(r < RmD) ter_f_C = 1.0_dp
    IF(r > RpD) ter_f_C = 0.0_dp
    IF( (r < RpD) .AND. (r > RmD) ) THEN
      ter_f_C = 0.5_dp *(1.0_dp - SIN(0.5_dp*PI*(r-bigR)/(bigD)))
    END IF
  END FUNCTION ter_f_C

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author I-Feng W. Kuo
! *****************************************************************************
  FUNCTION ter_f_C_d(tersoff,r)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_f_C_d

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_C_d', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: bigD, bigR, RmD, RpD

    bigR = tersoff%bigR
    bigD = tersoff%bigD
    RmD  = tersoff%bigR-tersoff%bigD
    RpD  = tersoff%bigR+tersoff%bigD
    ter_f_C_d = 0.0_dp
    IF(r < RmD) ter_f_C_d = 0.0_dp
    IF(r > RpD) ter_f_C_d = 0.0_dp
    IF( (r < RpD) .AND. (r > RmD) ) THEN
      ter_f_C_d = (0.25_dp*PI/bigD) * COS(0.5_dp*PI*(r-bigR)/(bigD)) / r
    END IF

  END FUNCTION ter_f_C_d

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author I-Feng W. Kuo
! *****************************************************************************
  FUNCTION ter_f_R(tersoff,r)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_f_R

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_R', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: A, lambda1

    A       = tersoff%A
    lambda1 = tersoff%lambda1
    ter_f_R = 0.0_dp
    ter_f_R = A * EXP(-lambda1*r)

  END FUNCTION ter_f_R

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author I-Feng W. Kuo
! *****************************************************************************
  FUNCTION ter_f_R_d(tersoff,r)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_f_R_d

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_R_d', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: A, f_R, lambda1

    A = tersoff%A
    lambda1 = tersoff%lambda1
    f_R = A * EXP(-lambda1*r)
    ter_f_R_d = 0.0_dp
    ter_f_R_d = lambda1 * f_R / r

  END FUNCTION ter_f_R_d

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author I-Feng W. Kuo
! *****************************************************************************
  FUNCTION ter_f_A(tersoff,r)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_f_A

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_A', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: B, lambda2

    B = tersoff%B
    lambda2 = tersoff%lambda2
    ter_f_A = 0.0_dp
    ter_f_A = - B * EXP(-lambda2*r)

  END FUNCTION ter_f_A

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author I-Feng W. Kuo
! *****************************************************************************
  FUNCTION ter_f_A_d(tersoff,r)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_f_A_d

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_A_d', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: B, lambda2

    B = tersoff%B
    lambda2 = tersoff%lambda2
    ter_f_A_d = 0.0_dp
    ter_f_A_d = - B *lambda2 * EXP(-lambda2*r) / r

  END FUNCTION ter_f_A_d

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author I-Feng W. Kuo
! *****************************************************************************
  FUNCTION ter_a_ij(tersoff,r)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_a_ij

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_a_ij', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: alpha, n

    n = tersoff%n
    alpha = tersoff%alpha
    ter_a_ij = 0.0_dp
    !Note alpha = 0.0_dp for the parameters in the paper so using simplified term
    !ter_a_ij = (1.0_dp+(alpha*ter_n_ij(tersoff,iparticle,jparticle,r))**n)**(-0.5_dp/n)
    ter_a_ij = 1.0_dp

  END FUNCTION ter_a_ij

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author I-Feng W. Kuo, Teodoro Laino
! *****************************************************************************
  FUNCTION ter_b_ij(tersoff,r_last_update_pbc,iparticle,jparticle, n_loc_size,&
       full_loc_list, loc_cell_v, cell_v,rcutsq)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    INTEGER, INTENT(IN)                      :: iparticle, jparticle, &
                                                n_loc_size
    INTEGER, DIMENSION(2, 1:n_loc_size)      :: full_loc_list
    REAL(KIND=dp), &
      DIMENSION(3, 1:n_loc_size)             :: loc_cell_v
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    REAL(KIND=dp), INTENT(IN)                :: rcutsq
    REAL(KIND=dp)                            :: ter_b_ij

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_b_ij', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: beta, n, zeta_ij

    n = tersoff%n
    beta = tersoff%beta
    ter_b_ij = 0.0_dp
    zeta_ij = ter_zeta_ij(tersoff,r_last_update_pbc,iparticle,jparticle,&
         n_loc_size,full_loc_list, loc_cell_v, cell_v,rcutsq)
    ter_b_ij = (1.0_dp+(beta*zeta_ij)**n)**(-0.5_dp/n)

  END FUNCTION ter_b_ij

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author I-Feng W. Kuo, Teodoro Laino
! *****************************************************************************
  FUNCTION ter_b_ij_d(tersoff,r_last_update_pbc,iparticle,jparticle,n_loc_size,&
       full_loc_list, loc_cell_v, cell_v,rcutsq)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    INTEGER, INTENT(IN)                      :: iparticle, jparticle, &
                                                n_loc_size
    INTEGER, DIMENSION(2, 1:n_loc_size)      :: full_loc_list
    REAL(KIND=dp), &
      DIMENSION(3, 1:n_loc_size)             :: loc_cell_v
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    REAL(KIND=dp), INTENT(IN)                :: rcutsq
    REAL(KIND=dp)                            :: ter_b_ij_d

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_b_ij_d', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: beta, beta_n, n, zeta_ij, &
                                                zeta_ij_n, zeta_ij_nm1

    n           = tersoff%n
    beta        = tersoff%beta
    beta_n      = beta**n
    zeta_ij     = ter_zeta_ij(tersoff,r_last_update_pbc,iparticle,jparticle,n_loc_size,&
         full_loc_list, loc_cell_v, cell_v,rcutsq)
    zeta_ij_nm1 = 0.0_dp
    IF (zeta_ij>0.0_dp) zeta_ij_nm1 = zeta_ij**(n-1.0_dp)
    zeta_ij_n   = zeta_ij**(n)

    ter_b_ij_d = 0.0_dp
    ter_b_ij_d = -0.5_dp * beta_n * zeta_ij_nm1 * &
                   ((1.0_dp + beta_n * zeta_ij_n)**((-0.5_dp/n) - 1.0_dp))

  END FUNCTION ter_b_ij_d

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      Using a local list of neighbors - [tlaino] 2007
!> \author I-Feng W. Kuo, Teodoro Laino
! *****************************************************************************
  FUNCTION ter_zeta_ij(tersoff,r_last_update_pbc,iparticle,jparticle,n_loc_size,&
       full_loc_list, loc_cell_v, cell_v,rcutsq)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    INTEGER, INTENT(IN)                      :: iparticle, jparticle, &
                                                n_loc_size
    INTEGER, DIMENSION(2, 1:n_loc_size)      :: full_loc_list
    REAL(KIND=dp), &
      DIMENSION(3, 1:n_loc_size)             :: loc_cell_v
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    REAL(KIND=dp), INTENT(IN)                :: rcutsq
    REAL(KIND=dp)                            :: ter_zeta_ij

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_zeta_ij', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilist, kparticle
    REAL(KIND=dp)                            :: cell_v_2(3), costheta, drij, &
                                                drik, expterm, f_C, gterm, &
                                                lambda3, n, rab2_max, rij(3), &
                                                rik(3)

    ter_zeta_ij = 0.0_dp
    n = tersoff%n
    lambda3 = tersoff%lambda3
    rab2_max = rcutsq
    rij(:) = r_last_update_pbc(jparticle)%r(:)-r_last_update_pbc(iparticle)%r(:)+cell_v
    drij = SQRT(DOT_PRODUCT(rij,rij))
    ter_zeta_ij = 0.0_dp
    DO ilist = 1, n_loc_size
       kparticle = full_loc_list(2,ilist)
       IF (kparticle==jparticle) CYCLE
       cell_v_2 = loc_cell_v(:,ilist)
       rik(:) = r_last_update_pbc(kparticle)%r(:)-r_last_update_pbc(iparticle)%r(:)+cell_v_2
       drik = DOT_PRODUCT(rik,rik)
       IF(drik>rab2_max) CYCLE
       drik = SQRT(drik)
       costheta = DOT_PRODUCT(rij,rik) / (drij*drik)
       IF(costheta<-1.0_dp) costheta=-1.0_dp
       IF(costheta>+1.0_dp) costheta=+1.0_dp
       f_C = ter_f_C(tersoff,drik)
       gterm = ter_g(tersoff,costheta)
       expterm = EXP((lambda3*(drij-drik))**3)
       ter_zeta_ij = ter_zeta_ij + f_C*gterm*expterm
    END DO

  END FUNCTION ter_zeta_ij

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>       Using a local list of neighbors - [tlaino] 2007
!> \author I-Feng W. Kuo, Teodoro Laino
! *****************************************************************************
  SUBROUTINE ter_zeta_ij_d(tersoff,r_last_update_pbc,iparticle,jparticle,f_nonbond,pv_nonbond,prefactor,&
       n_loc_size,full_loc_list, loc_cell_v, cell_v,rcutsq,use_virial)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    INTEGER, INTENT(IN)                      :: iparticle, jparticle
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_nonbond, pv_nonbond
    REAL(KIND=dp), INTENT(IN)                :: prefactor
    INTEGER, INTENT(IN)                      :: n_loc_size
    INTEGER, DIMENSION(2, 1:n_loc_size)      :: full_loc_list
    REAL(KIND=dp), &
      DIMENSION(3, 1:n_loc_size)             :: loc_cell_v
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    REAL(KIND=dp), INTENT(IN)                :: rcutsq
    LOGICAL, INTENT(IN)                      :: use_virial

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_zeta_ij_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilist, kparticle, nparticle
    REAL(KIND=dp)                            :: costheta, drij, drik, &
                                                expterm, expterm_d, f_C, &
                                                f_C_d, gterm, gterm_d, &
                                                lambda3, n, rab2_max
    REAL(KIND=dp), DIMENSION(3)              :: cell_v_2, dcosdri, dcosdrj, &
                                                dcosdrk, dri, drj, drk, rij, &
                                                rij_hat, rik, rik_hat

    n = tersoff%n
    lambda3 = tersoff%lambda3
    rab2_max = rcutsq

    rij(:) = r_last_update_pbc(jparticle)%r(:)-r_last_update_pbc(iparticle)%r(:)+cell_v
    drij = SQRT(DOT_PRODUCT(rij,rij))
    rij_hat(:) = rij(:) / drij

    nparticle = SIZE ( r_last_update_pbc )
    DO ilist = 1, n_loc_size
       kparticle = full_loc_list(2,ilist)
       IF (kparticle==jparticle) CYCLE
       cell_v_2 = loc_cell_v(:,ilist)
       rik(:) = r_last_update_pbc(kparticle)%r(:)-r_last_update_pbc(iparticle)%r(:)+cell_v_2
       drik = DOT_PRODUCT(rik,rik)

       IF(drik>rab2_max) CYCLE
       drik = SQRT(drik)
       rik_hat(:) = rik(:) / drik
       costheta = DOT_PRODUCT(rij,rik) / (drij*drik)
       IF(costheta<-1.0_dp) costheta=-1.0_dp
       IF(costheta>+1.0_dp) costheta=+1.0_dp

       dcosdrj(:) = (1.0_dp/(drij)) * ( rik_hat(:) - costheta*rij_hat(:) )
       dcosdrk(:) = (1.0_dp/(drik)) * ( rij_hat(:) - costheta*rik_hat(:) )
       dcosdri(:) = - ( dcosdrj(:) + dcosdrk(:) )

       f_C       = ter_f_C(tersoff,drik)
       f_C_d     = ter_f_C_d(tersoff,drik)
       gterm     = ter_g(tersoff,costheta)
       gterm_d   = ter_g_d(tersoff,costheta) !still need d(costheta)/dR term
       expterm   = EXP((lambda3*(drij-drik))**3)
       expterm_d = (3.0_dp)*(lambda3**3)*((drij-drik)**2)*expterm

       dri =  f_C_d * gterm   * expterm   * (rik) &
            + f_C   * gterm_d * expterm   * (dcosdri) &
            + f_C   * gterm   * expterm_d * (-rij_hat+rik_hat)

       !No f_C_d component for Rj
       drj =  f_C   * gterm_d * expterm   * (dcosdrj) &
            + f_C   * gterm   * expterm_d * (rij_hat)

       drk =  f_C_d * gterm   * expterm   * (-rik) &
            + f_C   * gterm_d * expterm   * (dcosdrk) &
            + f_C   * gterm   * expterm_d * (-rik_hat)

       f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + prefactor*dri(1)
       f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + prefactor*dri(2)
       f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + prefactor*dri(3)

       f_nonbond(1,jparticle) = f_nonbond(1,jparticle) + prefactor*drj(1)
       f_nonbond(2,jparticle) = f_nonbond(2,jparticle) + prefactor*drj(2)
       f_nonbond(3,jparticle) = f_nonbond(3,jparticle) + prefactor*drj(3)

       f_nonbond(1,kparticle) = f_nonbond(1,kparticle) + prefactor*drk(1)
       f_nonbond(2,kparticle) = f_nonbond(2,kparticle) + prefactor*drk(2)
       f_nonbond(3,kparticle) = f_nonbond(3,kparticle) + prefactor*drk(3)

       IF (use_virial) THEN
          pv_nonbond(1,1) = pv_nonbond(1,1) + prefactor*(rij(1)*drj(1) + rik(1)*drk(1))
          pv_nonbond(1,2) = pv_nonbond(1,2) + prefactor*(rij(1)*drj(2) + rik(1)*drk(2))
          pv_nonbond(1,3) = pv_nonbond(1,3) + prefactor*(rij(1)*drj(3) + rik(1)*drk(3))

          pv_nonbond(2,1) = pv_nonbond(2,1) + prefactor*(rij(2)*drj(1) + rik(2)*drk(1))
          pv_nonbond(2,2) = pv_nonbond(2,2) + prefactor*(rij(2)*drj(2) + rik(2)*drk(2))
          pv_nonbond(2,3) = pv_nonbond(2,3) + prefactor*(rij(2)*drj(3) + rik(2)*drk(3))

          pv_nonbond(3,1) = pv_nonbond(3,1) + prefactor*(rij(3)*drj(1) + rik(3)*drk(1))
          pv_nonbond(3,2) = pv_nonbond(3,2) + prefactor*(rij(3)*drj(2) + rik(3)*drk(2))
          pv_nonbond(3,3) = pv_nonbond(3,3) + prefactor*(rij(3)*drj(3) + rik(3)*drk(3))
       END IF
    END DO
  END SUBROUTINE ter_zeta_ij_d

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author I-Feng W. Kuo
! *****************************************************************************
  FUNCTION ter_g(tersoff,costheta)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    REAL(KIND=dp), INTENT(IN)                :: costheta
    REAL(KIND=dp)                            :: ter_g

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_g', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: c, c2, d, d2, h

    c = tersoff%c
    d = tersoff%d
    h = tersoff%h
    c2 = c*c
    d2 = d*d
    ter_g = 0.0_dp
    ter_g = 1.0_dp + (c2/d2) - (c2)/(d2+(h-costheta)**2)

  END FUNCTION ter_g

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author I-Feng W. Kuo
! *****************************************************************************
  FUNCTION ter_g_d(tersoff,costheta)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    REAL(KIND=dp), INTENT(IN)                :: costheta
    REAL(KIND=dp)                            :: ter_g_d

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_g_d', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: c, c2, d, d2, h, hc, sintheta

    c = tersoff%c
    d = tersoff%d
    h = tersoff%h
    c2 = c*c
    d2 = d*d
    hc = h-costheta

    sintheta = SQRT(1.0 - costheta**2)

    ter_g_d = 0.0_dp
    ! Still need d(costheta)/dR
    ter_g_d = (-2.0_dp * c2 * hc)/(d2 + hc**2)**2
  END FUNCTION ter_g_d

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>       Using a local list of neighbors - [tlaino] 2007
!> \author I-Feng W. Kuo, Teodoro Laino
! *****************************************************************************
  SUBROUTINE tersoff_forces(tersoff,r_last_update_pbc, cell_v, n_loc_size,&
       full_loc_list, loc_cell_v, iparticle,jparticle,f_nonbond,pv_nonbond,&
       use_virial, rcutsq)
    TYPE(tersoff_pot_type), POINTER          :: tersoff
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    INTEGER, INTENT(IN)                      :: n_loc_size
    INTEGER, DIMENSION(2, 1:n_loc_size)      :: full_loc_list
    REAL(KIND=dp), &
      DIMENSION(3, 1:n_loc_size)             :: loc_cell_v
    INTEGER, INTENT(IN)                      :: iparticle, jparticle
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_nonbond, pv_nonbond
    LOGICAL, INTENT(IN)                      :: use_virial
    REAL(KIND=dp), INTENT(IN)                :: rcutsq

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tersoff_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp) :: b_ij, b_ij_d, drij, f_A, f_A1, f_A2, f_A_d, f_C, f_C_d, &
      f_R, f_R1, f_R2, f_R_d, fac, prefactor, rij(3), rij_hat(3)

    CALL timeset ( routineN, handle )
    rij(:) = r_last_update_pbc(jparticle)%r(:)-r_last_update_pbc(iparticle)%r(:) + cell_v
    drij = SQRT(DOT_PRODUCT(rij,rij))
    rij_hat(:) = rij(:) / drij

    fac = - 0.5_dp
    b_ij  = ter_b_ij  (tersoff,r_last_update_pbc,iparticle,jparticle,n_loc_size,full_loc_list, loc_cell_v,cell_v,rcutsq)
    b_ij_d= ter_b_ij_d(tersoff,r_last_update_pbc,iparticle,jparticle,n_loc_size,full_loc_list, loc_cell_v,cell_v,rcutsq)
    f_A   = ter_f_A   (tersoff,drij)
    f_A_d = ter_f_A_d (tersoff,drij)
    f_C   = ter_f_C   (tersoff,drij)
    f_C_d = ter_f_C_d (tersoff,drij)
    f_R   = ter_f_R   (tersoff,drij)
    f_R_d = ter_f_R_d (tersoff,drij)

    ! Lets do the easy one first, the repulsive term
    ! Note a_ij = 1.0_dp so just going to ignore it...
    f_R1 = f_C_d * f_R * fac
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_R1*rij(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_R1*rij(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_R1*rij(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_R1*rij(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_R1*rij(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_R1*rij(3)

    IF (use_virial) THEN
       pv_nonbond(1,1) = pv_nonbond(1,1) - f_R1*rij(1)*rij(1)
       pv_nonbond(1,2) = pv_nonbond(1,2) - f_R1*rij(1)*rij(2)
       pv_nonbond(1,3) = pv_nonbond(1,3) - f_R1*rij(1)*rij(3)
       pv_nonbond(2,1) = pv_nonbond(2,1) - f_R1*rij(2)*rij(1)
       pv_nonbond(2,2) = pv_nonbond(2,2) - f_R1*rij(2)*rij(2)
       pv_nonbond(2,3) = pv_nonbond(2,3) - f_R1*rij(2)*rij(3)
       pv_nonbond(3,1) = pv_nonbond(3,1) - f_R1*rij(3)*rij(1)
       pv_nonbond(3,2) = pv_nonbond(3,2) - f_R1*rij(3)*rij(2)
       pv_nonbond(3,3) = pv_nonbond(3,3) - f_R1*rij(3)*rij(3)
    END IF

    f_R2 = f_C * f_R_d * fac
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_R2*rij(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_R2*rij(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_R2*rij(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_R2*rij(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_R2*rij(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_R2*rij(3)

    IF (use_virial) THEN
       pv_nonbond(1,1) = pv_nonbond(1,1) - f_R2*rij(1)*rij(1)
       pv_nonbond(1,2) = pv_nonbond(1,2) - f_R2*rij(1)*rij(2)
       pv_nonbond(1,3) = pv_nonbond(1,3) - f_R2*rij(1)*rij(3)
       pv_nonbond(2,1) = pv_nonbond(2,1) - f_R2*rij(2)*rij(1)
       pv_nonbond(2,2) = pv_nonbond(2,2) - f_R2*rij(2)*rij(2)
       pv_nonbond(2,3) = pv_nonbond(2,3) - f_R2*rij(2)*rij(3)
       pv_nonbond(3,1) = pv_nonbond(3,1) - f_R2*rij(3)*rij(1)
       pv_nonbond(3,2) = pv_nonbond(3,2) - f_R2*rij(3)*rij(2)
       pv_nonbond(3,3) = pv_nonbond(3,3) - f_R2*rij(3)*rij(3)
    END IF

    ! Lets do the f_A1 piece derivative of F_C
    f_A1 = f_C_d * b_ij * f_A * fac
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_A1*rij(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_A1*rij(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_A1*rij(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_A1*rij(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_A1*rij(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_A1*rij(3)

    IF (use_virial) THEN
       pv_nonbond(1,1) = pv_nonbond(1,1) - f_A1*rij(1)*rij(1)
       pv_nonbond(1,2) = pv_nonbond(1,2) - f_A1*rij(1)*rij(2)
       pv_nonbond(1,3) = pv_nonbond(1,3) - f_A1*rij(1)*rij(3)
       pv_nonbond(2,1) = pv_nonbond(2,1) - f_A1*rij(2)*rij(1)
       pv_nonbond(2,2) = pv_nonbond(2,2) - f_A1*rij(2)*rij(2)
       pv_nonbond(2,3) = pv_nonbond(2,3) - f_A1*rij(2)*rij(3)
       pv_nonbond(3,1) = pv_nonbond(3,1) - f_A1*rij(3)*rij(1)
       pv_nonbond(3,2) = pv_nonbond(3,2) - f_A1*rij(3)*rij(2)
       pv_nonbond(3,3) = pv_nonbond(3,3) - f_A1*rij(3)*rij(3)
    END IF

    ! Lets do the f_A2 piece derivative of F_A
    f_A2 = f_C * b_ij * f_A_d * fac
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_A2*rij(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_A2*rij(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_A2*rij(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_A2*rij(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_A2*rij(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_A2*rij(3)

    IF (use_virial) THEN
       pv_nonbond(1,1) = pv_nonbond(1,1) - f_A2*rij(1)*rij(1)
       pv_nonbond(1,2) = pv_nonbond(1,2) - f_A2*rij(1)*rij(2)
       pv_nonbond(1,3) = pv_nonbond(1,3) - f_A2*rij(1)*rij(3)
       pv_nonbond(2,1) = pv_nonbond(2,1) - f_A2*rij(2)*rij(1)
       pv_nonbond(2,2) = pv_nonbond(2,2) - f_A2*rij(2)*rij(2)
       pv_nonbond(2,3) = pv_nonbond(2,3) - f_A2*rij(2)*rij(3)
       pv_nonbond(3,1) = pv_nonbond(3,1) - f_A2*rij(3)*rij(1)
       pv_nonbond(3,2) = pv_nonbond(3,2) - f_A2*rij(3)*rij(2)
       pv_nonbond(3,3) = pv_nonbond(3,3) - f_A2*rij(3)*rij(3)
    END IF

    ! Lets do the f_A3 piece derivative of b_ij
    prefactor = f_C * b_ij_d * f_A * fac  ! Note need to do d(Zeta_ij)/dR
    CALL ter_zeta_ij_d(tersoff,r_last_update_pbc,iparticle,jparticle,f_nonbond,pv_nonbond,prefactor,&
         n_loc_size, full_loc_list, loc_cell_v,cell_v,rcutsq, use_virial)
    CALL timestop ( handle )
  END SUBROUTINE tersoff_forces

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      Fast implementation of the tersoff potential - [tlaino] 2007
!> \author Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE setup_tersoff_arrays(nonbonded, potparm, glob_loc_list, glob_cell_v, glob_loc_list_a, cell, error)
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    INTEGER, DIMENSION(:, :), POINTER        :: glob_loc_list
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: glob_cell_v
    INTEGER, DIMENSION(:), POINTER           :: glob_loc_list_a
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_tersoff_arrays', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iend, igrp, ikind, &
                                                ilist, ipair, istart, jkind, &
                                                nkinds, npairs, npairs_tot, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: work_list, work_list2
    INTEGER, DIMENSION(:, :), POINTER        :: list
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, cvi
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rwork_list
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot

    CPPostcondition(.NOT.ASSOCIATED(glob_loc_list),cp_failure_level,routineP,error,failure)
    CPPostcondition(.NOT.ASSOCIATED(glob_loc_list_a),cp_failure_level,routineP,error,failure)
    CPPostcondition(.NOT.ASSOCIATED(glob_cell_v),cp_failure_level,routineP,error,failure)
    failure = .FALSE.
    CALL timeset ( routineN, handle )
    npairs_tot = 0
    nkinds = SIZE(potparm%pot, 1)
    DO ilist=1,nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
       npairs=neighbor_kind_pair%npairs
       IF (npairs ==0) CYCLE
       Kind_Group_Loop1: DO igrp = 1, neighbor_kind_pair%ngrp_kind
          istart  = neighbor_kind_pair%grp_kind_start(igrp)
          iend    = neighbor_kind_pair%grp_kind_end(igrp)
          ikind   = neighbor_kind_pair%ij_kind(1,igrp)
          jkind   = neighbor_kind_pair%ij_kind(2,igrp)
          pot => potparm %pot( ikind, jkind )%pot
          npairs = iend-istart+1
          IF (pot%no_mb) CYCLE
          DO i = 1, SIZE(pot%type)
             IF (pot%type(i)==tersoff_type) npairs_tot = npairs_tot + npairs
          END DO
       END DO Kind_Group_Loop1
    END DO
    ALLOCATE(work_list(npairs_tot), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(work_list2(npairs_tot), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(glob_loc_list(2,npairs_tot), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(glob_cell_v(3,npairs_tot), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Fill arrays with data
    npairs_tot = 0
    DO ilist=1,nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
       npairs=neighbor_kind_pair%npairs
       IF (npairs ==0) CYCLE
       Kind_Group_Loop2: DO igrp = 1, neighbor_kind_pair%ngrp_kind
          istart  = neighbor_kind_pair%grp_kind_start(igrp)
          iend    = neighbor_kind_pair%grp_kind_end(igrp)
          ikind   = neighbor_kind_pair%ij_kind(1,igrp)
          jkind   = neighbor_kind_pair%ij_kind(2,igrp)
          list  => neighbor_kind_pair%list
          cvi   =  neighbor_kind_pair%cell_vector
          pot   => potparm %pot( ikind, jkind )%pot
          npairs = iend-istart+1
          IF (pot%no_mb) CYCLE
          CALL matvec_3x3(cell_v, cell%hmat, cvi)
          DO i = 1, SIZE(pot%type)
             ! TERSOFF
             IF (pot%type(i)==tersoff_type) THEN
                DO ipair = 1, npairs
                   glob_loc_list(:,npairs_tot+ipair) = list(:,istart-1+ipair)
                   glob_cell_v(1:3,npairs_tot+ipair) = cell_v(1:3)
                END DO
                npairs_tot = npairs_tot + npairs
             END IF
          END DO
       END DO Kind_Group_Loop2
    END DO
    ! Order the arrays w.r.t. the first index of glob_loc_list
    CALL sort(glob_loc_list(1,:),npairs_tot,work_list)
    DO ipair = 1, npairs_tot
       work_list2(ipair)=glob_loc_list(2,work_list(ipair))
    END DO
    glob_loc_list(2,:) = work_list2
    DEALLOCATE(work_list2, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(rwork_list(3,npairs_tot), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO ipair = 1, npairs_tot
       rwork_list(:,ipair)=glob_cell_v(:,work_list(ipair))
    END DO
    glob_cell_v = rwork_list
    DEALLOCATE(rwork_list, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(work_list, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(glob_loc_list_a(npairs_tot), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    glob_loc_list_a = glob_loc_list(1,:)
    CALL timestop ( handle )
  END SUBROUTINE setup_tersoff_arrays

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      Fast implementation of the tersoff potential - [tlaino] 2007
!> \author Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE destroy_tersoff_arrays(glob_loc_list, glob_cell_v, glob_loc_list_a, error)
    INTEGER, DIMENSION(:, :), POINTER        :: glob_loc_list
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: glob_cell_v
    INTEGER, DIMENSION(:), POINTER           :: glob_loc_list_a
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'destroy_tersoff_arrays', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(glob_loc_list)) THEN
       DEALLOCATE(glob_loc_list, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (ASSOCIATED(glob_loc_list_a)) THEN
       DEALLOCATE(glob_loc_list_a, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (ASSOCIATED(glob_cell_v)) THEN
       DEALLOCATE(glob_cell_v, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE destroy_tersoff_arrays

END MODULE manybody_tersoff

