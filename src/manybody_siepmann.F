!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!  \brief implementation of dipole and three-body part of Siepmann-Sprik 
!         potential
!>        dipole term: 3rd term in Eq. (1) in J. Chem. Phys., Vol. 102, p.511
!>        three-body term: Eq. (4) in J. Chem. Phys., Vol. 102, p. 511
!>        remaining terms of Siepmann-Sprik potential can be given via
!>        the GENPOT section
!> \par History
!>      12.2012 created
!> \author Dorothea Golze
! *****************************************************************************
MODULE manybody_siepmann

  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: pos_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: matvec_3x3
  USE message_passing,                 ONLY: mp_sum
  USE pair_potential_types,            ONLY: pair_potential_pp_type,&
                                             pair_potential_single_type,&
                                             siepmann_pot_type,&
                                             siepmann_type
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: setup_siepmann_arrays, destroy_siepmann_arrays,&
            siepmann_energy, siepmann_forces_v2, siepmann_forces_v3,&
            print_nr_oh_siepmann 
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_siepmann'

CONTAINS

! *****************************************************************************
!> \brief  energy of two-body dipole term and three-body term
!> \param pot_loc ...
!> \param siepmann ...
!> \param r_last_update_pbc ...
!> \param atom_a ...
!> \param atom_b ...
!> \param nloc_size ...
!> \param full_loc_list ...
!> \param loc_cell_v ...
!> \param cell_v ...
!> \param cell ...
!> \param drij ...
!> \param particle_set ...
!> \param nr_oh ...
!> \author Dorothea Golze - 11.2012 - University of Zurich
! *****************************************************************************
  SUBROUTINE siepmann_energy(pot_loc, siepmann, r_last_update_pbc, atom_a, atom_b,&
             nloc_size, full_loc_list,loc_cell_v, cell_v, cell, drij, particle_set,&
             nr_oh)

    REAL(KIND=dp), INTENT(OUT)               :: pot_loc
    TYPE(siepmann_pot_type), POINTER         :: siepmann
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    INTEGER, INTENT(IN)                      :: atom_a, atom_b, nloc_size
    INTEGER, DIMENSION(2, 1:nloc_size)       :: full_loc_list
    REAL(KIND=dp), DIMENSION(3, 1:nloc_size) :: loc_cell_v
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp)                            :: drij
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(INOUT)                   :: nr_oh

    CHARACTER(LEN=*), PARAMETER :: routineN = 'siepmann_energy', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: a_ij, D, E, f2, Phi_ij, &
                                                pot_loc_v2, pot_loc_v3

    a_ij   = siep_a_ij(siepmann,r_last_update_pbc,atom_a,atom_b,nloc_size,&
             full_loc_list, cell_v, siepmann%rcutsq,particle_set,&
             cell)
    Phi_ij  = siep_Phi_ij(siepmann,r_last_update_pbc,atom_a,atom_b,nloc_size,&
              loc_cell_v, cell_v, cell, siepmann%rcutsq,&
              particle_set, nr_oh)
    f2     = siep_f2 (siepmann,drij)
    D = siepmann%D
    E = siepmann%E

    !two-body part --> dipole term
    pot_loc_v2 = -D*f2*drij**(-3)*Phi_ij

    !three-body part
    pot_loc_v3 = E*f2*drij**(-siepmann%beta)*a_ij

    pot_loc = pot_loc_v2 + pot_loc_v3

  END SUBROUTINE siepmann_energy

! *****************************************************************************
!> \brief f2(r) corresponds to Equation (2) in J. Chem. Phys., Vol. 102, p.511
!> \param siepmann ...
!> \param r distance between oxygen and metal atom 
!> \retval siep_f2 ...
! *****************************************************************************
  FUNCTION siep_f2(siepmann,r)
    TYPE(siepmann_pot_type), POINTER         :: siepmann
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: siep_f2

    CHARACTER(LEN=*), PARAMETER :: routineN = 'siep_f2', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: rcut

    rcut=SQRT(siepmann%rcutsq)
    siep_f2 = 0.0_dp
    IF( r < rcut ) THEN
      siep_f2 = EXP(siepmann%B/(r-rcut))
    END IF
  END FUNCTION siep_f2

! *****************************************************************************
!> \brief f2(r)_d derivative of f2
!> \param siepmann ...
!> \param r distance between oxygen and metal atom 
!> \retval siep_f2_d ...
! *****************************************************************************
  FUNCTION siep_f2_d(siepmann,r)
    TYPE(siepmann_pot_type), POINTER         :: siepmann
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: siep_f2_d

    CHARACTER(LEN=*), PARAMETER :: routineN = 'siep_f2_d', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: B, rcut

    rcut=SQRT(siepmann%rcutsq)
    B=siepmann%B
    siep_f2_d = 0.0_dp
    IF( r < rcut ) THEN
      siep_f2_d = -B*EXP(B/(r-rcut))/(r-rcut)**2
    END IF

  END FUNCTION siep_f2_d

! *****************************************************************************
!> \brief exponential part of three-body term, see Equation (4) in J. Chem. 
!>        Phys., Vol. 102, p.511
!> \param siepmann ...
!> \param r_last_update_pbc ...
!> \param iparticle ...
!> \param jparticle ...
!> \param n_loc_size ...
!> \param full_loc_list ...
!> \param cell_v ...
!> \param rcutsq ...
!> \param particle_set ...
!> \param cell ...
!> \retval siep_a_ij ...
!> \par History
!>      Using a local list of neighbors
! *****************************************************************************
  FUNCTION siep_a_ij(siepmann,r_last_update_pbc,iparticle,jparticle,n_loc_size,&
       full_loc_list, cell_v,rcutsq,particle_set,cell)
    TYPE(siepmann_pot_type), POINTER         :: siepmann
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    INTEGER, INTENT(IN)                      :: iparticle, jparticle, &
                                                n_loc_size
    INTEGER, DIMENSION(2, 1:n_loc_size)      :: full_loc_list
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    REAL(KIND=dp), INTENT(IN)                :: rcutsq
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp)                            :: siep_a_ij

    CHARACTER(LEN=*), PARAMETER :: routineN = 'siep_a_ij', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: ilist, kparticle
    REAL(KIND=dp)                            :: costheta, drji, drjk, F, &
                                                rab2_max, rji(3), rjk(3), &
                                                theta

    siep_a_ij = 0.0_dp
    rab2_max = rcutsq
    F = siepmann%F
    CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,&
                        element_symbol=element_symbol)
    IF (element_symbol /= "O") RETURN
    rji(:) = -1.0_dp*(r_last_update_pbc(jparticle)%r(:)-r_last_update_pbc(iparticle)%r(:)+cell_v)
    drji = SQRT(DOT_PRODUCT(rji,rji))
    DO ilist = 1, n_loc_size
       kparticle = full_loc_list(2,ilist)
       IF (kparticle==jparticle) CYCLE
       rjk(:) = pbc(r_last_update_pbc(jparticle)%r(:),r_last_update_pbc(kparticle)%r(:),cell) 
       drjk = DOT_PRODUCT(rjk,rjk)
       IF(drjk>rab2_max) CYCLE
       drjk = SQRT(drjk)
       costheta = DOT_PRODUCT(rji,rjk) / (drji*drjk)
       IF(costheta<-1.0_dp) costheta=-1.0_dp 
       IF(costheta>+1.0_dp) costheta=+1.0_dp 
       theta = ACOS(costheta)
       siep_a_ij = siep_a_ij + EXP(F*(COS(theta/2.0_dp))**2)
    END DO
  END FUNCTION siep_a_ij

! *****************************************************************************
!> \brief derivative of a_ij
!> \param siepmann ...
!> \param r_last_update_pbc ...
!> \param iparticle ...
!> \param jparticle ...
!> \param f_nonbond ...
!> \param pv_nonbond ...
!> \param prefactor ...
!> \param n_loc_size ...
!> \param full_loc_list ...
!> \param loc_cell_v ...
!> \param cell_v ...
!> \param cell ...
!> \param rcutsq ...
!> \param use_virial ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>       Using a local list of neighbors
! *****************************************************************************
  SUBROUTINE siep_a_ij_d(siepmann,r_last_update_pbc,iparticle,jparticle,f_nonbond,&
                         pv_nonbond,prefactor,n_loc_size,full_loc_list,loc_cell_v,&
                         cell_v,cell,rcutsq,use_virial, error)
    TYPE(siepmann_pot_type), POINTER         :: siepmann
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
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(IN)                :: rcutsq
    LOGICAL, INTENT(IN)                      :: use_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'siep_a_ij_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilist, kparticle, nparticle
    REAL(KIND=dp)                            :: costheta, d_expterm, &
                                                dcos_thetahalf, drji, drjk, &
                                                F, rab2_max, theta
    REAL(KIND=dp), DIMENSION(3)              :: dcosdri, dcosdrj, dcosdrk, &
                                                dri, drj, drk, rji, rji_hat, &
                                                rjk, rjk_hat

    rab2_max = rcutsq
    F=siepmann%F
 
    rji(:) = -1.0_dp*(r_last_update_pbc(jparticle)%r(:)-r_last_update_pbc(iparticle)%r(:)+cell_v)
    drji = SQRT(DOT_PRODUCT(rji,rji))
    rji_hat(:) = rji(:) / drji
 
    nparticle = SIZE ( r_last_update_pbc )
    DO ilist = 1, n_loc_size
       kparticle = full_loc_list(2,ilist)
       IF (kparticle==jparticle) CYCLE
       rjk(:) = pbc(r_last_update_pbc(jparticle)%r(:),r_last_update_pbc(kparticle)%r(:),cell) 
       drjk = DOT_PRODUCT(rjk,rjk)
       IF(drjk>rab2_max) CYCLE
       drjk = SQRT(drjk)
       rjk_hat(:) = rjk(:) / drjk
       costheta = DOT_PRODUCT(rji,rjk) / (drji*drjk)
       IF(costheta<-1.0_dp) costheta=-1.0_dp
       IF(costheta>+1.0_dp) costheta=+1.0_dp

       dcosdri(:) = (1.0_dp/(drji)) * ( rjk_hat(:) - costheta*rji_hat(:) )
       dcosdrk(:) = (1.0_dp/(drjk)) * ( rji_hat(:) - costheta*rjk_hat(:) )
       dcosdrj(:) = - ( dcosdri(:) + dcosdrk(:) )

       theta  = ACOS(costheta)
       dcos_thetahalf = -1.0_dp/(2.0_dp*SQRT(1-costheta**2))
       d_expterm = -2.0_dp * F * COS(theta/2.0_dp) * SIN (theta/2.0_dp)&
                   * EXP(F *(COS(theta/2.0_dp))**2) 
       
        
       dri = d_expterm * dcos_thetahalf  * dcosdri 

       drj = d_expterm * dcos_thetahalf  * dcosdrj 

       drk = d_expterm * dcos_thetahalf  * dcosdrk
 
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
        CALL cp_unimplemented_error(fromWhere=routineP, &
             message="using virial with Siepmann-Sprik"//&
             " not implemented", &
             error=error, error_level=cp_failure_level)
       END IF
     END DO
  END SUBROUTINE siep_a_ij_d
! *****************************************************************************
!> \brief Phi_ij corresponds to Equation (3) in J. Chem. Phys., Vol. 102, p.511 
!> \param siepmann ...
!> \param r_last_update_pbc ...
!> \param iparticle ...
!> \param jparticle ...
!> \param n_loc_size ...
!> \param loc_cell_v ...
!> \param cell_v ...
!> \param cell ...
!> \param rcutsq ...
!> \param particle_set ...
!> \param nr_oh ...
!> \retval siep_Phi_ij ...
!> \par History
!>      Using a local list of neighbors
! *****************************************************************************
  FUNCTION siep_Phi_ij(siepmann,r_last_update_pbc,iparticle,jparticle,n_loc_size,&
       loc_cell_v, cell_v, cell, rcutsq, particle_set,nr_oh)
    TYPE(siepmann_pot_type), POINTER         :: siepmann
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    INTEGER, INTENT(IN)                      :: iparticle, jparticle, &
                                                n_loc_size
    REAL(KIND=dp), &
      DIMENSION(3, 1:n_loc_size)             :: loc_cell_v
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(IN)                :: rcutsq
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(INOUT), OPTIONAL         :: nr_oh
    REAL(KIND=dp)                            :: siep_Phi_ij

    CHARACTER(LEN=*), PARAMETER :: routineN = 'siep_Phi_ij', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: count_h, iatom, index_h1, &
                                                index_h2, natom
    REAL(KIND=dp)                            :: cosphi, drih, drix, drji, &
                                                h_max_dist, rab2_max, rih(3), &
                                                rih1(3), rih2(3), rix(3), &
                                                rji(3)

    siep_Phi_ij = 0.0_dp
    count_h=0
    index_h1=0
    index_h2=0
    rab2_max = rcutsq
    h_max_dist=2.27_dp ! 1.2 angstrom
    natom=SIZE(particle_set)
    CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,&
                        element_symbol=element_symbol)
    IF (element_symbol /= "O") RETURN
    rji(:) = -1.0_dp*(r_last_update_pbc(jparticle)%r(:)-r_last_update_pbc(iparticle)%r(:)+cell_v)
    drji = SQRT(DOT_PRODUCT(rji,rji))

    DO iatom=1,natom
       CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                          element_symbol=element_symbol)
       IF(element_symbol /= "H" ) CYCLE
       rih(:)= pbc(r_last_update_pbc(iparticle)%r(:),r_last_update_pbc(iatom)%r(:),cell)
       drih = SQRT(DOT_PRODUCT(rih,rih))
       IF(drih >= h_max_dist) CYCLE
       count_h=count_h+1 
       IF(count_h==1) THEN
        index_h1= iatom
       ELSEIF(count_h==2) THEN
        index_h2=iatom
       ELSE
        CALL stop_program(routineN,moduleN,__LINE__,&
             "Error in Siepmann-Sprik part: too many H atoms")
       ENDIF 
    ENDDO

    IF (count_h==0) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
             "No H atoms for O found")
    ELSEIF (count_h==1) THEN
        IF(siepmann%allow_oh_formation) THEN
           IF(PRESENT(nr_oh)) THEN
              nr_oh = nr_oh + 1
           ENDIF
           siep_Phi_ij = 0.0_dp
        ELSE
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Only one H atom of O atom found")
        ENDIF
    ENDIF
 
    IF(count_h==2) THEN
     !dipole vector rix of the H2O molecule 
     rih1(:)= pbc(r_last_update_pbc(iparticle)%r(:),r_last_update_pbc(index_h1)%r(:),cell)
     rih2(:)= pbc(r_last_update_pbc(iparticle)%r(:),r_last_update_pbc(index_h2)%r(:),cell)
     rix(:) = rih1(:) + rih2(:) 
     drix= SQRT(DOT_PRODUCT(rix,rix))
     cosphi=DOT_PRODUCT(rji,rix)/(drji*drix)
     IF(cosphi<-1.0_dp) cosphi=-1.0_dp 
     IF(cosphi>+1.0_dp) cosphi=+1.0_dp
     siep_Phi_ij=EXP(-8.0_dp*((cosphi-1)/4.0_dp)**4)
    ENDIF
  END FUNCTION siep_Phi_ij

! *****************************************************************************
!> \brief derivative of Phi_ij
!> \param siepmann ...
!> \param r_last_update_pbc ...
!> \param iparticle ...
!> \param jparticle ...
!> \param f_nonbond ...
!> \param prefactor ...
!> \param n_loc_size ...
!> \param full_loc_list ...
!> \param loc_cell_v ...
!> \param cell_v ...
!> \param cell ...
!> \param rcutsq ...
!> \param use_virial ...
!> \param particle_set ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>       Using a local list of neighbors
! *****************************************************************************
  SUBROUTINE siep_Phi_ij_d(siepmann,r_last_update_pbc,iparticle,jparticle,f_nonbond,&
                           prefactor,n_loc_size,full_loc_list, loc_cell_v,&
                           cell_v, cell, rcutsq,use_virial, particle_set, error)
    TYPE(siepmann_pot_type), POINTER         :: siepmann
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    INTEGER, INTENT(IN)                      :: iparticle, jparticle
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_nonbond
    REAL(KIND=dp), INTENT(IN)                :: prefactor
    INTEGER, INTENT(IN)                      :: n_loc_size
    INTEGER, DIMENSION(2, 1:n_loc_size)      :: full_loc_list
    REAL(KIND=dp), &
      DIMENSION(3, 1:n_loc_size)             :: loc_cell_v
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(IN)                :: rcutsq
    LOGICAL, INTENT(IN)                      :: use_virial
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'siep_Phi_ij_d', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: count_h, iatom, index_h1, &
                                                index_h2, natom
    REAL(KIND=dp)                            :: cosphi, dphi, drih, drix, &
                                                drji, h_max_dist, Phi_ij, &
                                                rab2_max
    REAL(KIND=dp), DIMENSION(3)              :: dcosdrh, dcosdri, dcosdrj, &
                                                drh, dri, drj, rih, rih1, &
                                                rih2, rix, rix_hat, rji, &
                                                rji_hat

    count_h=0
    index_h1=0
    index_h2=0
    rab2_max = rcutsq
    h_max_dist=2.27_dp ! 1.2 angstrom
    natom=SIZE(particle_set)
    Phi_ij = siep_Phi_ij(siepmann,r_last_update_pbc,iparticle,jparticle,&
             n_loc_size, loc_cell_v, cell_v, cell, rcutsq,&
             particle_set) 
    rji(:) = -1.0_dp*(r_last_update_pbc(jparticle)%r(:)-r_last_update_pbc(iparticle)%r(:)+cell_v)
    drji = SQRT(DOT_PRODUCT(rji,rji))
    rji_hat(:) = rji(:) / drji

    DO iatom=1,natom
       CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                          element_symbol=element_symbol)
       IF(element_symbol /= "H" ) CYCLE
       rih(:)= pbc(r_last_update_pbc(iparticle)%r(:),r_last_update_pbc(iatom)%r(:),cell)
       drih = SQRT(DOT_PRODUCT(rih,rih))
       IF(drih >= h_max_dist) CYCLE
       count_h=count_h+1 
       IF(count_h==1) THEN
        index_h1= iatom
       ELSEIF(count_h==2) THEN
        index_h2=iatom
       ELSE
        CALL stop_program(routineN,moduleN,__LINE__,&
             "Error in Siepmann-Sprik part: too many H atoms")
       ENDIF 
    ENDDO

    IF (count_h==0) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
             "No H atoms for O found")
    ELSEIF (count_h==1.AND..NOT.siepmann%allow_oh_formation) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
             "Only one H atom of O atom found")
    ENDIF 

    IF(count_h==2) THEN
     !dipole vector rix of the H2O molecule
     rih1(:)= pbc(r_last_update_pbc(iparticle)%r(:),r_last_update_pbc(index_h1)%r(:),cell)
     rih2(:)= pbc(r_last_update_pbc(iparticle)%r(:),r_last_update_pbc(index_h2)%r(:),cell)
     rix(:) = rih1(:) + rih2(:) 
     drix= SQRT(DOT_PRODUCT(rix,rix))
     rix_hat(:) = rix(:) / drix
     cosphi=DOT_PRODUCT(rji,rix)/(drji*drix)
     IF(cosphi<-1.0_dp) cosphi=-1.0_dp 
     IF(cosphi>+1.0_dp) cosphi=+1.0_dp
 
     dcosdrj(:) = (1.0_dp/(drji)) * (-rix_hat(:) + cosphi*rji_hat(:) )
     ! for H atoms:
     dcosdrh(:) = (1.0_dp/(drix)) * ( rji_hat(:) - cosphi*rix_hat(:) )
     dcosdri(:) = - dcosdrj - 2.0_dp * dcosdrh   
     
     dphi = Phi_ij * (-8.0_dp) * ( (cosphi -1)/4.0_dp)**3 
     
     dri = dphi * dcosdri 
     drj = dphi * dcosdrj 
     drh = dphi * dcosdrh 
  
     f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + prefactor*dri(1)
     f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + prefactor*dri(2)
     f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + prefactor*dri(3)
  
     f_nonbond(1,jparticle) = f_nonbond(1,jparticle) + prefactor*drj(1)
     f_nonbond(2,jparticle) = f_nonbond(2,jparticle) + prefactor*drj(2)
     f_nonbond(3,jparticle) = f_nonbond(3,jparticle) + prefactor*drj(3)
  
     f_nonbond(1,index_h1) = f_nonbond(1,index_h1) + prefactor*drh(1)
     f_nonbond(2,index_h1) = f_nonbond(2,index_h1) + prefactor*drh(2)
     f_nonbond(3,index_h1) = f_nonbond(3,index_h1) + prefactor*drh(3)

     f_nonbond(1,index_h2) = f_nonbond(1,index_h2) + prefactor*drh(1)
     f_nonbond(2,index_h2) = f_nonbond(2,index_h2) + prefactor*drh(2)
     f_nonbond(3,index_h2) = f_nonbond(3,index_h2) + prefactor*drh(3)

     IF (use_virial) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
           message="using virial with Siepmann-Sprik"//&
           " not implemented", &
           error=error, error_level=cp_failure_level)
     END IF

    ENDIF
  END SUBROUTINE siep_Phi_ij_d

! *****************************************************************************
!> \brief forces generated by the three-body term
!> \param siepmann ...
!> \param r_last_update_pbc ...
!> \param cell_v ...
!> \param n_loc_size ...
!> \param full_loc_list ...
!> \param loc_cell_v ...
!> \param iparticle ...
!> \param jparticle ...
!> \param f_nonbond ...
!> \param pv_nonbond ...
!> \param use_virial ...
!> \param rcutsq ...
!> \param cell ...
!> \param particle_set ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>       Using a local list of neighbors
! *****************************************************************************
  SUBROUTINE siepmann_forces_v3(siepmann,r_last_update_pbc, cell_v, n_loc_size,&
       full_loc_list, loc_cell_v, iparticle,jparticle,f_nonbond,pv_nonbond,&
       use_virial, rcutsq, cell, particle_set, error)
    TYPE(siepmann_pot_type), POINTER         :: siepmann
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
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'siepmann_forces_v3', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    REAL(KIND=dp)                            :: a_ij, beta, drji, E, f2, &
                                                f2_d, f_A1, f_A2, fac, &
                                                prefactor, rji(3), rji_hat(3)

    beta=siepmann%beta
    E=siepmann%E

    CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,&
                        element_symbol=element_symbol)
    IF (element_symbol /= "O") RETURN

    rji(:) = -1.0_dp*(r_last_update_pbc(jparticle)%r(:)-r_last_update_pbc(iparticle)%r(:)+cell_v)
    drji = SQRT(DOT_PRODUCT(rji,rji))
    rji_hat(:) = rji(:) / drji

    fac = - 1.0_dp  !gradient to force
    a_ij  = siep_a_ij (siepmann,r_last_update_pbc,iparticle,jparticle,n_loc_size,&
                       full_loc_list, cell_v,rcutsq,particle_set,cell)
    f2   = siep_f2   (siepmann,drji)
    f2_d = siep_f2_d (siepmann,drji)

    ! Lets do the f_A1 piece derivative of  f2
    f_A1 = E * f2_d * drji**(-beta) * a_ij * fac *(1.0_dp/drji) 
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_A1*rji(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_A1*rji(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_A1*rji(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_A1*rji(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_A1*rji(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_A1*rji(3)
 
    IF (use_virial) THEN
     CALL cp_unimplemented_error(fromWhere=routineP, &
          message="using virial with Siepmann-Sprik"//&
          " not implemented", error=error, error_level=cp_failure_level)
    END IF

    ! Lets do the f_A2 piece derivative of rji**(-beta)
    f_A2 = E * f2 *(-beta)*drji**(-beta-1) * a_ij * fac * (1.0_dp/drji)
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_A2*rji(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_A2*rji(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_A2*rji(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_A2*rji(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_A2*rji(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_A2*rji(3)
 
    IF (use_virial) THEN
     CALL cp_unimplemented_error(fromWhere=routineP, &
          message="using virial with Siepmann-Sprik"//&
          " not implemented", error=error, error_level=cp_failure_level)
    END IF

    ! Lets do the f_A3 piece derivative: of a_ij
    prefactor = E * f2 * drji**(-beta) * fac  
    CALL siep_a_ij_d(siepmann, r_last_update_pbc, iparticle, jparticle, f_nonbond,&
         pv_nonbond, prefactor, n_loc_size, full_loc_list,loc_cell_v, cell_v,&
         cell, rcutsq, use_virial, error)

  END SUBROUTINE siepmann_forces_v3

! *****************************************************************************
!> \brief forces generated by the dipole term
!> \param siepmann ...
!> \param r_last_update_pbc ...
!> \param cell_v ...
!> \param cell ...
!> \param n_loc_size ...
!> \param full_loc_list ...
!> \param loc_cell_v ...
!> \param iparticle ...
!> \param jparticle ...
!> \param f_nonbond ...
!> \param pv_nonbond ...
!> \param use_virial ...
!> \param rcutsq ...
!> \param particle_set ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>       Using a local list of neighbors
! *****************************************************************************
  SUBROUTINE siepmann_forces_v2(siepmann,r_last_update_pbc, cell_v, cell, n_loc_size,&
       full_loc_list, loc_cell_v, iparticle,jparticle,f_nonbond,pv_nonbond,&
       use_virial, rcutsq, particle_set, error)
    TYPE(siepmann_pot_type), POINTER         :: siepmann
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: n_loc_size
    INTEGER, DIMENSION(2, 1:n_loc_size)      :: full_loc_list
    REAL(KIND=dp), &
      DIMENSION(3, 1:n_loc_size)             :: loc_cell_v
    INTEGER, INTENT(IN)                      :: iparticle, jparticle
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_nonbond, pv_nonbond
    LOGICAL, INTENT(IN)                      :: use_virial
    REAL(KIND=dp), INTENT(IN)                :: rcutsq
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'siepmann_forces_v2', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    REAL(KIND=dp)                            :: D, drji, f2, f2_d, f_A1, &
                                                f_A2, fac, Phi_ij, prefactor, &
                                                rji(3)

    D=siepmann%D

    CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,&
                        element_symbol=element_symbol)
    IF (element_symbol /= "O") RETURN

    rji(:) = -1.0_dp*(r_last_update_pbc(jparticle)%r(:)-r_last_update_pbc(iparticle)%r(:) + cell_v)
    drji = SQRT(DOT_PRODUCT(rji,rji))

    fac = - 1.0_dp
    Phi_ij = siep_Phi_ij(siepmann,r_last_update_pbc,iparticle,jparticle,n_loc_size,&
             loc_cell_v, cell_v, cell, rcutsq,particle_set)
    f2   = siep_f2   (siepmann,drji)
    f2_d = siep_f2_d (siepmann,drji)

    ! Lets do the f_A1 piece derivative of  f2
    f_A1 = -D * f2_d * drji**(-3) * Phi_ij * fac *(1.0_dp/drji) 
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_A1*rji(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_A1*rji(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_A1*rji(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_A1*rji(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_A1*rji(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_A1*rji(3)
 
    IF (use_virial) THEN
     CALL cp_unimplemented_error(fromWhere=routineP, &
          message="using virial with Siepmann-Sprik"//&
          " not implemented", error=error, error_level=cp_failure_level)
    END IF

!   ! Lets do the f_A2 piece derivative of rji**(-3)
    f_A2 = -D * f2 *(-3.0_dp)*drji**(-4) * Phi_ij * fac * (1.0_dp/drji)
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_A2*rji(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_A2*rji(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_A2*rji(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_A2*rji(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_A2*rji(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_A2*rji(3)
 
    IF (use_virial) THEN
     CALL cp_unimplemented_error(fromWhere=routineP, &
          message="using virial with Siepmann-Sprik"//&
          " not implemented", error=error, error_level=cp_failure_level)
    END IF

    ! Lets do the f_A3 piece derivative: of Phi_ij
    prefactor = -D * f2 * drji**(-3) * fac  
    CALL siep_Phi_ij_d(siepmann, r_last_update_pbc, iparticle, jparticle, f_nonbond,&
         prefactor, n_loc_size, full_loc_list, loc_cell_v,cell_v, cell,&
         rcutsq, use_virial, particle_set, error)

  END SUBROUTINE siepmann_forces_v2

! *****************************************************************************
!> \brief ...
!> \param nonbonded ...
!> \param potparm ...
!> \param glob_loc_list ...
!> \param glob_cell_v ...
!> \param glob_loc_list_a ...
!> \param cell ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
! *****************************************************************************
  SUBROUTINE setup_siepmann_arrays(nonbonded, potparm, glob_loc_list, glob_cell_v,&
                                   glob_loc_list_a, cell, error)
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    INTEGER, DIMENSION(:, :), POINTER        :: glob_loc_list
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: glob_cell_v
    INTEGER, DIMENSION(:), POINTER           :: glob_loc_list_a
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_siepmann_arrays', &
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
             IF (pot%type(i)==siepmann_type) npairs_tot = npairs_tot + npairs
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
             ! SIEPMANN
             IF (pot%type(i)==siepmann_type) THEN
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
  END SUBROUTINE setup_siepmann_arrays

! *****************************************************************************
!> \brief ...
!> \param glob_loc_list ...
!> \param glob_cell_v ...
!> \param glob_loc_list_a ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE destroy_siepmann_arrays(glob_loc_list, glob_cell_v, glob_loc_list_a, error)
    INTEGER, DIMENSION(:, :), POINTER        :: glob_loc_list
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: glob_cell_v
    INTEGER, DIMENSION(:), POINTER           :: glob_loc_list_a
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'destroy_siepmann_arrays', &
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

  END SUBROUTINE destroy_siepmann_arrays

! *****************************************************************************
!> \brief prints the number of OH ions
!> \param nr_oh ...
!> \param mm_section ...
!> \param para_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE print_nr_oh_siepmann(nr_oh,mm_section,para_env,error)
    INTEGER, INTENT(INOUT)                   :: nr_oh
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_nr_oh_siepmann', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iw
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)

    CALL mp_sum(nr_oh,para_env%group)
    logger => cp_error_get_logger(error)

    iw= cp_print_key_unit_nr(logger,mm_section,"PRINT%PROGRAM_RUN_INFO",&
        extension=".mmLog",error=error)

    IF(iw > 0 .AND. nr_oh > 0) THEN
      WRITE(iw,'(/,A,T71,I10,/)') " SIEPMANN: number of OH ions at surface", nr_oh
    ENDIF

    CALL cp_print_key_finished_output(iw,logger,mm_section,"PRINT%PROGRAM_RUN_INFO",&
         error=error)

  END SUBROUTINE print_nr_oh_siepmann

END MODULE manybody_siepmann

