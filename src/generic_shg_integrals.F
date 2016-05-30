!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of contracted, spherical Gaussian integrals using the solid harmonic
!>        Gaussian (SHG) integral scheme. Routines for the following two-center integrals:
!>        i)  (a|O(r12)|b) where O(r12) is the overlap, coulomb operator etc.
!>        ii) (aba) and (abb) s-overlaps
!> \par Literature
!>      T.J. Giese and D. M. York, J. Chem. Phys, 128, 064104 (2008)
!>      T. Helgaker, P Joergensen, J. Olsen, Molecular Electronic-Structure
!>                                           Theory, Wiley
!> \par History
!>      created [05.2016]
!> \author Dorothea Golze
! **************************************************************************************************
MODULE generic_shg_integrals
   USE ai_construct_shg,                ONLY: &
        construct_dev_shg_ab, construct_int_shg_ab, construct_overlap_shg_aba, &
        construct_overlap_shg_abb, dev_overlap_shg_aba, dev_overlap_shg_abb, get_W_matrix, &
        get_dW_matrix, get_real_scaled_solid_harmonic
   USE ai_s_contract_shg,               ONLY: &
        contract_s_overlap_aba, contract_s_overlap_abb, contract_sint_ab_chigh, &
        contract_sint_ab_clow, s_coulomb_ab, s_gauss_ab, s_overlap_ab, s_overlap_abx, s_verf_ab, &
        s_verfc_ab, s_vgauss_ab
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE input_constants,                 ONLY: operator_coulomb,&
                                              operator_gauss,&
                                              operator_verf,&
                                              operator_verfc,&
                                              operator_vgauss
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: nsoset
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'generic_shg_integrals'

   PUBLIC :: two_center_shg_integrals, lri_precalc_angular_shg_part, &
             int_overlap_ab_shg, int_overlap_aba_shg, int_overlap_abb_shg

   ABSTRACT INTERFACE
! **************************************************************************************************
!> \brief Interface for the calculation of integrals over s-functions and their scalar derivatives
!>        with respect to rab2
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param omega ...
!> \param rab ...
!> \param v matrix storing the integrals and scalar derivatives
!> \param calculate_forces ...
! **************************************************************************************************
      SUBROUTINE ab_sint_shg(la_max, npgfa, zeta, lb_max, npgfb, zetb, omega, rab, v, calculate_forces)
      USE kinds,                        ONLY: dp
      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      LOGICAL, INTENT(IN)                                :: calculate_forces

      END SUBROUTINE ab_sint_shg
   END INTERFACE

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Calcululates the two-center integrals of the type (a|O(r12)|b) using the SHG scheme
!> \param r12_operator the integral operator, which depends on r12=|r1-r2|
!> \param vab integral matrix of spherical contracted Gaussian functions
!> \param dvab derivative of the integrals
!> \param rab distance vector between center A and B
!> \param fba basis at center A
!> \param fbb basis at center B
!> \param scona_shg SHG contraction matrix for A
!> \param sconb_shg SHG contraction matrix for B
!> \param omega parameter in the operator
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE two_center_shg_integrals(r12_operator, vab, dvab, rab, fba, fbb, scona_shg, sconb_shg, &
                                       omega, calculate_forces)

      INTEGER, INTENT(IN)                                :: r12_operator
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vab
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: dvab
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: fba, fbb
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scona_shg, sconb_shg
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: omega
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'two_center_shg_integrals', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: my_omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Waux_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: dWaux_mat

      PROCEDURE(ab_sint_shg), POINTER                    :: s_operator_ab => null()

      CALL timeset(routineN, handle)

      CALL precalc_angular_shg_part(fba, fbb, rab, Waux_mat, dWaux_mat, calculate_forces)
      my_omega = 1.0_dp

      SELECT CASE (r12_operator)
      CASE (operator_coulomb)
         s_operator_ab => s_coulomb_ab
      CASE (operator_verf)
         s_operator_ab => s_verf_ab
         IF (PRESENT(omega)) my_omega = omega
      CASE (operator_verfc)
         s_operator_ab => s_verfc_ab
         IF (PRESENT(omega)) my_omega = omega
      CASE (operator_vgauss)
         s_operator_ab => s_vgauss_ab
         IF (PRESENT(omega)) my_omega = omega
      CASE (operator_gauss)
         s_operator_ab => s_gauss_ab
         IF (PRESENT(omega)) my_omega = omega
      CASE DEFAULT
         CPABORT("Operator not available")
      END SELECT

      CALL int_operator_ab_shg(s_operator_ab, vab, dvab, rab, fba, fbb, scona_shg, sconb_shg, &
                               my_omega, Waux_mat, dWaux_mat, calculate_forces)

      DEALLOCATE (Waux_mat, dWaux_mat)

      CALL timestop(handle)

   END SUBROUTINE two_center_shg_integrals

! **************************************************************************************************
!> \brief precalculates the angular part of the SHG integrals
!> \param fba basis at center A
!> \param fbb basis at center B
!> \param rab distance vector between a and b
!> \param Waux_mat W matrix that contains the angular-dependent part
!> \param dWaux_mat derivative of the W matrix
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE precalc_angular_shg_part(fba, fbb, rab, Waux_mat, dWaux_mat, calculate_forces)

      TYPE(gto_basis_set_type), POINTER                  :: fba, fbb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: Waux_mat
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :), INTENT(OUT)              :: dWaux_mat
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'precalc_angular_shg_part', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, la_max, lb_max, lmax, mdim(3)
      INTEGER, DIMENSION(:), POINTER                     :: la_max_all
      REAL(KIND=dp)                                      :: rab2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Rc, Rs

      CALL timeset(routineN, handle)

      NULLIFY (la_max_all)
      rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)

      la_max = MAXVAL(fba%lmax)
      lb_max = MAXVAL(fbb%lmax)
      lmax = MAX(la_max, lb_max)

      ALLOCATE (la_max_all(0:lb_max))
      ALLOCATE (Rc(0:lmax, -2*lmax:2*lmax), Rs(0:lmax, -2*lmax:2*lmax))
      Rc = 0._dp
      Rs = 0._dp
      mdim(1) = MIN(la_max, lb_max)+1
      mdim(2) = nsoset(la_max)+1
      mdim(3) = nsoset(lb_max)+1
      ALLOCATE (Waux_mat(mdim(1), mdim(2), mdim(3)))
      ALLOCATE (dWaux_mat(3, mdim(1), mdim(2), mdim(3)))

      la_max_all(0:lb_max) = la_max
      !*** -rab, since Eq. in Ref. use Ra-Rb, not Rb-Ra
      CALL get_real_scaled_solid_harmonic(Rc, Rs, lmax, -rab, rab2)
      CALL get_W_matrix(la_max_all, lb_max, lmax, Rc, Rs, Waux_mat)
      IF (calculate_forces) THEN
         CALL get_dW_matrix(la_max_all, lb_max, Waux_mat, dWaux_mat)
      ENDIF

      DEALLOCATE (Rc, Rs, la_max_all)
      CALL timestop(handle)

   END SUBROUTINE precalc_angular_shg_part

! **************************************************************************************************
!> \brief calculate integrals (a|O(r12)|b)
!> \param s_operator_ab procedure pointer for the respective operator. The integral evaluation
!>        differs only in the calculation of the [s|O(r12)|s] integrals and their scalar
!>        derivatives.
!> \param vab integral matrix of spherical contracted Gaussian functions
!> \param dvab derivative of the integrals
!> \param rab distance vector between center A and B
!> \param fba basis at center A
!> \param fbb basis at center B
!> \param scona_shg SHG contraction matrix for A
!> \param sconb_shg SHG contraction matrix for B
!> \param omega parameter in the operator
!> \param Waux_mat W matrix that contains the angular-dependent part
!> \param dWaux_mat derivative of the W matrix
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE int_operator_ab_shg(s_operator_ab, vab, dvab, rab, fba, fbb, scona_shg, sconb_shg, &
                                  omega, Waux_mat, dWaux_mat, calculate_forces)

      PROCEDURE(ab_sint_shg), POINTER                    :: s_operator_ab
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vab
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: dvab
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: fba, fbb
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scona_shg, sconb_shg
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Waux_mat
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: dWaux_mat
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'int_operator_ab_shg', routineP = moduleN//':'//routineN

      INTEGER :: handle, iset, jset, la_max_set, lb_max_set, ndev, nds, nds_max, npgfa_set, &
                 npgfb_set, nseta, nsetb, nsgfa_set, nsgfb_set, nshella_set, nshellb_set, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, lb_max, npgfa, npgfb, nsgfa, &
                                                            nsgfb, nshella, nshellb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb, la, lb
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zeta, zetb
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: swork, swork_cont

      CALL timeset(routineN, handle)
      NULLIFY (la_max, lb_max, npgfa, npgfb, first_sgfa, first_sgfb, set_radius_a, &
               set_radius_b, zeta, zetb)

      ! basis ikind
      first_sgfa => fba%first_sgf
      la_max => fba%lmax
      la => fba%l
      npgfa => fba%npgf
      nsgfa => fba%nsgf_set
      nseta = fba%nset
      set_radius_a => fba%set_radius
      zeta => fba%zet
      nshella => fba%nshell
      ! basis jkind
      first_sgfb => fbb%first_sgf
      lb_max => fbb%lmax
      lb => fbb%l
      npgfb => fbb%npgf
      nsgfb => fbb%nsgf_set
      nsetb = fbb%nset
      set_radius_b => fbb%set_radius
      zetb => fbb%zet
      nshellb => fbb%nshell

      dab = SQRT(SUM(rab**2))

      la_max_set = MAXVAL(la_max)
      lb_max_set = MAXVAL(lb_max)

      ! allocate some work matrices
      npgfa_set = MAXVAL(npgfa)
      npgfb_set = MAXVAL(npgfb)
      nshella_set = MAXVAL(nshella)
      nshellb_set = MAXVAL(nshellb)
      nsgfa_set = MAXVAL(nsgfa)
      nsgfb_set = MAXVAL(nsgfb)
      ndev = 0
      IF (calculate_forces) ndev = 1
      nds_max = la_max_set+lb_max_set+ndev+1
      ALLOCATE (swork(npgfa_set, npgfb_set, nds_max))
      ALLOCATE (swork_cont(nds_max, nshella_set, nshellb_set))

      vab = 0.0_dp
      dvab = 0.0_dp

      DO iset = 1, nseta
         sgfa = first_sgfa(1, iset)

         DO jset = 1, nsetb

            sgfb = first_sgfb(1, jset)
            CALL s_operator_ab(la_max(iset), npgfa(iset), zeta(:, iset), &
                               lb_max(jset), npgfb(jset), zetb(:, jset), &
                               omega, rab, swork, calculate_forces)
            nds = la_max(iset)+lb_max(jset)+ndev+1
            CALL contract_sint_ab_chigh(npgfa(iset), nshella(iset), scona_shg(1:npgfa(iset), 1:nshella(iset), iset), &
                                        npgfb(jset), nshellb(jset), sconb_shg(1:npgfb(jset), 1:nshellb(jset), jset), &
                                        nds, swork(1:npgfa(iset), 1:npgfb(jset), 1:nds), &
                                        swork_cont(1:nds, 1:nshella(iset), 1:nshellb(jset)))
            CALL construct_int_shg_ab(la(:, iset), first_sgfa(:, iset), nshella(iset), &
                                      lb(:, jset), first_sgfb(:, jset), nshellb(jset), &
                                      swork_cont, Waux_mat, vab)
            IF (calculate_forces) THEN
               !*** -rab, since Eq. in Ref. use Ra-Rb, not Rb-Ra
               CALL construct_dev_shg_ab(la(:, iset), first_sgfa(:, iset), nshella(iset), &
                                         lb(:, jset), first_sgfb(:, jset), nshellb(jset), &
                                         -rab, swork_cont, Waux_mat, dWaux_mat, dvab)
            ENDIF
         END DO
      END DO

      DEALLOCATE (swork, swork_cont)

      CALL timestop(handle)
   END SUBROUTINE int_operator_ab_shg

! **************************************************************************************************
!> \brief calculate overlap integrals (a,b)
!> \param sab integral (a,b)
!> \param dsab derivative of sab
!> \param rab distance vector
!> \param fba basis at center A
!> \param fbb basis at center B
!> \param scona_shg contraction matrix A
!> \param sconb_shg contraxtion matrix B
!> \param Waux_mat W matrix that contains the angular-dependent part
!> \param dWaux_mat derivative of the W matrix
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE int_overlap_ab_shg(sab, dsab, rab, fba, fbb, scona_shg, sconb_shg, Waux_mat, dWaux_mat, &
                                 calculate_forces)

      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sab
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: dsab
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: fba, fbb
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scona_shg, sconb_shg
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Waux_mat
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: dWaux_mat
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'int_overlap_ab_shg', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iset, jset, la_max_set, &
                                                            lb_max_set, ndev, nds_max, npgfa_set, &
                                                            npgfb_set, nseta, nsetb, nshella_set, &
                                                            nshellb_set
      INTEGER, DIMENSION(:), POINTER                     :: la_max, lb_max, npgfa, npgfb, nshella, &
                                                            nshellb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb, la, lb
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: swork, swork_cont
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zeta, zetb

      CALL timeset(routineN, handle)
      NULLIFY (la_max, lb_max, npgfa, npgfb, first_sgfa, first_sgfb, set_radius_a, &
               set_radius_b, zeta, zetb)

      ! basis ikind
      first_sgfa => fba%first_sgf
      la_max => fba%lmax
      la => fba%l
      npgfa => fba%npgf
      nseta = fba%nset
      set_radius_a => fba%set_radius
      zeta => fba%zet
      nshella => fba%nshell
      ! basis jkind
      first_sgfb => fbb%first_sgf
      lb_max => fbb%lmax
      lb => fbb%l
      npgfb => fbb%npgf
      nsetb = fbb%nset
      set_radius_b => fbb%set_radius
      zetb => fbb%zet
      nshellb => fbb%nshell

      dab = SQRT(SUM(rab**2))

      la_max_set = MAXVAL(la_max)
      lb_max_set = MAXVAL(lb_max)

      ! allocate some work matrices
      npgfa_set = MAXVAL(npgfa)
      npgfb_set = MAXVAL(npgfb)
      nshella_set = MAXVAL(nshella)
      nshellb_set = MAXVAL(nshellb)
      ndev = 0
      IF (calculate_forces) ndev = 1
      nds_max = la_max_set+lb_max_set+ndev+1
      ALLOCATE (swork(npgfa_set, npgfb_set, nds_max))
      ALLOCATE (swork_cont(nds_max, nshella_set, nshellb_set))

      DO iset = 1, nseta

         DO jset = 1, nsetb

            IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

            CALL s_overlap_ab(la_max(iset), npgfa(iset), zeta(:, iset), &
                              lb_max(jset), npgfb(jset), zetb(:, jset), &
                              rab, swork, calculate_forces)
            CALL contract_sint_ab_clow(la(:, iset), npgfa(iset), nshella(iset), &
                                       zeta(:, iset), scona_shg(:, :, iset), &
                                       lb(:, jset), npgfb(jset), nshellb(jset), &
                                       zetb(:, jset), sconb_shg(:, :, jset), &
                                       swork, swork_cont, calculate_forces)
            CALL construct_int_shg_ab(la(:, iset), first_sgfa(:, iset), nshella(iset), &
                                      lb(:, jset), first_sgfb(:, jset), nshellb(jset), &
                                      swork_cont, Waux_mat, sab)
            IF (calculate_forces) THEN
               !*** -rab, since Eq. in Ref. use Ra-Rb, not Rb-Ra
               CALL construct_dev_shg_ab(la(:, iset), first_sgfa(:, iset), nshella(iset), &
                                         lb(:, jset), first_sgfb(:, jset), nshellb(jset), &
                                         -rab, swork_cont, Waux_mat, dWaux_mat, dsab)
            ENDIF
         END DO
      END DO

      DEALLOCATE (swork, swork_cont)

      CALL timestop(handle)

   END SUBROUTINE int_overlap_ab_shg

! **************************************************************************************************
!> \brief calculate integrals (a,b,fa)
!> \param abaint integral (a,b,fa)
!> \param dabdaint ...
!> \param rab distance vector between A and B
!> \param oba orbital basis at center A
!> \param obb orbital basis at center B
!> \param fba auxiliary basis set at center A
!> \param scona_shg contraction matrix for orb bas on A
!> \param sconb_shg contraction matrix for orb bas on B
!> \param sconca_shg contraction matrix for aux bas on A
!> \param cg_coeff Clebsch-Gordon coefficients
!> \param cg_none0_list list of none-zero Clebsch-Gordon coefficients
!> \param ncg_none0 number of non-zero Clebsch-Gordon coefficients
!> \param Waux_mat W matrix that contains the angular-dependent part
!> \param dWaux_mat derivative of the W matrix
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE int_overlap_aba_shg(abaint, dabdaint, rab, oba, obb, fba, scona_shg, &
                                  sconb_shg, sconca_shg, cg_coeff, cg_none0_list, ncg_none0, &
                                  Waux_mat, dWaux_mat, calculate_forces)

      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: abaint
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: dabdaint
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: oba, obb, fba
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: scona_shg, sconb_shg, sconca_shg, &
                                                            cg_coeff
      INTEGER, DIMENSION(:, :, :), INTENT(IN)            :: cg_none0_list
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: ncg_none0
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Waux_mat
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: dWaux_mat
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'int_overlap_aba_shg', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, iset, jset, kset, la_max_set, laa_max_set, lb_max_set, lca_max_set, ndev, &
         nds_max, nl_set, npgfa_set, npgfb_set, npgfca_set, nseta, nsetb, nsetca, nshella_set, &
         nshellb_set, nshellca_set
      INTEGER, DIMENSION(:), POINTER                     :: la_max, lb_max, lca_max, npgfa, npgfb, &
                                                            npgfca, nshella, nshellb, nshellca
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb, first_sgfca, la, &
                                                            lb, lca
      REAL(KIND=dp)                                      :: dab, rab2
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: swork, swork_cont
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b, set_radius_ca
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zeta, zetb, zetca

      CALL timeset(routineN, handle)
      NULLIFY (la_max, lb_max, lca_max, npgfa, npgfb, npgfca)
      NULLIFY (first_sgfa, first_sgfb, first_sgfca, set_radius_a, set_radius_b, &
               set_radius_ca, zeta, zetb, zetca)

      ! basis ikind
      first_sgfa => oba%first_sgf
      la_max => oba%lmax
      la => oba%l
      npgfa => oba%npgf
      nshella => oba%nshell
      nseta = oba%nset
      set_radius_a => oba%set_radius
      zeta => oba%zet
      ! basis jkind
      first_sgfb => obb%first_sgf
      lb_max => obb%lmax
      lb => obb%l
      npgfb => obb%npgf
      nshellb => obb%nshell
      nsetb = obb%nset
      set_radius_b => obb%set_radius
      zetb => obb%zet

      ! basis RI A
      first_sgfca => fba%first_sgf
      lca_max => fba%lmax
      lca => fba%l
      npgfca => fba%npgf
      nshellca => fba%nshell
      nsetca = fba%nset
      set_radius_ca => fba%set_radius
      zetca => fba%zet

      dab = SQRT(SUM(rab**2))
      rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)

      la_max_set = MAXVAL(la_max)
      lb_max_set = MAXVAL(lb_max)
      lca_max_set = MAXVAL(lca_max)
      npgfa_set = MAXVAL(npgfa)
      npgfb_set = MAXVAL(npgfb)
      npgfca_set = MAXVAL(npgfca)
      nshella_set = MAXVAL(nshella)
      nshellb_set = MAXVAL(nshellb)
      nshellca_set = MAXVAL(nshellca)
      !*** for forces: derivative+1 in auxiliary vector required
      ndev = 0
      IF (calculate_forces) ndev = 1

      laa_max_set = la_max_set+lca_max_set

      ! allocate some work storage....
      nds_max = laa_max_set+lb_max_set+ndev+1
      nl_set = INT((laa_max_set)/2)
      ALLOCATE (swork(npgfa_set, npgfb_set, npgfca_set, 0:nl_set, nds_max))
      ALLOCATE (swork_cont(nds_max, 0:nl_set, nshella_set, nshellb_set, nshellca_set))

      DO iset = 1, nseta

         DO jset = 1, nsetb

            IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

            DO kset = 1, nsetca

               IF (set_radius_b(jset)+set_radius_ca(kset) < dab) CYCLE

               CALL s_overlap_abx(la_max(iset), npgfa(iset), zeta(:, iset), &
                                  lb_max(jset), npgfb(jset), zetb(:, jset), &
                                  lca_max(kset), npgfca(kset), zetca(:, kset), &
                                  rab, swork, calculate_forces, calc_aba=.TRUE.)

               CALL contract_s_overlap_aba(la(:, iset), npgfa(iset), nshella(iset), zeta(:, iset), scona_shg(:, :, iset), &
                                           lb(:, jset), npgfb(jset), nshellb(jset), zetb(:, jset), sconb_shg(:, :, jset), &
                                           lca(:, kset), npgfca(kset), nshellca(kset), zetca(:, kset), sconca_shg(:, :, kset), &
                                           swork, swork_cont, calculate_forces)

               CALL construct_overlap_shg_aba(la(:, iset), first_sgfa(:, iset), nshella(iset), &
                                              lb(:, jset), first_sgfb(:, jset), nshellb(jset), &
                                              lca(:, kset), first_sgfca(:, kset), nshellca(kset), &
                                              cg_coeff, cg_none0_list, ncg_none0, &
                                              swork_cont, Waux_mat, abaint)
               IF (calculate_forces) THEN
                  !*** -rab, since Eq. in Ref. use Ra-Rb, not Rb-Ra
                  CALL dev_overlap_shg_aba(la(:, iset), first_sgfa(:, iset), nshella(iset), &
                                           lb(:, jset), first_sgfb(:, jset), nshellb(jset), &
                                           lca(:, kset), first_sgfca(:, kset), nshellca(kset), &
                                           cg_coeff, cg_none0_list, ncg_none0, &
                                           -rab, swork_cont, Waux_mat, dWaux_mat, dabdaint)
               ENDIF

            END DO
         END DO
      END DO

      DEALLOCATE (swork_cont)
      DEALLOCATE (swork)

      CALL timestop(handle)

   END SUBROUTINE int_overlap_aba_shg

! **************************************************************************************************
!> \brief calculate integrals (a,b,fb)
!> \param abbint integral (a,b,fb)
!> \param dabbint derivative of abbint
!> \param rab distance vector between A and B
!> \param oba orbital basis at center A
!> \param obb orbital basis at center B
!> \param fbb auxiliary basis set at center B
!> \param scona_shg contraction matrix for orb bas on A
!> \param sconb_shg contraction matrix for orb bas on B
!> \param sconcb_shg contraction matrix for aux bas on B
!> \param cg_coeff Clebsch-Gordon coefficients
!> \param cg_none0_list list of none-zero Clebsch-Gordon coefficients
!> \param ncg_none0 number of non-zero Clebsch-Gordon coefficients
!> \param Waux_mat W matrix that contains the angular-dependent part
!> \param dWaux_mat derivative of the W matrix
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE int_overlap_abb_shg(abbint, dabbint, rab, oba, obb, fbb, scona_shg, &
                                  sconb_shg, sconcb_shg, cg_coeff, cg_none0_list, ncg_none0, &
                                  Waux_mat, dWaux_mat, calculate_forces)

      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: abbint
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: dabbint
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: oba, obb, fbb
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: scona_shg, sconb_shg, sconcb_shg, &
                                                            cg_coeff
      INTEGER, DIMENSION(:, :, :), INTENT(IN)            :: cg_none0_list
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: ncg_none0
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Waux_mat
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: dWaux_mat
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'int_overlap_abb_shg', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, iset, jset, kset, la_max_set, lb_max_set, lbb_max_set, lcb_max_set, ndev, &
         nds_max, nl_set, npgfa_set, npgfb_set, npgfcb_set, nseta, nsetb, nsetcb, nshella_set, &
         nshellb_set, nshellcb_set
      INTEGER, DIMENSION(:), POINTER                     :: la_max, lb_max, lcb_max, npgfa, npgfb, &
                                                            npgfcb, nshella, nshellb, nshellcb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb, first_sgfcb, la, &
                                                            lb, lcb
      REAL(KIND=dp)                                      :: dab, rab2
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: swork, swork_cont
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b, set_radius_cb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zeta, zetb, zetcb

      CALL timeset(routineN, handle)
      NULLIFY (la_max, lb_max, lcb_max, npgfa, npgfb, npgfcb)
      NULLIFY (first_sgfa, first_sgfb, first_sgfcb, set_radius_a, set_radius_b, &
               set_radius_cb, zeta, zetb, zetcb)

      ! basis ikind
      first_sgfa => oba%first_sgf
      la_max => oba%lmax
      la => oba%l
      npgfa => oba%npgf
      nshella => oba%nshell
      nseta = oba%nset
      set_radius_a => oba%set_radius
      zeta => oba%zet
      ! basis jkind
      first_sgfb => obb%first_sgf
      lb_max => obb%lmax
      lb => obb%l
      npgfb => obb%npgf
      nshellb => obb%nshell
      nsetb = obb%nset
      set_radius_b => obb%set_radius
      zetb => obb%zet

      ! basis RI on B
      first_sgfcb => fbb%first_sgf
      lcb_max => fbb%lmax
      lcb => fbb%l
      npgfcb => fbb%npgf
      nshellcb => fbb%nshell
      nsetcb = fbb%nset
      set_radius_cb => fbb%set_radius
      zetcb => fbb%zet

      dab = SQRT(SUM(rab**2))
      rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)

      la_max_set = MAXVAL(la_max)
      lb_max_set = MAXVAL(lb_max)
      lcb_max_set = MAXVAL(lcb_max)
      npgfa_set = MAXVAL(npgfa)
      npgfb_set = MAXVAL(npgfb)
      npgfcb_set = MAXVAL(npgfcb)
      nshella_set = MAXVAL(nshella)
      nshellb_set = MAXVAL(nshellb)
      nshellcb_set = MAXVAL(nshellcb)
      !*** for forces: derivative+1 in auxiliary vector required
      ndev = 0
      IF (calculate_forces) ndev = 1

      lbb_max_set = lb_max_set+lcb_max_set

      ! allocate some work storage....
      nds_max = la_max_set+lbb_max_set+ndev+1
      nl_set = INT((lbb_max_set)/2)
      ALLOCATE (swork(1:npgfa_set, 1:npgfb_set, 1:npgfcb_set, 0:nl_set, 1:nds_max))
      ALLOCATE (swork_cont(nds_max, 0:nl_set, nshella_set, nshellb_set, nshellcb_set))

      DO iset = 1, nseta

         DO jset = 1, nsetb

            IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

            DO kset = 1, nsetcb

               IF (set_radius_a(iset)+set_radius_cb(kset) < dab) CYCLE

               CALL s_overlap_abx(la_max(iset), npgfa(iset), zeta(:, iset), &
                                  lb_max(jset), npgfb(jset), zetb(:, jset), &
                                  lcb_max(kset), npgfcb(kset), zetcb(:, kset), &
                                  rab, swork, calculate_forces, calc_aba=.FALSE.)

               CALL contract_s_overlap_abb(la(:, iset), npgfa(iset), nshella(iset), zeta(:, iset), scona_shg(:, :, iset), &
                                           lb(:, jset), npgfb(jset), nshellb(jset), zetb(:, jset), sconb_shg(:, :, jset), &
                                           lcb(:, kset), npgfcb(kset), nshellcb(kset), zetcb(:, kset), sconcb_shg(:, :, kset), &
                                           swork, swork_cont, calculate_forces)

               CALL construct_overlap_shg_abb(la(:, iset), first_sgfa(:, iset), nshella(iset), &
                                              lb(:, jset), first_sgfb(:, jset), nshellb(jset), &
                                              lcb(:, kset), first_sgfcb(:, kset), nshellcb(kset), &
                                              cg_coeff, cg_none0_list, &
                                              ncg_none0, swork_cont, Waux_mat, abbint)
               IF (calculate_forces) THEN
                  !*** -rab, since Eq. in Ref. use Ra-Rb, not Rb-Ra
                  CALL dev_overlap_shg_abb(la(:, iset), first_sgfa(:, iset), nshella(iset), &
                                           lb(:, jset), first_sgfb(:, jset), nshellb(jset), &
                                           lcb(:, kset), first_sgfcb(:, kset), nshellcb(kset), &
                                           cg_coeff, cg_none0_list, ncg_none0, -rab, swork_cont, &
                                           Waux_mat, dWaux_mat, dabbint)
               ENDIF

            END DO
         END DO
      END DO

      DEALLOCATE (swork_cont)
      DEALLOCATE (swork)

      CALL timestop(handle)

   END SUBROUTINE int_overlap_abb_shg

! **************************************************************************************************
!> \brief precalculates the angular part of the SHG integrals for the matrices
!>        (fa,fb), (a,b), (a,b,fa) and (b,fb,a); the same Waux_mat can then be used for all
!>        for integrals; specific for LRIGPW
!> \param oba orbital basis on a
!> \param obb orbital basis on b
!> \param fba aux basis on a
!> \param fbb aux basis on b
!> \param rab distance vector between a and b
!> \param Waux_mat W matrix that contains the angular-dependent part
!> \param dWaux_mat derivative of the W matrix
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE lri_precalc_angular_shg_part(oba, obb, fba, fbb, rab, Waux_mat, dWaux_mat, calculate_forces)

      TYPE(gto_basis_set_type), POINTER                  :: oba, obb, fba, fbb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: Waux_mat
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :), INTENT(OUT)              :: dWaux_mat
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'lri_precalc_angular_shg_part', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, isize, j, k, la_max, laa_max, &
                                                            lb_max, lbb_max, lca_max, lcb_max, &
                                                            li_max, lj_max, lmax, mdim(3), &
                                                            size_int(4, 2), temp
      INTEGER, DIMENSION(:), POINTER                     :: li_max_all
      REAL(KIND=dp)                                      :: rab2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Rc, Rs

      CALL timeset(routineN, handle)

      rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)

      !*** 1 Waux_mat of size (li_max,lj_max) for elements
      !                    i        j
      !    [aab]    --> (laa_max, lb_max)
      !    [bba]    --> (lbb_max, la_max) --> use for [abb]
      !    [ab] ri  --> (lca_max, lcb_max)
      !    [ab] orb --> (la_max , lb_max)

      la_max = MAXVAL(oba%lmax)
      lb_max = MAXVAL(obb%lmax)
      lca_max = MAXVAL(fba%lmax)
      lcb_max = MAXVAL(fbb%lmax)

      laa_max = la_max+lca_max
      lbb_max = lb_max+lcb_max
      li_max = MAX(laa_max, lbb_max)
      lj_max = MAX(la_max, lb_max, lcb_max)
      lmax = li_max

      ALLOCATE (li_max_all(0:lj_max))
      ALLOCATE (Rc(0:lmax, -2*lmax:2*lmax), Rs(0:lmax, -2*lmax:2*lmax))
      Rc = 0._dp
      Rs = 0._dp
      mdim(1) = li_max+lj_max+1
      mdim(2) = nsoset(li_max)+1
      mdim(3) = nsoset(lj_max)+1
      ALLOCATE (Waux_mat(mdim(1), mdim(2), mdim(3)))
      ALLOCATE (dWaux_mat(3, mdim(1), mdim(2), mdim(3)))
      !Waux_mat = 0._dp !.. takes time
      !dWaux_mat =0._dp !.. takes time

      !*** Waux_mat (li_max,lj_max) contains elements not needed,
      !*** make indixing so that only required ones are computed
      !*** li_max_all(j) --> li_max dependent on j
      size_int(1, :) = (/laa_max, lb_max/)
      size_int(2, :) = (/lbb_max, la_max/)
      size_int(3, :) = (/lca_max, lcb_max/)
      size_int(4, :) = (/la_max, lb_max/)

      li_max_all(:) = 0
      DO isize = 1, 4
         i = size_int(isize, 1)
         j = size_int(isize, 2)
         k = li_max_all(j)
         IF (k < i) li_max_all(j) = i
      ENDDO
      temp = li_max_all(lj_max)
      DO j = lj_max, 0, -1
         IF (li_max_all(j) < temp) THEN
            li_max_all(j) = temp
         ELSE
            temp = li_max_all(j)
         ENDIF
      ENDDO

      !*** -rab, since Eq. in Ref. use Ra-Rb, not Rb-Ra
      CALL get_real_scaled_solid_harmonic(Rc, Rs, lmax, -rab, rab2)
      CALL get_W_matrix(li_max_all, lj_max, lmax, Rc, Rs, Waux_mat)
      IF (calculate_forces) THEN
         CALL get_dW_matrix(li_max_all, lj_max, Waux_mat, dWaux_mat)
      ENDIF

      DEALLOCATE (Rc, Rs, li_max_all)

      CALL timestop(handle)

   END SUBROUTINE lri_precalc_angular_shg_part

END MODULE generic_shg_integrals
