!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of contracted, spherical Gaussian integrals using the (OS) integral
!>        scheme. Routines for the following two-center integrals:
!>        i)  (a|O(r12)|b) where O(r12) is the overlap, coulomb operator etc.
!>        ii) (aba) and (abb) s-overlaps
!> \par History
!>      created [06.2015]
!> \author Dorothea Golze
! **************************************************************************************************
MODULE generic_os_integrals
   USE ai_contraction_sphi,             ONLY: ab_contract,&
                                              abc_contract,&
                                              abcd_contract
   USE ai_derivatives,                  ONLY: dabdr_noscreen
   USE ai_operator_ra2m,                ONLY: operator_ra2m
   USE ai_operators_r12,                ONLY: ab_sint_os,&
                                              cps_coulomb2,&
                                              cps_gauss2,&
                                              cps_verf2,&
                                              cps_verfc2,&
                                              cps_vgauss2,&
                                              operator2
   USE ai_overlap,                      ONLY: overlap
   USE ai_overlap3,                     ONLY: overlap3
   USE ai_overlap_aabb,                 ONLY: overlap_aabb
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE debug_os_integrals,              ONLY: overlap_aabb_test,&
                                              overlap_ab_test,&
                                              overlap_abc_test
   USE input_constants,                 ONLY: operator_coulomb,&
                                              operator_gauss,&
                                              operator_verf,&
                                              operator_verfc,&
                                              operator_vgauss
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: ncoset
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'generic_os_integrals'

   PUBLIC :: int_operators_r12_ab_os, int_overlap_ab_os, int_ra2m_ab_os, int_overlap_aba_os, &
             int_overlap_abb_os, int_overlap_aabb_os

CONTAINS

! **************************************************************************************************
!> \brief Calcululates the two-center integrals of the type (a|O(r12)|b) using the OS scheme
!> \param r12_operator the integral operator, which depends on r12=|r1-r2|
!> \param vab integral matrix of spherical contracted Gaussian functions
!> \param dvab derivative of the integrals
!> \param rab distance vector between center A and B
!> \param fba basis at center A
!> \param fbb basis at center B
!> \param omega parameter in the operator
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE int_operators_r12_ab_os(r12_operator, vab, dvab, rab, fba, fbb, omega, calculate_forces)

      INTEGER, INTENT(IN)                                :: r12_operator
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: vab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dvab
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: fba, fbb
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: omega
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'int_operators_r12_ab_os', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: my_omega

      PROCEDURE(ab_sint_os), POINTER                     :: cps_operator2

      NULLIFY (cps_operator2)
      CALL timeset(routineN, handle)

      my_omega = 1.0_dp

      SELECT CASE (r12_operator)
      CASE (operator_coulomb)
         cps_operator2 => cps_coulomb2
      CASE (operator_verf)
         cps_operator2 => cps_verf2
         IF (PRESENT(omega)) my_omega = omega
      CASE (operator_verfc)
         cps_operator2 => cps_verfc2
         IF (PRESENT(omega)) my_omega = omega
      CASE (operator_vgauss)
         cps_operator2 => cps_vgauss2
         IF (PRESENT(omega)) my_omega = omega
      CASE (operator_gauss)
         cps_operator2 => cps_gauss2
         IF (PRESENT(omega)) my_omega = omega
      CASE DEFAULT
         CPABORT("Operator not available")
      END SELECT

      CALL int_operator_ab_os_low(cps_operator2, vab, dvab, rab, fba, fbb, my_omega, calculate_forces)

      CALL timestop(handle)

   END SUBROUTINE int_operators_r12_ab_os

! **************************************************************************************************
!> \brief calculate integrals (a|O(r12)|b)
!> \param cps_operator2 procedure pointer for the respective operator.
!> \param vab integral matrix of spherical contracted Gaussian functions
!> \param dvab derivative of the integrals
!> \param rab distance vector between center A and B
!> \param fba basis at center A
!> \param fbb basis at center B
!> \param omega parameter in the operator
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE int_operator_ab_os_low(cps_operator2, vab, dvab, rab, fba, fbb, omega, calculate_forces)

      PROCEDURE(ab_sint_os), POINTER                     :: cps_operator2
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: vab
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         INTENT(INOUT)                                   :: dvab
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: fba, fbb
      REAL(KIND=dp), INTENT(IN)                          :: omega
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'int_operator_ab_os_low', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iset, jset, lds, m1, m2, &
                                                            maxco, maxcoa, maxcob, maxl, maxla, &
                                                            maxlb, ncoa, ncoap, ncob, ncobp, &
                                                            nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      REAL(KIND=dp)                                      :: dab, rab2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: vac, vac_plus
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: devab, vwork
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb

      CALL timeset(routineN, handle)
      NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb, &
               first_sgfa, first_sgfb, set_radius_a, set_radius_b, rpgfa, rpgfb, &
               sphi_a, sphi_b, zeta, zetb)

      ! basis ikind
      first_sgfa => fba%first_sgf
      la_max => fba%lmax
      la_min => fba%lmin
      npgfa => fba%npgf
      nseta = fba%nset
      nsgfa => fba%nsgf_set
      rpgfa => fba%pgf_radius
      set_radius_a => fba%set_radius
      sphi_a => fba%sphi
      zeta => fba%zet
      ! basis jkind
      first_sgfb => fbb%first_sgf
      lb_max => fbb%lmax
      lb_min => fbb%lmin
      npgfb => fbb%npgf
      nsetb = fbb%nset
      nsgfb => fbb%nsgf_set
      rpgfb => fbb%pgf_radius
      set_radius_b => fbb%set_radius
      sphi_b => fbb%sphi
      zetb => fbb%zet

      CALL get_gto_basis_set(fba, maxco=maxcoa, maxl=maxla)
      CALL get_gto_basis_set(fbb, maxco=maxcob, maxl=maxlb)
      maxco = MAX(maxcoa, maxcob)
      IF (calculate_forces) THEN
         maxl = MAX(maxla+1, maxlb)
      ELSE
         maxl = MAX(maxla, maxlb)
      ENDIF
      lds = ncoset(maxl)

      rab2 = SUM(rab*rab)
      dab = SQRT(rab2)

      DO iset = 1, nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         ncoap = npgfa(iset)*ncoset(la_max(iset)+1)
         sgfa = first_sgfa(1, iset)

         DO jset = 1, nsetb

            ncob = npgfb(jset)*ncoset(lb_max(jset))
            ncobp = npgfb(jset)*ncoset(lb_max(jset)+1)
            sgfb = first_sgfb(1, jset)
            m1 = sgfa+nsgfa(iset)-1
            m2 = sgfb+nsgfb(jset)-1

            ! calculate integrals
            IF (calculate_forces) THEN
               ALLOCATE (vwork(ncoap, ncobp, la_max(iset)+lb_max(jset)+3), &
                         vac(ncoa, ncob), vac_plus(ncoap, ncobp), devab(ncoa, ncob, 3))
               devab = 0._dp
               vwork = 0.0_dp
               vac = 0.0_dp
               CALL operator2(cps_operator2, la_max(iset)+1, npgfa(iset), zeta(:, iset), la_min(iset), &
                              lb_max(jset)+1, npgfb(jset), zetb(:, jset), lb_min(jset), &
                              omega, rab, rab2, vac, vwork, maxder=1, vac_plus=vac_plus)
               CALL dabdr_noscreen(la_max(iset), npgfa(iset), zeta(:, iset), lb_max(jset), npgfb(jset), &
                                   vac_plus, devab(:, :, 1), devab(:, :, 2), devab(:, :, 3))
               DO i = 1, 3
                  CALL ab_contract(dvab(sgfa:m1, sgfb:m2, i), devab(:, :, i), sphi_a(:, sgfa:), &
                                   sphi_b(:, sgfb:), ncoa, ncob, nsgfa(iset), nsgfb(jset))
               ENDDO

            ELSE
               ALLOCATE (vwork(ncoa, ncob, la_max(iset)+lb_max(jset)+1), &
                         vac(ncoa, ncob), vac_plus(ncoap, ncobp), devab(ncoa, ncob, 3))
               vwork = 0.0_dp
               vac = 0.0_dp
               CALL operator2(cps_operator2, la_max(iset), npgfa(iset), zeta(:, iset), la_min(iset), &
                              lb_max(jset), npgfb(jset), zetb(:, jset), lb_min(jset), &
                              omega, rab, rab2, vac, vwork)
            ENDIF

            CALL ab_contract(vab(sgfa:m1, sgfb:m2), vac(1:ncoa, 1:ncob), sphi_a(:, sgfa:), sphi_b(:, sgfb:), &
                             ncoa, ncob, nsgfa(iset), nsgfb(jset))
            DEALLOCATE (vwork, vac, vac_plus, devab)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE int_operator_ab_os_low

! **************************************************************************************************
!> \brief calculate overlap integrals (a,b)
!> \param sab integral (a,b)
!> \param dsab derivative of sab with respect to A
!> \param ra position center A, for debug only
!> \param rb position center B, for debug only
!> \param rab distance vector between center A and B
!> \param fba basis at center A
!> \param fbb basis at center B
!> \param calculate_forces ...
!> \param debug integrals are debugged by recursive routines if requested
!> \param dmax maximal deviation between integrals when debugging
! **************************************************************************************************
   SUBROUTINE int_overlap_ab_os(sab, dsab, ra, rb, rab, fba, fbb, calculate_forces, debug, &
                                dmax)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: sab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dsab
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rb, rab
      TYPE(gto_basis_set_type), POINTER                  :: fba, fbb
      LOGICAL, INTENT(IN)                                :: calculate_forces, debug
      REAL(KIND=dp), INTENT(INOUT)                       :: dmax

      CHARACTER(LEN=*), PARAMETER :: routineN = 'int_overlap_ab_os', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iset, jset, lds, m1, m2, &
                                                            maxco, maxcoa, maxcob, maxl, maxla, &
                                                            maxlb, ncoa, ncob, nseta, nsetb, sgfa, &
                                                            sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: sint
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: devab, swork
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb

      CALL timeset(routineN, handle)
      NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb, &
               first_sgfa, first_sgfb, set_radius_a, set_radius_b, rpgfa, rpgfb, &
               sphi_a, sphi_b, zeta, zetb)

      ! basis ikind
      first_sgfa => fba%first_sgf
      la_max => fba%lmax
      la_min => fba%lmin
      npgfa => fba%npgf
      nseta = fba%nset
      nsgfa => fba%nsgf_set
      rpgfa => fba%pgf_radius
      set_radius_a => fba%set_radius
      sphi_a => fba%sphi
      zeta => fba%zet
      ! basis jkind
      first_sgfb => fbb%first_sgf
      lb_max => fbb%lmax
      lb_min => fbb%lmin
      npgfb => fbb%npgf
      nsetb = fbb%nset
      nsgfb => fbb%nsgf_set
      rpgfb => fbb%pgf_radius
      set_radius_b => fbb%set_radius
      sphi_b => fbb%sphi
      zetb => fbb%zet

      CALL get_gto_basis_set(fba, maxco=maxcoa, maxl=maxla)
      CALL get_gto_basis_set(fbb, maxco=maxcob, maxl=maxlb)
      maxco = MAX(maxcoa, maxcob)
      IF (calculate_forces) THEN
         maxl = MAX(maxla+1, maxlb)
      ELSE
         maxl = MAX(maxla, maxlb)
      ENDIF
      lds = ncoset(maxl)
      ALLOCATE (sint(maxco, maxco))
      IF (calculate_forces) THEN
         !derivatives will be stored in devab(:,:,2:4)
         ALLOCATE (swork(lds, lds, 4), devab(maxco, maxco, 4))
         devab = 0._dp
      ELSE
         ALLOCATE (swork(lds, lds, 1))
      ENDIF
      sint = 0._dp
      swork = 0._dp

      dab = SQRT(SUM(rab**2))

      sab = 0.0_dp
      dsab = 0.0_dp

      DO iset = 1, nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1, iset)

         DO jset = 1, nsetb

            IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

            ncob = npgfb(jset)*ncoset(lb_max(jset))
            sgfb = first_sgfb(1, jset)
            m1 = sgfa+nsgfa(iset)-1
            m2 = sgfb+nsgfb(jset)-1
            sint = 0._dp
            swork = 0._dp

            ! calculate integrals
            IF (calculate_forces) THEN
               devab = 0._dp
               CALL overlap(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                            lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                            rab, dab, sint, 0, .FALSE., swork, lds, sdab=devab)
               DO i = 1, 3
                  !NOTE: devab(:,:,2:4) contains all derivatives for lmin=0 to lmax=lmax
                  !      correct after contraction (multiply with zero for elements l < lmin)
                  CALL ab_contract(dsab(sgfa:m1, sgfb:m2, i), devab(:, :, i+1), sphi_a(:, sgfa:), &
                                   sphi_b(:, sgfb:), ncoa, ncob, nsgfa(iset), nsgfb(jset))
               ENDDO

            ELSE
               CALL overlap(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                            lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                            rab, dab, sint, 0, .FALSE., swork, lds)
            ENDIF
            ! debug if requested
            IF (debug) THEN
               CALL overlap_ab_test(la_max(iset), la_min(iset), npgfa(iset), zeta(:, iset), &
                                    lb_max(jset), lb_min(jset), npgfb(jset), zetb(:, jset), &
                                    ra, rb, sint, dmax)
            ENDIF

            CALL ab_contract(sab(sgfa:m1, sgfb:m2), sint, sphi_a(:, sgfa:), sphi_b(:, sgfb:), &
                             ncoa, ncob, nsgfa(iset), nsgfb(jset))
         END DO
      END DO

      IF (calculate_forces) THEN
         DEALLOCATE (devab)
      ENDIF

      DEALLOCATE (sint, swork)

      CALL timestop(handle)

   END SUBROUTINE int_overlap_ab_os

! **************************************************************************************************
!> \brief calculate integrals (a|(r-Ra)^(2m)|b)
!> \param sab integral (a|(r-Ra)^(2m)|b)
!> \param dsab derivative of sab with respect to A
!> \param rab distance vector between center A and B
!> \param fba fba basis at center A
!> \param fbb fbb basis at center B
!> \param m exponent m in operator (r-Ra)^(2m)
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE int_ra2m_ab_os(sab, dsab, rab, fba, fbb, m, calculate_forces)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: sab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dsab
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: fba, fbb
      INTEGER, INTENT(IN)                                :: m
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'int_ra2m_ab_os', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iset, jset, m1, m2, maxco, &
                                                            maxcoa, maxcob, maxl, maxla, maxlb, &
                                                            ncoa, ncob, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: failure
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: sint
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dsint
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: scon_a, scon_b, zeta, zetb

      failure = .FALSE.
      CALL timeset(routineN, handle)
      NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb, &
               first_sgfa, first_sgfb, scon_a, scon_b, zeta, zetb)

      ! basis ikind
      first_sgfa => fba%first_sgf
      la_max => fba%lmax
      la_min => fba%lmin
      npgfa => fba%npgf
      nseta = fba%nset
      nsgfa => fba%nsgf_set
      scon_a => fba%scon
      zeta => fba%zet
      ! basis jkind
      first_sgfb => fbb%first_sgf
      lb_max => fbb%lmax
      lb_min => fbb%lmin
      npgfb => fbb%npgf
      nsetb = fbb%nset
      nsgfb => fbb%nsgf_set
      scon_b => fbb%scon
      zetb => fbb%zet

      CALL get_gto_basis_set(fba, maxco=maxcoa, maxl=maxla)
      CALL get_gto_basis_set(fbb, maxco=maxcob, maxl=maxlb)
      maxco = MAX(maxcoa, maxcob)
      maxl = MAX(maxla, maxlb)
      ALLOCATE (sint(maxco, maxco))
      ALLOCATE (dsint(maxco, maxco, 3))

      dab = SQRT(SUM(rab**2))

      sab = 0.0_dp
      dsab = 0.0_dp

      DO iset = 1, nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1, iset)

         DO jset = 1, nsetb

            ncob = npgfb(jset)*ncoset(lb_max(jset))
            sgfb = first_sgfb(1, jset)
            m1 = sgfa+nsgfa(iset)-1
            m2 = sgfb+nsgfb(jset)-1
            sint(1:ncoa, 1:ncob) = 0._dp
            dsint(1:ncoa, 1:ncob, :) = 0._dp

            CALL operator_ra2m(la_max(iset), la_min(iset), npgfa(iset), zeta(:, iset), &
                               lb_max(jset), lb_min(jset), npgfb(jset), zetb(:, jset), &
                               m, rab, sint, dsint, calculate_forces)

            CALL ab_contract(sab(sgfa:m1, sgfb:m2), sint, scon_a(:, sgfa:), scon_b(:, sgfb:), &
                             ncoa, ncob, nsgfa(iset), nsgfb(jset))
            IF (calculate_forces) THEN
               DO i = 1, 3
                  CALL ab_contract(dsab(sgfa:m1, sgfb:m2, i), dsint(:, :, i), scon_a(:, sgfa:), &
                                   scon_b(:, sgfb:), ncoa, ncob, nsgfa(iset), nsgfb(jset))
               ENDDO
            ENDIF
         END DO
      END DO

      DEALLOCATE (sint, dsint)

      CALL timestop(handle)

   END SUBROUTINE int_ra2m_ab_os

! **************************************************************************************************
!> \brief calculate integrals (a,b,fa)
!> \param abaint integral (a,b,fa)
!> \param dabdaint derivative of abaint with respect to A
!> \param ra position center A, for debug only
!> \param rb position center B, for debug only
!> \param rab distance vector between center A and B
!> \param oba orbital basis at center A
!> \param obb orbital basis at center B
!> \param fba auxiliary basis set at center A
!> \param calculate_forces ...
!> \param debug integrals are debugged by recursive routines if requested
!> \param dmax maximal deviation between integrals when debugging
! **************************************************************************************************
   SUBROUTINE int_overlap_aba_os(abaint, dabdaint, ra, rb, rab, oba, obb, fba, &
                                 calculate_forces, debug, dmax)

      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: abaint
      REAL(KIND=dp), DIMENSION(:, :, :, :), OPTIONAL, &
         POINTER                                         :: dabdaint
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rb, rab
      TYPE(gto_basis_set_type), POINTER                  :: oba, obb, fba
      LOGICAL, INTENT(IN)                                :: calculate_forces, debug
      REAL(KIND=dp), INTENT(INOUT)                       :: dmax

      CHARACTER(LEN=*), PARAMETER :: routineN = 'int_overlap_aba_os', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iset, jset, kaset, m1, m2, &
                                                            m3, ncoa, ncob, ncoc, nseta, nsetb, &
                                                            nsetca, sgfa, sgfb, sgfc
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, lca_max, &
                                                            lca_min, npgfa, npgfb, npgfca, nsgfa, &
                                                            nsgfb, nsgfca
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb, first_sgfca
      REAL(KIND=dp)                                      :: dab, dac, dbc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: saba
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: sabda, sdaba, sdabda
      REAL(KIND=dp), DIMENSION(3)                        :: rac, rbc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b, set_radius_ca
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, rpgfca, sphi_a, sphi_b, &
                                                            sphi_ca, zeta, zetb, zetca

      CALL timeset(routineN, handle)
      NULLIFY (la_max, la_min, lb_max, lb_min, lca_max, lca_min, npgfa, npgfb, &
               npgfca, nsgfa, nsgfb, nsgfca)
      NULLIFY (first_sgfa, first_sgfb, first_sgfca, set_radius_a, set_radius_b, &
               set_radius_ca, rpgfa, rpgfb, rpgfca, sphi_a, sphi_b, sphi_ca, &
               zeta, zetb, zetca)

      ! basis ikind
      first_sgfa => oba%first_sgf
      la_max => oba%lmax
      la_min => oba%lmin
      npgfa => oba%npgf
      nseta = oba%nset
      nsgfa => oba%nsgf_set
      rpgfa => oba%pgf_radius
      set_radius_a => oba%set_radius
      sphi_a => oba%sphi
      zeta => oba%zet
      ! basis jkind
      first_sgfb => obb%first_sgf
      lb_max => obb%lmax
      lb_min => obb%lmin
      npgfb => obb%npgf
      nsetb = obb%nset
      nsgfb => obb%nsgf_set
      rpgfb => obb%pgf_radius
      set_radius_b => obb%set_radius
      sphi_b => obb%sphi
      zetb => obb%zet

      ! basis RI A
      first_sgfca => fba%first_sgf
      lca_max => fba%lmax
      lca_min => fba%lmin
      npgfca => fba%npgf
      nsetca = fba%nset
      nsgfca => fba%nsgf_set
      rpgfca => fba%pgf_radius
      set_radius_ca => fba%set_radius
      sphi_ca => fba%sphi
      zetca => fba%zet

      dab = SQRT(SUM(rab**2))

      DO iset = 1, nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1, iset)

         DO jset = 1, nsetb

            IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

            ncob = npgfb(jset)*ncoset(lb_max(jset))
            sgfb = first_sgfb(1, jset)
            m1 = sgfa+nsgfa(iset)-1
            m2 = sgfb+nsgfb(jset)-1

            ! calculate integrals abaint and derivative [d(a,b,a)/dA] dabdaint if requested
            rac = 0._dp
            dac = 0._dp
            rbc = -rab
            dbc = dab

            DO kaset = 1, nsetca

               IF (set_radius_b(jset)+set_radius_ca(kaset) < dab) CYCLE

               ncoc = npgfca(kaset)*ncoset(lca_max(kaset))
               sgfc = first_sgfca(1, kaset)
               m3 = sgfc+nsgfca(kaset)-1
               IF (ncoa*ncob*ncoc > 0) THEN
                  ALLOCATE (saba(ncoa, ncob, ncoc))
                  saba(:, :, :) = 0._dp
                  ! integrals
                  IF (calculate_forces) THEN
                     ALLOCATE (sdaba(ncoa, ncob, ncoc, 3))
                     ALLOCATE (sabda(ncoa, ncob, ncoc, 3))
                     ALLOCATE (sdabda(ncoa, ncob, ncoc, 3))
                     sdaba(:, :, :, :) = 0._dp
                     sabda(:, :, :, :) = 0._dp
                     sdabda(:, :, :, :) = 0._dp
                     CALL overlap3(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                   lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                   lca_max(kaset), npgfca(kaset), zetca(:, kaset), rpgfca(:, kaset), lca_min(kaset), &
                                   rab, dab, rac, dac, rbc, dbc, saba, sdaba, sabda)
                     !d(a,b,a)/dA = (da/dA,b,a) + (a,b,da/dA)
                     sdabda(:, :, :, :) = sdaba+sabda

                     DO i = 1, 3
                        CALL abc_contract(dabdaint(sgfa:m1, sgfb:m2, sgfc:m3, i), sdabda(:, :, :, i), &
                                          sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_ca(:, sgfc:), &
                                          ncoa, ncob, ncoc, nsgfa(iset), nsgfb(jset), nsgfca(kaset))
                     ENDDO

                     DEALLOCATE (sdaba, sabda, sdabda)
                  ELSE
                     CALL overlap3(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                   lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                   lca_max(kaset), npgfca(kaset), zetca(:, kaset), rpgfca(:, kaset), lca_min(kaset), &
                                   rab, dab, rac, dac, rbc, dbc, saba)
                  ENDIF
                  ! debug if requested
                  IF (debug) THEN
                     CALL overlap_abc_test(la_max(iset), npgfa(iset), zeta(:, iset), la_min(iset), &
                                           lb_max(jset), npgfb(jset), zetb(:, jset), lb_min(jset), &
                                           lca_max(kaset), npgfca(kaset), zetca(:, kaset), lca_min(kaset), &
                                           ra, rb, ra, saba, dmax)
                  ENDIF
                  CALL abc_contract(abaint(sgfa:m1, sgfb:m2, sgfc:m3), saba, &
                                    sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_ca(:, sgfc:), &
                                    ncoa, ncob, ncoc, nsgfa(iset), nsgfb(jset), nsgfca(kaset))
                  DEALLOCATE (saba)
               END IF
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE int_overlap_aba_os

! **************************************************************************************************
!> \brief calculate integrals (a,b,fb)
!> \param abbint integral (a,b,fb)
!> \param dabbint derivative of abbint with respect to A
!> \param ra position center A, for debug only
!> \param rb position center B, for debug only
!> \param rab distance vector between center A and B
!> \param oba orbital basis at center A
!> \param obb orbital basis at center B
!> \param fbb auxiliary basis set at center B
!> \param calculate_forces ...
!> \param debug integrals are debugged by recursive routines if requested
!> \param dmax maximal deviation between integrals when debugging
! **************************************************************************************************
   SUBROUTINE int_overlap_abb_os(abbint, dabbint, ra, rb, rab, oba, obb, fbb, calculate_forces, &
                                 debug, dmax)

      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: abbint
      REAL(KIND=dp), DIMENSION(:, :, :, :), OPTIONAL, &
         POINTER                                         :: dabbint
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rb, rab
      TYPE(gto_basis_set_type), POINTER                  :: oba, obb, fbb
      LOGICAL, INTENT(IN)                                :: calculate_forces, debug
      REAL(KIND=dp), INTENT(INOUT)                       :: dmax

      CHARACTER(LEN=*), PARAMETER :: routineN = 'int_overlap_abb_os', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iset, jset, kbset, m1, m2, &
                                                            m3, ncoa, ncob, ncoc, nseta, nsetb, &
                                                            nsetcb, sgfa, sgfb, sgfc
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, lcb_max, &
                                                            lcb_min, npgfa, npgfb, npgfcb, nsgfa, &
                                                            nsgfb, nsgfcb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb, first_sgfcb
      REAL(KIND=dp)                                      :: dab, dac, dbc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: sabb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: sdabb
      REAL(KIND=dp), DIMENSION(3)                        :: rac, rbc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b, set_radius_cb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, rpgfcb, sphi_a, sphi_b, &
                                                            sphi_cb, zeta, zetb, zetcb

      CALL timeset(routineN, handle)
      NULLIFY (la_max, la_min, lb_max, lb_min, lcb_max, lcb_min, npgfa, npgfb, &
               npgfcb, nsgfa, nsgfb, nsgfcb)
      NULLIFY (first_sgfa, first_sgfb, first_sgfcb, set_radius_a, set_radius_b, &
               set_radius_cb, rpgfa, rpgfb, rpgfcb, sphi_a, sphi_b, sphi_cb, &
               zeta, zetb, zetcb)

      ! basis ikind
      first_sgfa => oba%first_sgf
      la_max => oba%lmax
      la_min => oba%lmin
      npgfa => oba%npgf
      nseta = oba%nset
      nsgfa => oba%nsgf_set
      rpgfa => oba%pgf_radius
      set_radius_a => oba%set_radius
      sphi_a => oba%sphi
      zeta => oba%zet
      ! basis jkind
      first_sgfb => obb%first_sgf
      lb_max => obb%lmax
      lb_min => obb%lmin
      npgfb => obb%npgf
      nsetb = obb%nset
      nsgfb => obb%nsgf_set
      rpgfb => obb%pgf_radius
      set_radius_b => obb%set_radius
      sphi_b => obb%sphi
      zetb => obb%zet

      ! basis RI B
      first_sgfcb => fbb%first_sgf
      lcb_max => fbb%lmax
      lcb_min => fbb%lmin
      npgfcb => fbb%npgf
      nsetcb = fbb%nset
      nsgfcb => fbb%nsgf_set
      rpgfcb => fbb%pgf_radius
      set_radius_cb => fbb%set_radius
      sphi_cb => fbb%sphi
      zetcb => fbb%zet

      dab = SQRT(SUM(rab**2))

      DO iset = 1, nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1, iset)

         DO jset = 1, nsetb

            IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

            ncob = npgfb(jset)*ncoset(lb_max(jset))
            sgfb = first_sgfb(1, jset)
            m1 = sgfa+nsgfa(iset)-1
            m2 = sgfb+nsgfb(jset)-1

            ! calculate integrals abbint and derivative [d(a,b,b)/dA] dabbint if requested
            rac = rab
            dac = dab
            rbc = 0._dp
            dbc = 0._dp

            DO kbset = 1, nsetcb

               IF (set_radius_a(iset)+set_radius_cb(kbset) < dab) CYCLE

               ncoc = npgfcb(kbset)*ncoset(lcb_max(kbset))
               sgfc = first_sgfcb(1, kbset)
               m3 = sgfc+nsgfcb(kbset)-1
               IF (ncoa*ncob*ncoc > 0) THEN
                  ALLOCATE (sabb(ncoa, ncob, ncoc))
                  sabb(:, :, :) = 0._dp
                  IF (calculate_forces) THEN
                     ALLOCATE (sdabb(ncoa, ncob, ncoc, 3))
                     sdabb(:, :, :, :) = 0._dp
                     CALL overlap3(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                   lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                   lcb_max(kbset), npgfcb(kbset), zetcb(:, kbset), rpgfcb(:, kbset), lcb_min(kbset), &
                                   rab, dab, rac, dac, rbc, dbc, sabb, sdabc=sdabb)
                     DO i = 1, 3
                        CALL abc_contract(dabbint(sgfa:m1, sgfb:m2, sgfc:m3, i), sdabb(:, :, :, i), &
                                          sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_cb(:, sgfc:), &
                                          ncoa, ncob, ncoc, nsgfa(iset), nsgfb(jset), nsgfcb(kbset))
                     ENDDO
                     DEALLOCATE (sdabb)
                  ELSE
                     CALL overlap3(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                   lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                   lcb_max(kbset), npgfcb(kbset), zetcb(:, kbset), rpgfcb(:, kbset), lcb_min(kbset), &
                                   rab, dab, rac, dac, rbc, dbc, sabb)
                  ENDIF
                  ! debug if requested
                  IF (debug) THEN
                     CALL overlap_abc_test(la_max(iset), npgfa(iset), zeta(:, iset), la_min(iset), &
                                           lb_max(jset), npgfb(jset), zetb(:, jset), lb_min(jset), &
                                           lcb_max(kbset), npgfcb(kbset), zetcb(:, kbset), lcb_min(kbset), &
                                           ra, rb, rb, sabb, dmax)
                  ENDIF
                  CALL abc_contract(abbint(sgfa:m1, sgfb:m2, sgfc:m3), sabb, &
                                    sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_cb(:, sgfc:), &
                                    ncoa, ncob, ncoc, nsgfa(iset), nsgfb(jset), nsgfcb(kbset))
                  DEALLOCATE (sabb)
               ENDIF
            END DO

         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE int_overlap_abb_os

! **************************************************************************************************
!> \brief calculate overlap integrals (aa,bb)
!> \param saabb integral (aa,bb)
!> \param oba orbital basis at center A
!> \param obb orbital basis at center B
!> \param rab ...
!> \param ra ...
!> \param rb ...
!> \param debug integrals are debugged by recursive routines if requested
!> \param dmax maximal deviation between integrals when debugging
! **************************************************************************************************
   SUBROUTINE int_overlap_aabb_os(saabb, oba, obb, rab, ra, rb, debug, dmax)

      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: saabb
      TYPE(gto_basis_set_type), POINTER                  :: oba, obb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab, ra, rb
      LOGICAL, INTENT(IN)                                :: debug
      REAL(KIND=dp), INTENT(INOUT)                       :: dmax

      CHARACTER(LEN=*), PARAMETER :: routineN = 'int_overlap_aabb_os', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, iset, isgfa1, jset, jsgfa2, kset, ksgfb1, lds, lset, lsgfb2, m1, m2, m3, &
         m4, maxco, maxcoa, maxcob, maxl, maxla, maxlb, ncoa1, ncoa2, ncob1, ncob2, nseta, nsetb, &
         sgfa1, sgfa2, sgfb1, sgfb2
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: asets_equal, bsets_equal
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: swork
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: sint
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb

      CALL timeset(routineN, handle)
      NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb, &
               first_sgfa, first_sgfb, set_radius_a, set_radius_b, rpgfa, rpgfb, &
               sphi_a, sphi_b, zeta, zetb)

      ! basis ikind
      first_sgfa => oba%first_sgf
      la_max => oba%lmax
      la_min => oba%lmin
      npgfa => oba%npgf
      nseta = oba%nset
      nsgfa => oba%nsgf_set
      rpgfa => oba%pgf_radius
      set_radius_a => oba%set_radius
      sphi_a => oba%sphi
      zeta => oba%zet
      ! basis jkind
      first_sgfb => obb%first_sgf
      lb_max => obb%lmax
      lb_min => obb%lmin
      npgfb => obb%npgf
      nsetb = obb%nset
      nsgfb => obb%nsgf_set
      rpgfb => obb%pgf_radius
      set_radius_b => obb%set_radius
      sphi_b => obb%sphi
      zetb => obb%zet

      CALL get_gto_basis_set(oba, maxco=maxcoa, maxl=maxla)
      CALL get_gto_basis_set(obb, maxco=maxcob, maxl=maxlb)
      maxco = MAX(maxcoa, maxcob)
      maxla = 2*maxla
      maxlb = 2*maxlb
      maxl = MAX(maxla, maxlb)
      lds = ncoset(maxl)
      ALLOCATE (sint(maxco, maxco, maxco, maxco))
      ALLOCATE (swork(lds, lds))
      sint = 0._dp
      swork = 0._dp

      dab = SQRT(SUM(rab**2))

      DO iset = 1, nseta

         ncoa1 = npgfa(iset)*ncoset(la_max(iset))
         sgfa1 = first_sgfa(1, iset)
         m1 = sgfa1+nsgfa(iset)-1

         DO jset = iset, nseta

            ncoa2 = npgfa(jset)*ncoset(la_max(jset))
            sgfa2 = first_sgfa(1, jset)
            m2 = sgfa2+nsgfa(jset)-1

            DO kset = 1, nsetb

               ncob1 = npgfb(kset)*ncoset(lb_max(kset))
               sgfb1 = first_sgfb(1, kset)
               m3 = sgfb1+nsgfb(kset)-1

               DO lset = kset, nsetb

                  ncob2 = npgfb(lset)*ncoset(lb_max(lset))
                  sgfb2 = first_sgfb(1, lset)
                  m4 = sgfb2+nsgfb(lset)-1

                  ! check if sets are identical to spare some integral evaluation
                  asets_equal = .FALSE.
                  IF (iset == jset) asets_equal = .TRUE.
                  bsets_equal = .FALSE.
                  IF (kset == lset) bsets_equal = .TRUE.
                  ! calculate integrals
                  CALL overlap_aabb(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                    la_max(jset), la_min(jset), npgfa(jset), rpgfa(:, jset), zeta(:, jset), &
                                    lb_max(kset), lb_min(kset), npgfb(kset), rpgfb(:, kset), zetb(:, kset), &
                                    lb_max(lset), lb_min(lset), npgfb(lset), rpgfb(:, lset), zetb(:, lset), &
                                    asets_equal, bsets_equal, rab, dab, sint, swork, lds)
                  ! debug if requested
                  IF (debug) THEN
                     CALL overlap_aabb_test(la_max(iset), la_min(iset), npgfa(iset), zeta(:, iset), &
                                            la_max(jset), la_min(jset), npgfa(jset), zeta(:, jset), &
                                            lb_max(kset), lb_min(kset), npgfb(kset), zetb(:, kset), &
                                            lb_max(lset), lb_min(lset), npgfb(lset), zetb(:, lset), &
                                            ra, rb, sint, dmax)
                  ENDIF

                  CALL abcd_contract(saabb(sgfa1:m1, sgfa2:m2, sgfb1:m3, sgfb2:m4), sint, sphi_a(:, sgfa1:), &
                                     sphi_a(:, sgfa2:), sphi_b(:, sgfb1:), sphi_b(:, sgfb2:), ncoa1, ncoa2, &
                                     ncob1, ncob2, nsgfa(iset), nsgfa(jset), nsgfb(kset), nsgfb(lset))

                  ! account for the fact that some integrals are alike
                  DO isgfa1 = sgfa1, m1
                     DO jsgfa2 = sgfa2, m2
                        DO ksgfb1 = sgfb1, m3
                           DO lsgfb2 = sgfb2, m4
                              saabb(jsgfa2, isgfa1, ksgfb1, lsgfb2) = saabb(isgfa1, jsgfa2, ksgfb1, lsgfb2)
                              saabb(isgfa1, jsgfa2, lsgfb2, ksgfb1) = saabb(isgfa1, jsgfa2, ksgfb1, lsgfb2)
                              saabb(jsgfa2, isgfa1, lsgfb2, ksgfb1) = saabb(isgfa1, jsgfa2, ksgfb1, lsgfb2)
                           END DO
                        END DO
                     END DO
                  END DO

               END DO
            END DO
         END DO
      END DO

      DEALLOCATE (sint, swork)

      CALL timestop(handle)

   END SUBROUTINE int_overlap_aabb_os

END MODULE generic_os_integrals
