!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief initializes the environment for lri
!>        lri : local resolution of the identity
!> \par History
!>      created [06.2015]
!> \author Dorothea Golze
! **************************************************************************************************
MODULE lri_environment_init
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: copy_gto_basis_set,&
                                              gto_basis_set_type
   USE generic_shg_integrals,           ONLY: int_overlap_aba_shg
   USE generic_shg_integrals_init,      ONLY: contraction_matrix_shg,&
                                              contraction_matrix_shg_mix,&
                                              get_clebsch_gordon_coefficients
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: deallocate_bas_properties,&
                                              lri_env_create,&
                                              lri_environment_type
   USE mathconstants,                   ONLY: fac,&
                                              pi
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_environment_init'

   PUBLIC :: lri_env_init, lri_basis_init

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief initializes the lri env
!> \param qs_env ...
!> \param force_env_section ...
!> \param qs_kind_set ...
! **************************************************************************************************
   SUBROUTINE lri_env_init(qs_env, force_env_section, qs_kind_set)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: force_env_section
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CHARACTER(len=*), PARAMETER :: routineN = 'lri_env_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: ikind, lmax_ikind_orb, lmax_ikind_ri, &
                                                            maxl_orb, maxl_ri, nkind
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(gto_basis_set_type), POINTER                  :: lri_basis_set, orb_basis_set
      TYPE(lri_environment_type), POINTER                :: lri_env

      NULLIFY (atomic_kind_set, lri_basis_set, lri_env, orb_basis_set)
      CALL lri_env_create(lri_env)

      ! initialize the basic basis sets (orb and ri)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (lri_env%orb_basis(nkind), lri_env%ri_basis(nkind))
      maxl_orb = 0
      maxl_ri = 0
      DO ikind = 1, nkind
         NULLIFY (orb_basis_set, lri_basis_set)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type="ORB")
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=lri_basis_set, basis_type="LRI")
         NULLIFY (lri_env%orb_basis(ikind)%gto_basis_set)
         NULLIFY (lri_env%ri_basis(ikind)%gto_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            CALL copy_gto_basis_set(orb_basis_set, lri_env%orb_basis(ikind)%gto_basis_set)
            CALL copy_gto_basis_set(lri_basis_set, lri_env%ri_basis(ikind)%gto_basis_set)
         END IF
         lmax_ikind_orb = MAXVAL(lri_env%orb_basis(ikind)%gto_basis_set%lmax)
         lmax_ikind_ri = MAXVAL(lri_env%ri_basis(ikind)%gto_basis_set%lmax)
         maxl_orb = MAX(maxl_orb, lmax_ikind_orb)
         maxl_ri = MAX(maxl_ri, lmax_ikind_ri)
      END DO

      ! CG coefficients needed for lri integrals
      IF (ASSOCIATED(lri_env%cg_shg)) THEN
         CALL get_clebsch_gordon_coefficients(lri_env%cg_shg%cg_coeff, &
                                              lri_env%cg_shg%cg_none0_list, &
                                              lri_env%cg_shg%ncg_none0, &
                                              maxl_orb, maxl_ri)
      ENDIF
      CALL lri_basis_init(lri_env)
      ! init keywords

      ! check for debug (OS scheme)
      CALL section_vals_val_get(force_env_section, "DFT%QS%LRIGPW%DEBUG_LRI_INTEGRALS", &
                                l_val=lri_env%debug)
      ! integrals based on solid harmonic Gaussians
      CALL section_vals_val_get(force_env_section, "DFT%QS%LRIGPW%SHG_LRI_INTEGRALS", &
                                l_val=lri_env%use_shg_integrals)
      ! how to calculate inverse/pseuodinverse of overlap
      CALL section_vals_val_get(force_env_section, "DFT%QS%LRIGPW%LRI_OVERLAP_MATRIX", &
                                i_val=lri_env%lri_overlap_inv)
      CALL section_vals_val_get(force_env_section, "DFT%QS%LRIGPW%MAX_CONDITION_NUM", &
                                r_val=lri_env%cond_max)
      !
      CALL set_qs_env(qs_env, lri_env=lri_env)

   END SUBROUTINE lri_env_init

! **************************************************************************************************
!> \brief initializes the lri basis: calculates the norm, self-overlap
!>        and integral of the ri basis
!> \param lri_env ...
! **************************************************************************************************
   SUBROUTINE lri_basis_init(lri_env)
      TYPE(lri_environment_type), POINTER                :: lri_env

      CHARACTER(len=*), PARAMETER :: routineN = 'lri_basis_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: ikind, nkind
      INTEGER, DIMENSION(:, :, :), POINTER               :: orb_index, ri_index
      REAL(KIND=dp), DIMENSION(:), POINTER               :: orb_norm_r, ri_int_fbas, ri_norm_r, &
                                                            ri_norm_s
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: orb_ovlp, ri_ovlp
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: ovlp3, scon_orb, scon_ri
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: dovlp3, scon_mix
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis, ri_basis

      IF (ASSOCIATED(lri_env)) THEN
         IF (ASSOCIATED(lri_env%orb_basis)) THEN
            CPASSERT(ASSOCIATED(lri_env%ri_basis))
            nkind = SIZE(lri_env%orb_basis)
            CALL deallocate_bas_properties(lri_env)
            ALLOCATE (lri_env%bas_prop(nkind))
            DO ikind = 1, nkind
               NULLIFY (orb_basis, ri_basis)
               orb_basis => lri_env%orb_basis(ikind)%gto_basis_set
               IF (ASSOCIATED(orb_basis)) THEN
                  ri_basis => lri_env%ri_basis(ikind)%gto_basis_set
                  CPASSERT(ASSOCIATED(ri_basis))
                  NULLIFY (ri_norm_r)
                  CALL basis_norm_radial(ri_basis, ri_norm_r)
                  NULLIFY (orb_norm_r)
                  CALL basis_norm_radial(orb_basis, orb_norm_r)
                  NULLIFY (ri_norm_s)
                  CALL basis_norm_s_func(ri_basis, ri_norm_s)
                  NULLIFY (ri_int_fbas)
                  CALL basis_int(ri_basis, ri_int_fbas, ri_norm_s)
                  lri_env%bas_prop(ikind)%int_fbas => ri_int_fbas
                  NULLIFY (ri_ovlp)
                  CALL basis_ovlp(ri_basis, ri_ovlp, ri_norm_r)
                  lri_env%bas_prop(ikind)%ri_ovlp => ri_ovlp
                  NULLIFY (orb_ovlp)
                  CALL basis_ovlp(orb_basis, orb_ovlp, orb_norm_r)
                  lri_env%bas_prop(ikind)%orb_ovlp => orb_ovlp
                  NULLIFY (scon_ri)
                  CALL contraction_matrix_shg(ri_basis, scon_ri)
                  lri_env%bas_prop(ikind)%scon_ri => scon_ri
                  NULLIFY (scon_orb)
                  CALL contraction_matrix_shg(orb_basis, scon_orb)
                  lri_env%bas_prop(ikind)%scon_orb => scon_orb
                  NULLIFY (scon_mix)
                  CALL contraction_matrix_shg_mix(orb_basis, ri_basis, &
                                                  orb_index, ri_index, scon_mix)
                  lri_env%bas_prop(ikind)%scon_mix => scon_mix
                  lri_env%bas_prop(ikind)%orb_index => orb_index
                  lri_env%bas_prop(ikind)%ri_index => ri_index
                  NULLIFY (ovlp3)
                  ALLOCATE (ovlp3(orb_basis%nsgf, orb_basis%nsgf, ri_basis%nsgf))
                  ALLOCATE (dovlp3(orb_basis%nsgf, orb_basis%nsgf, ri_basis%nsgf, 3))
                  CALL int_overlap_aba_shg(ovlp3, dovlp3, (/0.0_dp, 0.0_dp, 0.0_dp/), &
                                           orb_basis, orb_basis, ri_basis, scon_orb, &
                                           scon_mix, orb_index, ri_index, &
                                           lri_env%cg_shg%cg_coeff, &
                                           lri_env%cg_shg%cg_none0_list, &
                                           lri_env%cg_shg%ncg_none0, &
                                           calculate_forces=.FALSE.)
                  lri_env%bas_prop(ikind)%ovlp3 => ovlp3
                  DEALLOCATE (orb_norm_r, ri_norm_r, ri_norm_s)
                  DEALLOCATE (dovlp3)
               END IF
            END DO
         END IF
      END IF

   END SUBROUTINE lri_basis_init

! **************************************************************************************************
!> \brief normalization for a contracted Gaussian s-function,
!>        spherical = cartesian Gaussian for s-functions
!> \param basis ...
!> \param norm ...
! **************************************************************************************************
   SUBROUTINE basis_norm_s_func(basis, norm)

      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), DIMENSION(:), POINTER                    :: norm

      CHARACTER(len=*), PARAMETER :: routineN = 'basis_norm_s_func', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ipgf, iset, isgf, ishell, jpgf, l, nbas
      REAL(KIND=dp)                                      :: aai, aaj, cci, ccj, expa, ppl

      NULLIFY (norm)

      nbas = basis%nsgf
      ALLOCATE (norm(nbas))
      norm = 0._dp

      DO iset = 1, basis%nset
         DO ishell = 1, basis%nshell(iset)
            l = basis%l(ishell, iset)
            IF (l /= 0) CYCLE
            expa = 0.5_dp*REAL(2*l+3, dp)
            ppl = pi**(3._dp/2._dp)
            DO isgf = basis%first_sgf(ishell, iset), basis%last_sgf(ishell, iset)
               DO ipgf = 1, basis%npgf(iset)
                  cci = basis%gcc(ipgf, ishell, iset)
                  aai = basis%zet(ipgf, iset)
                  DO jpgf = 1, basis%npgf(iset)
                     ccj = basis%gcc(jpgf, ishell, iset)
                     aaj = basis%zet(jpgf, iset)
                     norm(isgf) = norm(isgf)+cci*ccj*ppl/(aai+aaj)**expa
                  END DO
               END DO
               norm(isgf) = 1.0_dp/SQRT(norm(isgf))
            END DO
         END DO
      END DO

   END SUBROUTINE basis_norm_s_func

! **************************************************************************************************
!> \brief normalization for radial part of contracted spherical Gaussian
!>        functions
!> \param basis ...
!> \param norm ...
! **************************************************************************************************
   SUBROUTINE basis_norm_radial(basis, norm)

      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), DIMENSION(:), POINTER                    :: norm

      CHARACTER(len=*), PARAMETER :: routineN = 'basis_norm_radial', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ipgf, iset, isgf, ishell, jpgf, l, nbas
      REAL(KIND=dp)                                      :: aai, aaj, cci, ccj, expa, ppl

      NULLIFY (norm)

      nbas = basis%nsgf
      ALLOCATE (norm(nbas))
      norm = 0._dp

      DO iset = 1, basis%nset
         DO ishell = 1, basis%nshell(iset)
            l = basis%l(ishell, iset)
            expa = 0.5_dp*REAL(2*l+3, dp)
            ppl = fac(2*l+2)*SQRT(pi)/2._dp**REAL(2*l+3, dp)/fac(l+1)
            DO isgf = basis%first_sgf(ishell, iset), basis%last_sgf(ishell, iset)
               DO ipgf = 1, basis%npgf(iset)
                  cci = basis%gcc(ipgf, ishell, iset)
                  aai = basis%zet(ipgf, iset)
                  DO jpgf = 1, basis%npgf(iset)
                     ccj = basis%gcc(jpgf, ishell, iset)
                     aaj = basis%zet(jpgf, iset)
                     norm(isgf) = norm(isgf)+cci*ccj*ppl/(aai+aaj)**expa
                  END DO
               END DO
               norm(isgf) = 1.0_dp/SQRT(norm(isgf))
            END DO
         END DO
      END DO

   END SUBROUTINE basis_norm_radial

!*****************************************************************************
!> \brief integral over a single (contracted) lri auxiliary basis function,
!>        integral is zero for all but s-functions
!> \param basis ...
!> \param int_aux ...
!> \param norm ...
! **************************************************************************************************
   SUBROUTINE basis_int(basis, int_aux, norm)

      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), DIMENSION(:), POINTER                    :: int_aux, norm

      CHARACTER(len=*), PARAMETER :: routineN = 'basis_int', routineP = moduleN//':'//routineN

      INTEGER                                            :: ipgf, iset, isgf, ishell, l, nbas
      REAL(KIND=dp)                                      :: aa, cc, pp

      nbas = basis%nsgf
      ALLOCATE (int_aux(nbas))
      int_aux = 0._dp

      DO iset = 1, basis%nset
         DO ishell = 1, basis%nshell(iset)
            l = basis%l(ishell, iset)
            IF (l /= 0) CYCLE
            DO isgf = basis%first_sgf(ishell, iset), basis%last_sgf(ishell, iset)
               DO ipgf = 1, basis%npgf(iset)
                  cc = basis%gcc(ipgf, ishell, iset)
                  aa = basis%zet(ipgf, iset)
                  pp = (pi/aa)**(3._dp/2._dp)
                  int_aux(isgf) = int_aux(isgf)+norm(isgf)*cc*pp
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE basis_int

!*****************************************************************************
!> \brief self-overlap of lri basis for contracted spherical Gaussians.
!>        Overlap of radial part. Norm contains only normalization of radial
!>        part. Norm and overlap of spherical harmonics not explicitly
!>        calculated since this cancels for the self-overlap anyway.
!> \param basis ...
!> \param ovlp ...
!> \param norm ...
! **************************************************************************************************
   SUBROUTINE basis_ovlp(basis, ovlp, norm)

      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), DIMENSION(:, :), POINTER                 :: ovlp
      REAL(dp), DIMENSION(:), POINTER                    :: norm

      CHARACTER(len=*), PARAMETER :: routineN = 'basis_ovlp', routineP = moduleN//':'//routineN

      INTEGER                                            :: ipgf, iset, isgf, ishell, jpgf, jset, &
                                                            jsgf, jshell, l, li, lj, m_i, m_j, nbas
      REAL(KIND=dp)                                      :: aai, aaj, cci, ccj, expa, norm_i, &
                                                            norm_j, oo, ppl

      nbas = basis%nsgf
      ALLOCATE (ovlp(nbas, nbas))
      ovlp = 0._dp

      DO iset = 1, basis%nset
         DO ishell = 1, basis%nshell(iset)
            li = basis%l(ishell, iset)
            DO jset = 1, basis%nset
               DO jshell = 1, basis%nshell(jset)
                  lj = basis%l(jshell, jset)
                  IF (li == lj) THEN
                     l = li
                     expa = 0.5_dp*REAL(2*l+3, dp)
                     ppl = fac(2*l+2)*SQRT(pi)/2._dp**REAL(2*l+3, dp)/fac(l+1)
                     DO isgf = basis%first_sgf(ishell, iset), basis%last_sgf(ishell, iset)
                        m_i = basis%m(isgf)
                        DO jsgf = basis%first_sgf(jshell, jset), basis%last_sgf(jshell, jset)
                           m_j = basis%m(jsgf)
                           IF (m_i == m_j) THEN
                              DO ipgf = 1, basis%npgf(iset)
                                 cci = basis%gcc(ipgf, ishell, iset)
                                 aai = basis%zet(ipgf, iset)
                                 norm_i = norm(isgf)
                                 DO jpgf = 1, basis%npgf(jset)
                                    ccj = basis%gcc(jpgf, jshell, jset)
                                    aaj = basis%zet(jpgf, jset)
                                    oo = 1._dp/(aai+aaj)**expa
                                    norm_j = norm(jsgf)
                                    ovlp(isgf, jsgf) = ovlp(isgf, jsgf)+norm_i*norm_j*ppl*cci*ccj*oo
                                 END DO
                              END DO
                           ENDIF
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE basis_ovlp

END MODULE lri_environment_init
