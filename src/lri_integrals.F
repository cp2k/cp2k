!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates integrals for LRIGPW method
!>        lri : local resolution of the identity
!> \par History
!>      created JGH [08.2012]
!>      Dorothea Golze [02.2014] (1) extended, re-structured, cleaned
!>                               (2) heavily debugged
!>      split off JGH [11.2017]
!> \authors JGH
!>          Dorothea Golze
! **************************************************************************************************
MODULE lri_integrals
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE generic_os_integrals,            ONLY: int_overlap_ab_os,&
                                              int_overlap_aba_os,&
                                              int_overlap_abb_os
   USE generic_shg_integrals,           ONLY: get_abb_same_kind,&
                                              int_overlap_ab_shg_low,&
                                              int_overlap_aba_shg_low,&
                                              int_overlap_abb_shg_low,&
                                              lri_precalc_angular_shg_part
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: lri_environment_type,&
                                              lri_int_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
   TYPE int_type
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE             :: sabint
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE             :: sooint
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE          :: abaint
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE          :: abbint
   END TYPE int_type
   !
   TYPE dint_type
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE          :: dsabint
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE          :: dsooint
      REAL(KIND=dp), DIMENSION(:, :, :, :), ALLOCATABLE       :: dabaint
      REAL(KIND=dp), DIMENSION(:, :, :, :), ALLOCATABLE       :: dabbint
   END TYPE dint_type
! **************************************************************************************************

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_integrals'

   PUBLIC :: lri_int, lri_int2, lri_dint, lri_dint2, int_type, dint_type, allocate_int_type, deallocate_int_type

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief calcuates the lri integrals using solid harmonic Gaussians
!> \param lri_env ...
!> \param lrii ...
!> \param rab distance vector
!> \param obasa orb basis on A
!> \param obasb orb basis on B
!> \param fbasa aux basis on A
!> \param fbasb aux basis on B
!> \param iatom index atom A
!> \param jatom index atom B
!> \param ikind kind atom A
!> \param jkind kind atom B
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE lri_int(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                      iatom, jatom, ikind, jkind, calculate_forces)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_int_type), POINTER                        :: lrii
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: obasa, obasb, fbasa, fbasb
      INTEGER, INTENT(IN)                                :: iatom, jatom, ikind, jkind
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int', routineP = moduleN//':'//routineN

      IF (lri_env%use_shg_integrals) THEN
         CALL lri_int_shg(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                          iatom, jatom, ikind, jkind, calculate_forces)
      ELSE
         CALL lri_int_os(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                         iatom, jatom, ikind, calculate_forces)
      ENDIF

   END SUBROUTINE lri_int

! **************************************************************************************************
!> \brief calcuates the lri integrals using solid harmonic Gaussians
!> \param lri_env ...
!> \param lrii ...
!> \param lriint ...
!> \param rab distance vector
!> \param obasa orb basis on A
!> \param obasb orb basis on B
!> \param fbasa aux basis on A
!> \param fbasb aux basis on B
!> \param iatom index atom A
!> \param jatom index atom B
!> \param ikind kind atom A
!> \param jkind kind atom B
! **************************************************************************************************
   SUBROUTINE lri_int2(lri_env, lrii, lriint, rab, obasa, obasb, fbasa, fbasb, iatom, jatom, ikind, jkind)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_int_type), POINTER                        :: lrii
      TYPE(int_type)                                     :: lriint
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: obasa, obasb, fbasa, fbasb
      INTEGER, INTENT(IN)                                :: iatom, jatom, ikind, jkind

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int2', routineP = moduleN//':'//routineN

      INTEGER                                            :: nba, nfa
      INTEGER, DIMENSION(:, :, :), POINTER               :: fba_index, fbb_index, oba_index, &
                                                            obb_index
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dummy1, Waux_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: dummy2, dWaux_mat
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scon_fba, scon_fbb, scon_oba, scon_obb
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: scona_mix, sconb_mix

      dab = SQRT(SUM(rab*rab))
      IF (iatom == jatom .AND. dab < lri_env%delta) THEN
         CPASSERT(ALLOCATED(lriint%sabint))
         CPASSERT(ALLOCATED(lriint%sooint))
         CPASSERT(ALLOCATED(lriint%abaint))
         nba = obasa%nsgf
         nfa = fbasa%nsgf
         lriint%sabint(1:nfa, 1:nfa) = lri_env%bas_prop(ikind)%ri_ovlp(1:nfa, 1:nfa)
         lriint%sooint(1:nba, 1:nba) = lri_env%bas_prop(ikind)%orb_ovlp(1:nba, 1:nba)
         lriint%abaint(1:nba, 1:nba, 1:nfa) = lri_env%bas_prop(ikind)%ovlp3
      ELSE
         IF (lri_env%use_shg_integrals) THEN
            scon_oba => lri_env%bas_prop(ikind)%scon_orb
            scon_obb => lri_env%bas_prop(jkind)%scon_orb
            scon_fba => lri_env%bas_prop(ikind)%scon_ri
            scon_fbb => lri_env%bas_prop(jkind)%scon_ri
            scona_mix => lri_env%bas_prop(ikind)%scon_mix
            sconb_mix => lri_env%bas_prop(jkind)%scon_mix
            oba_index => lri_env%bas_prop(ikind)%orb_index
            fba_index => lri_env%bas_prop(ikind)%ri_index
            obb_index => lri_env%bas_prop(jkind)%orb_index
            fbb_index => lri_env%bas_prop(jkind)%ri_index
            !
            CALL lri_precalc_angular_shg_part(obasa, obasb, fbasa, fbasb, rab, Waux_mat, dWaux_mat, &
                                              .FALSE.)
            !*** (fa,fb)
            IF (ALLOCATED(lriint%sabint)) THEN
               CALL int_overlap_ab_shg_low(lriint%sabint, dummy1, rab, fbasa, fbasb, scon_fba, scon_fbb, &
                                           Waux_mat, dWaux_mat, .TRUE., .FALSE., contraction_high=.FALSE.)
            END IF
            !*** (a,b)
            IF (ALLOCATED(lriint%sooint)) THEN
               CALL int_overlap_ab_shg_low(lriint%sooint, dummy1, rab, obasa, obasb, scon_oba, scon_obb, &
                                           Waux_mat, dWaux_mat, .TRUE., .FALSE., contraction_high=.TRUE.)
            END IF
            !*** (a,b,fa)
            IF (ALLOCATED(lriint%abaint)) THEN
               CALL int_overlap_aba_shg_low(lriint%abaint, dummy2, rab, obasa, obasb, fbasa, &
                                            scon_obb, scona_mix, oba_index, fba_index, &
                                            lri_env%cg_shg%cg_coeff, lri_env%cg_shg%cg_none0_list, &
                                            lri_env%cg_shg%ncg_none0, &
                                            Waux_mat, dWaux_mat, .TRUE., .FALSE.)
            END IF
            !*** (a,b,fb)
            IF (ALLOCATED(lriint%abbint)) THEN
               IF (ikind == jkind) THEN
                  CPASSERT(ALLOCATED(lriint%abaint))
                  CALL get_abb_same_kind(lriint%abbint, dummy2, lriint%abaint, dummy2, &
                                         rab, obasa, fbasa, .TRUE., .FALSE.)
               ELSE
                  CALL int_overlap_abb_shg_low(lriint%abbint, dummy2, rab, obasa, obasb, fbasb, &
                                               scon_oba, sconb_mix, obb_index, fbb_index, &
                                               lri_env%cg_shg%cg_coeff, lri_env%cg_shg%cg_none0_list, &
                                               lri_env%cg_shg%ncg_none0, &
                                               Waux_mat, dWaux_mat, .TRUE., .FALSE.)
               ENDIF
            ENDIF
            DEALLOCATE (Waux_mat, dWaux_mat)
         ELSE
            !*** (fa,fb)
            IF (ALLOCATED(lriint%sabint)) THEN
               CALL int_overlap_ab_os(lriint%sabint, dummy1, rab, fbasa, fbasb, &
                                      .FALSE., lri_env%debug, lrii%dmax_ab)
            END IF
            !*** (a,b)
            IF (ALLOCATED(lriint%sooint)) THEN
               CALL int_overlap_ab_os(lriint%sooint, dummy1, rab, obasa, obasb, &
                                      .FALSE., lri_env%debug, lrii%dmax_oo)
            END IF
            !*** (a,b,fa)
            IF (ALLOCATED(lriint%abaint)) THEN
               CALL int_overlap_aba_os(lriint%abaint, dummy2, rab, obasa, obasb, fbasa, &
                                       .FALSE., lri_env%debug, lrii%dmax_aba)
            END IF
            !*** (a,b,fb)
            IF (ALLOCATED(lriint%abbint)) THEN
               CALL int_overlap_abb_os(lriint%abbint, dummy2, rab, obasa, obasb, fbasb, &
                                       .FALSE., lri_env%debug, lrii%dmax_abb)
            END IF
         ENDIF
      ENDIF

   END SUBROUTINE lri_int2

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param lri_env ...
!> \param lrii ...
!> \param rab ...
!> \param obasa ...
!> \param obasb ...
!> \param fbasa ...
!> \param fbasb ...
!> \param iatom ...
!> \param jatom ...
!> \param ikind ...
!> \param jkind ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE lri_dint(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                       iatom, jatom, ikind, jkind, calculate_forces)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_int_type), POINTER                        :: lrii
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: obasa, obasb, fbasa, fbasb
      INTEGER, INTENT(IN)                                :: iatom, jatom, ikind, jkind
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_dint', routineP = moduleN//':'//routineN

      IF (lri_env%use_shg_integrals) THEN
         CALL lri_int_shg(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                          iatom, jatom, ikind, jkind, calculate_forces)
      ELSE
         CALL lri_int_os(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                         iatom, jatom, ikind, calculate_forces)
      ENDIF

   END SUBROUTINE lri_dint

! **************************************************************************************************
!> \brief ...
!> \param lri_env ...
!> \param lrii ...
!> \param lridint ...
!> \param rab ...
!> \param obasa ...
!> \param obasb ...
!> \param fbasa ...
!> \param fbasb ...
!> \param iatom ...
!> \param jatom ...
!> \param ikind ...
!> \param jkind ...
! **************************************************************************************************
   SUBROUTINE lri_dint2(lri_env, lrii, lridint, rab, obasa, obasb, fbasa, fbasb, &
                        iatom, jatom, ikind, jkind)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_int_type), POINTER                        :: lrii
      TYPE(dint_type)                                    :: lridint
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: obasa, obasb, fbasa, fbasb
      INTEGER, INTENT(IN)                                :: iatom, jatom, ikind, jkind

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_dint2', routineP = moduleN//':'//routineN

      INTEGER                                            :: nba, nbb, nfa, nfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: fba_index, fbb_index, oba_index, &
                                                            obb_index
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dummy1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dummy2, Waux_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: dWaux_mat
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scon_fba, scon_fbb, scon_oba, scon_obb
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: scona_mix, sconb_mix

      dab = SQRT(SUM(rab*rab))
      nba = obasa%nsgf
      nbb = obasb%nsgf
      nfa = fbasa%nsgf
      nfb = fbasb%nsgf

      IF (iatom == jatom .AND. dab < lri_env%delta) THEN
         ! nothing to do
      ELSE
         IF (lrii%calc_force_pair) THEN
            IF (lri_env%use_shg_integrals) THEN
               scon_oba => lri_env%bas_prop(ikind)%scon_orb
               scon_obb => lri_env%bas_prop(jkind)%scon_orb
               scon_fba => lri_env%bas_prop(ikind)%scon_ri
               scon_fbb => lri_env%bas_prop(jkind)%scon_ri
               scona_mix => lri_env%bas_prop(ikind)%scon_mix
               sconb_mix => lri_env%bas_prop(jkind)%scon_mix
               oba_index => lri_env%bas_prop(ikind)%orb_index
               fba_index => lri_env%bas_prop(ikind)%ri_index
               obb_index => lri_env%bas_prop(jkind)%orb_index
               fbb_index => lri_env%bas_prop(jkind)%ri_index

               CALL lri_precalc_angular_shg_part(obasa, obasb, fbasa, fbasb, rab, &
                                                 Waux_mat, dWaux_mat, .TRUE.)
               !*** (fa,fb)
               IF (ALLOCATED(lridint%dsabint)) THEN
                  CALL int_overlap_ab_shg_low(dummy1, lridint%dsabint, rab, fbasa, fbasb, scon_fba, scon_fbb, &
                                              Waux_mat, dWaux_mat, .FALSE., .TRUE., contraction_high=.FALSE.)
               END IF
               !*** (a,b)
               IF (ALLOCATED(lridint%dsooint)) THEN
                  CALL int_overlap_ab_shg_low(dummy1, lridint%dsooint, rab, obasa, obasb, scon_oba, scon_obb, &
                                              Waux_mat, dWaux_mat, .FALSE., .TRUE., contraction_high=.TRUE.)
               END IF
               !*** (a,b,fa)
               IF (ALLOCATED(lridint%dabaint)) THEN
                  CALL int_overlap_aba_shg_low(dummy2, lridint%dabaint, rab, obasa, obasb, fbasa, &
                                               scon_obb, scona_mix, oba_index, fba_index, &
                                               lri_env%cg_shg%cg_coeff, lri_env%cg_shg%cg_none0_list, &
                                               lri_env%cg_shg%ncg_none0, &
                                               Waux_mat, dWaux_mat, .FALSE., .TRUE.)
               END IF
               !*** (a,b,fb)
               IF (ALLOCATED(lridint%dabbint)) THEN
                  IF (ikind == jkind) THEN
                     CPASSERT(ALLOCATED(lridint%dabaint))
                     CALL get_abb_same_kind(dummy2, lridint%dabbint, dummy2, lridint%dabaint, &
                                            rab, obasa, fbasa, .FALSE., .TRUE.)
                  ELSE
                     CALL int_overlap_abb_shg_low(dummy2, lridint%dabbint, rab, obasa, obasb, fbasb, &
                                                  scon_oba, sconb_mix, obb_index, fbb_index, &
                                                  lri_env%cg_shg%cg_coeff, lri_env%cg_shg%cg_none0_list, &
                                                  lri_env%cg_shg%ncg_none0, &
                                                  Waux_mat, dWaux_mat, .FALSE., .TRUE.)
                  ENDIF
               ENDIF
               DEALLOCATE (Waux_mat, dWaux_mat)

            ELSE

               !*** (fa,fb)
               IF (ALLOCATED(lridint%dsabint)) THEN
                  ALLOCATE (dummy1(nfa, nfb))
                  CALL int_overlap_ab_os(dummy1, lridint%dsabint, rab, fbasa, fbasb, &
                                         .TRUE., lri_env%debug, lrii%dmax_ab)
                  DEALLOCATE (dummy1)
               END IF
               !*** (a,b)
               IF (ALLOCATED(lridint%dsooint)) THEN
                  ALLOCATE (dummy1(nba, nbb))
                  CALL int_overlap_ab_os(dummy1, lridint%dsooint, rab, obasa, obasb, &
                                         .TRUE., lri_env%debug, lrii%dmax_oo)
                  DEALLOCATE (dummy1)
               END IF
               !*** (a,b,fa)
               IF (ALLOCATED(lridint%dabaint)) THEN
                  ALLOCATE (dummy2(nba, nbb, nfa))
                  CALL int_overlap_aba_os(dummy2, lridint%dabaint, rab, obasa, obasb, fbasa, &
                                          .TRUE., lri_env%debug, lrii%dmax_aba)
                  DEALLOCATE (dummy2)
               END IF
               !*** (a,b,fb)
               IF (ALLOCATED(lridint%dabbint)) THEN
                  ALLOCATE (dummy2(nba, nbb, nfb))
                  CALL int_overlap_abb_os(dummy2, lridint%dabbint, rab, obasa, obasb, fbasb, &
                                          .TRUE., lri_env%debug, lrii%dmax_abb)
                  DEALLOCATE (dummy2)
               END IF
            ENDIF
         ENDIF
      ENDIF

   END SUBROUTINE lri_dint2

! **************************************************************************************************
!> \brief calculates the lri intergrals according to the Obara-Saika (OS)
!>        scheme
!> \param lri_env ...
!> \param lrii ...
!> \param rab distance vector
!> \param obasa orb basis on center A
!> \param obasb orb basis on center B
!> \param fbasa aux basis on center A
!> \param fbasb aux basis on center B
!> \param iatom index atom A
!> \param jatom index atom B
!> \param ikind kind atom A
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE lri_int_os(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                         iatom, jatom, ikind, calculate_forces)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_int_type), POINTER                        :: lrii
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: obasa, obasb, fbasa, fbasb
      INTEGER, INTENT(IN)                                :: iatom, jatom, ikind
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_os', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nba, nbb, nfa, nfb
      LOGICAL                                            :: do_force
      REAL(KIND=dp)                                      :: dab

      CALL timeset(routineN, handle)

      dab = SQRT(SUM(rab*rab))
      nba = obasa%nsgf
      nbb = obasb%nsgf
      nfa = fbasa%nsgf
      nfb = fbasb%nsgf

      !*** calculate overlap integrals; for iatom=jatom this is the self-overlap
      IF (iatom == jatom .AND. dab < lri_env%delta) THEN
         !*** (fa,fa)
         lrii%sab(1:nfa, 1:nfa) = lri_env%bas_prop(ikind)%ri_ovlp(1:nfa, 1:nfa)
         !*** (a,a)
         lrii%soo(1:nba, 1:nba) = lri_env%bas_prop(ikind)%orb_ovlp(1:nba, 1:nba)
         !*** (a,a,fa)
         CALL int_overlap_aba_os(lrii%abaint, rab=rab, oba=obasa, obb=obasb, &
                                 fba=fbasa, calculate_forces=.FALSE., debug=lri_env%debug, &
                                 dmax=lrii%dmax_aba)
         IF (calculate_forces) THEN
            lrii%dsab = 0._dp
            lrii%dsoo = 0._dp
            lrii%dabdaint = 0.0_dp
            lrii%dabbint = 0.0_dp
         END IF
      ELSE
         IF (calculate_forces) THEN
            do_force = lrii%calc_force_pair
         ELSE
            do_force = .FALSE.
         END IF
         !*** (fa,fb)
         CALL int_overlap_ab_os(lrii%sab, lrii%dsab, rab, fbasa, fbasb, &
                                do_force, lri_env%debug, lrii%dmax_ab)
         !*** (a,b)
         CALL int_overlap_ab_os(lrii%soo, lrii%dsoo, rab, obasa, obasb, &
                                do_force, lri_env%debug, lrii%dmax_oo)
         !*** (a,b,fa)
         CALL int_overlap_aba_os(lrii%abaint, lrii%dabdaint, rab, obasa, obasb, fbasa, &
                                 do_force, lri_env%debug, lrii%dmax_aba)
         !*** (a,b,fb)
         CALL int_overlap_abb_os(lrii%abbint, lrii%dabbint, rab, obasa, obasb, fbasb, &
                                 do_force, lri_env%debug, lrii%dmax_abb)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE lri_int_os

! **************************************************************************************************
!> \brief calcuates the lri integrals using solid harmonic Gaussians
!> \param lri_env ...
!> \param lrii ...
!> \param rab distance vector
!> \param obasa orb basis on A
!> \param obasb orb basis on B
!> \param fbasa aux basis on A
!> \param fbasb aux basis on B
!> \param iatom index atom A
!> \param jatom index atom B
!> \param ikind kind atom A
!> \param jkind kind atom B
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE lri_int_shg(lri_env, lrii, rab, obasa, obasb, fbasa, fbasb, &
                          iatom, jatom, ikind, jkind, calculate_forces)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_int_type), POINTER                        :: lrii
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      TYPE(gto_basis_set_type), POINTER                  :: obasa, obasb, fbasa, fbasb
      INTEGER, INTENT(IN)                                :: iatom, jatom, ikind, jkind
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_shg', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nba, nbb, nfa, nfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: fba_index, fbb_index, oba_index, &
                                                            obb_index
      LOGICAL                                            :: do_force
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Waux_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: dWaux_mat
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scon_fba, scon_fbb, scon_oba, scon_obb
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: scona_mix, sconb_mix

      CALL timeset(routineN, handle)
      NULLIFY (scon_oba, scon_obb, scon_fba, scon_fbb, scona_mix, sconb_mix, &
               oba_index, obb_index, fba_index, fbb_index)
      dab = SQRT(SUM(rab*rab))
      nba = obasa%nsgf
      nbb = obasb%nsgf
      nfa = fbasa%nsgf
      nfb = fbasb%nsgf

      !*** calculate overlap integrals; for iatom=jatom this is the self-overlap
      IF (iatom == jatom .AND. dab < lri_env%delta) THEN
         !*** (fa,fa)
         lrii%sab(1:nfa, 1:nfa) = lri_env%bas_prop(ikind)%ri_ovlp(1:nfa, 1:nfa)
         !*** (a,a)
         lrii%soo(1:nba, 1:nba) = lri_env%bas_prop(ikind)%orb_ovlp(1:nba, 1:nba)
         !*** (a,a,fa)
         lrii%abaint(1:nba, 1:nba, 1:nfa) = lri_env%bas_prop(ikind)%ovlp3
         IF (calculate_forces) THEN
            lrii%dsab = 0._dp
            lrii%dsoo = 0._dp
            lrii%dabdaint = 0.0_dp
            lrii%dabbint = 0.0_dp
         END IF
      ELSE
         IF (calculate_forces) THEN
            do_force = lrii%calc_force_pair
         ELSE
            do_force = .FALSE.
         END IF
         scon_oba => lri_env%bas_prop(ikind)%scon_orb
         scon_obb => lri_env%bas_prop(jkind)%scon_orb
         scon_fba => lri_env%bas_prop(ikind)%scon_ri
         scon_fbb => lri_env%bas_prop(jkind)%scon_ri
         scona_mix => lri_env%bas_prop(ikind)%scon_mix
         sconb_mix => lri_env%bas_prop(jkind)%scon_mix
         oba_index => lri_env%bas_prop(ikind)%orb_index
         fba_index => lri_env%bas_prop(ikind)%ri_index
         obb_index => lri_env%bas_prop(jkind)%orb_index
         fbb_index => lri_env%bas_prop(jkind)%ri_index
         CALL lri_precalc_angular_shg_part(obasa, obasb, fbasa, fbasb, rab, Waux_mat, dWaux_mat, &
                                           do_force)
         !*** (fa,fb)
         CALL int_overlap_ab_shg_low(lrii%sab, lrii%dsab, rab, fbasa, fbasb, scon_fba, scon_fbb, &
                                     Waux_mat, dWaux_mat, .TRUE., do_force, contraction_high=.FALSE.)
         !*** (a,b)
         CALL int_overlap_ab_shg_low(lrii%soo, lrii%dsoo, rab, obasa, obasb, scon_oba, scon_obb, &
                                     Waux_mat, dWaux_mat, .TRUE., do_force, contraction_high=.TRUE.)
         !*** (a,b,fa)
         CALL int_overlap_aba_shg_low(lrii%abaint, lrii%dabdaint, rab, obasa, obasb, fbasa, &
                                      scon_obb, scona_mix, oba_index, fba_index, &
                                      lri_env%cg_shg%cg_coeff, lri_env%cg_shg%cg_none0_list, &
                                      lri_env%cg_shg%ncg_none0, &
                                      Waux_mat, dWaux_mat, .TRUE., do_force)
         !*** (a,b,fb)
         IF (ikind == jkind) THEN
            CALL get_abb_same_kind(lrii%abbint, lrii%dabbint, lrii%abaint, lrii%dabdaint, &
                                   rab, obasa, fbasa, .TRUE., do_force)
         ELSE
            CALL int_overlap_abb_shg_low(lrii%abbint, lrii%dabbint, rab, obasa, obasb, fbasb, &
                                         scon_oba, sconb_mix, obb_index, fbb_index, &
                                         lri_env%cg_shg%cg_coeff, lri_env%cg_shg%cg_none0_list, &
                                         lri_env%cg_shg%ncg_none0, &
                                         Waux_mat, dWaux_mat, .TRUE., do_force)
         ENDIF

         DEALLOCATE (Waux_mat, dWaux_mat)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE lri_int_shg

! **************************************************************************************************
!> \brief ...
!> \param lriint ...
!> \param lridint ...
!> \param nba ...
!> \param nbb ...
!> \param nfa ...
!> \param nfb ...
!> \param skip_sab ...
!> \param skip_soo ...
!> \param skip_aba ...
!> \param skip_abb ...
!> \param skip_dsab ...
!> \param skip_dsoo ...
!> \param skip_daba ...
!> \param skip_dabb ...
! **************************************************************************************************
   SUBROUTINE allocate_int_type(lriint, lridint, nba, nbb, nfa, nfb, &
                                skip_sab, skip_soo, skip_aba, skip_abb, &
                                skip_dsab, skip_dsoo, skip_daba, skip_dabb)

      TYPE(int_type), INTENT(INOUT), OPTIONAL            :: lriint
      TYPE(dint_type), INTENT(INOUT), OPTIONAL           :: lridint
      INTEGER, INTENT(IN)                                :: nba, nbb, nfa, nfb
      LOGICAL, INTENT(IN), OPTIONAL                      :: skip_sab, skip_soo, skip_aba, skip_abb, &
                                                            skip_dsab, skip_dsoo, skip_daba, &
                                                            skip_dabb

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_int_type', &
         routineP = moduleN//':'//routineN

      LOGICAL                                            :: do_aba, do_abb, do_daba, do_dabb, &
                                                            do_dsab, do_dsoo, do_sab, do_soo

      IF (PRESENT(lriint)) THEN
         do_sab = .TRUE.
         IF (PRESENT(skip_sab)) do_sab = .NOT. skip_sab
         do_soo = .TRUE.
         IF (PRESENT(skip_soo)) do_soo = .NOT. skip_soo
         do_aba = .TRUE.
         IF (PRESENT(skip_aba)) do_aba = .NOT. skip_aba
         do_abb = .TRUE.
         IF (PRESENT(skip_abb)) do_abb = .NOT. skip_abb
         !
         IF (do_sab) THEN
            IF (ALLOCATED(lriint%sabint)) DEALLOCATE (lriint%sabint)
            ALLOCATE (lriint%sabint(nfa, nfb))
            lriint%sabint = 0.0_dp
         END IF
         IF (do_soo) THEN
            IF (ALLOCATED(lriint%sooint)) DEALLOCATE (lriint%sooint)
            ALLOCATE (lriint%sooint(nba, nbb))
            lriint%sooint = 0.0_dp
         END IF
         IF (do_aba) THEN
            IF (ALLOCATED(lriint%abaint)) DEALLOCATE (lriint%abaint)
            ALLOCATE (lriint%abaint(nba, nbb, nfa))
            lriint%abaint = 0.0_dp
         END IF
         IF (do_abb) THEN
            IF (ALLOCATED(lriint%abbint)) DEALLOCATE (lriint%abbint)
            ALLOCATE (lriint%abbint(nba, nbb, nfb))
            lriint%abbint = 0.0_dp
         END IF
      END IF
      !
      IF (PRESENT(lridint)) THEN
         do_dsab = .TRUE.
         IF (PRESENT(skip_dsab)) do_dsab = .NOT. skip_dsab
         do_dsoo = .TRUE.
         IF (PRESENT(skip_dsoo)) do_dsoo = .NOT. skip_dsoo
         do_daba = .TRUE.
         IF (PRESENT(skip_daba)) do_daba = .NOT. skip_daba
         do_dabb = .TRUE.
         IF (PRESENT(skip_dabb)) do_dabb = .NOT. skip_dabb
         !
         IF (do_dsab) THEN
            IF (ALLOCATED(lridint%dsabint)) DEALLOCATE (lridint%dsabint)
            ALLOCATE (lridint%dsabint(nfa, nfb, 3))
            lridint%dsabint = 0.0_dp
         END IF
         IF (do_dsoo) THEN
            IF (ALLOCATED(lridint%dsooint)) DEALLOCATE (lridint%dsooint)
            ALLOCATE (lridint%dsooint(nba, nbb, 3))
            lridint%dsooint = 0.0_dp
         END IF
         IF (do_daba) THEN
            IF (ALLOCATED(lridint%dabaint)) DEALLOCATE (lridint%dabaint)
            ALLOCATE (lridint%dabaint(nba, nbb, nfa, 3))
            lridint%dabaint = 0.0_dp
         END IF
         IF (do_dabb) THEN
            IF (ALLOCATED(lridint%dabbint)) DEALLOCATE (lridint%dabbint)
            ALLOCATE (lridint%dabbint(nba, nbb, nfb, 3))
            lridint%dabbint = 0.0_dp
         END IF
      END IF

   END SUBROUTINE allocate_int_type

! **************************************************************************************************
!> \brief ...
!> \param lriint ...
!> \param lridint ...
! **************************************************************************************************
   SUBROUTINE deallocate_int_type(lriint, lridint)

      TYPE(int_type), INTENT(INOUT), OPTIONAL            :: lriint
      TYPE(dint_type), INTENT(INOUT), OPTIONAL           :: lridint

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_int_type', &
         routineP = moduleN//':'//routineN

      IF (PRESENT(lriint)) THEN
         IF (ALLOCATED(lriint%sabint)) DEALLOCATE (lriint%sabint)
         IF (ALLOCATED(lriint%sooint)) DEALLOCATE (lriint%sooint)
         IF (ALLOCATED(lriint%abaint)) DEALLOCATE (lriint%abaint)
         IF (ALLOCATED(lriint%abbint)) DEALLOCATE (lriint%abbint)
      END IF
      !
      IF (PRESENT(lridint)) THEN
         IF (ALLOCATED(lridint%dsabint)) DEALLOCATE (lridint%dsabint)
         IF (ALLOCATED(lridint%dsooint)) DEALLOCATE (lridint%dsooint)
         IF (ALLOCATED(lridint%dabaint)) DEALLOCATE (lridint%dabaint)
         IF (ALLOCATED(lridint%dabbint)) DEALLOCATE (lridint%dabbint)
      END IF

   END SUBROUTINE deallocate_int_type

END MODULE lri_integrals
