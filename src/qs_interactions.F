!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate the interaction radii for the operator matrix calculation.
!> \par History
!>      Joost VandeVondele : added exp_radius_very_extended
!>      24.09.2002 overloaded init_interaction_radii for KG use (gt)
!>      27.02.2004 integrated KG into QS version of init_interaction_radii (jgh)
!>      07.03.2006 new routine to extend kind interaction radius for semi-empiric (jgh)
!> \author MK (12.07.2000)
! *****************************************************************************
MODULE qs_interactions

  USE ao_util,                         ONLY: exp_radius
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: geminal_basis_set_type,&
                                             get_geminal_basis_set,&
                                             get_gto_basis_set,&
                                             gto_basis_set_type,&
                                             set_geminal_basis_set,&
                                             set_gto_basis_set
  USE cp_control_types,                ONLY: qs_control_type,&
                                             semi_empirical_control_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type,&
                                             set_potential
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE paw_proj_set_types,              ONLY: get_paw_proj_set,&
                                             paw_proj_set_type,&
                                             set_paw_proj_set
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_interactions'

! *** Public subroutines ***

  PUBLIC :: init_interaction_radii, init_se_nlradius

  PUBLIC :: write_paw_radii, write_ppnl_radii, write_ppl_radii, write_core_charge_radii,&
            write_pgf_orb_radii, write_geminal_radii

CONTAINS

! *****************************************************************************
!> \brief   Initialize all the atomic kind radii for a given threshold value.
!> \param qs_control ...
!> \param atomic_kind_set ...
!> \param error ...
!> \date    24.09.2002
!> \author  GT
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_interaction_radii(qs_control,atomic_kind_set,error)

    TYPE(qs_control_type), INTENT(IN)        :: qs_control
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_interaction_radii', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iexp_ppl, ikind, ip, ipgf, iprj_ppnl, iset, ishell, &
      j, l, lppl, lppnl, lprj, lprj_ppnl, maxl, maxnset, nexp_lpot, nexp_lsd, &
      nexp_ppl, nkind, nset
    INTEGER, DIMENSION(:), POINTER           :: ls, nct_lpot, nct_lsd, npgf, &
                                                nprj, nprj_ppnl, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: lshell
    LOGICAL                                  :: llsd, lpot, p_aux, p_ccr, &
                                                p_gem, p_lri, p_orb, p_paw, &
                                                p_ppl, p_ppnl, p_soft, &
                                                paw_atom
    REAL(KIND=dp) :: alpha_core_charge, alpha_ppl, ccore_charge, cerf_ppl, &
      core_charge_radius, gcca, kind_radius, ppl_radius, ppnl_radius, rcprj, &
      zab, zd, zeta
    REAL(KIND=dp), DIMENSION(:), POINTER     :: alpha_lpot, alpha_lsd, &
                                                alpha_ppnl, cexp_ppl, &
                                                set_radius
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cprj_ppnl, cval_lpot, &
                                                cval_lsd, pgf_radius, &
                                                rzetprj, zet
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: zetg
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set, &
                                                aux_fit_basis_set, lri_basis, &
                                                orb_basis_set, soft_basis
    TYPE(paw_proj_set_type), POINTER         :: paw_proj_set

    CALL timeset(routineN,handle)

    NULLIFY ( atomic_kind, all_potential, gth_potential, &
              orb_basis_set, aux_basis_set, soft_basis, lri_basis, &
              geminal_basis_set, paw_proj_set)
    NULLIFY (npgf, nprj_ppnl, nshell, lshell, nprj)
    NULLIFY (alpha_ppnl, cexp_ppl, set_radius, cprj_ppnl, pgf_radius, zet, gcc)

    p_orb  = .FALSE.
    p_aux  = .FALSE.
    p_ccr  = .FALSE.
    p_ppl  = .FALSE.
    p_ppnl = .FALSE.
    p_soft = .FALSE.
    p_paw  = .FALSE.
    p_gem  = .FALSE.
    p_lri  = .FALSE.

    nkind = SIZE(atomic_kind_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxnset=maxnset)

    nexp_ppl = 0

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           all_potential=all_potential,&
                           gth_potential=gth_potential,&
                           aux_basis_set=aux_basis_set,&
                           aux_fit_basis_set=aux_fit_basis_set,&
                           orb_basis_set=orb_basis_set,&
                           lri_basis_set=lri_basis,&
                           soft_basis_set=soft_basis,&
                           paw_proj_set=paw_proj_set,&
                           geminal_basis_set=geminal_basis_set,&
                           paw_atom=paw_atom)

      ! Calculate the orbital basis function radii ***
      ! For ALL electron this has to come before the calculation of the
      ! radii used to calculate the 3c lists.
      ! Indeed, there the radii of the GTO primitives are needed
      IF (ASSOCIATED(orb_basis_set)) THEN

        p_orb = .TRUE.
        CALL init_interaction_radii_orb_basis(orb_basis_set,qs_control%eps_pgf_orb,&
             qs_control%eps_kg_orb,error=error)

      END IF

      IF (ASSOCIATED(all_potential)) THEN

        p_ccr = .TRUE.
        CALL get_potential(potential=all_potential,&
                           alpha_core_charge=alpha_core_charge,&
                           ccore_charge=ccore_charge)

        ! Calculate the radius of the core charge distribution

        core_charge_radius = exp_radius(0,alpha_core_charge,&
                                        qs_control%eps_core_charge,&
                                        ccore_charge)

        CALL set_potential(potential=all_potential,&
                           core_charge_radius=core_charge_radius)

      ELSE IF (ASSOCIATED(gth_potential)) THEN

        p_ccr = .TRUE.
        p_ppl = .TRUE.
        CALL get_potential(potential=gth_potential,&
                           alpha_core_charge=alpha_core_charge,&
                           ccore_charge=ccore_charge,&
                           alpha_ppl=alpha_ppl,&
                           cerf_ppl=cerf_ppl,&
                           nexp_ppl=nexp_ppl,&
                           cexp_ppl=cexp_ppl,&
                           lpot_present=lpot,&
                           lsd_present=llsd,&
                           lppnl=lppnl,&
                           alpha_ppnl=alpha_ppnl,&
                           nprj_ppnl=nprj_ppnl,&
                           cprj_ppnl=cprj_ppnl)

        ! Calculate the radius of the core charge distribution
        core_charge_radius = exp_radius(0,alpha_core_charge,&
                                        qs_control%eps_core_charge,&
                                        ccore_charge)

        ! Calculate the radii of the local part
        ! of the Goedecker pseudopotential (GTH)
        ppl_radius = exp_radius(0,alpha_ppl,qs_control%eps_ppl,cerf_ppl)

        DO iexp_ppl=1,nexp_ppl
          lppl = 2*(iexp_ppl - 1)
          ppl_radius = MAX(ppl_radius,&
                           exp_radius(lppl,alpha_ppl,&
                                      qs_control%eps_ppl,&
                                      cexp_ppl(iexp_ppl)))
        END DO
        IF(lpot) THEN
          ! extended local potential
          CALL get_potential(potential=gth_potential,&
                             nexp_lpot = nexp_lpot,&
                             alpha_lpot = alpha_lpot,&
                             nct_lpot = nct_lpot,&
                             cval_lpot = cval_lpot)
          DO j=1,nexp_lpot
            DO i=1,nct_lpot(j)
              lppl = 2*(i-1)
              ppl_radius = MAX(ppl_radius,&
                           exp_radius(lppl,alpha_lpot(j),qs_control%eps_ppl,cval_lpot(i,j)))
            END DO
          END DO
        END IF
        IF(llsd) THEN
          ! lsd local potential
          CALL get_potential(potential=gth_potential,&
                             nexp_lsd = nexp_lsd,&
                             alpha_lsd = alpha_lsd,&
                             nct_lsd = nct_lsd,&
                             cval_lsd = cval_lsd)
          DO j=1,nexp_lsd
            DO i=1,nct_lsd(j)
              lppl = 2*(i-1)
              ppl_radius = MAX(ppl_radius,&
                           exp_radius(lppl,alpha_lsd(j),qs_control%eps_ppl,cval_lsd(i,j)))
            END DO
          END DO
        END IF

        ! Calculate the radii of the non-local part
        ! of the Goedecker pseudopotential (GTH)

        ppnl_radius = 0.0_dp

        p_ppnl    = (lppnl >= 0)

        DO l=0,lppnl
          DO iprj_ppnl=1,nprj_ppnl(l)
            lprj_ppnl = l + 2*(iprj_ppnl - 1)
            ppnl_radius = MAX(ppnl_radius,&
                              exp_radius(lprj_ppnl,alpha_ppnl(l),&
                                         qs_control%eps_pgf_orb,&
                                         cprj_ppnl(iprj_ppnl,l)))
          END DO
        END DO

        CALL set_potential(potential=gth_potential,&
                           core_charge_radius=core_charge_radius,&
                           ppl_radius=ppl_radius,&
                           ppnl_radius=ppnl_radius)

      END IF

      ! Calculate the aux fit orbital basis function radii

      IF (ASSOCIATED(aux_fit_basis_set)) THEN

        p_orb = .TRUE.
        CALL init_interaction_radii_orb_basis(aux_fit_basis_set,qs_control%eps_pgf_orb,&
             qs_control%eps_kg_orb,error=error)

      END IF

     ! Calculate the aux orbital basis function radii

      IF (ASSOCIATED(aux_basis_set)) THEN

        CALL init_interaction_radii_orb_basis(aux_basis_set,qs_control%eps_pgf_orb,error=error)

      END IF

     ! Calculate the geminal basis function radii

      IF (ASSOCIATED(geminal_basis_set)) THEN

        p_gem = .TRUE.
        CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,&
                                   nset=nset,&
                                   nshell=nshell,&
                                   npgf=npgf,&
                                   l=lshell,&
                                   ls=ls,&
                                   zet=zetg,&
                                   gcc=gcc,&
                                   pgf_radius=pgf_radius,&
                                   set_radius=set_radius)

        kind_radius = 0.0_dp

        DO iset=1,nset
          set_radius(iset) = 0.0_dp
          DO ipgf=1,npgf(iset)
            pgf_radius(ipgf,iset) = 0.0_dp
            DO ishell=1,nshell(iset)
              l = MAX(lshell(ishell,iset),ls(iset))
              gcca = gcc(ipgf,ishell,iset)
              zab = zetg(1,1,ipgf,iset) + zetg(2,2,ipgf,iset)
              zd  = (zetg(1,1,ipgf,iset) - zetg(2,2,ipgf,iset))**2 + 4._dp*zetg(1,2,ipgf,iset)
              zeta = 0.5_dp*MIN(zab+SQRT(zd),zab-SQRT(zd))
              pgf_radius(ipgf,iset) = MAX(pgf_radius(ipgf,iset),&
                                          exp_radius(l,zeta,qs_control%eps_pgf_orb,gcca))
            END DO
            set_radius(iset) = MAX(set_radius(iset),pgf_radius(ipgf,iset))
          END DO
          kind_radius = MAX(kind_radius,set_radius(iset))
        END DO

        CALL set_geminal_basis_set(geminal_basis_set=geminal_basis_set,&
                                   pgf_radius=pgf_radius,&
                                   set_radius=set_radius,&
                                   kind_radius=kind_radius)

      END IF

      ! Calculate the soft orbital basis function radii

      IF (ASSOCIATED(soft_basis)) THEN

        IF(paw_atom) THEN
          p_soft = .TRUE.
          CALL init_interaction_radii_orb_basis(soft_basis,qs_control%eps_pgf_orb,error=error)

        END IF
      END IF

      ! Calculate the lri basis function radii
      IF (ASSOCIATED(lri_basis)) THEN
          p_lri = .TRUE.
          CALL init_interaction_radii_orb_basis(lri_basis,qs_control%eps_pgf_orb,error=error)
      END IF

      ! Calculate the paw projector basis function radii

      IF (ASSOCIATED(paw_proj_set)) THEN

        p_paw = .TRUE.
        CALL get_paw_proj_set(paw_proj_set=paw_proj_set,&
                               nprj=nprj,&
                               maxl = maxl,&
                               zetprj=zet,&
                               rzetprj=rzetprj)
        rcprj = 0.0_dp
        rzetprj = 0.0_dp
        DO lprj = 0,maxl
           DO ip = 1,nprj(lprj)
             zeta = zet(ip,lprj)
             rzetprj(ip,lprj) = MAX(rzetprj(ip,lprj),&
                                      exp_radius(lprj,zeta,&
                                      qs_control%eps_pgf_orb,1.0_dp))
             rcprj = MAX(rcprj,rzetprj(ip,lprj))
           ENDDO
        ENDDO
        CALL set_paw_proj_set(paw_proj_set=paw_proj_set,&
                               rzetprj = rzetprj,&
                               rcprj=rcprj)

      END IF
    END DO  ! ikind

    CALL timestop(handle)

  END SUBROUTINE init_interaction_radii

! *****************************************************************************
!> \brief ...
!> \param orb_basis_set ...
!> \param eps_pgf_orb ...
!> \param eps_pgf_short ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE init_interaction_radii_orb_basis(orb_basis_set,eps_pgf_orb,eps_pgf_short,error)

    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    REAL(KIND=dp), INTENT(IN)                :: eps_pgf_orb
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: eps_pgf_short
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'init_interaction_radii_orb_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipgf, iset, ishell, l, nset
    INTEGER, DIMENSION(:), POINTER           :: npgf, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: lshell
    REAL(KIND=dp)                            :: eps_short, gcca, kind_radius, &
                                                short_radius, zeta
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pgf_radius, zet
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc

    IF (ASSOCIATED(orb_basis_set)) THEN

       IF(PRESENT(eps_pgf_short)) THEN
          eps_short = eps_pgf_short
       ELSE
          eps_short = eps_pgf_orb
       END IF
          
       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                              nset=nset,&
                              nshell=nshell,&
                              npgf=npgf,&
                              l=lshell,&
                              zet=zet,&
                              gcc=gcc,&
                              pgf_radius=pgf_radius,&
                              set_radius=set_radius)

       kind_radius = 0.0_dp
       short_radius = 0.0_dp

       DO iset=1,nset
         set_radius(iset) = 0.0_dp
         DO ipgf=1,npgf(iset)
           pgf_radius(ipgf,iset) = 0.0_dp
           DO ishell=1,nshell(iset)
             l = lshell(ishell,iset)
             gcca = gcc(ipgf,ishell,iset)
             zeta = zet(ipgf,iset)
             pgf_radius(ipgf,iset) = MAX(pgf_radius(ipgf,iset),&
                                         exp_radius(l,zeta,&
                                                    eps_pgf_orb,&
                                                    gcca))
             short_radius = MAX(short_radius,&
                             exp_radius(l,zeta,&
                             eps_short,&
                             gcca))
           END DO
           set_radius(iset) = MAX(set_radius(iset),pgf_radius(ipgf,iset))
         END DO
         kind_radius = MAX(kind_radius,set_radius(iset))
       END DO

       CALL set_gto_basis_set(gto_basis_set=orb_basis_set,&
                              pgf_radius=pgf_radius,&
                              set_radius=set_radius,&
                              kind_radius=kind_radius,&
                              short_kind_radius=short_radius)

     END IF

  END SUBROUTINE init_interaction_radii_orb_basis

! *****************************************************************************
!> \brief ...
!> \param se_control ...
!> \param atomic_kind_set ...
!> \param subsys_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE init_se_nlradius(se_control,atomic_kind_set,&
                   subsys_section,error)

    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_se_nlradius', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, nkind
    REAL(KIND=dp)                            :: kind_radius
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

    NULLIFY (orb_basis_set,atomic_kind)

    nkind = SIZE(atomic_kind_set)

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)

      IF (ASSOCIATED(orb_basis_set)) THEN

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               kind_radius=kind_radius)

        kind_radius = MAX (kind_radius,se_control%cutoff_exc)

        CALL set_gto_basis_set(gto_basis_set=orb_basis_set,&
                               kind_radius=kind_radius)

      END IF

    END DO

    CALL write_kind_radii(atomic_kind_set,subsys_section,error)

  END SUBROUTINE init_se_nlradius

! *****************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param subsys_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE write_kind_radii(atomic_kind_set,subsys_section,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=10)                        :: bas
    CHARACTER(LEN=default_string_length)     :: name, unit_str
    INTEGER                                  :: ikind, nkind, output_unit
    REAL(KIND=dp)                            :: conv, kind_radius
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    NULLIFY(orb_basis_set)
    bas="ORBITAL   "

    nkind = SIZE(atomic_kind_set)
!   *** Print the kind radii ***
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%RADII/KIND_RADII",extension=".Log",error=error)
    CALL section_vals_val_get(subsys_section,"PRINT%RADII%UNIT",c_val=unit_str,error=error)
    conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    IF (output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T56,A,T63,A,T75,A)")&
        "RADII: "//TRIM(bas)//" BASIS in "//TRIM(unit_str),&
        "Kind","Label","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,&
                             orb_basis_set=orb_basis_set)
        IF(ASSOCIATED(orb_basis_set)) THEN
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                 kind_radius=kind_radius)
          WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T69,F12.6)")&
            ikind,name,kind_radius*conv
        ELSE
           WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T72,A9)")&
              ikind,name,"no basis"
        END IF
      END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%RADII/KIND_RADII",error=error)

  END SUBROUTINE write_kind_radii

! *****************************************************************************
!> \brief  Write the radii of the core charge distributions to the output
!>         unit.
!> \param atomic_kind_set ...
!> \param subsys_section ...
!> \param error ...
!> \date    15.09.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_core_charge_radii(atomic_kind_set,subsys_section,error)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=default_string_length)     :: name, unit_str
    INTEGER                                  :: ikind, nkind, output_unit
    REAL(KIND=dp)                            :: conv, core_charge_radius
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gth_potential_type), POINTER        :: gth_potential

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%RADII/CORE_CHARGE_RADII",extension=".Log",error=error)
    CALL section_vals_val_get(subsys_section,"PRINT%RADII%UNIT",c_val=unit_str,error=error)
    conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    IF (output_unit>0) THEN
       nkind = SIZE(atomic_kind_set)
       WRITE (UNIT=output_unit,FMT="(/,T2,A,T56,A,T63,A,T75,A)")&
            "RADII: CORE CHARGE DISTRIBUTIONS in "//&
            TRIM(unit_str),"Kind","Label","Radius"
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name,&
               all_potential=all_potential,&
               gth_potential=gth_potential)
          IF (ASSOCIATED(all_potential)) THEN
             CALL get_potential(potential=all_potential,&
                  core_charge_radius=core_charge_radius)
             WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T69,F12.6)")&
                  ikind,name,core_charge_radius*conv
          ELSE IF (ASSOCIATED(gth_potential)) THEN
             CALL get_potential(potential=gth_potential,&
                  core_charge_radius=core_charge_radius)
             WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T69,F12.6)")&
                  ikind,name,core_charge_radius*conv
          END IF
       END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%RADII/CORE_CHARGE_RADII",error=error)
  END SUBROUTINE write_core_charge_radii

! *****************************************************************************
!> \brief   Write the orbital basis function radii to the output unit.
!> \param basis ...
!> \param atomic_kind_set ...
!> \param subsys_section ...
!> \param error ...
!> \date    05.06.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_pgf_orb_radii(basis,atomic_kind_set,subsys_section,error)

    CHARACTER(len=*)                         :: basis
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_pgf_orb_radii', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: bas
    CHARACTER(LEN=default_string_length)     :: name, unit_str
    INTEGER                                  :: ikind, ipgf, iset, nkind, &
                                                nset, output_unit
    INTEGER, DIMENSION(:), POINTER           :: npgf
    REAL(KIND=dp)                            :: conv, kind_radius, &
                                                short_kind_radius
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pgf_radius
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set, lri_basis_set, &
                                                orb_basis_set

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    NULLIFY(aux_basis_set,orb_basis_set,lri_basis_set)
    bas=" "
    bas(1:3) = basis(1:3)
    CALL uppercase ( bas )
    IF ( bas(1:3) == "ORB" ) THEN
      bas="ORBITAL   "
    ELSE IF ( bas(1:3) == "AUX" ) THEN
      bas="AUXILLIARY"
    ELSE IF ( bas(1:3) == "LRI" ) THEN
      bas="LOCAL RI"
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Undefined basis set type")
    END IF

    nkind = SIZE(atomic_kind_set)

!   *** Print the kind radii ***
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%RADII/KIND_RADII",extension=".Log",error=error)
    CALL section_vals_val_get(subsys_section,"PRINT%RADII%UNIT",c_val=unit_str,error=error)
    conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    IF (output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T46,A,T53,A,T63,A,T71,A)")&
        "RADII: "//TRIM(bas)//" BASIS in "//TRIM(unit_str),&
        "Kind","Label","Radius","OCE Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        IF ( bas(1:3) == "ORB" ) THEN
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name,&
                               orb_basis_set=orb_basis_set)
          IF(ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                   kind_radius=kind_radius,&
                                   short_kind_radius=short_kind_radius)
          END IF
        ELSE IF ( bas(1:3) == "AUX" ) THEN
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name,&
                               aux_basis_set=aux_basis_set)
          IF(ASSOCIATED(aux_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
                                   kind_radius=kind_radius)
          END IF
        ELSE IF ( bas(1:3) == "LOC" ) THEN
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name,&
                               lri_basis_set=lri_basis_set)
          IF(ASSOCIATED(lri_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=lri_basis_set,&
                                   kind_radius=kind_radius)
          END IF
        ELSE
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Undefined basis set type")
        END IF
        IF(ASSOCIATED(aux_basis_set).OR.ASSOCIATED(orb_basis_set)) THEN
          WRITE (UNIT=output_unit,FMT="(T45,I5,T53,A5,T57,F12.6,T69,F12.6)")&
            ikind,name,kind_radius*conv,short_kind_radius*conv
        ELSE
           WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T72,A9)")&
              ikind,name,"no basis"
        END IF
      END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%RADII/KIND_RADII",error=error)

!   *** Print the shell set radii ***
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%RADII/SET_RADII",extension=".Log",error=error)
    IF (output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T51,A,T57,A,T65,A,T75,A)")&
        "RADII: SHELL SETS OF "//TRIM(bas)//" BASIS in "//&
        TRIM(unit_str),"Kind","Label","Set","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        IF ( bas(1:3) == "ORB" ) THEN
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name,&
                               orb_basis_set=orb_basis_set)
          IF(ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                   nset=nset,&
                                   set_radius=set_radius)
          END IF
        ELSE IF ( bas(1:3) == "AUX" ) THEN
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name,&
                               aux_basis_set=aux_basis_set)
          IF(ASSOCIATED(aux_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
                                   nset=nset,&
                                   set_radius=set_radius)
          END IF
        ELSE IF ( bas(1:3) == "LOC" ) THEN
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name,&
                               lri_basis_set=lri_basis_set)
          IF(ASSOCIATED(lri_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=lri_basis_set,&
                                   nset=nset,&
                                   set_radius=set_radius)
          END IF
        ELSE
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Undefined basis set type")
        END IF
        IF(ASSOCIATED(aux_basis_set).OR.ASSOCIATED(orb_basis_set)) THEN
          WRITE (UNIT=output_unit,FMT="(T50,I5,T57,A5,(T63,I5,T69,F12.6))")&
            ikind,name,(iset,set_radius(iset)*conv,iset=1,nset)
        ELSE
           WRITE (UNIT=output_unit,FMT="(T50,I5,T58,A5,T73,A8)")&
              ikind,name,"no basis"
        END IF
      END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%RADII/SET_RADII",error=error)
!   *** Print the primitive Gaussian function radii ***
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%RADII/PGF_RADII",extension=".Log",error=error)
    IF (output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T51,A,T57,A,T65,A,T75,A)")&
        "RADII: PRIMITIVE GAUSSIANS OF "//TRIM(bas)//" BASIS in "//&
        TRIM(unit_str),"Kind","Label","Set","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        IF ( bas(1:3) == "ORB" ) THEN
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name,&
                               orb_basis_set=orb_basis_set)
          IF(ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                   nset=nset,&
                                   npgf=npgf,&
                                   pgf_radius=pgf_radius)
          END IF
        ELSE IF ( bas(1:3) == "AUX" ) THEN
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name,&
                               aux_basis_set=aux_basis_set)
          IF(ASSOCIATED(aux_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
                                   nset=nset,&
                                   npgf=npgf,&
                                   pgf_radius=pgf_radius)
          END IF
        ELSE IF ( bas(1:3) == "LOC" ) THEN
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name,&
                               lri_basis_set=lri_basis_set)
          IF(ASSOCIATED(lri_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=lri_basis_set,&
                                   nset=nset,&
                                   npgf=npgf,&
                                   pgf_radius=pgf_radius)
          END IF
        ELSE
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Undefined basis type")
        END IF
        IF(ASSOCIATED(aux_basis_set).OR.ASSOCIATED(orb_basis_set).OR. &
           ASSOCIATED(lri_basis_set)) THEN
          DO iset=1,nset
            WRITE (UNIT=output_unit,FMT="(T50,I5,T57,A5,T63,I5,(T69,F12.6))")&
              ikind,name,iset,&
              (pgf_radius(ipgf,iset)*conv,ipgf=1,npgf(iset))
          END DO
        ELSE
           WRITE (UNIT=output_unit,FMT="(T50,I5,T58,A5,T73,A8)")&
              ikind,name,"no basis"
        END IF
      END DO
   END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%RADII/PGF_RADII",error=error)
  END SUBROUTINE write_pgf_orb_radii

! *****************************************************************************
!> \brief   Write the geminal basis function radii to the output unit.
!> \param atomic_kind_set ...
!> \param subsys_section ...
!> \param error ...
!> \date    18.08.2009
!> \author  JGH
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_geminal_radii(atomic_kind_set,subsys_section,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=10)                        :: bas
    CHARACTER(LEN=default_string_length)     :: name, unit_str
    INTEGER                                  :: ikind, ipgf, iset, nkind, &
                                                nset, output_unit
    INTEGER, DIMENSION(:), POINTER           :: npgf
    REAL(KIND=dp)                            :: conv, kind_radius
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pgf_radius
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    NULLIFY(geminal_basis_set)

    nkind = SIZE(atomic_kind_set)

    ! Print the kind radii
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%RADII/KIND_RADII",extension=".Log",error=error)
    CALL section_vals_val_get(subsys_section,"PRINT%RADII%UNIT",c_val=unit_str,error=error)
    conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    IF (output_unit>0) THEN
      bas="GEMINAL"
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T55,A,T62,A,T74,A)")&
        "RADII: "//TRIM(bas)//" BASIS in "//TRIM(unit_str),"Kind","Label","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,&
                             geminal_basis_set=geminal_basis_set)
        IF(ASSOCIATED(geminal_basis_set)) THEN
          CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,&
                                 kind_radius=kind_radius)
          WRITE (UNIT=output_unit,FMT="(T54,I5,T62,A5,T69,F12.6)")&
            ikind,name,kind_radius*conv
        ELSE
           WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T72,A9)")&
              ikind,name,"no basis"
        END IF
      END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%RADII/KIND_RADII",error=error)

    ! Print the shell set radii
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%RADII/SET_RADII",extension=".Log",error=error)
    IF (output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T51,A,T57,A,T65,A,T75,A)")&
        "RADII: SHELL SETS OF "//TRIM(bas)//" BASIS in "//&
        TRIM(unit_str),"Kind","Label","Set","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,&
                             geminal_basis_set=geminal_basis_set)
        IF(ASSOCIATED(geminal_basis_set)) THEN
          CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,&
                                     nset=nset,&
                                     set_radius=set_radius)
          WRITE (UNIT=output_unit,FMT="(T50,I5,T57,A5,(T63,I5,T69,F12.6))")&
            ikind,name,(iset,set_radius(iset)*conv,iset=1,nset)
        ELSE
           WRITE (UNIT=output_unit,FMT="(T50,I5,T58,A5,T73,A8)")&
              ikind,name,"no basis"
        END IF
      END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%RADII/SET_RADII",error=error)
    ! Print the primitive Gaussian function radii
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%RADII/PGF_RADII",extension=".Log",error=error)
    IF (output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T51,A,T57,A,T65,A,T75,A)")&
        "RADII: PRIMITIVE GEMINALS OF "//TRIM(bas)//" BASIS in "//&
        TRIM(unit_str),"Kind","Label","Set","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,&
                             geminal_basis_set=geminal_basis_set)
        IF(ASSOCIATED(geminal_basis_set)) THEN
          CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,&
                                     nset=nset,&
                                     npgf=npgf,&
                                     pgf_radius=pgf_radius)
          DO iset=1,nset
            WRITE (UNIT=output_unit,FMT="(T50,I5,T57,A5,T63,I5,(T69,F12.6))")&
              ikind,name,iset,&
              (pgf_radius(ipgf,iset)*conv,ipgf=1,npgf(iset))
          END DO
        ELSE
           WRITE (UNIT=output_unit,FMT="(T50,I5,T58,A5,T73,A8)")&
              ikind,name,"no basis"
        END IF
      END DO
   END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%RADII/PGF_RADII",error=error)
  END SUBROUTINE write_geminal_radii

! *****************************************************************************
!> \brief   Write the radii of the exponential functions of the Goedecker
!>           pseudopotential (GTH, local part) to the logical unit number
!>          "output_unit".
!> \param atomic_kind_set ...
!> \param subsys_section ...
!> \param error ...
!> \date    06.10.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_ppl_radii(atomic_kind_set,subsys_section,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=default_string_length)     :: name, unit_str
    INTEGER                                  :: ikind, nkind, output_unit
    REAL(KIND=dp)                            :: conv, ppl_radius
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gth_potential_type), POINTER        :: gth_potential

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%RADII/GTH_PPL_RADII",extension=".Log",error=error)
    CALL section_vals_val_get(subsys_section,"PRINT%RADII%UNIT",c_val=unit_str,error=error)
    conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    IF (output_unit>0) THEN
      nkind = SIZE(atomic_kind_set)
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T56,A,T63,A,T75,A)")&
        "RADII: LOCAL PART OF GTH/ELP PP in "//&
        TRIM(unit_str),"Kind","Label","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,&
                             gth_potential=gth_potential)
        IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             ppl_radius=ppl_radius)
          WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T69,F12.6)")&
            ikind,name,ppl_radius*conv
        END IF
      END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%RADII/GTH_PPL_RADII",error=error)
  END SUBROUTINE write_ppl_radii

! *****************************************************************************
!> \brief  Write the radii of the projector functions of the Goedecker
!>          pseudopotential (GTH, non-local part) to the logical unit number
!>          "output_unit".
!> \param atomic_kind_set ...
!> \param subsys_section ...
!> \param error ...
!> \date    06.10.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_ppnl_radii(atomic_kind_set,subsys_section,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=default_string_length)     :: name, unit_str
    INTEGER                                  :: ikind, nkind, output_unit
    REAL(KIND=dp)                            :: conv, ppnl_radius
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gth_potential_type), POINTER        :: gth_potential

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%RADII/GTH_PPNL_RADII",extension=".Log",error=error)
    CALL section_vals_val_get(subsys_section,"PRINT%RADII%UNIT",c_val=unit_str,error=error)
    conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    IF (output_unit>0) THEN
      nkind = SIZE(atomic_kind_set)
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T56,A,T63,A,T75,A)")&
        "RADII: NON-LOCAL PART OF GTH PP in "//&
        TRIM(unit_str),"Kind","Label","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,&
                             gth_potential=gth_potential)
        IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             ppnl_radius=ppnl_radius)
          WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T69,F12.6)")&
            ikind,name,ppnl_radius*conv
        END IF
      END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%RADII/GTH_PPNL_RADII",error=error)
  END SUBROUTINE write_ppnl_radii

! *****************************************************************************
!> \brief   Write the radii of the one center projector
!> \param atomic_kind_set ...
!> \param subsys_section ...
!> \param error ...
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_paw_radii(atomic_kind_set,subsys_section,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=default_string_length)     :: name, unit_str
    INTEGER                                  :: ikind, nkind, output_unit
    LOGICAL                                  :: paw_atom
    REAL(KIND=dp)                            :: conv, rcprj
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(paw_proj_set_type), POINTER         :: paw_proj_set

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%RADII/GAPW_PRJ_RADII",extension=".Log",error=error)
    CALL section_vals_val_get(subsys_section,"PRINT%RADII%UNIT",c_val=unit_str,error=error)
    conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    IF (output_unit>0) THEN
       nkind = SIZE(atomic_kind_set)
       WRITE (UNIT=output_unit,FMT="(/,T2,A,T56,A,T63,A,T75,A)")&
            "RADII: ONE CENTER PROJECTORS in "//&
            TRIM(unit_str),"Kind","Label","Radius"
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name,&
               paw_atom=paw_atom,&
               paw_proj_set=paw_proj_set)
          IF (paw_atom .AND. ASSOCIATED(paw_proj_set)) THEN
             CALL get_paw_proj_set(paw_proj_set=paw_proj_set,&
                  rcprj=rcprj)
             WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T69,F12.6)")&
                  ikind,name,rcprj*conv
          END IF
       END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%RADII/GAPW_PRJ_RADII",error=error)
  END SUBROUTINE write_paw_radii

END MODULE qs_interactions
