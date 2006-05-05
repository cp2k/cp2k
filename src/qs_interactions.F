!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_interactions [1.0] *
!!
!!   NAME
!!     qs_interactions
!!
!!   FUNCTION
!!     Calculate the interaction radii for the operator matrix calculation.
!!
!!   AUTHOR
!!     MK (12.07.2000)
!!
!!   MODIFICATION HISTORY
!!     Joost VandeVondele : added exp_radius_very_extended
!!     24.09.2002 overloaded init_interaction_radii for KG use (gt)
!!     27.02.2004 integrated KG into QS version of init_interaction_radii (jgh)
!!     07.03.2006 new routine to extend kind interaction radius for semi-empiric (jgh)
!!
!!***
!******************************************************************************

MODULE qs_interactions

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type,&
                                             set_gto_basis_set
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell
  USE cp_control_types,                ONLY: qs_control_type,&
                                             semi_empirical_control_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE external_potential_types,        ONLY: all_potential_type,&
                                             elp_potential_type,&
                                             get_potential,&
                                             gth_potential_type,&
                                             set_potential
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE orbital_pointers,                ONLY: coset
  USE paw_proj_set_types,              ONLY: get_paw_proj_set,&
                                             paw_proj_set_type,&
                                             set_paw_proj_set
  USE qs_all_potential,                ONLY: radii_3c_orbxall
  USE qs_util,                         ONLY: exp_radius
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_interactions'


! *** Public subroutines ***

  PUBLIC :: init_interaction_radii, init_se_nlradius

! *** Public functions ***

  PUBLIC :: exp_radius_very_extended

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE init_interaction_radii(qs_control,cell,atomic_kind_set,force_env_section,error)

!   Purpose: Initialize all the atomic kind radii for a given threshold value.

!   History: - Creation (24.09.2002, gt)

!   ***************************************************************************

    TYPE(qs_control_type), INTENT(IN)        :: qs_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_interaction_radii'

    INTEGER :: i, iexp_ppl, ikind, ip, ipgf, iprj_ppnl, iset, ishell, istat, &
      l, lppl, lppnl, lprj, lprj_ppnl, maxl, maxnset, nexp_ppl, nkind, nset
    INTEGER, DIMENSION(:), POINTER           :: nexp_ppl_elp, npgf, nprj, &
                                                nprj_ppnl, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: lshell
    LOGICAL                                  :: p_aux, p_ccr, p_orb, p_paw, &
                                                p_ppl, p_ppnl, p_soft, &
                                                paw_atom
    REAL(KIND=dp) :: alpha_core_charge, alpha_ppl, ccore_charge, cerf_ppl, &
      core_charge_radius, eps_reduce, gcca, kind_radius, ppl_radius, &
      ppnl_radius, rcprj,  short_radius, zeta
    REAL(KIND=dp), DIMENSION(:), POINTER     :: alpha_ppl_elp, alpha_ppnl, &
                                                cexp_ppl, set_radius
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cexp_ppl_elp, cprj_ppnl, &
                                                pgf_radius, rzetprj, zet
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set, orb_basis_set, &
                                                soft_basis
    TYPE(paw_proj_set_type), POINTER         :: paw_proj_set

!   ---------------------------------------------------------------------------

    NULLIFY ( atomic_kind, all_potential, gth_potential, elp_potential, &
              orb_basis_set, aux_basis_set, soft_basis, paw_proj_set)
    NULLIFY ( npgf, nprj_ppnl, nshell, nexp_ppl_elp, lshell, nprj )
    NULLIFY ( alpha_ppnl, cexp_ppl, set_radius, alpha_ppl_elp, cprj_ppnl, &
             pgf_radius, zet, cexp_ppl_elp, gcc)

    p_orb  = .FALSE.
    p_aux  = .FALSE.
    p_ccr  = .FALSE.
    p_ppl  = .FALSE.
    p_ppnl = .FALSE.
    p_soft = .FALSE.
    p_paw  = .FALSE.

    eps_reduce = qs_control%gapw_control%eps_3c_reduce
    nkind = SIZE(atomic_kind_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxnset=maxnset)

    nexp_ppl = 0

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           all_potential=all_potential,&
                           gth_potential=gth_potential,&
                           elp_potential=elp_potential,&
                           aux_basis_set=aux_basis_set,&
                           orb_basis_set=orb_basis_set,&
                           soft_basis_set=soft_basis,&
                           paw_proj_set=paw_proj_set,&
                           paw_atom=paw_atom)

!     *** Calculate the orbital basis function radii ***
!         For ALL electron this has to come before the calculation of the 
!         radii used to calculate the 3c lists. 
!         Indeed, there the radii of the GTO primitives are needed
      IF (ASSOCIATED(orb_basis_set)) THEN

        p_orb = .TRUE.
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
                                                     qs_control%eps_pgf_orb,&
                                                     gcca))
              short_radius = MAX(short_radius,&
                             exp_radius(l,zeta,&
                             qs_control%eps_pgf_orb*eps_reduce,&
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


      IF (ASSOCIATED(all_potential)) THEN

        p_ccr = .TRUE.
        CALL get_potential(potential=all_potential,&
                           alpha_core_charge=alpha_core_charge,&
                           ccore_charge=ccore_charge)

!       *** Calculate the radius of the core charge distribution ***

        core_charge_radius = exp_radius(0,alpha_core_charge,&
                                        qs_control%eps_core_charge,&
                                        ccore_charge)

        CALL set_potential(potential=all_potential,&
                           core_charge_radius=core_charge_radius)

        ALLOCATE(all_potential%orb_radius_3c(nkind),STAT=istat)
        ALLOCATE(all_potential%set_radius_3c(maxnset,nkind),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "orb_radius_3c,set_radius_3c")

        CALL radii_3c_orbxall(atomic_kind_set,all_potential,qs_control%eps_core_charge)

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
                           lppnl=lppnl,&
                           alpha_ppnl=alpha_ppnl,&
                           nprj_ppnl=nprj_ppnl,&
                           cprj_ppnl=cprj_ppnl)

!       *** Calculate the radius of the core charge distribution ***

        core_charge_radius = exp_radius(0,alpha_core_charge,&
                                        qs_control%eps_core_charge,&
                                        ccore_charge)

!       *** Calculate the radii of the local part  ***
!       *** of the Goedecker pseudopotential (GTH) ***

        ppl_radius = exp_radius(0,alpha_ppl,qs_control%eps_ppl,cerf_ppl)

        DO iexp_ppl=1,nexp_ppl
          lppl = 2*(iexp_ppl - 1)
          ppl_radius = MAX(ppl_radius,&
                           exp_radius(lppl,alpha_ppl,&
                                      qs_control%eps_ppl,&
                                      cexp_ppl(iexp_ppl)))
        END DO

!       *** Calculate the radii of the non-local part ***
!       *** of the Goedecker pseudopotential (GTH)    ***

        ppnl_radius = 0.0_dp

        p_ppnl    = (lppnl >= 0)

        DO l=0,lppnl
          DO iprj_ppnl=1,nprj_ppnl(l)
            lprj_ppnl = l + 2*(iprj_ppnl - 1)
            ppnl_radius = MAX(ppnl_radius,&
                              exp_radius(lprj_ppnl,alpha_ppnl(l),&
                                         qs_control%eps_ppnl,&
                                         cprj_ppnl(iprj_ppnl,l)))
          END DO
        END DO

        CALL set_potential(potential=gth_potential,&
                           core_charge_radius=core_charge_radius,&
                           ppl_radius=ppl_radius,&
                           ppnl_radius=ppnl_radius)

      ELSE IF (ASSOCIATED(elp_potential)) THEN

        p_ccr = .TRUE.
        p_ppl = .TRUE.
        CALL get_potential(potential=elp_potential,&
                           alpha_core_charge=alpha_core_charge,&
                           ccore_charge=ccore_charge,&
                           alpha_ppl=alpha_ppl_elp,&
                           cerf_ppl=cerf_ppl,&
                           nexp_ppl=nexp_ppl_elp,&
                           cexp_ppl=cexp_ppl_elp)

!       *** Calculate the radius of the core charge distribution ***

        core_charge_radius = exp_radius(0,alpha_core_charge,&
                                        qs_control%eps_core_charge,&
                                        ccore_charge)

!       *** Calculate the radii of the  ***
!       *** effective local pseudopotential (elp) ***

        ppl_radius = exp_radius(0,alpha_ppl_elp(1),qs_control%eps_ppl,cerf_ppl)

        DO i=1,SIZE(nexp_ppl_elp)
          DO iexp_ppl=1,nexp_ppl_elp(i)
            lppl = 2*(iexp_ppl - 1)
            ppl_radius = MAX(ppl_radius,&
                           exp_radius(lppl,alpha_ppl_elp(i),&
                                      qs_control%eps_core_charge,&
                                      cexp_ppl_elp(i,iexp_ppl)))
          END DO
        END DO

        CALL set_potential(potential=elp_potential,&
                           core_charge_radius=core_charge_radius,&
                           ppl_radius=ppl_radius)

      END IF


!     *** Calculate the aux orbital basis function radii ***

      IF (ASSOCIATED(aux_basis_set)) THEN

        CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
                               nset=nset,&
                               nshell=nshell,&
                               npgf=npgf,&
                               l=lshell,&
                               zet=zet,&
                               gcc=gcc,&
                               pgf_radius=pgf_radius,&
                               set_radius=set_radius)

        kind_radius = 0.0_dp

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
                                                     qs_control%eps_pgf_orb,&
                                                     gcca))
            END DO
            set_radius(iset) = MAX(set_radius(iset),pgf_radius(ipgf,iset))
          END DO
          kind_radius = MAX(kind_radius,set_radius(iset))
        END DO

        CALL set_gto_basis_set(gto_basis_set=aux_basis_set,&
                               pgf_radius=pgf_radius,&
                               set_radius=set_radius,&
                               kind_radius=kind_radius)

      END IF

!     *** Calculate the soft orbital basis function radii ***

      IF (ASSOCIATED(soft_basis)) THEN

        IF(paw_atom) THEN
          p_soft = .TRUE.
          CALL get_gto_basis_set(gto_basis_set=soft_basis,&
                               nset=nset,&
                               nshell=nshell,&
                               npgf=npgf,&
                               l=lshell,&
                               zet=zet,&
                               gcc=gcc,&
                               pgf_radius=pgf_radius,&
                               set_radius=set_radius)

          kind_radius = 0.0_dp

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
                                                     qs_control%eps_pgf_orb,&
                                                     gcca))
              END DO
              set_radius(iset) = MAX(set_radius(iset),pgf_radius(ipgf,iset))
            END DO
            kind_radius = MAX(kind_radius,set_radius(iset))
          END DO

          CALL set_gto_basis_set(gto_basis_set=soft_basis,&
                               pgf_radius=pgf_radius,&
                               set_radius=set_radius,&
                               kind_radius=kind_radius)

        END IF
      END IF

!     *** Calculate the paw projector basis function radii ***

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

    IF (p_orb)  CALL write_pgf_orb_radii("orb",cell,atomic_kind_set,force_env_section,error)
    IF (p_aux)  CALL write_pgf_orb_radii("aux",cell,atomic_kind_set,force_env_section,error)
    IF (p_ccr)  CALL write_core_charge_radii(cell,atomic_kind_set,force_env_section,error)
    IF (p_ppl)  CALL write_ppl_radii(cell,atomic_kind_set,force_env_section,error)
    IF (p_ppnl) CALL write_ppnl_radii(cell,atomic_kind_set,force_env_section,error)
    IF (p_paw)  CALL write_paw_radii(cell,atomic_kind_set,force_env_section,error)

  END SUBROUTINE init_interaction_radii

! *****************************************************************************
  SUBROUTINE init_se_nlradius(se_control,cell,atomic_kind_set,&
                   force_env_section,error)

    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'init_se_nlradius'

    INTEGER                                  :: ikind, nkind
    REAL(KIND=dp)                            :: kind_radius
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

! ------------------------------------------------------------------------------

    NULLIFY (orb_basis_set,atomic_kind)

    nkind = SIZE(atomic_kind_set)

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)

      IF (ASSOCIATED(orb_basis_set)) THEN

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               kind_radius=kind_radius)

        kind_radius = MAX (kind_radius,se_control%rc_interaction)

        CALL set_gto_basis_set(gto_basis_set=orb_basis_set,&
                               kind_radius=kind_radius)

      END IF

    END DO

    CALL write_kind_radii(cell,atomic_kind_set,force_env_section,error)

  END SUBROUTINE init_se_nlradius

! *****************************************************************************

  SUBROUTINE write_kind_radii(cell,atomic_kind_set,force_env_section,error)

    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=10)                        :: bas
    CHARACTER(LEN=8)                         :: unit_of_length_name
    CHARACTER(LEN=default_string_length)     :: name
    INTEGER                                  :: ikind, nkind, output_unit
    REAL(KIND=dp)                            :: kind_radius, unit_of_length
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    NULLIFY(orb_basis_set)
    bas="ORBITAL   "

    CALL get_cell(cell=cell,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

    nkind = SIZE(atomic_kind_set)

!   *** Print the kind radii ***
    output_unit = cp_print_key_unit_nr(logger,force_env_section,&
         "SUBSYS%PRINT%RADII/KIND_RADII",extension=".Log",error=error)
    IF (output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T56,A,T63,A,T75,A)")&
        "RADII: "//TRIM(bas)//" BASIS in "//TRIM(unit_of_length_name),&
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
            ikind,name,kind_radius/unit_of_length
        ELSE
           WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T72,A9)")&
              ikind,name,"no basis"
        END IF
      END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
         "SUBSYS%PRINT%RADII/KIND_RADII",error=error)

  END SUBROUTINE write_kind_radii

! *****************************************************************************

  SUBROUTINE write_core_charge_radii(cell,atomic_kind_set,force_env_section,error)

!   Purpose: Write the radii of the core charge distributions to the output
!            unit.

!   History: - Creation (15.09.2000,MK)

!   ***************************************************************************

    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=8)                         :: unit_of_length_name
    CHARACTER(LEN=default_string_length)     :: name
    INTEGER                                  :: ikind, nkind, output_unit
    REAL(KIND=dp)                            :: core_charge_radius, &
                                                unit_of_length
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(gth_potential_type), POINTER        :: gth_potential

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,force_env_section,&
         "SUBSYS%PRINT%RADII/CORE_CHARGE_RADII",extension=".Log",error=error)
    IF (output_unit>0) THEN
       CALL get_cell(cell=cell,&
                     unit_of_length=unit_of_length,&
                     unit_of_length_name=unit_of_length_name)       
       nkind = SIZE(atomic_kind_set)
       WRITE (UNIT=output_unit,FMT="(/,T2,A,T56,A,T63,A,T75,A)")&
            "RADII: CORE CHARGE DISTRIBUTIONS in "//&
            TRIM(unit_of_length_name),"Kind","Label","Radius"
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name,&
               all_potential=all_potential,&
               elp_potential=elp_potential,&
               gth_potential=gth_potential)
          IF (ASSOCIATED(all_potential)) THEN
             CALL get_potential(potential=all_potential,&
                  core_charge_radius=core_charge_radius)
             WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T69,F12.6)")&
                  ikind,name,core_charge_radius/unit_of_length
          ELSE IF (ASSOCIATED(elp_potential)) THEN
             CALL get_potential(potential=elp_potential,&
                  core_charge_radius=core_charge_radius)
             WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T69,F12.6)")&
                  ikind,name,core_charge_radius/unit_of_length
          ELSE IF (ASSOCIATED(gth_potential)) THEN
             CALL get_potential(potential=gth_potential,&
                  core_charge_radius=core_charge_radius)
             WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T69,F12.6)")&
                  ikind,name,core_charge_radius/unit_of_length
          END IF
       END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
         "SUBSYS%PRINT%RADII/CORE_CHARGE_RADII",error=error)
  END SUBROUTINE write_core_charge_radii

! *****************************************************************************

  SUBROUTINE write_pgf_orb_radii(basis,cell,atomic_kind_set,force_env_section,error)

!   Purpose: Write the orbital basis function radii to the output unit.

!   History: - Creation (05.06.2000,MK)

!   ***************************************************************************

    CHARACTER(len=*)                         :: basis
    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=10)                        :: bas
    CHARACTER(LEN=8)                         :: unit_of_length_name
    CHARACTER(LEN=default_string_length)     :: name
    INTEGER                                  :: ikind, ipgf, iset, nkind, &
                                                nset, output_unit
    INTEGER, DIMENSION(:), POINTER           :: npgf
    REAL(KIND=dp)                            :: kind_radius, short_kind_radius,&
                                                unit_of_length
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pgf_radius
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set, orb_basis_set

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    NULLIFY(aux_basis_set,orb_basis_set)
    bas=" "
    bas(1:3) = basis(1:3)
    CALL uppercase ( bas )
    IF ( bas(1:3) == "ORB" ) THEN
      bas="ORBITAL   "
    ELSE IF ( bas(1:3) == "AUX" ) THEN
      bas="AUXILLIARY"
    ELSE 
      CALL stop_program("write_pgf_orb_radii","undefined basis type")
    END IF

    CALL get_cell(cell=cell,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

    nkind = SIZE(atomic_kind_set)

!   *** Print the kind radii ***
    output_unit = cp_print_key_unit_nr(logger,force_env_section,&
         "SUBSYS%PRINT%RADII/KIND_RADII",extension=".Log",error=error)
    IF (output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T46,A,T53,A,T63,A,T75,A)")&
        "RADII: "//TRIM(bas)//" BASIS in "//TRIM(unit_of_length_name),&
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
        ELSE
          CALL stop_program("write_pgf_orb_radii","undefined basis type")
        END IF
        IF(ASSOCIATED(aux_basis_set).OR.ASSOCIATED(orb_basis_set)) THEN
          WRITE (UNIT=output_unit,FMT="(T45,I5,T53,A5,T59,F12.6,T73,F12.6)")&
            ikind,name,kind_radius/unit_of_length,short_kind_radius/unit_of_length
        ELSE
           WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T72,A9)")&
              ikind,name,"no basis"
        END IF
      END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
         "SUBSYS%PRINT%RADII/KIND_RADII",error=error)

!   *** Print the shell set radii ***
    output_unit = cp_print_key_unit_nr(logger,force_env_section,&
         "SUBSYS%PRINT%RADII/SET_RADII",extension=".Log",error=error)
    IF (output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T51,A,T57,A,T65,A,T75,A)")&
        "RADII: SHELL SETS OF "//TRIM(bas)//" BASIS in "//&
        TRIM(unit_of_length_name),"Kind","Label","Set","Radius"
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
        ELSE
          CALL stop_program("write_pgf_orb_radii","undefined basis type")
        END IF
        IF(ASSOCIATED(aux_basis_set).OR.ASSOCIATED(orb_basis_set)) THEN
          WRITE (UNIT=output_unit,FMT="(T50,I5,T57,A5,(T63,I5,T69,F12.6))")&
            ikind,name,(iset,set_radius(iset)/unit_of_length,iset=1,nset)
        ELSE
           WRITE (UNIT=output_unit,FMT="(T50,I5,T58,A5,T73,A8)")&
              ikind,name,"no basis"
        END IF
      END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
         "SUBSYS%PRINT%RADII/SET_RADII",error=error)
!   *** Print the primitive Gaussian function radii ***
    output_unit = cp_print_key_unit_nr(logger,force_env_section,&
         "SUBSYS%PRINT%RADII/PGF_RADII",extension=".Log",error=error)
    IF (output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T51,A,T57,A,T65,A,T75,A)")&
        "RADII: PRIMITIVE GAUSSIANS OF "//TRIM(bas)//" BASIS in "//&
        TRIM(unit_of_length_name),"Kind","Label","Set","Radius"
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
        ELSE
          CALL stop_program("write_pgf_orb_radii","undefined basis type")
        END IF
        IF(ASSOCIATED(aux_basis_set).OR.ASSOCIATED(orb_basis_set)) THEN
          DO iset=1,nset
            WRITE (UNIT=output_unit,FMT="(T50,I5,T57,A5,T63,I5,(T69,F12.6))")&
              ikind,name,iset,&
              (pgf_radius(ipgf,iset)/unit_of_length,ipgf=1,npgf(iset))
          END DO
        ELSE
           WRITE (UNIT=output_unit,FMT="(T50,I5,T58,A5,T73,A8)")&
              ikind,name,"no basis"
        END IF
      END DO
   END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
         "SUBSYS%PRINT%RADII/PGF_RADII",error=error)
  END SUBROUTINE write_pgf_orb_radii

! *****************************************************************************

  SUBROUTINE write_ppl_radii(cell,atomic_kind_set,force_env_section,error)

!   Purpose: Write the radii of the exponential functions of the Goedecker
!            pseudopotential (GTH, local part) to the logical unit number
!            "output_unit".

!   History: - Creation (06.10.2000,MK)

!   ***************************************************************************

    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=8)                         :: unit_of_length_name
    CHARACTER(LEN=default_string_length)     :: name
    INTEGER                                  :: ikind, nkind, output_unit
    REAL(KIND=dp)                            :: ppl_radius, unit_of_length
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(gth_potential_type), POINTER        :: gth_potential

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,force_env_section,&
         "SUBSYS%PRINT%RADII/GTH_PPL_RADII",extension=".Log",error=error)
    IF (output_unit>0) THEN
      CALL get_cell(cell=cell,&
                    unit_of_length=unit_of_length,&
                    unit_of_length_name=unit_of_length_name)
      nkind = SIZE(atomic_kind_set)
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T56,A,T63,A,T75,A)")&
        "RADII: LOCAL PART OF GTH/ELP PP in "//&
        TRIM(unit_of_length_name),"Kind","Label","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,&
                             elp_potential=elp_potential,&
                             gth_potential=gth_potential)
        IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             ppl_radius=ppl_radius)
          WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T69,F12.6)")&
            ikind,name,ppl_radius/unit_of_length
        ELSE IF (ASSOCIATED(elp_potential)) THEN
          CALL get_potential(potential=elp_potential,&
                             ppl_radius=ppl_radius)
          WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T69,F12.6)")&
            ikind,name,ppl_radius/unit_of_length
        END IF
      END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
         "SUBSYS%PRINT%RADII/GTH_PPL_RADII",error=error)
  END SUBROUTINE write_ppl_radii

! *****************************************************************************

  SUBROUTINE write_ppnl_radii(cell,atomic_kind_set,force_env_section,error)

!   Purpose: Write the radii of the projector functions of the Goedecker
!            pseudopotential (GTH, non-local part) to the logical unit number
!            "output_unit".

!   History: - Creation (06.10.2000,MK)

!   ***************************************************************************

    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=8)                         :: unit_of_length_name
    CHARACTER(LEN=default_string_length)     :: name
    INTEGER                                  :: ikind, nkind, output_unit
    REAL(KIND=dp)                            :: ppnl_radius, unit_of_length
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gth_potential_type), POINTER        :: gth_potential

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,force_env_section,&
         "SUBSYS%PRINT%RADII/GTH_PPNL_RADII",extension=".Log",error=error)
    IF (output_unit>0) THEN
      CALL get_cell(cell=cell,&
                    unit_of_length=unit_of_length,&
                    unit_of_length_name=unit_of_length_name)
      nkind = SIZE(atomic_kind_set)
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T56,A,T63,A,T75,A)")&
        "RADII: NON-LOCAL PART OF GTH PP in "//&
        TRIM(unit_of_length_name),"Kind","Label","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,&
                             gth_potential=gth_potential)
        IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             ppnl_radius=ppnl_radius)
          WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T69,F12.6)")&
            ikind,name,ppnl_radius/unit_of_length
        END IF
      END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
         "SUBSYS%PRINT%RADII/GTH_PPNL_RADII",error=error)
  END SUBROUTINE write_ppnl_radii

! *****************************************************************************
  SUBROUTINE write_paw_radii(cell,atomic_kind_set,force_env_section,error)
  
!   Purpose: Write the radii of the one center projector 

!   ***************************************************************************

    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=8)                         :: unit_of_length_name
    CHARACTER(LEN=default_string_length)     :: name
    INTEGER                                  :: ikind, nkind, output_unit
    LOGICAL                                  :: paw_atom
    REAL(KIND=dp)                            :: rcprj, unit_of_length
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(paw_proj_set_type), POINTER         :: paw_proj_set

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,force_env_section,&
         "SUBSYS%PRINT%RADII/GAPW_PRJ_RADII",extension=".Log",error=error)
    IF (output_unit>0) THEN
       CALL get_cell(cell=cell,&
            unit_of_length=unit_of_length,&
            unit_of_length_name=unit_of_length_name)
       nkind = SIZE(atomic_kind_set)
       WRITE (UNIT=output_unit,FMT="(/,T2,A,T56,A,T63,A,T75,A)")&
            "RADII: ONE CENTER PROJECTORS in "//&
            TRIM(unit_of_length_name),"Kind","Label","Radius"
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=name,&
               paw_atom=paw_atom,&
               paw_proj_set=paw_proj_set)
          IF (paw_atom) THEN
             CALL get_paw_proj_set(paw_proj_set=paw_proj_set,&
                  rcprj=rcprj)
             WRITE (UNIT=output_unit,FMT="(T55,I5,T63,A5,T69,F12.6)")&
                  ikind,name,rcprj/unit_of_length
          END IF
       END DO
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
         "SUBSYS%PRINT%RADII/GAPW_PRJ_RADII",error=error)
  END SUBROUTINE write_paw_radii

! *****************************************************************************

! provides a radius that is keeping into account all the angular momentum
! components given the product of the two gaussians, and their original
! centers
!******
  FUNCTION exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,o1,o2,ra,rb,rp,&
                          zetp,eps,prefactor,cutoff) RESULT(radius)

    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, lb_max
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: pab
    INTEGER, OPTIONAL                        :: o1, o2
    REAL(KIND=dp), INTENT(IN)                :: ra(3), rb(3), rp(3), zetp, &
                                                eps, prefactor, cutoff
    REAL(KIND=dp)                            :: radius

    INTEGER                                  :: d, i, ico, j, jco, k, la(3), &
                                                lb(3), lxa, lxb, lya, lyb, &
                                                lza, lzb
    REAL(KIND=dp)                            :: bini, binj, coef(0:20,3), &
                                                polycoef(0:60), &
                                                prefactor_local, s1, s2

    radius=0.0_dp
    polycoef(:)=0.0_dp
    DO lxa=0,la_max
    DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
          DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
          DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
             la=(/lxa,lya,lza/)
             lb=(/lxb,lyb,lzb/)
             ico=coset(lxa,lya,lza)
             jco=coset(lxb,lyb,lzb)
             IF (PRESENT(pab)) THEN
               prefactor_local=prefactor*MAX(ABS(pab(o1+ico,o2+jco)),cutoff)
             ELSE
               prefactor_local=prefactor*MAX(1.0_dp,cutoff)
             ENDIF
             coef(:,:)=0.0_dp
             DO d=1,3
               bini=1.0_dp
               s1=1.0_dp
               DO i=0,la(d)
                  binj=1.0_dp
                  s2=1.0_dp
                  DO j=0,lb(d)
                    coef(la(d)+lb(d)-i-j,d)=coef(la(d)+lb(d)-i-j,d) +&
                                            bini*binj*s1*s2
                    binj=(binj*(lb(d)-j))/(j+1)
                    s2=s2*(rb(d)-rp(d))
                  ENDDO
                  bini=(bini*(la(d)-i))/(i+1)
                  s1=s1*(ra(d)-rp(d))
                ENDDO
             ENDDO
             DO i=0,la(1)+lb(1)
             DO j=0,la(2)+lb(2)
             DO k=0,la(3)+lb(3)
                polycoef(i+j+k)=MAX(polycoef(i+j+k),&
                                    ABS(prefactor_local*&
                                        coef(i,1)*coef(j,2)*coef(k,3)))
             ENDDO
             ENDDO
             ENDDO
  
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    DO i=0,la_max+lb_max
          radius=MAX(radius,exp_radius(i,zetp,eps,polycoef(i),1.0E-2_dp) )
    ENDDO

  END FUNCTION exp_radius_very_extended

! *****************************************************************************

END MODULE qs_interactions
