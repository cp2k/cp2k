!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
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
!!
!!   SOURCE
!******************************************************************************

MODULE qs_interactions

  USE kinds, ONLY: dp

  USE atomic_kind_types,        ONLY: atomic_kind_type,&
                                      get_atomic_kind
  USE basis_set_types,          ONLY: gto_basis_set_type,&
                                      get_gto_basis_set,&
                                      set_gto_basis_set
  USE dft_types,                ONLY: kg_control_type,&
                                      qs_control_type
  USE external_potential_types, ONLY: all_potential_type,&
                                      elp_potential_type,&
                                      get_potential,&
                                      gth_potential_type,&
                                      set_potential
  USE global_types,             ONLY: global_environment_type
  USE orbital_pointers,         ONLY: coset,ncoset
  USE simulation_cell,          ONLY: cell_type,&
                                      get_cell
  USE termination,              ONLY: stop_program



  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: init_interaction_radii

! *** Public functions ***

  PUBLIC :: exp_radius,&
            exp_radius_very_extended

! *****************************************************************************

  INTERFACE init_interaction_radii
    MODULE PROCEDURE init_interaction_radii_kg,init_interaction_radii_qs
  END INTERFACE

! *****************************************************************************
!!***

CONTAINS

! *****************************************************************************

  SUBROUTINE init_interaction_radii_kg(kg_control,cell,atomic_kind_set,globenv)

!   Purpose: Initialize all the atomic kind radii for a given threshold value.
!            (kim-gordon simulations)

!   History: - Creation (24.09.2002, gt)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(cell_type), POINTER                      :: cell
    TYPE(global_environment_type), INTENT(IN)     :: globenv
    TYPE(kg_control_type)                         :: kg_control

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(elp_potential_type), POINTER :: elp_potential
    TYPE(gth_potential_type), POINTER :: gth_potential
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    TYPE(gto_basis_set_type), POINTER :: aux_basis_set

    REAL(KIND = dp) :: alpha_core_charge,alpha_ppl,ccore_charge,cerf_ppl,&
                core_charge_radius,gcca,kind_radius,ppl_radius,zeta
    INTEGER  :: iexp_ppl,ikind,ipgf,iset,ishell,l,lppl,nexp_ppl,nkind,nset,i

    REAL(KIND = dp), DIMENSION(:), POINTER     :: set_radius,cexp_ppl,&
                                                  alpha_ppl_elp
    INTEGER, DIMENSION(:), POINTER      :: npgf,nshell,nexp_ppl_elp
    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: pgf_radius,zet,cexp_ppl_elp
    INTEGER, DIMENSION(:,:), POINTER    :: lshell
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: gcc

!   ---------------------------------------------------------------------------

    NULLIFY(atomic_kind,elp_potential,gth_potential,orb_basis_set,aux_basis_set)
    NULLIFY(alpha_ppl_elp,cexp_ppl,set_radius,npgf,nshell,nexp_ppl_elp)
    NULLIFY(pgf_radius,zet,cexp_ppl_elp,lshell,gcc)

    nkind = SIZE(atomic_kind_set)

    nexp_ppl = 0

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           elp_potential=elp_potential,&
                           gth_potential=gth_potential,&
                           aux_basis_set=aux_basis_set,&
                           orb_basis_set=orb_basis_set)

      IF (ASSOCIATED(elp_potential)) THEN

        CALL get_potential(potential=elp_potential,&
                           alpha_core_charge=alpha_core_charge,&
                           ccore_charge=ccore_charge,&
                           alpha_ppl=alpha_ppl_elp,&
                           cerf_ppl=cerf_ppl,&
                           nexp_ppl=nexp_ppl_elp,&
                           cexp_ppl=cexp_ppl_elp)

!       *** Calculate the radius of the core charge distribution ***

        core_charge_radius = exp_radius(0,alpha_core_charge,&
                                        kg_control%eps_core_charge,&
                                        ccore_charge)

!       *** Calculate the radii of the  ***
!       *** effective local pseudopotential (elp) ***

        ppl_radius = exp_radius(0,alpha_ppl_elp(1),kg_control%eps_ppl,cerf_ppl)

        DO i=1,size(nexp_ppl_elp)
          DO iexp_ppl=1,nexp_ppl_elp(i)
            lppl = 2*(iexp_ppl - 1)
            ppl_radius = MAX(ppl_radius,&
                           exp_radius(lppl,alpha_ppl_elp(i),&
                                      kg_control%eps_core_charge,&
                                      cexp_ppl_elp(i,iexp_ppl)))
          END DO
        END DO


        CALL set_potential(potential=elp_potential,&
                           core_charge_radius=core_charge_radius,&
                           ppl_radius=ppl_radius)


      ELSE IF (ASSOCIATED(gth_potential)) THEN

        CALL get_potential(potential=gth_potential,&
                           alpha_core_charge=alpha_core_charge,&
                           ccore_charge=ccore_charge,&
                           alpha_ppl=alpha_ppl,&
                           cerf_ppl=cerf_ppl,&
                           nexp_ppl=nexp_ppl,&
                           cexp_ppl=cexp_ppl) 

!       *** Calculate the radius of the core charge distribution ***

        core_charge_radius = exp_radius(0,alpha_core_charge,&
                                        kg_control%eps_core_charge,&
                                        ccore_charge)

!       *** Calculate the radii of the local part  ***
!       *** of the Goedecker pseudopotential (GTH) ***

        ppl_radius = exp_radius(0,alpha_ppl,kg_control%eps_ppl,cerf_ppl)

        DO iexp_ppl=1,nexp_ppl
          lppl = 2*(iexp_ppl - 1)
          ppl_radius = MAX(ppl_radius,&
                           exp_radius(lppl,alpha_ppl,&
                                      kg_control%eps_core_charge,&
                                      cexp_ppl(iexp_ppl)))
        END DO


        CALL set_potential(potential=gth_potential,&
                           core_charge_radius=core_charge_radius,&
                           ppl_radius=ppl_radius) 

      END IF


!     *** Calculate the orbital basis function radii (frozen density)***

      IF (ASSOCIATED(orb_basis_set)) THEN

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
                                                     kg_control%eps_pgf_orb,&
                                                     gcca))
            END DO
            set_radius(iset) = MAX(set_radius(iset),pgf_radius(ipgf,iset))
          END DO
          kind_radius = MAX(kind_radius,set_radius(iset))
        END DO

        CALL set_gto_basis_set(gto_basis_set=orb_basis_set,&
                               pgf_radius=pgf_radius,&
                               set_radius=set_radius,&
                               kind_radius=kind_radius)

      END IF

!     *** Calculate the orbital basis function radii (pol density) ***

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
                                                     kg_control%eps_pgf_orb,&
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



    END DO

    IF (ASSOCIATED(orb_basis_set)) THEN
      CALL write_pgf_orb_radii(cell,atomic_kind_set,globenv)
    END IF
    CALL write_core_charge_radii(cell,atomic_kind_set,globenv)
    IF (nexp_ppl > 0) CALL write_ppl_radii(cell,atomic_kind_set,globenv)

  END SUBROUTINE init_interaction_radii_kg

! *****************************************************************************

  SUBROUTINE init_interaction_radii_qs(qs_control,cell,atomic_kind_set,globenv)

!   Purpose: Initialize all the atomic kind radii for a given threshold value.

!   History: - Creation (24.06.1999, Matthias Krack)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(cell_type), POINTER                      :: cell
    TYPE(global_environment_type), INTENT(IN)     :: globenv
    TYPE(qs_control_type)                         :: qs_control

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(all_potential_type), POINTER :: all_potential
    TYPE(gth_potential_type), POINTER :: gth_potential
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    REAL(KIND = dp) :: alpha_core_charge,alpha_ppl,ccore_charge,cerf_ppl,&
                core_charge_radius,gcca,kind_radius,ppl_radius,ppnl_radius,zeta
    INTEGER  :: iexp_ppl,ikind,ipgf,iprj_ppnl,iset,ishell,l,lppl,lppnl,&
                lprj_ppnl,nexp_ppl,nkind,nset

    REAL(KIND = dp), DIMENSION(:), POINTER     :: alpha_ppnl,cexp_ppl,set_radius
    INTEGER, DIMENSION(:), POINTER      :: npgf,nprj_ppnl,nshell
    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: cprj_ppnl,pgf_radius,zet
    INTEGER, DIMENSION(:,:), POINTER    :: lshell
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: gcc

!   ---------------------------------------------------------------------------

    nkind = SIZE(atomic_kind_set)

    nexp_ppl = 0
    lppnl = -1

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           all_potential=all_potential,&
                           gth_potential=gth_potential,&
                           orb_basis_set=orb_basis_set)

      IF (ASSOCIATED(all_potential)) THEN

        CALL get_potential(potential=all_potential,&
                           alpha_core_charge=alpha_core_charge,&
                           ccore_charge=ccore_charge)

!       *** Calculate the radius of the core charge distribution ***

        core_charge_radius = exp_radius(0,alpha_core_charge,&
                                        qs_control%eps_core_charge,&
                                        ccore_charge)

        CALL set_potential(potential=all_potential,&
                           core_charge_radius=core_charge_radius)

      ELSE IF (ASSOCIATED(gth_potential)) THEN

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
      END IF

!     *** Calculate the orbital basis function radii ***

      IF (ASSOCIATED(orb_basis_set)) THEN

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

        CALL set_gto_basis_set(gto_basis_set=orb_basis_set,&
                               pgf_radius=pgf_radius,&
                               set_radius=set_radius,&
                               kind_radius=kind_radius)

      END IF

    END DO

    CALL write_pgf_orb_radii(cell,atomic_kind_set,globenv)
    CALL write_core_charge_radii(cell,atomic_kind_set,globenv)
    IF (nexp_ppl > 0) CALL write_ppl_radii(cell,atomic_kind_set,globenv)
    IF (lppnl > -1) CALL write_ppnl_radii(cell,atomic_kind_set,globenv)

  END SUBROUTINE init_interaction_radii_qs

! *****************************************************************************

  SUBROUTINE write_core_charge_radii(cell,atomic_kind_set,globenv)

!   Purpose: Write the radii of the core charge distributions to the output
!            unit.

!   History: - Creation (15.09.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(cell_type), POINTER                      :: cell
    TYPE(global_environment_type), INTENT(IN)     :: globenv

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(all_potential_type), POINTER :: all_potential
    TYPE(elp_potential_type), POINTER :: elp_potential
    TYPE(gth_potential_type), POINTER :: gth_potential
    CHARACTER(LEN=60)                 :: name
    CHARACTER(LEN=8)                  :: unit_of_length_name
    REAL(KIND = dp)                          :: core_charge_radius,unit_of_length
    INTEGER                           :: ikind,nkind,output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    CALL get_cell(cell=cell,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

    output_unit = globenv%scr

    nkind = SIZE(atomic_kind_set)

    IF (globenv%print%core_charge_radii) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T22,A,/)")&
        "RADII OF THE CORE CHARGE DISTRIBUTIONS IN "//&
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
          WRITE (UNIT=output_unit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
            ikind,name,core_charge_radius/unit_of_length
        ELSE IF (ASSOCIATED(elp_potential)) THEN
          CALL get_potential(potential=elp_potential,&
                             core_charge_radius=core_charge_radius)
          WRITE (UNIT=output_unit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
            ikind,name,core_charge_radius/unit_of_length
        ELSE IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             core_charge_radius=core_charge_radius)
          WRITE (UNIT=output_unit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
            ikind,name,core_charge_radius/unit_of_length
        END IF
      END DO
    END IF

  END SUBROUTINE write_core_charge_radii

! *****************************************************************************

  SUBROUTINE write_pgf_orb_radii(cell,atomic_kind_set,globenv)

!   Purpose: Write the orbital basis function radii to the output unit.

!   History: - Creation (05.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(cell_type), POINTER                      :: cell
    TYPE(global_environment_type), INTENT(IN)     :: globenv

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    CHARACTER(LEN=60)                 :: name
    CHARACTER(LEN=8)                  :: unit_of_length_name
    REAL(KIND = dp)                          :: kind_radius,unit_of_length
    INTEGER                           :: ikind,ipgf,iset,nkind,nset,output_unit
    REAL(KIND = dp), DIMENSION(:), POINTER   :: set_radius
    INTEGER, DIMENSION(:), POINTER    :: npgf
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: pgf_radius

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    CALL get_cell(cell=cell,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

    output_unit = globenv%scr

    nkind = SIZE(atomic_kind_set)

!   *** Print the kind radii ***

    IF (globenv%print%kind_radii) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T22,A,/)")&
        "RADII OF THE ORBITAL BASIS IN "//TRIM(unit_of_length_name),&
        "Kind","Label","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,&
                             orb_basis_set=orb_basis_set)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               kind_radius=kind_radius)
        WRITE (UNIT=output_unit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
          ikind,name,kind_radius/unit_of_length
      END DO
    END IF

!   *** Print the shell set radii ***

    IF (globenv%print%set_radii) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T18,A,T28,A)")&
        "RADII OF THE SHELL SETS OF THE ORBITAL BASIS IN "//&
        TRIM(unit_of_length_name),"Kind","Label","Set","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,&
                             orb_basis_set=orb_basis_set)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               nset=nset,&
                               set_radius=set_radius)
        WRITE (UNIT=output_unit,FMT="(/,T2,I5,3X,A5,(T16,I5,T23,F12.6))")&
          ikind,name,(iset,set_radius(iset)/unit_of_length,iset=1,nset)
      END DO
    END IF

!   *** Print the primitive Gaussian function radii ***

    IF (globenv%print%pgf_radii) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T18,A,T28,A)")&
        "RADII OF THE PRIMITIVE GAUSSIAN FUNCTIONS OF THE ORBITAL BASIS IN "//&
        TRIM(unit_of_length_name),"Kind","Label","Set","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,&
                             orb_basis_set=orb_basis_set)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               nset=nset,&
                               npgf=npgf,&
                               pgf_radius=pgf_radius)
        DO iset=1,nset
          WRITE (UNIT=output_unit,FMT="(/,T2,I5,3X,A5,T16,I5,(T23,F12.6))")&
            ikind,name,iset,&
            (pgf_radius(ipgf,iset)/unit_of_length,ipgf=1,npgf(iset))
        END DO
      END DO
    END IF

  END SUBROUTINE write_pgf_orb_radii

! *****************************************************************************

  SUBROUTINE write_ppl_radii(cell,atomic_kind_set,globenv)

!   Purpose: Write the radii of the exponential functions of the Goedecker
!            pseudopotential (GTH, local part) to the logical unit number
!            "output_unit".

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(cell_type), POINTER                      :: cell
    TYPE(global_environment_type), INTENT(IN)     :: globenv

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(gth_potential_type), POINTER :: gth_potential
    CHARACTER(LEN=60)                 :: name
    CHARACTER(LEN=8)                  :: unit_of_length_name
    REAL(KIND = dp)                          :: ppl_radius,unit_of_length
    INTEGER                           :: ikind,nkind,output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%ppl_radii) THEN
      CALL get_cell(cell=cell,&
                    unit_of_length=unit_of_length,&
                    unit_of_length_name=unit_of_length_name)
      output_unit = globenv%scr
      nkind = SIZE(atomic_kind_set)
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T22,A,/)")&
        "RADII OF THE LOCAL PART OF THE GTH PSEUDOPOTENTIAL IN "//&
        TRIM(unit_of_length_name),"Kind","Label","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,&
                             gth_potential=gth_potential)
        IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             ppl_radius=ppl_radius)
          WRITE (UNIT=output_unit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
            ikind,name,ppl_radius/unit_of_length
        END IF
      END DO
    END IF

  END SUBROUTINE write_ppl_radii

! *****************************************************************************

  SUBROUTINE write_ppnl_radii(cell,atomic_kind_set,globenv)

!   Purpose: Write the radii of the projector functions of the Goedecker
!            pseudopotential (GTH, non-local part) to the logical unit number
!            "output_unit".

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(cell_type), POINTER                      :: cell
    TYPE(global_environment_type), INTENT(IN)     :: globenv

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(gth_potential_type), POINTER :: gth_potential
    CHARACTER(LEN=60)                 :: name
    CHARACTER(LEN=8)                  :: unit_of_length_name
    REAL(KIND = dp)                          :: ppnl_radius,unit_of_length
    INTEGER                           :: ikind,nkind,output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%ppnl_radii) THEN
      CALL get_cell(cell=cell,&
                    unit_of_length=unit_of_length,&
                    unit_of_length_name=unit_of_length_name)
      output_unit = globenv%scr
      nkind = SIZE(atomic_kind_set)
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T22,A,/)")&
        "RADII OF THE NON-LOCAL PART OF THE GTH PSEUDOPOTENTIAL IN "//&
        TRIM(unit_of_length_name),"Kind","Label","Radius"
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,&
                             gth_potential=gth_potential)
        IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             ppnl_radius=ppnl_radius)
          WRITE (UNIT=output_unit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
            ikind,name,ppnl_radius/unit_of_length
        END IF
      END DO
    END IF

  END SUBROUTINE write_ppnl_radii


! *****************************************************************************

    FUNCTION exp_radius(l,alpha,threshold,prefactor,epsin) RESULT(radius)

!   Purpose: The radius of a primitive Gaussian function for a given threshold
!            is calculated.

!            g(r) = prefactor*r**l*exp(-alpha*r**2) - threshold = 0

!   History: - Creation (26.02.1999, Matthias Krack)
!            - Modified to bisection to avoid division by zero (01.2002, Joost)
!            - that means the old version contains bugs .......(04.2002, Joost)
!              so keep this one ...

!   ***************************************************************************

!   alpha    : Exponent of the primitive Gaussian function.
!   l        : Angular momentum quantum number l.
!   prefactor: Prefactor of the Gaussian function (e.g. a contraction
!              coefficient).
!   radius   : Calculated radius of the Gaussian function.
!   threshold: Threshold for radius.

!   epsiter: Convergence criterion.
!   g      : The function g(r).
!   maxiter: Maximum number of iterations.

!   ***************************************************************************

    USE termination, ONLY: stop_program

    REAL(KIND = dp), INTENT(IN) :: alpha,prefactor,threshold
    INTEGER, INTENT(IN)  :: l
    REAL(KIND = dp), INTENT(IN), OPTIONAL :: epsin

    REAL(KIND = dp) :: radius

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "FUNCTION exp_radius (MODULE mathlib)"

    REAL(KIND = dp) :: epsiter 
    INTEGER, PARAMETER  :: maxiter = 5000

!   *** Local variables ***

    REAL(KIND = dp) :: a,ar2,d,delta,dg,g,r,t,rlow,rhigh,rmid
    INTEGER  :: iter

!   ---------------------------------------------------------------------------
    if (present(epsin)) then
       epsiter=epsin
    else
       epsiter=1.0E-12_dp
    endif

!   *** Initialize function value ***

    radius = 0.0_dp

!   *** Load and check parameter values ***

    IF (l < 0) THEN
      CALL stop_program(routine,"The angular momentum quantum number is "//&
                                "negative")
    END IF

    IF (alpha == 0.0_dp) THEN
      CALL stop_program(routine,"The Gaussian function exponent is zero")
    ELSE
      a = ABS(alpha)
    END IF

    IF (threshold == 0.0_dp) THEN
      CALL stop_program(routine,"The requested threshold is zero")
    ELSE
      t = ABS(threshold)
    END IF

    IF (prefactor == 0.0_dp) THEN
      RETURN
    ELSE
      d = ABS(prefactor)
    END IF

!   *** Calculate the Maximum g(r) ***
    r = SQRT(0.5_dp*REAL(l,dp)/a)
    ar2 = a*r*r

    IF (l == 0) THEN
      g = d
    ELSE
      g = d*r**l*EXP(-ar2)
    END IF

    IF (t > g) THEN
      RETURN
    END IF

    rlow  = r
    rhigh = 2.0*rlow+1.0
    iter=0
    DO 
       iter=iter+1
       if (iter.gt.maxiter) CALL stop_program(routine,"Maximum number of iterations exceeded")
       g = d*rhigh**l*EXP(-a*rhigh**2)
       if (g < t) exit
       rlow = rhigh
       rhigh =  2.0*rlow+1.0 
    ENDDO

    DO iter=1,maxiter
       rmid= (rlow+rhigh)*0.5_dp
       ar2 = a*rmid*rmid
       g   = d*rmid**l*EXP(-ar2)
       if  (g.lt.t) then
           rhigh=rmid
       else
           rlow=rmid
       endif
       IF (ABS(rhigh-rlow).lt.epsiter) THEN
          radius=rhigh
          return
       ENDIF
    ENDDO 
    CALL stop_program(routine,"Maximum number of iterations exceeded")
     
  END FUNCTION exp_radius


! *****************************************************************************

! provides a radius that is keeping into account all the angular momentum
! components given the product of the two gaussians, and their original
! centers
!******
  FUNCTION exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,o1,o2,ra,rb,rp,&
                          zetp,eps,prefactor,cutoff) result(radius)

    integer ,INTENT(IN)     :: la_min,la_max,lb_min,lb_max
    real(KIND = dp) ,INTENT(IN)    :: ra(3),rb(3),rp(3),zetp,eps,prefactor
    real(KIND = dp) ,INTENT(IN)    :: cutoff
    real(KIND = dp) , DIMENSION(:,:), POINTER, optional   :: pab
    integer,optional        :: o1,o2
    real(KIND = dp)    :: radius
    integer     :: lxa,lxb,lya,lyb,lza,lzb,la(3),lb(3),ico,jco,handle
    real(KIND = dp)    :: prefactor_local
    real(KIND = dp) polycoef(0:60),coef(0:20,3),bini,binj,s1,s2
    integer :: d,i,j,k

    radius=0.0_dp
    polycoef(:)=0.0_dp
    DO lxa=0,la_max
    DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
          DO lza=max(la_min-lxa-lya,0),la_max-lxa-lya
          DO lzb=max(lb_min-lxb-lyb,0),lb_max-lxb-lyb
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
                polycoef(i+j+k)=max(polycoef(i+j+k),&
                                    abs(prefactor_local*&
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
          radius=max(radius,exp_radius(i,zetp,eps,polycoef(i),1.0E-2_dp) )
    ENDDO

  END FUNCTION exp_radius_very_extended

! *****************************************************************************

END MODULE qs_interactions
