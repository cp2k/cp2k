!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/interactions [1.0] *
!!
!!   NAME
!!     interactions
!!
!!   FUNCTION
!!     Calculate the interaction radii for the operator matrix calculation.
!!
!!   AUTHOR
!!     Matthias Krack (12.07.2000)
!!
!!   MODIFICATION HISTORY
!!     Joost VandeVondele : added exp_radius_very_extended
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE interactions

! *****************************************************************************

! Index:

! SUBROUTINE read_interaction_thresholds(start_section,end_section)
! SUBROUTINE write_interaction_thresholds(lunit)
! SUBROUTINE init_interaction_radii(globenv)
! SUBROUTINE calculate_core_charge_radii
! SUBROUTINE write_core_charge_radii(lunit)
! SUBROUTINE calculate_pgf_aux_radii
! SUBROUTINE write_pgf_aux_radii(lunit)
! SUBROUTINE calculate_pgf_orb_radii
! SUBROUTINE write_pgf_orb_radii(lunit)
! SUBROUTINE calculate_ppl_radii
! SUBROUTINE write_ppl_radii(lunit)
! SUBROUTINE calculate_ppnl_radii
! SUBROUTINE write_ppnl_radii(lunit)

! FUNCTION exp_radius(l,alpha,threshold,prefactor) RESULT(radius)
! FUNCTION exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,ra,rb,rp,zetp,eps,prefactor,cutoff) result(radius)

! *****************************************************************************

! eps_pgf_orb: Cutoff value for the interaction of the primitive Gaussian-type
!              functions (primitive basis functions).

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE atomic_kinds,          ONLY: all_potential_present,&
                                   gth_potential_present,&
                                   kind_info,&
                                   nkind
  USE global_types,          ONLY: global_environment_type
  USE method_specifications, ONLY: allchem,gpw
  USE orbital_pointers,      ONLY: coset,indco,ncoset
  USE timings,               ONLY: timeset,timestop

  IMPLICIT NONE

  PRIVATE

  REAL(wp) :: eps_core_charge,&
              eps_pgf_aux,&
              eps_pgf_orb,&
              eps_ppl,&
              eps_ppnl,&
              eps_rho_gspace,&
              eps_rho_rspace,&
              eps_gvg

! *** Public variables ***

  PUBLIC :: eps_rho_gspace,&
            eps_rho_rspace,&
            eps_gvg

! *** Public subroutines ***

  PUBLIC :: init_interaction_radii,&
            read_interaction_thresholds,&
            write_interaction_thresholds

! *** Public functions ***

  PUBLIC :: exp_radius
  PUBLIC :: exp_radius_very_extended

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE read_interaction_thresholds(start_section,end_section,globenv)

!   Purpose: Read the method specifications.

!   History: - Creation (13.06.2000, Matthias Krack)

!   ***************************************************************************

    USE input_utilities,  ONLY: finish_input_session,&
                                read_object,&
                                search,&
                                start_input_session
    USE string_utilities, ONLY: uppercase
    USE termination,      ONLY: stop_program

    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=*), INTENT(IN)              :: start_section,end_section

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_interaction_thresholds (MODULE interactions)"

!   *** Local variables ***

    CHARACTER(LEN=40) :: keyword
    REAL(wp)          :: eps_default
    LOGICAL           :: found

!   ---------------------------------------------------------------------------

!   *** Load the default values ***

    eps_core_charge = 1.0E-12_wp
    eps_pgf_aux = 1.0E-6_wp
    eps_pgf_orb = 1.0E-6_wp
    eps_ppl = 1.0E-12_wp
    eps_ppnl = 1.0E-12_wp
    eps_rho_gspace = 1.0E-8_wp
    eps_rho_rspace = 1.0E-8_wp
    eps_gvg= 1.0E-8_wp

    CALL start_input_session(globenv%input_file_name,globenv)

!   *** Search for the input section with the cell parameters ***

    CALL search(start_section,.TRUE.,found)

    IF (found) THEN

      DO

        CALL read_object(keyword,skip_lines=1)

        CALL uppercase(keyword)

        SELECT CASE (TRIM(keyword))
        CASE ("EPS_DEFAULT","DEFAULT")
          CALL read_object(eps_default)
          eps_core_charge = eps_default
          eps_pgf_aux = eps_default
          eps_pgf_orb = eps_default
          eps_ppl = eps_default
          eps_ppnl = eps_default
          eps_rho_gspace = eps_default
          eps_rho_rspace = eps_default
          eps_gvg = eps_default
        CASE ("EPS_CORE_CHARGE")
          CALL read_object(eps_core_charge)
        CASE ("EPS_PGF_AUX")
          CALL read_object(eps_pgf_aux)
        CASE ("EPS_PGF_ORB")
          CALL read_object(eps_pgf_orb)
        CASE ("EPS_PPL")
          CALL read_object(eps_ppl)
        CASE ("EPS_PPNL")
          CALL read_object(eps_ppnl)
        CASE ("EPS_RHO")
          CALL read_object(eps_rho_gspace)
          eps_rho_rspace = eps_rho_gspace
          eps_gvg = eps_rho_rspace
        CASE ("EPS_RHO_GSPACE")
          CALL read_object(eps_rho_gspace)
          eps_gvg = eps_rho_gspace
        CASE ("EPS_RHO_RSPACE")
          CALL read_object(eps_rho_rspace)
          eps_gvg = eps_rho_rspace
        CASE ("EPS_GVG")
          CALL read_object(eps_gvg)
        CASE DEFAULT
          IF (keyword == end_section) THEN
            EXIT
          ELSE
            CALL stop_program(routine,&
                              "Invalid keyword <"//TRIM(keyword)//&
                              "> found in the input section <"//&
                              TRIM(start_section)//">")
          END IF
        END SELECT

      END DO

    END IF

    CALL finish_input_session

  END SUBROUTINE read_interaction_thresholds

! *****************************************************************************

  SUBROUTINE write_interaction_thresholds(lunit)

!   Purpose: Write the interaction threshold values to the logical unit number
!            "lunit".

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: lunit

!   ---------------------------------------------------------------------------

    WRITE (UNIT=lunit,FMT="(/,/,T2,A,/)") "INTERACTION THRESHOLDS"

    WRITE (UNIT=lunit,FMT="(T3,A,T20,ES8.1)") "eps_pgf_orb:",eps_pgf_orb

    IF (allchem) THEN
      WRITE (UNIT=lunit,FMT="(T3,A,T20,ES8.1)")&
        "eps_pgf_aux:",eps_pgf_aux
    ELSE
      WRITE (UNIT=lunit,FMT="(T3,A,T20,ES8.1)")&
        "eps_core_charge:",eps_core_charge
    END IF

    IF (gth_potential_present) THEN
      WRITE (UNIT=lunit,FMT="(T3,A,T20,ES8.1)")&
        "eps_ppl:",eps_ppl,&
        "eps_ppnl:",eps_ppnl
    END IF

    WRITE (UNIT=lunit,FMT="(T3,A,ES8.1,/,T3,A,ES8.1,/,T3,A,ES8.1)")&
      "eps_rho_rspace:",eps_rho_rspace,&
      "eps_rho_gspace:",eps_rho_gspace,&
      "eps_gvg:",eps_gvg

  END SUBROUTINE write_interaction_thresholds

! *****************************************************************************

  SUBROUTINE init_interaction_radii(globenv)

!   Purpose: Initialize all the atomic kind radii for a given threshold value.

!   History: - Creation (24.06.1999, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: output_unit
    LOGICAL :: ionode

!   ---------------------------------------------------------------------------

    output_unit = globenv%scr
    ionode = globenv%ionode

    CALL calculate_pgf_orb_radii
    IF (ionode) CALL write_pgf_orb_radii(globenv)

    IF (allchem) THEN
      CALL calculate_pgf_aux_radii
      IF (ionode) CALL write_pgf_aux_radii(globenv)
    ELSE
      CALL calculate_core_charge_radii
      IF (ionode.AND.globenv%print%core_charge_radii) THEN
        CALL write_core_charge_radii(output_unit)
      END IF
    END IF

    IF (gth_potential_present) THEN
      CALL calculate_ppl_radii
      IF (ionode.AND.globenv%print%ppl_radii) THEN
        CALL write_ppl_radii(output_unit)
      END IF
      CALL calculate_ppnl_radii
      IF (ionode.AND.globenv%print%ppnl_radii) THEN
        CALL write_ppnl_radii(output_unit)
      END IF
    END IF

  END SUBROUTINE init_interaction_radii

! *****************************************************************************

  SUBROUTINE calculate_core_charge_radii

!   Purpose: Calculate the radii of the core charge distributions for a given
!            threshold.

!   History: - Creation (15.09.2000, Matthias Krack)

!   ***************************************************************************

!   *** Local variables ***

    INTEGER :: ikind

!   ---------------------------------------------------------------------------

    DO ikind=1,nkind
      kind_info(ikind)%core_charge_radius =&
        exp_radius(0,kind_info(ikind)%alpha_core_charge,eps_core_charge,&
                   kind_info(ikind)%ccore_charge)
    END DO

  END SUBROUTINE calculate_core_charge_radii

! *****************************************************************************

  SUBROUTINE write_core_charge_radii(lunit)

!   Purpose: Write the radii of the core charge distributions to the logical
!            unit number "lunit".

!   History: - Creation (15.09.2000, Matthias Krack)

!   ***************************************************************************

    USE cell_parameters, ONLY: unit_of_length,unit_of_length_name

    INTEGER, INTENT(IN) :: lunit

!   *** Local variables ***

    INTEGER :: ikind

!   ---------------------------------------------------------------------------

    WRITE (UNIT=lunit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T22,A,/)")&
      "RADII OF THE CORE CHARGE DISTRIBUTIONS IN "//&
      TRIM(unit_of_length_name),"Kind","Label","Radius"

    DO ikind=1,nkind
      WRITE (UNIT=lunit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
        ikind,kind_info(ikind)%label,&
        kind_info(ikind)%core_charge_radius/unit_of_length
    END DO

  END SUBROUTINE write_core_charge_radii

! *****************************************************************************

  SUBROUTINE calculate_pgf_aux_radii

!   Purpose: Calculate the radii of the auxiliary basis functions for a given
!            threshold.

!   History: - Creation (17.11.2000, Matthias Krack)

!   ***************************************************************************

    REAL(wp) :: gcca,kind_radius,pgf_radius,set_radius,zeta
    INTEGER  :: ipgf,ikind,iset,ishell,l

!   ---------------------------------------------------------------------------

    DO ikind=1,nkind
      kind_radius = 0.0_wp
      DO iset=1,kind_info(ikind)%aux_basis_set%nset
        set_radius = 0.0_wp
        DO ipgf=1,kind_info(ikind)%aux_basis_set%npgf(iset)
          pgf_radius = 0.0_wp
          DO ishell=1,kind_info(ikind)%aux_basis_set%nshell(iset)
            l = kind_info(ikind)%aux_basis_set%l(ishell,iset)
            gcca = kind_info(ikind)%aux_basis_set%gcc(ipgf,ishell,iset)
            zeta = kind_info(ikind)%aux_basis_set%zet(ipgf,iset)
            pgf_radius = MAX(pgf_radius,exp_radius(l,zeta,eps_pgf_aux,gcca))
          END DO
          kind_info(ikind)%aux_basis_set%pgf_radius(ipgf,iset) = pgf_radius
          set_radius = MAX(set_radius,pgf_radius)
        END DO
        kind_info(ikind)%aux_basis_set%set_radius(iset) = set_radius
        kind_radius = MAX(kind_radius,set_radius)
      END DO
      kind_info(ikind)%aux_basis_set%kind_radius = kind_radius
    END DO

  END SUBROUTINE calculate_pgf_aux_radii

! *****************************************************************************

  SUBROUTINE write_pgf_aux_radii(globenv)

!   Purpose: Write the auxiliary basis function radii to the output unit.

!   History: - Creation (17.11.2000, Matthias Krack)

!   ***************************************************************************

    USE cell_parameters, ONLY: unit_of_length,unit_of_length_name

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: ikind,ipgf,iset,output_unit

!   ---------------------------------------------------------------------------

    output_unit = globenv%scr

!   *** Print the kind radii ***

    IF (globenv%print%kind_radii) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T22,A,/)")&
        "RADII OF THE AUXILIARY BASIS IN "//TRIM(unit_of_length_name),&
        "Kind","Label","Radius"
      WRITE (UNIT=output_unit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
        (ikind,&
         kind_info(ikind)%label,&
         kind_info(ikind)%aux_basis_set%kind_radius/unit_of_length,&
         ikind=1,nkind)
    END IF

!   *** Print the shell set radii ***

    IF (globenv%print%set_radii) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T18,A,T28,A)")&
        "RADII OF THE SHELL SETS OF THE AUXILIARY BASIS IN "//&
       TRIM(unit_of_length_name),&
        "Kind","Label","Set","Radius"
      DO ikind=1,nkind
        WRITE (UNIT=output_unit,FMT="(/,T2,I5,3X,A5,(T16,I5,T23,F12.6))")&
          ikind,kind_info(ikind)%label,&
          (iset,&
           kind_info(ikind)%aux_basis_set%set_radius(iset)/unit_of_length,&
           iset=1,kind_info(ikind)%aux_basis_set%nset)
      END DO
    END IF

!   *** Print the primitive Gaussian function radii ***

    IF (globenv%print%pgf_radii) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T18,A,T28,A)")&
        "RADII OF THE PRIMITIVE GAUSSIAN FUNCTIONS OF THE AUXILIARY BASIS "//&
        "IN "//TRIM(unit_of_length_name),&
        "Kind","Label","Set","Radius"
      DO ikind=1,nkind
        DO iset=1,kind_info(ikind)%aux_basis_set%nset
          WRITE (UNIT=output_unit,FMT="(/,T2,I5,3X,A5,T16,I5,(T23,F12.6))")&
            ikind,kind_info(ikind)%label,iset,&
            (kind_info(ikind)%aux_basis_set%pgf_radius(ipgf,iset)/&
             unit_of_length,&
             ipgf=1,kind_info(ikind)%aux_basis_set%npgf(iset))
        END DO
      END DO
    END IF

  END SUBROUTINE write_pgf_aux_radii

! *****************************************************************************

  SUBROUTINE calculate_pgf_orb_radii

!   Purpose: Calculate the radii of the orbital basis functions for a given
!            threshold.

!   History: - Creation (24.06.1999, Matthias Krack)

!   ***************************************************************************

    REAL(wp) :: gcca,kind_radius,pgf_radius,set_radius,zeta
    INTEGER  :: ipgf,ikind,iset,ishell,l

!   ---------------------------------------------------------------------------

    DO ikind=1,nkind
      kind_radius = 0.0_wp
      DO iset=1,kind_info(ikind)%orb_basis_set%nset
        set_radius = 0.0_wp
        DO ipgf=1,kind_info(ikind)%orb_basis_set%npgf(iset)
          pgf_radius = 0.0_wp
          DO ishell=1,kind_info(ikind)%orb_basis_set%nshell(iset)
            l = kind_info(ikind)%orb_basis_set%l(ishell,iset)
            gcca = kind_info(ikind)%orb_basis_set%gcc(ipgf,ishell,iset)
            zeta = kind_info(ikind)%orb_basis_set%zet(ipgf,iset)
            pgf_radius = MAX(pgf_radius,exp_radius(l,zeta,eps_pgf_orb,gcca))
          END DO
          kind_info(ikind)%orb_basis_set%pgf_radius(ipgf,iset) = pgf_radius
          set_radius = MAX(set_radius,pgf_radius)
        END DO
        kind_info(ikind)%orb_basis_set%set_radius(iset) = set_radius
        kind_radius = MAX(kind_radius,set_radius)
      END DO
      kind_info(ikind)%orb_basis_set%kind_radius = kind_radius
    END DO

  END SUBROUTINE calculate_pgf_orb_radii

! *****************************************************************************

  SUBROUTINE write_pgf_orb_radii(globenv)

!   Purpose: Write the orbital basis function radii to the output unit.

!   History: - Creation (05.06.2000, Matthias Krack)

!   ***************************************************************************

    USE cell_parameters, ONLY: unit_of_length,unit_of_length_name

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: ikind,ipgf,iset,output_unit

!   ---------------------------------------------------------------------------

    output_unit = globenv%scr

!   *** Print the kind radii ***

    IF (globenv%print%kind_radii) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T22,A,/)")&
        "RADII OF THE ORBITAL BASIS IN "//TRIM(unit_of_length_name),&
        "Kind","Label","Radius"
      WRITE (UNIT=output_unit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
        (ikind,&
         kind_info(ikind)%label,&
         kind_info(ikind)%orb_basis_set%kind_radius/unit_of_length,&
         ikind=1,nkind)
    END IF

!   *** Print the shell set radii ***

    IF (globenv%print%set_radii) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T18,A,T28,A)")&
        "RADII OF THE SHELL SETS OF THE ORBITAL BASIS IN "//&
       TRIM(unit_of_length_name),&
        "Kind","Label","Set","Radius"
      DO ikind=1,nkind
        WRITE (UNIT=output_unit,FMT="(/,T2,I5,3X,A5,(T16,I5,T23,F12.6))")&
          ikind,kind_info(ikind)%label,&
          (iset,&
           kind_info(ikind)%orb_basis_set%set_radius(iset)/unit_of_length,&
           iset=1,kind_info(ikind)%orb_basis_set%nset)
      END DO
    END IF

!   *** Print the primitive Gaussian function radii ***

    IF (globenv%print%pgf_radii) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T18,A,T28,A)")&
        "RADII OF THE PRIMITIVE GAUSSIAN FUNCTIONS OF THE ORBITAL BASIS IN "//&
        TRIM(unit_of_length_name),&
        "Kind","Label","Set","Radius"
      DO ikind=1,nkind
        DO iset=1,kind_info(ikind)%orb_basis_set%nset
          WRITE (UNIT=output_unit,FMT="(/,T2,I5,3X,A5,T16,I5,(T23,F12.6))")&
            ikind,kind_info(ikind)%label,iset,&
            (kind_info(ikind)%orb_basis_set%pgf_radius(ipgf,iset)/&
             unit_of_length,&
             ipgf=1,kind_info(ikind)%orb_basis_set%npgf(iset))
        END DO
      END DO
    END IF

  END SUBROUTINE write_pgf_orb_radii

! *****************************************************************************

  SUBROUTINE calculate_ppl_radii

!   Purpose: Calculate the radii of the local part of the Goedecker pseudo-
!            potential (GTH).

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    REAL(wp) :: alpha_ppl,cerf_ppl,cexp_ppl,ppl_radius
    INTEGER  :: iexp_ppl,ikind,lppl

!   ---------------------------------------------------------------------------

    DO ikind=1,nkind
      IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN
        alpha_ppl = kind_info(ikind)%gth_potential%alpha_ppl
        cerf_ppl = kind_info(ikind)%gth_potential%cerf_ppl
        ppl_radius = exp_radius(0,alpha_ppl,eps_ppl,cerf_ppl)
        DO iexp_ppl=1,kind_info(ikind)%gth_potential%nexp_ppl
          lppl = 2*(iexp_ppl - 1)
          cexp_ppl = kind_info(ikind)%gth_potential%cexp_ppl(iexp_ppl)
          ppl_radius = MAX(ppl_radius,&
                           exp_radius(lppl,alpha_ppl,eps_ppl,cexp_ppl))
        END DO
        kind_info(ikind)%gth_potential%ppl_radius = ppl_radius
      END IF
    END DO

  END SUBROUTINE calculate_ppl_radii

! *****************************************************************************

  SUBROUTINE write_ppl_radii(lunit)

!   Purpose: Write the radii of the exponential functions of the Goedecker
!            pseudopotential (GTH, local part) to the logical unit number
!            "lunit".

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    USE cell_parameters, ONLY: unit_of_length,unit_of_length_name

    INTEGER, INTENT(IN) :: lunit

!   *** Local variables ***

    INTEGER :: ikind

!   ---------------------------------------------------------------------------

    WRITE (UNIT=lunit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T22,A,/)")&
      "RADII OF THE LOCAL PART OF THE GTH PSEUDOPOTENTIAL IN "//&
      TRIM(unit_of_length_name),"Kind","Label","Radius"

    DO ikind=1,nkind
      IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN
        WRITE (UNIT=lunit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
          ikind,kind_info(ikind)%label,&
          kind_info(ikind)%gth_potential%ppl_radius/unit_of_length
      END IF
    END DO

  END SUBROUTINE write_ppl_radii

! *****************************************************************************

  SUBROUTINE calculate_ppnl_radii

!   Purpose: Calculate the radii of the projector functions of the Goedecker
!            pseudopotential (GTH, non-local part).

!   History: - Creation (11.10.2000, Matthias Krack)

!   ***************************************************************************

    REAL(wp) :: alpha_ppnl,cprj_ppnl,ppnl_radius
    INTEGER  :: ikind,iprj_ppnl,l,lprj_ppnl,nprj_ppnl

!   ---------------------------------------------------------------------------

    DO ikind=1,nkind
      IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN
        ppnl_radius = 0.0_wp
        DO l=0,kind_info(ikind)%gth_potential%lppnl
          alpha_ppnl = kind_info(ikind)%gth_potential%alpha_ppnl(l)
          nprj_ppnl = kind_info(ikind)%gth_potential%nprj_ppnl(l)
          DO iprj_ppnl=1,nprj_ppnl
            lprj_ppnl = l + 2*(iprj_ppnl - 1)
            cprj_ppnl = kind_info(ikind)%gth_potential%cprj_ppnl(iprj_ppnl,l)
            ppnl_radius = MAX(ppnl_radius,&
                              exp_radius(lprj_ppnl,alpha_ppnl,eps_ppnl,&
                                         cprj_ppnl))
          END DO
        END DO
        kind_info(ikind)%gth_potential%ppnl_radius = ppnl_radius
      END IF
    END DO

  END SUBROUTINE calculate_ppnl_radii

! *****************************************************************************

  SUBROUTINE write_ppnl_radii(lunit)

!   Purpose: Write the radii of the projector functions of the Goedecker
!            pseudopotential (GTH, non-local part) to the logical unit number
!            "lunit".

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    USE cell_parameters, ONLY: unit_of_length,unit_of_length_name

    INTEGER, INTENT(IN) :: lunit

!   *** Local variables ***

    INTEGER :: ikind

!   ---------------------------------------------------------------------------

    WRITE (UNIT=lunit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T22,A,/)")&
      "RADII OF THE NON-LOCAL PART OF THE GTH PSEUDOPOTENTIAL IN "//&
      TRIM(unit_of_length_name),"Kind","Label","Radius"

    DO ikind=1,nkind
      IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN
        WRITE (UNIT=lunit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
          ikind,kind_info(ikind)%label,&
          kind_info(ikind)%gth_potential%ppnl_radius/unit_of_length
      END IF
    END DO

  END SUBROUTINE write_ppnl_radii

! *****************************************************************************

    FUNCTION exp_radius(l,alpha,threshold,prefactor) RESULT(radius)

!   Purpose: The radius of a primitive Gaussian function for a given threshold
!            is calculated.

!            g(r) = prefactor*r**l*exp(-alpha*r**2) - threshold = 0

!   History: - Creation (26.02.1999, Matthias Krack)
!            - Modified to bisection to avoid division by zero (01.2002, Joost)

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

    REAL(wp), INTENT(IN) :: alpha,prefactor,threshold
    INTEGER, INTENT(IN)  :: l

    REAL(wp) :: radius

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "FUNCTION exp_radius (MODULE mathlib)"

    REAL(wp), PARAMETER :: epsiter = 1.0E-2_wp
    INTEGER, PARAMETER  :: maxiter = 500

!   *** Local variables ***

    REAL(wp) :: a,ar2,d,delta,dg,g,r,t,rlow,rhigh,rmid
    INTEGER  :: iter

!   ---------------------------------------------------------------------------

!   *** Initialize function value ***

    radius = 0.0_wp

!   *** Load and check parameter values ***

    IF (l < 0) THEN
      CALL stop_program(routine,"The angular momentum quantum number is "//&
                                "negative")
    END IF

    IF (alpha == 0.0_wp) THEN
      CALL stop_program(routine,"The Gaussian function exponent is zero")
    ELSE
      a = ABS(alpha)
    END IF

    IF (threshold == 0.0_wp) THEN
      CALL stop_program(routine,"The requested threshold is zero")
    ELSE
      t = ABS(threshold)
    END IF

    IF (prefactor == 0.0_wp) THEN
      RETURN
    ELSE
      d = ABS(prefactor)
    END IF

!   *** Calculate the Maximum g(r) ***
    r = SQRT(0.5_wp*REAL(l,wp)/a)
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
       rmid= (rlow+rhigh)*0.5_wp
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
!***********************************************************************************
! provides a radius that is keeping into account all the angular momentum components
! given the product of the two gaussians, and their original centers
!******
  FUNCTION exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,ra,rb,rp, &
                                    zetp,eps,prefactor,cutoff) result(radius)
  implicit none

  integer ,INTENT(IN)     :: la_min,la_max,lb_min,lb_max
  real(wp) ,INTENT(IN)    :: ra(3),rb(3),rp(3),zetp,eps,prefactor
  real(wp) ,INTENT(IN)    :: cutoff,pab(ncoset(la_max),ncoset(lb_max))
  real(wp)    :: radius
  integer     :: lxa,lxb,lya,lyb,lza,lzb,la(3),lb(3),ico,jco,handle
  real(wp)    :: prefactor_local
  real(wp) polycoef(0:60),coef(0:20,3),bini,binj,s1,s2
  integer d,i,j,k

  CALL timeset("rad","I","",handle)
  radius=0.0_wp
  polycoef(:)=0.0_wp
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
           prefactor_local=prefactor*MAX(ABS(pab(ico,jco)),cutoff)
           coef(:,:)=0.0_wp
           DO d=1,3
             bini=1.0_wp
             s1=1.0_wp
             DO i=0,la(d)
                binj=1.0_wp
                s2=1.0_wp
                DO j=0,lb(d)
                  coef(la(d)+lb(d)-i-j,d)=coef(la(d)+lb(d)-i-j,d)+bini*binj*s1*s2
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
              polycoef(i+j+k)=max(polycoef(i+j+k),abs(prefactor_local*coef(i,1)*coef(j,2)*coef(k,3)))
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
        radius=max(radius,exp_radius(i,zetp,eps,polycoef(i)) )
  ENDDO
  CALL timestop(0.0_wp,handle)
  END FUNCTION

! *****************************************************************************



! *****************************************************************************

END MODULE interactions
