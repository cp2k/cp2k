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
! SUBROUTINE calculate_erfc_radii
! SUBROUTINE write_erfc_radii(lunit)
! SUBROUTINE calculate_erf_ppl_radii
! SUBROUTINE write_erf_ppl_radii(lunit)
! SUBROUTINE calculate_exp_ppl_radii
! SUBROUTINE write_exp_ppl_radii(lunit)
! SUBROUTINE calculate_pgf_aux_radii
! SUBROUTINE write_pgf_aux_radii(lunit)
! SUBROUTINE calculate_pgf_orb_radii
! SUBROUTINE write_pgf_orb_radii(lunit)
! SUBROUTINE calculate_prj_ppnl_radii
! SUBROUTINE write_prj_ppnl_radii(lunit)

! FUNCTION erfc_radius(alpha,threshold,prefactor) RESULT(radius)
! FUNCTION exp_radius(l,alpha,threshold,prefactor) RESULT(radius)

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

  IMPLICIT NONE

  PRIVATE

  REAL(wp) :: eps_erfc,&
              eps_erf_ppl,&
              eps_exp_ppl,&
              eps_pgf_aux,&
              eps_pgf_orb,&
              eps_prj_ppnl

! *** Public subroutines ***

  PUBLIC :: init_interaction_radii,&
            read_interaction_thresholds,&
            write_interaction_thresholds

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

    eps_erfc = 1.0E-12_wp
    eps_erf_ppl = 1.0E-12_wp
    eps_exp_ppl = 1.0E-12_wp
    eps_pgf_aux = 1.0E-6_wp
    eps_pgf_orb = 1.0E-6_wp
    eps_prj_ppnl = 1.0E-12_wp

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
          eps_erfc = eps_default
          eps_erf_ppl = eps_default
          eps_exp_ppl = eps_default
          eps_pgf_aux = eps_default
          eps_pgf_orb = eps_default
          eps_prj_ppnl = eps_default
        CASE ("EPS_ERFC")
          CALL read_object(eps_erfc)
        CASE ("EPS_ERF_PPL")
          CALL read_object(eps_erf_ppl)
        CASE ("EPS_EXP_PPL")
          CALL read_object(eps_exp_ppl)
        CASE ("EPS_PGF_AUX")
          CALL read_object(eps_pgf_aux)
        CASE ("EPS_PGF_ORB")
          CALL read_object(eps_pgf_orb)
        CASE ("EPS_PRJ_PPNL","EPS_PPNL")
          CALL read_object(eps_prj_ppnl)
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

    WRITE (UNIT=lunit,FMT="(T3,A,ES8.1)") "eps_pgf_orb:  ",eps_pgf_orb

    IF (allchem) THEN
      WRITE (UNIT=lunit,FMT="(T3,A,ES8.1)") "eps_pgf_aux:  ",eps_pgf_aux
    END IF

    IF (gpw.AND.all_potential_present) THEN
      WRITE (UNIT=lunit,FMT="(T3,A,ES8.1)") "eps_erfc:     ",eps_erfc
    END IF

    IF (gth_potential_present) THEN
      IF (allchem) THEN
        WRITE (UNIT=lunit,FMT="(T3,A,ES8.1)") "eps_erf_ppl:  ",eps_erf_ppl
      END IF
      WRITE (UNIT=lunit,FMT="(T3,A,ES8.1)") "eps_exp_ppl:  ",eps_exp_ppl,&
                                   "eps_prj_ppnl: ",eps_prj_ppnl
    END IF

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
    END IF

    IF (gpw.AND.all_potential_present) THEN
      CALL calculate_erfc_radii
      IF (ionode.AND.globenv%print%erfc_radii) THEN
        CALL write_erfc_radii(output_unit)
      END IF
    END IF

    IF (gth_potential_present) THEN
      IF (allchem) THEN
        CALL calculate_erf_ppl_radii
        IF (ionode.AND.globenv%print%erf_ppl_radii) THEN
          CALL write_erf_ppl_radii(output_unit)
        END IF
      END IF
      CALL calculate_exp_ppl_radii
      IF (ionode.AND.globenv%print%exp_ppl_radii) THEN
        CALL write_exp_ppl_radii(output_unit)
      END IF
      CALL calculate_prj_ppnl_radii
      IF (ionode.AND.globenv%print%prj_ppnl_radii) THEN
        CALL write_prj_ppnl_radii(output_unit)
      END IF
    END IF

  END SUBROUTINE init_interaction_radii

! *****************************************************************************

  SUBROUTINE calculate_erfc_radii

!   Purpose: Calculate the radii of the nuclear charge distributions for a
!            given threshold.

!   History: - Creation (15.09.2000, Matthias Krack)

!   ***************************************************************************

!   *** Local variables ***

    INTEGER :: ikind

!   ---------------------------------------------------------------------------

    DO ikind=1,nkind
      IF (ASSOCIATED(kind_info(ikind)%all_potential)) THEN
        kind_info(ikind)%all_potential%erfc_radius =&
          erfc_radius(kind_info(ikind)%alpha_erf,eps_erfc,&
                      REAL(kind_info(ikind)%z,wp))
      END IF
    END DO

  END SUBROUTINE calculate_erfc_radii

! *****************************************************************************

  SUBROUTINE write_erfc_radii(lunit)

!   Purpose: Write the radii of the nuclear charge distributions to the logical
!            unit number "lunit".

!   History: - Creation (15.09.2000, Matthias Krack)

!   ***************************************************************************

    USE cell_parameters, ONLY: unit_of_length,unit_of_length_name

    INTEGER, INTENT(IN) :: lunit

!   *** Local variables ***

    INTEGER :: ikind

!   ---------------------------------------------------------------------------

    WRITE (UNIT=lunit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T22,A,/)")&
      "RADII OF THE NUCLEAR CHARGE DISTRIBUTIONS IN "//&
      TRIM(unit_of_length_name),"Kind","Label","Radius"

    DO ikind=1,nkind
      IF (ASSOCIATED(kind_info(ikind)%all_potential)) THEN
        WRITE (UNIT=lunit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
          ikind,kind_info(ikind)%label,&
          kind_info(ikind)%all_potential%erfc_radius/unit_of_length
      END IF
    END DO

  END SUBROUTINE write_erfc_radii

! *****************************************************************************

  SUBROUTINE calculate_erf_ppl_radii

!   Purpose: Calculate the radii of the error functions of the Goedecker
!            pseudopotential (GTH, local part).

!   History: - Creation (07.11.2000, Matthias Krack)

!   ***************************************************************************

!   *** Local variables ***

    REAL(wp) :: alpha,c
    INTEGER  :: ikind

!   ---------------------------------------------------------------------------

    DO ikind=1,nkind
      IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN
        alpha = kind_info(ikind)%gth_potential%alpha_exp_ppl
        c = kind_info(ikind)%cerf
        kind_info(ikind)%gth_potential%erf_ppl_radius =&
          exp_radius(0,alpha,eps_erf_ppl,c)
      END IF
    END DO

  END SUBROUTINE calculate_erf_ppl_radii

! *****************************************************************************

  SUBROUTINE write_erf_ppl_radii(lunit)

!   Purpose: Write the radii of the error functions of the Goedecker
!            pseudopotential (GTH, local part) to the logical unit number
!            "lunit".

!   History: - Creation (07.11.2000, Matthias Krack)

!   ***************************************************************************

    USE cell_parameters, ONLY: unit_of_length,unit_of_length_name

    INTEGER, INTENT(IN) :: lunit

!   *** Local variables ***

    INTEGER :: ikind

!   ---------------------------------------------------------------------------

    WRITE (UNIT=lunit,FMT="(/,/,T2,A,/,/,T3,A,T10,A,T22,A,/)")&
      "RADII OF THE ERROR FUNCTIONS OF THE GTH PSEUDOPOTENTIAL IN "//&
      TRIM(unit_of_length_name),"Kind","Label","Radius"

    DO ikind=1,nkind
      IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN
        WRITE (UNIT=lunit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
          ikind,kind_info(ikind)%label,&
          kind_info(ikind)%gth_potential%erf_ppl_radius/unit_of_length
      END IF
    END DO

  END SUBROUTINE write_erf_ppl_radii

! *****************************************************************************

  SUBROUTINE calculate_exp_ppl_radii

!   Purpose: Calculate the radii of the exponential functions of the Goedecker
!            pseudopotential (GTH, local part).

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    REAL(wp) :: alpha_exp_ppl,cexp_ppl,exp_ppl_radius
    INTEGER  :: iexp_ppl,ikind,lppl

!   ---------------------------------------------------------------------------

    DO ikind=1,nkind
      IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN
        alpha_exp_ppl = kind_info(ikind)%gth_potential%alpha_exp_ppl
        exp_ppl_radius = 0.0_wp
        DO iexp_ppl=1,kind_info(ikind)%gth_potential%nexp_ppl
          lppl = 2*(iexp_ppl - 1)
          cexp_ppl = kind_info(ikind)%gth_potential%cexp_ppl(iexp_ppl)
          exp_ppl_radius = MAX(exp_ppl_radius,&
                               exp_radius(lppl,alpha_exp_ppl,eps_exp_ppl,&
                                          cexp_ppl))
        END DO
        kind_info(ikind)%gth_potential%exp_ppl_radius = exp_ppl_radius
      END IF
    END DO

  END SUBROUTINE calculate_exp_ppl_radii

! *****************************************************************************

  SUBROUTINE write_exp_ppl_radii(lunit)

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
      "RADII OF THE EXPONENTIAL FUNCTIONS OF THE GTH PSEUDOPOTENTIAL IN "//&
      TRIM(unit_of_length_name),"Kind","Label","Radius"

    DO ikind=1,nkind
      IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN
        WRITE (UNIT=lunit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
          ikind,kind_info(ikind)%label,&
          kind_info(ikind)%gth_potential%exp_ppl_radius/unit_of_length
      END IF
    END DO

  END SUBROUTINE write_exp_ppl_radii

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

  SUBROUTINE calculate_prj_ppnl_radii

!   Purpose: Calculate the radii of the projector functions of the Goedecker
!            pseudopotential (GTH, non-local part).

!   History: - Creation (11.10.2000, Matthias Krack)

!   ***************************************************************************

    REAL(wp) :: alpha_prj_ppnl,cprj_ppnl,prj_ppnl_radius
    INTEGER  :: ikind,iprj_ppnl,l,lprj_ppnl,nprj_ppnl

!   ---------------------------------------------------------------------------

    DO ikind=1,nkind
      IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN
        prj_ppnl_radius = 0.0_wp
        DO l=0,kind_info(ikind)%gth_potential%lppnl
          alpha_prj_ppnl = kind_info(ikind)%gth_potential%alpha_prj_ppnl(l)
          nprj_ppnl = kind_info(ikind)%gth_potential%nprj_ppnl(l)
          DO iprj_ppnl=1,nprj_ppnl
            lprj_ppnl = l + 2*(iprj_ppnl - 1)
            cprj_ppnl = kind_info(ikind)%gth_potential%cprj_ppnl(iprj_ppnl,l)
            prj_ppnl_radius = MAX(prj_ppnl_radius,&
                                  exp_radius(lprj_ppnl,alpha_prj_ppnl,&
                                             eps_prj_ppnl,cprj_ppnl))
          END DO
        END DO
        kind_info(ikind)%gth_potential%prj_ppnl_radius = prj_ppnl_radius
      END IF
    END DO

  END SUBROUTINE calculate_prj_ppnl_radii

! *****************************************************************************

  SUBROUTINE write_prj_ppnl_radii(lunit)

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
      "RADII OF THE PROJECTOR FUNCTIONS OF THE GTH PSEUDOPOTENTIAL IN "//&
      TRIM(unit_of_length_name),"Kind","Label","Radius"

    DO ikind=1,nkind
      IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN
        WRITE (UNIT=lunit,FMT="(T2,I5,3X,A5,2X,F12.6)")&
          ikind,kind_info(ikind)%label,&
          kind_info(ikind)%gth_potential%prj_ppnl_radius/unit_of_length
      END IF
    END DO

  END SUBROUTINE write_prj_ppnl_radii

! *****************************************************************************

  FUNCTION erfc_radius(alpha,threshold,prefactor) RESULT(radius)

!   Purpose: The radius of the function prefactor*erfc(alpha*r) for a given
!            threshold is calculated. The Newton-Raphson procedure is used to
!            locate the requested radius.

!            g(r) = prefactor*erfc(alpha*r) - threshold = 0

!   History: - Creation (04.08.1999, Matthias Krack)

!   ***************************************************************************

!   alpha    : Prefactor of r.
!   prefactor: Prefactor of the function.
!   radius   : Calculated radius of the function.
!   threshold: Threshold for the radius.

!   delta  : Newton-Raphson step length.
!   dg     : First derivative of the function g(r).
!   epsiter: Convergence criterion of the Newton-Raphson iteration.
!   g      : The function g(r).
!   maxiter: Maximum number of Newton-Raphson iterations.

!   ***************************************************************************

    USE mathconstants, ONLY: pi
    USE termination,   ONLY: stop_program

    REAL(wp), INTENT(IN) :: alpha,prefactor,threshold

    REAL(wp) :: radius

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "FUNCTION erfc_radius (MODULE mathlib)"

    REAL(wp), PARAMETER :: epsiter = 1.0E-12_wp
    INTEGER, PARAMETER  :: maxiter = 50

!   *** Local variables ***

    REAL(wp) :: a,ar,c,d,delta,dg,g,r,t
    INTEGER  :: iter

!   *** Functions ***

    REAL(wp), EXTERNAL :: erfc

!   ---------------------------------------------------------------------------

!   *** Initialize function value ***

    radius = 0.0_wp

!   *** Load and check parameter values ***

    IF (alpha == 0.0_wp) THEN
      CALL stop_program(routine,"The requested exponent is zero")
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

    c = -2.0_wp*a/SQRT(pi)

    t = t/d

!   *** Define an initial radius for the Newton-Raphson iteration ***

    r = SQRT(-LOG(t)/a**2)

!   *** Perform a Newton-Raphson procedure to find the requested radius ***

    DO iter=1,maxiter
      ar = a*r
      g = erfc(ar) - t
      dg = c*exp(-ar**2)
      delta = g/dg
      r = r - SIGN(MIN(ABS(delta),0.5_wp*r),delta)
      IF (ABS(delta) < epsiter) THEN
        radius = r
        RETURN
      END IF
    END DO

    CALL stop_program(routine,"Maximum number of iterations exceeded")

  END FUNCTION erfc_radius

! *****************************************************************************

  FUNCTION exp_radius(l,alpha,threshold,prefactor) RESULT(radius)

!   Purpose: The radius of a primitive Gaussian function for a given threshold
!            is calculated.

!            g(r) = prefactor*r**l*exp(-alpha*r**2) - threshold = 0

!   History: - Creation (26.02.1999, Matthias Krack)

!   ***************************************************************************

!   alpha    : Exponent of the primitive Gaussian function.
!   l        : Angular momentum quantum number l.
!   prefactor: Prefactor of the Gaussian function (e.g. a contraction
!              coefficient).
!   radius   : Calculated radius of the Gaussian function.
!   threshold: Threshold for radius.

!   delta  : Newton-Raphson step length.
!   dg     : Derivative of g.
!   epsiter: Convergence criterion of the Newton-Raphson iteration.
!   g      : The function g(r).
!   maxiter: Maximum number of Newton-Raphson iterations.

!   ***************************************************************************

    USE termination, ONLY: stop_program

    REAL(wp), INTENT(IN) :: alpha,prefactor,threshold
    INTEGER, INTENT(IN)  :: l

    REAL(wp) :: radius

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "FUNCTION exp_radius (MODULE mathlib)"

    REAL(wp), PARAMETER :: epsiter = 1.0E-12_wp
    INTEGER, PARAMETER  :: maxiter = 50

!   *** Local variables ***

    REAL(wp) :: a,ar2,d,delta,dg,g,r,t
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
      g = d*EXP(-ar2)
    ELSE
      g = d*r**l*EXP(-ar2)
    END IF

    IF (t < g) THEN

!     *** Ensure that the iteration starts on ***
!     *** the right-hand side of the maximum  ***

      r = 2.0_wp*r

    ELSE

!     *** The radius is zero, if threshold >= Maximum g(r) ***

      RETURN

    END IF

    t = t/d

    IF (l == 0) THEN

!     *** Analytical expression for l = 0 ***

      radius = SQRT(-LOG(t)/a)

    ELSE

!     *** Find outer radius for l > 0 by Newton-Raphson procedure ***

      DO iter=1,maxiter
        ar2 = a*r*r
        g = r**l*EXP(-ar2) - t
        dg = (REAL(l,wp) - 2.0_wp*ar2)*r**(l-1)*EXP(-ar2)
        delta = g/dg
        r = r - SIGN(MIN(ABS(delta),0.5_wp*r),delta)
        IF (ABS(delta) < epsiter) THEN
          radius = r
          RETURN
        END IF
      END DO

      CALL stop_program(routine,"Maximum number of iterations exceeded")

    END IF

  END FUNCTION exp_radius

! *****************************************************************************

END MODULE interactions
