!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/dft_types [1.0] *
!!
!!   NAME
!!     dft_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     - qs_control_type extended (18.01.2002,MK)
!!     - added pao (18.04.02,fawzi)
!!
!!   SOURCE
!******************************************************************************

MODULE dft_types

  USE kinds, ONLY: wp => dp,&
                   wp_size => dp_size

  USE band,         ONLY: band_structure_type
  USE global_types, ONLY: global_environment_type
  USE qs_parser,    ONLY: finish_parser,&
                          read_object,&
                          start_parser,&
                          stop_parser,&
                          test_object
  USE termination,  ONLY: stop_memory,&
                          stop_program

  IMPLICIT NONE

  PRIVATE

! *** Control parameters for a QUICKSTEP calculation ***

! eps_pgf_orb: Cutoff value for the interaction of the primitive Gaussian-type
!              functions (primitive basis functions).

  TYPE qs_control_type
!MK PRIVATE would be nice
    CHARACTER(LEN=10)               :: method
    REAL(wp)                        :: eps_core_charge,&
                                       eps_pgf_aux,&
                                       eps_pgf_orb,&
                                       eps_ppl,&
                                       eps_ppnl,&
                                       eps_rho_gspace,&
                                       eps_rho_rspace,&
                                       eps_gvg_rspace,&
                                       progression_factor,&
                                       relative_cutoff
    LOGICAL                         :: gapw,gpw,pao
    REAL(wp), DIMENSION(:), POINTER :: e_cutoff
  END TYPE qs_control_type

! *** Control parameters for a WAVE calculation ***

  TYPE wave_control_type
    REAL(wp) :: e_cutoff_dual,e_cutoff_wf
  END TYPE wave_control_type

! *** Control parameters for a DFT calculation ***

  TYPE dft_control_type
     TYPE(band_structure_type) :: band
     TYPE(qs_control_type)     :: qs_control
     TYPE(wave_control_type)   :: wave_control
     CHARACTER(LEN=40)         :: correlation_functional,&
                                  exchange_functional,&
                                  functional
     REAL(wp)                  :: density_cut,gradient_cut
     INTEGER                   :: nspins
     LOGICAL                   :: gradient_functional,smooth_gradient
     INTEGER, DIMENSION(3)     :: xc_mesh
  END TYPE dft_control_type

! *** Public data types ***

  PUBLIC :: dft_control_type,&
            qs_control_type,&
            wave_control_type

! *** Public subroutines ***

  PUBLIC :: read_dft_control,&
            read_qs_control,&
            read_wave_control,&
            write_dft_control,&
            write_qs_control,&
            write_wave_control

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE read_dft_control(dft_control,globenv)

!   Purpose: Read the DFT control parameters.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************

    TYPE(dft_control_type), POINTER           :: dft_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_dft_control (MODULE dft_types)"

!   *** Local variables ***

    CHARACTER(LEN=60) :: keyword,string
    INTEGER           :: istat
    LOGICAL           :: section_found

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(dft_control)) THEN
      ALLOCATE (dft_control,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"dft_control",0)
    END IF

!   *** Load the default values ***

    dft_control%correlation_functional = "None"
    dft_control%exchange_functional = "None"
    dft_control%functional = "Pade"
    dft_control%density_cut = 1.0E-10_wp
    dft_control%gradient_cut = 1.0E-8_wp
    dft_control%nspins = 1
    dft_control%gradient_functional = .FALSE.
    dft_control%smooth_gradient = .FALSE.
    dft_control%xc_mesh(1:3) = -1

!   *** Read the input section ***

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="DFT",&
                      section_found=section_found)

    DO WHILE (section_found.AND.(test_object(newline=.TRUE.) /= "EOS"))
      CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
      SELECT CASE (TRIM(keyword))
      CASE ("CORRELATION-FUNCTIONAL","C-FUNCTIONAL","C-FUN","CFUN")
        CALL read_object(string,lower_to_upper=.TRUE.)
        SELECT CASE (string)
        CASE ("LEE-YANG-PARR","LYP")
          dft_control%correlation_functional = "LYP"
          dft_control%functional = ""
          dft_control%gradient_functional = .TRUE.
        CASE ("NONE")
          dft_control%correlation_functional = "None"
          dft_control%functional = ""
        CASE ("PBE")
          dft_control%correlation_functional = "PBE"
          dft_control%functional = ""
          dft_control%gradient_functional = .TRUE.
        CASE ("PERDEW86","P86")
          dft_control%correlation_functional = "Perdew86"
          dft_control%functional = ""
          dft_control%gradient_functional = .TRUE.
        CASE ("VWN")
          dft_control%correlation_functional = "VWN"
          dft_control%functional = ""
        CASE DEFAULT
          CALL stop_parser(routine,&
                           "Invalid functional name <"//TRIM(string)//&
                           "> specified for the correlation functional")
        END SELECT
      CASE ("DENSITY_CUTOFF")
        CALL read_object(dft_control%density_cut)
      CASE ("EXCHANGE-CORRELATION-FUNCTIONAL","XC-FUNCTIONAL","XC-FUN",&
            "XCFUN","FUNCTIONAL","FUN")
        CALL read_object(string,lower_to_upper=.TRUE.)
        SELECT CASE (string)
        CASE ("BLYP")
          dft_control%exchange_functional = "Becke88"
          dft_control%correlation_functional = "LYP"
          dft_control%functional = ""
          dft_control%gradient_functional = .TRUE.
        CASE ("BP")
          dft_control%exchange_functional = "Becke88"
          dft_control%correlation_functional = "Perdew86"
          dft_control%functional = ""
          dft_control%gradient_functional = .TRUE.
        CASE ("HCTH/93","HCTH-93","HCTH93")
          dft_control%exchange_functional = "None"
          dft_control%correlation_functional = "None"
          dft_control%functional = "HCTH/93"
          dft_control%gradient_functional = .TRUE.
        CASE ("HCTH/120","HCTH-120","HCTH120")
          dft_control%exchange_functional = "None"
          dft_control%correlation_functional = "None"
          dft_control%functional = "HCTH/120"
          dft_control%gradient_functional = .TRUE.
        CASE ("HCTH/147","HCTH-147","HCTH147")
          dft_control%exchange_functional = "None"
          dft_control%correlation_functional = "None"
          dft_control%functional = "HCTH/147"
          dft_control%gradient_functional = .TRUE.
        CASE ("HCTH/407","HCTH-407","HCTH407")
          dft_control%exchange_functional = "None"
          dft_control%correlation_functional = "None"
          dft_control%functional = "HCTH/407"
          dft_control%gradient_functional = .TRUE.
        CASE ("NONE")
          dft_control%exchange_functional = "None"
          dft_control%correlation_functional = "None"
          dft_control%functional = "None"
        CASE ("PADE")
          dft_control%exchange_functional = "None"
          dft_control%correlation_functional = "None"
          dft_control%functional = "Pade"
        CASE ("PBE")
          dft_control%exchange_functional = "PBE"
          dft_control%correlation_functional = "PBE"
          dft_control%functional = ""
          dft_control%gradient_functional = .TRUE.
        CASE ("PERDEW86","P86")
          dft_control%exchange_functional = "Perdew86"
          dft_control%correlation_functional = "Perdew86"
          dft_control%functional = ""
          dft_control%gradient_functional = .TRUE.
        CASE ("SLATER","XALPHA")
          dft_control%exchange_functional = "Slater"
          dft_control%correlation_functional = "None"
          dft_control%functional = ""
        CASE ("VWN")
          dft_control%exchange_functional = "VWN"
          dft_control%correlation_functional = "VWN"
          dft_control%functional = ""
        CASE DEFAULT
          CALL stop_parser(routine,&
                           "Invalid functional name <"//TRIM(string)//&
                           "> specified for the exchange-correlation "//&
                           " functional")
        END SELECT
      CASE ("EXCHANGE-FUNCTIONAL","X-FUNCTIONAL","X-FUN","XFUN")
        CALL read_object(string,lower_to_upper=.TRUE.)
        SELECT CASE (string)
        CASE ("BECKE88","B88")
          dft_control%exchange_functional = "Becke88"
          dft_control%functional = ""
          dft_control%gradient_functional = .TRUE.
        CASE ("NONE")
          dft_control%exchange_functional = "None"
          dft_control%functional = ""
        CASE ("PBE")
          dft_control%exchange_functional = "PBE"
          dft_control%functional = ""
          dft_control%gradient_functional = .TRUE.
        CASE ("PERDEW86","P86")
          dft_control%exchange_functional = "Perdew86"
          dft_control%functional = ""
          dft_control%gradient_functional = .TRUE.
        CASE ("SLATER","XALPHA")
          dft_control%exchange_functional = "Slater"
          dft_control%functional = ""
        CASE ("VWN")
          dft_control%exchange_functional = "VWN"
          dft_control%functional = ""
        CASE DEFAULT
          CALL stop_parser(routine,&
                           "Invalid functional name <"//TRIM(string)//&
                           "> specified for the exchange functional")
        END SELECT
      CASE ("GRADIENT_CUTOFF")
        CALL read_object(dft_control%gradient_cut)
      CASE ("XC_MESH","MESH")
        CALL read_object(dft_control%xc_mesh(1))
        CALL read_object(dft_control%xc_mesh(2))
        CALL read_object(dft_control%xc_mesh(3))
      CASE ("NSPINS","SPINS","SPIN")
        CALL read_object(dft_control%nspins)
      CASE ("SPIN_POLARISATION","SPIN_POLARIZATION","LSD")
        dft_control%nspins = 2
      END SELECT
    END DO

    CALL finish_parser()

    IF (LEN_TRIM(dft_control%functional) == 0) THEN
      dft_control%functional = TRIM(dft_control%exchange_functional)//"-"//&
                               TRIM(dft_control%correlation_functional)
    END IF

  END SUBROUTINE read_dft_control

! *****************************************************************************

  SUBROUTINE read_qs_control(qs_control,globenv)

!   Purpose: Read the QUICKSTEP control parameters.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_control_type), INTENT(OUT)        :: qs_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_qs_control (MODULE dft_types)"

!   *** Local variables ***

    CHARACTER(LEN=60) :: keyword,string
    CHARACTER(LEN=3)  :: test_result
    REAL(wp)          :: value
    INTEGER           :: igrid_level,istat,ngrid_level
    LOGICAL           :: section_found

!   ---------------------------------------------------------------------------

!   *** Load the default values ***

    qs_control%method = "GPW"
    qs_control%eps_core_charge = 1.0E-12_wp
    qs_control%eps_pgf_aux = 1.0E-6_wp
    qs_control%eps_pgf_orb = 1.0E-6_wp
    qs_control%eps_ppl = 1.0E-12_wp
    qs_control%eps_ppnl = 1.0E-12_wp
    qs_control%eps_rho_gspace = 1.0E-8_wp
    qs_control%eps_rho_rspace = 1.0E-8_wp
    qs_control%eps_gvg_rspace = 1.0E-6_wp
    qs_control%progression_factor = 2.0_wp
    qs_control%relative_cutoff = 25.0_wp
    qs_control%gapw = .FALSE.
    qs_control%gpw = .TRUE.
    qs_control%pao = .FALSE.
    ngrid_level = 4
    ALLOCATE (qs_control%e_cutoff(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"qs_control%e_cutoff",&
                                     ngrid_level*wp_size)
    qs_control%e_cutoff(1) = 160.0_wp
    DO igrid_level=2,ngrid_level
      qs_control%e_cutoff(igrid_level) = qs_control%e_cutoff(igrid_level-1)/&
                                         qs_control%progression_factor
    END DO

!   *** Read the input section ***

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="QS",&
                      section_found=section_found)

    DO WHILE (section_found.AND.(test_object(newline=.TRUE.) /= "EOS"))
      CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
      SELECT CASE (TRIM(keyword))
      CASE ("DENSITY_CUTOFF","E_CUTOFF","ECUTOFF","CUTOFF","ECUT")
        CALL read_object(value)
        qs_control%e_cutoff(1) = MAX(0.01_wp,0.5_wp*value)
        DO igrid_level=2,ngrid_level
          test_result = test_object()
          IF ((test_result == "INT").OR.(test_result == "FLT")) THEN
            CALL read_object(value)
            qs_control%e_cutoff(igrid_level) = MAX(0.01_wp,0.5_wp*value)
          ELSE
            qs_control%e_cutoff(igrid_level) =&
              qs_control%e_cutoff(igrid_level-1)/qs_control%progression_factor
          END IF
        END DO
      CASE ("EPS_DEFAULT")
        CALL read_object(value)
        qs_control%eps_core_charge = value
        qs_control%eps_pgf_aux = value
        qs_control%eps_pgf_orb = value
        qs_control%eps_ppl = value
        qs_control%eps_ppnl = value
        qs_control%eps_rho_gspace = value
        qs_control%eps_rho_rspace = value
        qs_control%eps_gvg_rspace = value
      CASE ("EPS_CORE_CHARGE")
        CALL read_object(qs_control%eps_core_charge)
      CASE ("EPS_GVG_RSPACE","EPS_GVG")
        CALL read_object(qs_control%eps_gvg_rspace)
      CASE ("EPS_PGF_AUX")
        CALL read_object(qs_control%eps_pgf_aux)
      CASE ("EPS_PGF_ORB")
        CALL read_object(qs_control%eps_pgf_orb)
      CASE ("EPS_PPL")
        CALL read_object(qs_control%eps_ppl)
      CASE ("EPS_PPNL")
        CALL read_object(qs_control%eps_ppnl)
      CASE ("EPS_RHO")
        CALL read_object(value)
        qs_control%eps_rho_gspace = value
        qs_control%eps_rho_rspace = value
      CASE ("EPS_RHO_GSPACE")
        CALL read_object(qs_control%eps_rho_gspace)
      CASE ("EPS_RHO_RSPACE")
        CALL read_object(qs_control%eps_rho_rspace)
      CASE ("PROGRESSION_FACTOR","PROFAC")
        CALL read_object(value)
        qs_control%progression_factor = MAX(0.01_wp,value)
      CASE ("RELATIVE_CUTOFF","REL_CUTOFF")
        CALL read_object(value)
        qs_control%relative_cutoff = MAX(10.0_wp,value)
      CASE ("NGRID_LEVEL","NGRID")
        CALL read_object(ngrid_level)
        value = qs_control%e_cutoff(1)
        DEALLOCATE (qs_control%e_cutoff,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"qs_control%e_cutoff")
        ALLOCATE (qs_control%e_cutoff(ngrid_level),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"qs_control%e_cutoff",&
                                         ngrid_level*wp_size)
        qs_control%e_cutoff(1) = value
        DO igrid_level=2,ngrid_level
          qs_control%e_cutoff(igrid_level) =&
            qs_control%e_cutoff(igrid_level-1)/qs_control%progression_factor
        END DO
      CASE ("METHOD")
        CALL read_object(string,lower_to_upper=.TRUE.)
        SELECT CASE (string)
        CASE ("GAPW")
          qs_control%method = "GAPW"
          qs_control%gapw = .TRUE.
        CASE ("GPW")
          qs_control%method = "GPW"
          qs_control%gpw = .TRUE.
        CASE DEFAULT
          CALL stop_parser(routine,&
                           "Invalid method <"//TRIM(keyword)//&
                           "> specified for QUICKSTEP")
        END SELECT
      CASE ("PAO")
        qs_control%pao=.TRUE.
      END SELECT
    END DO

    CALL finish_parser()

  END SUBROUTINE read_qs_control

! *****************************************************************************

  SUBROUTINE read_wave_control(wave_control,globenv)

!   Purpose: Read the WAVE control parameters.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************

    TYPE(wave_control_type), INTENT(OUT)      :: wave_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_wave_control (MODULE dft_types)"

!   *** Local variables ***

    CHARACTER(LEN=60) :: keyword,string
    LOGICAL           :: section_found

!   ---------------------------------------------------------------------------

!   *** Load the default values ***

    wave_control%e_cutoff_wf = -1.0_wp
    wave_control%e_cutoff_dual = 4.0_wp

!   *** Read the input section ***

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="WAVE",&
                      section_found=section_found)

    DO WHILE (section_found.AND.(test_object(newline=.TRUE.) /= "EOS"))
      CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
      SELECT CASE (TRIM(keyword))
      CASE ("CUTOFF_DUAL","DUAL")
        CALL read_object(wave_control%e_cutoff_dual)
      CASE ("CUTOFF_WF","CUTOFF")
        CALL read_object(wave_control%e_cutoff_wf)
      END SELECT
    END DO

    CALL finish_parser()

  END SUBROUTINE read_wave_control

! *****************************************************************************

  SUBROUTINE write_dft_control(dft_control,globenv)

!   Purpose: Write the DFT control parameters to the output unit.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************

    TYPE(dft_control_type), POINTER           :: dft_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%dft_control_parameters) THEN
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T40,A)")&
        "DFT| Exchange-correlation functional",&
        ADJUSTR(dft_control%functional)
      WRITE (UNIT=output_unit,FMT="(T2,A,T70,I10)")&
        "DFT| Sets of spin states:",dft_control%nspins
      WRITE (UNIT=output_unit,FMT="(T2,A,T67,E13.5)")&
        "DFT| Cutoffs: density ",dft_control%density_cut,&
        "DFT|          gradient",dft_control%gradient_cut
      IF (dft_control%xc_mesh(1) == -1 ) THEN
        WRITE (UNIT=output_unit,FMT="(T2,A,T73,A)")&
          "DFT| XC mesh (PW density)","default"
      ELSE
        WRITE (UNIT=output_unit,FMT="(T2,A,T50,3I10)")&
          "DFT| XC mesh (PW density)",dft_control%xc_mesh(1:3)
      END IF
    END IF

  END SUBROUTINE write_dft_control

! *****************************************************************************

  SUBROUTINE write_qs_control(qs_control,globenv)

!   Purpose: Write the DFT control parameters to the output unit.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_control_type), INTENT(IN)         :: qs_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: igrid_level,ngrid_level,output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%dft_control_parameters) THEN
      ngrid_level = SIZE(qs_control%e_cutoff)
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T70,A)")&
        "QS| Method:",ADJUSTR(qs_control%method)
      WRITE (UNIT=output_unit,FMT="(T2,A,T70,I10)")&
        "QS| Number of grid levels:",SIZE(qs_control%e_cutoff)
      IF (ngrid_level == 1) THEN
        WRITE (UNIT=output_unit,FMT="(T2,A,T70,F10.1)")&
          "QS| Density cutoff [a.u.]:",qs_control%e_cutoff(1)
      ELSE
        WRITE (UNIT=output_unit,FMT="(T2,A,T70,F10.1)")&
          "QS| Density cutoff [a.u.]: 1) grid level",qs_control%e_cutoff(1)
        WRITE (UNIT=output_unit,FMT="(T2,A,I3,A,T70,F10.1)")&
          ("QS|                      ",igrid_level,") grid level",&
           qs_control%e_cutoff(igrid_level),&
           igrid_level=2,SIZE(qs_control%e_cutoff))
      END IF
      IF (qs_control%pao) THEN
         WRITE (UNIT=output_unit,FMT="(T2,A)") "QS| PAO active"
      END IF
      WRITE (UNIT=output_unit,FMT="(T2,A,T70,F10.1)")&
        "QS| Grid level progression factor:",qs_control%progression_factor
      WRITE (UNIT=output_unit,FMT="(T2,A,T70,F10.1)")&
        "QS| Relative density cutoff [a.u.]:",qs_control%relative_cutoff
      WRITE (UNIT=output_unit,FMT="(T2,A,T72,ES8.1)")&
        "QS| Interaction thresholds: eps_pgf_orb:",&
        qs_control%eps_pgf_orb,&
        "QS|                         eps_core_charge:",&
        qs_control%eps_core_charge,&
        "QS|                         eps_rho_gspace:",&
        qs_control%eps_rho_gspace,&
        "QS|                         eps_rho_rspace:",&
        qs_control%eps_rho_rspace,&
        "QS|                         eps_gvg_rspace:",&
        qs_control%eps_gvg_rspace,&
        "QS|                         eps_ppl:",&
        qs_control%eps_ppl,&
        "QS|                         eps_ppnl:",&
        qs_control%eps_ppnl
    END IF

  END SUBROUTINE write_qs_control

! *****************************************************************************

  SUBROUTINE write_wave_control(wave_control,globenv)

!   Purpose: Write the DFT control parameters to the output unit.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************

    TYPE(wave_control_type), INTENT(IN)       :: wave_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%dft_control_parameters) THEN
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T71,F10.3)")&
        "WAVE| Plane wave cutoff [a.u.]:",wave_control%e_cutoff_wf
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.2)")&
        "WAVE| Density cutoff (DUAL):",wave_control%e_cutoff_dual
     END IF

  END SUBROUTINE write_wave_control

! *****************************************************************************

END MODULE dft_types
