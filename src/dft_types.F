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
!!     - new xc code (31.05.02,jgh)
!!     - qs section revised (29.08.02,jgh)
!!     - kg section added (24.09.02,gt)
!!
!!   SOURCE
!******************************************************************************

MODULE dft_types

  USE kinds, ONLY: wp => dp,&
                   wp_size => dp_size

  USE band,           ONLY: band_structure_type
  USE global_types,   ONLY: global_environment_type
  USE pw_grid_types,  ONLY: pw_grid_type
  USE qs_parser,      ONLY: finish_parser,&
                            read_object,&
                            start_parser,&
                            stop_parser,&
                            test_object
  USE termination,    ONLY: stop_memory,&
                            stop_program
  USE xc_functionals, ONLY: xc_get_functional,&
                            xc_print_info, xc_set

  IMPLICIT NONE

  PRIVATE

! *** Control parameters for a KIM-GORDON calculation ***

  TYPE kg_control_type
    CHARACTER(LEN=10)               :: method
    REAL(wp)                        :: eps_core_charge,&
                                       eps_pgf_aux,&
                                       eps_pgf_orb,&
                                       eps_rho_gspace,&
                                       eps_rho_rspace,&
                                       eps_gvg_rspace,&
                                       progression_factor,&
                                       relative_cutoff
    LOGICAL                         :: polarization
    REAL(wp)                        :: cutoff
    REAL(wp), DIMENSION(:), POINTER :: e_cutoff
  END TYPE kg_control_type

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
    REAL(wp)                        :: cutoff
    REAL(wp), DIMENSION(:), POINTER :: e_cutoff
  END TYPE qs_control_type

! *** Control parameters for a WAVE calculation ***

  TYPE wave_control_type
     REAL(wp) :: e_cutoff_wf
     REAL(wp) :: e_cutoff_dual
  END TYPE wave_control_type

! *** Control parameters for a DFT calculation ***

  TYPE xc_grid_type
     LOGICAL                   :: pw_default
     TYPE ( pw_grid_type )     :: pw_grid
     INTEGER                   :: pw_mesh ( 3 )
     REAL ( wp )               :: pw_cutoff
  END TYPE xc_grid_type

  TYPE dft_control_type
     TYPE(band_structure_type)       :: band
     TYPE(kg_control_type)           :: kg_control
     TYPE(qs_control_type)           :: qs_control
     TYPE(wave_control_type)         :: wave_control
     TYPE(xc_grid_type)              :: integration_grid
     CHARACTER(LEN=40), DIMENSION(3) :: functional
     REAL(wp)                        :: density_cut,&
                                        gradient_cut
     INTEGER                         :: nspins,&
                                        charge
     LOGICAL                         :: forces,&
                                        gradient_functional,&
                                        smooth_gradient, &
                                        do_tddfpt_calculation
     LOGICAL, DIMENSION (2,3)        :: ftags
  END TYPE dft_control_type

! *** Public data types ***

  PUBLIC :: dft_control_type, xc_grid_type,&
            kg_control_type, qs_control_type,&
            wave_control_type

! *** Public subroutines ***

  PUBLIC :: read_dft_control,&
            read_kg_control,&
            read_qs_control,&
            read_wave_control,&
            write_dft_control,&
            write_kg_control,&
            write_qs_control,&
            write_wave_control

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE read_dft_control(dft_control,globenv)

!   Purpose: Read the DFT control parameters.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************
!!>----------------------------------------------------------------------------
!!  SECTION: &dft ... &end                                                    !
!!                                                                            !
!!  functional       [ ... ]                                                  !
!!  x-functional     [ ... ]                                                  !
!!  c-functional     [ ... ]                                                  !
!!  ke-functional    [ ... ]                                                  !
!!  density_cutoff   den_cut                                                  !
!!  gradient_cutoff  grad_cut                                                 !
!!  spin             nspins                                                   !
!!  lsd                                                                       !
!!  grid             plane_wave     [cutoff ecut, mesh nx ny nz]              !
!!  charge           charge                                                   !
!!                                                                            !
!!<----------------------------------------------------------------------------

    TYPE(dft_control_type), POINTER           :: dft_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_dft_control (MODULE dft_types)"

!   *** Local variables ***

    CHARACTER(LEN=60) :: keyword, string, str1, str2, str3
    INTEGER           :: istat
    LOGICAL           :: section_found
    REAL ( wp )       :: ecut, mesh(3)

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(dft_control)) THEN
      ALLOCATE (dft_control,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"dft_control",0)
    END IF

!   *** Load the default values ***

    dft_control % functional ( 1 ) = "SLATER"
    dft_control % functional ( 2 ) = "PZ"
    dft_control % functional ( 3 ) = "NONE"
    dft_control % density_cut = 1.0E-10_wp
    dft_control % gradient_cut = 1.0E-8_wp
    dft_control % nspins = 1
    dft_control % ftags = .FALSE.
    dft_control % gradient_functional = .FALSE.
    dft_control % smooth_gradient= .FALSE.
    dft_control % integration_grid % pw_default = .TRUE.
    dft_control % integration_grid % pw_cutoff  = -1._wp
    dft_control % integration_grid % pw_mesh    = -1_wp
    dft_control % do_tddfpt_calculation = .FALSE.
    dft_control % charge = 0

    SELECT CASE (globenv%run_type)
    CASE ("WAVEFUNCTION OPTIMIZATION")
      dft_control%forces = .FALSE.
    CASE DEFAULT
      dft_control%forces = .TRUE.
    END SELECT

!   *** Read the input section ***

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="DFT",&
                      section_found=section_found)

    DO WHILE (section_found.AND.(test_object(newline=.TRUE.) /= "EOS"))

      CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)

      SELECT CASE (TRIM(keyword))
      CASE ("EXCHANGE-CORRELATION-FUNCTIONAL","XC-FUNCTIONAL","XC-FUN",&
            "XCFUN","FUNCTIONAL","FUN")
        CALL read_object(string,lower_to_upper=.TRUE.)
        CALL xc_get_functional ( 0, string, dft_control % functional, &
                                    dft_control % ftags )
      CASE ("EXCHANGE-FUNCTIONAL","X-FUNCTIONAL","X-FUN","XFUN")
        CALL read_object(string,lower_to_upper=.TRUE.)
        CALL xc_get_functional ( 1, string, dft_control % functional, &
                                    dft_control % ftags )
      CASE ("CORRELATION-FUNCTIONAL","C-FUNCTIONAL","C-FUN","CFUN")
        CALL read_object(string,lower_to_upper=.TRUE.)
        CALL xc_get_functional ( 2, string, dft_control % functional, &
                                    dft_control % ftags )
      CASE ("KINETIC-ENERGY-FUNCTIONAL","KE-FUNCTIONAL","KE-FUN","KEFUN")
        CALL read_object(string,lower_to_upper=.TRUE.)
        CALL xc_get_functional ( 3, string, dft_control % functional, &
                                    dft_control % ftags )
      CASE ("DENSITY_CUTOFF")
        CALL read_object(dft_control%density_cut)
      CASE ("GRADIENT_CUTOFF")
        CALL read_object(dft_control%gradient_cut)
      CASE ("GRID")
        CALL read_object(str1,lower_to_upper=.TRUE.)
        SELECT CASE (TRIM(str1))
        CASE ( "PLANE_WAVE", "PW" )
          dft_control % integration_grid % pw_default = .FALSE.
          CALL read_object(str2,lower_to_upper=.TRUE.)
          SELECT CASE (TRIM(str2))
          CASE ( "MESH" )
            CALL read_object(dft_control % integration_grid % pw_mesh(1))
            CALL read_object(dft_control % integration_grid % pw_mesh(2))
            CALL read_object(dft_control % integration_grid % pw_mesh(3))
          CASE ( "CUTOFF" )
            CALL read_object(dft_control % integration_grid % pw_cutoff)
          END SELECT
        END SELECT
      CASE ("NSPINS","SPINS","SPIN")
        CALL read_object(dft_control%nspins)
      CASE ("SPIN_POLARISATION","SPIN_POLARIZATION","LSD")
        dft_control%nspins = 2
      CASE ("CHARGE")
         CALL read_object(dft_control%charge)
      CASE ("EXCITATIONS")
         CALL read_object(string, lower_to_upper=.TRUE.)
         SELECT CASE (TRIM(string))
         CASE ("TDLR","TDDFPT")
            dft_control % do_tddfpt_calculation = .TRUE.
         END SELECT
      END SELECT
    END DO

    CALL finish_parser()

    ! set derived information
    dft_control % gradient_functional = ANY ( dft_control % ftags ( 1, : ) )

  END SUBROUTINE read_dft_control

! *****************************************************************************
!!>----------------------------------------------------------------------------
!!  SECTION: &kg ... &end                                                     !
!!                                                                            !
!!  CUTOFF                 : ecut                                             !
!!  MULTI-GRID             : ngrid e1 e2 e3 ...                               !
!!  EPS_DEFAULT            : eps                                              !
!!  EPS_CORE_CHARGE        : eps                                              !
!!  EPS_GVG_RSPACE         : eps                                              !
!!  EPS_PGF_AUX            : eps                                              !
!!  EPS_PGF_ORB            : eps                                              !
!!  EPS_RHO                : eps                                              !
!!  EPS_RHO_GSPACE         : eps                                              !
!!  EPS_RHO_RSPACE         : eps                                              !
!!  PROGRESSION_FACTOR     : eprog                                            !
!!  RELATIVE_CUTOFF        : relcut                                           !
!!  NGRID_LEVEL            : ngrid                                            !
!!  METHOD                 : [KG_POL,KG_NOPOL]                                !
!!                                                                            !
!!<----------------------------------------------------------------------------

  SUBROUTINE read_kg_control(kg_control,globenv)

!   Purpose: Read the KIM-GORDON control parameters.

!   History: - Creation (24.9.2002)

!   ***************************************************************************

    TYPE(kg_control_type), INTENT(OUT)        :: kg_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_kg_control (MODULE dft_types)"

!   *** Local variables ***

    CHARACTER(LEN=60) :: keyword,string
    CHARACTER(LEN=3)  :: test_result
    REAL(wp)          :: value
    INTEGER           :: igrid_level,istat,ngrid_level
    LOGICAL           :: section_found
    LOGICAL           :: cutoff_set, multgrid_set, ngrid_set

!   ---------------------------------------------------------------------------

!   *** Load the default values ***

    kg_control%method = "KG_NOPOL"
    kg_control%eps_core_charge = 1.0E-12_wp
    kg_control%eps_pgf_aux = 1.0E-6_wp
    kg_control%eps_pgf_orb = 1.0E-6_wp
    kg_control%eps_rho_gspace = 1.0E-8_wp
    kg_control%eps_rho_rspace = 1.0E-8_wp
    kg_control%eps_gvg_rspace = 1.0E-6_wp
    kg_control%progression_factor = 2.0_wp
    kg_control%relative_cutoff = 25.0_wp
    kg_control%polarization = .FALSE.
!
    cutoff_set = .FALSE.
    multgrid_set = .FALSE.
    ngrid_set = .FALSE.

!   *** Read the input section ***

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="KG",&
                      section_found=section_found)

    DO WHILE (section_found.AND.(test_object(newline=.TRUE.) /= "EOS"))
      CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
      SELECT CASE (TRIM(keyword))
      CASE ("DENSITY_CUTOFF","E_CUTOFF","ECUTOFF","CUTOFF","ECUT")
        CALL read_object(value)
        kg_control%cutoff = MAX(0.01_wp,0.5_wp*value)
        cutoff_set = .TRUE.
      CASE ("MULTI_GRID","MULTGRID","MGRID")
        CALL read_object(igrid_level)
        IF ( ngrid_set ) THEN
           IF ( igrid_level /= ngrid_level ) CALL stop_parser(routine,&
              "Inconsistent values for number of multi grids")
        END IF
        ngrid_level = igrid_level
        ngrid_set = .TRUE.
        ALLOCATE (kg_control%e_cutoff(ngrid_level),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"kg_control%e_cutoff",&
                                         ngrid_level*wp_size)
        DO igrid_level=1,ngrid_level
          CALL read_object(value)
          kg_control%e_cutoff(igrid_level) = MAX(0.01_wp,0.5_wp*value)
        END DO
        multgrid_set = .TRUE.
      CASE ("EPS_DEFAULT")
        CALL read_object(value)
        kg_control%eps_core_charge = value
        kg_control%eps_pgf_aux = value
        kg_control%eps_pgf_orb = value
        kg_control%eps_rho_gspace = value
        kg_control%eps_rho_rspace = value
        kg_control%eps_gvg_rspace = value
      CASE ("EPS_CORE_CHARGE")
        CALL read_object(kg_control%eps_core_charge)
      CASE ("EPS_GVG_RSPACE","EPS_GVG")
        CALL read_object(kg_control%eps_gvg_rspace)
      CASE ("EPS_PGF_AUX")
        CALL read_object(kg_control%eps_pgf_aux)
      CASE ("EPS_PGF_ORB")
        CALL read_object(kg_control%eps_pgf_orb)
      CASE ("EPS_RHO")
        CALL read_object(value)
        kg_control%eps_rho_gspace = value
        kg_control%eps_rho_rspace = value
      CASE ("EPS_RHO_GSPACE")
        CALL read_object(kg_control%eps_rho_gspace)
      CASE ("EPS_RHO_RSPACE")
        CALL read_object(kg_control%eps_rho_rspace)
      CASE ("PROGRESSION_FACTOR","PROFAC")
        CALL read_object(value)
        kg_control%progression_factor = MAX(0.01_wp,value)
      CASE ("RELATIVE_CUTOFF","REL_CUTOFF")
        CALL read_object(value)
        kg_control%relative_cutoff = MAX(10.0_wp,value)
      CASE ("NGRID_LEVEL","NGRID")
        CALL read_object(ngrid_level)
        ngrid_set = .TRUE.
      CASE ("METHOD")
        CALL read_object(string,lower_to_upper=.TRUE.)
        SELECT CASE (string)
        CASE ("KG_POL")
          kg_control%method = "KG_POL"
          kg_control%polarization = .TRUE.
        CASE ("KG_NOPOL")
          kg_control%method = "KG_NOPOL"
          kg_control%polarization = .FALSE.
        CASE DEFAULT
          CALL stop_parser(routine,&
                           "Invalid method <"//TRIM(keyword)//&
                           "> specified for QUICKSTEP")
        END SELECT
      END SELECT
    END DO

    CALL finish_parser()

!   set defaults for cutoff and multigrids if not yet done
    IF ( .NOT. cutoff_set ) THEN
      IF ( .NOT. multgrid_set ) THEN
        kg_control%cutoff = MAX(0.01_wp,0.5_wp*20._wp)
      ELSE
        kg_control%cutoff = kg_control%e_cutoff(1)
      END IF
    END IF
    IF ( .NOT. ngrid_set ) ngrid_level = 1
    IF ( .NOT. multgrid_set ) THEN
      ALLOCATE (kg_control%e_cutoff(ngrid_level),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"kg_control%e_cutoff",&
                                       ngrid_level*wp_size)
      kg_control%e_cutoff(1) = kg_control%cutoff
      DO igrid_level=2,ngrid_level
        kg_control%e_cutoff(igrid_level) = kg_control%e_cutoff(igrid_level-1) &
            / kg_control%progression_factor
      END DO
    ELSE
      ! set cutoff to smallest value in multgrid available with >= cutoff
      DO igrid_level=ngrid_level,1,-1
        IF ( kg_control%cutoff <= kg_control%e_cutoff(igrid_level) ) THEN
          kg_control%cutoff = kg_control%e_cutoff(igrid_level)
          EXIT
        END IF
        ! set largest grid value to cutoff
        IF ( igrid_level == 1) THEN
          kg_control%cutoff = kg_control%e_cutoff(1)
        END IF
      END DO
    END IF
    ! check that multigrids are ordered
    DO igrid_level=2,ngrid_level
      IF ( kg_control%e_cutoff(igrid_level) > kg_control%e_cutoff(igrid_level-1) ) THEN
         CALL stop_program (routine,"Multi grids not ordered")
      END IF
    END DO

  END SUBROUTINE read_kg_control

! *****************************************************************************
!!>----------------------------------------------------------------------------
!!  SECTION: &qs ... &end                                                     !
!!                                                                            !
!!  CUTOFF                 : ecut                                             !
!!  MULTI-GRID             : ngrid e1 e2 e3 ...                               !
!!  EPS_DEFAULT            : eps                                              !
!!  EPS_CORE_CHARGE        : eps                                              !
!!  EPS_GVG_RSPACE         : eps                                              !
!!  EPS_PGF_AUX            : eps                                              !
!!  EPS_PGF_ORB            : eps                                              !
!!  EPS_PPL                : eps                                              !
!!  EPS_PPNL               : eps                                              !
!!  EPS_RHO                : eps                                              !
!!  EPS_RHO_GSPACE         : eps                                              !
!!  EPS_RHO_RSPACE         : eps                                              !
!!  PROGRESSION_FACTOR     : eprog                                            !
!!  RELATIVE_CUTOFF        : relcut                                           !
!!  NGRID_LEVEL            : ngrid                                            !
!!  METHOD                 : [GPW, GAPW]                                      !
!!  PAO                                                                       !
!!                                                                            !
!!<----------------------------------------------------------------------------

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
    LOGICAL           :: cutoff_set, multgrid_set, ngrid_set

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
!
    cutoff_set = .FALSE.
    multgrid_set = .FALSE.
    ngrid_set = .FALSE.

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
        qs_control%cutoff = MAX(0.01_wp,0.5_wp*value)
        cutoff_set = .TRUE.
      CASE ("MULTI_GRID","MULTGRID","MGRID")
        CALL read_object(igrid_level)
        IF ( ngrid_set ) THEN
           IF ( igrid_level /= ngrid_level ) CALL stop_parser(routine,&
              "Inconsistent values for number of multi grids")
        END IF
        ngrid_level = igrid_level
        ngrid_set = .TRUE.
        ALLOCATE (qs_control%e_cutoff(ngrid_level),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"qs_control%e_cutoff",&
                                         ngrid_level*wp_size)
        DO igrid_level=1,ngrid_level
          CALL read_object(value)
          qs_control%e_cutoff(igrid_level) = MAX(0.01_wp,0.5_wp*value)
        END DO
        multgrid_set = .TRUE.
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
        ngrid_set = .TRUE.
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

!   set defaults for cutoff and multigrids if not yet done
    IF ( .NOT. cutoff_set ) THEN
      IF ( .NOT. multgrid_set ) THEN
        qs_control%cutoff = MAX(0.01_wp,0.5_wp*320._wp)
      ELSE
        qs_control%cutoff = qs_control%e_cutoff(1)
      END IF
    END IF
    IF ( .NOT. ngrid_set ) ngrid_level = 4
    IF ( .NOT. multgrid_set ) THEN
      ALLOCATE (qs_control%e_cutoff(ngrid_level),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"qs_control%e_cutoff",&
                                       ngrid_level*wp_size)
      qs_control%e_cutoff(1) = qs_control%cutoff
      DO igrid_level=2,ngrid_level
        qs_control%e_cutoff(igrid_level) = qs_control%e_cutoff(igrid_level-1) &
            / qs_control%progression_factor
      END DO
    ELSE
      ! set cutoff to smallest value in multgrid available with >= cutoff
      DO igrid_level=ngrid_level,1,-1
        IF ( qs_control%cutoff <= qs_control%e_cutoff(igrid_level) ) THEN
          qs_control%cutoff = qs_control%e_cutoff(igrid_level)
          EXIT
        END IF
        ! set largest grid value to cutoff
        IF ( igrid_level == 1) THEN
          qs_control%cutoff = qs_control%e_cutoff(1)
        END IF
      END DO
    END IF
    ! check that multigrids are ordered
    DO igrid_level=2,ngrid_level
      IF ( qs_control%e_cutoff(igrid_level) > qs_control%e_cutoff(igrid_level-1) ) THEN
         CALL stop_program (routine,"Multi grids not ordered")
      END IF
    END DO

  END SUBROUTINE read_qs_control

! *****************************************************************************

  SUBROUTINE read_wave_control(wave_control,globenv)

!   Purpose: Read the WAVE control parameters.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************
!!>----------------------------------------------------------------------------
!!  SECTION: &wave ... &end                                                   !
!!                                                                            !
!!  ecut             e_cutoff_wf                                              !
!!  dual             e_cutoff_dual                                            !
!!                                                                            !
!!<----------------------------------------------------------------------------

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

      IF (dft_control%forces) THEN
        WRITE (UNIT=output_unit,FMT="(/,T2,A,T78,A)")&
          "DFT| Force calculation","YES"
      ELSE
        WRITE (UNIT=output_unit,FMT="(/,T2,A,T78,A)")&
          "DFT| Force calculation"," NO"
      END IF

      IF ( dft_control % nspins == 1 ) then
         !do nothing
      ELSE IF ( dft_control % nspins == 2 ) then
         WRITE ( output_unit, '( A, T60, A )' ) " DFT| ","spin polarized system"
      ELSE
         WRITE ( output_unit, '( A, T60, A, T77, I4 )' ) " DFT| ","sets of spin states:",&
            dft_control % nspins
      END IF

      IF ( dft_control % charge /= 0 ) then
         WRITE ( output_unit, '( A, T60, A, T77, I4 )' ) " DFT| ","charge:",&
              dft_control % charge
      END IF

      WRITE (UNIT=output_unit,FMT="(T2,A,T67,E14.5)")&
        "DFT| Cutoffs: density ",dft_control%density_cut,&
        "DFT|          gradient",dft_control%gradient_cut

      IF ( dft_control % smooth_gradient )  WRITE ( output_unit, '( A, T79, A )' ) &
           " DFT| Density gradient smoothing ","ON"

      IF ( dft_control % integration_grid % pw_default ) THEN
         WRITE ( output_unit, '( A, T74, A )' ) " DFT| XC mesh (PW density)", "default"
      ELSE
        IF ( dft_control % integration_grid % pw_mesh ( 1 ) > 0 ) &
              WRITE ( output_unit, '( A, T51, 3I10 )' ) &
              " DFT| XC mesh (PW density)", dft_control % integration_grid % pw_mesh
        IF ( dft_control % integration_grid % pw_cutoff > 0 ) &
              WRITE ( output_unit, '( A, T67, F10.1, "a.u." )' ) &
              " DFT| XC mesh (PW density)", dft_control % integration_grid % pw_cutoff
      END IF

      CALL xc_set ( dft_control % density_cut, dft_control % gradient_cut )
      CALL xc_print_info ( output_unit, dft_control % functional, &
                           dft_control % nspins /= 1 )

    END IF

  END SUBROUTINE write_dft_control

! *****************************************************************************

  SUBROUTINE write_kg_control(kg_control,globenv)

!   Purpose: Write the DFT control parameters to the output unit.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************

    TYPE(kg_control_type), INTENT(IN)         :: kg_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: igrid_level,ngrid_level,output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%dft_control_parameters) THEN
      ngrid_level = SIZE(kg_control%e_cutoff)
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T71,A)")&
        "KG| Method:",ADJUSTR(kg_control%method)
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
        "KG| Number of grid levels:",SIZE(kg_control%e_cutoff)
      IF (ngrid_level == 1) THEN
        WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
          "KG| Density cutoff [a.u.]:",kg_control%e_cutoff(1)
      ELSE
        WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
          "KG| Density cutoff [a.u.]:",kg_control%cutoff
        WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
          "KG| Multi grid cutoff [a.u.]: 1) grid level",kg_control%e_cutoff(1)
        WRITE (UNIT=output_unit,FMT="(T2,A,I3,A,T71,F10.1)")&
          ("KG|                         ",igrid_level,") grid level",&
           kg_control%e_cutoff(igrid_level),&
           igrid_level=2,SIZE(kg_control%e_cutoff))
      END IF
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
        "KG| Grid level progression factor:",kg_control%progression_factor
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
        "KG| Relative density cutoff [a.u.]:",kg_control%relative_cutoff
      WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
        "KG| Interaction thresholds: eps_pgf_orb:",&
        kg_control%eps_pgf_orb,&
        "KG|                         eps_core_charge:",&
        kg_control%eps_core_charge,&
        "KG|                         eps_rho_gspace:",&
        kg_control%eps_rho_gspace,&
        "KG|                         eps_rho_rspace:",&
        kg_control%eps_rho_rspace,&
        "KG|                         eps_gvg_rspace:",&
        kg_control%eps_gvg_rspace  
    END IF

  END SUBROUTINE write_kg_control

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
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T71,A)")&
        "QS| Method:",ADJUSTR(qs_control%method)
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
        "QS| Number of grid levels:",SIZE(qs_control%e_cutoff)
      IF (ngrid_level == 1) THEN
        WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
          "QS| Density cutoff [a.u.]:",qs_control%e_cutoff(1)
      ELSE
        WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
          "QS| Density cutoff [a.u.]:",qs_control%cutoff
        WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
          "QS| Multi grid cutoff [a.u.]: 1) grid level",qs_control%e_cutoff(1)
        WRITE (UNIT=output_unit,FMT="(T2,A,I3,A,T71,F10.1)")&
          ("QS|                         ",igrid_level,") grid level",&
           qs_control%e_cutoff(igrid_level),&
           igrid_level=2,SIZE(qs_control%e_cutoff))
      END IF
      IF (qs_control%pao) THEN
         WRITE (UNIT=output_unit,FMT="(T2,A)") "QS| PAO active"
      END IF
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
        "QS| Grid level progression factor:",qs_control%progression_factor
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
        "QS| Relative density cutoff [a.u.]:",qs_control%relative_cutoff
      WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
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
