!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Read xTB parameters.
!> \author JGH (10.2018)
! **************************************************************************************************
MODULE xtb_parameters

   USE basis_set_types,                 ONLY: allocate_sto_basis_set,&
                                              create_gto_from_sto_basis,&
                                              deallocate_sto_basis_set,&
                                              gto_basis_set_type,&
                                              set_sto_basis_set,&
                                              sto_basis_set_type
   USE cp_control_types,                ONLY: xtb_control_type
   USE cp_linked_list_input,            ONLY: cp_sll_val_next,&
                                              cp_sll_val_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                              parser_get_object
   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_list_get,&
                                              section_vals_type
   USE input_val_types,                 ONLY: val_get,&
                                              val_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE periodic_table,                  ONLY: get_ptable_info,&
                                              ptable
   USE physcon,                         ONLY: bohr,&
                                              evolt
   USE string_utilities,                ONLY: remove_word,&
                                              uppercase
   USE xtb_types,                       ONLY: xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER, PRIVATE :: nelem = 106
   !   H                                                                      He
   !   Li Be                                                 B  C  N  O  F    Ne
   !   Na Mg                                                 Al Si P  S  Cl   Ar
   !   K  Ca Sc                Ti V  Cr Mn Fe Co Ni Cu Zn    Ga Ge As Se Br   Kr
   !   Rb Sr Y                 Zr Nb Mo Tc Ru Rh Pd Ag Cd    In Sn Sb Te I    Xe
   !   Cs Ba La Ce-Lu          Hf Ta W  Re Os Ir Pt Au Hg    Tl Pb Bi Po At   Rn
   !   Fr Ra Ac Th Pa U        Np Pu Am Cm Bk Cf Es Fm Md    No Lr Rf Ha 106

!&<
   ! Element Valence
   INTEGER, DIMENSION(0:nelem), &
     PARAMETER, PRIVATE :: zval = (/-1, & !    0
                                     1, 2, & !    2
                                     1, 2, 3, 4, 5, 6, 7, 8, & !   10
                                     1, 2, 3, 4, 5, 6, 7, 8, & !   18
                                     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, & !   36
                                     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, & !   54
                                     1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
                                     4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, & !   86
                                    -1, -1, -1, 4, -1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1/)
!&>

!&<
   ! Element Pauling Electronegativity
   REAL(KIND=dp), DIMENSION(0:nelem), &
      PARAMETER, PRIVATE :: eneg = (/0.00_dp, & ! 0
                                     2.20_dp, 3.00_dp, & ! 2
                                     0.98_dp, 1.57_dp, 2.04_dp, 2.55_dp, 3.04_dp, 3.44_dp, 3.98_dp, 4.50_dp, & ! 10
                                     0.93_dp, 1.31_dp, 1.61_dp, 1.90_dp, 2.19_dp, 2.58_dp, 3.16_dp, 3.50_dp, & ! 18
                                     0.82_dp, 1.00_dp, 1.36_dp, 1.54_dp, 1.63_dp, 1.66_dp, 1.55_dp, 1.83_dp, &
                                     1.88_dp, 1.91_dp, 1.90_dp, 1.65_dp, 1.81_dp, 2.01_dp, 2.18_dp, 2.55_dp, 2.96_dp, 3.00_dp, & ! 36
                                     0.82_dp, 0.95_dp, 1.22_dp, 1.33_dp, 1.60_dp, 2.16_dp, 1.90_dp, 2.20_dp, &
                                     2.28_dp, 2.20_dp, 1.93_dp, 1.69_dp, 1.78_dp, 1.96_dp, 2.05_dp, 2.10_dp, 2.66_dp, 2.60_dp, & ! 54
                                     0.79_dp, 0.89_dp, 1.10_dp, &
                                     1.12_dp, 1.13_dp, 1.14_dp, 1.15_dp, 1.17_dp, 1.18_dp, 1.20_dp, 1.21_dp, &
                                     1.22_dp, 1.23_dp, 1.24_dp, 1.25_dp, 1.26_dp, 1.27_dp, & ! Lanthanides
                                     1.30_dp, 1.50_dp, 2.36_dp, 1.90_dp, 2.20_dp, 2.20_dp, 2.28_dp, 2.54_dp, &
                                     2.00_dp, 2.04_dp, 2.33_dp, 2.02_dp, 2.00_dp, 2.20_dp, 2.20_dp, & ! 86
                                     0.70_dp, 0.89_dp, 1.10_dp, &
                                     1.30_dp, 1.50_dp, 1.38_dp, 1.36_dp, 1.28_dp, 1.30_dp, 1.30_dp, 1.30_dp, &
                                     1.30_dp, 1.30_dp, 1.30_dp, 1.30_dp, 1.30_dp, 1.50_dp, & !  Actinides
                                     1.50_dp, 1.50_dp, 1.50_dp/)
!&>

!&<
   ! Shell occupation
   INTEGER, DIMENSION(1:5, 0:nelem) :: occupation = RESHAPE((/0,0,0,0,0, & ! 0
      1,0,0,0,0,  2,0,0,0,0, & ! 2
      1,0,0,0,0,  2,0,0,0,0,  2,1,0,0,0,  2,2,0,0,0,  2,3,0,0,0,  2,4,0,0,0,  2,5,0,0,0,  2,6,0,0,0, & ! 10
      1,0,0,0,0,  2,0,0,0,0,  2,1,0,0,0,  2,2,0,0,0,  2,3,0,0,0,  2,4,0,0,0,  2,5,0,0,0,  2,6,0,0,0, & ! 18
      1,0,0,0,0,  2,0,0,0,0,  2,0,1,0,0,  2,0,2,0,0,  2,0,3,0,0,  2,0,4,0,0,  2,0,5,0,0,  2,0,6,0,0, &
      2,0,7,0,0,  2,0,8,0,0,  2,0,9,0,0,  2,0,0,0,0, 2,1,0,0,0,  2,2,0,0,0,  2,3,0,0,0,  2,4,0,0,0,  2,5,0,0,0, 2,6,0,0,0, & ! 36
      1,0,0,0,0,  2,0,0,0,0,  2,0,1,0,0,  2,0,2,0,0,  2,0,3,0,0,  2,0,4,0,0,  2,0,5,0,0,  2,0,6,0,0, & !
      2,0,7,0,0,  2,0,8,0,0,  2,0,9,0,0,  2,0,0,0,0, 2,1,0,0,0,  2,2,0,0,0,  2,3,0,0,0,  2,4,0,0,0,  2,5,0,0,0, 2,6,0,0,0, & ! 54
      1,0,0,0,0,  2,0,0,0,0,  2,0,1,0,0, &
      2,0,1,0,0,  2,0,1,0,0,  2,0,1,0,0,  2,0,1,0,0,  2,0,1,0,0,  2,0,1,0,0,  2,0,1,0,0,  2,0,1,0,0, &
      2,0,1,0,0,  2,0,1,0,0,  2,0,1,0,0,  2,0,1,0,0,  2,0,1,0,0,  2,0,1,0,0, & ! Lanthanides
      2,0,2,0,0,  2,0,3,0,0,  2,0,4,0,0,  2,0,5,0,0,  2,0,6,0,0,  2,0,7,0,0,  2,0,8,0,0,  2,0,9,0,0, &
      2,0,0,0,0, 2,1,0,0,0,  2,2,0,0,0,  2,3,0,0,0,  2,4,0,0,0,  2,5,0,0,0,  2,6,0,0,0, & ! 86 (last element defined)
      0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0, & !
      0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0, &
      0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0, & ! Actinides
      0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0/), (/5, nelem+1/))
!&>

!&<
   ! COVALENT RADII
   ! based on "Atomic Radii of the Elements," M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
   ! in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
   ! edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
   ! corrected Nov. 17, 2010 for the 92nd edition.
   REAL(KIND=dp), DIMENSION(0:nelem), &
      PARAMETER, PRIVATE :: crad = (/0.00_dp, & ! 0
                                     0.32_dp, 0.37_dp, & ! 2
                                     1.30_dp, 0.99_dp, 0.84_dp, 0.75_dp, 0.71_dp, 0.64_dp, 0.60_dp, 0.62_dp, & ! 10
                                     1.60_dp, 1.40_dp, 1.24_dp, 1.14_dp, 1.09_dp, 1.04_dp, 1.00_dp, 1.01_dp, & ! 18
                                     2.00_dp, 1.74_dp, 1.59_dp, 1.48_dp, 1.44_dp, 1.30_dp, 1.29_dp, 1.24_dp, &
                                     1.18_dp, 1.17_dp, 1.22_dp, 1.20_dp, 1.23_dp, 1.20_dp, 1.20_dp, 1.18_dp, 1.17_dp, 1.16_dp, & ! 36
                                     2.15_dp, 1.90_dp, 1.76_dp, 1.64_dp, 1.56_dp, 1.46_dp, 1.38_dp, 1.36_dp, &
                                     1.34_dp, 1.30_dp, 1.36_dp, 1.40_dp, 1.42_dp, 1.40_dp, 1.40_dp, 1.37_dp, 1.36_dp, 1.36_dp, & ! 54
                                     2.38_dp, 2.06_dp, 1.94_dp, &
                                     1.84_dp, 1.90_dp, 1.88_dp, 1.86_dp, 1.85_dp, 1.83_dp, 1.82_dp, 1.81_dp, &
                                     1.80_dp, 1.79_dp, 1.77_dp, 1.77_dp, 1.78_dp, 1.74_dp, & ! Lanthanides
                                     1.64_dp, 1.58_dp, 1.50_dp, 1.41_dp, 1.36_dp, 1.32_dp, 1.30_dp, 1.30_dp, &
                                     1.32_dp, 1.44_dp, 1.45_dp, 1.50_dp, 1.42_dp, 1.48_dp, 1.46_dp, & ! 86
                                     2.42_dp, 2.11_dp, 2.01_dp, &
                                     1.90_dp, 1.84_dp, 1.83_dp, 1.80_dp, 1.80_dp, 1.51_dp, 0.96_dp, 1.54_dp, &
                                     1.83_dp, 1.50_dp, 1.50_dp, 1.50_dp, 1.50_dp, 1.50_dp, & !  Actinides
                                     1.50_dp, 1.50_dp, 1.50_dp/)
!&>

!&<
   ! Charge Limits (Mulliken)
   REAL(KIND=dp), DIMENSION(0:nelem), &
      PARAMETER, PRIVATE :: clmt = (/0.00_dp, & ! 0
                                     1.00_dp, 0.50_dp, & ! 2
                                     1.00_dp, 2.00_dp, 3.00_dp, 4.00_dp, 3.00_dp, 2.00_dp, 1.00_dp, 0.50_dp, & ! 10
                                     1.00_dp, 2.00_dp, 3.00_dp, 4.00_dp, 3.00_dp, 2.00_dp, 1.00_dp, 0.50_dp, & ! 18
                                     1.00_dp, 2.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, &
                                     3.00_dp, 3.00_dp, 3.00_dp, 2.00_dp, 2.03_dp, 2.00_dp, 2.00_dp, 2.00_dp, 1.00_dp, 0.50_dp, & ! 36
                                     1.00_dp, 2.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, &
                                     3.00_dp, 3.00_dp, 3.00_dp, 2.00_dp, 2.00_dp, 2.00_dp, 2.00_dp, 2.00_dp, 1.00_dp, 0.50_dp, & ! 54
                                     1.00_dp, 2.00_dp, 3.00_dp, &
                                     3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, &
                                     3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, & ! Lanthanides
                                     3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, &
                                     2.00_dp, 2.00_dp, 2.00_dp, 2.00_dp, 2.00_dp, 1.00_dp, 0.50_dp, & ! 86
                                     1.00_dp, 2.00_dp, 3.00_dp, &
                                     3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, &
                                     3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, & !  Actinides
                                     3.00_dp, 3.00_dp, 3.00_dp/)
!&>

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_parameters'

! *** Public data types ***

   PUBLIC :: xtb_parameters_init, xtb_parameters_read, xtb_parameters_set, init_xtb_basis, &
             xtb_set_kab

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param param ...
!> \param element_symbol ...
!> \param parameter_file_path ...
!> \param parameter_file_name ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE xtb_parameters_init(param, element_symbol, parameter_file_path, parameter_file_name, &
                                  para_env)

      TYPE(xtb_atom_type), POINTER                       :: param
      CHARACTER(LEN=2), INTENT(IN)                       :: element_symbol
      CHARACTER(LEN=*), INTENT(IN)                       :: parameter_file_path, parameter_file_name
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xtb_parameters_init', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=2)                                   :: enam, esym
      CHARACTER(len=default_string_length)               :: aname, filename
      INTEGER                                            :: i, ia, l
      LOGICAL                                            :: at_end, found
      TYPE(cp_parser_type), POINTER                      :: parser

      filename = ADJUSTL(TRIM(parameter_file_path))//ADJUSTL(TRIM(parameter_file_name))
      NULLIFY (parser)
      CALL parser_create(parser, filename, para_env=para_env)
      found = .FALSE.
      DO
         at_end = .FALSE.
         CALL parser_get_next_line(parser, 1, at_end)
         IF (at_end) EXIT
         CALL parser_get_object(parser, aname)
         enam = aname
         esym = element_symbol
         CALL uppercase(enam)
         CALL uppercase(esym)
         IF (enam == esym) THEN
            found = .TRUE.
            CALL parser_get_object(parser, param%eta)
            CALL parser_get_object(parser, param%xgamma)
            CALL parser_get_object(parser, param%alpha)
            CALL parser_get_object(parser, param%zneff)
            DO i = 1, 5
               CALL parser_get_object(parser, aname)
               ia = ICHAR(aname(1:1))
               IF (ia >= 49 .AND. ia <= 57) THEN
                  CALL parser_get_object(parser, param%kpoly(i))
                  CALL parser_get_object(parser, param%kappa(i))
                  CALL parser_get_object(parser, param%hen(i))
                  CALL parser_get_object(parser, param%zeta(i))
                  param%nshell = i
                  param%nval(i) = ia - 48
                  SELECT CASE (aname(2:2))
                  CASE ("s", "S")
                     param%lval(i) = 0
                  CASE ("p", "P")
                     param%lval(i) = 1
                  CASE ("d", "D")
                     param%lval(i) = 2
                  CASE ("f", "F")
                     param%lval(i) = 3
                  CASE DEFAULT
                     CPABORT("xTB PARAMETER ERROR")
                  END SELECT
                  CALL parser_get_next_line(parser, 1, at_end)
                  IF (at_end) EXIT
               ELSE
                  EXIT
               END IF
            END DO
            EXIT
         END IF
      END DO
      IF (found) THEN
         param%typ = "STANDARD"
         param%symbol = element_symbol
         param%defined = .TRUE.
         CALL get_ptable_info(element_symbol, number=ia)
         param%z = ia
         param%aname = ptable(ia)%name
         param%lmax = MAXVAL(param%lval(1:param%nshell))
         param%natorb = 0
         DO i = 1, param%nshell
            l = param%lval(i)
            param%natorb = param%natorb + (2*l + 1)
         END DO
         param%zeff = zval(ia)
      ELSE
         param%defined = .FALSE.
         CALL cp_warn(__LOCATION__, "xTB parameters for element "//element_symbol// &
                      " were not found in the parameter file "//ADJUSTL(TRIM(filename)))
      END IF
      CALL parser_release(parser)

   END SUBROUTINE xtb_parameters_init

! **************************************************************************************************
!> \brief Read atom parameters for xTB Hamiltonian from input file
!> \param param ...
!> \param element_symbol ...
!> \param xtb_section ...
! **************************************************************************************************
   SUBROUTINE xtb_parameters_read(param, element_symbol, xtb_section)

      TYPE(xtb_atom_type), POINTER                       :: param
      CHARACTER(LEN=2), INTENT(IN)                       :: element_symbol
      TYPE(section_vals_type), POINTER                   :: xtb_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xtb_parameters_read', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: label
      CHARACTER(len=20*default_string_length)            :: line_att
      INTEGER                                            :: i, ia, k, l, nshell
      LOGICAL                                            :: explicit, found, is_ok
      TYPE(cp_sll_val_type), POINTER                     :: list
      TYPE(section_vals_type), POINTER                   :: ap_section
      TYPE(val_type), POINTER                            :: val

      !
      ! This could probably be done nicer
      !
      NULLIFY (list, val)
      ap_section => section_vals_get_subs_vals(xtb_section, "ATOM_PARAMETER")
      CALL section_vals_get(ap_section, explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_list_get(ap_section, "_DEFAULT_KEYWORD_", list=list)
         found = .FALSE.
         nshell = 0
         DO
            is_ok = cp_sll_val_next(list, val)
            IF (.NOT. is_ok) EXIT
            CALL val_get(val, c_val=line_att)
            IF (found) THEN
               READ (line_att, *) label
               CALL remove_word(line_att)
               ia = ICHAR(label(1:1))
               IF (ia >= 49 .AND. ia <= 57) THEN
                  nshell = nshell + 1
                  k = nshell
                  param%nval(k) = ia - 48
                  SELECT CASE (label(2:2))
                  CASE ("s", "S")
                     param%lval(k) = 0
                  CASE ("p", "P")
                     param%lval(k) = 1
                  CASE ("d", "D")
                     param%lval(k) = 2
                  CASE ("f", "F")
                     param%lval(k) = 3
                  CASE DEFAULT
                     CPABORT("xTB PARAMETER ERROR")
                  END SELECT
                  !
                  READ (line_att, *) param%kpoly(k)
                  CALL remove_word(line_att)
                  READ (line_att, *) param%kappa(k)
                  CALL remove_word(line_att)
                  READ (line_att, *) param%hen(k)
                  CALL remove_word(line_att)
                  READ (line_att, *) param%zeta(k)
                  CALL remove_word(line_att)
               ELSE
                  EXIT
               END IF
            ELSE
               READ (line_att, *) label
               CALL remove_word(line_att)
               IF (label == element_symbol) THEN
                  found = .TRUE.
                  nshell = nshell + 1
                  k = nshell
                  READ (line_att, *) param%eta
                  CALL remove_word(line_att)
                  READ (line_att, *) param%xgamma
                  CALL remove_word(line_att)
                  READ (line_att, *) param%alpha
                  CALL remove_word(line_att)
                  READ (line_att, *) param%zneff
                  CALL remove_word(line_att)
                  READ (line_att, *) label
                  CALL remove_word(line_att)
                  ia = ICHAR(label(1:1))
                  CPASSERT((ia >= 49 .AND. ia <= 57))
                  param%nval(k) = ia - 48
                  SELECT CASE (label(2:2))
                  CASE ("s", "S")
                     param%lval(k) = 0
                  CASE ("p", "P")
                     param%lval(k) = 1
                  CASE ("d", "D")
                     param%lval(k) = 2
                  CASE ("f", "F")
                     param%lval(k) = 3
                  CASE DEFAULT
                     CPABORT("xTB PARAMETER ERROR")
                  END SELECT
                  !
                  READ (line_att, *) param%kpoly(k)
                  CALL remove_word(line_att)
                  READ (line_att, *) param%kappa(k)
                  CALL remove_word(line_att)
                  READ (line_att, *) param%hen(k)
                  CALL remove_word(line_att)
                  READ (line_att, *) param%zeta(k)
                  CALL remove_word(line_att)
               ENDIF
            END IF
         END DO
         IF (found) THEN
            param%typ = "STANDARD"
            param%symbol = element_symbol
            param%defined = .TRUE.
            CALL get_ptable_info(element_symbol, number=ia)
            param%z = ia
            param%aname = ptable(ia)%name
            param%lmax = MAXVAL(param%lval(1:param%nshell))
            param%natorb = 0
            param%nshell = nshell
            DO i = 1, param%nshell
               l = param%lval(i)
               param%natorb = param%natorb + (2*l + 1)
            END DO
            param%zeff = zval(ia)
         END IF
      END IF

   END SUBROUTINE xtb_parameters_read

! **************************************************************************************************
!> \brief Read atom parameters for xTB Hamiltonian from input file
!> \param param ...
! **************************************************************************************************
   SUBROUTINE xtb_parameters_set(param)

      TYPE(xtb_atom_type), POINTER                       :: param

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xtb_parameters_set', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, is, l, na
      REAL(KIND=dp), DIMENSION(5)                        :: kp

      IF (param%defined) THEN
         ! AO to shell pointer
         ! AO to l-qn pointer
         na = 0
         DO is = 1, param%nshell
            l = param%lval(is)
            DO i = 1, 2*l + 1
               na = na + 1
               param%nao(na) = is
               param%lao(na) = l
            END DO
         END DO
         !
         i = param%z
         ! Electronegativity
         param%electronegativity = eneg(i)
         ! covalent radius
         param%rcov = crad(i)*bohr
         ! shell occupances
         param%occupation(:) = occupation(:, i)
         ! check for consistency
         IF (ABS(param%zeff - SUM(param%occupation)) > 1.E-10_dp) THEN
            CALL cp_abort(__LOCATION__, "Element <"//TRIM(param%aname)//"> has inconsistent shell occupations")
         END IF
         ! orbital energies [evolt] -> [a.u.]
         param%hen = param%hen/evolt
         ! some forgotten scaling parameters (not in orig. paper)
         param%xgamma = 0.1_dp*param%xgamma
         param%kpoly(:) = 0.01_dp*param%kpoly(:)
         param%kappa(:) = 0.1_dp*param%kappa(:)
         ! we have 1/6 g * q**3 (not 1/3)
         param%xgamma = -2.0_dp*param%xgamma
         ! we need kappa l-indexed
         kp(:) = param%kappa(:)
         param%kappa(:) = 0.0_dp
         DO is = 1, param%nshell
            l = param%lval(is)
            IF (param%kappa(l + 1) == 0.0_dp) THEN
               param%kappa(l + 1) = kp(is)
            ELSE
               CPASSERT(ABS(param%kappa(l + 1) - kp(is)) < 1.e-10_dp)
            END IF
         END DO
         ! kx
         IF (param%kx < -10._dp) THEN
            ! use defaults
            SELECT CASE (param%z)
            CASE DEFAULT
               param%kx = 0.0_dp
            CASE (35) ! Br
               param%kx = 0.1_dp*0.381742_dp
            CASE (53) ! I
               param%kx = 0.1_dp*0.321944_dp
            CASE (85) ! At
               param%kx = 0.1_dp*0.220000_dp
            END SELECT
         END IF
         ! chmax
         param%chmax = clmt(i)
      END IF

   END SUBROUTINE xtb_parameters_set

! **************************************************************************************************
!> \brief ...
!> \param param ...
!> \param gto_basis_set ...
!> \param ngauss ...
! **************************************************************************************************
   SUBROUTINE init_xtb_basis(param, gto_basis_set, ngauss)

      TYPE(xtb_atom_type), POINTER                       :: param
      TYPE(gto_basis_set_type), POINTER                  :: gto_basis_set
      INTEGER, INTENT(IN)                                :: ngauss

      CHARACTER(LEN=*), PARAMETER :: routineN = 'init_xtb_basis', routineP = moduleN//':'//routineN

      CHARACTER(LEN=6), DIMENSION(:), POINTER            :: symbol
      INTEGER                                            :: i, nshell
      INTEGER, DIMENSION(:), POINTER                     :: lq, nq
      REAL(KIND=dp), DIMENSION(:), POINTER               :: zet
      TYPE(sto_basis_set_type), POINTER                  :: sto_basis_set

      IF (ASSOCIATED(param)) THEN
         NULLIFY (sto_basis_set)
         CALL allocate_sto_basis_set(sto_basis_set)
         nshell = param%nshell

         ALLOCATE (symbol(1:nshell))
         symbol = ""
         DO i = 1, nshell
            SELECT CASE (param%lval(i))
            CASE (0)
               WRITE (symbol(i), '(I1,A1)') param%nval(i), "S"
            CASE (1)
               WRITE (symbol(i), '(I1,A1)') param%nval(i), "P"
            CASE (2)
               WRITE (symbol(i), '(I1,A1)') param%nval(i), "D"
            CASE (3)
               WRITE (symbol(i), '(I1,A1)') param%nval(i), "F"
            CASE DEFAULT
               CPABORT('BASIS SET OUT OF RANGE (lval)')
            END SELECT
         END DO

         IF (nshell > 0) THEN
            ALLOCATE (nq(nshell), lq(nshell), zet(nshell))
            nq(1:nshell) = param%nval(1:nshell)
            lq(1:nshell) = param%lval(1:nshell)
            zet(1:nshell) = param%zeta(1:nshell)
            CALL set_sto_basis_set(sto_basis_set, name=param%aname, nshell=nshell, symbol=symbol, &
                                   nq=nq, lq=lq, zet=zet)
            CALL create_gto_from_sto_basis(sto_basis_set, gto_basis_set, ngauss=ngauss, ortho=.TRUE.)
         END IF

         ! this will remove the allocated arrays
         CALL deallocate_sto_basis_set(sto_basis_set)
         DEALLOCATE (symbol, nq, lq, zet)

      ELSE
         CPABORT("The pointer param is not associated")
      END IF

   END SUBROUTINE init_xtb_basis

! **************************************************************************************************
!> \brief ...
!> \param za ...
!> \param zb ...
!> \param xtb_control ...
!> \return ...
! **************************************************************************************************
   FUNCTION xtb_set_kab(za, zb, xtb_control) RESULT(kab)

      INTEGER, INTENT(IN)                                :: za, zb
      TYPE(xtb_control_type), INTENT(IN), POINTER        :: xtb_control
      REAL(KIND=dp)                                      :: kab

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xtb_set_kab', routineP = moduleN//':'//routineN

      INTEGER                                            :: j, z
      LOGICAL                                            :: custom

      kab = 1.0_dp
      custom = .FALSE.

      IF (xtb_control%kab_nval .GT. 0) THEN
         DO j = 1, xtb_control%kab_nval
            IF ((za == xtb_control%kab_types(1, j) .AND. &
                 zb == xtb_control%kab_types(2, j)) .OR. &
                (za == xtb_control%kab_types(2, j) .AND. &
                 zb == xtb_control%kab_types(1, j))) THEN
               custom = .TRUE.
               kab = xtb_control%kab_vals(j)
               EXIT
            END IF
         END DO
      END IF

      IF (.NOT. custom) THEN
         IF (za == 1 .OR. zb == 1) THEN
            ! hydrogen
            z = za + zb - 1
            SELECT CASE (z)
            CASE (1)
               kab = 0.96_dp
            CASE (5)
               kab = 0.95_dp
            CASE (7)
               kab = 1.04_dp
            CASE (28)
               kab = 0.90_dp
            CASE (75)
               kab = 0.80_dp
            CASE (78)
               kab = 0.80_dp
            END SELECT
         ELSEIF (za == 5 .OR. zb == 5) THEN
            ! Boron
            z = za + zb - 5
            SELECT CASE (z)
            CASE (15)
               kab = 0.97_dp
            END SELECT
         ELSEIF (za == 7 .OR. zb == 7) THEN
            ! Nitrogen
            z = za + zb - 7
            SELECT CASE (z)
            CASE (14)
               !xtb orig code parameter file
               ! in the paper this is Kab for B-Si
               kab = 1.01_dp
            END SELECT
         ELSEIF (za > 20 .AND. za < 30) THEN
            ! 3d
            IF (zb > 20 .AND. zb < 30) THEN
               ! 3d
               kab = 1.10_dp
            ELSEIF ((zb > 38 .AND. zb < 48) .OR. (zb > 56 .AND. zb < 80)) THEN
               ! 4d/5d/4f
               kab = 0.50_dp*(1.20_dp + 1.10_dp)
            END IF
         ELSEIF ((za > 38 .AND. za < 48) .OR. (za > 56 .AND. za < 80)) THEN
            ! 4d/5d/4f
            IF (zb > 20 .AND. zb < 30) THEN
               ! 3d
               kab = 0.50_dp*(1.20_dp + 1.10_dp)
            ELSEIF ((zb > 38 .AND. zb < 48) .OR. (zb > 56 .AND. zb < 80)) THEN
               ! 4d/5d/4f
               kab = 1.20_dp
            END IF
         END IF
      END IF

   END FUNCTION xtb_set_kab

END MODULE xtb_parameters

