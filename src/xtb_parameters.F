!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
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
   USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                              parser_get_object
   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE periodic_table,                  ONLY: get_ptable_info,&
                                              ptable
   USE physcon,                         ONLY: bohr,&
                                              evolt
   USE string_utilities,                ONLY: uppercase
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
                                     1.05_dp, 1.25_dp, & ! 2
                                     1.05_dp, 2.05_dp, 3.00_dp, 4.00_dp, 3.00_dp, 2.00_dp, 1.25_dp, 1.00_dp, & ! 10
                                     1.05_dp, 2.05_dp, 3.00_dp, 4.00_dp, 3.00_dp, 2.00_dp, 1.25_dp, 1.00_dp, & ! 18
                                     1.05_dp, 2.05_dp, 3.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, &
                                     3.50_dp, 3.50_dp, 3.50_dp, 2.50_dp, 2.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, 1.25_dp, 1.00_dp, & ! 36
                                     1.05_dp, 2.05_dp, 3.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, &
                                     3.50_dp, 3.50_dp, 3.50_dp, 2.50_dp, 2.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, 1.25_dp, 1.00_dp, & ! 54
                                     1.05_dp, 2.05_dp, 3.00_dp, &
                                     3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, &
                                     3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, & ! Lanthanides
                                     3.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, &
                                     2.50_dp, 2.50_dp, 3.50_dp, 3.50_dp, 3.50_dp, 1.25_dp, 1.00_dp, & ! 86
                                     1.05_dp, 2.05_dp, 3.00_dp, &
                                     3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, &
                                     3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, 3.00_dp, & !  Actinides
                                     3.00_dp, 3.00_dp, 3.00_dp/)
!&>

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_parameters'

! *** Public data types ***

   PUBLIC :: xtb_parameters_init, xtb_parameters_set, init_xtb_basis, xtb_set_kab
   PUBLIC :: metal, early3d, pp_gfn0

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param param ...
!> \param gfn_type ...
!> \param element_symbol ...
!> \param parameter_file_path ...
!> \param parameter_file_name ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE xtb_parameters_init(param, gfn_type, element_symbol, &
                                  parameter_file_path, parameter_file_name, &
                                  para_env)

      TYPE(xtb_atom_type), POINTER                       :: param
      INTEGER, INTENT(IN)                                :: gfn_type
      CHARACTER(LEN=2), INTENT(IN)                       :: element_symbol
      CHARACTER(LEN=*), INTENT(IN)                       :: parameter_file_path, parameter_file_name
      TYPE(mp_para_env_type), POINTER                    :: para_env

      SELECT CASE (gfn_type)
      CASE (0)
         CALL xtb0_parameters_init(param, element_symbol, parameter_file_path, &
                                   parameter_file_name, para_env)
      CASE (1)
         CALL xtb1_parameters_init(param, element_symbol, parameter_file_path, &
                                   parameter_file_name, para_env)
      CASE (2)
         CPABORT("gfn_type = 2 not yet supported")
      CASE DEFAULT
         CPABORT("Wrong gfn_type")
      END SELECT

   END SUBROUTINE xtb_parameters_init

! **************************************************************************************************
!> \brief ...
!> \param param ...
!> \param element_symbol ...
!> \param parameter_file_path ...
!> \param parameter_file_name ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE xtb0_parameters_init(param, element_symbol, parameter_file_path, parameter_file_name, &
                                   para_env)

      TYPE(xtb_atom_type), POINTER                       :: param
      CHARACTER(LEN=2), INTENT(IN)                       :: element_symbol
      CHARACTER(LEN=*), INTENT(IN)                       :: parameter_file_path, parameter_file_name
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=2)                                   :: esym
      CHARACTER(len=default_string_length)               :: aname, atag, filename
      INTEGER                                            :: i, l, zin, znum
      LOGICAL                                            :: at_end, found
      TYPE(cp_parser_type)                               :: parser

      filename = ADJUSTL(TRIM(parameter_file_path))//ADJUSTL(TRIM(parameter_file_name))
      CALL parser_create(parser, filename, apply_preprocessing=.FALSE., para_env=para_env)
      found = .FALSE.
      znum = 0
      CALL get_ptable_info(element_symbol, znum)
      DO
         at_end = .FALSE.
         CALL parser_get_next_line(parser, 1, at_end)
         IF (at_end) EXIT
         CALL parser_get_object(parser, aname)
         CALL uppercase(aname)
         IF (aname == "$Z") THEN
            CALL parser_get_object(parser, zin)
            IF (zin == znum) THEN
               found = .TRUE.
               DO
                  CALL parser_get_next_line(parser, 1, at_end)
                  IF (at_end) THEN
                     CPABORT("Incomplete xTB parameter file")
                  END IF
                  CALL parser_get_object(parser, aname)
                  CALL uppercase(aname)
                  SELECT CASE (aname)
                  CASE ("AO")
                     CALL parser_get_object(parser, atag)
                     CALL xtb_get_shells(atag, param%nshell, param%nval, param%lval)
                  CASE ("LEV")
                     DO i = 1, param%nshell
                        CALL parser_get_object(parser, param%hen(i))
                     END DO
                  CASE ("EXP")
                     DO i = 1, param%nshell
                        CALL parser_get_object(parser, param%zeta(i))
                     END DO
                  CASE ("EN")
                     CALL parser_get_object(parser, param%en)
                  CASE ("GAM")
                     CALL parser_get_object(parser, param%eta)
                  CASE ("KQAT2")
                     CALL parser_get_object(parser, param%kqat2)
                  CASE ("KCNS")
                     CALL parser_get_object(parser, param%kcn(1))
                     param%kcn(1) = param%kcn(1)*0.1_dp !from orig xtb code
                  CASE ("KCNP")
                     CALL parser_get_object(parser, param%kcn(2))
                     param%kcn(2) = param%kcn(2)*0.1_dp !from orig xtb code
                  CASE ("KCND")
                     CALL parser_get_object(parser, param%kcn(3))
                     param%kcn(3) = param%kcn(3)*0.1_dp !from orig xtb code
                  CASE ("REPA")
                     CALL parser_get_object(parser, param%alpha)
                  CASE ("REPB")
                     CALL parser_get_object(parser, param%zneff)
                  CASE ("POLYS")
                     CALL parser_get_object(parser, param%kpoly(1))
                  CASE ("POLYP")
                     CALL parser_get_object(parser, param%kpoly(2))
                  CASE ("POLYD")
                     CALL parser_get_object(parser, param%kpoly(3))
                  CASE ("KQS")
                     CALL parser_get_object(parser, param%kq(1))
                  CASE ("KQP")
                     CALL parser_get_object(parser, param%kq(2))
                  CASE ("KQD")
                     CALL parser_get_object(parser, param%kq(3))
                  CASE ("XI")
                     CALL parser_get_object(parser, param%xi)
                  CASE ("KAPPA")
                     CALL parser_get_object(parser, param%kappa0)
                  CASE ("ALPG")
                     CALL parser_get_object(parser, param%alpg)
                  CASE ("$END")
                     EXIT
                  CASE DEFAULT
                     CPABORT("Unknown parameter in xTB file")
                  END SELECT
               END DO
            ELSE
               CYCLE
            END IF
            EXIT
         END IF
      END DO
      IF (found) THEN
         param%typ = "STANDARD"
         param%symbol = element_symbol
         param%defined = .TRUE.
         param%z = znum
         param%aname = ptable(znum)%name
         param%lmax = MAXVAL(param%lval(1:param%nshell))
         param%natorb = 0
         DO i = 1, param%nshell
            l = param%lval(i)
            param%natorb = param%natorb + (2*l + 1)
         END DO
         param%zeff = zval(znum)
      ELSE
         esym = element_symbol
         CALL uppercase(esym)
         IF ("X " == esym) THEN
            param%typ = "GHOST"
            param%symbol = element_symbol
            param%defined = .FALSE.
            param%z = 0
            param%aname = "X "
            param%lmax = 0
            param%natorb = 0
            param%nshell = 0
            param%zeff = 0.0_dp
         ELSE
            param%defined = .FALSE.
            CALL cp_warn(__LOCATION__, "xTB parameters for element "//element_symbol// &
                         " were not found in the parameter file "//ADJUSTL(TRIM(filename)))
         END IF
      END IF
      CALL parser_release(parser)

   END SUBROUTINE xtb0_parameters_init

! **************************************************************************************************
!> \brief ...
!> \param param ...
!> \param element_symbol ...
!> \param parameter_file_path ...
!> \param parameter_file_name ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE xtb1_parameters_init(param, element_symbol, parameter_file_path, parameter_file_name, &
                                   para_env)

      TYPE(xtb_atom_type), POINTER                       :: param
      CHARACTER(LEN=2), INTENT(IN)                       :: element_symbol
      CHARACTER(LEN=*), INTENT(IN)                       :: parameter_file_path, parameter_file_name
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=2)                                   :: esym
      CHARACTER(len=default_string_length)               :: aname, atag, filename
      INTEGER                                            :: i, l, zin, znum
      LOGICAL                                            :: at_end, found
      TYPE(cp_parser_type)                               :: parser

      filename = ADJUSTL(TRIM(parameter_file_path))//ADJUSTL(TRIM(parameter_file_name))
      CALL parser_create(parser, filename, apply_preprocessing=.FALSE., para_env=para_env)
      found = .FALSE.
      znum = 0
      CALL get_ptable_info(element_symbol, znum)
      DO
         at_end = .FALSE.
         CALL parser_get_next_line(parser, 1, at_end)
         IF (at_end) EXIT
         CALL parser_get_object(parser, aname)
         CALL uppercase(aname)
         IF (aname == "$Z") THEN
            CALL parser_get_object(parser, zin)
            IF (zin == znum) THEN
               found = .TRUE.
               DO
                  CALL parser_get_next_line(parser, 1, at_end)
                  IF (at_end) THEN
                     CPABORT("Incomplete xTB parameter file")
                  END IF
                  CALL parser_get_object(parser, aname)
                  CALL uppercase(aname)
                  SELECT CASE (aname)
                  CASE ("AO")
                     CALL parser_get_object(parser, atag)
                     CALL xtb_get_shells(atag, param%nshell, param%nval, param%lval)
                  CASE ("LEV")
                     DO i = 1, param%nshell
                        CALL parser_get_object(parser, param%hen(i))
                     END DO
                  CASE ("EXP")
                     DO i = 1, param%nshell
                        CALL parser_get_object(parser, param%zeta(i))
                     END DO
                  CASE ("GAM")
                     CALL parser_get_object(parser, param%eta)
                  CASE ("GAM3")
                     CALL parser_get_object(parser, param%xgamma)
                  CASE ("CXB")
                     CALL parser_get_object(parser, param%kx)
                  CASE ("REPA")
                     CALL parser_get_object(parser, param%alpha)
                  CASE ("REPB")
                     CALL parser_get_object(parser, param%zneff)
                  CASE ("POLYS")
                     CALL parser_get_object(parser, param%kpoly(1))
                  CASE ("POLYP")
                     CALL parser_get_object(parser, param%kpoly(2))
                  CASE ("POLYD")
                     CALL parser_get_object(parser, param%kpoly(3))
                  CASE ("LPARP")
                     CALL parser_get_object(parser, param%kappa(2))
                  CASE ("LPARD")
                     CALL parser_get_object(parser, param%kappa(3))
                  CASE ("$END")
                     EXIT
                  CASE DEFAULT
                     CPABORT("Unknown parameter in xTB file")
                  END SELECT
               END DO
            ELSE
               CYCLE
            END IF
            EXIT
         END IF
      END DO
      IF (found) THEN
         param%typ = "STANDARD"
         param%symbol = element_symbol
         param%defined = .TRUE.
         param%z = znum
         param%aname = ptable(znum)%name
         param%lmax = MAXVAL(param%lval(1:param%nshell))
         param%natorb = 0
         DO i = 1, param%nshell
            l = param%lval(i)
            param%natorb = param%natorb + (2*l + 1)
         END DO
         param%zeff = zval(znum)
      ELSE
         esym = element_symbol
         CALL uppercase(esym)
         IF ("X " == esym) THEN
            param%typ = "GHOST"
            param%symbol = element_symbol
            param%defined = .FALSE.
            param%z = 0
            param%aname = "X "
            param%lmax = 0
            param%natorb = 0
            param%nshell = 0
            param%zeff = 0.0_dp
         ELSE
            param%defined = .FALSE.
            CALL cp_warn(__LOCATION__, "xTB parameters for element "//element_symbol// &
                         " were not found in the parameter file "//ADJUSTL(TRIM(filename)))
         END IF
      END IF
      CALL parser_release(parser)

   END SUBROUTINE xtb1_parameters_init

! **************************************************************************************************
!> \brief Read atom parameters for xTB Hamiltonian from input file
!> \param param ...
! **************************************************************************************************
   SUBROUTINE xtb_parameters_set(param)

      TYPE(xtb_atom_type), POINTER                       :: param

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
         IF (param%en == 0.0_dp) param%en = eneg(i)
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
         ! we need kpoly in shell order
         kp(:) = param%kpoly(:)
         param%kpoly(:) = 0.0_dp
         DO is = 1, param%nshell
            l = param%lval(is)
            param%kpoly(is) = kp(l + 1)
         END DO
         ! kx
         param%kx = 0.1_dp*param%kx
         IF (param%kx < -5._dp) THEN
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

      CHARACTER(LEN=6), DIMENSION(:), POINTER            :: symbol
      INTEGER                                            :: i, nshell
      INTEGER, DIMENSION(:), POINTER                     :: lq, nq
      REAL(KIND=dp), DIMENSION(:), POINTER               :: zet
      TYPE(sto_basis_set_type), POINTER                  :: sto_basis_set

      IF (ASSOCIATED(param)) THEN
         IF (param%defined) THEN
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
         END IF

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

! **************************************************************************************************
!> \brief ...
!> \param atag ...
!> \param nshell ...
!> \param nval ...
!> \param lval ...
!> \return ...
! **************************************************************************************************
   SUBROUTINE xtb_get_shells(atag, nshell, nval, lval)
      CHARACTER(len=*)                                   :: atag
      INTEGER                                            :: nshell
      INTEGER, DIMENSION(:)                              :: nval, lval

      CHARACTER(LEN=1)                                   :: ltag
      CHARACTER(LEN=10)                                  :: aotag
      INTEGER                                            :: i, j

      aotag = ADJUSTL(TRIM(atag))
      nshell = LEN(TRIM(aotag))/2
      DO i = 1, nshell
         j = (i - 1)*2 + 1
         READ (aotag(j:j), FMT="(i1)") nval(i)
         READ (aotag(j + 1:j + 1), FMT="(A1)") ltag
         CALL uppercase(ltag)
         SELECT CASE (ltag)
         CASE ("S")
            lval(i) = 0
         CASE ("P")
            lval(i) = 1
         CASE ("D")
            lval(i) = 2
         CASE DEFAULT
         END SELECT
      END DO

   END SUBROUTINE xtb_get_shells

! **************************************************************************************************
!> \brief ...
!> \param z ...
!> \return ...
! **************************************************************************************************
   FUNCTION metal(z) RESULT(ismetal)
      INTEGER                                            :: z
      LOGICAL                                            :: ismetal

      SELECT CASE (z)
      CASE DEFAULT
         ismetal = .TRUE.
      CASE (1:2, 6:10, 14:18, 32:36, 50:54, 82:86)
         ismetal = .FALSE.
      END SELECT

   END FUNCTION metal

! **************************************************************************************************
!> \brief ...
!> \param z ...
!> \return ...
! **************************************************************************************************
   FUNCTION early3d(z) RESULT(isearly3d)
      INTEGER                                            :: z
      LOGICAL                                            :: isearly3d

      isearly3d = .FALSE.
      IF (z >= 21 .AND. z <= 24) isearly3d = .TRUE.

   END FUNCTION early3d

! **************************************************************************************************
!> \brief ...
!> \param za ...
!> \param zb ...
!> \return ...
! **************************************************************************************************
   FUNCTION pp_gfn0(za, zb) RESULT(pparm)
      INTEGER                                            :: za, zb
      REAL(KIND=dp)                                      :: pparm

      pparm = 1.0_dp
      IF ((za > 20 .AND. za < 30) .OR. (za > 38 .AND. za < 48) .OR. (za > 56 .AND. za < 80)) THEN
         IF ((zb > 20 .AND. zb < 30) .OR. (zb > 38 .AND. zb < 48) .OR. (zb > 56 .AND. zb < 80)) THEN
            pparm = 1.1_dp
            IF (za == 29 .OR. za == 47 .OR. za == 79) THEN
               IF (za == 29 .OR. za == 47 .OR. za == 79) THEN
                  pparm = 0.9_dp
               END IF
            END IF
         END IF
      END IF

   END FUNCTION pp_gfn0

END MODULE xtb_parameters

