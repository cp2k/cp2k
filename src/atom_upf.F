!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE atom_upf
   USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                              parser_get_object,&
                                              parser_test_next_token
   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release
   USE kinds,                           ONLY: default_string_length,&
                                              dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   TYPE atom_upfpot_type
      CHARACTER(LEN=2)                             :: symbol
      CHARACTER(LEN=default_string_length)         :: pname
      INTEGER, DIMENSION(0:3)                       :: econf
      REAL(dp)                                      :: zion
      CHARACTER(LEN=default_string_length)         :: version
      ! <INFO>
      INTEGER                                       :: maxinfo = 100
      CHARACTER(LEN=default_string_length), DIMENSION(100) &
         :: info
      ! <HEADER>
      CHARACTER(LEN=default_string_length)         :: generated
      CHARACTER(LEN=default_string_length)         :: author
      CHARACTER(LEN=default_string_length)         :: date
      CHARACTER(LEN=default_string_length)         :: comment
      CHARACTER(LEN=4)                             :: pseudo_type
      CHARACTER(LEN=15)                            :: relativistic
      CHARACTER(LEN=default_string_length)         :: functional
      LOGICAL                                       :: is_ultrasoft = .FALSE.
      LOGICAL                                       :: is_paw = .FALSE.
      LOGICAL                                       :: is_coulomb = .FALSE.
      LOGICAL                                       :: has_so = .FALSE.
      LOGICAL                                       :: has_wfc = .FALSE.
      LOGICAL                                       :: has_gipaw = .FALSE.
      LOGICAL                                       :: paw_as_gipaw = .FALSE.
      LOGICAL                                       :: core_correction = .FALSE.
      REAL(dp)                                      :: total_psenergy
      REAL(dp)                                      :: wfc_cutoff
      REAL(dp)                                      :: rho_cutoff
      INTEGER                                       :: l_max
      INTEGER                                       :: l_max_rho
      INTEGER                                       :: l_local
      INTEGER                                       :: mesh_size
      INTEGER                                       :: number_of_wfc
      INTEGER                                       :: number_of_proj
      ! <MESH>
      REAL(dp)                                      :: dx
      REAL(dp)                                      :: xmin
      REAL(dp)                                      :: rmax
      REAL(dp)                                      :: zmesh
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: r, rab
      ! <NLCC>
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: rho_nlcc
      ! <LOCAL>
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: vlocal
      ! <NONLOCAL>
      REAL(dp), DIMENSION(:, :), ALLOCATABLE           :: dion
      REAL(dp), DIMENSION(:, :), ALLOCATABLE           :: beta
      INTEGER, DIMENSION(:), ALLOCATABLE              :: lbeta
      ! <SEMILOCAL>
      REAL(dp), DIMENSION(:, :), ALLOCATABLE           :: vsemi
   END TYPE atom_upfpot_type

   PRIVATE
   PUBLIC  :: atom_read_upf, atom_upfpot_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_upf'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param pot ...
!> \param upf_filename ...
! **************************************************************************************************
   SUBROUTINE atom_read_upf(pot, upf_filename)

      TYPE(atom_upfpot_type)                             :: pot
      CHARACTER(len=*), INTENT(IN)                       :: upf_filename

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_read_upf', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: nametag
      INTEGER                                            :: ib, ntag
      LOGICAL                                            :: at_end
      TYPE(cp_parser_type), POINTER                      :: parser

      ntag = 0
      NULLIFY (parser)
      CALL parser_create(parser, upf_filename)
      DO
         at_end = .FALSE.
         CALL parser_get_next_line(parser, 1, at_end)
         IF (at_end) EXIT
         CALL parser_get_object(parser, nametag, lower_to_upper=.TRUE.)
         IF (nametag(1:1) /= "<") CYCLE
         IF (ntag == 0) THEN
            ! we are looking for UPF tag
            IF (nametag(2:4) == "UPF") THEN
               CALL parser_get_object(parser, nametag, lower_to_upper=.TRUE.)
               ! read UPF file version
               CALL parser_get_object(parser, nametag, lower_to_upper=.TRUE.)
               pot%version = TRIM(nametag)
               CPASSERT(nametag(2:6) == "2.0.1")
               CALL parser_get_object(parser, nametag, lower_to_upper=.TRUE.)
               CPASSERT(nametag(1:1) == ">")
               ntag = 1
            ENDIF
         ELSE IF (ntag == 1) THEN
            ! we are looking for 1st level tags
            IF (nametag(2:8) == "PP_INFO") THEN
               CPASSERT(nametag(9:9) == ">")
               CALL upf_info_section(parser, pot)
            ELSEIF (nametag(2:10) == "PP_HEADER") THEN
               IF (.NOT. (nametag(11:11) == ">")) THEN
                  CALL upf_header_option(parser, pot)
               END IF
            ELSEIF (nametag(2:8) == "PP_MESH") THEN
               IF (.NOT. (nametag(9:9) == ">")) THEN
                  CALL upf_mesh_option(parser, pot)
               END IF
               CALL upf_mesh_section(parser, pot)
            ELSEIF (nametag(2:8) == "PP_NLCC") THEN
               IF (nametag(9:9) == ">") THEN
                  CALL upf_nlcc_section(parser, pot, .FALSE.)
               ELSE
                  CALL upf_nlcc_section(parser, pot, .TRUE.)
               END IF
            ELSEIF (nametag(2:9) == "PP_LOCAL") THEN
               IF (nametag(10:10) == ">") THEN
                  CALL upf_local_section(parser, pot, .FALSE.)
               ELSE
                  CALL upf_local_section(parser, pot, .TRUE.)
               END IF
            ELSEIF (nametag(2:12) == "PP_NONLOCAL") THEN
               CPASSERT(nametag(13:13) == ">")
               CALL upf_nonlocal_section(parser, pot)
            ELSEIF (nametag(2:13) == "PP_SEMILOCAL") THEN
               CALL upf_semilocal_section(parser, pot)
            ELSEIF (nametag(2:9) == "PP_PSWFC") THEN
               ! skip section for now
            ELSEIF (nametag(2:11) == "PP_RHOATOM") THEN
               ! skip section for now
            ELSEIF (nametag(2:7) == "PP_PAW") THEN
               ! skip section for now
            ELSEIF (nametag(2:6) == "/UPF>") THEN
               EXIT
            ENDIF
         END IF
      END DO
      CALL parser_release(parser)

      CPASSERT(ntag > 0)

      ! rescale projectors
      IF (ALLOCATED(pot%beta)) THEN
         DO ib = 1, pot%number_of_proj
            pot%beta(:, ib) = pot%beta(:, ib)/pot%r(:)
         END DO
      END IF

   END SUBROUTINE atom_read_upf
! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
! **************************************************************************************************
   SUBROUTINE upf_info_section(parser, pot)
      TYPE(cp_parser_type), POINTER                      :: parser
      TYPE(atom_upfpot_type)                             :: pot

      CHARACTER(LEN=default_string_length)               :: line, string
      INTEGER                                            :: icount, iline
      LOGICAL                                            :: at_end

      icount = 0
      DO
         CALL parser_get_next_line(parser, 1, at_end)
         CPASSERT(.NOT. at_end)
         iline = parser%buffer%present_line_number
         line = TRIM(parser%buffer%input_lines(iline))
         CALL parser_get_object(parser, string)
         IF (string(1:10) == "</PP_INFO>") EXIT
         icount = icount+1
         IF (icount > pot%maxinfo) CYCLE
         pot%info(icount) = line
      END DO
      pot%maxinfo = icount

   END SUBROUTINE upf_info_section

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
! **************************************************************************************************
   SUBROUTINE upf_header_option(parser, pot)
      TYPE(cp_parser_type), POINTER                      :: parser
      TYPE(atom_upfpot_type)                             :: pot

      CHARACTER(LEN=default_string_length)               :: line, string
      INTEGER                                            :: il
      LOGICAL                                            :: at_end

      DO
         IF (parser_test_next_token(parser) == "EOL") THEN
            CALL parser_get_next_line(parser, 1, at_end)
            CPASSERT(.NOT. at_end)
         END IF
         CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
         IF (string == "/>") EXIT
         SELECT CASE (string)
         CASE ("GENERATED")
            CALL parser_get_object(parser, pot%generated)
         CASE ("AUTHOR")
            CALL parser_get_object(parser, pot%author)
         CASE ("DATE")
            CALL parser_get_object(parser, pot%date)
         CASE ("COMMENT")
            CALL parser_get_object(parser, pot%comment)
         CASE ("ELEMENT")
            CALL parser_get_object(parser, line)
            pot%symbol = line(2:3)
         CASE ("PSEUDO_TYPE")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            pot%pseudo_type = line(2:il-1)
         CASE ("RELATIVISTIC")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            pot%relativistic = line(2:il-1)
         CASE ("IS_ULTRASOFT")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%is_ultrasoft
         CASE ("IS_PAW")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%is_paw
         CASE ("IS_COULOMB")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%is_coulomb
         CASE ("HAS_SO")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%has_so
         CASE ("HAS_WFC")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%has_wfc
         CASE ("HAS_GIPAW")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%has_gipaw
         CASE ("PAW_AS_GIPAW")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%paw_as_gipaw
         CASE ("CORE_CORRECTION")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%core_correction
         CASE ("FUNCTIONAL")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            pot%functional = line(2:il-1)
         CASE ("Z_VALENCE")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%zion
         CASE ("TOTAL_PSENERGY")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%total_psenergy
         CASE ("WFC_CUTOFF")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%wfc_cutoff
         CASE ("RHO_CUTOFF")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%rho_cutoff
         CASE ("L_MAX")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%l_max
         CASE ("L_MAX_RHO")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%l_max_rho
         CASE ("L_LOCAL")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%l_local
         CASE ("MESH_SIZE")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%mesh_size
         CASE ("NUMBER_OF_WFC")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%number_of_wfc
         CASE ("NUMBER_OF_PROJ")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%number_of_proj
         CASE DEFAULT
            CPASSERT(.FALSE.)
         END SELECT
      END DO

   END SUBROUTINE upf_header_option

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
! **************************************************************************************************
   SUBROUTINE upf_mesh_option(parser, pot)
      TYPE(cp_parser_type), POINTER                      :: parser
      TYPE(atom_upfpot_type)                             :: pot

      CHARACTER(LEN=default_string_length)               :: line, string
      INTEGER                                            :: il, jj
      LOGICAL                                            :: at_end

      DO
         IF (parser_test_next_token(parser) == "EOL") THEN
            CALL parser_get_next_line(parser, 1, at_end)
            CPASSERT(.NOT. at_end)
         END IF
         CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
         IF (string == ">") EXIT
         SELECT CASE (string)
         CASE ("DX")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%dx
         CASE ("XMIN")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%xmin
         CASE ("RMAX")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%rmax
         CASE ("MESH")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) jj
            CPASSERT(pot%mesh_size == jj)
         CASE ("ZMESH")
            CALL parser_get_object(parser, line)
            line = ADJUSTL(line)
            il = LEN_TRIM(line)
            READ (line(2:il-1), *) pot%zmesh
         CASE DEFAULT
            CPASSERT(.FALSE.)
         END SELECT

      END DO

   END SUBROUTINE upf_mesh_option

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
! **************************************************************************************************
   SUBROUTINE upf_mesh_section(parser, pot)
      TYPE(cp_parser_type), POINTER                      :: parser
      TYPE(atom_upfpot_type)                             :: pot

      CHARACTER(LEN=default_string_length)               :: line, string, string2
      INTEGER                                            :: icount, il, m, mc, ms
      LOGICAL                                            :: at_end

      DO
         CALL parser_get_next_line(parser, 1, at_end)
         CPASSERT(.NOT. at_end)
         CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
         SELECT CASE (string)
         CASE ("<PP_R")
            m = pot%mesh_size
            ms = pot%mesh_size
            mc = 1
            IF (string(6:6) /= ">") THEN
               ! options
               DO
                  IF (parser_test_next_token(parser) == "EOL") THEN
                     CALL parser_get_next_line(parser, 1, at_end)
                     CPASSERT(.NOT. at_end)
                  END IF
                  CALL parser_get_object(parser, string2, lower_to_upper=.TRUE.)
                  IF (string2 == ">") EXIT
                  SELECT CASE (string2)
                  CASE ("TYPE")
                     CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
                     CPASSERT(line == '"REAL"')
                  CASE ("SIZE")
                     CALL parser_get_object(parser, line)
                     line = ADJUSTL(line)
                     il = LEN_TRIM(line)
                     READ (line(2:il-1), *) ms
                     CPASSERT(ms <= m)
                  CASE ("COLUMNS")
                     CALL parser_get_object(parser, line)
                     line = ADJUSTL(line)
                     il = LEN_TRIM(line)
                     READ (line(2:il-1), *) mc
                  CASE DEFAULT
                     CPASSERT(.FALSE.)
                  END SELECT
               END DO
            END IF
            ALLOCATE (pot%r(m))
            pot%r = 0.0_dp
            icount = 1
            DO
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL parser_get_next_line(parser, 1, at_end)
                  CPASSERT(.NOT. at_end)
               ELSE IF (parser_test_next_token(parser) == "FLT") THEN
                  CALL parser_get_object(parser, pot%r(icount))
                  icount = icount+1
               END IF
               IF (icount > ms) EXIT
            END DO
         CASE ("<PP_RAB")
            IF (string(6:6) /= ">") THEN
               ! options
               DO
                  IF (parser_test_next_token(parser) == "EOL") THEN
                     CALL parser_get_next_line(parser, 1, at_end)
                     CPASSERT(.NOT. at_end)
                  END IF
                  CALL parser_get_object(parser, string2, lower_to_upper=.TRUE.)
                  IF (string2 == ">") EXIT
                  SELECT CASE (string2)
                  CASE ("TYPE")
                     CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
                     CPASSERT(line == '"REAL"')
                  CASE ("SIZE")
                     CALL parser_get_object(parser, line)
                     line = ADJUSTL(line)
                     il = LEN_TRIM(line)
                     READ (line(2:il-1), *) ms
                     CPASSERT(ms <= m)
                  CASE ("COLUMNS")
                     CALL parser_get_object(parser, line)
                     line = ADJUSTL(line)
                     il = LEN_TRIM(line)
                     READ (line(2:il-1), *) mc
                  CASE DEFAULT
                     CPASSERT(.FALSE.)
                  END SELECT
               END DO
            END IF
            ALLOCATE (pot%rab(m))
            pot%rab = 0.0_dp
            icount = 1
            DO
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL parser_get_next_line(parser, 1, at_end)
                  CPASSERT(.NOT. at_end)
               ELSE IF (parser_test_next_token(parser) == "FLT") THEN
                  CALL parser_get_object(parser, pot%rab(icount))
                  icount = icount+1
               END IF
               IF (icount > ms) EXIT
            END DO
         CASE ("</PP_MESH>")
            EXIT
         CASE DEFAULT
            !
         END SELECT
      END DO

   END SUBROUTINE upf_mesh_section

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
!> \param options ...
! **************************************************************************************************
   SUBROUTINE upf_nlcc_section(parser, pot, options)
      TYPE(cp_parser_type), POINTER                      :: parser
      TYPE(atom_upfpot_type)                             :: pot
      LOGICAL, INTENT(IN)                                :: options

      CHARACTER(LEN=default_string_length)               :: line, string
      INTEGER                                            :: icount, il, m, mc, ms
      LOGICAL                                            :: at_end

      m = pot%mesh_size
      ms = m
      mc = 1
      IF (options) THEN
         DO
            IF (parser_test_next_token(parser) == "EOL") THEN
               CALL parser_get_next_line(parser, 1, at_end)
               CPASSERT(.NOT. at_end)
            END IF
            CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
            IF (string == ">") EXIT
            SELECT CASE (string)
            CASE ("TYPE")
               CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
               CPASSERT(line == '"REAL"')
            CASE ("SIZE")
               CALL parser_get_object(parser, line)
               line = ADJUSTL(line)
               il = LEN_TRIM(line)
               READ (line(2:il-1), *) ms
               CPASSERT(ms <= m)
            CASE ("COLUMNS")
               CALL parser_get_object(parser, line)
               line = ADJUSTL(line)
               il = LEN_TRIM(line)
               READ (line(2:il-1), *) mc
            CASE DEFAULT
               CPASSERT(.FALSE.)
            END SELECT
         END DO
      END IF

      ALLOCATE (pot%rho_nlcc(m))
      pot%rho_nlcc = 0.0_dp
      icount = 1
      DO
         IF (parser_test_next_token(parser) == "EOL") THEN
            CALL parser_get_next_line(parser, 1, at_end)
            CPASSERT(.NOT. at_end)
         ELSE IF (parser_test_next_token(parser) == "FLT") THEN
            CALL parser_get_object(parser, pot%rho_nlcc(icount))
            icount = icount+1
         END IF
         IF (icount > ms) EXIT
      END DO

      CALL parser_get_next_line(parser, 1, at_end)
      CPASSERT(.NOT. at_end)
      CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
      CPASSERT(string == "</PP_NLCC>")

   END SUBROUTINE upf_nlcc_section

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
!> \param options ...
! **************************************************************************************************
   SUBROUTINE upf_local_section(parser, pot, options)
      TYPE(cp_parser_type), POINTER                      :: parser
      TYPE(atom_upfpot_type)                             :: pot
      LOGICAL, INTENT(IN)                                :: options

      CHARACTER(LEN=default_string_length)               :: line, string
      INTEGER                                            :: icount, il, m, mc, ms
      LOGICAL                                            :: at_end

      m = pot%mesh_size
      ms = m
      mc = 1
      IF (options) THEN
         DO
            IF (parser_test_next_token(parser) == "EOL") THEN
               CALL parser_get_next_line(parser, 1, at_end)
               CPASSERT(.NOT. at_end)
            END IF
            CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
            IF (string == ">") EXIT
            SELECT CASE (string)
            CASE ("TYPE")
               CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
               CPASSERT(line == '"REAL"')
            CASE ("SIZE")
               CALL parser_get_object(parser, line)
               line = ADJUSTL(line)
               il = LEN_TRIM(line)
               READ (line(2:il-1), *) ms
               CPASSERT(ms <= m)
            CASE ("COLUMNS")
               CALL parser_get_object(parser, line)
               line = ADJUSTL(line)
               il = LEN_TRIM(line)
               READ (line(2:il-1), *) mc
            CASE DEFAULT
               CPASSERT(.FALSE.)
            END SELECT
         END DO
      END IF

      ALLOCATE (pot%vlocal(m))
      pot%vlocal = 0.0_dp
      icount = 1
      DO
         IF (parser_test_next_token(parser) == "EOL") THEN
            CALL parser_get_next_line(parser, 1, at_end)
            CPASSERT(.NOT. at_end)
         ELSE IF (parser_test_next_token(parser) == "FLT") THEN
            CALL parser_get_object(parser, pot%vlocal(icount))
            icount = icount+1
         END IF
         IF (icount > ms) EXIT
      END DO

      ! Ry -> Hartree
      pot%vlocal = 0.5_dp*pot%vlocal

      CALL parser_get_next_line(parser, 1, at_end)
      CPASSERT(.NOT. at_end)
      CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
      CPASSERT(string == "</PP_LOCAL>")

   END SUBROUTINE upf_local_section

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
! **************************************************************************************************
   SUBROUTINE upf_nonlocal_section(parser, pot)
      TYPE(cp_parser_type), POINTER                      :: parser
      TYPE(atom_upfpot_type)                             :: pot

      CHARACTER(LEN=default_string_length)               :: line, string
      INTEGER                                            :: i1, i2, ibeta, icount, il, la, m, mc, &
                                                            ms, nbeta
      LOGICAL                                            :: at_end

      m = pot%mesh_size
      nbeta = pot%number_of_proj
      ALLOCATE (pot%dion(nbeta, nbeta), pot%beta(m, nbeta), pot%lbeta(nbeta))
      pot%dion = 0.0_dp
      pot%beta = 0.0_dp
      pot%lbeta = -1

      ibeta = 0
      DO
         CALL parser_get_next_line(parser, 1, at_end)
         CPASSERT(.NOT. at_end)
         CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
         IF (string(1:8) == "<PP_BETA") THEN
            ms = m
            ibeta = ibeta+1
            i1 = ibeta
            la = 0
            CPASSERT(ibeta <= nbeta)
            DO
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL parser_get_next_line(parser, 1, at_end)
                  CPASSERT(.NOT. at_end)
               END IF
               CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
               IF (string == ">") EXIT
               SELECT CASE (string)
               CASE ("TYPE")
                  CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
                  CPASSERT(line == '"REAL"')
               CASE ("SIZE")
                  CALL parser_get_object(parser, line)
                  line = ADJUSTL(line)
                  il = LEN_TRIM(line)
                  READ (line(2:il-1), *) ms
                  CPASSERT(ms <= m)
               CASE ("COLUMNS")
                  CALL parser_get_object(parser, line)
                  line = ADJUSTL(line)
                  il = LEN_TRIM(line)
                  READ (line(2:il-1), *) mc
               CASE ("INDEX")
                  CALL parser_get_object(parser, line)
                  line = ADJUSTL(line)
                  il = LEN_TRIM(line)
                  READ (line(2:il-1), *) i1
                  CPASSERT(i1 <= nbeta)
               CASE ("ANGULAR_MOMENTUM")
                  CALL parser_get_object(parser, line)
                  line = ADJUSTL(line)
                  il = LEN_TRIM(line)
                  READ (line(2:il-1), *) la
               CASE ("LABEL")
                  CALL parser_get_object(parser, line)
                  ! not used currently
               CASE ("CUTOFF_RADIUS_INDEX")
                  CALL parser_get_object(parser, line)
                  ! not used currently
               CASE ("CUTOFF_RADIUS")
                  CALL parser_get_object(parser, line)
                  ! not used currently
               CASE ("ULTRASOFT_CUTOFF_RADIUS")
                  CALL parser_get_object(parser, line)
                  ! not used currently
               CASE DEFAULT
                  CPASSERT(.FALSE.)
               END SELECT
            END DO
            pot%lbeta(i1) = la
            icount = 1
            DO
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL parser_get_next_line(parser, 1, at_end)
                  CPASSERT(.NOT. at_end)
               ELSE IF (parser_test_next_token(parser) == "FLT") THEN
                  CALL parser_get_object(parser, pot%beta(icount, i1))
                  icount = icount+1
               END IF
               IF (icount > ms) EXIT
            END DO
         ELSE IF (string(1:7) == "<PP_DIJ") THEN
            ms = nbeta*nbeta
            DO
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL parser_get_next_line(parser, 1, at_end)
                  CPASSERT(.NOT. at_end)
               END IF
               CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
               IF (string == ">") EXIT
               SELECT CASE (string)
               CASE ("TYPE")
                  CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
                  CPASSERT(line == '"REAL"')
               CASE ("SIZE")
                  CALL parser_get_object(parser, line)
                  line = ADJUSTL(line)
                  il = LEN_TRIM(line)
                  READ (line(2:il-1), *) ms
                  CPASSERT(ms <= m)
               CASE ("COLUMNS")
                  CALL parser_get_object(parser, line)
                  line = ADJUSTL(line)
                  il = LEN_TRIM(line)
                  READ (line(2:il-1), *) mc
               CASE DEFAULT
                  CPASSERT(.FALSE.)
               END SELECT
            END DO
            icount = 1
            DO
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL parser_get_next_line(parser, 1, at_end)
                  CPASSERT(.NOT. at_end)
               ELSE IF (parser_test_next_token(parser) == "FLT") THEN
                  i1 = (icount-1)/nbeta+1
                  i2 = MOD(icount-1, nbeta)+1
                  CALL parser_get_object(parser, pot%dion(i1, i2))
                  icount = icount+1
               END IF
               IF (icount > ms) EXIT
            END DO
         ELSE IF (string(1:7) == "<PP_QIJ") THEN
            ! skip this option
         ELSE IF (string(1:14) == "</PP_NONLOCAL>") THEN
            EXIT
         END IF
      END DO

      ! change units and scaling, beta is still r*beta
      pot%dion = 2.0_dp*pot%dion
      pot%beta = 0.5_dp*pot%beta

   END SUBROUTINE upf_nonlocal_section

! **************************************************************************************************
!> \brief ...
!> \param parser ...
!> \param pot ...
! **************************************************************************************************
   SUBROUTINE upf_semilocal_section(parser, pot)
      TYPE(cp_parser_type), POINTER                      :: parser
      TYPE(atom_upfpot_type)                             :: pot

      CHARACTER(LEN=default_string_length)               :: line, string
      INTEGER                                            :: i1, ib, icount, il, la, lmax, m, mc, ms
      LOGICAL                                            :: at_end

      m = pot%mesh_size
      lmax = pot%l_max
      ALLOCATE (pot%vsemi(m, lmax+1))
      pot%vsemi = 0.0_dp

      ib = 0
      DO
         CALL parser_get_next_line(parser, 1, at_end)
         CPASSERT(.NOT. at_end)
         CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
         SELECT CASE (string)
         CASE ("<PP_VNL")
            ms = m
            ib = ib+1
            i1 = ib
            la = 0
            CPASSERT(ib <= lmax+1)
            DO
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL parser_get_next_line(parser, 1, at_end)
                  CPASSERT(.NOT. at_end)
               END IF
               CALL parser_get_object(parser, string, lower_to_upper=.TRUE.)
               IF (string == ">") EXIT
               SELECT CASE (string)
               CASE ("TYPE")
                  CALL parser_get_object(parser, line, lower_to_upper=.TRUE.)
                  CPASSERT(line == '"REAL"')
               CASE ("SIZE")
                  CALL parser_get_object(parser, line)
                  line = ADJUSTL(line)
                  il = LEN_TRIM(line)
                  READ (line(2:il-1), *) ms
                  CPASSERT(ms <= m)
               CASE ("COLUMNS")
                  CALL parser_get_object(parser, line)
                  line = ADJUSTL(line)
                  il = LEN_TRIM(line)
                  READ (line(2:il-1), *) mc
               CASE ("L")
                  CALL parser_get_object(parser, line)
                  line = ADJUSTL(line)
                  il = LEN_TRIM(line)
                  READ (line(2:il-1), *) la
               CASE DEFAULT
                  CPASSERT(.FALSE.)
               END SELECT
            END DO
            i1 = la+1
            icount = 1
            DO
               IF (parser_test_next_token(parser) == "EOL") THEN
                  CALL parser_get_next_line(parser, 1, at_end)
                  CPASSERT(.NOT. at_end)
               ELSE IF (parser_test_next_token(parser) == "FLT") THEN
                  CALL parser_get_object(parser, pot%vsemi(icount, i1))
                  icount = icount+1
               END IF
               IF (icount > ms) EXIT
            END DO
         CASE ("</PP_SEMILOCAL>")
            EXIT
         CASE DEFAULT
            !
         END SELECT
      END DO
      ! Ry -> Hartree
      pot%vsemi = 0.5_dp*pot%vsemi

   END SUBROUTINE upf_semilocal_section

! **************************************************************************************************

END MODULE atom_upf
