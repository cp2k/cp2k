!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE fft_lib

   USE fft_kinds,                       ONLY: dp
   USE fft_plan,                        ONLY: BWFFT,&
                                              FWFFT,&
                                              fft_plan_type
   USE fftsg_lib,                       ONLY: fftsg_1d,&
                                              fftsg_2d,&
                                              fftsg_3d,&
                                              fftsg_do_cleanup,&
                                              fftsg_do_init,&
                                              fftsg_get_lengths
   USE fftw3_lib,                       ONLY: &
        fft_alloc => fftw_alloc, fft_dealloc => fftw_dealloc, fftw3_1d, fftw3_3d, &
        fftw3_create_plan_1d, fftw3_create_plan_3d, fftw3_destroy_plan, fftw3_do_cleanup, &
        fftw3_do_init, fftw3_get_lengths
   USE spfft_lib,                       ONLY: spfft_1d,&
                                              spfft_2d,&
                                              spfft_3d,&
                                              spfft_create_plan_1d,&
                                              spfft_create_plan_2d,&
                                              spfft_create_plan_3d,&
                                              spfft_destroy_plan
#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fft_lib'

   PUBLIC :: fft_do_cleanup, fft_do_init, fft_get_lengths, fft_create_plan_2d, fft_create_plan_3d
   PUBLIC :: fft_create_plan_1d, fft_1d, fft_2d, fft_3d, fft_destroy_plan
   PUBLIC :: fft_alloc, fft_dealloc
   PUBLIC :: FWFFT, BWFFT, fft_arbitrary_lengths_available

   LOGICAL, SAVE :: use_fftsg_sizes = .TRUE.
   INTEGER, SAVE :: fft_type = 0

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION fft_arbitrary_lengths_available()
      fft_arbitrary_lengths_available = (fft_type == 3) .OR. (fft_type == 4)
   END FUNCTION fft_arbitrary_lengths_available

! **************************************************************************************************
!> \brief Interface to FFT libraries
!> \param fftlib ...
!> \return ...
!> \par History
!>      IAB 09-Jan-2009 : Modified to use fft_plan_type
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
!> \author JGH
! **************************************************************************************************
   FUNCTION fft_library(fftlib) RESULT(flib)

      CHARACTER(len=*), INTENT(IN)                       :: fftlib
      INTEGER                                            :: flib

      SELECT CASE (fftlib)
      CASE DEFAULT
         CPABORT("Unknown FFT library: "//fftlib)
      CASE ("FFTSG")
         flib = 1
      CASE ("FFTW3")
         flib = 3
      CASE ("SPFFT")
         flib = 4
      END SELECT

   END FUNCTION fft_library

! **************************************************************************************************
!> \brief ...
!> \param wisdom_file ...
!> \param plan_style ...
!> \param fftsg_sizes ...
!> \param fftlib ...
! **************************************************************************************************
   SUBROUTINE fft_do_init(wisdom_file, plan_style, fftsg_sizes, fftlib)
      CHARACTER(LEN=*), INTENT(IN)                       :: wisdom_file
      INTEGER, INTENT(IN)                                :: plan_style
      LOGICAL, INTENT(IN)                                :: fftsg_sizes
      CHARACTER(LEN=*), INTENT(IN)                       :: fftlib

      use_fftsg_sizes = fftsg_sizes
      fft_type = fft_library(fftlib)

      SELECT CASE (fft_type)
      CASE DEFAULT
         CPABORT("fft_do_init")
      CASE (1)
         CALL fftsg_do_init()
      CASE (3, 4)
         CALL fftw3_do_init(wisdom_file, plan_style)
      END SELECT

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param wisdom_file ...
!> \param ionode ...
! **************************************************************************************************
   SUBROUTINE fft_do_cleanup(wisdom_file, ionode)
      CHARACTER(LEN=*), INTENT(IN)                       :: wisdom_file
      LOGICAL, INTENT(IN)                                :: ionode

      SELECT CASE (fft_type)
      CASE DEFAULT
         CPABORT("fft_do_cleanup")
      CASE (1)
         CALL fftsg_do_cleanup()
      CASE (3, 4)
         CALL fftw3_do_cleanup(wisdom_file, ionode)
      END SELECT

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param DATA ...
!> \param max_length ...
! **************************************************************************************************
   SUBROUTINE fft_get_lengths(DATA, max_length)
      INTEGER, INTENT(INOUT)                             :: max_length
      INTEGER, DIMENSION(max_length)                     :: DATA

      IF (use_fftsg_sizes) THEN
         CALL fftsg_get_lengths(DATA, max_length)
      ELSE
         SELECT CASE (fft_type)
         CASE DEFAULT
            CPABORT("fft_get_lengths")
         CASE (1)
            CALL fftsg_get_lengths(DATA, max_length)
         CASE (3, 4)
            CALL fftw3_get_lengths(DATA, max_length)
         END SELECT
      END IF

   END SUBROUTINE fft_get_lengths

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param plan ...
!> \param fft_in_place ...
!> \param n ...
!> \param zin ...
!> \param zout ...
! **************************************************************************************************
   SUBROUTINE fft_create_plan_3d(plan, fft_in_place, n, zin, zout)

      TYPE(fft_plan_type), INTENT(INOUT)                 :: plan
      LOGICAL, INTENT(IN)                                :: fft_in_place
      INTEGER, DIMENSION(3), INTENT(IN)                  :: n
      COMPLEX(KIND=dp), DIMENSION(PRODUCT(n)), &
         INTENT(INOUT)                                   :: zin, zout

      plan%fft_in_place = fft_in_place
      plan%n_3d = n
!$    plan%need_alt_plan = .FALSE.

      ! Planning only needed for FFTW3 and SpLA
      IF (fft_type == 3) THEN
         CALL fftw3_create_plan_3d(plan, zin, zout)
         plan%valid = .TRUE.
      ELSE IF (fft_type == 4) THEN
         CALL spfft_create_plan_3d(plan)
         plan%valid = .TRUE.
      END IF

   END SUBROUTINE fft_create_plan_3d

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param plan ...
!> \param trans ...
!> \param n1 ...
!> \param n2 ...
!> \param m ...
!> \param zin ...
!> \param zintermed ...
!> \param zout ...
! **************************************************************************************************
   SUBROUTINE fft_create_plan_2d(plan, trans, n1, n2, m, zin, zintermed, zout)

      TYPE(fft_plan_type), INTENT(INOUT)                 :: plan(2)
      LOGICAL, INTENT(IN)                                :: trans
      INTEGER, INTENT(IN)                                :: n1, n2, m
      COMPLEX(KIND=dp), DIMENSION(n1*n2*m), &
         INTENT(INOUT)                                   :: zin, zintermed, zout

      plan(1)%trans = trans
      plan(1)%n_3d = [n1, n2, m]
      plan(1)%n = n1
      plan(1)%m = m*n2
!$    plan(1)%need_alt_plan = .FALSE.

      plan(2)%trans = trans
      plan(2)%n = n2
      plan(2)%m = m*n1
!$    plan(2)%need_alt_plan = .FALSE.

      IF ((n1 > 0) .AND. (m > 0) .AND. (n2 > 0)) THEN
         ! Planning only needed for FFTW3
         IF (fft_type .EQ. 3) THEN
            CALL fftw3_create_plan_1d(plan(1), zin, zintermed)
            CALL fftw3_create_plan_1d(plan(2), zintermed, zout)
            plan(1)%valid = .TRUE.
            plan(2)%valid = .TRUE.
         ELSE IF (fft_type .EQ. 4) THEN
            CALL spfft_create_plan_2d(plan(1))
            plan(1)%valid = .TRUE.
            plan(2)%valid = .TRUE.
         ELSE
            plan(1)%valid = .FALSE.
            plan(2)%valid = .FALSE.
         END IF
      ELSE
         plan(1)%valid = .FALSE.
         plan(2)%valid = .FALSE.
      END IF

   END SUBROUTINE fft_create_plan_2d

! **************************************************************************************************
!> \brief ...
!> \param plan ...
!> \param fsign ...
!> \param scale ...
!> \param zin ...
!> \param zout ...
!> \param stat ...
! **************************************************************************************************
   SUBROUTINE fft_3d(plan, fsign, scale, zin, zout, stat)
      TYPE(fft_plan_type), INTENT(IN)                    :: plan
      INTEGER, INTENT(IN)                                :: fsign
      REAL(KIND=dp), INTENT(IN)                          :: scale
      COMPLEX(KIND=dp), DIMENSION(plan%n_3d(1), plan%&
         n_3d(2), plan%n_3d(3)), INTENT(INOUT)           :: zin
      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: zout
      INTEGER, INTENT(OUT)                               :: stat

      stat = fsign
      IF (plan%n_3d(1)*plan%n_3d(2)*plan%n_3d(3) > 0) THEN
         SELECT CASE (fft_type)
         CASE DEFAULT
            CPABORT("fft_3d")
         CASE (1)
            CALL fftsg_3d(plan%fft_in_place, fsign, scale, plan%n_3d, zin, zout)
         CASE (3)
            CALL fftw3_3d(plan, fsign, scale, zin, zout, stat)
         CASE (4)
            CALL spfft_3d(plan, fsign, scale, zin, zout, stat)
         END SELECT
      END IF
      ! stat is set to zero on error, -1,+1 are OK
      IF (stat .EQ. 0) THEN
         stat = 1
      ELSE
         stat = 0
      END IF

   END SUBROUTINE fft_3d

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param plan ...
!> \param trans ...
!> \param n ...
!> \param m ...
!> \param zin ...
!> \param zout ...
! **************************************************************************************************
   SUBROUTINE fft_create_plan_1d(plan, trans, n, m, zin, zout)
      TYPE(fft_plan_type), INTENT(INOUT)                 :: plan
      LOGICAL, INTENT(IN)                                :: trans
      INTEGER, INTENT(IN)                                :: n, m
      COMPLEX(KIND=dp), DIMENSION(n*m), INTENT(INOUT)    :: zin, zout

      plan%trans = trans
      plan%n = n
      plan%m = m
!$    plan%need_alt_plan = .FALSE.

      ! Planning only needed for FFTW3
      IF (fft_type .EQ. 3 .AND. (n*m .NE. 0)) THEN
         CALL fftw3_create_plan_1d(plan, zin, zout)
         plan%valid = .TRUE.
      ELSE IF (fft_type .EQ. 4) THEN
         CALL spfft_create_plan_1d(plan)
         plan%valid = .TRUE.
      ELSE
         plan%valid = .FALSE.
      END IF

   END SUBROUTINE fft_create_plan_1d

! **************************************************************************************************
!> \brief ...
!> \param plan ...
! **************************************************************************************************
   SUBROUTINE fft_destroy_plan(plan)
      TYPE(fft_plan_type), INTENT(INOUT)                 :: plan

! Planning only needed for FFTW3

      IF (plan%valid) THEN
         IF (fft_type == 3) THEN
            CALL fftw3_destroy_plan(plan)
            plan%valid = .FALSE.
         ELSE IF (fft_type == 4) THEN
            CALL spfft_destroy_plan(plan)
            plan%valid = .FALSE.
         END IF
      END IF

      plan = fft_plan_type()

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param plan ...
!> \param fsign ...
!> \param zin ...
!> \param zintermed ...
!> \param zout ...
!> \param scale ...
!> \param stat ...
! **************************************************************************************************
   SUBROUTINE fft_2d(plan, fsign, zin, zintermed, zout, scale, stat)
      TYPE(fft_plan_type), INTENT(IN)                    :: plan(2)
      INTEGER, INTENT(IN)                                :: fsign
      COMPLEX(KIND=dp), &
         DIMENSION(PRODUCT(plan(1)%n_3d)), &
         INTENT(INOUT), TARGET                           :: zin, zintermed, zout
      REAL(KIND=dp), INTENT(IN)                          :: scale
      INTEGER, INTENT(OUT)                               :: stat

      INTEGER                                            :: stat2

      IF (plan(1)%n*plan(1)%m > 0 .AND. plan(2)%n*plan(2)%m > 0) THEN
         SELECT CASE (fft_type)
         CASE DEFAULT
            CPABORT("fft_2d")
         CASE (1)
            CALL fftsg_2d(fsign, plan(1)%trans, scale, plan(1)%n_3d, zin, zintermed, zout)
            stat = fsign
         CASE (3)
            IF (fsign == FWFFT) THEN
               CALL fftw3_1d(plan(1), fsign, zin, zintermed, scale, stat)
               CALL fftw3_1d(plan(2), fsign, zintermed, zout, 1.0_dp, stat2)
            ELSE
               CALL fftw3_1d(plan(2), fsign, zin, zintermed, scale, stat)
               CALL fftw3_1d(plan(1), fsign, zintermed, zout, 1.0_dp, stat2)
            END IF
            stat = stat*stat2
         CASE (4)
            CALL spfft_2d(plan(1), fsign, zin, zintermed, zout, scale, stat)
         END SELECT
      END IF
      ! stat is set to zero on error, -1,+1 are OK
      IF (stat .EQ. 0) THEN
         stat = 1
      ELSE
         stat = 0
      END IF

   END SUBROUTINE fft_2d

! **************************************************************************************************
!> \brief ...
!> \param plan ...
!> \param fsign ...
!> \param zin ...
!> \param zout ...
!> \param scale ...
!> \param stat ...
! **************************************************************************************************
   SUBROUTINE fft_1d(plan, fsign, zin, zout, scale, stat)
      TYPE(fft_plan_type), INTENT(IN)                    :: plan
      INTEGER, INTENT(IN)                                :: fsign
      COMPLEX(KIND=dp), DIMENSION(plan%n, plan%m), &
         INTENT(INOUT), TARGET                           :: zin, zout
      REAL(KIND=dp), INTENT(IN)                          :: scale
      INTEGER, INTENT(OUT)                               :: stat

      stat = fsign
      IF (plan%n*plan%m > 0) THEN
         SELECT CASE (fft_type)
         CASE DEFAULT
            CPABORT("fft_1d")
         CASE (1)
            CALL fftsg_1d(fsign, plan%trans, plan%n, plan%m, zin, zout, scale)
         CASE (3)
            CALL fftw3_1d(plan, fsign, zin, zout, scale, stat)
         CASE (4)
            CALL spfft_1d(plan, fsign, zin, zout, scale, stat)
         END SELECT
      END IF
      ! stat is set to zero on error, -1,+1 are OK
      IF (stat .EQ. 0) THEN
         stat = 1
      ELSE
         stat = 0
      END IF

   END SUBROUTINE fft_1d

END MODULE

