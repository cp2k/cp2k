!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE fft_lib
   #:include "../pw_types.fypp"

   USE fft_kinds,                       ONLY: dp
   USE fft_plan,                        ONLY: BWFFT,&
                                              FWFFT,&
                                              fft_plan_type
   USE mathconstants, ONLY: z_zero
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
        USE cp_array_utils, ONLY: cp_1d_i_a_type
#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fft_lib'

   PUBLIC :: fft_do_cleanup, fft_do_init, fft_get_lengths, fft_create_plan_2d, fft_create_plan_3d
   PUBLIC :: fft_create_plan_1d, fft_1d, fft_2d, fft_3d, fft_destroy_plan
   PUBLIC :: fft_alloc, fft_dealloc
   PUBLIC :: FWFFT, BWFFT, fft_arbitrary_lengths_available
   PUBLIC :: fft_scatter, fft_gather

   LOGICAL, SAVE :: use_fftsg_sizes = .TRUE.
   INTEGER, SAVE :: fft_type = 0

   INTEGER, PARAMETER, PRIVATE :: do_fft_sg = 1
   INTEGER, PARAMETER, PRIVATE :: do_fft_fftw3 = 3

   INTERFACE fft_gather
      #:for kind in pw_kinds
         #:if kind[1]=="1"
            MODULE PROCEDURE fft_gather_p_${kind}$
         #:endif
      #:endfor
      #:for kind, kind2 in pw_kinds2
         #:if kind[1]=="1" and kind2[1]=="3"
            MODULE PROCEDURE fft_gather_s_${kind}$_${kind2}$
         #:endif
      #:endfor
   END INTERFACE

   INTERFACE fft_scatter
      #:for kind in pw_kinds
         #:if kind[1]=="1"
            MODULE PROCEDURE fft_scatter_p_${kind}$
         #:endif
      #:endfor
      #:for kind, kind2 in pw_kinds2
         #:if kind[1]=="1" and kind2[1]=="3"
            MODULE PROCEDURE fft_scatter_s_${kind}$_${kind2}$
         #:endif
      #:endfor
   END INTERFACE

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION fft_arbitrary_lengths_available()
      fft_arbitrary_lengths_available = (fft_type == do_fft_fftw3)
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
         CPABORT("Unknown FFT library")
      CASE ("FFTSG")
         flib = do_fft_sg
      CASE ("FFTW3")
         flib = do_fft_fftw3
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
         CPABORT("Unknown FFT library")
      CASE (do_fft_sg)
         CALL fftsg_do_init()
      CASE (do_fft_fftw3)
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
         CPABORT("Unknown FFT library")
      CASE (do_fft_sg)
         CALL fftsg_do_cleanup()
      CASE (do_fft_fftw3)
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
            CPABORT("Unknown FFT library")
         CASE (do_fft_sg)
            CALL fftsg_get_lengths(DATA, max_length)
         CASE (do_fft_fftw3)
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

      ! Planning only needed for FFTW3
      IF (fft_type == do_fft_fftw3) THEN
         CALL fftw3_create_plan_3d(plan, zin, zout)
         plan%valid = .TRUE.
      ELSE IF (fft_type == do_fft_sg) THEN
         ! Nothing to set up
      ELSE
         CPABORT("Unknown FFT library")
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
!> \param zout ...
! **************************************************************************************************
   SUBROUTINE fft_create_plan_2d(plan, trans, n1, n2, m, zin, zout)

      TYPE(fft_plan_type), INTENT(INOUT)                 :: plan(2)
      LOGICAL, INTENT(IN)                                :: trans
      INTEGER, INTENT(IN)                                :: n1, n2, m
      COMPLEX(KIND=dp), DIMENSION(n1*n2*m), &
         INTENT(INOUT)                                   :: zin, zout

      plan(1)%trans = trans
      plan(1)%n_3d = [n1, n2, m]
      plan(1)%n = n1
      plan(1)%m = m*n2
!$    plan(1)%need_alt_plan = .FALSE.

      plan(2)%trans = trans
      plan(2)%n_3d = [n1, n2, m]
      plan(2)%n = n2
      plan(2)%m = m*n1
!$    plan(2)%need_alt_plan = .FALSE.

      IF ((n1 > 0) .AND. (m > 0) .AND. (n2 > 0)) THEN
         ! Planning only needed for FFTW3
         IF (fft_type .EQ. do_fft_fftw3) THEN
            CALL fftw3_create_plan_1d(plan(1), zin, zout)
            CALL fftw3_create_plan_1d(plan(2), zin, zout)
            plan(1)%valid = .TRUE.
            plan(2)%valid = .TRUE.
         ELSE IF (fft_type == do_fft_sg) THEN
            plan(1)%valid = .FALSE.
            plan(2)%valid = .FALSE.
         ELSE
            CPABORT("Unknown FFT library")
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
            CPABORT("Unknown FFT library")
         CASE (do_fft_sg)
            CALL fftsg_3d(plan%fft_in_place, fsign, scale, plan%n_3d, zin, zout)
         CASE (do_fft_fftw3)
            CALL fftw3_3d(plan, fsign, scale, zin, zout, stat)
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
      IF (fft_type .EQ. do_fft_fftw3 .AND. (n*m .NE. 0)) THEN
         CALL fftw3_create_plan_1d(plan, zin, zout)
         plan%valid = .TRUE.
      ELSE IF (fft_type == do_fft_sg) THEN
         plan%valid = .FALSE.
      ELSE
         CPABORT("Unknown FFT library")
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
         IF (fft_type == do_fft_fftw3) THEN
            CALL fftw3_destroy_plan(plan)
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
            CPABORT("Unknown FFT library")
         CASE (do_fft_sg)
            CALL fftsg_2d(fsign, plan(1)%trans, scale, plan(1)%n_3d, zin, zintermed, zout)
            stat = fsign
         CASE (do_fft_fftw3)
            IF (fsign == FWFFT) THEN
               CALL fftw3_1d(plan(1), fsign, zin, zintermed, scale, stat)
               CALL fftw3_1d(plan(2), fsign, zintermed, zout, 1.0_dp, stat2)
            ELSE
               CALL fftw3_1d(plan(2), fsign, zin, zintermed, scale, stat)
               CALL fftw3_1d(plan(1), fsign, zintermed, zout, 1.0_dp, stat2)
            END IF
            stat = stat*stat2
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
            CPABORT("Unknown FFT library")
         CASE (do_fft_sg)
            CALL fftsg_1d(fsign, plan%trans, plan%n, plan%m, zin, zout, scale)
         CASE (do_fft_fftw3)
            CALL fftw3_1d(plan, fsign, zin, zout, scale, stat)
         END SELECT
      END IF
      ! stat is set to zero on error, -1,+1 are OK
      IF (stat .EQ. 0) THEN
         stat = 1
      ELSE
         stat = 0
      END IF

   END SUBROUTINE fft_1d
   

   #:for kind, type in pw_list
      #:if kind[1]=="1"
! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param c ...
!> \param scale ...
! **************************************************************************************************
         SUBROUTINE fft_gather_p_${kind}$ (cout, cin, map, ghat, yzq, scale)

            ${type}$, CONTIGUOUS, INTENT(INOUT)                       :: cout
            COMPLEX(KIND=dp), DIMENSION(:, :), CONTIGUOUS, INTENT(IN)      :: cin
            TYPE(cp_1d_i_a_type), DIMENSION(3, 2), INTENT(IN)          :: map
      INTEGER, DIMENSION(:, :), INTENT(IN) :: ghat
      INTEGER, DIMENSION(:, :), INTENT(IN) :: yzq
            REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

            CHARACTER(len=*), PARAMETER                        :: routineN = 'fft_gather_p'

            INTEGER                                            :: gpt, handle, l, m, mn, n

            CALL timeset(routineN, handle)

            ASSOCIATE (mapl => map(1, 1)%array, mapm => map(2, 1)%array, mapn => map(3, 1)%array)

               IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(l, m, mn, n) SHARED(cin, cout, scale, ghat, yzq)
                  DO gpt = 1, SIZE(cout)
                     l = mapl(ghat(1, gpt)) + 1
                     m = mapm(ghat(2, gpt)) + 1
                     n = mapn(ghat(3, gpt)) + 1
                     mn = yzq(m, n)
                     cout(gpt) = scale*${type2type("cin(l, mn)", "c3d", kind)}$
                  END DO
!$OMP END PARALLEL DO
               ELSE
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(l, m, mn, n) SHARED(cin, cout, ghat, yzq)
                  DO gpt = 1, SIZE(cout)
                     l = mapl(ghat(1, gpt)) + 1
                     m = mapm(ghat(2, gpt)) + 1
                     n = mapn(ghat(3, gpt)) + 1
                     mn = yzq(m, n)
                     cout(gpt) = ${type2type("cin(l, mn)", "c3d", kind)}$
                  END DO
!$OMP END PARALLEL DO
               END IF

            END ASSOCIATE

            CALL timestop(handle)

         END SUBROUTINE fft_gather_p_${kind}$

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param c ...
!> \param scale ...
! **************************************************************************************************
         SUBROUTINE fft_scatter_p_${kind}$ (cin, cout, map, ghat, yzq, scale)
            ${type}$, CONTIGUOUS, INTENT(IN)                          :: cin
            COMPLEX(KIND=dp), DIMENSION(:, :), CONTIGUOUS, INTENT(INOUT)   :: cout
            TYPE(cp_1d_i_a_type), DIMENSION(3, 2), INTENT(IN)          :: map
            INTEGER, DIMENSION(:, :), INTENT(IN) :: ghat
            INTEGER, DIMENSION(:, :), INTENT(IN) :: yzq
            REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

            CHARACTER(len=*), PARAMETER                        :: routineN = 'fft_scatter_p'

            INTEGER                                            :: gpt, handle, l, m, mn, n

            CALL timeset(routineN, handle)

               IF (.NOT. PRESENT(scale)) cout = z_zero

            ASSOCIATE (mapl => map(1, 1)%array, mapm => map(2, 1)%array, mapn => map(3, 1)%array)

               IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(l, m, mn, n) SHARED(cin, cout, scale, ghat, yzq)
                  DO gpt = 1, SIZE(cin)
                     l = mapl(ghat(1, gpt)) + 1
                     m = mapm(ghat(2, gpt)) + 1
                     n = mapn(ghat(3, gpt)) + 1
                     mn = yzq(m, n)
                     cout(l, mn) = ${type2type("scale*cin(gpt)", kind, "c3d")}$
                  END DO
!$OMP END PARALLEL DO
               ELSE
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(l, m, mn, n) SHARED(cin, cout, ghat, yzq)
                  DO gpt = 1, SIZE(cin)
                     l = mapl(ghat(1, gpt)) + 1
                     m = mapm(ghat(2, gpt)) + 1
                     n = mapn(ghat(3, gpt)) + 1
                     mn = yzq(m, n)
                     cout(l, mn) = ${type2type("cin(gpt)", kind, "c3d")}$
                  END DO
!$OMP END PARALLEL DO
               END IF

            END ASSOCIATE

            IF (ALLOCATED(map(1, 2)%array) .AND. ALLOCATED(map(2, 2)%array) .AND. ALLOCATED(map(3, 2)%array)) THEN

               ASSOCIATE (mapm => map(2, 2)%array, mapn => map(3, 2)%array, mapl => map(1, 2)%array)

                  IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(l, m, mn, n) SHARED(cin, cout, scale, ghat, yzq)
                     DO gpt = 1, SIZE(cin)
                        l = mapl(ghat(1, gpt)) + 1
                        m = mapm(ghat(2, gpt)) + 1
                        n = mapn(ghat(3, gpt)) + 1
                        mn = yzq(m, n)
                        cout(l, mn) = scale*#{if kind[0]=="c"}#CONJG#{endif}#(${type2type("cin(gpt)", kind, "c3d")}$)
                     END DO
!$OMP END PARALLEL DO
                  ELSE
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(l, m, mn, n) SHARED(cin, cout, ghat, yzq)
                     DO gpt = 1, SIZE(cin)
                        l = mapl(ghat(1, gpt)) + 1
                        m = mapm(ghat(2, gpt)) + 1
                        n = mapn(ghat(3, gpt)) + 1
                        mn = yzq(m, n)
                        cout(l, mn) = #{if kind[0]=="c"}#CONJG#{endif}#(${type2type("cin(gpt)", kind, "c3d")}$)
                     END DO
!$OMP END PARALLEL DO
                  END IF
               END ASSOCIATE
            END IF

            CALL timestop(handle)

         END SUBROUTINE fft_scatter_p_${kind}$
      #:endif
   #:endfor

   #:for kind, type, kind2, type2 in pw_list2
      #:if kind[1]=='1' and kind2[1]=='3'

! **************************************************************************************************
!> \brief Gathers the pw vector from a 3d data field
!> \param pw ...
!> \param c ...
!> \param scale ...
!> \par History
!>      none
!> \author JGH
! **************************************************************************************************
         SUBROUTINE fft_gather_s_${kind}$_${kind2}$ (cout, cin, map, ghat, scale)

            ${type}$, CONTIGUOUS, INTENT(INOUT)                       :: cout
            ${type2}$, CONTIGUOUS, INTENT(IN)   :: cin
            TYPE(cp_1d_i_a_type), DIMENSION(3, 2), INTENT(IN)          :: map
      INTEGER, DIMENSION(:, :), INTENT(IN) :: ghat
            REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

            CHARACTER(len=*), PARAMETER                        :: routineN = 'fft_gather_s'

            INTEGER                                            :: gpt, handle, l, m, n

            CALL timeset(routineN, handle)

            ASSOCIATE (mapl => map(1, 1)%array, mapm => map(2, 1)%array, mapn => map(3, 1)%array)

               IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(cin, cout, scale, ghat)
                  DO gpt = 1, SIZE(cout)
                     l = mapl(ghat(1, gpt)) + 1
                     m = mapm(ghat(2, gpt)) + 1
                     n = mapn(ghat(3, gpt)) + 1
                     cout(gpt) = scale*${type2type("cin(l, m, n)", kind2, kind)}$
                  END DO
!$OMP END PARALLEL DO
               ELSE
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(cin, cout, ghat)
                  DO gpt = 1, SIZE(cout)
                     l = mapl(ghat(1, gpt)) + 1
                     m = mapm(ghat(2, gpt)) + 1
                     n = mapn(ghat(3, gpt)) + 1
                     cout(gpt) = ${type2type("cin(l, m, n)", kind2, kind)}$
                  END DO
!$OMP END PARALLEL DO
               END IF

            END ASSOCIATE

            CALL timestop(handle)

         END SUBROUTINE fft_gather_s_${kind}$_${kind2}$

! **************************************************************************************************
!> \brief Scatters a pw vector to a 3d data field
!> \param pw ...
!> \param c ...
!> \param scale ...
!> \par History
!>      none
!> \author JGH
! **************************************************************************************************
         SUBROUTINE fft_scatter_s_${kind}$_${kind2}$ (cin, cout, map, ghat, scale)

            ${type}$, CONTIGUOUS, INTENT(IN)                 :: cin
            ${type2}$, CONTIGUOUS, INTENT(INOUT)               :: cout
            TYPE(cp_1d_i_a_type), DIMENSION(3, 2), INTENT(IN)          :: map
      INTEGER, DIMENSION(:, :), INTENT(IN) :: ghat
            REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

            CHARACTER(len=*), PARAMETER                        :: routineN = 'fft_scatter_s'

            INTEGER                                            :: gpt, handle, l, m, n

            CALL timeset(routineN, handle)

               ! should only zero the unused bits (but the zero is needed)
               IF (.NOT. PRESENT(scale)) cout = 0.0_dp

            ASSOCIATE (mapl => map(1, 1)%array, mapm => map(2, 1)%array, mapn => map(3, 1)%array)

               IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(cin, cout, scale, ghat)
                  DO gpt = 1, SIZE(cin)
                     l = mapl(ghat(1, gpt)) + 1
                     m = mapm(ghat(2, gpt)) + 1
                     n = mapn(ghat(3, gpt)) + 1
                     cout(l, m, n) = scale*${type2type("cin(gpt)", kind, kind2)}$
                  END DO
!$OMP END PARALLEL DO
               ELSE
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(cin, cout, ghat)
                  DO gpt = 1, SIZE(cin)
                     l = mapl(ghat(1, gpt)) + 1
                     m = mapm(ghat(2, gpt)) + 1
                     n = mapn(ghat(3, gpt)) + 1
                     cout(l, m, n) = ${type2type("cin(gpt)", kind, kind2)}$
                  END DO
!$OMP END PARALLEL DO
               END IF

            END ASSOCIATE

            IF (ALLOCATED(map(1, 2)%array) .AND. ALLOCATED(map(2, 2)%array) .AND. ALLOCATED(map(3, 2)%array)) THEN
               ASSOCIATE (mapl => map(1, 2)%array, mapm => map(2, 2)%array, mapn => map(3, 2)%array)

                  IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(cin, cout, scale, ghat)
                     DO gpt = 1, SIZE(cin)
                        l = mapl(ghat(1, gpt)) + 1
                        m = mapm(ghat(2, gpt)) + 1
                        n = mapn(ghat(3, gpt)) + 1
                        cout(l, m, n) = scale*#{if kind[0]=="c" and kind2[0]=="c"}#CONJG#{endif}#(${type2type("cin(gpt)", kind, kind2)}$)
                     END DO
!$OMP END PARALLEL DO
                  ELSE
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(cin, cout, ghat)
                     DO gpt = 1, SIZE(cin)
                        l = mapl(ghat(1, gpt)) + 1
                        m = mapm(ghat(2, gpt)) + 1
                        n = mapn(ghat(3, gpt)) + 1
                        cout(l, m, n) = #{if kind[0]=="c" and kind2[0]=="c"}#CONJG#{endif}#(${type2type("cin(gpt)", kind, kind2)}$)
                     END DO
!$OMP END PARALLEL DO
                  END IF

               END ASSOCIATE

            END IF

            CALL timestop(handle)

         END SUBROUTINE fft_scatter_s_${kind}$_${kind2}$
      #:endif
   #:endfor

END MODULE

