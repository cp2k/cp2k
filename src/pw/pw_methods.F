!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \note
!>      If parallel mode is distributed certain combination of
!>      "in_use" and "in_space" can not be used.
!>      For performance reasons it would be better to have the loops
!>      over g-vectors in the gather/scatter routines in new subprograms
!>      with the actual arrays (also the addressing) in the parameter list
!> \par History
!>      JGH (29-Dec-2000) : Changes for parallel use
!>      JGH (13-Mar-2001) : added timing calls
!>      JGH (26-Feb-2003) : OpenMP enabled
!>      JGH (17-Nov-2007) : Removed mass arrays
!>      JGH (01-Dec-2007) : Removed and renamed routines
!>      JGH (04-Jul-2019) : added pw_multiply routine
!>      03.2008 [tlaino] : Splitting pw_types into pw_types and pw_methods
!> \author apsi
! **************************************************************************************************
MODULE pw_methods
   #:include "pw_types.fypp"
   USE cp_log_handling, ONLY: cp_logger_get_default_io_unit, &
                              cp_to_string
   USE fft_tools, ONLY: BWFFT, &
                        FWFFT, &
                        fft3d
   USE kahan_sum, ONLY: accurate_dot_product, &
                        accurate_sum
   USE kinds, ONLY: dp
   USE machine, ONLY: m_memory
   USE mathconstants, ONLY: z_zero
   USE pw_copy_all, ONLY: pw_copy_match
   USE pw_fpga, ONLY: pw_fpga_c1dr3d_3d_dp, &
                      pw_fpga_c1dr3d_3d_sp, &
                      pw_fpga_init_bitstream, &
                      pw_fpga_r3dc1d_3d_dp, &
                      pw_fpga_r3dc1d_3d_sp
   USE pw_gpu, ONLY: pw_gpu_c1dr3d_3d, &
                     pw_gpu_c1dr3d_3d_ps, &
                     pw_gpu_r3dc1d_3d, &
                     pw_gpu_r3dc1d_3d_ps
   USE pw_grid_types, ONLY: HALFSPACE, &
                            PW_MODE_DISTRIBUTED, &
                            PW_MODE_LOCAL, &
                            pw_grid_type
   #:for space in pw_spaces
      #:for kind in pw_kinds
         USE pw_types, ONLY: pw_${kind}$_${space}$_type
      #:endfor
   #:endfor
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pw_zero, pw_structure_factor, pw_smoothing
   PUBLIC :: pw_copy, pw_axpy, pw_transfer, pw_scale
   PUBLIC :: pw_gauss_damp, pw_compl_gauss_damp, pw_derive, pw_laplace, pw_dr2, pw_write, pw_multiply
   PUBLIC :: pw_log_deriv_gauss, pw_log_deriv_compl_gauss, pw_log_deriv_mix_cl, pw_log_deriv_trunc
   PUBLIC :: pw_gauss_damp_mix, pw_multiply_with
   PUBLIC :: pw_integral_ab, pw_integral_a2b
   PUBLIC :: pw_dr2_gg, pw_integrate_function
   PUBLIC :: pw_set, pw_truncated
   PUBLIC :: pw_scatter, pw_gather
   PUBLIC :: pw_copy_to_array, pw_copy_from_array

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_methods'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.
   INTEGER, PARAMETER, PUBLIC  ::  do_accurate_sum = 0, &
                                  do_standard_sum = 1

   INTERFACE pw_zero
      #:for space in pw_spaces
         #:for kind in pw_kinds
            MODULE PROCEDURE pw_zero_${kind}$_${space}$
         #:endfor
      #:endfor
   END INTERFACE

   INTERFACE pw_scale
      #:for space in pw_spaces
         #:for kind in pw_kinds
            MODULE PROCEDURE pw_scale_${kind}$_${space}$
         #:endfor
      #:endfor
   END INTERFACE

   INTERFACE pw_write
      #:for space in pw_spaces
         #:for kind in pw_kinds
            MODULE PROCEDURE pw_write_${kind}$_${space}$
         #:endfor
      #:endfor
   END INTERFACE

   INTERFACE pw_integrate_function
      #:for space in pw_spaces
         #:for kind in pw_kinds
            MODULE PROCEDURE pw_integrate_function_${kind}$_${space}$
         #:endfor
      #:endfor
   END INTERFACE

   INTERFACE pw_set
      #:for space in pw_spaces
         #:for kind in pw_kinds
            MODULE PROCEDURE pw_set_value_${kind}$_${space}$
            MODULE PROCEDURE pw_zero_${kind}$_${space}$
         #:endfor
      #:endfor
   END INTERFACE

   INTERFACE pw_copy
      #:for space in pw_spaces
         #:for kind, kind2 in pw_kinds2_sameD
            MODULE PROCEDURE pw_copy_${kind}$_${kind2}$_${space}$
         #:endfor
      #:endfor
   END INTERFACE

   INTERFACE pw_axpy
      #:for space in pw_spaces
         #:for kind, kind2 in pw_kinds2_sameD
            MODULE PROCEDURE pw_axpy_${kind}$_${kind2}$_${space}$
         #:endfor
      #:endfor
   END INTERFACE

   INTERFACE pw_multiply
      #:for space in pw_spaces
         #:for kind, kind2 in pw_kinds2_sameD
            MODULE PROCEDURE pw_multiply_${kind}$_${kind2}$_${space}$
         #:endfor
      #:endfor
   END INTERFACE

   INTERFACE pw_multiply_with
      #:for space in pw_spaces
         #:for kind, kind2 in pw_kinds2_sameD
            MODULE PROCEDURE pw_multiply_with_${kind}$_${kind2}$_${space}$
         #:endfor
      #:endfor
   END INTERFACE

   INTERFACE pw_integral_ab
      #:for space in pw_spaces
         #:for kind, kind2 in pw_kinds2_sameD
            MODULE PROCEDURE pw_integral_ab_${kind}$_${kind2}$_${space}$
         #:endfor
      #:endfor
   END INTERFACE

   INTERFACE pw_integral_a2b
      #:for kind, kind2 in pw_kinds2_sameD
         #:if kind[1]=="1"
            MODULE PROCEDURE pw_integral_a2b_${kind}$_${kind2}$
         #:endif
      #:endfor
   END INTERFACE

   INTERFACE pw_gather
      #:for kind in pw_kinds
         #:if kind[1]=="1"
            MODULE PROCEDURE pw_gather_p_${kind}$
         #:endif
      #:endfor
      #:for kind, kind2 in pw_kinds2
         #:if kind[1]=="1" and kind2[1]=="3"
            MODULE PROCEDURE pw_gather_s_${kind}$_${kind2}$
         #:endif
      #:endfor
   END INTERFACE

   INTERFACE pw_scatter
      #:for kind in pw_kinds
         #:if kind[1]=="1"
            MODULE PROCEDURE pw_scatter_p_${kind}$
         #:endif
      #:endfor
      #:for kind, kind2 in pw_kinds2
         #:if kind[1]=="1" and kind2[1]=="3"
            MODULE PROCEDURE pw_scatter_s_${kind}$_${kind2}$
         #:endif
      #:endfor
   END INTERFACE

   INTERFACE pw_copy_to_array
      #:for space in pw_spaces
         #:for kind, kind2 in pw_kinds2_sameD
            MODULE PROCEDURE pw_copy_to_array_${kind}$_${kind2}$_${space}$
         #:endfor
      #:endfor
   END INTERFACE

   INTERFACE pw_copy_from_array
      #:for space in pw_spaces
         #:for kind, kind2 in pw_kinds2_sameD
            MODULE PROCEDURE pw_copy_from_array_${kind}$_${kind2}$_${space}$
         #:endfor
      #:endfor
   END INTERFACE

   INTERFACE pw_transfer
      #:for kind, kind2 in pw_kinds2
         #:if kind[1]=="1" and kind2[1]=="3"
            MODULE PROCEDURE pw_gather_s_${kind}$_${kind2}$_2
            MODULE PROCEDURE pw_scatter_s_${kind}$_${kind2}$_2
         #:endif
         #:for space in pw_spaces
            #:if kind[1]==kind2[1]
               MODULE PROCEDURE pw_copy_${kind}$_${kind2}$_${space}$
            #:endif
         #:endfor
         #:if kind2[0]=="c" and kind[1]=="3"
            MODULE PROCEDURE fft_wrap_pw1pw2_${kind}$_${kind2}$_rs_gs
         #:endif
         #:if kind[0]=="c" and kind2[1]=="3"
            MODULE PROCEDURE fft_wrap_pw1pw2_${kind}$_${kind2}$_gs_rs
         #:endif
      #:endfor
   END INTERFACE

CONTAINS
   #:for kind, type in pw_list
      #:for space in pw_spaces
! **************************************************************************************************
!> \brief Set values of a pw type to zero
!> \param pw ...
!> \par History
!>      none
!> \author apsi
! **************************************************************************************************
         SUBROUTINE pw_zero_${kind}$_${space}$ (pw)

            TYPE(pw_${kind}$_${space}$_type), INTENT(INOUT)                       :: pw

            CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_zero'

            INTEGER                                            :: handle

            CALL timeset(routineN, handle)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw)
            pw%array = ${type2type("0.0_dp", "r3d", kind)}$
!$OMP END PARALLEL WORKSHARE

            CALL timestop(handle)

         END SUBROUTINE pw_zero_${kind}$_${space}$

! **************************************************************************************************
!> \brief multiplies pw coeffs with a number
!> \param pw ...
!> \param a ...
!> \par History
!>      11.2004 created [Joost VandeVondele]
! **************************************************************************************************
         SUBROUTINE pw_scale_${kind}$_${space}$ (pw, a)

            TYPE(pw_${kind}$_${space}$_type), INTENT(INOUT)                       :: pw
            REAL(KIND=dp), INTENT(IN)                          :: a

            CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_scale'

            INTEGER                                            :: handle

            CALL timeset(routineN, handle)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(A, pw)
            pw%array = a*pw%array
!$OMP END PARALLEL WORKSHARE

            CALL timestop(handle)

         END SUBROUTINE pw_scale_${kind}$_${space}$

! **************************************************************************************************
!> \brief writes a small description of the actual grid
!>      (change to output the data as cube file, maybe with an
!>      optional long_description arg?)
!> \param pw the pw data to output
!> \param unit_nr the unit to output to
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
         SUBROUTINE pw_write_${kind}$_${space}$ (pw, unit_nr)

            TYPE(pw_${kind}$_${space}$_type), INTENT(in)                          :: pw
            INTEGER, INTENT(in)                                :: unit_nr

            INTEGER                                            :: iostatus

            WRITE (unit=unit_nr, fmt="('<pw>:{ ')", iostat=iostatus)

            WRITE (unit=unit_nr, fmt="(a)", iostat=iostatus) " in_use=${kind}$"
            IF (ASSOCIATED(pw%array)) THEN
               #:if kind[1]=='1'
                  WRITE (unit=unit_nr, fmt="(' array=<${kind[0]}$(',i8,':',i8,')>')") &
                     LBOUND(pw%array, 1), UBOUND(pw%array, 1)
               #:else
                  WRITE (unit=unit_nr, fmt="(' array=<${kind[0]}$(',i8,':',i8,',',i8,':',i8,',',i8,':',i8,')>')") &
                     LBOUND(pw%array, 1), UBOUND(pw%array, 1), LBOUND(pw%array, 2), UBOUND(pw%array, 2), &
                     LBOUND(pw%array, 3), UBOUND(pw%array, 3)
               #:endif
            ELSE
               WRITE (unit=unit_nr, fmt="(' array=*null*')")
            END IF

         END SUBROUTINE pw_write_${kind}$_${space}$

! **************************************************************************************************
!> \brief ...
!> \param fun ...
!> \param isign ...
!> \param oprt ...
!> \return ...
! **************************************************************************************************
         FUNCTION pw_integrate_function_${kind}$_${space}$ (fun, isign, oprt) RESULT(total_fun)

            TYPE(pw_${kind}$_${space}$_type), INTENT(IN)                          :: fun
            INTEGER, INTENT(IN), OPTIONAL                      :: isign
            CHARACTER(len=*), INTENT(IN), OPTIONAL             :: oprt
            REAL(KIND=dp)                                      :: total_fun

            INTEGER                                            :: iop

            iop = 0

            IF (PRESENT(oprt)) THEN
               SELECT CASE (oprt)
               CASE ("ABS", "abs")
                  iop = 1
               CASE DEFAULT
                  CPABORT("Unknown operator")
               END SELECT
            END IF

            total_fun = 0.0_dp

            #:if space == "rs"
               ! do reduction using maximum accuracy
               IF (iop == 1) THEN
                  total_fun = fun%pw_grid%dvol*accurate_sum(ABS(fun%array))
               ELSE
                  total_fun = fun%pw_grid%dvol*accurate_sum(${type2type("fun%array", kind, "r1d")}$)
               END IF
            #:elif space == "gs"
               IF (iop == 1) &
                  CPABORT("Operator ABS not implemented")
               #:if kind[1:]=="1d"
                  IF (fun%pw_grid%have_g0) total_fun = fun%pw_grid%vol*${type2type("fun%array(1)", kind, "r1d")}$
               #:else
                  CPABORT("Reciprocal space integration for 3D grids not implemented!")
               #:endif
            #:else
               CPABORT("No space defined")
            #:endif

            IF (fun%pw_grid%para%mode /= PW_MODE_LOCAL) THEN
               CALL fun%pw_grid%para%group%sum(total_fun)
            END IF

            IF (PRESENT(isign)) THEN
               total_fun = total_fun*SIGN(1._dp, REAL(isign, dp))
            END IF

         END FUNCTION pw_integrate_function_${kind}$_${space}$

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param value ...
! **************************************************************************************************
         SUBROUTINE pw_set_value_${kind}$_${space}$ (pw, value)
            TYPE(pw_${kind}$_${space}$_type), INTENT(IN)                          :: pw
            REAL(KIND=dp), INTENT(IN)                          :: value

            CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_set_value'

            INTEGER                                            :: handle

            CALL timeset(routineN, handle)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw,value)
            pw%array = ${type2type("value", "r3d", kind)}$
!$OMP END PARALLEL WORKSHARE

            CALL timestop(handle)

         END SUBROUTINE pw_set_value_${kind}$_${space}$
      #:endfor

      #:if kind[1]=="1"
! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param c ...
!> \param scale ...
! **************************************************************************************************
         SUBROUTINE pw_gather_p_${kind}$ (pw, c, scale)

            TYPE(pw_${kind}$_gs_type), INTENT(INOUT)                       :: pw
            COMPLEX(KIND=dp), DIMENSION(:, :), CONTIGUOUS, INTENT(IN)      :: c
            REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

            CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gather_p'

            INTEGER                                            :: gpt, handle, l, m, mn, n

            CALL timeset(routineN, handle)

            IF (pw%pw_grid%para%mode /= PW_MODE_DISTRIBUTED) THEN
               CPABORT("This grid type is not distributed")
            END IF

            ASSOCIATE (mapl => pw%pw_grid%mapl%pos, mapm => pw%pw_grid%mapm%pos, mapn => pw%pw_grid%mapn%pos, &
                       ngpts => SIZE(pw%pw_grid%gsq), ghat => pw%pw_grid%g_hat, yzq => pw%pw_grid%para%yzq)

               IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(l, m, mn, n), &
!$OMP             SHARED(c, pw, scale)
                  DO gpt = 1, ngpts
                     l = mapl(ghat(1, gpt)) + 1
                     m = mapm(ghat(2, gpt)) + 1
                     n = mapn(ghat(3, gpt)) + 1
                     mn = yzq(m, n)
                     pw%array(gpt) = scale*${type2type("c(l, mn)", "c3d", kind)}$
                  END DO
!$OMP END PARALLEL DO
               ELSE
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(l, m, mn, n), &
!$OMP             SHARED(c, pw)
                  DO gpt = 1, ngpts
                     l = mapl(ghat(1, gpt)) + 1
                     m = mapm(ghat(2, gpt)) + 1
                     n = mapn(ghat(3, gpt)) + 1
                     mn = yzq(m, n)
                     pw%array(gpt) = ${type2type("c(l, mn)", "c3d", kind)}$
                  END DO
!$OMP END PARALLEL DO
               END IF

            END ASSOCIATE

            CALL timestop(handle)

         END SUBROUTINE pw_gather_p_${kind}$

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param c ...
!> \param scale ...
! **************************************************************************************************
         SUBROUTINE pw_scatter_p_${kind}$ (pw, c, scale)
            TYPE(pw_${kind}$_gs_type), INTENT(IN)                          :: pw
            COMPLEX(KIND=dp), DIMENSION(:, :), CONTIGUOUS, INTENT(INOUT)   :: c
            REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

            CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_scatter_p'

            INTEGER                                            :: gpt, handle, l, m, mn, n

            CALL timeset(routineN, handle)

            IF (pw%pw_grid%para%mode /= PW_MODE_DISTRIBUTED) THEN
               CPABORT("This grid type is not distributed")
            END IF

            ASSOCIATE (mapl => pw%pw_grid%mapl%pos, mapm => pw%pw_grid%mapm%pos, mapn => pw%pw_grid%mapn%pos, &
                       ghat => pw%pw_grid%g_hat, yzq => pw%pw_grid%para%yzq, ngpts => SIZE(pw%pw_grid%gsq))

               IF (.NOT. PRESENT(scale)) c = z_zero

               IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(l, m, mn, n), &
!$OMP             SHARED(c, pw, scale)
                  DO gpt = 1, ngpts
                     l = mapl(ghat(1, gpt)) + 1
                     m = mapm(ghat(2, gpt)) + 1
                     n = mapn(ghat(3, gpt)) + 1
                     mn = yzq(m, n)
                     c(l, mn) = ${type2type("scale*pw%array(gpt)", kind, "c3d")}$
                  END DO
!$OMP END PARALLEL DO
               ELSE
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(l, m, mn, n), &
!$OMP             SHARED(c, pw)
                  DO gpt = 1, ngpts
                     l = mapl(ghat(1, gpt)) + 1
                     m = mapm(ghat(2, gpt)) + 1
                     n = mapn(ghat(3, gpt)) + 1
                     mn = yzq(m, n)
                     c(l, mn) = ${type2type("pw%array(gpt)", kind, "c3d")}$
                  END DO
!$OMP END PARALLEL DO
               END IF

            END ASSOCIATE

            IF (pw%pw_grid%grid_span == HALFSPACE) THEN

               ASSOCIATE (mapm => pw%pw_grid%mapm%neg, mapn => pw%pw_grid%mapn%neg, mapl => pw%pw_grid%mapl%neg, &
                          ghat => pw%pw_grid%g_hat, ngpts => SIZE(pw%pw_grid%gsq), yzq => pw%pw_grid%para%yzq)

                  IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(l, m, mn, n), &
!$OMP             SHARED(c, pw, scale)
                     DO gpt = 1, ngpts
                        l = mapl(ghat(1, gpt)) + 1
                        m = mapm(ghat(2, gpt)) + 1
                        n = mapn(ghat(3, gpt)) + 1
                        mn = yzq(m, n)
                        c(l, mn) = scale*#{if kind[0]=="c"}#CONJG#{endif}#(${type2type("pw%array(gpt)", kind, "c3d")}$)
                     END DO
!$OMP END PARALLEL DO
                  ELSE
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(l, m, mn, n) &
!$OMP             SHARED(c, pw)
                     DO gpt = 1, ngpts
                        l = mapl(ghat(1, gpt)) + 1
                        m = mapm(ghat(2, gpt)) + 1
                        n = mapn(ghat(3, gpt)) + 1
                        mn = yzq(m, n)
                        c(l, mn) = #{if kind[0]=="c"}#CONJG#{endif}#(${type2type("pw%array(gpt)", kind, "c3d")}$)
                     END DO
!$OMP END PARALLEL DO
                  END IF
               END ASSOCIATE
            END IF

            CALL timestop(handle)

         END SUBROUTINE pw_scatter_p_${kind}$
      #:endif
   #:endfor

   #:for kind, type, kind2, type2 in pw_list2_sameD
      #:for space in pw_spaces
! **************************************************************************************************
!> \brief copy a pw type variable
!> \param pw1 ...
!> \param pw2 ...
!> \par History
!>      JGH (7-Mar-2001) : check for pw_grid %id_nr, allow copy if
!>        in_use == COMPLEXDATA1D and in_space == RECIPROCALSPACE
!>      JGH (21-Feb-2003) : Code for generalized reference grids
!> \author apsi
!> \note
!>      Currently only copying of respective types allowed,
!>      in order to avoid errors
! **************************************************************************************************
         SUBROUTINE pw_copy_${kind}$_${kind2}$_${space}$ (pw1, pw2)

            TYPE(pw_${kind}$_${space}$_type), INTENT(IN)                          :: pw1
            TYPE(pw_${kind2}$_${space}$_type), INTENT(INOUT)                       :: pw2

            CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_copy'

            INTEGER                                            :: handle
            #:if kind[1:]=='1d'
               INTEGER :: i, j, ng, ng1, ng2, ns
            #:endif

            CALL timeset(routineN, handle)

            IF (pw1%pw_grid%spherical .neqv. pw2%pw_grid%spherical) &
               CPABORT("Both grids must be either spherical or non-spherical!")
            #:if space != "gs"
               IF (pw1%pw_grid%spherical) &
                  CPABORT("Spherical grids only exist in reciprocal space!")
            #:endif

            #:if kind[1]=='1'
               IF (.NOT. ASSOCIATED(pw1%pw_grid, pw2%pw_grid)) THEN
                  IF (pw1%pw_grid%spherical .AND. pw2%pw_grid%spherical) THEN
                     IF (pw_compatible(pw1%pw_grid, pw2%pw_grid)) THEN
                        ng1 = SIZE(pw1%array)
                        ng2 = SIZE(pw2%array)
                        ng = MIN(ng1, ng2)
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(ng, pw1, pw2)
                        pw2%array(1:ng) = ${type2type("pw1%array(1:ng)", kind, kind2)}$
!$OMP END PARALLEL WORKSHARE
                        IF (ng2 > ng) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(ng, ng2, pw2)
                           pw2%array(ng + 1:ng2) = ${type2type("0.0_dp", "r3d", kind2)}$
!$OMP END PARALLEL WORKSHARE
                        END IF
                     ELSE
                        CPABORT("Copies between spherical grids require compatible grids!")
                     END IF
                  ELSE
                     ng1 = SIZE(pw1%array)
                     ng2 = SIZE(pw2%array)
                     ns = 2*MAX(ng1, ng2)

                     IF ((pw1%pw_grid%id_nr == pw2%pw_grid%reference)) THEN
                        IF (ng1 >= ng2) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng2, pw1, pw2)
                           DO i = 1, ng2
                              j = pw2%pw_grid%gidx(i)
                              pw2%array(i) = ${type2type("pw1%array(j)", kind, kind2)}$
                           END DO
!$OMP END PARALLEL DO
                        ELSE
                           CALL pw_zero(pw2)
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng1, pw1, pw2)
                           DO i = 1, ng1
                              j = pw2%pw_grid%gidx(i)
                              pw2%array(j) = ${type2type("pw1%array(i)", kind, kind2)}$
                           END DO
!$OMP END PARALLEL DO
                        END IF
                     ELSE IF ((pw2%pw_grid%id_nr == pw1%pw_grid%reference)) THEN
                        IF (ng1 >= ng2) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng2, pw1, pw2)
                           DO i = 1, ng2
                              j = pw1%pw_grid%gidx(i)
                              pw2%array(i) = ${type2type("pw1%array(j)", kind, kind2)}$
                           END DO
!$OMP END PARALLEL DO
                        ELSE
                           CALL pw_zero(pw2)
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng1, pw1, pw2)
                           DO i = 1, ng1
                              j = pw1%pw_grid%gidx(i)
                              pw2%array(j) = ${type2type("pw1%array(i)", kind, kind2)}$
                           END DO
!$OMP END PARALLEL DO
                        END IF
                     ELSE
                        #:if kind=='c1d' and kind2=='c1d' and space=="gs"
                           CALL pw_copy_match(pw1, pw2)
                        #:else
                           CPABORT("Copy not implemented!")
                        #:endif
                     END IF

                  END IF

               ELSE
!$OMP PARALLEL WORKSHARE PRIVATE(i) DEFAULT(NONE) SHARED(pw1, pw2)
                  pw2%array = ${type2type("pw1%array", kind, kind2)}$
!$OMP END PARALLEL WORKSHARE
               END IF
            #:else
               IF (ANY(SHAPE(pw2%array) /= SHAPE(pw1%array))) &
                  CPABORT("3D grids must be compatible!")
               IF (pw1%pw_grid%spherical) &
                  CPABORT("3D grids must not be spherical!")
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,pw2)
               pw2%array(:, :, :) = ${type2type("pw1%array(:, :, :)", kind, kind2)}$
!$OMP END PARALLEL WORKSHARE
            #:endif

            CALL timestop(handle)

         END SUBROUTINE pw_copy_${kind}$_${kind2}$_${space}$

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param array ...
! **************************************************************************************************
         SUBROUTINE pw_copy_to_array_${kind}$_${kind2}$_${space}$ (pw, array)
            TYPE(pw_${kind}$_${space}$_type), INTENT(IN)                          :: pw
            ${type2}$, INTENT(INOUT)   :: array

            CHARACTER(len=*), PARAMETER :: routineN = 'pw_copy_to_array'

            INTEGER                                            :: handle

            CALL timeset(routineN, handle)

            #:if kind[1]=="1"
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw, array)
               array(:) = ${type2type("pw%array(:)", kind, kind2)}$
!$OMP END PARALLEL WORKSHARE
            #:else
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw, array)
               array(:, :, :) = ${type2type("pw%array(:, :, :)", kind, kind2)}$
!$OMP END PARALLEL WORKSHARE
            #:endif

            CALL timestop(handle)
         END SUBROUTINE pw_copy_to_array_${kind}$_${kind2}$_${space}$

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param array ...
! **************************************************************************************************
         SUBROUTINE pw_copy_from_array_${kind}$_${kind2}$_${space}$ (pw, array)
            TYPE(pw_${kind}$_${space}$_type), INTENT(IN)                          :: pw
            ${type2}$, INTENT(IN)      :: array

            CHARACTER(len=*), PARAMETER :: routineN = 'pw_copy_from_array'

            INTEGER                                            :: handle

            CALL timeset(routineN, handle)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw, array)
            pw%array = ${type2type("array", kind2, kind)}$
!$OMP END PARALLEL WORKSHARE

            CALL timestop(handle)
         END SUBROUTINE pw_copy_from_array_${kind}$_${kind2}$_${space}$

! **************************************************************************************************
!> \brief pw2 = alpha*pw1 + beta*pw2
!>      alpha defaults to 1, beta defaults to 1
!> \param pw1 ...
!> \param pw2 ...
!> \param alpha ...
!> \param beta ...
!> \param allow_noncompatible_grids ...
!> \par History
!>      JGH (21-Feb-2003) : added reference grid functionality
!>      JGH (01-Dec-2007) : rename and remove complex alpha
!> \author apsi
!> \note
!>      Currently only summing up of respective types allowed,
!>      in order to avoid errors
! **************************************************************************************************
         SUBROUTINE pw_axpy_${kind}$_${kind2}$_${space}$ (pw1, pw2, alpha, beta, allow_noncompatible_grids)

            TYPE(pw_${kind}$_${space}$_type), INTENT(IN)                          :: pw1
            TYPE(pw_${kind2}$_${space}$_type), INTENT(INOUT)                       :: pw2
            REAL(KIND=dp), INTENT(in), OPTIONAL                :: alpha, beta
            LOGICAL, INTENT(IN), OPTIONAL                      :: allow_noncompatible_grids

            CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_axpy'

            INTEGER                                            :: handle
            LOGICAL                                            :: my_allow_noncompatible_grids
            REAL(KIND=dp)                                      :: my_alpha, my_beta
            #:if kind[1]=='1'
               INTEGER                                            :: i, j, ng, ng1, ng2
            #:endif

            CALL timeset(routineN, handle)

            my_alpha = 1.0_dp
            IF (PRESENT(alpha)) my_alpha = alpha

            my_beta = 1.0_dp
            IF (PRESENT(beta)) my_beta = beta

            my_allow_noncompatible_grids = .FALSE.
            IF (PRESENT(allow_noncompatible_grids)) my_allow_noncompatible_grids = allow_noncompatible_grids

            IF (my_beta /= 1.0_dp) THEN
               IF (my_beta == 0.0_dp) THEN
                  CALL pw_zero(pw2)
               ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw2,my_beta)
                  pw2%array = pw2%array*my_beta
!$OMP END PARALLEL WORKSHARE
               END IF
            END IF

            #:if kind[1]=='1'
               IF (ASSOCIATED(pw1%pw_grid, pw2%pw_grid)) THEN

                  IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL WORKSHARE PRIVATE(i) DEFAULT(NONE) SHARED(pw1, pw2)
                     pw2%array = pw2%array + ${type2type("pw1%array", kind, kind2)}$
!$OMP END PARALLEL WORKSHARE
                  ELSE IF (my_alpha /= 0.0_dp) THEN
!$OMP PARALLEL WORKSHARE PRIVATE(i) DEFAULT(NONE) SHARED(pw2,pw1,my_alpha)
                     pw2%array = pw2%array + my_alpha*${type2type("pw1%array", kind, kind2)}$
!$OMP END PARALLEL WORKSHARE
                  END IF

               ELSE IF (pw_compatible(pw1%pw_grid, pw2%pw_grid) .OR. my_allow_noncompatible_grids) THEN

                  ng1 = SIZE(pw1%array)
                  ng2 = SIZE(pw2%array)
                  ng = MIN(ng1, ng2)

                  IF (pw1%pw_grid%spherical) THEN
                     IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL DO PRIVATE(i) DEFAULT(NONE) SHARED(ng, pw1, pw2)
                        DO i = 1, ng
                           pw2%array(i) = pw2%array(i) + ${type2type("pw1%array(i)", kind, kind2)}$
                        END DO
!$OMP END PARALLEL DO
                     ELSE IF (my_alpha /= 0.0_dp) THEN
!$OMP PARALLEL DO PRIVATE(i) DEFAULT(NONE) SHARED(pw2,pw1,my_alpha,ng)
                        DO i = 1, ng
                           pw2%array(i) = pw2%array(i) + my_alpha*${type2type("pw1%array(i)", kind, kind2)}$
                        END DO
!$OMP END PARALLEL DO
                     END IF
                  ELSE IF ((pw1%pw_grid%id_nr == pw2%pw_grid%reference)) THEN
                     IF (ng1 >= ng2) THEN
                        IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng, pw1, pw2)
                           DO i = 1, ng
                              j = pw2%pw_grid%gidx(i)
                              pw2%array(i) = pw2%array(i) + ${type2type("pw1%array(j)", kind, kind2)}$
                           END DO
!$OMP END PARALLEL DO
                        ELSE IF (my_alpha /= 0.0_dp) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(my_alpha, ng, pw1, pw2)
                           DO i = 1, ng
                              j = pw2%pw_grid%gidx(i)
                              pw2%array(i) = pw2%array(i) + my_alpha*${type2type("pw1%array(j)", kind, kind2)}$
                           END DO
!$OMP END PARALLEL DO
                        END IF
                     ELSE
                        IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng, pw1, pw2)
                           DO i = 1, ng
                              j = pw2%pw_grid%gidx(i)
                              pw2%array(j) = pw2%array(j) + ${type2type("pw1%array(i)", kind, kind2)}$
                           END DO
!$OMP END PARALLEL DO
                        ELSE IF (my_alpha /= 0.0_dp) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(my_alpha, ng, pw1, pw2)
                           DO i = 1, ng
                              j = pw2%pw_grid%gidx(i)
                              pw2%array(j) = pw2%array(j) + my_alpha*${type2type("pw1%array(i)", kind, kind2)}$
                           END DO
!$OMP END PARALLEL DO
                        END IF
                     END IF
                  ELSE IF ((pw2%pw_grid%id_nr == pw1%pw_grid%reference)) THEN
                     IF (ng1 >= ng2) THEN
                        IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng, pw1, pw2)
                           DO i = 1, ng
                              j = pw1%pw_grid%gidx(i)
                              pw2%array(i) = pw2%array(i) + ${type2type("pw1%array(j)", kind, kind2)}$
                           END DO
!$OMP END PARALLEL DO
                        ELSE IF (my_alpha /= 0.0_dp) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(my_alpha, ng, pw1, pw2)
                           DO i = 1, ng
                              j = pw1%pw_grid%gidx(i)
                              pw2%array(i) = pw2%array(i) + my_alpha*${type2type("pw1%array(j)", kind, kind2)}$
                           END DO
!$OMP END PARALLEL DO
                        END IF
                     ELSE
                        IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(ng, pw1, pw2)
                           DO i = 1, ng
                              j = pw1%pw_grid%gidx(i)
                              pw2%array(j) = pw2%array(j) + ${type2type("pw1%array(i)", kind, kind2)}$
                           END DO
!$OMP END PARALLEL DO
                        ELSE IF (my_alpha /= 0.0_dp) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE) SHARED(my_alpha, ng, pw1, pw2)
                           DO i = 1, ng
                              j = pw1%pw_grid%gidx(i)
                              pw2%array(j) = pw2%array(j) + my_alpha*${type2type("pw1%array(i)", kind, kind2)}$
                           END DO
!$OMP END PARALLEL DO
                        END IF
                     END IF
                  ELSE
                     CPABORT("Grids not compatible")
                  END IF
                  #:else
                  IF (ASSOCIATED(pw1%pw_grid, pw2%pw_grid)) THEN
                     IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1, pw2)
                        pw2%array = pw2%array + ${type2type("pw1%array", kind, kind2)}$
!$OMP END PARALLEL WORKSHARE
                     ELSE IF (my_alpha /= 0.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1, pw2, my_alpha)
                        pw2%array = pw2%array + my_alpha*${type2type("pw1%array", kind, kind2)}$
!$OMP END PARALLEL WORKSHARE
                     END IF

                  ELSE IF (pw_compatible(pw1%pw_grid, pw2%pw_grid) .OR. my_allow_noncompatible_grids) THEN

                     IF (ANY(SHAPE(pw1%array) /= SHAPE(pw2%array))) &
                        CPABORT("Noncommensurate grids not implemented for 3D grids!")

                     IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,pw2)
                        pw2%array = pw2%array + ${type2type("pw1%array", kind, kind2)}$
!$OMP END PARALLEL WORKSHARE
                     ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,pw2,my_alpha)
                        pw2%array = pw2%array + my_alpha*${type2type("pw1%array", kind, kind2)}$
!$OMP END PARALLEL WORKSHARE
                     END IF
                     #:endif

                  ELSE

                     CPABORT("Grids not compatible")

                  END IF

                  CALL timestop(handle)

                  END SUBROUTINE pw_axpy_${kind}$_${kind2}$_${space}$

! **************************************************************************************************
!> \brief pw_out = pw_out + alpha * pw1 * pw2
!>      alpha defaults to 1
!> \param pw_out ...
!> \param pw1 ...
!> \param pw2 ...
!> \param alpha ...
!> \author JGH
! **************************************************************************************************
                  SUBROUTINE pw_multiply_${kind}$_${kind2}$_${space}$ (pw_out, pw1, pw2, alpha)

                     TYPE(pw_${kind}$_${space}$_type), INTENT(INOUT)                       :: pw_out
                     TYPE(pw_${kind}$_${space}$_type), INTENT(IN)                          :: pw1
                     TYPE(pw_${kind2}$_${space}$_type), INTENT(IN) :: pw2
                     REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha

                     CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_multiply'

                     INTEGER                                            :: handle
                     REAL(KIND=dp)                                      :: my_alpha

                     CALL timeset(routineN, handle)

                     my_alpha = 1.0_dp
                     IF (PRESENT(alpha)) my_alpha = alpha

                     IF (.NOT. ASSOCIATED(pw_out%pw_grid, pw1%pw_grid) .OR. .NOT. ASSOCIATED(pw_out%pw_grid, pw2%pw_grid)) &
                        CPABORT("pw_multiply not implemented for non-identical grids!")

                     IF (my_alpha == 1.0_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw_out, pw1, pw2)
                        pw_out%array = pw_out%array + pw1%array*${type2type("pw2%array", kind2, kind)}$
!$OMP END PARALLEL WORKSHARE
                     ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(my_alpha, pw_out, pw1, pw2)
                        pw_out%array = pw_out%array + my_alpha*pw1%array*${type2type("pw2%array", kind2, kind)}$
!$OMP END PARALLEL WORKSHARE
                     END IF

                     CALL timestop(handle)

                  END SUBROUTINE pw_multiply_${kind}$_${kind2}$_${space}$

! **************************************************************************************************
!> \brief ...
!> \param pw1 ...
!> \param pw2 ...
! **************************************************************************************************
                  SUBROUTINE pw_multiply_with_${kind}$_${kind2}$_${space}$ (pw1, pw2)
                     TYPE(pw_${kind}$_${space}$_type), INTENT(INOUT)                       :: pw1
                     TYPE(pw_${kind2}$_${space}$_type), INTENT(IN)                          :: pw2

                     CHARACTER(LEN=*), PARAMETER                        :: routineN = 'pw_multiply_with'

                     INTEGER                                            :: handle

                     CALL timeset(routineN, handle)

                     IF (.NOT. ASSOCIATED(pw1%pw_grid, pw2%pw_grid)) &
                        CPABORT("Incompatible grids!")

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw1,pw2)
                     pw1%array = pw1%array*${type2type("pw2%array", kind2, kind)}$
!$OMP END PARALLEL WORKSHARE

                     CALL timestop(handle)

                  END SUBROUTINE pw_multiply_with_${kind}$_${kind2}$_${space}$

! **************************************************************************************************
!> \brief Calculate integral over unit cell for functions in plane wave basis
!>      only returns the real part of it ......
!> \param pw1 ...
!> \param pw2 ...
!> \param sumtype ...
!> \param just_sum ...
!> \param local_only ...
!> \return ...
!> \par History
!>      JGH (14-Mar-2001) : Parallel sum and some tests, HALFSPACE case
!> \author apsi
! **************************************************************************************************
               FUNCTION pw_integral_ab_${kind}$_${kind2}$_${space}$ (pw1, pw2, sumtype, just_sum, local_only) RESULT(integral_value)

                     TYPE(pw_${kind}$_${space}$_type), INTENT(IN)                          :: pw1
                     TYPE(pw_${kind2}$_${space}$_type), INTENT(IN)                          :: pw2
                     INTEGER, INTENT(IN), OPTIONAL                      :: sumtype
                     LOGICAL, INTENT(IN), OPTIONAL                      :: just_sum, local_only
                     REAL(KIND=dp)                                      :: integral_value

                     CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_integral_ab'

                     INTEGER                                            :: handle, loc_sumtype
                     LOGICAL                                            :: my_just_sum, my_local_only

                     CALL timeset(routineN, handle)

                     loc_sumtype = do_accurate_sum
                     IF (PRESENT(sumtype)) loc_sumtype = sumtype

                     my_local_only = .FALSE.
                     IF (PRESENT(local_only)) my_local_only = local_only

                     IF (.NOT. ASSOCIATED(pw1%pw_grid, pw2%pw_grid)) THEN
                        CPABORT("Grids incompatible")
                     END IF

                     my_just_sum = .FALSE.
                     IF (PRESENT(just_sum)) my_just_sum = just_sum

                     ! do standard sum
                     IF (loc_sumtype == do_standard_sum) THEN

                        ! Do standard sum

                        #:if kind=="r1d" and kind2=="r1d"
                           integral_value = DOT_PRODUCT(pw1%array, pw2%array)
                        #:elif kind=="r3d" and kind2=="r3d"
                           integral_value = SUM(pw1%array*pw2%array)
                        #:elif kind[0]=="r"
                           integral_value = SUM(pw1%array*REAL(pw2%array, KIND=dp)) !? complex bit
                        #:elif kind2[0]=="r"
                           integral_value = SUM(REAL(pw1%array, KIND=dp)*pw2%array) !? complex bit
                        #:else
                           integral_value = SUM(REAL(CONJG(pw1%array) &
                                                     *pw2%array, KIND=dp)) !? complex bit
                        #:endif

                     ELSE

                        ! Do accurate sum
                        #:if kind[0]=="r" and kind2[0]=="r"
                           integral_value = accurate_dot_product(pw1%array, pw2%array)
                        #:elif kind[0]=="r"
                           integral_value = accurate_sum(pw1%array*REAL(pw2%array, KIND=dp)) !? complex bit
                        #:elif kind2[0]=="r"
                           integral_value = accurate_sum(REAL(pw1%array, KIND=dp)*pw2%array) !? complex bit
                        #:else
                           integral_value = accurate_sum(REAL(CONJG(pw1%array) &
                                                              *pw2%array, KIND=dp)) !? complex bit
                        #:endif

                     END IF

                     IF (.NOT. my_just_sum) THEN
                        #:if space == "rs"
                           integral_value = integral_value*pw1%pw_grid%dvol
                        #:elif space == "gs"
                           integral_value = integral_value*pw1%pw_grid%vol
                        #:else
                           #:stop "Unknown space: "+space
                        #:endif
                     END IF

                     #:if kind[1]=="1"
                        IF (pw1%pw_grid%grid_span == HALFSPACE) THEN
                           integral_value = 2.0_dp*integral_value
                           IF (pw1%pw_grid%have_g0) integral_value = integral_value - &
                           #:if kind[0]=="c"
                                                                     REAL(CONJG(pw1%array(1))*pw2%array(1), KIND=dp)
                              #:elif kind2[0]=="r"
                              pw1%array(1)*pw2%array(1)
                              #:else
                              pw1%array(1)*REAL(pw2%array(1), KIND=dp)
                              #:endif
                           END IF
                        #:endif

                        IF (.NOT. my_local_only .AND. pw1%pw_grid%para%mode == PW_MODE_DISTRIBUTED) &
                           CALL pw1%pw_grid%para%group%sum(integral_value)

                        CALL timestop(handle)

                     END FUNCTION pw_integral_ab_${kind}$_${kind2}$_${space}$
                     #:endfor

                     #:if kind[1]=="1"
! **************************************************************************************************
!> \brief ...
!> \param pw1 ...
!> \param pw2 ...
!> \return ...
! **************************************************************************************************
                        FUNCTION pw_integral_a2b_${kind}$_${kind2}$ (pw1, pw2) RESULT(integral_value)

                           TYPE(pw_${kind}$_gs_type), INTENT(IN)                          :: pw1
                           TYPE(pw_${kind2}$_gs_type), INTENT(IN) :: pw2
                           REAL(KIND=dp)                                      :: integral_value

                           CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_integral_a2b'

                           INTEGER                                            :: handle

                           CALL timeset(routineN, handle)

                           IF (.NOT. ASSOCIATED(pw1%pw_grid, pw2%pw_grid)) THEN
                              CPABORT("Grids incompatible")
                           END IF

                           #:if kind[0]=="c"
                              #:if kind2[0]=="c"
                                 integral_value = accurate_sum(REAL(CONJG(pw1%array)*pw2%array, KIND=dp)*pw1%pw_grid%gsq)
                              #:else
                                 integral_value = accurate_sum(REAL(CONJG(pw1%array), KIND=dp)*pw2%array*pw1%pw_grid%gsq)
                              #:endif
                           #:elif kind2[0]=="c"
                              integral_value = accurate_sum(pw1%array*REAL(pw2%array, KIND=dp)*pw1%pw_grid%gsq)
                           #:else
                              integral_value = accurate_sum(pw1%array*pw2%array*pw1%pw_grid%gsq)
                           #:endif
                           IF (pw1%pw_grid%grid_span == HALFSPACE) THEN
                              integral_value = 2.0_dp*integral_value
                           END IF

                           integral_value = integral_value*pw1%pw_grid%vol

                           IF (pw1%pw_grid%para%mode == PW_MODE_DISTRIBUTED) &
                              CALL pw1%pw_grid%para%group%sum(integral_value)
                           CALL timestop(handle)

                        END FUNCTION pw_integral_a2b_${kind}$_${kind2}$
                     #:endif
                  #:endfor

                  #:for kind, type, kind2, type2 in pw_list2
                     #:for space in pw_spaces
                        #:for space2 in pw_spaces

                           #:if space != space2 and ((space=="rs" and kind[1]=="3" and kind2[0]=="c") or (space=="gs" and kind2[1]=="3" and kind[0]=="c"))
! **************************************************************************************************
!> \brief Generic function for 3d FFT of a coefficient_type or pw_r3d_rs_type
!> \param pw1 ...
!> \param pw2 ...
!> \param debug ...
!> \par History
!>      JGH (30-12-2000): New setup of functions and adaptation to parallelism
!>      JGH (04-01-2001): Moved routine from pws to this module, only covers
!>                        pw_types, no more coefficient types
!> \author apsi
!> \note
!>       fft_wrap_pw1pw2
! **************************************************************************************************
                              SUBROUTINE fft_wrap_pw1pw2_${kind}$_${kind2}$_${space}$_${space2}$ (pw1, pw2, debug)

                                 TYPE(pw_${kind}$_${space}$_type), INTENT(IN)                  :: pw1
                                 TYPE(pw_${kind2}$_${space2}$_type), INTENT(INOUT)               :: pw2
                                 LOGICAL, INTENT(IN), OPTIONAL                      :: debug

                                 CHARACTER(len=*), PARAMETER                        :: routineN = 'fft_wrap_pw1pw2'

                                 COMPLEX(KIND=dp), DIMENSION(:, :), CONTIGUOUS, POINTER         :: grays
                                 COMPLEX(KIND=dp), DIMENSION(:, :, :), CONTIGUOUS, POINTER      :: c_in, c_out
                                 INTEGER                                            ::  handle, handle2, my_pos, nrays, &
                                                                                       out_unit
                                 #:if (space=="rs" and kind=="r3d" and kind2=="c1d") or (space=="gs" and kind=="c1d" and kind2=="r3d")
                                    INTEGER, DIMENSION(3)                              :: nloc
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                                    LOGICAL                                            :: use_pw_gpu
#endif
                                 #:endif
                                 INTEGER, DIMENSION(:), POINTER                     :: n
                                 LOGICAL                                            :: test
                                 REAL(KIND=dp)                                      :: norm

                                 CALL timeset(routineN, handle2)
                                 out_unit = cp_logger_get_default_io_unit()
                                 CALL timeset(routineN//"_"//TRIM(ADJUSTL(cp_to_string( &
                                                                          CEILING(pw1%pw_grid%cutoff/10)*10))), handle)

                                 NULLIFY (c_in)
                                 NULLIFY (c_out)

                                 IF (PRESENT(debug)) THEN
                                    test = debug
                                 ELSE
                                    test = .FALSE.
                                 END IF

                                 !..check if grids are compatible
                                 IF (.NOT. ASSOCIATED(pw1%pw_grid, pw2%pw_grid)) THEN
                                    IF (pw1%pw_grid%dvol /= pw2%pw_grid%dvol) THEN
                                       CPABORT("PW grids not compatible")
                                    END IF
                                    IF (pw1%pw_grid%para%group /= pw2%pw_grid%para%group) THEN
                                       CPABORT("PW grids have not compatible MPI groups")
                                    END IF
                                 END IF

                                 !..prepare input
                                 #:if space == "rs"
                                    norm = 1.0_dp/pw1%pw_grid%ngpts
                                 #:elif space == "gs"
                                    norm = 1.0_dp
                                 #:else
                                    #:stop "Unknown space tag: "+space
                                 #:endif

                                 n => pw1%pw_grid%npts

                                 IF (pw1%pw_grid%para%mode == PW_MODE_LOCAL) THEN

                                    !
                                    !..replicated data, use local FFT
                                    !

                                    IF (test .AND. out_unit > 0) THEN
                                       WRITE (out_unit, '(A)') " FFT Protocol "
                                       #:if space=="rs"
                                          WRITE (out_unit, '(A,T76,A)') "  Transform direction ", "FWFFT"
                                          WRITE (out_unit, '(A,T72,A)') "  in space ", "REALSPACE"
                                          WRITE (out_unit, '(A,T72,A)') "  out space ", "REALSPACE"
                                       #:else
                                          WRITE (out_unit, '(A,T76,A)') "  Transform direction ", "BWFFT"
                                          WRITE (out_unit, '(A,T66,A)') "  in space ", "RECIPROCALSPACE"
                                          WRITE (out_unit, '(A,T66,A)') "  out space ", "RECIPROCALSPACE"
                                       #:endif
                                       WRITE (out_unit, '(A,T66,E15.6)') "  scale factor", norm
                                    END IF

                                    #:if space=="rs"
                                       #:if kind==kind2=="c3d"
                                          c_in => pw1%array
                                          c_out => pw2%array
                                          CALL fft3d(FWFFT, n, c_in, c_out, scale=norm, debug=test)
                                       #:elif kind=="r3d" and kind2=="c3d"
                                          pw2%array = CMPLX(pw1%array, 0.0_dp, KIND=dp)
                                          c_out => pw2%array
                                          CALL fft3d(FWFFT, n, c_out, scale=norm, debug=test)
                                       #:elif kind=="c3d" and kind2=="c1d"
                                          c_in => pw1%array
                                          ALLOCATE (c_out(n(1), n(2), n(3)))
                                          ! transform
                                          CALL fft3d(FWFFT, n, c_in, c_out, scale=norm, debug=test)
                                          ! gather results
                                          IF (test .AND. out_unit > 0) WRITE (out_unit, '(A)') "  PW_GATHER : 3d -> 1d "
                                          CALL pw_gather_s_${kind2}$_c3d(pw2, c_out)
                                          DEALLOCATE (c_out)
                                       #:elif kind=="r3d" and kind2=="c1d"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                                          CALL pw_gpu_r3dc1d_3d(pw1, pw2, scale=norm)
#elif defined (__PW_FPGA)
                                          ALLOCATE (c_out(n(1), n(2), n(3)))
                                          ! check if bitstream for the fft size is present
                                          ! if not, perform fft3d in CPU
                                          IF (pw_fpga_init_bitstream(n) == 1) THEN
                                             CALL pw_copy_to_array(pw1, c_out)
#if (__PW_FPGA_SP && __PW_FPGA)
                                             CALL pw_fpga_r3dc1d_3d_sp(n, c_out)
#else
                                             CALL pw_fpga_r3dc1d_3d_dp(n, c_out)
#endif
                                             CALL zdscal(n(1)*n(2)*n(3), norm, c_out, 1)
                                             CALL pw_gather_s_${kind2}$_c3d(pw2, c_out)
                                          ELSE
                                             CALL pw_copy_to_array(pw1, c_out)
                                             CALL fft3d(FWFFT, n, c_out, scale=norm, debug=test)
                                             CALL pw_gather_s_${kind2}$_c3d(pw2, c_out)
                                          END IF
                                          DEALLOCATE (c_out)
#else
                                          ALLOCATE (c_out(n(1), n(2), n(3)))
                                          c_out = 0.0_dp
                                          CALL pw_copy_to_array(pw1, c_out)
                                          CALL fft3d(FWFFT, n, c_out, scale=norm, debug=test)
                                          CALL pw_gather_s_${kind2}$_c3d(pw2, c_out)
                                          DEALLOCATE (c_out)
#endif
                                       #:endif
                                    #:else
                                       #:if kind=="c3d" and kind2=="c3d"
                                          c_in => pw1%array
                                          c_out => pw2%array
                                          CALL fft3d(BWFFT, n, c_in, c_out, scale=norm, debug=test)
                                       #:elif kind=="c3d" and kind2=="r3d"
                                          c_in => pw1%array
                                          ALLOCATE (c_out(n(1), n(2), n(3)))
                                          CALL fft3d(BWFFT, n, c_in, c_out, scale=norm, debug=test)
                                          ! use real part only
                                          IF (test .AND. out_unit > 0) WRITE (out_unit, '(A)') "  REAL part "
                                          pw2%array = REAL(c_out, KIND=dp)
                                          DEALLOCATE (c_out)
                                       #:elif kind=="c1d" and kind2=="c3d"
                                          c_out => pw2%array
                                          IF (test .AND. out_unit > 0) WRITE (out_unit, '(A)') "  PW_SCATTER : 3d -> 1d "
                                          CALL pw_scatter_s_${kind}$_c3d(pw1, c_out)
                                          CALL fft3d(BWFFT, n, c_out, scale=norm, debug=test)
                                       #:elif kind=="c1d" and kind2=="r3d"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                                          CALL pw_gpu_c1dr3d_3d(pw1, pw2, scale=norm)
#elif defined (__PW_FPGA)
                                          ALLOCATE (c_out(n(1), n(2), n(3)))
                                          ! check if bitstream for the fft size is present
                                          ! if not, perform fft3d in CPU
                                          IF (pw_fpga_init_bitstream(n) == 1) THEN
                                             CALL pw_scatter_s_${kind}$_c3d(pw1, c_out)
                                             ! transform using FPGA
#if (__PW_FPGA_SP && __PW_FPGA)
                                             CALL pw_fpga_c1dr3d_3d_sp(n, c_out)
#else
                                             CALL pw_fpga_c1dr3d_3d_dp(n, c_out)
#endif
                                             CALL zdscal(n(1)*n(2)*n(3), norm, c_out, 1)
                                             ! use real part only
                                             CALL pw_copy_from_array(pw2, c_out)
                                          ELSE
                                             IF (test .AND. out_unit > 0) WRITE (out_unit, '(A)') "  PW_SCATTER : 3d -> 1d "
                                             CALL pw_scatter_s_${kind}$_c3d(pw1, c_out)
                                             ! transform
                                             CALL fft3d(BWFFT, n, c_out, scale=norm, debug=test)
                                             ! use real part only
                                             IF (test .AND. out_unit > 0) WRITE (out_unit, '(A)') "  REAL part "
                                             CALL pw_copy_from_array(pw2, c_out)
                                          END IF
                                          DEALLOCATE (c_out)
#else
                                          ALLOCATE (c_out(n(1), n(2), n(3)))
                                          IF (test .AND. out_unit > 0) WRITE (out_unit, '(A)') "  PW_SCATTER : 3d -> 1d "
                                          CALL pw_scatter_s_${kind}$_c3d(pw1, c_out)
                                          ! transform
                                          CALL fft3d(BWFFT, n, c_out, scale=norm, debug=test)
                                          ! use real part only
                                          IF (test .AND. out_unit > 0) WRITE (out_unit, '(A)') "  REAL part "
                                          CALL pw_copy_from_array(pw2, c_out)
                                          DEALLOCATE (c_out)
#endif
                                       #:endif
                                    #:endif

                                    IF (test .AND. out_unit > 0) WRITE (out_unit, '(A)') " End of FFT Protocol "

                                 ELSE

                                    !
                                    !..parallel FFT
                                    !

                                    IF (test .AND. out_unit > 0) THEN
                                       WRITE (out_unit, '(A)') " FFT Protocol "
                                       #:if space=="rs"
                                          WRITE (out_unit, '(A,T76,A)') "  Transform direction ", "FWFFT"
                                          WRITE (out_unit, '(A,T72,A)') "  in space ", "REALSPACE"
                                          WRITE (out_unit, '(A,T72,A)') "  out space ", "REALSPACE"
                                       #:else
                                          WRITE (out_unit, '(A,T76,A)') "  Transform direction ", "BWFFT"
                                          WRITE (out_unit, '(A,T66,A)') "  in space ", "RECIPROCALSPACE"
                                          WRITE (out_unit, '(A,T66,A)') "  out space ", "RECIPROCALSPACE"
                                       #:endif
                                       WRITE (out_unit, '(A,T66,E15.6)') "  scale factor", norm
                                    END IF

                                    my_pos = pw1%pw_grid%para%group%mepos
                                    nrays = pw1%pw_grid%para%nyzray(my_pos)
                                    grays => pw1%pw_grid%grays

                                    #:if space=="rs"
                                       #:if kind=="c3d" and kind2=="c1d"
                                          !..prepare input
                                          c_in => pw1%array
                                          grays = z_zero
                                          !..transform
                                          IF (pw1%pw_grid%para%ray_distribution) THEN
                                             CALL fft3d(FWFFT, n, c_in, grays, pw1%pw_grid%para%group, &
                                                        pw1%pw_grid%para%yzp, pw1%pw_grid%para%nyzray, &
                                                        pw1%pw_grid%para%bo, scale=norm, debug=test)
                                          ELSE
                                             CALL fft3d(FWFFT, n, c_in, grays, pw1%pw_grid%para%group, &
                                                        pw1%pw_grid%para%bo, scale=norm, debug=test)
                                          END IF
                                          !..prepare output
                                          IF (test .AND. out_unit > 0) &
                                             WRITE (out_unit, '(A)') "  PW_GATHER : 2d -> 1d "
                                          CALL pw_gather_p_${kind2}$ (pw2, grays)
                                       #:elif kind=="r3d" and kind2=="c1d"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                                          ! (no ray dist. is not efficient in CUDA)
                                          use_pw_gpu = pw1%pw_grid%para%ray_distribution
                                          IF (use_pw_gpu) THEN
                                             CALL pw_gpu_r3dc1d_3d_ps(pw1, pw2, scale=norm)
                                          ELSE
#endif
!..   prepare input
                                             nloc = pw1%pw_grid%npts_local
                                             ALLOCATE (c_in(nloc(1), nloc(2), nloc(3)))
                                             CALL pw_copy_to_array(pw1, c_in)
                                             grays = z_zero
                                             !..transform
                                             IF (pw1%pw_grid%para%ray_distribution) THEN
                                                CALL fft3d(FWFFT, n, c_in, grays, pw1%pw_grid%para%group, &
                                                           pw1%pw_grid%para%yzp, pw1%pw_grid%para%nyzray, &
                                                           pw1%pw_grid%para%bo, scale=norm, debug=test)
                                             ELSE
                                                CALL fft3d(FWFFT, n, c_in, grays, pw1%pw_grid%para%group, &
                                                           pw1%pw_grid%para%bo, scale=norm, debug=test)
                                             END IF
                                             !..prepare output
                                             IF (test .AND. out_unit > 0) &
                                                WRITE (out_unit, '(A)') "  PW_GATHER : 2d -> 1d "
                                             CALL pw_gather_p_${kind2}$ (pw2, grays)
                                             DEALLOCATE (c_in)

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                                          END IF
#endif
                                       #:endif
                                    #:else
                                       #:if kind=="c1d" and kind2=="c3d"
                                          !..prepare input
                                          IF (test .AND. out_unit > 0) &
                                             WRITE (out_unit, '(A)') "  PW_SCATTER : 2d -> 1d "
                                          grays = z_zero
                                          CALL pw_scatter_p_${kind}$ (pw1, grays)
                                          c_in => pw2%array
                                          !..transform
                                          IF (pw1%pw_grid%para%ray_distribution) THEN
                                             CALL fft3d(BWFFT, n, c_in, grays, pw1%pw_grid%para%group, &
                                                        pw1%pw_grid%para%yzp, pw1%pw_grid%para%nyzray, &
                                                        pw1%pw_grid%para%bo, scale=norm, debug=test)
                                          ELSE
                                             CALL fft3d(BWFFT, n, c_in, grays, pw1%pw_grid%para%group, &
                                                        pw1%pw_grid%para%bo, scale=norm, debug=test)
                                          END IF
                                          !..prepare output (nothing to do)
                                       #:elif kind=="c1d" and kind2=="r3d"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                                          ! (no ray dist. is not efficient in CUDA)
                                          use_pw_gpu = pw1%pw_grid%para%ray_distribution
                                          IF (use_pw_gpu) THEN
                                             CALL pw_gpu_c1dr3d_3d_ps(pw1, pw2, scale=norm)
                                          ELSE
#endif
!..   prepare input
                                             IF (test .AND. out_unit > 0) &
                                                WRITE (out_unit, '(A)') "  PW_SCATTER : 2d -> 1d "
                                             grays = z_zero
                                             CALL pw_scatter_p_${kind}$ (pw1, grays)
                                             nloc = pw2%pw_grid%npts_local
                                             ALLOCATE (c_in(nloc(1), nloc(2), nloc(3)))
                                             !..transform
                                             IF (pw1%pw_grid%para%ray_distribution) THEN
                                                CALL fft3d(BWFFT, n, c_in, grays, pw1%pw_grid%para%group, &
                                                           pw1%pw_grid%para%yzp, pw1%pw_grid%para%nyzray, &
                                                           pw1%pw_grid%para%bo, scale=norm, debug=test)
                                             ELSE
                                                CALL fft3d(BWFFT, n, c_in, grays, pw1%pw_grid%para%group, &
                                                           pw1%pw_grid%para%bo, scale=norm, debug=test)
                                             END IF
                                             !..prepare output
                                             IF (test .AND. out_unit > 0) &
                                                WRITE (out_unit, '(A)') "  Real part "
                                             CALL pw_copy_from_array(pw2, c_in)
                                             DEALLOCATE (c_in)
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
                                          END IF
#endif
                                       #:endif
                                    #:endif
                                 END IF

                                 IF (test .AND. out_unit > 0) THEN
                                    WRITE (out_unit, '(A)') " End of FFT Protocol "
                                 END IF

                                 CALL timestop(handle)
                                 CALL timestop(handle2)

                              END SUBROUTINE fft_wrap_pw1pw2_${kind}$_${kind2}$_${space}$_${space2}$
                           #:endif
                        #:endfor
                     #:endfor

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
                        SUBROUTINE pw_gather_s_${kind}$_${kind2}$_2(pw1, pw2, scale)

                           TYPE(pw_${kind2}$_gs_type), INTENT(IN)                       :: pw1
                           TYPE(pw_${kind}$_gs_type), INTENT(INOUT)   :: pw2
                           REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

                           CALL pw_gather_s_${kind}$_${kind2}$ (pw2, pw1%array, scale)

                        END SUBROUTINE pw_gather_s_${kind}$_${kind2}$_2

! **************************************************************************************************
!> \brief Gathers the pw vector from a 3d data field
!> \param pw ...
!> \param c ...
!> \param scale ...
!> \par History
!>      none
!> \author JGH
! **************************************************************************************************
                        SUBROUTINE pw_gather_s_${kind}$_${kind2}$ (pw, c, scale)

                           TYPE(pw_${kind}$_gs_type), INTENT(INOUT)                       :: pw
                           ${type2}$, CONTIGUOUS, INTENT(IN)   :: c
                           REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

                           CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gather_s'

                           INTEGER                                            :: gpt, handle, l, m, n

                           CALL timeset(routineN, handle)

                           ASSOCIATE (mapl => pw%pw_grid%mapl%pos, mapm => pw%pw_grid%mapm%pos, mapn => pw%pw_grid%mapn%pos, &
                                      ngpts => SIZE(pw%pw_grid%gsq), ghat => pw%pw_grid%g_hat)

                              IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(c, pw, scale)
                                 DO gpt = 1, ngpts
                                    l = mapl(ghat(1, gpt)) + 1
                                    m = mapm(ghat(2, gpt)) + 1
                                    n = mapn(ghat(3, gpt)) + 1
                                    pw%array(gpt) = scale*${type2type("c(l, m, n)", kind2, kind)}$
                                 END DO
!$OMP END PARALLEL DO
                              ELSE
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(c, pw)
                                 DO gpt = 1, ngpts
                                    l = mapl(ghat(1, gpt)) + 1
                                    m = mapm(ghat(2, gpt)) + 1
                                    n = mapn(ghat(3, gpt)) + 1
                                    pw%array(gpt) = ${type2type("c(l, m, n)", kind2, kind)}$
                                 END DO
!$OMP END PARALLEL DO
                              END IF

                           END ASSOCIATE

                           CALL timestop(handle)

                        END SUBROUTINE pw_gather_s_${kind}$_${kind2}$

! **************************************************************************************************
!> \brief Scatters a pw vector to a 3d data field
!> \param pw ...
!> \param c ...
!> \param scale ...
!> \par History
!>      none
!> \author JGH
! **************************************************************************************************
                        SUBROUTINE pw_scatter_s_${kind}$_${kind2}$_2(pw1, pw2, scale)

                           TYPE(pw_${kind}$_gs_type), INTENT(IN)                 :: pw1
                           TYPE(pw_${kind2}$_gs_type), INTENT(INOUT)               :: pw2
                           REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

                           CALL pw_scatter_s_${kind}$_${kind2}$ (pw1, pw2%array, scale)

                        END SUBROUTINE pw_scatter_s_${kind}$_${kind2}$_2

! **************************************************************************************************
!> \brief Scatters a pw vector to a 3d data field
!> \param pw ...
!> \param c ...
!> \param scale ...
!> \par History
!>      none
!> \author JGH
! **************************************************************************************************
                        SUBROUTINE pw_scatter_s_${kind}$_${kind2}$ (pw, c, scale)

                           TYPE(pw_${kind}$_gs_type), INTENT(IN)                 :: pw
                           ${type2}$, CONTIGUOUS, INTENT(INOUT)               :: c
                           REAL(KIND=dp), INTENT(IN), OPTIONAL                :: scale

                           CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_scatter_s'

                           INTEGER                                            :: gpt, handle, l, m, n

                           CALL timeset(routineN, handle)

                           ASSOCIATE (mapl => pw%pw_grid%mapl%pos, mapm => pw%pw_grid%mapm%pos, mapn => pw%pw_grid%mapn%pos, &
                                      ghat => pw%pw_grid%g_hat, ngpts => SIZE(pw%pw_grid%gsq))

                              ! should only zero the unused bits (but the zero is needed)
                              IF (.NOT. PRESENT(scale)) c = 0.0_dp

                              IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(c, pw, scale)
                                 DO gpt = 1, ngpts
                                    l = mapl(ghat(1, gpt)) + 1
                                    m = mapm(ghat(2, gpt)) + 1
                                    n = mapn(ghat(3, gpt)) + 1
                                    c(l, m, n) = scale*${type2type("pw%array(gpt)", kind, kind2)}$
                                 END DO
!$OMP END PARALLEL DO
                              ELSE
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(c, pw)
                                 DO gpt = 1, ngpts
                                    l = mapl(ghat(1, gpt)) + 1
                                    m = mapm(ghat(2, gpt)) + 1
                                    n = mapn(ghat(3, gpt)) + 1
                                    c(l, m, n) = ${type2type("pw%array(gpt)", kind, kind2)}$
                                 END DO
!$OMP END PARALLEL DO
                              END IF

                           END ASSOCIATE

                           IF (pw%pw_grid%grid_span == HALFSPACE) THEN

                              ASSOCIATE (mapl => pw%pw_grid%mapl%neg, mapm => pw%pw_grid%mapm%neg, mapn => pw%pw_grid%mapn%neg, &
                                         ghat => pw%pw_grid%g_hat, ngpts => SIZE(pw%pw_grid%gsq))

                                 IF (PRESENT(scale)) THEN
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(c, pw, scale)
                                    DO gpt = 1, ngpts
                                       l = mapl(ghat(1, gpt)) + 1
                                       m = mapm(ghat(2, gpt)) + 1
                                       n = mapn(ghat(3, gpt)) + 1
                                       c(l, m, n) = scale*#{if kind[0]=="c" and kind2[0]=="c"}#CONJG#{endif}#(${type2type("pw%array(gpt)", kind, kind2)}$)
                                    END DO
!$OMP END PARALLEL DO
                                 ELSE
!$OMP PARALLEL DO PRIVATE(gpt, l, m, n) DEFAULT(NONE) SHARED(c, pw)
                                    DO gpt = 1, ngpts
                                       l = mapl(ghat(1, gpt)) + 1
                                       m = mapm(ghat(2, gpt)) + 1
                                       n = mapn(ghat(3, gpt)) + 1
                                       c(l, m, n) = #{if kind[0]=="c" and kind2[0]=="c"}#CONJG#{endif}#(${type2type("pw%array(gpt)", kind, kind2)}$)
                                    END DO
!$OMP END PARALLEL DO
                                 END IF

                              END ASSOCIATE

                           END IF

                           CALL timestop(handle)

                        END SUBROUTINE pw_scatter_s_${kind}$_${kind2}$
                     #:endif
                  #:endfor

! **************************************************************************************************
!> \brief Multiply all data points with a Gaussian damping factor
!>        Needed for longrange Coulomb potential
!>        V(\vec r)=erf(omega*r)/r
!>        V(\vec g)=\frac{4*\pi}{g**2}*exp(-g**2/omega**2)
!> \param pw ...
!> \param omega ...
!> \par History
!>      Frederick Stein (12-04-2019) created
!> \author Frederick Stein (12-Apr-2019)
!> \note
!>      Performs a Gaussian damping
! **************************************************************************************************
                  SUBROUTINE pw_gauss_damp(pw, omega)

                     TYPE(pw_c1d_gs_type), INTENT(INOUT)                       :: pw
                     REAL(KIND=dp), INTENT(IN)                          :: omega

                     CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gauss_damp'

                     INTEGER                                            :: handle
                     REAL(KIND=dp)                                      :: omega_2

                     CALL timeset(routineN, handle)
                     CPASSERT(omega >= 0)

                     omega_2 = omega*omega
                     omega_2 = 0.25_dp/omega_2

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw,omega_2)
                     pw%array = pw%array*EXP(-pw%pw_grid%gsq*omega_2)
!$OMP END PARALLEL WORKSHARE

                     CALL timestop(handle)

                  END SUBROUTINE pw_gauss_damp

! **************************************************************************************************
!> \brief Multiply all data points with the logarithmic derivative of a Gaussian
!> \param pw ...
!> \param omega ...
!> \note
!>      Performs a Gaussian damping
! **************************************************************************************************
                  SUBROUTINE pw_log_deriv_gauss(pw, omega)

                     TYPE(pw_c1d_gs_type), INTENT(IN)                          :: pw
                     REAL(KIND=dp), INTENT(IN)                          :: omega

                     CHARACTER(len=*), PARAMETER :: routineN = 'pw_log_deriv_gauss'

                     INTEGER                                            :: handle
                     REAL(KIND=dp)                                      :: omega_2

                     CALL timeset(routineN, handle)
                     CPASSERT(omega >= 0)

                     omega_2 = omega*omega
                     omega_2 = 0.25_dp/omega_2

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw,omega_2)
                     pw%array = pw%array*(1.0_dp + omega_2*pw%pw_grid%gsq)
!$OMP END PARALLEL WORKSHARE

                     CALL timestop(handle)
                  END SUBROUTINE pw_log_deriv_gauss

! **************************************************************************************************
!> \brief Multiply all data points with a Gaussian damping factor
!>        Needed for longrange Coulomb potential
!>        V(\vec r)=erf(omega*r)/r
!>        V(\vec g)=\frac{4*\pi}{g**2}*exp(-g**2/omega**2)
!> \param pw ...
!> \param omega ...
!> \par History
!>      Frederick Stein (12-04-2019) created
!> \author Frederick Stein (12-Apr-2019)
!> \note
!>      Performs a Gaussian damping
! **************************************************************************************************
                  SUBROUTINE pw_compl_gauss_damp(pw, omega)

                     TYPE(pw_c1d_gs_type), INTENT(INOUT)                       :: pw
                     REAL(KIND=dp), INTENT(IN)                          :: omega

                     CHARACTER(len=*), PARAMETER :: routineN = 'pw_compl_gauss_damp'

                     INTEGER                                            :: cnt, handle, i
                     REAL(KIND=dp)                                      :: omega_2, tmp

                     CALL timeset(routineN, handle)

                     omega_2 = omega*omega
                     omega_2 = 0.25_dp/omega_2

                     cnt = SIZE(pw%array)

!$OMP PARALLEL DO PRIVATE(i, tmp) DEFAULT(NONE) SHARED(cnt, pw,omega_2)
                     DO i = 1, cnt
                        tmp = -omega_2*pw%pw_grid%gsq(i)
                        IF (ABS(tmp) > 1.0E-5_dp) THEN
                           pw%array(i) = pw%array(i)*(1.0_dp - EXP(tmp))
                        ELSE
                           pw%array(i) = pw%array(i)*(tmp + 0.5_dp*tmp*(tmp + (1.0_dp/3.0_dp)*tmp**2))
                        END IF
                     END DO
!$OMP END PARALLEL DO

                     CALL timestop(handle)

                  END SUBROUTINE pw_compl_gauss_damp

! **************************************************************************************************
!> \brief Multiply all data points with the logarithmic derivative of the complementary Gaussian damping factor
!> \param pw ...
!> \param omega ...
!> \note
! **************************************************************************************************
                  SUBROUTINE pw_log_deriv_compl_gauss(pw, omega)

                     TYPE(pw_c1d_gs_type), INTENT(IN)                          :: pw
                     REAL(KIND=dp), INTENT(IN)                          :: omega

                     CHARACTER(len=*), PARAMETER :: routineN = 'pw_log_deriv_compl_gauss'

                     INTEGER                                            :: handle, i
                     REAL(KIND=dp)                                      :: omega_2, tmp

                     CALL timeset(routineN, handle)

                     omega_2 = omega*omega
                     omega_2 = 0.25_dp/omega_2

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,tmp) &
!$OMP             SHARED(pw,omega_2)
                     DO i = 1, SIZE(pw%array)
                        tmp = omega_2*pw%pw_grid%gsq(i)
                        ! For too small arguments, use the Taylor polynomial to prevent division by zero
                        IF (ABS(tmp) >= 0.003_dp) THEN
                           pw%array(i) = pw%array(i)*(1.0_dp - tmp*EXP(-tmp)/(1.0_dp - EXP(-tmp)))
                        ELSE
                           pw%array(i) = pw%array(i)*(0.5_dp*tmp - tmp**2/12.0_dp)
                        END IF
                     END DO
!$OMP END PARALLEL DO

                     CALL timestop(handle)

                  END SUBROUTINE pw_log_deriv_compl_gauss

! **************************************************************************************************
!> \brief Multiply all data points with a Gaussian damping factor and mixes it with the original function
!>        Needed for mixed longrange/Coulomb potential
!>        V(\vec r)=(a+b*erf(omega*r))/r
!>        V(\vec g)=\frac{4*\pi}{g**2}*(a+b*exp(-g**2/omega**2))
!> \param pw ...
!> \param omega ...
!> \param scale_coul ...
!> \param scale_long ...
!> \par History
!>      Frederick Stein (16-Dec-2021) created
!> \author Frederick Stein (16-Dec-2021)
!> \note
!>      Performs a Gaussian damping
! **************************************************************************************************
                  SUBROUTINE pw_gauss_damp_mix(pw, omega, scale_coul, scale_long)

                     TYPE(pw_c1d_gs_type), INTENT(INOUT)                       :: pw
                     REAL(KIND=dp), INTENT(IN)                          :: omega, scale_coul, scale_long

                     CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gauss_damp_mix'

                     INTEGER                                            :: handle
                     REAL(KIND=dp)                                      :: omega_2

                     CALL timeset(routineN, handle)

                     omega_2 = omega*omega
                     omega_2 = 0.25_dp/omega_2

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw, omega_2, scale_coul, scale_long)
                     pw%array = pw%array*(scale_coul + scale_long*EXP(-pw%pw_grid%gsq*omega_2))
!$OMP END PARALLEL WORKSHARE

                     CALL timestop(handle)

                  END SUBROUTINE pw_gauss_damp_mix

! **************************************************************************************************
!> \brief Multiply all data points with the logarithmic derivative of the mixed longrange/Coulomb potential
!>        Needed for mixed longrange/Coulomb potential
!> \param pw ...
!> \param omega ...
!> \param scale_coul ...
!> \param scale_long ...
!> \note
! **************************************************************************************************
                  SUBROUTINE pw_log_deriv_mix_cl(pw, omega, scale_coul, scale_long)

                     TYPE(pw_c1d_gs_type), INTENT(IN)                          :: pw
                     REAL(KIND=dp), INTENT(IN)                          :: omega, scale_coul, scale_long

                     CHARACTER(len=*), PARAMETER :: routineN = 'pw_log_deriv_mix_cl'

                     INTEGER                                            :: handle, i
                     REAL(KIND=dp)                                      :: omega_2, tmp

                     CALL timeset(routineN, handle)

                     omega_2 = omega*omega
                     omega_2 = 0.25_dp/omega_2

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,tmp) &
!$OMP             SHARED(pw,omega_2,scale_long,scale_coul)
                     DO i = 1, SIZE(pw%array)
                        tmp = omega_2*pw%pw_grid%gsq(i)
                        pw%array(i) = pw%array(i)*(1.0_dp + scale_long*tmp*EXP(-tmp)/(scale_coul + scale_long*EXP(-tmp)))
                     END DO
!$OMP END PARALLEL DO

                     CALL timestop(handle)

                  END SUBROUTINE pw_log_deriv_mix_cl

! **************************************************************************************************
!> \brief Multiply all data points with a complementary cosine
!>        Needed for truncated Coulomb potential
!>        V(\vec r)=1/r if r<rc, 0 else
!>        V(\vec g)=\frac{4*\pi}{g**2}*(1-cos(g*rc))
!> \param pw ...
!> \param rcutoff ...
!> \par History
!>      Frederick Stein (07-06-2021) created
!> \author Frederick Stein (07-Jun-2021)
!> \note
!>      Multiplies by complementary cosine
! **************************************************************************************************
                  SUBROUTINE pw_truncated(pw, rcutoff)

                     TYPE(pw_c1d_gs_type), INTENT(INOUT)                       :: pw
                     REAL(KIND=dp), INTENT(IN)                          :: rcutoff

                     CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_truncated'

                     INTEGER                                            :: handle, i
                     REAL(KIND=dp)                                      :: tmp

                     CALL timeset(routineN, handle)
                     CPASSERT(rcutoff >= 0)

!$OMP PARALLEL DO PRIVATE(i,tmp) DEFAULT(NONE) SHARED(pw, rcutoff)
                     DO i = 1, SIZE(pw%array)
                        tmp = SQRT(pw%pw_grid%gsq(i))*rcutoff
                        IF (tmp >= 0.005_dp) THEN
                           pw%array(i) = pw%array(i)*(1.0_dp - COS(tmp))
                        ELSE
                           pw%array(i) = pw%array(i)*tmp**2/2.0_dp*(1.0 - tmp**2/12.0_dp)
                        END IF
                     END DO
!$OMP END PARALLEL DO

                     CALL timestop(handle)

                  END SUBROUTINE pw_truncated

! **************************************************************************************************
!> \brief Multiply all data points with the logarithmic derivative of the complementary cosine
!>        This function is needed for virials using truncated Coulomb potentials
!> \param pw ...
!> \param rcutoff ...
!> \note
! **************************************************************************************************
                  SUBROUTINE pw_log_deriv_trunc(pw, rcutoff)

                     TYPE(pw_c1d_gs_type), INTENT(IN)                          :: pw
                     REAL(KIND=dp), INTENT(IN)                          :: rcutoff

                     CHARACTER(len=*), PARAMETER :: routineN = 'pw_log_deriv_trunc'

                     INTEGER                                            :: handle, i
                     REAL(KIND=dp)                                      :: rchalf, tmp

                     CALL timeset(routineN, handle)
                     CPASSERT(rcutoff >= 0)

                     rchalf = 0.5_dp*rcutoff
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,tmp) &
!$OMP             SHARED(pw,rchalf)
                     DO i = 1, SIZE(pw%array)
                        tmp = rchalf*SQRT(pw%pw_grid%gsq(i))
                        ! For too small arguments, use the Taylor polynomial to prevent division by zero
                        IF (ABS(tmp) >= 0.0001_dp) THEN
                           pw%array(i) = pw%array(i)*(1.0_dp - tmp/TAN(tmp))
                        ELSE
                           pw%array(i) = pw%array(i)*tmp**2*(1.0_dp + tmp**2/15.0_dp)/3.0_dp
                        END IF
                     END DO
!$OMP END PARALLEL DO

                     CALL timestop(handle)

                  END SUBROUTINE pw_log_deriv_trunc

! **************************************************************************************************
!> \brief Calculate the derivative of a plane wave vector
!> \param pw ...
!> \param n ...
!> \par History
!>      JGH (06-10-2002) allow only for inplace derivatives
!> \author JGH (25-Feb-2001)
!> \note
!>      Calculate the derivative dx^n(1) dy^n(2) dz^n(3) PW
! **************************************************************************************************
                  SUBROUTINE pw_derive(pw, n)

                     TYPE(pw_c1d_gs_type), INTENT(INOUT)                       :: pw
                     INTEGER, DIMENSION(3), INTENT(IN)                  :: n

                     CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_derive'

                     COMPLEX(KIND=dp)                                   :: im
                     INTEGER                                            :: handle, m

                     CALL timeset(routineN, handle)

                     IF (ANY(n < 0)) &
                        CPABORT("Nonnegative exponents are not supported!")

                     m = SUM(n)
                     im = CMPLX(0.0_dp, 1.0_dp, KIND=dp)**m

                     IF (n(1) == 1) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw)
                        pw%array(:) = pw%array(:)*pw%pw_grid%g(1, :)
!$OMP END PARALLEL WORKSHARE
                     ELSE IF (n(1) > 1) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(n, pw)
                        pw%array(:) = pw%array(:)*(pw%pw_grid%g(1, :)**n(1))
!$OMP END PARALLEL WORKSHARE
                     END IF

                     IF (n(2) == 1) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw)
                        pw%array(:) = pw%array(:)*pw%pw_grid%g(2, :)
!$OMP END PARALLEL WORKSHARE
                     ELSE IF (n(2) > 1) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(n, pw)
                        pw%array(:) = pw%array(:)*(pw%pw_grid%g(2, :)**n(2))
!$OMP END PARALLEL WORKSHARE
                     END IF

                     IF (n(3) == 1) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw)
                        pw%array(:) = pw%array(:)*pw%pw_grid%g(3, :)
!$OMP END PARALLEL WORKSHARE
                     ELSE IF (n(3) > 1) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(n, pw)
                        pw%array(:) = pw%array(:)*(pw%pw_grid%g(3, :)**n(3))
!$OMP END PARALLEL WORKSHARE
                     END IF

                     ! im can take the values 1, -1, i, -i
                     ! skip this if im == 1
                     IF (ABS(REAL(im, KIND=dp) - 1.0_dp) > 1.0E-10_dp) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(im, pw)
                        pw%array(:) = im*pw%array(:)
!$OMP END PARALLEL WORKSHARE
                     END IF

                     CALL timestop(handle)

                  END SUBROUTINE pw_derive

! **************************************************************************************************
!> \brief Calculate the Laplacian of a plane wave vector
!> \param pw ...
!> \par History
!>      Frederick Stein (01-02-2022) created
!> \author JGH (25-Feb-2001)
!> \note
!>      Calculate the derivative DELTA PW
! **************************************************************************************************
                  SUBROUTINE pw_laplace(pw)

                     TYPE(pw_c1d_gs_type), INTENT(INOUT)                       :: pw

                     CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_laplace'

                     INTEGER                                            :: handle

                     CALL timeset(routineN, handle)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw)
                     pw%array(:) = -pw%array(:)*pw%pw_grid%gsq(:)
!$OMP END PARALLEL WORKSHARE

                     CALL timestop(handle)

                  END SUBROUTINE pw_laplace

! **************************************************************************************************
!> \brief Calculate the tensorial 2nd derivative of a plane wave vector
!> \param pw ...
!> \param pwdr2 ...
!> \param i ...
!> \param j ...
!> \par History
!>      none
!> \author JGH (05-May-2006)
!> \note
! **************************************************************************************************
                  SUBROUTINE pw_dr2(pw, pwdr2, i, j)

                     TYPE(pw_c1d_gs_type), INTENT(INOUT)                       :: pw, pwdr2
                     INTEGER, INTENT(IN)                                :: i, j

                     CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_dr2'

                     INTEGER                                            :: cnt, handle, ig
                     REAL(KIND=dp)                                      :: gg, o3

                     CALL timeset(routineN, handle)

                     o3 = 1.0_dp/3.0_dp

                     cnt = SIZE(pw%array)

                     IF (i == j) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ig,gg) SHARED(cnt, i, o3, pw, pwdr2)
                        DO ig = 1, cnt
                           gg = pw%pw_grid%g(i, ig)**2 - o3*pw%pw_grid%gsq(ig)
                           pwdr2%array(ig) = gg*pw%array(ig)
                        END DO
!$OMP END PARALLEL DO
                     ELSE
!$OMP PARALLEL DO PRIVATE (ig) DEFAULT(NONE) SHARED(cnt, i, j, pw, pwdr2)
                        DO ig = 1, cnt
                           pwdr2%array(ig) = pw%array(ig)*(pw%pw_grid%g(i, ig)*pw%pw_grid%g(j, ig))
                        END DO
!$OMP END PARALLEL DO
                     END IF

                     CALL timestop(handle)

                  END SUBROUTINE pw_dr2

! **************************************************************************************************
!> \brief Calculate the tensorial 2nd derivative of a plane wave vector
!>      and divides by |G|^2. pwdr2_gg(G=0) is put to zero.
!> \param pw ...
!> \param pwdr2_gg ...
!> \param i ...
!> \param j ...
!> \par History
!>      none
!> \author RD (20-Nov-2006)
!> \note
!>      Adapted from pw_dr2
! **************************************************************************************************
                  SUBROUTINE pw_dr2_gg(pw, pwdr2_gg, i, j)

                     TYPE(pw_c1d_gs_type), INTENT(INOUT)                       :: pw, pwdr2_gg
                     INTEGER, INTENT(IN)                                :: i, j

                     INTEGER                                            :: cnt, handle, ig
                     REAL(KIND=dp)                                      :: gg, o3
                     CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_dr2_gg'

                     CALL timeset(routineN, handle)

                     o3 = 1.0_dp/3.0_dp

                     cnt = SIZE(pw%array)

                     IF (i == j) THEN
!$OMP PARALLEL DO PRIVATE (ig) DEFAULT(NONE) PRIVATE(gg) SHARED(cnt, i, o3, pw, pwdr2_gg)
                        DO ig = pw%pw_grid%first_gne0, cnt
                           gg = pw%pw_grid%g(i, ig)**2 - o3*pw%pw_grid%gsq(ig)
                           pwdr2_gg%array(ig) = gg*pw%array(ig)/pw%pw_grid%gsq(ig)
                        END DO
!$OMP END PARALLEL DO
                     ELSE
!$OMP PARALLEL DO PRIVATE (ig) DEFAULT(NONE) SHARED(cnt, i, j, pw, pwdr2_gg)
                        DO ig = pw%pw_grid%first_gne0, cnt
                           pwdr2_gg%array(ig) = pw%array(ig)*(pw%pw_grid%g(i, ig)*pw%pw_grid%g(j, ig)) &
                                                /pw%pw_grid%gsq(ig)
                        END DO
!$OMP END PARALLEL DO
                     END IF

                     IF (pw%pw_grid%have_g0) pwdr2_gg%array(1) = 0.0_dp

                     CALL timestop(handle)

                  END SUBROUTINE pw_dr2_gg

! **************************************************************************************************
!> \brief Multiplies a G-space function with a smoothing factor of the form
!>      f(|G|) = exp((ecut - G^2)/sigma)/(1+exp((ecut - G^2)/sigma))
!> \param pw ...
!> \param ecut ...
!> \param sigma ...
!> \par History
!>      none
!> \author JGH (09-June-2006)
!> \note
! **************************************************************************************************
                  SUBROUTINE pw_smoothing(pw, ecut, sigma)

                     TYPE(pw_c1d_gs_type), INTENT(INOUT)                       :: pw
                     REAL(KIND=dp), INTENT(IN)                          :: ecut, sigma

                     CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_smoothing'

                     INTEGER                                            :: cnt, handle, ig
                     REAL(KIND=dp)                                      :: arg, f

                     CALL timeset(routineN, handle)

                     cnt = SIZE(pw%array)

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ig, arg, f) SHARED(cnt, ecut, pw, sigma)
                     DO ig = 1, cnt
                        arg = (ecut - pw%pw_grid%gsq(ig))/sigma
                        f = EXP(arg)/(1 + EXP(arg))
                        pw%array(ig) = f*pw%array(ig)
                     END DO
!$OMP END PARALLEL DO

                     CALL timestop(handle)

                  END SUBROUTINE pw_smoothing

! **************************************************************************************************
!> \brief ...
!> \param grida ...
!> \param gridb ...
!> \return ...
! **************************************************************************************************
                  ELEMENTAL FUNCTION pw_compatible(grida, gridb) RESULT(compat)

                     TYPE(pw_grid_type), INTENT(IN)                     :: grida, gridb
                     LOGICAL                                            :: compat

                     compat = .FALSE.

                     IF (grida%id_nr == gridb%id_nr) THEN
                        compat = .TRUE.
                     ELSE IF (grida%reference == gridb%id_nr) THEN
                        compat = .TRUE.
                     ELSE IF (gridb%reference == grida%id_nr) THEN
                        compat = .TRUE.
                     END IF

                  END FUNCTION pw_compatible

! **************************************************************************************************
!> \brief Calculate the structure factor for point r
!> \param sf ...
!> \param r ...
!> \par History
!>      none
!> \author JGH (05-May-2006)
!> \note
! **************************************************************************************************
                  SUBROUTINE pw_structure_factor(sf, r)

                     TYPE(pw_c1d_gs_type), INTENT(INOUT)                       :: sf
                     REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: r

                     CHARACTER(len=*), PARAMETER :: routineN = 'pw_structure_factor'

                     INTEGER                                            :: cnt, handle, ig
                     REAL(KIND=dp)                                      :: arg

                     CALL timeset(routineN, handle)

                     cnt = SIZE(sf%array)

!$OMP PARALLEL DO PRIVATE (ig, arg) DEFAULT(NONE) SHARED(cnt, r, sf)
                     DO ig = 1, cnt
                        arg = DOT_PRODUCT(sf%pw_grid%g(:, ig), r)
                        sf%array(ig) = CMPLX(COS(arg), -SIN(arg), KIND=dp)
                     END DO
!$OMP END PARALLEL DO

                     CALL timestop(handle)

                  END SUBROUTINE pw_structure_factor

                  END MODULE pw_methods
